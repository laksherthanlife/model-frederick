# Source audit — is the citation culture real or decorative?

**Audited 2026-08-26** against `src/ystwin/` (53 modules), `scripts/` (12), and every PMID, PMCID
and DOI in the repository. Identifiers were resolved through NCBI E-utilities and the Crossref
API; claims were checked against full text (PMC, publisher, or the repo's own vendored
BioNumbers export) wherever the text was reachable. Anything I could not reach is marked
**UNVERIFIED** with what I tried.

---

## 1. Verdict

**Both, and the split is clean enough to name.** Where a *document* was written — `docs/research/UPR_ANCHOR.md`, `docs/research/G4_STATISTICS.md` — the sourcing is genuinely excellent: **34 of 34** PMIDs there resolve to exactly the paper the surrounding prose names, the ValidPrime 60 % ceiling is quoted correctly and really is stated three times in the paper, and the doc even flags a "bibliographic trap" and which PMC records are image-only. Where identifiers were *typed inline beside code*, the culture is decorative: in `generator/stress_panel.py` — the repository's single most-imported module — **14 of 31 PMIDs (45 %) resolve to entirely unrelated papers**, including a butterfly-wing paper standing in for the DNA-damage regulon and a hyperspectral dermoscopy paper standing in for the ATP sensor. The named authors and years are right and the biology they assert is right; the numeric identifiers are wrong. The single test guarding this asserts `bool(module.source)` — that a citation *exists*, which a wrong PMID satisfies.

Genuinely sourced: of 97 scalar module-level numeric constants, **24 carry an attached docstring and 12 name a source or derivation**; I verified 11 of those against primary literature and 9 held exactly. But `stress_panel.py` alone contains **136 numeric parameters** (76 stressor target weights, 25 EC50s, 25 lethal doses, 6 crosstalk weights, 4 cascade weights) of which **two** — the DTT and H2O2 lethal doses — are measured, and **none** carries a per-number citation.

---

## 2. Scoreboard

### 2.1 Constant census (mechanically enumerated)

| Population | Count |
| --- | ---: |
| Scalar module-level numeric constants (`src/ystwin` + `scripts`) | **97** |
| …with an attached docstring | 24 |
| …whose docstring names a source or states a derivation | 12 |
| …bare, no docstring, no comment | **73** |
| Numeric function-default arguments | **133** |
| Numeric parameters inside `stress_panel.py` data structures | **136** |
| Distinct PMIDs in the repository | **80** |
| Distinct PMCIDs | 3 |
| Distinct DOIs | 36 |

### 2.2 Classification of the constants I checked in depth (35)

| Category | Count | Examples |
| --- | ---: | --- |
| **Cited and correct** | **13** | `_MAX_CORRECTABLE_GDNA_SHARE=0.60`, `ROGFP2_MIDPOINT_MV=-280`, `CONCENTRATION_FLOOR/CEILING_M`, `thermodynamic.CYTOSOLIC_PH=7.5`, `_RESTING_MAINTENANCE=0.7`, `stress_panel.CYTOSOLIC_PH=(6.4,7.4)`, QUEEN Kd/pH/Mg window |
| **Cited but WRONG** | **5** | 14 broken PMIDs in `stress_panel.py` (counted as one finding per distinct claim below), `PANEL["mOrange2"].maturation_h=2.3`, `PANEL["YFP"].maturation_h=0.3`, `PANEL["CFP"]=0.4`, `_LOG_K_MAGNESIUM` ionic-strength attribution |
| **Cited vaguely** | **4** | `GDCW_PER_OD=0.3` ("a literature scale"), `_GLUCOSE_QMAX=22.0`, `PANEL` maturation block ("order-of-magnitude literature figures"), `literature.LIBRARY` van Hoek entries (`verified=False`) |
| **Asserted, honestly** | **6** | `g1_optical.od_linear_max=1.0`, `MIN_PLATES_FOR_INTERVAL=3`, `DEFAULT_SIMULATIONS=200`, `MAINTENANCE_PER_ACTIVITY=1200` ("a plausible scale, not a fitted one"), `_DEFAULT_BASAL=0.9` (fitted to own plates), `MEASURED_*` panel constants |
| **Asserted but DRESSED UP** | **7** | `DEFAULT_WELL_CV=0.10`, `_MAINTENANCE_PER_ACTIVITY=300.0`, `stress_panel` module docstring's topology claim, `REPORTERS.also_reads` weights, `STRESSORS` EC50s, `_POOL_POTENCY`/`_GENERAL_POTENCY`, `ThermodynamicData` magnesium docstring |

### 2.3 PMID health by location

| Location | PMIDs | Correct | Wrong | Wrong % |
| --- | ---: | ---: | ---: | ---: |
| `docs/research/UPR_ANCHOR.md` + qPCR docs | 34 | 34 | 0 | **0 %** |
| `docs/`, `tests/`, `bridge/`, `fba/` (other) | 14 | 14 | 0 | 0 % |
| `generator/context.py` | 1 | 0 | 1 | 100 % |
| **`generator/stress_panel.py`** | **31** | **17** | **14** | **45 %** |
| `docs/research/G4_STATISTICS.md` (PMCID) | 1 | 0 | 1 | 100 % |
| **Total** | **81** | **65** | **16** | **20 %** |

---

## 3. Cited-WRONG and dressed-up findings

### 3.1 CRITICAL — fourteen PMIDs in the #1 hub resolve to unrelated papers

`generator/stress_panel.py` has fan-in 8 in-package / 40 total (`docs/ARCHITECTURE.md` §2.4 ranks
it first). `docs/research/IDENTIFIABILITY.md` and `docs/FINDINGS.md:181` both justify treating the
loading matrix `L` as *known* on the grounds that it is "transcribed from the literature **with
PMIDs**". That justification is what these identifiers carry.

Every named author/year below is correct and the intended paper does support the biological
claim. Only the numeric identifier is wrong — which is the worst possible failure mode, because
it is invisible to a reader who trusts the string and fatal to anyone who follows it.

| Line | Code claims | PMID cited | What that PMID actually is | Correct PMID | Severity |
| --- | --- | --- | --- | --- | --- |
| 197 `cell_wall` | Jung & Levin 1999 — Rlm1 carries CWI transcription | 10361303 | Tovar 1999, *Mol Microbiol* — "The mitosome… in *Entamoeba histolytica*" | **10594829** | HIGH |
| 200 `dna_damage` | Huang 1998 — Crt1 represses RNR, released on checkpoint | 9663389 | Galant 1998, *Curr Biol* — butterfly achaete-scute homolog / wing scales | **9741624** | HIGH |
| 203 `calcium` | Stathopoulos & Cyert 1997 — Crz1 drives CDRE transcription | 9420334 | Sellers 1998, *Genes Dev* — retinoblastoma protein / E2F | **9407035** | HIGH |
| 206 `carbon` | Young 2003 — Snf1 relieves Mig1, activates Adr1/Cat8 | 12482873 | Despouy 2003, *JBC* — cyclin D3 / retinoic acid receptors | **12676948** | HIGH |
| 209 `nitrogen` | Beck & Hall 1999 — TOR inhibition releases Gln3 | 10604474 | Liang 1999, *Nature* — beclin 1 and autophagy (same issue, off by 4) | **10604478** | HIGH |
| 212 `hypoxia` | Kwast 1998 — heme-dependent Rox1 repression | 9515917 | Maupin-Furlow 1998, *J Bacteriol* — 20S proteasome of *Methanosarcina* | **9510529** | HIGH |
| 215 `copper` | Thiele 1988 — Ace1 activates CUP1 | 3062374 | Robinson 1988, *MCB* — vacuolar protein sorting mutants | **3043194** | HIGH |
| 218 `zinc` | Zhao & Eide 1997 — Zap1 autoregulates, drives ZRT1/2 | 9223284 | Minvielle-Sebastia 1997, *PNAS* — poly(A)-binding protein | **9271382** | HIGH |
| 221 `sulfur` | Thomas & Surdin-Kerjan 1997 — Met4 sulfur regulon | 9409149 | Brown 1997, *MMBR* — "Archaea and the prokaryote-to-eukaryote transition" (**off by one**; the intended paper is the very next article in the same issue) | **9409150** | HIGH |
| 228 `retrograde` | Liao & Butow 1993 — RTG1/RTG3 → CIT2 | 8500171 | Fütterer 1993, *Cell* — ribosome migration on cauliflower mosaic virus RNA | **8422683** | HIGH |
| 231 `alkaline_ph` | Lamb & Mitchell 2003 — Rim101 cleaved/activated at alkaline pH | 12588995 | Kren 2003, *MCB* — War1p and weak-acid stress | **12509465** ⚠ | HIGH |
| 237, 415 `peroxide` / `HyPer7` | Pak 2020 — HyPer7 reports peroxide, pH-insensitive | 32160542 | Laver 2020, *Cell Rep* — RNA-binding protein Rasputin/G3BP | **32130885** | HIGH |
| 240, 420 `atp` / `QUEEN-2m` | Yaginuma 2014 — QUEEN, ratiometric ATP sensor | 24815987 | Vasefi 2014, *Sci Rep* — polarization-sensitive hyperspectral **dermoscopy** | **25283467** | **CRITICAL** (see 3.2) |
| 246, 434 `nadh` / `Peredox` | Hung 2011 — Peredox, cytosolic NADH:NAD | 21982710 | Krahmer 2011, *Cell Metab* — phosphatidylcholine synthesis / lipid droplets (same issue, off by 4) | **21982714** | HIGH |

⚠ `alkaline_ph` is doubly weak: even the *intended* Lamb & Mitchell 2003 (MCB 23:677) is about
Rim101 repressing `NRG1`/`SMP1`, not about its proteolytic activation at alkaline pH. The
proteolysis result is Li et al. 2004 / Xu & Mitchell 2001. Fix the identifier **and** the claim.

**Two further broken identifiers outside `stress_panel.py`:**

| Location | Code claims | Cited | Actually | Severity |
| --- | --- | --- | --- | --- |
| `generator/context.py:35` | Schüller & Entian 1996 — ICL1 moves >200-fold only on transfer to a non-fermentative source | PMID 8710508 | Esposito 1996, *NAR* — "The complete nucleotide sequence of bacteriophage HP1 DNA". I searched PubMed for every Schüller HJ + Entian KD co-authored paper (1330335, 2002006, 3049255, 2823078) — **none is from 1996**, so the intended reference could not be identified either | MEDIUM |
| `docs/research/G4_STATISTICS.md:110` | ValidPrime (Laurell 2012) | PMC3326325 | Sharifulin 2011, *NAR* — ribosomal protein S26 / mRNA binding. Correct ID is **PMC3326333** (PMID 22228834) | MEDIUM |

**Why the test suite does not catch this.** `tests/test_network_wiring.py:54` is
`test_every_module_cites_a_source` → `assert all(module.source for module in MODULES.values())`,
and `tests/test_stress_panel.py:30` is `assert module.source, name`. Both assert that a citation
string is *non-empty*. A PMID pointing at a butterfly paper passes both. This is the mechanism by
which a citation culture becomes decorative: the invariant enforced is presence, not truth.

### 3.2 CRITICAL — the headline sensor build rests on a claim its source does not make

`generator/stress_panel.py:417–425`:

```python
"QUEEN-2m": Reporter(..., source="Yaginuma 2014, PMID 24815987 -- ratiometric ATP sensor,
                     shown in yeast; Kd 4.5 mM ATP at 25 C, flat against Mg2+ over 1-2 mM,
                     and flat against pH over 7.3-8.8, ...")
```

I read Yaginuma et al. 2014, *Sci Rep* **4**:6522 (PMC4185378), the real QUEEN paper. Findings:

- **Kd 4.5 mM at 25 °C for QUEEN-2m** — VERIFIED verbatim.
- **"nearly unaffected at pH 7.3–8.8 and also by physiological Mg2+ concentrations in the 1–2 mM range"** — VERIFIED verbatim. The transcription of these three numbers is exact.
- **"shown in yeast" — FALSE.** The paper's title is *"Diversity in ATP concentrations in a single **bacterial** cell population"*. QUEEN was developed and demonstrated in *E. coli*. There is no *S. cerevisiae* work in it. QUEEN *was* later ported to budding and fission yeast — by **Takaine et al. 2019, *J Cell Sci*, doi:10.1242/jcs.230649** — which the repository does not cite anywhere.

This is load-bearing, not cosmetic:

- `Reporter.demonstrated_in_yeast` defaults to `True` and is overridden to `False` only for
  Peredox, so QUEEN-2m inherits `True` on the strength of this source string.
- `demonstrated_in_yeast()` is documented as "a build restricted to these is one the wet lab can
  start on; the rest would each need a sensor ported first, which is a different project."
- `RECOMMENDED_BUILD` resolves to `['STRE-general', 'UPRE-ER', 'TRX2-oxidative', 'QUEEN-2m']`
  and `THREE_SENSOR_BUILD` to `['UPRE-ER', 'TRX2-oxidative', 'QUEEN-2m']`, whose docstring says
  *"all three demonstrated in S. cerevisiae."*
- `tests/test_sensor_availability.py:73` asserts
  `all(REPORTERS[r].demonstrated_in_yeast for r in RECOMMENDED_BUILD.stress_reporters)`.
- `analysis/sensor_selection.py:181` records that naming the axes "swaps HyPer7 for QUEEN… lifts
  ATP to 0.67 and held-out transfer from five modules to seven."

So the repository's flagship design recommendation — and the specific 5→7-module transfer gain
attributed to including an ATP channel — is gated on a boolean whose only stated evidence is a
paper about *E. coli*, carrying a PMID for a paper about skin imaging. The *conclusion* survives
(Takaine 2019 does establish QUEEN in yeast), but nothing in the repository supports it.

**HyPer7 has the same structural hole, unflagged.** `REPORTERS["HyPer7"]` also inherits
`demonstrated_in_yeast=True`, and its source string (Pak 2020) makes no yeast claim at all — Pak
et al. characterised HyPer7 in HeLa/zebrafish. Published HyPer7-in-yeast work I could find is in
*Komagataella phaffii*, not *S. cerevisiae*. **UNVERIFIED for *S. cerevisiae*.**

### 3.3 HIGH — `analysis/power.py::DEFAULT_WELL_CV = 0.10`: the honesty statement points the wrong way

The docstring (`power.py:129–146`) says, in full:

> "The value is not measured for these constructs… **published well-to-well spread for yeast plate
> reporter assays is generally larger, not smaller, so this is more likely optimistic than
> pessimistic**… Until then, **replicate counts from this module are lower bounds**."

**The directional claim gestures at "published spread" and names nothing.** It is unresolvable as
written, so it is a *dressed-up* assertion, not a cited one. Worse, the evidence I could find
points the other way:

- Kensy et al. 2009, *Microb Cell Fact* (PMC2700080) — the standard validation of online
  biomass + fluorescence monitoring in shaken microtiter plates, i.e. exactly this assay class —
  reports **"data reproducibility of smaller than 5 % standard deviation, when cultivating the
  same clone in the same medium on a microtiter plate."** That is *half* of 0.10.
- Plate-assay QC convention treats **CV < 10 % as the acceptance ceiling** for a high-quality
  assay, not as a typical biological floor.

If the true within-plate well CV is nearer 0.05, then 0.10 is *pessimistic*, and the replicate
counts this module returns are **upper bounds, not lower bounds** — the opposite of what the
docstring tells a reader to conclude. `replicates_needed()` is a number somebody then goes and
grows cultures for, and `docs/ARCHITECTURE.md` §5.1 calls this module's interval "the
best-behaved uncertainty in the repo". A wrong safety direction on the best-behaved number is
worth more than a wrong constant.

The neighbouring `biological_cv: float = 0.12` (between-plate, `power.py:199/370/445`) is a bare
default with **no docstring at all** and no source, though the repo *has* a measured between-plate
figure — `generator/panel_experiment.OBSERVED_ACTIVITY_CV = 0.146`, from 28 matched conditions on
two real plates. Two plate-level spread constants, one measured and one invented, live in the
same package and are not reconciled.

*(Nothing here disputes the variance-decomposition argument itself — separating plate-level from
well-level noise is correct and the docstring's reasoning about why the plate term cancels out of
a within-plate dose slope is sound.)*

### 3.4 HIGH — the fluorophore maturation panel disagrees with the repo's own vendored data

`generator/unmixing.py:45–52`, guarded only by the comment *"Maturation half-times are
order-of-magnitude literature figures, not measured here."* No source is named for any entry.

| Fluorophore | Code | Published | Source |
| --- | ---: | --- | --- |
| **mOrange2** | **2.3 h** | **4.5 h** | Shaner et al. 2008, *Nat Methods* — mOrange2 maturation t½ **4.5 h**, up from mOrange's 2.5 h; the slowdown was the acknowledged cost of the 25× photostability gain |
| **YFP** | **0.3 h** | **39 min = 0.65 h** | **BNID 102974 / 108677, organism = *Saccharomyces cerevisiae*** — in the repo's own `data/kaggle/BioNumbers_Nov2024.csv` |
| **CFP** | **0.4 h** | **49 min = 0.82 h** | **BNID 106883 / 108678, organism = *S. cerevisiae*** — same file |
| eGFP | 0.4 h | 14 min–4 h across BNID 107000/107002/106884 | within a very wide published spread |
| mApple | 0.6 h | ~0.5 h (Shaner 2008) | plausible |

Two things make this worse than a rounding error:

1. **The correct values for YFP and CFP, measured *in this organism*, are already in the
   repository** — in the BioNumbers export that `tests/test_constants_against_bionumbers.py`
   already opens to check four *other* constants. The test file's own preamble says BioNumbers
   "was assembled by people with no interest in this model, which is what makes it worth checking
   against." The fluorophore panel was never put through that check.
2. **2.3 h is not mOrange2's number and is not mOrange's either.** It sits between them, which is
   the signature of a half-remembered figure rather than a transcription slip.

**Blast radius.** `archive/tests/test_deconvolve.py:14` states the design conclusion: *"mOrange2 is
therefore not a liability. It is the best mu sensor in the panel,"* and `RECOMMENDED_BUILD`
carries `reference_fluorophore="mOrange2"`. The conclusion is *robust*: the load-bearing quantity
is the ratio (7.7× as coded, 6.9× using the published values), and the true mOrange2/YFP contrast
is if anything slightly *less* extreme but in the same direction. So the recommendation survives —
but any quantitative statement about maturation-mismatch magnitude (`_steady_state`,
`ratio_series`, the spectral-feasibility work) is running on numbers that are 2× off in two
channels and ~2× off in a third.

### 3.5 HIGH — `_MAINTENANCE_PER_ACTIVITY`: a 4× library/script split, neither grounded

| Location | Value | Provenance |
| --- | ---: | --- |
| `bridge/latent_bridge.py:26` | `300.0` | inline comment `# mmol ATP/gDW/h per unit promoter activity`. **No docstring, no derivation, no source.** |
| `scripts/run_scenarios.py:35` | `1200.0` | `# ATP cost per unit activity: a plausible scale, not a fitted one.` |

`run_scenarios.py` is the **only** caller and passes `1200.0` explicitly, so the library default
of `300.0` is dead — but it is a public keyword argument (`maintenance_per_activity=`), so any new
caller silently gets a quarter of the effect the one working pipeline uses. Nothing in the
repository explains where either number came from, and neither is fitted.

The magnitude is not negligible. `atp_maintenance = 0.7 + k × excess` on a Yeast9 NGAM of 0.7
(verified: `r_4046` bounds are exactly `(0.7, 0.7)` in the actual model file). At an excess of
~1e-3 activity units, `k = 1200` raises maintenance by ~170 % and `k = 300` by ~43 %. The
`NGAM B` column in `outputs/scenario_predictions.csv` is that number.

**Mitigating, and it is substantial:** `latent_constraints()` *refuses to run* without
`acknowledge_unvalidated=True`, and stamps every result with the G4 refusal string. The branch is
labelled research-only in the dataclass, in the printed report, and in `compare_branches`'s
`load_bearing=False`. So the number is arbitrary but it is quarantined. The defect is the silent
4× fork, not the use.

Related, and unsourced: `latent_constraints(basal_activity: float = 1.0e-3)` — the activity
treated as unstressed, which sets the zero point of the whole effect. No docstring, no source.

### 3.6 MEDIUM — `stress_panel.py`'s opening docstring contradicts the file and the tests

Lines 3–6, the first thing any reader of the #1 hub sees:

> "…Msn2/4 sits under most of them, **Rpn4 answers to Msn2/4, Yap1 and Hsf1** together, and Hog1
> routes its transcriptional output through Msn2."

But line 187–190 of the same file:

```python
"proteasome": Module("Rpn4", "PACE", "... RPN4 has no STRE",
                     driven_by={"oxidative": 0.45, "heat": 0.40}),
```

There is **no** Msn2/4 → Rpn4 edge, and the module's own source string says `RPN4 has no STRE`.
`tests/test_network_wiring.py` locks in the removals (`test_the_general_stress_module_does_not_
drive_the_upr`, `…_the_heat_module`, `…_the_oxidative_module`, `test_the_osmotic_arrow_is_not_
drawn_from_the_general_stress_module`). The data structure and the test suite are right — the
literature has RPN4 under Yap1/Hsf1/Pdr1-3 with no canonical STRE, and Hahn 2006 confirmed the
YRE and HSE by promoter mutagenesis (**VERIFIED**). The module docstring is the stale artefact of
the topology that was corrected, and it is the summary a reader will quote.

`docs/ARCHITECTURE.md` §6.5 already keeps a list of docstrings that disagree with behaviour. This
one is not on it and belongs there.

### 3.7 MEDIUM — `bridge/thermodynamic.py`: the docstring says magnesium is not modelled; the code models it

`ThermodynamicData`'s class docstring (line ~127):

> "**Magnesium binding is NOT modelled here**, and it matters most for exactly the adenylates."

Lines 187–190 of the same file:

```python
log_k = _LOG_K_MAGNESIUM.get(key)
if log_k is not None and free_magnesium_m > 0:
    bound = 10.0**log_k * free_magnesium_m
    value -= GAS_CONSTANT_KJ * STANDARD_TEMPERATURE_K * np.log1p(bound)
```

Magnesium binding *is* modelled, for exactly the three adenylate/phosphate species the docstring
says it matters most for, and `free_magnesium_m` defaults to 1 mM. Anyone reading the docstring
would double-correct.

Separately, `_LOG_K_MAGNESIUM = {"cpd00002": 4.19, "cpd00008": 3.17, "cpd00009": 1.88}` is
documented as *"from Alberty's tabulations **at I = 0.25 M**."* **UNVERIFIED** — I could not reach
Alberty's tables (the 2003 monograph is not open access; eQuilibrator's FAQ host presented an
invalid certificate; a web search returned only a paraphrase of my own query). The attribution
should be checked against a primary table before it is relied on: Mg association constants fall
substantially with ionic strength, so *which* ionic strength these three values are quoted at
changes the transformed ΔG_f of every adenylate. Flag, do not yet assert.

### 3.8 MEDIUM — the loading matrix is not "transcribed from the literature with PMIDs"

`docs/research/IDENTIFIABILITY.md:17` and `docs/FINDINGS.md:181` both defend treating `L` as a
known operator — the entire Path A identifiability result, and the README's "2 → 2 of 7, 4 → 4 of
7, 7 → 7 of 7" table — on the grounds that `L` is *"transcribed from the literature with PMIDs."*

Measured against `reporter_loadings()`, which is what actually builds `L`:

| | With a resolvable reference | Total |
| --- | ---: | ---: |
| `MODULES` (module identity + regulon) | 24 | 24 |
| `REPORTERS` (which promoter reads which module) | **5** | 24 |
| `REPORTERS.also_reads` crosstalk weights — **the off-diagonal of `L`** | **0** | 6 |
| `STRESSORS` | **7** | 25 |
| `STRESSORS.targets` weights | **0** | 76 |

The PMIDs are attached to the *module list* — which regulons exist. `L`'s entries are the
`also_reads` weights, and every one of those six numbers (0.25, 0.30, 0.35, 0.40, 0.25, 0.30) is
a bare literal with a mechanism sentence and no reference: *"HSP104 carries both HSE and STRE"* is
true and gives no basis for 0.35 rather than 0.2 or 0.5. Same for the 76 stressor target weights
and the 25 EC50s that drive `module_ec50`, `healthy_ladder` and every simulated panel.

Two of the 136 numeric parameters in the file are measured (`DTT.lethal_dose=1.55` and
`H2O2.lethal_dose=1.0`, both from the repo's own plates, both documented at length and honestly —
including that the EC50s beside them are *not* tested). The provenance sentence in the docs
should be narrowed to what is true: the *topology* is literature-transcribed; the *weights* are
expert priors.

### 3.9 LOW — duplicated thresholds that can drift apart

`viz/figures.py:1003` re-declares `CORRECTABLE_GDNA_CEILING = 0.60` and `:1011`
`PASS_MARGIN_CYCLES = 3.0` as fresh literals instead of importing `qpcr._MAX_CORRECTABLE_GDNA_
SHARE` and `qpcr._MIN_RT_MINUS_MARGIN`. The figure and the gate can silently disagree. Same
pattern for `run_power.py::DOSES`/`FOLDS` vs `viz/figures.py::DESIGN_DOSES`/`POWER_FOLDS`.

### 3.10 LOW — the one constant with a real provenance note has it detached from the constant

`stress_panel.py:37` is `_DEFAULT_BASAL = 0.9`, immediately followed by `CYTOSOLIC_PH` and *its*
docstring; the string at lines 48–54 that explains `_DEFAULT_BASAL`'s calibration is an orphan
expression statement eleven lines below the assignment, attached to nothing. No documentation
tool will associate them. (`docs/ARCHITECTURE.md` §6.5 item 5 already records this — credit
where due.)

I checked the arithmetic the orphan claims: with `basal = 0.9`, `UPRE-ER` across
`healthy_ladder("DTT", 5)` gives folds **1.056, 1.107, 1.194, 1.319, 1.418** — matching the stated
"measured induction across the whole DTT ladder is 1.38 to 1.47 fold" at the top rung.
**Constant VERIFIED against the repo's own plates.** It is honestly derived; only its
documentation is orphaned.

---

## 4. PMID / DOI / PMCID resolution table

### 4.1 `generator/stress_panel.py` — the 17 that are correct

| PMID | Resolves to | Claim attached | Verdict |
| --- | --- | --- | --- |
| 8641288 | Martínez-Pastor 1996, *EMBO J* 15:2227 — Msn2p/Msn4p required for STRE induction | ESR = Msn2/4 + STRE | **Correct** |
| 8898193 | Cox JS 1996, *Cell* 87:391 — regulating a TF that controls the UPR | Cox & Walter; Ire1 splices HAC1, so UPR is post-transcriptional | **Correct** |
| 10411744 | Boy-Marcotte 1999, *Mol Microbiol* 33:274 | (a) "Hsf1 and Msn2/4 regulons overlap only slightly"; (b) "25→38 °C, half of induced proteins Msn2/4-dependent" | **Correct — verified verbatim.** Paper: "Among 52 proteins induced by a shift from 25 °C to 38 °C, half… dependent upon Msn2p and/or Msn4p… The two sets of proteins overlapped only slightly." |
| 37467033 | Ciccarelli 2023, *Mol Biol Cell* 34:ar101 | "coupling is compensatory" | **Correct — verified.** "Hsf1 and Msn2/4 display functional compensatory induction following HS." |
| 9472026 | Görner 1998, *Genes Dev* 12:586 — Msn2p nuclear localisation | "nuclear entry survives hog1 deletion" | **Correct paper**; specific hog1Δ result not re-read (UNVERIFIED detail, consistent with the paper's scope) |
| 11918814 | Owsianik 2002, *Mol Microbiol* 43:1295 — control of 26S proteasome expression by MDR TFs | RPN4 promoter elements | **Correct** |
| 16556235 | Hahn 2006, *Mol Microbiol* 60:240 — stress regulatory network for proteasome expression via Hsf1 | "YRE and HSE in the RPN4 promoter, both confirmed by mutagenesis" | **Correct — verified.** (The separate MMS-0.02 % figure at line 316: **UNVERIFIED**, paywalled.) |
| 8670839 | Yamaguchi-Iwai 1996, *EMBO J* 15:3377 — iron-regulated DNA binding by Aft1 | iron = Aft1/Aft2 + FeRE | **Correct** |
| 18627600 | Salin 2008, *BMC Genomics* 9:333 — selenite stress networks | (a) "a weak YRE upstream of AFT2"; (b) "Rpn4 and Pdr1 form a positive loop" | **Correct — both verified.** (a) "one YRE (TTAGTCA) conserved in *Saccharomyces sensu stricto*… 147 base pairs upstream from the ATG"; (b) "a positive transcriptional loop connects *RPN4* and *PDR1*". Note the code's hedge "single-source and never replicated" *understates* it — the paper also confirms Yap1 binding at the AFT2 promoter by ChIP-chip |
| 18469822 | Gutscher 2008, *Nat Methods* 5:553 — real-time imaging of glutathione redox potential | roGFP2-Grx1 reports E_GSH as an equilibrium ratio | **Correct** |
| 22705944 | Kojer 2012, *EMBO J* 31:3169 | "cytosolic E_GSH at −306 mV and matrix at −301 mV" | **Correct — verified exactly.** "an EGSH[cytosol] of approximately −306±1.3 mV"; matrix and IMS "−301±2 and −301±5". Probe: Grx1-roGFP2 |
| 23762325 | Ayer 2013, *PLoS ONE* 8:e65240 | "−350 mV with a measured rather than assumed cytosolic pH" (`stress_panel`); "Ayer 2013 reports 7.5 in the yeast cytosol using pHluorin" (`thermodynamic.CYTOSOLIC_PH`) | **Correct — both verified.** "E_GSH between −340 to −350 mV"; "pH 7.5±0.2 in the cytosol"; "The pH value obtained using pHluorin was used for calculation of E_GSH." Minor: Ayer used roGFP2, not Grx1-roGFP2 |
| 9671304 | Miesenböck 1998, *Nature* 394:192 — pH-sensitive GFPs | ratiometric pHluorin reports cytosolic pH | **Correct** |
| 20581803 | Dechant 2010, *EMBO J* 29:2515 — cytosolic pH as a second messenger for glucose | `CYTOSOLIC_PH = (6.4, 7.4)`; prose "Resting is about 7.2, and glucose removal drops it to 6.4" | **Correct** — glucose starvation drops cytosolic pH to **~6.4** from **~7.4**, attributed to Pma1 inhibition, all verified. **Minor inconsistency:** the prose says resting is 7.2 while the tuple it documents uses 7.4; the literature supports the tuple, so fix the prose |
| 9393686 | Larsson 1997, *J Bacteriol* 179:7243 — glycolytic flux conditionally correlated with ATP | "glycolytic flux correlates NEGATIVELY with intracellular ATP and not at all with the ATP/ADP ratio" | **Correct** |
| 25955212 | Zhao 2015, *Cell Metab* 21:777 — SoNar NAD+/NADH sensor | SoNar never applied in yeast | **Correct paper**; the negative ("never applied in yeast") is unfalsifiable from one PMID but consistent |
| 25401080 | Knudsen 2014, *AMB Express* 4:81 — NADH-dependent biosensor in *S. cerevisiae* | "GPD2-promoter fusion… a transcriptional proxy and not a ratio" | **Correct** |

### 4.2 `generator/stress_panel.py` — the 14 that are wrong

See §3.1. Every row is a HIGH-or-above severity mismatch.

### 4.3 Elsewhere in the repository — all correct

| Location | PMID / ID | Resolves to | Verdict |
| --- | --- | --- | --- |
| `bridge/tmfa.py:8`, `tests/test_tmfa.py:8` | 17172310 | Henry 2007, *Biophys J* 92:1792 — TMFA | **Correct.** `CONCENTRATION_FLOOR_M = 1e-5` / `CEILING = 2e-2` **verified verbatim**: "The concentrations of all intracellular species were restricted within the ranges observed in the cell (between 10−5 M and 0.02 M)" |
| `bridge/tmfa.py:182` | 18383140 | Canelas 2008, *Biotechnol Bioeng* 100:734 — cytosolic free NAD/NADH ratio | **Correct** |
| `bridge/thermodynamic.py:10` | 16788595 | Kümmel 2006, *Mol Syst Biol* — network-embedded thermodynamic analysis | **Correct** |
| `bridge/thermodynamic.py:11` | 25028891 | Martínez 2014, *Biophys J* 107:493 — thermodynamic curation of human and yeast GEMs | **Correct** |
| `tests/test_thermodynamic_bridge.py` | 17287356 | Vemuri 2007, *PNAS* 104:2402 — increasing NADH oxidation reduces overflow | **Correct** |
| `tests/test_thermodynamic_bridge.py` | 21335394 | Agrimi 2011, *AEM* 77:2239 — mitochondrial NAD carriers | **Correct** |
| `docs/FINDINGS.md:1392` | 26098102 | Su 2015, *PLoS ONE* 10:e0130840 — redox imbalance / 7-DHC | **Correct** |
| `tests/test_physiology_validation.py:27` | 2566299 | Postma 1989, *AEM* 55:468 — Crabtree effect, CBS 8066 | **Correct** |
| `tests/test_physiology_validation.py:49` | 27317316 | Vos 2016, *Microb Cell Fact* 15:111 — maintenance energy at near-zero growth | **Correct** |
| `tests/test_atp_sensor.py:30`, `docs/ATP_SENSOR.md:51` | 9427394 | Kratzer & Schüller 1997, *Mol Microbiol* 26:631 — ACS1 control by CAT8/ADR1/UME6 | **Correct** |
| `docs/research/G4_STATISTICS.md` | 19246619 | Bustin 2009, *Clin Chem* 55:611 — MIQE | **Correct paper.** The attached claim ("mandates no numeric margin") is **PARTIALLY VERIFIED** — consistent with the MIQE checklist as I know it, but the full text is paywalled and I could not confirm verbatim |
| `docs/research/G4_STATISTICS.md` | 40272429 | Bustin 2025, *Clin Chem* 71:634 — MIQE 2.0 | **Correct** |
| `docs/research/G4_STATISTICS.md` | 16899212 | Nordgård 2006, *Anal Biochem* 356:182 — error propagation in relative RT-PCR | **Correct** |
| `docs/research/UPR_ANCHOR.md` | **34 PMIDs** — 8898194 Sidrauski 1996; 9323131 Sidrauski 1997; 9348528 Kawahara 1997; 9430730 Kawahara 1998; 9382810 Chapman 1997; 11595189 Rüegsegger 2001; 27692069 Di Santo 2016; 15314654 Leber 2004; 35163590 Hata 2022; 36284099 Matabishi-Bibi 2022; 39995043 Matabishi-Bibi 2025; 10847680 Travers 2000; 16371132 Kimata 2006; 8423809 Kohno 1993; 40959222 Geronimo 2025; 40772766 Ishiwata-Kimata 2025; 29165698 Sarkar 2018; 20625545 Pincus 2010; 21444684 Rubio 2011; 19874630 Teste 2009; 21452013 Vaudano 2011; 37454173 Đermić 2023; 10948426 Wiame 2000; 1503777 Bickler 1992; 8313910 Kuge 1994; 10347154 Lee 1999; 9712873 Godon 1998; 15135069 / 15225652 Tsuzi 2004; 11102521 Gasch 2000; 17959824 Brauer 2008; 12184808 Vandesompele 2002; 22228834 Laurell 2012 | **34 / 34 correct.** Every identifier matches the paper the surrounding prose names |

### 4.4 PMCIDs and DOIs

| ID | Resolves to | Verdict |
| --- | --- | --- |
| **PMC3326325** (`G4_STATISTICS.md:110`) | Sharifulin 2011, *NAR* 40:3056 — ribosomal protein S26 | **WRONG.** Should be PMC3326333 |
| PMC7646510 (`literature.py`, `stress_panel.py`) | Phosphoproteome response to DTT | Not independently re-read — **UNVERIFIED**, plausible |
| PMC3847955 (`literature.py`) | Dilution / growth-rate-dependent gene expression | **UNVERIFIED**, plausible |
| **All 36 DOIs** (`G4_STATISTICS.md`, `docs/research/*`, `docs/FINDINGS.md`) | Resolved via the Crossref API; every one returns the author, year, journal and title the citing text states — Laurell 2012, Bustin 2009/2010/2025, dMIQE 2020, Li 2022, Đermić 2023, Padhi 2016, Hellemans 2007, Nordgård 2006, Bilgrau 2016, Ruijter 2009/2021, Walther & Schüller 2001, Cranmer 2020, Argelaguet 2018/2020, Nolan 2006, Kapteyn 2021, Modrák 2025, Frazier 2020, Säilynoja 2022, Eriksson 2004, Scheuerer 2015, Gneiting 2007 ×2, Matheson & Winkler 1976, Dawid 1984, Schölkopf 2021, Muratore 2021/2022, Ledermann 1937, Jakobi 1995, Goyal 2014, Koos 2013 | **36 / 36 correct** |

**Verified against the ValidPrime full text (PMC3326333), the single most load-bearing external
number in the qPCR path:** the docstring at `qpcr.py:77–88` claims the paper "states three times
that correction holds only 'as long as the DNA contribution to the total signal is <60 %', above
which its own software refuses." All three statements are present —

> "the correction was less precise when the gDNA background exceeded 60 % of the total signal";
> "efficient correction… is possible, as long as the DNA contribution to the total signal is <60 %";
> "it is not advisable to perform correction on samples where the DNA-derived signal exceeds 60 %"

— and the software does grade such samples `HIGHDNA` and decline to correct. **Cited and correct,
including the count.** The derived arithmetic is correct too: 0.60 share ↔ 0.737 cycles ("about
0.74"), 3 cycles ↔ 12.5 %, 5 cycles ↔ 3.1 %, 0.6 cycles ↔ 66 %, 1.0 ↔ 50 %; and
`G4_STATISTICS.md`'s variance-inflation column `√(1+r²)/(1−r)` reproduces to 15.75 / 2.92 / 2.24 /
1.15 against the table's 15.72 / 2.90 / 2.24 / 1.15.

---

## 5. Method-correctness findings

Each implementation was read against its canonical reference. Six are correct as prescribed; four
carry a caveat.

### Correct as prescribed

| Method | Location | Canonical reference | Finding |
| --- | --- | --- | --- |
| **Systematic resampling** | `estimator.py::_resample` | Kitagawa 1996; Doucet & Johansen 2009 | **Correct.** `positions = (rng.random() + arange(N))/N` gives u ~ U[0, 1/N) with the required deterministic stride; `searchsorted(cumsum(w), positions)` is the correct inverse-CDF map (side='left' returns the first index with cumsum ≥ u); weights are reset to uniform after. Lower variance than multinomial, exactly as the docstring says |
| **Phase randomisation** | `analysis/nulls.py::_phase_randomise_1d` | Theiler et al. 1992, *Physica D* 58:77 | **Correct FT surrogate.** Magnitudes preserved exactly; phases i.i.d. U[0,2π); DC phase zeroed; Nyquist phase zeroed for even *n*; `irfft` enforces Hermitian symmetry so the output is real by construction; mean restored. This is the linear-Gaussian-process surrogate the docstring claims, not the amplitude-adjusted (AAFT) variant — and it does not claim to be. The docstring's "approximately the same autocorrelation" is conservative: the periodogram, hence the circular autocorrelation, is preserved exactly |
| **BIC** | `analysis/latent.py:123` | Schwarz 1978 | **Correct.** `n·ln(σ̂²) + k·ln(n)` is `−2·loglik + k·ln(n)` for a Gaussian up to an additive constant in *n*, which cancels in model comparison at fixed *n*. `select_dimension` correctly declines to use it (the fit is exact at *k* = channels) and uses held-out error instead — the right call, and stated |
| **Fisher information** | `analysis/design.py::fisher_information` | Standard linear-Gaussian result | **Correct.** `Lᵀ S⁻¹ L` with per-channel variance scaled to `‖row‖`; the pseudo-inverse diagonal for standard errors; and — a good catch — an explicit null-space leakage test so a module inside the null space reports `inf` rather than a finite-looking diagonal. No canonical reference is named, but the derivation is written out in the module docstring, so this is *derived*, not *asserted* |
| **NIS / χ² acceptance band** | `scripts/run_calibration_nis.py::_acceptance_band` | Bar-Shalom, Li & Kirubarajan 2001, §5.4 | **Correct.** Scalar innovations → NIS ~ χ²(1); the *K*-step average ~ χ²(K)/K; band = `[χ²_{α/2}(K)/K, χ²_{1−α/2}(K)/K]`. The docstring's warning that quoting the single-sample band [0.001, 5.02] for a 25-step average "would make almost anything pass" is exactly right. The whiteness test (lag-1 autocorrelation against ±2/√K) is the standard companion and the script is correct that NIS alone can be bought back by retuning σ |
| **ddCq / Livak** | `qpcr.py::delta_delta_cq` | Livak & Schmittgen 2001; Pfaffl 2001 | **Correct**, and better than Livak: `fold = (1+E)^(−ΔΔCq)` generalises to E < 1 as documented. It refuses rather than guesses on a missing reference gene, a missing per-sample reference reading, or a missing control dose — three failure modes that normally produce a silent wrong answer. Technical replicates are averaged on the Cq (log) scale, which is the standard convention |

### Carrying a caveat

| Method | Location | Finding | Severity |
| --- | --- | --- | --- |
| **TOST / equivalence** | `analysis/uncertainty.py::equivalence` | The **interval level does not match the canonical TOST**. Schuirmann (1987) equivalence at level α uses the **100(1−2α) %** interval — a 90 % CI for α = 0.05. `equivalence()` consumes `FoldChange`, whose `confidence` defaults to **0.95**, so requiring the whole 95 % CI inside the margin is a TOST at **α = 0.025 per side**, not 0.05. The direction is conservative (it will declare EQUIVALENT less often than nominal), so no false claim can come out of it — but the docstring says "Two one-sided tests, in interval form" without saying which α that interval delivers, and a reader will assume 0.05. Separately, the DIFFERENT arm tests the interval against **1.0**, which is a plain difference test, not part of TOST. The *design* is right — pre-naming the margin is the load-bearing part and the docstring is excellent on why | MEDIUM |
| **Cluster bootstrap** | `analysis/uncertainty.py::fold_change` | The two-stage resample (clusters with replacement, then units within each drawn cluster) is a recognised construction for clustered data (Davison & Hinkley 1997 §3.8; Field & Welsh 2007, *JRSS-B*), and the per-plate-then-combine ordering with a geometric mean on the log ratio is right and well argued. **The caveat is the t-inflation**: `lo/hi` are rescaled about the centre by `t_{1−α/2}(G−1) / z_{1−α/2}` and the docstring calls this *"the standard small-sample correction."* It is not. Cameron, Gelbach & Miller (2008) recommend T(G−1) critical values for the **cluster-robust t-statistic** and the wild cluster bootstrap-*t*; there is no standard correction of this form for a **percentile** bootstrap interval. It is a reasonable ad-hoc widening — the docstring's "and is honest about being approximate" is fair — but "standard" dresses it up. Call it what it is | MEDIUM |
| **Minimum cluster count** | `analysis/uncertainty.py::MIN_PLATES_FOR_INTERVAL = 3` | **Honestly asserted, and the combinatorics check out.** For G = 3, distinct resample multisets = C(5,3) = **10** ✓ (docstring says 10); for G = 2, C(3,2) = **3**, of which two are degenerate ✓. Against the literature, 3 is far below anything anyone recommends — Cameron & Miller (2015) document severe over-rejection of cluster-robust inference below G ≈ 30–50, and Bell & McCaffrey (2002) exists precisely because small-G behaviour is bad. The docstring says exactly this ("Three is already poor… not as a count anyone should design to"), so this is a model of an honestly-asserted constant, not a defect. Refusing at G = 2 with a reason string rather than emitting a tight-looking number is the right behaviour | — |
| **Savitzky-Golay** | `growth.py`, `reporter.py`, `analysis/uncertainty.py` | Savitzky & Golay 1964 assumes **uniform sampling**, and `growth.py::_regular_grid` correctly resamples onto a uniform grid when the clock drifts and interpolates the derivative back — a subtlety most implementations miss. Window/polyorder guards are correct (odd length, `window > polyorder`, capped at series length). **Caveat:** `reporter.py::promoter_activity` uses `polyorder=3` for the derivative while `growth.py` uses `polyorder=2`, and neither default is justified anywhere; `default_activity_window_h` picks `duration/6` clipped to `[4·dt, duration/3]` with no stated basis. `tests/test_reporter_noise.py` does document the real cost (double differentiation through a maturation step amplifies noise), which is the right thing to have measured | LOW |
| **Wilson score interval** | `analysis/power.py::PowerEstimate` | **Correct.** Centre `(p + z²/2n)/(1 + z²/n)` and half-width `(z/(1+z²/n))·√(p(1−p)/n + z²/4n²)` is Wilson (1927) exactly, clipped to [0,1]. The docstring's reason for choosing it over Wald — "these estimates live at 0.0 and 1.0, where a Wald interval has zero width and covers nothing" — is precisely why Wilson is the right choice. `replicates_needed(require_confidence=True)` requiring the *lower* bound to clear the target, with the stated reason (a walk that stops at the first point estimate over the line stops on the first *lucky* draw), is a genuinely sophisticated correction | — |
| **β-carotene pathway** | `fba/carotenoid.py` | **Verified empirically against the real Yeast9 model** (`yeast-GEM.xml`, v9.0.2). All four reactions return `check_mass_balance() == {}`. I initially suspected the `proton: −4` in CRTI was spurious; it is not — Yeast9's FAD is `C27H30N9O15P2` (−3) and FADH2 is `C27H33N9O15P2` (−2), so FADH2 = FAD + 3 H, +1 charge, and four reductions genuinely consume 4 H⁺. Gene assignments are correct for *Xanthophyllomyces dendrorhous*: crtE = GGPP synthase, crtI = phytoene desaturase, and the bifunctional crtYB split across a phytoene-synthase (crtB) and lycopene-cyclase (crtY) reaction sharing one gene rule. Stoichiometries are right (2 GGPP → phytoene + 2 PPi; four desaturations to lycopene). **Caveat:** the module names no source for any of it, and its balance test is `pytest.mark.integration` and needs a GSMM the repo does not ship — so on a machine without `YSTWIN_YEAST_GEM` the whole check silently skips | — |

---

## 6. Unsourced but load-bearing, ranked by blast radius

Ranked using `docs/ARCHITECTURE.md` §2.4's fan-in table.

| # | Constant(s) | Where | Fan-in | What it drives | Backing |
| --: | --- | --- | ---: | --- | --- |
| **1** | **76 stressor target weights + 25 EC50s + 6 `also_reads` crosstalk weights** | `generator/stress_panel.py` | **8 / 40** | Every simulated panel; `reporter_loadings` → `fisher_information` → `module_standard_errors` → `RECOMMENDED_BUILD` → the README's identifiability table and the 5→7-module transfer claim | **None.** A mechanism sentence per stressor; no number carries a reference. Documented in the repo *as if* sourced ("transcribed from the literature with PMIDs") — see §3.8 |
| **2** | `_POOL_POTENCY = 0.35`, `_GENERAL_POTENCY = 2.5`, `_LETHAL_MULTIPLE = 3.0`, `_LETHAL_HILL = 2.5`, `_MU_MAX = 0.40` | `stress_panel.py:250–255` | **8 / 40** | `module_ec50` — the entire separation between two regulons driven by one agent, which is what makes any dose ladder identifiable at all. Also `viability`, `growth_rate`, `healthy_ladder` | **None.** No docstrings. The `Stressor` class docstring rationalises the *ordering* ("pools lead… the general stress response trails") but nothing supports 0.35 vs 0.5, or 2.5 vs 2.0. **Dressed up** |
| **3** | `DEFAULT_WELL_CV = 0.10`, `biological_cv = 0.12` | `analysis/power.py` | 3 / 12 | `replicates_needed()` — the number of cultures a person grows. §5.1 of ARCHITECTURE calls this "the best-behaved uncertainty in the repo" | 0.10 is dressed up with an unsupported and probably inverted directional claim (§3.3). 0.12 has **no docstring at all**, and is unreconciled with the measured `OBSERVED_ACTIVITY_CV = 0.146` |
| **4** | `PANEL` maturation half-times | `generator/unmixing.py:45–52` | 0 / 6 | Ratio-cancellation analysis; the "mOrange2 is the best μ sensor" recommendation; `RECOMMENDED_BUILD.reference_fluorophore` | Wrong for 3 of 5 entries, and two of the correct values are in the repo's own BioNumbers export (§3.4) |
| **5** | `GDCW_PER_OD = 0.3` vs `gdcw_per_od = 0.42` | `run_calibration_nis.py:55` vs `generator/plate.py:84` | — | Biomass units for the filter-calibration result and for every synthetic plate | Two values, 40 % apart, neither cited. Literature range for *S. cerevisiae* is ~0.3–0.66 g DCW/L per OD600, commonly 0.4–0.5; deriving from BNID 100986 (3×10⁷ cells/mL at OD 1) with ~15 pg/cell gives ~0.45. **0.42 is well-supported and uncited; 0.3 sits at the bottom of the range.** *Mitigating and verified:* the NIS docstring's claim that it "cancels out of the NIS, which is a ratio" is **true** — priors and observation model both scale by it, so it drops out of both channels. The `0.42` in `generator/plate.py` has no such excuse |
| **6** | `_MAINTENANCE_PER_ACTIVITY = 300.0` / `1200.0`; `basal_activity = 1.0e-3` | `bridge/latent_bridge.py`, `scripts/run_scenarios.py` | 1 / 4 | The `NGAM B` and Branch-B flux columns of `outputs/scenario_predictions.csv` | Neither cited, 4× apart (§3.5). Heavily quarantined by `acknowledge_unvalidated` and the G4 refusal string — the best-mitigated unsourced number in the repo |
| **7** | `_HEAT_PER_DEGREE = 0.06`, `_HYPOXIA_THRESHOLD = 0.21`, `_CARBON` growth rates | `generator/context.py:50–52, 22` | 1 / 7 | Baseline module activity and growth rate under every non-default culture context; the context-transfer results | `_CARBON` and `_PHASE` **do** carry good sourced docstrings (Walther & Schüller 2001 DOI **verified**; Gasch 2000) — but the Schüller & Entian PMID beside them is broken (§3.1), and the three scalars below them are bare |
| **8** | `panel_calibration.py`'s seven thresholds (`_AGREEMENT_FOLD = 3.0`, `_HEALTHY_FRACTION = 0.7`, `_SATURATION_FRACTION = 0.15`, `_MIN_FIT_QUALITY = 0.5`, …) | `generator/panel_calibration.py:39–45` | 1 / 4 | Whether a design-side constant is judged to agree with a real-plate ladder — i.e. the pass/fail of the repo's own reality check | **None.** Seven bare literals, the largest undocumented cluster in the package, and they set the bar the generator is graded against |
| **9** | `g1_optical` thresholds (`min_od_above_blank = 0.02`, `max_decline_fraction = 0.15`, `min_od_fold = 1.15`, `min_reporter_above_background = 50.0`) | `gates/g1_optical.py:38–42` | 0 / 6 | Which wells are called quantitative — upstream of every growth rate and therefore of every activity | Argued, not sourced. **`od_linear_max = 1.0` is exemplary**: it says outright it is "a placeholder until a dilution series is run: **it is the single most important number to measure for the platform**," and gives the physical reason (0.3–0.5 cm path in a 96-well plate vs the 1 cm cuvette figure). **Honestly asserted — confirmed** |
| **10** | `growth._DEFAULT_DETECTION_OD = 0.02`, `_DEFAULT_RATE_QUANTILE = 0.9`, `_MIN_WINDOW_POINTS = 5`; `reporter` window/polyorder defaults | `growth.py:20, 79–80`; `reporter.py` | 5 / 10, 6 / 23 | The growth window and the smoothing that every μ and every activity passes through | Bare literals in the two lowest-level physics modules. `growth_rate_uncertainty` beside them is a model of the opposite: fully derived, every term measurable on the trace, with the reason it must be derived rather than fitted spelled out |

**Not on this list, and worth saying so.** `OBSERVED_ACTIVITY_CV = 0.146`, `MEASURED_GROWTH_RATE_SE
= 0.0117`, `MEASURED_ACTIVITY_CV = 0.14`, `REFERENCE_GROWTH_RATE = 0.22`, `_DEFAULT_BASAL = 0.9`
and the DTT/H2O2 lethal doses are all measured on this project's own plates, with the sample size,
the spread and the derivation in the docstring. `MEASURED_GROWTH_RATE_SE` even records that an
earlier grid-searched version came out 4× too large and why that was unidentifiable. That is the
standard the rest of the constants should be held to, and it already exists in the codebase.

---

## 7. What to fix first

1. **Correct the 14 PMIDs in `generator/stress_panel.py`, the 1 in `generator/context.py`, and the PMCID in `G4_STATISTICS.md`.** Replacements for 13 of the 15 are in §3.1; `context.py`'s Schüller & Entian 1996 needs the reference identified from scratch or removed. Two hours of work that restores the provenance claim the identifiability result is defended with. — *highest ratio of credibility restored to effort in the repository.*

2. **Make the citation test check truth, not presence.** Replace `assert module.source` with a test that extracts every `PMID \d+` from `MODULES`, `STRESSORS` and `REPORTERS`, resolves it through E-utilities (network-marked, skipped offline), and asserts the returned `sortfirstauthor` and year match the author/year written in the same string. That single test converts the whole class of defect from invisible to impossible. Run it over `docs/` too.

3. **Fix the QUEEN-2m and HyPer7 provenance, because a build recommendation depends on it.** Cite **Takaine et al. 2019, *J Cell Sci*, doi:10.1242/jcs.230649** for QUEEN in yeast; strike "shown in yeast" from the Yaginuma attribution and keep Yaginuma for the Kd/pH/Mg window, which is correct. Give HyPer7 an explicit `demonstrated_in_yeast` value with its own citation rather than letting it inherit the default. Consider making `demonstrated_in_yeast` a **required** field so no reporter can acquire the flag by omission.

4. **Retract or repair `DEFAULT_WELL_CV`'s directional claim.** Either name the publications that show larger spread, or replace the sentence with "the direction of this error is unknown" and stop saying replicate counts are lower bounds. Better: measure it. The docstring already says how — the zero-dose wells of any plate, three technical replicates per construct per plate — and the repo already has those plates.

5. **Fix the fluorophore maturation panel from the repo's own BioNumbers export.** mOrange2 → 4.5 h (Shaner 2008); YFP → 0.65 h (BNID 102974/108677, *S. cerevisiae*); CFP → 0.82 h (BNID 106883/108678, *S. cerevisiae*). Add them to `tests/test_constants_against_bionumbers.py`, which already opens that file. Re-check anything quantitative that depended on the old ratio.

6. **Resolve the `_MAINTENANCE_PER_ACTIVITY` fork.** Delete the unused `300.0` default and make the parameter required, so the one number in play is the one in the script and it cannot silently fork again. Same treatment for `gdcw_per_od` (0.3 vs 0.42) — pick one, source it, put it in `stress_panel.py` or `paths.py` and import it in both places.

7. **Repair the two docstrings that contradict their own code**, and add both to `ARCHITECTURE.md` §6.5: `stress_panel.py`'s opening topology paragraph (§3.6, it contradicts the file and four tests) and `ThermodynamicData`'s "Magnesium binding is NOT modelled here" (§3.7, it is).

8. **Narrow the provenance sentence in `IDENTIFIABILITY.md` and `FINDINGS.md`.** "Transcribed from the literature with PMIDs" is true of the module list and false of the matrix entries. Something like: *"the topology is literature-transcribed and PMID-cited; the crosstalk weights and stressor potencies are expert priors with no per-number source, and the identifiability result is conditional on them."* Then run the `shuffled_loadings` null that `analysis/nulls.py` already provides against the real matrix — the module was built for exactly this question and would convert a provenance argument into evidence.

9. **Document `_POOL_POTENCY`, `_GENERAL_POTENCY`, `_LETHAL_MULTIPLE`, `_LETHAL_HILL`** — even just "chosen so the pool/specific/general orderings separate at plausible dose spacings; not fitted, not sourced" moves them from *dressed up* to *honestly asserted*, which is the whole difference. Do the same for `panel_calibration.py`'s seven thresholds.

10. **Housekeeping:** move `_DEFAULT_BASAL`'s orphan docstring to line 38 where it attaches; import `CORRECTABLE_GDNA_CEILING` / `PASS_MARGIN_CYCLES` in `viz/figures.py` from `qpcr.py` instead of re-declaring; change `equivalence()`'s docstring to state the α its interval actually delivers (or take `confidence=0.90` for a 0.05 TOST); and change the cluster bootstrap's "the standard small-sample correction" to "an approximate widening, not a standard procedure."

11. **Verify the Alberty magnesium constants** (§3.7). `4.19 / 3.17 / 1.88` attributed to "I = 0.25 M" is UNVERIFIED and I could not reach a primary table. Check it against Alberty (2003) or eQuilibrator directly before any adenylate ΔG is quoted.

---

## Appendix — what I could not verify, and what I tried

| Claim | Attempted | Status |
| --- | --- | --- |
| Alberty Mg log K at I = 0.25 M | Web search; eQuilibrator FAQ (TLS certificate failure); Crossref | **UNVERIFIED** |
| MIQE 2009 "mandates no numeric margin" verbatim | Oxford Academic (paywalled to abstract) | **PARTIALLY VERIFIED** — consistent with the checklist as known, not confirmed verbatim |
| Hahn 2006 "RPN4 induction at 0.02 percent MMS" | Web search; Wiley (paywalled) | **UNVERIFIED** — the paper and the YRE/HSE claims are confirmed; this specific dose is not |
| Görner 1998 "nuclear entry survives hog1 deletion" | PMID resolved; full text not re-read | **UNVERIFIED detail**, correct paper |
| `EXPECTATIONS["hac1_total_fold"] = 1.0` ("relative HAC1 abundance unchanged over at least 12 h") | Both PMIDs (34985696, 34985697) resolve correctly to *Methods Mol Biol* chapters on HAC1 splicing detection; chapters are paywalled | **UNVERIFIED**, and broadly consistent with the standard account (HAC1 is regulated by Ire1-mediated unconventional splicing, not by transcription) |
| PMC7646510, PMC3847955 | Not re-read | **UNVERIFIED**, plausible |
| HyPer7 in *S. cerevisiae* | Web search — found *Komagataella phaffii* only | **UNVERIFIED for *S. cerevisiae*** |
| `REPORTERS["TRX2-oxidative"]` "keeps roughly two thirds of its H2O2 induction in yap1"; `REPORTERS["STRE-general"]` "lost in msn2 msn4 and retained in yap1"; "one pH unit shifts a glutathione probe by 60 mV against a genuine peroxide response of 40 to 50" | No reference given in code to check against; the 60 mV/pH half is independently **VERIFIED** (2.303RT/F at 303.15 K = 60.1 mV, which `generator/redox.py::NERNST_MV_PER_PH` computes correctly) | **The 40–50 mV response figure is UNVERIFIED** and is the load-bearing half of the `interpretable()` gate |
