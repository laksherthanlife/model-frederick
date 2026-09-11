# Public data this model can be tested against

Every quantitative claim in this repository is either simulation-derived or measured on two
plates from one lab. [`CLAIM_BOUNDARY.md`](../CLAIM_BOUNDARY.md) records that nothing is
Tier 3. This file is the search for data that would change that.

The twin is host- and product-agnostic: it takes a fluorescence and an optical-density time
series from one well and returns dilution-corrected promoter activity. Nothing in that
contract is about yeast, and nothing is about beta-carotene. The search was run accordingly.

---

## 1. Verdict

**External validation is reachable, and the best dataset is not in E. coli.** It is the
Swain lab's **omniplate** deposit (Edinburgh DataShare, `10.7488/ds/3263`, CC BY 4.0): 819
budding-yeast wells across **nine biological replicate plates**, each carrying OD, two
fluorescence channels, an **untagged wild-type strain** and a media blank on the same plate,
at ~100 timepoints. It is a 24 MB tidy TSV — download and go. It supplies, from outside this
lab, the three things this project cannot supply itself: a **measured per-biomass
autofluorescence term** (Tier 0 → measured), **n ≥ 3 biological replicates** (INCONCLUSIVE →
an interval), and a large external hold-out for the forecast that
[`PREDICTION.md`](PREDICTION.md) already scored and lost. Second is Zaslaver 2009 (PLoS Comput
Biol): 1,920 promoters × 6 media × 52 timepoints, promoter activity **and OD** on the same
wells.

---

## 2. Ranked table

Ranked by what each unlocks, not by size. A dataset that enables a registered
predicted-then-measured claim outranks one that corroborates a constant.

| # | Dataset | What it contains | Licence | URL | Effort | What it unlocks |
|---|---|---|---|---|---|---|
| **1** | **omniplate — HXT deletions** (Montaño-Gutierrez, Manzanaro Moreno & Swain 2021) | *S. cerevisiae*. 9 experiments (separate plates/days), **819 wells**, 69–120 timepoints each. Columns: `OD`, `GFP60`, `GFP80`, `AutoFL`. 35 GFP-tagged strains + **`WT` untagged** + **`Null` media blank**. 2 glucose levels (0.2 % / 2 %) → a real growth-rate contrast. Media blank already subtracted; the per-biomass autofluorescence term is *not*, which is exactly the decomposition in `observation.py`. 13 MB TSV (24 MB zip incl. a raffinose OD-only set) | **CC BY 4.0** (`license_text` inside the zip) | [datashare.ed.ac.uk/handle/10283/4192](https://datashare.ed.ac.uk/handle/10283/4192) · zip: `https://datashare.ed.ac.uk/download/DS_10283_4192.zip` | **Download-and-go.** Tidy long-format TSV, columns already named `experiment / well / condition / strain / time / OD / …` | **Tier 3 is reachable here.** Register a forecast, then score it on 9 held-out plates. Also: the first *measured* autofluorescence anywhere in this project's evidence base; n=9 where the repo has n=2; an external test of the dilution correction across a 10× growth-rate contrast |
| **2** | **Zaslaver et al. 2009 Dataset S1** — genome-wide promoter activity | *E. coli*. **1,920 promoters × 6 media × 52 timepoints**, 16-min resolution over 14 h, exponential → stationary. Each sheet holds **two stacked blocks: `PA` then `OD`**, same rows, same wells. Media: glucose, glycerol, no-AA, phosphate-limited, nitrogen-limited, 4 % ethanol. Includes `U66` empty-vector and `Empty` control rows. 21.5 MB XLS | **CC BY** (PLoS) | [`10.1371/journal.pcbi.1000545.s021`](https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1000545.s021) (verified 21,487,616 B) · mirror: [weizmann.ac.il/…/zaslaver_itzkovitch_excel_output.xls](https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/DownloadableData/zaslaver_itzkovitch_excel_output.xls) | **Light.** Needs `xlrd` for the legacy `.xls` and a 2-block parser (`PA` at row 0, `OD` at row 1922). ~20 lines | An external forecast target at n≈11,500 growth curves, in the exact regime the twin fails in (decelerating culture, `_propagate` gives µ no drift). Also the **prior-art anchor**: `PA = dGFP/dt/OD` is the Alon-lab formalism since 2002, which the README currently presents as rederived |
| **3** | **iGEM InterLab 2018 / Beal et al. 2020** — OD calibration across 244 labs | Per-lab JSON: `microsphere_dilution_calibration`, `CFU_calibration`, `LUDOX_water_calibration`, `fluorescein_dilution_calibration`, `plate_data_sets`. 8 constitutive-GFP *E. coli* strains, 144 OD + 144 fluorescence readings per team. 789 KB zip | **CC BY 4.0** | [`42003_2020_1127_MOESM3_ESM.zip`](https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-020-01127-5/MediaObjects/42003_2020_1127_MOESM3_ESM.zip) (verified, 10 files) | **Light.** JSON, already structured per lab | `od_linear_max` (protocol P4) and the OD→cell-count factor from an external, 244-lab protocol, with a measured *between-instrument* spread. Turns "a placeholder" into "a placeholder with a published distribution around it" |
| **4** | **Guerra et al. 2022** — FP maturation *in budding yeast* | Measured maturation half-times in *S. cerevisiae* by optogenetic pulse + microfluidic single-cell imaging: **mCitrine 10.4 min**, mVenus 20.8, mNeonGreen 11.6, sfGFP 6.9, Cerulean 9.7, mTurquoise2 66.6, mTFP1 76.5, mScarlet-I 31.2, mCherry 52.3, tdTomato 93, mKate2 129 | **CC BY 4.0** | [pubs.acs.org/doi/10.1021/acssynbio.1c00387](https://pubs.acs.org/doi/10.1021/acssynbio.1c00387) | **Table transcription.** Supporting Information is free-access | The repo's reporter is **mCitrine** and its vendored BioNumbers row says *YFP, 39 min* (BNID 102974). That row is Gordon et al. 2007, cycloheximide-chase on generic YFP. Guerra measures **mCitrine specifically, in yeast**: 10.4 min. `ReporterKinetics.k_mat` becomes a measurement, and `archive/code/deconvolve.py`'s maturation-mismatch trick gets real numbers |
| **5** | **Keren et al. 2013 Table S3/S4** | ~900 *S. cerevisiae* + ~1,800 *E. coli* promoter activities across 10 growth conditions incl. **NaCl 1 M, sorbitol 1 M, 39 °C, amino-acid starvation**, with **mean doubling time per condition**. Scalar per promoter × condition (not time-resolved). 586 KB XLSX | **CC BY 3.0** | [europepmc.org/…/PMC3817408/supplementaryFiles](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3817408/supplementaryFiles) → `msb201359-s2.xlsx` | The headline result — 60–90 % of promoters change between conditions by a **global scaling factor set by the condition, not the promoter** — is the growth-rate confound this repo exists to remove, measured genome-wide in yeast under stress. Directly corroborates or refutes the "most of the sensor response was dilution" claim at n≈900 instead of n=4 |
| **6** | **Gerosa et al. 2013 Dataset S2/S3** | *E. coli*, 12 constitutive promoter-GFP reporters × 18 nutrient conditions + arginine-depletion and diauxic-shift time courses, 10-min interval, OD and GFP on the same wells. Promoter activity as `dGFP/(dt·OD)`. 1 MB XLSX (15 time-course sheets) | **CC BY-NC-SA 3.0** — *non-commercial, share-alike* | [europepmc.org/…/PMC3658269/supplementaryFiles](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3658269/supplementaryFiles) → `msb201314-s3.xlsx` | A small, clean, *dynamic* external test of the same formalism: does the twin reproduce their promoter activity from their conditions? Note the licence forbids commercial reuse — check before it lands in a shared artefact |
| **7** | **Gasch et al. 2000 — GEO GSE18** | 156 arrays, *S. cerevisiae* ESR. Curated per-stress subsets are individually downloadable: **GDS113 / GDS31 (DTT time course)**, **GDS17 (H₂O₂)**, GDS30 (diamide), GDS20 (hyper-osmotic), GDS16 (heat shock), GDS108 (menadione), GDS19 (nitrogen depletion) | Public domain (NCBI GEO) | [ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18) | **Light** via GEO SOFT/matrix | The crosstalk topology in `generator/stress_panel.py::reporter_loadings` is Tier 0 "transcribed from the literature". GSE18 is the literature. GDS113/GDS31 and GDS17 are **the same two stressors the plates used**, so the loadings can be re-derived from data rather than transcribed. Moves a Tier 0 block toward Tier 1(b) |
| **8** | **Granados et al. 2018 PNAS** — single-cell TF dynamics | *S. cerevisiae*, single-cell nuclear-translocation time series for **10 transcription factors** (Msn2, Msn4, Dot6, Tod6, Sfp1, Maf1, Mig1, Mig2, Hog1, Yap1) across transitions into **carbon, osmotic and oxidative stress**. 7 JSON files, ~110 MB | Edinburgh DataShare end-user licence (check terms) | [`10.7488/ds/2214`](https://doi.org/10.7488/ds/2214) → [handle/10283/2901](https://datashare.ed.ac.uk/handle/10283/2901) · zip 117 MB verified | The one dataset that tests [`IDENTIFIABILITY.md`](IDENTIFIABILITY.md) §4b — that the **time axis** carries more modules than the channel count. Ten simultaneous latent readouts with genuinely different kinetics (specialists fast/narrow, generalists slow/broad) is the ground truth the Ledermann argument needs. Companion: Granados et al. 2017 eLife, [handle/10283/2708](https://datashare.ed.ac.uk/handle/10283/2708), osmotic-stress survival, 1.9 MB |
| **9** | **yeast-GEM `data/physiology/`** | `chemostatData_Tobias2013.tsv` — glucose/O₂/NH₃ uptake vs experimental growth rate, 11+ steady states, same 4-column format as the `chemostatData_VanHoek1998.tsv` already vendored here. Plus `flux_data_anaerobic.tsv` (Jouhten 2008), `biomassComposition_Lahtvee2017.tsv` | **CC BY 4.0** | [raw.githubusercontent.com/SysBioChalmers/yeast-GEM/main/data/physiology/](https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/main/data/physiology/chemostatData_Tobias2013.tsv) | **Zero effort** — drop-in replacement format. A second held-out physiology check alongside van Hoek and Vos, at no parsing cost. Caveat: it is *in* yeast-GEM, so it is plausibly part of what yeast-GEM was fitted to — check before calling it held out |
| **10** | **Xia et al. 2022** (Nat Commun 13:2819) | Glucose-limited chemostats at **9 dilution rates, 0.025 → ~0.4 /h**, absolute transcriptome + proteome + phosphoproteome + metabolome | CC BY 4.0 | [doi.org/10.1038/s41467-022-30513-2](https://doi.org/10.1038/s41467-022-30513-2) | Medium — multi-omics tables | A genuinely held-out critical-growth-rate test. `CLAIM_BOUNDARY.md` records the van Hoek failure (critical µ ≈ 0.28 vs 0.45 predicted) and names proteome allocation as the missing mechanism. This dataset **measures that allocation** across the transition |
| **11** | **Hackett et al. 2016** (Science 354:aaf2786) | **25 steady-state chemostats** = 5 nutrient limitations × 5 growth rates, with fluxes for 56 reactions, >100 metabolites, 370 enzymes | Check publisher terms (Science) | [science.org/doi/10.1126/science.aaf2786](https://www.science.org/doi/10.1126/science.aaf2786) | Medium–heavy; supplementary tables behind a publisher page | The widest nutrient × growth-rate physiology grid in yeast. Held-out checks across *limitation type*, not just dilution rate |
| **12** | **FlopR example data** (Fedorec et al. 2020) | Tecan Spark and BioTek Cytation raw exports + parsed CSVs: OD, GFP, mCherry time series on the same wells, plus fluorescein/microsphere calibration plates and plate layouts | **GPL-3.0** (package) | [github.com/ucl-cssb/flopr](https://github.com/ucl-cssb/flopr) → `examples/plate_reader/` · Zenodo [10.5281/zenodo.3977408](https://zenodo.org/records/3977408) (CC BY 4.0) | **Download-and-go**, but small | A parser regression corpus: two instrument export formats the repo does not yet read, with the correct parse already committed next to them. Cheap insurance for `plate/synergy.py` |
| **13** | **Beal et al. 2018** (InterLab 2016) and **Beal et al. 2022** (multicolor) | 2018: 92 institutions, plate reader + flow, fluorescein/FITC + colloidal silica, units of **µM FITC/OD** and **MEFL**. 2022: extends calibration to red (sulforhodamine-101) and blue (Cascade Blue) alongside green (fluorescein, 488/530 ± 30) | CC BY 4.0 / CC BY-NC 4.0 | [10.1371/journal.pone.0199432](https://doi.org/10.1371/journal.pone.0199432) · [10.1093/synbio/ysac010](https://doi.org/10.1093/synbio/ysac010) | Light | The absolute-units route. `gain` and `gdcw_per_od` are unidentifiable together in `generator/calibrate.py`; a fluorescein standard curve fixes the fluorescence half of that product in units another lab can reproduce. Note: **6-hour endpoint, not a time series** — this is calibration data, not validation data |
| **14** | **Zaslaver et al. 2006 supplementary** | `All_strains_forWEB.xls` (443 KB, the strain list) and `ClusteredData_Figure3.mat` (765 KB, clustered diauxic-shift promoter activity). Plus **Ronen et al. 2002 SOS data** — 8 operons × 50 timepoints × 4 UV doses, promoter activity only, 37 KB | Not stated on the page — assume all-rights-reserved unless the paper says otherwise | [weizmann.ac.il/mcb/alon/EcoliTranscLibrary/](https://www.weizmann.ac.il/mcb/alon/EcoliTranscLibrary/Zaslaver_supplMat.html) · [sosdata.zip](https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/DownloadableData/sosdata.zip) | Light (`.mat` v5) | Provenance and citation, mostly. The 2006 page carries **no OD channel**; the SOS data is promoter activity only. Both are superseded for validation purposes by #2 |

### Methods worth reading even though they are not data

- **Pavlou, Cinquemani, Geiselmann & de Jong 2022** (*Biophys J* 121:4179) — "Maturation models
  of fluorescent proteins are necessary for unbiased estimates of promoter activity." An
  independent group proving the exact correction `reporter.py::promoter_activity` implements
  is *necessary*, with the size of the bias when it is skipped: RFP promoter activity delayed
  by 39–53 min uncorrected vs ~3–7 min corrected. This is the closest thing to a peer-reviewed
  warrant for the repo's central move, and it is not cited anywhere in the tree.
- **Lichten, White, Clark & Swain 2014** (*BMC Biotechnol* 14:11) — two-wavelength spectral
  unmixing to remove **yeast** autofluorescence in a plate reader. The method that produced
  the `GFP` / `AutoFL` channel pair in dataset #1, and the answer to `observation.py`'s
  standing instruction to measure autofluorescence rather than guess it.
- **Stevenson, Chen, Chew, Dunlop, Swain & Pilizota 2016** (*Sci Rep* 6:38828) — "General
  calibration of microbial growth in microplate readers." The protocol behind P1/P4. Its
  supporting data exists but is 5–17 GB of images ([`10.7488/ds/1980`](https://doi.org/10.7488/ds/1980),
  [`10.7488/ds/1447`](https://doi.org/10.7488/ds/1447)) — read the paper, not the deposit.
- **Swain, Stevenson, Leary, Montaño-Gutierrez, Clark, Vogel & Pilizota 2016**
  (*Nat Commun* 7:13766) — Gaussian-process inference of growth rate *with error bars* from
  OD. A direct comparator for `growth.py`, with 4 small CSVs at [`10.7488/ds/1405`](https://doi.org/10.7488/ds/1405).

---

## 3. The one to do first

**Dataset:** `HXTdeletions_r.tsv` from Edinburgh DataShare `10.7488/ds/3263`.

```bash
curl -L -o /tmp/omni.zip \
  "https://datashare.ed.ac.uk/download/DS_10283_4192.zip"
unzip /tmp/omni.zip -d data/external/omniplate/
# CC BY 4.0 — the licence text ships inside the zip as `license_text`
```

What is in it, read off the file rather than the abstract:

| | |
|---|---|
| Rows | 82,145 |
| Columns | `experiment, condition, strain, time, well, OD, GFP60, AutoFL, GFP80, commontime, hxt` |
| Experiments (plates) | 9 — `20180430_…` through `20180730_…` |
| Wells | 819, median 97 timepoints each (range 69–120) |
| Conditions | `Glu 0.2 %`, `Glu 2 %` |
| Strains | 35 GFP-tagged (HXT1–7 × {WT, *mth1Δ*, *rgt2Δ*, *snf3Δ*, *std1Δ*}), plus `WT` (untagged) and `Null` |
| `Null` | media blank: OD 0.0001 ± 0.0009, GFP60 0.10 ± 1.28 — i.e. **the media blank is already subtracted** |
| `WT` | untagged, real growth (mean OD 1.32), **GFP60 mean 42.1, AutoFL mean 155** — the per-biomass autofluorescence term, measured |
| Tagged strains | GFP60 mean 50 (Hxt4 *snf3Δ*) to 983 (Hxt6 *mth1Δ*) |

Two facts make this the right first job. First, the file's blank convention **is already the
repo's convention** — media subtracted, autofluorescence not — so `observation.py`'s
decomposition maps onto it without reinterpretation. Second, `WT` GFP60 (42) is **up to 85 %
of the dimmest tagged strains' signal**. That is not a rounding term. It is the same
correction `DATA_INVENTORY.md` sweeps blind between 0 % and 30 %, except here it was measured.

### Step 1 — measure the thing that has never been measured (½ day, no registration needed)

Fit `ReporterOptics.autofluorescence` by regressing `GFP60` on `OD` over the `strain == "WT"`
wells, per experiment. Then rerun `reporter.py::promoter_activity` over the 35 tagged strains
twice — once with the fitted term, once with it set to zero — and report how much the
recovered activity ranking moves.

This is the external analogue of the sweep in `DATA_INVENTORY.md` §"Autofluorescence has
never been measured", run on data where the answer is known. It answers a question the repo
currently cannot: *at a realistic autofluorescence share, does omitting the term change which
conclusions you draw?* Publishable as a Tier 2 result — held out from the model, on a
measurement it never saw — the same row as the Vos 2016 retentostat check.

### Step 2 — the registered forecast (this is the Tier 3 route)

`CLAIM_BOUNDARY.md` §"What would move things up a tier" item 4 says a registered forecast is
the only route to Tier 3 and needs no bench work. It has been blocked on the fact that the
only data available is the data the model was built on. That is no longer true.

**Script.** New scripts/run_external_validation.py (proposed). Reuses `analysis/splits.py`,
`analysis/nulls.py::skill_score`, `estimator.py`, and the CRPS/WIS scoring specified in
[`CONTRACT.md`](../CONTRACT.md) Level 2.

**Split.** `heldout_replicate` with `group_key=("strain", "condition")` and the held-out axis
`experiment`. This is the plate-as-unit-of-replication rule
[`SPLITS.md`](../SPLITS.md) §1 already enforces, and the dataset satisfies it honestly: the
nine experiments are nine separate days. Write the manifest with `scripts/make_splits.py` and
**record the SHA-256 before looking at any held-out plate.**

**Prediction.** Condition `estimator.py` on the first 2 h of each held-out well; forecast OD
and GFP60 for the remainder. Same protocol as `PREDICTION.md` §0, so the two results are
directly comparable — and that is the point.

**Score.** CRPS and the weighted interval score against the baseline ladder from
`CONTRACT.md`: mean-predictor, naive uncorrected RFU/OD, and a six-point log-linear
growth-only extrapolator. Report **skill scores**, so losing to the trivial baseline shows as
a negative number. Report 95 % interval coverage against nominal.

**Register first.** Fill in the `PREDICTION.md` §5.7 checklist, commit it, *then* run. The
existing forecast lost — 23 % of level on OD at 2 h, 40–49 % coverage on a nominal 95 %
interval — and the cause is named: `_propagate` gives growth rate no drift, so a decelerating
culture is extrapolated flat. The 0.2 % glucose arm of this dataset decelerates hard. **The
honest expectation is that it loses again**, and a registered, external, negative Tier 3
result is worth more than an unregistered positive one.

**Then the diagnostic.** Two glucose levels means two growth-rate regimes. Score separately.
If the twin loses only in the slow arm, the deficit is exactly the drift term, and the fix is
a growth-rate random walk in `_propagate` rather than wider error bars — which
`PREDICTION.md` already argues would hide the defect rather than fix it.

**Then the confound.** For each of the 35 tagged strains, compare naive `GFP60/OD` fold change
between `Glu 2 %` and `Glu 0.2 %` against the dilution-corrected fold. The README's central
table — "most of the sensor response was dilution" — has never been checked outside this lab.
Here it can be, at n=35 strains × 9 plates, against a growth-rate contrast the experimenters
imposed on purpose.

**Reality check on Step 2.** These are *constitutive-ish* metabolic promoters under a nutrient
shift, not stress reporters under a dose ladder. It tests the growth-dilution machinery and
the forecast, which is the part that has failed. It does **not** test the stress-attribution
layer, and the report must say so.

---

## 4. Gaps — things no public dataset provides

Searched and not found. Each is a reason to run a plate, not a reason to look harder.

1. **A yeast stress-reporter dose-response with OD and fluorescence time series.** This is the
   central gap. Keren 2013 has yeast stress conditions but scalar activities; Gasch has
   yeast stress dynamics but by microarray, no reporter; omniplate has yeast reporter time
   series but nutrient shifts, not stressor doses. **Nothing published combines all three.**
   The NewProtocol plates are, as far as this search reaches, the only dataset of their kind —
   which is an argument for finishing them properly (`DATA_INVENTORY.md` items 1–4), not for
   substituting someone else's.

2. **gDCW per OD₆₀₀ for *S. cerevisiae* with a citable measured range.** Not in the vendored
   BioNumbers export — it has cells/mL per OD (BNID 100986: 3 × 10⁷) and dry-weight-normalised
   quantities, but no OD→gDCW factor. Published values scatter widely because the factor is
   strain-, medium- and reader-specific, which is the honest answer: **it cannot be borrowed.**
   `generator/plate.py` currently carries `gdcw_per_od = 0.42` with no source, and
   `observation.py` refuses without one. Protocol P1's filter-and-dry step is the only fix.

3. **mCitrine autofluorescence at this excitation/emission on this reader.** Guerra 2022 gives
   maturation, not autofluorescence. Lichten 2014 gives the *method*. omniplate gives an
   external magnitude on a different instrument. None of these is a number for this reader —
   and the sweep in `DATA_INVENTORY.md` shows three of eight dose-response readings turn on it.
   **One plate of BY4741 in the mCitrine channel.**

4. **Raw GFP for the Alon promoter library.** Every large promoter-library deposit checked —
   Zaslaver 2006, Zaslaver 2009, Keren 2013, Gerosa 2013 — publishes **promoter activity**,
   already `dGFP/dt/OD`. Zaslaver 2009 alone also publishes OD. Raw fluorescence can be
   reconstructed as `G(t) = G(0) + ∫ PA·OD dt`, but then re-deriving PA from the reconstruction
   is **circular by construction** and cannot test the inversion. Use the OD block for the
   forecast, and the PA block as a target — not as a way to manufacture raw fluorescence.

5. **Flapjack has no public data corpus.** The instance named in the papers,
   `flapjack.rudge-lab.org`, is **NXDOMAIN** as of this search — no DNS record, not merely
   down. The software is alive and MIT-licensed (`flapjacksynbio/pyFlapjack`,
   `flapjacksynbio/flapjack2`, both pushed in 2026), and it is a reasonable *storage* choice.
   It is not a source of validation data today.

6. **Causton et al. 2001 is not in GEO** and the MBoC page returns 403 to automated fetches.
   Its topology contribution is largely redundant with GSE18. Cite it; do not plan on
   downloading it.

7. **The SGD mirror of Gasch 2000 is dead.** `sgd-archive.yeastgenome.org/expression/microarray/…`
   404s at every path tried. **Use GEO GSE18**, which is live and has the per-stress GDS
   subsets already curated.

8. **No public dataset pairs a fluorescent reporter with a qPCR anchor** on the same cultures.
   The G4 anchor problem in [`G4_STATISTICS.md`](G4_STATISTICS.md) has no external check.

9. **Nothing supplies the reward `R`.** `CONTRACT.md` records it as undefined. That is a design
   decision, not a data gap, and no dataset will fill it.

---

## 5. Annotated bibliography

All citations verified against PubMed (esummary) or Crossref on 2026-08-26. All URLs fetched
and confirmed live, with byte counts where a file was downloaded. Nothing here is transcribed
from a search snippet.

### Reporter / promoter-activity time series

**Montaño-Gutierrez LF, Manzanaro Moreno N, Farquhar IL, Huo Y, Bandiera L, Swain PS.**
"Analysing and meta-analysing time-series data of microbial growth and gene expression from
plate readers." *PLoS Comput Biol* 18(5):e1010138, 2022. doi:10.1371/journal.pcbi.1010138.
PMID 35617352. Data: Montaño-Gutierrez, Manzanaro Moreno & Swain (2021), Edinburgh DataShare,
doi:[10.7488/ds/3263](https://doi.org/10.7488/ds/3263), CC BY 4.0. Software:
[pypi.org/project/omniplate](https://pypi.org/project/omniplate/).
— **The best find in this search.** The paper's two autofluorescence corrections (GFP-specific
linear unmixing, and a general untagged-control interpolation) are precisely what
`observation.py` asks for and has never had. Contents verified by download: 82,145 rows, 819
wells, 9 plates, `OD/GFP60/GFP80/AutoFL`, `WT` untagged and `Null` blank present. `GFP80`
tracks `GFP60` at a near-constant ratio ≈ 7.8–8.0 (same emission at a second detector gain);
`AutoFL` is a second, distinct wavelength — the pair is what the Lichten unmixing consumes.

**Zaslaver A, Kaplan S, Bren A, Jinich A, Mayo A, Dekel E, Alon U, Itzkovitz S.** "Invariant
distribution of promoter activities in *Escherichia coli*." *PLoS Comput Biol* 5(10):e1000545,
2009. doi:10.1371/journal.pcbi.1000545. PMID 19851443. CC BY.
— Dataset S1 is described in the paper as "All promoter activities and OD at each condition",
and that is exactly what it is: verified 21,487,616 bytes, 6 sheets (Glucose, no Glucose,
no AA, Phosphate limited, Nitrogen limited, Ethanol), each with a `PA` block at row 0 and an
`OD` block at row 1922, 1,920 rows each. Methods, verbatim: *"Promoter activity was calculated
as the temporal derivative of the background subtracted GFP intensity divided by the OD,
PA = dGFP/dt/OD"*, at 16-min resolution over 14 h, 52 timepoints, Tecan Infinite F200,
384-well, 30 °C. **This sentence is the formalism `reporter.py` rederives**, published in 2009
and traceable to Ronen 2002. The README should cite it rather than imply independence.

**Zaslaver A, Bren A, Ronen M, Itzkovitz S, Kikoin I, Shavit S, Liebermeister W, Surette MG,
Alon U.** "A comprehensive library of fluorescent transcriptional reporters for *Escherichia
coli*." *Nat Methods* 3(8):623–8, 2006. doi:10.1038/nmeth895. PMID 16862137.
— The library. Its supplementary page carries the strain list (443 KB `.xls`) and clustered
diauxic-shift activity (765 KB `.mat`); **no OD channel and no raw GFP.** Licence unstated on
the page.

**Ronen M, Rosenberg R, Shraiman BI, Alon U.** "Assigning numbers to the arrows:
parameterizing a gene regulation network by using accurate expression kinetics." *PNAS*
99(16):10555–60, 2002. doi:10.1073/pnas.152046799.
— The origin of the promoter-activity formalism. Data downloaded and inspected: 8 SOS operons
× 50 timepoints (0–294 min, 6-min interval) × 4 experiments (UV 5 and 20 J m⁻²), **promoter
activity in GFP/min/OD only** — the raw channels were not deposited. Small, but the canonical
citation for the equation.

**Keren L, Zackay O, Lotan-Pompan M, Barenholz U, Dekel E, Sasson V, Aidelberg G, Bren A,
Zeevi D, Weinberger A, Alon U, Milo R, Segal E.** "Promoters maintain their relative activity
levels under different growth conditions." *Mol Syst Biol* 9:701, 2013. doi:10.1038/msb.2013.59.
PMID 24169404. CC BY 3.0.
— Yeast: Tecan Freedom EVO + Infinite F500, **20-min interval, OD600 + YFP + RFP on the same
wells**. But the deposit (Table S3) is **one scalar promoter activity per promoter per
condition**, not the time series — verified by opening `msb201359-s2.xlsx`. Table S4 gives
mean doubling time per condition, which is what makes the growth-rate comparison possible
anyway.

**Gerosa L, Kochanowski K, Heinemann M, Sauer U.** "Dissecting specific and global
transcriptional regulation of bacterial gene expression." *Mol Syst Biol* 9:658, 2013.
doi:10.1038/msb.2013.14. PMID 23591774. **CC BY-NC-SA 3.0.**
— 12 promoter-GFP reporters × 18 conditions, TECAN Infinite M200, 10-min interval, OD and
fluorescence on the same wells. Datasets verified by download: S2 = steady-state activities
across 19 conditions; **S3 = 15 time-course sheets**; S4 = arginine depletion; S5 = simulated
diauxic shift. All are `Promoter activity (GFP/(OD*h))` — **no raw GFP or OD column.**

### Calibration and absolute units

**Beal J, Farny NG, Haddock-Angelli T, Selvarajah V, Baldwin GS, Buckley-Taylor R, Gershater M,
Kiga D, Marken J, Sanchania V, Sison A, Workman CT; iGEM Interlab Study Contributors.**
"Robust estimation of bacterial cell count from optical density." *Commun Biol* 3:512, 2020.
doi:10.1038/s42003-020-01127-5. PMID 32943734. CC BY 4.0. (Author correction PMID 33110148.)
— 244 labs, three protocols (CFU, LUDOX/water, **silica microsphere serial dilution** —
recommended), 8 constitutive-GFP *E. coli* strains. Complete per-lab data verified by
download: 789 KB zip, 10 JSON files including `microsphere_dilution_calibration.json` and
`plate_data_sets.json`. Also gives the **effective linear range** of an instrument, which is
`od_linear_max` by another name.

**Beal J, Haddock-Angelli T, Baldwin G, Gershater M, Dwijayanti A, Storch M, de Mora K,
Lizarazo M, Rettberg R; iGEM Interlab Study Contributors.** "Quantification of bacterial
fluorescence using independent calibrants." *PLoS ONE* 13(6):e0199432, 2018.
doi:10.1371/journal.pone.0199432. PMID 29928012. CC BY 4.0.
— 92 institutions. Fluorescein/FITC + colloidal silica → **µM FITC/OD**; SpheroTech RCP-30-5A
→ **MEFL**. **Endpoint at 6 h, not a time course** — calibration, not validation.

**Beal J, Telmer CA, Vignoni A, Boada Y, Baldwin GS, Shin J, Lee B, Wang R, ... Haddock-Angelli
T.** "Multicolor plate reader fluorescence calibration." *Synth Biol (Oxf)* 7(1):ysac010, 2022.
doi:10.1093/synbio/ysac010. PMID 35949424. **CC BY-NC 4.0.**
— Extends the fluorescein protocol (488 ex / 530 ± 30 em — mCitrine's window) to red
(sulforhodamine-101) and blue (Cascade Blue). The route to calibrating a multi-fluorophore
panel in units another lab can reproduce, which `generator/unmixing.py` will need.
Companion dataset: González-Cebrián A, Borràs-Ferrís J, Boada Y, Vignoni A, Ferrer A, Picó J,
"PLATERO dataset", Zenodo [10.5281/zenodo.7071949](https://zenodo.org/records/7071949),
CC BY 4.0 — a single 85 KB fluorescein calibration workbook.

**Stevenson K, McVey AF, Clark IBN, Swain PS, Pilizota T.** "General calibration of microbial
growth in microplate readers." *Sci Rep* 6:38828, 2016. doi:10.1038/srep38828. PMID 27958314.
— The protocol behind P1/P4: when OD is proportional to cell number and when it is not.
Supporting data is 5.4 GB + 2.5 GB ([10.7488/ds/1980](https://doi.org/10.7488/ds/1980)) and
16.7 GB ([10.7488/ds/1447](https://doi.org/10.7488/ds/1447)) — the paper is the deliverable.

**Fedorec AJH, Robinson CM, Wen KY, Barnes CP.** "FlopR: An Open Source Software Package for
Calibration and Normalization of Plate Reader and Flow Cytometry Data." *ACS Synth Biol*
9(9):2258–2266, 2020. doi:10.1021/acssynbio.0c00296. PMID 32854500. Package GPL-3.0 at
[github.com/ucl-cssb/flopr](https://github.com/ucl-cssb/flopr); Zenodo release
[10.5281/zenodo.3977408](https://zenodo.org/records/3977408), CC BY 4.0.
— `examples/plate_reader/tecan_spark/200228_example_data.csv` and its `_parsed.csv` are a raw
export and its correct parse, side by side, with OD/GFP/mCherry and a layout file. Same for
BioTek Cytation and Neo2.

**iGEM Measurement Tools, "Fluorescence-Tutorials."**
[github.com/iGEM-Measurement-Tools/Fluorescence-Tutorials](https://github.com/iGEM-Measurement-Tools/Fluorescence-Tutorials).
No licence declared on the repository. Contains `Plate Reader Examples/iGEM 2019 Plate Reader
Fluorescence Calibration - Example.xlsx` and a worked MATLAB analysis.

### Reporter kinetics

**Guerra P, Vuillemenot L-A, Rae B, Ladyhina V, Milias-Argeitis A.** "Systematic *In Vivo*
Characterization of Fluorescent Protein Maturation in Budding Yeast." *ACS Synth Biol*
11(3):1129–1141, 2022. doi:10.1021/acssynbio.1c00387. PMID 35180343. CC BY 4.0.
— **The correction the repo needs.** Optogenetic (EL222) transcription pulse + microfluidic
single-cell imaging + delay-differential model with asymmetric-division volume dynamics, in
*S. cerevisiae*. **mCitrine 10.4 min.** mVenus 20.8, mNeonGreen 11.6, sfGFP 6.9, Cerulean 9.7,
mTurquoise2 66.6, mTFP1 76.5, mScarlet-I 31.2, mCherry 52.3, tdTomato 93, mKate2 129 (the
red/cyan ones two-step).

**Gordon A, Colman-Lerner A, Chin TE, Benjamin KR, Yu RC, Brent R.** "Single-cell quantification
of molecules and rates using open-source microscope-based cytometry." *Nat Methods* 4(2):175–81,
2007. doi:10.1038/nmeth1008. PMID 17237792.
— **The primary source behind BNID 102974 and 106883**, which the vendored
`data/kaggle/BioNumbers_Nov2024.csv` carries as *YFP 39 min* and *ECFP 49 min* in
*S. cerevisiae*. Method: pheromone-induced YFP, 30 min of synthesis, then cycloheximide block,
then time to maximum fluorescence at **25 °C**. Reported as 39 ± 7 and 49 ± 9 min, CV < 0.1.
The BNID is faithful; the number is generic YFP at 25 °C, not mCitrine at 30 °C, which is why
Guerra 2022 supersedes it for this repo. (BNIDs 108677/108678 are duplicates of the same pair.)

**Balleza E, Kim JM, Cluzel P.** "Systematic characterization of maturation time of fluorescent
proteins in living cells." *Nat Methods* 15(1):47–51, 2018. doi:10.1038/nmeth.4509.
PMID 29320486.
— 50 cyan-to-far-red FPs by high-precision time-lapse in *E. coli*, two temperatures.
Supplementary Tables 1–2 carry the half-times. Use for any fluorophore Guerra 2022 did not
measure, and flag the host difference when you do.

**Pavlou A, Cinquemani E, Geiselmann J, de Jong H.** "Maturation models of fluorescent proteins
are necessary for unbiased estimates of promoter activity." *Biophys J* 121(21):4179–4188, 2022.
doi:10.1016/j.bpj.2022.09.021. PMID 36146937.
— Independent proof that skipping the maturation term biases promoter activity, with sizes:
GFPmut2 maturation ~9 min; slow-maturing RFP promoter activity delayed **39–53 min** without
correction versus ~3–7 min with it. Data and MATLAB code in the article's supporting material;
no repository DOI.

**Lichten CA, White R, Clark IBN, Swain PS.** "Unmixing of fluorescence spectra to resolve
quantitative time-series measurements of gene expression in plate readers." *BMC Biotechnol*
14:11, 2014. doi:10.1186/1472-6750-14-11. (Verified via Crossref; not returned by a PubMed
title search.)
— The two-wavelength yeast autofluorescence correction, and the reason dataset #1 has both a
`GFP` and an `AutoFL` column.

### Stress-response topology

**Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D, Brown PO.**
"Genomic expression programs in the response of yeast cells to environmental changes."
*Mol Biol Cell* 11(12):4241–57, 2000. doi:10.1091/mbc.11.12.4241. PMID 11102521.
GEO **GSE18**, 156 samples. Curated subsets: **GDS113 / GDS31 (DTT)**, **GDS17 (H₂O₂)**,
GDS30 (diamide), GDS20 (hyper-osmotic), GDS33 (hypo-osmotic), GDS16/GDS15/GDS34/GDS36 (heat
shock), GDS108 (menadione), GDS19 (nitrogen depletion), GDS18 (stationary phase), GDS21
(carbon sources), GDS115 (amino-acid/adenine starvation). Also GSM1080/GSM1081 — Msn2 and
Msn4 overexpression arrays, directly relevant to the ESR channel in
`IDENTIFIABILITY.md` §4a.

**Causton HC, Ren B, Koh SS, Harbison CT, Kanin E, Jennings EG, Lee TI, True HL, Lander ES,
Young RA.** "Remodeling of yeast genome expression in response to environmental changes."
*Mol Biol Cell* 12(2):323–37, 2001. doi:10.1091/mbc.12.2.323. PMID 11179418.
— **Not deposited in GEO** (searched by PMID and by title). The MBoC page returns 403 to
automated fetches. Cite for the topology; do not plan a download.

### Single-cell dynamics

**Granados AA, Pietsch JMJ, Cepeda-Humerez SA, Farquhar IL, Tkačik G, Swain PS.** "Distributed
and dynamic intracellular organization of extracellular information." *PNAS* 115(23):6088–6093,
2018. doi:10.1073/pnas.1716659115. PMID 29784812. Data: Edinburgh DataShare
doi:[10.7488/ds/2214](https://doi.org/10.7488/ds/2214), 8 files (README + 7 JSON, 3.87–48.85 MB),
zip verified at 117,075,062 bytes.
— Ten TFs, three stress transitions, single-cell time series. The generalist/specialist split
(Msn2/4, Tod6, Dot6, Maf1, Sfp1 versus Hog1, Yap1, Mig1/2) is an empirical version of the
crosstalk topology `stress_panel.py` transcribes, and the **timing differences are the signal
`IDENTIFIABILITY.md` §4b argues is untapped.**

**Granados AA, Crane MM, Montaño-Gutierrez LF, Tanaka RJ, Voliotis M, Swain PS.** "Distributing
tasks via multiple input pathways increases cellular survival in stress." *eLife* 6:e21415,
2017. doi:10.7554/eLife.21415. PMID 28513433. Data: doi:10.7488/ds/2043,
[handle/10283/2708](https://datashare.ed.ac.uk/handle/10283/2708), README + 1.89 MB zip.

**Bergen AC, Olsen RA, Sedlackova H, Gasch AP, ...** "Modeling single-cell phenotypes links
yeast stress acclimation to transcriptional repression and pre-stress cellular states."
*eLife* 11:e82017, 2022. doi:10.7554/eLife.82017. PMID 36350693.
— Single-cell yeast stress acclimation from the Gasch lab. Data availability not opened in
this pass — **UNVERIFIED whether a public deposit exists.** Worth ten minutes if the time-axis
route in §4b is pursued.

### Physiology

**SysBioChalmers/yeast-GEM**, `data/physiology/`, CC BY 4.0.
`chemostatData_Tobias2013.tsv` (glucose/O₂/NH₃ uptake vs growth rate, same 4-column layout as
the vendored `chemostatData_VanHoek1998.tsv`), `flux_data_anaerobic.tsv` (Jouhten 2008),
`biomassComposition_Lahtvee2017.tsv`, `aminoAcid_Bjorkeroth2020.tsv`.
— Zero-effort additions. **Caveat:** these ship inside yeast-GEM, so treat "held out" as
unproven until you check the model's fitting history — the same care `CLAIM_BOUNDARY.md`
already takes over the Vos 2016 entry.

**Xia J, Sánchez BJ, Chen Y, Campbell K, Kasvandik S, Nielsen J.** "Proteome allocations change
linearly with the specific growth rate of *Saccharomyces cerevisiae* under glucose limitation."
*Nat Commun* 13:2819, 2022. doi:10.1038/s41467-022-30513-2. PMID 35595797. CC BY 4.0.
— Nine dilution rates, 0.025 → ~0.4 /h, absolute transcriptome + proteome + phosphoproteome +
metabolome. The mechanism the van Hoek failure blames, measured.

**Hackett SR, Zanotelli VRT, Xu W, Goya J, Park JO, Perlman DH, Gibney PA, Botstein D,
Storey JD, Rabinowitz JD.** "Systems-level analysis of mechanisms regulating yeast metabolic
flux." *Science* 354(6311):aaf2786, 2016. doi:10.1126/science.aaf2786. PMID 27789812.
— 25 steady-state chemostats = 5 limitations × 5 growth rates. Widest grid available. Check
Science's reuse terms before redistribution.

**Boender LGM, de Hulster EAF, van Maris AJA, Daran-Lapujade PAS, Pronk JT.** "Quantitative
physiology of *Saccharomyces cerevisiae* at near-zero specific growth rates." *Appl Environ
Microbiol* 75(17):5607–14, 2009. doi:10.1128/AEM.00429-09. PMID 19592533.
— The retentostat series that Vos 2016 (already used here) extends.

**Canelas AB, Harrison N, Fazio A, Zhang J, Pitkänen J-P, van den Brink J, Bakker BM,
Bogner L, Bouwman J, Castrillo JI, ... Nielsen J.** "Integrated multilaboratory systems biology
reveals differences in protein metabolism between two reference yeast strains." *Nat Commun*
1:145, 2010. doi:10.1038/ncomms1150. PMID 21266995.
— Multi-laboratory chemostat physiology for CEN.PK113-7D and S288c. The closest thing to an
inter-lab reproducibility baseline on the physiology side, and a useful sanity bound on how
much of a discrepancy is strain and how much is model.

---

## Provenance

Searched 2026-08-26. Citations verified via NCBI E-utilities `esummary` (PubMed) and the
Crossref REST API; GEO accessions via E-utilities `esearch`/`esummary` on `db=gds`. Every URL
in the tables was issued a request and its status recorded; where a file is described by its
contents (Zaslaver 2009 Dataset S1, the omniplate TSVs, the Gerosa and Keren workbooks, the
Beal 2020 data zip, the Ronen SOS archive) the file was downloaded and opened, not inferred
from the paper. Three negative results — the dead SGD Gasch mirror, the dead Flapjack instance,
and Causton's absence from GEO — are recorded in §4 because a broken link found late costs
more than one found now.

Not searched: PRIDE (proteomics, out of scope for a fluorescence twin), ArrayExpress/BioStudies
beyond what GEO mirrors, and non-English deposits.
