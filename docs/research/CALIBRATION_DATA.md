# Calibration data for COSMIC-dFBA task ordering

## Verdict

**A sufficient dataset exists.** Elizondo & Saa 2025 (PMID 40891387, doi
`10.1021/acssynbio.5c00256`, ACS Synth Biol 14:3457–3472, open access at PMC12455641)
reports, for **three β-carotene-producing CEN.PK2-1c strains at two dilution rates in a
defined-medium carbon-limited aerobic chemostat**, every quantity the §4.3 lexicographic
ordering needs except CO2 and O2:

- **β-carotene and lycopene production rates** separately, in nmol/gDCW/h, with 95% CIs (Table 2)
- **growth rate**, measured as the dilution rate, with 95% CIs (Table S8)
- **specific glucose uptake rate**, in mmol/gDCW/h, with 95% CIs (Table S4)
- **ethanol, acetate and glycerol secretion rates**, same units and table, same CIs
- **biomass**, as steady-state OD600 (Table S9) times a per-strain measured OD→gDCW factor (Table S7)

Six independent (strain, µ) conditions, all on one carbon source, all in the paper's own
units — no OD factor assumed, no unit conversion invented. This is the paper the hunt was
for, and it did not exist when `run_d1.py` was written.

The missing CO2 is the one gap that closes honestly by difference, because it is the
**large** term being derived from **measured small** ones, not the reverse (§4). The
remaining caveats are real but do not block the method: the data are six **steady states**,
not a time course, so they license the task ordering and the lifted product floor — which
is exactly D1's blocker — but not a validated dynamic forecast; and the contents are those
of the López β-car strains (4–32 mg/gDCW), not of the 79 mg/gDCW record holder.

Second-best is **Bubphasawan et al. 2025** (PMID 40944773) for a figure-only but complete
β-carotene shake-flask set, and **Scalcinati et al. 2012** (PMID 22938570) for the only
sufficient set that also carries CO2, O2 and a closed carbon balance — on α-santalene, not
a carotenoid. One paywalled paper, **Verwaal et al. 2010** (PMID 20632327), is the only
other carotenoid chemostat study in this organism and could be a second sufficient dataset
*with* CO2; it could not be obtained (§8). Recommended route: **use Elizondo 2025
directly** (§6).

---

## 1. The sufficient dataset, in full

All identifiers verified against NCBI E-utilities; titles quoted back in §7. Numbers below
were read directly from the article and from the supplementary PDF
(`sb5c00256_si_001.pdf`), not from a secondary source.

**Cultivation.** YNB 1× without amino acids or carbohydrate, supplemented with histidine,
leucine, tryptophan, methionine and uracil for the CEN.PK2-1c auxotrophies; **12 g/L
glucose**; 5 g/L ammonium sulfate. Carbon-limited aerobic stirred-tank bioreactors,
continuous mode, 500 mL working volume, D = 0.1 and 0.25 /h, 200 rpm, 30 ± 0.5 °C, pH 5.0
(1 M NaOH), airflow 0.3 L/min. A **defined medium** — which is why the carbon accounting in
§4 is meaningful at all.

**Table S4 — extracellular fluxes, mmol/gDCW/h, mean (95% CI).** Negative is uptake.

| strain | µ (1/h) | glucose | ethanol | acetate | glycerol |
| --- | ---: | ---: | ---: | ---: | ---: |
| β-car2 | 0.101 | −6.68 [−7.223, −6.176] | 6.40 [5.873, 6.961] | 0.639 [0.582, 0.700] | 0.453 [0.424, 0.484] |
| β-car2 | 0.254 | −14.40 [−15.54, −13.32] | 16.12 [15.21, 17.21] | 0.000 | 1.353 [1.186, 1.526] |
| β-car3 | 0.101 | −7.13 [−7.694, −6.580] | 6.01 [5.725, 6.300] | 0.623 [0.550, 0.698] | 0.592 [0.528, 0.658] |
| β-car3 | 0.254 | −11.72 [−17.79, −5.66] | 15.84 [15.06, 16.64] | 0.000 | 0.000 |
| β-car4 | 0.101 | −5.54 [−5.949, −5.157] | 2.90 [2.697, 3.122] | 0.518 [0.342, 0.697] | 0.573 [0.512, 0.636] |
| β-car4 | 0.254 | −14.28 [−15.69, −12.92] | 14.55 [11.72, 17.68] | 0.000 | 0.000 |

**Table 2 — product rates, nmol/gDCW/h, mean (95% CI).**

| strain | µ (1/h) | lycopene | β-carotene |
| --- | ---: | ---: | ---: |
| β-car2 | 0.101 | 94.4 [85.4, 104.4] | 151.3 [140.3, 163.4] |
| β-car2 | 0.254 | 84.3 [74.1, 94.9] | 184.2 [138.8, 230.5] |
| β-car3 | 0.101 | 472.0 [430.7, 518.5] | 227.1 [207.6, 240.9] |
| β-car3 | 0.254 | 126.5 [123.2, 129.9] | 294.8 [198.6, 391.4] |
| β-car4 | 0.101 | 760.7 [647.5, 877.1] | 185.3 [166.4, 201.8] |
| β-car4 | 0.254 | 392.9 [352.6, 434.6] | 453.6 [361.2, 547.4] |

**Tables S7–S9 — biomass.** OD→gDCW factors, (gDCW/L)/OD600: β-car2 0.533 ± 0.017,
β-car3 0.541 ± 0.004, β-car4 0.498 ± 0.012, each with SE and degrees of freedom. Dilution
rates 0.101 ± 0.00157 and 0.254 ± 0.00149 /h (n = 8 each). Steady-state OD600 per strain
per rate, n = 12–19. So biomass concentration is measured, not assumed, and the factor is
**per strain** — the thing `EXTERNAL_CAROTENOID_BOUND.md` §8 refused to borrow across papers.

**It lifts the floor, and here is the mechanism.** Reading the D1 frontier from
`EXTERNAL_CAROTENOID_BOUND.md` §2, `q_max(µ) = 0.207469 (1 − µ/0.376829)` mmol/gDW/h, the
β-carotene task efficiency is:

| strain | µ | qP measured (mmol/gDCW/h) | qP max, ec model | task efficiency |
| --- | ---: | ---: | ---: | ---: |
| β-car2 | 0.101 | 1.513e-04 | 0.1519 | 0.100 % |
| β-car3 | 0.101 | 2.271e-04 | 0.1519 | 0.150 % |
| β-car4 | 0.101 | 1.853e-04 | 0.1519 | 0.122 % |
| β-car2 | 0.254 | 1.842e-04 | 0.0676 | 0.272 % |
| β-car3 | 0.254 | 2.948e-04 | 0.0676 | 0.436 % |
| β-car4 | 0.254 | 4.536e-04 | 0.0676 | 0.671 % |

β-carotene is the **lowest-efficiency task by three orders of magnitude**, so the ordering
ranks it last, pins biomass and ethanol at their measured values first, and the product
floor rises to whatever the pinned tasks leave. That is precisely the mechanism §4.3
describes. The `qP max` column is indicative only — it is the ec model's frontier at a
protein-pool calibration fitted elsewhere, not an LP run under Elizondo's own fixed uptake.
Running that LP is the first implementation step, not a literature question.

---

## 2. Candidate table

Ranked by sufficiency. Every PMID and DOI resolved against NCBI E-utilities or Crossref;
titles in §7. "DERIVABLE" states what from. Requirements: (a) titre or content, (b) growth
or biomass over time, (c) carbon uptake rate, (d) secreted byproduct rates.

| # | Paper | strain / product / mode | (a) | (b) | (c) | (d) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **1** | **Elizondo & Saa 2025** · PMID 40891387 | β-car2/3/4 (CEN.PK2-1c) · β-carotene + lycopene · **carbon-limited aerobic chemostat**, 12 g/L glucose, D = 0.101 & 0.254 /h | **REPORTED** 151–454 nmol/gDCW/h β-car, 84–761 lycopene, 95% CI (Table 2) | **REPORTED** µ = D, 0.101 ± 0.0016 and 0.254 ± 0.0015 /h (S8); OD600 (S9) × OD→gDCW 0.498–0.541 (S7) | **REPORTED** 5.54–14.40 mmol/gDCW/h, 95% CI (S4) | **REPORTED** ethanol 2.90–16.12, acetate 0–0.639, glycerol 0–1.353 mmol/gDCW/h (S4). **CO2, O2 ABSENT** | **sufficient** |
| 2 | Scalcinati et al. 2012 · PMID 22938570 | CEN.PK113-5D derivatives · **α-santalene** (sesquiterpene, not a carotenoid) · two-phase ISPR glucose-limited chemostat, D = 0.05 & 0.10 /h | REPORTED total santalene mg/24h/L; rate in Cmmol/g biomass/h; yield Cmmol/Cmmol glucose | REPORTED D = µ 0.048–0.10 /h; Y_sx g/g | REPORTED r_s mmol/g/h (Table 1) | REPORTED r_CO2, r_O2, r_etoh, r_acet mmol/g/h, RQ, residual glucose, **carbon recovery %**. Glycerol measured by HPLC (Aminex HPX-87H) but **not tabulated** | **sufficient, wrong product** — the only set found that carries CO2 and O2 |
| 3 | Bubphasawan et al. 2025 · PMID 40944773 | Sp_Bc Δgal80 · β-carotene · shake flask, 2 % sucrose sole C source, triplicate; plus 5-L fed-batch | REPORTED content mg/gDCW **and** titre mg/L over time (Fig 2A,B); 2.29 ± 0.16 mg/gDCW | REPORTED DCW g/L over time (Fig 2C) | **DERIVABLE** from residual sucrose/fructose/glucose over time (Fig 2D) ÷ DCW. Sucrose, not glucose | **REPORTED** ethanol g/L (Fig 2E) **and glycerol g/L (Fig 2F)** over time, by HPLC (Shodex SH1011, RID). Acetate, CO2 ABSENT | **partially sufficient** — everything present but **only as figure panels**; medium is YP-based so total carbon in is not accounted |
| 4 | Nunta/Watcharawipas group 2024 · PMID 38921419 | Sp_Bc-Can001 Δgal80 · **canthaxanthin** · 5-L fed-batch, semi-defined, 20 g/L sucrose batch + feed, 120 h, duplicate | REPORTED 60.36 ± 1.51 mg/L @ 96 h; titre over time (Fig 7a) | REPORTED DCW g/L over time (Fig 7a) | **DERIVABLE** from residual sucrose/fructose/glucose (Fig 7a) ÷ DCW | **REPORTED** ethanol peak 26.71 ± 0.00 g/L @ 72 h; glycerol 7.94 ± 0.52 g/L at end; both over time (Fig 7a). Acetate, CO2 ABSENT | **partially sufficient** — figure-only; complex N (30 g/L yeast extract + peptone) |
| 5 | López et al. 2019 · PMID 31380362 | SM14 (S288c) & βcar1.2 (CEN.PK2-1c) · β-carotene · 250 mL YPD shake flask **and** 1-L defined-medium fed-batch, 20 g/L glucose batch + 450 g/L feed | REPORTED 21 & 5.8 mg/gDCW flask; 6.48 & 6.91 mg/gDCW and 209.0 & 739.6 mg/L fed-batch (Table 2); content over time (Figs 1C,D) | REPORTED OD600 every 2–3 h, 1 OD600 = 0.4 g/L measured; µ_max 0.43 & 0.40 /h; fed-batch µ_set 0.13 → 0.03 /h | **DERIVABLE** from HPLC glucose time course ÷ biomass (Figs 1C,D, 2A,B). Fed-batch feed is analytic but needs Y_sx, which is **not given numerically** | **PARTIAL** — glucose, **ethanol** and acetic acid all by HPLC over time (ethanol 12.3 g/L @ 68 h, >18 @ 77 h). **Glycerol ABSENT (0 mentions). CO2 ABSENT** | **partially sufficient** — glycerol and CO2 missing; no supplementary data, so all rates need figure digitisation |
| 6 | López et al. 2020 · PMID 33102463 | CEN.PK2 β-Car/β-iono series · β-carotene + β-ionone · YPD shake flask, 72 h | REPORTED 4 → 16 → 32 mg/gDCW total carotenoids | REPORTED growth curves (Fig 3B, Table S3) | **ABSENT as a measurement** — q_glucose is *computed* as the optimal minimum by the GSMM and compared against **Van Hoek 1998**, a non-carotenoid strain | **ABSENT** — no ethanol (0 measurements), no glycerol (0 mentions), no acetate. The 80–100 % "carbon recovery" is **model-computed, not measured** | **insufficient** for (c) and (d) — but see §5, it is the published precedent for the composition route |
| 7 | Reyes et al. 2014 · PMID 24262517 | SM14 lineage · total carotenoids · batch | REPORTED 6 mg/gDCW glucose phase, 15 mg/gDCW ethanol phase | not verified (paywalled, no PMC) | UNVERIFIED | UNVERIFIED — the phase split implies glucose and ethanol were followed | **UNVERIFIED** — worth a paywalled look; it is SM14's origin |
| 8 | Olson et al. 2016 · PMID 27423881 | SM14 and YLH2 · β-carotene · batch bench-top bioreactor, YNB glucose, C:N 8.8 and 50 | REPORTED 16.8 ± 1.8 and 25.52 ± 2.15 mg/gDCW (SM14); 5.68 ± 1.24 and 22.58 ± 0.11 (YLH2) | REPORTED biomass traces (Figs 1, 2b, 5) | **DERIVABLE** from the glucose time course (Figs 1, 5) ÷ biomass; no q_S tabulated | **PARTIAL** — ethanol and acetate time courses REPORTED. **Glycerol ABSENT. CO2 ABSENT** | **partially sufficient** — read from the publisher landing page, not the PDF; treat the panel-level claims as provisional |
| 8b | **Verwaal et al. 2010 · PMID 20632327** | two carotenogenic *S. cerevisiae* strains · carotenoids · **carbon-limited chemostat** | UNVERIFIED — "different absolute carotenoid levels", no units in the abstract | UNVERIFIED | UNVERIFIED | UNVERIFIED | **UNVERIFIED — the highest-value unresolved lead.** The Daran/Delft-affiliated group conventionally tabulates full chemostat physiology (q_S, q_EtOH, q_glycerol, q_acetate, q_CO2, q_O2), so this is *plausibly* a second sufficient dataset. Wiley paywalled; the green-OA record (OpenAlex, Semantic Scholar) resolves to a WUR landing page carrying **no file**. Do not assume it |
| 9 | Bu et al. 2020 · PMID 33062054 | YBX-01→20 · β-carotene · YPD + 10 g/L acetate, flask, dodecane overlay | REPORTED 149.8 mg/L intracellular | REPORTED growth curve (Fig 2b) | **ABSENT as a rate** — glucose consumption plotted, no biomass-normalised rate | **ABSENT** — no glycerol (0), no ethanol rate, no acetate consumption rate | insufficient |
| 10 | Fathi et al. 2021 · PMID 33332529 | *S. cerevisiae* + LIP2/7/8 · β-carotene 46.5 mg/gDCW · YNB + 1 % olive oil, flask | REPORTED | REPORTED final DCW only (0.33–0.35 g/L for the control) | **ABSENT** — 0 mentions of glucose, uptake or consumption; substrate is triacylglycerol | **ABSENT** — glycerol appears only as a hydrolysis product in discussion; no rates | insufficient |
| 11 | Su et al. 2020 · PMID 32236768; Sun et al. 2020 · PMID 33616900; Arhar et al. 2024 · PMID 39215465; ~~Bubphasawan et al. 2024 · PMID 38710418~~ (**WITHDRAWN** by the publisher); Verwaal 2007 · PMID 17496128; Yamano 1994 · PMID 7765036; Li 2013 · PMID 23718229; Xie 2015 · PMID 25475893; Rabeharindranto 2019 · PMID 30723675; Zhao 2021 · PMID 33605428 | various · β-carotene or lycopene · flask, bioreactor, fed-batch | REPORTED in all cases | partial to absent | **ABSENT in all** | **ABSENT in all** | insufficient — titre-only or titre-plus-biomass; see `EXTERNAL_CAROTENOID_BOUND.md` §3 for the full 30-paper survey, which was checked and not repeated here |

**The search space is small, which is why this is a short list.** PubMed
`(chemostat OR "continuous culture") AND (carotenoid OR carotene OR lycopene) AND
Saccharomyces` returns **exactly two records**: Elizondo 2025 and Verwaal 2010. The axis is
exhausted rather than under-searched, and two independent searches — one weighted to
chemostat and fed-batch, one to flux-analysis and GEM-coupled papers — converged on the
same single sufficient paper, including on the carbon-limitation caveat recorded in §8.

**Three useful negatives.** No 13C-MFA of a carotenoid-producing *S. cerevisiae* exists;
~15 Europe PMC full-text query variants over carotene/lycopene/isoprenoid/farnesene/
squalene/amorphadiene/santalene/bisabolene returned only reviews and other organisms.
**COSMIC-dFBA itself contains no yeast data**: PMID 38387677 was calibrated and validated
exclusively on CHO-cell fed-batch bioreactor runs producing a monoclonal antibody, and
`github.com/LewisLabUCSD/COSMIC-dFBA` ships `COSMIC_dFBA.m` and a licence with no data
files — there is no worked yeast example to copy. And **PMID 39537637** (Nat Commun 15:9844,
isopentenol utilization pathway in yeast) looks promising on title but its only
`mmol/gDCW/h` figures are **FBA bounds**, not measurements.

**Named non-*cerevisiae* alternatives**, for completeness and not as candidates: Larroude
et al. 2018 (PMID 28986998), 6.5 g/L and 90 mg/gDCW β-carotene in *Yarrowia lipolytica*
fed-batch at 0.048 g/g glucose — a partial (c), no (d); two *Y. lipolytica* chemostat papers
(PMIDs 35955650, 37742215) have q_glucose and byproducts but the product is lipid. None may
be compared against a *S. cerevisiae* GSMM, for the reason
`EXTERNAL_CAROTENOID_BOUND.md` §8 already gives.

---

## 3. Method note, because the next search will need it

Uptake rates are never in abstracts, so abstract search cannot find them. What worked was
**Europe PMC REST full-text search with bare quoted phrases** (the `FULL_TEXT:` field prefix
returns zero hits and is unusable):

```
.../rest/search?query=%22carotene%22%20AND%20%22specific%20glucose%20uptake%20rate%22&format=json
```

The single most productive move was searching for the **instrument** rather than the
analyte: `"Aminex"`, `"Shodex SH1011"`, `"refractive index"` find papers that measured
glucose, glycerol, ethanol and acetate in one injection, which predicts a usable dataset far
better than any word about carotenoids.

---

## 4. The derivation arithmetic: can the missing rates be inferred?

**The answer depends entirely on which direction you derive in, and the asymmetry is
brutal.** Basis: β-carotene 536.87 g/mol, 40 C; ethanol 46.068, 2 C; glycerol 92.094, 3 C;
acetate 60.052, 2 C; glucose 180.156, 6 C; biomass CH1.8O0.5N0.2 at 24.6 g per C-mol.

### 4a. Deriving the product by difference — hopeless

Worked on the best candidate that *needs* it, López 2019's SM14 fed-batch at 68 h
(PMID 31380362, Table 2 plus the ethanol and residual-glucose figures quoted in its text):

| term | as reported | mmol C/L | vs product C |
| --- | ---: | ---: | ---: |
| biomass | 32.3 g/L | 1313.0 | 84× |
| ethanol | 12.3 g/L | 534.0 | **34.3×** |
| residual glucose | 2.4 g/L | 79.9 | 5.1× |
| β-carotene | 209.0 mg/L | **15.57** | 1× |
| measured non-CO2 total | | 1942.5 | |

β-carotene is **0.80 %** of the measured non-CO2 carbon, and 0.65–0.72 % of the glucose fed
(taking Y_xs in the plausible aerobic window 0.45–0.50 g/g, since the feed integral is not
reported). The error budget follows immediately:

| relative error on ethanol | mmol C/L | as a multiple of the whole product carbon |
| ---: | ---: | ---: |
| 2 % | 10.7 | 0.69× |
| 5 % | 26.7 | **1.71×** |
| 10 % | 53.4 | 3.43× |
| 20 % | 106.8 | 6.86× |

**A 2 % error on ethanol alone already consumes two-thirds of the entire product carbon.**
To recover the product flux to 10 % by difference you would need ethanol to **0.29 %
relative** — better than HPLC does on a clean standard, and far better than any
difference-based inference. The residual does not merely swallow the product flux; it is
34× larger than it. Route 1, as posed, is dead. Say so and stop.

The same arithmetic kills **inferring glycerol** by difference. On López 2019's βcar1.2
fed-batch, a glycerol yield anywhere in the ordinary aerobic window gives:

| Y_glycerol (g/g) | g/L | mmol C/L | vs product C |
| ---: | ---: | ---: | ---: |
| 0.01 | 2.25 | 73.4 | 1.3× |
| 0.03 | 6.76 | 220.3 | 4.0× |
| 0.05 | 11.27 | 367.2 | 6.7× |

The window itself spans five times the product carbon, and CO2 — the largest sink — is
unmeasured, so there is nothing to difference against.

### 4b. Deriving CO2 by difference — sound, and it closes the one gap in the sufficient dataset

Reverse the direction and the arithmetic reverses with it. On Elizondo's defined medium
every other term is measured, so CO2 is the residual, and a residual only needs to be
known roughly for the ordering to rank it:

| strain | µ | C in (glucose) | C out measured | **CO2 by difference** | product C as % of C in |
| --- | ---: | ---: | ---: | ---: | ---: |
| β-car2 | 0.101 | 40.08 | 19.55 | **20.53** | 0.0245 % |
| β-car2 | 0.254 | 86.40 | 46.63 | 39.77 | 0.0124 % |
| β-car3 | 0.101 | 42.78 | 19.18 | 23.60 | 0.0654 % |
| β-car3 | 0.254 | 70.32 | 42.02 | 28.30 | 0.0240 % |
| β-car4 | 0.101 | 33.24 | 12.70 | 20.54 | 0.1138 % |
| β-car4 | 0.254 | 85.68 | 39.46 | 46.22 | 0.0395 % |

All in mmol C/gDCW/h. Precision: the glucose 95% CI on β-car2 at µ = 0.101 is ±7.8 %,
which is ±3.14 mmol C/gDCW/h on the input and therefore roughly **±15 % on derived CO2** —
useless for a product flux, entirely adequate for a term whose only job is to be ranked
against biomass and ethanol. Sanity check: β-car2 at µ = 0.101 gives fermentative CO2 equal
to ethanol (6.40) plus respiratory 14.1, and Van Hoek's most fermentative state
(D = 0.40 /h) measured q_CO2 = 18.9 mmol/gDCW/h, so 20.5 is the right order for a culture
more fermentative still.

**The general rule this establishes: never derive the small term from the large residual;
derive the large residual from the measured small ones.** Route 1 fails for the product and
for glycerol, and succeeds for CO2 — and CO2 is the only term the sufficient dataset lacks.

---

## 5. Could a non-carotenoid dataset calibrate the non-product tasks?

Yes, it is defensible, it has a published precedent, and it is now unnecessary.

**The precedent.** López 2020 (PMID 33102463) did exactly this composition: it constrained
a contextualised *iMM904* GSMM with its own measured β-carotene, β-ionone and biomass
yields, and took **the specific growth and CO2 rates from Van Hoek, van Dijken & Pronk
1998** — a non-carotenoid strain — then computed the optimal minimum specific glucose
uptake and found it agreed with Van Hoek's measurement to ~10 % worst case, 3.3 % average,
with 80–100 % carbon-balance closure. So a referee has already accepted the move in print.

**The repo's own copy of that reference is a truncated subset.**
`data/physiology/chemostatData_VanHoek1998.tsv` holds four columns and four rows:

```
Drate  GlucoseUptake  O2uptake  CO2production
0.025  0.30  0.80  0.80   …  0.15  1.70  3.90  4.20
```

D ≤ 0.15 /h — entirely below the critical dilution rate, so **every row has zero ethanol and
zero glycerol by construction**, and neither column is even present. The source paper's
Table 1 (PMID 9797269) runs to D = 0.40 /h and tabulates q_glucose, q_O2, q_CO2,
**q_ethanol, q_glycerol, q_acetate and q_pyruvate**, with ethanol onset at D_crit = 0.28 /h
(0.11 at 0.28, 2.3 at 0.30, 13.9 at 0.40 mmol/gDCW/h) and glycerol appearing only at
D ≥ 0.35. Extending the vendored file from its own cited source costs an hour and no wet
lab. The `yeast-GEM` sibling files do not help: `chemostatData_Tobias2013.tsv` carries only
`GLCxtI / O2xtI / NH3xtI / growth`, with no ethanol or glycerol column either.

**Where the composition is defensible, and where it is not.** A glucose-limited exponential
fed-batch at µ = 0.1 /h is metabolically close to a glucose-limited chemostat at
D = 0.1 /h, so borrowing is legitimate *when the regimes match*. They frequently do not:

- López 2019's **βcar1.2** at µ_set 0.13 → 0.03 /h behaves as Van Hoek's CBS 8066 does —
  respiratory, no ethanol accumulation. Composition is defensible.
- López 2019's **SM14** produced **12.3 g/L ethanol at µ ≈ 0.1 /h**, where Van Hoek's strain
  produces exactly zero below 0.28 /h. The paper names the cause: an "exacerbated Crabtree
  effect". Composition is refuted by construction for this strain.
- Elizondo's own strains are worse still: **q_ethanol = 6.40 mmol/gDCW/h at µ = 0.101**,
  against Van Hoek's zero, and q_glucose 6.68 against Van Hoek's 1.1 at D = 0.10 — a
  **6× higher uptake at the same growth rate**, i.e. a biomass yield roughly six times
  lower.

That last line is the finding that settles the question. **The non-product tasks are not
strain-independent in the strains that matter.** Ethanol is the second-highest-priority
task in the ordering, and a borrowed ethanol efficiency of zero versus a measured 6.40
inverts the ranking outright. Composition mixes strains and conditions **incoherently for
exactly the carotenoid strains this project cares about**, and would have quietly produced
a wrong ordering. Do the Van Hoek extension anyway — it is cheap, it fixes a truncated
vendored file, and it gives a respiratory reference to contrast against — but do not
calibrate the ordering from it.

---

## 6. The three routes, ranked, and the recommendation

| route | cost | what it actually licenses | verdict |
| --- | --- | --- | --- |
| **1. Derive by carbon balance** | zero | **Nothing** for the product or glycerol: the residual is 34× the product carbon and 0.29 % ethanol precision is unobtainable (§4a). **CO2 only**, to ±15 %, which is enough (§4b) | **rejected as a substitute; adopted for CO2 alone** |
| **2. Compose — non-product tasks from a standard physiology dataset** | ~1 h to extend the vendored Van Hoek file from its own source | A method demonstration on a **chimera**. Refuted for SM14 and for all three Elizondo strains, whose ethanol and glucose rates differ from the reference by ∞ and 6× at the same µ (§5) | **do the file extension; do not calibrate from it** |
| **3. Minimum wet-lab experiment** | non-product half: ~1 day, reagents only. Product half: **blocked** — no strain in this lab carries the carotenoid pathway (`docs/PARKED.md`, `docs/ARCHITECTURE.md`) | Everything, for this project's own strain, medium and instruments — but not until a producer exists | **the eventual answer, not the current one** |

### Route 3 in concrete form, since it will be needed eventually

**Non-product half, runnable now.** One shake flask of BY4741 in the standard medium.
Sample at 6 timepoints spanning exponential growth, the diauxic shift and post-diauxic
consumption — the shift is not optional, because Reyes 2014 and López 2019 both show
carotenoid content quadrupling *after* glucose is exhausted, so a single exponential-phase
ordering would be calibrated in the wrong phase. At each timepoint: OD600 on the plate
reader, and **one HPLC injection of the supernatant**. An Aminex HPX-87H or Shodex SH1011
column with dilute acid eluent and a refractive-index detector returns **glucose, glycerol,
ethanol and acetate in a single ~20 min run** — this is the method in PMID 41376837 and
PMID 40944773, and it is why "HPLC on the supernatant at N timepoints" genuinely does cover
requirements (c) and (d) at once. Three biological replicates: **18 injections plus
standards, one day, one instrument**.

That licenses: measured task efficiencies for biomass, ethanol, glycerol and acetate on
*this* strain in *this* medium on *this* instrument, with replicate error bars, plus CO2 by
difference per §4b. It does **not** license anything about the product task, any other
medium, or the claim that the ordering is condition-independent — for that, run the same
flask at a second glucose concentration and you have two orderings to compare.

**Product half.** The identical protocol on a carotenoid producer, adding a pellet
extraction and the carotenoid HPLC run already scoped in `docs/PARKED.md` (β-carotene
450 nm, plus phytoene 285 nm and lycopene 470 nm from the same injection). It cannot be
scheduled until a strain exists, and strain acquisition — not measurement — is the cost.
One gDCW-per-OD600 calibration on that strain is mandatory and cheap: pigmented cells
scatter light differently from wild type, which is why Elizondo measured the factor per
strain and why `EXTERNAL_CAROTENOID_BOUND.md` §8 refused to borrow one.

### Recommendation

**Use Elizondo & Saa 2025 (PMID 40891387).** It is free, it is open access, its
supplementary PDF is one HTTP request, and it is the only source found that reports all four
required quantities as **biomass-normalised rates with confidence intervals** on
β-carotene-producing *S. cerevisiae* in a defined medium — no digitising, no borrowed OD
factor, no unit conversion of ours between the paper and the model. Concretely:

1. Vendor Tables 2, S4, S7, S8 and S9 as a TSV under `data/physiology/`, alongside the
   Van Hoek file and in the same shape.
2. Close CO2 by difference per §4b and record it as derived, not measured.
3. Fix the six uptake vectors, run the six LPs, and compute the task efficiencies for real
   instead of against the borrowed ec frontier used in §1.
4. Extend `chemostatData_VanHoek1998.tsv` to D = 0.40 /h with the ethanol, glycerol,
   acetate and pyruvate columns its own source already publishes. Use it as the
   **respiratory contrast**, not as calibration.
5. Book the BY4741 flask from §Route 3 anyway. It is one day, and it is the only work here
   whose output is still valid after the carotenoid strain lands.

Do not spend money on Route 2 or wait on Route 3 to start. The dataset exists.

---

## 7. Identifier verification

Every identifier cited above was resolved and the returned title checked against the claim.
Resolved via NCBI E-utilities `esummary`; DOIs as returned by the same records.

| PMID | DOI | title as returned |
| --- | --- | --- |
| 40891387 | 10.1021/acssynbio.5c00256 | Complex Kinetic Models Predict β-Carotene Production and Reveal Flux Limitations in Recombinant *Saccharomyces cerevisiae* Strains |
| 38387677 | 10.1016/j.ymben.2024.02.012 | COSMIC-dFBA: A novel multi-scale hybrid framework for bioprocess modeling |
| 22938570 | 10.1186/1475-2859-11-117 | Combined metabolic engineering of precursor and co-factor supply to increase α-santalene production by *Saccharomyces cerevisiae* |
| 40944773 | 10.1186/s40643-025-00936-y | Sustainable β-carotene production by engineered *S. cerevisiae* using sucrose and agricultural by-products |
| 38921419 | 10.3390/jof10060433 | Metabolic Engineering of *Saccharomyces cerevisiae* for Production of Canthaxanthin, Zeaxanthin, and Astaxanthin |
| 31380362 | 10.3389/fbioe.2019.00171 | Build Your Bioprocess on a Solid Strain—β-Carotene Production in Recombinant *Saccharomyces cerevisiae* |
| 33102463 | 10.3389/fbioe.2020.578793 | Engineering *Saccharomyces cerevisiae* for the Overproduction of β-Ionone and Its Precursor β-Carotene |
| 9797269 | 10.1128/AEM.64.11.4226-4233.1998 | Effect of specific growth rate on fermentative capacity of baker's yeast |
| 24262517 | 10.1016/j.ymben.2013.11.002 | Improving carotenoids production in yeast via adaptive laboratory evolution |
| 27423881 | 10.1007/s10295-016-1808-9 | Characterization of an evolved carotenoids hyper-producer of *Saccharomyces cerevisiae* through bioreactor parameter optimization and Raman spectroscopy |
| 33062054 | 10.1186/s13068-020-01809-6 | Engineering endogenous ABC transporter with improving ATP supply and membrane flexibility enhances the secretion of β-carotene in *Saccharomyces cerevisiae* |
| 33332529 | 10.1093/femsyr/foaa068 | Metabolic engineering of *Saccharomyces cerevisiae* for production of β-carotene from hydrophobic substrates |
| 41376837 | 10.1016/j.mec.2025.e00267 | Model-guided chemical environment and metabolic network design to couple pathways with cell fitness |
| 20632327 | 10.1002/yea.1807 | Heterologous carotenoid production in *Saccharomyces cerevisiae* induces the pleiotropic drug resistance stress response |
| 39537637 | 10.1038/s41467-024-54298-8 | Yeast metabolism adaptation for efficient terpenoids synthesis via isopentenol utilization |
| 28986998 | 10.1002/bit.26473 | A synthetic biology approach to transform *Yarrowia lipolytica* into a competitive biotechnological producer of β-carotene |

The remaining PMIDs in row 11 of §2 (17496128, 7765036, 23718229, 25475893, 30723675,
32236768, 33605428, 33616900, 38710418 — **WITHDRAWN**, see below — and 39215465) were
verified in the same run and are already recorded in
`docs/EXTERNAL_CAROTENOID_BOUND.md` §3; they are not re-tabulated here.

**PMID 38710418 (Bubphasawan 2024) was retracted by Elsevier after this survey was made.**
It appears in row 11's "insufficient" bucket, so nothing here was ever calibrated on it, but
it is struck rather than silently left in a citation list.

---

## 8. What could not be established

- **Verwaal et al. 2010 (PMID 20632327) is the outstanding lead and is unresolved.** It is
  the only carotenoid **carbon-limited chemostat** study in *S. cerevisiae* other than
  Elizondo, and its group conventionally publishes a full chemostat physiology table. If it
  does, it is a second sufficient dataset **with CO2 and O2** — the one thing Elizondo
  lacks. Wiley paywalled; the green-OA record resolves to a Wageningen landing page with no
  file attached, and academia.edu returns 403. One institutional or interlibrary pull would
  settle it, and it is the single highest-value paywalled request on this list.
- **Olson et al. 2016 (PMID 27423881) is resolved only to panel level.** Its contents,
  biomass traces, glucose time course and ethanol/acetate traces were read from the
  publisher landing page rather than the PDF, so the figure-level claims in §2 row 8 are
  provisional. Glycerol and CO2 are absent either way, so it cannot become sufficient.
- **Reyes et al. 2014 (PMID 24262517) is unresolved.** Its glucose-phase versus
  ethanol-phase content split (6 vs 15 mg/gDCW) implies both were followed; the underlying
  traces may exist in it. Paywalled, no PMC. **UNVERIFIED.**
- **Van Hoek 1998 Table 1's ethanol, glycerol, acetate and pyruvate values quoted in §5
  were read out of the PMC page by an automated fetch, not transcribed from the PDF by
  hand.** The confidence check that they are right: the fetched D = 0.025 row
  (q_glucose 0.3, q_O2 0.8, q_CO2 0.8) matches the vendored TSV's first row exactly on all
  three columns. That is corroboration, not verification. Transcribe from the PDF before
  vendoring the extension.
- **Elizondo's chemostats give steady states, not a trajectory.** Six (strain, µ) points
  license the task-efficiency ordering and the lifted product floor. They do not license a
  validated *dynamic* forecast, which is what COSMIC-dFBA ultimately produces. Whether a
  static ordering transfers to a fed-batch trajectory is untested here and is the next open
  question, not a settled one.
- **The `qP max` and task-efficiency columns in §1 are indicative.** They use the ec-model
  frontier from `EXTERNAL_CAROTENOID_BOUND.md` §2, computed at a different growth regime and
  with the ec model's own protein-pool calibration, rather than an LP run under Elizondo's
  measured uptake bounds. The conclusion — β-carotene ranks last by three orders of
  magnitude — is far too large to be overturned by that mismatch, but the numbers are not
  the ones the implementation should record.
- **Elizondo's chemostat may not have been carbon-limited as the paper describes it, and two
  independent checks agree on this.** Taking β-car2 at µ = 0.101, OD600 0.812 × 0.533
  gDCW/OD gives 0.433 gDCW/L; q_glucose 6.68 mmol/gDCW/h then implies 2.89 mmol/L/h consumed,
  which at D = 0.101 /h is **5.2 g/L of the 12 g/L fed — about 57 % of the glucose
  unaccounted for**. Two yields follow and they answer different questions: on *consumed*
  glucose, Y = µ/(q_S · 0.180) = **0.084 g/g**; on *fed* glucose, X/S_in = **0.036 g/g**.
  Both are far below the ~0.5 g/g of a genuinely glucose-limited aerobic chemostat, and
  ethanol at 6.40 mmol/gDCW/h at µ = 0.1 is itself inconsistent with carbon limitation. No
  carbon balance is reported in the paper. This does not affect the usability of the measured
  rates, which stand on their own and are what the ordering consumes — but the regime must
  **not** be described as glucose-limited in this repo's prose, and the derived CO2 of §4b
  inherits the same unclosed input, so it should be recorded as derived-under-an-assumption.
- **No CO2 or O2 measurement exists for any carotenoid *S. cerevisiae* dataset found.**
  Scalcinati 2012 has both, on α-santalene. Every carotenoid paper lacks them.
- **No 13C-MFA of a carotenoid-producing *S. cerevisiae* exists**, as far as ~15 Europe PMC
  full-text query variants can establish. Absence of evidence over open-access full text
  only; paywalled full text was not searchable.
</content>
