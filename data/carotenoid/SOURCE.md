# Carotenoid calibration data — provenance

Two files, one source: **Elizondo & Saa 2025**, PMID `40891387`, doi
`10.1021/acssynbio.5c00256`, ACS Synth Biol 14(9):3457–3472, open access at PMC12455641.
Title as returned by NCBI E-utilities `esummary`:

> Complex Kinetic Models Predict β-Carotene Production and Reveal Flux Limitations in
> Recombinant Saccharomyces cerevisiae Strains.

`docs/research/CALIBRATION_DATA.md` establishes why this is the only sufficient dataset.
This file records where the numbers physically came from and what they mean.

## Where the numbers came from

Not from the article PDF and not from the supplementary PDF. From the **authors' own
released data repository**, `github.com/SysBioengLab/BcarGRASP`, which the article links
from its Methods. That repository holds the pre-processing output the article's tables were
rendered from, so the values here are at full precision and did not pass through a
transcription step.

| this file | upstream path in `SysBioengLab/BcarGRASP` | the article's rendering |
| --- | --- | --- |
| `elizondo2025_steady_states.tsv`, `q_*` columns | `data_pre-processing/fluxes/output/FluxesBiomassStatistics.xlsx` | Table S4 (glucose, ethanol, acetate, glycerol); Table 2 (lycopene, β-carotene) |
| `elizondo2025_steady_states.tsv`, `mu_*` | `.../DilutionRatesStatistics.xlsx` | Table S8 |
| `elizondo2025_steady_states.tsv`, `od600*` | `.../BiomassMeasurementsStatistics.xlsx` | Table S9 |
| `elizondo2025_steady_states.tsv`, `od_to_gdcw*` | `.../BiomassConversionStatistics.xlsx` | Table S7 |
| `elizondo2025_relative_mrna.tsv` | `data_pre-processing/transcripts/output/ConditionsRelativeLimits.xlsx`, sheet `4A` | Figure 5 heatmap |

**Cross-check performed.** Table 2 was re-extracted independently from the JATS XML that
NCBI E-utilities returns for PMC12455641 and compared against the repository file. All
twelve product rates agree: β-car2 94.4 / 151.3, β-car3 472.0 / 227.1, β-car4 760.7 / 185.3
nmol/gDCW/h at µ = 0.101, and 84.3 / 184.2, 126.5 / 294.8, 392.9 / 453.6 at µ = 0.254
(lycopene / β-carotene). Table S4 could not be cross-checked the same way — the
supplementary PDF is served by ACS and returns 403, and PMC does not mirror it — but the
repository file is the upstream of that table, not a copy of it.

**Column order in the mRNA sheet was verified, not assumed.** `ConditionsRelativeLimits.xlsx`
labels its columns `Var2_1 … Var2_18` with no condition names: eighteen columns are six
conditions × (lower 95, mean, upper 95). The ordering `2D01, 2D025, 3D01, 3D025, 4D01,
4D025` is pinned by the article's own text, which quotes β-car3 at low growth rate as
CrtE ≈ 61 %, CrtI ≈ 89 %, CrtYB ≈ 68 % and β-car2 as ≈ 27 %, ≈ 51 %, ≈ 33 %. Triple 3 reads
0.6147 / 0.8907 / 0.6848 and triple 1 reads 0.2710 / 0.5106 / 0.3276. Triple 5 is exactly
1.000 for every gene, which is β-car4 at low growth rate, the normalising condition.

## Units and meaning

`elizondo2025_steady_states.tsv`

| column | unit | note |
| --- | --- | --- |
| `mu_per_h` | 1/h | the dilution rate, measured, not the setpoint. µ = D at steady state |
| `q_glucose` | mmol/gDCW/h | **negative is uptake**, as the authors write it |
| `q_ethanol`, `q_acetate`, `q_glycerol` | mmol/gDCW/h | positive is secretion |
| `q_lycopene`, `q_betacarotene` | mmol/gDCW/h | Table 2 reports these in nmol/gDCW/h; this file keeps the repository's mmol so every flux column shares one unit. Multiply by 1e6 to recover the article's numbers |
| `od600`, `od_to_gdcw_per_l` | — , (gDCW/L)/OD600 | the conversion factor is **per strain**, measured, and must not be borrowed across strains |
| `biomass_g_per_l` | g/L | derived here as `od600 × od_to_gdcw_per_l`; the only computed column |
| `*_lo95`, `*_hi95` | same as parent | 95 % confidence interval as released |

`elizondo2025_relative_mrna.tsv` — RT-qPCR, **relative** expression, dimensionless, every
condition normalised to β-car4 at µ = 0.101.

## What the product columns are, and are not

`q_lycopene` and `q_betacarotene` are **accumulation** rates, not reaction fluxes.
Carotenoids are not secreted; they sit in the membrane and leave the vessel inside the
biomass. At chemostat steady state that makes each one the growth dilution of an
intracellular pool:

    q_lycopene      = µ · [lycopene]
    q_betacarotene  = µ · [beta-carotene]

so the intracellular **content** in mmol/gDCW is `q / µ`, and it is a measured quantity, not
a fitted one. The cyclase reaction flux equals `q_betacarotene`; the desaturase flux equals
`q_lycopene + q_betacarotene`. That identity is what makes this dataset usable for the
kinetic layer at all, and it is derived in `docs/research/KINETIC_FIT.md` §2.

## Two caveats that travel with the data

**The vessel was probably not carbon-limited**, whatever the article calls it. Taking β-car2
at µ = 0.101, biomass 0.4326 g/L × q_glucose 6.683 mmol/gDCW/h = 2.89 mmol/L/h consumed,
against 12 g/L × 0.101 /h = 6.73 mmol/L/h fed. About 57 % of the glucose is unaccounted for,
and a genuinely glucose-limited aerobic chemostat does not secrete 6.4 mmol/gDCW/h of
ethanol at µ = 0.1. The measured rates stand on their own — they are what the fit consumes —
but nothing here may be described as glucose-limited. Recorded at greater length in
`docs/research/CALIBRATION_DATA.md` §8.

**Relative mRNA is not enzyme per gDCW.** RT-qPCR relative expression is normalised to
reference genes, and total RNA per gDCW itself rises with growth rate in yeast, so the
across-growth-rate ratios in `elizondo2025_relative_mrna.tsv` carry a confound that the
across-strain ratios at one growth rate do not. `KINETIC_FIT.md` §4 uses it only where that
confound is absent, and says so where it is not.
