# The van Hoek 1998 chemostat reference

## What changed

`data/physiology/chemostatData_VanHoek1998.tsv` held 4 rows and 4 columns. It now holds
**10 rows and 10 columns** — the whole of Table 1 of the source paper. The four rows that
were already there agree with the source exactly, and their original four fields are
byte-identical; the extension only appends.

The columns that were missing are the ones `bridge/regulation.py` needs. Ethanol and
glycerol carry the overflow carbon in yeast, and without them a task-efficiency ordering
over this dataset can only see biomass and CO2.

## Source

| | |
| --- | --- |
| Paper | van Hoek P, van Dijken JP, Pronk JT |
| Title, as returned by PubMed | **"Effect of specific growth rate on fermentative capacity of baker's yeast."** |
| Journal | *Applied and Environmental Microbiology* 64(11):4226–4233, November 1998 |
| PMID | **9797269** |
| DOI | **10.1128/aem.64.11.4226-4233.1998** |
| Full text | Free to read at PMC106631 |
| Table | Table 1, "Cell yields, metabolic fluxes, and carbon recovery as a function of the dilution rate in aerobic, glucose-limited chemostat cultures of *S. cerevisiae* DS28911" |

Verified by resolving PMID 9797269 through NCBI E-utilities `esummary` and quoting the
returned title back. The DOI is recorded in the lowercase form the vendored citation table
carries; PubMed returns the same string with `AEM` capitalised, and both resolve.

**Access.** Free to read, **not open access.** Europe PMC reports `isOpenAccess: N` and no
licence, and PMC refuses full text over `efetch` ("the publisher of this article does not
allow downloading of the full text in XML form"). The rendered PMC article page does serve
Table 1 as real HTML table markup, which is where these numbers come from — not from a
figure, and not from a paywalled PDF. Nothing here was digitised off a plot.

**Strain and regime.** *S. cerevisiae* DS28911, an industrial baker's yeast. Aerobic
glucose-limited chemostat, defined mineral medium, 1.0 L working volume, 30 °C, 800 rpm,
air at 0.5 L/min, dissolved oxygen above 60 % of air saturation, no detectable
oscillations. Batch µ_max on glucose is quoted as 0.42 /h, and the paper cites that to
another work rather than measuring it here.

**This is not CBS 8066.** `docs/research/CALIBRATION_DATA.md` §5 calls it "Van Hoek's
CBS 8066" twice. The strain is DS28911. CBS 8066 appears in this paper once, in a
reference title. That file is not this lane's to edit.

## The table

Every value below is from Table 1. Fluxes are **mmol/gDW/h** and no conversion was needed —
footnote *a* of the table reads "Fluxes (q) are expressed as millimoles per gram of dry
yeast biomass per hour", and the Fig. 1 legend repeats the same units for q_ethanol, q_O2
and q_CO2 independently. `Drate` is /h. `BiomassYield_g_per_g` is grams dry biomass per gram
glucose. `CarbonRecovery_pct` is a percentage.

| Drate | Glucose | O2 | CO2 | Ethanol | Glycerol | Acetate | Pyruvate | Yield | C recovery |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.025 | 0.30 | 0.80 | 0.80 | 0 | 0 | 0 | 0 | 0.45 | 98.9 |
| 0.05 | 0.60 | 1.30 | 1.40 | 0 | 0 | 0 | 0 | 0.47 | 95.0 |
| 0.10 | 1.10 | 2.50 | 2.70 | 0 | 0 | 0 | 0 | 0.48 | 96.0 |
| 0.15 | 1.70 | 3.90 | 4.20 | 0 | 0 | 0 | 0 | 0.49 | 102.4 |
| 0.20 | 2.30 | 5.30 | 5.70 | 0 | 0 | 0 | 0 | 0.48 | 100.9 |
| 0.25 | 2.80 | 7.00 | 7.50 | 0 | 0 | 0 | 0 | 0.48 | 102.6 |
| 0.28 | 3.40 | 7.40 | 8.00 | 0.11 | 0 | 0.08 | 0.01 | 0.46 | 97.0 |
| 0.30 | 4.50 | 6.10 | 8.80 | 2.30 | 0 | 0.41 | 0.01 | 0.37 | 99.1 |
| 0.35 | 8.60 | 5.10 | 14.90 | 9.50 | 0.05 | 0.62 | 0.03 | 0.23 | 99.4 |
| 0.40 | 11.10 | 3.70 | 18.90 | 13.90 | 0.15 | 0.60 | 0.05 | 0.20 | 97.9 |

**Every zero is a non-detect, not a measured zero.** Table 1 footnote *c* reads "0, below
detection limit". The TSV writes them as `0.00` so the file parses as numeric; treat them as
"below the detection limit of the assay named below", which is a bound, not a point value.
Stated detection limits: glucose ca. 5 µM by glucose oxidase kit, acetic acid ca. 0.2 mM by
enzymic kit. Ethanol, glycerol and pyruvate were HPLC (Aminex HPX-87H, 60 °C, 5 mM H2SO4 at
0.6 mL/min; pyruvate by UV at 214 nm, ethanol and glycerol by refractive index) and the
paper states no detection limit for them.

**The critical dilution rate is 0.28 /h.** Below it, metabolism is fully respiratory: no
fermentation products, RQ within 8 % of unity, yield flat at 0.45–0.49 g/g. Ethanol is
first detectable at D = 0.28 and rises to 13.9 by D = 0.40, where the yield has collapsed
to 0.20. Glycerol appears only at D ≥ 0.35. `tests/test_physiology_validation.py` asserts
the critical rate is "near 0.28 /h" on this paper's authority, and Table 1 supports that
exactly.

**Oxygen uptake is not monotonic.** q_O2 peaks at 7.4 at D = 0.28 and then *falls* to 3.7 at
D = 0.40 as the cells switch to fermentation. Any code that treats oxygen uptake as
increasing with growth rate is wrong on this dataset above D_crit.

### The transcription is checked by arithmetic, not by eye

The prior literature pass (`CALIBRATION_DATA.md` §8) recorded that these values had been
read off the PMC page by an automated fetch and asked for a hand check before vendoring.
Rather than re-read them, the check here is that the table's own two derived columns are
**recomputed from the transcribed fluxes and compared against what the paper printed**. A
single mis-transcribed digit in glucose, CO2 or ethanol breaks the carbon balance
visibly — that is the point of doing it this way.

Carbon in is `q_glucose × 6`. Carbon out is `q_CO2 × 1 + q_ethanol × 2 + q_acetate × 2 +
q_pyruvate × 3 + q_glycerol × 3`, plus biomass carbon `D × 0.48 / 12.011 × 1000` C-mmol/gDW/h,
using the 48 % biomass carbon content that Table 1 footnote *b* names. Yield is
`D / (q_glucose × 0.180156)`.

| Drate | C recovery recomputed | printed | Δ | yield recomputed | printed | Δ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.025 | 99.9 | 98.9 | +1.0 | 0.463 | 0.45 | +2.8 % |
| 0.05 | 94.4 | 95.0 | −0.6 | 0.463 | 0.47 | −1.6 % |
| 0.10 | 101.5 | 96.0 | +5.5 | 0.505 | 0.48 | +5.1 % |
| 0.15 | 99.9 | 102.4 | −2.5 | 0.490 | 0.49 | −0.0 % |
| 0.20 | 99.2 | 100.9 | −1.7 | 0.483 | 0.48 | +0.6 % |
| 0.25 | 104.1 | 102.6 | +1.5 | 0.496 | 0.48 | +3.3 % |
| 0.28 | 96.1 | 97.0 | −0.9 | 0.457 | 0.46 | −0.6 % |
| 0.30 | 97.2 | 99.1 | −1.9 | 0.370 | 0.37 | 0.0 % |
| 0.35 | 95.7 | 99.4 | −3.7 | 0.226 | 0.23 | −1.8 % |
| 0.40 | 96.8 | 97.9 | −1.1 | 0.200 | 0.20 | 0.0 % |

Carbon recovery agrees within 5.5 percentage points worst case and under 2 on seven of ten
rows; yield agrees within 5.1 % worst case and under 1 % on six of ten. Both residuals are
what rounding the q columns to two significant figures produces — at D = 0.10, a printed
yield of 0.48 implies q_glucose = 1.156, which prints as 1.1 or 1.2 depending on the true
value. The transcription is sound. It does not recover the significant figures the authors
rounded away, and no arithmetic here can.

## What each column unblocks

`bridge/regulation.py` implements the COSMIC-dFBA §4.3 lexicographic task-efficiency
ordering (PMID 38387677, doi 10.1016/j.ymben.2024.02.012). Each task's efficiency is
measured flux ÷ its individually maximised FBA flux, and `MetabolicTask` refuses to yield a
flux it was not given — an unmeasured task raises `MissingMeasurement`, and a measured one
without a `source` string raises too. So each column is one task that can now be ranked
instead of refused.

| column | task | what it unblocks |
| --- | --- | --- |
| `Drate` | biomass | µ = D. Already present. |
| `GlucoseUptake` | — | the uptake bound the LP runs under, not a task. Already present. |
| `O2uptake` | — | the respiratory ceiling that makes overflow appear at all. Already present. |
| `CO2production` | CO2 | already present, and the only secretion task the truncated file could rank. |
| **`Ethanol`** | ethanol | **the largest secretion task above D_crit** — 13.9 mmol/gDW/h at D = 0.40, larger than glucose uptake. Nothing could rank it before. |
| **`Glycerol`** | glycerol | the redox-overflow task. Non-zero only at D ≥ 0.35, so it also fixes where the ordering has 3 rankable secretion tasks and where it has 2. |
| **`Acetate`** | acetate | third secretion task, non-zero from D = 0.28. |
| **`Pyruvate`** | pyruvate | fourth, and small — 0.01 to 0.05. Useful as the low-efficiency end of the ordering. |
| `BiomassYield_g_per_g` | — | not a flux. Consistency check on glucose uptake against µ, as used above. |
| `CarbonRecovery_pct` | — | not a flux. The paper's own closure, 95–103 %, which is what licenses deriving a missing term by difference on this dataset. |

**What this dataset is for, and what it is not for.** `CALIBRATION_DATA.md` §5 settles this
and the extension does not reopen it: the non-product task efficiencies are **not**
strain-independent in the strains this project cares about, so this table is the
**respiratory contrast**, not the calibration set. DS28911 secretes no ethanol at all below
D = 0.28, where Elizondo's β-carotene strains secrete 6.40 mmol/gDW/h at µ = 0.101 and take
up glucose six times faster. Ranking tasks off this table and applying the order to a
carotenoid strain would invert the ethanol rank outright.

## Disagreement with `fba/physiology.py`

**`REFERENCE_AEROBIC_BATCH` cites this paper for five numbers that are not in it.** Do not
edit `src/ystwin/fba/physiology.py` on the strength of this section — it is another lane's
file, and this is a report.

```
REFERENCE_AEROBIC_BATCH = PhysiologyReference(
    name="aerobic glucose batch",
    growth_rate=0.40, glucose_uptake=21.3, oxygen_uptake=7.8,
    ethanol_secretion=27.4, co2_secretion=20.4,
    source="van Hoek, van Dijken & Pronk 1998, Appl Environ Microbiol 64:4226",
)
```

The strings `21.3`, `27.4` and `20.4` **do not occur anywhere in PMID 9797269**, in the body
text or in Table 1. The paper contains exactly one table and three figures; there is no
batch-culture flux table in it, and no batch culture was run — every steady state is a
glucose-limited chemostat. `tests/test_physiology.py` additionally attributes the reference
to "CEN.PK113-7D, aerobic glucose-excess batch"; this paper's strain is DS28911 and
CEN.PK113-7D is not mentioned.

The nearest thing the cited source contains is its highest dilution rate, and it disagrees
on every measured flux:

| quantity | `REFERENCE_AEROBIC_BATCH` | van Hoek Table 1 at D = 0.40 /h | disagreement |
| --- | ---: | ---: | ---: |
| growth rate, /h | 0.40 | 0.40 | matches |
| glucose uptake | 21.3 | **11.1** | 1.92× |
| oxygen uptake | 7.8 | **3.7** | 2.11× |
| ethanol secretion | 27.4 | **13.9** | 1.97× |
| CO2 production | 20.4 | **18.9** | 1.08× |

Three of the four fluxes are close to **twice** the source's value at the same growth rate,
which reads more like a different culture regime — glucose-excess batch is a real and
different phenotype — than like transcription drift. The likely reading is that the
constants are genuine aerobic batch measurements of CEN.PK113-7D from **some other source**,
and that only the citation is wrong. That other source was not found.

The other van Hoek 1998 paper was checked and ruled out: PMID 9603825, "Effects of pyruvate
decarboxylase overproduction on flux distribution at the pyruvate branch point in
*Saccharomyces cerevisiae*", AEM 64(6):2133–2140. It does use CEN.PK113-7D, but its two
tables carry µ_max and glycolytic enzyme activities only — no extracellular fluxes — and
`21.3`, `27.4` and `20.4` do not occur in it either. Those are the only two van Hoek 1998
records in PubMed.

`src/ystwin/generator/literature.py` already carries 21.3 and 27.4 as
`verified=False`, and `docs/SOURCE_AUDIT.md` §2.2 already classifies the `literature.LIBRARY`
van Hoek entries as "cited vaguely". This section upgrades that: the attribution is not
vague, it is **wrong**, and `fba/physiology.py` states the same numbers with no
`verified` flag at all.

### A second lead, weaker

`src/ystwin/bridge/physiology_bridge.py:30` reads
`_GLUCOSE_QMAX = 22.0 # mmol/gDW/h, van Hoek et al. 1998 aerobic batch`.

**22.0 mmol/gDW/h does occur in PMID 9797269, and it is not a glucose uptake rate.** It is
the *fermentative capacity* — an offline assay in which chemostat-grown cells are moved to
anaerobic conditions under CO2 with 2 % glucose, and the resulting **ethanol** production
rate is measured. The paper: "the fermentative capacity showed only a small further
increase, up to 22.0 mmol of ethanol · g of dry yeast biomass⁻¹ · h⁻¹ at D = 0.40 h⁻¹". The
highest in-situ glucose uptake this paper reports is 11.1.

So a glucose Vmax of 22.0 attributed to "van Hoek 1998 aerobic batch" is plausibly an
anaerobic ethanol capacity read as an aerobic glucose uptake. This is inference from a
matching number and a matching citation, not proof — 22.0 is also almost exactly 2 × 11.1,
and the constant may have an unrelated origin. Worth one check by whoever owns that file.

## What could not be obtained

- **Error bars on Table 1.** The table prints single values. Fig. 1 shows means and standard
  deviations of duplicate assays at different timepoints within the same steady state, but
  as a plot. The replication is therefore analytical, within one chemostat run per dilution
  rate, and the TSV carries no uncertainty column because the source publishes none in
  numeric form. Digitising Fig. 1 would give the spread and is not done here.
- **Significant figures beyond what was printed.** q values are given to two significant
  figures. The 2–5 % residuals in the arithmetic check above are that rounding and cannot be
  recovered.
- **Biomass concentration, residual glucose, and protein content per row.** Protein content
  is Fig. 1C only. Residual glucose is not tabulated. Neither is needed for the task
  ordering.
- **The true source of `REFERENCE_AEROBIC_BATCH`'s five constants.** Both van Hoek 1998
  papers are ruled out. Unresolved.
- **A machine-readable copy.** PMC serves no XML full text for this publisher and Europe
  PMC's repository endpoint returns 403, so the table was parsed out of the rendered PMC
  page's HTML table markup. That markup is structured, not OCR, and the arithmetic check
  above is the guard against a parse error.

## Provenance

The table is published measurement on an industrial baker's yeast strain, by another
laboratory, in 1998. It is not this lab's data and nothing in it is a claim of this
repository. Under `docs/CLAIM_BOUNDARY.md` it is external literature: a model result
agreeing with it is held-out agreement at best, and it never becomes evidence about the
strains in `data/plates/`.
