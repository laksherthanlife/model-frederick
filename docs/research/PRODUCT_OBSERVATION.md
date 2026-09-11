# Reading the product on the plate reader

The observation layer read biomass (OD600) and reporter (mCitrine RFU) and did not read
the product. That left the last link of the chain unmeasurable on the instrument this lab
already owns, and it is the link the twin exists to forecast.

Beta-carotene is a pigment. It absorbs near 450 nm, lycopene near 470, phytoene near 285.
So the same reader, in the same well, on the same run, can give all three quantities:
scatter at 600 nm for biomass, fluorescence for the stress reporters, absorbance at 450 nm
for titre. `src/ystwin/observation.py` now does that.

**This channel is not free, and the thing that is not free is the cells.** A whole-cell
reading at 450 nm is pigment absorbance *plus* light scattering, and at the contents this
organism reaches the two terms are the same order of magnitude. Attributing the whole
reading to pigment is the same class of error as the naive-versus-dilution-corrected fold
already recorded here: a number that looks like a titre and is mostly turbidity.

## 1. The model

Whole cells:

```
A(lambda) = blank(lambda)
          + scatter_per_od600 * OD600_true
          + flattening * eps(lambda) * (product * biomass / 1000) * path_length
```

clipped at `detector_max` when one is declared. Units: `product` in mmol/gDCW — the basis
`kinetic/carotenoid.py` carries its pools in — `biomass` in gDCW/L, `eps` in M⁻¹cm⁻¹,
`path_length` in cm. `OD600_true` is `biomass / gdcw_per_od`, so the dry-weight factor is
required here as it is in `observe_od`.

Extracted sample, `observe_absorbance_extract`:

```
A(lambda) = blank(lambda) + eps(lambda) * (product * harvested / volume / 1000) * path_length
```

No cells, so no scattering term. No packaging, so the solution coefficient is the right
one. This path consults neither `scatter_per_od600` nor `flattening`, and it is the
reference the whole-cell path is calibrated against.

`product_from_absorbance` inverts the whole-cell form. That is the direction the instrument
is used in, and it is where the scattering term earns its keep: subtracting a measured
scattering floor is the difference between a titre and a turbidity.

## 2. Scattering: modelled explicitly, and refused without a measurement

The design follows Myers, Curtis & Curtis 2013 (PMID 24499615), which is directly on this
problem — separating chromophore absorbance from cell scattering in an optical-density
reading. Their statement of the split is the one this module implements: OD at a *robust*
wavelength "is primarily the result of light scattering and does not vary with culture
conditions", whereas OD at a *sensitive* wavelength "is additionally dependent on light
absorption by the organism's pigments", and their remedy is correlation against off-peak
light attenuation.

So `scatter_per_od600` is **absorbance at 450 nm per unit linearised OD600, contributed by
cells alone**. It is measured, not derived:

> Grow the isogenic **non-producing** strain. Run a dilution series. Read every well at
> both 600 nm and 450 nm. Linearise the 600 nm readings through `calib/od.py`. The slope of
> A450 on linearised OD600 is the ratio. One plate, one strain, and it is then the same for
> every producing well on the reader.

What the module does **not** do is extrapolate the ratio from a wavelength power law.
Yeast cells are several microns across, far larger than 450 nm, so the scattering is deep
in the Mie regime where the wavelength dependence is weak and not a clean exponent. An
asserted `lambda^-n` would be exactly the plausible substitute number this repository
refuses. `observe_absorbance` raises with the missing measurement named, and the message
points at `observe_absorbance_extract` as the route that needs no such number — a refusal
that does not say what to do instead gets worked around.

Two tests pin this: `test_a_well_with_no_product_still_reads_well_above_the_blank`, and
`test_ignoring_the_scattering_floor_would_have_inflated_the_titre`, which shows the naive
inversion overstating content by more than twofold at ordinary numbers.

## 3. Packaging: the second confound, also refused by default

Pigment inside cells is not pigment in solution. Confining a chromophore to particles
flattens its apparent absorption spectrum relative to the same quantity dissolved — the
sieve or package effect, Duysens 1956,
doi [10.1016/0006-3002(56)90380-8](https://doi.org/10.1016/0006-3002(56)90380-8). A
solution extinction coefficient applied to intact cells therefore over-reads, by a factor
nobody in this project has measured.

`PigmentOptics.flattening` carries it, in (0, 1]. `None` refuses. `1.0` is permitted and is
an explicit claim that there is no packaging effect, which for a membrane-localised
carotenoid is a strong claim and has to be typed deliberately — the same idiom
`inner_filter_coeff=0.0` already uses. A value above 1 is rejected at construction:
packaging can only reduce apparent absorbance, so a fitted value above 1 is a scattering
ratio that is too small, not an optical effect.

Measuring it needs no new equipment: read one culture whole-cell, extract it, read the
extract, and take the ratio at matched pigment quantity.

> **Note on the Crossref record.** Crossref returns the Duysens title as "The flattering of
> the absorption spectrum of suspensions…" and the author as "Duyens". Both are typos in
> the registered metadata, not in the citation. The journal, year and subject are right, so
> the refreshed `doi_titles.csv` will carry the typo; that is the registry's text, kept
> verbatim rather than corrected, so the table stays a record of what resolves.

## 4. The extinction coefficient, with its solvent

A carotenoid extinction coefficient is not a property of the carotenoid. It is a property
of the carotenoid **in a solvent**, at a **wavelength**: beta-carotene's maximum sits near
450 nm in petroleum ether and near 466 nm in chloroform. That solvent dependence is the
subject of Craft & Soares 1992, doi [10.1021/jf00015a013](https://doi.org/10.1021/jf00015a013),
whose title is "Relative solubility, stability, and absorptivity of lutein and
β-carotene in organic solvents". So `PigmentExtinction` makes `wavelength_nm` and `solvent`
required fields and refuses construction without them.

The value encoded:

| | |
| --- | --- |
| Quantity | Molar extinction coefficient of **beta-carotene** |
| Solvent | **petroleum ether** |
| Wavelength | **450 nm** (the maximum in that solvent) |
| Tabulated form | A(1%, 1 cm) = **2592** |
| Molar form | 2592 × 536.87 / 10 = **1.3916 × 10⁵ M⁻¹cm⁻¹** |
| Primary tabulation | Britton, Liaaen-Jensen & Pfander, *Carotenoids: Handbook* (2004), doi [10.1007/978-3-0348-7836-4](https://doi.org/10.1007/978-3-0348-7836-4) |
| Read directly, with solvent and wavelength stated | doi [10.1038/s41598-026-45956-6](https://doi.org/10.1038/s41598-026-45956-6): "E 1%1cm is the extinction coefficient of β-carotene in petroleum ether (2592)", with A read at 450 nm |

**What was verified and how.** The handbook DOI was resolved through Crossref and returns
*Carotenoids: Handbook*, Birkhäuser Basel 2004, editors Britton, Liaaen-Jensen and Pfander.
The handbook itself is not open, so the numeric value was not read from it; it was read
from the open-access record above, which states the value, the solvent and the wavelength
together. The conversion to molar units is arithmetic on the molecular weight 536.87 g/mol,
the same figure `EXTERNAL_CAROTENOID_BOUND.md` and `scripts/parked/run_d1.py` use, and the
test pins the derivation rather than the number.

**Corroboration, and its defect.** Three carotenoid-engineering papers quote the same
handbook as **138,900 M⁻¹cm⁻¹ at 450 nm** — PMC4510654, PMC6516660, PMC10148934 — which is
0.2% from the value above. None of the three names a solvent. That is precisely the failure
`PigmentExtinction` makes impossible in this repository, and it is recorded here rather
than smoothed over.

**What is still open.** Whether the handbook's 450 nm entry is petroleum ether or hexane is
not settled by what was read: one open-access source gives 2592 at 453 nm and attributes it
to "ethanol or petroleum ether". Since 453 nm is the hexane maximum and 450 nm the
petroleum-ether maximum, the wobble is a wavelength-solvent pairing, not a disagreement
about the number. The module turns that ambiguity into an enforced check instead of a
silent assumption: `PigmentOptics` refuses to pair a coefficient with a filter at a
different wavelength.

## 5. The inner filter, and its sign

mCitrine is excited near 516 nm and emits near 529. Beta-carotene absorbs across 400–500 nm
with a shoulder past that. In a strain that both produces carotenoid and carries a
fluorescent reporter, the product eats the reporter's excitation light. This is guaranteed,
not hypothetical, and `ReporterOptics.inner_filter_coeff` already anticipated it.

**The sign convention: product REDUCES apparent fluorescence.** `observe_rfu` multiplies
the per-cell signal by `exp(-inner_filter_coeff * carotenoid)`, so RFU falls as titre rises
at fixed promoter activity. `correct_inner_filter` divides by the same factor, so the
correction can only move a reading **up**.

Getting that backwards is the most expensive single error available in this module. The
inferred promoter activity would be too high in exactly the producing strains the twin
exists to describe, and the error would grow as the strain improves — an artefact that
reads as product-induced stress and would survive every shuffled-channel control, because
the coupling is real signal.

Pinned by three tests in `tests/test_product_absorbance.py::TestTheInnerFilterSign`:
product reduces the reading; the attenuation deepens monotonically with titre; the
correction never lowers a reading. The refusal path is pinned too — with
`inner_filter_coeff=None` the correction raises rather than treating an uncalibrated filter
as an absent one.

**Making the coupling one quantity instead of two.**
`inner_filter_coeff_from_extinction` computes the coefficient the product channel implies:
transmitted fraction is `10^-A`, so with `A = flattening * eps * (q * X / 1000) * l` the
natural-log coefficient is `ln(10) * flattening * eps * X * l / 1000` per mmol/gDCW. It
refuses an extinction coefficient measured at any wavelength other than the excitation
wavelength — reusing the 450 nm peak for a 516 nm excitation would overstate the coupling
several-fold and therefore over-correct, which is the sign error again, only quieter.

The derived value is a **lower bound** on the measured spike-in coefficient, because it
accounts for the excitation path only and emission at 529 nm is reabsorbed as well. A
measured coefficient below this bound is a result, not a rounding difference.

The coefficient carries a biomass, because `attenuation()` is a function of intracellular
concentration while absorbance is a function of well concentration. It is therefore valid
near the density it was calibrated at. That is a real limitation of the existing signature
and is left as it stands rather than changed under another lane's tests.

## 6. Dynamic range: this channel saturates long before OD600 does

Worked at the content `EXTERNAL_CAROTENOID_BOUND.md` computes on the D1 frontier — 32.84
mg/gDCW, which is 0.0612 mmol/gDCW — with `flattening = 1` and a 0.5 cm path:

| biomass | pigment absorbance at 450 nm, solution basis |
| ---: | ---: |
| 1 gDCW/L | **4.26 AU** |
| 0.5 gDCW/L | 2.13 AU |
| 0.1 gDCW/L | 0.43 AU |

The whole-cell figure is lower by the flattening factor, and the scattering term adds to
it. Even so: a good producer at ordinary plate densities runs the 450 nm channel past any
reader's linear range while the 600 nm channel is still comfortable. Two consequences.

1. **The product channel needs its own dilution and its own linear-range gate.** G1's
   `od_linear_max` is a threshold on the 600 nm channel; it says nothing about 450 nm.
2. **A clipped reading must not be inverted.** `product_from_absorbance` raises when the
   reading is at or above a declared `detector_max`, because a clipped reading carries no
   concentration information and would return the same content for every well past the
   ceiling.

Detector compression *below* the ceiling is not modelled. `ODCalibration.saturation_k` is
fitted at 600 nm and reusing it at 450 nm would be an assumption, so the module applies a
hard ceiling only and leaves the curvature to a dilution series at 450 nm. Plate-reader
absorbance is in any case not cuvette absorbance — the path length is set by fill volume
and the geometry differs from a cuvette (Stevenson et al. 2016, PMID 27958314) — which is
why `path_length_cm` is a per-experiment argument rather than a plate constant.

## 7. What this channel cannot see

**Phytoene is colourless at 450 nm.** It absorbs near 285 nm. A strain blocked at crtI
accumulates phytoene and reads as *zero product* on this channel while carrying a full
carbon flux into the branch. The absorbance channel measures coloured carotenoid, not
pathway flux, and the two separate exactly where `kinetic/carotenoid.py` says the branch
can be limited. HPLC remains the only way to resolve phytoene, lycopene and beta-carotene,
and it is one run for all three.

**Lycopene reads at 470 nm and overlaps beta-carotene at 450.** A strain accumulating both
gives a single 450 nm number that is a mixture. Nothing here deconvolves it; the module
refuses the wrong-wavelength pairing, which is as far as one filter can go.

## 8. Claim tier

Under [`../CLAIM_BOUNDARY.md`](../CLAIM_BOUNDARY.md):

| Item | Tier |
| --- | --- |
| The absorbance model and its inverse | **Tier 1**, self-consistent in simulation. Round-trip and Beer-Lambert arithmetic, correctly carried out on a model we adopted. Nothing here has met a plate. |
| eps = 1.3916 × 10⁵ M⁻¹cm⁻¹, petroleum ether, 450 nm | **Literature, with its solvent and wavelength verified.** Not measured on this reader, and a solution coefficient is not a whole-cell coefficient. |
| `scatter_per_od600` | **Not asserted at all.** No default exists; the module refuses. This is the single most valuable number to measure for the channel. |
| `flattening` | **Not asserted at all.** Same treatment. |
| `path_length_cm`, `blank`, `detector_max` | Per-experiment arguments with no defaults beyond a zero blank. |
| The inner-filter sign | **Structural**, and tested. It follows from the direction of Beer-Lambert absorption, not from any fitted number. |

No constant in this file has been given a plausible value in place of a missing
measurement. Three quantities the channel needs are missing, and all three refuse.

## 9. What it would cost to make this a measurement

In cost order.

1. **One dilution series of the non-producing strain, read at 600 and 450 nm.** Gives
   `scatter_per_od600` and, from the same plate, the 450 nm linear range. No new strain, no
   new reagent.
2. **One whole-cell reading beside one extract of the same culture.** Gives `flattening`.
   Needs a producing strain, which does not yet exist here.
3. **A purified-carotenoid spike-in into non-producing wells, read in the mCitrine
   channel.** Gives `inner_filter_coeff` empirically, and checks it against the lower bound
   `inner_filter_coeff_from_extinction` computes.

Item 1 is runnable now and is the one that turns a refusal into a channel.

## Provenance

Written 2026-08-26 against `src/ystwin/observation.py` and
`tests/test_product_absorbance.py`. Every identifier here was resolved before it was
written down: the two DOIs and the handbook through Crossref, PMID 24499615 and PMID
27958314 through PubMed, and the three PMC records by reading the sentence that carries the
coefficient. `scripts/refresh_citations.py` vendors what each one resolves to.
