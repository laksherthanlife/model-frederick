# Triacetic acid lactone in *Kluyveromyces marxianus* — provenance, and why no spec ships

**One source:** Lang X, Besada-Lombana PB, Li M, Da Silva NA, Wheeldon I. 2020.
"Developing a broad-range promoter set for metabolic engineering in the thermotolerant
yeast *Kluyveromyces marxianus*." *Metab Eng Commun* 11:e00145. PMID `32995271`,
PMC `PMC7508702`, doi `10.1016/j.mec.2020.e00145`. Open access, CC-BY.

Title and author list as returned by NCBI E-utilities `esummary` for PMID 32995271, not
taken on anyone's word. The full text was read as the Europe PMC JATS XML
(`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7508702/fullTextXML`, 120,315 bytes)
and the supplementary material as `mmc1.docx` from the Europe PMC `supplementaryFiles`
endpoint. Every number below came out of one of those two fetches.

**`docs/SECOND_DATASET_HUNT.md` ranks this dataset first of four for fit to the pathway
chain. This directory is the record of the attempt to use it, and the attempt failed.**
It ships the two axes that *are* numeric so nobody repeats the retrieval, and it states
exactly which column is missing and where it physically is.

## The verdict, first

**A `data/pathways/tal.toml` cannot be written honestly from this paper.** Three
independent blockers, any one of which is sufficient:

| # | blocker | severity |
| --- | --- | --- |
| 1 | The product axis — per-strain, per-temperature TAL — exists **only as bar heights in Fig. 4B and Fig. S7**. No table, no source data, no repository. | fatal, and not fixable by re-reading |
| 2 | **TAL is assayed in the supernatant.** It is an extracellular product, so its terminal node's `fate` is `secreted`, which `pathway/spec.py` refuses at load. | fatal, and structural |
| 3 | Expression is numeric at **30 °C only**. At 37 and 41 °C it is a log2-fold bar chart (Fig. 3C), and the authors state the expression/product rank correspondence is **lost** at those two temperatures. | fatal for 12 of the 18 states |

No fit was performed, no calibration is shipped, and no number in this directory was read
off a chart.

## Blocker 1 — the product is a bitmap

Figure 4B is captioned "Specific TAL production at late exponential phase at 30, 37 and
41 °C", and its y-axis reads **`TAL (mg/L)/OD600`**. Eighteen bars, six promoters ×
three temperatures, error bars, and **no printed values**. Figure S7 is the same panel at
late stationary phase. Reading the axis *label* off the image is how the units above were
established; reading the bar *heights* off it is digitisation, and this directory does not
do it.

What the article states numerically about that panel, in full:

> "The six promoters tested resulted in a wide range of TAL specific titers, covering a
> **17.8-fold change** between the highest and the lowest measured across all temperatures
> and promoters."

One ratio. That is one constraint on eighteen unknowns.

Everything else the article quantifies about TAL is a **different experiment** and does not
supply the panel:

- 58–64 mg/L (late exponential) and 94–114 mg/L (stationary) — Fig. S8, two strains
  (P*NC1* vs P*ScADH2*), 37 °C, **2 % SD-His glucose**, not the xylose panel.
- 1.2 g/L and 0.82 g/L — Fig. S9, strain **KM1 ΔURA3**, multi-copy plasmid, 1 % SXCA,
  48 h. A different host strain and a different plasmid copy number.

**Searched and not found.** The article carries **no data-availability statement**, and
the strings `github`, `zenodo`, `figshare` and `deposited` do not occur anywhere in the
full text. The supplementary file `mmc1.docx` was downloaded and unpacked: it contains
**Table S1 (primers)**, **Table S2 (promoter sequences)** and ten figure images. There is
no numeric results table in it. The Elsevier vector PDF was not retrievable through PMC,
Europe PMC or the NCBI OA service; and it would not have helped, because bar geometry is
still a plot and not a reported value.

An **erratum exists** — *Metab Eng Commun* 13:e00186 (PMID 34765440, PMC8569586),
"Erratum regarding previously published articles in volumes 9, 10 and 11" — and it is
recorded here because anyone using this paper must check it. It is a blanket
declaration-of-interest notice covering many articles across three volumes, not a data
correction to this one.

## Blocker 2 — TAL leaves the cell, and the spec loader knows it

From Materials and methods, "Determination of triacetic acid lactone levels", verbatim:

> "For the HPLC-UV assay, samples were **centrifuged at 2500× g for 5 min and 1 mL of the
> supernatant was collected** and stored at 4 °C for further analysis. […] For the
> spectrophotometric assay, **the supernatant was diluted 20-fold** and absorbance was
> measured […] at 277 nm."

Both assays measure the **medium**. Nothing intracellular is quantified anywhere in the
paper. So the terminal node's fate is `Fate.SECRETED`, and `PathwaySpec.__post_init__`
raises on it:

> "this module models intracellular products only. A secreted species has no `mu*[X]`
> dilution term — it leaves through a transporter at a rate the cell sets, not one growth
> sets"

That is the correct refusal and not an obstacle to route around. Declaring the node
`diluted` instead would make the file load, make the solver return a content in mg/gDCW,
and make every one of those numbers meaningless — which is precisely the failure
`Fate.SECRETED` exists to prevent. `tests/test_tal_pathway.py` pins both halves: that the
secreted declaration is refused, and that the *only* thing standing between the refusal and
a wrong number is which fate the author types.

## Blocker 3 — the expression axis stops at 30 °C, and the authors say so

Promoter strength is numeric **once**, in the caption to Fig. 3B, for the six promoters
that drive 2-PS in the TAL panel. At 37 and 41 °C, Fig. 3C reports only
`log2(RFU/RFU@30 °C)` as unlabelled bars. The article's own reading of the consequence,
verbatim:

> "relative promoter strengths for 2-PS expression […] **matched** the relative EGFP
> expression from the same promoter set **at 30 °C** (Fig. S10A). The consistency between
> the rank order of promoter strength judged by TAL production and separately by EGFP
> fluorescence **was lost at 37 and 41 °C** (Fig. S10)."

So the pairing assumption `docs/SECOND_DATASET_HUNT.md` flagged — expression paired by
promoter identity from a separate EGFP strain — is **stated by the authors to fail** on
two of the three temperatures. Twelve of the eighteen states have no usable expression
value even in principle, and this is a fact from the paper rather than a caution invented
here.

Their proposed mechanism is worth carrying, because it is a direct statement that the
one-enzyme flux law is insufficient for this pathway:

> "TAL biosynthesis does not rely solely on the expression of 2-PS, but also on the
> expression of upstream enzymes […] acetyl-CoA and malonyl-CoA pools also are likely
> affected by temperature and TAL titers may not necessarily correlate exclusively with
> levels of 2-PS protein."

Compare `data/pathways/phb.toml`, whose whole lesson is the same one from the other
direction: precursor supply, not entry-enzyme dosage.

## What IS numeric, and is shipped here

### `lang2020_promoter_strength_xylose_30c.tsv`

From the Fig. 3B caption, verbatim: "The absolute value of RFU/OD for each promoter is
25 ± 18 (ADH1), 47 ± 7 (HHF1), 95 ± 14 (NC1), 10 ± 1 (PGK), 15 ± 11 (SSA3) and 86 ± 3
(TEF3) at 30 °C."

| column | unit | note |
| --- | --- | --- |
| `rfu_per_od600` | RFU/OD600 | background-subtracted (Methods: "All relative fluorescence intensity (RFU) shown in this study are background subtracted", blank plasmid pIW578 as background) |
| `rfu_per_od600_sd` | RFU/OD600 | biological triplicates, standard deviation |
| `relative_to_tef3` | — | **derived here**, `rfu_per_od600 / 86`. The only computed column |

**The caption does not say which carbon source those absolutes belong to**, and the panel
compares two. The assignment to **xylose** is not assumed — it is pinned by the article's
own prose:

> "In glucose, P*PGK* was found to be a medium level promoter, reaching **28 % of P*TEF3***,
> but in xylose expression was reduced to **less than 12 % of P*TEF3***. Growth in xylose
> had the opposite effect on P*ADH1*, increasing expression to **28 % of P*TEF3***."

10/86 = 11.6 %, which is "less than 12 %", and 25/86 = 29.1 %, which is the ~28 % quoted
for ADH1 **in xylose**. If the absolutes were the glucose set, PGK would have to read 28 %
of TEF3 and it reads 11.6 %. The check is reproduced as an assertion in
`tests/test_tal_pathway.py`.

**Two of the six values are barely resolved.** SSA3 is 15 ± 11 and ADH1 is 25 ± 18 —
relative standard deviations of 73 % and 72 % on biological triplicates. Any fit using this
axis inherits that, and a log-space law inherits it badly, because the lower 1-SD edge of
SSA3 is 4 RFU/OD.

**These are EGFP reporter strains, not the TAL producers.** Separate transformants
carrying `promoter–EGFP` on the same low-copy backbone. Expression is paired to the TAL
panel by **promoter identity only**; it was never measured on a cell that makes TAL.

### `lang2020_growth_rates_xylose.tsv`

From the Results, verbatim: "The growth rates on xylose at 30, 37, and 41 °C were
significantly higher than at 45 °C (**0.28 h⁻¹ at 30 °C, 0.34 h⁻¹ at 37 °C, and 0.35 h⁻¹
at 41 °C**)", with "0.14 h⁻¹ with xylose at 45 °C" from the sentence before.

| column | unit | note |
| --- | --- | --- |
| `temperature_c` | °C | |
| `mu_per_h` | 1/h | quoted from the Results text, not read off Fig. 3A |
| `in_scope_for_the_2ps_panel` | yes/no | whether a TAL panel was run at that temperature |

**These are per TEMPERATURE, not per strain**, and they are measured on
CBS6556 ΔHIS3 ΔURA3 carrying the **blank** vector pIW578 — not on any of the six 2-PS
strains, and not on the EGFP strains. Three distinct values for six genotypes: µ carries
**zero** genotype information in this dataset. Within a temperature, dividing by µ is
dividing every strain by the same constant, so it cannot separate strains — which is the
opposite of the β-carotene chemostat, where µ is the setpoint that defines the state.
The 45 °C row is marked out of scope: no TAL panel was run there.

### The product's identity, for whenever the numbers arrive

Triacetic acid lactone is 4-hydroxy-6-methyl-2H-pyran-2-one. **PubChem CID 54675757**,
`MolecularFormula` C6H6O3, `MolecularWeight` **126.11**, `IUPACName`
"4-hydroxy-6-methylpyran-2-one" — fetched from the PUG REST property endpoint by name, not
computed and not remembered. Recorded here because it is the one field of a would-be
`tal.toml` that *is* establishable from a source, and `tests/test_tal_pathway.py` uses it to
build the chain it then demonstrates the loader refusing.

## Two corrections to `docs/SECOND_DATASET_HUNT.md`

Record these; the entry as written is more favourable than the paper supports.

1. **"Product is specific, biomass-normalised — the only candidate whose product is in the
   same class as `q`" is not right.** Figure 4B's axis is `TAL (mg/L)/OD600`: an
   **extracellular concentration divided by an optical density**. It has no time dimension,
   so it is not a specific productivity. It is a specific *titre* — time-integrated, and
   in the same class as the titres the hunt document rejects the other three candidates
   for. Recovering `q` from it needs `q ≈ µ · (P/X)`, which holds only if `q` was constant
   from inoculation through balanced exponential growth, and that is an assumption nobody
   has checked here. Since µ differs by temperature (1.25-fold from 30 to 41 °C), the
   assumption is not neutral: it would rescale the temperature axis.

2. **The catch is worse than "paired by promoter identity".** The hunt document records
   that expression was measured in separate EGFP strains. True, and additionally: the
   authors report that the pairing's rank order **fails** at 37 and 41 °C, and expression is
   not numeric there in any case.

Neither correction changes the hunt's headline — no second chemostat expression series
exists — and TAL was ranked first on the honest grounds that it was the only candidate with
any µ at all. It still is. It is the *product* column that turns out to be unreachable.

## What would make this dataset usable

Exactly one thing, and it is a request to the authors rather than an analysis:
**the source values behind Fig. 4B** — eighteen numbers with their standard deviations.
With those, the fit is still not the β-carotene one (TAL is secreted, so `spec.py` refuses
it, and only the six 30 °C states have a paired expression), but the *statistics* would
for the first time be better resolved than the calibration this repository actually uses:

    beta-carotene   3 strains   permutation floor 1/3!  = 0.167   -- cannot reach 0.05
    TAL             6 strains   permutation floor 1/6!  = 0.00139 -- can reach 0.05

Six genotypes is the first panel in this repository that could return a significant
leave-one-strain-out result rather than a best-of-six one. That is why it was worth the
retrieval, and it is why the missing eighteen numbers are worth asking for.
