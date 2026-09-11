# Glycogen and trehalose contents at near-zero growth — provenance

Provenance for **`boender2009_storage_carbohydrates.tsv`** in this directory. One source,
three measured states, and a unit conversion that has two defensible answers.

**Boender LGM, de Hulster EAF, van Maris AJA, Daran-Lapujade PAS, Pronk JT. 2009.
Quantitative physiology of *Saccharomyces cerevisiae* at near-zero specific growth rates.**
*Appl Environ Microbiol* 75(17):5607–5614. PMID `19592533`, PMC`2737911`,
doi `10.1128/AEM.00429-09`. © 2009 American Society for Microbiology.

**There is a published correction and it does not touch anything used here.** *Appl Environ
Microbiol* 2009 Dec;75(23):7578, doi `10.1128/aem.02344-09`, PMC`2786401`. Fetched and read
this session. It corrects the y-axis unit prefix of **Figure 4C only** — `q_lactate`,
`q_acetate` and `q_succinate` from `mmol · g⁻¹ · h⁻¹` to `μmol · g⁻¹ · h⁻¹`. No glycogen,
trehalose, growth-rate, viability or maintenance value is amended. Every number below stands
as first published.

## Where the numbers physically came from

From the **PMC full text HTML** at `https://pmc.ncbi.nlm.nih.gov/articles/PMC2737911/`,
downloaded and stripped to text this session. Not from the PDF, not from the abstract, not
from a figure. The article is **not open access** — Europe PMC returns no full text for
PMC2737911 — but PMC itself serves the HTML.

**Every value below is from the Results or Materials and Methods prose. Nothing here is
digitised from a figure.** The retentostat time course (Fig. 2), glucose-consumption
kinetics (Fig. 3) and product rates (Fig. 4) exist only as plots and are therefore
**unusable** by the standing rule; they are not in the TSV.

### The one paragraph the storage-carbohydrate numbers come from

Verbatim, Results, final paragraph before `DISCUSSION`:

> In chemostat cultures of *S. cerevisiae*, accumulation levels of glycogen exhibit a
> negative correlation with the specific growth rate (17, 32, 43). To investigate whether
> increased accumulation of storage carbohydrates also occurred at near-zero specific growth
> rates, glycogen and trehalose contents were measured in biomass samples from the
> steady-state chemostat cultures (*D* of 0.025 h⁻¹) and from 22-day retentostat cultures.
> Indeed, glycogen contents increased from 4.3% ± 0.8% (glucose equivalents/biomass) in the
> chemostat cultures to 9.1% ± 0.6% in the retentostat cultures. The glycogen content of the
> retentostat-grown cells was close to that of cells from nitrogen-starved, glucose-grown
> shake-flask cultures (13% ± 1%). The trehalose contents of the chemostat and retentostat
> samples did not differ significantly (1.0% ± 0.4% under both conditions).

**The `±` is unattributed.** The paper never says what it means for these four values. It is
not safe to call it a standard deviation: Fig. 1's legend says "error bars indicate the
standard deviations" of three independent chemostat cultures, while Figs. 2–5 say "average ±
mean deviation of measurements on two independent cultures". The Results text gives no
legend. The TSV columns are therefore named `*_pct_pm` — "plus or minus, as printed" — and
not `*_sd`. **Do not weight anything by these.**

### The assay, verbatim, Materials and Methods

> For analysis of trehalose and glycogen, exactly 20 ml of culture broth was centrifuged (at
> 4°C for 5 min at 10,000 × *g*), washed with cold demineralized water, and finally
> suspended in an exact volume to obtain samples with 5.00 g·liter⁻¹ biomass and stored at
> −20°C. Trehalose and glycogen measurements were performed as described by Parrou and
> Francois (32). Glucose released by glycogen and trehalose conversion was determined using
> the UV method based on Roche kit number 0716251 (Almere, The Netherlands). Trehalose and
> glycogen amounts were determined in triplicate measurements for each sample.

So the instrument reading is **released glucose**, in both assays. That is what makes the
unit conversion below have two readings rather than one. Parrou & Francois 1997 (*Anal
Biochem* 248:186–188), reference 32, was **not fetched** — it is behind Elsevier — so how
*they* express their results is not established here, and the ambiguity is left open rather
than closed by assumption.

### Culture conditions, verbatim, Materials and Methods

> The prototrophic laboratory strain *S. cerevisiae* CEN.PK113-7D (*MAT* a *MAL* 2-8ᶜ *SUC* 2)
> was used in the present study.

> Triplicate anaerobic chemostat cultivations were performed at a dilution rate of 0.025 h⁻¹
> in 2-liter fermentors (Applikon, Schiedam, The Netherlands) with a stirrer speed of 800
> rpm. The working volume was kept at 1.4 liters by means of an electrical level sensor. To
> maintain anaerobic conditions, both the fermentor and the medium vessel were sparged with
> N₂ (5.0; Linde Gas Benelux, The Netherlands) at a flow rate of 0.70 liters·min⁻¹ and ca. 5
> ml·min⁻¹, respectively. … Temperature was kept constant at 30°C, and the pH was controlled
> at 5.0 by automatic addition of 2 M potassium hydroxide.

> Retentostat cultures were operated at a dilution rate of 0.025 h⁻¹ and sparged with
> nitrogen at the same rate as the chemostat cultures to maintain anaerobic conditions.

Anaerobic growth factors were supplied — "ergosterol (final concentration, 10 mg·liter⁻¹)
and Tween-80 (final concentration, 420 mg·liter⁻¹)". This matters for anyone reading
`docs/EXTERNAL_PRODUCT_VALIDATION.md` §2: an unsupplemented anaerobic glucose culture does
not grow in this model, and these cultures were supplemented, so there is no conflict.

**The feed glucose concentration is NOT STATED in the paper.** The medium is cited as
reference 52 (Verduyn synthetic medium) and the 2% glucose figure that does appear belongs
to the *precultures*, not to the chemostat feed; the separate "20 g·liter⁻¹ glucose" belongs
to purity-check *agar plates*. `C_s,in` is defined symbolically in the equations and never
given a value in the text. It is therefore absent from the TSV rather than borrowed.

### Growth rates, verbatim

Chemostat: `D` = 0.025 h⁻¹, and at chemostat steady state µ = D. Used as `mu_per_h = 0.025`.

Retentostat — **the paper gives a bound, not a value, and gives two different doubling
times for it.** Abstract:

> After 22 days of cultivation, specific growth rates had decreased below 0.001 h⁻¹
> (doubling time of >700 h).

Results:

> When specific growth rates were calculated based on the viable fraction of the biomass in
> the retentostat cultures (Fig. 2A), they decreased to below 0.001 h⁻¹, corresponding to a
> doubling time of over 650 h (Fig. 2B).

Both say `< 0.001 h⁻¹`; they disagree on the doubling time (>700 h vs >650 h). Neither is a
measured µ. The TSV records `mu_per_h = NA` and `mu_upper_bound_per_h = 0.001` for this row.
**Neither doubling time is carried into the TSV**, because they conflict and nothing here
needs them.

### Two numbers read but deliberately NOT used

Recorded so a reader knows they were seen and rejected, not missed.

- **Maintenance coefficient, 0.50 mmol glucose·g⁻¹·h⁻¹.** Abstract: "The viable biomass
  concentration in the retentostats could be accurately predicted by a maintenance
  coefficient of 0.50 mmol of glucose g⁻¹ of biomass h⁻¹ calculated from anaerobic,
  glucose-limited chemostat cultures grown at dilution rates of 0.025 to 0.20 h⁻¹." Results
  adds: "this is equivalent to a maintenance requirement for ATP (*m*ATP) of 1 mmol of ATP
  g⁻¹ of biomass h⁻¹." Not in the TSV: it is a whole-cell energetic parameter, it belongs to
  no single row, and `docs/superseded/maintenance-coefficient.md` already records what
  happened the last time a maintenance number entered this repository through a side door.
- **Viability, 91% ± 8% → 79% ± 6%.** Results: "The fluorescence microscopy kit indicated a
  decrease in the viability from 91% ± 8% in the chemostat cultures to 79% ± 6% after 22
  days of retentostat cultivation (Fig. 2C). CFU counts confirmed a decrease of culture
  viability although the loss of viability indicated by this method was more pronounced (ca.
  60% viability after 22 days)." Not in the TSV: the two assays disagree by 19 percentage
  points and the product layer has no viability term to receive either.

## The unit conversion, with the arithmetic, and why there are two answers

The paper reports **percentages**. The solver wants **mmol·gDCW⁻¹**. Getting between them
requires knowing what the numerator of the percentage is, and the paper's own parenthetical
says "**glucose equivalents**/biomass" while the assay measures **released glucose**. Two
readings survive that, they differ by 11%, and neither is refuted by anything I read. Both
are in the TSV, in separate columns, so no downstream user inherits a silent choice.

Molar masses, computed from IUPAC conventional atomic weights C = 12.011, H = 1.008,
O = 15.999:

| species | formula | arithmetic | g·mol⁻¹ |
| --- | --- | --- | ---: |
| glucosyl (anhydroglucose) residue | C₆H₁₀O₅ | 6(12.011) + 10(1.008) + 5(15.999) = 72.066 + 10.080 + 79.995 | **162.141** |
| D-glucose | C₆H₁₂O₆ | 6(12.011) + 12(1.008) + 6(15.999) = 72.066 + 12.096 + 95.994 | **180.156** |
| trehalose | C₁₂H₂₂O₁₁ | 12(12.011) + 22(1.008) + 11(15.999) = 144.132 + 22.176 + 175.989 | **342.297** |

Cross-checked against PubChem this session (`pug/compound/cid/5793,7427/property/…`), which
returns C₆H₁₂O₆ / 180.16 for CID 5793 (glucose) and C₁₂H₂₂O₁₁ / 342.30 for CID 7427
(trehalose). C₆H₁₀O₅ is not looked up: it is the formula **yeast-GEM v9.0.2 itself assigns
to `s_0773`, glycogen**, read out of `data/gem/yeast-GEM.xml.gz` this session.

A percentage *p* (% w/w) is *10p* mg per gDCW.

**Reading A — the literal one.** The percentage is the mass of *glucose released*, per mass
of biomass. Hydrolysis is 1 glucose per glucosyl residue (yeast-GEM `r_0511` and `r_0463`
each release one C6 unit per residue), and 2 glucose per trehalose (`r_0194`:
`trehalose + H₂O → 2 D-glucose`), so:

    glycogen   mmol glucosyl/gDCW  = 10p / 180.156
    trehalose  mmol trehalose/gDCW = 10p / 180.156 / 2

**Reading B — the percentage is already the polymer's own mass.**

    glycogen   mmol glucosyl/gDCW  = 10p / 162.141
    trehalose  mmol trehalose/gDCW = 10p / 342.297

Worked, to the digits in the TSV:

| state | *p* | A: 10*p*/180.156 (÷2 for trehalose) | B |
| --- | ---: | ---: | ---: |
| glycogen, chemostat *D* = 0.025 | 4.3 % | 43 / 180.156 = **0.238682** | 43 / 162.141 = **0.265201** |
| glycogen, retentostat 22 d | 9.1 % | 91 / 180.156 = **0.505118** | 91 / 162.141 = **0.561240** |
| glycogen, N-starved flask | 13.0 % | 130 / 180.156 = **0.721597** | 130 / 162.141 = **0.801771** |
| trehalose, both states | 1.0 % | 10 / 180.156 / 2 = **0.027754** | 10 / 342.297 = **0.029214** |

All in mmol·gDCW⁻¹. Under reading A the chemostat glycogen is `0.238682 × 162.141 = 38.70`
mg of polymer per gDCW — **not** 43, which is the mass of the glucose the assay saw. That
7.4% gap between "what was weighed out" and "what is in the cell" is the whole reason both
columns exist.

The `±` values are **not** propagated into the molar columns. Since it is not established
whether they are SDs, a propagated interval would look like a statistic and would not be
one.

## Row-by-row notes on `boender2009_storage_carbohydrates.tsv`

`NA` means the paper does not state it. It never means zero.

- **`chemostat_D0025`** — the only row `solve_pathway` can accept. µ = D = 0.025 h⁻¹,
  anaerobic, glucose-limited, 30 °C, pH 5.0, CEN.PK113-7D.
- **`retentostat_22d`** — `in_scope_for_solve_pathway = no`, for two independent reasons.
  (i) µ is a bound, not a value. (ii) Even given a value, the terminal solve is `X = v_in/µ`,
  which diverges as µ → 0 and has no degradation term; at near-zero growth the glycogen pool
  is set by synthesis-minus-degradation, not by washout. `dilution_rate_per_h` is still 0.025
  because the pump was left there — in a retentostat with full cell retention D is no longer
  µ, which is the point of the experiment.
- **`n_starved_shake_flask`** — `13% ± 1%`. Almost everything about this culture is `NA`, and
  that is not an oversight in the transcription: **the paper describes this culture nowhere**.
  It appears once, in the Results sentence quoted above, with no Materials and Methods entry,
  no strain, no medium, no growth rate, no temperature, and — checked against the sentence
  character by character — **no citation attached to it**. It is a comparison the authors
  offer without provenance. Kept in the file because it was measured and hiding it would be
  worse; flagged out of scope because a state with no growth rate is not a state this model
  can be given.

## What this dataset can and cannot test

**Can:** the solver's arithmetic, on one row. Given the chemostat content and µ = 0.025 h⁻¹,
the steady-state identity `q = µ·[X]` gives an implied **net** glycogen synthesis flux of
`0.025 × 0.238682 = 5.967 × 10⁻³` mmol glucosyl·gDCW⁻¹·h⁻¹ under reading A
(`0.025 × 0.265201 = 6.630 × 10⁻³` under reading B). Feeding that back through
`solve_pathway` returns 0.238682 (resp. 0.265201) mmol·gDCW⁻¹ with zero round-trip error and
`carbon_closes = 0`. That is a **round trip**, not a prediction — it exercises the solver and
nothing else.

**Cannot:** anything about `pathway/flux.py`. That law needs a relative expression number for
the entry enzyme, and **Boender 2009 measures no transcript, no protein and no enzyme
activity of any gene.** There is no expression column in this file because there is no
expression measurement in the paper. With one in-scope state there is also nothing to hold
out, so no fit and no leave-one-out score is possible even in principle.

Nor is the missing input the only problem: glycogen is a **native** pathway, and
`pathway/flux.py`'s own validation reports every native yeast gene scoring **−0.46 to −0.86**
leave-one-strain-out. The law is calibrated on heterologous cassette dosage. Applying it to
GSY2 would be running it in the regime its authors demonstrated it fails in.

**No calibration file is shipped in this directory.** There is nothing to fit and nothing to
fit it on.

---

# Second file: `boender2011_transcriptome_index.tsv`

Everything above is Boender 2009. This section is a **separate source and a separate file**,
added because the "Cannot" paragraph immediately above is true of Boender 2009 and **not**
true of the companion paper.

**Boender LGM, van Maris AJA, de Hulster EAF, Almering MJH, van der Klei IJ, Veenhuis M,
de Winde JH, Pronk JT, Daran-Lapujade P. 2011. Cellular responses of *Saccharomyces
cerevisiae* at near-zero growth rates: transcriptome analysis of anaerobic retentostat
cultures.** *FEMS Yeast Res* 11:603–620. PMID `22093745`, PMC`3498732`,
doi `10.1111/j.1567-1364.2011.00750.x`.

Read this session from the Europe PMC JATS XML,
`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3498732/fullTextXML`. Open access, so
unlike the 2009 paper the machine-readable full text is available.

## Why this file contains no contents

**It cannot.** Boender 2011's glycogen contents exist only in Fig. 6c. Its sole table is
Table 1, *"Enrichement for MIPS categories (primary categories indicated in capitals), KEGG
pathways (indicated in italics) and TFs (indicated in bold)"* (the misspelling is the
source's) — an enrichment table with no physiology in it. Every glycogen mention in the body
text is qualitative. Figure legend, verbatim, typo included:

> Fig 6 (a) Electron micrographs of *Saccharomyces cerevisiae* before starting the
> retentostat (pictures I and II) and after 22 days in retentostat (pictures III and IV).
> Glycogen was stained with phosphotungstic acid in pictures II and IV. GG, glycogen
> granules; LD, lipid droplets. The size bar respresents 1 μm. (b) Mean-normalized expression
> of genes involved in glycogen synthesis and degradation according to growth rate in the
> combined chemostat/retentostat dataset. (c) Intracellular glycogen contents as a function
> of the specific growth rate in the combined chemostat/retentostat dataset.

`glycogen_content_source = figure_only` on every row, and **nothing was digitised**. So this
file is an **index**: it records which growth rate each transcriptome array set belongs to,
and where the arrays are. Its value is the accessions, not contents.

## The gap this file exists to make visible

The 2009 paper has the contents and only a **bound** on µ. The 2011 paper has a **point** µ
and only a figure for the contents. **No single publication supplies a matched
(µ, glycogen content) pair at near-zero growth**, and the two papers' retentostats are not
the same cultures — 2011, verbatim: *"Duplicate retentostat cultures were started from
independent anaerobic, glucose-limited chemostat cultures grown at the same dilution rate."*

They also disagree about trehalose. Boender 2009 reports 1.0 % ± 0.4 % in both states.
Boender 2011 says, verbatim:

> Expression levels of key genes in trehalose accumulation metabolism also changed at
> near-zero growth rates, but intracellular trehalose levels were below detection level
> throughout the experiments.

No detection limit is given anywhere in the 2011 paper, so `trehalose_status` is a flag and
never a number. **That disagreement is the evidence that a 2009 content must not be paired
with a 2011 µ**, whatever the shared protocol.

## Every number in `boender2011_transcriptome_index.tsv`

**`mu_per_h = 0.0006`, `mu_pm_per_h = 0.0001`** — Results, *"Cultivation at near-zero
specific growth rates in retentostats"*:

> Over 22 days of retentostat cultivation, the estimated specific growth rate progressively
> decreased to 0.0006 ± 0.0001 h⁻¹ ( Fig. 1a ) and the budding index decreased to 15%
> ( Fig. 1d ), which is a typical value for non-growing *S. cerevisiae* ( Lewis *et al*.,
> 1993 ).

The `±` is unattributed here exactly as it is in the 2009 paper — the Fig. 1 legend says
"Data points represent average ± mean deviation of measurements on two independent
cultures", but this value is body text, not a legend. Column named `_pm`, as above. The
budding index (15 %) is **not** carried into the TSV: it is not a state variable this model
consumes.

**`day`, `n_arrays`, `geo_accession`** — Materials and methods, *"Microarrays and
transcriptome analysis"*:

> Independent duplicate retentostat cultures were subjected to microarray analysis at four
> time points after switching the effluent line to the filter unit (2, 9, 16 and 22 days).
> Microarray analysis of independent, triplicate anaerobic glucose-limited chemostat cultures
> grown at a specific growth rate of 0.025 h⁻¹ ( *t* = 0) were also performed as part of this
> study, resulting in a dataset of 11 arrays. These array data can be retrieved from Genome
> Expression Omnibus (GEO, http://www.ncbi.nlm.nih.gov/geo/ ) with series number GSE22574 .

**`n_arrays = 2` on the four retentostat rows is arithmetic, not a printed number, and the
arithmetic is here so it can be checked.** "Independent duplicate … at four time points"
gives 4 × 2 = 8; the sentence's own total is 11; "independent, triplicate" chemostat *t* = 0
supplies 3; 8 + 3 = 11 ✓. So 2 per timepoint is forced, not assumed.

**`mu_per_h = 0.03, 0.05, 0.1, 0.2`, their `n_arrays`, and the ArrayExpress accession** —
same section:

> These data were combined with previously published microarray datasets obtained from
> chemostat cultures grown under identical conditions, but at specific growth rates of 0.03,
> 0.05, 0.1 and 0.2 h⁻¹

> For specific growth rates of 0.03, 0.1 and 0.2 h⁻¹, microarray data were derived from three
> independent replicates, whereas cultures at 0.05 h⁻¹ were performed in duplicate, thereby
> resulting in a final dataset of 22 microarrays. These additional chemostat-based microarray
> data are available from ArrayExpress with the accession number E-MTAB-78
> ( http://www.ebi.ac.uk/microarray-as/ae/ ) and from GEO with the accession number GSE11452 .

Second cross-check: 3 + 2 + 3 + 3 = 11, and 11 + 11 = 22 ✓ matches the stated "final dataset
of 22 microarrays". Both totals close.

**`feed_glucose_g_per_l = 50` / `25`** — same section:

> The only notable difference was the glucose concentration in the feed of 25 g L⁻¹ for
> chemostat cultures at 0.03, 0.05, 0.1 and 0.2 h⁻¹ and of 50 g L⁻¹ for chemostat cultures at
> 0.025 h⁻¹ and retentostat cultures.

This is the feed concentration the 2009 paper never states (see above), and it is recorded
**only** on 2011 rows. It must not be back-filled into the 2009 file: the 2009 paper prepared
its medium in 120-litre batches and 2011 in 40-litre batches, and nothing establishes the two
media were identical.

**`mu_per_h = 0.025` on `chemostat_D0025_t0`** — the dilution rate, at steady state µ = D:
*"Anaerobic glucose-limited retentostats with a dilution rate of 0.025 h⁻¹ and a working
volume of 1.4 L were operated as described previously ( Boender et al ., 2009 )."*

## What this second file changes about what the dataset can test

The "Cannot" paragraph above is correct **for Boender 2009**. Boender 2011 changes its scope
in one specific way and in no other.

**It supplies a genome-wide expression measurement matched to a growth-rate series** — 22
arrays across µ = 0.0006 to 0.2 h⁻¹, at GSE22574 and GSE11452 / E-MTAB-78. That is the input
`pathway/flux.py` needs and the 2009 paper lacks. Neither accession was downloaded or checked
in this session; they are recorded because the paper names them.

**And the paper already reports the answer, in words.** Results, verbatim:

> Several genes involved in glycogen synthesis and degradation were strongly up-regulated at
> near-zero growth rates ( Fig. 6b ). While expression of *GSY1* and *GSY2* , which encode
> glycogen synthases, was not affected, *GAC1* , whose gene product activates the glycogen
> synthases via phosphorylation, was up-regulated at near-zero growth rate.

Read against the 2009 contents — 4.3 % → 9.1 %, a rough doubling across the same transition —
**the entry-enzyme transcript does not move while the pool doubles.** A flux law reading GSY2
mRNA predicts no change where the measurement shows a doubling. Control is post-translational,
through the PP1 targeting subunit GAC1 (`YOR178C`, "protein phosphatase regulator GAC1" in the
vendored SGD annotation), which catalyses no step in the chain and therefore cannot be named
in `entry_enzyme` at all.

That is a **falsification of the entry-enzyme premise on a native pathway**, stated by the
source in prose. It is not a scored result: scoring it would require pulling GSE22574,
matching probes to GSY1/GSY2/GAC1, and pairing them with contents that exist only in a figure.
The first two steps are possible; the third is not, which is why every row here says
`in_scope_for_solve_pathway = no`.

## Relationship between the two files in this directory

They are **complementary and non-overlapping by construction**.

| file | source | what it holds | rows the solver can take |
| --- | --- | --- | --- |
| `boender2009_storage_carbohydrates.tsv` | Boender 2009 | the measured glycogen and trehalose contents | 1 |
| `boender2011_transcriptome_index.tsv` | Boender 2011 | growth rates, array counts, accessions; **no contents** | 0 |

No number appears in both. The 2009 file has no accession column because the 2009 paper
deposited no arrays; the 2011 file has no content column because the 2011 contents are
figure-only. Anyone tempted to join them on `µ` should read "The gap this file exists to make
visible" above first.
