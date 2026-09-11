# PHB (poly-3-hydroxybutyrate) — provenance

Three files, one source: **Kocharin K, Chen Y, Siewers V, Nielsen J. 2012**, *Engineering of
acetyl-CoA metabolism for the improved production of polyhydroxybutyrate in Saccharomyces
cerevisiae*, **AMB Express 2:52**. PMID `23009357`, PMC `PMC3519744`, doi
`10.1186/2191-0855-2-52`. Open Access, CC-BY 2.0.

## Where the numbers physically came from

The **Europe PMC full-text JATS XML**, fetched this session:

```
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3519744/fullTextXML
```

That is the publisher's own marked-up text and tables, not a PDF and not a summary. Table 3
and Table 2 were extracted from their `<table-wrap>` elements; every text number below is
quoted verbatim from the same document. **No number in these files was digitised from a
figure**, and the one place where that would have been necessary is recorded as a hole
rather than filled.

| file | what is in it |
| --- | --- |
| `kocharin2012_states.tsv` | one row per measured state. Six rows: SCKK005 and SCKK006 in the aerobic batch bioreactor, and all four strains in shake flask. **Measured values only.** |
| `kocharin2012_text_ratios.tsv` | strain-pair ratios the Results state in words. Kept separate because for two of them **no absolute value exists anywhere in the paper**, so they are not reducible to the states table. |
| `kocharin2012_derived.tsv` | every quantity obtained by arithmetic, with its formula and its inputs. Nothing here is a measurement and the file says so per row. |

---

## 1. What the paper measures, verbatim

### Table 3 — "Yields and kinetic parameters obtained from batch cultivations"

Footnote, verbatim: *"The values were calculated from at least triplicate fermentations
(n ≥ 3) and represent as mean ± SD."* and *"* The values are significantly difference at
p-value ≤ 0.05."*

| | Unit | SCKK005 | SCKK006 |
| --- | --- | --- | --- |
| Maximum specific growth rate | h-1 | 0.27 ± 0.02 | 0.28 ± 0.00 |
| Glucose consumption rate | g (g DW h)-1 | 1.80 ± 0.09 | 2.24 ± 0.33 |
| Biomass yield on glucose | g (g glc)-1 | 0.15 ± 0.01 | 0.13 ± 0.02 |
| Ethanol yield on glucose | g (g glc)-1 | 0.35 ± 0.05 | 0.35 ± 0.07 |
| Glycerol yield on glucose* | g (g glc)-1 | 0.05 ± 0.00 | 0.07 ± 0.00 |
| Acetate yield on glucose* | g (g glc)-1 | 0.02 ± 0.00 | 0 |
| PHB yield on glucose* | mg (g glc)-1 | 0.02 ± 0.01 | 0.13 ± 0.02 |
| Biomass yield on ethanol | g (g EtOH)-1 | 0.46 ± 0.27 | 0.45 ± 0.08 |
| PHB yield on ethanol* | mg (g EtOH)-1 | 0.22 ± 0.04 | 6.09 ± 1.44 |

### Numbers stated only in the running text

Each is quoted exactly as the XML renders it.

- **PHB specific productivity, glucose phase.** *"Finally, as a result of employing the
  acetyl-CoA boost plasmid to channel carbon from ethanol to the PHB pathway, the specific
  productivity of PHB in the glucose phase in the strain carrying the acetyl-CoA boost
  plasmid (SCKK006) was 99.3 ± 4 μmole (g DW-1 h-1), 16.5 times higher compared to the
  strain carrying the empty acetyl-CoA plasmid (SCKK005)."* **See §4 — this value does not
  reconcile with Table 3 and must not be used unresolved.**
- **PHB yield on ethanol, and its ratio.** *"In the ethanol phase, the PHB yield on ethanol
  in SCKK006 was 6.09 ± 1.44 mg (g EtOH)-1, which was approximately 25-fold higher than for
  SCKK005 that had a yield of 0.22 ± 0.04 mg (g EtOH)-1."*
- **Maximum PHB titre.** *"The maximum PHB titer detected during the ethanol phase in
  SCKK005 was 1.85 mg⋅L-1 while SCKK006, which contained both the PHB biosynthesis and the
  acetyl-CoA boost plasmid, reached a titer of 43.11 mg⋅L-1 after 36 h of batch
  fermentation. In SCKK005, the PHB level remained at the same concentration until the end
  of fermentation while the PHB titer in SCKK006 tended to decrease after 50 h."*
- **Maximum acetate.** *"the maximum acetate concentration of 0.42 ± 0 g⋅L-1 detected in
  SCKK005 was 3 times higher than the acetate concentration of 0.12 ± 0.01 g⋅L-1 detected
  in SCKK006."*
- **Shake-flask growth rate, panel-level.** *"Thus, the maximum specific growth rate in
  shake flasks was 0.18-0.20 h-1 while the maximum specific growth rate in the bioreactors
  was 0.27-0.28 h-1."* **This is a range over the strains, not a strain-resolved value.**
  It is carried in `mu_max_panel_lo_per_h` / `mu_max_panel_hi_per_h`, identical on all four
  flask rows precisely so that nobody reads it as strain-specific.
- **Shake-flask genotype ratio.** *"The recombinant strain with both the acetyl-CoA boost
  plasmid and the PHB plasmid (SCKK006) produced an 18 times higher final concentration of
  PHB (at 120 h) compared to the reference strain, SCKK005 (Figure 2 b)."* The **only**
  quantitative statement attached to the four-strain flask panel.
- **Acetate in the cit2Δ strain.** *"more than 6.5 g L-1 of acetate was detected in the
  medium after 120 of fermentation of SCKK009"*. A **lower bound**, and the column is named
  `acetate_120h_lower_bound_g_per_l` for that reason.
- **Medium.** *"Glucose was added at a concentration of 20 g L-1."* (Bioreactor cultivation),
  and for the flasks *"45 mL of defined minimal medium in a 100 mL unbaffled flask"* with
  *"20 g L-1glucose as carbon source"* in the Figure 2 legend.

### Strains — Table 2, and the Methods sentence that constructs them

Table 2 gives, verbatim: `SCKK005 | MATa SUC2 MAL2-8cura3-52 his3-Δ1 | pIYC04/pKK01`,
`SCKK006 | ... | pIYC08/pKK01`, `SCKK009 | ... cit2Δ | pIYC08/pKK01`,
`SCKK010 | ... mls1Δ | pIYC08/pKK01`.

The construction sentence is from **Methods**, not the figure legend: *"Strain SCKK005 was
constructed by transforming plasmids pKK01 and pIYC04 into strain CEN.PK113-11C. Plasmids
pKK01 and pIYC08 were co-transformed into strain CEN.PK113-11C for the construction of
SCKK006. Strain SCKK009 and SCKK010 were constructed by co-transformation of plasmids pKK01
and pIYC08 into SIYC32 and SCIYC33, respectively."*

---

## 2. The chemistry the spec encodes, verbatim

*"The first enzyme of the pathway is acetyl-CoA C-acetyltransferase [EC 2.3.1.9], encoded by
phaA, which catalyzes the condensation of two acetyl-CoA molecules to form acetoacetyl-CoA
(Peoples and Sinskey 1989b). The next step is the reduction of acetoacetyl-CoA to
(R)-3-hydroxybutyryl-CoA, which is catalyzed by NADPH-dependent acetoacetyl-CoA reductase
[EC 1.1.1.36] encoded by phaB (Peoples and Sinskey 1989b). Finally, PHA synthase [EC 2.3.1.-]
encoded by phaC, catalyzes the polymerization of (R)-3-hydroxybutyryl-CoA monomers to PHB
(Peoples and Sinskey 1989a)."*

That single passage supplies the whole chain and the precursor stoichiometry of 2.

Cytosolic, not mitochondrial or peroxisomal: *"only the cytoplasmic thiolase participates in
PHB biosynthesis."*

Promoters, from Methods: *"PhaA was cloned into pSP-GM2 into the SpeI/SacI sites between the
PGK1 promoter and the ADH1 terminator. Then, PhaB was cloned into the BamHI/SalI sites
between the TEF1 promoter and the CYC1 terminator of the same vector to yield pSP-GM2-phaAB.
PhaC was cloned into the MCS of pSP-GM2 vector the TEF1 promoter and the CYC1 terminator."*

### Molar masses, and where each came from

| species | value | source, fetched this session |
| --- | ---: | --- |
| acetoacetyl-CoA | 851.6 g/mol | PubChem PUG-REST, CID 92153, `MolecularFormula` C25H40N7O18P3S, `MolecularWeight` 851.6 |
| (R)-3-hydroxybutanoyl-CoA | 853.6 g/mol | PubChem PUG-REST, CID 11966146, C25H42N7O18P3S, `MolecularWeight` 853.6 |
| PHB, **per repeat unit** | 86.09 g/mol | KEGG `C06143` gives `FORMULA (C4H6O2)n`. PubChem returns `MolecularWeight` **86.09** for C4H6O2 twice independently — CID 7302 (gamma-butyrolactone) and CID 4093 (methacrylic acid) |

Cross-check on the repeat unit by the other route: 3-hydroxybutanoic acid (PubChem CID 441,
C4H8O3, 104.10) minus water (PubChem CID 962, H2O, 18.015) for the ester bond gives 86.085.
Agrees to the rounding.

Cross-check on acetoacetyl-CoA against the vendored genome-scale model: `yeast-GEM.xml.gz`
(v9.0.2, `data/gem/`) carries `s_0367`, `name="acetoacetyl-CoA"`, `compartment="c"`,
`fbc:charge="-4"`, `fbc:chemicalFormula="C25H36N7O18P3S"` — the same species as the −4 anion,
four protons lighter than PubChem's neutral acid. Consistent.

The precursor id in the spec, `s_0373`, is read from the same file:
`<species metaid="s_0373" id="s_0373" name="acetyl-CoA" compartment="c" fbc:charge="-4"
fbc:chemicalFormula="C23H34N7O17P3S">`, with `<compartment id="c" name="cytoplasm"/>`.
**yeast-GEM v9.0.2 has no (R)-3-hydroxybutyryl-CoA metabolite at all** — grepping the
uncompressed SBML for `hydroxybutyr` returns nothing — so the FBA audit layer can anchor on
acetyl-CoA and on acetoacetyl-CoA, and on nothing downstream of those.

---

## 3. What is NOT here, and why

**There is no PHB content in mg per gDW anywhere in the text or the tables.** Grepping the
full XML for `mg/gDW` and for `mg (g DW` returns zero hits. The four-strain content is
plotted in **Figure 2b and only there**, so under this repository's rule it is not a
measurement and it is not in these files. `phb_content_mg_per_gdw` is therefore `NA` on
every row — the column exists so the hole is visible in the data rather than only in prose.

**There is no expression measurement of any kind.** No RT-qPCR, no proteomics, no copy
number. This matters more than it first looks, because `pathway/flux.py` predicts flux from
relative expression of the entry enzyme and there is nothing here to feed it.

**The genotype axis is a precursor axis.** All four strains carry the *same* PHB plasmid
pKK01, so relative PhaA expression is 1.0 across the panel by construction. What varies is
pIYC08 (ADH2 / ALD6 / acsL641P / ERG10) against the empty pIYC04, plus the cit2Δ and mls1Δ
backgrounds. Constant entry-enzyme dosage against a stated 16.5-fold change in specific
productivity is precisely the case the flux law cannot describe. **So this dataset cannot
fit, and cannot score, `pathway/flux.py`** — and no PHB entry belongs in
`pathway/calibrations.py`. Read the other way it is a finding: for this pathway, flux is set
by precursor supply.

**Every state is batch.** Shake flask or batch bioreactor, never a chemostat. The solver's
`d[X]/dt = 0` holds at best during balanced exponential growth on glucose, which is why the
derived quantities below are all glucose-phase.

**PHB is removed by something the balance does not model.** *"the PHB titer in SCKK006 tended
to decrease after 50 h"*, against the Introduction's own *"the lack of enzymes for PHB
depolymerization"*.

**PHB is not secreted, and this was checked rather than assumed.** It is assayed out of the
biomass — *"10–20 mg of dried cells were weighed and boiled in 1 mL of concentrated sulfuric
acid for 60 min"* — and the string `secret` occurs zero times in the whole XML. So
`fate = "diluted"` is right and `PathwaySpec` accepts the spec, correctly.

---

## 4. The arithmetic, and one contradiction inside the paper

All of it is in `kocharin2012_derived.tsv`; this section shows the working.

### 4a. Content per biomass, by unit cancellation

`Y_PHB/glc ÷ Y_X/glc` = (mg PHB / g glc) ÷ (g DW / g glc) = **mg PHB / g DW**. No model and
no molar mass is involved — it is PHB made per biomass made during the glucose phase.

    SCKK005:  0.02 / 0.15 = 0.1333 mg (g DW)-1
    SCKK006:  0.13 / 0.13 = 1.0000 mg (g DW)-1

It equals the cellular content only if PHB is neither degraded nor carried over, which the
falling titre above warns against late in the fermentation.

### 4b. Specific productivity, by the paper's own equation

The paper states its equations verbatim: *"The specific product formation rate was
calculated by using the equation: rp=μmax⋅YspYsx and the specific glucose consumption rate
was calculated by using the equation: rs=μmaxYsx."*

    SCKK005:  rp = 0.27 × 0.02 / 0.15 = 0.0360 mg (g DW h)-1
    SCKK006:  rp = 0.28 × 0.13 / 0.13 = 0.2800 mg (g DW h)-1

In moles of repeat unit, dividing by 86.09 mg/mmol:

    SCKK005:  0.0360 / 86.09 = 4.1817e-4 mmol (g DW h)-1 = 0.4182 μmol (g DW h)-1
    SCKK006:  0.2800 / 86.09 = 3.2524e-3 mmol (g DW h)-1 = 3.2524 μmol (g DW h)-1

**The equation checks out against Table 3's own glucose rate**, which is the reason to trust
it: `rs = μmax / Ysx` gives 0.27/0.15 = **1.80** against Table 3's tabulated 1.80 ± 0.09
(exact), and 0.28/0.13 = **2.154** against 2.24 ± 0.33 (inside one SD).

### 4c. The Table 3 yields reproduce both measured titres to 4.9 %

Independent check, using only Table 3 and the stated 20 g/L glucose:

    ethanol formed        = Y_EtOH/glc × 20 g/L        = 0.35 × 20 = 7.0 g/L   (both strains)
    predicted PHB titre   = Y_PHB/glc × 20  +  Y_PHB/EtOH × 7.0

    SCKK005:  0.02×20 + 0.22×7.0  =  0.40 + 1.54  =  1.94 mg/L   vs measured  1.85  → 1.049×
    SCKK006:  0.13×20 + 6.09×7.0  =  2.60 + 42.63 = 45.23 mg/L   vs measured 43.11  → 1.049×

Both 4.9 % high, by the *same* factor. Table 3 and the reported titres are one coherent
dataset.

### 4d. The contradiction: 99.3 μmol (g DW)-1 h-1 is ~30× too large

The text's own figure for the quantity computed in §4b, for SCKK006, is
**99.3 ± 4 μmol (g DW)-1 h-1**. Section 4b gives **3.2524**. The ratio is **30.5×**.

Three further checks, all of which put the text value on the wrong side:

1. **Against the yield ratio.** The text says 16.5× between the strains. Table 3's PHB yield
   on glucose gives 0.13/0.02 = **6.5×**, and the §4b productivities give 0.28/0.036 =
   **7.78×**. Neither is 16.5.
2. **Against the titre.** 99.3 μmol (g DW)-1 h-1 at μ = 0.28 /h implies a steady-state
   content of 0.0993/0.28 × 86.09 = **30.53 mg/gDW**, i.e. 3.05 % of dry weight. The glucose
   phase alone makes 0.13 × 20 = 2.6 g/L of biomass, so that content implies **79.4 mg/L of
   PHB in the glucose phase** — 1.84× the paper's measured *maximum* titre of 43.11 mg/L,
   which it says was reached in the *ethanol* phase.
3. **Reading it as nmol instead does not rescue it** — that would be 0.0993 μmol (g DW)-1
   h-1, 33× *smaller* than §4b rather than larger.

**Conclusion, recorded and not resolved:** the tabulated yields, the titres and the paper's
own equations agree with each other; the single text value 99.3 does not agree with any of
them. `kocharin2012_states.tsv` carries it in a column named
`q_phb_glucose_phase_reported_umol_per_gdw_h` with a pointer to this section, and
`kocharin2012_derived.tsv` carries the content it implies purely to show the size of the gap.
**Do not use 99.3 as an entry flux.** Use §4b, which is derived from numbers that
cross-check three ways, and label it derived.

---

## 5. What this dataset can and cannot test

**Can:** that `load_pathway("phb")` accepts a three-node polymer pathway and that
`solve_pathway` walks it. Verified: the spec loads, `calibratable` is `False` (one measurable
node — correct, and the spec says why), and feeding the §4b fluxes at the Table 3 growth
rates returns 0.1333 and 1.0000 mg/gDCW with `carbon_closes` exactly 0.

**Cannot:** score anything. Two independent reasons, and each alone is sufficient.

1. *The flux law has no input here.* No expression was measured, and entry-enzyme dosage is
   constant across the whole strain panel by construction.
2. *The solver check is an identity, not a test.* The only route from this paper to a content
   number is `Y_PHB/glc ÷ Y_X/glc`, and at balanced growth the solver's terminal relation
   `content = v_in/μ` reduces to exactly that same ratio when `v_in` comes from
   `rp = μ·Ysp/Ysx`. The solver reproduces the number because it is the same arithmetic, not
   because it predicted anything.

An independent test needs a **measured** PHB content per gDW paired with a growth rate — the
quantity locked in Figure 2b here. Kocharin & Nielsen 2013 (PMID 23514405), already named in
`docs/WHAT_IS_LEFT.md`, reports content in mg/gDW across eleven chemostat steady states and
is the dataset that would actually score this spec.

---

## 6. Note on the shared directory and the shared spec

`data/phb/` is written by two contributors covering two papers, and `data/pathways/phb.toml`
serves both. This document covers **only** the `kocharin2012_*.tsv` files; the provenance of
`kocharin2013_chemostat_states.tsv` is not mine and is not asserted here.

What was checked directly, and is worth having on the record: the spec written for the 2012
genotype axis reproduces the 2013 chemostat table exactly. Feeding that table's eleven
`(q_phb_mmol_per_gdcw_h, mu_per_h)` pairs through `solve_pathway` with
`data/pathways/phb.toml` returns its eleven `phb_mg_per_gdw` values to within 4e-6 relative,
with `carbon_closes` at 1e-18. The two tables were built independently and land on the same
86.09 g/mol repeat-unit basis, so the spec is not tuned to either one.
