# The beta-carotene ceiling, tested against the published literature

No strain in this lab carries the carotenoid pathway, so there is no measurement of ours
to predict. What is testable without wet lab is whether the ceiling `fba/fva.py` computes
actually *bounds* the beta-carotene contents that have been published in
*Saccharomyces cerevisiae*. A bound no published measurement can violate is a decoration.

**§§1–9 test the FBA ceiling (32.84 mg/gDCW, `fba/fva.py` at `growth_fraction=0.90`).
§§10–13, added 2026-09-11, test a DIFFERENT and later assertion against the same
literature: the kinetic content ceiling `pathway/capacity.py` computes, 1.2483 mg/gDCW,
which is asserted to hold "for every genotype at every growth rate", together with the
phytoene `passthrough` node in `data/pathways/beta_carotene.toml`.** Two asserted constants,
one literature. Nothing in §§10–13 changes a claim record, a gate or any code.

**Outcome: the bound as `FINDINGS.md` stated it is refuted.** Four published
beta-carotene contents exceed 32.84 mg/gDCW, the largest by 2.41x, and they do so in the
papers' own units with no unit conversion involved. The cause is not the model and not
the arithmetic — it is the **90% growth fraction**, which is a free parameter in
`scripts/parked/run_d1.py` and not a measurement. Restated at the growth rate each paper
actually implies, nothing in the literature violates the model.

## 1. The ceiling, reproduced

The number comes from `fba/fva.py::product_flux_range`, called by
`scripts/parked/run_d1.py` on the model returned by
`fba/carotenoid.py::add_beta_carotene_pathway`.

| | |
| --- | --- |
| Model | **GECKO `ecYeastGEM_batch.xml`**, built on yeast-GEM 8.3.4 — resolved by `paths.ec_yeast_gem()` |
| Not | plain Yeast9 v9.0.2 (`paths.yeast_gem()`) |
| Growth reaction | `r_2111` |
| Product reaction | `DM_betacarotene_c`, a demand — beta-carotene accumulates in membranes, so there is no exchange |
| Glucose bound | `None`; the ec model's own protein-pool constraint sets uptake |
| Growth fraction | **0.90** of unconstrained maximum |
| Units the code produces | **mmol/gDW/h**, a flux — not mg/gDCW |

**Which model matters, and the choice is right.** `fba/physiology.py` validated against
van Hoek, van Dijken & Pronk 1998 (aerobic glucose batch; PMID 9797269) gives:

| model | growth /h | measured | rel. err |
| --- | ---: | ---: | ---: |
| Yeast9 v9.0.2, measured uptakes applied | 0.346 | 0.400 | **−13.5%** |
| ecYeastGEM 8.3.4 batch, no hand-set uptakes | 0.377 | 0.400 | 5.8% |

A ceiling computed on plain Yeast9 would describe an organism growing 13.5% too slow, so it
would be worthless as a bound. The ec model reproduces the fermentative subset (growth,
glucose, ethanol) to 6–16% and is the only one either document may bound anything with. It
still misses oxygen by 27% and CO2 by 75% when left unconstrained; the ceiling below does not depend on either. The earlier 73%/65% figures were artefacts of an unsourced reference phenotype and of cap_uptake being a no-op on the split ec model.

**What the code returns**, on this machine, both GSMMs present:

| growth % | mu (1/h) | min | max (mmol/gDW/h) | rel. width |
| ---: | ---: | ---: | ---: | ---: |
| 100% | 0.3768 | 0 | ~0 | — |
| 99% | 0.3731 | 0 | 2.075e-03 | 1.000 |
| 95% | 0.3580 | 0 | 1.037e-02 | 1.000 |
| **90%** | **0.3391** | **0** | **2.0747e-02** | **1.000** |
| 75% | 0.2826 | 0 | 5.187e-02 | 1.000 |
| 50% | 0.1884 | 0 | 1.037e-01 | 1.000 |
| 25% | 0.0942 | 0 | 1.556e-01 | 1.000 |

**The conversion to mg/gDCW is not in the module.** It is `run_d1.py`'s
`flux_from_content` inverted. That function is `q = C * mu / MW` with `MW = 536.87` g/mol,
so:

```
C  =  q * MW / mu
   =  0.0207469 mmol/gDW/h  x  536.87 g/mol  /  0.339146 /h
   =  32.84 mg/gDCW
```

That is where "about 33 mg/gDCW" comes from. It is arithmetic on the flux, correctly done.

## 2. The frontier is a straight line, so "33 mg/gDCW" is a statement about the growth fraction

The growth–product frontier is exactly linear — max abs. residual 7.1e-05 mmol/gDW/h over
twelve growth fractions:

```
q_max(mu)  =  0.207469  x  (1 - mu / 0.376829)     mmol/gDW/h
```

with 0.207469 mmol/gDW/h the product flux at zero enforced growth, confirmed directly by
LP. Substituting into `C = q * MW / mu` gives a closed form for the content ceiling:

```
C(mu)  =  111.38 / mu  -  295.58        mg/gDCW
```

| growth fraction f | mu (1/h) | content ceiling (mg/gDCW) |
| ---: | ---: | ---: |
| 0.99 | 0.3731 | 2.99 |
| 0.95 | 0.3580 | 15.56 |
| **0.90** | **0.3391** | **32.84** |
| 0.75 | 0.2826 | 98.53 |
| 0.50 | 0.1884 | 295.58 |
| 0.25 | 0.0942 | 886.75 |
| 0.10 | 0.0377 | 2660.2 |

**The mechanism.** The ec model's protein pool is a fixed budget. Carbon and enzyme
capacity not spent on biomass are available for the product, and the substitution is
one-for-one, which is why the frontier is a line rather than a curve. Expressing the
ceiling as *specific content* then divides by mu, and that division is what makes the
number explode as growth slows. As `mu -> 0` the content ceiling diverges: **FBA places no
finite bound on specific product content in a non-growing cell.**

So `32.84 mg/gDCW` is not a property of the model. It is a property of the model *and* the
choice `growth_fraction=0.90`, which is `product_flux_range`'s default and appears nowhere
in any measurement.

## 3. Published beta-carotene in *S. cerevisiae*

Every PMID below was resolved against NCBI E-utilities and the returned title checked
against the claim. Titres are **as reported, in the paper's own units**.

| PMID | first author, year | journal | strain / genotype | carbon source | mode | titre as reported |
| --- | --- | --- | --- | --- | --- | --- |
| 7765036 | Yamano 1994 | Biosci Biotechnol Biochem 58:1112 | R7; *Erwinia uredovora* crtE/crtB/crtI/crtY on 2u episomal vector, 2 copies crtY | 2% galactose, YNB | 3 d, 30 °C, stationary | 103 ug/g dry wt |
| 17496128 | Verwaal 2007 | Appl Environ Microbiol 73:4342 | CEN.PK113-6B; *X. dendrorhous* crtYB + crtI x2 + crtE + ScTHMG1 truncated, all genomically integrated, TDH3p | 2% glucose, YNB | shake flask, 72 h, 225 rpm | 0.5 -> **5.9 mg/g (dw)** |
| 19801484 | Ukibe 2009 | Appl Environ Microbiol 75:7205 | INVSc1; crtI + crtYB on 2u plasmids, +/- BTS1 | 2% glucose, SC | 20 mL in 100 mL flask, 20 °C, 5 d | ~18 ug/g dw; ~390 ug/g dw with BTS1 |
| 21573686 | Lange 2011 | Appl Microbiol Biotechnol 91:1611 | G175 (YEplac-CaroSA); *X. dendrorhous* genes + *S. aureus* mvaK1 under adh1p, episomal | glucose; YPD and YNB | not stated in abstract | 14.3 mg/L at 9 g CDM/L (YPD); **3,897 ug/g CDM** at 1 g CDM/L (YNB) |
| 22080204 | Yan 2012 | Curr Microbiol 64:152 | industrial wine yeast T73-63 | grape juice | not stated | 5.89 mg/g |
| 22086347 | Yan 2012 | Curr Microbiol 64:159 | recombinant; HMG-CoA reductase overexpression + ergosterol inhibitor | not stated | not stated | 2.05 mg/g dw control -> 6.29 mg/g dw with 100 mg/L ketoconazole |
| 23718229 | Li 2013 | FEMS Microbiol Lett 345:94 | codon optimisation of crtI (5 codons) and crtYB (8); tHMG1 vs *S. aureus* mva | UNVERIFIED | UNVERIFIED | ~200% improvement; **absolute titre UNVERIFIED** — hard paywall, no accessible mirror |
| 23860829 | Xie 2014 | Biotechnol Bioeng 111:125 | decentralized assembly, marker-recyclable integrative pMRI, GAL10-GAL1 bidirectional promoters | galactose (GAL system) | shake flask | 11 mg/g DCW total carotenoids (72.57 mg/L); **7.41 mg/g DCW beta-carotene** |
| 24262517 | Reyes 2014 | Metab Eng 21:26 | crtE/crtYB/crtI producer, adaptive laboratory evolution under periodic H2O2 | not verified | not verified | 6 -> up to 18 mg/g DCW **carotenoids** (not disambiguated) |
| 25423750 | Wang 2014 | Sheng Wu Gong Cheng Xue Bao 30:1204 | BY4742-derived BW02 | not stated | not stated | 1.56 mg/g DCW |
| 25475893 | Xie 2015 | Metab Eng 28:8 | sequential glucose-regulated control of the FPP downstream / upstream / competitive branches | glucose | high-cell-density fermentation | 1156 mg/L (**20.79 mg/g DCW**) carotenoid |
| 27423881 | Olson 2016 | J Ind Microbiol Biotechnol 43:1355 | SM14, evolved from the Reyes 2014 lineage; PTDH3-crtYB/crtI/crtE, dCTT1 | YNB, C:N 8.8 -> 50 | bench-top bioreactor | **25.52 +/- 2.15 mg/g DCW**; parental 5.68 -> 22.58 mg/g dcw |
| 30138874 | Yamada 2018 | Bioresour Technol 268:616 | YPH499; cocktail delta-integration | not stated | 96 h | 52.3 mg/L |
| 30723675 | Rabeharindranto 2019 | Metab Eng Commun 8:e00086 | BY4741; bidomain and tridomain CrtYB–CrtB–CrtI fusions, integrative pMRI, GAL1-GAL10 | 20 g/L galactose, YPG | 50 mL flask, 28 °C, 140 rpm, 72 h | "two times more beta-carotene" than the natural configuration; **absolute values in a figure only — UNVERIFIED** |
| 31380362 | López 2019 | Front Bioeng Biotechnol 7:171 | SM14 and βcar1.2 (integrated cassettes incl. tHMG1) | 20 g/L glucose | 50 mL in 250 mL baffled flask, 72 h; and 1-L fed-batch, 77–80 h | flask **21 mg/gDCW** (SM14, 159.6 mg/L), **5.8 mg/gDCW** (βcar1.2); fed-batch 210 mg/L at 32.3 gDCW/L (**6.48 mg/gDCW**) and 750 mg/L at 107.1 gDCW/L (**6.91 mg/gDCW**). States **1 OD600 = 0.4 g/L** |
| 32236768 | Su 2020 | J Ind Microbiol Biotechnol 47:383 | growth-phase-dependent control of the pathway, MVA step, squalene branch | not stated | 5-L bioreactor, one-stage | 0.034 -> **33.1 mg/g CDW** and 1.48 g/L — **lycopene, not beta-carotene** |
| 32272391 | Cheng 2020 | Bioresour Technol 308:123275 | SR8B, engineered for xylose use and beta-carotene | sorghum xylose hydrolysate, 17.4 / 32 / 66 g/L xylose | not stated | 82.50 / 93.56 / 114.50 mg/L; **8.29 / 8.10 / 7.32 mg/g DCW** |
| 33062054 | Bu 2020 | Biotechnol Biofuels 13:168 | YBX-01 -> YBX-20; Snq2p ABC transporter + ATP supply + membrane fluidity | YPD, 2% glucose + 10 g/L acetate | 50 mL in 250 mL flask, 72 h, 10% dodecane overlay | 149.8 mg/L intracellular + 10.1 mg/L secreted at OD600 22.5 (full text); **no gDCW/L, so unconvertible** |
| 33102463 | López 2020 | Front Bioeng Biotechnol 8:578793 | CEN.PK2; single then multi-copy carotenogenic genes, HMGR1 bottleneck relieved; PhCCD1 for beta-ionone | not stated | shake flask, 72 h | 4 -> 16 -> **32 mg/gDCW total carotenoids**; 33 mg/L beta-ionone |
| 33332529 | Fathi 2021 | FEMS Yeast Res 21:foaa068 | *Y. lipolytica* LIP2/LIP7/LIP8 lipases + crtI/crtYB/crtE integrated | **1% (v/v) olive oil** in YPD or YNB | shake flask, 4 d | 477.9 mg/L (YPD + oil); highest content **46.5 mg/g DCW** (YNB + oil, at only 1.4–1.5 g DCW/L) |
| 33605428 | Zhao 2021 | Biotechnol Bioeng 118:2043 | ARE1/ARE2 overexpression; PAH1/DPP1/LPP1 deletion | not stated | not stated | fold changes only (1.5x, 2x, 2.4x); **no absolute titre in the abstract** |
| 33616900 | Sun 2020 | Biotechnol Bioeng 117:3522 | xylose-fermenting host; crtYB, crtI, crtE, no tHMG1 | **xylose** | fed-batch bioreactor | 772.8 mg/L |
| 34983533 | Bu 2022 | Microb Cell Fact 21:3 | lipid-droplet TAG metabolism + IZH1p-driven ERG9 downregulation + extra constitutive pathway | + 2 mM oleic acid | not stated | **11.4 mg/g DCW** and 142 mg/L |
| 35024023 | Yamada 2022 | Eng Life Sci 22:4 | ultrasound-irradiated two-phase extractive fermentation, isopropyl myristate | not stated | not stated | 264 mg/L total carotenoid |
| 38213763 | Fan 2024 | Biodes Res 6:0026 | multidimensional optimisation | not stated | shake flask | 166.79 +/- 10.43 mg/L |
| 38710418 **WITHDRAWN** | ~~Bubphasawan 2024~~ | Bioresour Technol 2024:130799 — **retracted by the publisher; Crossref carries an explicit retraction and the title is prefixed "WITHDRAWN:". Not evidence. Superseded by PMID 40944773, which reports 2.29 and 2.90 mg/gDCW by HPLC — ~31× and ~21× lower.** | GAL80 deletion, removing the galactose-induction requirement | **sucrose**; and molasses + fish meal | not stated in abstract | 727.8 +/- 68.0 mg/L at **71.8 +/- 0.4 mg/g DCW**; 354.9 +/- 8.2 mg/L at **60.5 +/- 4.3 mg/g DCW** |
| 38829459 | Yamada 2024 | World J Microbiol Biotechnol 40:230 | HP100_74, point and structural mutagenesis | not stated | not stated | 37.6 mg/L (parent 20.1 mg/L) |
| 39215465 | Arhar 2024 | J Appl Microbiol 135:lxae224 | BCY2073; integrated BTS1 + ER-localised CrtYB, episomal high-copy **C-terminally truncated cytosolic CrtI**, tHMG1 | "an optimized medium" (abstract); full text reports 20 g/L glucose + lactate + acetate — UNVERIFIED, paywalled, no PMC | shake flask (full text: 96 h, 20 °C — UNVERIFIED) | **79 mg/g DCW**, a 76-fold improvement; authors call it the highest content for *S. cerevisiae* to date |
| 40168627 | Lin 2025 | J Agric Food Chem 73:9187 | platform strain via systematic engineering + ARTP mutagenesis | not stated in abstract | fed-batch, 5-L bioreactor | **2.09 g/L**, stated as the highest *S. cerevisiae* titre to date; **no content reported** |
| 40944773 | Bubphasawan 2025 | Bioresour Bioprocess 12:96 | sucrose + agricultural by-products | sucrose; by-products | fed-batch mentioned | 23.30 +/- 4.22 mg/L at **2.29 +/- 0.16 mg/g DCW**; 17.02 +/- 0.40 mg/L at **2.90 +/- 0.21 mg/g DCW** |

**Two corrections to the starting list this task was given.** Ukibe 2009 is not a
beta-carotene paper in sake yeast; it is an **astaxanthin** paper in the laboratory diploid
INVSc1, and beta-carotene appears only as a pathway intermediate in partial constructs.
"Xie 2015" and "Xie 2014" are two different papers (PMID 25475893 and 23860829) and both
are included. A "Zhou 2017" beta-carotene paper in *S. cerevisiae* could not be located and
is reported as not found rather than guessed at.

## 4. Unit conversion

Our ceiling is mg/gDCW. Published titres arrive as mg/gDCW, mg/g dry weight, ug/g CDM,
mg/L, and g/L. The rule applied: **convert only where the paper itself supplies the
biomass, and show the arithmetic.** No OD-to-gDCW factor was assumed anywhere.

**Only one conversion is load-bearing**, and it uses a biomass figure printed in the same
paper:

| row | arithmetic | result |
| --- | --- | ---: |
| Lange 2011, YPD (PMID 21573686) | 14.3 mg/L / 9 g CDM/L | 1.59 mg/g CDM |

**Two more conversions were performed and then found to be unnecessary**, which is the
best available check on the method. López 2019 reports fed-batch titres *and* biomass *and*
its own specific contents. Converting from the first two reproduces the third:

| row | our arithmetic | our result | paper's own figure |
| --- | --- | ---: | ---: |
| López 2019, SM14 fed-batch | 210 mg/L / 32.3 gDCW/L | 6.50 mg/gDCW | **6.48 mg/gDCW** |
| López 2019, βcar1.2 fed-batch | 750 mg/L / 107.1 gDCW/L | 7.00 mg/gDCW | **6.91 mg/gDCW** |

Agreement to 0.3% and 1.3%. The paper's own numbers are used in the test table; ours are
shown only to demonstrate that the division does what it should when the inputs exist.

**Unit rescalings, not conversions** (no biomass needed): 103 ug/g -> 0.103 mg/g;
3,897 ug/g CDM -> 3.897 mg/g CDM; 18 and 390 ug/g -> 0.018 and 0.390 mg/g;
79 mg g-1 cell dry weight -> 79 mg/gDCW.

**Rows left unconverted**, because the paper gives a volumetric titre and no biomass at the
time of measurement: PMID 30138874 (52.3 mg/L), 33616900 (772.8 mg/L), 35024023
(264 mg/L), 38213763 (166.79 mg/L), 38829459 (37.6 mg/L), **40168627 (2.09 g/L, the
*S. cerevisiae* record titre)**, and Fathi 2021's YPD arm (477.9 mg/L — the 46.5 mg/g DCW
in that paper is the YNB arm, a different medium, so it cannot be transferred). PMID
33062054 reports 149.8 mg/L intracellular at OD600 22.5 but no gDCW/L and no OD factor;
33605428 reports only fold changes. **No conversion was invented for any of these.** An
OD-to-gDCW factor would convert three of them, and is deliberately not applied — see §8.

**Internal consistency check.** Where a paper reports titre *and* content, dividing them
must return a plausible biomass. It does, in every case, which is the strongest available
evidence that these numbers are not being misread:

| paper | titre / content | implied biomass |
| --- | --- | ---: |
| Bubphasawan 2024, sucrose | 727.8 / 71.8 | 10.14 gDCW/L |
| Bubphasawan 2024, molasses | 354.9 / 60.5 | 5.87 gDCW/L |
| Bubphasawan 2025 | 23.30 / 2.29 | 10.17 gDCW/L |
| Bu 2022 | 142.0 / 11.4 | 12.46 gDCW/L |
| Xie 2014 | 72.57 / 7.41 | 9.79 gDCW/L |
| Xie 2015 | 1156 / 20.79 | 55.6 gDCW/L |
| Cheng 2020, HCB / MCB / NCB | 114.50 / 7.32, 93.56 / 8.10, 82.50 / 8.29 | 15.64, 11.55, 9.95 gDCW/L |

Arhar 2024 divides to 8.51 gDCW/L (672 / 79, its titre UNVERIFIED — full text only), and
the two Bubphasawan papers land on the
same two biomass values (10.1 and 5.9 gDCW/L) from independently reported titre–content
pairs three years apart while differing 30-fold in content. Both are internally consistent;
the 2024 strain is simply much better.

One cross-paper check falls out of this. Bu 2022 reports 142 mg/L, 11.4 mg/gDCW and
OD600 30.30, which implies **0.411 g/L per OD600** — against the **0.40 g/L per OD600**
López 2019 states outright. Two unrelated groups, agreeing to 3%. That is reassuring about
the factor and is still not a licence to apply it to papers that report neither: see §8.

## 5. The test

**28 rows convert.** The bound under test is `C <= 32.84 mg/gDCW`. The last two columns
invert the closed form of §2: `mu_cap = 111.38 / (C + 295.58)` is the highest growth rate
at which the model permits that content at all.

| strain / condition | PMID | mg/gDCW | analyte | ceiling / obs | verdict | implied mu cap (1/h) | % of mu_max |
| --- | --- | ---: | --- | ---: | :---: | ---: | ---: |
| **Arhar 2024, retargeted CrtI** | 39215465 | **79.00** | **beta-car, HPLC** | 0.42 | **OVER — REFUTES** | 0.2974 | 78.9% |
| ~~Bubphasawan 2024, sucrose~~ | ~~38710418~~ | ~~71.80~~ | **WITHDRAWN** | — | **STRUCK** | — | — |
| ~~Bubphasawan 2024, molasses~~ | ~~38710418~~ | ~~60.50~~ | **WITHDRAWN** | — | **STRUCK** | — | — |
| Fathi 2021, YNB + olive oil | 33332529 | 46.50 | **beta-car, HPLC** | 0.71 | **OVER — REFUTES** | 0.3256 | 86.4% |
| Su 2020, 5-L bioreactor | 32236768 | 33.10 | lycopene | 0.99 | **OVER** | 0.3389 | 89.9% |
| López 2020, HMGR1 relieved | 33102463 | 32.00 | total | 1.03 | below | 0.3400 | 90.2% |
| Olson 2016, SM14 | 27423881 | 25.52 | **A453-sum** | 1.29 | below | 0.3469 | 92.1% |
| Olson 2016, parental, high C:N | 27423881 | 22.58 | **A453-sum** | 1.45 | below | 0.3501 | 92.9% |
| López 2019, SM14 flask | 31380362 | 21.00 | **A453-sum** | 1.56 | below | 0.3518 | 93.4% |
| Xie 2015, high-cell-density | 25475893 | 20.79 | total | 1.58 | below | 0.3521 | 93.4% |
| Reyes 2014, evolved | 24262517 | 18.00 | total | 1.82 | below | 0.3552 | 94.3% |
| Bu 2022, final strain | 34983533 | 11.40 | beta-car | 2.88 | below | 0.3628 | 96.3% |
| Cheng 2020, NCB | 32272391 | 8.29 | beta-car | 3.96 | below | 0.3665 | 97.3% |
| Cheng 2020, MCB | 32272391 | 8.10 | beta-car | 4.05 | below | 0.3668 | 97.3% |
| Xie 2014, flask | 23860829 | 7.41 | beta-car | 4.43 | below | 0.3676 | 97.6% |
| Cheng 2020, HCB | 32272391 | 7.32 | beta-car | 4.49 | below | 0.3677 | 97.6% |
| López 2019, βcar1.2 fed-batch | 31380362 | 6.91 | **A453-sum** | 4.75 | below | 0.3682 | 97.7% |
| López 2019, SM14 fed-batch | 31380362 | 6.48 | **A453-sum** | 5.07 | below | 0.3688 | 97.9% |
| Yan 2012, + ketoconazole | 22086347 | 6.29 | beta-car | 5.22 | below | 0.3690 | 97.9% |
| Verwaal 2007 | 17496128 | 5.90 | beta-car | 5.57 | below | 0.3695 | 98.0% |
| Yan 2012, grape juice | 22080204 | 5.89 | beta-car | 5.58 | below | 0.3695 | 98.0% |
| López 2019, βcar1.2 flask | 31380362 | 5.80 | **A453-sum** | 5.66 | below | 0.3696 | 98.1% |
| Lange 2011, YNB + glucose | 21573686 | 3.90 | beta-car | 8.43 | below | 0.3719 | 98.7% |
| Bubphasawan 2025, by-products | 40944773 | 2.90 | beta-car | 11.32 | below | 0.3732 | 99.0% |
| Bubphasawan 2025, sucrose | 40944773 | 2.29 | beta-car | 14.34 | below | 0.3739 | 99.2% |
| Yan 2012, control | 22086347 | 2.05 | beta-car | 16.02 | below | 0.3742 | 99.3% |
| Lange 2011, YPD | 21573686 | 1.59 | beta-car | 20.67 | below | 0.3748 | 99.5% |
| Wang 2014, BW02 | 25423750 | 1.56 | beta-car | 21.05 | below | 0.3749 | 99.5% |

### Correction, 2026-08-28: the four López 2019 rows are NOT β-carotene

They were recorded here as analyte `beta-car`, and `HARD_TESTS.md` §2 called the same number
"total carotenoid". **`HARD_TESTS.md` was right and this table was wrong.** The methods
section of PMID 31380362 settles it:

> "Carotenoid quantification was performed by measuring the absorbance at **453 nm** of the
> hexane extracts, and then converted into concentrations using a standard curve of
> β-carotene ranging from 0.5 to 10 mg/L."

**No chromatography at all.** A single-wavelength absorbance reading of a crude hexane
extract reports the sum of every *coloured* carotenoid present — principally β-carotene
**and lycopene**, plus γ-carotene and neurosporene. Phytoene is colourless (λmax ≈ 286 nm)
and does not interfere. Calibrating against a pure β-carotene standard does **not** make the
measurement specific: it only fixes the extinction coefficient used to convert a summed
absorbance. The paper reports no carotenoid separately from any other.

That distinction is not pedantic here, because **lycopene is the direct precursor of
β-carotene in this pathway and the model predicts it accumulating to more than the product**.
At D = 0.18 /h, `outputs/registered_prediction_D018.csv` gives:

| strain | β-carotene mg/gDCW | lycopene mg/gDCW | what an A453 assay reads | vs the 1.2483 ceiling |
| --- | ---: | ---: | ---: | :---: |
| b-car2 | 0.5354 | 0.2408 | 0.7762 | below (0.62×) |
| b-car3 | 0.9095 | 0.8609 | 1.7703 | **ABOVE (1.42×)** |
| b-car4 | 1.0223 | 1.4503 | 2.4726 | **ABOVE (1.98×)** |

So on the model's own numbers, **an absorbance assay reads above the ceiling in two of three
strains while β-carotene itself stays below it in all three.** An A453 number above 1.2483
is therefore not evidence against the ceiling — it is the expected reading. This is the trap
López fell into, and this table repeated it by relabelling their number.

**Consequence for M5:** absorbance at a single wavelength **cannot** test the ceiling, and
`MEASUREMENTS_NEEDED.md` M5 step 6 is right to require chromatographic separation with an
authentic standard, reporting β-carotene, lycopene and phytoene from the same injection.

**A second, unresolved caveat on the same rows.** López 2019 states `1 OD600 = 0.4 g/L`
"determined experimentally" and the flask biomass is derived that way rather than weighed.
The fed-batch rows are ambiguous: 72.6 OD600 × 0.4 implies 29.0 gDCW/L against the 32.3
gDCW/L in the paper's own Table 2, so those two may be gravimetric. §4 above uses these rows
to cross-check the OD-to-gDCW factor; that check should be read as weaker than it appears
until which biomass is which has been confirmed from the paper.

Yamano 1994 (0.103 mg/g) and Ukibe 2009 (0.018 and 0.390 mg/g) sit two to three orders of
magnitude below and are omitted from the table only for length.

**Verdict: REFUTED.** **Three** rows exceed the ceiling, **two of them beta-carotene
specifically** — Arhar 2024 (79.00 mg/gDCW, HPLC) and Fathi 2021 (46.50, HPLC) — with Su
2020 exceeding it on lycopene. The tightest *surviving* margin is López 2020 at 1.03x, so
the literature has produced a content within 3% of the bound as well as three that break it.
The bound as written was not conservative; it sits inside the published range, and has done
since at least 2021.

*Corrected 2026-08-30.* This paragraph read "five of 28 rows ... four of them beta-carotene"
and named four groups including Pathumthani. Two of those five rows are the Bubphasawan
2024 pair, struck from the table above when PMID 38710418 was found to be **withdrawn by
the publisher** — but the verdict's counts were never restated, so a retracted paper went
on carrying the refutation's headline arithmetic.

The refutation is not a single outlier, and it does not need the withdrawn paper. **Two
independent groups** — Graz (PMID 39215465) and Tehran/Lyngby (33332529) — break it on
beta-carotene by HPLC, and Guangzhou (32236768) matches it on lycopene at 1.008x.

## 6. Which of the three causes it is

**Not the conversion.** The decisive rows — 79 mg/g DCW (PMID 39215465), 71.8 and
60.5 mg/g DCW (38710418 — **since WITHDRAWN, see above; not evidence**) and 46.5 mg/g DCW
(33332529) — are printed by their authors in
mg/g DCW (or `mg g-1 cell dry weight`), the same unit as the ceiling. No biomass estimate,
no OD factor, no arithmetic of ours stands between the
paper and the comparison. All three pass the titre/content consistency check of §4. A
conversion artefact cannot explain this refutation.

**Not the model.** The absolute maximum product flux is 0.207469 mmol/gDW/h at zero
enforced growth, so by §2 the model permits even 79 mg/gDCW at any `mu <= 0.2974 /h`, which
is 78.9% of mu_max. Every refuting row sits comfortably inside the feasible region of the
*same* ec model once growth is allowed to be what it actually was. The stoichiometry and
the protein pool are not contradicted by anything published.

**It is the growth fraction.** `growth_fraction=0.90` is `product_flux_range`'s default. It
was never measured, and it is a bad assumption for exactly the cultures that produce high
contents: carotenoid strains are harvested at 72–120 h, in late-exponential or post-diauxic
phase, at growth rates far below 90% of maximum. The one paper here that reports enough to
see this directly, López 2019, reaches 107.1 gDCW/L in an 80 h fed-batch — a regime in
which mu is small by construction. Choosing f = 0.90 and then comparing against
stationary-phase contents compares two different physiological states.

**The general form of the error.** Because `C(mu) = 111.38/mu - 295.58` diverges as
`mu -> 0`, *any* finite mg/gDCW ceiling is a claim about growth rate wearing the costume of
a claim about metabolism. Quoting one without the other cannot be right.

## 7. What the bound is when stated correctly

Two forms survive.

**(a) The frontier, as a line.** `q_max(mu) = 0.207469 (1 - mu/0.376829)` mmol/gDW/h. This
is the honest statement and it is falsifiable — but only by a paper reporting product flux
and growth rate together, which none of the 30 papers above does. **No published
measurement can currently test it.**

**(b) The growth-free yield bound.** Maximising the demand with growth free on the ec
model gives **0.1984 g beta-carotene per g glucose** (0.207469 mmol/gDW/h against
3.116 mmol glucose/gDW/h), against a carbon-only limit of 0.447 g/g. This needs no growth
rate. Only one paper reports sugar loading and titre together, and it is a xylose paper, so
the comparison is indicative rather than clean:

| PMID 32272391 | sugar fed | titre | yield | vs. 0.1984 g/g |
| --- | ---: | ---: | ---: | ---: |
| HCB | 66.0 g/L xylose | 114.50 mg/L | 0.00173 g/g | 114x below |
| MCB | 32.0 g/L xylose | 93.56 mg/L | 0.00292 g/g | 68x below |
| NCB | 17.4 g/L xylose | 82.50 mg/L | 0.00474 g/g | 42x below |

The yield bound survives by 42–114x. It is a true bound and a very loose one.

## 8. Limitations

- **Abstract-level evidence for most rows.** Full text was read for Yamano 1994, Verwaal
  2007, Ukibe 2009, Rabeharindranto 2019, López 2019, Fathi 2021, Bu 2020, Bu 2022 and
  Fan 2024. Every other row rests on the abstract, which is why "not stated" appears so
  often under carbon source and cultivation mode. A titre in an abstract is the paper's own
  headline number and is the right thing to test, but the conditions attached to it are
  under-reported here.
- **Li 2013 (PMID 23718229) contributes nothing quantitative.** The identifier is verified;
  the titre is behind a paywall with no accessible mirror and is recorded as UNVERIFIED
  rather than estimated.
- **Rabeharindranto 2019 (PMID 30723675) is unscored.** Its absolute contents exist only in
  a bar chart. Reading values off a figure would have added two rows of invented precision.
- **Analyte is not always beta-carotene.** Four rows report *total carotenoids* (24262517,
  25475893, 31380362 fed-batch, 33102463) and one reports *lycopene* (32236768). Lycopene
  is C40H56 with the same 536.87 g/mol, so it is stoichiometrically interchangeable for this
  test; total carotenoids is an upper bound on the beta-carotene fraction, which makes those
  rows conservative in the direction of *not* refuting.
- **Carbon source differs from the model's on most refuting rows.** The ceiling is computed
  on glucose. PMID 38710418 (**withdrawn**) used sucrose and molasses; PMID 33332529 uses olive oil, a much
  more reduced substrate that should raise a carbon-limited ceiling. The bound was stated
  with no carbon-source qualifier, so it is answerable to these, but a glucose-only test is
  narrower. The top row helps here: Arhar 2024's medium is reported (in full text this task
  could not access) as glucose plus lactate and acetate, which would make the largest
  violation a broadly glucose-grown one. **That medium is UNVERIFIED**, so the strictly
  glucose-only version of the test rests on PMID 32236768's 33.1 mg/g CDW — a 1.01x margin,
  too thin to call on its own.
- **No OD-to-gDCW factor was applied, on purpose.** Two papers independently imply
  0.40–0.41 g/L per OD600 (§4), which is tempting to reuse. It is not reused: the factor is
  strain-, medium- and instrument-specific, carotenoid-laden cells scatter light differently
  from the wild type, and three of the unconvertible rows would swing the test if it were
  wrong. Their absence is a reported gap, not a silent normalisation.
- **A 39.5 g/L beta-carotene figure circulating in reviews is *Yarrowia lipolytica*, not
  *S. cerevisiae*.** It is not in this table and must not be compared against a ceiling
  computed on a *S. cerevisiae* GSMM. The *S. cerevisiae* record titre is 2.09 g/L
  (PMID 40168627).
- **The ec model is yeast-GEM 8.3.4; `paths.yeast_gem()` is 9.0.2.** The version mismatch is
  documented in `paths.py` and is unresolved. The ceiling inherits 8.3.4's biomass
  composition and protein-pool calibration.
- **mu_max = 0.376829 /h is the model's, not a measurement.** The measured value in
  `fba/physiology.py` is 0.400 /h. Using the measured mu_max would move the 90% ceiling but
  not the conclusion, since the refutation is about the choice of fraction, not its
  denominator.
- **This is not a test of the twin.** It tests one number the twin's documentation asserts.
  D1's structural result — feasible range `[0, ceiling]`, relative width 1.000, so FBA
  bounds the product and never predicts it — is untouched and is confirmed again by every
  row of the §1 table.

## 9. What would make this a sharp test

One measurement, absent from all 30 papers: **the specific growth rate at the moment of
harvest**, or equivalently the biomass time course near the sampling point. With it, each
published content becomes a point in `(mu, C)` and the frontier of §2 becomes directly
falsifiable instead of vacuous. Without it, the mg/gDCW form of the bound cannot be tested
and should not be quoted.

Second best, and cheaper to obtain from existing papers: **sugar consumed alongside titre**,
which makes form (b) testable on glucose. Only one of the 30 papers reports both, and on
xylose.

For our own strain, when it exists, this is nearly free: record OD at the carotenoid
sampling point and the one before it. That is what turns the ceiling from a decoration into
a bound.

A recent kinetic-model study of recombinant *S. cerevisiae* beta-carotene strains reaches
the same structural conclusion this repository reached from FVA — that the flux limitation
is not stoichiometric (Elizondo 2025, ACS Synth Biol, PMID 40891387).

## Reproduce

```bash
python3 scripts/parked/run_d1.py        # the flux ceiling: 2.0747e-02 mmol/gDW/h at mu=0.3391
```

Needs both GSMMs; set `YSTWIN_EC_YEAST_GEM` and `YSTWIN_YEAST_GEM` if they are not at the
documented defaults. The mg/gDCW figures in §1–2 are `q * 536.87 / mu` on that output, and
the frontier constants come from the same call at other `growth_fraction` values. The
`mu_cap` column of §5 is `111.38 / (C + 295.58)`.

---

# Part II — the kinetic ceiling and the phytoene node (added 2026-09-11)

Part I above tests `fba/fva.py`'s 32.84 mg/gDCW. Part II tests two *different* asserted
constants, both of which are live in the shipped package:

| # | assertion | where it lives | what it says |
| --- | --- | --- | --- |
| **A** | **content ≤ 1.2483 mg/gDCW** | `src/ystwin/pathway/capacity.py` | "content ≤ vmax_per_growth = 2.3252e-3 mmol/gDCW × 536.87 g/mol = 1.2483 mg/gDCW **for every genotype at every growth rate**" |
| **B** | **phytoene is a `passthrough` node** | `data/pathways/beta_carotene.toml`, node `phytoene` | `rate_law = "passthrough"` with `measurable = false`. `pathway/solve.py`'s own rate-law table: `passthrough   v_out = v_in   ->   X = 0`. The phytoene pool is asserted to be **exactly zero** |

Assertion A reproduced on this machine: `2.3252e-3 × 536.87 = 1.2483` mg/gDCW. **A is not a
growth-fraction artefact like Part I's.** Its derivation cancels μ exactly, so "those were
stationary-phase harvests" is not available as a defence — the claim is that no genotype
reaches 1.2483 mg/gDCW at *any* growth rate.

## 10. The ceiling table

### The inclusion rule, applied before anything was counted

A published value enters this table only if **both** hold:

1. **Species-resolved.** β-carotene separated from lycopene and phytoene by chromatography,
   with identity by retention time against an authentic standard. A single-wavelength
   absorbance total is not β-carotene — it sums every coloured carotenoid and reads above any
   one of them (§5 of Part I, and `tests/test_absorbance_cannot_test_the_ceiling.py`).
2. **A real denominator.** mg/gDCW against dry cell mass. **No mg/L, per-OD600 or per-cell
   value was converted anywhere in this section**, and no `gdcw_per_od` was assumed for any
   strain. Papers rejected on this rule are counted in §12, not hidden.

Denominator classes, because they are not equally strong and pretending otherwise would be
the same error as converting an OD:

| class | meaning |
| :---: | --- |
| **W** | the paper states its **weighing procedure verbatim** (drying temperature, time, balance or filter) |
| **S** | the paper reports mg/gDCW and a DCW determination, but **no procedure is recorded in this repository's notes and the full text was not re-opened in this pass** |
| **C** | biomass from a **measured, strain-specific OD→gDCW factor** (not a borrowed one) |

### THE TABLE — every usable published β-carotene content, sorted

| mg/gDCW | × 1.2483 | | PMID | paper, strain | genotype | cultivation | assay — what makes it specific | den. |
| ---: | ---: | :---: | --- | --- | --- | --- | --- | :---: |
| 0.0180 | 0.014× | below | 19801484 | Ukibe 2009 | INVSc1; crtI + crtYB on 2µ, **no** BTS1 | 20 mL in 100 mL flask, 20 °C, 250 rpm, 5 d | HPLC; β-carotene at RT 25.0 min, resolved from astaxanthin, zeaxanthin, canthaxanthin, echinenone | S |
| 0.0478 | 0.038× | below | 30723675 | Rabeharindranto 2019, **Y/B/I** | BY4741 Δgal80 ddp1::CrtE-tHMG1; CrtY, CrtB, CrtI as three **separated** domains | 50 mL YPG, 28 °C, 140 rpm, 72 h | HPLC-PDA 200–600 nm; β-carotene integrated at 454 nm, phytoene at 282 nm, eight carotenoids quantified off Sigma + Carotenature standards, apocarotenal internal standard | **W** |
| 0.3900 | 0.312× | below | 23718229 | Li 2013 | codon-optimised crtI/crtYB + *S. aureus* mva | 72 h | HPLC LC-20AT + photodiode array | S |
| **1.2517** | **1.003×** | **ABOVE** | 30723675 | Rabeharindranto 2019, **Y(yb)B/I** | native bifunctional CrtYB + separate CrtI | as above | as above | **W** |
| 1.9800 | 1.586× | ABOVE | 35547038 | Jiao 2018, **JX8** | BY4741; *B. trispora* GGPPS/CARB/CARRP polycistronic via P2A, episomal | 50 mL minimal, 30 °C, 200 rpm, 72 h | HPLC, YMC **carotenoid C30** 250×4.6 mm 3 µm — the carotenoid-isomer column; β-carotene and lycopene reported separately | **W** |
| 2.0200 | 1.618× | ABOVE | 35547038 | Jiao 2018, **JX5** | as above, CARB first in the operon | as above | as above | **W** |
| 2.1000 | 1.682× | ABOVE | 35547038 | Jiao 2018, JX5 fed-batch | as above | 3 L bioreactor, 2 L working volume, 750 g/L glucose feed, 32 h | as above | **W** |
| 2.1500 | 1.722× | ABOVE | 40944773 | Bubphasawan 2025 | *S. pararoseus* crtE/crtYB/crtI on PGAL1 / PGAL10 / PGAL7 | 5 L fed-batch on sucrose, 120 h | HPLC Vanquish, Hypersil GOLD C18, DAD 450 nm, Sigma PHR1239; identity by RT **and** UV-vis spectrum | S |
| 2.7221 | 2.181× | ABOVE | 30723675 | Rabeharindranto 2019, **Y(yb)B(ek)I** | tridomain CrtYB–CrtB–CrtI fusion | as above | as above | **W** |
| 3.8970 | 3.122× | ABOVE | 21573686 | Lange 2011 | G175 (YEplac-CaroSA); *X. dendrorhous* genes + *S. aureus* mvaK1, episomal | YNB + 2% glucose, 30 °C, 120 rpm, 96 h | HPLC; Table 2 separates β-carotene from lycopene | S |
| 5.1050 | 4.089× | ABOVE | 29468122 | Wadhwa 2016 | ABC 276 (S288c); *R. toruloides* GGPPS + PSY1 + CRTI on centromeric TEF plasmids | 100 mL SD + amino acids, 30 °C, 250 rpm, 5 d | HPLC C18 + photodiode array; β-carotene, lycopene (95 ± 37) **and phytoene (2727 ± 1421 µg/gDCW)** each quantified off its own standard, phytoene from CaroteNature | **W** |
| 5.9180 | 4.741× | ABOVE | 17496128 | Verwaal 2007, **YB/I/E+tHMG1+I** | CEN.PK113-6B; *X. dendrorhous* crtYB + crtI ×2 + crtE + tHMG1, integrated | YNB + 2% glucose, 30 °C, 225 rpm, 72 h | HPLC, Kontron 440 **diode array**, spectra recorded online; five carotenoids resolved in one run | **W** |
| 7.4100 | 5.936× | ABOVE | 23860829 | Xie 2014 | marker-recyclable integrative pMRI, GAL10-GAL1 bidirectional | YPG + 2% galactose, 30 °C, 230 rpm, 72 h | HPLC UV/VIS 450 nm. The same paper's **spectrophotometric total is 11 mg/gDCW**, 1.5× higher — that number is not used here | S |
| 11.400 | 9.132× | ABOVE | 34983533 | Bu 2022 | FY1679-01B YBX-41; lipid-droplet TAG metabolism + IZH1p-driven ERG9 down | YPD + 2 mM oleic acid, 72 h | HPLC Agilent 1200 C18, 450 nm, Sigma β-carotene standard | S |
| 46.500 | 37.250× | ABOVE | 33332529 | Fathi 2021 | CEN.PK113-7D SC-LIP278β; *Y. lipolytica* lipases + crtE/crtYB/crtI, integrated. **Host is *S. cerevisiae*** | YNB + 1% olive oil, 30 °C, 150 rpm, 72 h | HPLC Thermo Discovery HS F5; β-carotene at RT 7.6 min, 450 nm, Sigma C4582 standard | **W** |
| 79.000 | 63.285× | ABOVE | 39215465 | Arhar 2024 | BCY2073; BTS1 + ER-targeted CrtYB integrated, **truncated cytosolic CrtI** episomal, tHMG1 | shake flask, 96 h, 20 °C | HPLC Waters e2695, **286 and 448 nm**, β-carotene *and* phytoene authentic standards (CaroteNature); lycopene never above detection | **W** |

**16 values, 12 papers, eight named strain backgrounds** — INVSc1, BY4741, S288c/ABC 276,
CEN.PK113-6B, CEN.PK113-7D, FY1679-01B, G175 and BCY2073 (the Bubphasawan, Xie and Li hosts
are not recorded here, so eight is a floor) — and **four distinct heterologous gene sources**:
*X. dendrorhous*, *R. toruloides*, *B. trispora*, *S. pararoseus*. The assertion under test is **genotype-independence**, so a spread of
genotypes is exactly what it has to survive.

### The count

> **13 of 16 values exceed 1.2483 mg/gDCW, by 1.003× to 63.3×. 10 of the 12 papers produce at
> least one value above it.** The three values below are a strain with the cyclase but no
> BTS1 (0.018), a strain whose three enzyme domains are deliberately unfused and which is the
> worst construct in its own paper (0.0478), and one 2013 cassette (0.390).
>
> **Every usable published β-carotene content above 0.4 mg/gDCW exceeds the ceiling. There is
> no measurement in the band between 0.390 and 1.2483.** The ceiling does not sit at the top
> of the published distribution, or in it; it sits under all of it bar three trace values.

**Restricting to the strictest denominator class changes nothing.** Of the 10 class-**W**
values — the ones whose paper prints its weighing procedure verbatim — **9 exceed the
ceiling**, from 6 independent papers (Rabeharindranto, Jiao, Wadhwa, Verwaal, Fathi, Arhar),
spanning 1.003× to 63.3×. The verdict does not rest on the rows whose denominator this
repository has not audited.

### The only things the ceiling does bound: the data it was fitted to

| condition | strain | μ (1/h) | β-carotene mg/gDCW | × 1.2483 | lycopene mg/gDCW | × 1.2483 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 2D025 | b-car2 | 0.2543 | 0.3889 | 0.312× | 0.1779 | 0.14× |
| 3D025 | b-car3 | 0.2543 | 0.6223 | 0.498× | 0.2670 | 0.21× |
| 4D025 | b-car4 | 0.2543 | 0.9576 | 0.767× | 0.8294 | 0.66× |
| 2D01 | b-car2 | 0.1010 | 0.8044 | 0.644× | 0.5017 | 0.40× |
| 4D01 | b-car4 | 0.1010 | 0.9850 | 0.789× | **4.0443** | **3.24×** |
| 3D01 | b-car3 | 0.1010 | **1.2073** | **0.967×** | **2.5094** | **2.01×** |

Contents are `q / μ × 536.87` on `data/carotenoid/elizondo2025_steady_states.tsv`, the
repository's own algebra on its own vendored data. The 3D01 figure reproduces the 96.7%
`docs/HARD_TESTS.md` records independently, which is the check on the arithmetic. Denominator
class **C**: Elizondo's biomass comes from a per-strain measured OD→gDCW factor (Table S7),
which is a measured conversion for that strain — not a borrowed one, and not a weighing.

**Two things follow from the right-hand columns, and they are separate.**

* The **cyclase-capacity arithmetic** is about the terminal node only, so lycopene sitting
  above 1.2483 does not contradict it directly.
* The **"membrane-holding limit" reading** of the fitted ceiling — that a yeast cell cannot
  physically carry more than ~1.25 mg/gDCW of a C40 carotenoid, which is the reading
  `data/carotenoid_batch/SOURCE.md` weighs against the storage-engineered rows — **is dead on
  the calibration set's own numbers.** In 4D01 the *same cells* hold 4.0443 mg/gDCW of
  lycopene, 3.24× the asserted ceiling, at a genuine chemostat steady state where balanced
  growth holds and no batch-harvest caveat applies. The C40 that the ceiling says cannot be
  held is being held, in the dataset the ceiling was fitted to.

### What does not rescue assertion A

* **Not "those are stationary-phase harvests where μ → 0".** Assertion A states that μ
  cancels. A bound advertised as μ-independent cannot be defended by the growth rate. *(The
  real structural repair is different and worth naming: the derivation assumes balanced
  growth, and a 72–120 h batch harvest accumulates product with the dilution term switched
  off. That is an argument for restating the claim as a steady-state-only bound — not for
  keeping it as written.)*
* **Not the genotype input.** `pathway/capacity.py` already records that a 100× rise in
  `entry_expression` moves the content by 0.06%, from 1.1065 to 1.2476 mg/gDCW. Nothing a
  caller can turn reaches 1.2517, let alone 79.
* **Not the assay.** Nine of the 13 exceedances come from runs that resolve β-carotene from
  lycopene *and* phytoene in the same injection; three of them (Rabeharindranto, Wadhwa,
  Arhar) quantify phytoene against an authentic CaroteNature or Sigma standard.
* **Not the denominator.** See the class-**W** subtotal above, and §12.

### One record correction this section forces

`data/carotenoid_batch/published_batch_titres.tsv` carries PMID 30723675 with an empty
`content_mg_per_gdcw`, `dfba_ready = no`, and the note *"absolute values in a figure
only — UNVERIFIED"*. **That is wrong.** The absolute contents are in Supplementary Table S3
(`mmc2.docx`), captioned *"Carotenoid content in the different strains (µg/g DCW) after 72
hours of growth in YPG at 28ºC"*, machine-readable, eight carotenoids × sixteen strains, and
the rows reconcile to their own totals (4064.9 + 1251.7 + 252.1 + 39.4 + 374.1 = 5982.2
exactly, the printed total). **This document does not change that file.** Promoting the row
is a data edit and belongs in a separate reviewed act; it is named here so it is not lost.

## 11. The phytoene table

### The assertion, and what it is worth

`solve.py` gives `passthrough  →  X = 0`. So `beta_carotene.toml` asserts that **the phytoene
pool is exactly zero** in a chain whose own comment block cites Verwaal 2007 as the source of
"the bifunctional crtYB and **the phytoene accumulation**". The paper the chemistry comes from
is the paper that measured the pool the node sets to zero.

### The method that can see phytoene — and it is not "HPLC"

Phytoene has three conjugated double bonds and is **colourless**, λmax ≈ **286 nm**. It is
invisible at 450–453 nm. **The discriminator is the detection wavelength, not the
instrument**: an HPLC method specified at a single 450 nm channel is exactly as blind to
phytoene as a plate reader. Fathi 2021 and Bu 2022 run HPLC against authentic β-carotene
standards and still cannot see phytoene, because 450 nm is the only channel. Only a
sub-300 nm channel or a diode-array detector recording the full spectrum can.

Every fraction below comes from a run that has one. **A phytoene fraction needs no
denominator at all** — it is a ratio of two peaks in one injection, so gravimetric dry
weight, OD conversions and cell counts all cancel. The denominator rule of §10 blocks
*contents*; it does not block *fractions*, which is why Chen 2016 appears here and not there.

### THE TABLE — every usable published phytoene fraction

**Group 1 — the repository's own architecture** (a bifunctional synthase/cyclase plus a
separate desaturase, i.e. a complete route to β-carotene):

| phytoene % of quantified carotenoid | strain | PMID | paper | detection | β-carotene in the same injection |
| ---: | --- | --- | --- | --- | --- |
| 23.2% | Y/B/I (three separated domains) | 30723675 | Rabeharindranto 2019 | PDA 200–600 nm, integrated at 282 nm, Sigma phytoene standard | 47.8 µg/gDCW |
| 29% | YB/I/E+I (second crtI copy) | 17496128 | Verwaal 2007 | Kontron 440 diode array, spectra online | 1,627 µg/g dw (68%) |
| 33% | YB/I (episomal) | 17496128 | Verwaal 2007 | as above | 2 µg/g dw (67%) |
| 36.9% | Y(yb)B(ek)I (tridomain fusion) | 30723675 | Rabeharindranto 2019 | PDA, 282 nm, Sigma standard | 2,722.1 µg/gDCW |
| **48%** | **YB/I/E+tHMG1+I** | 17496128 | Verwaal 2007 | Kontron 440 diode array | **5,918 µg/g dw (52%)** |
| 64% | YB/I/E (episomal) | 17496128 | Verwaal 2007 | as above | 85 µg/g dw (17%); lycopene 79 (16%) |
| 67.9% | Y(yb)B/I (native configuration) | 30723675 | Rabeharindranto 2019 | PDA, 282 nm, Sigma standard | 1,251.7 µg/gDCW |
| 73% | YB/I/BTS1 (episomal) | 17496128 | Verwaal 2007 | Kontron 440 diode array | 13 µg/g dw (16%) |
| 86% | YB/I/E (integrated) | 17496128 | Verwaal 2007 | as above | 141 µg/g dw (9%) |
| 92% | YB/I/E+tHMG1 | 17496128 | Verwaal 2007 | as above | 501 µg/g dw (4.5%) |
| 94% | YB/I/BTS1 (integrated) | 17496128 | Verwaal 2007 | as above | 15 µg/g dw (3%) |
| 100% | Y(yb)B (**no CrtI at all**) | 30723675 | Rabeharindranto 2019 | PDA, 282 nm, Sigma standard | none — 9,665.3 µg/gDCW is all phytoene |

**n = 11 strains carrying both a synthase and a desaturase, two independent laboratories
(Wageningen 2007, Toulouse 2019). Range 23.2%–94%. Median 64%.** *(The `Y(yb)B` row is
listed for completeness and excluded from that statistic: with no desaturase, 100% is
arithmetic, not biology.)*

**Group 2 — other engineered *S. cerevisiae* carotenoid strains, same method quality:**

| phytoene % | strain / note | PMID | paper | detection |
| ---: | --- | --- | --- | --- |
| 3.2%–89.0% | twelve CrtB+CrtI strains with **no cyclase**: B(ek)I 3.2, B(gs)I 23.2, B/I 44.7, B(cpr)I 50.5, BI 56.7, B(yb)I 62.3, IB 79.7, I(cpr)B 87.3, I(yb)B 87.4, I(gs)B 89.0, I(ek)B 93.0, B 100 | 30723675 | Rabeharindranto 2019 | 282 nm, Sigma standard |
| 34.4% | ABC 276, *R. toruloides* genes; 2727 ± 1421 of 7927 µg/gDCW quantified | 29468122 | Wadhwa 2016 | 285 nm, A(1%,1cm) = 750, CaroteNature standard |
| 26% | SyBE_Sc14C35, a **lycopene** strain (no cyclase) | 27329233 | Chen 2016 | 287 nm, Sigma phytoene standard |
| 3.99% | SyBE_Sc14D14, the final fed-batch lycopene strain **after deliberate CrtI fine-tuning** | 27329233 | Chen 2016 | as above |

**Refused, and the refusals matter:**

* **Arhar 2024 (39215465) — the best method in the entire set, and no extractable fraction.**
  286 nm, CaroteNature phytoene standard, gravimetric filtration CDW, and it is the strain
  that breaks the ceiling by 63×. Its phytoene values are plotted in figures only; the text
  gives relative changes ("dropped by 56%"). **The figures were not digitised and no number is
  reported.** This is the single highest-value missing measurement in Part II.
* **Yamano 1994 (7765036) — 11% phytoene, secondary and unverified.** The figure is Verwaal
  2007's report of it; the primary is paywalled. Recorded as secondary, not entered as a
  measurement.
* **Elizondo 2025 (40891387) — no phytoene measurement exists, and `measurable = false` is
  verified rather than assumed.** The full text contains "phytoene" exactly once, inside a
  cited reference title. The authors' released repository (`SysBioengLab/BcarGRASP`) carries
  `MetabolitesMeasurements.xlsx` with sheets for glucose, ethanol, acetate, glycerol,
  **lycopene and β-carotene — and no phytoene sheet.**

### Scoring assertion B

**Refuted as a description of the pathway.** In 11 complete β-carotene strains from two
laboratories, phytoene is **23%–94% of quantified carotenoid, median 64%** — normally the
*largest single species* in the chain, not a trace. The starkest single comparison is inside
the repository's own citation: in Verwaal's `YB/I/E+tHMG1+I`, the strain whose 5.918 mg/gDCW
β-carotene is row 12 of §10, phytoene is **5,380 µg/g dw** — a pool 0.91× the size of the
product, asserted by `beta_carotene.toml` to be **zero**.

**What this does NOT license.** None of these strains is an Elizondo strain, and phytoene was
not measured there. So:

* **No numeric correction to `BETA_CAROTENE_FLUX.alpha` is proposed, and none should be made
  from this table.** `alpha` was fitted against `q_lycopene + q_betacarotene`; if the Elizondo
  strains sat anywhere in the 23–94% band their architecture produces elsewhere, the total
  post-synthase flux would be low by roughly 1.3× to 17×. **That is a range imported from
  other genotypes, not a measurement of these ones.**
* What the table *does* support is a documentation repair: `docs/MEASUREMENTS_NEEDED.md`
  currently anchors the expectation on Chen 2016 as *"a phytoene content near 4% of the
  carotenoid total is expected"*. **4% is the best case of a lycopene paper with no cyclase,
  reached only after CrtI was deliberately fine-tuned, and the same paper says "approximately
  26% of the total carotenoid was phytoene" for its earlier strain and that "phytoene was one
  of the major components ... in most of our engineered strains".** The expectation is
  anchored on the wrong pathway and the wrong end of that paper's own distribution.
* **`measurable = false` on the phytoene node is correct and now verified** against the
  authors' released data, not merely asserted. The reachable act is the one
  `docs/MEASUREMENTS_NEEDED.md` already names — ask the Elizondo/Saa group whether raw
  chromatograms were retained **and at what wavelengths they were acquired.** If a single
  450–478 nm channel was used, no re-integration is possible and only new injections can
  close it.

### One more constant that exists, in units this repository cannot use

`beta_carotene.toml` says changing the phytoene node to `saturating` "requires CrtI kinetics
and pool measurements that this dataset does not supply". **The CrtI kinetics now exist, in
*S. cerevisiae*, in vivo**: Fournié & Truan et al., EMBO J 2024 (PMID 39322757) varied
intracellular phytoene 430-fold and fitted an irreversible Michaelis–Menten model, giving for
*X. dendrorhous* CrtI — the repository's own enzyme — Vmax_cell 168 ± 7 µM/h, K½_cell
380 ± 49 µM, kcat_cell 0.030 ± 0.006 s⁻¹.

**They must still be refused, and for a stated reason rather than an absent one.** These are
intracellular **µM**. Converting µM to mmol/gDCW needs a measured cell volume per gDCW that
this repository does not have, and assuming one is precisely the move this codebase refuses
everywhere else. The honest form of the sentence in the TOML is therefore **"exist, in units
this repository cannot convert without an unmeasured factor"**, not "do not exist".

## 12. The honest denominator count

The rules in §10 are strict on purpose. Counting what they threw out is the real finding
about this literature — and it is a **larger** number than the number that survived.

Scope: *S. cerevisiae* papers carrying a β-carotene or carotenoid quantity that reached the
question *"does this yield a usable mg/gDCW β-carotene content?"* Papers excluded earlier on
species or analyte are listed separately so nothing is hidden.

| outcome | n papers | PMIDs |
| --- | ---: | --- |
| **PASSED — species-resolved assay + a real denominator** | **12** | 17496128, 19801484, 21573686, 23718229, 23860829, 29468122, 30723675, 33332529, 34983533, 35547038, 39215465, 40944773 |
| **REJECTED — assay is an absorbance sum, not β-carotene** | **13** | 24262517, 25475893, 27423881, 29789567, 31380362, 32272391, 32478054, 33102463, 33498600, 33616900, 35024023, 36814715, 38650275 |
| **REJECTED — assay never described in the paper** | **2** | 22080204, 22086347 (both PDFs read; neither states a method) |
| **REJECTED — mg/L or per-OD only, no denominator at all** | **5** | 30138874, 33062054, 38213763, 38829459, 40168627 |
| **REJECTED — mg/gDCW printed, weighing step described NOWHERE** | **5** | 27329233, 29161329, 39581972, 40615910, 42662202 |
| **no absolute number published (fold changes only)** | **1** | 33605428 |
| *blocked on access, not on the standard* | *13* | *7765036, 20632327, 23861041, 24486029, 25423750, 26179622, 26749524, 29170968, 31595456, 34865010, 37953664, 38324606, 38412934* |
| *out of scope — wrong host (*Y. lipolytica*)* | *5* | *32636824, 33746920, 34690947, 36094200, 37055369* |
| *out of scope — wrong analyte* | *5* | *26063466 (β-ionone), 32236768 (lycopene), 34268298 (ergosterol), 36369967 (lutein), 36408203 (astaxanthin)* |

> **25 papers were rejected on the measurement standard against 12 that passed.**
> **15 of the 25 fail on the assay** (13 absorbance sums + 2 with no stated method);
> **10 fail on the denominator** (5 with no denominator at all, 5 that print one and never
> say how the gram was obtained).

### The subtler failure mode is the important one

The expected blocker was mg/L-and-per-OD reporting. It is real — five papers, including two
(33062054 Bu 2020, 38213763 Fan 2024) with a clean specific HPLC assay that a single
supplementary sentence would unblock. **But it is not the dominant mode among the papers this
pass newly surfaced.** The dominant mode is worse: **modern, open-access, impeccable
chromatography, a stated mg/gDCW denominator, and no weighing step anywhere in the text.**

| PMID | paper | value refused | would have been | why refused |
| --- | --- | ---: | ---: | --- |
| 40615910 | Li 2025 (Car09, three XdCrt copies) | 23.6 mg/gDCW | **18.9×** | "dry cell weight" 0 occurrences, "dried" 0, "lyophi" 0, "weigh" 0; DCW appears only inside the conversions themselves |
| 42662202 | Khamwachirapithak 2026 | 9.38 mg/gDCW | **7.5×** | "dry cell weight" appears once — in the Glossary of Abbreviations |
| 40615910 | Li 2025 (Car05) | 9.3 mg/gDCW | **7.5×** | as above |
| 39581972 | Liu 2024 (residual β-carotene, vitamin-A strain) | 4.88 mg/gDCW | **3.9×** | "dry cell weight" 0 occurrences; growth is OD600 throughout |
| 29161329 | Bu 2017 (industrial wine yeast T73-63) | 3.72 mg/gDCW | **3.0×** | HPLC procedure by reference [7] and dry-weight procedure by reference [23], both paywalled |
| 27329233 | Chen 2016 (total carotenoid, lycopene strain) | 60.94 mg/gDCW | — | no dry-weight determination in the full text; its **fractions** are used in §11, which need no denominator |

**An unaudited denominator is indistinguishable from an OD conversion, and an OD conversion is
the one thing this repository will not accept.** That is why the ceiling survived as long as
it did — not because the literature is silent, but because its denominators are mostly
unaudited.

**And note which way the strictness cuts.** Every refused value above that carries a number
**also exceeds 1.2483 mg/gDCW**, by 3.0× to 18.9×. The rule is costing the refutation
evidence, not manufacturing it. A looser rule would make assertion A look *worse*, never
better.

### The clean exception, stated for consistency

`data/carotenoid_batch/SOURCE.md` records an explicit weighing procedure for only two or three
of the nine `dfba_ready` rows. The strict bar was applied to the **new** candidates; it was not
retroactively applied to Ukibe 2009, Lange 2011, Li 2013, Xie 2014 and Bubphasawan 2025, whose
full texts were not re-opened in this pass. They are marked class **S** in §10 and the
class-**W** subtotal is given separately precisely so a reader can discount them without
recomputing anything. **Both papers added in this pass carry the weighing sentence verbatim**,
so they are at least as well-evidenced as any row already in that table.

## 13. What this settles, and what it does not

### Settled

1. **Assertion A — `content ≤ 1.2483 mg/gDCW` for every genotype at every growth rate — is
   false.** 13 of 16 usable published values exceed it, from 10 of 12 papers, across at least
   8 host backgrounds and 4 heterologous gene sources, by 1.003× to 63.3×. Restricted to
   papers that print their weighing procedure verbatim it is still 9 of 10 values from 6
   papers. There is no published β-carotene content between 0.390 and 1.2483 mg/gDCW: the
   ceiling is not near the top of the distribution, it is under all of it except three trace
   values and the six calibration states it was fitted to.
2. **The "membrane-holding limit" reading of that ceiling is refuted by the calibration set
   itself** — b-car4 at μ = 0.101 /h holds 4.0443 mg/gDCW of lycopene, the same C40, 3.24× the
   ceiling, in steady state.
3. **Assertion B — phytoene pool = 0 — is false as a description of this pathway.** Phytoene is
   23%–94% (median 64%) of quantified carotenoid in 11 complete strains from two laboratories,
   and is the single largest carotenoid species in 6 of the 11.
4. **`measurable = false` on the phytoene node is correct**, and now verified against the
   Elizondo authors' own released data rather than asserted.
5. **The literature's real defect is its denominators**, quantified in §12: 25 papers rejected
   against 12 passed, and the growth area is papers that print mg/gDCW without ever saying how
   the gram was obtained.

### NOT settled — and this must not be read as validation of the engine

**This is a test of two asserted constants. It is not a validation of the mechanistic route,
and it cannot become one.**

* **A one-sided bound test is not a prediction test.** Refuting a ceiling shows the model
  permits too little. It says nothing about whether the model gets any *particular* strain's
  content right, and no row above was predicted by anything — each is a number the model is
  compared against after the fact.
* **Nothing here matches environment, genotype and measurement on one strain.** Every row is a
  different genotype, in a different medium, in a different vessel, harvested at a different
  time, measured on a different instrument. Validating the chain
  `environment → stress state → constrained flux → titre` needs all four fixed on the *same*
  strain and varied one at a time. This survey holds none of them fixed.
* **Not one paper in the set reports the specific growth rate at harvest.** That is the single
  missing measurement `docs/HARD_TESTS.md` already names, and it is why the μ-dependent form
  of any bound stays untestable from the literature. (Assertion A's μ-independence is the
  reason it *could* be tested here at all.)
* **The mechanistic layer is untouched by all of this.** `mech/engine.py`, the stress state,
  the environment channel and the flux law are neither confirmed nor refuted by a ceiling
  comparison. Reading §§10–13 as evidence that the twin works would be a category error:
  what has been shown is that one of its stated limits is wrong and one of its stated
  simplifications is contradicted by the pathway's own literature.
* **No claim record was changed.** `data/current_claims.json` is untouched; so is every file
  under `data/frozen_evidence/`, every gate, every ratchet and every line of code. Promoting
  the Rabeharindranto row in `published_batch_titres.tsv` (§10), restating the
  `MEASUREMENTS_NEEDED.md` phytoene expectation (§11) and narrowing the CrtI-kinetics sentence
  in `beta_carotene.toml` (§11) are each a separate reviewed act.

### The three measurements that would close what is left open

1. **Arhar 2024's phytoene fraction** (39215465) — figure source data or an author request.
   It would give a phytoene fraction for the strain that breaks the ceiling by 63×.
2. **Elizondo's acquisition wavelengths** (40891387) — one question to the authors decides
   whether M8 is closable by re-integration or only by new injections.
3. **Verwaal 2010** (20632327, *Yeast* 27:983–998) — the only carbon-limited **chemostat**
   carotenoid paper besides the calibration set, reporting contents at stated dilution rates.
   Wiley returns 403 and there is no mirror; one institutional retrieval would supply the
   `(μ, content)` pairs nothing else in this literature has. Caveat: its abstract's 211 and
   1752 µg/g figures are *total* carotenoids, so it may land in the absorbance tier even once
   obtained.
