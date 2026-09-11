# Stress rules: what the literature actually establishes

**Claim under test.** *"Stress in a yeast cell still follows certain rules which could be
modelled."* Compiled 2026-08-26. Every PMID below was resolved against PubMed and its title
read back before it was written down; the resolved titles are vendored in
`data/citations/pmid_titles.csv`. Citation lists that already exist are not repeated —
`docs/research/PARAMETER_SOURCES.md` §§ on the ESR and on growth-vs-defence carries the
bibliography, and this note adds the numbers and the verdicts.

---

## 0. The four verdicts, first

**Q1 — Are there quantitative rules?** **Yes, but not the rules the claim assumes.** The
established rules are about *magnitude and membership* — how many genes, what fraction, how
they scale with growth rate — and they are strong. What is **not** established is a rule
mapping *dose* to *response*. Dose is encoded temporally, not in a single saturating curve
(PMID 22179789); the osmotic branch's steady-state output is provably **independent** of dose
(PMID 19596242); a single reporter carries only **1.2–1.3 bits** about the signal, enough for
*which* stress but explicitly **not** *how much* (PMID 25985085); and the response to the same
dose depends on the cell's stress history through non-overlapping genetic routes (PMID
22102822). So "rules" survives; "a dose-response law" does not. Put plainly: **a latent stress
state for yeast is well-founded as a classifier of stress *identity* and poorly founded as an
estimator of stress *magnitude*.** That distinction should be written into the contract before
anything is fitted.

**Q1a — the growth confound, which matters most.** Both extreme positions are refuted.
Growth rate is not a nuisance to be dismissed: **>25%** of all yeast genes are linearly
correlated with growth rate independent of the limiting nutrient (PMID 17959824), **~half**
the genome moves with growth rate with **>80% overlap** with the stress-response genes (PMID
17105650), and slow growth alone confers heat resistance (PMID 19056679) and stress survival
(PMID 21965291). The strongest form of the sceptical case says the ESR is *mostly* an artefact
of cells redistributing into G1 (PMID 24952590) or of metabolic-cycle phase (PMID 22456505).
But the ESR is **not reducible** to it: when division arrest and stress defence were
experimentally decoupled, "the ESR cannot be explained by changes in growth rate or cell-cycle
phase" (PMID 30078561). **A latent stress state is therefore identifiable in
principle and confounded in practice, and the confound is not removable by regressing out
growth — it is a second, partly shared cause.**

**Q2 — Does a latent stress state model already exist?** **Effectively yes; the novelty claim
does not survive as stated.** No paper is titled that, but a ten-reporter yeast panel decoded to
an internal stress representation exists (PMID 29784812), a rank-1 latent state fitted to ~900
yeast fluorescent promoter reporters exists (PMID 24169404), a one-dimensional latent stress
state already coupled to metabolism exists (PMID 22456505), and a five-biosensor yeast toolbox
for reading intracellular state in a bioprocess exists and is genome-integrated (PMID 35069506).
What is still new is the **coupling** — reporter panel → identified state → yeast GEM bounds —
and the identifiability argument for it. §2.

**Q3 — Is regulation-to-metabolism coupling established for yeast?** **The machinery yes, the
benefit no.** "All the methods have a lower overall predictive capability compared to pFBA", and
for yeast specifically every method **except E-Flux** had a median error above pFBA (PMID
24762745). The one quantified head-to-head over nine yeast conditions puts **pFBA at r = 0.8337
and plain E-Flux at 0.7829** — below the no-data control; only E-Flux2 edges ahead, by Δr = 0.035
(PMID 27327084). E-Flux's own scale factor is unidentifiable: refitting it between anaerobic and
aerobic yeast moves it **40.06 → 86.56**, and above ~96 the constraints stop binding at all
(PMID 34839828). Nobody has coupled a stress-module model to yeast bounds; what the yeast stress
literature supports is **raised maintenance ATP** (PMID 27307591, 28779005). And the reason
upper-bound scaling underperforms is measurable: **>50% of capacity sits in reserve in seven
metabolic superpathways, >80% in glucose metabolism** (PMID 32312967), so the bounds being
tightened were mostly not binding. §3.

**Q4 — Beta-carotene parameters.** **The configuration is documented, the kinetics do not
exist, and the bridge is rich.** No Km and no kcat for *X. dendrorhous* crtE, crtYB or crtI in
any host, and **no kcat for phytoene synthase in any organism** — so all six parameters in
`kinetic/carotenoid.py` would be fitted, not sourced. Which step limits is **contested**: the
only published kinetic model of this pathway in yeast puts flux control on **CrtYB**, not crtI
(PMID 40891387), and a third line of evidence says the phytoene block is **compartmentalisation**
rather than turnover (PMID 39215465). **Six** papers report a stress readout and a titre in the
same strain. Two fix the sign of the coupling in *opposite* directions — β-carotene accumulation
makes yeast **more** H₂O₂-sensitive (PMID 19801484) while raised carotenoid rescues `yap1Δ`
(PMID 10673398) — and the best dataset ran in chemostats, which removes the Q1a confound by
design (PMID 20632327, raw arrays public as GEO GSE8451). §4.

---

## 1. Q1: the quantitative rules of the yeast stress response

### 1.1 Size and membership of the ESR — solid

| Claim | Number | Organism / strain | Citation | How solid |
| --- | --- | --- | --- | --- |
| Genes with a common drastic response across ~12 stressors | **~900** | *S. cerevisiae*, DNA microarray | Gasch 2000, PMID 11102521 | **Solid.** The anchor. Replicated in substance by every later study. |
| Genome fraction involved in *any* environmental response | **>50%** | *S. cerevisiae* | Causton 2001, PMID 11179418 | Solid. Note this is not the ESR — it is everything. |
| Common signature across conditions | **~10% of yeast genes** | *S. cerevisiae* | Causton 2001, PMID 11179418 | Solid; the same order as Gasch's ~900 of ~6,000 genes, though not the same number. |
| ESR partition: induced arm | **~300 transcripts (iESR)** | *S. cerevisiae* | Ho 2018, PMID 30078561 | **Solid.** The modern canonical split. |
| ESR partition: repressed arm | **~600 transcripts (rESR)**, ribosomal protein + RiBi | *S. cerevisiae* | Ho 2018, PMID 30078561 | **Solid.** Two-thirds of the ESR is *repression*. |
| ESR conserved beyond *S. cerevisiae* | Present in *S. pombe*, *C. albicans*; presence/absence tracks ecological niche | ascomycetes | Gasch 2007, PMID 17605132 | Qualitative, review-level. |
| The 25-year retrospective on all of the above | — | — | Gasch 2026, PMID 42165638 | The current statement of what is settled and what is not. |

**Consequence for the model.** `src/ystwin/generator/stress_panel.py` models the ESR as one
scalar, Msn2/4 acting through STRE, induction only. That is the **~300-gene third** of the
ESR. The ~600-gene repressed arm runs through Dot6/Tod6 (PMID 19901341) and Sfp1 (PMID
15353587), not Msn2/4, and it is the arm that couples to growth. See §5.

### 1.2 Growth-rate coupling — the strongest set of numbers in this note

| Claim | Number | Organism / condition | Citation | How solid |
| --- | --- | --- | --- | --- |
| Genes linearly correlated with growth rate, independent of limiting nutrient | **>25% of all yeast genes** | 36 steady-state chemostats, 6 limiting nutrients | Brauer 2008, PMID 17959824 | **Solid.** Six nutrients, one linear law. |
| Positively correlated genes | ribosomal function | as above | Brauer 2008, PMID 17959824 | Solid |
| Negatively correlated genes | peroxisomal function | as above | Brauer 2008, PMID 17959824 | Solid |
| G0/G1 population fraction vs growth rate | **linear**, independent of limiting nutrient | as above | Brauer 2008, PMID 17959824 | Solid |
| Genes affected by specific growth rate | **~half of all yeast genes** | chemostats, generation time 2–35 h | Regenberg 2006, PMID 17105650 | Solid |
| Overlap between growth-rate-regulated genes and stress-response genes | **>80%** | as above | Regenberg 2006, PMID 17105650 | **Solid, and this is the confound in one number.** |
| Growth rate predictable from expression of a small gene set | high accuracy, robust across platform, condition, and to *S. bayanus* and *S. pombe* | *S. cerevisiae* + 2 spp. | Airoldi 2009, PMID 19119411 | Solid. An "instantaneous growth rate" is a *measurable latent variable*. |
| Slow growth alone → heat-shock resistance | yes, and independent of respiratory activity | chemostats + heat pulse | Lu 2009, PMID 19056679 | Solid. Causality runs growth → tolerance. |
| Mutant growth rate vs severe-stress survival | **inverse correlation**, confirmed in chemostats | deletion collection; heat, acid, oxidative | Zakrzewska 2011, PMID 21965291 | Solid direction; effect size not a single coefficient. |
| Much of the ESR attributable to metabolic-cycle phase transitions (HOC→LOC) | "much of the common ESR" + cross-protection | *S. cerevisiae*, *S. pombe*, human | Slavov 2012, PMID 22456505 | **Contested.** The strongest form of the growth-explains-ESR argument. |
| **The ESR is not explained by growth rate or cell-cycle phase** | direct test: division arrested, stress applied | transcriptome + proteome + polysome | Ho 2018, PMID 30078561 | **Solid, and it is the refutation of the above.** |
| Functionally unrelated deletion strains share one expression signature; it is **highly similar to the ESR**, and grows stronger the more slowly the strain grows | — | *S. cerevisiae*, deletion strains | O'Duibhir 2014, PMID 24952590 | Solid |
| **The slow-growth signature and the ESR "mainly reflect a redistribution of cells over different cell cycle phases", chiefly more G1 — not a direct single-cell response** | — | as above | O'Duibhir 2014, PMID 24952590 | **Contested.** This is the strongest statement of the confound, and it is what PMID 30078561 was designed to test. |

**Verdict on Q1a.** The confound is real, large (>80% overlap, PMID 17105650), and
bidirectional (PMID 19056679, 21965291: slow growth *causes* stress tolerance). It is also
**not the whole story** (PMID 30078561). A latent stress state that does not carry growth rate
as a co-varying state is measuring growth rate and calling it stress. A latent stress state
that *regresses growth out* deletes two-thirds of the real ESR signal, because the rESR is
where the shared variance lives. **The only defensible structure is a joint state — stress and
growth as two coupled latents, not one plus a correction.** That is what `stress_panel.py`'s
`growth_rate()` already gestures at and does not yet do: it makes growth a deterministic
function of dose, so growth carries no information of its own.

### 1.3 Proteome allocation and the cost of the stress response

Scott/Hwa allocation laws are established for bacteria (PMID 21097934: ribosomal-protein
fraction rises linearly with growth rate, with the linear coefficient set by nutrient quality
and translational capacity). The yeast equivalents exist but are **weaker and messier**:

| Claim | Number | Organism | Citation | How solid |
| --- | --- | --- | --- | --- |
| Bacterial allocation law: empirical constraints on resource allocation to protein synthesis, quantitatively predicting the effect of translation-inhibiting antibiotics on expression and of gratuitous protein expression on growth | phenomenological laws | *E. coli* | Scott 2010, PMID 21097934 | Solid, and **not yeast**. This is the template the project would need a yeast analogue of. |
| Functional protein groups linear in specific growth rate | translation proteins "perfect linear increase"; glycolysis and **chaperones linearly decrease** under respiratory conditions | *S. cerevisiae*, glucose-limited chemostats | Xia 2022, PMID 35595797 | **Solid, and the closest thing to a yeast allocation law.** |
| Excess (non-translating) ribosomal protein | **constant ≈8% of the proteome**; **≈25%** of RPs in fast-growing cells do not translate; the excess fraction **increases as growth rate decreases** | *S. cerevisiae*, many media | Metzl-Raz 2017, PMID 28857745 | **Solid.** A quantified reserve, i.e. cells are *not* at the allocation optimum. |
| Cost of producing unneeded protein is not a single ribosome budget | transcription-limited in low phosphate, translation-limited in low nitrogen, both in rich media; cells adapt by growing larger | *S. cerevisiae* | Kafri 2016, PMID 26725116 | Solid, and it **refutes** the simple "burden = ribosome competition" model. |
| Enzyme usage controls flux in amino-acid biosynthesis | via enzyme-constrained GEM | *S. cerevisiae* | Xia 2022, PMID 35595797 | Model-derived. |
| Proteome mass reallocated from amino-acid biosynthesis into translation on supplementation, raising growth rate | — | *S. cerevisiae* | Björkeroth 2020, PMID 32817546 | Solid. Reallocation, measured. |
| Protein efficiency and allocation determine the fermentative/respiratory switch | — | *S. cerevisiae* | Chen 2019, PMID 31405984 | Solid for Crabtree; not a stress law. |
| Growth-rate-dependent *and* nutrient-specific allocation, i.e. one growth axis does not suffice | — | *S. pombe* | Kleijn 2022, PMID 35228260 | Solid, wrong species; a useful contrast. |
| Under stress, translational capacity is **redistributed**, not uniformly scaled | rESR repression frees translation factors; failure to repress rESR *delays* induced-protein production | *S. cerevisiae* | Ho 2018, PMID 30078561 | **Solid, and it contradicts a uniform capacity multiplier.** |

**Verdict.** There is a yeast allocation law for *growth* (PMID 35595797, 28857745). There is
**no** published Scott/Hwa-equivalent partition law for *stress* — no fitted φ_stress(dose)
trading against φ_ribosome. The mechanism is established (PMID 30078561) and the coefficient
is not. This is the single largest genuinely-open quantitative gap in Q1.

### 1.4 Dose-response form — no canonical functional form exists

This is where the claim is weakest, and the finding is not "nobody has measured it" but
"the measurements say a single curve is the wrong object".

| Claim | Number | Organism | Citation | How solid |
| --- | --- | --- | --- | --- |
| Stress identity *and* intensity are encoded in the **amplitude, duration or frequency** of Msn2 nuclear translocation | three distinct dynamic modes; glucose limitation → pulsatile, dose-dependent **frequency**; oxidative → sustained, dose-dependent **amplitude** | *S. cerevisiae*, single cell | Hao 2011, PMID 22179789 | **Solid.** One agent's dose axis and another's are different variables. |
| Msn2 pulses are frequency-modulated | discrete nuclear localisation bursts | *S. cerevisiae* | Cai 2008, PMID 18818649 | Solid |
| Promoters decode TF dynamics differently; decoding trades off noise against control | promoter-specific thresholds and timescales | *S. cerevisiae* | Hansen 2013, PMID 24189399 | Solid |
| **Information a single Msn2 target gene can carry** | **1.2–1.3 bits** (HXK1, SIP18, amplitude modulation); **1.11 bits** (HXK1, frequency modulation) | *S. cerevisiae*, single cell | Hansen 2015, PMID 25985085 | **Solid.** ~2–2.5 distinguishable levels per reporter. |
| Intrinsic upper bound after removing extrinsic noise | **~1.5–1.6 bits** (AM); **1.36 bits** (FM) | as above | Hansen 2015, PMID 25985085 | Solid |
| Two genes jointly | **1.67 bits** (1× diploid); **1.83 bits** (2× diploid) | as above | Hansen 2015, PMID 25985085 | Solid — **integration buys very little.** |
| Conclusion the authors draw | transduction is "limited to error-free transduction of signal **identity**, but **not** signal intensity information" | as above | Hansen 2015, PMID 25985085 | **Solid, and directly load-bearing.** |
| Four distinct expression programs can be induced by TF dynamics alone | 4 | *S. cerevisiae* | Hansen 2016, PMID 27046808 | Solid (dispatch/commentary format). |
| A panel of TFs encodes environment: **generalists** (Msn2/4, Tod6, Dot6, Maf1, Sfp1) encode *which* of several stresses but **only if stress is high**; **specialists** (Hog1, Yap1, Mig1/2) encode one stress, faster and over a **wider range of magnitudes** | mutual information across 10 TFs | *S. cerevisiae*, single cell | Granados 2018, PMID 29784812 | **Solid, and it is the design rule for a panel.** |
| Dot6 encodes almost as much information as Msn2 | — | as above | Granados 2018, PMID 29784812 | Solid |
| Hog1 nuclear enrichment **perfectly adapts**: steady-state output is *independent* of steady-state input | integral feedback, one integrating mechanism, low noise | *S. cerevisiae*, single cell | Muzzey 2009, PMID 19596242 | **Solid, and fatal to a steady-state Hill for the osmotic arm.** |
| Osmo-adaptation dynamics are dominated by a **fast** Hog1 negative feedback needing no protein synthesis; a much slower feedback through gene expression appears only after large shocks | system identification from periodic stimuli | *S. cerevisiae* | Mettetal 2008, PMID 18218902 | Solid. The transcriptional arm a promoter fusion reads is the *slow* one. |
| A fitted, published dynamic HOG dose-response model exists | biochemical network + thermodynamic volume model | *S. cerevisiae* | Klipp 2005, PMID 16025103 | Solid, and it is a **kinetic** model, not a Hill in dose. |
| Glycerol accumulation, glycolytic flux and growth quantified together under hyperosmotic dose | quantitative time courses | *S. cerevisiae* | Petelenz-Kurdziel 2013, PMID 23762021 | Solid; a candidate parameter source. |
| Expression change under acute stress correlates **poorly** with the genes required to survive that stress | — | *S. cerevisiae* deletion collection | Berry 2008, PMID 18753408 (citing Giaever) | **Solid, and it decouples "induction" from "fitness".** |
| Msn2 and Msn4 play **non-redundant, condition-specific** roles — "arguing against a generic general-stress function" | — | *S. cerevisiae* | Berry 2008, PMID 18753408 | **Solid, and it contradicts a single scalar ESR.** |
| Acquired tolerance to the *same* severe stress proceeds by **non-overlapping** gene sets depending on the prior mild stress | little overlap across salt / DTT / heat pretreatment → H₂O₂ | *S. cerevisiae* deletion libraries | Berry 2011, PMID 22102822 | **Solid.** The response is history-dependent. |
| Cross-protection variation maps to gene-expression variation as a heritable trait | linkage mapping | *S. cerevisiae* segregants | Stuecker 2018, PMID 29649251 | Solid; genotype enters the dose-response. |
| Single-cell ESR is heterogeneous with both intrinsic and extrinsic components | — | *S. cerevisiae*, scRNA-seq | Gasch 2017, PMID 29240790 | Solid; a population mean hides a distribution. |

**Verdict on Q1's dose-response sub-question.** **There is no canonical functional form, and
the literature explains why rather than merely failing to supply one.** Hill-type saturation
in dose is the wrong object for at least three of the panel's stressors: osmotic (perfect
adaptation, PMID 19596242), glucose limitation (frequency-coded, PMID 22179789), and any
stressor applied to a cell with a stress history (PMID 22102822). Where a fitted form exists
it is a *kinetic ODE model of one pathway* (PMID 16025103), not a static dose curve.

---

## 2. Q2: does a latent stress state model already exist?

**Verdict. Yes — every constituent piece is published, and two combinations land close enough
that the novelty claim as stated does not survive.** No paper is *titled* "a latent stress state
model of *S. cerevisiae*", so the phrase is free. But: a ten-reporter yeast panel decoded to a
multidimensional internal stress representation exists (PMID 29784812); a rank-1 latent state
fitted to ~900 yeast *fluorescent promoter reporters* exists, with a resource-allocation model
explaining it (PMID 24169404); the argument that the whole ESR reduces to **one** latent
variable that is itself metabolic exists (PMID 22456505); pushing an inferred regulatory state
into genome-scale flux bounds exists, in *E. coli* (PMID 20876091); and a five-biosensor yeast
toolbox for reading the intracellular state during a bioprocess exists and is genome-integrated
(PMID 35069506). There is also a standing **counter-result**: yeast stress reporters transmit
stress identity but not intensity (PMID 25985085). Claiming "first low-dimensional
representation of the yeast stress response" would be claiming 2000–2013 work.

### 2.1 Prior art, closest first

| What it is | Number that matters | Citation | Distance from this project |
| --- | --- | --- | --- |
| **Ten fluorescent TF reporters (Msn2/4, Dot6, Tod6, Maf1, Sfp1, Mig1/2, Hog1, Yap1) decoded per single cell to extracellular stress identity and magnitude by mutual information** | 10 reporters; generalist/specialist split | Granados 2018, PMID 29784812 | **Closest on the panel side.** Does the panel-design reasoning too. Not coupled to metabolism. |
| **~900 *S. cerevisiae* promoters measured with fluorescent reporters; 60–90% change between conditions by a single global scaling factor depending only on condition, not promoter identity; "only several scaling factors suffice"** | **60–90% rank-1**; ~900 promoters | Keren 2013, PMID 24169404 | **Closest on the latent-state side, and the most dangerous result for a 24-module panel.** See §5. |
| **Five genetically encoded biosensors (pH, ATP, oxidative stress, glycolytic flux, ribosome production) integrated marker-free at one chr X site, in *S. cerevisiae* and *S. boulardii*, for real-time intracellular state in bioprocesses** | 5 sensors, 1 locus | Torello Pianale 2021, PMID 35069506 | **The panel, already built.** Includes the neutrality checks. |
| Same toolbox run under dynamic pH and glucose with quantified robustness metrics | — | Torello Pianale 2025, PMID 40219637 | The dynamic-perturbation use case. |
| **Two-reporter single-cell state (Msn2 + Dot6) with size, growth rate and cell-cycle phase in microfluidics; Msn2/Dot6 discordance defines subpopulations; post-stress growth recovery correlates *positively* with Dot6 activation** | 2 reporters + 4 phenotypes | Bergen 2022, PMID 36350693 | **A multi-reporter latent-state analysis with predictive validity, from the Gasch lab.** |
| The ESR reduces to one latent variable (HOC→LOC metabolic-cycle phase fraction), already coupled to metabolism via O₂ consumption | 1 dimension | Slavov 2012, PMID 22456505 | **The most damaging single hit.** A 1-D latent stress state coupled to metabolism, published 2012. |
| Latent TF-activity state from expression + connectivity, **with explicit identifiability criteria** | — | Liao 2003 (NCA), PMID 14673099 | Owns the identifiability argument for latent regulatory states. |
| Literal linear state-space model, hidden state = TF activity, yeast expression + ChIP | — | Li 2006, PMID 16403793 | This *is* a hidden-state dynamical model of yeast regulation. |
| Variational Bayesian state-space model with explicit hidden factors | — | Beal 2005, PMID 15353451 | Method precedent. |
| ~50 latent module activities fitted **directly to the Gasch stress compendium** | ~50 modules | Segal 2003, PMID 12740579 | A low-dimensional latent description of the yeast stress response, 2003. |
| SVD "eigengenes"/"eigenarrays", arrays classified by cellular state | — | Alter 2000, PMID 10963673 | The canonical early one. Cell cycle, not stress. |
| Fuzzy k-means soft decomposition of the author's own stress compendium | — | Gasch 2002, PMID 12429058 | Overlapping low-dimensional ESR decomposition. |
| ICA latent components on yeast arrays, beating PCA/k-means on functional coherence | — | Lee 2003, PMID 14611662 | Method precedent. |
| iModulons: 92 independently regulated latent signals with a condition × signal **activity matrix** | 92 | Sastry 2019, PMID 31797920 | Prokaryote, but productised. The template exists. |
| Growth rate as a scalar latent physiological variable inverted from expression | >25% of genome linear in it | Brauer 2008, PMID 17959824; Airoldi 2009, PMID 19119411 | A working "physiological index" model. The thing your state must be shown not to be. |
| The slow-growth/ESR signature as one latent confound (G1 fraction) | — | O'Duibhir 2014, PMID 24952590 | The confound the identifiability argument must handle. |
| **Information ceiling: Msn2-type reporters transmit signal identity, not intensity** | 1.2–1.3 bits/gene; 1.67–1.83 bits for two | Hansen 2015, PMID 25985085 | **A published negative result on exactly this identifiability question.** |
| Bayesian inference of inaccessible molecular states from yeast time-lapse reporter data, intrinsic/extrinsic noise separated, Bayes-factor model selection | — | Zechner 2014, PMID 24412977 | Latent-state inference from a yeast reporter, done rigorously. |
| Reporter *dynamics*, not level, carry the state (Crz1, frequency-modulated bursts) | — | Cai 2008, PMID 18818649 | Structural constraint on any static-level state. |
| How promoters decode the latent TF signal | — | Hansen 2013, PMID 24189399; Sweeney 2023, PMID 37087734 | The observation model. |
| **Inferred regulatory state → FBA flux bounds (PROM)** | — | Chandrasekaran 2010, PMID 20876091 | *E. coli* / *M. tuberculosis*. **The pattern is 16 years old.** |
| Kalman-filter state estimation of baker's yeast cultivation from a gas sensor array | — | Yousefi-Darani 2021, PMID 33716616; review PMID 33174065 | Bioprocess soft-sensor precedent, no stress state. |
| Hybrid PLSR + unscented Kalman filter soft sensor reporting intracellular states | — | Hermann 2025, PMID 40564470 | *E. coli*, not yeast. |
| Yeast hybrid fermentation twin, modular for twin platforms | — | Valencia-Velásquez 2025, PMID 40849368 | **No stress state.** |
| Digital twins and model-assisted DoE | — | Kuchemüller 2021, PMID 32797268 | Framing only. |

### 2.2 What would actually be new

Narrow and mechanical, not conceptual. Three candidates, in descending order of defensibility:

1. **Closing the loop reporter panel → identified low-dimensional stress state → time-varying
   constraints on a yeast GEM.** PROM (PMID 20876091) did regulatory-state-to-bounds in
   *E. coli* from transcriptomes, not from reporters. No yeast equivalent driven by an
   engineered reporter panel was found.
2. **A formal identifiability/observability treatment of that state from that panel**, which
   is what `docs/research/IDENTIFIABILITY.md` is for. Caveat: NCA (PMID 14673099) already owns
   identifiability criteria for latent regulatory states, and Keren 2013 (PMID 24169404) plus
   Hansen 2015 (PMID 25985085) mean the honest result may be a *negative* one.
3. **Running it online as a twin.** No published yeast digital twin carries an explicit latent
   stress state.

Any write-up must cite Granados, Keren, Slavov, Hansen and PROM as prior art **up front**, and
must prove the inferred state is not growth rate (PMID 19119411) or metabolic-cycle phase
(PMID 22456505) in disguise.

---

## 3. Q3: coupling regulation to metabolism in yeast

**Verdict. The machinery is established; its value is not, and on yeast the record is
unflattering.** Only three papers validate expression-to-flux methods against yeast ¹³C fluxes.
Machado & Herrgård found "**all the methods have a lower overall predictive capability compared
to pFBA**", and specifically for yeast that every method **except E-Flux** had a median error
above pFBA (PMID 24762745). The one quantified head-to-head puts pFBA at r = **0.8337** and
plain **E-Flux at 0.7829 — below the no-data control** (PMID 27327084). Enzyme-constrained
yeast models are mature, and their budget is a product of three numbers of which two are
assumed. **Nobody has coupled a stress-module activity model to yeast reaction bounds**; what
the yeast stress literature actually supports is a **raised ATP maintenance** mechanism (PMID
27307591, 28779005). The operational consequence for `src/ystwin/bridge/regulation.py`: score
the E-Flux layer against **both** a pFBA control and a raised-NGAM control, or the result means
nothing.

### 3.1 Benchmarks with yeast flux validation — the numbers

| Benchmark | Yeast content | Result | Citation |
| --- | --- | --- | --- |
| 7 methods (GIMME, iMAT, E-Flux, Lee-12, MADE, RELATCH*, GX-FBA) vs a **pFBA control** | *S. cerevisiae* Rintala 5 O₂ levels (+ 2 *E. coli* studies) | "**all of the methods (with the exception of E–Flux) present a median prediction error above that of pFBA**"; overall "**All the methods have a lower overall predictive capability compared to pFBA**"; and "**no improvement in predictive ability when proteomics data is used instead of transcriptomics**". Errors are published as box-plots of a normalised Euclidean error, so *by how much* is not quantified in the paper. | Machado & Herrgård 2014, PMID 24762745 (erratum *PLoS Comput Biol* 2014;10(10):e1003989) |
| Uncentred Pearson vs ¹³C-MFA over 20 conditions, **9 of them yeast** (Yeast5 / iMM904 / iND750) | Rintala O₂ series + Celton acetoin | **pFBA 0.8337 ± 0.1800**; E-Flux2 (with min-L2) **0.8683 ± 0.0964**; **plain E-Flux 0.7829 ± 0.1307, range [0.3506, 0.9223]**; standard FBA 0.7952; SPOT 0.8030; Lee-12 0.5792. **The best expression method beats pFBA by Δr = +0.035; plain E-Flux loses to it.** | Kim & Lun 2016, PMID 27327084 |
| LBFBA (linear-bound FBA) vs pFBA | 5 yeast O₂ conditions, iMM904 | **The only method that cleanly beats pFBA**: lower normalised error in **all 5** yeast conditions, with transcriptomics and proteomics. *E. coli* median error ~**20% vs pFBA ~40%**. Cost: it must be **trained on a paired ¹³C dataset**, and its fitted expression→flux slopes can be **negative**. Also reports that the yeast best-case error was high in 4/5 conditions, i.e. **iMM904 itself disagrees with the measured fluxes**. | Tian & Reed 2018, PMID 29878053 |
| 21 conditions, **no yeast** (*E. coli*, *B. subtilis*, *Synechocystis*, *Synechococcus*) | — | **Uptake rates known: pFBA 0.768 > E-Flux2 0.725 > SPOT 0.650.** Uptake rates unknown: **SPOT 0.583 > E-Flux2 0.385 > pFBA 0.237.** The sharpest statement of when expression data earns its place. | Bhadra-Lobo 2020, PMID 32903284 |
| Fitting E-Flux's implicit proportionality constant (PC) | yeast (Rintala) + *E. coli* | **PC = 40.06 anaerobic vs 86.56 aerobic — a 2.2× swing between two conditions of the same organism.** Above **PC ≈ 96** the aerobic-yeast E-Flux bounds become entirely non-binding. Predicted aerobic ethanol **2.8 ± 1.12 vs 0 measured**; anaerobic **5.86 ± 3.81 vs 9.46**. PC improves flux prediction in ~70% of cases. | Sinha 2021, PMID 34839828 |
| Threshold and method sensitivity of context-specific extraction | human, not yeast | Model size swung from **<600 to >1800 reactions**; **gene-expression threshold is the largest contributor (p < 1.5×10⁻¹³)**. | Opdam 2017, PMID 28215528 |
| Decision sensitivity across 32 human tissues | human | **Thresholding contributes 38.5% of PC1 in active-reaction content; tissue identity only 15.8%.** The preprocessing choice matters more than the biology. | Richelle 2019, PMID 31323017; companion PMID 30986217; guidelines PMID 36566974 |

Framework and review context, for completeness: PMID 25285097 (unifying framework), PMID
22934050 (the taxonomy review), PMID 21943338 (TIGER toolbox), PMID 19714220 (E-Flux's origin
paper — *M. tuberculosis*; **E-Flux has never had a first-party yeast validation**).

### 3.2 Stress coupled to a yeast GEM — several routes, none of them a stress-module model

| Route | What was done | Numbers | Citation |
| --- | --- | --- | --- |
| **Maintenance ATP** — the mechanism the yeast stress literature actually supports | Ethanol, NaCl-osmotic and heat, **three graded levels each, in chemostats fixed at D = 0.1 h⁻¹ so growth rate is controlled**; RNA-seq → condition-specific GEMs → FBA fitted **additional maintenance ATP**. All three stresses converge on raised maintenance ATP → respiration saturates → overflow metabolism. | Respiration linear in glucose uptake only to **~3.4 mmol/gDW/h**; transcriptome PC1 = **63%** of variance, PC2 = **13%** | Lahtvee 2016, PMID 27307591 |
| Same, inside an ecModel | 38 °C glucose-limited chemostat at 0.1 h⁻¹ reproduced **by raising NGAM**, not by expression-mapped bounds | — | Sánchez 2017, PMID 28779005 |
| **Measured proteome under stress on ecYeastGEM** — the closest quantitative precedent | ecYeastGEM + measured absolute proteomes at D = 0.1 h⁻¹ under **high temperature, low pH and KCl osmotic** stress; condition-responsive enzymes defined at relative usage ≥ 0.95 | Conclusion, verbatim: "upregulation and high saturation of enzymes in amino acid metabolism are common across organisms and conditions, suggesting the relevance of **metabolic robustness in contrast to optimal protein utilization** as a cellular objective for microbial growth under stress" | Domenzain 2022, PMID 35773252 |
| **A genuine stress model on an ecYeast backbone** | etcYeast: per-enzyme melting temperature, T_opt and ΔCp, fitted by SMC-ABC | 100 posterior models each **R² > 0.9**; only **9 enzymes (1%)** have mean Tm < 42 °C; **ERG1** predicted most rate-limiting above T_opt and validated by ortholog swap. Enzyme-T_opt prior R² = 0.5, RMSE 13 °C; posterior vs BRENDA T_opt only **r = 0.49 (p = 0.075, n = 14)** | Li 2021, PMID 33420025 |
| **Signalling → ecModel, the closest analogue to this project's design** | Boolean nutrient-signalling module (Snf1, TORC1, PKA) coupled through a regulatory network onto an enzyme-constrained central-carbon yeast model | Mean \|log₁₀\| protein-abundance error **2.62 → 1.55** (respiration) and **3.56 → 2.32** (fermentation); proteins within one order of magnitude **40.8% / 65.5%**. **Read the absolute values: even improved, it is 1.5–2.3 orders of magnitude out.** Gain came from isoenzyme usage, which pure ecModels collapse. | Österberg 2021, PMID 33836000 |
| **Stress transcriptomes → yeast GEM at scale** | **163 single-cell GEMs built with GIMME** from the Gasch 2017 osmotic-stress scRNA-seq (80 high-osmotic, 83 reference) | Random forest on the resulting flux profiles separates stressed from reference cells at **78% accuracy** | Zhang 2024 (Yeast9), PMID 39134886 |
| Regulatory-network integration in yeast | IDREAM = EGRIN + PROM on yeast. States plainly that **PROM "underperforms in eukaryotes"** | IDREAM **2,588 influences (91 TFs → 794 metabolic genes)** vs PROM **31,075 associations (177 TFs)**; TF-deletion growth-defect **AUC 0.67 (IDREAM) vs 0.58 (PROM)** | Wang 2017, PMID 28520713 |
| The original yeast rFBA | **55 TFs regulating 750 metabolic genes**; 10 TF-deletion strains × 12 carbon sources | **13** significant growth over-predictions, reduced to **8** after adding inferred edges. Nutrient regulation, not stress. | Herrgård 2006, PMID 16606697 |
| rFBA and PROM origins | *E. coli*; *E. coli* + *M. tuberculosis* | PROM was **never applied to yeast by its authors** | Covert 2004, PMID 15129285; Chandrasekaran 2010, PMID 20876091 |

### 3.3 The enzyme budget, and how much of it is assumed

| Quantity | Value | Source |
| --- | --- | --- |
| The constraint itself | pool bound = **P_total × f × σ** | Sánchez 2017, PMID 28779005 |
| f, enzyme-in-model per total cellular protein | **0.4461** (measured, PaxDb / Yeast 7.6) | Sánchez 2017, PMID 28779005 |
| σ, saturation | **fitted per regime: 51% / 46% / 44%** (Van Hoek chemostat / CEN.PK chemostat / batch) | Sánchez 2017, PMID 28779005 |
| P_total measured for the proteomics condition | **0.448 g/gDW**, of which measured enzymes **0.283**, leaving **0.165 g/gDW** unmeasured | Sánchez 2017, PMID 28779005 |
| **What `ecYeastGEM_batch` actually uses** | `sigma = 0.5`, `Ptot = 0.5` (source comment: **"Assumed constant"**), `f = 0.5`, `gR_exp = 0.41` → **protein pool = 0.125 g protein/gDW** | **Verified directly in GECKO source**: `geckomat/getModelParameters.m` (v2.0.0) and `tutorials/full_ecModel/YeastGEMAdapter.m` (main). Not a literature value. |
| ecYeast7 scale | 6,741 reactions, 3,388 metabolites, 764 enzymes; median kcat **70.9 s⁻¹**, median MW **48.2 kDa**; central-carbon enzymes faster and smaller (**120 s⁻¹, 41.7 kDa**) than secondary metabolism (**45.1 s⁻¹, 53.6 kDa**) | Sánchez 2017, PMID 28779005 |
| What the enzyme constraint buys | µmax across 12 carbon sources × 3 media: **average relative error 8%**; proteomics integration **reduced flux variability in >60% of reactions**; reproduces the Crabtree switch that Yeast7 cannot | Sánchez 2017, PMID 28779005 |
| **A single pool is the wrong model** | Constraints are **compartment-specific** (cytosol, plasma membrane, mitochondrial matrix and inner membrane) plus ribosome capacity plus a cytosolic density equality. Minimal unspecified-protein fraction **0.245 g/g on glucose, 0.49 on galactose, 0.34 on maltose**. Which constraint binds moves with condition: **glucose transporter below 0.28 h⁻¹**, a **mitochondrial** constraint at ethanol onset, a **cytosolic volume** constraint under sugar excess | Elsemman 2022, PMID 35145105 |
| Apparent turnover numbers are condition-dependent | kapp and ATP maintenance regressed from flux + protein data **vary with oxygen and nutrient availability**; Crabtree is limited by **mitochondrial** proteome capacity and secondarily ribosomes, **not** overall proteome capacity | Dinh 2023, PMID 37080482 |
| σ is not a constant | Yeast allocates **~0.1 g/gCDW** to energy metabolism; apparent saturation rises with glucose uptake and then asymptotes | Chen 2019, PMID 31405984 |
| **The reserve, and it is enormous** | **>50% of capacity held in reserve in seven metabolic superpathways; >80% in glucose metabolism.** Translational reserve >50% for **2,490 of 3,361 expressed genes (74%)**. Ribosomes carry up to **30%** sub-stoichiometric proteins | Yu 2020, PMID 32312967 |
| Model lineage | Yeast8 / ecYeast8 ecosystem; Yeast9 = 1,162 genes, 2,805 metabolites, 4,130 reactions, MEMOTE +27%, growth R² = 0.842, but synthetic-lethality **recall 13.1% / precision 10.0%** | Lu 2019, PMID 31395883; Zhang 2024, PMID 39134886 |

**Two conclusions the project should adopt.** (a) A stress layer that scales one global enzyme
pool is scaling the quantity shown not to bind (PMID 35145105, 37080482). (b) Yu 2020 (PMID
32312967) is the deeper problem: with >50% reserve in most superpathways, **an expression-scaled
upper-bound layer is mostly tightening constraints that were never binding**, which is a
mechanistic explanation for why E-Flux fails to beat pFBA.

### 3.4 Failure modes, with the numbers to design against

| Failure mode | The number | Citation |
| --- | --- | --- |
| Does not beat a no-data control | "All the methods have a lower overall predictive capability compared to pFBA"; yeast: all but E-Flux worse than pFBA | PMID 24762745 |
| …and the same when uptake rates are known | pFBA **0.768** vs E-Flux2 **0.725** vs SPOT **0.650** | PMID 32903284 |
| **E-Flux's scale factor is unidentifiable** | PC **40.06 → 86.56** between two conditions; **PC > 96 makes the constraints vanish** | PMID 34839828 |
| **E-Flux's flux solution is non-unique** | The same input admits solutions scoring **r = 0.35 to 0.92**. A reported correlation with no stated tie-break is not a result. Fixed by adding min-L2 (E-Flux2) | PMID 27327084 |
| Cannot predict growth-rate *increases* | E-Flux only ever scales bounds **down**, so growth is set by the objective; methods with no growth constraint (iMAT, Lee-12, RELATCH*, GX-FBA) predict **no growth at all** | PMID 24762745 |
| Proteomics does not rescue it | "we did not observe any improvements … when we used proteomic rather than transcriptomic data" | PMID 24762745 |
| **mRNA→protein is weak in yeast specifically** | 5,354 mRNAs + 2,198 proteins × 10 conditions: overall correlation **0.46** (0.88 for the 202 differentially expressed); translation efficiency varies **>400-fold**; mRNA explains **61%** of protein-synthesis variance. **And: "only mitochondrial fluxes were positively associated with changes at the transcript level."** | Lahtvee 2017, PMID 28365149 |
| The classic yeast measurement | Equal mRNA → protein varied **>20-fold**; invariant protein with mRNA varying **up to 30-fold** | Gygi 1999, PMID 10022859 |
| In principle | "transcript levels by themselves are not sufficient to predict protein levels in many scenarios" | Liu 2016, PMID 27104977 |
| Hard kcat bounds still admit infeasible solutions | "the true solution space is even smaller" — regulation of activity and substrate under-saturation are not accounted for | PMID 28779005 |
| **The optimality objective itself fails under stress** | Stress proteomes support "metabolic robustness in contrast to optimal protein utilization as a cellular objective" | PMID 35773252 |
| The GEM is the ceiling | Best-achievable error high in 4/5 yeast conditions → iMM904 disagrees with measured fluxes; Yeast9 synthetic-lethality recall 13.1% | PMID 29878053; PMID 39134886 |
| Even signalling+ecModel hybrids are order-of-magnitude at best | mean \|log₁₀\| protein error **1.55 / 2.32** after a one-order improvement | PMID 33836000 |

Three of these argue **for** what `bridge/regulation.py` already does. Module-level rather than
gene-level activity is the right response to the mRNA→protein evidence (PMID 28365149, 10022859,
27104977) and to PMID 35595797's finding that functional-group totals correlate well where
individual genes do not. E-Flux needing no threshold sidesteps PMID 31323017's 38.5%. And the
file's own note that E-Flux cannot lift a lower bound off zero is exactly PMID 24762745's
"cannot predict growth-rate increases", arrived at independently.

## 4. Q4: beta-carotene in *S. cerevisiae* — model parameters

Titres are **not** repeated here; `docs/EXTERNAL_CAROTENOID_BOUND.md` already carries 28
converted rows and the refutation of the 32.84 mg/gDCW ceiling. This section adds what a kinetic
or constraint model needs, and the stress bridge.

**Verdict.** The pathway configuration is well documented. **The enzyme kinetics do not exist**:
no Km and no kcat for *X. dendrorhous* crtE, crtYB or crtI in any host, and **no kcat for
phytoene synthase in any organism at all**. Which step limits is **contested**, and the most
recent yeast-specific kinetic model says **CrtYB**, not crtI. No measured FPP/GGPP split and no
GGPP pool exists for a carotenogenic yeast. What *does* exist, richly, is the bridge: **six
papers report a stress readout and a titre in the same strain**, and they fix the *sign* of the
coupling in both directions.

### 4.1 The stress-and-titre validation candidates — read this first

| Rank | Paper | What it reports together | Why it matters |
| --- | --- | --- | --- |
| **1** | **Bu 2020, PMID 33062054** — "Engineering endogenous ABC transporter with improving ATP supply and membrane flexibility enhances the secretion of β-carotene in *Saccharomyces cerevisiae*" | Titre **+ quantitative proteome + elevated intracellular ROS** in one strain against a matched non-producer. 67.8 → 149.8 mg/L intracellular + 10.1 mg/L secreted. 85 proteins up / 20 down: TCA (Sdh1-3, Lsc1/2, Cit2/3), ATP synthase, **ergosterol (Erg1/4/6/20/24/25)**, lipid (Acc1, Ino1), cell wall (Mpg1, Crh1, Pst1); **Pdr5p 1.20×, Pdr10p 1.26×**. ROS elevated at 36 h. SNQ2 overexpression cost ~27.9% of OD₆₀₀. | **Everything a stress-coupled model needs to fit, in one paper.** |
| **2** | **Verwaal 2010, PMID 20632327** — "Heterologous carotenoid production in *S. cerevisiae* induces the pleiotropic drug resistance stress response" | Genome-wide transcriptome + carotenoid level in **two strains at two production levels**, in glucose-limited chemostats at **D = 0.10 h⁻¹**, 30 °C, pH 5.0. High producer specifically induces **PDR** (ABC + MFS transporters); the low producer shows **no clear genome-wide response** — i.e. a **stress threshold**. `pdr10Δ` transformants had decreased growth rate *and* lower carotenoid. Conclusion: membrane stress plus a secretion response. **Raw arrays are public: GEO GSE8451, 14 Affymetrix arrays (GPL90), full genotypes.** | **The only item where the primary data is obtainable.** Chemostat means mu is fixed, which removes the §1.2 confound *by design*. Reported chemostat contents of 211 and 1752 µg/g come from a secondary summary — **UNVERIFIED**. |
| **3** | **Yan 2011, PMID 21120656** — "Important role of catalase in the production of β-carotene by recombinant *S. cerevisiae* under H₂O₂ stress" | **The cleanest oxidant-dose ↔ titre pair available** (two doses, not a curve). 0.5 and 1.0 mM H₂O₂ significantly *stimulated* β-carotene; in `ctt1Δ` the titre was **16.7% and 36.7% lower** at those doses; no growth difference between `ctt1Δ` and parent. | `healthy_ladder("H2O2")` runs 0.034–0.55 mM, so **0.5 mM is on the ladder and 1.0 mM is the panel's own measured lethal dose** — i.e. this paper reports a *rising* titre at a dose the panel calls lethal, which is itself a testable disagreement. CTT1 is the panel's `STRE-general` marker gene, here used as a genetic control rather than a reporter. |
| **4** | **Ukibe 2009, PMID 19801484** — "Metabolic engineering of *S. cerevisiae* for astaxanthin production and oxidative stress tolerance" | Survival on **1.7 mM H₂O₂**. Ketocarotenoid (astaxanthin) strains **more** tolerant; **β-carotene-only strains were *more sensitive* than non-carotenogenic cells.** β-carotene ≈390 µg/gDW; astaxanthin ≈29 µg/g. Growth rate matched controls. | **Fixes the sign: β-carotene accumulation is pro-oxidant in yeast.** This is the result that stops the model assuming the product is protective. |
| **5** | **Méndez-Alvarez 2000, PMID 10673398** | The opposite-sign control: raised carotenoid rescues `yap1Δ` against H₂O₂, paraquat, menadione and UV, and **diphenylamine inhibition of carotenogenesis removes the protection** — a causality test. | Pair with #4 to bracket the sign. The two together say the sign depends on the carotenoid species, not on carotenoids generally. |
| **6** | **Bu 2017, PMID 29161329** | Whole-metabolome **burden** measurement, producer vs parent: acetate, glycerol, citrate, pyruvate, succinate, fatty acids, **ergosterol** and energy metabolites all lower in the recombinant. 10 g/L acetate in exponential phase → **+39.3%** β-carotene concentration, restoring acetyl-CoA/FA/ergosterol/ATP. | A burden measurement **plus a rescue that reverses it** — the closest thing to a burden coefficient. |
| 7 | **Reyes 2014, PMID 24262517** (already a repo PMID) | Periodic H₂O₂ shocking; **6 → 18 mg/gDCW**; transcriptome of the hyper-producers (lipid and mevalonate biosynthesis up). | Selection under a panel stressor. |
| 8 | **Godara 2019, PMID 31595456** | ALE under stress; YMRCTy1-3 mutation raises carotenoid; **positive correlation between lipid content and β-carotene**. | Tolerance and production co-selected. |
| 9 | UPR ↔ growth ↔ carotenoid: **Nguyen 2022, PMID 36255243**; **Ishiwata-Kimata 2025, PMID 40772766**; **Hao 2024, PMID 39536088** | Constitutive HAC1i expands the ER and raises TAG and carotenoid. Nguyen: HAC1i cells "grew considerably slowly", and grew *faster* under **weak** ER stress. Ishiwata-Kimata: `hda3Δ` restores growth with UPR targets still elevated. Hao: `CLN3Δ` gives +83% cell volume, **+76.9% DCW/cell**, +82% protein, **+41% carotenoid**, with UPR up and framed as *beneficial*. | Directly relevant: the repo's own plate reporters are **UPRE1/UPRE2**. These three disagree about whether UPR activation is a cost, which is itself the finding. |
| 10 | **Gosselin-Monplaisir 2023, PMID 37511557** | Not a dataset — **the assay**. The AOPY method measures ROS/radical scavenging in *living* yeast and was shown to read β-carotene inside living producer cells without lysis. | The way to generate a matched ROS-and-titre sample in-house. |
| — | Adjacent host only: **Park 2022, PMID 35671979** (*Y. lipolytica*) | BHT antioxidant → **>20-fold** retinol; GSH2 overexpression improves further; **<5% of consumed β-carotene recoverable** as retinal/retinol, i.e. severe *in vivo* oxidative destruction of the product. | Mechanism, not a fit target. The same 1% BHT intervention in *S. cerevisiae* is **Lin 2024, PMID 39581972** (readout is vitamin A). |

**A confound that must be read before interpreting any PDR signal.** **Chen 2024, PMID
38324606** — "Transcription Factor Pdr3p Promotes Carotenoid Biosynthesis by Activating GAL
Promoters in *Saccharomyces cerevisiae*" — shows Pdr3p raises carotenoid **not** by secretion but
by directly activating GAL promoter UASs (EMSA + qRT-PCR + transcriptomics). **In a GAL-driven
strain, a PDR signal is therefore not unambiguously stress**; it may be the expression system
talking. Verwaal's strains are TDH3-driven, so item 2 survives this, but any GAL-driven
replication does not.

### 4.2 The one existing kinetic model of this pathway in yeast

**Elizondo 2025, PMID 40891387** — "Complex Kinetic Models Predict β-Carotene Production and
Reveal Flux Limitations in Recombinant *Saccharomyces cerevisiae* Strains". Already a repo PMID.

| Item | Value |
| --- | --- |
| Strains | β-car2/3/4 in CEN.PK2-1c; CrtE/CrtYB/CrtI from *Phaffia rhodozyma* (= *X. dendrorhous*) at loci **XI-5, XI-3, X-2**; promoters **TEF1, PGK1, TDH3**; terminators ADH1, CYC1 |
| Cultivation | Carbon-limited aerobic chemostats at **D = 0.1 and 0.25 h⁻¹** |
| Fluxes, β-car4 @ 0.101 h⁻¹ | β-carotene **185.3 nmol/gDCW/h** [166.4, 201.8]; lycopene **760.7** [647.5, 877.1] |
| Fluxes, β-car4 @ 0.254 h⁻¹ | β-carotene **453.6 nmol/gDCW/h** [361.2, 547.4]; lycopene **392.9** [352.6, 434.6] |
| Method | ABC-GRASP kinetic ensembles; thermodynamics via eQuilibrator/Component Contribution. **Rate constants are Bayesian-sampled, not taken from literature kcat/Km — so the paper contains no reusable constants.** |
| **Flux control** | **CrtYB exerts the highest control at both growth rates.** 67% CrtYB up → **+145.6%** (low mu) / **+84.2%** (high mu) β-carotene flux. "Upregulation of CrtI and CrtE in the simulations had no relevant effect on the pathway." ERG13 also high control; ERG10 discarded. |
| Not reported | intracellular GGPP or FPP; growth penalty; biomass-yield comparison |

These are **the only published flux values for a carotenogenic *S. cerevisiae***, and they come
with mu attached. Note the units: **nmol**/gDCW/h, four orders below the mmol/gDW/h that
`kinetic/carotenoid.py` and `fba/carotenoid.py` work in.

### 4.3 Which step limits — the literature disagrees, and this changes the model

| Position | Evidence | Citation |
| --- | --- | --- |
| **crtI (phytoene desaturase)** | The primary claim, in cyanobacteria: "phytoene desaturation is a rate-limiting step in carotenoid biosynthesis" | Chamovitz 1993, PMID 8349618 |
| crtI, yeast-specific | Phytoene was **86–94% of total carotenoid** in low-crtI strains, and an extra crtI copy was needed. **But qPCR showed crtI was the highest-expressed carotenogenic gene** — so the block is not transcriptional; the authors proposed an unsuitable membrane environment or missing electron carriers | Verwaal 2007, PMID 17496128 |
| crtI, in planta | — | Al-Babili 2006, PMID 16488912 |
| **CrtYB**, and this is the yeast kinetic-model answer | Highest flux control at both dilution rates; CrtI and CrtE upregulation had no relevant effect | Elizondo 2025, PMID 40891387 |
| **crtE** | Named as one of two bottleneck enzymes | Song 2026, PMID 41645597 |
| **Neither — it is compartmentalisation, not turnover** | **Phytoene is stored in lipid droplets; CrtI and CrtYB sit on the ER.** Mitochondrial CrtI was non-functional; mitochondrial CrtYB improved titre; **cytosolic CrtI with its C-terminal membrane anchor deleted increased β-carotene**. 79 mg/gDCW | Arhar 2024, PMID 39215465 |
| Lipid droplets are the major store | Direct LD separation, quantification and disruption; 2 mM oleate +36.4%; IZH1 promoter on ERG9 +31.7% "without adversely affecting cell growth"; final 11.4 mg/gDCW, 142 mg/L | Bu 2022, PMID 34983533 |

**This is a correction to make explicitly.** An earlier reading of Verwaal 2007 alone would put
`vmax_crti` as the binding parameter in `kinetic/carotenoid.py`. The two most recent lines of
evidence say otherwise: flux control sits on the **bifunctional CrtYB** (PMID 40891387), and the
phytoene block is a **localisation** problem (PMID 39215465) that a well-mixed Michaelis-Menten
branch cannot represent at all. **The `crtyb_pool` shared-capacity feature the file already has
is the right structure and is where the evidence now points.** Its most influential unknown —
crtI turnover — has never been measured in any organism.

### 4.4 Enzyme kinetics: what exists

| Enzyme | Km | kcat | Verdict |
| --- | --- | --- | --- |
| **BTS1** GGPPS, *S. cerevisiae* | FPP 0.0032–0.0037 mM; IPP 0.0008–0.0017 mM | **2.5 *or* 0.025 *or* 0.037 s⁻¹** | The only enzyme with a full {Km, kcat} pair. **kcat differs ~100× across three papers from the same lab at identical conditions** (PMID 16554305, 19245203, 23534508) — pick one, document it, never average |
| **crtE, *X. dendrorhous*** | — | — | **Nothing. This is the gene in the strain.** |
| **crtYB, *X. dendrorhous*** (synthase + cyclase) | — | — | **Nothing for either domain — and this is the step PMID 40891387 says controls the flux** |
| **crtI, *X. dendrorhous*** | — | — | **Nothing** |
| ERG20 FPPS, *S. cerevisiae* | IPP 4 µM, DMAPP 8 µM, GPP 14 µM; spec. act. 2.33 µmol/min/mg (PMID 234442) | — | No yeast kcat. Would have to be back-calculated |
| ERG20 F96C / F96W | — | — | No *in vitro* kinetics anywhere |
| Bacterial/plant proxies | crtE *E. uredovora* Km FPP 0.011, GPP 0.009, IPP 0.036 mM (PMID 8215396); IdsA *C. glutamicum* full pair (PMID 25181035); crtB *E. herbicola* S₀.₅ GGPP ≈35 µM with substrate inhibition >100 µM (PMID 12641468); crtI *P. ananatis* Km phytoene **0.0176 mM** (PMID 22745782); crtI *R. gelatinosus* 14.8 mM (PMID 20887710) — **~1000× the *P. ananatis* value, almost certainly mol% in liposome rather than bulk aqueous; do not put these on one axis**; crtY *P. haeundaensis* Km lycopene 3.5 µM (PMID 23412054) — **BRENDA transcribes 35 µM, 10× the source; use the paper** | — | Proxies only. **No kcat exists for phytoene synthase in any organism.** |

So all six parameters in `kinetic/carotenoid.py` would be fitted, not sourced. Taking them from
GECKO's automated kcat matching (PMID 35773252) would put a database default behind a claim about
this pathway, which is what `docs/research/CIRCULARITY.md` exists to prevent. **The correct action is to
keep the module parked and record why.**

### 4.5 Flux partitioning, pools and growth — the honest state

- **No measured FPP/GGPP split** between ergosterol and the carotenoid branch in a carotenogenic
  *S. cerevisiae*, and **no measured GGPP pool in any carotenogenic yeast**. The engineering
  literature manipulates the split (ERG9 promoter swaps, PEST/ERAD degrons) and reports titre:
  +31.7% (PMID 34983533), +42.3% lycopene (PMID 38621758), +86% nerolidol (PMID 27939849).
- **Nearest pool substitute:** Rubat 2017, PMID 28854674 measures intracellular DMAPP/IPP, GPP
  and FPP (nM/OD) across ERG20 A99X variants with mu for each — WT BY4741 **mu_max 0.17 h⁻¹**,
  variants 0.02–0.14, WT split ≈**38% GPP / 62% FPP**. **GGPP was not measured, and the strain is
  not carotenogenic.** Values came from one automated full-text read — **check the figures before
  use**.
- **Best-quantified branch-point cost:** Peng 2017, PMID 27939849 — over-supplying FPP pushed
  **squalene to 1% of biomass and the specific growth rate declined**; destabilising Erg9p
  returned squalene to wild type **with no growth penalty**. Sesquiterpene, not carotenoid, but
  the cleanest partition ↔ growth-cost pair in yeast.
- **BTS1** raises the GGPP pool but Jiang 1995, PMID 7665600 is genetic — **no kinetic constants
  in it; do not cite it for any**. GGPP-node dose-response libraries: PMID 29118396 (11 promoters
  × 10 combinations, GGOH 18.45 → 161.82 mg/L), PMID 31689467 (ERG20 F96C as a GGPP synthase).
- **Growth and burden.** Only three published mu values are usable for producers, all dilution
  rates: **D = 0.10 h⁻¹** (PMID 20632327) and **D = 0.10 and 0.25 h⁻¹** (PMID 40891387). Verwaal
  2007 (PMID 17496128) and Bu 2022 (PMID 34983533) say growth was unaffected; **Yan 2012, PMID
  22080204 states outright that "the cell growth was inhibited by the metabolic burden induced by
  the production of heterogeneous β-carotene"**; Verwaal 2010 reports reduced growth only in
  `pdr10Δ`. **No paper reports a measured mu for a high-titre carotenogenic strain against its
  parent.** That is a gap, not a search failure.

**Confirmed absent by searching:** no paper reports HSP/heat-shock gene expression alongside a
carotenoid titre in *S. cerevisiae*; no glutathione (GSH/GSSG) measurement in a carotenogenic
*S. cerevisiae*; no measured GGPP pool in any carotenogenic yeast.

---

## 5. What this means for the model we are building

### 5.1 Findings that CONTRADICT a parameter or structure in `generator/stress_panel.py`

These are the high-value items. Line references are pointers, not identifiers.

| # | What the file does | What the literature says | Citation | Severity |
| --- | --- | --- | --- | --- |
| 1 | `MODULES["ESR"]` is Msn2/4 → STRE, **induction only**. There is no repressive module, no PAC/RRPE element, and no RP/RiBi reporter anywhere in the panel. | The ESR is **~300 induced + ~600 repressed**. Two-thirds of it is the repressed arm, which runs through **Dot6/Tod6** and **Sfp1**, not Msn2/4. Dot6 carries almost as much information as Msn2. Dot6 activation *positively* predicts post-stress growth recovery. | Ho 2018, PMID 30078561; Lippman 2009, PMID 19901341; Marion 2004, PMID 15353587; Granados 2018, PMID 29784812; Bergen 2022, PMID 36350693 | **Highest.** The missing arm is the one that couples to growth, which is the confound this repo already carries. |
| 2 | One scalar `ESR` module, driven by every stressor through a weight. | "Msn2p and Msn4p play nonredundant and condition-specific roles in gene-expression regulation, **arguing against a generic general-stress function**." | Berry 2008, PMID 18753408 | High. A single shared scalar is the structure this refutes. |
| 3 | `module_response` applies **one static saturating Hill (n = 1) in dose** per module. | Hog1 nuclear enrichment **perfectly adapts** — steady-state output is independent of steady-state input — and the transcriptional arm is the *slow* feedback. Dose is encoded in pulse **frequency** under glucose limitation and in **amplitude** under oxidative stress; one functional form cannot serve both. | Muzzey 2009, PMID 19596242; Mettetal 2008, PMID 18218902; Hao 2011, PMID 22179789; Cai 2008, PMID 18818649 | **High for the osmotic and carbon arms specifically.** `NaCl`, `sorbitol` and `glucose_starvation` are the affected entries. |
| 4 | `viability()` multiplies every module alike — "a cell losing the capacity to transcribe cannot spare one regulon." | Under stress translational capacity is **redistributed** toward induced transcripts, not uniformly scaled; blocking rESR repression *delays* induced-protein production. | Ho 2018, PMID 30078561; Bergen 2022, PMID 36350693 | High. A uniform multiplier is the opposite of the measured mechanism. |
| 5 | `growth_rate(stressor, dose) = mu_max * viability(dose)` — growth is a deterministic function of dose, so it carries no information of its own. | Causality runs **both** ways: slow growth alone induces heat resistance and stress tolerance, and >25–50% of the transcriptome co-varies with growth rate with >80% overlap with the stress genes. | Lu 2009, PMID 19056679; Zakrzewska 2011, PMID 21965291; Brauer 2008, PMID 17959824; Regenberg 2006, PMID 17105650; O'Duibhir 2014, PMID 24952590 | High. Growth must be a co-varying latent, not a derived quantity. |
| 6 | `reporter_loadings` assigns 24 modules distinct loading vectors over 24 promoter fusions. | **60–90% of promoter activity change across conditions is a single global scaling factor** that depends only on the condition, not the promoter — measured on ~900 yeast fluorescent reporters. A single Msn2 target gene carries **1.2–1.3 bits**; two carry **1.67–1.83**. | Keren 2013, PMID 24169404; Hansen 2015, PMID 25985085 | **Highest, jointly with #1.** The panel's effective rank is near 1, not 24. This belongs in `docs/research/IDENTIFIABILITY.md`. |
| 7 | One `ec50` per stressor-module pair, fixed. | Response to the same dose depends on the cell's stress **history** through non-overlapping genetic routes, and cross-protection capacity varies heritably between strains. | Berry 2011, PMID 22102822; Stuecker 2018, PMID 29649251 | Medium. Already noted in `PARAMETER_SOURCES.md`; the new part is that the routes do not overlap at all. |
| 8 | Dose flows one way: agent → module. There is **no path by which the product changes a stressor's effective dose.** | β-carotene accumulation makes yeast **more** H₂O₂-sensitive; raised carotenoid **rescues** `yap1Δ` against H₂O₂, paraquat, menadione and UV, and inhibiting carotenogenesis removes the protection. The sign depends on the carotenoid species. | Ukibe 2009, PMID 19801484; Méndez-Alvarez 2000, PMID 10673398 | **High for the product half of the project.** The two halves cannot be joined until a product can modify a dose. |
| 9 | Nothing in the panel represents **where** a molecule is. | Phytoene sits in lipid droplets while CrtI and CrtYB sit on the ER; deleting CrtI's membrane anchor to make it cytosolic *increased* β-carotene. The phytoene block is a localisation problem, not a turnover one. | Arhar 2024, PMID 39215465; Bu 2022, PMID 34983533 | Medium for the panel, **high for `kinetic/carotenoid.py`**, whose well-mixed Michaelis-Menten branch cannot express it. |

### 5.2 Findings that SUPPORT a parameter already in the file

- `_GENERAL_POTENCY = 2.5` — the ESR trailing the specific regulon. Granados 2018 (PMID
  29784812) independently supports the *ordering*: generalists (Msn2/4, Dot6, Tod6, Maf1,
  Sfp1) encode stress identity "**only if stress is high**", while specialists (Hog1, Yap1,
  Mig1/2) encode over "a **wider range of magnitudes**" and faster. The ordering is now
  evidenced; the factor 2.5 remains asserted.
- `_MU_MAX = 0.40` /h. Within the band the file already brackets; nothing found against it.
- The `Kind.RATIOMETRIC` / `Kind.TRANSCRIPTIONAL` split, and the growth-dilution time constant
  `1/(mu + k_deg)`. Mettetal 2008 (PMID 18218902) reinforces it: the fast feedback needs no
  protein synthesis and a promoter fusion cannot see it at all.
- **`bridge/latent_bridge.py`'s entire mechanism.** The file maps latent stress activity onto the
  **ATP maintenance** reaction (`r_4046`), a choice it makes without a literature citation.
  Lahtvee 2016 (PMID 27307591) is that citation: ethanol, salt and temperature stress, three
  graded levels each, in chemostats "**in which specific growth rate-dependent changes are
  eliminated**", concluding that "**increased ATP demand for cellular maintenance underpins a
  general stress response and is responsible for the onset of overflow metabolism**". Sánchez
  2017 (PMID 28779005) reproduced 38 °C growth in ecYeast7 purely by raising NGAM. **Raised
  maintenance ATP is the mechanism the yeast stress literature actually supports, and it is the
  one this repo already picked.** `_MAINTENANCE_PER_ACTIVITY = 300.0` remains asserted, but it
  now has a published dataset it could be calibrated against rather than guessed.

### 5.3 Rules that could become constraints or tests

| Rule | Number | Use it as |
| --- | --- | --- |
| ESR partition 300 iESR / 600 rESR | 1:2 | A **structural prior** on a two-arm latent state: one induced, one repressed, opposite signs, with the repressed arm carrying the growth coupling. |
| Growth-rate-linear transcriptome | >25% of genes, linear, nutrient-independent | A **test**: regress the inferred latent state on measured mu. If R² is high, the state is growth rate. This is the falsification test the model currently lacks. |
| Overlap of growth-regulated and stress-regulated genes | >80% of *genes* | A **budget, with a caveat**: only ~20% of ESR genes move independently of growth rate. Converting a gene-count overlap into a variance budget is an inference, not a measurement — state it as one. |
| Excess ribosomal protein | ≈8% of proteome, rising as mu falls | An **allocation constraint** for the enzyme-budget layer, and the closest thing to a yeast φ_ribosome law. |
| Information per reporter | 1.2–1.3 bits; 1.67–1.83 for two | A **hard identifiability bound**. Convert to distinguishable dose levels (~2–3) and compare against the number of ladder rungs `healthy_ladder` spends. If the ladder has 5 rungs and the reporter carries 2.4 levels, three rungs are wasted. |
| Global scaling factor across promoters | 60–90% | A **null model for the panel**: fit rank-1 first, and only claim module resolution for the residual. `analysis/nulls.py` is where this goes. |
| Perfect adaptation of Hog1 | steady-state output independent of dose | A **refusal**: the osmotic module must not be given a sustained dose-graded activity. |
| Stress → raised maintenance ATP → respiration saturates → overflow | respiration linear in glucose uptake only to **~3.4 mmol/gDW/h** | A **constraint and a test** for `latent_bridge.py`: the twin should place the respirofermentative switch where Lahtvee 2016 (PMID 27307591) places it, and `_MAINTENANCE_PER_ACTIVITY` should be fitted to that dataset rather than asserted. |
| Metabolic and translational reserve | **>50%** in 7 superpathways, **>80%** in glucose metabolism, **74%** of genes with >50% translational reserve | A **refusal**: do not expect an upper-bound-scaling layer to change the flux solution. Report which bounds actually became active; if none did, say so. |
| E-Flux's proportionality constant | **40.06 → 86.56** between two yeast conditions; **> ~96 = inert** | A **required diagnostic**: report the fitted scale and whether the bounds bind, every time. |
| pFBA and raised-NGAM controls | pFBA r = **0.8337** on 9 yeast conditions | **Mandatory comparison arms.** A regulation layer scored without them has not been shown to add anything. |

### 5.4 The one structural recommendation

Make the latent state **two-dimensional and jointly identified**: a stress axis and a growth
axis, with the rESR arm loading on both. The literature makes every single-axis option
indefensible — a stress-only state is confounded by >80% overlapping growth signal, and a
growth-corrected state throws away the two-thirds of the ESR where the shared variance lives.
Then add one repressive reporter (an RP or RiBi promoter, or the ribosome-production biosensor
already built and genome-integrated in Torello Pianale 2021, PMID 35069506) so the second axis
is measured rather than inferred.

### 5.5 The validation route, which now exists

The repo's stated blocker is that "no strain produces the product", so Level 3 is unreachable
(`docs/CONTRACT.md`, `docs/FINDINGS.md`). Published datasets partly lift that, because they
contain a stress readout and a titre in the same strain:

1. **Verwaal 2010 (PMID 20632327)** — chemostat at D = 0.10 h⁻¹, so mu is fixed rather than
   inferred, and the arrays are public (**GEO GSE8451**, 14 arrays, GPL90, full genotypes).
   A **held-out test of the regulation layer**: given a high-producing configuration, the twin
   should predict elevated `xenobiotic` (Pdr1/Pdr3 → PDRE) activity, and given the low producer
   should predict **no genome-wide response** — the paper reports exactly that threshold. That is
   a **binary, pre-registerable prediction**, which is what Tier 3 requires and what
   `docs/CLAIM_BOUNDARY.md` says nothing currently reaches. Verwaal's strains are TDH3-driven, so
   the Pdr3p/GAL-promoter confound (PMID 38324606) does not apply — but it would to any
   GAL-driven replication.
2. **Bu 2020 (PMID 33062054)** — titre, quantitative proteome and intracellular ROS in one
   strain against a matched non-producer. The richest single fit target, and the one that can
   score module activity against measured protein rather than against transcript.
3. **Yan 2011 (PMID 21120656)** — two H₂O₂ doses (0.5, 1.0 mM) with titre and a `ctt1Δ`
   genetic control. `healthy_ladder("H2O2")` tops out at 0.55 mM, so 0.5 mM is on the ladder
   and **1.0 mM is where the panel puts the lethal dose while this paper reports a rising
   titre** — a direct, testable disagreement with `STRESSORS["H2O2"].lethal_dose`.
4. **Ukibe 2009 (PMID 19801484) and Méndez-Alvarez 2000 (PMID 10673398)** — the sign test.
   β-carotene accumulation makes yeast *more* H₂O₂-sensitive; raised carotenoid rescues `yap1Δ`.
   The model must be able to express **both**, and `stress_panel.py` currently expresses neither:
   there is no path by which a product changes a stressor's effective dose.
5. **Pan 2018 (PMID 29643937)** — acetic acid at 91.6 mM, inside
   `healthy_ladder("acetic_acid")`, with CTT1, a ROS level, a growth readout and a titre.
6. **Reyes 2014 (PMID 24262517)** — H₂O₂ dose, titre 6 → 18 mg/gDCW, transcriptome.

Item 1 is the cheapest genuine Tier 3 opportunity in this document.

---

## 6. What nobody has established

This is the section that decides what is worth attempting. Each item is an absence found by
looking, not an absence assumed.

| # | The gap | Why it is genuinely open | Worth attempting? |
| --- | --- | --- | --- |
| 1 | **A yeast proteome-allocation law for stress.** There is a fitted allocation law for *growth* (PMID 35595797, 28857745, 32817546) and a bacterial template (PMID 21097934), but no fitted φ_stress(dose) trading against φ_ribosome in *S. cerevisiae*. The mechanism is established (PMID 30078561) and the coefficient is not. | The measurement exists in pieces — Xia 2022 has absolute proteomes across growth rates, and chaperone fraction *decreases* linearly with growth rate in it. Nobody has done the stress axis. | **Yes.** This is the highest-value gap. It is a regression on published data, not a new experiment. |
| 2 | **A dose-to-response law of any functional form.** No canonical Hill, no fitted threshold, no biphasic parameterisation for yeast stress-gene induction against dose. What exists is per-pathway kinetic ODE models (PMID 16025103) and single-cell dynamic encoding (PMID 22179789). | The literature explains the absence structurally: dose lives in temporal features, the osmotic branch discards it (PMID 19596242), and reporters carry identity but not intensity (PMID 25985085). | **Only as a negative result.** Attempting a dose law and reporting that it is unidentifiable is publishable and honest; attempting it and reporting a fitted EC50 is not. |
| 3 | **Reporter-panel-driven latent state → time-varying yeast GEM constraints.** PROM did state→bounds in *E. coli* from transcriptomes (PMID 20876091); Granados did panel→state in yeast without metabolism (PMID 29784812); nobody has joined them in yeast. | Both halves are published; the join is not. | **Yes — this is the defensible novelty claim.** Frame it as the coupling, not the latent state. |
| 4 | **A formal identifiability treatment of a stress state from a promoter-fusion panel.** | NCA (PMID 14673099) owns the general criteria; Keren 2013 (PMID 24169404) supplies the 60–90% rank-1 result that makes the yeast case hard; nobody has combined them into a bound on how many modules a panel of *n* fusions can resolve. | **Yes,** and it is cheap: it is arithmetic on published numbers plus this repo's own plates. |
| 5 | **A yeast digital twin carrying an explicit latent stress state.** Yeast hybrid fermentation twins exist without one (PMID 40849368); soft sensors estimate physical, not stress, state (PMID 33716616). | Nothing found. | Yes, but it is the *last* step, not the first. |
| 6 | **Whether the ESR is stress-specific or growth-derived, settled.** PMID 24952590 and 22456505 say mostly growth/cell-cycle; PMID 30078561 says not reducible to it. Both are careful papers. | The field has not converged. | **Do not attempt to settle it.** Build a model that is correct under either, i.e. a joint stress-and-growth state (§5.4). |
| 7 | **Enzyme kinetics for the carotenogenic enzymes.** No Km and no kcat for *X. dendrorhous* crtE, crtYB or crtI in any host, and **no kcat for phytoene synthase in any organism**. | Which step limits is itself contested (PMID 8349618 / 17496128 vs 40891387 vs 41645597 vs 39215465), so the missing constant is also the decisive one. | **As a wet-lab item, yes** — but it is enzymology, not modelling. Until then `kinetic/carotenoid.py` stays parked, and stating that is the deliverable. |
| 8 | **A quantitative burden law for the carotenoid pathway.** The product both *imposes* stress (membrane/PDR, PMID 20632327) and *relieves* it (antioxidant, PMID 24262517, 29643937). No paper puts a coefficient on either arm. | Direction established both ways; magnitude neither. | **Yes,** and it is the one place where this project's two halves are genuinely one problem rather than two glued together. |
| 9 | **A measured growth rate for a high-titre carotenogenic strain against its parent.** Only three usable mu values exist for producers, all chemostat dilution rates (PMID 20632327, 40891387). Papers split on whether burden costs growth at all: PMID 17496128 and 34983533 say unaffected, PMID 22080204 says growth was inhibited by the burden. | Nobody has run the matched comparison. | **Yes, and it is the cheapest wet-lab item in this document** — one growth curve pair. |
| 10 | **A measured FPP/GGPP split, or any GGPP pool, in a carotenogenic yeast.** | The engineering literature manipulates the split and reports titre. The nearest pool measurement (PMID 28854674) is not carotenogenic and did not measure GGPP. | Yes; it is what would unblock `fba/carotenoid.py`'s branch-point assumption. |
| 11 | **A stress-conditional enzyme budget.** Enzyme-constrained yeast models establish the budget under nutrient conditions (PMID 28779005, 35773252, 35145105, 37080482); none parameterises how the budget or its compartment split moves under a stressor. | Dinh 2023 shows the parameters *do* vary with oxygen and nutrients, so they should vary with stress; nobody has measured it. | Yes, and it is the natural home for finding #1. |

### 6.1 What is NOT worth attempting

- **Claiming a first low-dimensional representation of the yeast stress response.** That is
  2000–2013 work (PMID 10963673, 12429058, 12740579, 14611662, 24169404).
- **Claiming a first latent-state-to-flux-bounds coupling.** PROM, 2010 (PMID 20876091).
- **Fitting a single EC50 per stressor-module pair and treating it as a parameter of yeast.**
  The response is history-dependent (PMID 22102822), strain-dependent (PMID 29649251), and
  cell-to-cell heterogeneous (PMID 29240790, 36350693).
- **A 24-module panel from 24 promoter fusions.** Keren 2013 (PMID 24169404) and Hansen 2015
  (PMID 25985085) together put the achievable rank far below 24.
