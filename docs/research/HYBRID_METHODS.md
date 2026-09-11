# Hybrid ML + dFBA: the six admitted limits, and what to take from each

The state of the art in hybrid machine-learning / dynamic-FBA bioprocess modelling is COSMIC-dFBA
(Gopalakrishnan et al., *Metabolic Engineering* 82:183–192, 2024, doi 10.1016/j.ymben.2024.02.012,
PMID 38387677). A classifier predicts the fraction `f` of cells in a production rather than a growth
state from bioreactor metabolite concentrations; two state-specific models with hand-assigned
objective priorities are blended as `v_net = f·v_p + (1−f)·v_g`; the result is 90 % and 72 % better
than standard dFBA on cell density and antibody titre, and within 10 % of measurement on both.
Those four numbers are verified from the abstract. The paper body is paywalled and both bioRxiv
HTML and PDF return 403, so the classifier metrics (F1 0.731, Matthews 0.454, 94/130 timepoints
within 0.1) and the Fig 5B amino-acid signs are taken from the brief and are marked **paywalled**
where they carry an argument.

Compiled 2026-08-26. Every identifier below was resolved against Crossref, NCBI E-utilities or
Europe PMC and its title read back; anything I could not resolve is marked UNVERIFIED and is not
relied on. Repo numbers quoted here come from `docs/FINDINGS.md`, `docs/CLAIM_BOUNDARY.md`,
`docs/research/PREDICTION.md` and `src/ystwin/fba/solver.py`.

---

## Verdict table

| # | The admitted limit | Does the literature already solve it? | Verdict |
| --- | --- | --- | --- |
| L1 | Cannot determine the *cause* of the state shift | **Partly. Built, once, in a GEM — for a soil bacterium, on simulated labels.** Cause classifiers from omics are mature and accurate. Nothing closes the loop. | **Adapt** — the cheap version is free: read the cause off the LP dual |
| L2 | Objectives are specified in advance, not learned | **Yes, and the answer is negative.** A COSMIC-dFBA co-author published the refutation in 2026: objectives are non-identifiable and the constraints do the work. | **Ignore** |
| L3 | The state classifier is weak | **Yes, and cheaply.** Put the state in the filter instead of in a classifier. No published head-to-head exists, but this repo's own diagnostics already name the defect. | **Adopt** |
| L4 | Worse than baseline on nine amino acids | **Yes, fully explained, with numbers.** Six of the nine are the documented overflow-catabolism set. This is not an ML failure. | **Adopt as a guard**, not as a fix |
| L5 | Eight conditions, one CHO line, one product | **No. The whole field validates on one cell line.** Genuine transfer exists only without a GEM, or only in bacteria. | **Ignore** as method, **adopt** as claim discipline |
| L6 | No proteome / enzyme allocation | **Largely yes, for the shift that matters.** A fixed growth objective plus a compartment protein budget already predicts Crabtree and overflow. | **Adopt** — this is the architecture change |

---

## L1. It cannot determine the cause of the state shift

**Verdict.** Cause-resolved state inference exists and is accurate, but not as one thing. Classifying
which environmental cause a cell is under, from expression or from a reporter panel, is a mature
problem solved to 70–98 % balanced accuracy. Making a metabolic model *respond* to a named cause is
older still — regulatory FBA has done it since 2001. What nobody has built is the join: an inference
that reads which cause is acting from observable data and then changes the flux prediction
accordingly. One 2026 paper comes closest and gets three quarters of the way, in a GEM, citing
COSMIC-dFBA as its motivation — but its cause labels come from the model's own shadow prices over a
sampled design space rather than from measured data, the organism is a soil bacterium, and the
diagnosis does not feed back into the simulation. So the gap is real, but it is narrower and more
specific than "nobody has done this". State it as: **no cause-conditioned state mixture has been
built for any organism, and no cause inference for any bioprocess runs on measured data.**

**The one paper that does it inside a GEM.** Kim, Choi & Kim, "Causal AI digital twin for bioprocess
bottleneck diagnosis via metabolic flexibility and rigidification maps", *iScience* 29:116820, 2026,
doi 10.1016/j.isci.2026.116820, PMID 42502377. Verified from the record and abstract: regimes are
labelled by the **active-constraint set** over a Latin-hypercube design space, intracellular state is
encoded as targeted FVA *widths* across 30 reactions/modules, and XGBoost + SHAP classifies the
regime, on *Stenotrophomonas maltophilia* SO-1 anchored only by endpoint OD600 at 32 h. The
classifier scores (macro-F1 0.991; cross-organism transfer to *E. coli* iML1515 at macro-F1 0.972)
and the sentence naming COSMIC-dFBA are reported but the article is not open access here —
**UNVERIFIED**. The mechanism is what matters and the abstract confirms it: the cause is read off
which constraint is binding, which is information the LP already produces.

**Cause attribution from LP duals is 34 years old.** Savinell & Palsson gave the shadow-price
formalism and applied it to hybridoma metabolism in 1992 (PMIDs 1593896 and 1593897, *J Theor Biol*
154:421–454 and 455–473), and Chen et al. use dual prices in a CHO GEM explicitly —
"the corresponding dual prices provide fruitful information concerning coupling relationships
between nutrients" (doi 10.1038/s41540-019-0103-6, *npj Syst Biol Appl* 5:25, 2019). Nobody has
turned that into a state classifier.

**Which-stress classification, where the numbers are.**

| Work | Input | Causes separated | Score |
| --- | --- | --- | --- |
| Kim, Zorraquino & Tagkopoulos 2015, doi 10.1371/journal.pcbi.1004127 | *E. coli* expression | strain, growth phase, medium, oxygen, antibiotic, carbon source | balanced accuracy 70.0 ± 3.5 % to 98.3 ± 2.3 %; +10.6 ± 1.0 % from a joint model |
| Elad et al. 2008, doi 10.1021/es801489a | five stress-promoter *lux* reporters | 5 toxicants + control | causative agent at ≤3 % error in 30 min, zero false negatives, in tap water and wastewater |
| Sastry et al. 2019, doi 10.1038/s41467-019-13483-w | *E. coli* expression | ICA components (iModulons) whose named factors *are* candidate causes | — |

Elad 2008 is the important one for this project: a small panel of stress-promoter reporters, read
online, identifying which agent is acting. That is the architecture already in `analysis/latent.py`
and `analysis/stress_model.py` — a latent basis from readings alone, then a supervised map to named
modules.

**The mathematical shape is proven elsewhere.** Attributing an observed response to competing named
causes by non-negative factorisation is field-proven at scale: Alexandrov et al., "Signatures of
mutational processes in human cancer", *Nature* 500:415–421, 2013, doi 10.1038/nature12477 — >20
signatures over 4,938,362 mutations, each attributed to a named process. Nobody has ported it to
bioprocess stress attribution.

**Two pieces of counter-evidence that should temper the claim.**

- Lahtvee et al., doi 10.1091/mbc.E16-03-0187 (*Mol Biol Cell* 27:2505–2514, 2016). Chemostats, so
  the growth-rate confound is removed: ethanol, salt and temperature stress in yeast **converge on
  mitochondrial metabolism**, and the shared response is an increased ATP demand. Different causes
  look the same at the maintenance-ATP and overflow level. So cause is poorly identifiable from
  fluxes alone — COSMIC-dFBA's single-`f` assumption is not obviously wrong on the channel it reads.
- O'Duibhir et al., doi 10.15252/msb.20145172 (*Mol Syst Biol* 10:732, 2014). A large share of an
  apparently stress-specific expression signature is a cell-cycle-distribution artefact of slower
  growth. A cause classifier that does not remove growth rate first learns "slow", not "why slow".
  `analysis/latent.py` already takes dilution-corrected activities, which is the right input.

**The reverse direction — cause into the model — is old and works.** rFBA (Covert, Schilling &
Palsson, doi 10.1006/jtbi.2001.2405) and PROM (doi 10.1073/pnas.1005139107) both key the model's
constraints to which condition is present. Ghodba et al. (doi 10.1016/j.jbiotec.2025.08.010,
*J Biotechnol* 408:61–71, 2025) do it in a CHO dFBA: kinetic constraints as explicit functions of pH
and temperature, on 20 Ambr250 fed-batches, R² ≥ 0.97 on growth and titre. Cause-*conditioned*, not
cause-*inferring*. It is the complementary half.

**And COSMIC-dFBA's own group took the first step without finishing it.** Gopalakrishnan et al., doi
10.1016/j.ymben.2024.07.009 (*Metab Eng* 85:94–104, 2024), five months later: multi-omic
characterisation found "topological changes in peripheral metabolic pathway expression associated
with phase shifts" and nutrient-uptake bottlenecks invisible to exometabolomics alone. They
characterised the phase shifts with omics and did not turn it into a cause classifier or feed it back.

**What we should do: adapt, and take the free version first.** Do not build a cause classifier from
new data. Two things, in order.

1. **Read the cause off the dual.** Every solve this project already performs produces shadow prices.
   Recording which exchange bound has the largest positive shadow price at each timestep costs
   nothing, is not a heuristic, and is the most authoritative available signal for which constraint
   is binding. It gives a cause label with no new measurement and no new model.
2. **Keep the latent-to-named-module map, but hold it to the existing gate.** `analysis/stress_model.py`
   already does the Elad-2008 architecture. `bridge/latent_bridge.py` is the join to a constraint, and
   G4 returned INCONCLUSIVE for all four constructs, so the branch does not yet earn its place. The
   novelty claim survives — a *cause-conditioned* state mixture has been built for no organism — but
   it must be stated against Kim 2026, not against silence, and it must not be claimed until the
   latent branch predicts an excluded anchor.

---

## L2. Objectives are defined in advance, not learned

**Verdict.** Do not learn the objective. This is not an open opportunity; it is a closed question with
a negative answer, and the closing paper was written by a COSMIC-dFBA co-author. Objective learning
has a 23-year literature with methods that work, but every method that recovered a known objective
needed a near-complete intracellular flux map from ¹³C MFA, and the 2026 state of the art concludes
that objectives are substantially non-identifiable and that predictive accuracy is driven by the
constraints rather than by the objective. Independently, the yeast fermentation literature ran the
experiment that settles it for a project like this one: replacing five hand-specified phase
objectives with a single continuous one changed nothing.

**What each method needed, and whether it worked.**

| Method | Organism, network | Data required | Recovered a known objective? |
| --- | --- | --- | --- |
| ObjFind, doi 10.1002/bit.10617 (2003) | *E. coli* | ¹³C isotopomer flux map, 2 conditions | Yes; biomass coefficient 9× (aerobic) and 15× (anaerobic) the next largest |
| BOSS, doi 10.1186/1471-2105-9-43 (2008) | yeast central metabolism, 62 reactions | full flux map over all 62 reactions | Yes; tolerated ~15 % flux noise. Bilevel duality made the problem non-convex; needed multi-restart |
| invFBA, doi 10.1186/s13059-016-0968-2 (2016) | *E. coli*, *S. oneidensis* | **complete** flux distribution | Yes on noise-free simulated data. **No** on real ¹³C data from six evolved *E. coli* strains — biomass maximisation was not in the inferred objective space for any of them. Recovery is "highly reduced" between 1 % and 10 % noise |
| CellTarget, doi 10.1038/s41540-026-00700-8 (2026) | CHO, CHOmpact (144 rxns) and iCHO2441 | ¹³C MFA, 4 conditions | Good forward flux prediction; beat an inverse-FBA benchmark. Objective **non-identifiable** |

**CellTarget is the paper that decides this.** Monteiro et al., *npj Systems Biology and Applications*
12:79, 2026, PMID 41932911. Senior author Cleo Kontoravdi, who is a co-author of COSMIC-dFBA. Direct
quotes from the full text:

> "The results reveal substantial non-identifiability of cellular objectives, with predictive accuracy
> driven more strongly by network constraints than by objective specificity."

> "Despite comparable flux predictions, the inferred cellular objectives differ between model sizes for
> the same experimental data. This highlights a lack of identifiability of the cellular objective that
> persists across network scales."

> "The lack of convergence toward a unique coefficient structure further suggests that cellular
> behaviour cannot be explained by optimisation of a single dominant objective. Moreover, the flux
> constraints appear to play a more decisive role than the objective itself."

> "When expanding the dimension of cellular objective, the limiting factor appears to be the
> availability of data capable of resolving metabolic priorities beyond a low-dimensional set."

**And the yeast head-to-head says the objective machinery is not load-bearing.** Two mSystems papers
on the same *Saccharomyces* batch fermentations. MPMO (doi 10.1128/msystems.00260-21, 2021) uses five
phases with phase-specific objectives — the direct yeast analogue of COSMIC-dFBA's state-specific
priorities. IMC (doi 10.1128/msystems.01615-24, 2025) replaces all of it with one continuous
time-varying objective. Result: mean R² **0.93 for IMC versus 0.91 for MPMO**, and a
Kolmogorov–Smirnov test "did not reveal statistical evidence that both R-squared value distributions
differed at any reasonable confidence level (P-value ≈ 0.7)". Removing the hand-specified phase
objectives entirely made no measurable difference. Neither model is enzyme-constrained — which is
where L6 picks this up.

Static single objectives *do* fail on yeast fermentation, so this is not an argument for one fixed
objective either. IMC tested them: "maximization of biomass resulted in no growth during the
exponential phase, maximizing ATP resulted in lack of growth". Objective selection still matters;
objective *learning* does not pay. For yeast specifically, Schnitzer et al.
(doi 10.1371/journal.pone.0276112) found max-growth essential and the parsimonious solution or an
added growth-independent energy cost the improvement that helps — which is the pFBA-plus-maintenance
default, not a learned objective. Schuetz et al. remain the practical selection guide: 11 objectives
× 8 constraints against ¹³C fluxes over 6 conditions, no single objective works everywhere
(doi 10.1038/msb4100162), with metabolism sitting near a three-dimensional Pareto surface across nine
bacteria (doi 10.1126/science.1216882).

**What we should do: ignore.** This project has ~a dozen extracellular metabolites and tens of
timepoints, no ¹³C flux map, and assay noise far above the 1 % where invFBA still recovers a known
objective. Learning an objective here would produce a fitted coefficient vector that is not
identifiable and would be indistinguishable from tuning. Keep max-growth plus a growth-independent
ATP maintenance term, and spend the effort on the constraint side. Cite CellTarget as the honest
answer to the stated limitation, including that it is a negative result from inside the same group —
that is a stronger point in a survey than presenting L2 as open.

---

## L3. The state classifier is weak

**Verdict.** The classifier is weak because it is a classifier. F1 0.731 and Matthews 0.454
(**paywalled**) is a per-timepoint discriminative label fitted at n = 130, on a quantity that is
continuous, latent, temporally smooth, and already governed by an ODE the model solves anyway.
Nothing about that problem calls for a better classifier; it calls for the state to live in the
filter. This project already has the filter. The right move is not to add a model but to move the
state into one that exists — and this repo's own diagnostics already say so, by an independent route,
before any of this literature is consulted.

**What the repo already measured.** From `docs/research/PREDICTION.md` and `docs/CLAIM_BOUNDARY.md`:
conditioning `estimator.py` on the first 2 h of 84 real wells and forecasting the remaining 12
timepoints, the twin **loses to a six-point log-linear extrapolation** on OD by 23 % of level, and
covers 40–49 % of the time with a nominal 95 % interval. The cause is already named:
`_propagate` gives growth rate no drift, so a decelerating culture is extrapolated flat.

The innovation-calibration test this section used to cite alongside that loss is not currently
available. The 167-well medians of 1.76 / 4.44 against [0.52, 1.63] belong to an unarchived
2000-particle research probe — `data/current_claims.json` carries it as `nis.research.shipped`,
role `historical`, binding `pending`. Run on the four committed plates, the shipped CLI refuses
instead: `outputs/nis_channel_summary.csv` reports `tested_well_channels: 0` and verdict
INCONCLUSIVE, every well-channel below the ESS floor. So the forecast loss is the live evidence
here, and the "no setting of the random-walk scales satisfies both tests" argument is suspended,
not available for reuse.

A decelerating culture *is* the state shift. COSMIC-dFBA's `f` and this project's missing drift term
are the same missing quantity.

**What would beat LDA plus a logistic curve at this sample size.**

| Option | What it is | Fit here |
| --- | --- | --- |
| Continuous latent state inside the existing particle filter | add a bounded state (logit of `f`, or a drift on growth rate) to the state vector in `estimator.py::_propagate`; the mechanistic ODE is the transition model | **Best.** ~20 lines. Non-negativity and boundedness are already handled by the multiplicative-walk pattern. Gives `f` a posterior instead of a point label, and inherits the NIS and whiteness tests for free |
| Interacting multiple model / Rao–Blackwellised particle filter | keep a discrete mode (which cause) and marginalise it analytically; continuous state by SMC. doi 10.1109/9.1299 (Blom & Bar-Shalom 1988); arXiv:1301.3853 (Doucet, de Freitas, Murphy & Russell; the arXiv posting is 2013, the paper is UAI 2000) | **Second.** This is the natural home for L1's cause label: one mode per cause, marginalised out |
| Recurrent switching linear dynamical system | the switch probability depends on the continuous state, so the model says which coordinate drove the switch. arXiv:1610.08466 — note the preprint is titled "Recurrent switching linear dynamical systems" while the AISTATS 2017 version is "Bayesian Learning and Inference in Recurrent Switching Linear Dynamical Systems" | **Read for the idea, do not adopt the model.** State-dependent switching is exactly the L1 mechanism, but the machinery is built for thousands of trials |
| Particle MCMC for the parameters | doi 10.1111/j.1467-9868.2009.00736.x | Later. Useful once the dynamics are right; it will not fix wrong dynamics |
| GP state-space model | arXiv:1406.4905 | **Ignore.** Non-parametric dynamics is the wrong tool when a validated mechanistic ODE is in hand and n is in the hundreds |
| Deep learning | — | **Ignore.** 130 timepoints, or 4342 scored points here, supports a state-space model with a handful of parameters, not a learned transition function |

**No published head-to-head exists.** I searched for a comparison of a switching state-space model
against a discriminative classifier for bioprocess phase inference and found none. That absence is
itself worth stating: the field has not tested the substitution, so the argument here has to rest on
the sample-size and structure argument plus this repo's own diagnostics, not on a benchmark.

**What we should do: adopt.** Give `_propagate` a growth-rate drift, driven by a bounded latent state.
Score it with the NIS and whiteness tests that already exist and the forecast task that already ran.
It is the smallest change in this document and the only one whose failure mode is already measured.

---

## L4. Their method is worse than baseline for about nine amino acids

**Verdict.** Fig 5B is not a bug in the ML layer and it is not the lexicographic ordering starving
low-priority tasks. It is the expected, documented, quantified behaviour of a mammalian GEM asked for
amino-acid exchange. CHO cells consume the branched-chain and aromatic amino acids far in excess of
stoichiometric demand and catabolise the surplus to secreted by-products the model has no reactions
for. A growth-plus-product objective has no term that can generate excess uptake, so it clamps uptake
at the biomass and product requirement. Sharpening the growth and titre predictions sharpens that
clamp — and moves the amino-acid predictions further from a measurement the clamp cannot reach. The
two results in Fig 5B are causally linked, not merely coincident. This is the failure mode this
project inherits directly, because Yeast9 and ecYeastGEM pin amino-acid demand the same way.

**Six of the nine are the documented overflow set.** Gonzalez et al., "Comprehensive stable-isotope
tracing of glucose and amino acids identifies metabolic by-products and their sources in CHO cell
culture", *PNAS* 121:e2403033121, 2024, doi 10.1073/pnas.2403033121: comprehensive ¹³C tracing in two
IgG-producing CHO lines verified **45 secreted by-products**, the majority derived from glucose,
**leucine, isoleucine, valine, tyrosine, tryptophan, methionine and phenylalanine**.

| COSMIC-dFBA Fig 5B, worse than dFBA (**paywalled**) | In the ¹³C by-product source set? |
| --- | --- |
| leucine, isoleucine, valine | yes (BCAA catabolism) |
| phenylalanine, tryptophan | yes (aromatic catabolism) |
| methionine | yes |
| threonine, histidine, lysine | no — but all three are essential in mammalian GEMs |

Six of nine overlap exactly. Mulukutla et al. name the two pathways: Phe–Tyr, giving 3-phenyllactate
and 4-hydroxyphenyllactate, and BCAA, giving isovalerate, isobutyrate and 2-methylbutyrate
(doi 10.1016/j.ymben.2019.03.001, *Metab Eng* 54:54–68, 2019; BCAT1 knockout eliminated all BCAA
by-products).

**How badly a CHO GEM does on this, in numbers.** Széliová et al., doi 10.1371/journal.pcbi.1009022
(*PLoS Comput Biol* 17:e1009022, 2021), pFBA on iCHO1766 against 20 ¹³C-MFA datasets from 6
publications, predicting minimum uptake rates:

| Metabolite | R² | median relative error |
| --- | ---: | ---: |
| glucose | 0.92 | 4.1 % |
| glutamine | 0.75 | 50.4 % |
| asparagine | 0.12 | 24.2 % |
| all remaining amino acids | **≤ 0.06** | **65.8–176 %** |

Verified verbatim from the full text. In the same work, measured uptake was sometimes *below* the
model's stoichiometric minimum (tyrosine and cysteine, three datasets each) and had to be overwritten
— the model could not even accept the data. Essential amino-acid exchanges are conventionally inputs
to a mammalian GEM, not outputs of it.

**Testing the brief's four hypotheses against the evidence.**

| Hypothesis | Verdict |
| --- | --- |
| Transporter-limited, not objective-driven | **Strongly supported, and dominant.** Sreejan et al. built the alternative: a CHO multi-amino-acid multi-transporter kinetic model giving uptake as a function of medium concentrations and transporter levels, with no cellular objective anywhere, validated on held-out uptake data (doi 10.1016/j.jbiotec.2021.12.003). Ley et al. is the cleanest refutation of objective-driven uptake: Hpd and Gad2 disruptions gave **unchanged amino-acid uptake rates** while growth rose up to 19 % and IVCD up to 50 % (doi 10.1016/j.ymben.2019.09.005) |
| Protein turnover dominates net exchange | **Partly, at the transport level.** Nicolae et al. found all exchanged amino acids except asparagine exchange *reversibly* with the medium, reversibility varying in time while net flux stayed constant, and named "simultaneous synthesis and catabolism" (doi 10.1186/1752-0509-8-50). Net exchange genuinely is a small difference of large bidirectional fluxes. No paper quantifies proteolysis' share in CHO — the weakest-documented of the four |
| Lexicographic ordering starves low-priority tasks | **Real mechanism, thinly documented, and not the main cause here.** Lexicographic optimisation is the standard fix for non-unique exchange fluxes in dFBA — DFBAlab exists because "few existing simulators address ... nonunique exchange fluxes" (doi 10.1002/bit.24748). So in any lexicographic dFBA the amino-acid exchange values are set by where they sit in the ordering, and reporting them as predictions reports an implementation choice. Whether COSMIC-dFBA uses lexicographic ordering is in the paywalled Methods — **UNVERIFIED** |
| Degeneracy: poorly constrained by the biomass objective alone | **Strongly supported and quantified.** Mahadevan & Schilling's alternate-optima result (doi 10.1016/j.ymben.2003.09.002) is the general statement; Széliová's R² ≤ 0.06 is the CHO measurement. Chen et al. name it: "overly-restrictive constraints, including essential amino acid exchange fluxes that can lead to improper predictions" (doi 10.1038/s41540-019-0103-6) |

**The direction of the error is known and matches.** Schinn et al., doi 10.1002/bit.27714, state it in
one sentence: "the models assume an optimal and highly efficient metabolism, and therefore tend to
**underestimate amino acid consumption**". Their fix is to bolt ML onto the GEM — the same genre as
COSMIC-dFBA. They needed ML *because* the GEM under-predicts amino-acid consumption. That is the same
defect reappearing as negative bars.

**What we should do: adopt as a guard, not chase a fix.** Three actions, all cheap.

1. **Never report an amino-acid exchange flux as a prediction.** In this project's vocabulary they are
   inputs or unconstrained outputs, not claims. `docs/CLAIM_BOUNDARY.md` is the right place for this
   and the boundary should say so explicitly, the way it already refuses oxygen- and OUR-linked
   claims because the GEM misses oxygen by 65 %.
2. **Expect the same coupling here.** Any improvement to growth or product prediction will tighten the
   stoichiometric clamp on amino-acid uptake. If a future result shows growth improving while
   amino-acid agreement degrades, that is the predicted behaviour and not evidence of a bug.
3. **If amino acids ever have to be predicted, add the sink, not a better objective.** Either a
   transporter-kinetics layer in the Sreejan style, or data-driven bounds in the NEXT-FBA style
   (doi 10.1016/j.ymben.2025.03.010 — ANNs trained on exometabolomics against ¹³C fluxomics emit
   bounds that constrain the GEM). Both fix the constraint. Neither touches the objective.

---

## L5. Eight media conditions, one CHO line, one product

**Verdict.** The field does validate on single-cell-line datasets, and eight conditions is
representative rather than unusually small. No hybrid ML + genome-scale *dynamic bioprocess* model
published to date has demonstrated cross-clone or cross-product transfer. Two things do generalise,
and neither is the thing COSMIC-dFBA is: bacterial hybrids generalise across media and across
hundreds of designs, and one Gaussian-process hybrid genuinely transfers across products — without a
GEM. Say this plainly; it is a fact about the field, not a criticism of one paper.

| Work | Scope | Validation size | Transfer shown |
| --- | --- | --- | --- |
| COSMIC-dFBA, doi 10.1016/j.ymben.2024.02.012 | CHO, dynamic, GEM | 8 media conditions, 1 line, 1 antibody | none |
| doi 10.34133/csbj.0078 (2026) | CHO, ODE + constraint-based + ML | **23 fed-batch cultures** | none claimed |
| Moreno-Paz et al., doi 10.1111/1751-7915.13995 | yeast, ecYeast8 inside dFBA | 9 bioreactor datasets: chemostat, batch, fed-batch, a Δpdc producer | across wild-type-like strains; **no accuracy metric reported**, trajectory agreement only |
| dAMN, doi 10.1093/bioinformatics/btag230 | *E. coli* + *P. putida*, neural + dFBA | 280 and 81 media | generalises to **unseen media**, mean R² ≥ 0.9; separate model per organism, no cross-organism transfer |
| AMN, doi 10.1038/s41467-023-40380-0 | *E. coli* | 73 media | cross-validation R² = 0.78 |
| Oyetunde et al., doi 10.1371/journal.pone.0210558 | *E. coli* cell factories, GEM features + ML | ~1200 designs from ~100 papers | **genuinely cross-strain and cross-product**, Pearson r 0.8–0.93 on unseen data — but a static titre/yield/rate regressor, not dynamics |
| Hutter et al., doi 10.1002/bit.27907 | hybrid GP with entity-embedding vectors | multiple products | **genuinely cross-product**, with interpretable product similarity — **no GEM** |

**And the GEM layer itself may not carry clone-specific information.** Strain et al.,
doi 10.1002/bit.28366 (*Biotechnol Bioeng* 120:2460–2478, 2023): CHO cell-line-specific GEMs "better
capture extracellular phenotypes but failed to improve intracellular reaction rate predictions". That
is part of why cross-clone GEM transfer is unproven — there may be less clone information in the
network than the effort assumes.

**What we should do: ignore as a method, adopt as claim discipline.** There is no published technique
to import. What is importable is the standard: `docs/research/GENERALIZATION.md` already argues this
project should stop being one team's project, and this section supplies the external benchmark for
what "demonstrated" means in the field. The honest bar is low, which means a modest cross-condition
demonstration here would be competitive — and also means nobody gets to claim generalisation from
one strain.

---

## L6. Proteome and enzyme allocation

**Verdict.** Yes, largely — and this is the finding that should change the architecture. For the class
of shift that motivates state-specific objectives in the first place, a proteome-constrained model
gets there with **one fixed objective**: maximise growth, subject to compartment-specific protein
pools. The state shift is not a change of objective; it is a change in *which constraint is binding*.
That is a strictly better mechanism, because it is mechanistic, it is falsifiable, and it produces the
cause of the shift as a by-product — which is L1's answer as well. Enzyme constraints do not subsume
*everything*: they miss signalling-driven regulation, they miss stationary and decay phase, and they
trade a hand-specified objective for a hand-calibrated protein pool. But they subsume the part
COSMIC-dFBA built machinery for.

**The single-objective result.** Elsemman et al., "Whole-cell modeling in yeast predicts
compartment-specific proteome constraints that drive metabolic strategies", *Nature Communications*
13:801, 2022, doi 10.1038/s41467-022-28467-6. Verified verbatim: the model "predicts metabolic fluxes
and corresponding protein expression by constraining compartment-specific protein pools and
**maximising growth rate**". Under glucose limitation a mitochondrial constraint limits growth at the
onset of ethanol formation — the Crabtree effect. Under sugar excess, a constraint on total cytosolic
volume dictates overflow metabolism. No objective switching anywhere; the *binding constraint* switches.

Supporting the same conclusion from three directions:

- Sánchez et al., doi 10.15252/msb.20167411 (GECKO, *Mol Syst Biol* 13, 2017), state the problem
  plainly: traditional GEMs cannot show overflow metabolism unless ad hoc constraints or objective
  functions are imposed. ecYeast7 reproduces the Crabtree transition, carbon-source specificity and
  knockout shifts that Yeast7 misses.
- Salvy & Hatzimanikatis, doi 10.1073/pnas.2013836118 (*PNAS* 118:e2013836118, 2021): **diauxie itself
  emerges** as an optimal growth strategy under resource-allocation constraints. No regulatory rules,
  no phase objectives.
- Molenaar et al., doi 10.1038/msb.2009.82 (2009): the coarse-grained argument that the growth-strategy
  shift is a tradeoff in cellular economics, from growth-rate optimisation alone.

**The comparison you asked for does not exist — and that is the novelty claim.** No paper compares an
enzyme-constrained model against an objective-priority or state-switching model on the same dataset.
The nearest thing is a gap with a specific shape: MPMO and IMC (§L2) are the yeast twin of
COSMIC-dFBA's phase-objective machinery, run on *Saccharomyces* batch fermentations, and **neither
uses an ecModel anywhere**. Meanwhile Moreno-Paz et al. (doi 10.1111/1751-7915.13995) already ran
ecYeast8 inside dFBA across chemostat, batch and fed-batch and got the critical dilution rate
(0.27 h⁻¹) and ethanol onset that Yeast8 misses — but reported no error metric. **Running
ecYeastGEM against MPMO/IMC-style phase objectives on the same fermentation has not been done.**

**What this repo has already measured, which is the strongest single argument.**

| Repo result | Source | What it means |
| --- | --- | --- |
| Plain Yeast9 v9.0.2 given measured uptake rates overpredicts growth by **73 %** (0.69 vs 0.40 /h) | `docs/FINDINGS.md` | the unconstrained model is not usable for physiology |
| GECKO ec model reproduces growth, glucose and ethanol to **6–16 % from its protein pool alone** — but misses oxygen by **65 %** and CO₂ by **62 %** | `docs/FINDINGS.md` | the protein budget carries most of the phenotype; the gas exchange is still wrong |
| van Hoek 1998 critical growth rate ≈ 0.28 /h (PMID 9797269); an oxygen ceiling alone puts it near 0.45 /h. Recorded as a held-out **failure** | `docs/CLAIM_BOUNDARY.md` | "the real limit is proteome allocation, which ecYeastGEM models and this does not" — already written down, and pcYeast is the paper that confirms it |
| D1 measured FVA relative width **1.000** for product flux: the feasible interval runs from zero to the ceiling at every growth level | `src/ystwin/fba/solver.py` | with the current constraints, any single product flux is a property of the solver. Enzyme constraints are the mechanism that narrows this |

The repository has independently reproduced the L6 argument as a held-out failure and written the
diagnosis in the right words. It has not acted on it.

**Where enzyme constraints are not enough — state these, or the recommendation is dishonest.**

- **Ribosome and TOR signalling.** pcYeast fails here, and says so: under translation inhibition
  "the ribosomal proteome fraction increased much less with inhibition than the model predicted",
  concluding that "the expression of ribosomes in yeast is dominantly regulated by environmental
  nutrient signalling and less by internal cues" (verified verbatim, doi 10.1038/s41467-022-28467-6).
- **Glucose repression and isoenzyme usage** needed a Boolean signalling network bolted onto an
  ecModel (Österberg et al., doi 10.1371/journal.pcbi.1008891).
- **The transient of a shift** is not set by the pool alone: decFBA gets the diauxic shift but not the
  lag; an explicit rate-of-enzyme-change constraint was needed (Karlsen et al.,
  doi 10.1371/journal.pone.0280077).
- **Nitrogen limitation breaks the tight-budget assumption.** Yu et al., doi 10.1038/s41467-020-15749-0:
  seven metabolic superpathways hold over 50 % capacity in reserve, glucose metabolism over 80 %, and
  over 50 % translational reserve for 2490 of 3361 expressed genes. If this project's fermentation runs
  into N limitation, a saturated protein pool is the wrong assumption.
- **The pool is fitted.** GECKO assumes an average enzyme saturation and calibrates it; GECKO 2.0
  reports that a large share of kcat values come from other organisms or non-specific mechanisms
  (PMID 35773252; GECKO 3.0 protocol, PMID 38238583, makes ML-predicted kcats the default). Resource
  balance analysis of yeast needed condition-specific kapp and maintenance ATP, with "estimated
  parameter values ... found to vary with oxygen and nutrient availability"
  (doi 10.1016/j.ymben.2023.04.009). You trade hand-specified objective priorities for one or two
  hand-calibrated allocation parameters. That is a better trade, not a free one.
- **Stationary and decay phase are essentially untouched** by ecModels; the yeast fermentation work
  that models those phases uses non-enzyme-constrained models with phase objectives.

**One directly relevant precedent for the product side.** pcSecYeast applies proteome constraints to
recombinant protein secretion (doi 10.1038/s41467-022-30689-7), which is the closest published
analogue to a yeast product-burden project.

**What we should do: adopt.** Make ecYeastGEM_batch the model of record for any physiology or product
claim, and stop treating the growth/production state as an objective choice. Two specifics.

1. **The protein pool is the state variable, not `f`.** A production state is a reallocation of the
   same budget, so it belongs in the constraint. This also removes the need to specify two models and
   blend them.
2. **The plumbing is already correct and already tested.** `fba/physiology.py::cap_uptake` reads the
   irreversible split from the model rather than the filename, after measuring that setting
   `lower_bound` on `r_1714` capped nothing while glucose flowed at 17.93 through `r_1714_REV`
   against a requested cap of 21.3. That was the hard part of running an ecModel in a dFBA loop and
   it is done.

---

## Ranking: how much the fix buys, per unit of effort

| Rank | Limit | Effort | What adopting buys | Why here |
| --- | --- | --- | --- | --- |
| 1 | **L3** state inference | ~20 lines in `estimator.py::_propagate` | turns the project's largest measured failure into a scored improvement | The defect is diagnosed by the scored forecast loss of 23 %, and the scoring machinery already runs. The innovation-calibration route is suspended (`nis.current.*` refuse; `tested_well_channels: 0`), so this rests on one measurement, not three. Nothing else in this document can be scored the day it is written |
| 2 | **L6** enzyme constraints | days: wire the ecModel into the prediction path | replaces the entire state-objective layer with a mechanism, and narrows an FVA width of 1.000 | The model file, the manifest hash and the irreversible-split handling all exist. The repo has already reproduced the argument as a held-out failure |
| 3 | **L4** amino acids | hours: a paragraph in `CLAIM_BOUNDARY.md` | prevents an unfalsifiable claim and pre-explains a result that will otherwise look like a bug | Cheapest item here. Buys honesty, not accuracy — but it is the failure mode this project inherits and it is fully explained by literature |
| 4 | **L1** cause resolution | small for the dual, large for the classifier | the novelty claim, and a cause label for free | Shadow prices cost nothing and are the authoritative signal. The classifier half is gated on G4, which is INCONCLUSIVE, so it cannot be claimed yet |
| 5 | **L5** validation scope | none available | nothing to import | The field's bar is low. Useful as a benchmark for claims, not as a method |
| 6 | **L2** learned objectives | weeks, for a negative result | nothing | CellTarget already ran this experiment with better data than this project has and found the objective non-identifiable. IMC and MPMO tie at P ≈ 0.7 |

## The one change

**Give `estimator.py::_propagate` a growth-rate drift, driven by a bounded latent state.** Three
reasons it goes first. It is the project's own most important open defect, named in its own words —
the filter is over-confident, its innovations are not white, no setting of the random-walk scales
fixes both, and the diagnosis is that a decelerating culture is extrapolated flat. It is L3's fix,
and it is L1's substrate: once the phase fraction is a filtered state with a posterior rather than a
per-timepoint label, a discrete cause mode can be attached to it and marginalised out. And it is the
only change here that arrives already falsifiable, because the NIS test, the whiteness test and the
84-well forecast task all exist and all currently fail in the same direction, so the change either
moves three numbers or it does not.

L6 is the larger architectural win and should be the second change, not the first. There is no point
blending two genome-scale states by a fraction `f` while the filter that would supply `f` cannot
track a culture slowing down.
