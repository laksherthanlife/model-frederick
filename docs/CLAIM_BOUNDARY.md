# The claim boundary

Every number in the README is one of four kinds of thing, and the difference between them
is larger than the difference between any two numbers within one kind. This file assigns a
tier to each claim so that a reader who spot-checks one knows what they have checked.

The tiers are not invented here. They are the vocabulary of
[`research/CIRCULARITY.md`](research/CIRCULARITY.md) §5, each with a citable precedent.

| Tier | Name | What it means | Precedent |
| --- | --- | --- | --- |
| **0** | **Asserted** | A constant chosen for plausibility. Nothing measures it; nothing yet contradicts it. | — |
| **1** | **Self-consistent in simulation** | The fitting procedure correctly inverts the generator it was given, or a deterministic computation was carried out correctly on a model we adopted. Necessary, and much weaker than it looks. | Simulation-Based Calibration — Talts et al., arXiv:1804.06788 |
| **2** | **Robust across generator families** | Holds when the generator's *structure*, not just its parameters, is perturbed — or when it holds on a measurement the model never saw. Reported with a UCBOG where it comes from simulation. | Domain randomization (Peng et al., ICRA 2018); SOB/UCBOG (Muratore et al., TPAMI 2021) |
| **3** | **Predicted then measured** | A prediction registered *before* the measurement and scored against it. | Sim-vs-Real Correlation Coefficient — Kadian et al., *IEEE RA-L* 2020 |

## Scope: intracellular products, and a terminal node that may also leave

Stated before the tiers because it decides which of them can apply at all.

The product layer solves ``d[X]/dt = v_in - v_out - mu*[X] = 0`` at every node. The
``mu*[X]`` term is growth dilution. For a product that accumulates inside the cell —
beta-carotene, lycopene, PHB, squalene — it is the only thing besides the next enzyme that
removes a pool, and the balance closes on it alone.

A **terminal** node may declare two further outlets, and each is refused at SOLVE time
until its constants are supplied, by name:

| fate | outlet | constants required |
|---|---|---|
| `diluted` | — | none |
| `degraded` | `k_deg*[X]` | `degradation_rate_per_h` |
| `secreted` | `vmax*[X]/(km + [X])` | `secretion_vmax_mmol_per_gdcw_h`, `secretion_km_mmol_per_gdcw`, `growth_rate_range` |

An **intermediate** declaring either is still refused at LOAD, and no measurement lifts it:
the walk carries one flux forward from each node to the next, so an intermediate leaking
sideways breaks the walk rather than a coefficient in it.

**Where the line moved, and why.** Both non-dilution fates were once refused at load, which
made every storage compound and every secreted product *unrepresentable* rather than
*under-parameterised*. That confused two different things. A pathway whose chemistry this
walk cannot describe should be refused at load; a pathway one measured number short should
be refused where the number is needed, so the message can name the missing measurement and
the experiment that would produce it. `DEGRADED` moved on 2026-09-01 and `SECRETED`
followed. A model that answers a question it cannot answer is worse than one that declines
— and a model that declines when it is only missing a constant tells the reader less than
it could.

**Nothing was unlocked by this.** No shipped pathway declares a secreted node, and none
honestly can yet: the dataset that motivated the outlet (Kastberg 2025) reports no absolute
quantity of product anywhere, only relative label-free fold changes. Farnesene, resveratrol
and TAL remain unanswerable, and now say so with the name of the missing measurement.

**Why the export outlet saturates and the degradation outlet does not.** Not symmetry, and
not preference. Kastberg 2025 measured the intracellular proteome and the secretome on the
same K. phaffii chemostat samples: the human insulin precursor rises 4.12 log2FC (17.4×)
inside the cell while the secretome does not move, and Mambalgin-1 rises 4.15 log2FC inside
against 1.8 log2FC secreted. A first-order outlet secretes in proportion to the pool by
construction, so no rate constant reproduces a 17× pool beside a flat supernatant. For
GPH1/SGA1 turnover no such curvature has been observed, so approximating it with one
constant is still the better trade.

## The product prediction was circular until August 2026, and is not now

`Environment.pathway_flux` was a caller-supplied constant, and both entry points obtained it
as ``q_lycopene + q_betacarotene`` from the very chemostat state being predicted. The model
measured the product in order to predict the product. `scripts/predict_product.py` printed
"this is a residual and not a prediction" in its own output, which was honest, and no
non-circular source of that number existed anywhere in the repository.

**Capacity-only FVA is not an allocation law.** In the previously tested formulation,
product upper bounds left zero product flux feasible. This is a property of those
constraints, not a theorem that all FBA task policies or state mixtures cannot supply a
heterologous flux. Allocation assumptions and matched host/uptake conditions must be
stated and tested separately.

**Expression-to-flux evidence.** `scripts/fit_pathway_flux.py` compares measured genes
using leave-one-strain-out predictions and a geometric-mean baseline fitted only on the
training strains. Earlier skill values used the full-data mean and are withdrawn.
The corrected descriptive ranking still favors the heterologous genes, but ERG9 has
slightly positive skill. “Every native gene is negative” and the inference that the
ranking excludes normalization or growth confounding are not supported.

**Whole-chain validation.** There are three independent held-out strains, each with two
condition predictions. Gene selection, model selection, entry-flux fitting and branch
kinetics must occur within training folds. β-carotene and lycopene are reported jointly,
with constant-rate, constant-content and entry-only comparisons. The earlier fixed-CrtE
scores (14.2% median error on β-carotene and 26.3% on lycopene) are historical scores of a
retrospectively selected model, not nested-selection performance. Current claims require
the regenerated fold-level record, including failed predictions and measured interval
coverage; they do not establish prospective or general biological validation.

**Nothing in this repository is Tier 3.** Forward predictions have now been scored on
held-out real data, and the result splits on dose: the fitted dose-response beats
carrying the nearest measured dose forward at every interior rung at or below 0.5 mM
(skill +0.09 to +0.37) and loses above it (−0.18, −0.68), with the sign separating on
`dose / EC50`. That is recorded in
[`WHY_INTERPOLATION_WORKS.md`](WHY_INTERPOLATION_WORKS.md) and
[`research/PREDICTION.md`](research/PREDICTION.md). It is a result, and the winning half
has a mechanism behind it that predicts the direction before the data is seen — but every
one of those scores was computed after the measurement rather than registered before it,
so none of it clears the bar.

**One thing is now positioned to.** `outputs/registered_prediction_D018.csv` predicts
beta-carotene for three strains at D = 0.18 /h, with intervals, written down before any
measurement exists at that dilution rate. It is Tier 0 today and stays Tier 0 if nobody
runs the chemostat. It becomes the first Tier 3 entry in this file if somebody does — and
it can be scored rather than rationalised either way, which is the property the tier is
actually about.

## Three things this tiering is meant to stop

**Tier 1 numbers reading as Tier 2 numbers.** Same-teacher recovery does not validate shared
biological assumptions. Simulation Optimization Bias is non-negative **in expectation** for
an optimized finite-sample objective under the assumptions in `research/CIRCULARITY.md` §1;
it is not a pointwise guarantee that every held-out simulator metric overstates real-plate
performance. Independent held-out evaluation can be valid within the simulator distribution
while the biological reality gap remains unmeasured.

**Tier 0 constants reading as measurements.** A number in the source with three significant
figures looks like it came from an instrument. Several did not.

**A result measured on real data reading as a prediction.** G1, D2, G4, the dose-response
correction and the calibration work are all measured on real plates, and none of them
forecasts anything. They are *diagnostic*: they say what a measurement is worth, which is a
different and more modest job than saying what the next one will be. On the DNV capability
scale that is level 2 reaching for level 3 — see [`CONTRACT.md`](CONTRACT.md).

---

## Tier 0 — Asserted

A constant chosen for plausibility. Listed so that no reader mistakes one for a
measurement, and so that measuring one is visibly worth doing.

| Constant | Where | Status |
| --- | --- | --- |
| latent-state → ATP maintenance scale | `bridge/latent_bridge.py::_MAX_STRESS_MAINTENANCE` = 6.5 mmol ATP/gDW/h at full scale | **An envelope, and the value it replaced was wrong by 10⁴.** `_MAINTENANCE_PER_ACTIVITY = 300.0` multiplied a promoter activity in RFU/OD/h and called the product mmol ATP/gDW/h; on the measured constructs that demanded an NGAM of 316,000–391,000 against a resting 0.7, where Yeast9 goes infeasible at 19.19. The replacement is bounded by reproducing Lahtvee 2016's (PMID 27307591) glucose-limited chemostat design *in the model* — not read off their figure, which is a bar chart — so it is still **Tier 0**: an envelope consistent with a measured design, not a measured value. It also inherits the observed activity range, which is why `activity_full_scale` is an argument rather than a hidden constant. Still refused by G4 at every call site. [`superseded/maintenance-coefficient.md`](superseded/maintenance-coefficient.md). |
| reporter autofluorescence | `observation.py` | **Never measured on any plate.** Two reporter-free BY4741 controls were plated and both were read in OD600 only. See [`DATA_INVENTORY.md`](DATA_INVENTORY.md). |
| `od_linear_max` | `gates/g1_optical.py` | Placeholder. `calib/od_linearity.py` turns a dilution series into it; the series has not been run (protocol P4). **A linear range has since been measured for the FLUORESCENCE detector** (`calib/gain_linearity.py`, from the second PMT gain in every `+75100` protocol) — it is a photomultiplier and this is a photodiode, so it is a measurement of the other channel and this constant stays a placeholder. `tests/test_gain_linearity.py` holds it at 1.0 so the two cannot be confused. |
| H2O2 EC50, 0.5 mM | `generator/stress_panel.py` | Literature. The shape-free crossing agrees (0.498, 0.396) but the biphasic fit is unidentifiable because every dose that would show the falling limb killed the culture — so this is a literature value that survived a check, not a measured one. |
| reporter loadings and crosstalk topology | `generator/stress_panel.py::reporter_loadings` | Transcribed from the literature with PMIDs, and corrected against it once (see the README's wiring table). Transcribed is not measured. **Everything in the known-loadings path below rests on these.** |
| the generator's single configuration `ξ̄` | `generator/stress_panel.py` | One crosstalk topology, one set of EC50s, one loading matrix. `research/CIRCULARITY.md` §7 Fix 1 is the change that would move the whole Tier 1 block toward Tier 2. |

Two Tier 0 constants have been **independently corroborated**, which is worth distinguishing
from measured. Six generator constants were chosen from the literature while building the
panel and only afterwards compared against BioNumbers — the Milo lab's curated database,
assembled by people with no interest in this model. Five agreed; the sixth disagreement is a
real split in the published record and is recorded as such. That raises confidence in the
constants. It does not make them measurements of *these* cultures.

## Tier 1 — Self-consistent in simulation

Two distinct things live here and it is worth keeping them apart, because only the first
carries the Simulation Optimization Bias.

**(a) The fitter inverts its own generator.** Data was produced by `stress_panel.py` /
`panel_experiment.py` / `culture.py`, a method was fitted to it, and the method recovered
what was put in. Optimising and scoring on the same generator draw is exactly the situation
the SOB theorem describes, so every number here is an overestimate of the plate value.

| Claim | Where in the README |
| --- | --- |
| Detection vs attribution power (0% / 4% / 46% / 98% and the nutrient-axis rows) | *What the design has to be* |
| Collinearity 0.998 dose-only vs 0.570 crossed | *Simulation-first* |
| Latent states found: one-at-a-time 2 states / 70% spurious vs crossed 2 / 8% | *What the design has to be* |
| Transfer R²: 0.064 blocked, 0.353 all pairs, 0.539 six selected | *Transfer* |
| Combination lift 0.038 → 0.429 (11×), all-pairs 0.161 | *Where the numbers stand after both audits* |
| Pair design across 25 stressors: 0.029 / 0.236 / **0.295** / 0.281 / 0.121 | *The pair design* |
| Latent width vs redox recovery and held-out module score | *Latent width* |
| Wide vs buildable panel, modules recovered and held-out score | *What it reports* |
| Deployable model applied to an unseen plate: r = 0.999 redox … 0.812 oxidative | *What it reports* |
| ESR across contexts (0.06 / 0.35 / 0.63 narrow vs 0.21 / 0.11 / 0.32 broad) | *Training across culture conditions* |
| Estimator recovers a known input to 1.1% median relative error | *Measured estimator behaviour* |
| Module recovery 0.93–0.99 against a shuffled control near 0.02 | *Transfer* |

That last row is the one Tier 1 claim in the block that has been attacked by a null.
`analysis/recovery.py::module_recovery(shuffle=True)` is a real null and it was already
there. It is a null on *module recovery*, not on transfer or power — see
[`NULL_RESULTS.md`](NULL_RESULTS.md).

**(b) A deterministic computation on a model we adopted.** No fitting, no scoring on the
same draw, so no SOB. These are facts about Yeast9 / ecYeastGEM / the ModelSEED and
eQuilibrator tables, computed correctly. They are only as good as those models.

| Claim | Note |
| --- | --- |
| D1: feasible product-flux range `[0, ceiling]`, relative width 1.000 at every growth level | Structural, and the most robust Tier 1 result here — a capacity bound moves the top of the interval and cannot lift the floor off zero. |
| D3: the regulation layer does **not** narrow the product range, and provably cannot | **Tier 1 negative, and sharper than D1's.** Per-reaction E-Flux regulation replaces D1's single ATP-maintenance scalar; relative width is 1.000 before and 1.000 after, at all seven growth fractions and for all four stressor probes — **0 of 28 rows narrowed, every regulated floor exactly zero**. The reason is structural, not empirical: E-Flux (Colijn 2009, PMID 19714220) constrains `|v| ≤ c·e`, a bound on how much flux is *permitted*, and no upper bound can make a flux mandatory. A repression sweep down to 1% of reference capacity leaves the floor at zero with growth already at 1.8% of maximum. This retires the objection that D1's width was an artefact of the ATP-scalar coupling. |
| ~~The ceiling corresponds to ~33 mg/gDCW at 90% of max growth~~ **REFUTED as a bound** | The arithmetic is right (32.84 mg/gDCW) and the bound is not: three published *S. cerevisiae* contents exceed it, up to 2.41x, two of them beta-carotene by HPLC from two independent groups. (Read "five ... across four groups" until 2026-08-30, a count that still included the two Bubphasawan rows struck when PMID 38710418 was found withdrawn.) The cause is `growth_fraction=0.90`, which is nowhere measured — the frontier is exactly linear, so the content ceiling is `111.38/mu - 295.58` mg/gDCW and **diverges as mu -> 0**. FBA places no finite bound on specific content in a non-growing cell. See [`EXTERNAL_CAROTENOID_BOUND.md`](EXTERNAL_CAROTENOID_BOUND.md). |
| TMFA: growth 0.08584 preserved; no flux range narrowed on 86 reactions | An LP fact, and the swing table that explains it is arithmetic. |
| The Mg²⁺ and membrane-potential corrections; the eQuilibrator vs group-contribution table | The inputs are Tier 0 tabulated energies. |
| Hours to 90% of steady state (5.8 / 7.7 / 11.5 / 23.0 / 115) | Exact consequence of `dR/dt = k − (μ + k_deg)R`. Arithmetic, given the model. |
| Surrogate 161× per call, 148,000× batched, median error 0.000% | A software benchmark on this code, not a claim about yeast. |

### The identifiability tables are Tier 1 *and* they are about the wrong path

This is the correction that most changes how the README reads, and it comes from
[`research/IDENTIFIABILITY.md`](research/IDENTIFIABILITY.md) §1–2.

Two modules in this package answer to the word "identifiable", and they are not the same
problem:

- **Path A — loadings known.** `analysis/design.py` and `analysis/sensor_selection.py` build
  `L' S⁻¹ L` from `reporter_loadings`, which is *given*. Identifiability is the rank of `L`,
  so *k* channels give ≤ *k* modules. **The README's `2 → 2 of 7`, `4 → 4 of 7`, `7 → 7 of 7`
  table is correct within this framing**, and so is the later, stronger statement that
  N channels resolve N of the 24 state variables.
- **Path B — loadings fitted.** `analysis/latent.py::_factorise` runs an SVD and *estimates*
  the loadings. This is exploratory factor analysis, and it is bound by Ledermann (1937):
  `m ≤ (2p + 1 − √(8p+1))/2`. At *p* = 4 that is **one** factor. Not four.

`run_training.py` and `run_transfer.py` — every transfer and every "general stress state"
result — go through **Path B**. The identifiability tables come from **Path A**. *No line of
code connects them, and the Path A tables are not a warrant for the Path B results.*

Two consequences the README must carry rather than imply:

1. The recommended four-stress-channel build supports **one** fitted factor on Path B. The
   README's latent-width table selects width 4 for that build and remarks approvingly that
   "nothing is compressed away when width equals channels". That configuration is the
   **full-rank degeneracy**: a rank-*p* model of a *p*-column matrix reproduces any masked
   entry, so the held-out error is not held out. `select_dimension` lands on it 1–2 times in
   30 at *p* = 4–5, and the Ledermann ceiling is not enforced anywhere in the code.
2. The recoverable ground is to stop doing factor analysis at all on Path B, since `L` is
   already known: solve `y = Lx + ε` under non-negativity and sparsity. Its specific,
   serious catch is that sparsity fails for ESR — a module that is on almost always is not
   sparse — which is why the dedicated `STRE-general` channel is the only route to it.

## Tier 2 — Robust across generator families, or held out from the model

Cross-family comparisons and changed-equation teacher checks now exist; their scope and
negative results are recorded in `CROSS_FAMILY.md` and architecture section 10. They must
not be pooled with same-teacher recovery or presented as biological validation. The
assumptions behind simulation optimism are separated in `research/CIRCULARITY.md`.

The other kind of evidence is agreement with a measurement the model never saw. Two
historical physiology examples follow, including a recorded failure.

| Claim | Status |
| --- | --- |
| Vos 2016 retentostat: measured 0.039 ± 0.003 mmol glucose/gDW/h, model predicts **0.036** | Genuinely held out — the maintenance parameter was set before that experiment was published, and the chemostat series yeast-GEM *was* fitted to is excluded deliberately. **Not Tier 3**: the measurement predates this project and nothing was registered in advance, and the parameter being vindicated is yeast-GEM's, not this repository's. |
| van Hoek 1998: critical growth rate near 0.28 /h; an oxygen ceiling alone puts it near 0.45 | A held-out **failure**, recorded as a test that asserts the discrepancy. The real limit is proteome allocation, which ecYeastGEM models and this does not. |

The forward-prediction score in `research/PREDICTION.md` §0 belongs in this row rather than
in Tier 3, and it is the most important entry in this file. Conditioning `estimator.py` on
the first 2 h of 84 real wells and forecasting the remaining 12 timepoints, **the twin loses
to a six-point log-linear extrapolation** on OD by 23% of level at 2 h, ties on mCitrine, and
covers 40–49% of the time with a nominal 95% interval. The loss is diagnostic rather than
diffuse — `|bias|/sd > 1`, and the cause is named: `_propagate` gives growth rate no drift,
so a decelerating culture is extrapolated flat. A χ² test on the filter's own innovations
finds the same defect by an independent route. It was scored after the fact, so it is not
Tier 3; it is a measured, held-out, honestly negative result, which is the best thing in the
repository.

### Published HOG: a narrower retrospective parameter-fit holdout

The source-backed HOG workflow fits measured phospho-Hog1 observations, not generated latent
labels. Its preliminary Western holdouts are explicitly excluded from the source publication's
parameter fit. They were nevertheless available during original model development, so this is
narrower evidence than an experiment wholly unseen by the model and is not prospective Tier 3.

The estimated activation balance is conditional on the remaining published parameters and an
imposed culture-density trajectory. The absolute Western scale includes literature assumptions;
separate fast rates and the full reaction structure are not learned. The refit does not improve
the stronger-dose waveform over the published parameters. Increasing modeled Gpd1 availability
also does not force the GEM to produce more glycerol in the reported capacity scenario.
`data/hog2013/fit_protocol.json` and each run's machine-generated scores state the exact scope.
Neither all-stress dynamics nor realized product titres are validated by this benchmark.

### The phenotype every physiology claim was scored against was not in the paper it cited

This is the correction that most changes the Tier 2 row above, and it changes the *sign* of
two published findings. `REFERENCE_AEROBIC_BATCH` carried glucose 21.3, oxygen 7.8, ethanol
27.4 and CO2 20.4 mmol/gDW/h, attributed to van Hoek, van Dijken & Pronk 1998. Checked
against the full text of **PMID 9797269**, the strings `21.3`, `27.4` and `20.4` occur
**zero times**, the paper runs no batch culture at all — every steady state in it is a
glucose-limited chemostat — and its strain is DS28911, not the CEN.PK113-7D the tests
claimed. The reference now carries that paper's own Table 1 at D = 0.40 /h: glucose 11.1,
oxygen 3.7, ethanol 13.9, CO2 18.9.

Two findings reverse rather than move:

| Was | Is |
| --- | --- |
| "Plain Yeast9 overpredicts growth by **73%**" | It **under**predicts by **13.5%**. Against the corrected reference both models pass at 35% tolerance. |
| "The ec model misses oxygen by **65%**" | **27%**. Two artefacts inflated it: the reference oxygen was 7.8 against a measured 3.7, and `fba/physiology.py::cap_uptake` was a no-op on the ec model because that model is irreversibly split, so `r_1992` was never actually constrained. |

The old numbers were not preserved as an alternative reading: a number nobody can source is
not a measurement, and there is no honest way to keep it while saying so.
[`superseded/reference-phenotype.md`](superseded/reference-phenotype.md).

**A second mis-citation on the same paper, and it is the more instructive one.**
`bridge/physiology_bridge.py` carried `_GLUCOSE_QMAX = 22.0` as an aerobic glucose uptake
ceiling from van Hoek 1998. That number *is* in the paper, as a different quantity: "the
fermentative capacity showed only a small further increase, up to **22.0 mmol of ethanol** ·
g of dry yeast biomass⁻¹ · h⁻¹", measured in an offline **anaerobic** assay under CO2. Wrong
metabolite, wrong gas regime. The paper's highest in-situ glucose uptake is 11.1.

**What this costs the tier table.** A citation guard that checks "does this identifier
resolve" passes both errors, and so does one that checks "does this figure appear in that
paper". Only reading the sentence around the number catches the second. So the Tier 2 rows
above are only as good as the hand-checks in [`CITATION_AUDIT.md`](CITATION_AUDIT.md), and
`tests/test_citations.py` says so in its own docstring — *"whether a paper supports the
sentence next to it is a judgement, and no test makes it"*.

## Measured on real plates — diagnostic, not predictive

These are not simulation. They are also not predictions. They say what a measurement is
worth. Everything here rests on the four NewProtocol exports, whose actual coverage is in
[`DATA_INVENTORY.md`](DATA_INVENTORY.md) and is smaller than the README used to say.

| Claim | Status and correction |
| --- | --- |
| Plate map recovered, 1.8% and 3.4% median error, both plates independently, logbook-confirmed | The strongest measured result in the file. Note it is a *recovery*, so the map is an output, not an input. |
| G1 pass rates 80% / 100% / 86%, zero linear-range failures on the NewProtocol plates | Measured. |
| D2: zero negative activity, implied loss rate 0.0000 /h on the good plates | Measured, and it **withdrew** an earlier alarm. |
| D2: median R² of `log(RFU/OD)` on growth = **0.82** over 8 G1-passing cultures | Measured — and now attacked. A phase-randomised null on a *superset* of 146 ungated wells gives observed 0.522 against null 0.117, p median 0.022. The coupling survives. **~0.12 of any reported variance share is autocorrelation artefact**, and **43% of wells do not individually clear their own null** — which matters because the correction is applied well by well. The 0.522 is not the 0.82: different well set, no G1 gate. Not a contradiction. See [`NULL_RESULTS.md`](NULL_RESULTS.md). |
| Noise: activity CV 0.146, growth CV 0.33, σ_μ = 0.0117 /h across 210 wells | Measured, and the two independently obtained terms combine to within 3% of the measured total. |
| Lethality: DTT halves growth at 1.45 and 1.65 mM; H2O2 at 0.79 and 1.24 mM | Measured directly from growth, which is why it survived when the reporter-inferred EC50 did not. |
| DTT EC50 "REFUTED at 0.304 mM" | **Withdrawn by the first audit.** Biased low by up to 87%, misspecified fit, peak measured where growth was already at 71% of control. Neither EC50 was ever identifiable. |
| Dose-response, naive vs dilution-corrected (eight rows) | Measured, and **all eight are INCONCLUSIVE at n = 2** — see below. |
| G4 verdicts, all four INCONCLUSIVE | For the ER pair this understates the problem — see below. |

### The dose-response table rests on n = 2, not n = 4

From `DATA_INVENTORY.md`, read off the exports rather than the logbook. The README's
"4 biological replicates" is true of the **OD600** channel and false of the **fluorescence**
channel, which is the one every reported result rests on. Replicates 2 and 4 read
fluorescence on columns 1–3 rows A–G only, and columns 1–3 are UPRE1. So:

| Construct | Biological replicates *with fluorescence* |
| --- | ---: |
| UPRE1 | 4 (two unusable as they stand) |
| UPRE2 | **2** |
| NativeYap1 | **2** |
| AlteredYap1 | **2** |

`analysis/uncertainty.py::fold_change` refuses an interval below three replicates. A third
plate has since been recovered — it was being dropped by a reader that stopped after 21
wells of 87 on a column-stacked export, not by anything wrong with the data — so the table
now has intervals where it previously had eight refusals. **Six calls resolve**: UPRE1 and
UPRE2 each clear "no change" at 0.2, 0.5 and 1.0 mM DTT. The rest stay INCONCLUSIVE, which
is unresolved and not refuted.

**This section previously said all eight rows were INCONCLUSIVE and was left saying it after
the third plate arrived.** That is recorded here rather than quietly corrected, because a
claim boundary that drifts in the *conservative* direction is still a boundary nobody can
rely on, and this document is the one place in the repository whose whole value is being
current.

One point estimate also moves. Computing the fold **within each plate and then combining**,
rather than as a ratio of pooled means, cancels a plate-level shift that multiplies dosed
and control wells alike. `AlteredYap1 @ 1.0 mM H2O2` corrects from **1.50 to 1.44** on that
change alone, and `NativeYap1 @ 1.0 mM` from 1.34 to 1.27. Every other row is unchanged to
two decimals.

### The autofluorescence sweep moves the high-dose rows and resolves none of them

Because autofluorescence has never been measured (Tier 0, above), the per-cell signal
`(RFU − media_blank)/(OD − od_blank)` still contains it: the media blank removes the
`background` term, and autofluorescence is a *per-biomass* term that survives that
subtraction entirely. In the ODE inversion an unremoved constant `a` contributes exactly
`a·μ` to the recovered activity, so it can be swept exactly. Sweeping `a` from 0% to 30% of
each construct's own 0 mM signal, over `outputs/autofluorescence_sensitivity.csv`:

- **The low-dose point estimates barely move.** `NativeYap1 @ 0.1 mM` runs 0.92 → 0.91 and
  `AlteredYap1 @ 0.1 mM` runs 1.01 → 1.02. Every one of those rows carries an interval
  containing one at every value of `a`, and its verdict stays INCONCLUSIVE throughout.
- **The 2–5 mM DTT point estimates move.** `UPRE1 @ 5 mM` runs 0.93 → 1.06 and `UPRE2 @
  5 mM` runs 1.06 → 1.22, and three of those rows cross from INCONCLUSIVE to DIFFERENT
  inside the sweep: `UPRE2 @ 2 mM` at a = 20%, `UPRE1 @ 2 mM` and `UPRE2 @ 5 mM` at a = 30%.
  The same correction enlarges the resolved ones: `UPRE2 @ 1.0 mM` goes 1.56 → 1.71.

**What this sweep does not license, and what the registry says instead.** This section
previously read *"The 'no induction at all' calls at 0.1–0.2 mM are robust"* and *"The
'entirely dilution' calls at 2–5 mM DTT are not"*, treating insensitivity to `a` as evidence
that there was nothing to detect. Both readings are formally withdrawn.
`data/current_claims.json` carries `biosensor.no_induction.native_low`, `.native_second`,
`.altered_low` and `biosensor.entirely_dilution.upre1` at role `historical`, each with
`"margins": null` and `"margin_source": null`, and each requiring that justified
practical-equivalence margins be declared *before* the interval is interpreted at all. The
registry states the reason in its own words: **"no SESOI is invented by this registry"**, and
none is invented here either.

So the correct reading of every one of those rows is UNRESOLVED. An interval containing one
is a failure to resolve a direction, not a measurement of no response; a point estimate that
holds still under an assumed background is a stable *unresolved* number; and a nonsignificant
fold does not become equivalence by being repeated. Equivalence needs the complete interval
to sit strictly inside two margins that were justified in advance, and no such margins exist
for this panel.

<!-- audit:retracted "no induction at all" calls at 0.1–0.2 mM are robust -->
<!-- audit:retracted "entirely dilution" calls at 2–5 mM DTT are not -->

Those markers are the enforcement, not the decoration. Both sentences survived here for weeks
after the registry withdrew what they assert, because unmarked prose is invisible to
`audit_claims.py` — the scope policy says so directly: *"Unmarked prose is not discovered by
decimal, citation-year, test-count or keyword heuristics."* `check_retracted_claims` now fails
if either phrase reappears in any document but this one and `docs/superseded/`.

One plate of BY4741 read in the mCitrine channel would constrain `a`, which is worth having
for the 2–5 mM rows. It would not resolve the low-dose rows, because what blocks those is a
missing margin, not a missing background.

### The ER anchor did not merely fail to conclude; it failed as a measurement

From [`research/G4_STATISTICS.md`](research/G4_STATISTICS.md). The right unit is not cycles
but the fraction of the +RT signal that is genomic DNA, `2^(−ΔCq)`:

| Construct | Anchor | RT− margin | %gDNA | Revised verdict |
| --- | --- | ---: | ---: | --- |
| UPRE1 | HAC1 | +0.13 cyc | **91.4%** | **REFUTED as a measurement.** Transcript is 8.6% of the signal. Variance inflation 15.7×, and at literature-typical Cq noise a correction returns a *negative* transcript quantity about a third of the time. Past the 60% ceiling of ValidPrime, the only validated correction, by a wide margin. |
| UPRE2 | HAC1 | +0.60 cyc | **66.0%** | **REFUTED as a measurement.** Past the ceiling, though only just. |
| NativeYap1 | TRX2 | +8.6 cyc | 0.3% | INCONCLUSIVE, genuinely. Clean anchor, underpowered. |
| AlteredYap1 | TRX2 | +8.6 cyc | 0.3% | INCONCLUSIVE, genuinely. Closest to decidable. |

The distinction is not pedantic. For the ER pair **nothing was learned about the
biosensors** — the anchor assay did not measure transcript, which is a protocol finding, not
a biosensor finding. Only the oxidative pair is evidence about the constructs at all, and
there the required replicate count goes as the **square** of the inverse SNR, so the
productive move is a larger effect (a better anchor gene, or a dose inside the 0.5–1 mM
responsive window) rather than 80 more plates.

A useful correction to received wisdom while we are here: **MIQE never mandated a 5-cycle
RT− margin.** That number is a software default whose real content is "≤3% background", and
it originated as a *no-template-control* rule.

## G-d and G-e: what is unreachable, what is unbuilt, and what already works

These three are different states and the file previously ran them together. Being blocked by
a design is not the same as being blocked by a method, and neither is the same as not yet
having been written.

**G-d — generalising the state to an unseen stressor — is unreachable on the plates that
exist, and the obstruction is the design, not the estimator.** The wet-lab data carries two
stressors, and each is read only by its own dedicated construct pair: UPRE1/UPRE2 read DTT,
NativeYap1/AlteredYap1 read H2O2. Holding either stressor out removes every channel that
reads its private modules, leaving at most the two channels dosed with the other agent —
the zero-peer case the redundancy law says recovers nothing. So the cross-stressor number
cannot be *computed* on these plates. It does not evaluate to zero; there is nothing to
evaluate. `SPLITS.md` records the same fact as `heldout_stressor` being available on the
simulated panel and not on the real one, and the two datasets therefore support disjoint
claims.

Where the test *could* be run it was: 18 (held-out stressor, own target) pairs on Gasch
GSE18, arrays nobody in this project measured. The fitted basis **fails** its
rotated-subspace null there — +0.688 against a null median of +0.695, p = 0.990, skill
−0.025 — the same verdict the simulated version gets and by a wider margin
([`EXTERNAL_VALIDATION.md`](EXTERNAL_VALIDATION.md) §1). That is a real, held-out negative
result and it stands.

**G-e — forecasting the product — needs a kinetic or expression layer, and that layer is
unbuilt rather than blocked.** Three independent lines say the forecast cannot come from the
stoichiometry: D1 measured the feasible range at `[0, ceiling]`; the regulation layer showed
structurally that no capacity bound can lift the floor (Tier 1(b) above); and the field puts
its own forecasts elsewhere — COSMIC-dFBA (PMID 38387677) forecasts in an ML layer with dFBA
only propagating it, and Elizondo & Saa 2025 (PMID 40891387) is titled *"Complex Kinetic
Models Predict β-Carotene Production and Reveal Flux Limitations in Recombinant
Saccharomyces cerevisiae Strains"*. What has changed is that **the calibration set now exists**:
Elizondo & Saa report β-carotene and lycopene production rates, growth rate, glucose uptake
and ethanol/acetate/glycerol secretion for three producing CEN.PK2-1c strains at two
dilution rates, six independent (strain, µ) conditions in the paper's own units
([`research/CALIBRATION_DATA.md`](research/CALIBRATION_DATA.md)). It did not exist when
`run_d1.py` was written. Two limits remain and are real: those are six *steady states*, not
a time course, so they license the task ordering and a lifted product floor but not a
validated dynamic forecast; and no strain in this project makes the product
([`PARKED.md`](PARKED.md)).

**Held-out dose prediction, by contrast, works and is scored.** It is the one forward
prediction in the repository that has been run on real data at all, and it splits on dose:

| held out | dose / EC50 | test rows | skill vs. carrying the nearest measured dose forward |
| ---: | ---: | ---: | ---: |
| 0.1 mM | 0.7 | 18 | **+0.086** |
| 0.2 mM | 1.5 | 36 | **+0.368** |
| 0.5 mM | 3.7 | 18 | **+0.101** |
| 1.0 mM | 7.4 | 36 | **−0.683** |
| 2.0 mM | 14.8 | 18 | −0.178 |

Positive at every interior rung at or below 0.5 mM, negative above it, with the sign
separating on `dose / EC50` and a mechanism that predicts the direction before the data is
seen ([`WHY_INTERPOLATION_WORKS.md`](WHY_INTERPOLATION_WORKS.md)). It is **not Tier 3** —
scored after the measurement, not registered before it — and it is a modest win, not a
decisive one: an earlier version read +0.164/+0.429/+0.236 on two replicates and halved
when a third was recovered.

## What would move things up a tier

In cost order, cheapest first.

1. **One plate of BY4741 read in the mCitrine channel.** Decides three of the eight
   dose-response rows. Turns a Tier 0 constant into a measurement.
2. **Configure the reader to fluoresce all 12 columns, and put blanks in every channel
   read.** Costs three wells. Stops the n = 2 problem recurring.
3. **Recover UPRE1's replicates 2 and 4 under a background sweep.** Analysis only; takes
   UPRE1 to n = 4 with a stated bound rather than a silently assumed background.
4. **Register a forecast before running it.** `research/PREDICTION.md` §5.7 is a
   pre-registration checklist and §5 specifies the protocol. This is the only route to
   Tier 3, and it needs no new bench work — the data is already on disk.
5. **Challenge the declared generator distribution and report the objective a UCBOG bounds.**
   `research/CIRCULARITY.md` §7–8. Family checks and a dimension/design-choice bound now exist
   in `CROSS_FAMILY.md`; they are not a bound on all fitted weights or unknown biology.
   Mak/Morton/Wood's monotonicity result concerns expected sample-optimum bias under iid
   sampling, not every finite-run score or bootstrap bound.
6. **Two more full biological replicates**, which puts every row of the dose-response table
   above an interval.

## Provenance

Assignments made 2026-08-26 against the README as it then stood, the source tree under
`src/ystwin`, and the five research notes linked above. The audit that checks the README's
*factual* claims against the tree — module names, script paths, test counts — is
`scripts/audit_claims.py`; this file is the other half, and it is not automated.

Revised the same day for four results that landed after the first pass: the regulation
layer's Tier 1 negative on the product range; the stress-to-ATP coefficient, which was
wrong by 10⁴ and is now an asserted envelope; the reference phenotype, which was not in the
paper it cited and which reverses two published findings; and the G-d / G-e / held-out-dose
split above, added because "blocked", "unbuilt" and "already works" were being reported as
one state.
