# ystwin — reading gene activity out of glowing cells

Engineers put **fluorescent reporters** into cells: wire a gene of interest to a protein
that glows, and the brightness tells you when that gene switches on. It's one of the most
common measurements in biology.

It's also, done naively, wrong — for a reason this whole project exists to handle.

**One sentence for what this does:**

> Give it one well's brightness and cloudiness readings over time, and it tells you how
> hard that gene was actually switched on — with an honest error bar, and a refusal when
> the data isn't good enough to answer.

Nothing in that sentence is specific to an organism or a product. The correction it
performs applies to any dividing cell carrying any reporter, and that generality is the
point: the same maths reads a stress sensor in yeast, a promoter library in *E. coli*, or
an inducible circuit in anything else that grows.

What *is* specific is named as such and kept behind the general part — the yeast stress
catalogue, the yeast metabolic model, one plate reader's file format, and product workflows
with their own calibration and validation limits. Those are worked examples, not the subject.
A catalogue entry is not automatically a mechanistic pathway or a learned biological state.

We test it on baker's yeast (*Saccharomyces cerevisiae*) engineered to glow under stress,
because that's the data we have.

---

## The problem, in plain terms

We put a **fluorescent reporter** into yeast: a stress gene is wired to a protein that
glows yellow. More stress → gene switches on → more glow. Add a stressor (we use DTT and
hydrogen peroxide), read the plate every half hour, and watch the glow rise.

That's the idea. Here's the catch that this project exists to deal with.

**The glow builds up when cells grow slowly, even with no extra stress.**

Think of it like adding food colouring to a glass of water while someone is also topping
the glass up. If they top it up quickly, the colour stays weak — it's being *diluted*. If
they slow down, the colour looks stronger. Nothing changed about how fast you added the
colouring. The glass just stopped filling.

Cells work the same way. Every time a cell divides, it splits its glowing protein between
two daughters. Fast growth = strong dilution = dim. Slow growth = weak dilution = bright.

And here's the trap: **stress makes cells grow slower.** So a stressor makes the glow go
up *twice* — once because it really switched the gene on, and once just because the cells
slowed down. Read the raw brightness and you can't tell those apart.

This model separates them. The equation is:

```
change in glow  =  how hard the gene is on  −  (growth rate + decay) × current glow
```

Rearrange it and you can back out "how hard the gene is on" from readings you can
actually take. That number — not raw brightness — is what everything else is built on.

---

## What we found

Six things, in order of how much they matter.

### 1. Most of the "sensor response" was dilution

When we corrected for growth, a lot of the dose response disappeared.

| Sensor | Dose | Looks like | Actually is | Resolvable? |
|---|---|---|---|---|
| NativeYap1 | 0.1 mM H₂O₂ | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->1.28× brighter | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->0.92× — no resolved induction; equivalence untested** | no |
| UPRE1 | 2.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE1;stressor=DTT;dose_mM=2.0;late_fraction=0.75" -->2.02× brighter | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=2.0;late_fraction=0.75" -->1.06× — compatible with dilution; equivalence untested** | no |
| UPRE2 | 1.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.93× brighter | <!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.56× — real, and this one resolves | **yes** |
| AlteredYap1 | 0.5 mM H₂O₂ | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.5;late_fraction=0.75" -->1.69× brighter | <!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.5;late_fraction=0.75" -->1.37× — real, and new at the fourth repeat | **yes** |

Growth correction separates apparent brightening from inferred promoter activity.
An interval containing unity is unresolved, not evidence of “no response.” That claim
requires a predeclared equivalence margin and an interval lying entirely inside it.
The displayed intervals are conditional on the observation and bootstrap assumptions.

The last column is the part it is easiest to skip. UPRE2 at 1.0 mM is a call — its
interval is [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.30, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.83] and clears 1.0. A different sensor, AlteredYap1 at 1.0 mM H₂O₂, has a point estimate of <!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=1.0;late_fraction=0.75" -->1.00, but its interval is now **refused** because the bootstrap produces too many invalid ratios.
The old interval endpoints are withdrawn. Missing inferential support is not a precisely
flat response or evidence of practical equivalence.

**Seven of these are provable, and it was five before the fourth repeat.** The fourth plate
was excluded because every spreadsheet export of it had already had its blank subtracted
away; the instrument's own file still had the background, and reading it directly put the
plate back. UPRE1 clears "no change" at 0.2, 0.5 and 1.0 mM; UPRE2 now clears it at 0.2,
0.5 **and** 1.0 mM; and **AlteredYap1 at 0.5 mM H₂O₂ is the first per-dose call either
oxidative sensor has ever produced.** Nothing that was a call stopped being one.

**UPRE2 at 0.5 mM is the one to watch, because it has now been a call twice and not a call
once.** It cleared at n=2, lost the call when the estimator's smoothing window was measured
rather than inherited (its lower bound had been sitting at 1.014, fourteen thousandths clear
of nothing), and clears again at n=4 with a lower bound of
<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE2;stressor=DTT;dose_mM=0.5;late_fraction=0.75" -->1.14 — and now at every window we sweep rather than at one. That history is the
argument for reading the interval and not the fold. Details in
[`docs/FINDINGS.md`](docs/FINDINGS.md) and
[`outputs/estimator_accuracy.csv`](outputs/estimator_accuracy.csv).

**Multiplicity and resampling require separate checks.** The family contains twenty-four
requested folds, including unestimable cases. Multiplicity cannot be reduced by dropping
failed cases. The minimum two-sided **sign-test** p-value with four independent plates is
**0.125**; this bound applies to that test, not to every possible inferential method.

The actual bootstrap resamples both plates and wells. Its attainable quantiles and coverage
must therefore be evaluated under that two-stage procedure, not inferred from plate
multisets alone. Nominal intervals, multiplicity-adjusted decisions and equivalence tests
are reported separately. A tiny bootstrap tail estimate does not by itself establish
calibrated coverage or turn unresolved induction into “no response.”
[`outputs/fold_multiplicity.csv`](outputs/fold_multiplicity.csv).

**What survives most easily is the question asked once per sensor instead of six times.**
"Does this reporter respond to dose at all" is one test, and it can be answered by permuting
dose labels *within* each plate — which holds the plate fixed rather than estimating it from
four. All four sensors pass at p ≤ 0.0004, Bonferroni-corrected, **including both oxidative
sensors.** So the honest summary of the biosensor work is that the evidence these reporters
track dose is strong, and the evidence that any particular dose gives any particular fold is
weaker — though no longer absent, which is what the fourth plate bought.

How well the rest survives being poked is worth stating exactly. All of them keep their
direction across every measurement window we tried, and all survive assuming up to 30% of
the glow is the cells' own background. The stronger claim — that the error bar stays clear
of "no change" at *every* window, not just the default one — now holds for **all seven**,
where at three plates it held for four of five. The row that used to be the exception,
UPRE2 at 0.2 mM, reaches
<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE2;stressor=DTT;dose_mM=0.2;late_fraction=0.9" -->1.046 at the widest window instead of falling through it. The high doses still come
back **"can't tell"**, which fits: that is where the chemical is killing the cells.

### 2. The dataset is the size the design says, and three times we were wrong about why not

The design says four repeats. **All four are now usable, for every sensor.** Getting there
took three corrections, and every one of them was a fault in how we read the files rather
than anything wrong with the cultures.

| Sensor | Plates with brightness | Plates actually usable |
|---|---|---|
| UPRE1, UPRE2, NativeYap1, AlteredYap1 | 4 | **4** |

This section has said **two** and then **three**, with these explanations, all wrong:

- *"Two plates read only columns 1–3."* They didn't. All four carry 87 wells. Two of them
  write the plate **down** the page instead of across it — one block per group of columns,
  only the first block carrying the Time column — and our parser stopped after 21 wells.
- *"Two plates have no blank well in the glow channel."* They do. The blanks live in a
  later block, which is the same fault seen from the other end.
- *"The fourth plate's background is gone and cannot be recovered."* It is gone from every
  **spreadsheet** of it: the operator exported through a blank-subtraction transform, so
  20260804's H1–H3 read ±0.002 where the rest of the plate spans 0.09 to 0.61. It is not
  gone from the plate reader's own experiment file, which stores the untransformed
  readings. There H1–H3 read **0.087/0.090/0.088 OD and 347/349/354 RFU**.

Recovering 20260728 took the panel from two usable repeats to three. Reading 20260804's
instrument file takes it to four, and the arithmetic between the two copies is checked
rather than trusted: subtracting the per-timepoint mean of H1:H3 from the instrument file
reproduces the blank-subtracted spreadsheet **exactly** in the glow channel — all 2,175
readings, every disagreement either zero or the third of a count that rounding three
integers to one must produce — and to 0.0005 in the cloudiness channel, which is the most a
three-decimal export can say about a four-decimal reading. Both copies are committed, and
[`tests/test_xpt_export.py`](tests/test_xpt_export.py) re-runs that comparison from the
committed text alone, with no instrument file present.

**The n=2 table stood on this page after the repository's own
[`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md) had retracted it** — which is a
documentation failure of exactly the kind [`scripts/audit_claims.py`](scripts/audit_claims.py)
exists to prevent, and could not see, because that table carried no marker. It does now.

**One plate is one repeat, and that is now enforced rather than assumed.** The fourth plate
is committed twice — as the instrument file and as the blank-subtracted spreadsheet, because
the second is the evidence of what the first still has — so something has to decide which
copy is the biological replicate. `run_sensor_characterisation.one_file_per_plate` does, by
the logbook plate the two share, and refuses rather than ranking when two copies of one
plate are both spreadsheets.

Nothing was ever wrong with any of those cultures.

### 3. A control measurement was made, and the export threw it away

To know how much of the glow is the *reporter* and how much is the cell glowing on its
own, you measure a plain yeast strain with no reporter in it. **That plate was run, on
2026-08-07, and it was read in the glow channel** — the raw instrument file carries all 96
wells at the same gain, interval and optics as the biosensor replicates. The `.xlsx` export
dropped the channel, so for three weeks every document here said the measurement had never
been taken. There is **one** such plate, not two: the second file is the same run exported
again.

This section previously read *"Two such plates exist. Both were read for cloudiness only —
never in the glow channel."* Both halves were wrong, and what found them was reading the
instrument file instead of the spreadsheet.

It matters because it changes conclusions. Sweeping across plausible values for that
background:

- the low-dose rows stay put (they move by ≤0.02) and stay unresolved, equivalence untested ✅
- three high-dose rows **change verdict** inside the sweep — UPRE1 at 5 mM runs 0.93 → 1.06 ⚠️

A first pass over the recovered channel puts the reporter-free per-cell signal near **249
RFU/OD**, about 6% of the median reporter well — which, if it survives being done properly,
puts the background at the *bottom* of the swept range, where those three rows stay
unresolved rather than crossing. It is not done properly yet: it splices two plates read eleven days apart, and
the control plate has no medium-only well that fluoresces, so the subtraction can only make
249 an over-estimate. **Turning it into a per-construct `a` is now a data-reduction task,
not a plate to run** — `docs/research/XPT_INVENTORY.md` §2.2 carries the layout and the
evidence.

### 4. The ER sensor's independent check was measuring the wrong thing

To confirm a glowing sensor really tracks stress, you check it against a separate method —
here, qPCR, which counts how many copies of a gene's message the cell made.

Every qPCR run includes a control that should come back **empty**. Ours didn't. It came
back almost as strong as the real sample, which means what we measured was mostly
**leftover DNA**, not gene activity.

| Sensor | Contamination in the gene being measured | …and in the gene it is measured *against* |
|---|---|---|
| UPRE1 | **91.4%** genomic DNA | 14.2% |
| UPRE2 | **66.0%** | 6.9% |
| NativeYap1 | 0.3% | **9.8%** |
| AlteredYap1 | 0.3% | **9.5%** |

The only published method for correcting this works up to 60%. We're past it. **Nothing
was learned about the ER sensors** — that's a protocol failure, not a sensor failure.

**The second column is the part an earlier version of this page left out, and it changes the
word "clean".** Every one of these numbers is a *ratio* to a reference gene, UBC — and the
UBC channel carries 7–14% leftover DNA in all four constructs. So the oxidative sensors are
clean in the gene we care about and are being divided by something that is not. The effect is
small — [`docs/research/G4_ANCHOR.md`](docs/research/G4_ANCHOR.md) puts it at roughly a 15%
phantom fold — but "clean", unqualified, was the wrong word for it. The oxidative sensors are
still the usable pair; they just need a bigger effect to measure than we had power for.

### 5. Four colours do not identify an unrestricted stress network

The original design compared four channels with seven selected stress programmes, not an
exhaustive inventory of yeast stress biology. The broader generator now catalogues 24
entries: 19 transcriptional programmes and five metabolite pools, across 25 stressors.
Those overlapping entries are not all mechanistic or independently identifiable states.

For the fitted, static factor-analysis problem considered here, the *Ledermann bound*
limits four channels to **one** generic factor under its assumptions. This is not a universal
bound on supervised state estimation with supplied labels, loadings or dynamics. The newer
teacher/student loop deliberately supplies informative virtual reporters for three synthetic
aggregates; recovering them does not demonstrate that four real lab reporters resolve yeast's
stress network.

Two ways forward that *aren't* blocked, both in
[`docs/research/IDENTIFIABILITY.md`](docs/research/IDENTIFIABILITY.md).

### 6. The "hidden stress state" idea is not ours, and we say so

Reading several stress sensors at once and boiling them down to a small hidden state has
already been done in yeast, and an earlier version of this page did not say so.

- **2018.** Ten fluorescently tagged stress regulators, decoded cell by cell into an
  internal picture of what the cell is being hit with (Granados et al., PNAS,
  PMID 29784812).
- **2013.** About 900 yeast promoters measured with fluorescent reporters. Between 60% and
  90% of them shift from one condition to the next by a **single shared scaling factor**
  that depends on the condition and not on the promoter. Most of what a big sensor panel
  sees is one knob, not many (Keren et al., PMID 24169404).
- **2015.** A number on the ceiling: one such sensor carries a little over one bit —
  enough to say *which* stress, not *how much* (Hansen & O'Shea, PMID 25985085).

So the honest version is: **a hidden stress state is well founded as a way of telling
stresses apart, and poorly founded as a way of measuring how bad one is.** What would be
new here is the *joining* — sensor panel to hidden state to a metabolic model — and the
argument about how many hidden causes a panel this size can honestly resolve.

**Our own results agree with theirs, and that is the strongest part of it.** When we threw
away our fitted axes and used random ones of the same rank, the random ones did just as
well — twice: once in simulation, and once on public data we did not collect. The reason is
the one the 2013 paper predicts. A single component carries **51%** of the variation, and it
is the general stress response, which every stressor sets off. Any random rotation keeps
that component, so it scores the same.

Full prior-art list and the nulls: [`docs/FINDINGS.md`](docs/FINDINGS.md) and
[`docs/research/STRESS_RULES.md`](docs/research/STRESS_RULES.md) §2.

---

## What the model can and can't do

We grade every claim so nobody has to guess:

| Tier | Meaning | Example |
|---|---|---|
| **0** | A number we picked because it seemed reasonable | ATP cost per unit of gene activity — which was recently found to be wrong by a factor of ten thousand, and is now an envelope rather than a guess |
| **1** | Works in simulation, tested against our own simulator | All the power and transfer tables |
| **2** | Survives when we change the simulator's structure | Two physiology checks |
| **3** | Predicted first, measured after | **Nothing yet** — but [`outputs/registered_prediction_D018.csv`](outputs/registered_prediction_D018.csv) is now on the table: three strains at D = 0.18 /h, written down before any measurement exists there. It becomes Tier 3 if someone runs the chemostat, and stays Tier 0 if nobody does. |

That last row is the honest headline. **No claim here has been registered as a prediction
and then checked against a new measurement.**

We also tested the model on itself and it failed. Asked "how confident are you?", it says
it's far more certain than it turns out to be — and the pattern of its errors shows the
*equations* are missing something, not just the error bars. Retuning the error bars would
hide it rather than fix it. Details in
[`docs/research/PREDICTION.md`](docs/research/PREDICTION.md).

**A note on scope.** Our own yeast don't make beta-carotene yet — the wet-lab side of this
project is the biosensor work above. The product side is calibrated and scored against
someone else's published chemostat strains (Elizondo 2025), and it is clearly labelled as
that wherever it appears.

---

## Dynamic physical simulation

`predict_protocol` runs one amount-based reactor clock with signalling, expression,
metabolic pools, growth, death, washout and CrtE/CrtI/CrtYB product dynamics. Protocol
adapters in the culture, dynamic-FBA and hybrid modules use this same engine. Optical
observation changes do not change physical truth.

This example deliberately opts into **unvalidated structural and kinetic priors**:

```python
from ystwin.predict import (
    Control, EngineParameters, Protocol, ProtocolGenotype,
    initialize, predict_protocol,
)

parameters = EngineParameters.prior()
genotype = ProtocolGenotype.prior(product=True, reporter=True)
initial = initialize(
    parameters, genotype, volume_l=1.0, biomass_gdw_l=0.05,
    medium_mM={"glucose": 111.0, "nitrogen": 10.0, "oxygen": 0.2, "osmolyte": 100.0},
)
protocol = Protocol(
    (0.0, 0.5, 1.0),
    (Control(0.0, oxygen_transfer_per_h=10.0, oxygen_saturation_mM=0.2),),
)
result = predict_protocol(protocol, genotype, initial, parameters)
print(result.validity)
```

The result exposes units, compartments, inventories, assumptions and unsupported scope.
It has no flat genotype-independent carotenoid ceiling. The optional paired
acetate/proton-export prior is not a fully identified Pma1 regulation model. Supported
thermodynamic rate constraints require explicit activities and usable energies; they do
not certify uncovered reactions. Complete charge/electron, phosphate/CoA and calibrated
whole-cell prediction remain unresolved. Parameter/source evidence is versioned in
`data/parameter_evidence.json`; software conservation tests are not biological validation.

For state inference, `mech.inference.PhysicalParticleFilter` propagates complete physical
states through this same engine and observer. It requires an explicit ensemble, weights,
parameters, genotype, calibration and measurement-noise standard deviations. Resampling
preserves event/history and ancestor identities; forecasts do not mutate the filter.
Its uncertainty is conditional on those supplied assumptions, not independently calibrated
biological coverage. The older four-state estimator remains a separate comparison model.

## Empirical product comparison

The older expression-to-flux pathway model remains an explicitly empirical comparison,
not the physical engine's default law.

```
genotype: one RT-qPCR number   ->  pathway flux         PREDICTS
environment                    ->  growth rate          BOUNDS
stressor and dose              ->  does a steady state exist at all
flux + growth rate             ->  every pool in the chain   PREDICTS
flux + genome-scale model      ->  feasible? what does it cost?   AUDITS
```

**This used to be circular and now isn't.** The pathway flux was a number you had to hand
in, and both entry points obtained it by adding up the measured product — the model measured
beta-carotene in order to predict beta-carotene. It now comes from relative expression of
the pathway's first heterologous enzyme.

The product evidence comprises **three held-out strains and six condition predictions**,
not six independent strains. β-carotene and lycopene must be scored together. Gene choice,
model choice, entry-flux calibration and branch kinetics belong inside the training folds;
constant-content and constant-rate predictors are necessary comparisons.

The earlier fixed-CrtE analysis reported median errors of 14.2% for β-carotene and 26.3%
for lycopene. Those are **historical, retrospectively selected-model scores**, not the
performance of a nested gene-selection procedure. Current results must come from the
regenerated fold records produced by `scripts/predict_product.py`.

The expression-law baseline previously used the full-data geometric mean despite being
labelled training-only. Correcting it invalidates the statement that every native gene
has negative skill: ERG9 is slightly positive. The heterologous genes still rank highest
in the descriptive comparison, but that ranking does not identify cassette dosage as the
causal explanation or exclude measurement artefacts. Gene selection is evaluated within
folds rather than justified by the same outcomes used to rank genes.

**Capacity and allocation are different questions.** In the tested capacity-only
formulation, pathway upper bounds leave zero product flux feasible. This does not establish
that every FBA formulation is incapable of predicting allocation. Metabolic task policies,
state mixtures and mechanistic allocation add distinct assumptions that require their own
validation. Host-model identity and condition-matched molecular uptake constraints must be
recorded before comparing feasible product ranges or growth costs.

**Scope: products that stay inside the cell.** The maths balances what the pathway makes
against what growth dilutes. A secreted product leaves through a transporter instead, and the
same equations return a number that means nothing — so a pathway declared as secreted is
**refused at load time** rather than answered. Beta-carotene, lycopene, PHB and squalene are
in; farnesene and resveratrol are not.

**The pathway structure is declared as data.** `data/pathways/<product>.toml` supplies the
chain, rate-law choices and molar masses. Calibration and GEM installation can still require
product-specific code; the general pathway solver is not the whole product workflow.

---

## Synthetic teacher/student laboratory and stress-map scope

`scripts/run_in_silico_loop.py` fits a student to generated state, derivative and control
labels. It learns **three synthetic aggregates** (UPR, oxidative stress and burden) driven
by **two imposed inputs**, observed through **four virtual reporters** plus cell density.
The encoder, dynamics and controller weights are fitted; product amounts are not fitting
inputs. The state definitions, informative sensor design, control-law family, GEM and
native enzyme priors are supplied assumptions, not learned biological mechanisms.

This is narrower than the broader stress catalogue above. Only DTT/UPR and H2O2/oxidative
have optional mechanistic reporter routes in `generator/panel_experiment.py`; those routes
are not the new student's virtual sensor panel. The **15-branch conceptual atlas** in
`data/gem/stress_response_map.json` adds no kinetics. Its audit queries actual ORF/GPR
associations in each supplied GEM without optimization. Metabolic chemistry or GPR presence
is not evidence that signaling, repression or a full stress response is implemented.

The assessment separates `teacher_recovery_passed`,
`cross_equation_generalization_passed`, and `biological_validation` (always false here).
The frozen `closed_loop_04` review passed same-teacher recovery but failed the changed-equation
challenge under the same 0.01 forecast-RMSE criterion; both errors are recorded in the review.
Small teacher/student product errors cancel shared model errors; they are not biological
validation. The earlier local review is not included in a clean checkout and is not
public numerical evidence. Reproduce an assessment from a newly generated run using the
commands below; keep its model/configuration identity with the report.

```bash
python3 scripts/run_in_silico_loop.py --output-dir /tmp/ystwin-learning-example --skip-products
python3 scripts/audit_stress_map.py --output-dir /tmp/ystwin-stress-map-example
python3 scripts/run_in_silico_loop.py --assess-run /tmp/ystwin-learning-example \
  --assessment-output /tmp/ystwin-stress-map-example/teacher_recovery_review.json
```

The audit needs both GEM files; the review needs the existing local run and a new output
file. See [`docs/REPRODUCING.md`](docs/REPRODUCING.md#10-in-silico-teacherstudent-loop)
for artifacts and execution semantics, and
[`docs/research/CIRCULARITY.md`](docs/research/CIRCULARITY.md) for leakage, shared-model
assumptions and the limits of simulation optimization bias.

## Published HOG model and measured-data fitting

A separate route, `scripts/run_hog_learning.py`, executes the published
Petelenz-Kurdziel et al. 2013 osmoadaptation SBML: **29 dynamic states and 58 reactions**,
including Hog1, Gpd1, glycerol, Fps1 regulation and volume feedback. Original model/data
files, checksums, the 2014 mutant-caption correction, and source assumptions are recorded
in `data/hog2013/`. This does not silently replace the synthetic teacher above.

The loader preserves **407 observations** with source cells and measurement provenance.
The frozen protocol fits **13 WT phospho-Hog1 observations at the calibration NaCl dose**
and reserves **90 preliminary Western observations** that the source explicitly says were
not used for parameter fitting. They are not independent of the paper's model-development
process. Exact doses and source cells are in `data/hog2013/fit_protocol.json`.

The data constrain an **activation/deactivation balance**, not two separately resolved fast
rates. The common fast timescale remains a prior. The positive observation scale is estimated
from training data, not mistaken for a directly measured absolute phosphorylation fraction.

The refit **does not improve the stronger-dose prediction** over the published parameters,
although both dynamic models outperform persistent activation. Exact parameter estimates and
measurement-block waveform errors are emitted in `fit.json` and `validation.json`, alongside
per-trace `heldout_scores.csv` in the requested output directory. The earlier local run
labelled `hog_learning_run_03` is not included in a clean checkout; use the public source
and explicit reproduction command below rather than treating that local path as evidence.

Gpd1 abundance updates the existing GEM enzyme budget without inventing a flux: the demonstrated
bound increases, while the chosen glucose-only glycerol-capacity probe remains unchanged.
Retention and secretion have distinct product tasks. This is not validated titre prediction or
discovery of all stress-response laws; most kinetic parameters and the imposed OD trajectory
remain priors.

```bash
python3 scripts/run_hog_learning.py --output-dir /tmp/ystwin-hog-example --gem-check --plot
```

Use a new or empty destination. See [reproduction section 12](docs/REPRODUCING.md#12-published-hog-model-and-measured-data-fitting)
for source acquisition, numerical checks, artifacts and remaining limitations.

## Product-blind native training and transfer

`scripts/run_native_training.py` fits native glycerol synthesis/diffusion and empirical
native growth/exchange curves, then freezes source files, model parameters and the executable
prediction dependencies. It does not import product-specific pathways or fitted product
constants. A scoped read audit records the native inputs actually used.

`scripts/evaluate_frozen_transfer.py` accepts new chemistry and experimental conditions only
after that freeze. Its generic chemical installer checks reaction balance and explicit
compartment roles; it never switches behavior on a product name. Native glycerol commitments
consume actual carbon/cofactor resources, and learned molecular uptake caps cannot reopen
closed exchanges or relax hard limits.

The lycopene transfer study is deliberately qualified: the selected chemostats measure an
intermediate in a carotenoid pathway, fungal desaturase redox chemistry is not independently
resolved, and missing catalytic/consumption information does not identify a point forecast.
The strict verified-chemistry lane remains blocked. Separately registered, explicitly unverified
redox hypotheses use the same frozen native model and produce sealed feasible envelopes before
outcome scoring. All variants and infeasible cases are retained; no best chemistry is selected
using the outcomes. The completed comparison is recorded in
`outputs/heldout_transfer_summary.json`: some conditions fail native feasibility, and the
remaining envelopes are too broad to establish accurate lycopene prediction.

The model remained outcome-blind, but one data curator accidentally encountered product-rate
values after model freezing. That incident is retained in the audit trail, not erased by the
independent condition-only verification. This is not claimed as fully blinded operator knowledge.
Source/evaluation qualifications live in `data/holdout_transfer/`; the native checkpoint is in
`outputs/native_training_run_02`. See [reproduction section 13](docs/REPRODUCING.md#13-native-training-before-an-unseen-product-test)
for the freeze, prediction and scoring contracts and reporting limitations.

## Biology-learning audit and native consistency

The native learner uses real measurements, but that does not establish that its supplied
kinetic laws, cross-model coupling or allocation rules were learned. The deeper review is in
`outputs/biology_learning_review.json`.

`scripts/run_native_reconciliation.py` now separates baseline host/data inconsistency from
learned glycerol obligations. It tests primary printed precision, simultaneous native
exchanges and explicitly labelled energy/protein counterfactuals without changing the source
model. Rounding is not a confidence interval, and a feasible witness is not a realized flux.
No new physiological multiplier was adopted.

`scripts/audit_biology_learning.py` binds measured source cells and quantities, inventories
fixed versus estimated structure, and runs retrospective null/baseline and identifiability
checks. Its verified claim gate re-collects authoritative evidence; serialized success flags,
reference-shaped metadata and hashes alone cannot authorize biological claims. Stronger
unimplemented evidence grades fail closed.

`scripts/check_biology_readiness.py` combines fresh native consistency solves with that verified
claim gate. The current result is **not ready to claim biological-law learning**. Independent
enzyme evidence and an unopened external signaling-data catalogue are recorded under
`data/biochemical_evidence/` and `data/validation_candidates/`. Mixed-extract activity is not
intrinsic turnover, and nuclear localization is not phosphorylation concentration.
See [reproduction section 14](docs/REPRODUCING.md#14-native-reconciliation-and-biology-learning-readiness)
for executable checks and the distinction between expected claim refusal and a software failure.

## Relative orders and measured-cohort development

`analysis/partial_orders.py` represents supported relative orders without assigning invented
magnitudes. Quantity, calibration, assay context and source are explicit. Incomparability is
not a tie, temporal precedence is not causal dependence, and an earlier event is not necessarily
more important. Declared ties form quotient classes; inconsistent cycles are reported rather
than silently broken.

Order-only queries can preserve a direction while leaving its magnitude unbounded. Explicit
range assumptions enable validated affine LP bounds. Nonlinear forward-model sweeps report
sampled robustness only; a missing or zero-probability conditioning domain does not acquire an
imaginary uniform distribution. The native example in `outputs/native_order_summary.csv` shows
that oxygen-uptake order need not follow growth-rate order. Its rounding bands are assumptions
about reported precision, not biological confidence intervals.

The source-cohort development pipeline also completed its approved fit and comparison. It
uses actual published intensity ratios, a fixed earlier-frame cohort, training-only fits, and
sealed predictions before development-response release. The transient family was selected,
with its parameter-bound warning retained. This is **retrospective processed-cohort development**,
not independent biological-law validation, dose generalization or a certified prospective forecast.
The result is recorded in
`outputs/native_population_development/native_training_01_phase_c_01_scoring/development_report.json`.

The consolidated utilities include `scripts/ecmodel_isoprenoid_prior.py`, which exports model
turnover priors and a declared capacity sweep without treating them as enzyme measurements,
and the measurand-aware `published_cassettes.rows_above_content` query. Neither restores the
old unsupported dosage-scaling or localization claims. See reproduction sections 15 and 16.

The public evidence bundle in `data/frozen_evidence/native_v1/` supports replay without the
original machine-specific records. Its manifest has a separate pinned identity; declared storage
paths are neutralized while scientific content is retained. Original byte integrity is not
claimed when those local originals are unavailable. The replay checks source/code bindings and
recomputes recorded native predictions and score arithmetic; it neither fits again nor upgrades
the biological claim. See reproduction section 17 and
`outputs/consolidation_checks/public_replay_final.json`.

## Running it

```bash
python3 -m pip install -e ".[dev]"
python3 -m pytest tests/ -q          # the full suite
```

<!-- audit:test_count --> The suite is **10142 tests**.

That number is not decoration: `scripts/audit_claims.py` re-collects it and fails if this
line drifts. The marker is how it finds the claim -- an earlier version scanned the prose
for "N tests" and flagged its own changelog.

The analysis scripts need real plate files, which live outside this repo. Point them at
your copy:

```bash
export YSTWIN_PLATES="/path/to/Biosensor Testing"
python3 scripts/run_gates.py                    # which wells are trustworthy
python3 scripts/run_sensor_characterisation.py  # dose response, corrected
python3 scripts/run_calibration_nis.py          # is the model's confidence honest?
```

Without the data, everything skips cleanly instead of failing -- and that is enforced, not just intended: `tests/conftest.py` carries the skip for every
fixture backed by real data. Three files used to bypass it and put 49 tests
into ERROR on a clean checkout.
Setup detail: [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

**Checking our own work** — these are meant to fail when something drifts:

```bash
python3 scripts/audit_claims.py           # does this README match the code?
python3 scripts/audit_reproducibility.py  # could a stranger rerun this?
python3 scripts/audit_determinism.py      # same seed, same answer?
```

---

## Where things are

| | |
|---|---|
| [`docs/FINDINGS.md`](docs/FINDINGS.md) | The full research log — every result, every audit, every number we got wrong and fixed |
| [`docs/CONTRACT.md`](docs/CONTRACT.md) | Exactly what goes in and what comes out |
| [`docs/CLAIM_BOUNDARY.md`](docs/CLAIM_BOUNDARY.md) | Every claim, with its tier |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | How the code fits together — start here to contribute |
| [`docs/MECHANISTIC_LAYER.md`](docs/MECHANISTIC_LAYER.md) | `src/ystwin/mech/` — what it measures, what it refuses, and the three headline results that are negatives |
| [`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md) | What the plates actually measured |
| [`docs/PRODUCT_ARCHITECTURE.md`](docs/PRODUCT_ARCHITECTURE.md) | Predicting product X through pathway Y under environment Z — the chain end to end |
| [`docs/PROTOCOLS.md`](docs/PROTOCOLS.md) | Wet-lab protocols, each naming the function that consumes its result |
| [`docs/FIGURES.md`](docs/FIGURES.md) | The four figures that carry the argument, two of them negatively |
| [`docs/research/`](docs/research/) | Deep dives with citations |

---

## The idea behind all of it

**The model is allowed to say "I don't know."**

Most of this code is arithmetic. The part worth copying is that it refuses. A gate marks a
well unusable rather than analysing it anyway. A correction won't run without a
measurement it needs. An error bar comes back empty at two repeats instead of pretending
to be narrow.

That's why the results here are smaller than they first looked — and why what's left can
be trusted.

## How far this is from what it is meant to be

[`docs/DISTANCE_TO_THE_VISION.md`](docs/DISTANCE_TO_THE_VISION.md), and it is not
flattering. It scores the distance at **≈ 30%**, split three ways: the engineering is at
~85%, the scientific scope at ~30%, and the evidence at ~5%.

All five layers the goal names — GEM, FBA, regulation, latent stress state, thermodynamics —
are now **reachable from the prediction chain**, and `ProductPrediction.layers` reports on
every call which of them ran and what it did. That is a smaller change than it sounds, and
the third column is why: **only two of the seven can move the answer.** The GEM audit and
the thermodynamic gate can refute a number and never set one; E-Flux is provably inert
because every module a stressor touches is *induced* and E-Flux applies upper bounds only;
and the latent stress branch is computed, attached, and labelled with G4's refusal rather
than fed in. Assembling them is what demonstrated that, not what fixed it.

Of eight questions a user would really ask, **two get a real answer**. Three are correctly
refused. Three — ethanol instead of glucose, 37 °C, a stressor — return the *identical
number*, because those inputs never reach it.

What exists is **one product, one pathway, from genotype, at one dilution rate, one host,
one paper**. The environment axis is not partially working; it is refuted on the only
dataset that can test it, and no second dataset exists to rescue it.

**Nothing here is Tier 3**, and that is the row that dominates the score. One chemostat at
D = 0.18, against a prediction already written down in
[`outputs/registered_prediction_D018.csv`](outputs/registered_prediction_D018.csv), would be
the first — and would move the number further than any amount of code, because the
engineering column has under 15 points left in it.

## What happened when the model was attacked

[`docs/HARD_TESTS.md`](docs/HARD_TESTS.md). Three passes were told to make it produce a
confidently wrong answer. They found eight things, six were real, and **two of the six were
false claims in this repository's own documentation**.

The two that matter most: the model would happily report a pool weighing **36× the cell that
holds it**, and it turned out to assert a **hard ceiling of 1.25 mg/gDCW on β-carotene that
no genotype at any growth rate can pass** — which nobody knew it was making.

**They got different remedies, and this paragraph used to claim they got the same one.** It
read "Both are now refusals rather than numbers", and only the first is. An impossible pool
raises `ImplausibleContent` and returns nothing:
`predict_product` at an entry expression of 10⁴ refuses rather than reporting 30 g of
lycopene per gram of cell. The ceiling is **a note, not a refusal** — a prediction within
10% of it says so and hands the number over anyway, which is what
[`docs/HARD_TESTS.md`](docs/HARD_TESTS.md) §2 has always described. Asked for β-carotene at
D = 0.18 the chain still answers 1.247 mg/gDCW at an entry expression of 100, capped at a
bound the next paragraph says is wrong by 63×.

That is the defensible design — a note is what a model owes you when its own bound is known
to be too low, since refusing would withhold the only number it has — but it is not a
refusal, and calling it one overstated the repair in the one direction that flatters the
model. Found by the 2026-08-31 audit, which ran the chain rather than reading the sentence.

**The ceiling has since been refuted outright, by about 63×.** Arhar 2024 (PMID 39215465)
measured 79 mg/gDCW by HPLC on gravimetric dry weight with lycopene below detection — a
specific assay, not a total-carotenoid absorbance sum. This paragraph used to quote a "17×
below a published measurement" comparison; that number came from a total-carotenoid reading
this repository has since reclassified as non-specific, so it is withdrawn in favour of the
measurement that actually settles it. The repair is **structural** rather than a bigger
constant: capacity is set by crtYB dosage, which `Genotype` does not carry. See
[`docs/HARD_TESTS.md`](docs/HARD_TESTS.md) and
[`docs/EXTERNAL_CAROTENOID_BOUND.md`](docs/EXTERNAL_CAROTENOID_BOUND.md).

None of the failures were arithmetic. They were all at the **boundaries** — what the model
does outside its calibration, what it refuses, and what it claims about itself — which is
exactly where none of the other 2,400 tests were looking.

## What an audit found that the tests did not

[`docs/AUDIT_2026_08.md`](docs/AUDIT_2026_08.md). Every defect in it was found by *running*
the code. The suite stood at 2205 tests when the audit began and passed throughout it;
neither integrity gate flagged any of what follows.

The short version: one script had been unable to start since a refactor months ago, one
comparison branch was arithmetically incapable of producing a difference, two scripts
quietly rebuilt tracked tables from half their inputs and exited 0, both integrity gates
raised false alarms of their own, and the six numbers in the table at the top of this file
were a replicate behind. A green suite says the code does what its tests say. It says
nothing about a script that cannot run or a number that stopped being true.
