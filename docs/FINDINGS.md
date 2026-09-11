# ystwin — yeast stress digital twin, model v2

> **Given a dose and a well's OD and fluorescence time series, this twin returns a
> posterior over that well's growth rate and promoter activity — dilution removed,
> uncertainty attached, and refused outright where the optical channel is not
> quantitative.**

That is the contract. `docs/CONTRACT.md` derives it: the three ambition levels, which one
is committed to and which two are not, and the S / D / O / U / Q / R breakdown that
Kapteyn et al. (*Nat. Comput. Sci.* 1:337–347, 2021) use to say what a twin is. The short
version of the gap: **R, the reward, is undefined**, so this is a predictive twin and not a
prescriptive one; and Level 3 — a distribution over the β-carotene titre curve — is blocked
twice over, because no strain produces the product and because FBA cannot supply it.

Instrumented, tested implementation of the decision-critical layers of the
beta-carotene digital-twin plan, plus the diagnostics that decide what the rest of
the architecture may claim.

## yeast-GEM cannot reach van Hoek's growth under van Hoek's own uptakes, on nine rows of ten

**2026-09-09.** **Tier 1 — deterministic computation on an adopted model.** Nothing is
fitted here; it is only as good as yeast-GEM 9.0.2 and the transcribed Table 1.

Every dilution rate in `data/physiology/chemostatData_VanHoek1998.tsv` was put to plain
yeast-GEM under **that row's own measured glucose and oxygen**, both as equalities
(`constrain_uptake(..., mode="fixed")`), biomass `r_2111` maximised, GLPK. Efficiency is the
COSMIC-style ratio the task ordering computes: measured &divide; model maximum. Above 1 it is
not a low priority, it is a contradiction, and `task_priority_order` raises rather than rank it.

| D, /h | glucose | O2 | model max | efficiency | verdict |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0.025 | 0.30 | 0.80 | 0.02337 | 1.070 | refuses |
| 0.050 | 0.60 | 1.30 | 0.04554 | 1.098 | refuses |
| 0.100 | 1.10 | 2.50 | 0.09337 | 1.071 | refuses |
| 0.150 | 1.70 | 3.90 | 0.14720 | 1.019 | refuses |
| 0.200 | 2.30 | 5.30 | 0.19861 | 1.007 | refuses |
| 0.250 | 2.80 | 7.00 | 0.22898 | 1.092 | refuses |
| **0.280** | **3.40** | **7.40** | **0.29427** | **0.9515** | **REACHABLE** |
| 0.300 | 4.50 | 6.10 | 0.27783 | 1.080 | refuses |
| 0.350 | 8.60 | 5.10 | 0.33533 | 1.044 | refuses |
| 0.400 | 11.10 | 3.70 | 0.34598 | 1.156 | refuses — **the row `physiology.py` uses** |

Reproduced here independently of the pass that first reported it; every figure agrees to the
digits printed, `0.09338` excepted, which is `0.09337` at this solver's tolerance.

**This is not a data-selection bug.** There is no better row to pick. Nine of ten steady
states in the only chemostat series this repository has vendored are outside what the model
can do, so the question is not "which row" but whether yeast-GEM reproduces this dataset at
all. It does not.

### The refusal has two different causes, on either side of the critical rate

Lifting the two constraints one at a time separates them. Glucose stays fixed at the measured
value throughout; "yield" is `D / (q_glucose x 0.180156)`, in g biomass per g glucose.

| D, /h | as computed | O2 equality lifted | + maintenance deleted | measured yield |
| ---: | ---: | ---: | ---: | ---: |
| 0.025 | 1.070 | 1.065 | 0.935 | 0.463 |
| 0.050 | 1.098 | 0.996 | 0.935 | 0.463 |
| 0.100 | 1.071 | 1.055 | **1.020** | **0.505** |
| 0.150 | 1.019 | 1.012 | 0.990 | 0.490 |
| 0.200 | 1.007 | 0.992 | 0.976 | 0.483 |
| 0.250 | 1.092 | 1.015 | **1.002** | **0.496** |
| 0.280 | 0.952 | 0.934 | 0.924 | 0.457 |
| 0.300 | 1.080 | 0.754 | 0.748 | 0.370 |
| 0.350 | 1.044 | 0.459 | 0.457 | 0.226 |
| 0.400 | 1.156 | 0.406 | 0.405 | 0.200 |

**Above D = 0.28 the cause is oxygen, and only oxygen.** Lift the oxygen equality and the three
super-critical rows clear by a wide margin — 0.754, 0.459, 0.406 — while deleting maintenance
on top of that moves them by under one part in a hundred. van Hoek's q_O2 **falls** from 7.4 to
3.7 as the cells ferment, and yeast-GEM cannot make that much biomass on that little oxygen.
It is a structural disagreement about biomass per oxygen in a fermenting Crabtree-positive cell,
and no parameter in the maintenance account touches it.

**Below D = 0.28 the cause is the yield, and two rows are unreachable outright.** With oxygen
free and `r_4046` deleted, yeast-GEM's absolute glucose-to-biomass ceiling is **0.495 g/g** —
flat across the series, because at that point it is pure stoichiometry. Two of Table 1's own
implied yields are **above** it: 0.505 at D = 0.100 and 0.496 at D = 0.250. No setting of any
parameter in this model reaches those two rows. The other four sub-critical rows sit below
0.495 and are refused by the fixed-oxygen equality and the 0.7 mmol ATP/gDW/h maintenance
charge acting together; lifting both flips all four. So the sub-critical refusals are a mix —
two stoichiometric and unfixable, four the joint work of an equality and a fitted constant.

### D = 0.280 is not a coincidence, and it is not a code coupling either

The single reachable row is exactly `generator/culture.py`'s `CRITICAL_GROWTH_RATE_PER_H = 0.28`.
Nothing on this path reads that constant — it is imported only by `predict.py`,
`generator/context.py` and `run_environment_sweep.py`, never by `fba/` or `bridge/regulation.py`.
The match is a **shared cause**, and both effects are downstream of the same physiology:

- D = 0.28 is where van Hoek's q_O2 **peaks at 7.4**, the maximum anywhere in Table 1. The
  respiratory supply is at its largest, so the oxygen limit is at its weakest.
- D = 0.28 is simultaneously the first row whose yield **leaves the plateau** — 0.457 g/g
  against 0.483–0.505 on the four rows above it — because ethanol and acetate first appear
  here. It is the only row in the series where yeast-GEM's own yield on the row's uptakes
  (0.480 g/g) is **larger** than the yield the measurement asks for.

Those two pressures are minimal at the same row and nowhere else, and the reason they coincide
is the definition of a critical dilution rate: maximum respiratory capacity, at the onset of
overflow. The constant is 0.28 because van Hoek's cells switch there; the model reaches 0.28
because van Hoek's cells switch there. Causal, through the biology, not through the code.

### Five of the nine refusals survive the source's own printing precision

Stated because it bounds the claim. Table 1 prints fluxes to one decimal, so a printed `2.3`
is anything in [2.25, 2.35). Re-run at the **most favourable** reading — every q raised by
0.05 and both uptakes relaxed from equalities to caps — the refusals at D = 0.025, 0.150,
0.200 and 0.250 fall below 1 and are **inside the source's rounding**. Five survive:

| D, /h | efficiency, as printed | efficiency, best case | verdict |
| ---: | ---: | ---: | --- |
| 0.050 | 1.098 | 1.038 | refuses still |
| 0.100 | 1.071 | 1.042 | refuses still |
| 0.300 | 1.080 | 1.070 | refuses still |
| 0.350 | 1.044 | 1.036 | refuses still |
| **0.400** | **1.156** | **1.147** | **refuses still** |

The headline is nine of ten as computed; the part that no reading of the paper can explain
away is **five of ten**, and the row this repository actually uses is the worst of them. Under
the milder `mode="cap"` alone, without the rounding allowance, it is eight of ten.

### What it means for the COSMIC task ordering

`scripts/run_regulation.py` builds three measured tasks — biomass, ethanol, CO2 — from
`REFERENCE_AEROBIC_BATCH`, which is Table 1 at D = 0.40. `task_priority_order` refuses at the
biomass task with efficiency 1.156, and `regulation_matched_tasks.json` records a refusal in
place of an order. **That refusal is correct and should not be relaxed.** An efficiency above 1
means the denominator is smaller than the numerator it is supposed to normalise; ranking it
would be ranking a number that does not exist.

What this finding changes is the diagnosis. The refusal was readable as a bad row choice, and
it is not: it is yeast-GEM disagreeing with the whole dataset. Producing a task ordering is
therefore not a matter of repairing the experiment — it is a matter of deciding which of two
things to publish.

**Option A — move the reference to D = 0.280.** The only row the model can reach, at
efficiency 0.9515 with 4.9% of headroom. **This was run, not assumed.** Feeding
`measured_task_report` a D = 0.280 reference (glucose 3.40, O2 7.40, ethanol 0.11, CO2 8.00)
returns `status="optimal"` on both models and a complete order where D = 0.40 returns a
refusal:

| model | rank 1 | rank 2 | rank 3 | constrained growth |
| --- | --- | --- | --- | ---: |
| plain yeast-GEM 9.0.2 | biomass, 0.9515 | CO2, 0.9001 | ethanol, 0.3246 | 0.29427 |
| ecYeastGEM_batch 8.3.4 | biomass, 0.9368 | CO2, 0.9411 | ethanol, 0.3360 | — |

(The rank-2 efficiency can exceed rank 1's because each round recomputes the maxima under the
floors already imposed; the round-1 ranking is what fixes the order.) So the answer to "is the
ordering then produced" is yes, on both models, and it is the same order on both. The costs
are real and must be stated with it: ethanol at D = 0.280 is **0.11 mmol/gDW/h**, 73 times
below the CO2 task at 8.00 and 126 times below the D = 0.40 value, and van Hoek's footnote
*c* makes the sub-critical zeros non-detects rather than measurements — so the ordering is
computed at the very edge of the assay, where the ethanol rank is decided by a number near the
detection limit. It also moves the reference off the respiro-fermentative regime that the
β-carotene work cares about and onto the boundary of it. And it changes the five committed
constants: `bridge/maintenance_calibration.py:313` and `scripts/maintenance_scale.py:96` read
`glucose_uptake` and `oxygen_uptake` straight out of the reference and would silently produce
different numbers, while `bridge/stress_energetics.py`, `fba/stress_ph.py` and `mech/ph.py`
quote the D = 0.40 operating point in prose and recorded results that would go stale. The
honest reading of a 4.9% margin is that the model *just barely* clears one row out of ten, not
that it reproduces the physiology — the ordering above would rest on that margin.

**Option B — keep D = 0.400 and publish the discrepancy.** Nothing moves, and the refusal
becomes a reported result rather than a missing output: yeast-GEM, given the measured glucose
and oxygen of a published chemostat, cannot make the growth that was measured with them, on
nine rows of ten and on five of ten under the most generous reading of the source. The cost is
that the COSMIC §4.3 task-ordering experiment produces **no output at all** on this dataset,
and that stays true until a model or a dataset changes. The gain is that the strongest number
in the finding is the one that is currently blocking it — 1.156 at D = 0.40 — so the block and
the result are the same fact, and the experiment's absence is explained rather than unexplained.

**Unchanged pending that decision.** `REFERENCE_AEROBIC_BATCH` still carries D = 0.40 and its
five constants exactly as committed. Only the label moved: the reference is now
`REFERENCE_GLUCOSE_LIMITED_CHEMOSTAT`, with the old name kept as an alias because twenty-odd
modules and tests import it, and `aerobic_batch_constraints` is now
`chemostat_uptake_constraints` on the same terms. The paper ran no batch culture.

## Autofluorescence measured at last, and it is too small to matter

**2026-08-30.** `observation.py` insists autofluorescence be measured on an isogenic
reporter-free strain rather than guessed, and three documents recorded that as never done.
It was done. BY4741 was plated on 2026-08-07 under the same protocol, gain, interval and
optics as biosensor replicates 2-4, and read in both channels; the `.xlsx` export dropped
the mCitrine channel and `plate/gen5.py` recovers it.

**It is a slope, not a ratio.** With no reporter the model is a straight line in biomass,
`RFU = background + a*OD`. Taking the median of `RFU/OD` -- the obvious move, and what a
first pass did -- measures `a + background/OD`, which at the optical densities a four-hour
plate reaches is mostly background. That returned **249 RFU/OD**. The slope returns **18**.

One slope per WELL, then a t interval over the thirty wells that grew. Twenty-five readings
of one well are one well.

| | |
| --- | ---: |
| a | **+18.2 RFU/OD** |
| 95% CI over wells | **[-6.1, +42.4]** |
| includes zero | **yes** |
| as a fraction of a reporter well (4354 RFU/OD) | **0.42%**, upper bound **0.97%** |

**And that closes the sensitivity it was blocking.**
`outputs/autofluorescence_sensitivity.csv` sweeps this fraction at 0, 5, 10, 20 and 30% and
four rows change verdict somewhere in that sweep. The lowest fraction at which **any** verdict
changes is **5%** -- five times the upper bound measured here. **Every call in the panel
survives the measured autofluorescence.** That is the sweep's question answered, and it is not
the same as claiming the term is zero: what is established is a bound, and the bound is below
where anything moves.

The 2026-08-29 report that found the channel and then concluded the sweep had to stay is in
[`superseded/autofluorescence-first-pass.md`](superseded/autofluorescence-first-pass.md). It
is not restated here: its plate reading is right, its median `RFU/OD` is the method this
section replaces, and the 239 RFU/OD it prints is that ratio over a different well set than
the 249 quoted above.

`scripts/measure_autofluorescence.py` -> `outputs/autofluorescence_measured.csv`.

## A resolved identifier is not the same as the right molecule

**2026-08-30.** Running a second product through the chain to check the machinery is not
carotenoid-specific found a defect in the estimator this session had just promoted to primary.

`bridge/equilibrator.py` resolved a metabolite by KEGG identifier and used whatever came back.
yeast-GEM writes a polymer as one representative repeat unit; eQuilibrator holds a specific
oligomer. So `s_0773` glycogen -- **C6H10O5** in the model -- resolved to a **tetrasaccharide,
C24H42O21**, and glycogen synthase came back at **-849.2 kJ/mol** for a glycosyl transfer
textbooks put near -13. A mass-unbalanced reaction, scored without complaint.

**41 of 1,720 resolved metabolites mismatch**, and they are one failure mode rather than
scattered noise: polymers and isoprenoids, the model's repeat unit against a fixed oligomer.
`s_0641` dodecaprenyl diphosphate is C60 in the model and matched a **C10**.

The resolution is now checked against the formula, on heavy atoms only -- protonation state
legitimately differs at pH, and refusing UDP for one hydrogen would throw away most of the
coverage. A mismatch is dropped from the compound table entirely, so a reaction containing one
is **uncovered** rather than mis-scored. It costs 46 of 1,289 single-compartment reactions.

**The check that had to run before trusting anything else this session:** every metabolite of
the carotenoid pathway and the thiolase matches. Those results stand.

## CrtI was written with the wrong electron acceptor, and it cost the model its only oxygen dependence

**2026-08-30.** `fba/carotenoid.py` wrote phytoene desaturase as

    phytoene + 4 FAD + 4 H+ -> lycopene + 4 FADH2

with the comment "electrons to FAD, four protons balance charge" and no citation. It balances
mass and charge. It is also **thermodynamically impossible**: component contribution scores it
at **+166.2 +- 12.8 kJ/mol**, uphill by thirteen standard errors, on a step Verwaal 2007
(PMID 17496128) measured running to completion in this organism.

**The gate had been reporting that as uncertainty in the literature.** All three carotenoid
steps returned `cannot_say`, and the reason recorded in `data/thermo/refuted_energies.tsv` was
that the two estimators straddle zero on CRTI -- "a disagreement about which way the step
goes". It was a disagreement about which reaction we had written down.

**The defect is that FAD in CrtI is a prosthetic group, not a substrate.** A bacterial-type
phytoene desaturase holds one flavin and reoxidises it; four free FADH2 as products makes the
cell pay for four flavin reductions it never performs. Free FAD is also a poor oxidant,
E'0 = -0.219 V, which is the whole of the +166.

**That eQuilibrator is right about flavins was checked, not assumed.** yeast-GEM carries
exactly one cytosolic FAD reaction, soluble fumarate reductase `r_0455`. From
E'0(fumarate/succinate) = +0.031 V and E'0(FAD/FADH2) = -0.219 V the textbook value is
-2 x 96.485 x 0.250 = **-48.2 kJ/mol**:

| estimator | r_0455 | vs textbook |
| --- | ---: | ---: |
| component contribution (eQuilibrator) | **-43.5** | 4.7 |
| ModelSEED group contribution | -12.6 | 35.6 |

**Written as the flavin oxidase it is,** with O2 as terminal acceptor and hydrogen peroxide as
the two-electron product:

| CRTI stoichiometry | dGr'0 kJ/mol |
| --- | ---: |
| `+ 4 free FAD -> + 4 FADH2` (as shipped) | **+166.2 +- 12.8** |
| `+ 4 O2 -> + 4 H2O2` (flavin oxidase; adopted) | **-292.7 +- 17.4** |
| `+ 2 O2 -> + 4 H2O` (4-electron; not used) | -676.5 +- 8.4 |

The four-electron form is not used: flavoproteins reduce O2 by two electrons at a single site.
The choice of acceptor is a modelling decision and is recorded as one; that the shipped form
was impossible is not.

**What it fixed downstream.** The sign disagreement is gone, so the refutation now refuses a
MAGNITUDE and not a verdict, and the whole pathway gates as `runs` when both estimators are
available (`PreferredEnergies.both_dg0`, and the sign check now runs on EVERY reaction). The remaining disagreement is 22 to 124 kJ/mol
on values of 57 to 292 -- and it is **not specific to C40 polyenes**, which is what the
refutation used to claim: `CRTE`, which touches no C40 species at all, disagrees by the
largest ratio of the four at 8.6x.

**And it gave the model back an oxygen dependence it should always have had.** Carotenoid
synthesis in yeast is strictly aerobic. Before this the pathway consumed no O2, so the model
would have predicted beta-carotene in an anaerobic chemostat:

| O2 uptake limit, mmol/gDCW/h | max beta-carotene flux | mg/gDCW/h |
| ---: | ---: | ---: |
| 0 | **0** | **0** |
| 0.1 | 0.0212 | 11.4 |
| 1.0 | 0.0692 | 37.2 |
| >= 2.0 | 0.0863 | 46.4 |

The FBA ceiling itself barely moved, 0.0873 -> 0.0863 mmol/gDCW/h.

**This does not make the prediction environment-sensitive, and the gap is now specific.**
`predict.py` still records the culture context and discards it: twelve environments spanning
three carbon sources, two oxygen fractions and two temperatures return **one** content,
1.106538 mg/gDCW, at D = 0.15 /h. What moves the number is the pump (1.150 -> 0.982 mg/gDCW
across D = 0.11 to 0.25) and the genotype (1.107 -> 1.219 across 1x to 4x crtE). What is new
is that a MECHANISM now exists for oxygen to matter, in the FBA layer, where before there was
none. Wiring it needs a gas fraction converted to an uptake bound, which needs a kLa nobody
here has measured -- so it is named as the missing measurement rather than invented.

## The thiolase lead was a unit error, and the estimator that refutes it was already here

**2026-08-30.** The one mechanism that ever explained Kocharin's 4.24× PHB flux difference
between carbon feeds is refuted. The count of mechanisms that fail to explain it is six.

**The defect.** `bridge/thermodynamic.py` built its formation-energy table by calling
`pytfa`'s `MetaboliteThermo` without passing `thermo_unit`. It defaults to `'kJ/mol'`; the
vendored ModelSEED database is `kcal/mol` and declares so in its own `units` field. The pH
and ionic-strength transform was therefore computed in kJ on top of a tabulated energy in
kcal, and every reported energy was on a scale that was neither.

**What exposed it** is a value from outside both estimators. Transformed to pH 7.0 and
I = 0.25 M, water's standard formation energy is −155.66 kJ/mol (Alberty 2003, table 4.2).
The table read in its own unit and converted once gives −155.70. Read under the default it
gives **+24.85**, and no scale on which the formation energy of water is positive is a scale.

**Three independent checks agree the correction is right, not merely different:**

| check | mixed scale | corrected | reference |
| --- | ---: | ---: | ---: |
| water, ΔG′f at pH 7, I = 0.25 | +24.85 | **−155.70** | −155.66 |
| ATP hydrolysis | −37.32 | **−28.65** | −30.5 |
| vs eQuilibrator, 1,098 shared reactions | RMSE 320.1, r 0.435 | **RMSE 116.6, r 0.942** | — |

**What it costs.** The thiolase `r_0103` goes from +15.6 to **+38.07 kJ/mol**, and its
acetyl-CoA threshold from 221 µM to **19 mM**. Measured cytosolic acetyl-CoA is ~10 µM on
glucose and ~425 µM on ethanol, and the lead was that these straddle the threshold. They do
not. Nothing straddles 19 mM.

**And it was checkable all along, which is the part worth keeping.** eQuilibrator never
touched the unit bug. Its energy for the same reaction is **+24.96 kJ/mol** — the literature
value for a thiolase condensation, ~+26 — and it puts the threshold at **1.4 mM**, so ethanol
misses on that reading too:

| reading | dG′r0 | threshold | glucose 10 µM | ethanol 425 µM |
| --- | ---: | ---: | --- | --- |
| ModelSEED, unit bug (published) | +15.61 | 221 µM | below | **ABOVE** |
| ModelSEED, unit corrected | +38.07 | 19,042 µM | below | below |
| eQuilibrator (matches literature) | +24.96 | 1,414 µM | below | below |

The refutation does not depend on the fix. `bridge/equilibrator.py` has been in this
repository since the ribose-5-phosphate isomerase work, is documented there as the better of
the two estimators, and would have said this whenever it was asked. It was never asked about
the reaction the project's biggest metabolic claim rested on.

**What survives.** The step is uphill and has a threshold; that was never in doubt. The gate
is unaffected — it returned `cannot_run` on the glucose state and now returns it on ethanol
too, and refusing an infeasible step does not require the step to be feasible somewhere. The
wiring, the guards, the three states and the direction check are untouched. `36.2 kJ/mol` for
the naive quadrature is also untouched: the tabulated formation *errors* are in kJ/mol in the
source and were never on the mixed scale.

**Collateral.** The two estimators' disagreement on the carotenoid steps shrank from 91, 275
and 252 kJ/mol to **22, 225 and 124** — part of the old gap was the unit, not the chemistry.
It did not shrink enough to change `refuted_uncertainty`, and `CRTI` still disagrees in sign:
−59 against +166 kJ/mol, which needs no threshold to refuse.

**Pinned by** `tests/test_thermodynamic_physics.py::TestTheTableIsReadInItsOwnUnit` (water,
the declared unit, and the eQuilibrator comparison) and `tests/test_thiolase_threshold.py`,
whose `test_ethanol_also_sits_below_it` asserted the opposite until today.

## The fourth plate was on disk the whole time, and it is worth two calls and a correction

**2026-08-30.** The biosensor panel is **n=4**. Plate `20260804` had been excluded because
every `.xlsx` export of it was blank-subtracted; the instrument's own `.xpt` still carries
the background, `plate/gen5.py` reads it, and it is now committed as text in `data/plates`
beside the subtracted export it was exported to.

**The delta, at the canonical late window, is two calls gained and none lost:**

| view | at n=3 | at n=4 |
| --- | ---: | ---: |
| folds estimable of 24 | 20 | 20 |
| clearing 1.0 unadjusted | **5** | **7** |
| Bonferroni family-wise over the family of 24 | **unreachable** | **unreachable** |
| holding at all six late windows | 4 of 5 | **7 of 7** |
| exact sign-test floor | `2*(1/2)**3` = 0.25 | `2*(1/2)**4` = **0.125** |

The two new calls, pinned to the table that produced them:

| construct | dose | fold | 95% CI at n=4 | 95% CI at n=3 |
| --- | --- | ---: | --- | --- |
| UPRE2 | 0.5 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE2;stressor=DTT;dose_mM=0.5;late_fraction=0.75" -->1.42 | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE2;stressor=DTT;dose_mM=0.5;late_fraction=0.75" -->1.14, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE2;stressor=DTT;dose_mM=0.5;late_fraction=0.75" -->1.64] | [0.99, 1.73] |
| AlteredYap1 | 0.5 mM H2O2 | <!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.5;late_fraction=0.75" -->1.37 | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.5;late_fraction=0.75" -->1.06, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.5;late_fraction=0.75" -->1.70] | [0.90, 1.96] |

The second is **the first per-dose call either oxidative sensor has ever produced.** Every
call that existed at n=3 is still a call at n=4.

**A family-wise claim was made here and then withdrawn the same day.** This section briefly
reported that `UPRE1` at 0.5 mM DTT survived correction for its whole family, at
[1.009, 2.068] — "the first family-wise-corrected biosensor claim in this repository's
history". It was not one. The 2026-08-30 audit found two defects underneath it, and fixing
either alone removes it:

- `scripts/fold_multiplicity.py` set the Bonferroni divisor from `len(table)`, the **20
  estimable** folds, rather than from the **24** construct × dose cells the panel asks
  about. A fold that could not be estimated was making the correction on every other fold
  weaker.
- `analysis/adjusted_interval` answered any level it was given. A cluster bootstrap on four
  plates is supported on plate multisets, so its smallest atom is `(1/4)**4` and the finest
  two-sided level it can express is **0.0078**. The family-wise level is `0.05/24` =
  **0.0021**. The interval was a quantile of a tail the distribution has no mass in.

So the honest row is not "nothing survives the correction" but **the correction cannot be
applied at four plates at all**, and `adjusted_interval` now refuses instead of returning a
number. **Five plates is the first count at which the question can be asked.** That is a
sharper statement than the one it replaces, and it is the second time in this file that
counting the resolution of the design mattered more than the measurement.

**What did not improve.** The two-sided exact-sign test still has an attainable-p floor
above 0.05 at four plates; six plates are needed for that test, not for every possible
inferential method. Bootstrap tail fractions are not permutation p-values. The per-sensor
comparison below separates historical Bonferroni results from the current Holm adjustment:

| construct | stressor | historical Spearman at n=3 | current at n=4 | historical Bonferroni x4 at n=3 | current Holm at n=4 |
| --- | --- | ---: | ---: | ---: | ---: |
| UPRE1 | DTT | +0.946 | <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=UPRE1;measurement=dose_response_permutation" -->+0.947 | 2e-04 | <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=UPRE1;measurement=dose_response_permutation" -->0.0002 |
| UPRE2 | DTT | +0.906 | <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=UPRE2;measurement=dose_response_permutation" -->+0.917 | 2e-04 | <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=UPRE2;measurement=dose_response_permutation" -->0.0002 |
| NativeYap1 | H2O2 | +0.593 | <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=NativeYap1;measurement=dose_response_permutation" -->+0.469 | 2e-04 | <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=NativeYap1;measurement=dose_response_permutation" -->0.0005 |
| AlteredYap1 | H2O2 | +0.607 | <!-- audit:value table=outputs/fold_multiplicity.csv column=statistic row="construct=AlteredYap1;measurement=dose_response_permutation" -->+0.434 | 4e-04 | <!-- audit:value table=outputs/fold_multiplicity.csv column=p_holm row="construct=AlteredYap1;measurement=dose_response_permutation" -->0.0010 |

All four current Holm-adjusted tests reject the declared within-plate permutation null.
The smaller oxidative correlations remain visible, but processing and adjustment changes
prevent attributing the p-value differences solely to adding one plate. In the historical
comparison, `AlteredYap1` at 1.0 mM H2O2 moved from 1.49 to **1.00**, interval [0.24, 2.50].
That historical non-call is retained, not promoted to a current point estimate.

**What the plate cost to recover, and what was checked before using it.** Subtracting the
per-timepoint mean of H1:H3 from the archive reproduces the subtracted export exactly in
fluorescence — 2,175 readings, every disagreement zero or exactly a third, which is what
rounding the mean of three integers back to an integer must produce — and to 0.0005 in
absorbance, the arithmetic maximum for a three-decimal export of a four-decimal reading. That
comparison now runs from the committed text with no `.xpt` present
(`tests/test_xpt_export.py`). Three further things had to be right and one of them was not:

- **The logbook's blank wells for this plate were wrong**, recorded as H4-H6 the way
  20260803's entry was before it was corrected. The archive carries all 96 wells, and in the
  reporter channel H1-H3 read 341 RFU against 56 for H4-H12 — medium fluoresces at 480/530
  and water does not. Gen5 agrees: its own subtraction on this plate is the mean of H1:H3.
  Committing the archive made this load-bearing rather than cosmetic, because the archive
  *has* H4-H6 and the stale entry would have selected them.
- **The dose ladder and the layout are the recorded ones.** The `.xpt` holds plate reads and
  no per-construct plot sheets, so the ladder comes from the logbook — and the subtracted
  export's own sheets agree with it dose for dose, while recovering the layout from the
  archive returns the recorded columns at 3.4% fit.
- **One plate is one biological replicate.** Both copies stay committed, so
  `run_sensor_characterisation.one_file_per_plate` decides which is the replicate, by the
  logbook plate they share, and refuses rather than ranking when neither is an instrument
  file. Without it the panel would have reported n=5 for four cultures.

`scripts/plate_readings.py --write --window`, then `scripts/fold_multiplicity.py`.
`outputs/sensor_characterisation.csv` (335 well rows, was 251),
`outputs/late_window_sensitivity.csv`, `outputs/autofluorescence_sensitivity.csv`,
`outputs/fold_multiplicity.csv`.

## Twenty-four folds, corrected for each other, and none of the five survives

**2026-08-29. Superseded above on the plate count: this entry is the n=3 state.** The
reasoning holds and the counts moved — 5 calls became 7 and the sign-test floor 0.25 became
0.125. **The family-wise column did not move**, and this entry's own sentence below —
"at three plates a family-wise correction is not conservative, it is unreachable" — turned
out to be true at four as well, for a reason it did not yet name: the smallest atom of a
four-plate cluster bootstrap is `(1/4)**4`. It briefly read "0 to 1" on 2026-08-30; see the
withdrawal above.

`outputs/late_window_sensitivity.csv` tests twenty-four folds against 1.0,
twenty of them estimable, each at a nominal 95%. Nothing had divided by twenty.

Doing so is not the interesting part. **At three plates a family-wise correction is not
conservative, it is unreachable.**

| view | calls |
| --- | ---: |
| unadjusted 95% intervals | **5** |
| Bonferroni family-wise over m = 20 | **0** |
| exact sign test over 3 plates | p floor 0.25 — nothing clears at any conventional alpha |

The middle row is arithmetic, not evidence: the interval is built from a t quantile with two
degrees of freedom, and `t.ppf(0.975, 2) = 4.30` against `t.ppf(1 - 0.05/40, 2) = 17.3`, so
every width more than triples for reasons unrelated to the measurements. The bottom row says
why. With three plates all on one side, the smallest two-sided p an exact distribution-free
test can return is `2 * (1/2)**3 = 0.25`. **Six plates is the first count that reaches 0.05
on a single fold; ten is the first that survives correcting for twenty.**

**The bootstrap tail probabilities are therefore not p-values.** They read 0.0003 for six
different folds -- a thousandfold below what the design supports. A cluster bootstrap on
three clusters enumerates draws from a three-point empirical distribution: it can say the
three plates agree and it cannot say a fourth would. `analysis/uncertainty.py` already knew
this about its *intervals*, which is why it t-inflates them and refuses below three plates.
The same caution never reached the tail. Holm and BH are computed and reported and change
nothing, because they are corrections applied to numbers that cannot be believed.

### What survives: one test per sensor instead of six

"Does this reporter respond to dose at all" is one test, and it is not bounded by the plate
count, because permuting dose labels **within** each plate holds the plate's inoculum, medium
batch, reader gain and incubator position fixed rather than estimating them from three.

| construct | stressor | plate-averaged Spearman | permutation p | Bonferroni x4 |
| --- | --- | ---: | ---: | ---: |
| UPRE1 | DTT | +0.946 | 5e-05 | 2e-04 |
| UPRE2 | DTT | +0.906 | 5e-05 | 2e-04 |
| AlteredYap1 | H2O2 | +0.607 | 1e-04 | 4e-04 |
| NativeYap1 | H2O2 | +0.593 | 5e-05 | 2e-04 |

*These are the n=3 values and `outputs/fold_multiplicity.csv` now holds the n=4 ones; the
2026-08-30 entry above has both columns side by side.*

All four, comfortably, **including both oxidative sensors, which produced zero per-dose
calls.** So the biosensor result is better than the per-dose table said in one direction and
worse in the other: the evidence that these reporters track dose is strong, and the evidence
that any particular dose gives any particular fold is weak. Only the second claim was being
made. `scripts/fold_multiplicity.py`, `outputs/fold_multiplicity.csv`.

What it does not license is a number. "Activity tracks dose" is a within-plate statement; the
size of a response at one dose is exactly the between-plate quantity three plates cannot pin.

## The crosstalk table is one replicate, and its own SOURCE.md said three

**2026-08-29.** `data/crosstalk/erox_2026-08_endpoint.tsv` is described by
`data/crosstalk/SOURCE.md` as "the workbook's own biological-replicate mean across the three
plates". It is **replicate 1 alone**, and the `std/2` columns are the technical triplicate
spread within that single plate. The `Summary` sheet of `ER&OX-summary.xlsx` lays the three
replicates side by side at columns B, O and AB; only the first was transcribed.

Checked cell by cell. UPRE1 at 0 mM reads `0.320333 / 1774.0` in the TSV and in block 1,
against `0.517333 / 2691` in block 2 and `0.487667 / 2373.67` in block 3.

Two consequences that are not cosmetic:

- **The numbers move.** The workbook's own n=3 aggregate gives UPRE2 + H2O2 at 1 mM as
  0.840 ± 0.113 against the 0.982 carried here; NativeYap1 + DTT at 2 mM as 1.158 ± 0.130
  against 1.084.
- **The viability gate flips.** `score_sensor_crosstalk.py` reads UPRE1 as 90% viable at
  1 mM H2O2 on this plate; on the other two it is **40% and 31%**, already dying at a dose the
  script treats as readable. Replicate 1 is also the plate with the lowest control OD.

The verdict's *direction* is probably unaffected — the workbook's own n=3 off-target folds are
still near 1.0 — but this repository refuses an interval below three biological replicates,
and the crosstalk claim was resting on one. **Re-deriving the TSV from all three blocks is the
fix**, and it has to carry the workbook's own caveat, verbatim from `Summary!B2`: *"Note that
for the second replicate, the original values are not blanked."*


## The estimator's smoothing window was chosen by a rule tuned on a different experiment

**2026-08-29.** `reporter.py::default_activity_window_h` sets the window every
dilution-corrected number here depends on. It justified that with a five-cell table in its
own docstring -- `0.5 h -> 14.1% ... 12 h -> 1.6%` -- that **no script produced and no test
read**, measured over a 24-hour run at 145 points. The committed exports are 25 points over
4.00 hours, where the rule returned 0.667 h: below the smallest cell in the table justifying
it. The only test of the rule checked it at 12 h and 48 h, six and twelve times the plates'
density. It was tested only where it works.

`analysis/estimator_accuracy.py` is the producer, and every parameter in it is read off the
plates rather than invented -- 25 points over 4.00 h, 1.54 doublings, mu falling 0.563 to
0.217 /h, 1.1% multiplicative reader noise taken as the residual of each channel about a
smooth fit in log space. `outputs/estimator_accuracy.csv` is the table.

**Three findings, and the second is the one that matters.**

**1. The error is set by the window's duration in hours, not by its share of the run and not
by its point count.** Two geometries six times apart in density agree at matched duration and
differ fourfold at matched point count:

| window | 145 pts / 24 h | 25 pts / 4 h |
| --- | ---: | ---: |
| 0.33 h | 9.8% | 8.2% |
| 1.17 h | 5.4% | 4.5% |
| 2.17 h | 2.3% | 2.0% |
| **3 points** | **9.8%** (145 pts / 24 h) | **2.1%** (25 pts / 24 h) |

A rule expressed as a fraction of the run is therefore correct only at the run length it was
tuned on. Fixed by adding an absolute floor, `TARGET_ACTIVITY_WINDOW_H = 2.0`, which the
existing third-of-the-run cap then binds to 1.333 h on these plates. Every longer geometry is
unchanged: a 24 h run's sixth is already 4 h and outranks it.

**2. What the window buys is precision, not accuracy -- so the point estimates did not move
and one call did.** A fold is a ratio of two wells on one plate, and a smoother treats
numerator and denominator alike, so the recovered fold is under 0.6% low at every window from
0.667 h to 2.0 h. What changes is the spread: sd 0.116 at 0.667 h against 0.067 at 2.0 h.
Re-deriving the sweep moved twenty-one folds by at most a few thousandths and took the calls
at the default window from **six to five**. UPRE2 at 0.5 mM DTT had been clearing 1.0 with a
lower bound of 1.014 and now straddles it. A result that survives only at one setting of an
arbitrary knob was never worth the sixth place it held.

**3. Ignoring chromophore maturation costs ~10% on a single well and ~2% on a fold.** The
pipeline runs `ReporterKinetics(k_deg=0.0)` against a YFP that takes tens of minutes to
mature, and `promoter_activity_from_total` refuses maturation outright -- there is no code
path by which a correction could reach a committed number. Generated with a 20-minute
maturation and inverted without one, a single well's late activity is biased +10.8%; the same
bias on a fold is +2.5%, because a ratio cancels what both wells share. Costed rather than
assumed, in either direction.

### And a defect the re-derivation exposed

`analysis/uncertainty.py::fold_change` could return `estimable=True` with a **NaN point
estimate and a real-looking interval**. `_combine` takes a geometric mean over positive
per-plate folds, so a condition where every plate's recovered activity is negative returns
NaN -- and the bootstrap guard only requires half the *resamples* to be positive, so it
passed. `excludes_unity` reads `low > 1.0 or high < 1.0` and knows nothing of the point, so
AlteredYap1 at 2.0 mM H2O2 -- naive fold **-0.22**, a culture losing signal because the
peroxide is killing it -- came back as `[0.003, 0.585]`, entirely below 1.0, which reads as a
confident finding of strong repression. Now refused, with the reason named. The invariant is
tested over the whole shipped table: an estimable fold must have a fold.


## Claim boundary — read this before any table

Every number below is one of four kinds of thing, and the difference between the kinds is
larger than the difference between any two numbers within one kind. Tiers, with their
precedents, are defined in `docs/research/CIRCULARITY.md` §5; the per-claim assignment is
`docs/CLAIM_BOUNDARY.md`.

| Tier | Name | What it means |
| --- | --- | --- |
| **0** | **Asserted** | A constant chosen for plausibility. Nothing measures it. |
| **1** | **Self-consistent in simulation** | The fitter correctly inverts the generator it was given, or a deterministic computation was carried out correctly on a model we adopted. Necessary, and much weaker than it looks. |
| **2** | **Robust across generator families** | Holds when the generator's *structure* is perturbed — or when it holds on a measurement the model never saw. |
| **3** | **Predicted then measured** | Registered *before* the measurement, scored against it. |

Where this repository actually stands, and it is uncomfortable:

- **Nothing is Tier 3.** One forward prediction has now been scored on held-out real data,
  and **the twin lost** — beaten by a six-point log-linear extrapolation on OD by 23% of
  level at a 2 h horizon, with nominal 95% intervals covering 40–49% of the time. The loss
  names its own cause: `estimator.py::_propagate` gives growth rate no drift, so a
  decelerating culture is extrapolated flat. It was scored after the fact rather than
  registered before it, so it is not Tier 3. `docs/research/PREDICTION.md` §0.
- **The historical transfer, power and identifiability tables are Tier 1.** They evaluate
  a declared simulator, not biological transfer. Optimization bias is an expectation-level
  result for selected finite-sample objectives, not a sign guarantee for each metric
  (`docs/research/CIRCULARITY.md` §1). The later `docs/CROSS_FAMILY.md` records structural
  checks and a configuration-choice UCBOG; these do not retrospectively validate every table.
- **G1, D2, G4, the dose-response correction and the calibration work are measured on real
  plates — and they are diagnostic, not predictive.** They say what a measurement is worth,
  which is a different and more modest job than saying what the next one will be.
- **Besides that loss, Tier 2 evidence exists in exactly two places, and one of them is a
  failure**: the Vos 2016 retentostat prediction holds on a genuinely held-out measurement,
  and the van Hoek 1998 critical growth rate does not. Both are on the physiology side.

Three corrections that this labelling pass forced into the prose below, flagged here because
they change what earlier versions of this file said:

1. **The dataset is not 4 biological replicates for fluorescence.** Replicates 2 and 4 read
   columns 1–3 only, which is UPRE1. UPRE2, NativeYap1 and AlteredYap1 are at **n = 2**, so
   **all eight dose-response readings are INCONCLUSIVE**, and `AlteredYap1 @ 1.0 mM`
   corrects from 1.43 to **1.31** under a plate-paired estimator. `docs/DATA_INVENTORY.md`.
2. **The ER anchor did not fail to conclude; it failed as a measurement.** UPRE1's
   +0.13-cycle RT− margin is **91.4% genomic DNA**, far past the 60% ceiling of the only
   validated correction. Nothing was learned about the ER biosensors.
   `docs/research/G4_STATISTICS.md`.
3. **"One channel buys one module" is a known-loadings result and does not license the
   transfer work.** The transfer and training results go through a *fitted*-loadings path,
   where the Ledermann bound gives four channels **one** factor, not four.
   `docs/research/IDENTIFIABILITY.md`.

Where a section below reports a table, it carries a tier badge. Where it does not, assume
Tier 1 and check `docs/CLAIM_BOUNDARY.md`.

## How this file is arranged

Design work first (what the experiment has to be), then what the real plates measured, then
the audits that corrected earlier numbers, then the metabolic and thermodynamic layers.

| | |
| --- | --- |
| **Design, simulated** | [Simulation-first](#simulation-first) · [What the design has to be](#what-the-design-has-to-be-measured-not-asserted) · [What is new and what is not](#what-is-new-here-and-what-is-not) · [Maturation](#maturation-is-information-not-a-defect) · [Co-expression](#co-expression-changes-the-answer) · [Fisher information](#design-on-fisher-information-not-on-a-spanning-score) · [Transfer](#transfer-train-on-some-stressors-predict-one-never-seen) · [The trained model](#the-trained-model) · [Contexts](#training-across-culture-conditions) · [The time axis](#the-time-axis) |
| **Measured, real plates** | [Scope](#scope-two-separate-things-live-here) · [Findings that constrain the plan](#findings-that-constrain-the-plan) · [Calibration](#calibration-what-the-uploaded-plates-did-to-the-panel) · [Optical linear range](#measuring-the-optical-linear-range) · [Physiology validation](#validating-against-published-physiology) · [BioNumbers](#kaggle-and-what-it-was-and-was-not-good-for) |
| **Corrections** | [Audit](#audit-where-the-earlier-numbers-were-wrong) · [Second audit](#second-audit-three-structural-gaps-closed) · [Third audit](#third-audit-labelling-what-each-number-is) · [Where the numbers stand](#where-the-numbers-stand-after-both-audits) · [Removed](#removed) |
| **Metabolism** | [The metabolic link](#the-metabolic-link-and-what-the-literature-would-not-support) · [TMFA](#thermodynamic-flux-analysis) · [E_GSH](#the-e_gsh-disagreement-resolved) · [Thermodynamic data](#better-thermodynamic-data) |
| **Deeper detail** | `docs/CONTRACT.md` · `docs/CLAIM_BOUNDARY.md` · `docs/DATA_INVENTORY.md` · `docs/NULL_RESULTS.md` · `docs/ARCHITECTURE.md` · `docs/REPRODUCING.md` · `docs/PROTOCOLS.md` · `docs/PARKED.md` · `docs/ATP_SENSOR.md` · `docs/research/` |

## Simulation-first

The generator exists so the pipeline can be trained and stress-tested before the wet lab
catches up. It is mechanistic, not fitted to curves: three equations produce a plate, and
the plate is written as a Synergy export that re-enters through the same readers and gates
as a real file. Uncalibrated it lands at OD 0.096-0.886 against the real 0.094-0.800 and
passes G1 at 80%, matching the 2026-07-22 plate.

**The decoupling grid is the point.** Every perturbation run so far slows growth, so stress
and `1/mu` are one axis and no latent stress state is separable from growth at any sample
size. `generator/design.py` crosses the dose ladder with nutrient limitation, which slows
growth *without* stressing, and `collinearity()` scores a planned design before it is run:
0.993 for the dose-only series already performed, 0.694 for the crossed grid. On simulated
data where truth is known, only the decoupled design recovers the true coefficients
(+1, -1) stably; the collinear one is unbiased but three times noisier.

## What the design has to be, measured not asserted

Two conditional simulated analyses compare designs; neither establishes biological state identifiability.

**Detecting a dose slope is easy; attributing it to dose is not.** At six simulated biological
replicates, with growth rate also in the model:

**Tier 1 — self-consistent in simulation.** These current estimates use the declared generator, seed and 120 Monte Carlo campaigns, not independent biological validation (`docs/research/CIRCULARITY.md` §1). The registered table retains Wilson uncertainty; 100% observed success is not certainty about true power.

| Design | collinearity | 1.2x | 1.5x | 2.0x | 3.0x |
| --- | --- | --- | --- | --- | --- |
| dose only (as run) | <!-- audit:value table=outputs/power_analysis.csv column=collinearity row="design=dose only (as run);question=attribute;fold=1.2" -->0.993 | <!-- audit:value table=outputs/power_analysis.csv column=power_n6 row="design=dose only (as run);question=attribute;fold=1.2" scale=100 -->**7.5%** | <!-- audit:value table=outputs/power_analysis.csv column=power_n6 row="design=dose only (as run);question=attribute;fold=1.5" scale=100 -->**16.7%** | <!-- audit:value table=outputs/power_analysis.csv column=power_n6 row="design=dose only (as run);question=attribute;fold=2.0" scale=100 -->52.5% | <!-- audit:value table=outputs/power_analysis.csv column=power_n6 row="design=dose only (as run);question=attribute;fold=3.0" scale=100 -->96.7% |
| dose x 2 nutrients | 0.587 | 100% | 100% | 100% | 100% |
| dose x 3 nutrients | <!-- audit:value table=outputs/power_analysis.csv column=collinearity row="design=dose x 3 nutrients;question=attribute;fold=1.2" -->0.507 | <!-- audit:value table=outputs/power_analysis.csv column=power_n6 row="design=dose x 3 nutrients;question=attribute;fold=1.2" scale=100 -->100% | 100% | 100% | 100% |

The within-plate variance component that table is run at is measured, not guessed. From the
zero-dose wells of the two usable plates, three wells per construct per plate differing only
in position, the lognormal sigma of recovered activity has a median of <!-- audit:value table=outputs/power_analysis.csv column=well_cv row="design=dose only (as run);question=attribute;fold=1.5" -->**0.052** across
eight construct-by-plate groups. Kensy et al. 2009 report under 5% between replicate wells of
one clone, which agrees with it. `outputs/power_analysis.csv` carries Wilson bounds per cell
and records the `well_cv` it was run at, so a change to that constant shows up in the
artefact rather than only in the code.

That constant reached its current value through two earlier states — no within-plate term at
all, then a guessed 0.10 that was unreferenced and backwards — and the full record of both,
with why the alarm about the second was itself wrong, is in
[`superseded/power-attribution.md`](superseded/power-attribution.md). It is not restated
here: this document carries the current row, the archive carries the versions.

The measured corrected effects are 1.2-1.5x. **At 1.5x the current design has <!-- audit:value table=outputs/power_analysis.csv column=power_n6 row="design=dose only (as run);question=attribute;fold=1.5" scale=100 -->16.7% power to
attribute the response to dose rather than to the growth the dose caused; adding one
nutrient level takes it to 100%.** Detection alone needs 3 replicates at 1.5x and above and
5 at 1.2x, which is why detection is the wrong question to power for.

A chemostat is *not* required. Estimating mu from a plate OD trace adds negligible error
here -- the growth coefficient recovers to -0.998 against a true -1.0, with no attenuation
-- so the nutrient axis is what matters, not the vessel.

**Stressors must be crossed, not run one at a time.** With one stressor on at a time the
other is off, so the two axes are anti-correlated by construction and a single latent
factor explains 70% of every channel without any shared biology:

**Tier 1 — self-consistent in simulation.** Scored within the declared generator configuration, not validated against biology (`docs/research/CIRCULARITY.md` §1). No null or UCBOG accompanied this historical table.

| Design | states found | variance one state explains |
| --- | --- | --- |
| one stressor at a time | 2 | 70% (spurious) |
| crossed stressors | 2 | 8% |
| crossed, with a real shared burden | 1 | 98% |

Only the crossed design can tell a genuine general stress axis from an artefact of the
dosing schedule. `analysis/latent.py` selects its dimension by held-out reconstruction
error, not by counting sensors, and an in-sample criterion picks the wrong dimension here
because a fit with as many factors as channels is exact.

## The whole stress landscape, and what five channels can reach

> **Superseded in two places, both below.** The panel is no longer seven regulons — see
> [The whole stress landscape, and two kinds of sensor](#the-whole-stress-landscape-and-two-kinds-of-sensor)
> — and the claim that ESR, osmotic and proteasome are *intrinsically* out of reach was the
> old 0.35 projection threshold talking, not the algebra; see
> [Design on Fisher information](#design-on-fisher-information-not-on-a-spanning-score).
> The rank statement itself survives both. Kept because the reasoning about ESR's ubiquity
> is what motivated the dedicated channel, and that conclusion did not change.

`generator/stress_panel.py` carries seven yeast stress regulons with their literature
crosstalk -- ESR (Msn2/4), UPR (Hac1), oxidative (Yap1/Skn7), heat (Hsf1), osmotic (Hog1
via Msn2), proteasome (Rpn4, driven by all three of Msn2/4, Yap1 and Hsf1) and iron
(Aft1/2, the only one with no shared driver) -- plus eight stressors that hit them in
known patterns. Nothing here is limited to ER and oxidative.

`analysis/sensor_selection.py` chooses which reporters to build under a channel cap, by
how much of the module space a set spans rather than by how many sensors it contains.

**Tier 1 — known-loadings path only.** This table is a rank statement about a *given*
loading matrix, and it does not license the transfer or training results further down.

| Stress channels | Modules identifiable |
| --- | --- |
| 2 | 2 of 7 |
| 4 | 4 of 7 |
| 7 | 7 of 7 |

One channel buys one module -- crosstalk gives nothing away free. With five fluorophores
and one spent on the growth reference, four stress channels reach **UPR, oxidative, heat
and iron**, leaving ESR, osmotic and proteasome out of reach.

### Two questions answer to the word "identifiable", and they have different ceilings

The distinction is from `docs/research/IDENTIFIABILITY.md` §1–2, and it matters because this
file previously reported both as one.

**Path A — loadings known.** `analysis/design.py` builds `L' S^-1 L` from
`reporter_loadings`, and `analysis/sensor_selection.py` scores candidate sets the same way.
`L` is *given*, transcribed from the literature with PMIDs, so this is a linear inverse
problem with a known operator and identifiability is the rank of `L`. **The table above is
correct within this framing**, and so is the stronger `N channels resolve N of 24` further
down.

**Path B — loadings fitted.** `analysis/latent.py` runs an SVD and *estimates* the loadings.
That is exploratory factor analysis, bound by Ledermann (1937): `m <= (2p + 1 - sqrt(8p+1))/2`.

| observed channels *p* | 2 | 3 | **4** | **5** | 6 | 7 | 8 | 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| max integer fitted factors | 0 | 1 | **1** | **2** | 3 | 3 | 4 | 6 |

`scripts/run_training.py` and `scripts/run_transfer.py` — every transfer result and every
"general stress state" result in this file — go through **Path B**. The table above comes
from **Path A**. *No line of code connects them.* Four fitted channels support **one**
factor, not four.

Two consequences, stated rather than implied:

- The latent-width table further down selects width 4 for the four-channel build and remarks
  that "nothing is compressed away when width equals channels". That configuration is the
  **full-rank degeneracy** — a rank-*p* model of a *p*-column matrix reproduces any masked
  entry, so the held-out error is not held out. `select_dimension` lands on it 1–2 times in
  30 at *p* = 4–5, and the Ledermann ceiling is enforced nowhere in the code.
- The recoverable ground is to stop doing factor analysis on Path B at all, because `L` is
  already known: solve `y = Lx + eps` under non-negativity and sparsity, which is a real
  identifiability constraint rather than a regulariser. Its catch is specific and serious —
  **sparsity fails for ESR**, because a module that is on almost always is not sparse. Which
  is the same conclusion the paragraph below reaches by a different route.

**The general stress state is the hard one, and the fix is nearly free.** ESR is read by
six of the seven reporters as crosstalk and drives five of them. That ubiquity is what
makes it unidentifiable: no reporter is uniquely ESR unless one is dedicated to it, so a
pure spanning criterion skips it. Forcing `STRE-general` into the set costs **0.4%** of the
score and makes ESR identifiable. Build it.

## What is new here, and what is not

A latent stress state read off a yeast reporter panel is **not** new, and earlier versions
of this file left that unsaid — which reads as a claim to have invented it. Four papers
bound the claim, and three of them bound it downward.

| Prior art | What it already did | What it leaves |
| --- | --- | --- |
| Keren 2013, *Mol Syst Biol*, **PMID 24169404** — "Promoters maintain their relative activity levels under different growth conditions" | Measured ~900 *S. cerevisiae* promoters with fluorescent reporters. **60–90% of promoters** change between conditions by one constant **global scaling factor that depends only on the condition, not on the promoter's identity**; "only several scaling factors suffice". | A ceiling, not a gap. If most of what a promoter panel does across conditions is one scale factor, the rank a 24-module panel can reach is near 1, not 24. |
| Granados 2018, *PNAS*, **PMID 29784812** | Decoded **ten** yeast fluorescent TF reporters — nuclear-translocation dynamics, per single cell, by mutual information — into "a precise, distributed internal representation of extracellular change", and derived the panel-design rule with it: generalists (Msn2/4, Dot6, Tod6, Maf1, Sfp1) encode *which* of several stresses but only when stress is high; specialists (Hog1, Yap1, Mig1/2) encode one stress, faster and over a wider range of magnitudes. | Panel → state is done, in this organism, eight years ago. It is not coupled to metabolism, and it reads TF *localisation* rather than promoter fusions. |
| Hansen & O'Shea 2015, *eLife*, **PMID 25985085** | Put an information ceiling on the readout. A single Msn2 target gene carries **1.2–1.3 bits**, two carry 1.67–1.83, and the conclusion the authors draw is that transduction is "limited to error-free transduction of signal **identity**, but **not** signal intensity information". | The division that matters: the state is well founded as a classifier and poorly founded as a magnitude estimator. |
| Kim, Choi & Kim 2026, *iScience*, **PMID 42502377** — "Causal AI digital twin for bioprocess bottleneck diagnosis via metabolic flexibility and rigidification maps" | A digital twin that labels cause-resolved regimes and diagnoses bottlenecks inside a GEM — three quarters of the join this project is attempting. | Its regime labels come from the **model's own active-constraint set** over a Latin-hypercube design space, not from any measurement, and the organism is *Stenotrophomonas maltophilia*, not yeast. It is reported to cite COSMIC-dFBA as its motivation; the article was not obtainable here, so that sentence is **unverified** (`docs/research/HYBRID_METHODS.md` §2). |

**The honest restatement, and it is narrower than what stood here before.** A latent stress
state for yeast is **well founded as a classifier of which stress and poorly founded as a
magnitude estimator**. What is new in this repository is not the latent state. It is two
other things:

1. **The coupling.** Reporter panel → identified state → time-varying constraints on a
   yeast GEM. PROM (Chandrasekaran & Price 2010, PMID 20876091) did regulatory-state-to-bounds
   in *E. coli* and *M. tuberculosis* from transcriptomes, sixteen years ago; Granados did
   panel-to-state in yeast without metabolism. Nobody has joined them in yeast.
2. **The identifiability argument.** A bound on how many modules a panel of *n* fusions can
   resolve. NCA (Liao 2003, PMID 14673099) owns the general criteria for latent regulatory
   states and Keren supplies the rank result that makes the yeast case hard; the two have
   not been combined. This is the [Path A / Path B](#two-questions-answer-to-the-word-identifiable-and-they-have-different-ceilings)
   distinction above, and it is cheap — arithmetic on published numbers plus this
   repository's own plates.

The full prior-art sweep, twenty-odd hits ranked by distance, is
`docs/research/STRESS_RULES.md` §2. Any write-up must carry Granados, Keren, Slavov
(PMID 22456505) and Hansen up front, and must show the inferred state is not growth rate in
disguise.

**The repository's own nulls agree with the prior art, and that agreement is the strongest
part of the argument.** The fitted basis fails its rotated-subspace null twice — in
simulation (+0.0499 observed against a null of +0.0362, p = 0.268,
[`NULL_RESULTS.md`](NULL_RESULTS.md)) and again on Gasch arrays nobody here measured
(+0.688 against +0.695, p = 0.990, skill **−0.025**,
[`EXTERNAL_VALIDATION.md`](EXTERNAL_VALIDATION.md) §6). The measured reason has the shape
Keren's result predicts: **one component carries 51% of the channel variance, and that
component is the general stress response**, which every one of the eight stressors drives. A random rotation
of correlated channels retains that axis and predicts the regulon means just as well, so
the identity of the fitted axes contributes nothing. Two independent lines — one published
in 2013, one measured here in 2026 on data the model never saw — reach the same conclusion
by different routes.

What passes its null is the **known-loadings** readout, `x = y (L⁺)′`: median Spearman
0.465 against a shuffled-loadings null of 0.091, p = 0.005. That is the Path A / Path B
split restated in nulls rather than in algebra, and it is the reason the recoverable route
is to stop fitting loadings that are already known.

## Maturation is information, not a defect

An earlier version of this README said to keep mOrange2 out of a ratio. That was wrong.

Steady state with a maturation step is `R = k*m / [(mu+m+kdeg)(mu+kdeg)]`. A *matched*
pair cancels mu, so its ratio reads relative activity. A *mismatched* pair does not: its
ratio is monotonic in mu, from 1 to `m_A/m_B`. So a deliberately mismatched pair is a
growth-rate sensor, and the wider the mismatch the more sensitive it is:

**Tier 1 — arithmetic, given the model.** An exact consequence of `dR/dt = k - (mu + k_deg) R`, not a measurement.

| Reference pair | Growth sensitivity at mu = 0.25 |
| --- | --- |
| **YFP / mOrange2** | **0.295** |
| eGFP / mOrange2 | 0.274 |
| mApple / mOrange2 | 0.235 |
| CFP / eGFP | 0.000 |

mOrange2 is the slowest-maturing fluorophore in the panel and therefore the *only* one
that makes a usable growth sensor. Every pair without it that was tested scores lower, and
matched pairs score exactly zero. It is not a liability; it is the reference channel.

There are two routes to growth, with different calibration burdens.
`growth_from_reference` needs one channel's absolute activity and works whatever its
maturation. `growth_from_pair` needs only a relative calibration but does require the
mismatch. Once growth is known, `solve_activities` divides each channel's maturation out
exactly -- so maturation is calibrated, never avoided.

## Co-expression changes the answer

Reporters in one cell share `mu` exactly, so a ratio between two of them cancels the
dilution term outright -- structurally, not as a correction, and with no growth estimate
needed:

    R_stress / R_reference  ->  k_stress / k_reference

Simulated, the ratio is invariant to eight decimal places across mu from 0.40 to 0.05.
That rescues the design already run. Attribution power at six biological replicates:

**Tier 1 — self-consistent in simulation.** Scored within the declared generator configuration, not validated against biology (`docs/research/CIRCULARITY.md` §1). No null or UCBOG accompanied this historical table.

| Effect | 1 reporter, dose only | 1 reporter + nutrient axis | + co-expressed reference |
| --- | --- | --- | --- |
| 1.2x | 0% | 99% | **100%** |
| 1.5x | 3% | 100% | **100%** |
| 2.0x | 44% | 100% | **100%** |

**Growth must then be left out of the model.** It has no true effect on the ratio, and on a
collinear design putting it back inflates the dose standard error and throws away the
power the ratio just bought -- 12% instead of 100% at 1.5x. Verify the cancellation once
with `ratio_growth_invariance`, which does need growth to vary independently, rather than
carrying the term in every analysis.

**The fluorophore pairing matters, and one of the planned five breaks it.** Maturation
must match or `mu` does not cancel:

**Tier 1 — arithmetic, given the model.** An exact consequence of `dR/dt = k - (mu + k_deg) R`, not a measurement.

| Pair | Departure from invariance |
| --- | --- |
| YFP / YFP | 0.0000 |
| eGFP / YFP | 0.017 |
| mApple / YFP | 0.048 |
| **mOrange2 / YFP** | **0.226** |

mOrange2 matures in ~2.3 h against ~0.3 h for YFP, leaving a **23% growth-dependent
residual** -- larger than the 1.2-1.5x effects being measured. Pair YFP with eGFP or
mApple; keep mOrange2 out of a ratio unless its maturation is calibrated and corrected.

Two further conditions: the reference must be genuinely constitutive (a stress-responsive
reference divides the signal out, and power falls below 30%), and spectral unmixing error
adds directly to the ratio noise.

## Scope: two separate things live here

**Measured, against real plates.** The wet lab is characterising biosensors — four
mCitrine reporters, UPRE1/UPRE2 for ER stress (challenged with DTT) and
NativeYap1/AlteredYap1 for oxidative stress (challenged with H2O2). **No carotenoid
pathway is in these strains.** G1, D2 and G4 bear on those sensors as they stand.

**Designed ahead of the strain.** D1, `kinetic/carotenoid.py` and the product
inner-filter correction concern the product layer. No product data exists yet, so
nothing in the current data constrains them and nothing in them constrains the
current data. They exist so the product layer is settled before the pathway lands,
not after.

```bash
python3 -m pip install -e ".[dev]"
python3 -m pytest tests/ -q                        # the full suite
python3 -m pytest tests/ -q -m "not integration"   # skips the tests that touch real
                                                   # data or the GSMM
python3 -m pytest tests/ -q --cov=ystwin           # line coverage
python3 scripts/run_gates.py                       # G1 then D2 on the real plates
python3 scripts/parked/run_d1.py                   # FVA width on Yeast9 + crt pathway
python3 scripts/audit_claims.py                    # does this file describe the code that exists?
```

The counts these commands print are deliberately **not** written here. The one audited test
count lives behind the `audit:test_count` marker in the README, where
`scripts/audit_claims.py` re-collects it and fails on drift; three unaudited copies of it in
this block had already gone stale by several hundred.

Run the last one before trusting any of the rest. It checks the prose against the tree —
module names, script paths, test counts, pinned dependencies — and writes
`outputs/audit_claims.csv`. It exists because this file put the test count at a hundred and
forty-nine while the suite had grown past twelve hundred, and documented a module that had
been deleted, and both survived because nothing checked. A reader who spot-checks one number
and finds it wrong has no reason to trust the next one. `docs/REPRODUCING.md` has the full
run order.

## What is here

| Module | Purpose |
| --- | --- |
| `plate/synergy.py` | Synergy H1 kinetic export reader. Finds channel blocks by content, trims protocol padding, keeps repeat reads of one fluorophore distinct, aligns channels onto one clock. |
| `growth.py` | Specific growth rate from `d(ln OD)/dt` via Savitzky-Golay. |
| `reporter.py` | Reporter as a **dynamic state**: `dR/dt = k_synth - (mu + k_deg) R`, optional maturation. `promoter_activity` inverts dilution and maturation. |
| `gates/g1_optical.py` | G1: is an optical channel quantitative for this well? Linear range, blank separation, sustained decline, dynamic range. |
| `diagnostics/dilution_confound.py` | D2: how much of a reporter's dynamic range is growth rate rather than promoter activity. |
| `fba/carotenoid.py` | crtE / crtYB(PSY) / crtI / crtYB(LCY) installed into Yeast9 or GECKO, mass balanced, named by the real construct genes. |
| `fba/physiology.py` | Validates a GSMM against a measured phenotype before any bound it produces is used. |
| `fba/fva.py` | D1: feasible product-flux range and effective-capacity sweep. |
| `fba/surrogate.py` | Interpolated stand-in for repeated LP solves: 161x faster per call, 148,000x batched, median error 0.000%. Refuses to extrapolate. |
| `kinetic/carotenoid.py` | GGPP -> phytoene -> lycopene -> beta-carotene with a shared bifunctional crtYB pool and growth dilution. Intermediate ratios name the limiting enzyme. |
| `calib/od.py` | Dilution-series fit for the reader's saturation; yields `od_linear_max` for G1. |
| `calib/kdeg.py` | Chase fit, photobleaching/degradation split, and the free lower bound on loss rate. |
| `observation.py` | State -> instrument reading. Includes the beta-carotene inner-filter effect, which is inert until a carotenoid-producing strain exists and refuses to correct without a measured coefficient. |
| `estimator.py` | Particle filter over biomass / reporter / promoter activity / growth rate. Missing channels are skipped, never zero-filled. |
| `archive/code/module_admission.py` | Refuses a latent module whose anchor flux the model cannot reproduce. |
| `gates/g4_anchor.py` | G4: PASS / REFUTED / INCONCLUSIVE against an independent anchor. A failing anchor QC forces INCONCLUSIVE, and an unmeasured growth rate is reported as untested rather than as cleared. |
| `qpcr.py` | ddCq relative expression, no-RT contamination QC, and a reader for the annotated Cq layout. |
| `plate/dose_response.py` | Per-construct dose-response sheets: the reporter signal as the team reads it. |
| `plate/layout.py` | Recovers the plate map by matching the derived sheets, and reports the margin over the runner-up so a weak recovery is visible. |
| `generator/culture.py` | S0 mechanistic culture: growth inhibition, promoter induction, reporter dilution. `nutrient_factor` slows growth without stressing. |
| `generator/plate.py` / `export.py` | Whole synthetic plates written as Synergy exports that round-trip through the real readers. |
| `generator/calibrate.py` | S1: fits an uploaded plate; reports gain x gdcw_per_od as unidentifiable. |
| `generator/panel_gates.py` | Physical, posterior-predictive and diversity gates on the path every result runs through. `scripts/run_training.py` refuses to train on data that fails them. |
| `generator/design.py` | The decoupling grid, and `collinearity` as an explicit diagnostic on a planned design. |
| `generator/literature.py` | Grounded parameters with sources, and published expectations kept separate from them. |
| `analysis/power.py` | Detection vs attribution power; replicates needed per design; closed-form reporter so a sweep runs in seconds. |
| `analysis/latent.py` | Multi-channel latent state, dimension by held-out error, missing readings masked not zero-filled. |
| `generator/unmixing.py` | Spectral mixing and least-squares unmixing for the five planned fluorophores. |
| `generator/stress_panel.py` | Seven yeast stress regulons with literature crosstalk, eight stressors, seven reporters. |
| `archive/code/deconvolve.py` | Joint N-channel solve: growth from a reference or a mismatched pair, then every activity with maturation divided out. |
| `analysis/sensor_selection.py` | Which reporters to build under a channel cap, and which modules that set can actually separate. |
| `bridge/physiology_bridge.py` | Branch A: measured physiology -> supply-limited bounds -> flux. Load-bearing. |
| `bridge/latent_bridge.py` | Branch B: latent state -> ATP maintenance. Refused by G4 at every call site. |
| `bridge/regulation.py` | D3: module activity -> regulon (SGD/MacIsaac binding motifs) -> the GEM's own gene-reaction rules -> per-reaction bounds. E-Flux, plus the lexicographic task-efficiency ordering that refuses without a measured flux. |

## Findings that constrain the plan

**D1 (design, no product data yet) — FBA bounds beta-carotene, it does not predict it.** The feasible range of
the product flux is `[0, ceiling]` at every growth level tested; relative width is
`1.000` throughout, with or without a capacity bound. A capacity ceiling moves the
top of the interval and never lifts the floor off zero. Product forecasting has to
come from a kinetic or expression layer; FBA's role is bounding, cofactor and
precursor accounting, and the growth-product trade-off curve.

**The ceiling was quoted as a bound and it does not hold.** This paragraph used to say
that at 90% of maximum growth the ceiling corresponds to about 33 mg/gDCW, that "reported
engineered titres of 5–20 mg/gDCW" sit 1.6–6.6x below it, and that the bound was therefore
informative. The 5–20 mg/gDCW range was not a literature figure at all — it is the
hardcoded grid `(5, 10, 20, 50)` in `scripts/parked/run_d1.py`, and no source was ever
cited. Scored against the actual literature, **the bound is refuted**: **two** published
beta-carotene contents in *S. cerevisiae* exceed the exact ceiling of 32.84 mg/gDCW —
79 mg/g DCW (PMID 39215465) and 46.5 mg/g DCW (PMID 33332529) — in the papers' own units,
with no unit conversion of ours in between. A third row, lycopene at 33.1 mg/g CDW
(PMID 32236768), matches it. The tightest surviving margin is 1.03x (PMID 33102463).

*Corrected 2026-08-30.* This read "four ... exceed" and cited 71.8 and 60.5 mg/g DCW from
**PMID 38710418, which was withdrawn by the publisher** and is struck from
`docs/EXTERNAL_CAROTENOID_BOUND.md`'s table. The retraction was recorded there and in
`MEASUREMENTS_NEEDED.md` and never propagated here, because the guard that exists to catch
exactly this scanned only that one file's table rows. The refutation does not need the
withdrawn paper: two independent groups break the bound on beta-carotene by HPLC.

The cause is the **growth fraction**, not the model and not the arithmetic. The
growth–product frontier is exactly linear, `q_max(mu) = 0.207469 (1 − mu/0.376829)`
mmol/gDW/h, so the content ceiling is `C(mu) = 111.38/mu − 295.58` mg/gDCW and **diverges
as `mu -> 0`**: FBA places no finite bound on specific content in a non-growing cell. Every
refuting row is feasible on the same ec model at the growth rate its 72–120 h harvest
implies. `growth_fraction=0.90` is `product_flux_range`'s default, not a measurement, so
any single mg/gDCW number is a claim about growth rate wearing the costume of a claim about
metabolism. What survives is the frontier itself, and a growth-free yield bound of
0.1984 g beta-carotene per g glucose, which the literature clears by 42–114x. Full scoring
of 28 convertible rows against 30 verified papers, the conversion arithmetic, and what one
extra measurement would make falsifiable:
[`EXTERNAL_CAROTENOID_BOUND.md`](EXTERNAL_CAROTENOID_BOUND.md).

Forcing 500 mg/gDCW-equivalent is still infeasible, and that statement is robust: the flux
it implies at mu = 0.3391 /h is 0.316 mmol/gDW/h, above the absolute maximum product flux
of 0.207 mmol/gDW/h at *any* growth rate.

**D3 (design, simulation only) — D1's conclusion survives, and for a stronger reason than
D1 gave.** D1 collapsed the whole 24-module latent state into one ATP maintenance scalar
(`bridge/latent_bridge.py`, `_MAINTENANCE_PER_ACTIVITY = 300.0`, which this repository's own
audit calls "a plausible scale, not a fitted one") and concluded that FBA cannot forecast
product. The obvious objection was that the conclusion might be about the *coupling* rather
than about FBA: a bound spanning `[0, ceiling]` is what you get when the only thing you tell
the model is "the cell is stressed, charge it ATP". `bridge/regulation.py` replaces the
scalar with per-reaction regulation and `scripts/run_regulation.py` tests it. **The
relative width does not move: 1.000 before, 1.000 after, at all seven of D1's growth
fractions, for all four stressor probes — 0 of 28 rows narrowed, and every regulated floor
is exactly zero.** The objection is answered and D1 stands.

The reason is structural rather than empirical, which is what makes it a sharper result
than D1's. `relative_width` is `(max − min) / max`, and `min` is zero because nothing in
the network *requires* product formation. E-Flux (Colijn 2009, PMID 19714220) constrains
`|v| ≤ c·e` — a bound on how much flux is *permitted*. No such constraint can make any flux
mandatory, so no capacity-based regulation layer, however specific, can lift that floor.
`EFLUX_LIFTS_LOWER_BOUND is False` is asserted in code, not argued in prose. A repression
depth sweep confirms it empirically down to 1% of reference capacity — at which point
growth has fallen to 1.8% of maximum and the floor is *still* zero.

What the mapping actually reached, so the null is not a null about nothing:

| | |
| --- | --- |
| module → gene | SGD curated regulation records restricted to MacIsaac 2006 (PMID 16522208) conserved binding motifs — binding evidence only, independent of any expression dataset. Same restriction, same reason, as `run_external_validation.py`. |
| gene → reaction | the GEM's own `gene_reaction_rule`. Not curated here. |
| reach | **19 of 24 modules have a factor; 15 reach ≥1 reaction; 154 distinct regulon genes are in the models.** ecYeastGEM: 680 reactions, 124 of them enzyme draws. Plain Yeast9: 334 reactions. |
| gaps, named | **UPR reaches zero reactions** — Hac1's conserved-motif regulon is 7 genes and the GSMM carries none of them, so the module two of the four wet-lab constructs report on cannot constrain metabolism at all. `calcium` has 0 records at this evidence level; `zinc` and `xenobiotic` have regulons the models do not carry; the 5 metabolite pools have no factor by construction. |

Three further findings, each of which was a surprise:

1. **On today's panel the layer can only *relax* bounds, never tighten them.** Every
   negatively-weighted arm in `STRESSORS[...].targets` points at a metabolite pool —
   `DTT → redox −0.90`, `menadione → nadh −0.40` — and a pool has no transcription factor,
   so no regulon. Across the four probes, **0 of the repressed module-arms have a regulon**
   and 0 of 48 candidate reactions were tightened. The signed direction is
   implemented and honoured (a negative weight lowers a bound; absolute values are never
   taken), but the panel supplies nothing for it to act on. That is a fact about
   `stress_panel.py`, not about the mechanism — and it is fixable, because the ESR is
   famously bidirectional and its repressed arm is simply not encoded yet.
2. **A single global protein pool is the wrong thing to scale.** The reference capacity is
   the load-bearing choice in E-Flux, and the honest one read off a GECKO model is
   `prot_pool_budget / enzyme_MW`. Measured, that ceiling is loose by more than an order of
   magnitude: at maximum growth the busiest enzyme in the H2O2 regulons draws **5.8%** of
   its own implied ceiling and the median draws none of it, so a scale of 0.45 is still
   non-binding. The literature predicts exactly this — Elsemman 2022 (PMID 35145105) and
   Dinh 2023 (PMID 37080482) both find the binding proteome constraint compartment-specific
   rather than global. Against a tighter *computed* reference (each reaction's own FVA
   maximum at 90% growth) the same repression does bite — growth falls to 0.45×— and the
   relative width is *still* 1.000.
3. **The ceiling moves the wrong way.** Under the tight reference, repression lowers growth
   and therefore *raises* the product ceiling 5.3–5.9×, because the growth–product frontier
   is linear and inverse. Any claim that a regulation layer "tightens the product bound"
   needs to check which end moved and why.

**What would lift the floor, and the measurement it needs.** COSMIC-dFBA (Gopalakrishnan
2024, PMID 38387677) §4.3 orders metabolic tasks lexicographically by *task efficiency* =
measured flux / individually maximised FBA flux, and fixes each leader at its measured flux
before re-ranking. That fixing is a lower bound, and it is the only mechanism here that can
lift a floor. It is implemented (`task_priority_order`, `apply_task_bounds`) and it
**refuses rather than defaults** when a task has no measured flux, because substituting a
simulated one puts the FBA maximum in both numerator and denominator, makes every efficiency
1.0, and turns the resulting narrow interval into an artefact of the substitution. Run on
this repository:

- On the tasks van Hoek 1998 measured, it **refuses on biomass**: the measured 0.40 /h
  exceeds ecYeastGEM's own maximum of 0.3768 /h, so biomass efficiency exceeds 1 and no
  priority can be formed. That is the D1.0 model-validity finding restated as a hard stop.
- On the two secretion tasks the model can support it works: ethanol efficiency **0.199**,
  CO2 **0.123**, ethanol ranked first. Two lower bounds lifted.
- And the product range is **still `[0, 0.01786]`, relative width 1.000**. A lower bound on
  other tasks does not put one on the product.
- On beta-carotene it refuses, as it must. **The missing measurement, named: the specific
  beta-carotene production rate `q_p` in mmol/gDW/h for a strain carrying
  crtE/crtYB/crtI**, obtainable from content in mg/gDCW plus the specific growth rate over
  the same interval via `q_p = content/MW · mu`. No such strain exists (`docs/PARKED.md`),
  so the numerator does not exist.

The honest reading of §4.3 is worth stating plainly, because it changes what the mechanism
promises: fixing a measured product rate collapses the interval **onto the measurement**
(relative width 0.000 at any `q_p`). It makes the model *consistent with* a measured product
rate; it does not forecast an unmeasured one. For COSMIC-dFBA that is sufficient — the
forecast comes from its ML layer and dFBA only propagates it. For this project's G-e it
means what D1 meant: **the product number has to come from outside the stoichiometry.** The
regulation layer is worth having for cofactor and precursor accounting and for the
growth–product frontier, and it retires the hypothesis that D1's width was an artefact of
the ATP-scalar coupling.

Precedent, so the novelty is not overstated: PROM (Chandrasekaran & Price 2010, PMID
20876091) is this exact pattern in *E. coli* and *M. tuberculosis* and is sixteen years old;
the field's own systematic comparison (Machado & Herrgard 2014, PMID 24762745) found
expression-integration methods frequently fail to beat plain pFBA. Everything above is
**Tier 1** — a deterministic computation on models and a binding map we adopted. The
magnitudes relating regulon activity to enzyme capacity are **asserted** and marked as such
in the source. `scripts/run_regulation.py`; `outputs/regulation_*.csv`.

**Model validity gates which model may be used — and the previous version of this
paragraph was wrong in both directions.** It said plain Yeast9 overpredicts growth by
**73%** and the ec model misses oxygen by **65%**. Both figures were artefacts.

The reference phenotype carried glucose 21.3, ethanol 27.4 and CO2 20.4 attributed to van
Hoek 1998. Those strings occur **zero times** in PMID 9797269, which ran no batch culture
— every steady state in it is a glucose-limited chemostat — and whose strain is DS28911,
not the CEN.PK113-7D the tests claimed. The reference now carries that paper's own Table 1
at D = 0.40: glucose 11.1, oxygen 3.7, ethanol 13.9, CO2 18.9. Separately, `cap_uptake`
was a no-op on the ec model because it is irreversibly split, so oxygen had never actually
been constrained there.

Corrected, against van Hoek's measured D = 0.40 phenotype:

| | growth | ethanol | CO2 | passes at 35% |
| --- | ---: | ---: | ---: | :---: |
| plain Yeast9, uptakes applied | **−13.5%** | +15.3% | +8.1% | yes |
| GECKO ec, uptakes applied | −12.9% | +13.2% | +2.7% | yes |

So plain Yeast9 does not overpredict; it *under*predicts by 13.5%, and both models
reproduce the phenotype to within about 15%. The ec model is better, marginally, not
dramatically.

What survives, and is the real reason to prefer the ec model: left unconstrained it puts
growth within **5.8%** of the measured rate from its protein pool alone. But it reaches
that growth on **61% more glucose** and **113% more ethanol** — a glucose-excess optimum
measured against a glucose-limited chemostat. The growth agreement is real and it is not
evidence that the flux distribution is. Oxygen-, OUR- and respiration-linked claims remain
unsupported: the ec model still misses oxygen by 27% and CO2 by 75% when left to itself.

**D2 — the reporter needs its own state, but the loss-rate alarm was an artefact of
the old plates.** Treating fluorescence as an instantaneous readout of state,
`y_t ~ p(y | z_t, ...)`, omits that reporter concentration relaxes to
`k_synth/(mu + k_deg)`. That structural correction stands regardless.

The *quantitative* alarm did not survive better data. On the two 2026-07 plates
(high inoculum, OD to 1.9, OD falling in most wells) inverting with `k_deg = 0` gave
impossible negative promoter activity in 52-72% of timepoints, implying a loss rate
of 0.049-0.074 /h. On the three **NewProtocol replicates** (2026-07-22, 07-28, 08-03)
the same analysis gives **zero negative activity and an implied loss rate of 0.0000 /h
in every replicate**. The apparent excess loss was driven by the corrupted OD
denominator, not by the reporter. mCitrine behaves as a stable YFP: loss is dilution.

What survives on the good data is the growth coupling, and it is *stronger* than pure
dilution predicts — fitted slope +2.4 to +5.6 against a predicted +1.0, on 3-20
estimable wells per replicate. Something couples reporter output to growth rate beyond
dilution. That is now the open question, not the loss rate.

Among cultures that pass G1 and admit the statistic, the median share of
`log(RFU/OD)` variance explained by growth rate alone is **0.82**, and the fitted slope is
positive in every well (+0.49 to +0.84 against the `+1.0` pure dilution predicts) — the sign
and roughly the magnitude the mechanism calls for. That rests on 8 cultures, and the
quasi-steady-state limit it depends on holds for only ~14% of each trace.

The preliminary experiment supports no conclusion either way: quasi-steady state holds there
only in wells with high, falling OD, which is exactly what G1 rejects, so the statistic is
estimable in none of its 53 G1-passing cultures.

**That coupling has now been attacked with a null, and it survives.** Both sides of the
regression are smooth autocorrelated series derived from the *same* optical density trace,
so they correlate spuriously as a matter of course. The right null is phase randomisation of
the growth-rate regressor — same power spectrum, same autocorrelation, phases destroyed —
and not a plain shuffle, which produces white noise that any smoothness-sensitive statistic
rejects for the wrong reason. `tests/test_nulls.py` demonstrates that failure mode directly.
Over 146 wells on the two NewProtocol plates that admit the statistic:

| | median R2 | IQR |
| --- | ---: | --- |
| observed | **0.522** | [0.293, 0.869] |
| phase-randomised null | **0.117** | [0.073, 0.175] |

The observed median is more than four times the null, so the relationship survives the
strongest like-for-like null available. Two qualifications, both real, and both belong next
to the 0.82 wherever it is quoted:

- **About 0.12 of any reported variance share is autocorrelation artefact.** Any well whose
  R2 sits near 0.12 is indistinguishable from smoothness. The honest statistic is the
  *excess* over the null, not the raw share.
- **It is a majority effect, not a universal one.** 83 of 146 wells (57%) beat their own
  null at p < 0.05, median p = 0.022 — which only just clears. So "growth explains most of
  the reporter signal" is supported on aggregate and is **not** established well by well,
  and the correction is applied well by well.

The 0.522 is not the 0.82 and the two are not in contradiction: this sweep did not apply the
G1 optical-quality gate, so it runs over a superset of wells including poorer ones. The
apples-to-apples comparison is observed against null on the *same* well set. Re-running the
sweep behind G1 is a small job and would give the gated equivalent. `docs/NULL_RESULTS.md`.

**G1 — the NewProtocol fixed it.** On the 2026-07 plates 53/90 and **8/42** wells pass,
the dominant failure being `linear_range`: cultures run to raw OD 1.8-1.9 in a 96-well
plate. On the NewProtocol replicates, inoculated at OD ~0.11, pass rates are
**80%, 100% and 86%** with **zero linear-range failures**. Remaining failures are
`sustained_decline` and low `dynamic_range`.

**G4 — the biosensors were tested against an independent anchor, and only one
construct is currently testable.** Anchors: Hac1 for the ER reporters, TRX2 for the
oxidative ones, both against UBC, over two biological replicates.

Replicate 1 (2026-07-24) was initially unusable: its export is a zip whose entries are
lowercased, which no OOXML reader will open, and once repaired its Target and Sample
columns are blank. The logbook records the qPCR plate layout and it is identical across
all three runs, so `read_positional_cq` restores the annotation from well position,
refusing to do so unless the data agree with the layout. That recovered 48 readings and
took the anchor from n=2 to **n=3** with no re-run and no re-export.

**Measured on real qPCR plates.** The gate's own verdicts, and then what the arithmetic of
`docs/research/G4_STATISTICS.md` does to two of them:

| Construct | Anchor | Stressor | Gate verdict | %gDNA of the +RT signal | Revised |
| --- | --- | --- | --- | ---: | --- |
| UPRE1 | Hac1 | DTT | INCONCLUSIVE (RT- margin +0.13 cyc, clears 0%) | **91.4%** | **REFUTED as a measurement** |
| UPRE2 | Hac1 | DTT | INCONCLUSIVE (RT- margin +0.60 cyc, clears 0%) | **66.0%** | **REFUTED as a measurement** |
| NativeYap1 | TRX2 | H2O2 | INCONCLUSIVE, anchor SNR 0.39 | 0.3% | INCONCLUSIVE, genuinely |
| AlteredYap1 | TRX2 | H2O2 | INCONCLUSIVE, anchor SNR 1.19 | 0.3% | INCONCLUSIVE, genuinely |

**The ER anchor did not fail to conclude. It failed as a measurement, and the difference
decides how the team should talk about it.** Cycles are a log-scale proxy; the quantity that
decides whether a reading means anything is the contamination fraction,
`%gDNA = 2^(-dCq) x 100`. At UPRE1's +0.13-cycle margin, genomic DNA is **91.4%** of the +RT
signal and transcript is 8.6%. Variance inflation on a subtraction is **15.7x** — a factor
independent of sigma, so it cannot be beaten by better pipetting — and at literature-typical
Cq noise (sigma 0.1–0.3 cycles) the correction returns a *negative* transcript quantity about
a third of the time. The only published gDNA correction, ValidPrime (Laurell 2012), is
validated to **60% gDNA** and emits `HIGHDNA` above it. Both ER constructs sit past that
ceiling. **Do not attempt a gDNA correction on this data.**

So for the ER pair, **nothing has been learned about the biosensors**. That is a protocol
finding, not a biosensor finding, and it is the honest thing to report. Only the oxidative
pair is evidence about the constructs at all.

**And the oxidative pair does not want 80 more plates.** `gates/g4_anchor.py` derives the
replicate count as `n x (target_snr / snr)^2` — about **80** for NativeYap1 at SNR 0.39 and
about **9** for AlteredYap1 at SNR 1.19 — and the quadratic is the whole story:
**doubling the effect size is worth four times as many replicates.** Two of the three anchor
RNA doses (0, 0.2, 0.5 mM) sit at or below the bottom of the 0.5–1 mM responsive window the
dose-response analysis identified, which is a plausible and cheap explanation for the small
effect. Sampling the anchor at 1.0 mM, or picking a better-behaved anchor gene, may recover
more SNR than any number of extra replicates.

One correction to received wisdom, since it is what the 3-cycle gate is judged against:
**MIQE never mandated a 5-cycle RT- margin.** The number is a software default whose real
content is "<=3% background", and it originated as a *no-template-control* rule. MIQE 2009,
2010 and 2.0 (2025) require only that the +RT/-RT comparison be reported and that the
experimenter state their own tolerance — and MIQE 2.0 explicitly endorses correcting rather
than discarding, which is why the arithmetic above matters more than the convention.

**The third replicate withdrew the one PASS.** At n=2 `AlteredYap1` cleared the bar with
an anchor effect of +0.75 fold at SNR 1.6. Replicate 1 puts its 0.5 mM fold at 1.18
against 2.09 and 1.41, which widens the between-replicate spread and drops the SNR to
1.19. The n=2 result was optimistic; no construct currently passes G4.

The two sensor pairs are challenged with **different agents**, so their dose axes are
not comparable and `assert_single_stressor` refuses to pool them.

The ER anchor is not measuring transcript, in **all three** replicates. Its
no-reverse-transcriptase control amplifies at essentially the same cycle as the +RT
sample, so the Hac1 signal in the UPRE samples is largely genomic DNA. The logbook names
the primer pair — **"Hac1 Exon"** — which explains both failures at once: exon primers
amplify the genomic locus readily, and they measure *total* HAC1. Two fixes, independent of each other: DNase-treat
the RNA and use intron-spanning primers, and pick a better UPR readout. `HAC1` is
unconventionally spliced by Ire1 during the UPR, so *total* HAC1 transcript barely
moves even when the pathway is fully on — which is what the data show. **KAR2** (BiP)
is the conventional downstream UPR target and would give a real dose response.

The oxidative anchor is clean (RT- margin 8.6 cycles, 100% pass). It is simply
underpowered at n = 2.

The one PASS carries a caveat the gate states itself: growth rate per dose is not
available without a well-to-dose plate map, so the growth confound is untested for it.

**Sensor characterisation — most of the apparent dose response is growth dilution.**
The plate map is not in the exports, so it was *recovered*: searching every
construct-to-column assignment for the one that reproduces the team's own per-construct
sheets from the raw wells. Both usable plates independently return the same answer —
UPRE1 cols 1-3, UPRE2 4-6, NativeYap1 7-9, AlteredYap1 10-12, rows A-G ascending dose —
at 1.8% and 3.4% median error.

**The wet-lab logbook independently confirms that map**, and supplies things no
recovery could reach: blank positions per plate, the OD 0.75 inoculation target, and
the fact that RNA for the qPCR anchor was taken from rows A, C and D (0, 0.2, 0.5 mM)
of the same plates after each read — so the anchor and the reporter are the *same
cultures*, not merely the same conditions. Layouts are transcribed in
`plate/layout.py::RECORDED_PLATES`.

Two discrepancies between logbook and export, both recorded rather than silently
resolved. The 2026-08-03 entry lists blanks at H4-H6, but that export contains no H4-H6
wells and its H1-H3 behave like medium (OD drifts 1.10x while the reporter stays flat at
433 -> 426, against a culture's 639 -> 1705); the entry looks copied from the 07-28
template, and H1-H3 is used. And the H2O2 stock preparation tops out at 4 mM while the
plate map writes 5 mM in row G; the exports carry 4 mM, which is what is used.

With growth available, the reporter's apparent induction can be split from the growth
slowdown the stressor causes. Fold change at the end of the run, relative to 0 mM.

**Measured on real plates, at four biological replicates.** The table below has been
recomputed twice: at n=3 when a parsing fault that was dropping a plate was found, and at
n=4 when `20260804`'s instrument file was read in place of its blank-subtracted export.
Both times every magnitude moved, and both times the direction of travel was the same —
more rows resolve. The n=2 version said none did; seven do now.

| Construct | Dose | Naive fold | Corrected, plate-paired | 95% CI | Verdict | Reading |
| --- | --- | ---: | ---: | --- | --- | --- |
| NativeYap1 | 0.1 mM H2O2 | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->1.28 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->0.92** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->0.74, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->1.08] | INCONCLUSIVE | induction unresolved; practical equivalence not established |
| NativeYap1 | 0.2 mM H2O2 | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.2;late_fraction=0.75" -->1.43 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.2;late_fraction=0.75" -->1.03** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=NativeYap1;stressor=H2O2;dose_mM=0.2;late_fraction=0.75" -->0.80, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=NativeYap1;stressor=H2O2;dose_mM=0.2;late_fraction=0.75" -->1.27] | INCONCLUSIVE | induction unresolved; practical equivalence not established |
| AlteredYap1 | 0.1 mM H2O2 | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->1.28 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->1.01** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->0.94, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.1;late_fraction=0.75" -->1.10] | INCONCLUSIVE | induction unresolved; practical equivalence not established |
| AlteredYap1 | 1.0 mM H2O2 | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=1.0;late_fraction=0.75" -->1.52 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=1.0;late_fraction=0.75" -->1.00** | REFUSED: invalid bootstrap ratios | INCONCLUSIVE | point estimate retained; old interval endpoints withdrawn |
| UPRE1 | 1.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE1;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.83 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.44** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE1;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.22, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE1;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.68] | **DIFFERENT** | real, ~25% smaller than it looks |
| UPRE2 | 1.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.93 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.56** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.30, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE2;stressor=DTT;dose_mM=1.0;late_fraction=0.75" -->1.83] | **DIFFERENT** | real, ~24% smaller than it looks |
| UPRE1 | 2.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE1;stressor=DTT;dose_mM=2.0;late_fraction=0.75" -->2.02 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=2.0;late_fraction=0.75" -->1.06** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE1;stressor=DTT;dose_mM=2.0;late_fraction=0.75" -->0.92, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE1;stressor=DTT;dose_mM=2.0;late_fraction=0.75" -->1.21] | INCONCLUSIVE | dilution-sensitive; no equivalence established; see the sweep |
| UPRE1 | 5.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE1;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->1.76 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->0.93** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE1;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->0.80, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE1;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->1.09] | INCONCLUSIVE | dilution-sensitive; no equivalence established; see the sweep |
| UPRE2 | 5.0 mM DTT | <!-- audit:value table=outputs/late_window_sensitivity.csv column=naive_fold row="construct=UPRE2;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->2.02 | **<!-- audit:value table=outputs/late_window_sensitivity.csv column=fold row="construct=UPRE2;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->1.06** | [<!-- audit:value table=outputs/late_window_sensitivity.csv column=low row="construct=UPRE2;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->0.93, <!-- audit:value table=outputs/late_window_sensitivity.csv column=high row="construct=UPRE2;stressor=DTT;dose_mM=5.0;late_fraction=0.75" -->1.23] | INCONCLUSIVE | dilution-sensitive; no equivalence established; see the sweep |

Every cell above is pinned to `outputs/late_window_sensitivity.csv` at the canonical late
window, so this table cannot drift from the file that produced it without
`scripts/audit_claims.py` failing. It drifted once, which is why they are there.

**The resolvable set is 0.2–1.0 mM DTT on both ER sensors plus 0.5 mM H2O2 on
AlteredYap1**, seven calls in all: UPRE1 and UPRE2 each clear "no change" at 0.2, 0.5 and
1.0 mM, and AlteredYap1 clears it at 0.5 mM. At n=3 it was five, and neither oxidative
sensor produced a resolvable induction at any dose.

**No row resolves in the other direction at any window any more, and three did.** The four
H2O2 rows above 1.0 mM — both oxidative sensors at 2.0 and 4.0 mM — have recovered activity
below background, so fewer than three plates return a positive per-plate fold and the
geometric mean rests on too few clusters to bootstrap. They are now inestimable at all six
windows, which makes the sweep exactly 20 of 24 estimable everywhere rather than 20 at the
narrow windows and 22–23 at the wide ones.

Two of the three losses are an estimator fix and one is the plate, and they are worth
keeping apart:

| row | window | at n=3, as shipped | at n=4 | why it went |
| --- | --- | --- | --- | --- |
| AlteredYap1 @ 4.0 mM H2O2 | 0.8 | 0.007 [0.006, 0.883] | not estimable | estimator |
| NativeYap1 @ 2.0 mM H2O2 | 0.9 | 0.038 [0.024, 0.318] | not estimable | estimator |
| UPRE2 @ 2.0 mM DTT | 0.9 | 0.831 [0.728, 0.965] | 0.916 [0.734, 1.249] | the plate |

The first two were dying cultures reported as confident findings of strong repression.
`analysis/uncertainty.py::fold_change` now counts the plates that actually entered the
geometric mean rather than the plates that carried the dose, and a fold resting on fewer
than three positive plates — these rest on one or two — is refused instead of being widened
as though it rested on four. **The
fourth plate is what exposed it**: the point estimate does not move when a plate that
contributes nothing is added, but the t quantile does, and `AlteredYap1` at 2.0 mM went from
[0.010, 1.349] to [0.013, 0.473] on that alone.

**That fix changes nothing at the canonical window** — every fold, bound and estimable flag
at `late_fraction=0.75` is bit-identical between the three-plate table under the old rule and
under the new one — so all of 5 → 7 is the fourth plate and none of it is the estimator.

The AlteredYap1 rows are the ones worth dwelling on, because between them they are the trap
and the escape from it. At 1.0 mM the point estimate was 1.49 at n=3, close enough to
UPRE2's 1.49 that a reader scanning the fold column would have called them the same result;
the intervals said otherwise, [0.82, 2.73] against [1.20, 1.74]. The fourth plate settled
which was which — UPRE2 went to 1.56 [1.30, 1.83] and AlteredYap1 to **1.00** [0.24, 2.50].
Meanwhile the same sensor's 0.5 mM row, which nobody was quoting, became a call. The lesson
is the one the intervals were already carrying: read the width, not the point.

**The fold is computed within each plate and then combined**, rather than as a ratio of
pooled means. A plate-level shift that multiplies a plate's dosed and control wells alike
cancels in the first and does not cancel in the second. At four replicates the two rows that
still move on that change alone are `AlteredYap1` at 1.0 mM (1.26 -> 1.00) and `NativeYap1`
at 1.0 mM (1.31 -> 1.18); every other row is unchanged to two decimals. Both are H2O2 rows
at the top of the usable range, which is where the plates disagree most.

**Every row rests on four biological replicates.** That took three separate corrections to
the reading of the files, none to the cultures: a parser that stopped after 21 wells of 87
on a column-stacked export (which recovered 20260728), an operator's blank-subtraction
transform that had to be undone by reading the instrument file rather than its export (which
recovered 20260804), and a logbook blank position that the instrument file contradicted.

Three was the threshold `analysis/uncertainty.py::fold_change` requires; four is the first
count at which a family-wise correction over the twenty estimable folds leaves anything
standing. INCONCLUSIVE is still the majority verdict and still means unresolved rather than
refuted. `docs/DATA_INVENTORY.md` has the per-plate coverage.

### Autofluorescence has never been measured, and it splits this table in two

> **Superseded in part, 2026-08-30.** The claim below that both reporter-free plates were
> read in OD600 only is **false**, and so is "two such strains". The raw
> `20260807_BY4741_ER_Oxidative.xpt` carries `mCitrine:480,530` for all 96 wells under the
> same protocol as replicates 2–4; the `.xlsx` export dropped the channel. The 2026-08-14
> file is the same plate exported again, not a second plating. So `a` is a data-reduction
> task, not a plate to run. Everything below about what an unremoved `a` *does* — the
> `a * mu` contribution and the sweep — is unaffected and still stands.
> See `docs/research/XPT_INVENTORY.md` §2.2.

`observation.py` writes the measurement model as
`RFU = background + [autofluorescence + gain * R] * biomass * exp(-eps * carotenoid)`, and
its docstring insists autofluorescence be measured on an isogenic reporter-free strain
rather than guessed. **A reporter-free strain was plated and its `.xlsx` export was OD600
only** — which is what this entry originally recorded as "never measured". The instrument
file kept the mCitrine channel; see the note above.

What the analysis subtracts instead is the media blank, which is the `background` term.
Autofluorescence is a *per-biomass* term and survives that subtraction entirely, so the
per-cell signal still contains it. The ODE inversion makes this exactly sweepable: an
unremoved constant `a` contributes exactly `a * mu` to the recovered activity, since a
constant has zero derivative. Sweeping `a` as a fraction of each construct's own 0 mM signal:

| Construct | Dose | a=0% | a=5% | a=10% | a=20% | a=30% |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| NativeYap1 | 0.1 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=NativeYap1;stressor=H2O2;dose_mM=0.1;autofluorescence_fraction=0.0" -->0.92 | 0.92 | 0.92 | 0.92 | 0.91 |
| AlteredYap1 | 0.1 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=0.1;autofluorescence_fraction=0.0" -->1.01 | 1.01 | 1.02 | 1.02 | 1.02 |
| AlteredYap1 | 1.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=AlteredYap1;stressor=H2O2;dose_mM=1.0;autofluorescence_fraction=0.0" -->1.00 | 1.01 | 1.03 | 1.06 | 1.10 |
| UPRE1 | 2.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=2.0;autofluorescence_fraction=0.0" -->1.06 | 1.08 | 1.10 | 1.16 | **1.22** |
| UPRE1 | 5.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=UPRE1;stressor=DTT;dose_mM=5.0;autofluorescence_fraction=0.0" -->0.93 | 0.95 | 0.97 | 1.01 | **1.06** |
| UPRE2 | 5.0 | <!-- audit:value table=outputs/autofluorescence_sensitivity.csv column=fold row="construct=UPRE2;stressor=DTT;dose_mM=5.0;autofluorescence_fraction=0.0" -->1.06 | 1.09 | 1.11 | 1.16 | **1.22** |

At four biological replicates, and the a=0% column is pinned cell by cell so it cannot drift
from the file that produced it. Every row moved when the fourth plate landed; the two that
moved furthest are AlteredYap1 at 1.0 mM, whose whole sweep dropped by about half a fold,
and the high-DTT rows, which all rose by 0.03 to 0.06.

**The low-dose point estimates are relatively insensitive to this assumed background sweep.**
That sensitivity result does not establish "no induction at all": no justified practical-equivalence
margins were supplied, and an interval containing one is not evidence of equivalence.

**The higher-DTT point estimates are more background-sensitive.** Changes across an assumed
background grid do not by themselves establish either "entirely dilution" or "a small real
induction". Both earlier interpretations are withdrawn; the intervals, unavailable intervals,
and independent plate counts must remain part of the conclusion.

The size of `a` is not identified by this sweep. Matched non-reporter controls can constrain
background, but one additional plate does not by itself establish biological equivalence or
resolve every inferential uncertainty.

Growth correction substantially changes the dose-response interpretation. It separates a
reporter's per-cell dilution term from the naive readout under the stated observation model;
it does not prove that every residual change is a specific stress response. The lower-dose
non-detections and higher-dose background sensitivity therefore remain distinct limitations,
not proof of a universal responsive window.

Caveats, in order of size. Two plates, which for three of the four constructs is also two
biological replicates in the fluorescence channel. `NativeYap1` was inoculated 3.7x lower than
`UPRE1` (blank-corrected OD 0.044 vs 0.162) — **the logbook explains why**: its overnight
culture came up at OD 0.9 against 3.7 for the others, so it was made up in 10 mL rather
than 3 mL, a different dilution from every other well on the plate. That leaves it close
enough to the detection floor that its growth estimate stays above the physiological
maximum for yeast, so its correction is over-applied and its corrected folds are a lower
bound. Both Yap1 constructs go below background above 1 mM, so those rows are
uninterpretable.

## Consequences for the plan

1. `y_t ~ p(y | z_t, m_t, x_t)` in §8.1 needs the reporter as its own ODE state.
   Built: `reporter.py`, `observation.py`, and the filter in `estimator.py`, which
   infers promoter activity and growth rate jointly with uncertainty.
2. Measure `k_deg` for the actual construct (`docs/PROTOCOLS.md` P2). On the good
   plates the free lower bound is already 0, so the chase is now a confirmation
   rather than a rescue.
3. Measure the OD linear range (`docs/PROTOCOLS.md` P1) and set
   `OpticalQualityGate.od_linear_max` from it. The NewProtocol inoculation density
   should be kept.
4. Any oxygen-linked latent module is refused by `archive/code/module_admission.py` until a
   model reproduces respiration or the anchor is measured directly.
5. Product forecasting moved off FBA onto `kinetic/carotenoid.py`; FBA keeps bounds,
   cofactor accounting and the precursor supply that feeds the kinetic branch.

## Measured estimator behaviour

**Tier 1 — self-consistent in simulation.** The estimator recovers a known input from data
it generated itself. Read the next subsection before quoting any of it.

The promoter-activity inversion recovers a known input to **1.1% median relative
error at 1% reader CV** (correlation 0.999) with 10-minute sampling. Error scales
about linearly with noise (5% CV → 5.4%). The smoothing window is the sensitive
knob: 0.5 h gives 14% error, 4–8 h gives ~1%, 12 h begins to bias. The default
scales with run length rather than being a constant.

`implied_min_k_deg` uses `-d(ln R)/dt - mu`, not `-(dR/dt)/R - mu`. On a known
0.30 /h decay the ratio form drifts to 4.84 as the window widens to 8 h; the log
form stays within 0.30–0.32.

### On real traces, asked to forecast, the twin loses

That 1.1% is the estimator inverting its own generator. Asked instead to predict readings it
had not seen, on real plates, it is beaten by trivial extrapolation — and this is the only
scored forward prediction the project has. Full protocol and baselines in
`docs/research/PREDICTION.md` §0.

Conditioning `ParticleFilter` on the raw OD600 and mCitrine readings for `t <= 2.0 h`
(13 of 25 timepoints) on 84 real culture wells, then forecasting the remaining 12 timepoints
with no further observations:

| horizon | OD RMS rel. err. | best fair baseline | FL RMS rel. err. | best fair baseline |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 h | 16.5% | **7.7%** (causal climatology) | 12.0% | **4.7%** (log-linear) |
| 1.0 h | 26.5% | **12.5%** | 14.6% | **11.3%** |
| 1.5 h | 38.2% | **16.7%** | 20.5% | 21.7% |
| 2.0 h | 53.5% | **21.7%** | **28.4%** | 36.3% |

Empirical coverage of the nominal 95% forecast interval is **0.40–0.49** on both channels.
Paired per well against a six-point log-linear extrapolation — the fairest single comparison,
since both see exactly the same 13 points from that one well — the baseline wins on OD by
23.0% of level at 2 h and ties on mCitrine.

Two things make this worth having rather than merely embarrassing, and one thing this section
used to count as a third has been withdrawn.

**The loss is diagnostic, not diffuse.** `|bias|/sd` is 2.17 on OD and 1.33 on FL, so the
error is systematic. The bias is positive at every horizon: the twin over-predicts. The cause
is in `estimator.py::_propagate` — growth rate is a driftless random walk, so the forecast
carries the last inferred `mu` forward, and these cultures decelerate hard inside a four-hour
run (`mu_max` 0.557 falling to 0.277 on the 2026-07-22 plate). A model with no deceleration
term extrapolating a decelerating culture must overshoot, and does.

**A second, independent statistic was claimed here, and it is withdrawn.** This paragraph
read that a chi-squared test on the filter's own innovations reported the posterior as
over-confident, "median time-averaged NIS 1.76 on OD, 4.44 on mCitrine, against an
acceptance band of [0.52, 1.63]", and that the two tests together trapped the defect because
no setting of the random-walk scales satisfied both. That pair is a 2000-particle research
probe on 167 wells over 2 plates at K = 25 — `data/current_claims.json` files it as
`nis.research.shipped`, role `historical`, binding `pending`, because its configuration and
innovations were never archived. The shipped CLI, regenerated, does not reproduce it and does
not replace it: `outputs/nis_channel_summary.csv` covers 319 well-channels per channel over
4 plates at K = 22 and reports **`tested_well_channels: 0`** with verdict **INCONCLUSIVE**,
every well-channel refused below the ESS floor and 208 of 638 below the particle-diversity
floor as well. The registry states the consequence as three refusals (`nis.current.od`,
`nis.current.rfu`, `nis.current.channel_refusal`).

The descriptive medians on that run are 45.31 (OD) and 19.32 (mCitrine) against the K = 22
band [0.499, 1.672], but `outputs/nis_manifest.json` says in its own `aggregation_caveat`
that they have no chi-squared null, so they are not a verdict of over-confidence — they are
what an untested channel happens to look like. **The calibration gate currently returns no
test.** The forecast loss above is therefore the only live measurement of this defect, not
one of two, and the "no setting of the walks satisfies both" argument has no second test to
stand on until a configuration that clears the ESS and particle-diversity floors is run and
archived.

**The reporter ODE is the one part that earns its keep.** On mCitrine the twin is beaten at
short horizon and pulls level at 1.5–2 h, exactly where naive log-linear extrapolation starts
overshooting because an exponential fit cannot know that `dR/dt = k - (mu + k_deg) R`
saturates. It is not yet enough to win.

And one methodological warning bought for free: on aggregate RMS the twin *wins* on
fluorescence at 2 h; on the paired per-well mean-absolute comparison it *ties*. The
disagreement is real — the baseline has fat tails that RMS punishes and MAE does not, while
the twin is uniformly mediocre. **Which scoring rule you pick changes the verdict**, which is
the single strongest argument for pre-registering the rule before looking. Doing that is the
only route to a Tier 3 claim, and it needs no new bench work.

## Not yet built

dFBA integration over the full culture and gates G3-G5. The estimator, kinetic branch and
surrogate they depend on are now in place. The augmentation generator is not "not yet built"
but **removed** — its gates live on the panel path in `generator/panel_gates.py`; see
[Removed](#removed).

**Three states, kept apart, because this file used to run them together.**

| | Status | Why |
| --- | --- | --- |
| **G-d**, generalising the state to an unseen stressor | **Unreachable on these plates** | Two stressors, and each read only by its own construct pair — UPRE1/UPRE2 read DTT, NativeYap1/AlteredYap1 read H2O2. Hold either out and every channel reading its private modules goes with it, leaving at most the two channels dosed with the other agent. The number cannot be computed here; it does not evaluate to zero. Where it *could* be run — 18 pairs on Gasch GSE18 — the fitted basis **fails** its rotated null, +0.688 against +0.695, p = 0.990. |
| **G-e**, forecasting the product | **Unbuilt, not blocked** | The forecast has to come from a kinetic or expression layer, and the calibration set for one now exists (Elizondo & Saa 2025, PMID 40891387 — six steady states, product rate, growth rate and uptakes in the paper's own units). Steady states, not a time course, and no strain here makes the product. `docs/research/CALIBRATION_DATA.md`. |
| Held-out **dose** prediction | **Works, and is scored** | Positive skill at every interior rung at or below 0.5 mM, negative above it, sign separating on `dose / EC50`. Not Tier 3 — scored after the measurement. `docs/WHY_INTERPOLATION_WORKS.md`. |

`docs/CLAIM_BOUNDARY.md` states each of the three at length.

## Data

Analysis reads from sibling directories, not from this package. The
`Biosensor Testing/` NewProtocol replicates are the usable set: 7 DTT doses x 3
technical replicates x 4 constructs (UPRE1, UPRE2, NativeYap1, AlteredYap1), across four
plates.

**"Four biological replicates" is now true of the fluorescence channel too**, which is the
one every reported result rests on, and it took two corrections to get there. This table read
`21` and `UPRE1 only` for two of the four plates; both readings were faults in this
repository's readers rather than in the plates. What the files actually hold:

| Plate | Source read | OD wells | Fluorescence wells | Blanks in fluorescence |
| --- | --- | ---: | ---: | --- |
| `20260722` | `.xlsx` export | 87 | **87** | H1-H3 |
| `20260728` | `.xlsx` export | 87 | **87** | H4-H6 |
| `20260803` | `.xlsx` export | 87 | **87** | H1-H3 |
| `20260804` | **`.xpt` archive** | **96** | **96** | H1-H3 |

Replicates 2 and 4 looked like 21-well UPRE1-only reads because they write the plate down
the page instead of across it — one block per group of columns, only the first carrying the
Time column — and the reader stopped after the first block. Replicate 4 then looked
unusable a second time because every `.xlsx` of it had been exported through a
blank-subtraction transform; `plate/gen5.py` reads the instrument file, which never had one
applied. Nothing was ever wrong with those cultures. `docs/DATA_INVENTORY.md`.

Several workbooks stack a raw block and a blank-subtracted block in one sheet;
`SynergyRun.raw_channel` selects the unblanked one, detected structurally by the presence
of negative values. That heuristic is why `20260804` could not be rescued from its export:
the workbook holds *only* the subtracted block, so there is no unblanked one to select.

**No plate has ever been read for autofluorescence.** `20260807_BY4741_ER_Oxidative.xlsx`
(33 wells) and `20260814_BY4741_0mM_OD600 Data.xlsx` (6 wells) are the reporter-free
controls, and both are OD600 only.


## The whole stress landscape, and two kinds of sensor

The panel is no longer seven regulons. `generator/stress_panel.py` now carries 24 state
variables and 25 stressors: 19 transcriptional regulons (ESR, UPR, oxidative, heat,
osmotic, proteasome, iron, cell wall, DNA damage, calcium, carbon, nitrogen, hypoxia,
copper, zinc, sulfur, xenobiotic, retrograde, alkaline pH) and 5 metabolite pools
(glutathione redox, peroxide, ATP:ADP, cytosolic pH, NADH:NAD).

The pools matter because they are read by a different kind of sensor. A promoter fusion
has to be transcribed, translated and matured, accumulates over hours, and is diluted by
growth, so reading activity off it means dividing out a growth rate estimated from noisy
optical density. A ratiometric biochemical sensor -- roGFP2-Grx1, HyPer7, QUEEN, pHluorin,
Peredox -- reports an equilibrium between two spectral forms of one molecule: seconds to
settle, reversible, and the ratio cancels concentration so growth never enters. The whole
dilution problem this repository is built around simply does not apply to them.

That makes them look strictly better on information, and an unconstrained search fills
every channel with them. It should not: they are all GFP- or cpYFP-derived
excitation-ratio sensors emitting near 510 nm, so two of them in one strain give one
measurement. `Reporter.spectral_slot` encodes that, and selection refuses builds that
cannot be made. A promoter fusion claims no slot, because the promoter is separable from
the protein and can be built onto whichever fluorophore is free.

The recommended build is therefore a mix: **STRE-general, UPRE-ER, PACE-proteasome,
roGFP2-Grx1, and mOrange2 as the constitutive reference**. Four measurements plus growth,
one of them immune to the growth correction entirely.

## Design on Fisher information, not on a spanning score

Three earlier criteria were ad hoc and all had the same root cause -- no noise model, so
"identifiable" had no statistical meaning. `analysis/design.py` replaces them:

- information is `L' S^-1 L`, and a module's standard error is the matching diagonal of its
  pseudo-inverse, so identifiability is a precision claim against an effect worth detecting
  rather than a projection weight clearing 0.35;
- noise is multiplicative, so a reporter that reads more loudly carries proportionally more
  noise and gains nothing. Only orthogonality survives. Row-normalising was the right
  instinct in the old score; here it falls out of the noise model instead of being asserted;
- a build is chosen by weighting the modules that matter, not by forcing a reporter in;
- the score is a log pseudo-determinant, so no ridge silently prices one more estimable
  module against better conditioning.

The rank ceiling is exact and unforgiving: **N channels resolve N state variables**, whichever
N you choose out of 24. No module is intrinsically unreachable -- the earlier claim that
ESR, osmotic and proteasome were out of reach was the 0.35 threshold talking.

## The wiring, corrected against the literature

The first version pushed a general-stress arm into nearly every regulon. Most of those
edges are not real, and removing them changed the conclusions:

**Tier 0 — transcribed from the literature, then corrected against it.** These edges are asserted inputs to everything on the known-loadings path.

| edge | verdict |
| --- | --- |
| ESR to UPR | do not encode -- Ire1 splices HAC1 mRNA, and HAC1 has no STRE |
| ESR to heat | do not encode -- Hsf1 and Msn2/4 act in parallel with compensatory, negative coupling |
| ESR to oxidative | unsupported |
| ESR to osmotic | arrow reversed -- Hog1 acts on Msn2/4 activity, and that is signalling |
| ESR to proteasome | weak and redundant -- RPN4 has no STRE, and reporter and RNA-seq disagree |
| oxidative to proteasome | solid -- Yap1 to RPN4 through a YRE, mutagenesis-confirmed |
| heat to proteasome | solid -- Hsf1 to RPN4 through an HSE, mutagenesis and ChIP |
| iron independent of ESR | correct, with one weak single-source exception, Yap1 to AFT2 |

Promoter crosstalk is now separate from the cascade. `Module.driven_by` is one factor
transcribing the gene for another; `Reporter.also_reads` is a second response element in a
promoter. One weight was standing in for two mechanisms and double-counting the first.

## Transfer: train on some stressors, predict one never seen

This is the test of whether a stress state is general rather than a lookup table. Withhold
a stressor entirely, fit on the rest, reveal some of its channels and predict the others,
against an oracle allowed to see it. `analysis/transfer.py`.

Two things had to be fixed before the numbers meant anything. Reading a state off unseen
data by least squares is an ill-posed inverse problem, and it produced R2 of -1e7; the fix
is a MAP state under the prior the training states already define, with the noise taken
from the cross-validated error rather than the in-sample residual, which collapses to zero
once the states saturate the channels. And the design question had to stop being proxied:
three analytical measures of "does this design span the space" failed, so designs are now
scored by simulating the transfer they deliver.

**The experiment decides whether transfer works, not the biology.** One stressor at a time
gives one ray per stressor in module space, so a held-out stressor is reachable only if
another happens to point the same way:

**Tier 1 — self-consistent in simulation**, and on the *fitted*-loadings path, where Ledermann caps four channels at one factor. See the two-paths note above.

| design | wells | median transfer R2 |
| --- | --- | --- |
| blocked, one stressor at a time | 144 | 0.064 |
| all 66 pairs | 936 | 0.353 |
| **6 selected pairs** | **216** | **0.539** |

Six chosen pairs beat all sixty-six at a quarter of the wells, verified on held-out noise
seeds. More combinations do not keep adding information -- the fitted subspace follows
whatever the design samples most, so redundant mixtures crowd it with covered directions.

Transfer tracks how much of a held-out stressor lies in the fitted subspace, which is the
mechanism and the honest limit. Well-aligned stressors reach 92-95% of oracle; stressors
whose modules nothing else excites do not transfer, and no amount of modelling fixes that.

The transferable state wants **3 dimensions** out of 24. Entry-wise cross-validation says 1,
and is not wrong -- it asks which rank best fills a gap in a well whose other channels are
known. Selection has to follow the task.

Latent states do carry the modules rather than merely compressing them: recovery is 0.93 to
0.99 for regulons the panel reads and excites, against a shuffled control near 0.02.

Run `python scripts/run_transfer.py` for the whole report.


## The trained model

`analysis/stress_model.py` is the artefact that gets used: reporter readings in, named
module activities out. `scripts/run_training.py` fits it and writes it to `outputs/`.

Two pieces are learned, and they are learned differently. The latent basis comes from the
readings alone and is unsupervised, so it transfers to any plate read on the same channels.
The map from that basis to named modules is supervised, and simulation is the only place
the labels exist to fit it. That is what a mechanistic generator buys: the readout is
fitted where truth is known and then applied to real readings, which carry no labels.

Inference is MAP under the prior the training states define, so a plate outside the
training range is pulled toward the training mean rather than producing an enormous state
that happens to match the channels it was shown. The model records which channels it
expects and refuses anything else, so it cannot be handed the wrong columns silently.

### Latent width has to be chosen on the job the model will do

Selecting width on channel prediction picked **one** state for the four-channel build, which
then reported almost nothing for the glutathione pool -- despite the build carrying a
glutathione sensor. One state cannot hold the general stress response and a metabolite pool
at once. Selecting on module recovery instead (`select_width_for_modules`) picks the full
width and roughly doubles what the model reports:

**Tier 1 — self-consistent in simulation**, and on the *fitted*-loadings path, where Ledermann caps four channels at one factor. See the two-paths note above.

| latent width | redox recovered | held-out mean over driven modules |
| --- | --- | --- |
| 1 | 0.07 | 0.275 |
| 2 | 0.92 | 0.416 |
| 3 | 0.94 | 0.489 |
| 4 | 1.00 | 0.573 |

The held-out column here predates the `module_transfer` correction and is restricted to
*driven* modules, which the correction affects least; the width ordering it selects on is
unchanged, and `select_width_for_modules` still picks the full width.

On a narrow panel the answer is the full width, and that is not a failure. There is nothing
to compress out of four channels; what the latent layer buys there is the readout that
names modules, not the reduction. Compression earns its place only when channels outnumber
the structure behind them.

### What it reports

Scored leave-one-stressor-out in module units -- trained without a stressor, does it report
the right biology when that stressor arrives? Both stages are blind to it.

**Tier 1 — self-consistent in simulation**, and on the *fitted*-loadings path, where Ledermann caps four channels at one factor. See the two-paths note above.

| model | channels | latent width | in-sample | median modules recovered, held out | mean score on the modules it drives |
| --- | --- | --- | --- | --- | --- |
| every promoter fusion plus one ratiometric sensor | 20 | 5 | 10 / 24 | 5 / 24 | 0.448 |
| **the five-channel build that can be made** | **4 + reference** | **5** | **11 / 24** | **1 / 24** | **0.249** |
| the three-sensor build | 3 + growth | 4 | 9 / 24 | 0 / 24 | 0.111 |

**In-sample and held-out are different numbers and the gap is the whole story.** The
buildable panel describes 11 of 24 modules on plates it was fitted on and a median of
**one** on a stressor it has never seen. The three-sensor build reaches zero, and the
20-channel panel — twice the width — reaches five. Note the buildable panel now describes
*more* in-sample than the wide one (11 against 10): nothing is compressed away when
latent width equals channel count. That is not
a failure of fitting; it is what one ray per stressor in module space buys.

An earlier reading of this table gave 15 / 24 and 10 / 24 as *held-out* counts. Those came
from scoring modules whose held-out truth was identically zero, where predicting "off"
beat the training mean for free. The retraction, the mechanism and the corrected numbers
are in [`superseded/module-transfer-inflation.md`](superseded/module-transfer-inflation.md).

What transfers is set by redundancy, not by what the panel reads. Across all 76 (held-out
stressor, own target) pairs, a module with no other stressor driving it is recovered 0 times
out of 6; with six or more peers, 84% of the time (Spearman +0.59). The three-sensor build
carries a dedicated ATP channel and still recovers nothing for `antimycin_A` or
`glucose_starvation`, because each is the panel's only real driver of its own axis.

On the modules it does read the buildable panel is the *better* model in-sample -- UPR 0.393
against the wide panel's 0.134, ATP 0.518 against 0.034 -- because nothing is compressed
away when width equals channels. What it cannot reach is what its promoters do not read:
iron, cell wall, copper, sulfur and the rest come back near zero, which is the rank ceiling
doing exactly what it says.

Loaded from disk and applied to a plate it never saw, the deployable model predicts redox
at r = 0.999, proteasome 0.995, ESR 0.984, UPR 0.938 and oxidative 0.812, with an
unstressed control reading near zero.

### The pair design across all 25 stressors

`select_combinations` searched 300 candidate pairs, ranking on a representative ten held
out and verifying the winner on all twenty-five and on fresh noise seeds:

**Tier 1 — self-consistent in simulation**, and on the *fitted*-loadings path, where Ledermann caps four channels at one factor. See the two-paths note above.

| pairs | median transfer, scored on all 25 |
| --- | --- |
| 0, one agent at a time | 0.029 |
| 6 | 0.236 |
| **7** | **0.295** |
| 8 | 0.281 |
| all 300 | 0.121 |

Seven is the peak; the eighth pair lowers the score. Crossing every pair is two and a half
times worse than seven chosen ones at forty times the wells.

That last comparison is conditional on the setting it was derived in -- 25 stressors on 12
channels with six revealed -- and can reverse on a narrower panel, because a design is only
optimal for the readout it was optimised against. What holds everywhere tried is the margin
over dosing one agent at a time. Re-run the search when the panel changes.


## Calibration: what the uploaded plates did to the panel

The panel's dose parameters were literature estimates. Two have measured ladders behind
them, so `scripts/run_calibration.py` fits them and reports whether the encoded value
survived. Two things came out of it, and one is a change to the generator itself.

### The dose response is biphasic, and a saturating Hill could not say so

Corrected activity rises to a peak near 1 mM and then falls -- UPRE1 runs 1.00, 1.15, 1.27,
1.38 and back to 0.91 at 5 mM DTT. A dose high enough to induce strongly is also high
enough to start stopping the cell. `module_response` now multiplies induction by a
viability term, so the curve peaks and turns down instead of approaching a ceiling.

The two agents fail differently and the difference is not cosmetic:

- **At 2 mM peroxide the wells are dead.** Growth is near zero and blank-corrected
  fluorescence goes negative. That is an absent measurement, not a small one, and
  `usable_points` drops it rather than fitting it.
- **Under DTT the cells are alive but slowed, and the naive readout inverts the answer.**
  Raw fluorescence per OD climbs to 2.17-fold at 2 mM while true activity is falling,
  because reporter accumulates in cells that have stopped diluting it. This is the growth
  confound the whole repository is built around, visible on your own plates.

### What the verdicts were

**Measured on real plates — diagnostic, not predictive.**

| agent | encoded | fitted | shape-free half-max | verdict |
| --- | --- | --- | --- | --- |
| DTT | 1.0 mM | 0.11, 0.14 | **0.304, 0.304** | **REFUTED** |
| H2O2 | 0.5 mM | unidentifiable | 0.498, 0.396 | INCONCLUSIVE |

DTT's encoded EC50 was wrong by more than threefold, and both replicate constructs cross
half-maximal at the same 0.304 mM, so the refutation does not rest on one construct or one
plate. The panel now carries the measured value and its source says `measured, not
assumed`. DTT also tolerates a far wider window than the default assumes -- lethal is
around 7 mM against an inducing 0.3, a 24-fold range where the default guesses three.

Peroxide is inconclusive for a reason worth stating: every dose that would have shown the
falling limb killed the cells, so the biphasic fit is unidentifiable and the optimiser ran
to its bounds. Reporting that fit as a refutation would have corrected the panel against
numbers the data never contained, so `fit_dose_response` refuses a parameter pinned at a
bound or a lethal dose below the inducing one. The shape-free crossing still works, and it
**confirms** the literature 0.5 mM. To settle the lethal dose, dose between 1 and 2 mM --
the ladder currently jumps straight from a live culture to a dead one.

### Where a weak fit needed a second opinion

The biphasic fit reaches R2 of only 0.53 on seven DTT points with four free parameters,
which is too little to set a value from alone. `half_maximal_dose` interpolates the
crossing from basal to peak and assumes nothing about the shape; both replicates land on
0.304. The refutation rests on that, not on the fit.

## Measuring the optical linear range

`calib/od_linearity.py` turns a dilution series into `od_linear_max`, the number the
optical gate flags as most worth measuring and still carries as a placeholder. Protocol P4
in `docs/PROTOCOLS.md` has the bench steps. It reports whether saturation was reached at
all, so a series that never left the linear range gives a lower bound rather than a limit.


## Audit: where the earlier numbers were wrong

Everything above was re-examined against the plates rather than against itself. Four things
came out of it, and three were corrections to claims made earlier in this file.

### The calibration verdict was overturned

The first pass reported DTT's EC50 as measured at 0.304 mM and refuted the literature 1.0.
That was wrong three times over:

- the shape-free crossing it rested on is **biased low by up to 87%**, because the peak it
  measures against is itself pulled down by the toxicity that turns the curve over;
- the biphasic fit that disagreed with it was **misspecified**, leaving systematic residuals
  of +12% at the peak and -22% past it, with the lethality exponent pinned at its bound in
  every replicate;
- both estimates were anchored on a peak measured where **growth was already at 71% of
  control**, and for peroxide at 23%.

Restricted to doses where the culture stays healthy, both ladders are still climbing.
Neither EC50 was ever identifiable, both keep their literature values, and
`induction_is_identifiable` now says so rather than fitting one anyway.

What the plates do measure is **lethality**, because growth measures it directly instead of
inferring it from a reporter through a model of transcription. DTT halves growth at 1.45
and 1.65 mM; the reporter route had said 7.2. Peroxide halves at 0.79 and 1.24 mM.

DTT's window turns out to be under twofold -- lethal at 1.55 against an inducing 1.0 -- and
that is why the ladder cannot characterise the promoter. The dose ladder in
`panel_dataset` now spans a fraction of each agent's measured lethal dose instead of a
multiple of its EC50, so no simulated well is dosed past the point where a culture reports.

### The generator had no growth in it

`culture.py` states the standard in its own docstring: growth sits in the reporter's
balance, so a culture whose promoter never moves still shows an apparent dose response, and
data without that term validates the pipeline against a strawman. `panel_experiment.py` was
that strawman -- readings were loadings times activity times a noise draw, with no growth,
no dilution and no maturation anywhere in it. Every transfer and training number was
produced on data lacking the confound this repository exists to handle.

It is derivable rather than invented. Activity is recovered as `k = dR/dt + mu*R`, so an
error in the estimated growth rate propagates as `k_hat = k * (1 + eps/(mu + k_deg))`.
Growth falls with dose, so the same absolute error becomes a larger relative one at every
rung -- a bias with dose structure, which no amount of measurement noise reproduces.

### The noise was three to seven times too low

Measured from your own replicates across 28 matched conditions:

**Measured on real plates — diagnostic, not predictive.**

| | measured | what the design work had used |
| --- | --- | --- |
| activity CV, plate to plate | **0.146** | 0.02-0.05 |
| growth rate CV | **0.33** (IQR 0.19-0.46) | 0.20, assumed |

Both are now named constants in `panel_experiment.py`, and the scripts use them.

### What survived, and what shrank

**Tier 1 — self-consistent in simulation.** Scored within the declared generator configuration, not validated against biology (`docs/research/CIRCULARITY.md` §1). No null or UCBOG accompanied this historical table.

| claim | as first reported | at measured noise and confound |
| --- | --- | --- |
| combinations lift transfer | 8-10x | **2.4x** |
| alignment predicts transfer | r = 0.49 | **r = 0.85** |
| wide model, held-out modules | 0.551 | **0.358** |
| buildable panel, held-out modules | 0.442 | **0.282** |

The design finding survives and the mechanism claim is **stronger** than first measured, but
the headline lift was roughly double what the data supports. One earlier finding reversed
outright: arbitrary pairs scoring below no pairs at all was an artifact of a ladder running
past the lethal dose, where the extra wells bought nothing because nothing was reporting.
Dosed inside the healthy range, arbitrary pairs help -- and the chosen design still helps
more.

### What is still not modelled

- One saturation per stressor drives every module it targets, so a single agent cannot hit
  two regulons at different concentrations. Real dose responses do.
- The panel path has no acceptance gates. The older augment path has three -- physical,
  posterior-predictive and diversity -- and nothing imports it.
  **Closed since:** `generator/panel_gates.py` carries all three on the path every result
  runs through, and the augment module was deleted rather than left as a second
  implementation. See [Gates on the panel path](#gates-on-the-panel-path) and
  [Removed](#removed).
- No plate position effects, no edge effects, no temporal autocorrelation, no maturation
  delay on the transcriptional channels.


## Second audit: three structural gaps closed

### One agent can now reach two regulons at two concentrations

A single saturation drove every module a stressor targets, so peroxide engaged Yap1 and
Msn2/4 at exactly the same concentration and a dose series traced one ray. Two orderings
are defensible enough to encode as defaults: a metabolite pool moves at the concentration
of the chemistry itself, with no threshold to cross, and the general stress response
trails the specific sensor because it answers substantial damage rather than a signal.
`module_ec50` applies them; explicit overrides win where something better is known.

Unflattening the dose axis **improved** the findings rather than costing anything: the
design lift went from 2.4x to 4.0x and the mechanism correlation from 0.85 to 0.87, because
a stressor is no longer a single direction however finely it is dosed.

### Every promoter now sits on a constitutive floor

Readings rose from zero, so a fully induced module looked like an infinite fold change. The
calibration fit puts UPRE's basal at 1847 against an induced 713 at peak, and the measured
induction across the whole DTT ladder is 1.38 to 1.47 fold. Starting from zero inflated the
contrast by roughly an order of magnitude.

The posterior-predictive gate is what caught it -- simulated spread came out at 1.5 against
0.17 measured, and no amount of added noise closes a gap that is a missing floor.

### The noise model was double-counting, and the first fix for it was a bandage

Passing the measured activity CV and the measured growth CV together was wrong: the
observed activity CV **already contains** the growth-correction error, because it is
measured on corrected activity. Adding the growth term on top produced ten times the
dynamic range of a real plate.

The first repair was worse than the problem. Two free parameters were grid-searched until
a simulated plate reproduced one measured summary statistic -- which is a ridge of equally
good answers, not an identification, and the point taken off it carried no information.
Trying to recover the split from the plates directly fails outright: with two plates per
condition each CV is a one-degree-of-freedom estimate, and regressing CV squared on
1/mu squared returns an R2 of **0.07**.

**Neither term needed fitting.** A growth rate is the slope of log optical density against
time, so the error in it is the standard error of that slope -- residual spread over the
spread of the time points and the root of their number, every term measurable on the trace
itself. `growth_rate_uncertainty` computes it, and across 210 real wells it is

    sigma_mu = 0.0117 1/h    (IQR 0.0079 to 0.0200, median window 25 points)

That is an **absolute rate error**, which was the parameterisation error underneath
everything else. The old parameter was a CV multiplied by the control's growth rate, worth
about 0.048 -- four times too large. Expressed properly the dose structure stops being
something to dial in: the same 0.0117 is 3% of a healthy growth rate and 12% of one slowed
to 0.10, and it becomes so on its own.

That leaves one quantity bounded rather than two searched. The residual measurement CV is
set at 0.14 under the observed total of 0.146, the growth term is measured independently,
and what the pair produce is then **checked** -- 0.133 against 0.146, with shape features of
[0.182, 1.373, 0.458] against [0.173, 1.308, 0.523] measured. The gate accepts it.

The same measurement grounds the Fisher penalty in `sensor_selection.py`, which had carried
0.20 as a stated assumption. A promoter fusion pays the relative growth error it inherits;
a ratiometric sensor pays none of it. That is `growth_penalty_cv()` -- the measured absolute
error over `REFERENCE_GROWTH_RATE`, the 0.22 1/h a stressed culture actually runs at, which
comes to 0.053.

It is **computed, not written down**. A derived number frozen as a literal agrees with its
inputs only until one of them is remeasured, which is the grid-search failure one step
removed, so the quotient is taken in code and tested against its two sources.

The recommended build is unchanged under a penalty nearly four times smaller, which is worth
knowing: that conclusion never depended on the assumed value being large -- and that could
not have been claimed while the number was invented.

### Plate geometry

Edge wells evaporate, the medium concentrates, and the culture behaves as though dosed
higher than the pipette says. It cannot be folded into the noise term: measurement noise is
independent per well and averages away over replicates, while a position effect belongs to
the position and every replicate in that ring inherits it. A generator with noise but no
geometry overstates what replication buys. `plate_layout.py` assigns wells deterministically,
so a design that puts every replicate of one treatment down a single column can be seen to
have confounded treatment with position.

## Gates on the panel path

`panel_gates.py` carries the three the augmentation path already had, and
`scripts/run_training.py` now refuses to train on data that fails them.

**physical** -- a reading below the constitutive floor, past the reader's ceiling, or a
plate with no induction anywhere.

**posterior-predictive** -- the one that earns its keep, because it is the only check that
compares the generator against measured data rather than against its own assumptions. It
runs on a matched slice: the reference is one promoter under one agent, so holding a
twenty-channel panel to that spread would refuse it for covering more biology.

**diversity** -- a batch of near-copies, which passes everything else while carrying the
information of a single plate.

## Where the numbers stand after both audits

**Tier 1 — self-consistent in simulation.** Scored within the declared generator configuration, not validated against biology (`docs/research/CIRCULARITY.md` §1). No null or UCBOG accompanied this historical table.

| claim | first reported | after the audits |
| --- | --- | --- |
| combinations lift transfer | 8-10x | **11x** (0.038 to 0.429) |
| all pairs versus a chosen few | 2.4x worse | **2.7x worse** (0.161 versus 0.429) |
| wide model, held-out modules | 0.551 | **0.247** |
| buildable panel, held-out modules | 0.442 | **0.141** |

The design findings held or strengthened; the absolute recovery numbers roughly halved once
the generator stopped flattering itself. What did not survive is the claim that combinations
rescue a stressor with a private module. They cannot: holding out BPS withdraws every
treatment containing BPS, so its pairs leave with it and the iron axis is as absent as
before. The fix for a private module is a second agent that reaches it -- cobalt chloride
perturbs iron handling, and with it in the panel the axis survives BPS being withdrawn.


## A better route to the same quantity: mu cancels

`reporter.promoter_activity` inverts the per-cell concentration: smooth `R = F/X`,
differentiate it, add `(mu + k_deg) R`. Writing the same equation in terms of total
signal removes the growth rate entirely. From `R = F/X`,

    dR/dt = (dF/dt)/X - R mu

so `k_synth = dR/dt + (mu + k_deg) R` collapses to

    k_synth = (dF/dt)/X + k_deg (F/X)

and `mu` has cancelled exactly, for any `k_deg`. This is the form the wider
genetic-circuit characterisation literature uses; the per-cell route is the one this
project derived independently.

Two things follow. It takes one numerical derivative rather than two -- differentiating
a raw total instead of a quotient of two noisy channels. And more importantly it does
not use `mu` at all, so the growth-rate estimate's error, itself a smoothed derivative
of a noisy optical trace, is removed rather than propagated.

On simulated data with a known constant activity and a decelerating culture at 2% reader
noise, RMSE against truth falls from **8.0% to 5.7%** -- a 28% reduction.

On the 168 real NewProtocol wells:

| | per-cell route | total-signal route |
| --- | ---: | ---: |
| late-window activity, median | 1637.3 | 1630.6 |
| step-to-step roughness, median | 594.9 | **358.6** |
| wells with negative late activity | 11% | 11% |

The two agree on the answer, which is the check that the algebra is right rather than
merely different. The total route is **40% smoother**. And the negative-activity fraction
is unchanged, which matters: that signal is what refutes the stable-reporter assumption
on the July plates, and a reformulation that quietly smoothed it away would be hiding the
finding rather than improving the estimate.

`promoter_activity` is kept, not replaced. It is the only route that handles a maturation
step, because cancelling `mu` also loses the immature pool -- so the total-signal form
refuses a maturing reporter rather than returning the mature-only answer.

## Fourth audit: the gate was silently skipping the plates it was written for

`scripts/run_gates.py` resolved its channels by exact name -- `"OD600"`, `"mCitrine"` --
and returned `None` without a message when neither matched. A plate read more than once
names its blocks `OD600[1]`, `mCitrine[2]` and so on, which is every NewProtocol export.
So the two replicates carrying a usable blank were dropped in silence, and the
improvement this gate exists to demonstrate had **no committed table behind it** while
being quoted throughout this document. Both failures are now fixed: channels resolve
through `raw_channel`, which handles the suffix, and an unusable plate returns a note
rather than nothing, because an absent plate and a plate that passed nothing look
identical in a report that says neither.

With the plates actually gated, three of the four documented D2 findings reproduce
exactly, and one G1 figure does not.

| Plate | G1 pass | Negative-activity fraction | Implied min k_deg |
| --- | --- | ---: | ---: |
| 2026-07-01 | 53/90 (59%) | 0.72 | 0.0490 /h |
| 2026-07-09 | 8/42 (19%) | 0.52 | 0.0742 /h |
| 2026-07-22 | **67/84 (80%)** | **0.00** | **0.0000 /h** |
| 2026-08-03 | **52/86 (60%)** | -- | -- |

The July plates land inside the documented 52-72% negative-activity range and the
0.049-0.074 /h implied loss rate. The 2026-07-22 replicate reproduces its documented 80%
pass rate exactly and shows zero negative activity with an implied loss rate of exactly
zero -- which is the central D2 claim, that the apparent excess loss was a property of
the old plates and not of the reporter, now standing on a committed artefact.

**The 2026-08-03 pass rate does not reproduce.** This document reports 86%; the gate
returns **60%** on the 86 wells shared between that plate's two 87-well reads, which is
the correct pairing. The 21-well subset cannot be the source either -- it carries no
blank, so the gate cannot run on it at all. The discrepancy is unresolved and recorded
rather than reconciled; the earlier figure should not be quoted until someone establishes
what well set produced it.

## Third audit: labelling what each number is

The first two audits re-computed numbers. This one re-computed almost nothing and asked a
different question: **for each claim in this file, what kind of evidence stands behind it?**
The answer is `docs/CLAIM_BOUNDARY.md`, and going through the exercise changed six things
that had nothing to do with arithmetic.

### The document was asserting things about the code that were false

`scripts/audit_claims.py` now compares the prose against the tree on every claim a machine
can check. When it was first run it failed six ways:

| what the README said | what was true |
| --- | --- |
| a hundred and forty-nine tests, 97% line coverage | **1268 collected and passing**, and coverage regenerated at **81%** — the 97% was wrong too, and by more than the test count was |
| generator/augment.py documented in the contents table | deleted; the orphaned `.pyc` is what gave it away |
| bridge/directionality.py and generator/combination_design.py referenced | both deleted |
| the D1 quickstart command pointed at scripts/run_d1.py | the script is parked, at `scripts/parked/run_d1.py` |

None of these is a scientific error and that is exactly why they mattered. A reader who
spot-checks one number and finds it wrong has no reason to trust the next one, which is
expensive in a document where most of the numbers are careful. The audit runs in CI-shaped
form now: `python3 scripts/audit_claims.py`, non-zero exit on any failure.

### Four scientific claims changed, and each is corrected where it is made

| what changed | where it landed |
| --- | --- |
| **The dataset is smaller than this file said** — not four biological replicates in fluorescence but **two**, for three of the four constructs. All eight dose-response readings are INCONCLUSIVE at that replication, and `AlteredYap1 @ 1.0 mM` moves 1.43 -> **1.31** under a plate-paired estimator. Point estimates and mechanism unchanged; the evidence that they differ from 1.0 is what is missing. | [Data](#data), [dose response](#findings-that-constrain-the-plan) |
| **One measurement did not fail to conclude, it failed.** The ER anchor's RT- margins put genomic DNA at 91.4% and 66.0% of the +RT signal, past the ceiling of the only validated correction. Nothing was learned about the ER biosensors — a protocol finding, reported as one. | [G4](#findings-that-constrain-the-plan) |
| **Two identifiability questions were being reported as one.** The `2 -> 2 of 7` table is a known-loadings rank statement; the transfer and training results go through a *fitted*-loadings path where Ledermann caps four channels at **one** factor. No line of code connected them. | [the two-paths note](#the-whole-stress-landscape-and-what-five-channels-can-reach) |
| **The strongest correlational claim now has a null under it.** Observed median R2 0.522 against a phase-randomised null of 0.117 over 146 wells: it survives, which strengthens it — but ~0.12 of any variance share is autocorrelation artefact and 43% of wells do not individually clear their own null. | [D2](#findings-that-constrain-the-plan) |

### And nothing here is Tier 3

The historical transfer, power and recovery tables are simulation-scoped. The earlier
claim here that every such number was *provably* optimistic over real biology overstated
the theorem: it concerns the expectation of an optimized finite-sample objective under an
assumed distribution, not each independent held-out score. Family checks and a limited
configuration-choice UCBOG have since been run (`docs/CROSS_FAMILY.md`); neither validates
the shared biological assumptions or retrospectively bounds every table. See
`docs/research/CIRCULARITY.md` §1 and §4 for the corrected scope.

The one scored forward prediction on real data exists and it is a **loss**
([above](#on-real-traces-asked-to-forecast-the-twin-loses)). Registering the next one before
running it is the only route to Tier 3, and it costs no bench time.

### What this audit could not resolve

- **Orphaned bytecode remains in `src/ystwin/**/__pycache__`** for all three deleted modules.
  It is harmless at runtime and the check that reports it is doing its job; clearing the
  caches clears it.
- **Coverage is not checked by anything.** It was regenerated for this pass, but
  `audit_claims.py` compares only the test *count* against the tree, so the coverage figure
  can drift silently the way the old 97% did. Re-run `pytest --cov=ystwin` before quoting it.
- **The autofluorescence constant is still unmeasured**, and three of the eight dose-response
  rows turn on it.
- **No UCBOG accompanied these historical tables.** A later configuration-choice estimate
  appears in `docs/CROSS_FAMILY.md`; it does not bound each table's biological reality gap.

## Training across culture conditions

Every simulated well was exponential-phase glucose at 30 degrees. Two of the twenty-four
modules are set by the baseline before any agent is added -- the carbon regulon is
derepressed the moment glucose runs out, and the general stress response is high in
stationary phase whatever else is happening. A model trained in that one corner has never
seen a raised ESR without a stressor behind it, which is precisely the reading a plate gives
at the end of a run.

`generator/context.py` carries the baseline: carbon source, growth phase, temperature,
aeration and strain. Each shifts what the stressor is added to and how fast the culture is
diluting its reporter while being measured -- not what the stressor does.

Asked to report the ESR of an unstressed culture read at three points in a run:

**Tier 1 — self-consistent in simulation.** Scored within the declared generator configuration, not validated against biology (`docs/research/CIRCULARITY.md` §1). No null or UCBOG accompanied this historical table.

| read at | true ESR | trained in one corner | trained across seven contexts |
| --- | --- | --- | --- |
| exponential | 0.00 | 0.06 | 0.21 |
| diauxic | 0.45 | 0.35 | 0.11 |
| stationary | 0.75 | **0.63** | **0.32** |

The narrow model reports 0.12 against a true 0.75. It is not slightly wrong; it cannot
represent the case at all. Breadth halves the error on the conditions that matter and costs
accuracy on the one corner it used to have to itself, which is the right trade for a model
that will meet plates read at whatever time the run ended.

### A well too slow to correct is dropped, not fitted

Activity is recovered by dividing growth out, so a stationary culture divides by nearly
nothing and a growth estimate landing slightly negative sends the recovered activity below
zero -- an arithmetic artifact, not a weak promoter, and exactly what the real plates showed
at 2 mM peroxide. The generator now drops those wells and reports how many, as an
experimenter would. Keeping them would teach a model that promoters run backwards in
stationary phase. The physical gate caught this the first time contexts were switched on.


## The time axis

Every simulated reading was a steady state; a real plate is read on a timetable. Reporter
accumulates as `dR/dt = k - (mu + k_deg) R`, so it approaches `k/(mu + k_deg)` with a time
constant of `1/(mu + k_deg)`, and reading before that reads a fraction of the response.

The bias runs the wrong way, which is why it has to be in the generator rather than noted in
a caveat. **The relaxation rate is the growth rate**, so the stressed culture -- slow-growing,
the one the experiment is about -- takes longest to arrive. An early read understates stress
most in exactly the wells where stress is greatest.

**Tier 1 — arithmetic, given the model.** An exact consequence of `dR/dt = k - (mu + k_deg) R`, not a measurement.

| culture | mu | hours to 90% of steady state |
| --- | --- | --- |
| unstressed exponential | 0.40 | 5.8 |
| mild stress | 0.30 | 7.7 |
| the dose that halves growth | 0.20 | 11.5 |
| heavily stressed | 0.10 | 23.0 |
| stationary | 0.02 | 115 |

A ratiometric sensor carries none of this. It reports an equilibrium between two forms of one
molecule and is at its final value whenever the plate is read -- a second, independent reason
to carry one, and only visible once time exists in the model.

### What this says about the uploaded runs

Every plate ran **4.1 hours**. The fastest-growing control needs 5.5, and NativeYap1 at
1 mM DTT needs 23.7. **No well on any plate reached steady state**, and the shortfall is
differential: across the DTT ladder a 2 h read captures 55% of the response at the bottom
dose and 48% at the top, so the measured dose response is compressed by roughly a tenth at
its most informative end.

The corrected activity estimator does not formally need steady state -- it recovers
`k = dR/dt + mu*R` at any time -- so this does not invalidate the activity numbers. What it
does mean is that the naive fluorescence-per-OD readout on these plates is uninterpretable
on its own, and that the high-dose wells, which are the ones the ladder exists to resolve,
had the least time to say anything.

### A longer run would not fix the noise, and it is worth knowing why

The obvious recommendation is a longer run, and on the growth-rate uncertainty it is right:
the standard error of a log-OD slope falls as roughly the window length to the three-halves,
so doubling the run improves it about threefold. But that term is not what limits the
measurement:

    residual        0.140
    growth term     0.053     (sigma_mu 0.0117 / mu 0.22)
    predicted total 0.150     against a measured 0.146

The growth-correction error is **12.6% of the variance**. Everything else is residual
well-to-well and plate-to-plate variation, which a longer run does nothing about. Run longer
to reach steady state and to give slow wells time to respond -- not to reduce noise.

That the two independently obtained terms combine to within 3% of the measured total is the
strongest check the noise model has: `sigma_mu` came from the optical density traces and the
residual from replicate spread, and neither was fitted to the total they reproduce.


## The metabolic link, and what the literature would not support

Seven literature sweeps were run before any of this was encoded. Three findings changed the
plan, and two of them killed the coupling that was about to be built.

### A sensor reads a pool, and a pool does not set a flux

- **ATP.** Larsson 1997 (PMID 9393686), yeast chemostats: glycolytic flux correlates
  **negatively** with intracellular ATP, and **not at all** with the ATP/ADP ratio. Mapping
  a measured ATP level onto a flux has the opposite sign to intuition, and the ratio carries
  no signal whatsoever.
- **NADH.** Vemuri 2007 (PMID 17287356) drained cytosolic NADH hard enough to abolish 80% of
  glycerol production and the critical dilution rate **did not move** (0.27 +/- 0.02 against
  0.29 +/- 0.01, not significant); only the mitochondrial sink moved it. Agrimi 2011
  (PMID 21335394) perturbed the mitochondrial NAD pool, shifted D_crit, and left the
  critical **glucose flux** invariant, concluding that q_glucose is the governing variable.
  Su 2015 (PMID 26098102) cut the free ratio by 78% and moved ethanol by 7.9%.
- There is **no published quantitative relation** between cytosolic NADH:NAD and ethanol
  yield, respiratory quotient, q_glucose or dilution rate in yeast. The free cytosolic ratio
  has been measured at exactly two points, ever, by one lab.

Cytosolic redox is a **consequence** of the flux, not its control variable. A twin that
drives fermentation from a cytosolic redox reading would be encoding a correlate.

### What does have precedent, and it is narrower

Kummel 2006 (PMID 16788595) placed the adenylate energy charge and the NADH/NAD ratio on
**iND750**, and Martinez 2014 (PMID 25028891) did the same on **Yeast 5** -- both yeast
genome-scale models, both as bounds on Gibbs energies constraining reaction
**directionality**, which way a reaction may run rather than how fast. Both papers set those
bounds deliberately loose so they would not bind.

So `bridge/thermodynamic.py` supplies transformed formation energies and `dG' = dG'0 + RT ln Q`,
and nothing above it claims a measured pool sets a rate.

### Nobody has put a *metabolite* biosensor on a GSMM, and that is not entirely encouraging

Across 35 distinct query formulations, **no paper constrains a genome-scale metabolic model
with a genetically encoded metabolite biosensor, in any organism**. The nearest prior art is
Takhaveev 2023 (reporter-derived rates as biomass stoichiometry in a non-genome-scale
model), Monteiro 2019 (the inversion -- a constraint-based model calibrating the biosensor),
Labhsetwar 2013/2017 (fluorescent-protein **abundance** as enzyme-capacity bounds) and
Ishchuk 2022 (a biosensor validating GSMM predictions). The niche is open, but two of the
three mechanisms that would justify filling it are the ones refuted above.

**This is a claim about metabolite biosensors and it does not extend to the stress state.**
A latent stress state read off a yeast *promoter* panel has substantial prior art —
Granados 2018 decoded ten of them in this organism — and the honest version of that
claim is [What is new here, and what is not](#what-is-new-here-and-what-is-not).

### Two things the tables got wrong until they were checked

**The energies need transforming.** Untransformed, the tabulated formation energies give
**-6.8 kJ/mol** for ATP hydrolysis against a textbook -30. Transformed to pH 7.5 and
0.25 M ionic strength they give **-39.5**, which is what the textbook figure becomes once
the magnesium it assumes is removed. Magnesium binding is not modelled, and it matters most
for exactly the adenylates.

**pH is a parameter, not a constant.** Much of the 40 mV disagreement between the two
published E_GSH literatures for the yeast cytosol -- Kojer 2012 at -306 mV against Ayer 2013
at -350 mV -- comes from one line assuming pH 7.0 where the other measured 7.5.

### A sensor in the panel that has never worked in this organism

**Peredox has no published application in *S. cerevisiae*.** Neither has SoNar. The only
genetically encoded NADH readout in yeast is a GPD2-promoter fusion (Knudsen 2014,
PMID 25401080), which reports a transcriptional proxy and cannot give a ratio at all.
`Reporter.demonstrated_in_yeast` now records this, and `demonstrated_in_yeast()` filters the
catalogue down to what a wet lab could start on today.

The restriction turns out to cost **nothing**: from four channels to nine the selector picks
the same set either way. So no conclusion rested on it -- the guard is against a future
mistake rather than a correction to a past one. roGFP2-Grx1, which the recommended build
does use, is thoroughly established in yeast.


## Thermodynamic flux analysis

`bridge/tmfa.py`. Two bugs were found by auditing this before trusting it, and both
invalidated earlier conclusions.

### A herbicide's formation energy was being used for ribose-5-phosphate

The ModelSEED alias tables map one identifier to several compounds and give no way to
choose between them. Taking whichever came last mapped **ribose-5-phosphate** (via KEGG
C03736) onto `cpd19028`, **Chloramben** -- a chlorinated benzoic acid herbicide.

The consequence was not subtle. Ribose-5-phosphate isomerase came out at **+75.4 kJ/mol**,
an isomerase between two pentose phosphates whose real value is near +0.5, and the database
reported an error of only 1.7 on it. Being confidently wrong, it made an essential reaction
thermodynamically impossible and the whole pentose phosphate pathway blocked growth.

Every mapping is now checked against the metabolite's chemical formula, ignoring hydrogen so
protonation state does not matter. **88 mappings were rejected**, coverage went 56.1% to
54.7%, and growth stopped collapsing.

### One binary per reaction over-constrains; it needs one per direction

A single direction binary forces every constrained reaction to carry a signed dG' whether or
not it is running. Reactions at zero flux still had to pick a side, which squeezed the
concentration variables hard enough that the model needed a window spanning thirteen orders
of magnitude before it would grow at all.

Two binaries with `zf + zr <= 1` leave an inactive reaction's dG' free, which is what Henry
2007 does.

### With both fixed, the answer to the original question is still no

Growth is preserved exactly -- 0.08584 plain, 0.08584 constrained -- across glycolysis, the
TCA cycle, pyruvate metabolism, the pentose phosphate pathway, oxidative phosphorylation and
mitochondrial transport. On 86 reactions sharing 132 concentration variables, applying the
measured ATP/ADP and NADH/NAD windows narrows **no** flux range.

The reason is arithmetic, and it survives the audit. Every metabolite in the concentration
quotient gets a vote weighted by its range: the generic window is 7.6 ln-units, about
19 kJ/mol per metabolite, where a pinned cofactor ratio contributes around 5.5.

**Tier 1 — deterministic computation on an adopted model.** Nothing is fitted here, so no simulation-optimisation bias; it is only as good as the GSMM and the thermodynamic tables.

| reaction | swing from the measured couple | swing from the unmeasured substrates |
| --- | --- | --- |
| malate dehydrogenase, cytosolic | 5.5 kJ/mol | 57.5 |
| glycerol-3-phosphate dehydrogenase | 5.5 | 57.5 |
| pyruvate carboxylase | 5.5 | 95.8 |

The unmeasured substrates outvote the measured couple seven to seventeen fold and move to
cancel whatever it says. **A single biosensor is the wrong shape of measurement for this
job** -- Kummel 2006 constrained iND750 with a measured metabolome, dozens of species at
once. That is a statement about the method rather than about the sensor, and it says what
would work: broad metabolomics, not one couple.

### Both missing physical terms are now in

**Magnesium.** ATP, ADP and phosphate bind Mg2+, and at the roughly 1 mM free Mg2+ of a
yeast cytosol most ATP is MgATP. Binding stabilises the more highly charged species
preferentially, so it costs hydrolysis some drive. Association constants from Alberty
(log K 4.19, 3.17, 1.88 at I = 0.25 M) move ATP hydrolysis from **-41.9 to -37.3 kJ/mol** --
the +4.6 shift those constants predict. The remaining gap to the textbook -30.5 is the
SEED table's own offset, not a missing term.

**Membrane potential.** A charge crossing a membrane does `zF*dPsi` whatever the chemistry
does. `COMPARTMENTS` carries a pH and a potential for each, with the cytosol as the
electrical reference; the matrix sits at -160 mV, so **a proton entering it releases
-15.4 kJ/mol** that a chemistry-only energy cannot see.

### Which exposed a real limitation of the model's transport reactions

Applying the electrical term everywhere broke mitochondrial transport outright. The cause is
not the physics: yeast-GEM writes many carriers as a bare uniport of a charged species, so
charging the full `zF*dPsi` gives **+77 kJ/mol for PRPP transport** and **+124 for
propionyl-CoA**, which no carrier does. The real mechanisms are symports and antiports whose
counter-ion the stoichiometry omits.

The term is therefore applied where the net charge moved is one or less -- what a real
electrogenic carrier does, and what the ADP/ATP exchanger genuinely is, swapping ATP(4-) out
for ADP(3-) in. Above that the counter-ion is missing from the equation and the energy
**cannot be estimated from the stoichiometry given**, so it is refused rather than
fabricated. That is what Henry 2007 does with any unreliable energy, and it is the
difference between declining to answer and answering wrongly.

Growth is preserved at 0.08584 across glycolysis, the TCA cycle, pyruvate metabolism, the
pentose phosphate pathway, oxidative phosphorylation and mitochondrial transport.

### The assumed concentration window excludes a measured value

Free cytosolic NADH in yeast runs **3 to 30 uM** (Canelas 2008, PMID 18383140) against the
**10 uM floor** Henry 2007 assumes, so most of the measured range sits underneath it. A
measurement now overrides the window rather than being clipped to it: the window exists
because concentrations are unknown, and where one is known the assumption has no standing.

## Removed

Three modules are gone from the tree. They are named here without backticks, because
`scripts/audit_claims.py` treats a backticked module name as an assertion that the file
exists — which is the check that caught the first two still being documented as live.

**generator/augment.py** — S2 augmentation with physical, posterior-predictive and
diversity gates. It duplicated gates that `generator/panel_gates.py` now applies on the path
every result runs through, and nothing imported it. The findings are kept as tests.

**bridge/directionality.py** — produced an infeasible model and was superseded by
`bridge/tmfa.py`, whose audit found the two bugs recorded above. The findings are kept as
tests.

**generator/combination_design.py** — combination selection, now
`analysis/experiment_design.py`. Its results (the pair-design tables above) survived the
move; only the module did not.

All three left orphaned bytecode behind in `src/ystwin/**/__pycache__`, which is how the
first two were found still sitting in the contents table. `audit_claims.py` now checks for
that directly: a `.pyc` with no matching source is a deleted module, and a deleted module is
a documentation claim about to go stale. Clearing the caches (`find src -name __pycache__
-exec rm -rf {} +`) clears the check.


## Validating against published physiology

`data/physiology/` and `tests/test_physiology_validation.py`. The chemostat series shipped
with yeast-GEM is the data its own growth-associated maintenance was **fitted to**, so
reproducing it to 5% confirms the pipeline runs and validates nothing. These are
measurements the model never saw.

### The maintenance prediction holds, and it is genuinely held out

**Tier 2 — held out.** The one place in this file where a number agrees with a measurement
the model never saw. **Still not Tier 3**, for two reasons worth stating: the measurement
predates this project, so nothing was registered in advance; and the parameter being
vindicated is yeast-GEM's, not this repository's.

Vos 2016 (PMID 27317316) ran an aerobic retentostat at near-zero growth and measured glucose
consumption of **0.039 +/- 0.003 mmol/gDW/h**. The model predicts **0.036** -- inside the
measurement error, from a maintenance parameter set before that experiment was published.

### The model does not ferment at all

Nothing in yeast-GEM limits respiration, so glucose is respired however fast it arrives:

**Tier 1 — deterministic computation on an adopted model.** Nothing is fitted here, so no simulation-optimisation bias; it is only as good as the GSMM and the thermodynamic tables.

| glucose uptake | growth | ethanol | oxygen |
| --- | --- | --- | --- |
| 4.0 | 0.353 | 0.00 | 8.9 |
| 8.0 | 0.709 | 0.00 | 17.6 |
| 12.0 | 1.066 | 0.00 | 26.4 |

**There is no Crabtree effect**, and the model asks for more than double the highest oxygen
uptake ever measured in this organism.

### A measured ceiling restores it

Postma 1989 (PMID 2566299) puts maximum oxygen uptake at **12 mmol O2/gDW/h** in CBS 8066.
Imposing it:

**Tier 1 — deterministic computation on an adopted model.** Nothing is fitted here, so no simulation-optimisation bias; it is only as good as the GSMM and the thermodynamic tables.

| glucose uptake | growth | ethanol | oxygen |
| --- | --- | --- | --- |
| 5.0 | 0.442 | 0.00 | 11.1 |
| 6.0 | 0.492 | 0.98 | 12.0 |
| 8.0 | 0.535 | 4.36 | 12.0 |

Ethanol appears, and the cap binds once it does.

**It is still wrong about where — and that is the other half of Tier 2, the failing half.**
van Hoek 1998 puts the critical growth rate near 0.28 /h;
an oxygen ceiling alone puts it near 0.45. The real limit is proteome allocation rather than
gas exchange, which is what ecYeastGEM models and what this does not. Recorded as a test
that asserts the discrepancy rather than hiding it.


## Kaggle, and what it was and was not good for

Credentials were already working -- the earlier failure was `kagglehub.whoami()`, which in
1.0.2 authenticates by a different path than the classic API key.

**It is thin for yeast physiology, and one result shows why "reputable" has to be checked.**
The top hit for "saccharomyces cerevisiae" is `usharengaraju/floodsdamageindia`, a flood
damage dataset carrying a yeast title. The best-voted yeast transcriptomics set turns out to
be genotype comparisons -- itc1, swr1, INO80, rrp6, all chromatin remodellers -- with no
stress conditions in it at all, so it cannot test the stress panel.

**One dataset was worth having.** BioNumbers (`haotiannnnn/bionumbers-dataset`) is the Milo
lab's curated database, 10k+ values each with a citation, assembled by people with no
interest in this model. Every constant below was chosen from the literature while building
the generator and only afterwards compared:

**Tier 0 — asserted, then independently corroborated.** Corroboration raises confidence in a
constant. It does not make it a measurement of *these* cultures.

| quantity | the generator | BioNumbers |
| --- | --- | --- |
| cytosolic free NAD/NADH | 101-320 | **101(+-14) to 320(+-45)** (BNID 108145) |
| ATP concentration | 0.9-4.4 mM, mean 1.93 | **1.9 mM** (BNID 106020) |
| doubling time, rich medium | mu_max 0.40 /h = 104 min | **~100 min** (BNID 100270) |
| cytosolic GSSG | 4 uM | **4.0 uM** (BNID 103548) |
| cytosolic E_GSH | -300 to -320 mV | **-289 mV** (BNID 103543) |

The last row is the one disagreement, and it is a real split in the published record rather
than an error: BNID 103543 is the rxYFP measurement, the -300 to -320 range is the
Grx1-roGFP2 line, and about 40 mV of the difference is one pH unit in the Nernst term. The
panel records the roGFP2 range and the reason for the gap rather than picking one silently.

`tests/test_constants_against_bionumbers.py` pins all six so they cannot drift, and skips
cleanly if the export is absent.


## The E_GSH disagreement, resolved

It was never a disagreement. The couple is GSSG + 2H+ + 2e- -> 2GSH, one proton per
electron, so its potential moves **60 mV per pH unit**. Kojer 2012 reports -306 mV assuming
pH 7.0; Ayer 2013 reports -350 measured alongside a pHluorin that read 7.5. Half a pH unit
is 30 mV, which is most of the 44 mV between them. `generator/redox.py` computes E_GSH from
oxidation degree and pH, and the two literatures land on top of each other.

### Which exposed a flaw in the recommended build

Ayer measured roughly **one pH unit of cytosolic acidification under H2O2** -- 60 mV of
apparent shift -- against a genuine E_GSH response to 1 mM H2O2 of **40 to 50 mV**. The
artifact is larger than the signal and has the same sign, and their own text says the pH
change "contributed significantly to the observed changes in redox state".

The build carried roGFP2 with no way to measure pH, and **could not have added one**:
pHluorin is also a green excitation-ratio sensor and shares the slot. Under any acidifying
stressor the redox channel would have returned a number nobody could attribute.

`interpretable()` now refuses a build carrying a pH-sensitive sensor without a pH reading,
and the selector resolved it by swapping roGFP2 for **HyPer7**, which Pak 2020 engineered to
be pH-insensitive, occupies the same slot, and reads peroxide instead of glutathione. The
build is now **STRE-general, UPRE-ER, PACE-proteasome, HyPer7** plus the mOrange2 reference.

## Better thermodynamic data

`bridge/equilibrator.py`. The pytfa table is group contribution, which fails where substrate
and product are structurally similar because the errors do not cancel. eQuilibrator is
component contribution from the group that curates BioNumbers, handles pH, ionic strength
and magnesium natively, and reports an uncertainty per reaction.

**Tier 1 — deterministic computation on an adopted model.** Nothing is fitted here, so no simulation-optimisation bias; it is only as good as the GSMM and the thermodynamic tables.

| reaction | group contribution | component contribution | truth |
| --- | --- | --- | --- |
| ribose-5-phosphate isomerase | **+75.4** | **+2.1 +/- 0.8** | near +0.5 |
| ATP hydrolysis | -37.3 (with hand-added Alberty constants) | **-30.5 +/- 0.3** | **-30.5** |
| phosphopentomutase | **+77.4** | -11.5 +/- 1.5 | a mutase, so near zero |

It gets ATP hydrolysis exactly right without the magnesium correction having to be bolted
on, and covers 61% of the model against 54.7%.

### A bug the swap introduced and a test caught

The loader caches one underlying object, so two instances built at different magnesium
silently returned the same energies -- the second load overwrote the first's settings and
nothing errored. Conditions are now applied per query, with two tests holding them apart.
