# The mechanistic layer

`src/ystwin/mech/` and the module that joins it to everything else, `mech/chain.py`.

Read this before reading the code. It says what the layer measures, what it asserts, what it
refuses, and what would change each — and it states the results that argue against the layer
as plainly as the ones that argue for it, because three of the four headline results are
negatives.

The legacy chain results below were produced by running this repository and are checked by
`tests/test_mech_chain.py`. The separate published-HOG extension at the end has its own
source, observation and numerical verification contracts.

---

## 0. Which entry point runs which physics

**Read this before running anything here.** The module you import decides which physics you
get. The routes below do not share a state space, they are not scored against the same
thing, and only the ones *below* the shared engine produced any number this repository
publishes.

| Entry point | Physics it runs | What it is validated against |
|---|---|---|
| `predict.py::predict_protocol`, and `simulate_protocol` re-exported unchanged from `mech/adapters.py` by `generator/culture.py`, `fba/dynamic.py` and `hybrid.py` | **the shared engine**, `mech/engine.py` — one stiff clock, unit-bearing amounts, carbon/nitrogen ledgers, an explicit `validity` | **nothing.** No dataset here has been scored against it, and no script under `scripts/` calls it at all |
| `predict.py::predict_product` | the empirical route: a fitted flux law times relative entry-enzyme expression, then the closed-form steady state in `pathway/solve.py`, with `generator/context.py` bounding mu | the Elizondo 2025 chemostat states, leave-one-strain-out (`pathway/calibrations.py::BETA_CAROTENE_FLUX`) |
| `predict.py::predict_product(mech=...)` | the same, plus this document's chain (`mech/chain.py`) as a BOUND on mu | nothing further. The chain reports bands, and §5.1's assembled gate refuses a caller who asks to *fit* it |
| `generator/culture.py::simulate_empirical_comparison` (legacy name `simulate_culture`) | the synthetic generator's own three equations; the dose enters as two constants (§2(2)) | nothing — it is the model the synthetic cultures come *from*, not a fit to them |
| `fba/dynamic.py::simulate_fixed_volume_comparison` (legacy name `simulate_batch`) | fixed-volume dynamic FBA, integrating rates it is handed | nothing. "It has never been scored" is that module's own docstring |
| `hybrid.py::simulate_hybrid_comparison` (legacy name `simulate_hybrid_case`) | a calibrated sensor estimate reaching a GSMM allocation solve, coupled only by the empirical growth-retention factor | nothing. Its summary reports `biological_validation` as `False` and its own `prediction_kind` as an allocation-conditioned scenario, "not validated realized production" |
| `in_silico.py::simulate_controlled_comparison` (legacy name `simulate_controlled_product`) | an explicit four-control GSMM teacher over that same fixed-volume integration | itself. Same-teacher recovery inside a synthetic world; the summary's `truth_scope` says it is not measured biological production |

The legacy name in each row is retained and still works. It is an alias for the comparison
route in its own row and never for the shared engine, and each module says so at the alias.

**The shared engine is the documented physical entry, and it is not validated.** Those two
statements are not in tension. `EngineParameters` refuses to build at all without an
explicit `allow_prior=True`, and the refusal says why: the connected model carries
unvalidated structural priors even where its rate constants are measured. What the engine
earns by being the documented entry is that there is now ONE implementation of the physics
rather than several that could quietly disagree — `tests/test_engine_adapters.py` pins that
each protocol entry above is exactly one call into `mech/engine.py` with the caller's own
objects, rebuilding nothing. It earns nothing whatever about being right.

```python
from ystwin.mech import engine
from ystwin.mech.contracts import Control, Genotype, Protocol
from ystwin.predict import predict_protocol

parameters = engine.EngineParameters.prior()      # ScientificRefusal without allow_prior
genotype = Genotype.prior()
state = engine.initialize(parameters, genotype, volume_l=0.01, biomass_gdw_l=0.2,
                          medium_mM={"glucose": 30, "nitrogen": 10, "oxygen": 0.2})
protocol = Protocol((0.0, 0.015, 0.03),
                    (Control(0.0, oxygen_transfer_per_h=20, oxygen_saturation_mM=0.2),))
result = predict_protocol(protocol, genotype, state, parameters)
result.validity.biological_validation   # False, and nothing in this repository moves it
```

**Every leave-one-strain-out figure this package quotes was measured with the mechanism
ABSENT.** `mech=` defaults to `None` on `predict_product` and no caller in `src/` or
`scripts/` passes anything else; nothing under `scripts/` calls the shared engine. So every
committed table in `outputs/` came from a comparison route and not from the engine: the flux
calibration and its leave-one-strain-out error from `predict_product`, the hybrid benchmark
from `simulate_hybrid_comparison`, the closed-loop runs from
`simulate_controlled_comparison`. **The empirical route remains the one the published
numbers came from.** Switching the mechanism on can move them (§2(1)), which is exactly why
it is off by default rather than quietly on, and why no figure in this document may be
re-attributed to the engine.

**What would change that.** A scoring run of `predict_protocol` against a dataset the engine
did not see. There is none here, and the missing measurement is named rather than
approximated: a time-resolved trajectory — biomass, substrate, product and reporter on one
clock — for a strain and medium outside the calibration set. Until one exists, a number out
of the shared engine is a simulation, not a prediction.

---

## 1. What it is

Eight modules make the chain this document is about. Seven landed in Phase 0 and none of
them was connected to anything; the eighth is the join. `src/ystwin/mech/` now holds
**seventeen** code modules beside `__init__.py` — it held twelve when this section was
written, and the five that arrived since are the shared engine and its contracts, adapters,
inference and signalling, which section 0 maps and this section does not describe.
`mech/oxidative.py` and `mech/upr.py` drive a reporter channel through
`generator/panel_experiment.py`'s own mechanism seam, which §3 and §5.4 take account of.
`mech/kinetic_sbml.py` and `mech/hog.py` supply the separate published-model workflow
summarized at the end; they do not silently enter this legacy chain:

| Module | What it computes | States |
|---|---|---|
| `mech/params.py` | The `Param` grade and the `REFUSED` sentinel; the free-scalar gate | — |
| `mech/state.py` | `Window`, `StateVar`, `MechState`; the timescale catalogue | — |
| `mech/integrate.py` | The pinned stiff driver (BDF) and the tau/T reduction audit | — |
| `mech/ph.py` | Weak-acid entry, repaired charge balance, the Pma1 ATP bill and Citrine quench; not a panel promoter-response route | `A_i`, charge-linked `pH_c` (algebraically reducible on the audited windows) |
| `mech/burden.py` | The declared-genotype growth tax and the product pool | `product_fraction` |
| `mech/population.py` | Plasmid segregation and the generation clock | `plasmid_bearing`, `generations` |
| `mech/ablation.py` | The scoring protocol every result here is scored under | — |
| **`mech/chain.py`** | **The join: environment in, existing product and reporter out** | — |

`mech/chain.py` runs them in order and hands the result to code that already existed —
`pathway/solve.py` for the product, `generator/culture.py` and `reporter.py` for
the reporter. It writes no second copy of either.

```
context (carbon, temperature, oxygen, phase)  ->  mu_env             environment   IN_CHAIN
stress panel's own dose curve, where no mechanism applies             stress panel  IN_CHAIN
medium pH + total acid  ->  [AH]_o                                    weak acid     IN_CHAIN
[AH]_o  ->  (A_i, pH_c) integrated over the declared window           pH states     IN_CHAIN
[A-]_i  ->  Pma1 ATP bill  ->  mu_max                                 proton bill   IN_CHAIN
declared cassette copies  ->  mu_max                                  burden        IN_CHAIN
mu = min(setpoint, mu_max)                                            vessel        IN_CHAIN
mu  ->  (F, generations)                                              population    IN_CHAIN
mu  ->  solve_pathway  ->  content;  x mean F  ->  culture content    metabolism    IN_CHAIN
mu, dose  ->  simulate_culture  ->  per-cell reporter                 reporter      IN_CHAIN
pH_c(t)  ->  Citrine quench  ->  what the instrument sees             photophysics  IN_CHAIN
[AH]_o  ->  Yeast9 pfba delta                                         GEM depth     AUDITS
proteome allocation arm; fed-batch ethanol accumulation               allocation,
                                                                      ethanol       REPORTED
```

The layer vocabulary is `predict.py`'s own `LayerState`, and `ChainResult.layer_report()`
prints the same shape as `ProductPrediction.layer_report()`. A layer that ran and provably
could not have changed anything says `INERT`; one that needs an argument says `NOT_RUN` and
names the argument.

---

## 2. The four questions, answered with numbers

### (1) Does the environment now move the product?

**Yes, from one distinct value to six.** Ten environments at a held fed-batch setpoint of
mu_set = 0.2543 /h, beta-carotene, one genotype:

| Condition | incumbent `predict.py` | chain, per bearing cell | chain, per gram of culture |
|---|---|---|---|
| reference | 0.97656355 | 0.97656355 | 0.31700026 |
| ethanol carbon | 0.97656355 | 1.11791951 | 0.58884694 |
| galactose | 0.97656355 | 0.97656355 | 0.31700026 |
| 37 C | 0.97656355 | 0.97656355 | 0.31700026 |
| 20 C | 0.97656355 | 0.99072116 | 0.33486299 |
| pH 4 medium | 0.97656355 | 0.97656355 | 0.31700026 |
| anoxia | 0.97656355 | 0.97656355 | 0.31700026 |
| H2O2 1 mM | **REFUSED** | 1.04627953 | 0.42076401 |
| acetic 175 mM, pH 4 | **REFUSED** | 0.98640259 | 0.32926914 |
| ethanol, 37 C | 0.97656355 | 1.13097195 | 0.62977149 |

All units mg/gDCW. The incumbent's spread across the eight it accepts is **1.0000x**, one
number to eight decimals. The chain's is **1.9867x** across six distinct values, and it
answers the two the incumbent refuses.

The same comparison at mu_set = 0.18 /h is the one the build brief names: ethanol, 37 C,
pH 6, anoxia, H2O2, acetic acid and DTT all returned **1.07098788 mg/gDCW** to eight
decimals. The chain reproduces that number exactly for the reference condition — it must,
because nothing there has a mechanism — and separates ethanol (1.11791951 per cell,
0.58884694 per gram of culture, +25.4%) and ethanol at 37 C (1.13097195 / 0.62977149,
+34.1%) from it.

**The size of this answer depends on the setpoint, and the table above is quoted at the
setpoint where it is largest.** That is a real caveat, not a presentational one, because the
only route the environment has to the product is to make the setpoint UNREACHABLE — so the
number of conditions that separate is just the number whose mu_max falls below mu_set. Both
ends, measured:

| mu_set | incumbent | conditions that MOVE | distinct culture contents | culture spread |
|---|---|---|---|---|
| **0.18** (the brief's own baseline) | 1 value, 0 refused | **2 of 10** | 3 | **1.34x** |
| 0.2543 (the top of the fitted range) | 1 value, 2 refused | **5 of 10** | 6 | **1.99x** |

At mu_set = 0.18 the only conditions that move are the two ethanol contexts; 37 C, 20 C,
galactose, pH 4, anoxia, H2O2 and 175 mM acetic acid all still return the reference number,
because 0.18 /h is below what every one of them can reach. The honest one-line answer to
"does environment now move product" is therefore: **yes, but only by refusing a setpoint, and
at the brief's own setpoint it moves 2 conditions of 10 rather than 5.**

**And here is the honest part: the route that carries this is not the pH block.** It is two
things, in this order of size.

- **The vessel rule.** A feed cannot hold a growth rate the strain cannot reach.
  `predict.py::_growth_rate` checks that when a stressor is present and does not check it at
  all for a context, so an ethanol fed-batch at mu_set = 0.18 /h is accepted although
  `context_growth_rate` puts that culture's maximum at 0.14 /h. `mech/chain.py` applies
  `mu = min(mu_set, mu_max)` on both branches, which is the arithmetic `fba/fedbatch.py`
  already enforces at design time. Ablation on the growth channel, at the ethanol/37 C
  condition: **|d mu| = 0.1261 /h = 10.8x** `MEASURED_GROWTH_RATE_SE` (0.0117 /h). That is
  a correction to the incumbent, not new mechanism, and §6 carries the one-line diff.
- **The population layer.** Content per bearing cell moves little because the pathway sits
  near its own fitted ceiling; content per gram of CULTURE moves a lot, because a slower
  culture spends fewer generations in 120 h and keeps more of its plasmid. Time-averaged
  bearing fraction runs 0.19–0.80 across the declared sweep at the acid condition.

The acid's own growth arm reaches the product only above 161.43 mM acetate at medium pH 4
(the dose at which mu_max crosses the setpoint), against a declared lethal dose of 180 mM —
and §3 explains why that number should not be quoted as a prediction.

### (2) How many ODE states now exist?

**Five mechanistic states, against zero in the old stress path**, plus the three the
incumbent generator already integrated:

```
mechanism            5   A_i, pH_c, product_fraction, plasmid_bearing, generations
incumbent generator  3   biomass, immature reporter, mature reporter
total                8
```

The old stress path integrated the biomass logistic and the reporter and nothing else.
**The dose did reach both of those states — but only as a constant scalar.** In
`generator/culture.py::simulate_culture` it enters as `params.growth_rate_at(dose)`, one
number fixed for the whole run, and as `params.promoter_activity_at(dose)`, an amplitude held
constant in time (`np.full_like(t, k_synth)`). So the precise claim is not that the
environment reached nothing; it is that **no state of the environment's own was integrated,
and the dose carried no timescale** — it set two rates and then stopped being a variable.

What the five add is that two of them (`A_i`, `pH_c`) are driven BY the environment and have
their own measured time constants, so the dose now has a shape in time rather than only an
amplitude. That is the difference the geometry test in §3 is built to detect — and §3 is also
where it is reported that most of what the test detects is the fluorophore, not the cell.
`ode_state_census()` computes this table.

### (3) What does the tau/T audit say?

Run on the assembled state, in both declared windows, at 38.72 mM undissociated acetate:

**PLATE_READ_4H** (T = 4.14 h, mu = 0.245–0.363 /h)

| state | tau/T | bias | verdict |
|---|---|---|---|
| `A_i` | 0.00109 – 0.0814 | 0.1% – 8.1% | eliminate |
| `pH_c` | 0.000381 – 0.00064 | 0.0% – 0.1% | eliminate |
| `product_fraction` | 0.665 – 0.986 | 51.7% – 62.8% | keep |
| `plasmid_bearing` | 4.14 – 22.5 | 2.2% – 11.2% | straddles |
| `generations` | 0.461 – 0.683 | 40.8% – 52.5% | keep |

**FEDBATCH_5D** (T = 120 h, mu = 0.101 /h)

| state | tau/T | bias | verdict |
|---|---|---|---|
| `A_i` | 3.75e-05 – 0.00295 | 0.0% – 0.3% | eliminate |
| `pH_c` | 1.32e-05 – 2.21e-05 | 0.0% – 0.0% | eliminate |
| `product_fraction` | 0.0825 | 8.3% | eliminate |
| `plasmid_bearing` | 0.513 – 1.88 | 44.0% – 77.6% | keep |
| `generations` | 0.0572 | 5.7% | keep (INTEGRATING; the ratio test does not apply) |

**Three** states integrated on the plate, **two** in the vessel — and the two are exactly
the two `REVISED_BUILD_LIST.md` predicted would survive there. The criterion is
`tau/T < 0.146`, which is `panel_experiment.OBSERVED_ACTIVITY_CV`, and
`ReducedSystem.audit_table()` prints the caveat that a systematic bias and a random CV are
not the same kind of quantity.

**The `pH_c` rows are the repaired ones, and what they replace is this block's own
criterion-4 headline.**
<!-- audit:retracted pH_c is a state whose verdict FLIPS between the windows -->
The plate row used to read **keep**, at a tau/T above the criterion and a bias over half,
and that difference between the windows *was* the headline: irreducible on the plate,
algebra in the vessel. It was an artefact. `weak_acid_rhs` did not conserve charge at
mu > 0 — it diluted the anion and not the protons — so the error was mu-proportional and its
decay was read as a growth-driven slow mode. Repaired, charge balance leaves the coupled
system one-dimensional, `pH_c` is slaved to `A_i` **exactly**, and it eliminates in BOTH
windows. Three claims are withdrawn with it. Each withdrawn value sits beside the
measurement that replaced it in `mech/ph.py::ph_states` and in
`tests/test_mech_chain.py::test_both_acid_states_reduce_in_both_windows`, which is where
they are quoted rather than here — `audit_claims.py`'s unmarked-number ratchet stands at its
measured count, and a withdrawn value has no `outputs/` cell to be marked against:

- *"The two windows' reduction sets are not nested."* They are. Both acid states eliminate
  in both, and the vessel's two integrated states are a strict subset of the plate's three.
- *"The slow mode is growth-driven, tau of order one over mu."* That tau was the
  charge-balance error's own relaxation. The surviving mode is set by entry rather than by
  growth, so its tau is the SAME band in both windows and the two tau/T rows above differ
  only because the windows differ in length.
- *"`reporter_state_ablation` clears its floor."* Post-repair it does not, by two to three
  orders of magnitude — it was scoring the error itself. What survives on this arm is the
  Citrine quench in §6, which is algebra rather than a state.

### (4) How many reactions does pH now reach, and at what depth?

The baseline was **one** — the Citrine brightness multiplier, in one script, at depth 1.

Measured on the vendored yeast-GEM v9.0.2 with the acetic-acid uncoupling installed
(`fba/stress_ph.py::add_weak_acid_uncoupling`) and the proton budget closed, at a load of
3 mmol/gDCW/h, comparing parsimonious FBA before and after:

```
477 reactions move in EVERY run; a given run reports 484-487   across all 8 compartments

depth from the acid's own entry step      0:    3
(breadth-first on the metabolite-sharing  1:  141
 graph, restricted to reactions that      2:  254
 actually moved)                          3:   71
                                          4:   16
                                          5:    1
                                          6:    1
```

**The exact count is not reproducible and the histogram above is one run's.** pFBA has
alternate optima on this model, so the moved set is not a function of the constraints alone.
Measured over 8 repeats in one process: the count is 484 on the first call and 487 on every
later one, the union is 494, the intersection is 477, and **17 reactions (3.4% of the union)
flip between runs**. What IS stable across all 8 repeats: the maximum depth (6), the number
of compartments (8), and the shape of the histogram — the mass sits at depths 1 and 2 and
decays. Quote "about 480 reactions over seven depths", never a single integer.

`reachable_reactions()` and `reaction_depths()` compute it; the layer is `AUDITS`, because a
flux redistribution cannot make this chain's prediction move, only disagree with it. The
tests pin this as inequalities (`> 400`, `max >= 4`) for exactly this reason, and
`test_the_moved_set_is_degenerate_and_the_count_must_not_be_pinned` measures the instability
so that nobody later "tightens" those bounds into a number the solver does not owe them.

---

## 3. The geometry test, which is the sharpest check here — and it is half a negative

In the old generator each stressor's dose x time surface is a straight line through the
origin: every dose row is a scalar multiple of one common time profile, so a fitted latent
state of rank *k* cannot beat a random rotation of the same rank.

`ray_share()` is `sigma_1^2 / sum(sigma_i^2)` of the (dose, time) matrix, uncentred. On the
one channel the mechanism reaches — acetic acid at medium pH 4.5, six rungs from 5 to 90 mM,
19 transcriptional reporters, the 4.14 h read:

| | median | min | max |
|---|---|---|---|
| incumbent | 0.999516 | 0.995574 | 0.999516 |
| mechanistic | 0.975763 | 0.975763 | 0.980943 |

**The ray share falls, and the off-ray energy — the part a latent state could actually use —
rises by a factor of 50.1.** The mechanism that produces it is a dose-dependent *shape*
rather than a dose-dependent amplitude: the anion pool fills on its own time constant, the
growth rate moves over the read as it does, and the Citrine quench multiplies the trace by a
curve whose normalised shape differs between doses.

**And it is still 0.9758, not 0.5.** The surface remains overwhelmingly one-dimensional. On
this evidence a fitted rank-2 latent state would still be close to unidentifiable, and the
mechanism has NOT yet bought the thing it was built for. Saying otherwise from a 50x on a
2% residual would be reading a ratio as a magnitude.

**Worse, the decomposition says the whole effect is an instrument artefact.**
`dose_time_surface` takes `quench` and `growth_coupling` separately so the two halves can be
run apart. On `UPRE-ER`, the same ladder:

| what is switched on | ray share |
|---|---|
| incumbent (neither) | 0.99951607 |
| the growth coupling alone | 0.99950939 |
| the Citrine quench alone | 0.97572351 |
| both | 0.97576315 |

The pH-dependent brightness does all of it; the cell's own dynamics — an anion pool filling
on its own time constant and moving the growth rate over the read — move the ray share by
**6.7e-6**, a 1.001x change in off-ray energy. So what the mechanism bent is the shape of
what the READER records, not the shape of what the promoter does. That is a real and useful
result — it says a dose-response fitted on this instrument at low medium pH carries a
systematic pH slope — but it is not the mechanistic identifiability the block was built for,
and calling it that would be the error this document exists to avoid.

**The panel-wide number has not moved at all.** Over all 25 stressors and 19 transcriptional
reporters (475 surfaces, six-rung healthy ladders):

| | median | min | surfaces changed |
|---|---|---|---|
| incumbent | 0.999248 | 0.985960 | — |
| mechanistic | 0.999248 | 0.978482 | **19 of 475 (4.0%)** |

Identical to nine decimal places. This surface is `mech/chain.py`'s, and the chain routes
only `acetic_acid`, so twenty-four of twenty-five surfaces are the incumbent's, unchanged.

**The second seam does route two stressors, and it does not rescue the median either.**
`generator/panel_experiment.py::panel_dataset(mechanism=...)` drives `TRX2-oxidative` from
`mech/oxidative.py` under H2O2 and `UPRE-ER` from `mech/upr.py` under DTT —
`MECHANISED_CHANNELS`, two of the panel's 25 stressors and four of its 600 reporter by
stressor cells. Measured there, only H2O2's ray share collapses; DTT's barely moves, the
other 23 stressors are unchanged to machine zero, and the panel median moves in the fourth
decimal. `tests/test_panel_mechanistic.py` pins every figure and states in its own docstring
that quoting that median as evidence of a mechanistic dataset would be wrong.

That is the plain statement the build brief asked for: **on the panel as a whole the
mechanism has not bought what it was built for, on either seam, and it cannot until many
more than two stressors have an entry point.** One stressor's geometry is what was bought.

---

## 4. What is MEASURED, what is ASSERTED, what is REFUSED

### Measured — 23 MEASURED constants are registered on the in-chain path, plus 2 DERIVED

Counted by running `ParamRegistry.by_tag()` over `IN_CHAIN_REGISTRIES`: **23 MEASURED,
2 DERIVED, 4 SWEPT, 2 ASSERTED, 1 REFUSED — and zero BORROWED, zero uncited.** The last two
are the ones worth checking, and `uncited()` and `borrowed()` both return empty.

All of them live in `mech/ph.py` and `mech/burden.py`, each with its identifier at the point
of use. The four that do the most work:

- **Acetic-acid permeability across the *yeast* plasma membrane**, 1.4e-5 cm/s (Gabba 2020),
  with the warning in the parameter's own source string that the paper's prose transposes
  acetic and formic and that the 990e-5 figure in the same table is the *vesicle* number.
- **Resting cytosolic pH**, 7.08 ± 0.05 (Orij 2012, 2,088 biological replicates), and the
  same paper's measured invariance of pH_c to medium pH — which is why the medium-pH knob
  does nothing on its own in this model and only acts through weak-acid partitioning.
- **Lactic-acid permeability, zero** (Gabba 2020). A measured negative, and the block's
  sharpest unfitted test.
- **Growth loss per cassette copy**, 0.01 (Kafri 2016), linear with no threshold inside the
  range measured.

### Asserted — 2, and both belong to the assembly rather than to the biology

`mech/chain.py::CHAIN_PARAMS` holds exactly two: the initial plasmid-bearing fraction
(0.95, the same value `mech/population.py` declares, carried so the two cannot drift) and
the reference copy number (1, a declaration about the panel strains and not a measurement of
them). Both are tagged `ASSERTED` and both count as free scalars in the gate below.

### Swept — 4 declared sensitivity axes, and this is why `run_chain` has no default

`beta_cytosolic` (90–200 mM/pH), `cytosolic_volume` (1.0–2.7 mL/gDCW),
`p_loss_per_division` (0.0145–0.0847 per division), `b_burden` (0–0.234). `SweepPoint` has
no default constructor and `run_chain` requires one, so no number can leave this layer
without the corner that produced it. `chain_band()` runs all 16 corners; at the 175 mM acid
condition the culture content spans **0.186–0.823 mg/gDCW**, which is a wider interval than
any effect reported above and is the correct thing to quote.

### Refused — 33 named holes, each with the measurement that would close it

`mech/burden.py` refuses `m_stress`, the burden-by-stress interaction, which is the loop
this repository cannot close: even the *sign* of an environmental input to that block is
unmeasured for all 25 stressors. `mech/ph.py::NOT_BUILT` holds 18 (Pma1's Vmax and pKa, the anion export kinetics, the basal
proton conductance, the membrane potential, the vacuolar flux, the Rim101/Msn2/Haa1 gains,
and the permeabilities of propionic, benzoic, sorbic and CO2 — `acid()` raises
`AcidUnmeasured` for those by name). `mech/population.py::NOT_BUILT` holds 10, including all
five constants of the ethanol growth-inhibition law, and `mech/burden.py` holds 5 across its
two registries. The chain never reads any of them, and `float()` on one raises with the
reason and the missing measurement attached — verified: `float(NOT_BUILT["jmax_pump"])`
raises `RefusedValue`, not a number.

---

## 5. The four results that argue against this layer

Stated here, not buried, because each is pinned by a test.

### 5.1 The assembled free-scalar gate FAILS: 7 free scalars against 3 independent targets

Every block passes criterion (e) on its own — pH 2 against 3, population 2 against 2, burden
1 against 1. Assembled, the free scalars add and the targets do not, because
`GROWTH_RATE_FLOOR` is *one* growth assay that four registries each count as their own and
`REPORTER_ACTIVITY_FLOOR` is *one* plate CV that two count. Summing the per-block gates gives
7 against 8 and passes; de-duplicating the targets by name gives **7 against 3** and does
not.

```
free: beta_cytosolic, cytosolic_volume, p_loss_per_division, b_burden, m_stress,
      F0_bearing, n_copies_reference
distinct targets: half_yield_transfer, growth_rate, reporter_activity
```

`assembled_gate()` computes it; `require_assembled_gate()` raises with the count in the
message. The consequence is enforced in the API rather than written in a comment: this layer
reports bands, and a caller who asks to *fit* it is refused.

**What would change it.** Any one of: a replica-plating theta on this team's own construct
(closes `p_loss_per_division`); an OD600 growth curve of the bearing strain against an
isogenic empty-vector control in the same medium (closes `b_burden`); a wet/dry determination
with an inulin space (closes `cytosolic_volume`); a pH-clamp titration against a
simultaneous ratiometric pHluorin read (closes `beta_cytosolic` — and needs dual-excitation
optics the manifest does not have). Three of the four are cheap.

### 5.2 The acid's growth arm is below its floor everywhere it is physical

The pump-holds bill computes what holding pH_c at rest costs. At the resting pH the trapped
anion is `[AH]_o x 10^(7.08 - 4.757)` = 210x the outside undissociated concentration, so the
implied intracellular pool passes `mech/ph.py`'s own 500 mM osmotic flag above about **3.7
mM total acetate at medium pH 4.5** (2.8 mM at pH 4.0, 15.5 mM at pH 5.5). At the largest
dose where the state is still one a cell could hold:

```
d mu = 0.00395 /h = 0.34x MEASURED_GROWTH_RATE_SE
```

Below the floor of the assay that would measure it. Every larger number this arm reports —
including the 0.0894–0.1749 /h at 150 mM and the 161.43 mM crossover in §2 — is arithmetic
on a state a cell cannot occupy. `ChainResult.bill_is_physical` says so per call, the chain
emits a note, and `mech/ph.py` flags rather than clips because clipping would hide it.

The reason the arm is a bound and not a prediction is named: the cell exports the anion, and
`vmax_anion_export` and `k_anion_export` are both `NOT_BUILT`. The chain therefore carries
the two pH arms as a **bracket** — the bill (pH held, ATP paid, no quench) and the no-pump
trajectory (nothing held, no ATP, pH falls) — and reports the growth cost from the first and
the reporter artefact from the second, saying in a note that neither is the cell.

### 5.3 The geometry moved because of the optics, not because of the cell

Section 3's decomposition. The Citrine quench carries the whole 50x; the anion pool and the
growth coupling it drives carry 1.001x. `mech/ph.py` says as much in `citrine_quench`'s own
docstring — it calls the effect an artefact and gives its size against measured induction
folds of ~1.5 — and this is that sentence measured on a surface.

### 5.4 The chain reaches one stressor out of twenty-five, and the panel seam two

See §3. In this chain the pH block is the only one with an environmental entry point that
`stress_panel.STRESSORS` can express, and only for `acetic_acid`, because it is the only
weak acid on the panel with a measured *S. cerevisiae* permeability. Counting the other seam
as well takes it to **two of twenty-five** — H2O2 and DTT, through `mech/oxidative.py` and
`mech/upr.py` — and each remaining omission has a reason rather than a schedule, listed in
`panel_experiment.MECHANISED_CHANNELS`. Tunicamycin shares the UPR module but its entry is
REFUSED; menadione, diamide, cobalt and copper share the oxidative module but a redox cycler
and a metal do not reach Yap1 by the peroxide dose law.

---

## 6. What clears its floor

Scored under `docs/ABLATION_PROTOCOL.md`. Only the two observables with a measured floor in
this repository are scored; the product content is quoted as an effect size and NOT scored,
because it is an HPLC quantity and no floor for that assay has been measured here. Borrowing
the plate CV for it would be exactly the error the protocol exists to prevent.

| Piece | Observable | Floor | Effect | vs floor |
|---|---|---|---|---|
| The vessel rule (ethanol, 37 C) | growth rate | 0.0117 /h | 0.1261 /h | **10.8x** |
| The vessel rule (ethanol) | growth rate | 0.0117 /h | 0.1143 /h | **9.8x** |
| Citrine quench, 20 mM acetate pH 4.5 | reporter activity | 0.146 CV | 0.475–0.698 | **3.3–4.8x** |
| Citrine quench, 60 mM acetate pH 4.5 | reporter activity | 0.146 CV | 0.690–0.855 | **4.7–5.9x** |
| Declared burden, 8 copies | growth rate | 0.0117 /h | 0.0283 /h | **2.4x** |
| **The acid ATP bill, where physical** | growth rate | 0.0117 /h | 0.00395 /h | **0.34x — FAILS** |
| The E-Flux regulation layer | — | — | 0 bounds tightened | **0.00x — INERT** |

**`reporter_state_ablation` is absent from this table because it no longer belongs in it.**
It used to pass, and §2(3) is why it does not: it was scoring the charge-balance error. With
`pH_c` slaved to `A_i` exactly the residual is the anion pool's own dilution offset and not a
pH state, so this block has no state whose removal is measurable on an instrument it owns.
The two quench rows survive because they score the *algebra*. `mech/ph.py` carries the
withdrawn ratio beside the measured one.

The two quench rows are banded over all 16 declared sweep corners AND both ends of the
measured plate growth band, on the 4.14 h read. Swept over the buffering corners alone at one
growth rate the band is narrower (3.8–4.9x and 5.2–5.9x); the wider figure is the one to
quote, because the growth rate is measured with a spread and the quench depends on it.

The burden row is a *declared genotype axis*, not an environmental channel: nothing in
`L0`–`L5` writes a copy number, and `mech/burden.py` says so in the layer detail. It is
listed because it clears a floor, not because it closes criterion 1.

---

## 7. The change to `predict.py`, which has landed

**This section used to propose two diffs. Both are applied**, verified by reading
`predict.py` in the working tree rather than carried over from the earlier audit. That
module is owned separately, so what guarantees the behaviour is its own test suite and not
this paragraph — treat what follows as a record of the current state.

The defect was in `_growth_rate`: the reachability check sat inside `if environment.stressor
is not None`, so a setpoint with no stressor was returned unchecked and an ethanol fed-batch
at mu_set = 0.18 /h was accepted although `context_growth_rate` puts that culture's maximum
at 0.14 /h. It now reads the context ceiling on **every** call, multiplies it by the retained
fraction only when a stressor is present, and raises `SetpointUnreachable` outside the
stressor branch — the arithmetic `mech/chain.py` already applies in its own `vessel` layer
and `fba/fedbatch.py` enforces at design time. Re-measured against the current code, the
refusal fires for a stressor-free ethanol context at 0.18 /h and names 0.140 /h.

The behaviour change is therefore live, not hypothetical. Four contexts refuse where they
once answered, over the setpoints this repository actually uses (0.101, 0.18, 0.2543 /h):
`ethanol` (mu_max 0.1400, refuses at 0.18 and 0.2543), `ethanol` + 37 C (0.1282, the same
two), `growth_phase="stationary"` (0.0200, refuses at all three) and 20 C (0.2434, refuses
at 0.2543) — eight (context, setpoint) pairs, each ceiling re-derived from
`context_growth_rate` at this checkpoint. `scripts/product_environment_sweep.py` and
`scripts/env_to_product.py` both build ethanol contexts across that setpoint range; both
already catch `SetpointUnreachable` and record the refusal as a row, so those tables carry
refusals rather than failing.

The second, smaller diff has landed too. `predict_product` takes a `mech` argument and
`_mech_layers` imports `mech.chain` **inside the function**, the way `_gem_audit` imports
`fba.audit`, because `mech/chain.py` imports `predict.py` at module level for the
`LayerStatus` vocabulary and a module-level import the other way would be a cycle. It
defaults to `None`, for the reason section 0 gives.

---

## 8. What would change each conclusion

| Conclusion | What would overturn it |
|---|---|
| The environment moves the product | Nothing cheap. It rests on `context_growth_rate`'s measured carbon and temperature rates and on the vessel identity; if `context_growth_rate` is wrong the whole incumbent is |
| The acid growth arm is below its floor | A measured `vmax_anion_export`. The bill is an upper bound with the anion trapped; a real export rate lowers the pool into the physical range and the effect with it. Measuring it makes the arm *smaller*, not larger |
| The geometry has barely moved | More stressors with a mechanistic entry point, or a shorter read that resolves the pH transient. At 10-minute cadence over 4.14 h the fast arm is under-sampled by design. Two now exist on the panel seam and the median still does not move, so a handful more will not do it either |
| The assembled gate fails | One replica-plating plate (`p_loss_per_division`), or one paired growth curve against an empty vector (`b_burden`). Either takes the count to 6 against 3 — still failing. Both, plus an inulin-space volume, takes it to 4 against 3. It needs three of the four |
| Five states is the right number | The tau/T audit itself: on the fed-batch three of the five already reduce to algebra, so the vessel that makes the titre runs a two-state model, and on the plate `pH_c` reduces too. A measured `vmax_anion_export` is what would make `pH_c` a state again — export breaks the one-to-one between anion made and proton delivered that makes charge balance an algebraic constraint here |
| pH reaches 484 reactions | The proton budget. With `r_4527` open, Yeast9 disposes of a cytosolic proton for a third of an ATP through pyrophosphate excretion and the load is under-priced; `constrain_proton_budget` is what makes the number mean anything |

---

## 9. Running it

This runs the **chain** of sections 1–8, which is a bound on mu and a band, not the shared
physical engine. For the engine, use the `predict_protocol` snippet in section 0; the two
are different physics and neither substitutes for the other.

```python
from ystwin.mech.chain import run_chain, chain_band, SweepPoint
from ystwin.mech.state import FEDBATCH_5D
from ystwin.predict import Environment, Genotype
from ystwin.generator.context import CultureContext
from ystwin.pathway import calibrations
from ystwin.pathway.spec import load_pathway

spec = load_pathway("beta_carotene")
environment = Environment(context=CultureContext(carbon_source="ethanol"),
                          growth_rate_setpoint_per_h=0.2543)

# One corner. The SweepPoint is required, and that is the gate of section 5.1 in the API.
result = run_chain(environment, sweep=SweepPoint.midpoint(), window=FEDBATCH_5D,
                   genotype=Genotype(1.0), spec=spec,
                   flux_calibration=calibrations.BETA_CAROTENE_FLUX,
                   kinetics=calibrations.BETA_CAROTENE_KINETICS,
                   construct="UPRE1")
print(result.layer_report())
print(result.summary())

# The shape a number should actually be quoted in: 16 corners of the declared sweep.
band = [r.product.population_content_mg_per_gdcw for r in chain_band(environment, ...)]
```

`mech.chain.provenance()` regenerates the parameter tables of section 4, including every
source string, straight from the registries — so this document cannot drift from the code
without the regeneration showing it.

---

## 10. What an independent re-run confirmed, and the four numbers it corrected

Every headline above was recomputed from the public API in a separate process, at the
integration pass, rather than being carried over from the block reports. Reproduced exactly:
the incumbent's `1.07098788` across all eight named environments; the chain's six distinct
contents and 1.9867x culture spread at mu_set = 0.2543; the 0.186–0.823 mg/gDCW declared
band; the ODE census (5 + 3) and both tau/T verdict sets **as they stood before the
`weak_acid_rhs` repair, which then changed the plate's** — see §2(3); the panel-wide geometry
(median 0.999248 unmoved, 19 of 475 surfaces changed) and the acetic-channel fall to 0.9758
with the quench carrying 50.1x of the 50.0x; the assembled gate at 7 free against 3 targets;
and every row of section 6 except the two that are re-banded there.

Four corrections came out of that pass, and they are in the text above rather than in a
footnote:

1. **The GEM depth count was a first-run artefact.** 484 is what the first pFBA call in a
   process returns and 487 what every later one does; 17 reactions flip between alternate
   optima. Section 2(4) now quotes the stable core and the band, and two tests pin the
   degeneracy so the inequalities are never tightened into a false integer.
2. **"Zero of its integrated states were driven by the environment" was wrong.** The dose
   reaches both incumbent states through `growth_rate_at(dose)` and `promoter_activity_at(dose)`
   — as constants. The defensible claim, now in section 2(2), is that the dose carried no
   *timescale*, not that it carried nothing.
3. **The provenance counts were understated**: 23 MEASURED (not 13) and 33 refusals (not 27),
   both read off the registries by `by_tag()` and `refusals()`.
4. **Criterion 1 is setpoint-dependent and was quoted at its best setpoint.** At the brief's
   own mu_set = 0.18 the chain separates 3 of 10 conditions, not 6, because the only route the
   environment has to the product is to make the setpoint unreachable. Section 2(1) now carries
   both ends.

The `A_i` and `pH_c` rows of the tau/T tables were also re-derived at the dose the section
states. **That sentence used to end "which moved them slightly and changed no verdict", and
that is now wrong**: the charge-balance repair in `weak_acid_rhs` moved `pH_c` on the plate
by orders of magnitude and its verdict from keep to eliminate. §2(3) carries the retraction
and the three claims that went with it.

## 11. Published HOG extension: a different evidence source

`mech/hog.py` executes the original Petelenz-Kurdziel osmoadaptation model through the
explicitly supported SBML subset in `mech/kinetic_sbml.py`. `analysis/hog_data.py` retains
source identities, measured versus inferred quantities, assay normalization and whole-assay
splits. `analysis/hog_learning.py` estimates a conditional activation balance from measured
phospho-Hog1 data, rather than fitting generated latent-state labels.

This route does not satisfy every old-chain biological claim merely by living in `mech/`.
Most reaction parameters remain published priors, the source OD trajectory is imposed,
protein concentrations are literature-scaled, and the supplied reaction structure is not
learned. The measured preliminary Western holdouts were excluded from the source parameter
fit, but not from the original model-development process. Separate fast rates are not
resolved by the available sampling; their balance is the fitted quantity.

Only the volume-corrected Gpd1 protein-amount ratio is mapped to available enzyme capacity.
The SBML Fps1 regulator remains a native diagnostic: it is neither enzyme abundance nor an
identified GEM permeability law. Source dilution and measurement-volume corrections must
not be applied twice, and near-zero numerical states are not clipped to manufacture a
physical control. Increasing an upper capacity also does not force a glycerol flux; the
reported GEM probe explicitly shows that limitation.

Use `scripts/run_hog_learning.py`, not the legacy chain entry point. Source artifacts and
the frozen fitting protocol live in `data/hog2013/`; numerical results live in each run's
machine-generated files. Architecture and reproduction section 12 describe the workflow,
its independent numerical checks, its measured-data holdouts and its remaining gaps.
