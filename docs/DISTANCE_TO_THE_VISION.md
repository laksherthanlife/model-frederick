# How far is this from the thing it is supposed to be

The vision, as stated: *a model of GEM + FBA + metabolism + regulation + latent stress
state, such that given a pathway and a set of environment conditions it predicts how much
product you get.*

This is a historical accounting of the vision, with corrections where claims outgrew their
evidence. The dated scores, coverage tallies and engineering percentages below are not a
current CI report. The integration corrections distinguish source-cohort descriptions from
independent biology, and missing/refused requests from invariant predictions. Final retained
artifacts still require governed regeneration; externally inspected candidates are not
silently installed or promoted by this document.

---

## 1. The headline: five named layers, all of them now reachable — and one of them predicts

`predict.py` is the prediction chain. **Every layer the vision names is reachable from
it**, and a returned `ProductPrediction.layers` reports which layers ran and what they did.
Reachability is not a promise that every request reaches every layer: a setpoint refusal
can precede a product result. The comparison runner retains that request with missing
content and an environment layer marked `not-run`, rather than inventing a reported layer.

That is the change. It is *not* the same as "the vision is built", and the third column is
why: reachable is not the same as load-bearing.

| the vision names | built? | tested? | reached from `predict.py`? | **sets the number?** |
| --- | :---: | :---: | :---: | :---: |
| GEM / FBA audit | yes | yes | yes, opt-in (`audit_model=`) | **no** — audits |
| regulation (E-Flux) | yes | yes | yes, inside the audit | **no** — inert |
| latent stress state | yes | yes | yes (`latent=`) | **no** — reported, refused by G4 |
| thermodynamics | yes | yes | yes (`thermo=`) | **no** — audits |
| metabolism (pathway solve) | yes | yes | yes, every call | **yes** |
| flux law | yes | yes | yes, every call | **yes** |
| stress panel | yes | yes | yes, every call | only in batch |

**One layer sets the number, and that is the honest picture rather than a defect.** The
earlier version of this section counted 3,069 lines under `bridge/` and `fba/` that the
prediction never touched and called them unassembled. They are assembled now. What
assembling them proved is what the reasons below always said: four of the five *cannot*
move the answer, and running them makes that visible per call instead of arguable per
document.

- **FBA cannot supply a flux.** Cap every pathway reaction at the measured magnitude and the
  FVA floor stays at exactly zero. An upper bound cannot make a flux mandatory. It is opt-in
  because it costs an 18-second SBML parse and can only agree or refute; the returned content
  is bit-for-bit identical with and without it, and
  `tests/test_everything_is_wired.py` asserts that with `==` rather than a tolerance.
- **E-Flux tightens nothing**, for a structural reason: every transcriptional module a
  stressor touches is *induced*, so its scale is ≥ 1, and E-Flux applies upper bounds only.
  Sixty scales compute, zero bind. It now runs where a ceiling is the only thing that could
  matter — inside the audit — and reports `no reaction reached` there.
- **The latent stress branch is refused by G4**, so it is computed, attached and *labelled
  with that refusal* rather than fed in. Deleting it would hide the one route by which stress
  could reach the product; using it would launder an unvalidated branch into a headline
  number. Reporting it is the third option, and it is the one taken.
- **Thermodynamics was the exception, and it is now joined.** `pathway/thermo_gate.py`
  gates every step the spec declares chemistry for, at the concentrations that solve itself
  produced, swept across the cytosolic-volume convention. It cannot move the number — a
  driving force says whether a step may run, never how fast — so what it adds is the one
  thing the solver could not check about its own output: whether carbon was moved through a
  step pointing the other way.

**What the join actually bought, on real chemistry — and what it then cost.** `phb`'s entry
step is yeast-GEM's cytosolic thiolase `r_0103`, and the gate reaches it from inside a
prediction. The table this section carried until 2026-08-30 read:

| feed | cytosolic acetyl-CoA | verdict | dG |
| --- | ---: | --- | ---: |
| glucose | 10 µM | **cannot run** | +15.61 kJ/mol |
| ethanol | 425 µM | ~~**runs**~~ | ~~−3.29 kJ/mol~~ |

The second row is withdrawn. `bridge/thermodynamic.py` was reading the vendored ModelSEED
table — which declares itself kcal/mol — through `pytfa` under its default `thermo_unit`
of kJ/mol, so the pH transform landed in kJ on top of an energy in kcal. Water came out at
**+24.85 kJ/mol** on that scale, which is what exposed it. Corrected:

| feed | cytosolic acetyl-CoA | verdict | dG |
| --- | ---: | --- | ---: |
| glucose | 10 µM | **cannot run** | +38.07 kJ/mol |
| ethanol | 425 µM | **cannot run** | +19.17 kJ/mol |

The threshold moves from 221 µM to 19 mM and **both feeds are below it**. The straddle was
the whole of the thiolase lead, and it existed only in the erroneous number. This is not an
artefact of the correction either: eQuilibrator never touched that bug, its +24.96 kJ/mol is
the literature value for a thiolase condensation, and it puts the threshold at 1.4 mM —
ethanol misses on that reading too. The check was available from inside this repository from
the day the eQuilibrator backend landed, and nobody ran it.

What the chain can do is unchanged: it *asks* the question, on a solved pathway, and item (4)
below is still what answers it. What is gone is the answer we thought we had.

**And what it did buy, after two of our own errors were cleared out of the way.**
`beta_carotene`'s three steps all came back `cannot_say` until 2026-08-30, and neither reason
was a fact about carotenoid chemistry. The first was a naming gap, closed by
`data/thermo/heterologous_metabolites.tsv`. The second was `fba/carotenoid.py` writing CrtI
with four free FAD as terminal acceptor — thermodynamically impossible at +166.2 kJ/mol, on a
step Verwaal 2007 measured running to completion. Written as the flavin oxidase it is, the
sign disagreement between the two estimators is gone.

**This paragraph used to end "all three steps gate as `runs` at −116.5, −327.1 and −294.2
kJ/mol with nothing ungated." That is false and was false when written.** Running the gate
returns `0/3 declared steps carry an energy; 3 cannot say`, and
`tests/test_everything_is_wired.py:200` asserts exactly that —
`{s.feasibility for s in report.steps} == {Feasibility.CANNOT_SAY}`. `git show` puts the test
and the claim in the same commit, so the contradiction was born rather than drifted into.
Corrected 2026-09-02. The magnitude disagreement between the estimators refuses a NUMBER, and
without a number there is no verdict — `cannot_say` is the honest output, not `runs`.

### How far along is this, as a number

A single percentage needs a denominator, so here is the one used, with every row measured
rather than judged. **The weighting across the three groups is a judgement and is stated as
one**; the rows inside them are not.

| dimension | measured | score |
| --- | --- | ---: |
| layers reachable from `predict.py` | 7 of 7 | 100% |
| lines under `bridge/` + `fba/` the chain reaches | 2,655 of 4,436 (8 modules of 14) | 60% |
| integrity gates passing | claims 175/175, determinism 72/72, reproducibility 284 pass / 1 skip / **0 fail** | 100% |
| suite | 3,241 collected, 0 failing | ~100% |
| **engineering** | | **~85%** |

*The four rows above are measured; the group score is not a mean of them and never was.
The rows average ~94% while engineering is scored ~85%, because the score is a judgement
about what remains rather than an aggregate — and the 2026-08-30 audit named nine open
engineering items, so a rows-mean would be the wrong number to quote. Read the rows as
evidence and the group score as an opinion informed by them. `equilibrator.py`, `tmfa.py`,
`carotenoid.py`, `surrogate.py`, `insulin.py` and `secretion.py` are the six modules the
chain does not reach, and the first of those is the estimator that refuted the thiolase
lead. The line ratio fell from 78% to 60% without the chain losing anything: `insulin.py`
and `secretion.py` were written after that row was last measured, and both are secreted-
product chemistry that `PathwaySpec` refuses at load time by design.*
| layers that can move the answer | 2 of 7 | 29% |
| user questions with a real answer (§2) | 2 of 8 | 25% |
| historical axes-working tally | genotype yes, environment labelled **refuted** — the universal inference is withdrawn in §5 | 50% (historical) |
| pathway specs that load and calibrate | 1 of 3 | 33% |
| declared steps the energy tables cover | 1 of 4 | 25% |
| **scientific scope** | | **~30%** |
| claims at Tier 3 — predicted first, measured after | **0** | 0% |
| biosensor claims that survive multiplicity | **0 of 7** per-dose — the family-wise level is unreachable at four plates, so `adjusted_interval` refuses it rather than quoting resampling noise; **4 of 4** per-sensor | — |
| resolution of the headline claim | permutation *p* = 1/6 = 0.167, **the floor at three strains** | — |
| **evidence** | | **~5%** |

Weighting engineering 20 / scope 40 / evidence 40 — because the vision is a *predictive*
claim, so what it predicts and whether anyone checked should outweigh whether the code is
tidy — gives:

> ## **≈ 30%**

**The wiring moved the first group and almost nothing else.** Engineering went from roughly
70% to 85%; the headline moved about three points. That is the honest accounting of it, and
it is worth stating plainly because assembling five layers *feels* like the largest possible
change and is not: four of the five cannot move a number, which is what wiring them
demonstrated rather than what it fixed.

**All three integrity gates now exit 0, which they have never done before.** The last red row
was `igem_results()` resolving to a sibling directory by convention, so every July plate
result was reproducible on exactly one machine. It is fixed the way this repository already
knew how: `plate_readings.py --export` now writes both wet-lab plate sets as text and
verifies both value-for-value, the resolver's fallback is gone, and `run_gates.py` and
`run_d2.py` read the committed text. Neither set of workbooks can be redistributed — both
name a private individual in `docProps/core.xml` — and neither needs to be. Fifteen checks
report SKIP, which is not a pass: they are the tables whose raw inputs are not on this
machine, and the staleness is still printed in hours.

**What is left in the engineering column is now under 15 points with nothing named in it.**
That is the ceiling, not an invitation: the column is worth 20% of a score dominated by the
two rows below, and both of those are chemostats.

### What 2026-08-30 changed: the instrument, not the twin

**The grand-vision number did not move.** Every scope row is identical — 2 of 8 questions,
one calibratable pathway, the environment axis refuted, 1 of 4 pathway steps with an energy —
and evidence is still zero at Tier 3.

**But there are two claims in this repository and only one of them is the twin.** The README's
own one-sentence contract is *"give it one well's brightness and cloudiness readings over time,
and it tells you how hard that gene was actually switched on"*. That is the **instrument**, and
it is what the whole package is built on. It just got its best evidence, from a dataset this
project did not collect for the purpose and had never opened.

`plate/gen5.py` reads the Gen5 `.xpt` binary — every plate this project has, at one more
decimal, all 96 wells instead of the 45–93 the exports kept. That opened the auto feedback loop
plates: a galactose-inducible circuit carrying its own repressor, 121 timepoints over 20 hours,
three biological replicates. And on them the naive and corrected readouts **disagree about the
direction of the response**:

| strain | naive `RFU/OD` 0 → 1.8% | corrected activity |
| --- | --- | --- |
| No-LacI | 3390 → 1354 (**0.40×**) | −151 → **+85** |
| No-Gal | 2928 → 2270 (0.78×) | −119 → **+263** |

On the biosensor plates the correction *removes* an apparent response, which always invites
"the correction is eating signal". Here it **recovers** a monotone, triplicated dose response
that `RFU/OD` was inverting. Same thesis, and the direction that is far harder to argue with.

**With one rival explanation that this experiment cannot exclude.** Galactose slows these
cultures, so a dose ladder is also a growth ladder: pooled Spearman(galactose, `mu_late`) is
+0.891, and conditioning the dose association on growth rate takes it from **+0.839 to
+0.494**. Within a single strain on a single plate it survives almost intact (median partial
ρ +1.000, 15 of 16 blocks positive), so growth rate does not account for the ranking — but
roughly half the pooled signal is shared with it. `mu_late` is a column of
`outputs/afl_circuit.csv` so this is checkable rather than assumed, and separating the two
needs a galactose-insensitive promoter on the same ladder, which nobody has built. Counted
on the instrument's side of the line, not the twin's, and counted smaller than it first
looked.

These plates are also the only ones in the project measured where the estimator is accurate:
20 hours gives a 3.33 h window and ~1.3% error, against the four-hour plates' floor and ~8%.

**And the same reader moved a number, which the rest of this section did not.** Plate
`20260804` was outside the biosensor panel because every `.xlsx` of it had been exported
through a blank-subtraction transform; the `.xpt` still carries the background, and reading
it takes the panel from **n=3 to n=4** for all four constructs. Seven per-dose folds now
clear "no change" against five, and all seven hold across the whole late-window sweep against
four of five.

**No fold survives a family-wise correction, and the reason is sharper than that.** This
paragraph briefly claimed one did — UPRE1 at 0.5 mM DTT. The 2026-08-30 audit found the
divisor was the 20 *estimable* folds rather than the 24 *asked about*, and, underneath that,
that `adjusted_interval` was answering a level four plates cannot express: a cluster
bootstrap on four plates has a smallest atom of `(1/4)**4`, so the finest two-sided level is
0.0078 against a family-wise 0.0021. **Five plates is the first count at which the question
can be asked at all.** The exact sign-test floor moves with the plate count, from 0.25 to
0.125, so it is still above 0.05: six plates remains the first count at which a single fold
could be significant on its own terms.

That is a real gain and it is on the instrument's side of the line, not the twin's. It also
cost the project one number in the other direction: `AlteredYap1` at 1.0 mM H2O2 went from
1.49 to 1.00 once a fourth plate disagreed, and the per-sensor permutation test weakened on
both oxidative sensors.

**So: the instrument is in much better shape and the twin is exactly where it was.** That is
worth separating rather than averaging, because the score below weights the twin.

### What 2026-08-29 changed, and what it did not

**The number did not move, and that is the finding.** Engineering is a point or two higher —
all three gates green, one more layer of self-checking. Scope is identical: still 2 of 8
questions, still one calibratable pathway, still the environment axis refuted. Evidence is
still zero at Tier 3.

What changed is **whether the claims already there can be believed**, and in every case the
answer got smaller and firmer:

*This table is the 2026-08-29 state. The fourth plate moved the first row again the next
day — 7 unadjusted, and family-wise unreachable rather than merely unmet — as the section
above records.*

| claim | before | after |
| --- | --- | --- |
| resolvable biosensor doses | 6 | **5**, then **0** under family-wise correction |
| what replaces them | — | **4 of 4** sensors show a dose response at *p* ≤ 1e-4 |
| crosstalk specificity | n = **1**, mislabelled as n = 3 | n = **3**, with the on-target half |
| activity estimator | window from a rule tuned on 24 h runs | measured at the 4 h geometry |
| `fold_change` on a dying culture | reported a call | refused |

A day spent removing claims and adding one better one. That is not the same as progress
toward the vision and should not be counted as it: **nothing here moved the product model,
because nothing here could.** What it did is make the biosensor half of the project defensible
— the half that has real data — so that when the wet lab does produce a held-growth-rate run,
the thing it lands on is not resting on a single plate, an uncorrected family of twenty, and a
smoother nobody had measured.

> **Updated 2026-09-02: that run is no longer a chemostat.** `fba/fedbatch.py` replaces it
> with an exponential fed-batch, and every "chemostat at D = x" below should be read as "a
> vessel held at mu = x" — see [`PROTOCOLS.md` §P3](PROTOCOLS.md). The substitution is exact
> for the pools law, because the settled fed-batch residual substrate IS the chemostat's Monod
> residual at the same growth rate. It is NOT exact at the top of the range: `D = 0.35` was
> already refused here, and `mu = 0.254` is now refused too, because it is the strain's
> measured uptake ceiling.

**The cheapest point per unit of work is not code.** One chemostat at D = 0.18 against
`outputs/registered_prediction_D018.csv` converts the evidence row from 0% to a real Tier 3
entry — the first in the project's history — and moves the headline further than any
refactor can.

**And the second cheapest is one strain, which is also not code.** The genotype axis — the
one axis scored as working — rests on three strains, and three strains cap the permutation
*p* at 1/6 = 0.167 no matter what analysis is applied. A fourth takes the floor to
1/24 = 0.042, the first value that can cross 0.05.
[`SECOND_DATASET_HUNT.md`](SECOND_DATASET_HUNT.md) records that the authors' own release
does not contain one: `SysBioengLab/BcarGRASP` carries the same six conditions and three
strains, and its `ERG9` tables are an estimated squalene flux rather than a fourth strain.
So it is a chemostat run, and it is on none of the five specifications in
[`MEASUREMENTS_NEEDED.md`](MEASUREMENTS_NEEDED.md).

## 2. What a user can actually ask

Eight questions from the **historical** interface check, with the answers then reported.
These values and the old two-of-eight tally are retained for context, not current coverage:
notably the ethanol row below was an unreachable held setpoint, not a valid invariant
prediction. The corrected six-request comparison follows the historical notes.

| question | historically reported answer |
| --- | --- |
| chemostat at D = 0.18, this strain | **1.071 mg/gDCW** |
| a strain with 10× the cassette dosage | **1.235 mg/gDCW** |
| chemostat at D = 0.05 | refused |
| chemostat at D = 0.35 | refused |
| shake flask / batch | **1.117 mg/gDCW at 20 °C on glucose** — see correction below |
| the same strain fed ethanol | 1.071 — *identical* |
| the same strain at 37 °C | 1.071 — *identical* |
| the same strain under 1 mM DTT | 1.071 — *identical* |

> **Correction (2026-09-03): the flask row was wrong, and it was wrong about the vessel.**
> The refusal is about `mu` leaving the fitted window `[0.101, 0.2543]`, not about batch
> culture. Calling `predict_product` with this repo's own carbon × temperature model returns
> real numbers for **6 of 12 flask contexts**: glucose at 20 °C 0.9907, glucose at 15 °C
> 1.1167, galactose at 25 °C 0.9989, galactose at 20 °C 1.0826, ethanol at 30 °C 1.1179,
> ethanol at 25 °C 1.1414 mg/gDCW. A shake flask at 20 °C on glucose is already inside the
> model's answerable set.
>
> The deeper point is `scripts/batch_sufficiency.py`: the balance `solve_pathway` closes is
> **per gram of biomass**, so the vessel does not appear in it. Integrating it forward at
> constant `mu` converges on exactly the steady-state answer, in **3.7–4.0 generations to
> 10%, 4.8–5.1 to 5%, 7.2–7.5 to 1%** — the same count in a flask as in a chemostat.
> **`mu` must be KNOWN and roughly CONSTANT; it does not have to be HELD.** Every "chemostat
> at D = x" ask in this document should be read as "a culture whose `mu` you measured and
> which held it long enough", and `docs/PROTOCOLS.md` §P4 is how to run that in a flask.

The old summary said **two of eight real answers**, three refusals and three identical
returns. It had already corrected an earlier "one of eight / four identical" tally, and
the dated flask correction changed the coverage again. None of those counts is a current
coverage guarantee. In particular, an environment that cannot sustain the requested μ
must refuse rather than return the same number and be counted as evidence of invariance.

**Corrected fixed-gene comparison, inspected externally while retained artifacts await
governed regeneration:** `scripts/env_to_product.py` retains all six β-carotene requests.

| requested μ (/h) | glucose | ethanol | answered / refused requests |
| ---: | --- | --- | --- |
| 0.101 | answered | answered | 2 / 0 |
| 0.15 | answered | `SetpointUnreachable` | 1 / 1 |
| 0.2543 | answered | `SetpointUnreachable` | 1 / 1 |

That is **four answered plus two refused**, not six answers. The ethanol context maximum
is 0.140 /h. Refused rows preserve carbon source and requested μ, the refusal class and
reason, missing returned growth/content and layer flag, and an environment layer state of
`not-run`. Answered rows have a `reported` environment layer which does not set the number.
Only the low-rate pair has two finite answers to compare. Their equality follows from the
fixed-gene μ-only prediction, not biological invariance; higher-rate matched-carbon
comparisons are unavailable. `tests/test_env_to_product.py` requires this full population
and missingness rather than allowing `nunique` to turn one surviving answer into a finding.

## 3. Pathways

Three specs exist. **One loads and can be calibrated.**

| spec | loads | calibratable | why not |
| --- | :---: | :---: | --- |
| `beta_carotene` | yes | **yes** | — |
| `phb` | yes | no | one measurable node, and no expression measurement exists |
| `glycogen` | **no** | — | refused: GPH1/SGA1 degrade it and the solver has no term for that |

The generality claim is *structural* — nothing in `src/` names a product, a new one is a
TOML file — and it has been demonstrated exactly once, on a pathway that shares the host,
the promoter set and the paper with the original.

## 4. Every claim, and how it was earned

| claim | number | how | tier |
| --- | --- | --- | --- |
| genotype → product | 14.2% median, worst 1.33× | leave-one-**strain**-out | **2** |
| expression → pathway flux | skill +0.63, 1.22× typical | leave-one-**strain**-out | **2** |
| growth rate, external paper | +5.5% vs Torello Pianale | fully held out | **2** |
| environment → product | corrected log-flux constant matches the comparator; three nonconstant candidates remain negative | retrospective state-LOO against fold-training mean; random-design percentiles are not candidate nulls | inference pending, not a universal refutation |
| registered D = 0.18 | written, unmeasured | — | 0 |
| β-carotene ceiling 1.25 mg/gDCW | **REFUTED ~63×** — Arhar 2024 PMID 39215465 measured 79 mg/gDCW by HPLC on gravimetric DCW with lycopene below detection | `Genotype.cassette` now carries per-gene dosage, but `capacity_for` refuses to scale the anchor at another crtYB dosage; `published_cassettes.py` reports only two distinct numeric crtYB dosages across different study contexts, not an identified dosage-to-capacity law | **refuted** |
| thiolase threshold mechanism | historical ordering held 8 of 9 within each convention and across their cross product, failing at 2.7 mL/gDW with [CoA] = 300 µM; the subsequent unit correction leaves both feeds below threshold (§1) | **two labs spliced**; the former favorable ordering depended on the erroneous energy | original threshold lead refuted |

The former "three Tier 2 claims, one refutation" tally is historical, not a current evidence
inventory. The fixed-CrtE strain-holdout numbers above were selected after examining the
same cohort and are not the corrected nested selected-pipeline scores; see
[`EXTERNAL_PRODUCT_VALIDATION.md` §1](EXTERNAL_PRODUCT_VALIDATION.md). The universal
environment-refutation claim is withdrawn while the poor candidate scores and failed
thiolase lead remain. Nothing here is prospective Tier 3 evidence.

## 5. The distance, stated plainly

The vision is **product X, pathway Y, environment conditions**. What exists is:

> **one product**, through **one pathway**, from **genotype**, at **one dilution rate**,
> in **one host**, calibrated on **one paper**.

The environment axis — the reason to build a twin at all rather than a lookup table — is
not independently validated by Kocharin's eleven condition summaries. The old stronger
verdict, "refuted" because every candidate lay inside a permutation null, is withdrawn.
It mixed a full-data-mean comparator with a random-design demonstration and inferred
unidentifiability from neither a valid candidate null nor independent experimental units.

Under the corrected fold-training-mean **log-flux** comparison the constant has zero skill
to roundoff, while all three nonconstant candidates remain negative. Their
`pending_exchangeability` status withholds inference. Separately, a retrospective **content**
fit improves the median within the same cohort but is more than 2× wrong on two of eleven
states. These observations neither erase the poor predictions nor establish an environment
law on new biology. [`EXTERNAL_PRODUCT_VALIDATION.md` §3](EXTERNAL_PRODUCT_VALIDATION.md)
keeps the old and corrected comparisons distinct; §7 states what is still missing.

---

## 6. What to do, ranked by what it buys per unit of work

### Cheap, and worth doing regardless

1. ~~**Wire the thermodynamics layer into the chain.**~~ — **done, and here is exactly how
   far it got.** `predict.py` takes `thermo=` and gates every step the spec declares
   chemistry for, across the cytosolic-volume convention, at the concentrations that solve
   produced. `refuse_thermodynamically_blocked=True` turns the report into a refusal for a
   caller who wants one; the default reports, because the tables cover 44.5% of yeast-GEM's
   reactions and a refusal resting mostly on `cannot_say` would fire unevenly.

   What made the join possible was the missing half this entry named: **a spec that declares
   its own chemistry.** `Node` now carries `metabolite`, `reaction` and `stoichiometry`, and
   both shipped TOMLs declare all three. The stoichiometry is what lets the gate run with no
   SBML parse; `PathwaySpec.check_steps_against(model)` compares the two copies whenever a
   model is at hand, so the duplication is a cross-check on `fba/carotenoid.py` rather than a
   second copy of it.

   **The `background_m` caveat this entry raised has NOT gone away, and it should not be read
   as having.** `phb`'s thioesters are still `passthrough` and still solve to zero, and
   acetyl-CoA is still not a node, so the acetyl-CoA level in the thiolase verdict is still
   supplied rather than solved. What changed is that the verdict is now produced inside a
   prediction and labelled, instead of only inside `scripts/thiolase_threshold.py`. What
   would close it properly is measurement (4).

   *What it bought:* `predict_product` on `phb` returns **cannot run** at 10 µM acetyl-CoA
   and, since the 2026-08-30 unit correction, **cannot run** at 425 µM as well — so what the
   chain can state about a solved pathway is a refutation rather than the mechanism §4 used
   to describe. `beta_carotene` returns **`cannot_say` on all three steps** — the naming gap
   and the impossible CrtI stoichiometry were both cleared, which removed the sign
   disagreement, but the estimators still differ in magnitude by more than the gate will
   accept, so no step carries a usable energy. (This line read "`runs` on all three steps, at
   −116.5, −327.1 and −294.2 kJ/mol" until 2026-09-02; the code and its test always said
   otherwise.) See §4 and `docs/FINDINGS.md`.
2. ~~**Make the FBA audit reachable from `predict_product`**~~ — **done.** It takes an
   already-loaded model, so a caller pays the 18-second SBML parse once and reuses it, and a
   caller who does not want the audit pays nothing. The number is bit-for-bit unchanged with
   and without it, which is the property that makes it an audit: an upper bound cannot make a
   flux mandatory, so it can only ever say whether the network could carry what the kinetics
   asked for. `audit_reaction` and `audit_glucose_uptake` are required alongside the model —
   the first because this function cannot know which of a GSMM's reactions is the product,
   the second because the growth cost swings sevenfold with it.
3. **Give `phb` and any future spec an expression axis or mark it non-calibratable at load.**
   It is currently only discoverable by reading `spec.calibratable`.

### The measurements, in order of what they settle

4. **Cytosolic acetyl-CoA, one lab, one quench protocol, glucose- and ethanol-limited at
   matched D.** *Addresses:* the 43-fold cross-laboratory pool splice and a specified
   thermodynamic condition. The original ranking called this the highest-value experiment
   and the sole route to an environment model. That claim no longer stands: the energy-unit
   correction refuted the proposed threshold ordering (§1), and the invalid flux-null bar
   did not establish unidentifiability. Same-laboratory concentrations would not by
   themselves identify a flux law or validate product prediction.
5. **β-carotene content at a *known growth rate*.** *Settles:* the one part of the ceiling
   that is still open. The level is not: §4 of this file records the ceiling as **REFUTED by
   ~63×** (Arhar 2024, PMID 39215465, 79 mg/gDCW by HPLC). This item previously read
   "nothing has ever tested it", which the same document refutes eleven lines earlier. What
   no refuting paper reports is a growth rate — they are all flask cultures — so the
   μ-dependence the model asserts is untested even though the level is settled.

### The structural work, which needs the measurements first

9. **A concentration-based flux law.** `v = kcat·E·[S]ⁿ/(Kⁿ + [S]ⁿ)` with n fixed by
   chemistry. This is the only candidate that explains *why* the growth exponent flips by
   feed rather than noting that it does. **Blocked on (4)** — and hard-blocked on the
   carotenoid axis, where intracellular FPP and GGPP concentrations do not exist in usable
   form in any published work.
10. **A second host or promoter set for the flux law.** Three strains from one paper is one
    experiment. Whether the *shape* transfers is untested.
11. **A second product that is genuinely unrelated.** PHB was the right choice and it cannot
    calibrate anything. What is needed is an intracellular, non-secreted, non-degraded product
    with both a content series *and* an expression measurement. Nothing found so far has both.

### What not to do

- **Do not promote an environment term on the old null-bar argument.** The former rule
  requiring skill near 0.08 at three parameters is withdrawn: it was an artificial
  random-design percentile, not a candidate significance or identifiability threshold.
  `scripts/flux_null_baseline.py` now separates that demonstration from four descriptive
  candidates with inference pending. A new term needs a declared target, independent-unit
  evaluation and any justified selection-aware inference, not retuning to a percentile.
- **Do not wire E-Flux into the prediction.** It is provably inert, for a structural reason.
- **Do not extend to secreted or degraded products** without adding the missing terms first.
  Both are refused at load time and both refusals are correct.

## 7. The honest summary

The genotype comparisons are retrospective within one source cohort, not proof of
transport to new biology. On the environment axis, three fixed log-flux candidates still
score poorly, the separate content fit has a bad tail, and two β-carotene requests now
correctly refuse. None of that licenses the former universal null-bar/identifiability
claim, or treating missing answers as invariance.

The proposed thiolase ordering was refuted by the energy-unit correction, not identified
as the mechanism; the two-laboratory concentration splice remains unresolved. Item (4)
would measure that splice under one protocol, not by itself close the distance to a model.
Independent cultivation identities, a declared validation/selection procedure and new
measurements are still needed. Governed regeneration can correct the retained artifacts;
it cannot supply independent biological evidence.
