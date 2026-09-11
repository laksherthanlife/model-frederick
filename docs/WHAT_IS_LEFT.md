# What is left

Measured against the stated goal: *a model of GEM + FBA + metabolism + regulation + latent
stress state, such that given a pathway and a set of environment conditions it predicts how
much product you get.*

Written after the August 2026 rebuild. The layer inventory and dated progress reports
retain that history, not a current claim that every "done" item has independent biological
validation. The environment and conditional-kinetic corrections below narrow those claims;
final retained-artifact regeneration/adoption remains separate from inspecting a candidate.

## Against the five named pieces

| Piece | State | Evidence |
| --- | --- | --- |
| **GEM + FBA** | **built, reached opt-in, and cannot move the number** | `fba/audit.py` checks feasibility, the precursor budget and the growth cost, and cannot supply a flux -- capping every pathway reaction at the measured magnitude leaves the FVA floor at exactly zero. Two earlier versions of this row were wrong in opposite directions: one said it was in the chain when `predict.py` did not import it, the next said it was not reachable at all. It is reachable, through `audit_model=`, and `ProductPrediction.layers` reports it as `audits` -- ran, can refute, cannot move. `tests/test_everything_is_wired.py` asserts the content is identical with and without it using `==`. |
| **Metabolism** | **done** | `pathway/solve.py` walks any declared chain in closed form. Reproduces the hand-written β-carotene solver bitwise; carbon closes to 1e-19. |
| **Product from environment + genotype** | **historical fixed-candidate result, not general validation** | Leave-one-strain-out: median 14.2% error, worst 1.33×, from expression and dilution rate, excluding the held-out strain's product from those fits. CrtE/model choice was made after examining the same cohort; this is not the corrected nested selected-pipeline score or validation of other intracellular products. |
| **Regulation** | **wired, and provably inert — which is the finding** | E-Flux now reaches `fba/audit.py`, which was the right place: scaling a ceiling cannot set a flux, but it can narrow the envelope a prediction is checked against. Having asked the question, the answer is that it tightens **zero** bounds for any stressor at any dose. Reasons below. |
| **Latent stress state → product** | **computed, attached, and deliberately not used** | In a chemostat the pump fixes μ, and μ was stress's only channel, so a stressed prediction is bit-for-bit the unstressed one. `predict.py` emits a note saying "no route modelled", not "no effect" — because silence there reads as a finding. Passing `latent=` now runs `bridge/latent_bridge.py` and returns its maintenance demand on the result under `LayerState.REPORTED`, carrying G4's refusal in the label. Reporting rather than using is the whole point: it keeps the one candidate route visible without letting an unvalidated branch reach a headline number, and the `==` assertion in `tests/test_everything_is_wired.py` is what stops that changing by accident. |

| **Thermodynamics** | **built, joined, and limited by table coverage rather than by wiring** | `pathway/thermo_gate.py` is reached from `predict.py` through `thermo=`. On `phb` it reproduces the thiolase threshold inside a prediction — `cannot run` at 10 µM cytosolic acetyl-CoA, and since the 2026-08-30 kcal/kJ correction in `bridge/thermodynamic.py`, `cannot run` at 425 µM too: the threshold is 19 mM, not 221 µM, and both feeds are below it. That refutes the thiolase lead, and the gate is unaffected — refusing an infeasible step does not require the step to be feasible somewhere. On `beta_carotene` all three steps come back **`cannot_say`** — verified by running it, and asserted by `tests/test_everything_is_wired.py:200`. (This cell claimed **`runs`** at −116.5, −327.1 and −294.2 kJ/mol until 2026-09-02, in contradiction of the test written in the same commit.) They were `cannot_say` for two days and the reason was ours, not the literature's: `fba/carotenoid.py` wrote CrtI with four free FAD as terminal acceptor, which is thermodynamically impossible (+166.2 kJ/mol on a step Verwaal 2007 measured running to completion) and made the two estimators straddle zero. Written as the flavin oxidase it is, the sign disagreement is gone; the remaining magnitude disagreement refuses a NUMBER and not a VERDICT, because both estimators agree in sign by more than they differ. |

## Why regulation cannot bite, exactly

Wiring E-Flux into the audit took about fifty lines, as estimated. What it bought was not a
narrower envelope but a precise statement of why there will not be one, and that turns out to
be the more useful thing.

`gene_scales` computes `scale = max(0, 1 + activity)` and E-Flux applies **upper bounds
only**. So a bound tightens if and only if a module is **repressed**. Two facts, both
measured rather than argued:

1. **Every transcriptional module a stressor touches is induced.** Positive activity gives a
   scale at or above 1, which raises a ceiling. Raising a ceiling of 1000 to 1500 does
   nothing to a flux of 10.
2. **Every repressing signal in the whole panel sits on a metabolite pool** — `redox` under
   DTT, `nadh` under menadione, `atp` under glucose starvation and antimycin A, `ph` under
   acetic acid. Pools are not transcriptional, carry no regulon, and `gene_scales` **refuses**
   them rather than pretending. The one direction that could tighten a bound is precisely
   the direction with no route to a gene.

So sixty scales compute and zero bind — not a wiring fault, a structural one.
`tests/test_fba_audit.py::TestRegulationIsWiredAndProvablyInert` pins both facts, so if a
repressed transcriptional module is ever added, or a **measured** reference capacity replaces
the model-implied one, the tests fail and that is the signal to look again.

## The stress gap, stated precisely

This is the largest hole and it is worth being exact about, because it is easy to describe
as smaller than it is.

Two routes by which stress could plausibly move product, and why neither is wired:

1. **Maintenance-ATP burden.** Stress raises non-growth ATP demand, which diverts carbon.
   `bridge/latent_bridge.py` implements exactly this. It is **refused by G4** — the anchor
   that would validate it failed as a measurement, not merely failed to conclude. Wiring it
   in would launder an unvalidated branch into a headline number.
2. **Flux redirection through regulation.** E-Flux scales reaction *upper bounds*. An upper
   bound cannot lower a flux that is not already at its ceiling, and the product flux here
   sits ~100× below its ceiling. So E-Flux cannot move this prediction even in principle.

So stress reaching product needs either a validated latent→maintenance anchor, or a
mechanism nobody in this repository has yet proposed. It is not a wiring job.

## Ranked, with what each would actually buy

## The environment axis needs independent evaluation, not a universal null bar

Kocharin provides eleven condition summaries from one reported genotype/study: four
glucose, three ethanol and four mixed-feed states across four dilution rates. The 4.24×
feed contrast is descriptive. Neither eleven rows nor three feed labels supplies the
cultivation/biological-replicate identifiers or randomization scheme needed for inference.

**The old identifiability argument is withdrawn.** This section previously quoted skills
−0.100, −0.153, −0.140 and −0.162 and said a three-parameter model had to clear a roughly
0.08 null threshold. The comparator used the full-data mean, including the held-out
response; the purported null redrew Gaussian covariates as well as shuffling the response.
It was not a null for these fixed candidates and established neither unidentifiability nor
a required sample size. The full historical table remains in
[`EXTERNAL_PRODUCT_VALIDATION.md` §3](EXTERNAL_PRODUCT_VALIDATION.md).

The corrected, not-yet-adopted `scripts/flux_null_baseline.py` candidate uses each fold's
**training mean of log flux**. The constant matches that comparator (skill zero to
roundoff); log(μ), carbon source and their combination still have negative RMSE skills
(approximately −0.0485, −0.0364 and −0.0561). These are source-cohort state-LOO descriptions,
not independent validation. All four candidates remain visible with
`inference_status=pending_exchangeability`; no candidate p-value or percentile is issued.
The five-row random-design demonstration is separate and supplies no replacement bar.

The PHB **content** fit in `scripts/env_to_product.py` answers another question: within this
same cohort its median fold error improves from 1.505× to 1.320×, while the count wrong by
more than 2× increases from one to two. The shipped law is reported separately, not scored
as a held-out refit. For β-carotene, all six fixed-gene requests remain: four answered and
two ethanol setpoint refusals, with missing content and `not-run` layers. Only μ = 0.101 /h
has two answered carbon sources; missing higher-rate ethanol cases do not establish
invariance. None of these facts validates an environment law on new biology.

What remains is a declared target/candidate, traceable independent cultivations and a
justified validation or exchangeability design. If model selection is part of a future
claim, evaluate that selection too; neither beating a constant nor a demonstration
percentile substitutes for it.

1. **A second product.** The generality claim is structural — `spec.py` and `solve.py` name
   no product, and a new one is a TOML file — but structure is not evidence. **Kocharin &
   Nielsen 2013** (PMID 23514405, CC-BY) is the right first test: PHB, eleven chemostat
   steady states, three carbon sources × four dilution rates, **one genotype**, content in
   mg/gDW. One genotype across eleven environments is what makes "flux is a strain constant"
   falsifiable rather than fitted. Needs a PHB spec and one expression number.
2. **A third dilution rate, sited outside [0.101, 0.254].** You chose D = 0.05 or 0.35, and
   that is the right call: a point *between* the anchors moves the 95% width on the growth
   exponent from 0.46 to 0.47, and a point outside halves it. Note the consequence — the
   model currently **refuses** outside that window, so the measurement is informative
   precisely because the model cannot answer there yet.
3. ~~**Regulation into the audit.**~~ **Done, and it buys nothing — which is the finding,
   and it is now measured twice.** The fifty lines were written; the table at the top of
   this file and the section below carry the result. Re-measured 2026-09-05 across **all 25
   stressors, each driven to three times its own EC50**, the number of GEM bounds E-Flux
   tightens is still **zero**. Do not list it again as a route: an upper bound cannot make a
   flux mandatory, and the one direction that could bind — repression — sits on metabolite
   pools that carry no regulon.
4. **A second host or promoter set for the flux law.** Three strains from one paper is one
   experiment. The law reads relative expression, so a different normalisation would need
   its own scalar — whether the *shape* transfers is untested.
5. **Batch mode.** μ_max ≈ 0.4 /h is outside the flux calibration's range, so the model
   refuses there. A shake-flask user asks exactly that question. Needs a batch calibration
   state, not an architecture change.

## What the adversarial pass changed about this list

[`HARD_TESTS.md`](HARD_TESTS.md) has the eight findings. Three of them belong here because
they change what "left to do" means:

- **The model asserts a hard ceiling on β-carotene content that no genotype can pass** —
  1.25 mg/gDCW — **refuted by ~63×** (Arhar 2024, PMID 39215465, 79 mg/gDCW by HPLC on
  gravimetric DCW). The older "17× below a published measurement" framing rested on a
  total-carotenoid absorbance sum this repository has since reclassified as non-specific.
  Closing that gap is now the single
  highest-value experiment available, above the third dilution rate, because it tests a
  claim the model makes about *every strain that will ever be built* rather than about six
  that were measured.
- ~~**The FBA audit is not in the chain**, and this file said it was.~~ **Done.** It is
  reachable through `audit_model=`, opt-in because it costs an 18-second SBML parse and can
  only agree or refute — and that design decision is now stated rather than hidden behind
  either of the two false claims this row has carried, first that it was in the chain and
  then that it was unreachable. `tests/test_everything_is_wired.py` asserts the returned
  content is identical with and without it, using `==`.
- **Phytoene's `passthrough` is refuted**, not merely unverified (Chen 2016, PMID 27329233:
  3.99% of total carotenoid). Fixing it needs a CrtI vmax and Km that nobody has measured
  for these strains, so it is a measurement rather than a code change.

## Standing limits that no amount of code fixes

- **Nothing is Tier 3.** `outputs/registered_prediction_D018.csv` is registered and
  unmeasured. It becomes Tier 3 when someone runs the chemostat and not before.
- **Phytoene is declared `passthrough` and that is an assumption.** Verwaal 2007 reports it
  accumulating. Until it is measured, what the flux layer predicts is the *post-CrtI* flux
  and its relation to true pathway flux is unknown. Declared in the spec so it is greppable
  rather than implicit.
- **crtYB is bifunctional** — one protein running both the synthase step that sets the flux
  and the cyclase step that splits it. The flux law reads crtE and the branch capacity is
  fitted flat across strains whose crtYB copy number is 1:2:3. This entry said "both cannot
  be true of one protein's abundance. Flagged, not modelled" until 2026-09-04, then promoted
  a large F-test p-value to support for equal capacities. That inference is withdrawn:
  these are full-data conditional fit diagnostics, not measured protein abundances or a
  test establishing absence of a dosage effect. Fitting one capacity per strain alongside
  the growth term gives
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=capacity spread across strains" -->1.1494x between the largest and
  smallest across a 3× copy span, non-monotone in copy number, correlating with the measured
  CrtYB transcript at Pearson <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=pearson log mrna vs log capacity" -->+0.0593,
  and the two extra parameters do not survive a nested F test (F(2, 2) =
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F statistic" -->0.5666, p =
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F p" -->0.6383). Non-rejection on six states and two denominator degrees of freedom does not
  establish equal capacities, the sign of a dosage effect or adequate power. See
  `docs/MEASUREMENTS_NEEDED.md` for the proposed dosage series; these conditional fits are
  not the nested forward new-strain validation.
- **Secreted products are out of scope**, refused at load time. That was a decision, not an
  oversight, and it excludes farnesene and resveratrol.
