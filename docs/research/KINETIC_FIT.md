# Fitting the carotenoid branch

`kinetic/carotenoid.py` was parked because "no strain carries the pathway, so nothing here
can be fitted". That reason expired. Elizondo & Saa 2025 (PMID `40891387`,
doi `10.1021/acssynbio.5c00256`) published six chemostat steady states on three
β-carotene CEN.PK2-1c strains with lycopene and β-carotene reported **separately**, and
`docs/research/CALIBRATION_DATA.md` established it as the only sufficient dataset in the
literature. This file is what happened when the branch was actually fitted to it.

## Verdict

**One of the three steps is identifiable and it now carries a fitted, held-out law. The
other two are not identifiable and stay unfitted.** Specifically:

- The **cyclase** step has both its substrate concentration and its rate measured, six
  times. A two-parameter law fitted to it removes **82 %** of the baseline's squared
  out-of-sample error (leave-one-out RMSE in log rate 0.186 against 0.437 for predicting
  the training mean; worst held-out state 33 % off, median 12 %).
- The **synthase** and **desaturase** steps have a pinned flux and an unmeasured substrate,
  so no (vmax, km) pair is determined. `calibrated_kinetics()` raises rather than returning
  four plausible numbers.
- The form the module was written around — a Michaelis-Menten cyclase whose vmax does not
  depend on growth rate — is **refuted**, in all three strains, in the same direction, with
  no fitting involved.
- Flux control lands on the cyclase, which agrees with what Elizondo report for CrtYB. The
  agreement is corroboration and not independent confirmation, for a reason §7 gives.

Everything below is **Tier 1** under [`CLAIM_BOUNDARY.md`](../CLAIM_BOUNDARY.md): a
deterministic fit to real measurements, scored on states it did not see, but computed after
the measurement rather than registered before it. Nothing here is Tier 3.

Reproduce with `python scripts/fit_carotenoid_kinetics.py`.

---

## 1. The extracted table

Six steady states. Rates nmol/gDCW/h, contents nmol/gDCW. Provenance, units and the
cross-check against the printed Table 2 are in
[`data/carotenoid/SOURCE.md`](../../data/carotenoid/SOURCE.md); the numbers came from the
authors' own released pre-processing output, not from a PDF transcription.

| condition | strain | µ (1/h) | q_glucose | q_ethanol | q_lycopene | q_β-carotene | [lycopene] | [β-carotene] | desaturase flux | cyclase share | CrtYB mRNA |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2D01 | β-car2 | 0.101 | −6.68 | 6.40 | 94.4 | 151.3 | 935 | 1498 | 245.7 | 0.616 | 0.328 |
| 2D025 | β-car2 | 0.254 | −14.40 | 16.12 | 84.3 | 184.2 | 331 | 724 | 268.5 | 0.686 | 0.292 |
| 3D01 | β-car3 | 0.101 | −7.13 | 6.01 | 472.0 | 227.1 | 4674 | 2249 | 699.1 | 0.325 | 0.685 |
| 3D025 | β-car3 | 0.254 | −11.72 | 15.84 | 126.5 | 294.8 | 497 | 1159 | 421.3 | 0.700 | 0.420 |
| 4D01 | β-car4 | 0.101 | −5.54 | 2.90 | 760.7 | 185.3 | 7533 | 1835 | 946.0 | 0.196 | 1.000 |
| 4D025 | β-car4 | 0.254 | −14.28 | 14.55 | 392.9 | 453.6 | 1545 | 1784 | 846.5 | 0.536 | 0.620 |

`q_glucose`/`q_ethanol` are mmol/gDCW/h and negative is uptake. The last four columns are
derived here, not reported: contents are `q / µ`, desaturase flux is the two product rates
summed, cyclase share is `q_β-carotene / desaturase flux`, and CrtYB mRNA is relative
expression normalised to β-car4 at µ = 0.101.

## 2. Why the contents are measurements

The single fact that makes this dataset usable is that **carotenoids are not secreted**.
They sit in the membrane and leave the vessel inside the cells, so a chemostat's only
outlet for them is washout with the biomass. At steady state that makes each reported rate
the growth dilution of an intracellular pool:

```
d[lycopene]/dt  = 0 = v_crtI − v_LCY − µ·[lycopene]      ⇒  q_lycopene     = µ·[lycopene]
d[β-car]/dt     = 0 = v_LCY          − µ·[β-carotene]    ⇒  q_β-carotene   = µ·[β-carotene]
```

Three consequences, and the first two are the whole reason a fit is possible:

1. **The lycopene content is a measurement**, `q_lycopene / µ`, with no parameter in it.
   The cyclase therefore has its substrate concentration *and* its rate observed, six
   times — which is a textbook Michaelis-Menten dataset.
2. **The cyclase flux is `q_β-carotene` exactly**, and the desaturase flux is
   `q_lycopene + q_β-carotene`.
3. **The reported lycopene rate is not the desaturase flux.** It is only the part that did
   *not* go on to the product. Reading it as crtI flux understates the pathway by up to 5×
   here — in β-car2 the cyclase carries 62 % of the carbon, so washout is 38 % of what crtI
   made. `BranchSteadyState.lycopene_accumulation` exists to keep the two apart and a test
   asserts they differ.

## 3. Free parameters against observations

The honest count, before any fitting.

| | count |
| --- | ---: |
| steady states | 6 |
| product rates measured per state | 2 (lycopene, β-carotene) |
| **independent observations for the cyclase** | **6** |
| kinetic parameters in the full branch | 6 (three vmax, three km) |
| of those, constrained by any measurement | **2** |

The second product is not a seventh and eighth observation for the cyclase: given the
desaturase flux, `q_lycopene` and `q_β-carotene` sum to it, so one determines the other.
Six observations is the number, and it is what every candidate in §4 is scored against.

**Four of the six parameters are structurally unidentifiable, not merely poorly
constrained.** Phytoene and GGPP were never measured. The desaturase flux is pinned at
`q_lycopene + q_β-carotene`, but for *any* (vmax, km) with vmax above that flux there is a
phytoene pool `p = km·v/(vmax − v)` delivering it exactly. The script tabulates six such
parameter sets spanning three orders of magnitude in vmax; all six reproduce all six states
to machine precision, differing only in a pool nobody measured. Reporting one of them would
be reporting the optimiser's starting point.

So `calibrated_kinetics()` raises `NotImplementedError` naming the two steps and the one
measurement that would lift it: a phytoene peak, which the same HPLC injection already
resolves at 285 nm. A model with six parameters where two are fitted and four are asserted
is not a fitted model, and it must not be shipped looking like one.

## 4. What the data refutes, before what it supports

**A Michaelis-Menten cyclase with a growth-rate-independent vmax is refuted.** No fitting
is involved. MM says the rate is non-decreasing in substrate concentration for *any* vmax
and *any* km. Each strain was run at two dilution rates, which is a paired two-point test of
exactly that with the enzyme complement held as fixed as an experiment can hold it, and all
three go the wrong way:

| strain | [lycopene] ratio (high/low) | q_β-carotene ratio | MM violated |
| --- | ---: | ---: | --- |
| β-car2 | 2.82 | 0.82 | yes |
| β-car3 | 9.40 | 0.77 | yes |
| β-car4 | 4.88 | 0.41 | yes |

Nine times the substrate buys three-quarters of the rate. No parameter choice rescues that.

The refutation survives every attempt to absorb it into something other than growth rate.
Below, leave-one-out over all six states, target `q_β-carotene`, predicted from lycopene
content and growth rate only. Skill is against the better of three baselines — predicting
the training mean rate, the training mean *content* times µ, and the training mean
product ratio — the strongest of which is the mean content at RMSE(log) 0.437.

| model | free params | RMSE(log) train | RMSE(log) LOO | MAPE LOO | worst state | skill |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline — training mean rate | 1 | — | 0.454 | 41.0 % | 78 % | −0.08 |
| baseline — training mean content | 1 | — | **0.437** | 41.9 % | 135 % | 0.00 |
| baseline — training mean ratio | 1 | — | 0.977 | 147.7 % | 537 % | −3.99 |
| vmax shared *(the module's original form)* | 2 | 0.361 | 0.468 | 41.2 % | 68 % | −0.15 |
| vmax ∝ CrtYB mRNA | 2 | 0.488 | 0.586 | 59.9 % | 196 % | −0.79 |
| vmax per strain, growth-independent | 4 | 0.275 | <!-- audit:value table=outputs/carotenoid_fit.csv column="RMSE(log) LOO" row="model=vmax per strain, growth-independent" -->0.550 | 49.4 % | 145 % | <!-- audit:value table=outputs/carotenoid_fit.csv column="skill vs best baseline" row="model=vmax per strain, growth-independent" -->-0.5845 |
| capacity per strain × µ | 4 | 0.085 | <!-- audit:value table=outputs/carotenoid_fit.csv column="RMSE(log) LOO" row="model=capacity per strain x mu" -->0.205 | 17.4 % | 30 % | <!-- audit:value table=outputs/carotenoid_fit.csv column="skill vs best baseline" row="model=capacity per strain x mu" -->+0.78 |
| vmax ∝ mRNA × µ^β | 3 | 0.271 | 0.420 | 33.1 % | 112 % | 0.08 |
| vmax ∝ µ, no saturation | 1 | 0.242 | 0.283 | 25.6 % | 46 % | 0.58 |
| vmax ∝ µ^β | 3 | 0.105 | 0.205 | 20.2 % | 33 % | 0.78 |
| **vmax ∝ µ** | **2** | **0.106** | **0.186** | **17.4 %** | **33 %** | **0.82** |

Every row varies exactly one thing against `vmax shared`, and its name says which. That rule
was added on 2026-09-04 after it turned out to have been broken. The row now called `vmax per
strain, growth-independent` was called `vmax per strain` and varied two things — per-strain
capacities and the loss of the growth term — so its −0.58 was quoted here, in
`docs/MEASUREMENTS_NEEDED.md`, in `docs/WHAT_IS_LEFT.md` and in
`tests/test_kinetic_carotenoid.py` as evidence about GENE DOSAGE, when the factor carrying
the loss was growth-rate-independence. The missing row is the one above it.

Four things this table says that a single fitted number would not:

- **Every growth-rate-independent form loses to the baseline.** Shared vmax, vmax tracking
  measured CrtYB mRNA, and vmax free per strain with no growth term all score negative. The
  last one is four parameters on six points. Its training residual is **sixth of eight**, not
  second: the ordering is 0.0845 (per strain × mu), 0.1054 (mu^beta), 0.1058 (mu), 0.2424
  (mu, no saturation), 0.2706 (mRNA x mu^beta), 0.2751 (per strain, growth-independent),
  0.3611 (shared), 0.4880 (mRNA). Even so it is what over-parameterisation looks like: four
  parameters buying the best fit available to a family that has no predictive skill at all.
- **Per-strain capacity is NOT that family, and the corrected row says so.** Given the growth
  term, three capacities on six points score <!-- audit:value table=outputs/carotenoid_fit.csv column="RMSE(log) LOO" row="model=capacity per strain x mu" -->0.205 held out — above
  the baseline and a hair behind the pooled fit's 0.186. So the reason to pool is not that a
  per-strain capacity fails; it is that the three capacities come out
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=capacity spread across strains" -->1.1494x apart and non-monotone
  across strains whose crtYB copy number is 1:2:3, correlating with the measured CrtYB
  transcript at Pearson <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=pearson log mrna vs log capacity" -->+0.0593,
  and that the extra two parameters do not pay for themselves on a nested F test:
  F(2, 2) = <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F statistic" -->0.5666, p =
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F p" -->0.6383, RSS
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=rss restricted" -->0.067167 against
  <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=rss full" -->0.042874. Two degrees of freedom in the denominator make
  that a weak test, and the weakness is the honest reason to pool.
- **The measured enzyme proxy makes things worse, not better.** Scaling vmax by the RT-qPCR
  CrtYB level is the mechanistically obvious move and it scores −0.79. Using it as an
  authoritative signal was the right instinct and the data rejected it; §6 says why that
  matters more than the score.
- **Freeing the growth exponent does not help.** `vmax ∝ µ^β` fits marginally better in
  sample and generalises worse, and its fitted β lands at 0.91–1.18 across folds. β = 1 is
  what the data wants, so the two-parameter form is the one to keep.
- **Saturation is real.** Dropping km costs 0.24 of skill, so the cyclase is not first-order
  in lycopene over this range.

## 5. The fitted law, with its uncertainty

Winner, in the space it is actually identified in — a saturating relation between the two
**contents**, with no explicit growth-rate term:

```
[β-carotene] = capacity · [lycopene] / (km + [lycopene])        mmol/gDCW
q_β-carotene = µ · [β-carotene]                                 mmol/gDCW/h
```

| parameter | fitted | fit scatter, 95 % | measurement error, 95 % | across the six LOO refits |
| --- | ---: | --- | --- | --- |
| capacity | **2325 nmol/gDCW** (1.25 mg/gDCW) | [1756, 3078] | [2112, 2599] | [2192, 2695] |
| km | **597 nmol/gDCW** (0.32 mg/gDCW) | [296, 1205] | [408, 910] | [426, 786] |

Residual sd in log rate 0.130. The two log parameters correlate +0.85, so they are not
independently pinned — a larger capacity with a larger km fits nearly as well, which is why
the km interval spans a factor of four and why the point value alone would mislead.

**Both numbers are FITTED. Neither is a literature value and neither may be presented as
one.** There was nothing to look up: PubMed returns **zero** records for
`"lycopene cyclase" AND (kcat OR "turnover number")` and zero for
`"phytoene desaturase" AND (…)`; the two hits for `"phytoene synthase" AND (kcat OR
"turnover number" OR "catalytic constant")` are PMID 39322757, a methods review on in vivo
enzymology, and PMID 35699140, whose kcat belongs to a regulator's disulfide reductase
activity and not to the synthase. That is absence of evidence over the searched surface,
not proof of absence — BRENDA's pages are JavaScript-rendered and its API needs
credentials, and SABIO-RK's REST endpoint now returns its single-page app, so neither
database could be queried directly. The corroborating fact is that **Elizondo's own team
hit the same wall**: they sampled a thermodynamically-constrained ensemble (ABC-GRASP)
rather than parameterising from a database, which is what you do when there is nothing to
look up.

**The two uncertainty estimates disagree by about 2×, and the direction matters.** The
asymptotic interval propagates the *scatter of the fit* — how badly a two-parameter law
misses six real points — and the Monte Carlo propagates the *measurement error* Elizondo
reports. Fit scatter is the wider one, which says the residual is model error rather than
measurement error. More replicates would not tighten this; a better functional form, or a
third dilution rate, would.

## 6. The held-out score, and the uncomfortable part

Leave-one-out over all six states. Each prediction comes from a fit that never saw the
state it predicts, and the baseline column is the mean of the same five training states.

| held-out state | measured | predicted | error | baseline |
| --- | ---: | ---: | ---: | ---: |
| 2D01 | 151.3 | 141.6 | −6.4 % | 269.0 |
| 2D025 | 184.2 | 244.1 | +32.5 % | 262.4 |
| 3D01 | 227.1 | 198.6 | −12.5 % | 253.8 |
| 3D025 | 294.8 | 257.3 | −12.7 % | 240.3 |
| 4D01 | 185.3 | 246.5 | +33.0 % | 262.2 |
| 4D025 | 453.6 | 420.7 | −7.3 % | 208.5 |

**Against the best baseline the law wins four of six, not five, and there is no tie.**
The table above scores against *training-mean rate* (LOO RMSE 0.4537). The stronger
baseline is *training-mean content* (0.4372) — and it is the one the headline skill of
0.82 is computed from, so the two must be quoted together. Against it the fit loses 2D01
(fit −6.4 % against baseline +3.5 %) and 4D01 (+33.0 % against −19.2 %), and 3D01 is a
loss rather than a tie: 12.5 % against 11.8 %.

Reporting the weaker baseline in the table while computing the skill from the stronger one
is the error this document exists to avoid, and it was in the first draft of it.

The fitted capacity moves from 2192.1 to 2695.4 across the six leave-one-out refits against
a point estimate of 2325.2 — that is −5.7 % / +15.9 %, not the ±11 % first written, which
was the spread about the range midpoint rather than about the value reported beside it. No
single state carries the result, but the spread is asymmetric and 4D01 pulls it upward.

**Now the part that should not be smoothed over.** A growth-rate-independent content
relation is *arithmetically the same statement* as a cyclase vmax proportional to growth
rate. That is not what a fixed complement of Michaelis-Menten enzyme does, and Elizondo's
own RT-qPCR has CrtYB mRNA **falling** by a third from the low dilution rate to the high one
(β-car4: 1.00 → 0.62; β-car3: 0.68 → 0.42; β-car2: 0.33 → 0.29). The enzyme is going the
wrong way to explain a vmax that rises 2.5-fold. Two candidate mechanisms, and this dataset
does not separate them:

- **A capacity ceiling rather than a rate limit.** `capacity` is a maximum β-carotene
  *content*, ~1.25 mg/gDCW. Carbon arriving beyond it has to stay as lycopene. This reading
  needs no growth-rate-dependent enzyme at all and is why the law is stated in content space.

  **The "plausibly set by how much the membrane will hold" gloss is struck, as of
  2026-09-04.** It was offered here as a live, undecided mechanism while this repository's own
  `data/carotenoid_batch/published_batch_titres.tsv` already closed it. Four rows report an
  HPLC β-carotene content above 1.2483 mg/gDCW in strains with **no storage engineering** —
  Lange 2011 at 3.90, Verwaal 2007 at 5.90, Xie 2014 at 7.41, Arhar 2024 at 79.00 mg/gDCW —
  which is 3.1× to 63.3× the fitted ceiling. A membrane that holds 79 mg/gDCW in one ordinary
  strain is not what limits another to 1.25. The 63× refutation is recorded at
  `pathway/capacity.py` and in `docs/DISTANCE_TO_THE_VISION.md`; what is new here is only that
  it also disposes of this candidate mechanism, which nothing had connected.

  **One caveat travels with that and a fixer must keep it.** All four are batch or stationary
  harvests. There μ → 0, which removes the growth dilution the steady-state law depends on, so
  they do not refute the ceiling *as a steady-state quantity* by simple comparison — they
  refute it as a physical holding capacity, which is exactly the reading being struck.
- **A genuinely growth-coupled vmax**, from something the mRNA does not capture —
  translation, folding, or membrane insertion of a protein whose substrate is in the lipid
  phase.

Relative mRNA is normalised to reference genes and total RNA per gDCW itself rises with
growth rate in yeast, so the *across-growth-rate* mRNA ratios carry a confound the
across-strain ones do not — the mRNA argument against a growth-coupled vmax is weaker than
it looks.

### The exponent is identifiable, and this section used to say it was not

This paragraph previously read "**two dilution rates cannot establish a functional form**;
`vmax ∝ µ` and any other monotone curve through the same pair are indistinguishable here."
That is correct for a fit whose only leverage is the growth-rate axis — two points, two
parameters, nothing left over — and it is the wrong description of this fit.

The calibration set is **three strains** at two dilution rates, and the three strains carry
pathway fluxes spanning 3.5×. The cyclase is Michaelis–Menten, so the map from flux to
product rate is *curved*, and the same µ observed at three flux levels on a curved response
is not the same measurement three times. **The leverage is on the flux axis, not the growth
axis**, which is why it was missed.

Fitted with the exponent free, on all six states:

| | capacity | Km | exponent | RMSE(log) |
| --- | ---: | ---: | ---: | ---: |
| free | 2.458e-3 | 5.345e-4 | **1.041** | 0.0777 |

with three residual degrees of freedom. Every other exponent fits worse on the point
estimate, and not narrowly — µ^0.5 gives 0.170 and µ^1.5 gives 0.124 against 0.078 for µ^1.
Dropping any single state moves the estimate to 0.90–1.12; dropping any whole **strain** —
which removes an entire flux level, the axis the leverage comes from — moves it to
1.01–1.09.

**None of those numbers is a confidence interval, and an earlier version of this section
quoted them as though they were.** A leave-one-out spread measures how stable a fit is under
resampling; each refit still sees most of the data, so the estimates cluster whether or not
the parameter is well determined. The right quantity is the profile likelihood:

| | value |
| --- | ---: |
| point estimate | 1.04 |
| **95% profile interval** (F-test, n=6, k=3) | **[0.53, 1.78]** |
| leave-one-out spread | [0.85, 1.20] |

The interval is **3.6× the resampling spread**. So the honest statement has three parts, and
this file has at various times asserted only one of them:

- the exponent is **bounded** — 0 and 3 are outside, so the data is not indifferent to the
  law, and the original "two dilution rates cannot establish a functional form" is wrong;
- the exponent is **not pinned to 1** — µ^0.6 and µ^1.7 are both inside, so "every other
  exponent is ruled out" is also wrong;
- µ^0.5 both **fits worse by a wide margin** and **sits inside the 95% interval**. Both are
  true. Only the second bears on what has been excluded.

`tests/test_growth_exponent_identifiability.py` pins all of that, including the
counterfactual: collapse the three strains to their mean flux and µ¹ and µ² *do* tie to
within 5%, which is the degenerate case the old paragraph was describing.

**What this does not rescue is the mechanism.** CrtYB mRNA still falls by about a third from
the low dilution rate to the high one, which is the wrong direction for "capacity rises
because there is more enzyme". The exponent is measured; why it is 1 is not settled, and §6
lays out the two candidates.

`calibrated_cyclase()` still **refuses** outside [0.101, 0.254] /h rather than
extrapolating, which is exactly where a batch culture at µ_max ≈ 0.4 /h would want it. That
refusal is about the *range* the law was measured over, which a third dilution rate would
widen; it is no longer about whether the law itself is pinned down.

## 7. Does flux control land on CrtYB?

**Elizondo say yes, and this fit is consistent with it without independently establishing
it.**

What the paper says, quoted from the PMC full text (PMC12455641): "the promiscuous CrtYB
enzyme exerts the highest control over β-carotene production at different growth rates in
the best producer"; "CrtYB is a bottleneck downstream of lycopene in the β-car4 strain, and
increasing its expression is necessary for increasing β-carotene production"; a 67 % CrtYB
upregulation raises the β-carotene sink flux ~146 % and ~84 % at the two dilution rates.
They put entry-flux control on **ERG13**, not on any crt gene, and report that
**upregulating CrtI and CrtE had no relevant effect** — a *decrease* in CrtI is what
redirects flux out of lycopene in the reference strain.

**Our model's structure agrees, on the half of it the data can see.**

- The bifunctional crtYB with a shared enzyme pool is the right shape: "promiscuous CrtYB"
  is Elizondo's own term, and their detailed structure splits it into four sub-reactions
  (`CrtBa, CrtBb, CrtYa, CrtYb`) precisely because one protein runs both domains.
- The **cyclase** domain, not crtI, is the step that saturates. The fitted km, 597
  nmol/gDCW, sits inside the measured range of lycopene contents (331–7533), so the cyclase
  runs between 36 % and 93 % saturated across the six states. A step running far below
  saturation passes on whatever arrives and carries no control; one running near saturation
  does.
- The pellet says the same thing directly. β-car4 at the low dilution rate carries four
  times more lycopene than β-carotene, which is the accumulation signature
  `limiting_step()` already reports as `crtYB_LCY`.

**Where the agreement stops.** Elizondo rank CrtYB against ERG13, CrtI, CrtE, ERG10 and the
rest of the mevalonate pathway using measured transcripts and a model ensemble. This fit
cannot rank anything against the cyclase, because phytoene and GGPP were not measured and
the two upstream steps are unidentifiable (§3). What it establishes is that the cyclase is
partly saturated and that a law built on that predicts held-out states — which is
*consistent with* CrtYB holding control, and is not a second independent measurement of it.
Calling it confirmation would be double-counting one dataset.

One structural caveat travels with this. Elizondo's MCA is reported for **β-car4**, the best
producer, at both dilution rates. Our fit pools all three strains, which the data supports —
but not for the reason this paragraph gave until 2026-09-04.
<!-- audit:retracted a per-strain capacity scores −0.58 out of sample -->
It said "a per-strain capacity
scores −0.58 out of sample", which was the score of a row that ALSO dropped the growth term
(see §4). Given the growth term, a per-strain capacity scores
<!-- audit:value table=outputs/carotenoid_fit.csv column="RMSE(log) LOO" row="model=capacity per strain x mu" -->0.205 held out and beats the baseline. What supports pooling is
that the three fitted capacities sit <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=capacity spread across strains" -->1.1494x
apart against a 3× copy-number span, and that the two extra parameters do not survive a
nested F test (p = <!-- audit:value table=outputs/carotenoid_parsimony.csv column="value" row="quantity=nested F p" -->0.6383). That pooling is a finding
about these three strains, not a claim that the cyclase capacity is strain-independent in
general.

## 8. What this does and does not license

**Licensed.** The product forecast can now come off a kinetic layer rather than off FBA,
which is what D1 and the regulation layer said it had to do. Given a lycopene content and a
growth rate inside the calibrated window, the cyclase step returns a number with a held-out
error of 17 % and a stated interval. `fba/carotenoid.py` remains what it was — the
accounting layer — and this is the layer that forecasts, exactly the split COSMIC-dFBA uses
(doi `10.1016/j.ymben.2024.02.012`).

**Not licensed, and each one names its blocker.**

- **A full-branch forecast from GGPP supply.** Four of six parameters are unfitted. The twin
  cannot yet go GSMM → GGPP → product end to end without asserting them. *Blocker: one
  phytoene HPLC peak, 285 nm, same injection.*
- **Anything outside 0.101–0.254 /h.** Two dilution rates interpolate; they do not establish
  a form. The function refuses. *Blocker: a third dilution rate.*
- **A dynamic forecast.** These are six steady states, not a time course. Whether a law
  fitted to steady states transfers to a fed-batch trajectory is untested — the same open
  question `CALIBRATION_DATA.md` §8 already records.
- **Any claim about this lab's strain.** No strain here carries the pathway. These are
  López-lineage CEN.PK2-1c strains at 1–5 mg/gDCW total carotenoid, not the 79 mg/gDCW
  record holder, and the capacity constant is theirs.
- **A mechanism for the growth-rate term.** §6. The fit is good; the mechanism is not
  established, and the two candidates make different predictions outside the window.

Two facts from the source travel with all of the above: the vessel was probably **not**
carbon-limited despite the paper's description (57 % of the fed glucose is unaccounted for,
and 6.4 mmol/gDCW/h of ethanol at µ = 0.1 is not a carbon-limited phenotype), and no CO₂ or
O₂ was measured. Neither affects the product rates the fit consumes; both bound how far the
strains' physiology may be generalised. `CALIBRATION_DATA.md` §8 carries the arithmetic.

## 9. Files

| file | what it is |
| --- | --- |
| `data/carotenoid/elizondo2025_steady_states.tsv` | the six states, all fluxes, 95 % CIs |
| `data/carotenoid/elizondo2025_relative_mrna.tsv` | RT-qPCR relative expression, 14 genes × 6 conditions |
| `data/carotenoid/SOURCE.md` | provenance, units, the two caveats |
| `scripts/fit_carotenoid_kinetics.py` | refutation, identifiability demonstration, candidates, LOO, uncertainty |
| `src/ystwin/kinetic/carotenoid.py` | `ELIZONDO2025`, `calibrated_cyclase`, the refusals |
| `tests/test_kinetic_carotenoid.py` | 32 tests; the refutation is asserted so it cannot be reverted quietly |
