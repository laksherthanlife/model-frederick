# From zero scored predictions to one

Historical research notes on forecast contracts and exploratory experiments. This is not
the current validation specification. Current scientific conclusions are bound through
`data/current_claims.json` to reviewed reproduction receipts.

The NIS values 1.76/4.44 below came from a different, incompletely recorded research
configuration. Its exact selection, seed and prior manifest were not recovered, so those
values are not substituted for the current diagnostic; `data/current_claims.json` carries
them as `nis.research.shipped`, role `historical`, binding `pending`.

This paragraph also claimed that "the archived 500-particle CLI configuration separately
reproduced 19.79/54.45 at the starting code snapshot". **That is retracted.** No output, no
manifest and no registry entry anywhere in this repository records that pair, so there is
nothing to reproduce it from and no configuration it can be attributed to. The current
diagnostic is the one in `outputs/nis_channel_summary.csv`: 319 well-channels per channel
over 4 plates at K = 22, descriptive medians 45.31 (OD) and 19.32 (mCitrine) against the
band [0.499, 1.672], `tested_well_channels: 0`, verdict INCONCLUSIVE — which is a refusal,
not a magnitude.

<!-- audit:retracted reproduced 19.79/54.45 at the starting code snapshot -->

The repaired CLI uses opening-only priors, mean-preserving physical-state noise, explicit recorded
blanks, scored-read exclusions and particle/ancestor diagnostics. Its default observation
model is now a normalized conditional Gaussian: each particle's expected reading sets its
noise scale, and predictive measurement variance is the prior-weighted mean of the squared
scales. Finite signed readings are retained; positive physical-prior initialization is a
separate requirement. Configuration and population are recorded in `outputs/nis_manifest.json`.
An ESS/diversity refusal remains inconclusive, not a zero NIS or a successful calibration.
The explicitly named legacy observation mode preserves the older outcome-dependent
pseudo-likelihood for historical comparisons, not current calibration claims.

The sections below retain the historical hypotheses and partial experimental descriptions.
Different configurations must not be conflated, and a source citation or a remembered
number is not a complete reproduction record.

---

## The contract, in one sentence

> **ystwin maps one well's raw OD600 and mCitrine readings over the first two hours, plus
> pre-registered reporter and growth priors, to a joint predictive distribution over that
> well's raw readings for the following two hours; it is never given the well's later
> readings, its construct, or its dose; and it exists to decide, before a plate finishes,
> whether a well is still behaving like its cohort.**

Two more ambitious candidates are stated in §1.2. This is the one to commit to, and §1.3
says why. The template that generated it — *inputs → distribution-valued outputs over a
horizon, never given X, supporting decision D* — is §1.1.

## Verdict — the fastest route from zero scored predictions to one

**Run the NIS calibration test on `estimator.py` against the real plate traces.** No ground
truth, no baseline, no rival model, no producing strain, no new plate: the null is a χ²
distribution and the statistic is built from the filter's own innovations. It took under an
hour here and already answers, on all 167 eligible wells. The posterior is **over-confident**
(median time-averaged NIS 1.76 on OD and 4.44 on mCitrine, against an acceptance band of
[0.52, 1.63]) and the innovations are almost never white (median lag-1 autocorrelation
+0.62 and +0.90 against a ±0.40 bound; 4% of wells pass on mCitrine). Run both tests
together and they trap the defect: **no setting of the random-walk scales satisfies both at
once** (§6.5b), so the deterministic dynamics are wrong, not their variance — `_propagate`
gives growth rate no drift, so the propagated mean cannot follow a decelerating culture.
The forecast-ahead task of §0 finds the same missing term independently, as 40–49% coverage
of a nominal 95% interval and a bias-dominated +37% overshoot.

Both are pre-registerable today, and both are losses. A scored loss that names its own
cause is the deliverable.

---

## 0. The headline, before the theory

The audit was right that no scored forward prediction exists. It is also, it turns out,
about forty lines of glue away — the pieces are all in `estimator.py`. So the first thing
this document did was run it.

**Protocol.** For every culture well on the two plates that carry fluorescence in all four
construct blocks, condition `ParticleFilter` on the raw OD600 and raw mCitrine readings for
`t <= 2.0 h` (13 of 25 timepoints), then forecast the remaining 12 timepoints with no
further observations, and score the forecast against the withheld readings. Priors as
stated in `TwinPriors`, with the observation sigmas set to the measured white-noise
component (§6.2). 84 wells, 2 plates, 4342 scored (well, channel, horizon) points.

**Result.**

| horizon | OD RMS rel. err. | best fair baseline | FL RMS rel. err. | best fair baseline |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 h | 16.5% | **7.7%** (causal clim.) | 12.0% | **4.7%** (log-linear) |
| 1.0 h | 26.5% | **12.5%** (causal clim.) | 14.6% | **11.3%** (log-linear) |
| 1.5 h | 38.2% | **16.7%** (causal clim.) | 20.5% | **21.7%** (log-linear) |
| 2.0 h | 53.5% | **21.7%** (causal clim.) | **28.4%** | 36.3% (log-linear) |

"Causal clim." is the leave-one-plate-out shape climatology of §4 — a fair baseline, not the
same-plate oracle. Baselines are deterministic, so their RMS relative error is directly
comparable to the twin's; the distributional comparison is below.

| horizon | OD | FL |
| ---: | ---: | ---: |
| empirical coverage of the nominal 95% forecast interval, 0.5/1.0/1.5/2.0 h | 0.40 / 0.40 / 0.47 / 0.47 | 0.49 / 0.46 / 0.46 / 0.49 |
| \|bias\| / predictive sd at 2.0 h | 2.17 | 1.33 |

**Verdict on the first scored prediction: the twin loses.** Paired per-well against a
six-point log-linear extrapolation — the fairest single comparison, since both see exactly
the same 13 training points from that one well — the baseline wins on OD by 23.0% of level
at 2 h (95% well-bootstrap CI [18.2, 28.3], twin better in 10% of 84 wells) and ties on FL
at 2 h (baseline better by 2.0% of level, CI [−1.6, +5.7], twin better in 48% of wells).
The 95% forecast intervals cover 40–49% of the time instead of 95%.

Two things follow, and they are the reason this is worth having.

**First, the loss is diagnostic, not diffuse.** `|bias|/sd > 1` says the error is
systematic, not noise. The bias is positive at every horizon on both channels (+37% OD,
+14% FL at 2 h) — the twin over-predicts. The cause is in `estimator.py::_propagate`:
growth rate is a driftless random walk, so the forecast carries the last inferred `mu`
forward, and these cultures decelerate hard inside the four-hour run (`mu_max` 0.557
falling to `mu_late` 0.277 on the 2026-07-22 plate — a factor of two). A model with no
deceleration term extrapolating a decelerating culture must overshoot, and does.

**Second, the two channels fail differently, and that is informative.** On OD the twin is
beaten by trivial extrapolation everywhere and by a wide margin: the biomass model, as
written, adds nothing and subtracts a lot. On mCitrine the twin is beaten at short horizon
and pulls level at 1.5–2 h — where the naive log-linear extrapolation starts overshooting
badly (+15.5% bias at 2 h) because an exponential fit cannot know that
`dR/dt = k_syn - (mu + k_deg) R` saturates. That is the reporter ODE doing exactly the work
it was built to do, and it is the one place in the whole pipeline where the mechanism
demonstrably earns its keep. It is not yet enough to win.

**Third, the intervals fail the same way the filter's do.** Coverage of 0.40–0.49 at
nominal 0.95 says the forecast intervals are roughly half the width they need to be. §6.5
tests the *filtering* posterior by an entirely different route — a χ² test on the filter's
own innovations, no forecast and no baseline involved — and finds the same sign:
over-confident, median normalised innovation squared 1.76 on OD and 4.44 on mCitrine
against an acceptance band of [0.52, 1.63]. Two independent statistics, one conclusion.

**And a methodological warning bought for free.** On the aggregate RMS table above, the
twin *wins* on FL at 2 h (28.4% vs 36.3%). On the paired per-well mean-absolute comparison
it *ties*. The disagreement is real: the log-linear baseline has fat tails (a handful of
wells where the exponential fit explodes) which RMS punishes and MAE does not, while the
twin is uniformly mediocre. **Which scoring rule you pick changes the verdict.** That is
the single strongest argument in this document for pre-registering the rule before looking
(§4, §7).

---

## 1. What a digital twin's contract actually is

### 1.1 The template

A twin's contract is the sentence that makes it falsifiable. Five slots, and the third is
the one nobody writes down:

> **`<TWIN>` maps `<INPUTS, fixed at t=0>` to `<OUTPUTS, distribution-valued>` over
> `<HORIZON>`, is never given `<WITHHELD>`, and exists to support `<DECISION>`.**

The `WITHHELD` slot is what turns a fitting exercise into a prediction. Everything else can
be satisfied by a curve fit. The competing team's one-liner — *"fixed environment
e=[I,O,S] → complete product curve P(t); the model never receives an environmental time
series"* — is strong precisely because that clause is explicit and load-bearing: it forbids
the model from being fed the thing that would make the task trivial.

The literature converges on this shape from three directions.

**The provenance direction (Grieves; Glaessgen & Stargel).** Grieves' original framing is
three-part — physical product, virtual product, and the data connection between them — and
Glaessgen & Stargel's NASA/USAF definition is explicitly about *as-built, as-maintained*
state plus a multi-physics simulation used to *predict remaining life*. Both make the
digital twin a claim about a **specific individual asset over its future**, not a claim
about a population. That is the part this project has been getting wrong at the level of
framing: the reported results are population-level audits (pass rates, variance shares,
verdicts), and a twin's contract is per-instance and forward.

**The capability-ladder direction (San, Rasheed & Kvamsdal 2021; originally DNVGL-RP-A204).**
Note the attribution, which is commonly got wrong: the 0-5 ladder is *not* in Rasheed, San &
Kvamsdal's *IEEE Access* 2020 paper -- that paper lists eight "values" of a digital twin and
no capability levels at all. The ladder appears in San, Rasheed & Kvamsdal,
*GAMM-Mitteilungen* 44(4):e202100007, 2021 (San first author), which credits it onward to
DNVGL-RP-A204 (2020). Its levels run standalone
→ descriptive → diagnostic → predictive → prescriptive → autonomous. This project's
real-data output sits squarely at **diagnostic**: G1 pass rates, D2 variance shares, G4
verdicts, dilution-corrected dose responses, plate-map recovery. Excellent diagnostics.
Zero predictive. The ladder is useful because it names the gap without disparaging what
exists — diagnostic is a real rung, it is just not the one the word "twin" is usually
claiming.

**The formal direction (Kapteyn, Pretorius & Willcox).** This is the vocabulary the project
should adopt, because it is a probabilistic graphical model rather than a diagram. Their
abstraction has exactly six quantities — **physical state `S`, digital state `D`,
observational data `O`, control inputs `U`, quantities of interest `Q`, reward `R`** — and
the asset–twin system is a **dynamic decision network** (a dynamic Bayesian network plus
decision nodes), unrolled from `t = 0` through the current time `t_c` to a prediction
horizon `t_p`. The belief state factorises as

```
p(D_0..D_tc, Q_0..Q_tc, R_0..R_tc | o_0..o_tc, u_0..u_tc)
    = Π_{t=0..tc} [ φ_t^update · φ_t^QoI · φ_t^evaluation ]

φ_t^update     = p(D_t | D_{t-1}, U_{t-1} = u_{t-1}, O_t = o_t)
φ_t^QoI        = p(Q_t | D_t)
φ_t^evaluation = p(R_t | D_t, Q_t, U_t = u_t, O_t = o_t)
```

with prediction achieved by extending the same belief state past `t_c` to `t_p`. Verified
against arXiv:2012.05841v2, equations (1)–(4).

Mapped onto this project, the fit is close enough to be worth writing down:

| Kapteyn quantity | this twin |
| --- | --- |
| `S` physical state | true biomass, per-cell reporter, promoter activity, instantaneous μ — **never observed** |
| `D` digital state | `estimator.py`'s four particle dimensions plus `TwinPriors` |
| `O` observational data | raw OD600 and mCitrine, 25 points, 10 min apart |
| `U` control inputs | dose and construct at `t = 0`; nothing thereafter — **a fixed, single-shot control** |
| `Q` quantities of interest | dilution-corrected promoter activity; fold change; (parked: product flux) |
| `R` reward | **absent.** Nothing in the repo evaluates a decision. |

Four consequences, and the last two are the diagnosis.

1. **Digital state is not physical state**, and the paper is explicit that "the digital
   state space will generally be only a subset of the physical state space." `estimator.py`
   already has this right — `promoter_activity` is *a state with a posterior*, and the
   docstring says so. Kapteyn gives that choice a name and a justification.
2. **`φ_t^update` is exactly what `ParticleFilter.update` computes**, and `φ_t^QoI` is
   exactly `reporter.py::promoter_activity`. The repo has implemented two of the three
   factors without naming them.
3. **`U` is degenerate here, and that is the honest reason this is a modest twin.** Dose is
   set once and never revised; there is no control loop. The competing team's contract has
   the same degeneracy (`e = [I, O, S]` fixed at `t = 0`), so this is a property of
   plate-based biology, not a defect peculiar to this project. Say so rather than implying
   a control loop exists.
4. **`R` is missing, and that is a defect.** `run_scenarios.py` produces flux tables and no
   decision. That is why it reads as a report rather than a twin: nothing in it can be
   *acted on*, so nothing in it can be *wrong in a way that costs something*. The
   Candidate A contract in §1.2 supplies a minimal but real `R` — the cost of continuing to
   read a well that has already left the quantitative range.

**The bioprocess direction (hybrid semi-parametric modelling; Oliveira; von Stosch et al.;
Glassey).** This is the school this project actually belongs to — mechanistic skeleton
(the reporter ODE, the dilution term, the GSMM) with data-driven parts bolted where
mechanism is unavailable. Its standard validation idiom is exactly the one missing here:
train on some conditions, predict a *withheld condition*, report the error. Note that
`analysis/transfer.py::leave_one_stressor_out` already implements that shape, on simulated
data.

### 1.2 Three candidate contracts for this twin

**Candidate A — Trace forecaster (the one to commit to).**

> ystwin maps a single well's raw OD600 and mCitrine readings over the first two hours,
> plus the pre-registered reporter and growth priors, to a joint predictive distribution
> over that well's raw OD600 and mCitrine readings for the following two hours; it is never
> given the well's later readings, its construct identity, or its dose; and it exists to
> decide, before a plate finishes, whether a well is behaving like its cohort or has
> already left the quantitative range.

- Inputs: 13 timepoints × 2 channels, one well.
- Outputs: distribution over 12 × 2 withheld readings.
- Withheld: everything after 2 h; the construct label; the dose. **Construct and dose
  withheld is what makes this hard and honest** — the filter is not allowed to know it is
  looking at 5 mM DTT.
- Horizon: 2 h, which is the whole remaining run.
- Decision: early termination / early flagging. Modest, real, and already implicitly
  desired (G1 exists to answer "is this well quantitative?" *after* the fact; this answers
  it during).

**Candidate B — Dose-response interpolator.**

> ystwin maps the dilution-corrected promoter activity measured at a subset of the dose
> ladder for one construct on one plate to a predictive distribution over that construct's
> activity at a held-out dose on the same plate; it is never given the held-out dose's
> readings; and it exists to decide which doses the next plate need not repeat.

Cheaper to state, much weaker: it is interpolation within one plate, and §6.4 shows the
between-plate spread swamps it.

**Candidate C — Cross-construct transfer (aspirational, do not commit).**

> ystwin maps the full dose-response of UPRE1 under DTT to a predictive distribution over
> UPRE2's dose-response under DTT on the same plates; never given UPRE2's readings; to
> decide whether a new promoter needs its own full ladder.

This is the one that would impress. It is also the one the data cannot support: n=2 plates
for UPRE2 and a between-plate fold reproducibility of 8.7% GCV on a corrected effect of
1.1–1.5× (§6.4). Park it.

### 1.3 Which to commit to, and why

**Commit to Candidate A.** Four reasons, in order of weight.

1. **It needs no producing strain, no anchor, no qPCR, no autofluorescence measurement.**
   It predicts the two channels that were actually recorded, in their raw units, which is
   the only quantity in this project with no unmeasured calibration constant standing
   between model and data. `observation.py`'s inner-filter term stays inert;
   `gates/g4_anchor.py` is not consulted; `fba/` stays parked. `kinetic/carotenoid.py` is no longer parked -- it is fitted to six measured chemostat steady states, see docs/research/KINETIC_FIT.md.
2. **It has thousands of scored points and needs no new bench time.** 4342 today, and
   ~8600 if the two UPRE1-only plates are brought in under a background sweep.
3. **It is the purest available test of the one piece of mechanism the project owns** —
   `reporter.py`'s `dR/dt = k_synth - (mu + k_deg) R`. Everything else in the twin is
   either a diagnostic or parked.
4. **It is comparable in kind to the competing team's contract.** Both are "fixed
   information at t=0 → distribution over a future curve, no environmental time series
   fed in". Theirs is a product curve; this one is a reporter curve. That is a difference
   in what strain exists, not a difference in rigour.

The honest caveat, stated in the contract's own terms: Candidate A licenses a claim about
**wells on plates like these**, not about new plates. There are four plates and two of them
are usable, so the unit of independent replication for any "this generalises" claim is 2.
§7.4 says exactly what to do about that.

---

## 2. How to score a probabilistic forward prediction

### 2.1 The rule: CRPS, per channel, on the raw reading

A scoring rule `S(F, y)` takes a predictive distribution `F` and the realised value `y`.
It is **proper** if the forecaster's expected score is optimised by reporting their true
belief, and **strictly proper** if that optimum is unique. Gneiting & Raftery is the
reference that organises this; the practical upshot is short.

- **Log score** `-log f(y)`. Strictly proper, and the right thing when a density is
  available. Wrong here: it is infinite whenever a particle ensemble puts no mass at `y`,
  which will happen (the twin's coverage is 0.47, so 53% of points sit in a tail, and some
  sit outside the ensemble's support entirely). Do not use the log score on a
  1200-particle ensemble at a 2 h horizon.
- **CRPS** `∫ (F(z) - 1{y ≤ z})² dz`. Strictly proper, finite for any `F`, in the units of
  the observable, and it degenerates to the absolute error when `F` is a point mass — which
  means **a deterministic baseline and a probabilistic twin can be scored on one scale
  without dressing the baseline in a fake distribution.** That last property is what makes
  it the right choice for this project, where the baselines are genuinely deterministic.
  For an equally weighted ensemble `x_1..x_n` the exact form is

  ```
  CRPS = (1/n) Σ_i |x_i - y|  -  (1/(2n²)) Σ_i Σ_j |x_i - x_j|
  ```

  which is what the §0 numbers use.
- **Energy score** generalises CRPS to vectors, and the **variogram score** to
  correlation-sensitive multivariate cases.

### 2.2 Trajectory-valued predictions: what to actually do

The prediction here is a *trajectory* — 12 timepoints × 2 channels, strongly
autocorrelated. Three options, and the recommendation is to do the first two and skip the
third.

1. **Per-horizon CRPS, reported as a curve against lead time, then averaged.** The
   standard in weather and epidemic verification. It is *not* a proper scoring rule for the
   joint trajectory — it ignores the dependence between horizons — but it is proper for
   each marginal, it is interpretable, and it localises failure (§0's table is only useful
   because it is per-horizon: it shows the twin's problem grows with lead time, which
   names the missing deceleration term).
2. **Energy score on the whole 24-vector**, as one number that is proper for the joint.
   Report it as the single headline figure so the per-horizon curve cannot be cherry-picked.
   Do not use it as the only number — it is uninterpretable on its own.
3. **A trajectory-level distance** (dynamic time warping, functional depth). Skip: not
   proper, and there is no established convention to compare against.

**Do not average CRPS across channels.** OD is ~0.3 units and mCitrine is ~1200 RFU; the
sum is the fluorescence score with rounding error. Score each channel separately, or
score both after dividing by the channel's own training-period mean level (which is what
the "% of level" columns in §0 do) and say which normalisation was used.

### 2.3 Calibration and sharpness are separate, and calibration comes first

The organising principle (Gneiting, Balabdaoui & Raftery) is **maximise sharpness subject
to calibration**. A proper scoring rule mixes the two into one number, which is what makes
it decision-relevant and also what makes it useless for diagnosis. So report both:

- **PIT histogram.** For each scored point, `u = F(y)` — the predictive CDF evaluated at
  the realisation. Under calibration `u ~ Uniform(0,1)`. For an ensemble the discrete
  analogue is the **rank histogram** (a.k.a. Talagrand diagram): the rank of `y` among the
  `n` ensemble members, uniform over `1..n+1` under calibration. U-shaped means
  under-dispersed (intervals too narrow); dome-shaped means over-dispersed; sloped means
  biased. **This project will get a violent U with a strong slope** — coverage 0.47 at
  nominal 0.95 is the U, and `|bias|/sd = 2.17` on OD is the slope.
- **Coverage of central intervals** at the nominal levels you intend to quote. Reported in
  §0. This is the number to put in front of a judge, because "our 95% interval contains
  the truth 47% of the time" needs no statistical training to evaluate.
- **Sharpness**, unconditional on the outcome: mean predictive interval width. Only
  meaningful once calibrated, and quoting it before then is how over-confident models get
  praised.

### 2.4 Prequential / forward-chaining validation

Dawid's prequential principle: a forecasting system should be assessed only on the
sequence of predictions it actually issued, each conditioned solely on the data available
when it was issued — never on a refit that saw the whole series. Two operational
consequences here, both currently violated somewhere in the repo:

- **The split must be temporal, not random.** Holding out random timepoints from a 10-minute
  trace leaks almost everything, because a neighbour 10 minutes away is nearly the same
  measurement. `estimator.py::update` already refuses a backwards observation and says why
  ("the filter does not run backwards"; "which is what the firewall requires of a
  forecast"). That firewall is correct and has never been exercised on real data.
- **Any hyperparameter chosen by looking at the test region is not a prediction.** The
  smoothing window in `reporter.py::default_activity_window_h` was tuned on simulated
  induction pulses, which is fine and is documented as such. `TwinPriors`' walk scales are
  documented as "stated rather than tuned against the outcome being tested" — which is
  exactly the prequential requirement, written into a docstring, and it means the §0
  numbers are legitimate as they stand. Preserve that. The moment `activity_walk` is tuned
  to improve the forecast score, the score stops being a prediction. If tuning is needed,
  tune on the two UPRE1-only plates and score on the two full plates, and say so.

### 2.5 Visual predictive checks: the right diagnostic idiom for this data shape

VPC is the pharmacometrics standard for exactly this data shape — a population model, many
individuals, each a noisy trajectory. The procedure: simulate many replicates of the whole
dataset from the model, and overlay the simulated 5th/50th/95th percentile bands on the
observed 5th/50th/95th percentiles, binned by time. Mis-specification shows up as the
observed percentile line leaving the simulated band.

**Prediction-corrected VPC** (Bergstrand, Hooker, Wallin & Karlsson) is the variant this
project needs, and the reason is precise: a plain VPC is invalid when individuals differ in
their independent variables, because binning across heterogeneous individuals inflates the
apparent spread. Here the heterogeneity is enormous and structural — 7 doses × 4 constructs
in one plate, and doses 2 and 4 mM H₂O₂ *kill* the culture (§6.3). A plain VPC pooled over
a plate would be uninterpretable. Prediction correction normalises each observation by the
ratio of the bin's typical prediction to that individual's own prediction before pooling,
which is precisely the fix.

**Use VPC as the diagnostic beside the scoring, not instead of it.** It has no test
statistic and no null; it is a picture. Its value here is that it will show *where in the
trace* the model breaks, which the aggregate CRPS cannot. Recommendation: 500–1000
simulation replicates, binned on the 10-minute grid (no binning decision needed — the grid
is already discrete and shared), stratified by construct, prediction-corrected within
dose.

### 2.6 The epidemic forecast hubs: the governance model to copy

The COVID-19 Forecast Hub and CDC FluSight are the best available template for
*pre-registered, scored, multi-model* prediction, and their credibility comes from
governance rather than statistics:

- **A fixed quantile grid.** Every model submits the same quantiles of the same target, so
  no model can win by predicting a different thing.
- **Submission deadlines before the truth exists.** The forecast is timestamped before the
  observation.
- **A frozen truth source, named in advance**, with a stated policy for revisions.
- **A published baseline model** that every submission is scored against.
- **Weighted Interval Score (WIS)** as the metric — a weighted average of interval scores
  across a quantile grid, which approximates CRPS while being computable from quantiles
  alone, and which **decomposes into sharpness (interval width) + an overprediction
  penalty + an underprediction penalty.**

  A naming caution, because the wrong form is near-universal: Bracher et al. describe the
  components as *sharpness* and *calibration*, and the word "dispersion" appears exactly
  once in their paper -- in the PIT-histogram section, not as a component name. The
  familiar three-way "dispersion / underprediction / overprediction" labelling is
  downstream vocabulary from the `scoringutils` R package and later hub papers. Use it if
  you like, but do not attribute it to Bracher et al. 2021.

That decomposition is the single most useful import for this project. §0's finding —
"the error is bias, not spread, and the bias is positive" — is exactly a WIS
overprediction-component finding, and having it fall out of the primary metric rather than
a side calculation is worth the small extra implementation cost.

**Recommendation:** score with CRPS on the ensemble as primary (it is exact, and the
ensemble is available), and additionally report WIS on a fixed quantile grid with its
three-way decomposition. They will agree to within a few percent; the decomposition is
what you are buying.

### 2.7 Skill scores: the baseline-relative convention

Weather verification's convention, which makes a score interpretable:

```
SS = (S - S_ref) / (S_perfect - S_ref)
```

For a negatively-oriented score with `S_perfect = 0` this reduces to
`SS = 1 - S/S_ref`. Positive is better than reference, zero is no better, negative is
worse. Report **CRPSS against the named baseline**, per channel per horizon, with a
bootstrap interval. From §0, the twin's CRPSS against log-linear extrapolation is negative
on OD at every horizon and approximately zero on FL at 2 h.

State the reference in the same breath as the number, always. "CRPSS = 0.31" is
meaningless; "CRPSS = 0.31 against persistence" is a weak claim; "CRPSS = 0.31 against
leave-one-plate-out shape climatology" is a strong one.

### 2.8 Conformal prediction: the escape hatch, and why n=2 nearly closes it

Conformal prediction turns any point predictor into interval predictions with **finite-sample
marginal coverage**, assuming only exchangeability of the calibration and test points — no
correctness of the model at all. Given `n` calibration residuals and target miscoverage
`α`, the interval is `prediction ± q` where `q` is the `⌈(n+1)(1-α)⌉`-th smallest absolute
residual. This is the right tool for a twin whose mechanistic model is known to be wrong,
which is this one.

The arithmetic is the constraint. `⌈(n+1)(1-α)⌉ ≤ n` requires **`n ≥ ⌈1/α⌉ - 1`**: 19
calibration points for 95% coverage, 9 for 90%. Below that the conformal quantile is
`+∞` and the honest interval is the whole line.

So the feasibility depends entirely on what you are willing to call exchangeable:

| calibration unit | n available | valid at 95%? | what the interval means |
| --- | ---: | --- | --- |
| biological replicate (plate) | 2 (4 for UPRE1) | **no** | nothing |
| construct × dose cell | 28 | yes | coverage over cells, one plate |
| well | ~84 per plate | yes | coverage over wells |
| (well, horizon) point | 4342 | yes, but | points are not exchangeable — heavily autocorrelated |

The recommendation: **conformalise at the well level, within horizon.** For each horizon
`h`, calibrate on the wells of the held-out plates and apply to the target plate. n ≈ 84
per horizon is comfortably above 19, wells within a plate are a defensible exchangeability
unit, and doing it per horizon avoids pretending that a 30-minute and a 2-hour forecast
share an error distribution (they do not — 12% vs 28% on FL).

Two refinements worth knowing about and probably not worth implementing yet:

- **Conformalised quantile regression** gives locally adaptive widths rather than one
  constant band, which matters because the error scales with level here.
- **Adaptive conformal inference** and **EnbPI** handle time series and distribution shift
  by updating the effective `α` online. Relevant if this ever runs live during a plate;
  overkill for retrospective scoring of a 4-hour run.

The honest framing for a judge: **conformal gives you intervals with a coverage guarantee
you can state without believing the model.** It does not give you a better mean. §0's
problem is the mean.

---

## 3. Ranked menu of prediction tasks on current data

Ranked by (defensibility × cost) — what to do first, not what is most impressive.

### Rank 1 — Forecast-ahead within a trace (Candidate A)

**Already run in §0.** Condition on `t ≤ 2 h`, forecast to 4 h, score CRPS and WIS per
channel per horizon.

- **Baselines:** persistence; 6-point log-linear extrapolation; leave-one-plate-out shape
  climatology (§4).
- **Scoring:** CRPS from the particle ensemble, per channel, per horizon; energy score on
  the joint 24-vector as the headline; rank histogram and 95% coverage as the calibration
  companion.
- **Power:** not the binding constraint. 4342 scored points, 84 wells, but only **2
  independent plates**. The per-well paired bootstrap in §0 already resolves differences of
  ~2–3% of level; a plate-level claim rests on n=2 and should not be made.
- **Licenses:** "the reporter ODE forecasts an unseen half of a real trace better/worse
  than trivial extrapolation, by this much, with calibrated/uncalibrated intervals." That
  is a genuine predictive claim about a real measurement.
- **What could go wrong:** (a) the 2 h split is a choice, and moving it changes the answer
  — pre-register 2.0 h and report a sensitivity sweep separately, clearly labelled as
  post-hoc; (b) 12 of the 25 timepoints sit in the split's training half and the trace is
  only 4 h, so the longest honest horizon is 2 h — do not extrapolate the finding to
  overnight cultures; (c) the top H₂O₂ doses kill the culture and OD *declines*, which no
  term in `_propagate` can represent — either exclude them with a pre-stated rule or accept
  that they will dominate the error.

### Rank 2 — Filter calibration by NIS on real traces

**Already run in §6.5.** Needs no ground truth at all, which is its whole point.

- **Scoring:** χ² test on the time-averaged normalised innovation squared, plus the
  innovation-whiteness autocorrelation test. Full definitions in §6.
- **Power:** `K` = 25 steps per well × 167 wells across 2 plates. The per-well 95% χ²
  acceptance band at `K` = 25 is [0.52, 1.63]; the observed medians are 1.76 (OD) and 4.44
  (mCitrine), with 32% and 12% of wells inside, and whiteness passing in 28% and 4%. This
  is decisive at the available n — no test in this document is better powered.
- **Licenses:** "the twin's stated posterior uncertainty is/is not honest, tested against
  a χ² null." A falsifiable claim about the uncertainty rather than the mean — and the only
  claim in this document that requires no baseline, because the null is a distribution
  rather than a rival model.
- **What could go wrong:** NIS tests the *whole* stated prior specification, walk scales
  included, so on its own it cannot say *which* prior is wrong. That is why the whiteness
  test must run alongside rather than instead: NIS constrains the width, whiteness
  constrains the structure, and the pair is jointly identifying. §6.5b shows why that
  matters — the two statistics move in opposite directions as the walk scales widen, so no
  scale satisfies both, which localises the defect to the deterministic dynamics. Running
  NIS alone would have yielded a retuned sigma and a false pass.

### Rank 3 — Leave-one-biological-replicate-out on the dose response

Fit on plates `1..k-1`, predict plate `k`'s dilution-corrected fold at each construct × dose.

- **Baselines:** the geometric mean of the training plates' folds; the naive uncorrected
  RFU/OD fold; a growth-only predictor (`mu_max` ratio).
- **Power: bad, and the direction of the answer is close to knowable in advance.** §6.4
  measures the between-plate log-sd of the corrected fold at 0.083 (8.7% GCV) against
  effects of 1.1–1.5×, and — the finding that should stop this task from being the
  headline — **the naive uncorrected fold is more reproducible between plates than the
  corrected one, log-sd 0.037 vs 0.083.** The dilution correction involves a numerical
  derivative and a growth-rate estimate, and `mu_max` itself has a between-plate log-sd of
  0.48 (62% GCV). So on a cross-replicate reproducibility metric the mechanistic
  correction will very likely *lose to the raw ratio it was built to replace*.
- **Licenses:** with n=2 usable plates for three of four constructs, leave-one-out means
  train on one plate and predict the other. That licenses almost nothing, and
  `analysis/uncertainty.py::MIN_PLATES_FOR_INTERVAL = 3` is the module's own correct
  refusal of exactly this.
- **Do it anyway, and publish the loss.** It costs nothing, it is the honest answer to
  "does your correction help?", and the correct rebuttal is already available and is worth
  making explicitly: *reproducibility and unconfoundedness are different properties.* The
  naive fold is a better-conditioned statistic of the same data and a worse estimate of
  promoter activity. Losing a reproducibility contest is not evidence against the
  correction — but pretending the contest was not lost is what an audit would catch.

### Rank 4 — Leave-one-dose-out interpolation

Hold out one interior dose per construct; predict its corrected activity from the rest of
that plate's ladder.

- **Baselines:** linear interpolation in log-dose between the bracketing doses (very
  strong); the construct's plate mean; a fitted Hill curve on the remaining doses.
- **Power:** 5 usable doses per construct (0–1.0 mM for the oxidative pair, since 2 and 4
  mM give negative recovered activity), 4 constructs, 2 plates → 3 interior held-out doses
  × 4 constructs × 2 plates = 24 scored points. Thin but not empty.
- **Licenses:** "the model interpolates the ladder", which is a weak claim because
  log-dose linear interpolation on a smooth monotone curve is nearly unbeatable.
- **Verdict:** low value. The interesting version is Rank 5.

### Rank 5 — Leave-the-top-dose-out extrapolation

Hold out the highest dose; predict it from the rest.

- **This is a genuinely hard test and it will be failed, for a reason worth stating.** The
  DTT dose response is **non-monotone**: corrected activity rises to a peak at 1.0 mM and
  falls at 2.0 and 5.0 mM (UPRE1: 1.29 → 0.93 → 0.85 on plate 1). Extrapolating a
  monotone-looking rise into a fall is not possible without a mechanism for the fall, and
  the twin has none. The oxidative pair is worse: at 2 and 4 mM H₂O₂ the culture dies and
  the recovered activity goes *negative*.
- **Licenses:** if passed, a lot — extrapolation is the claim everyone wants. If failed,
  it licenses a precise statement of the model's domain of validity, which is itself
  publishable and is more useful than a passed interpolation.
- **Recommendation:** run it, expect failure, and report the failure as a **stated
  applicability boundary**: "the twin is calibrated for doses at or below the growth-inhibition
  threshold and refuses to extrapolate past it." That mirrors
  `fba/surrogate.py`'s existing behaviour, which "refuses to extrapolate" — the same
  discipline, applied to dose instead of to interpolation grid.

### Rank 6 — Cross-construct transfer (UPRE1 → UPRE2)

- **Baselines:** predict UPRE2 = UPRE1's own curve; predict UPRE2 = plate mean.
- **Power:** n=2 plates for UPRE2. Between-construct fold differences at 1.0 mM DTT are
  1.39 (UPRE1) vs 1.48 (UPRE2) against a between-plate log-sd of ~0.11. The effect is
  smaller than the noise. **Not detectable.**
- **Recommendation:** do not run. It would produce a number and the number would mean
  nothing. `analysis/transfer.py` already has the machinery; point it at simulated data,
  where it currently lives, and say so.

---

## 4. Baselines — the right null

A prediction with no baseline is a magic trick. The ladder below is ordered by how hard it
is to beat, and **all four rungs must be reported**, because a model that beats persistence
and loses to climatology has learned nothing about biology and a great deal about
arithmetic.

### The ladder

**Rung 0 — Persistence.** Predict the last observed value. Bias-dominated and easy to
beat (§0: −8.9% to −28.6% bias on OD). Its only job is to prove the target has dynamics
worth predicting. It does: 32.3% error at a 2 h horizon on OD.

**Rung 1 — Extrapolation of the same well.** Six-point log-linear fit to the training
window, extrapolated. **This is the baseline that matters**, because it sees exactly the
same information as the twin — one well, first 13 timepoints — and applies no biology at
all. Any claim that the mechanism helps is a claim about this comparison. It is also
strong: 26.7% (OD) and 36.3% (FL) at 2 h, and near the noise floor at 0.5 h on FL (4.7%
against a 3.3% technical CV).

**Rung 2 — Climatology, in two flavours that must not be confused.**
- *Causal / leave-one-plate-out:* mean log-trajectory shape from the **other** plates,
  anchored at the target well's last training value. Fair. 21.7% (OD) and 37.6% (FL) at 2 h.
- *Oracle / same-plate:* shape from the target plate's own 84 wells including future
  timepoints. **Not a baseline — a ceiling.** 21.4% (OD) and 26.9% (FL). It uses future
  information the twin is forbidden, and it is a useful upper bound on what any
  single-well model could achieve. Label it as a ceiling wherever it appears.

**Rung 3 — Permutation and shuffled-label nulls.** For any claim of the form "the
prediction depends on dose/construct", refit with dose labels shuffled *within plate* and
confirm the score collapses. Two cautions specific to this codebase:
- `observation.py`'s docstring already flags the general principle — the inner-filter term
  is "a path from the quantity being predicted back into the channel predicting it, so it
  is real signal that shuffled-channel controls cannot detect." A shuffle test cannot find
  leakage that runs through a physical mechanism. Shuffles bound one failure mode, not all.
- Shuffle *within* plate. Shuffling across plates destroys the plate effect too, and since
  the plate effect is enormous here (§6.4) that would make the null artificially easy.

### The reporting rule

For each task, one table: rows = baselines in ladder order plus the twin plus the oracle
ceiling; columns = horizon; cells = CRPS and CRPSS-against-Rung-1 with a bootstrap
interval. If the twin does not appear between Rung 2 and the ceiling, say so in the first
sentence of the results section rather than the last.

---

## 5. The scoring protocol, specified for implementation

Precise enough to build. Written as a pre-registration document, which is how it should be
committed — **to git, before the scoring script is run.**

### 5.1 Frozen inputs

- **Truth source:** the four `.xlsx` exports resolved by `paths.biosensor_plates()`, read
  with `plate.synergy.read_synergy_kinetic(path, prefer_raw=True)`, channels
  `raw_channel("OD600")` and `raw_channel("mCitrine")`. Record the SHA-256 of each file in
  the pre-registration. **Note the trap in §6.1: the second block of each channel is the
  blank-subtracted copy, not a repeat read. Use `raw_channel`, which already picks
  correctly.**
- **Eligible wells:** wells present in both channels, excluding blanks identified by
  `diagnostics.dilution_confound.detect_blank_wells` or by
  `plate.layout.blanks_for_export(path.name)` where the logbook records them. Blanks are
  excluded, not predicted. On the two plates carrying all four blocks this is 83–84 wells;
  the two UPRE1-only plates have no in-channel blank and are excluded from the primary
  analysis (see §5.6).
- **Time grid:** as exported, 25 points, ~10 min apart, 0.138–4.138 h (OD) and
  0.152–4.152 h (mCitrine). Align with `SynergyRun.aligned(reference)`; do not resample.
- **Split:** `t_split = 2.0 h`. Training = indices with `t <= t_split` (13 points).
  Test = the remaining 12. Fixed in advance, not swept.

### 5.2 Frozen model specification

Every number below is registered before scoring and not touched afterwards.

```
gdcw_per_od        = 0.40 g/L per linearised OD unit     # literature; unidentifiable with gain
optics.gain        = 1.0                                  # sets reporter units; absorbed by R0
optics.background  = mean blank mCitrine reading, per plate
optics.autofluorescence = 0.0                             # NEVER MEASURED - see docs/DATA_INVENTORY.md
optics.inner_filter_coeff = None                          # no producing strain; stays inert
od_blank           = mean blank OD reading, per plate

TwinPriors.biomass            = (X0, 0.05*X0)   where X0 = (OD[0] - od_blank) * gdcw_per_od
TwinPriors.reporter           = (R0, 0.10*R0)   where R0 = (RFU[0] - background) / X0
TwinPriors.promoter_activity  = (0.35*R0, 0.30*0.35*R0)
TwinPriors.growth_rate        = (0.35, 0.10)
TwinPriors.k_deg              = 0.0             # D2: mCitrine loss is pure dilution
TwinPriors.activity_walk      = 0.15            # as shipped, stated not tuned
TwinPriors.growth_walk        = 0.06            # as shipped
TwinPriors.biomass_walk       = 0.02            # as shipped
TwinPriors.od_rel_sigma       = 0.028           # MEASURED white component, section 6.2
TwinPriors.rfu_rel_sigma      = 0.013           # MEASURED white component, section 6.2
n_particles = 5000, seed = 20260826, resample_threshold = 0.5
```

Note `gain × gdcw_per_od` is unidentifiable (`generator/calibrate.py` says so). It does
not matter here: `R0` is derived from the data and scales inversely with `gain`, so the
predicted RFU `background + gain·R·X` is invariant to the split. Say this in the
pre-registration so nobody later mistakes the arbitrary `gain = 1.0` for a fitted value.

### 5.3 The forecast

**One change to `estimator.py` is required and it is a genuine bug fix, not an
accommodation.** `ParticleFilter.forecast(horizon_h)` calls `_propagate(x, horizon_h)`
once. Because `mu` is a random walk with `sd ∝ sqrt(dt)` and biomass integrates it as
`X·exp(mu_next·dt)`, a single step of length `T` gives `Var(log X)` contribution
`σ_μ² T³`, whereas iterating `n` steps of `h = T/n` gives
`σ_μ² h³ · Σ_{k,l=1..n} min(k,l) = σ_μ² T³ (n+1)(2n+1)/(6n²) → σ_μ² T³/3`.
**One step is asymptotically 3× too wide in log-biomass variance.** Verified both ways at
`T = 2 h`, `n = 12` (the data grid), 40 000 particles: the exact formula predicts
`sd(log X) = 0.1079` for the iterated path and `0.1721` for one step; measured 0.1084 and
0.1720. The predictive-sd ratio on biomass is 1.60 (the asymptotic 1.73 for the μ-walk term
alone, diluted by `biomass_walk` and the prior spread on μ). The reporter state inherits
the same error, ratio 1.56.

The irony is worth recording: that excess width partially compensates for the missing
deceleration term, so the shipped one-step forecast would show *better* interval coverage
than the correct iterated one — for entirely the wrong reason, and without touching the
+37% mean bias that dominates the error. Fix the propagation; fix the bias separately.

Required addition:

```python
def forecast_path(self, times_h: Sequence[float]) -> list[Posterior]:
    """Iterated forecast on a given grid, with no further observations."""
    # resample once to an equally weighted ensemble, then step _propagate
    # along the supplied grid, yielding a Posterior at each time.

def forecast_ensemble(self, times_h: Sequence[float]) -> np.ndarray:
    """The raw particle ensemble at each forecast time: shape (len(times), n_particles, 4).

    Mean and std are not enough. CRPS, WIS and the rank histogram all need the
    predictive distribution, and a Gaussian summary of a bounded, right-skewed
    ensemble is a different forecast from the one the filter actually made.
    """
```

Then, for each forecast time, map the state ensemble through the observation model —
`observe_od(b, 0.0, calibration, gdcw_per_od)` and
`observe_rfu(r, b, 0.0, optics)` — to get the **predictive ensemble in reading units**.
Score in reading units. Never in state units: the state carries the unidentifiable
`gain × gdcw_per_od` and the reading does not.

### 5.4 The scores

For each (well, channel, horizon) with predictive ensemble `x_1..x_n` (equally weighted)
and realisation `y`:

1. **CRPS**, exact ensemble form of §2.1.
2. **WIS** on the quantile grid `α/2` for
   `α ∈ {0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90}` plus the
   median — i.e. the 23-quantile grid the COVID Hub used, taken from the ensemble's
   empirical quantiles. Report the sharpness / overprediction / underprediction
   decomposition (see the naming caution above).
3. **Absolute error of the predictive median**, so the point-forecast and the
   distributional verdict can be compared. §0 shows they disagree; that disagreement must
   be visible rather than resolved by choosing one.
4. **PIT value** `u = (#{x_i < y} + 0.5·#{x_i = y}) / n`, for the histogram.
5. **Interval indicators** at 50%, 80% and 95% central levels, for the coverage table.

Aggregate: mean CRPS per (channel, horizon); mean over horizons as the headline per
channel; **energy score** on the joint 24-vector per well as the single joint number.

### 5.5 Baselines, precisely

- **B0 persistence:** `pred(t) = y(t_split)`.
- **B1 log-linear-6:** OLS of `log y` on `t` over the last 6 training points;
  `pred(t) = exp(intercept + slope·t)`.
- **B2 climatology-causal:** `shape(t) = mean over other plates of mean over wells of log y(t)`;
  `pred(t) = y(t_split)·exp(shape(t) - shape(t_split))`.
- **B3 climatology-oracle (CEILING, label as such):** as B2 but `shape` from the target
  plate's own wells.

All four are deterministic. **CRPS of a point forecast is its absolute error** — no
dressing required, which is the property that makes CRPS the right primary rule here.
Optionally also report a *dressed* variant of B1, with a Gaussian whose sd is the
horizon-specific RMS training residual, to check the twin is not merely losing to a point
forecast on a sharpness technicality. It is not: §0's coverage of 0.47 is a
calibration failure, not a sharpness one.

### 5.6 Uncertainty on the score difference — and the honest denominator

This is where most of the defensibility lives.

- **Primary:** paired per-well difference in mean score, bootstrapped over **wells**
  (n ≈ 84), stratified by plate. Licenses: "on plates like these, per well." §0 uses this.
- **Secondary:** cluster bootstrap over **plates** (n = 2). Report it and report that
  n = 2 makes it uninformative, exactly as `analysis/uncertainty.py` refuses an interval
  below `MIN_PLATES_FOR_INTERVAL = 3`. **Do not compute a plate-level interval at n = 2.**
  The module already contains the correct refusal and the reasoning behind it ("two
  clusters admit three distinct resamples, of which two are degenerate"); reuse it rather
  than making an exception for the flagship result.
- **Never bootstrap over (well, horizon) points.** 4342 points, but a well's 12 horizons
  are one realisation of one trajectory. Bootstrapping points would divide the true
  standard error by ~√12 and produce the tightest, most wrong interval in the report.

The two UPRE1-only plates (2026-07-28, 2026-08-04) have no blank inside the fluorescence
channel, so `optics.background` is unmeasured on them. They can be added under the
bounded-background sweep `0 ≤ b ≤ min(F)` that `docs/DATA_INVENTORY.md` proposes, reported
as a range of scores rather than a score. That raises the plate count for UPRE1 to 4 —
still below `MIN_PLATES_FOR_INTERVAL + 1`, but it is the only route to n > 2 that costs no
bench time.

### 5.7 Pre-registration checklist

Commit before running the scoring script:

- [ ] SHA-256 of every input export.
- [ ] `t_split`, the horizon set, the eligible-well rule.
- [ ] Every value in §5.2, with a source for each (measured / literature / stated).
- [ ] Primary scoring rule, primary channel, primary horizon, primary baseline.
- [ ] The bootstrap unit and the claim it licenses.
- [ ] The pass criterion, as a number: e.g. *"the twin passes if its mean CRPS on mCitrine
      at the 2 h horizon is below B1's, with the paired per-well 95% bootstrap interval
      excluding zero."* Write it down and then let it fail if it fails. On §0's numbers it
      fails (CI [−1.6, +5.7] straddles zero).
- [ ] The `git rev-parse HEAD` of the code that will produce the numbers.

---

## 6. Filter calibration: NIS and NEES, and the test to run on `estimator.py`

This is the section to implement first, because it needs no ground truth, no baseline, no
new measurement and no producing strain — and because it has already produced a decisive
answer (§6.5).

### 6.1 Why NIS and not NEES

Both come from target tracking, where the question "is my filter lying about its own
uncertainty?" is operationally urgent.

**NEES — Normalized Estimation Error Squared.** With true state `x_k`, posterior mean
`x̂_k` and posterior covariance `P_k`:

```
ε_x(k) = (x_k - x̂_k)ᵀ P_k⁻¹ (x_k - x̂_k)
```

If the filter is consistent, `ε_x(k) ~ χ²_n` with `n = dim(x)`, so `E[ε_x] = n`. Averaged
over `N` independent runs, `ε̄_x = (1/N) Σ ε_x^(i)` and `N·ε̄_x ~ χ²_{nN}`, giving the
two-sided acceptance region `[χ²_{nN,α/2}/N, χ²_{nN,1-α/2}/N]`. The normalised form
`ANEES = ε̄_x / n` has expectation 1 under consistency, which is the convention worth
adopting because it makes the target value 1 regardless of state dimension.

**NEES is unavailable here.** It requires the true state, and the true state is biomass,
per-cell mature reporter concentration, promoter activity and instantaneous growth rate —
none of which is measured on these plates. NEES is available only against the synthetic
generator, where truth is known. **That is a legitimate use and should be run** — it
validates the estimator's implementation before the estimator's honesty is tested on real
data — but it is a simulation result, which is what the audit correctly objected to.

**NIS — Normalized Innovation Squared.** With measurement `z_k`, predicted measurement
`ẑ_k|k-1`, innovation `z̃_k = z_k - ẑ_k|k-1` and innovation covariance `S_k`:

```
ε_z(k) = z̃_kᵀ S_k⁻¹ z̃_k ,     S_k = H_k P_k|k-1 H_kᵀ + R_k
```

Under consistency `ε_z(k) ~ χ²_m` with `m = dim(z)`. **NIS needs no ground truth** — only
the measurement and the filter's own prediction of it. That is the property that makes it
the right test for this project. Time-averaged over `K` steps,
`ε̄_z = (1/K) Σ_k ε_z(k)`, and `K·ε̄_z ~ χ²_{mK}` if the innovations are independent, so
the acceptance band is `[χ²_{mK,α/2}/K, χ²_{mK,1-α/2}/K]`.

Reading the result:
- `ε̄_z / m ≈ 1` — the claimed uncertainty is honest.
- `> 1` — over-confident; `S` too small; intervals too narrow.
- `< 1` — over-dispersed; `S` too large; intervals too wide.

**The whiteness test, which is not optional.** Time-averaged NIS assumes the innovations
are serially independent — that is where the `χ²_{mK}` comes from. So the whiteness of the
standardised innovations `s_k = z̃_k / √S_k` must be tested alongside, not instead:

```
ρ̂(1) = corr(s_k, s_{k+1}),     |ρ̂(1)| ≲ 2/√K under whiteness
```

A filter can pass NIS and fail whiteness, and that combination means the *structure* is
wrong even though the *width* happens to be right. §6.5 finds exactly the inverse-ish case,
and the whiteness failure is the more robust one.

### 6.2 Adapting `S_k` to a particle filter

There is no `H` and no `P` in `estimator.py` — the filter is nonlinear and
sample-based. The ensemble analogue is exact and easy. At each update, **before**
reweighting, with propagated particles `x_i` and *prior* weights `w_i`:

```
h_i     = h(x_i)                        # observe_od(...) or observe_rfu(...)
ẑ       = Σ_i w_i h_i                    # predictive mean of the measurement
P_z     = Σ_i w_i (h_i - ẑ)²             # state-induced predictive variance
R       = (rel_sigma · |z|)²             # the filter's own measurement variance
S       = P_z + R
ε_z     = (z - ẑ)² / S                   # dof = 1 per scalar channel
```

For the two channels jointly, `m = 2`; treat them separately when the aim is to find out
*which* channel is mis-specified, which it is here.

Two subtleties worth writing into the docstring.

- `_log_likelihood` sets `sigma = rel_sigma · max(|obs|, 1e-6)` — scaled by the
  **observed** value, not the predicted one. Harmless for filtering, mildly circular for
  NIS (`R` depends on `z`). It biases `S` upward when `z` overshoots, which flatters the
  test. Prefer `sigma = rel_sigma · ẑ` for the NIS computation and note the change, or
  report both. The §6.5 numbers use the observed-value form, i.e. the flattering one, and
  the filter fails anyway.
- Particle-filter consistency has a second failure mode with no Kalman analogue:
  **degeneracy**. If `ESS` collapses, `P_z` is computed from a handful of distinct
  particles and is underestimated, so NIS is inflated for a reason that has nothing to do
  with model correctness. `Posterior.effective_sample_size` already exists; report ESS
  alongside every NIS value and drop steps below a pre-registered floor (`ESS < 100`).
  The classical ESS-based degeneracy criterion for sequential importance sampling is the
  right citation for that floor.

### 6.3 Setting `R` honestly — which is the precondition for the whole test

NIS is only a test of the *model* if `R` is a measurement of the *instrument*. Otherwise it
is a test of a guess. `TwinPriors` ships `od_rel_sigma = 0.02` and `rfu_rel_sigma = 0.03`,
neither documented as measured. So they were measured.

**A trap first, because it looks like the answer and is not.** Each export contains two
blocks per channel (`OD600[1]`/`OD600[2]`, `mCitrine[1]`/`mCitrine[2]`), which look like
independent repeat reads and would give a clean instrument-noise estimate. **They are not.**
Block 2 is block 1 minus a constant: OD `0.257→0.160, 0.263→0.166, 0.278→0.180` (−0.097
throughout), mCitrine `870→461, 895→484, 942→530` (−409 to −412, the variation being
integer rounding). `KineticBlock.blank_subtracted` is `False` for the first and `True` for
the second, and `SynergyRun.raw_channel` already picks correctly — but any hand-rolled
"repeat read" noise estimate from these would report the reader as ~10× more precise than
it is. **There are no independent repeat reads in this dataset.**

What is available is the three technical replicates per (construct, dose), and they
decompose usefully. On the log scale, splitting the n=3 spread into a persistent per-well
multiplicative offset, the shared time course, and a white residual, over all four plates
and all construct × dose cells:

| channel | persistent per-well sd | white residual sd | total sd | persistent share of variance |
| --- | ---: | ---: | ---: | ---: |
| OD600 | 0.039 | **0.028** | 0.052 | 66% |
| mCitrine | 0.033 | **0.013** | 0.032 | **87%** |

*Computed from all four exports via `plate.synergy.read_synergy_kinetic` and
`plate.layout.NEWPROTOCOL_LAYOUT`; 3 wells × 25 times per cell.*

Three conclusions, and the third is the important one.

1. **The correct `R` for a single-well filter is the white component:**
   `od_rel_sigma ≈ 0.028`, `rfu_rel_sigma ≈ 0.013`. The shipped `0.02` is slightly too
   tight on OD; the shipped `0.03` is **2.3× too wide** on mCitrine. §6.5 shows that
   substituting the honest values makes the filter score *worse*, because the shipped
   `0.03` was quietly absorbing dynamic model error as instrument noise. That is the
   argument for measuring `R`: not to improve the filter, but to stop it hiding.
2. **If the filter is run on the mean of 3 wells**, `R` is the white component over 3 plus
   the *shared* component, which does not average down. That is `analysis/uncertainty.py`'s
   "the unit of replication is the plate, not the well" argument, one level down. Filter
   single wells; it is cleaner and there are 84 of them.
3. **87% of the mCitrine technical spread is a persistent per-well offset, and the
   observation model has no term for it.** `observation.py` writes
   `RFU = background + [autofluorescence + gain·R]·biomass·exp(-eps·carotenoid)` — every
   term shared across wells. A real per-well multiplicative offset (pipetting volume,
   optical path, cell-density-independent well-to-well gain) is therefore unmodelled.
   *This looks like it should be the cause of the calibration failure in §6.5, and it is
   not* — §6.5b prototypes the offset as a nuisance state and the whiteness violation does
   not budge, because the filter's free biomass and reporter states already absorb a
   constant. The decomposition is a correct statement about the data and a wrong inference
   about the model. Kept here because it is the right `R`, and because the failed inference
   is instructive.

### 6.4 What the plate-level spread implies for every cross-replicate task

Same data, computed for context because it decides Ranks 3–6 above. Between-plate log-sd
of the fold change over each construct's own 0 mM control, across the two plates with all
four constructs:

| quantity | median between-plate log-sd | geometric CV |
| --- | ---: | ---: |
| naive uncorrected RFU/OD fold | **0.037** | 3.8% |
| dilution-corrected activity fold | 0.083 | 8.7% |
| `mu_max` ratio | 0.481 | 61.8% |

*From `outputs/sensor_characterisation.csv`, folds recomputed within plate then compared.*

The corrected fold is **2.2× less reproducible in log-sd (5× in variance)** than the naive
ratio it was built to replace, and `mu_max` — an input to the correction — varies by 62%
GCV between plates. Individual cells are worse than the median suggests: AlteredYap1 at
1.0 mM H₂O₂ gives folds of **1.99 and 0.87** on the two plates. Those two replicates
disagree about the *sign* of the effect.

This is not an argument against the dilution correction, and it must not be reported as
one. The correction targets an unconfounded estimate of promoter activity; the naive ratio
targets a reproducible summary statistic. A numerical derivative plus a growth-rate
estimate is inherently noisier than a raw ratio, and that is the price of removing the
`1/μ` gain. But it does mean **any prediction task scored on cross-replicate reproducibility
of the corrected quantity is a task the twin is set up to lose**, and it is why §3 ranks
forecast-ahead above leave-one-replicate-out.

It also means the top two H₂O₂ doses are outside the model's domain entirely: at 2 and 4 mM
the recovered activity is **negative** (−38, −53, −6.8, −172 …) because the culture dies,
OD declines, and the inversion `dR/dt + (μ + k_deg)R` has no branch for it. Pre-register
their exclusion, or pre-register that they are included and will dominate.

### 6.5 The test, run

The above was all preparation. Here is the result on **the full eligible population**: both
plates that carry an in-channel blank, every culture well, `K = 25` steps, both channels.
167 wells. Acceptance band for the time-averaged NIS at `K = 25`, `m = 1`, `α = 0.05`:
**[0.52, 1.63]**. Whiteness bound `2/√25 = 0.40`. Minimum ESS over all wells and steps was
1000 of 2000, i.e. the resampling floor — **no degeneracy anywhere**, so §6.2's caveat does
not bite and every well counts.

*A methodological note that is itself a finding: an earlier version of this section reported
the opposite sign, from a 36-well subset. `cult[:36]` is not a random sample of a plate —
it is the low-dose rows, where the model is nearly right. The subset said "over-dispersed";
the population says "over-confident". **Never characterise a filter on a well-ordered
subset of a plate.** The numbers below are the population.*

| configuration | OD NIS | OD in-band | OD ac | OD white | mCitrine NIS | mCitrine in-band | mCitrine ac | mCitrine white |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **shipped sigmas (0.020/0.030), shipped walks** | **1.76** | 32% | **+0.62** | 28% | **4.44** | 12% | **+0.90** | **4%** |
| measured white sigmas (0.028/0.013), shipped walks | 3.17 | 22% | +0.74 | 17% | 5.35 | 14% | +0.84 | 14% |
| + per-well offset state (sd 0.039/0.033) | 2.73 | 23% | +0.76 | 17% | 4.90 | 16% | +0.84 | 14% |
| walks ×3 (`gw` 0.18 / `aw` 0.45) | 0.65 | 18% | +0.48 | 41% | 1.76 | 20% | +0.52 | 37% |
| walks ×8 (`gw` 0.48 / `aw` 1.20 / `bw` 0.10) | **0.32** | 17% | **+0.06** | **75%** | **0.47** | 32% | **+0.15** | **67%** |

*"in-band" = fraction of wells whose time-averaged NIS lies in [0.52, 1.63]; "ac" = median
lag-1 autocorrelation of standardised innovations; "white" = fraction of wells with
|ac| < 0.40.*

**The verdict.**

1. **With the shipped priors the filter is over-confident, on both channels.** Median
   time-averaged NIS 1.76 (OD) and 4.44 (mCitrine), both above the band's upper bound of
   1.63, with only 32% and 12% of wells inside. The filter's stated posterior is narrower
   than its own innovations justify — by a factor of about 4 in variance on the reporter
   channel. **This agrees in direction with §0's forecast finding** (nominal 95% intervals
   covering 40–49%), which is a useful internal consistency check: filtering and
   forecasting fail the same way.
2. **The failure is grossly heterogeneous across the plate, and the heterogeneity is
   biological.** Median mCitrine NIS by construct, at the measured-white setting:
   **AlteredYap1 58.8, NativeYap1 30.7, UPRE1 2.15, UPRE2 2.22.** The oxidative pair is an
   order of magnitude worse than the ER pair, and within each the failure grows with dose
   (mCitrine NIS by dose row: 2.1, 4.4, 6.7, 8.1, 4.9, 42.4, 24.2). The two top H₂O₂ doses
   are where the culture is dying — the regime §6.4 shows produces negative recovered
   activity. **Stratify by construct and dose, or report a median that is really a
   statement about H₂O₂ toxicity.**
3. **The whiteness test fails everywhere and under every configuration that keeps NIS
   near 1.** Median lag-1 autocorrelation +0.62 to +0.90 against a ±0.40 bound; on
   mCitrine with the shipped priors only **4% of wells** have white innovations.
4. **Using the correctly measured `R` makes NIS worse, and that is the right outcome.**
   Moving from the shipped `rfu_rel_sigma = 0.030` to the measured white value 0.013 pushes
   the mCitrine NIS from 4.44 to 5.35. The shipped value was partly absorbing model error
   as though it were instrument noise. Measuring `R` honestly does not improve the filter;
   it stops the filter from hiding. That is the whole reason §6.3 was necessary.

### 6.5b Two fixes tested, both rejected — and what that identifies

The same table, read as an experiment on two candidate fixes.

**Fix 1 rejected: the per-well offset state changes nothing.** Suggested by §6.3's finding
that 87% of the fluorescence technical variance is a persistent per-well offset. Prototyped
as an extra nuisance dimension (log-normal, sd 0.039 OD / 0.033 mCitrine, zero walk,
carried through resampling). Median NIS moves 3.17 → 2.73 and 5.35 → 4.90; the
autocorrelation moves **+0.74 → +0.76 and +0.84 → +0.84**, i.e. not at all. In hindsight
the reason should have been anticipated: the filter's biomass and reporter states are free
and initialised from that well's own first reading, so a *constant* multiplicative per-well
gain is already absorbed by the existing state. A constant cannot produce a mismatch that
grows with time, and the mismatch here grows with time. §6.3's decomposition is correct
about the data and wrong about the cause.

**Fix 2 rejected: no random-walk scale satisfies both tests.** Widening the walks is the
other obvious knob, and it works on exactly one of the two statistics at a time:

- shipped walks — NIS **too high** (1.76 / 4.44) *and* innovations strongly correlated
  (+0.62 / +0.90). Failing both.
- walks ×3 — NIS lands near the band (0.65 / 1.76), autocorrelation still far outside it
  (+0.48 / +0.52). Whiteness in only 41% / 37% of wells.
- walks ×8 — autocorrelation collapses into the bound (+0.06 / +0.15, white in 75% / 67% of
  wells) and NIS falls **below** the band (0.32 / 0.47). The posterior is now so wide it
  claims almost nothing.

The in-band fraction never exceeds 32%, at any setting, on either channel.

**That trade-off is the identification.** A one-dimensional scale parameter cannot fix a
two-dimensional failure, which means **the deficiency is in the deterministic part of the
dynamics, not the stochastic part.** `_propagate` gives growth rate no drift at all — `mu`
is a driftless random walk — so the propagated *mean* cannot follow a decelerating culture,
and the only way to get the observations inside the predictive distribution is to inflate
its width until the mean stops mattering.

**That trade-off is the identification.** A one-dimensional scale parameter cannot fix a
two-dimensional failure, which means **the deficiency is in the deterministic part of the
dynamics, not the stochastic part.** `_propagate` gives growth rate no drift at all — `mu`
is a driftless random walk — so the propagated *mean* cannot follow a decelerating culture,
and the only way to keep the observations inside the predictive distribution is to inflate
its width until the mean stops mattering.

**And this is the same defect §0 found from the other end.** The forecast overshoots by
+37% on OD with `|bias|/sd = 2.17`: bias-dominated, growing with lead time, positive.
Independent test, independent statistic, same missing term. That convergence is worth more
than either result alone, and it names the next piece of work precisely:

> **Give `mu` a deterministic deceleration term.** A substrate-depletion or logistic
> capacity term in `_propagate` — `dmu/dt` driven by an explicit resource state, or a
> mean-reverting `mu` with a declining target — so that the propagated mean tracks the
> culture and the walk carries only genuine uncertainty. `generator/culture.py` already
> contains growth inhibition and `nutrient_factor`, so the mechanism exists in this
> codebase on the simulation side and has simply never been given to the estimator.

After that change, re-run this table. The claim "this twin's stated uncertainty is honest"
becomes available exactly when the NIS medians sit in [0.52, 1.63] **and** the
autocorrelation sits inside ±0.40 at the same time. Neither test alone would have found
this; running one of them and stopping would have produced a tuned sigma and a false pass.

### 6.6 Implementation, concretely

The §6.5 numbers were produced by reaching into `ParticleFilter` internals (`_propagate`,
`_log_likelihood`, `_w`, `_x`) because the pre-update predictive ensemble is not exposed.
That is fine for a probe and unacceptable for a reported result. The minimal, honest
addition to `estimator.py`:

```python
@dataclass(frozen=True)
class Innovation:
    """One channel's innovation and the filter's own claim about its spread."""
    time_h: float
    channel: str
    observed: float
    predicted: float          # z-hat, the prior predictive mean
    state_variance: float     # P_z, from the propagated ensemble
    measurement_variance: float  # R
    effective_sample_size: float

    @property
    def total_variance(self) -> float:
        return self.state_variance + self.measurement_variance

    @property
    def nis(self) -> float:
        """Normalized innovation squared. One degree of freedom."""
        return (self.observed - self.predicted) ** 2 / self.total_variance

    @property
    def standardised(self) -> float:
        return (self.observed - self.predicted) / math.sqrt(self.total_variance)
```

`ParticleFilter.update` returns `Posterior`; give it a sibling
`update_with_innovations(obs) -> tuple[Posterior, list[Innovation]]`, or accumulate onto
`self.innovations`. The computation happens between `_propagate` and the reweighting in
`update`, reusing the same `observe_od` / `observe_rfu` calls `_log_likelihood` already
makes — so it costs one extra pass over the particles per channel and no extra observation
evaluations if the predicted values are cached.

Then a script — call it `scripts/run_calibration_nis.py`, beside the existing
`scripts/run_calibration.py` — that:

1. resolves the plates via `paths.biosensor_plates()` and skips cleanly when absent, as the
   other scripts do;
2. runs one filter per eligible well over all 25 timepoints;
3. writes per-(well, channel, step) innovations to `outputs/nis_innovations.csv`;
4. reports, per channel: median and IQR of time-averaged NIS, the χ² band for the actual
   `K` used, the fraction of wells inside it, the median lag-1 autocorrelation with its
   `2/√K` bound, and the ESS distribution;
5. **and refuses to print a verdict for any well whose minimum ESS fell below the
   pre-registered floor**, in the style of `g4_anchor.py`'s INCONCLUSIVE — a filter that
   degenerated has not been tested, and saying so is the whole house style.

Companion run on the generator, where truth exists, for **ANEES**: `generate_plate` returns
`plate.truth`, so `ε_x = (x - x̂)ᵀ P⁻¹ (x - x̂)` is computable with the sample covariance of
the weighted particles. Report `ANEES = ε̄_x / 4` with the `χ²_{4N}/(4N)` band. This
validates the estimator against a known truth; the real-data NIS validates its honesty.
Report both, and never let the simulated one stand in for the real one.

---

## 7. What this all licenses, and what it does not

**Available today, no new data:**
- The NIS/whiteness verdict of §6.5–6.5b. Falsifiable, scored against a χ² null, no
  baseline and no truth required, and it identifies its own defect rather than merely
  reporting a failure. **Strongest single result in this document.**
- The forecast-ahead scores of §0. A scored forward prediction with a named baseline, which
  the twin currently loses on OD and ties on mCitrine at 2 h.
- The reproducibility ladder of §6.4, which is what makes Ranks 3–6 honest.
- The one concrete piece of work both tests name: a deterministic deceleration term for
  `mu` in `estimator.py::_propagate`, borrowed from `generator/culture.py`, which already
  has growth inhibition on the simulation side. After that, both tests re-run unchanged and
  the pass criteria in §5.7 and §6.5b decide it. That is a closed loop, which is what the
  project has been missing.
- All computable from `outputs/` and the four exports in an afternoon of glue.

**Available for one plate of bench time:** the BY4741 reporter-free strain read in the
mCitrine channel. `docs/DATA_INVENTORY.md` already argues for it on other grounds.
`autofluorescence` currently enters the pre-registration as `0.0`, and it is the only
frozen input in §5.2 with no measurement and no literature value behind it.

**Available for two plates:** n = 4 for all four constructs, which crosses
`MIN_PLATES_FOR_INTERVAL` and makes every cross-replicate task in §3 estimable rather than
refused. Combined with reading fluorescence on all 12 columns and putting blanks inside
every channel, this is the difference between a report and a result.

**Not available and should not be claimed:** anything about the product. No carotenoid
pathway exists in these strains, so `fba/fva.py`, `kinetic/carotenoid.py`,
`observation.py`'s inner-filter term and the entire latent branch of `run_scenarios.py`
stay parked. The right thing to say about them is what `docs/PARKED.md` and the module
docstrings already say — designed ahead of the strain, constrained by nothing in the
current data, and not quoted as results.

**And one thing to stop doing.** `scripts/run_scenarios.py` is the only end-to-end chain
and it runs on `generate_plate(...)`. It should keep existing — it is the integration test
for the chain — but it should be renamed or clearly labelled so it cannot be read as a
result, and `MAINTENANCE_PER_ACTIVITY = 1200.0` ("a plausible scale, not a fitted one")
should stay exactly as annotated. It is not a prediction and the file should never imply
it is. The prediction is §0, and §0 is a loss, and a scored loss is worth more than an
unscored chain.

---

## 8. Annotated bibliography

Every entry below was checked against a page fetched during the writing of this document.
Anything that could not be confirmed is marked `UNVERIFIED` and says what was and was not
confirmed. Bibliographic traps are collected in §8.7.

### 8.1 Proper scoring rules

**Gneiting, T. & Raftery, A. E. (2007). "Strictly Proper Scoring Rules, Prediction, and
Estimation." *Journal of the American Statistical Association* 102(477), 359–378.**
DOI [10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437).
Author copy: <https://sites.stat.washington.edu/raftery/Research/PDF/Gneiting2007jasa.pdf>
— The organising reference for §2.1. Defines properness and strict properness, and gives
the CRPS in both the integral form and the kernel/expectation form that yields the exact
ensemble estimator used in §0. Also introduces the **energy score** as the multivariate
generalisation, which is the §2.2 recommendation for the joint 24-vector.

**Matheson, J. E. & Winkler, R. L. (1976). "Scoring Rules for Continuous Probability
Distributions." *Management Science* 22(10), 1087–1096.**
DOI [10.1287/mnsc.22.10.1087](https://doi.org/10.1287/mnsc.22.10.1087)
— The origin of the continuous ranked probability score. Cite this, not only Gneiting &
Raftery, when the point is CRPS's provenance rather than its properties.

**Scheuerer, M. & Hamill, T. M. (2015). "Variogram-Based Proper Scoring Rules for
Probabilistic Forecasts of Multivariate Quantities." *Monthly Weather Review* 143(4),
1321–1334.** <https://journals.ametsoc.org/view/journals/mwre/143/4/mwr-d-14-00269.1.xml>
— Relevant if the trajectory's *correlation structure* ever becomes the object of the
claim. The energy score is relatively insensitive to mis-specified correlations; the
variogram score is built to catch them. Not needed for the first scored prediction. DOI
prefix `10.1175/MWR-D-14-00269.1` is stated on the AMS landing page.

### 8.2 Calibration and sharpness

**Gneiting, T., Balabdaoui, F. & Raftery, A. E. (2007). "Probabilistic Forecasts,
Calibration and Sharpness." *Journal of the Royal Statistical Society: Series B
(Statistical Methodology)* 69(2), 243–268.**
DOI [10.1111/j.1467-9868.2007.00587.x](https://doi.org/10.1111/j.1467-9868.2007.00587.x)
— The source of the "maximise sharpness subject to calibration" paradigm and of the PIT
histogram as the calibration diagnostic. §2.3 is this paper applied to a particle filter.
The key move for this project is the insistence that a single proper score cannot
substitute for a calibration check — which is exactly why §0 reports coverage beside CRPS.

### 8.3 Prequential validation

**Dawid, A. P. (1984). "Present Position and Potential Developments: Some Personal Views.
Statistical Theory: The Prequential Approach." *Journal of the Royal Statistical Society:
Series A (General)* 147(2), 278–290.**
DOI [10.2307/2981683](https://doi.org/10.2307/2981683).
Copy: <https://www.cs.ubc.ca/~murphyk/MLRG/dawid84Prequential.pdf>
— The principle behind §2.4 and behind §5.7's pre-registration checklist: a forecasting
system is to be judged only on the forecasts it actually issued, each using only the data
available at issue time. `estimator.py::update`'s refusal to accept a backwards
observation is this principle implemented in a type check.

### 8.4 Nulls and permutation tests

**Ojala, M. & Garriga, G. C. (2010). "Permutation Tests for Studying Classifier
Performance." *Journal of Machine Learning Research* 11, 1833–1863.**
<https://www.jmlr.org/papers/v11/ojala10a.html>
— The reference for §4 Rung 3. Its useful contribution here is the distinction between
permuting *labels* (tests whether there is any label–feature relation) and permuting
*within features* (tests whether a specific feature carries the relation). For this project
the analogue is shuffling dose *within plate* rather than across, which is the §4 caution.

### 8.5 Digital twin definitions and formalism

**Kapteyn, M. G., Pretorius, J. V. R. & Willcox, K. E. (2021). "A probabilistic graphical
model foundation for enabling predictive digital twins at scale." *Nature Computational
Science* 1(5), 337–347.**
DOI [10.1038/s43588-021-00069-0](https://doi.org/10.1038/s43588-021-00069-0).
Preprint: [arXiv:2012.05841](https://arxiv.org/abs/2012.05841) (v2, 9 Feb 2021) — same
title, same three authors, and the formalism in §1.1 was read off the v2 PDF directly
(equations 1–4, Figures 2–3), so the notation quoted there is exact rather than
paraphrased.
— **The single most useful reference in this bibliography.** Supplies the six-quantity
abstraction (`S`, `D`, `O`, `U`, `Q`, `R`), the dynamic-decision-network formalism, the
belief-state factorisation into `φ^update`, `φ^QoI` and `φ^evaluation`, and the explicit
separation of an assimilation phase from a prediction phase extending to a horizon `t_p`.
Demonstrated on a structural digital twin of a UAV, calibrated on experimental data. The
statement that "the digital state space will generally be only a subset of the physical
state space" is the formal licence for `estimator.py`'s choice to make promoter activity a
latent state rather than a derived number. Affiliations at publication: Kapteyn at MIT
Aeronautics & Astronautics, Pretorius at The Jessara Group, Willcox at the Oden Institute,
UT Austin.

### 8.7 Traps, including three in this project's own data

**Data traps found while writing this document — these matter more than the bibliographic
ones, because they would have silently corrupted the scoring.**

1. **The second block of each channel is not a repeat read.** Each export exposes
   `OD600[1]`/`OD600[2]` and `mCitrine[1]`/`mCitrine[2]`. They look like duplicate reads and
   would give a beautiful instrument-noise estimate. Block 2 is block 1 minus a constant
   (OD −0.097 exactly; mCitrine −409 to −412, the variation being integer rounding).
   `KineticBlock.blank_subtracted` distinguishes them and `SynergyRun.raw_channel` already
   picks correctly — but a hand-rolled noise estimate from the pair would report the reader
   as roughly 10× more precise than it is, and would then make every NIS value in §6.5
   wrong in the flattering direction. **There are no independent repeat reads in this
   dataset.**
2. **The plate count differs by channel and by plate, and `README.md` states the OD figure.**
   `docs/DATA_INVENTORY.md` already documents this. Confirmed independently here: 4
   biological replicates with fluorescence for UPRE1, **2** for UPRE2 / NativeYap1 /
   AlteredYap1, because replicates 2 and 4 read fluorescence on columns 1–3 only. Any power
   calculation that uses "4 replicates" is wrong for three of four constructs.
3. **The top two H₂O₂ doses are outside the model's domain, and fail silently.** At 2 and 4
   mM the recovered `activity_late` is **negative** (−38, −53, −6.8, −172 in
   `outputs/sensor_characterisation.csv`) because the culture dies, OD declines and the
   inversion `dR/dt + (μ + k_deg)R` has no branch for a shrinking population. A fold change
   computed over these is a ratio of negative numbers and will look finite and plausible.
   Exclude them by a pre-registered rule, or state that they are included and will dominate
   every aggregate.

**A methodological trap, demonstrated rather than warned about.** §0 finds that the
aggregate RMS comparison and the paired per-well mean-absolute comparison give *opposite
verdicts* on the same forecasts (twin wins on FL at 2 h by RMS; ties or loses by paired
MAE), because the log-linear baseline has fat tails that RMS punishes and MAE does not.
Pre-register the rule.

**A reasoning trap, also demonstrated.** §6.3 measures that 87% of the fluorescence
technical variance is a persistent per-well offset, and §6.5's whiteness failure looks like
its obvious consequence. §6.5b prototypes the fix and it does nothing, because a constant
offset is already absorbed by free states initialised from the data. A correct measurement
plus a plausible mechanism is still not a diagnosis until the fix is tested.

**A sampling trap, also demonstrated, and the most embarrassing one.** An earlier draft of
§6.5 characterised the filter on `cult[:36]` — the first 36 culture wells of one plate. On
a plate laid out with dose along rows, the first 36 wells are the three lowest dose rows,
where the model is nearly right. That subset reported median NIS 0.42–0.72 (**over-dispersed**);
the full 167-well population reports 1.76–4.44 (**over-confident**). Same code, same
priors, opposite verdict, because the subset was the easy third of the plate. On a plate
with structure along both axes, *any* prefix of the well list is a stratified sample of
something. Enumerate the whole eligible set or sample it at random.

