# What survives a null

This records what happens when nulls are applied to the claims here. The machinery is
`analysis/nulls.py`.

An early version wrongly said the project had no nulls anywhere:
`analysis/recovery.py::module_recovery` already accepted `shuffle=True`, and the transfer
runner reported a shuffled module-recovery control. Transfer comparisons were added later,
as recorded below; the power/null gaps remain separately listed at the end. A null-test
non-rejection must not be turned into equivalence or a claim of no biological effect.

## The growth-dilution coupling is real, and about 12 points of it was never real

**Claim under test.** "Among cultures that pass G1 and admit the statistic, the median
share of `log(RFU/OD)` variance explained by growth rate alone is 0.82."

**Why it needed a null.** The regression is `log(specific)` on `−log(mu + k_deg)`. Both
sides are smooth, autocorrelated time series, and both are derived from the *same*
optical density trace — `specific = (RFU − background)/(OD − blank)` and
`mu = d(ln OD)/dt`. Two smooth autocorrelated series correlate spuriously as a matter of
course, because the effective number of independent observations in a 25-point trace is
far below 25. Here they are not even independent by construction.

**The null.** Phase randomisation of the growth-rate regressor: same power spectrum,
same autocorrelation, phases destroyed. This is the standard surrogate for testing
against a linear Gaussian process with the observed spectrum, and it is the *right* null
here precisely because a plain shuffle is not — a shuffle produces white noise with the
observed histogram, which any smoothness-sensitive statistic rejects for the wrong
reason. `tests/test_nulls.py` demonstrates that failure mode directly: on independent
smooth traces the shuffle null over-rejects and phase randomisation stays near nominal.

**Result**, 146 wells across the two NewProtocol plates that admit the statistic:

| | median R² | IQR |
| --- | ---: | --- |
| observed | **0.522** | [0.293, 0.869] |
| phase-randomised null | **0.117** | [0.073, 0.175] |

Wells beating their own null at p < 0.05: **83 / 146 (57%)**. Median p-value **0.022**.

**Reading.** The coupling is real. Autocorrelation alone buys R² ≈ 0.12, and the observed
median is more than four times that, so the relationship survives the strongest
like-for-like null available. That is a genuine strengthening of the finding: it was
previously a number, and it is now a number that has been attacked.

Two qualifications, both real.

**About 0.12 of the reported variance share is artefact.** Any well whose R² sits near
0.12 is indistinguishable from smoothness. The honest statistic is the *excess* over the
null, not the raw share.

**It is a majority effect, not a universal one.** Forty-three percent of wells do not
individually clear their own null, and the median p of 0.022 only just clears. So
"growth explains most of the reporter signal" is supported *on aggregate* and is not
established well-by-well — which matters, because the correction is applied well-by-well.

**Note on comparability.** The 0.522 median here is not the README's 0.82. This sweep did
not apply the G1 optical-quality gate, so it runs over a superset of wells including
poorer ones. The apples-to-apples comparison is observed against null on the *same* well
set — 0.522 against 0.117 — not either figure against the other. Re-running this sweep
behind G1 is a small job and would give the gated equivalent.

## Transfer: no detected superiority over the random-subspace comparison

**Historical numerical comparison, not the current nested-width pipeline.**
`scripts/run_transfer.py` gained its own null and the then-headline configuration did not
reject it. Its current scorer uses outer-training channel means and nested width selection.
The old table is kept to record the unsuccessful result, not silently presented as a
fresh run under the revised protocol:

| | median transfer R² | p |
| --- | ---: | ---: |
| **observed** | **+0.0499** | -- |
| null: random axes, same rank (`rotated_subspace`) | **+0.0362** | 0.268 |
| null: no cross-channel structure (`matched_marginals`) | −0.0012 | 0.024 |

The fitted model's score was above the rotated-null median, with p = 0.268: no detected
superiority under that test. The former "skill over it +0.051" is also withdrawn as a
summary of this table; it does not follow from the displayed observed/null values under
`(observed - null) / (1 - null)`. No replacement scalar is needed to interpret the
non-rejection.

The earlier recorded comparison was observed +0.0298 against rotated null +0.0351 at
p = 0.634. H2O2 EC50 corrections changed the generator. The old attribution of later
movement to metabolite-pool arms was unsupported for these readers: `CROSS_FAMILY.md`
records that those arms are not reader-visible. Both historical directions are retained;
neither establishes equality or absence of transfer information.

The historical marginals-null rejection (p = 0.024) addresses a different, weaker
comparison about cross-channel structure. It does not validate the fitted latent axes.
**The old "arbitrary axes carry as much transfer information" inference is withdrawn.**
Failure to detect superiority is not an equivalence test, and no equivalence margin or
power calculation was supplied. The unsuccessful rotated-null comparison remains a
limitation, not evidence that latent structure has no role.

Two things stop this being the end of the transfer story, and neither rescues the
headline number.

The historical **oracle** range, 0.016 to 0.208, came from a model allowed to see the
held-out stressor. It is an in-sample diagnostic, not held-out transfer or a measured
ceiling on every estimator. Its low values motivate examining the channel design without
establishing an impossibility result.

This was one configuration. The historical eleven-fold combination lift, 0.038 to 0.429,
came from `select_combinations` choosing pairs on transfer over a different reporter set
and dose ladder. Separate probes over eight stressors gave +0.0005 rising to +0.0164 with
four arbitrary co-dosed pairs. Those poor probes remain, but neither arbitrary pairs nor
simply choosing a `--scale` budget evaluates the full pair-selection pipeline independently.
The older lift is not confirmed by the fixed-configuration null comparison.

`CROSS_FAMILY.md` keeps three scopes separate: fixed **three-state** transfer/nulls,
named-family dimension/design selection (three states in the corrected default candidate),
and a separate **two-state SPOTA candidate** sampled independently of reference sets. Their
JSON configurations are authoritative. The named-family score difference is descriptive;
the approximate UCBOG targets only the separate candidate's expected configuration gap
under the uniform synthetic sampler. **Selected-pipeline optimism remains pending.** None
of the historical non-rejections, the candidate bound or a freshly inspected external
artifact establishes equivalence or independent biological validation.

## Seed replication: one build comparison does not survive it

`run_training.py` now trains at three dataset seeds rather than one. The seed draws a
different simulated experiment, so anything that moves across seeds was a property of the
draw.

| build | latent width | modules recovered, per seed | mean score on driven modules |
| --- | --- | --- | --- |
| five-channel | 5, 5, 5 | 5, 5, 5 -> **5 [5, 5]** | 0.255, 0.244, 0.249 -> **0.249 [0.244, 0.255]** |
| three-sensor | 4, 4, 4 | 5, 4, 6 -> **5 [4, 6]** | 0.120, 0.094, 0.130 -> **0.120 [0.094, 0.130]** |

**These are full-run numbers.** An earlier version of this table reported
`run_training.py --quick` output -- four doses and three replicates instead of six and
six -- without saying so, and drew two conclusions from it that the full run does not
support. Corrected below, and the script now stamps its mode into every table it writes
so the two cannot be confused again.

**Module count still does not separate the builds**, but not for the reason first given.
Both medians are 5. What differs is stability: five-channel returns 5 at every seed, while
three-sensor spans 4 to 6. A single-seed comparison could have drawn either build ahead of
the other purely by which seed it ran.

**The score separates them by about a factor of two**, not by the third the `--quick` run
suggested. Ranges 0.244-0.255 against 0.094-0.130, comfortably non-overlapping. The
`--quick` figures understated the gap by flattering the three-sensor build, which is the
one with fewer channels and therefore the one a smaller dose ladder helps most.

Latent width was stable at both builds in this run. A review pass reported instability
(5, 3, 5) at full size which I could not reproduce over three seeds; recorded as
unresolved rather than asserted either way.

## The first held-out prediction on real plates: wins inside the range, loses outside it

`scripts/run_heldout_score.py` reads `outputs/split_manifest.csv`, fits a dose-response on
the training groups, and predicts the held-out ones. Two functional forms are scored: the
repo's own **biphasic** model from `generator/panel_calibration.py` -- basal plus
saturating induction over a viability term, four parameters -- and a **log-linear** fit as
a second baseline.

Interpolation's partition space has only five interior rungs, so it is **exhausted** rather
than sampled -- one dose held out at a time, which makes the partitions comparable:

| held out | test wells | biphasic | nearest-dose | skill |
| ---: | ---: | ---: | ---: | ---: |
| 0.1 mM | 12 | 309.9 | 370.7 | **+0.164** |
| 0.2 mM | 12 | 222.0 | 388.5 | **+0.429** |
| 0.5 mM | 12 | 296.1 | 387.6 | **+0.236** |
| 1.0 mM | 12 | 532.6 | 294.8 | **-0.807** |
| 2.0 mM | 12 | 896.0 | 760.9 | -0.178 |

Wins on three of five, and the sign separates on dose relative to the **fitted EC50** rather
than on dose itself -- crossover between 3.7 and 7.4 multiples of it. Every scored prediction
comes from UPRE1 or UPRE2, the only fits `DoseFit.identifiable` vouches for, whose EC50s are
0.121 and 0.148 mM. The mechanism, the two hypotheses it beat, and what it does not
establish are in [`WHY_INTERPOLATION_WORKS.md`](WHY_INTERPOLATION_WORKS.md).

Claimable: **on an interior dose at or below 0.5 mM the biphasic fit beats carrying the
nearest measured dose forward on every partition tested, by 0.16 to 0.43.**

**Extrapolation has nothing to score.** `DoseFit.identifiable` is false for all four
constructs there, with the reason attached: the fitted lethal dose sits below the fitted
EC50, so the fit is saying the culture dies before it induces. That is not a dose-response
and cannot be extrapolated from. The -1.913 once reported was the score of predictions made
from fits the module had already refused -- a refusal, not a failure, and the same
distinction the gates draw elsewhere.

Three earlier versions of the interpolation number were reported and replaced. They are in
[`superseded/heldout-interpolation.md`](superseded/heldout-interpolation.md) rather than
here, with what each claimed and why it moved.

What survives from all of them: **"use the nearest dose you measured" is a strong baseline**
for a smooth response on a fine ladder -- strong enough to beat a wrong functional form
everywhere and a right one outside its range. It belongs beside every future prediction here.

The two unscoreable splits are informative rather than missing:

`heldout_construct` returns no scoreable test row, because it withholds whole constructs
and the model is fitted per construct. That is the honest exposure of a limitation --
this model cannot transfer across constructs at all, and a cross-construct claim needs a
model that shares structure between them.

`heldout_replicate` has zero test rows, because the manifest's only feasible replicate
split is UPRE1 on plates 20260728 and 20260804 -- the two with no reporter-channel blank,
which is why they are absent from the characterisation table. The same missing blank that
blocks the fold change blocks the replicate-level prediction, from the other direction.

## Still without a null

Recorded so the gaps are visible rather than implied:
- **"These channels identify these modules."** Needs `shuffled_loadings`: does a random
  loading matrix of the same shape and scale produce comparable Fisher information?
- **Detection and attribution power tables.** These are simulation-derived and now carry
  Wilson intervals, but no null: a random design's power at the same replicate count is
  not reported.
- **Latent dimension selection.** See `docs/research/IDENTIFIABILITY.md` — at four to
  five channels the cross-validated choice occasionally lands on the full-rank
  degeneracy, and the Ledermann ceiling is not enforced at all.
