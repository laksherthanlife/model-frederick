# G4 — the anchor problem: what the existing data can and cannot settle

Research note. Written against the data in `outputs/g4_rt_minus_qc.csv`,
`outputs/g4_anchor_fold_change.csv` and `outputs/g4_reporter_response.csv`, and against
the code in `src/ystwin/qpcr.py`, `src/ystwin/gates/g4_anchor.py` and `scripts/run_g4.py`.

Every number below was recomputed from those files; none is copied from the README.

---

## 1. Verdict

**The ER arm is unsalvageable and must be re-run. The oxidative arm can be pushed from
INCONCLUSIVE to a bounded, defensible statement today, but not to PASS.**

The Hac1 assay is 91% genomic DNA (median RT- margin +0.13 cycles, UPRE1). A gDNA
subtraction is *median-unbiased* there but its 95% interval on a true 2-fold change spans
eleven orders of magnitude, and 7 of 36 Hac1 readings have a **negative** margin, which
falsifies the additive model the correction rests on. Total HAC1 is also the wrong
transcript: HAC1 is regulated by Ire1-dependent splicing, not by transcription. Two
independent fatal faults, so re-run with **KAR2** and DNase-treated RNA.

The TRX2 arm is clean but was measured at the wrong doses against an unstable normaliser.
At n=3 the minimum detectable effect is **2.77-fold**; the anchor moved ~1.5-fold. Pooled
TOST does yield one real result (below). New wet-lab work is unavoidable for a PASS.

---

## 2. The gDNA correction — the algebra, and where it dies

### 2.1 The model and the estimator

Write `N = E^(-Cq)` for the relative starting quantity at a fixed threshold, with `E = 2`
for a perfect doubling per cycle. On one RNA sample split into a +RT and a −RT reaction:

```
N_plus  = N_cDNA + N_gDNA        (reverse transcriptase present)
N_minus =          N_gDNA        (reverse transcriptase omitted)
=>  N_cDNA = N_plus - N_minus
```

Define the **margin** `m = Cq_minus − Cq_plus`. Then the gDNA share of the +RT signal is

```
f = N_minus / N_plus = 2^(-m)
```

and the cDNA-only Cq is recovered without ever leaving Cq space:

```
Cq_cDNA = Cq_plus − log2(1 − 2^(−m))
```

This is a strictly *additive* correction in linear space and a *subtractive* one in Cq
space. It is well defined for any `m > 0` and undefined for `m ≤ 0`.

### 2.2 What the correction costs: the variance inflation factor

Let `s` be the SD of a single Cq measurement. A relative quantity has
`SD(N)/N = ln(2)·s`. Treating the +RT and −RT wells as independent,

```
Var(N_cDNA) = (ln2·s)² (N_plus² + N_minus²)

SD(N_cDNA)/N_cDNA   sqrt(1 + 2^(−2m))
------------------ = ----------------- ≡ VIF(m)
    ln2·s              1 − 2^(−m)
```

`VIF(m)` is the factor by which correcting inflates the noise, expressed as an effective
Cq SD `s_eff = s · VIF(m)`. In this dataset the **technical-replicate Cq SD is s = 0.11
cycles** (median over 72 duplicate pairs; UBC 0.10, TRX2 0.10, Hac1 0.11 — the Hac1 *mean*
is 0.46 only because one pair differs by 5 cycles).

| margin `m` | gDNA share `f` | `VIF(m)` | effective Cq SD | uncorrected bias on a true 2× |
| ---: | ---: | ---: | ---: | ---: |
| **0.13** (UPRE1 Hac1) | **91.4 %** | **15.7×** | 1.73 cyc | reads 1.09× |
| 0.27 | 82.9 % | 7.6× | 0.84 cyc | reads 1.17× |
| **0.60** (UPRE2 Hac1) | **66.0 %** | **3.5×** | 0.39 cyc | reads 1.34× |
| 1.00 | 50.0 % | 2.24× | 0.25 cyc | reads 1.50× |
| **1.72** | 30.4 % | **1.50×** | 0.17 cyc | reads 1.70× |
| 2.00 | 25.0 % | 1.37× | 0.15 cyc | reads 1.75× |
| **2.68** | 15.6 % | **1.20×** | 0.13 cyc | reads 1.84× |
| **3.00** (repo's bar) | 12.5 % | **1.15×** | 0.13 cyc | reads 1.88× |
| **3.52** | 8.7 % | **1.10×** | 0.12 cyc | reads 1.91× |
| 8.09 (TRX2) | 0.4 % | 1.00× | 0.11 cyc | reads 1.99× |

Solving `VIF(m) ≤ c` for the margin gives the thresholds in bold:

```
VIF ≤ 1.50  ⟺  m ≥ 1.72 cycles
VIF ≤ 1.20  ⟺  m ≥ 2.68 cycles
VIF ≤ 1.10  ⟺  m ≥ 3.52 cycles
```

**The repo's existing 3-cycle bar (`_MIN_RT_MINUS_MARGIN = 3.0` in `qpcr.py`) is not a
folk rule — it is almost exactly the margin at which the gDNA correction costs 15 % extra
noise and leaves ≤ 6 % bias on a 2-fold effect.** Keep it, and now it has a derivation.

### 2.3 The identifiability floor

Two separate conditions must hold before the correction means anything.

**(a) `m` must be significantly positive.** The margin is a difference of two Cq values.
With two +RT wells and one −RT well, `SE(m) = s·sqrt(1/2 + 1) = 0.135 cycles`. Requiring
`m > 2·SE(m)` gives an absolute floor of

```
m_min ≈ 0.27 cycles
```

UPRE1's median Hac1 margin is **+0.13 cycles — below the floor**. It is not statistically
distinguishable from "the −RT and +RT reactions contain the same thing".

**(b) The additive model must not be falsified.** `N_plus ≥ N_minus` is a structural
requirement: adding cDNA cannot reduce the template. **7 of 36 Hac1 readings have `m < 0`**
(minimum **−5.03** cycles, 13 Aug), and **28 of 36 have `m < 1`**. A negative margin means
the +RT and −RT wells are not a matched pair — RT-buffer inhibition, unequal template, or a
mis-assigned well. Whatever the cause, `N_plus − N_minus` is then negative and the
correction returns a complex or clipped number. There is nothing to correct.

### 2.4 Monte Carlo: what actually happens if you do it anyway

20 000 simulations per row. True cDNA fold change = 2.00, `s = 0.11` cycles, correction
applied at both doses.

| margin | gDNA share | recovered fold (median) | 95 % range |
| ---: | ---: | ---: | :--- |
| **0.13** | 91.4 % | 1.91 | **[0.00, 3.0 × 10¹¹]** |
| 0.60 | 66.0 % | 2.01 | [1.12, 4.34] |
| 1.00 | 50.0 % | 2.01 | [1.36, 3.13] |
| 2.00 | 25.0 % | 2.00 | [1.53, 2.63] |
| 3.00 | 12.5 % | 2.00 | [1.58, 2.55] |
| 8.60 | 0.3 % | 2.00 | [1.62, 2.48] |

At `m = 0.13` the estimator is **unbiased and useless**. That is the whole answer to
"is the correction defensible at a 0.13-cycle margin": it is *correct* and it *carries no
information*. At `m ≥ 3` the interval is indistinguishable from the contamination-free
case, which is the same threshold the VIF algebra gives.

### 2.5 Running it on the real data confirms the algebra

Applying `Cq_corr = Cq − log2(1 − 2^(−m))` to every reading (target **and** reference) and
recomputing `delta_delta_cq`:

| construct | anchor | raw fold at 0.5 mM (3 reps) | gDNA-corrected |
| --- | --- | --- | --- |
| AlteredYap1 | TRX2 | 2.09, 1.41, 1.18 → SNR 1.19 | 2.37, 1.41, 1.16 → SNR **1.01** |
| NativeYap1 | TRX2 | 1.63, 1.17, 0.73 → SNR 0.39 | 1.93, 1.15, 0.74 → SNR **0.45** |
| UPRE2 | Hac1 | 1.33, 1.45, 1.29 | 1.54, 3.09, 1.62 |
| UPRE1 | Hac1 | 0.90, 2.44, 1.01 | **5.1 × 10⁸, 1.8 × 10⁹, 2.9** |

The UPRE1 row is the correction detonating exactly as predicted. On TRX2 the correction is
a wash — it removes a real but tiny bias (gDNA share 0.4 %) and adds its own noise via the
contaminated **UBC** denominator.

### 2.6 The finding nobody has looked at: the reference gene is contaminated too

`scripts/run_g4.py` line 125 filters the QC frame to the *anchor* gene:

```python
margin = qc[(qc.construct == construct) & (qc.target == target)]   # target is Hac1 or TRX2
result = anchor_agreement(frame, anchor_qc_pass_rate=float(margin.passed.mean()))
```

**UBC's own RT- margins never reach a verdict.** They are printed (line 90 iterates all
targets) but never gate. Recomputed:

| target | median margin | pass rate @ 3 cyc | gDNA share |
| --- | ---: | ---: | ---: |
| TRX2 | 7.75 – 8.09 | 100 % | 0.4 % |
| **UBC** | **2.82 – 3.85** | **44 – 67 %** | **7 – 14 %** |
| Hac1 | 0.13 – 0.60 | 0 % | 66 – 91 % |

Per replicate, UBC's mean margin is **2.41 (24 Jul)**, 3.37 (11 Aug), 4.08 (13 Aug) — the
24 Jul run fails outright (max margin 3.55). **Overall 76 of 144 readings clear the
3-cycle bar; the reference gene is half of the failures.**

Because `ΔCq = Cq_target − Cq_reference`, a reference biased low by `b` cycles biases every
`ΔCq` high by `b`, and `ΔΔCq` carries `b(dose) − b(0)`. Observed within-replicate swings in
UBC margin reach 1.3 cycles, i.e. ±0.2 cycles of spurious ΔΔCq — a 15 % phantom fold
change. Small, but it is pure artefact and it is currently invisible to the gate.

**Implementation (small, do it now):** in `run_g4.py`, gate on the *minimum* of the anchor's
and the reference's pass rate, not the anchor's alone —

```python
anchor_qc = qc[(qc.construct == construct) & (qc.target == target)].passed.mean()
ref_qc    = qc[(qc.construct == construct) & (qc.target == REFERENCE_TARGET)].passed.mean()
result = anchor_agreement(frame, anchor_qc_pass_rate=min(anchor_qc, ref_qc))
```

and widen `anchor_agreement`'s failure message so it names *which* gene failed. This turns
a silent contamination of the denominator into a stated one.

### 2.7 A second, independent diagnostic the −RT wells already give you

The −RT reaction is a *gDNA assay*. `HAC1`, `UBC` and `TRX2` are all single-copy loci in a
haploid genome, so in the same −RT reaction they should amplify at comparable Cq for
comparable primer efficiency. They do not:

```
mean RT- Cq:   UBC 22.94    TRX2 24.17    Hac1 26.66
Hac1 − UBC = +3.72 cycles  ->  13x lower apparent template on identical DNA
```

The "Hac1 Exon" primer pair is **13-fold less efficient than the UBC pair on the same
template**. So the Hac1 assay is not merely measuring the wrong thing — it is a poor assay
on top of that. Any replacement primer pair should be validated by this same test before it
is trusted: run it on the −RT wells and check it lands within ~1 cycle of the reference
pair.

---

## 3. Statistical rescue routes, ranked

Baseline for every route below: per-observation SD of `ΔΔCq` on the oxidative arm is
**σ = 0.45 log2 units** (residual variance around the dose means; 0.38 AlteredYap1,
0.55 NativeYap1). Between-replicate SD of the top-dose log2 fold is 0.42 (AlteredYap1)
and 0.58 (NativeYap1).

### R1. Equivalence testing (TOST) on the pooled oxidative anchor — **the one route that works today**

**Idea.** The gate's third verdict conflates two different states: "we did not detect an
effect" and "we have positive evidence the effect is small". TOST separates them. Pool the
two Yap1 constructs (they share TRX2, UBC, the stressor and the dose ladder), giving n = 6
construct-replicate units, and test both one-sided nulls against a smallest effect size of
interest (SESOI).

**Result on the real data** (`y` = log2 fold at 0.5 mM H₂O₂, all six units):

```
mean +0.379  sd 0.516  95% CI [-0.163, +0.921]   (fold 1.30x, CI [0.89x, 1.89x])

SESOI ±0.5 log2 (1.4x):  TOST p = 0.295   not equivalent
SESOI ±1.0 log2 (2.0x):  TOST p = 0.016   EQUIVALENT
SESOI ±1.6 log2 (3.0x):  TOST p = 0.001   EQUIVALENT
```

**This is a real, publishable result extractable from the existing three replicates:**
*at 0.5 mM H₂O₂ the TRX2 anchor's response is bounded below 2-fold (equivalence test,
p = 0.016, n = 6).* Per construct at n = 3 only the 3-fold bound is reachable
(AlteredYap1 p = 0.028, NativeYap1 p = 0.026).

**Why it matters for the gate.** The reporters move 1.68× (AlteredYap1) and 1.89×
(NativeYap1) at that dose. An anchor bounded below 2-fold cannot certify or refute a
1.7-fold reporter — the two are the same size. The honest verdict is not "the anchor is
broken" but **"the anchor and the reporter agree to within the anchor's resolution, and
that resolution is too coarse to be evidence."**

**Implementation.** Add to `gates/g4_anchor.py`:

```python
def anchor_equivalence(effects, sesoi_log2: float = 1.0, alpha: float = 0.05):
    """Two one-sided tests: is the anchor's effect provably smaller than sesoi?"""
    n = len(effects); m = np.mean(effects); se = np.std(effects, ddof=1)/np.sqrt(n)
    p = max(stats.t.sf((m + sesoi_log2)/se, n-1), stats.t.cdf((m - sesoi_log2)/se, n-1))
    return p < alpha, p
```

and give `AnchorResult` a fourth verdict, `BOUNDED` — "the anchor moved by less than the
SESOI, so no reporter claim of that size can be tested against it". That is strictly more
informative than INCONCLUSIVE and it is what the data support.

**Failure mode.** TOST is only as defensible as the SESOI. Setting it *post hoc* to
whatever the data will clear is p-hacking in reverse. Fix the SESOI from the biology
*before* looking: for a stress anchor, 2-fold is the conventional line for "a real
transcriptional response", and the reporter's own dynamic range (1.7–2.6×) is the natural
alternative. State it in the code as a default and never tune it per construct.

### R2. Partial pooling across the two constructs sharing an anchor — **cheap, modest, real**

**Idea.** NativeYap1 and AlteredYap1 measure the *same gene* (TRX2), normalised to the
*same gene* (UBC), under the *same stressor*, in the *same qPCR runs*. Their anchor
responses are two draws from one distribution, so the between-replicate variance component
should be estimated from six units, not three.

**Result.** Method-of-moments decomposition of the top-dose effect:

```
per-construct means:  AlteredYap1 +0.601   NativeYap1 +0.157
within-construct var  0.259      between-construct var (MoM)  0.012
shrinkage weight  lambda = tau2/(tau2 + sigma2/n) = 0.122
  AlteredYap1  +0.601 -> +0.406      NativeYap1  +0.157 -> +0.352
SE of the pooled effect  0.208   (vs 0.294 per construct)   ->  t = +1.82
```

λ = 0.12 means the data support **near-complete pooling**: the two constructs' anchors are
statistically indistinguishable, which is exactly what should happen when they are the same
assay on different strains. Pooling cuts the SE by `sqrt(2)` — equivalent to doubling the
replicates for free.

**Expected gain.** Minimum detectable effect at 80 % power falls from **2.77-fold (n = 3,
per construct)** to **1.56-fold (n = 6, pooled)**. That is the single largest free gain
available.

**It still does not reach significance:** t = 1.82, p ≈ 0.13.

**Implementation.** `anchor_agreement` currently takes one construct's frame. Add an
optional `pool_key` column (default: the anchor gene) and estimate `anchor_spread` from all
rows sharing that key while reporting `anchor_effect` per construct. A full Bayesian version
(`brms`: `y ~ dose + (dose | replicate) + (1 | construct)`) buys little beyond this at
k = 3 — see R6 — but it does give an honest posterior instead of a point SNR.

**Failure mode.** Pooling is only licensed because the anchor gene, reference gene, stressor
and dose ladder are shared. `qpcr.py::assert_single_stressor` already refuses to pool across
DTT and H₂O₂ — extend the same discipline: refuse to pool across anchor genes.

### R3. Regress the noisy anchor **on** the precise reporter, not the reporter on the anchor

**Idea.** The current gate computes `_correlation(per_dose.reporter, per_dose.anchor_mean)`
— a symmetric Pearson r between two per-dose means. Both are measured with error, so r is
attenuated by `sqrt(rel_x · rel_y)` (Spearman's correction for attenuation). Recomputing
reliabilities from the data:

```
anchor reliability (ICC of the dose signal):  NativeYap1 0.00   AlteredYap1 0.33
  -> reliability of a 3-replicate dose mean:  NativeYap1 0.00   AlteredYap1 0.60
reporter reliability across 3 independent plates: ~0.98
   (CV at 0.5 mM: AlteredYap1 3.6%, NativeYap1 4.1%, UPRE1 3.4%, UPRE2 4.9%)
```

Maximum attainable `|r_obs| = sqrt(rel_anchor · 0.98)`:

| anchor reliability | max attainable r | gate needs 0.70 |
| ---: | ---: | :--- |
| 0.10 | 0.31 | impossible |
| 0.25 | 0.49 | impossible |
| **0.50** | **0.70** | exactly at the bar |
| 0.75 | 0.86 | reachable |

**So the gate's `min_correlation = 0.7` is unreachable unless the anchor's per-dose
reliability exceeds 0.50.** NativeYap1's is 0.00. This is not a threshold that needs
loosening — it is a specification on the *anchor*, and it should be stated as one.

**The fix.** OLS is unbiased when the *predictor* carries no error. The reporter is 25×
more precise than the anchor, so it is the predictor. Regress
`anchor_log2 ~ reporter_response`, replicate-centred (which makes it a within-replicate,
paired contrast), using all `doses × replicates` observations rather than 3 dose means:

```
NativeYap1                slope +0.184 +-1.015   t=+0.47  df=5   p=0.66
AlteredYap1               slope +0.916 +-1.178   t=+2.00  df=5   p=0.10
POOLED Yap1 (common slope) slope +0.455 +-0.665  t=+1.51  df=11  p=0.16
UPRE1                     slope +1.360 +-7.134   t=+0.49  df=5   p=0.64
UPRE2                     slope +0.890 +-2.089   t=+1.10  df=5   p=0.32
```

Nothing reaches significance, but the estimator is now unbiased, has a confidence interval,
and uses 9 points per construct instead of 3.

**The hard structural problem this exposes.** With `n_doses = 3`, Fisher's z has
`SE = 1/sqrt(n − 3)` — **division by zero.** *The gate's correlation criterion has no
confidence interval at 3 doses at all.* It is a point estimate with one degree of freedom.

| n_doses | 95 % lower bound if r = 0.90 | if r = 0.95 |
| ---: | ---: | ---: |
| 3 | **undefined** | **undefined** |
| 4 | −0.45 | −0.13 |
| 5 | +0.09 | +0.42 |
| 6 | +0.33 | +0.60 |
| 7 | +0.46 | +0.69 |
| **8** | +0.53 | **+0.74** ✓ |
| 10 | +0.62 | +0.80 ✓ |

**To certify `r ≥ 0.70` you need ≥ 8 shared doses and a near-perfect correlation.** The
experiment supplied 3. This is the cheapest thing to fix in the whole report — see §5.

The same arithmetic condemns the growth test. Recomputed `growth R²` over the 3 shared
doses versus all 7:

| construct | growth R², 3 shared doses | growth R², all 7 doses |
| --- | ---: | ---: |
| AlteredYap1 | 0.27 | 0.81 |
| NativeYap1 | 0.79 | 0.33 |
| UPRE1 | 0.70 | 0.02 |
| **UPRE2** | **0.98 (would REFUTE)** | 0.78 |

UPRE2 crosses the `growth_confound_r2 = 0.9` refutation line on 3 points and falls well
below it on 7. **A gate that can refute a construct on an R² with one degree of freedom is
not a gate.** Require `n_doses ≥ 5` before the correlation and growth tests run at all, and
return INCONCLUSIVE with that reason otherwise.

### R4. Better normalisation — the largest single variance sink, and it is not the anchor

**This is the most surprising result in the dataset.** Decomposing the between-replicate
variance of the top-dose ΔΔCq into its target and reference parts:

```
NativeYap1:  Var(dCq_TRX2) = 0.202   Var(dCq_UBC) = 0.859   Cov = 0.361
             -> Var(ddCq)  = 0.202 + 0.859 - 2(0.361) = 0.339
AlteredYap1: Var(dCq_TRX2) = 2.305   Var(dCq_UBC) = 3.190   Cov = 2.659
             -> Var(ddCq)  = 2.305 + 3.190 - 2(2.659) = 0.177
```

**The reference gene contributes 4× more variance than the anchor for NativeYap1
(0.859 vs 0.202), and for NativeYap1 normalising to UBC makes the measurement *noisier*
than not normalising at all** (Var(ΔΔCq) 0.339 > Var(ΔCq_TRX2) 0.202).

UBC is not stable here. Its Cq range across the three doses within one
construct-replicate unit:

```
median 0.81 cycles, max 2.23 cycles
exceeds 0.5 cycles in 11 of 12 construct-replicate units
```

A normaliser is supposed to be flat to well under half a cycle. UBC swings up to 2.2 cycles
(4.6-fold) — **larger than the anchor effect it is meant to normalise** (~0.6 log2).

**And the sign of the whole result depends on the normaliser.** Recomputing the pooled
oxidative anchor four ways:

| normalisation | pooled mean log2FC at 0.5 mM | SD | SNR |
| --- | ---: | ---: | ---: |
| TRX2 / UBC (current) | **+0.379** | 0.516 | 0.73 |
| TRX2, no normalisation | **−0.604** | 1.005 | 0.60 |
| TRX2 / UBC(−RT), gDNA as loading control | −0.149 | 0.755 | 0.20 |
| TRX2 / geomean(UBC, UBC(−RT)) | +0.115 | 0.384 | 0.30 |

**TRX2 goes up if you normalise to UBC and down if you do not.** That is not an
underpowered measurement — it is a measurement whose *direction* is set by an analysis
choice. No amount of statistics fixes it; only a validated multi-gene normaliser does.

**Expected gain from more reference genes.** With independent reference errors,
`s_R → s_R/sqrt(k)`:

| target SD | reference SD | k = 1 | k = 2 | k = 3 |
| ---: | ---: | ---: | ---: | ---: |
| 0.45 | 0.80 | 0.918 | 0.723 | **0.645** |
| 0.45 | 1.14 | 1.226 | 0.923 | **0.797** |

Going from one reference gene to three cuts SD(ΔΔCq) by ~1.4×, **worth twice the
replicates and costing two extra rows on a plate you are already running.**

**Implementation.** `qpcr.py::delta_delta_cq` takes `reference_target: str`. Change it to
accept `reference_targets: Sequence[str]` and normalise to the **geometric mean of their
Cq values** (equivalently, the arithmetic mean in Cq space), which is the geNorm/qBase
convention. The single-string form stays supported for the existing runs.

**Failure mode.** Averaging in a *second unstable* reference makes things worse, not better.
Validate the candidate set with geNorm/NormFinder on the actual samples before adopting it —
which needs ≥ 3 candidates measured, so it is a plate decision, not an analysis decision.

### R5. Dose-response slope instead of a two-group contrast — **less than you would hope**

**Idea.** `anchor_agreement` computes the effect as `per_replicate[top_dose] −
per_replicate[control_dose]` — a two-point contrast that throws away the middle dose.
Fitting a slope over all doses uses everything.

**Algebra.** `SE(slope) = σ / sqrt(R · Sxx)` with `Sxx = Σ(x − x̄)²`.

| dose ladder | Sxx | SE(slope), R = 3 | SE, R = 4 | SE, R = 6 |
| --- | ---: | ---: | ---: | ---: |
| 0, 0.2, 0.5 (**current**) | 0.127 | 0.730 | 0.632 | 0.516 |
| 0, 0.2, 0.5, 1.0 | 0.567 | 0.345 | 0.299 | 0.244 |
| 0, 0.1, 0.2, 0.5, 1.0 | 0.652 | 0.322 | 0.279 | 0.228 |
| 0, 0.1, 0.2, 0.5, 1.0, 2.0 | 2.893 | 0.153 | 0.132 | 0.108 |

**Result on the real data:** the slope SNR is essentially the same as the two-point SNR
(AlteredYap1 1.36 vs 1.42; NativeYap1 0.29 vs 0.27; pooled 0.75 vs 0.73). Adding the middle
dose buys nothing here **because the between-replicate variance is a replicate-level
offset, not per-observation noise** — every dose in a replicate is shifted together, so
averaging across doses does not average it away.

The slope route pays only when the dose *range* widens. The gain from `Sxx` is real
(0.127 → 2.893, a 4.8× SE reduction) but it is a gain from **going to higher doses**, not
from adding doses in the middle.

**Failure mode.** A linear slope over a saturating dose-response is mis-specified, and both
reporters here are non-monotone above 1 mM H₂O₂ (response goes *negative* at 2 and 4 mM —
the cells are dying). Fit on log-dose over the monotone range, or use the measured reporter
response as the regressor (R3), which is shape-agnostic by construction.

### R6. Bayesian hierarchical modelling — honest, not powerful

**Idea.** `y[rep, dose, construct] ~ Normal(a + b·dose + u[rep], sigma)` with
`u[rep] ~ Normal(0, tau)` and a partially-pooled `tau` across constructs, fitted in Stan/brms.

**Expected gain: near zero in power, large in honesty.** The reason is `k = 3`. Between-
replicate variance estimation at k = 3 is notoriously unstable — the method-of-moments
estimate is frequently zero (as it nearly is here: `tau² = 0.012` against `sigma² = 0.259`),
and a REML/ML fit will often return `tau = 0` exactly, which then reports a fake-precise
interval. The standard fixes are (i) a weakly informative prior that keeps the variance
away from zero, and (ii) a small-sample interval adjustment.

The `t` penalty at small k is the blunt version of the same point:

```
k = 2:  t(0.975, 1)  = 12.706   (6.5x wider than the normal 1.96)
k = 3:  t(0.975, 2)  =  4.303   (2.2x wider)
k = 4:  t(0.975, 3)  =  3.182   (1.6x wider)
k = 6:  t(0.975, 5)  =  2.571   (1.3x wider)
```

**Going from 3 to 4 replicates shrinks every interval by 26 % before a single extra data
point is considered** — the fourth replicate is worth more than the third or the fifth.

**Failure mode, and it is fatal to one common idea:** biological replicate is **perfectly
confounded with qPCR run** in this design (one culture set per plate: 24 Jul, 11 Aug,
13 Aug). No hierarchical model can separate a run effect from a biological effect when they
are the same three levels. If the between-replicate spread is really plate-to-plate
technical variation, the only fix is a design fix — run two biological replicates on one
plate, or include an inter-run calibrator sample on every plate.

### R7. What does **not** work — recorded so nobody retries it

- **Bootstrap CIs at n = 3.** Percentile bootstrap on the AlteredYap1 effect gives
  fold 1.52 [1.18, 2.09] against the t-interval's 1.52 [0.73, 3.14] — **2.5× too narrow**.
  With only 3³ = 27 distinct resamples the bootstrap cannot reach outside the observed
  range and badly undercovers. Use the t-interval or a Bayesian posterior; do not bootstrap
  three numbers.
- **gDNA correction as a rescue for the oxidative arm.** Applied end-to-end it changes
  AlteredYap1's SNR from 1.19 to 1.01 and NativeYap1's from 0.39 to 0.45. Net zero. TRX2's
  gDNA share is 0.4 %; there is nothing there to recover.
- **Using the −RT signal as a loading control.** Tested (row 3 of the R4 table): SNR drops
  to 0.20. gDNA carryover tracks DNase efficiency, not input RNA mass, so it is not a
  loading measure.
- **More replicates on the current design.** The gate already computes this
  (`replicates_needed`): 80 for NativeYap1, 9 for AlteredYap1. Nine is conceivable; eighty
  is not. And both numbers are estimated from an SNR that is itself based on 3 points, so
  they carry enormous uncertainty. Changing the *anchor gene* is worth far more than
  changing `n` — see §4.

---

## 5. The confirmatory experiment

### 5.1 Why the current design could never have worked

Two numbers, side by side.

**Where the RNA was taken.** `plate/layout.py::RECORDED_PLATES["20260722"]` records:
*"RNA extracted from rows A, C and D (0, 0.2, 0.5 mM) after the read."* Rows A–G carry
7 doses. **RNA was taken from the bottom three.**

**Where the reporter actually responds** (mean of 3 plates, fold vs 0 mM):

| dose (mM) | 0 | 0.1 | 0.2 | 0.5 | **1.0** | **2.0** | 4.0 / 5.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AlteredYap1 (H₂O₂) | 1.00 | 1.25 | 1.38 | 1.68 | **1.69** | −0.10 | −0.14 |
| NativeYap1 (H₂O₂) | 1.00 | 1.33 | 1.48 | 1.89 | **1.75** | −0.26 | −0.23 |
| UPRE1 (DTT) | 1.00 | 1.07 | 1.14 | 1.41 | 1.86 | **2.11** | 1.86 |
| UPRE2 (DTT) | 1.00 | 1.09 | 1.15 | 1.40 | 1.89 | **2.18** | 2.00 |

The three sampled doses cover the bottom **40 %** of the reporter's H₂O₂ range and the
bottom **19 %** of its DTT range. The anchor was asked to confirm a reporter signal at the
doses where the reporter itself has barely moved.

**Minimum detectable effect** at σ = 0.45 log2, 80 % power, α = 0.05 two-sided:

| replicates | per construct | pooled across the 2 constructs sharing the anchor |
| ---: | ---: | ---: |
| 3 | **2.77-fold** | 1.56-fold |
| 4 | 1.94-fold | 1.43-fold |
| 5 | 1.69-fold | 1.36-fold |
| 6 | 1.56-fold | 1.32-fold |
| 8 | 1.43-fold | 1.26-fold |
| 12 | 1.32-fold | 1.20-fold |

The anchor moved ~1.5-fold. **At n = 3 the design's floor was 2.77-fold.** It was
underpowered by construction, and no third replicate was ever going to fix that.

Power as a function of the *true* anchor effect:

| true fold at top dose | R = 3 | R = 4 | R = 5 | R = 6 | R = 8 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1.5× | 0.26 | 0.43 | 0.59 | 0.72 | 0.88 |
| 2× | 0.54 | 0.83 | 0.95 | 0.99 | 1.00 |
| 3× | 0.85 | 0.99 | 1.00 | 1.00 | 1.00 |
| 5× | 0.98 | 1.00 | 1.00 | 1.00 | 1.00 |
| 10× | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

**Read this table sideways.** Detecting a 1.5-fold anchor at 80 % power needs R ≈ 8.
Detecting a 5-fold anchor needs R = 3. **Choosing an anchor gene with a large induction is
worth about three times as many replicates as choosing more replicates.** That is the whole
argument of §4 in one number.

### 5.2 The design

Exploit the one genuine strength of the existing protocol, which the README already names:
RNA was taken from the *same wells* as the fluorescence read. A paired design has
`Var(paired difference) = 2σ²(1 − ρ)`:

| ρ (reporter–anchor correlation within a culture) | SD ratio vs unpaired | effective-n multiplier |
| ---: | ---: | ---: |
| 0.3 | 0.84 | 1.4× |
| 0.5 | 0.71 | 2.0× |
| 0.7 | 0.55 | 3.3× |
| 0.9 | 0.32 | 10× |

Keep that. Change four things: the anchor gene, the doses, the reference set, and the DNase
step.

