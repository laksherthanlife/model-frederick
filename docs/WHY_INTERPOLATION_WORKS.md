# Why the fit beats the baseline below 0.5 mM and not above

The measured result, from `scripts/run_heldout_score.py` over every interior rung:

| held out | dose / EC50 | test rows | skill |
| ---: | ---: | ---: | ---: |
| 0.1 mM | 0.7 | 18 | **+0.086** |
| 0.2 mM | 1.5 | 36 | **+0.368** |
| 0.5 mM | 3.7 | 18 | **+0.101** |
| 1.0 mM | 7.4 | 36 | **−0.683** |
| 2.0 mM | 14.8 | 18 | −0.178 |

**These magnitudes are smaller than the first published set, and that is the honest
direction.** The table was first computed on two biological replicates. A third was
recovered when a plate that had been dropped for a parsing fault was read correctly, and
re-running on it took the test set from 12 rows per dose to 18-36 and cut the skill at
0.1 and 0.5 mM by roughly half. The earlier figures were +0.164, +0.429 and +0.236.

What survived is the structure, which is the part the mechanism predicts: skill is
positive at every interior rung at or below 0.5 mM and negative above it, the sign still
separates on `dose / EC50`, and the crossover still sits between 3.7 and 7.4. What did
not survive is any description of this as a decisive win. It is a modest one.

The sign separates on dose relative to the fitted EC50, not on dose itself. Every scored
prediction comes from UPRE1 or UPRE2 — the only two constructs whose fits
`DoseFit.identifiable` vouches for — and their fitted EC50s are 0.121 and 0.148 mM. So the
crossover sits between **3.7 and 7.4 multiples of the model's own EC50.**

## The mechanism

The induction term is Michaelis–Menten in form, `d / (K + d)`. Its sensitivity to dose is

```
d/dd [ d / (K + d) ]  =  K / (K + d)²
```

which peaks at `d = 0` and falls as `1/d²` once `d ≫ K`. At eight times K the model's
prediction moves 6% as much per unit dose as it does at K; at fifteen times, 2%.

That is the whole of it. Deep in saturation the model **cannot move its prediction with
dose**, while the measured activity still varies — so the fitted curve contributes almost
nothing there and the residual is whatever the plateau height happens to be. Nearer K the
same parameter is strongly constrained by every training point, and pooling beats the local
answer.

This is the standard result behind optimal design for saturating models: the Fisher
information for `K` is concentrated at doses of order `K`, which is why D-optimal designs
for a Michaelis–Menten curve place their support there rather than spreading it up the
plateau (Matthews & Allcock, *Optimal designs for Michaelis–Menten kinetic studies*,
**PMID 14748040**; the general nonlinear-design framework is Box & Lucas, *Design of
Experiments in Non-Linear Situations*, **Biometrika 1959, doi:10.2307/2332810**).

The corollary that matters here: **a dose ladder that reaches far past its own EC50 buys
almost no information for the doses at the top.** The ladder runs to 5 mM against an EC50
of 0.13 — a factor of 37 — so most of its rungs sit where the model is least able to say
anything.

## Why the baseline is hard to beat at all

"Carry the nearest measured dose forward" is persistence, the standard hard baseline in
forecast verification for exactly this reason: on a smooth response measured on a fine
ladder its error is bounded by the local slope times the gap, which is small wherever the
response is flat. Since a saturating response *is* flat at high dose, **persistence gets
better in the same region where the model gets worse.** The two effects compound, which is
why the crossover is sharp rather than gradual.

## What this does and does not establish

**Established.** The sign of the skill separates perfectly on `dose / EC50` across all
five partitions, with the boundary between 3.7 and 7.4. The mechanism is a property of the
functional form, not of this dataset, and it predicts the direction before the data is
seen.

**Not established.** The *magnitude* does not track sensitivity — 0.2 mM scores better than
0.1 mM, where sensitivity is higher. Correlation of skill with log-sensitivity is +0.64
across five points, which is suggestive and no more. Five partitions cannot distinguish
this mechanism from others that predict the same ordering.

**Two mechanisms I tested and rejected**, recorded because a story that survives one test
is worth less than one that survived three:

*Leverage.* Prediction variance for a fitted model grows as `(x − x̄)² / Sxx`, so the
obvious explanation is that held-out doses far from the design centre are predicted worse.
It **fails**: the leverage term is *lowest* at 1.0 mM (0.006), where skill is worst, and
higher at 0.1 and 0.2 mM, where the model wins. Four of five partitions come out the wrong
way round.

*Lethal-parameter instability.* Withholding a high dose removes the data constraining the
viability denominator, so `lethal_dose` should become unidentifiable. **Partly supported
and not sufficient**: withholding 1.0 mM does give the largest spread in fitted lethal dose
(4.04 against 1.96–3.26 elsewhere) and the worst skill, and withholding 2.0 mM distorts the
fitted EC50 most (1.564 against 0.55–0.77). But the ordering across the other three
partitions does not follow, so it is at best a contributing term.

## The testable prediction

If saturation is the mechanism, the crossover should **move with the EC50** rather than
sitting at a fixed dose. A construct with a higher EC50 should keep winning to higher
doses. The two identifiable fits here differ by only 20% in EC50, which is not enough
spread to test it — so this is a prediction the current data cannot check, and the cheapest
way to check it is a construct with a genuinely different EC50 rather than more replicates
of these.

## Reproduce

```bash
python3 scripts/run_heldout_score.py     # the table, and outputs/heldout_interpolation_by_dose.csv
```

Superseded versions of this result, and why each moved, are in
[`superseded/heldout-interpolation.md`](superseded/heldout-interpolation.md).
