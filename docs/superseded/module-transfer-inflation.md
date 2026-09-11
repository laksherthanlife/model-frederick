# Superseded: how many modules transfer to an unseen stressor

**Current result:** the model recovers a held-out stressor's modules only where *other*
stressors in the panel drive the same module. Where the held-out stressor is a module's
only driver, recovery is zero — 0 of 6 such pairs, mean score 0.016.

Everything below was measured with a scoring bug in `module_transfer`, and every count it
produced was too high.

---

## The bug

`module_transfer` scored each module as a skill score against predicting the training mean:

```python
total = np.sum((truth - baseline) ** 2)          # baseline = training mean
left  = np.sum((truth - predicted) ** 2)
score = 0.0 if total <= 1e-12 else 1 - left / total
```

The guard catches the case where the truth *equals* the baseline. It does not catch the
case that actually dominates: the held-out stressor leaves a module **untouched**, so its
truth is a constant zero while the training mean sits well above zero because other
stressors drive that module hard.

Then any prediction closer to zero than the training mean scores well, without tracking
anything. Held out caffeine — which drives `cell_wall`, `nitrogen` and `ESR`, and nothing
else — the metric reported:

| module | truth (all held-out wells) | training mean | prediction | score |
| --- | ---: | ---: | ---: | ---: |
| peroxide | 0.0000 | 0.0770 | 0.0368 | **0.641** |
| redox | 0.0000 | 0.1299 | 0.0674 | **0.623** |
| oxidative | 0.0000 | 0.1264 | 0.0685 | **0.604** |
| proteasome | 0.0000 | 0.0825 | 0.0590 | **0.440** |
| copper | 0.0000 | 0.0187 | 0.0156 | **0.286** |

Every one of those has a truth span of exactly zero. Shrinking toward zero is free, and
the metric paid for it — in proportion to how strongly the *rest of the panel* drove the
module. The score was partly a measure of the panel, not of the model.

## What the numbers were, and are

Modules scoring above the 0.25 recovery threshold, three-sensor build, 8 replicates pooled:

| held out | own targets | reported | actual | inflation |
| --- | ---: | ---: | ---: | --- |
| tunicamycin | 3 | 8 | **1** | copper, iron, nadh, oxidative, peroxide, ph, redox |
| menadione | 5 | 9 | **6** | UPR, alkaline_ph, dna_damage |
| antimycin_A | 4 | 9 | **0** | all nine |
| DTT | 4 | 7 | **2** | copper, iron, nadh, oxidative, peroxide |
| H2O2 | 4 | 15 | **7** | UPR, alkaline_ph, calcium, cell_wall, dna_damage, ph, sulfur, zinc |
| glucose_starvation | 3 | 9 | **0** | all nine |
| sorbitol | 2 | 7 | **1** | copper, nadh, oxidative, peroxide, proteasome, redox |
| copper_sulfate | 3 | 0 | **0** | — |
| caffeine | 3 | 7 | **1** | copper, nadh, oxidative, peroxide, proteasome, redox |

Two stressors go to zero. `antimycin_A` and `glucose_starvation` both drive ATP, the build
has a dedicated ATP channel (QUEEN), and neither recovers it.

## The claims this retracts

- **"A stressor on one of its axes transfers well — 8 or more modules."** Withdrawn.
  Measured 1, 6 and 0 for the three stressors it was asserted on.
- **"The canonical three transfer — 5 or more modules."** Withdrawn. Measured 2, 7, 0.
- The two off-axis honesty tests (`sorbitol`, `caffeine`) had been **failing
  permissively** — recovering 6–7 modules against a bound of 3. They pass now at 1 and 1.
  The bound was right and the metric was wrong, which is the reverse of how it looked.

## What survived, and what replaced it

`StressModel.recovery` is **not** affected. It centres each module's truth on its own mean,
so a constant module has zero total and is guarded correctly. The in-sample headline —
9 of 24 modules for the three-sensor build, 5 without pooling — stands.

The replacement claim is graded and mechanistic. Across all 76 (held-out stressor, own
target) pairs in the panel, scored against how many other stressors drive the same module:

| peers driving the module | pairs | mean score | recovered |
| ---: | ---: | ---: | ---: |
| 0 | 6 | 0.016 | **0%** |
| 1–2 | 28 | 0.081 | 7% |
| 3–5 | 17 | 0.288 | 47% |
| 6+ | 25 | 0.291 | 84% |

Spearman +0.591. The zero-peer row is the falsifiable edge and is now asserted directly:
`TPEN→zinc`, `fluconazole→xenobiotic`, `methionine_starvation→sulfur`,
`sodium_hydroxide→alkaline_ph`, `copper_sulfate→copper`, `antimycin_A→retrograde` — six
pairs, none recovered.

So the honest description of what transfers is **interpolation across redundant coverage**,
not generalisation to unseen biology. A module has to appear twice in the panel, not once.

## The lesson

An R² or skill score against a constant truth is not a measure of fit — there is no
variance to explain, and whatever the baseline's offset happens to be becomes free credit.
The guard has to be on the **variance of the truth**, not on the distance between truth
and baseline. `_RECOVERABLE_SPAN` in `analysis/stress_model.py` is that guard.

It is worth noting how this was found: not by the metric looking wrong, but by a
falsification test failing in the permissive direction and being taken seriously instead of
having its bound relaxed.
