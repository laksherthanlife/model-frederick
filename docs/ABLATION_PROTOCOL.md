# ABLATION_PROTOCOL.md — how criterion (a) is scored, and why the rule had to be written down

Implementation: `src/ystwin/mech/ablation.py`. Tests: `tests/test_ablation.py` (39, green).
Every number in this document was produced by that code and can be reproduced with
`python3 -m pytest tests/test_ablation.py -q`; nothing here is quoted from a design document.

---

## 0. The rule, in one paragraph

**An ablation compares two FITTED models at stated parameter counts, never a mechanistic
trace against a frozen incumbent.** Both models are refit, in the same call, on the same
rows. Both parameter counts are reported. The effect is scored against a floor that carries
the name of the assay that measured it, and the floor's assay must be the observable's own.
Reporter activity is scored against **0.146** and growth against **0.0117 /h**; for any
other observable the caller names the assay and its floor, and if they cannot, the call
refuses.

## 1. Why this exists

Four of the five block rescues in the 2026-09-05 adjudication failed in the same way, and
`REVISED_BUILD_LIST.md` §5 records it: each scored a mechanistic trace against an incumbent
that had been **frozen, stripped of a term, or denied its own free parameter**. The sharpest
case is measured rather than argued — refitting the *real* incumbent
(`generator/culture.py::simulate_culture`, logistic, free K) to the TORC1 advocate's own
lagged trace gives a relative RMSE of 0.0021–0.0068, i.e. **0.015x–0.047x the 0.146 floor**,
where the advocate had reported a pass after deleting the incumbent's capacity term and
comparing against `X0·exp(mu·t)`.

The error is one-directional. Every instance of it **inflates the mechanistic block's
apparent contribution**. Ten blocks are scheduled to be scored against one incumbent, so an
architecture that does not make the rule executable will repeat it against itself.

The second error the protocol prevents is the wrong floor. Route C2 scored an **FBA product
ceiling** against the **plate-reader activity CV** — the wrong floor for the wrong
observable, in a repository whose own memory records the product ceiling being refuted by
63x. That is not a slip anyone would repeat knowingly; it is what happens when a floor is a
float.

## 2. What is enforced, and how

| Rule | Enforcement | Fails as |
|---|---|---|
| Both models are refit here | `ablate()` calls `model.refit(data)` on each side itself; the interface method is named `refit` because a caller may not hand in a prepared fit | — |
| Neither model is frozen | Each side is refit a second time on a **deterministically perturbed** copy of the same rows, and its parameters must move | `FrozenIncumbent` |
| Both counts are reported | `AblationResult.report()` prints both `n_free`, always | — |
| Removing a piece cannot add a parameter | reduced `parameter_names` may not be longer than full | `ValueError` |
| The reduced model is a restriction of the full one | `nested_out_of_sample()` requires strict subset of `parameter_names` | `NotNested` |
| The floor is named | the floor is a `mech.params.Target` — observable, assay, value, units, source — and a bare float is refused | `FloorNotNamed` |
| The floor came from this observable's assay | `Target.assay` is compared with `Observable.assay` | `FloorNotNamed` |
| A dimensional floor cannot bound a dimensionless effect | `scored="absolute"` requires identical units; `scored="relative"` requires a floor in `RELATIVE_FLOOR_UNITS` | `ValueError` |
| Free scalars ≤ independent targets | `require_free_scalar_gate()`, and `nested_out_of_sample()` runs it with one target per held-out fold | `FreeScalarGateFailed` |

The frozen-model check is the one worth dwelling on, because it is the difference between a
rule and a paragraph. A model is not asked whether it is frozen; it is **measured**. It is
refit on the same rows with the response column tilted by a fixed 5% cosine, and if no
fitted parameter moves by more than 1e-6 of itself, it is frozen however it describes
itself. A model that declares **zero** free scalars is a different thing and is allowed: it
is reported as parameter-free, which is a strength.

## 3. The floors

Two, and only two, are measured on this platform. Both are re-exported from
`generator/panel_experiment.py` rather than re-typed.

| Registered name | Floor | Units | Assay | Source |
|---|---|---|---|---|
| `reporter_activity` | **0.146** | relative CV | 96-well plate reader, mCitrine 480/530, OD600-corrected | `OBSERVED_ACTIVITY_CV` — total plate-to-plate CV across 28 matched conditions on the 2026-07-22 and 2026-08-03 plates |
| `growth_rate` | **0.0117** | 1/h — **absolute** | slope of log OD600 against time | `MEASURED_GROWTH_RATE_SE` — median across 210 real wells with a fittable window, IQR 0.0079–0.0200 |

`0.0117 /h` is absolute on purpose: it is 3% of a healthy growth rate and 12% of one slowed
to 0.10, so turning it into a CV at the point of use changes what it means. `ablate()`
refuses to score a relative effect against it.

### 3.1 The observables that have no floor here, and the refusal each returns

`REFUSED_FLOORS` in `mech/ablation.py`. Each entry is the message the caller gets, and each
names the missing measurement rather than a substitute.

| Observable | Why there is no floor |
|---|---|
| `product_ceiling` | an FBA ceiling is an LP bound, not a reading. Route C2 scored one against the plate CV; name the analytical method that would measure the product and its replicate CV |
| `titre` | no titre assay has been run on these strains. The 22.2% that gets reached for is the entry-flux law's leave-one-**strain**-out generalisation error (`bridge/stress_diversion.py`), which is a cross-strain prediction error and not any instrument's repeatability |
| `product_content` | absorbance on a plate reader sums lycopene with β-carotene and reads above the stoichiometric ceiling. An HPLC method and its own CV would close it |
| `flux` | a GEM flux is solved, not read. 13C-MFA confidence intervals would be the floor and none has been run here |

## 4. Using it

```python
from ystwin.mech import ablation as ab

block = ab.reporter_block("20260803_ER&oxidativestress_Replicate3.xlsx", "UPRE1")

result = ab.ablate(
    full_model=ab.MaturingReporter(),        # with the piece
    reduced_model=ab.SaturatingReporter(),   # without it
    observable=ab.endpoint_fold_induction(block),
)                                            # floor resolved from the registry, or refused
print(result.report())
```

```
ABLATION  reporter_activity [fold over the zero-dose control] -- endpoint fold induction at 5 mM, read at 4 h
  full     MaturingReporter         6 free scalars   value 2.70103
  reduced  SaturatingReporter       5 free scalars   value 2.64348
  both refit on the same 175 rows; neither is frozen
  effect   0.0217716 (relative)
  floor    0.146 relative CV -- assay: 96-well plate reader, mCitrine 480/530, OD600-corrected
           generator/panel_experiment.py::OBSERVED_ACTIVITY_CV -- total plate-to-plate CV across 28 matched conditions on the 2026-07-22 and 2026-08-03 plates
  verdict  0.149x the floor -- DOES NOT CLEAR THE FLOOR
```

A block writing its own ablation supplies three things: a model **with** the piece, its
restriction **without**, and an `Observable` naming what moves, the rows both are fitted to,
and the instrument that would measure it. If the observable is not `reporter_activity` or
`growth_rate`, it must arrive with its own `Target`.

### 4.1 Criterion (e), callable

```python
gate = ab.free_scalar_gate(model, targets)     # reports
ab.require_free_scalar_gate(model, targets)    # refuses: FreeScalarGateFailed
```

Every **fitted** scalar counts as free — that is what fitting means — so the gate counts
`parameter_names` rather than asking a model to grade itself. A target the model was allowed
to search for does **not** count, for `ARCHITECTURE_GAPS.md` 0.3's reason: scoring on
whether *some point in a sweep* reproduces a number is fitting with extra steps, and such a
target in the denominator lets a sweep pay for itself. `Target(..., fitted=True)` declares
one and it drops out of the count.

`nested_out_of_sample()` runs the gate itself, with **one target per held-out fold**, so a
model with more free scalars than the ladder has rungs is refused before it is scored.

## 5. Tier 0.3, run: does the maturation state beat the incumbent's algebra?

`REVISED_BUILD_LIST.md` Tier 0.3 asks for a nested out-of-sample comparison on the existing
data, before Phase 1: *current algebra vs the mech model with all sweeps free,
leave-one-dose-out, parameters counted.* This is that comparison, on the nested pair that
exists today.

**The pair.** `SaturatingReporter` is the incumbent's own form written out so it can be
refit: `k(d) = k_basal + k_max·d/(K_dose + d)` is `stress_panel.module_response`'s saturating
dose term, and `dR/dt = k − λR` is `kinetics.reporter_at_time`'s balance. Five free scalars.
`MaturingReporter` adds the **L8 maturation state** — an immature pool feeding the observed
mature one — parameterised by `tau_mat` so that `tau_mat = 0` is a corner of the fitted box
and reproduces the reduced model exactly. Six free scalars, strictly nested.

**The data.** All twelve committed (plate, construct) blocks: three NewProtocol plates that
carry a per-construct dose ladder × four built strains, 175 rows each, read from
`data/plates/` through `plate/replay.py`. Leave-one-dose-out, both models refit on every
fold. Errors are RMSE on the held-out rows over the mean absolute signal of that fold's
**training** rows, so the denominator carries nothing the fit did not have. The two end
rungs of a ladder are extrapolation by construction, so the table reports the five interior
folds.

| Plate | Construct | `tau_mat` fitted (h) | LODO interior: algebra | + maturation | delta | ×floor | Ablation effect | ×floor |
|---|---|---|---|---|---|---|---|---|
| 20260722 | UPRE1 | 0.614 | 0.1380 | 0.1373 | +0.0007 | +0.005 | 0.0035 | 0.024 |
| 20260722 | UPRE2 | 0.768 | 0.1447 | 0.1426 | +0.0021 | +0.014 | 0.0116 | 0.079 |
| 20260722 | NativeYap1 | 0.272 | 0.8692 | 0.8686 | +0.0005 | +0.004 | 0.0033 | 0.023 |
| 20260722 | AlteredYap1 | 0.188 | 0.8127 | 0.8124 | +0.0003 | +0.002 | 0.0021 | 0.015 |
| 20260803 | UPRE1 | 0.783 | 0.1248 | 0.1178 | +0.0070 | +0.048 | 0.0218 | 0.149 |
| 20260803 | UPRE2 | 0.776 | 0.1040 | 0.0994 | +0.0046 | +0.032 | 0.0151 | 0.104 |
| 20260803 | NativeYap1 | 0.233 | 0.9002 | 0.9004 | −0.0002 | −0.002 | 0.0000 | 0.000 |
| 20260803 | AlteredYap1 | 0.108 | 0.8539 | 0.8541 | −0.0002 | −0.001 | 0.0021 | 0.014 |
| 20260804 | UPRE1 | 3.404 | 0.1193 | 0.1052 | +0.0141 | +0.096 | 0.0297 | 0.204 |
| 20260804 | UPRE2 | 4.150 | 0.1280 | 0.1083 | +0.0197 | **+0.135** | 0.0283 | 0.194 |
| 20260804 | NativeYap1 | 2.070 | 0.6937 | 0.6897 | +0.0040 | +0.028 | 0.0000 | 0.000 |
| 20260804 | AlteredYap1 | 0.655 | 0.6961 | 0.6949 | +0.0012 | +0.008 | 0.0000 | 0.000 |

**The result.**

1. **The ablation clears the floor on none of the twelve blocks.** Its largest effect is
   **0.0297 = 0.204x** the 0.146 floor, on 20260804/UPRE1.
2. **Out of sample the extra state never earns its parameter.** It wins on ten of twelve
   blocks and loses on two, and its best win is **0.0197 = 0.135x** the floor. Seven times
   too small to be seen by the assay that would have to see it.
3. **`tau_mat` is not identified by these traces.** It lands between **0.108 h and 4.150 h**
   — a **38x** spread — and **5.5x** across the three replicates of UPRE1 alone. That is the
   measurement standing behind the revised list's experiment 1.2 (*measure `k_mat` on this
   instrument*), and the reason `MaturingReporter` never presents its fitted value as a
   maturation time.

### 5.1 What criterion (d) says about the state that was added

The reduction rule is `tau/T < 0.146`: on the 4 h plate read these blocks come from, that is
`tau < 0.584 h`. The fitted `tau_mat` straddles it. Four of the twelve blocks (0.108, 0.188,
0.233, 0.272 h) sit **below** the threshold, so the maturation state would be reduced to
algebra there and the ablation is asking whether a state that should not exist helps; eight
sit above it, up to 4.150 h, so the state survives the rule. **All four that reduce are H2O2
blocks and all six DTT blocks keep the state**, which is a second reason the H2O2 rows carry
no weight.

On a 5-day fed-batch (`tau < 17.5 h`) every one of the twelve reduces. The maturation state
does not exist in the vessel that makes the titre under this architecture's own rule.

**What this result is not.**

- It is **not experiment 0.2 / A4.** That is the lag-1 innovation autocorrelation of +0.954
  in `estimator.py`, a different measurement on a different object, and it remains to run.
- The four H2O2 rows are **a comparison between two models that both miss**, not a finding
  about the maturation state. Their interior LODO error is 0.69–0.90 of block scale against
  0.10–0.14 on the DTT blocks, because the 2 and 4 mM peroxide wells go **negative** after
  the dilution correction — `panel_experiment.PanelDataset.n_unusable` documents exactly
  this — and neither a saturating nor a maturing model describes a negative reporter.
- It does **not** say the L8 maturation state is wrong. It says that on the reporter
  channel, on the data in hand, at equal fitting effort, it is below the floor of the assay
  that would have to see it. A block is deleted from the build for failing its ablation
  test, not tuned; a block whose observable does not exist yet is a different verdict.

## 6. Reading a verdict honestly

Three things the protocol deliberately does not do.

**It does not adjudicate a tie.** `clears_floor` is strictly `ratio > 1.0`. An effect equal
to the floor is not above it.

**It does not punish a parameter.** There is no AIC, BIC or F-test here. The out-of-sample
comparison is the answer to "did the extra scalar pay for itself", and it answers it by
holding data out, which needs no assumption about the noise model. The parameter counts are
reported so a reader can see what was spent; they do not enter the verdict.

**It does not know whether an effect is real.** It knows whether it is larger than the noise
of a named instrument. Those are different claims, and conflating them is how a 22.2%
cross-strain generalisation error came to be used as a within-strain per-condition floor.

## 7. The protocol's own constants

Three, all numerical, all `ASSERTED` and registered in `ABLATION_PARAMS`. None of them enters
a prediction, so none is a degree of freedom criterion (e) is about — a check that fires on
5% and on 8% alike has not spent a parameter.

| Name | Value | What it is |
|---|---|---|
| `probe_amplitude` | 0.05 | size of the deterministic tilt used to prove a model responds to its data |
| `probe_moved` | 1e-6 relative | how far a fitted parameter must move to count as having moved |
| `fit_max_nfev` | 20000 | evaluation budget for one least-squares fit; reached by none of the twelve blocks |

The fitted models' bounds and starting points are derived from each frame's own scale rather
than asserted as magnitudes, for the reason `pathway/flux.py` gives about its extrapolation
window: a guard whose bounds are unrelated to the data cannot fire where it matters.
