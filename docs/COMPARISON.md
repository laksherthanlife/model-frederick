# A head-to-head, re-measured

**A ("ours")** — `IGEM/model-v2`, package `ystwin`. A yeast stress digital twin built around
biosensor characterisation on real plates.

**B ("theirs")** — `IGEM/Validation-Notebooks`, a clone of `github.com/IGEM-NUS-2026/Validation-Notebooks`.
A synthetic DBTL benchmark that treats exact Yeast9 dynamic pFBA as a *hidden wet lab*.

Measured 2026-08-26 on this machine. Every number below has the command that produced it.

**A was under active development throughout the measurement**, moving from commit `9f8c953`
(10 commits) to `fe66743` (15 commits) while this was being written. Static counts below are at
**`fe66743`**; the two full test-suite runs were executed at **`9f8c953`**, where the suite
collected 1,445 rather than the current 1,452. Nothing in the conclusions turns on those seven
tests. B was static at `e8682ca` (5 commits).

---

## 1. Verdict

The earlier comparison is stale in both directions, and B has moved too — it is at 5 commits,
not 1, and only 2 of its 11 requirements are pinned.

A closed most of the process gap: it has history, pins, four nulls, a hashed split manifest,
three repo-level audits, and a scored held-out forward prediction **on real plates that it
loses** — the single most credible artefact in either repository. B still has 2.3× the source,
11 acceptance-gate audits, 69,120 LP solves, 3×3 seed replication, world-separated deltas and
runtime leakage guards that A has no equivalent of.

But B ships none of its evidence: `results/` is empty and gitignored, and **18 of its 126 tests
fail on a fresh checkout**, including the leakage gates that are its main methodological claim.
A's suite is 1,445/1,445 green with data, and 1,265 passed / 180 skipped / **0 failed** with
every data path pointed at nothing.

A is still plainly behind on three things: single-seed results, one generator configuration,
and — worst — none of its new hardening is imported by the scripts that produce its headline
tables.

---

## 2. The re-measured table

| | A (ours) | B (theirs) | Earlier claim |
| --- | ---: | ---: | --- |
| Source LOC, non-test | **18,152** (src 13,090 + scripts 5,062) | **41,982** (all in `scripts/`) | A 9,284 / B 41,982 |
| Source files | 87 (67 src + 20 scripts) | 88 | — |
| Test LOC / files | **14,591** / 89 | **2,376** / 23 | A 11,875 / B 2,376 |
| Tests collected | **1,452** | **126** | A 1,204 / B 126 |
| Suite, this machine | **1,445 passed, 0 failed, 0 skipped** (at `9f8c953`) | **95 passed, 13 skipped, 18 failed** | not measured |
| Suite, no external data | **1,265 passed, 180 skipped, 0 failed** | n/a — B has no offline mode | not measured |
| Statement coverage | **93%** (4,972 stmts, 335 missed, 67 modules) | no coverage config | "81%" (was 97%) |
| Real wet-lab data | **Yes** — 4 NewProtocol Synergy exports, 2 ER plates, Bio-Rad Cq, 2 ATP plates | **None** | unchanged |
| Held-out predictive validation | **Yes, on real plates** — 84 wells, 4,342 scored points, twin loses | Yes, on its own simulator | A "None" — **overturned** |
| Baselines / nulls | **4 surrogates + skill score** (`analysis/nulls.py`); 3 forecast baselines | 3 outcome-blind selectors + 5 workflow tracks | A "None" — **overturned** |
| Figures (distinct / files) | **5 / 10** + `figures/report.html` | **16 / 26** | A 0 / B 34 |
| Notebooks | 0 | 9 (mean 2.4 code cells, **0 stored outputs**) | A 0 / B 8 |
| Acceptance-gate scripts | 3 (`audit_claims`, `audit_reproducibility`, `audit_determinism`) | 11 `audit_*.py` | A "No" — **partly overturned** |
| Pinned deps | **9 pins, all `==`**, each checked against `importlib.metadata` | **2 of 11** (`cobra`, `swiglpk`) | A "No" / B "Yes" — **both wrong** |
| Commits | **15** | **5** | A 0 / B 1 — **both stale** |
| Tracked files | 253 | 157 | — |
| Results committed | 34 CSVs in `outputs/` | 5 CSVs; `results/` empty and gitignored | — |
| Prose docs | 8,337 lines, 19 files | 7,080 lines, 5 files (66–88% duplicated between three of them) | — |

### Commands

```bash
# LOC and file counts
find src scripts -name '*.py' | xargs cat | wc -l          # A: 18,152
find scripts     -name '*.py' | xargs cat | wc -l          # B: 41,982
find tests -name '*.py' | xargs cat | wc -l                # A 14,591 / B 2,376

# tests
python3 -m pytest --collect-only -q | awk '{s+=$2} END {print s}'   # A: 1,452
python3 -m pytest --collect-only -q | tail -1                       # B: 126 collected
python3 -m pytest -q -p no:cacheprovider                            # both

# A with every asset path pointed at nothing (the fresh-clone condition)
env YSTWIN_IGEM_RESULTS=/nowhere/x YSTWIN_PLATES=/nowhere/x YSTWIN_ATP_SENSOR=/nowhere/x \
    YSTWIN_THERMO=/nowhere/x YSTWIN_QPCR=/nowhere/x YSTWIN_QPCR_RAW=/nowhere/x \
    YSTWIN_YEAST_GEM=/nowhere/x YSTWIN_EC_YEAST_GEM=/nowhere/x \
    python3 -m pytest -q -p no:cacheprovider

# coverage
python3 -m pytest -q --cov=src/ystwin --cov-report=term; coverage report | tail -1

# history, tracking, figures
git log --oneline | wc -l          # A: 15 / B: 5
git ls-files | wc -l               # A: 253 / B: 157
ls figures/ | sed 's/\.[a-z]*$//' | sort -u | wc -l   # A: 6 (5 figs + report) / B: 16
ls scripts/audit_*.py | wc -l      # A: 3 / B: 11
grep -c '==' requirements.txt      # B: 2 of 11 lines
```

Three rows deserve a caveat rather than a victory lap.

**A's 1,452 vs B's 126 is not a like-for-like ratio.** A has 1,609 `assert` statements across
1,452 tests — **1.1 asserts per test**, with only 15 `parametrize` decorators, so the count is
genuine but the unit is tiny. B has 428 asserts across 126 tests — **3.4 per test**, at 18.9
test-LOC per test against A's 10.0. A's style is one-assertion-per-behaviour; B's is fewer,
fatter integration gates. A has roughly **3.8×** B's assertions, not 11×. Several of B's are
integrity gates with no counterpart anywhere in A —
`test_preprocessing_uses_training_statistics_only`,
`test_training_only_pca_ignores_extreme_test_curve`,
`test_regime_labels_are_posthoc_not_curve_generators`,
`test_minimum_baseline_selectors_are_outcome_blind_and_unique`. On *methodological* assertions
per test, B is ahead.

**A's coverage is 93%, not the 81% previously claimed and not the 97% claimed before that.**
The 81% figure came from a stale `.coverage` file that predated `viz/figures.py` being tested.
The weakest module is `qpcr.py` at 57%.

**B's figure count was inflated.** 26 files, but they are 16 distinct figures each emitted as
`.png` and `.svg`; 10 png + 16 svg. The "34" in the earlier table counted something else.

---

## 3. What A does better

**It has real measurements, and B has none.** Four NewProtocol Synergy exports, two 2026-07 ER
plates, Bio-Rad Cq exports and two ATP-sensor plates. B contains no experimental data of any
kind and no comparison to a published measurement — `grep -riE "literature value|measured in
vivo|chemostat|retentostat"` over B's five markdown files returns nothing about yeast
physiology. Its own claim boundary (`README.md:127–132`, four lines of prose) concedes "This is
a synthetic computational validation only."

**It has a scored forward prediction on held-out real data, and it lost.** `docs/research/PREDICTION.md`
§0: condition `estimator.py`'s particle filter on the first 2 h (13 of 25 timepoints) of 84
real wells across 2 plates, forecast the remaining 12, score 4,342 (well, channel, horizon)
points. Paired per-well against a six-point log-linear extrapolation, the twin **loses on OD by
23.0% of level at 2 h** (95% well-bootstrap CI [18.2, 28.3], twin better in 10% of 84 wells) and
ties on mCitrine. Nominal 95% forecast intervals cover **40–49%**. The document then names the
cause — `estimator.py::_propagate` gives growth rate no drift, so a decelerating culture is
extrapolated flat. It does **not** get a second, independent confirmation. This paragraph
used to cite `scripts/run_calibration_nis.py` for "median time-averaged NIS 1.76 on OD, 4.44
on mCitrine against a χ² band of [0.52, 1.63]"; that pair comes from a 2000-particle
research probe on 167 wells over 2 plates at K = 25, which `data/current_claims.json` carries
as `nis.research.shipped` at role `historical` with its binding `pending` — the archived
configuration and innovations were never recovered, so it cannot be re-run. Regenerating the
script gives a different answer: `outputs/nis_channel_summary.csv` reports 319 well-channels
per channel over 4 plates at K = 22, **`tested_well_channels: 0`** and verdict
**INCONCLUSIVE**, every well-channel excluded below the ESS floor (208 of 638 also below the
particle-diversity floor). Its descriptive medians are 45.31 and 19.32 against the K = 22
band [0.499, 1.672], and `outputs/nis_manifest.json` warns in its own field that these are
"per-well K bands; descriptive channel aggregates have no chi-squared null". So the honest
statement is that the calibration gate returns **no test at all**, not a second diagnosis.
The forecast loss stands on its own. What this document still records correctly is that the
aggregate-RMS table and the paired-MAE table *disagree on the winner*, and the right lesson
about pre-registering the scoring rule. **B has no comparable self-inflicted loss on anything
it holds out.**

**It works on a machine that has never seen a plate.** 1,265 passed / 180 skipped / **0 failed**
with every `YSTWIN_*` variable pointed at a nonexistent directory. `scripts/audit_reproducibility.py:704`
(`check_offline_suite`) automates exactly that check, and `--offline-suite` runs the suite twice
and requires the second run to have zero failures and at least as many skips. B has no
equivalent, and does not survive the condition (§4, §7).

**Its splits are structural refusals, not conventions.** `src/ystwin/analysis/splits.py` (1,047
lines) makes a split that straddles a batch axis *impossible to construct*: `make_split(frame,
"heldout_replicate", group_key=("construct","dose_mM"))` raises `SplitNotPossible` with the
reason. It separates `interpolation` from `extrapolation` as different claims, makes
`extrapolation` deterministic (the top rung is the top rung), takes every combination containing
a held-out stressor with it, and refuses to pool constructs challenged with different agents by
delegating to `qpcr.STRESSOR_FOR_CONSTRUCT` rather than trusting the frame's own labels.
`partition_hash` (`splits.py:281`) puts a SHA-256 of the whole assignment in every manifest row.
B's split integrity is asserted in tests
(`test_split_manifest_is_environment_level_and_has_no_overlap`); A's is enforced at the
constructor.

**It refuses rather than reports.** `analysis/uncertainty.py:52` sets `MIN_PLATES_FOR_INTERVAL = 3`
and `fold_change` returns `INCONCLUSIVE` below it — which is why **all eight rows** of A's
headline dose-response table are currently INCONCLUSIVE (`docs/DATA_INVENTORY.md`). `gates/g1_optical.py`
refuses wells outside the linear range. That is a repository that argues itself down.

**Its claim tiering is real and cited.** `docs/CLAIM_BOUNDARY.md` (20 KB) assigns every README
number to Tier 0–3 with a published precedent each (Talts SBC; Peng domain randomization;
Muratore SOB/UCBOG; Kadian SRCC), states plainly that **nothing is Tier 3**, and lists twelve
Tier 1 claims by name. Its earlier description of all such scores as "provably optimistic"
was overbroad: `research/CIRCULARITY.md` now distinguishes expected optimization bias from
unmeasured biological transfer. B's claim boundary in this comparison was four lines.

**Its audits report their own failures.** `outputs/audit_reproducibility.csv`: 102 PASS, **33
FAIL**, 1 SKIP — stale tables, five asset resolvers found by sibling-directory convention, five
outputs produced but never committed. `outputs/audit_determinism.csv`: 44 PASS, 2 FAIL, one of
which is `seed changes the answer: module_recovery — every seed gives the same answer, so the
seed is accepted and not used`. That check — asserting seed *sensitivity*, in both directions,
with a required written reason to declare an entry point insensitive — is something B does not
have at all.

**Its dependency discipline is stricter.** 9 `==` pins, and `check_dependency_pins`
(`audit_reproducibility.py:244`) fails if the pinned version differs from what
`importlib.metadata` reports as installed. It also asserts the solver is GLPK, because a
degenerate LP has a face of optima and each solver picks a different vertex. B pins 2 of 11
requirements and records the versions it happened to run with in prose (`README.md:399` —
"dependency versions recorded in the run were COBRApy `0.31.1`, optlang `1.9.1`, NumPy `2.5.1`,
pandas `2.3.3`, and Python `3.12.13`"). Recorded is not pinned: none of numpy, pandas or optlang
is constrained in `requirements.txt`.

---

## 4. What B does better

### 4.1 Runtime leakage guards, not just test-time ones — **A should steal this**

`scripts/open_ended_baselines.py:22` declares `FORBIDDEN_OUTCOME_COLUMNS` (objective_value,
final_product, feasible, hidden_regime_id, …) and `assert_no_hidden_inputs` (`:232`) *raises* if
any appears in a proposal frame. `deployment_benchmark.py:1130` scans every outgoing bundle for
eight leak terms and marks it unapproved. `audit_open_ended_benchmark_readiness.py:400` writes a
whole `leakage_matrix()` CSV grading each campaign's leak risk in prose. A has structural refusal
inside `splits.py` but nothing that scans an *artefact* for contamination on the way out.

**Cost to A:** low. A `forbidden_columns` check in `analysis/splits.py` and a leak scan in
`build_report.py` is a day. **Verdict: adopt.**

### 4.2 Seed replication in every result-producing run — **A should steal this**

`run_fixed_environment_validation.py:42–43` (`dataset_seeds = (101, 202, 303)`,
`model_seeds = (11, 22, 33)`) and `run_dfba_state_machine_validation.py:49–50`
(`(701, 702, 703)` × `(17, 29, 43)`) each declare a 3×3 grid, with a separate `fast_*` pair for
tests so the CI path is one seed and the reported path is nine. Paired per-split differences are
reported across them (`README.md:359–362`: −0.004959 extrapolation, +0.001970 interpolation,
+0.011622 held-out) — and B says in the same paragraph that one paired comparison per split makes
its uncertainty estimate "a limitation rather than a stable model-ranking claim".

A does not do this. `run_training.py:81,107` runs at `seed=0`. `run_power.py:52,66` at `seed=0`,
`:82,84` at `seed=1`. Only `run_transfer.py:59` sweeps (`range(5) if scale else range(2)`).
**Every one of A's twelve Tier 1 headline numbers is a single draw with no reported spread.**

**Cost to A:** low-to-moderate. The generator already takes a seed; the change is a loop and a
column. `audit_determinism.py` already proves the seeds are live. **Verdict: adopt — this is the
cheapest real improvement available.**

### 4.3 More than one generator configuration — **A should steal this**

B runs two hidden worlds, `baseline_world` and `strong_oxidative_burden_world`
(`run_prospective_dbtl_benchmark.py:45–47`), and reports **world-separated deltas** (+0.976554
baseline, +0.597597 strong oxidative burden) alongside the pooled effect. A has exactly one
generator configuration ξ̄ and says so — `CLAIM_BOUNDARY.md` Tier 0, last row, and
`research/CIRCULARITY.md` §7 Fix 1.

Be honest about how much this buys B: `design_benchmark_exact.py:100` `regime_configs()` builds
its worlds by `replace(base, psy_synthesis_base=…, des_decay_base=…, …)` — **same
`capacity_mode` (`dynamic_congestion_feedback`), same structure, perturbed parameters**. That is
the weak form of domain randomization, over two hand-chosen points rather than a sampled
distribution, and B reports no UCBOG either. Two is nonetheless strictly more than one.

**Cost to A:** moderate. `generator/stress_panel.py` would need its constants lifted into a
sampled configuration. Days, not hours. **Verdict: adopt, after 4.2.**

### 4.4 Hash-frozen protocols and assets

B hashes almost everything: `gem_backend.py:134` SHA-256s the SBML asset and publishes the digest
(`README.md:393–396`); `deterministic_candidate_hash` (`deployment_benchmark.py:265`),
`effective_model_hash` (`final_dbtl_benchmark_freeze.py:124`), `phenotype_fingerprint` (`:162`),
`sha256_frame` over trajectory/constraint/flux tables (`audit_open_ended_benchmark_readiness.py:63`),
and a world hash whose test asserts it changes with hidden parameters
(`test_world_hash_changes_with_hidden_world_parameters`).

A hashes one thing: the split partition (`splits.py:281`). It does not hash the Yeast9 asset it
loads, the generator configuration, or any output table.

**Cost to A:** low for the asset digest (ten lines in `paths.py` plus an audit row); moderate for
a generator-config hash. **Verdict: adopt the asset digest; the rest is lower priority than 4.1–4.3.**

### 4.5 Machine-readable acceptance gates with declared criteria

B writes gate tables — `prospective_pilot_acceptance.csv` (11/11 pass),
`dfba_surrogate_decision_gate.csv` (**3/7 pass, verdict `stop_surrogate_not_solved`**),
`gem_dynamic_capacity_acceptance.csv` (verdict `gem_capacity_feedback_too_weak`, PC1 did not fall
below the requested 99%). Criteria are stated with thresholds, the verdict is a string in a CSV,
and **two of the three named above are refusals**. That is a genuinely good instrument and B used
it to stop itself.

A's nearest equivalents are its biological gates (`g1_optical.py`, `g4_anchor.py`) and its three
repo audits. It has no *methodological* gate — nothing that says "this analysis is not licensed
until criteria 1–7 pass". `run_calibration_nis.py:49` pre-registers `MIN_ESS = 100.0` "before
looking at the result", which is the right instinct at the scale of one constant.

**Cost to A:** moderate. **Verdict: adopt for the two claims A most wants to make (transfer lift,
attribution power), not repository-wide.**

### 4.6 Scale, and whether it matters

B ran 69,120 Yeast9 LP solves in the full prospective benchmark, 18,000 in the dynamic-capacity
grid, 13,824 in the finite-pool run. A has run nothing at that scale.

**Mostly decoration, with one exception.** Scale did not make B's headline effect strong: the
paired hybrid-minus-conventional mean is +0.787075 with bootstrap CI **[+0.083642, +1.565145]**
across **20** world-seed pairs — a lower bound one-ninth of the point estimate, from 20 paired
units. 69,120 LP solves bought 20 independent comparisons. What scale *did* buy B is the thing in
4.2 and 4.3: enough runs to report spread across seeds and worlds at all. That is the part worth
copying, and it is cheap; the LP count is not.

### 4.7 Notebooks with a declared reproducibility contract

scripts/notebook_repro.py (in the other repository) gives each notebook an `ExperimentSpec` naming its required output
files and its regeneration command, and `ensure_outputs` fails with that command printed.
`notebooks/README.md` tabulates all nine. A has no notebooks — for a judge, that is a gap in the
narrative layer regardless of what the code does.

Two honest deductions. B's notebooks average **2.4 code cells** and carry **zero stored outputs
and zero execution counts** (`json.load(nb)['cells']`) — most of their content is markdown with
result tables *transcribed by hand*, so the numbers a reader sees are not regenerated when the
notebook runs. And on a fresh clone they cannot run at all, because the CSVs they declare are
gitignored.

**Cost to A:** low. A's `figures/report.html` (built by `scripts/build_report.py`, every number
read from a committed `outputs/` table, no timestamps, SVGs inlined with per-figure id prefixes)
is arguably the better artefact already. **Verdict: A should add 2–3 notebooks for the judge, and
keep report.html as the source of truth.**

### 4.8 A rival's contract, stated in one sentence

B's `README.md:19–23` — *"fixed environment `e = [I, O, S]` → complete product curve `P(t)` …
The model never receives an environmental time series in the canonical pipeline"* — is a
falsifiable contract with an explicit withheld clause, and A's `CONTRACT.md` says outright that
this sentence is what A could not previously produce. A now has one
(`docs/CONTRACT.md`, borrowed from Kapteyn et al. 2021), and
it is more careful — it fills in all six of S/D/O/U/Q/R and declares `R` **undefined**. This row
is now roughly even; it was not before.

### 4.9 Hand-rolled numerics as a portability guarantee

`tests/test_ode_regime_models.py:26–32` asserts the validation script's source contains no
`scipy`, `sklearn` or `torch`. B's MLPs, PCA and Spearman (`design_benchmark_exact.py:38`,
`spearman_without_scipy`) are numpy by hand. It removes a whole class of version drift. A depends
on scipy (pinned 1.17.0) and
is fine — but B's *test* is a nice pattern: an assertion about the dependency surface, not just
about behaviour.

---

## 5. Where both are weak

**Neither has a Tier 3 claim.** Neither registered a prediction before making the measurement.
A says so explicitly (`CLAIM_BOUNDARY.md`: "Nothing in this repository is Tier 3"); B does not
raise the question.

**Neither reports a UCBOG, and both are circular in the same way.** A's transfer, power and
latent-recovery numbers are fitted to and scored on `generator/stress_panel.py`. B's Experiments
1–3 are fitted to and scored on its own reduced surrogate; the Simulation Optimization Bias is
non-negative in both cases.

**Is A's circularity actually better than B's? No — and B's is marginally better.** B's dFBA
oracle sits on Yeast9, a third-party published model neither team wrote and neither team fits,
so the stoichiometric core is genuinely exogenous. But the *dynamics* on top of it — the
seven-state machine, the seven `gamma_*` and seven `pathway_cap_*` constants, the sensitivities
and decay rates in `GEMCultureConfig` (`gem_backend.py:70–125`, ~55 hand-set floats) — are
asserted by B exactly as `stress_panel.py`'s EC50s and loadings are asserted by A. B has two
settings of them; A has one.
That is the whole of B's advantage on this axis. A's Tier 1b block (D1's `[0, ceiling]` structural
result, the TMFA swing table, the ecYeastGEM sweep) is the mirror image: deterministic
computation on adopted third-party models, no fitting, no SOB. **Neither repository has a claim
that is robust across generator *families*.**

**Both have a headline that their own deeper analysis undercuts.** B's README reports held-out
R² of 1.0000 / 0.9999 / 0.9959 / 0.9982 for Experiments 1A/1B/1C/2B (`README.md:119–124`) — while
its own surrogate audit finds training-only PCA explains **99.45% of trajectory variance with one
component** (`README.md:364–366`), that simple baselines are competitive with the neural models
(`:372–376`), and that the gate reads `stop_surrogate_not_solved` at **3/7 criteria passed**
(`:378–383`). An R² of 1.0000 on a one-dimensional trajectory family is not evidence of a model;
it is evidence of an easy task, and B's own audit says so 250 lines further down. A's README
quoted a 0.82 variance share and "4 biological replicates" that `DATA_INVENTORY.md` and
`NULL_RESULTS.md` have since cut to a null-corrected excess over 0.117 and to n=2 INCONCLUSIVE
for every construct. Both repositories bury the correction below the claim — but A has since
committed the correction to the README itself (`ce08b10`, "Correct the replicate count: every
construct is n=2, not just three of four"), and B has not.

**Both keep results outside git that their prose depends on.** B's `results/` is 8 empty
directories and `data/*.csv` is gitignored with 5 of 12 force-added; A's audit flags 5 tables
"produced but never committed" and 3 PNGs that no script writes.

**Neither has CI.** No `.github/workflows` in either. A's three audits are built to be CI steps
(`docs/AUDITS.md` says so) and are not wired to one; B's 11 audits cannot run in CI at all,
because their inputs are gitignored.

**Both have machine-specific paths.** B hardcodes `/Users/julianedberthartono/Jelly/igem/data/raw`
in four places (`gem_backend.py:25`, `prepare_prospective_dbtl_hpc_bundle.py:18`,
`prepare_final_dbtl_hpc_bundle.py:17`, `run_harder_shift_validation.py:1527`) with an env override
on only the first. A resolves five assets by sibling-directory convention. The difference is that
A's `audit_reproducibility.py` **detects and prints all five as FAIL**; B's failure is silent
until a test dies.

**Neither reports a proper scoring rule on a probabilistic forecast.** A's PREDICTION.md §1
specifies CRPS / weighted interval score and a baseline ladder, and §0 reports RMS relative error
and coverage instead. B scores deterministic RMSE throughout.

---

## 6. The three things A should steal next, ranked

**1. Seed replication in every result-producing script.** (§4.2) Highest ratio of credibility to
cost in this document. `run_training.py`, `run_power.py`, `run_transfer.py` and `run_g4.py` all
run at a single seed today. Loop over 3–5, report median and spread, and the twelve Tier 1 rows
in `CLAIM_BOUNDARY.md` stop being point estimates. B's `DFBAConfig` (3×3) is the model. Cost:
hours. Blocks nothing.

**2. A second generator configuration, and world-separated reporting.** (§4.3) Lift
`stress_panel.py`'s EC50s, loadings and crosstalk into a sampled configuration, run the transfer
and power sweeps under two or three of them, and report per-configuration deltas beside the
pooled number. This is `research/CIRCULARITY.md` §7 Fix 1, which A has written down and not run,
and it is the only route from Tier 1 to Tier 2 that needs no bench time. Cost: days.

**3. Runtime leakage guards on artefacts.** (§4.1) `assert_no_hidden_inputs`-style column
refusal in `analysis/splits.py`, and a forbidden-term scan over anything `build_report.py` emits.
A's structural split refusal is the better *design*; what it lacks is the check that runs on the
table after the fact. Cost: a day.

Two runners-up, in case the top three land early: SHA-256 the Yeast9 asset in `paths.py` and put
the digest in an audit row (§4.4); and add 2–3 notebooks over the committed `outputs/` CSVs so a
judge has a narrative entry point (§4.7).

---

## 7. The single sharpest criticism of each

**Of B:** *its methodology is unfalsifiable by anyone but its author.* `results/` is empty,
`data/*.csv` is gitignored with 5 of 12 force-added, and `python3 -m pytest` on a clean clone
gives **18 failed, 95 passed, 13 skipped**. The failures are not logic bugs — they are missing
artefacts and a hardcoded path (`test_asset_discovery_finds_raw_yeast_gem`,
`tests/test_yeast_gem_backend.py:20–22`, asserts a file exists at
`/Users/julianedberthartono/...`, two lines after a fixture in the same file that correctly
*skips* on the same condition). Critically, **the failures include the integrity gates**:
`test_no_culture_leakage_and_real_gem_zero_surrogate`,
`test_calibration_manifest_has_no_exact_test_leakage_labels`, and
`test_public_scientist_tools_do_not_expose_hidden_verifier`. B's central methodological claim —
that its hidden wet lab never leaked into its workflows — rests on three tests that no one else
can run. 69,120 LP solves and 11 acceptance gates are worth nothing a reviewer can verify.

**Of A:** *the hardening is beside the pipeline, not inside it.* Grep the scripts that produce
A's headline tables:

```
run_training.py   : 0 references to nulls, uncertainty, splits, estimator
run_transfer.py   : 0
run_power.py      : 0
run_g4.py         : 0
run_scenarios.py  : 0
```

`analysis/nulls.py` is imported by exactly one non-test site — `audit_determinism.py:346`, as a
determinism specimen. `analysis/splits.py` is imported only by `make_splits.py`, which writes the
manifest; **nothing reads it**, so `SPLITS.md`'s "a result cites a row" is aspirational. And
`estimator.py`, the repository's only posterior, is imported by `run_calibration_nis.py` and
`audit_determinism.py` and by no results script at all — exactly what `CONTRACT.md` admits when
it says Level 1 is "one integration away".

One module is a genuine exception and should be named: `analysis/uncertainty.py` **is** wired
in — `run_d2.py`, `run_sensor_characterisation.py` and `build_report.py` all import it, which is
why the dose-response table reads INCONCLUSIVE rather than reporting eight tight-looking numbers.
That is the pattern the other three modules need, and it is proof the wiring is tractable.

The consequence is concrete: A built four null models and applied them to **one** claim
(growth-dilution, in `NULL_RESULTS.md`); the transfer claim's null is recorded as "an open
question, not yet a finding" and the power tables have none. A built a split module that refuses
leaky partitions and scored nothing through it. It built a filter with a calibration test that
found the filter over-confident, and no reported number uses the filter. **A's methodology has
been demonstrated on itself and not yet turned on its own results** — which, unlike B's problem,
is a wiring job of days rather than a rebuild.

---

## 8. Which repository wins, and with whom

**A judge, in a 20-minute read, would probably still find B more convincing.** It has nine
notebooks, 16 figures, a 1,781-line README with an experiment ladder, 69,120 LP solves, "hidden
wet lab", 11 acceptance gates, and a completed prospective benchmark with a positive paired
effect. A has five figures, no notebooks, a headline dose-response table whose every row reads
INCONCLUSIVE, and a flagship prediction that loses.

One caveat on A's side, and it landed mid-measurement: at `9f8c953` A had **no root README** —
`git status` showed `R README.md -> docs/FINDINGS.md`, so the entry point was a 104 KB research
log. Commit `08f975f` ("Make the README readable, and move the research log behind it") fixed
that. A's README is now 220 lines of plain language with a "What the model can and can't do"
section, and it is a better front door than B's 1,781-line wall. That was the single highest-value
judge-facing change available and it is already done; §6 is what remains.

**A scientist would trust A more, and it is not close.** A is the only one of the two that
measured anything, the only one whose predictions were scored against something it did not
generate, the only one that reports its own audit failures (33 of 136) in a committed CSV, the
only one that runs clean on a machine with no data, and the only one that assigns every claim a
tier with a citable precedent and then says none of them reaches the top. B's most valuable
single act — the `stop_surrogate_not_solved` gate — is exactly the kind of thing A does routinely.

The synthesis is unflattering to both: **B has the better benchmark design and no evidence a
reader can check; A has the better epistemics and has not yet pointed them at its own headline
numbers.** §6 is the shorter list of the two.

---

## 9. Addendum — what changed within the hour, and it matters

This report was committed as part of A's `32d0f5a`, which landed minutes after it was written
and acted on §7's criticism of A directly. Recording it here rather than editing the body,
so the snapshot above stays honest about what was true when it was measured.

**`run_transfer.py` now runs its own null, and A's transfer claim does not survive it.**
From `32d0f5a` and the rewritten `docs/NULL_RESULTS.md`:

| | median transfer R² | p |
| --- | ---: | ---: |
| observed | **+0.0298** | — |
| `rotated_subspace` — random axes, same rank | **+0.0351** | 0.634 |
| `matched_marginals` — no cross-channel structure | −0.0007 | 0.024 |

**The random-rotation null scores higher than the fitted model.** Skill over it is +0.031 on a
scale where 1.0 closes the gap to perfect. The fitted latent axes carry no more transfer
information than arbitrary axes of the same rank; the only null the model clears is the one
saying the channels are correlated at all, which any low-rank projection of correlated data
exploits. Consistent with the Ledermann result in `research/IDENTIFIABILITY.md` — four channels
do not span seven regulons — the oracle column (allowed to see the held-out stressor) tops out
at 0.208, so the ceiling is low for reasons that are not about the fitting.

Three consequences for this comparison.

**§7's criticism of A is now half-answered, and the answer cost A a headline number.** One of
the five unwired scripts is wired; the other four (`run_training`, `run_power`, `run_g4`,
`run_scenarios`) are not. The fix was wired *into the script* rather than run as a one-off probe,
so it re-runs whenever the claim is regenerated — which is the correct form and the one B's
acceptance gates also take.

**§6's ranking is unchanged but its first item just got cheaper to justify.** A null applied to
one claim killed it. The power tables and the identifiability tables still have none. Seed
replication (§6.1) and a second generator configuration (§6.2) remain the top two, and the case
for them is now stronger, not weaker: a single-seed number that has never met a null is exactly
what just failed.

**The README's eleven-fold combination lift is neither confirmed nor refuted.** It comes from a
configuration that still has not been reproduced (`NULL_RESULTS.md`, "Transfer: an open question").
The Tier 1 row in `CLAIM_BOUNDARY.md` for `0.038 → 0.429` should be read as unverified rather
than as surviving anything.

**On the head-to-head, this cuts both ways and net favours A.** A lost a headline number — the
transfer result was one of the twelve Tier 1 claims and one of the more quotable ones. But it
lost it *to its own instrument, in public, in a commit that says so in its title*. B has no
episode of that kind anywhere in its history: its nearest equivalent, the
`stop_surrogate_not_solved` gate, stopped work before a claim was made rather than withdrawing
one already published. A repository that deletes its own best number an hour after building the
tool that could is making the strongest possible argument about how to read its remaining ones.
