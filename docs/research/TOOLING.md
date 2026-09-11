# What to build on, and what to keep

**Question asked.** ~18,000 lines of hand-written source (13,110 under `src/ystwin`, 5,101
under `scripts/`, plus 14,591 lines of tests). Module by module: is there a well-maintained
library that does this better, and what would adopting it cost?

Compiled 2026-08-26. Every maintenance fact below was checked against the GitHub/GitLab API
or PyPI on that date; every capability claim was checked against library source or DeepWiki,
not against a README. Where DeepWiki had no index the source was read directly, and that is
marked. Two claims were verified by running code, and are marked ✅.

---

## 1. Verdict

Roughly **8% of the source is genuine reinvention** — about 1,500 lines across
`bridge/tmfa.py`, `bridge/thermodynamic.py`'s transport energetics, `analysis/latent.py`,
and `estimator.py`'s filter core. The rest is either problem-specific (gates, refusal,
provenance, plate dialects) or thin enough that a dependency costs more than it saves.

The single highest-value adoption is **`pytfa`'s `ThermoModel`** — already an installed
dependency, used for 2% of what it offers, while `bridge/tmfa.py` hand-rolls the other 98%
*less completely* (no proton/water handling, no transport ΔG, no per-compartment
concentration ranges, ΔG° treated as exact, no relaxation for the infeasibility the module's
own docstring describes).

The most *uncomfortable* finding is not a library at all: **Flapjack already computes the
core inversion**, and its formulation is cleaner than this repo's. See §3.

---

## 2. Module by module

Maintenance column: latest PyPI release date / last repo push / open issues. "Verdict" is
adopt · partially adopt · do not adopt · keep.

### 2.1 State estimation

| Module | LOC | Candidate | Maintenance | Verdict | Migration cost |
|---|---:|---|---|---|---|
| `estimator.py` filter core | 407 | **`particles`** (Chopin) | MIT · PyPI `0.4` 2023-11-06 · GH push 2026-02-19 · 13 issues | **Partially adopt** | 2–3 days |
| ″ `_resample` (~10 LOC) | — | `filterpy.monte_carlo` | MIT · PyPI `1.4.5` **2018-10-10** · GH push 2024-02-07 · 85 issues | **Do not adopt** | — |
| ″ | — | `dynamax` | MIT · `1.0.2` 2026-06-25 · active | **Do not adopt** | — |
| ″ | — | `numpyro` / `pyro` | Apache-2.0 · both active | **Do not adopt** | — |
| ″ | — | `stonesoup` (Dstl) | MIT · `1.9.1` 2026-06-24 · 638★ · active | **Do not adopt** | — |
| `Innovation` / NIS | ~60 | *(nothing)* | — | **Keep** | — |

**`particles` — partially adopt.** DeepWiki-verified against the source: it ships
`backward_sampling_ON2`, `backward_sampling_mcmc` (O(N)), `backward_sampling_reject`,
`backward_sampling_qmc`, and two-filter smoothing — i.e. the smoother the contract says is
missing. Resampling schemes: `multinomial`, `residual`, `stratified`, `systematic`, `ssp`,
`killing` — systematic included. Missing channels are handled *natively and better than
here*: `dists.IndepProd` composed with `dists.FlatNormal` (fixed-time missingness) or
`dists.MixMissing` (random missingness) gives exactly the "a missing channel contributes no
likelihood term" semantics the module docstring insists on, without the `if obs.od is not
None` branching.

**The non-negativity justification does not survive contact with the library.** The
docstring says *"a particle filter rather than a Kalman variant because activity cannot be
negative"*. That argues against a **Kalman filter**, not against `particles`. `particles`
ships `TruncNormal`, `Gamma`, `LogNormal`, and `LogD` (a log-transformed distribution) — the
multiplicative random walk in `_propagate` is `LogD(Normal(...))` in one line. The custom
filter is justified by the *diagnostics*, not by the state constraint.

**No library ships NIS.** Checked: `filterpy.stats` has a four-line `NEES` and nothing else;
`particles` has none (DeepWiki: *"explicit filter calibration diagnostics like NIS or coverage
checks are not mentioned"*); Stone Soup ships OSPA and SIAP but not NIS/NEES. `Innovation` —
with its prior-predictive discipline, its ESS travelling alongside the residual, and its
docstring admitting that scaling σ by the *observed* value biases the test toward passing —
is the most defensible original code in the repository. Keep it verbatim.

**Recommended shape:** define the twin as a `particles.state_space_models.StateSpaceModel`
subclass; keep `Innovation` and `_predictive` as a `particles` collector. Net change is
roughly −180 lines of filter machinery, +FFBS smoothing, +an off-the-shelf PMMH/SMC² route
to the reporter kinetics.

### 2.2 Latent state and factorisation

| Module | LOC | Candidate | Maintenance | Verdict | Migration cost |
|---|---:|---|---|---|---|
| `analysis/latent.py` | 150 | `fancyimpute.IterativeSVD` | Apache-2.0 · PyPI `0.7.0` **2021-10-25** · GH push 2023-10-25 | **Do not adopt** | — |
| ″ | — | `sklearn` `FactorAnalysis` / `NMF` / `SparsePCA` | BSD-3 · `1.9.0` 2026-06-02 | **Do not adopt** as written | — |
| ″ | — | `mofapy2` (+MEFISTO) | LGPL-3.0 · `0.7.4` 2026-04-23 · active | **Evaluate, do not adopt yet** | ~1 week |
| ″ → replacement | — | **`scipy.optimize.nnls` + `sklearn.linear_model.Lasso(positive=True)`** | BSD-3 · both active | **Adopt — delete the FA path** | 1–2 days |
| `analysis/transfer.py`, `stress_model.py` | 375 | — | — | **Keep** (downstream of the above) | — |

`fancyimpute/iterative_svd.py` is literally the same algorithm as `latent.py::_factorise`.
It is also **currently broken**: open issue from 2026-02-23, *"argument force_all_finite
incompatible with newer scikit-learn"*, no release in five years. Recommending it would be
recommending a dead library.

`sklearn.decomposition.FactorAnalysis` does **not** accept NaN (DeepWiki-verified against
`_assert_all_finite_element_wise`); neither does `NMF`. So the one thing `latent.py` does
that a library does not is the thing sklearn cannot do. That is not an argument for keeping
it — it is an argument for **not doing factor analysis at all.**

`docs/research/IDENTIFIABILITY.md` §4a already reached this conclusion: `reporter_loadings`
is a *known* operator, so the problem is `y = Lx + ε` with `L` known, `x ≥ 0`, `x` sparse —
non-negative least squares with an L1 penalty, not exploratory factor analysis. That is
`scipy.optimize.nnls` and `Lasso(positive=True)`, both two-line calls in maintained BSD-3
libraries. It replaces 150 lines of iterated SVD *and* escapes the Ledermann bound that binds
only the fitted-loadings path. **This is the highest-value change in the whole document that
also deletes code.**

MOFA's **MEFISTO** variant (GP priors over a covariate, per DeepWiki) is the only library
that speaks to IDENTIFIABILITY §4b — factor analysis with the time axis kept rather than
collapsed. It is worth an experiment. It does not break the Ledermann bound (MOFA relies on
ARD/spike-and-slab for *interpretability*, not identification — the identifiability doc
already records that "identifiab" appears zero times in either MOFA paper), and LGPL-3.0 is a
licence question for a repo that has not stated one.

### 2.3 Uncertainty, nulls, power, design

| Module | LOC | Candidate | Maintenance | Verdict | Migration cost |
|---|---:|---|---|---|---|
| `analysis/uncertainty.py::fold_change` | 523 | `scipy.stats.bootstrap` | BSD-3 · `1.18.1` 2026-08-21 | **Do not adopt** | — |
| ″ | — | `arch.bootstrap` | NCSA · `8.0.0` 2025-10-21 · active | **Do not adopt** | — |
| ″ | — | `pingouin.compute_bootci` | GPL-3.0 · `0.6.1` 2026-03-28 | **Do not adopt** | — |
| ″ `equivalence` | — | `statsmodels.stats.weightstats.ttost_*` | BSD-3 · `0.14.6` 2025-12-05 | **Do not adopt**; cite | — |
| ″ `activity_uncertainty` | — | — | — | **Keep** | — |
| `analysis/nulls.py` | 321 | **`pyunicorn.timeseries.Surrogates`** | BSD-3 · `1.0.0` 2026-07-17 · active | **Partially adopt** (cross-check) | 1 day |
| ″ | — | `neurokit2.signal_surrogate` | MIT · `0.2.13` 2026-03-02 | **Do not adopt** | — |
| ″ | — | `nolds` | MIT · `0.6.3` 2025-11-30 | **Do not adopt — ships no surrogates** | — |
| `analysis/power.py` | 513 | `statsmodels.stats.power` | BSD-3 · active | **Do not adopt** | — |
| `analysis/design.py` + `sensor_selection.py` | 393 | `Pyomo.DoE` | BSD · active | **Do not adopt now** | — |
| ″ | — | `pyDOE3` | BSD-3 · 2026-02-09 | **Do not adopt** (classical DOE, not model-based) | — |
| `analysis/splits.py` | 1,051 | `sklearn.model_selection` group splitters | BSD-3 · active | **Partially adopt** | 2–3 days |

**No maintained Python library does a cluster bootstrap.** Verified three ways.
`scipy.stats.bootstrap` supports percentile/basic/BCa, `paired=`, and `batch=`, and
DeepWiki confirms *"no explicit parameter or mechanism … to directly support resampling of
clusters"*. `pingouin.compute_bootci` handles paired/unpaired only, and lacks BCa.

⚠️ **A DeepWiki answer I had to correct.** Asked about `arch`, DeepWiki replied that
`IndependentSamplesBootstrap` *"specifically supports clustered resampling where the
resampling unit is a group"*. I read `arch/bootstrap/base.py` directly: the classes are
`IIDBootstrap`, `IndependentSamplesBootstrap`, `CircularBlockBootstrap`,
`StationaryBootstrap`, `MovingBlockBootstrap`, `MOONBootstrap`. `IndependentSamplesBootstrap`
is for independent samples of *unequal length* (two-sample comparisons), not for clustered
data. There is no cluster bootstrap in `arch`.

So `fold_change`'s two-stage resample — plates with replacement, then wells within each drawn
plate — has no off-the-shelf equivalent. Its `MIN_PLATES_FOR_INTERVAL = 3` refusal has no
equivalent anywhere at all. **Keep it.**

`statsmodels` *does* ship TOST (`ttost_mean`, `ttost_ind`, `ttost_paired`, `ztost_*`,
`binom_tost`, `tost_proportions_2indep`, `tost_poisson_2indep`). None of them takes a
confidence interval and a margin, which is the form `equivalence()` needs, because the
interval it must test is the cluster-bootstrap one. Keep the function; cite
`statsmodels.stats.weightstats` in the docstring so a reader knows the convention is standard
and not invented here.

`pyunicorn.timeseries.Surrogates` is the real overlap on `nulls.py`. Source-read (not
indexed on DeepWiki): `correlated_noise_surrogates` (FT phase randomisation — the exact
Theiler 1992 construction `_phase_randomise_1d` implements), `AAFT_surrogates`,
`refined_AAFT_surrogates` (IAAFT), `twin_surrogates`, plus `test_pearson_correlation` and
`test_threshold_significance`. BSD-3, v1.0.0, Python 3.14 supported. **But** it pulls
`igraph`, `h5netcdf[h5py]`, `tqdm` and a Cython build to get one function this repo has in 12
lines. Recommendation: *validate* `phase_randomised` against `correlated_noise_surrogates`
in a test, and adopt IAAFT from pyunicorn only if a null ever needs the amplitude
distribution preserved as well as the spectrum. Do not take the dependency now.

The brief's premise that **`nolds` ships surrogate methods is wrong** — its only modules are
measures.py (external), datasets.py (external), examples.py (external), and there is no surrogate function in any of
them. NeuroKit2's `signal_surrogate` offers IAAFT and random shuffle only — no plain FT
phase randomisation, no `rotated_subspace`, no `matched_marginals`, no `shuffled_loadings`.
Three of the four surrogates in `nulls.py` are specific to *this* problem's null hypotheses
and exist nowhere else.

`statsmodels.stats.power` covers t-tests, ANOVA, chi-square and normal approximations —
DeepWiki confirms **none of them cover clustered or hierarchical designs**. `power.py` fits
one slope per plate and tests the slopes, and its `DEFAULT_WELL_CV = 0.052` is *measured*
from the zero-dose wells. No closed-form power library reproduces that. Keep. (The
`PowerEstimate` Wilson interval is a nice touch that closed-form power functions cannot give
you at all.)

`Pyomo.DoE` is real and maintained: FIM by finite-difference sensitivities, D/A/E/modified-E
optimality, `Experiment` objects with `unknown_parameters`/`experiment_inputs`/
`experiment_outputs`/`measurement_error` suffixes, defaults to IPOPT+MA57. That is a heavy
apparatus for `design.py`, whose FIM is `L' S^-1 L` for a *known* `L` — three lines of numpy.
**But** it is the right tool for the design question the repo has not asked yet:
IDENTIFIABILITY §6 wants a densely-sampled early time course, and "which timepoints" is a
model-based optimal-design problem over a dynamic model with sensitivities. File it as the
tool for that, not as a replacement for `design.py`.

`splits.py` at 1,051 lines is the largest single reimplementation candidate in `analysis/`.
sklearn's `GroupKFold`, `StratifiedGroupKFold`, `LeaveOneGroupOut`, `GroupShuffleSplit` cover
the mechanics, and DeepWiki confirms they raise `ValueError` when a split is infeasible. What
they do **not** have: `SplitNotPossible` carrying a *reason*, `_confounding_reason`,
`partition_hash`, `split_manifest`, or the refusal to pool constructs. Delegate the index
arithmetic; keep the manifest and the refusal vocabulary on top. Realistic saving: 300–400
lines.

### 2.4 Thermodynamics and metabolism

| Module | LOC | Candidate | Maintenance | Verdict | Migration cost |
|---|---:|---|---|---|---|
| `bridge/tmfa.py` | 270 | **`pytfa.ThermoModel`** | Apache-2.0 · PyPI `0.9.4` 2022-06-27 · GH push 2026-07-27 · 12 issues | **Adopt** | 3–5 days |
| `bridge/thermodynamic.py::electrical_work` + `dg_prime` | ~130 of 369 | **`equilibrator_api.ComponentContribution.multicompartmental_standard_dg_prime`** | MIT · `0.7.0` 2026-04-16 · GitLab active 2026-08-15 | **Adopt** | 1–2 days |
| `bridge/thermodynamic.py` SEED alias mapping | ~120 | `pytfa.io.base` + `annotate_from_lexicon` | as above | **Partially adopt** | 1 day |
| `bridge/equilibrator.py` | 140 | — | — | **Keep** (thin, correct wrapper) | — |
| `fba/*` | 582 | `cobrapy` (already used) | LGPL/GPL · `0.32.1` 2026-08-11 | **Keep** | — |

**This is the clearest reinvention in the repository.** `bridge/tmfa.py`'s docstring derives
the Henry 2007 formulation from scratch and builds it with raw `optlang` variables. pytfa —
*already an installed dependency, already pinned in `pyproject.toml`* — implements exactly
that formulation and more. Source-read from the installed 0.9.1:

| `tmfa.py` builds by hand | pytfa provides | pytfa additionally provides |
|---|---|---|
| `lnC_{met}` variables | `LogConcentration` | LC fixed to `log(10**-pH)` for H⁺; LC = 0 for H₂O |
| `zf_`/`zr_` binaries | `ForwardUseVariable` / `BackwardUseVariable` | `SimultaneousUse` constraint |
| `dG_{rxn}` variable | `DeltaG` | `DeltaGstd` **bounded by ±deltaGRerr**, `DeltaGErr` |
| `dGdef_`, `fwd_`, `rev_`, `sgnf_`, `sgnr_` | `NegativeDeltaG`, `Forward/BackwardDeltaGCoupling`, `Forward/BackwardDirectionCoupling` | `ThermoDisplacement` (Γ), `ForbiddenProfile` |
| one global `[1e-5, 2e-2]` window | per-compartment `c_min`/`c_max` | — |
| — | — | transport ΔG via `find_transported_mets` |
| — | `analysis.variability_analysis` (TVA) | `find_directionality_profiles`, `calculate_dissipation` |
| — | `optim.relaxation.relax_dgo` / `relax_lc` | — |

That last row matters most. `tmfa.py`'s opening sentence is *"Thermodynamic flux analysis,
because fixing concentrations made the model infeasible"*, and pytfa ships two functions
whose entire purpose is diagnosing and relaxing exactly that infeasibility — which reaction's
ΔG° or which metabolite's concentration window is responsible. The repo currently answers
that question by widening a global window and reporting a coverage fraction.

The one thing `tmfa.py` does that pytfa does not is **one binary per direction** rather than
one per reaction, so an inactive reaction's ΔG′ stays free. pytfa's `SimultaneousUse`
constraint (`FU + BU ≤ 1`) is the same idea; verify the equivalence on a small model before
migrating rather than assuming it.

**On the pytfa pin.** `pyproject.toml` pins `pytfa==0.9.1` because *"0.9.1 clears
`MetaboliteThermo.__dict__` inside `__init__`, which is why `thermodynamic.py` reads
`deltaGf_tr` rather than calling `calcDGis()`, and a release that fixes that would silently
return different formation energies."* Verified in the installed source: `__init__` calls
`calcDGis()` once, stores it as `deltaGf_tr`, then **replaces** `__dict__` with a dict that no
longer carries `self.RT`, `MAX_pH`, `MIN_pH`, `Debye_Huckel_B` or `Adjustment` — so a second
`calcDGis()` call would raise. The workaround is correct. **And no such release exists:**
`pytfa/thermo/metabolite.py` has not been touched since 2019-01-22, three releases before the
current 0.9.4. Moving 0.9.1 → 0.9.4 is safe on this specific point and should be done as part
of the migration.

**equilibrator.** `ComponentContribution` exposes `multicompartmental_standard_dg_prime` —
transport across a membrane with a potential and a pH difference, the thing
`thermodynamic.py::electrical_work` computes from `FARADAY_KJ_PER_V` and a hand-built
`moves` dict. Worse, `bridge/equilibrator.py::reaction_energy` currently **refuses** any
multi-compartment reaction moving more than one net charge (`if len(compartments) > 1 and
abs(net) > 1.0: return None`) — a refusal that exists only because the hand-rolled electrical
term is not trusted, and that the library function makes unnecessary. `standard_dg_prime_multi`
additionally returns a **covariance matrix** across reactions, which is what a TFA model
should use for `DeltaGstd` bounds instead of treating ΔG° as exact.

### 2.5 I/O, physics, reporting

| Module | LOC | Candidate | Maintenance | Verdict | Migration cost |
|---|---:|---|---|---|---|
| `plate/synergy.py` | 298 | `Plateo` | MIT · `0.3.1` 2025-05-29 | **Do not adopt** | — |
| ″ | — | `pyphotometrics` | **does not exist on PyPI** | **N/A** | — |
| ″ | — | Flapjack `registry/upload.py` | MIT · see §3 | **Do not adopt; read it** | — |
| `generator/synergy.py` | 57 | `plate/synergy.py` (in-repo) | — | **Merge — this is internal duplication** | 1 day |
| `plate/dose_response.py` | 120 | — | — | **Keep** | — |
| `growth.py` | 171 | `croissance` (biosustain) | Apache-2.0 · `1.2.0` 2021-08-09 · GH 2026-02-04 | **Do not adopt** | — |
| `reporter.py` | 208 | Flapjack / `wellfare` | see §3 | **Keep** | — |
| `qpcr.py` | 506 | `pyQPCR` and friends | all dead (2013, 2021, 0–2★) | **Keep** | — |
| `viz/figures.py` | 1,179 | matplotlib (already used) | — | **Keep** | — |
| `gates/*`, `diagnostics/*`, `calib/*` | 1,038 | *(nothing)* | — | **Keep** | — |
| `scripts/` pipeline | 5,101 | `snakemake` | MIT · `9.25.2` 2026-08-18 | **Partially adopt** | 3–5 days |
| ″ | — | `nextflow` | Apache-2.0 · `26.04.6` 2026-07-09 | **Do not adopt** (JVM/Groovy) | — |
| `outputs/` versioning | — | `dvc` | Apache-2.0 · `3.67.1` 2026-03-31 · **repo now `treeverse/dvc`** | **Do not adopt yet** | — |

**There is no maintained plate-reader kinetic parser in Python.** Checked exhaustively:

- **Plateo** has a `parsers/` directory, and every parser in it reads *one value per well* or
  a picklist: Labcyte Echo logfiles and picklists, Tecan EVO picklists, AATI fragment
  analyzer, NanoQuant reads, Roche LightCycler qPCR, generic tables. Its `Plate`/`Well`/
  `WellContent` model has no time axis. It cannot represent a kinetic run.
- **`pyphotometrics`** returns 404 on PyPI. The only GitHub match is a 2016 repo with 0 stars
  and an unrelated purpose. The brief's candidate does not exist.
- **Opentrons** is liquid handling; **Benchling** is an ELN with a proprietary API. Neither
  parses instrument exports.
- **Flapjack** has a Synergy parser (§3), but it keys off an `End Kinetic` marker and a
  paired metadata sheet, and refuses anything else.

`plate/synergy.py` locates its header rows *by content*, disambiguates repeat reads of one
fluorophore as `mCitrine[1]`/`[2]` because they differ in gain, refuses a bare fluorophore
name when a run holds several, and detects blank-subtracted blocks from negative values. None
of that exists anywhere else. **Keep it, and treat it as an asset worth publishing
separately** — it is the most reusable single file in the repo for anyone outside this
project.

Its `blank_subtracted` heuristic is the weakest point (GENERALIZATION.md §5.6 already says
so): "negative values mean a blank was already subtracted" is a string-free heuristic where
the export's own block name (`Blank ...`, which `generator/synergy.py` *does* key off) is the
authoritative signal. That is an internal fix, not a library adoption.

**Two Synergy readers in one repo is the duplication that actually matters here.**
`plate/synergy.py` (298 LOC) and `generator/synergy.py` (57 LOC) parse the same instrument's
exports with different assumptions, and `plate/dose_response.py` parses a third dialect.
Before importing anyone else's parser, merge these.

**Snakemake, not Nextflow, not DVC (yet).** The three audit scripts total 1,806 lines and
are genuinely good — `audit_reproducibility.py` checks pins against installed versions, path
resolution outside the checkout, solver identity, and table-vs-code staleness;
`audit_determinism.py` discovers seeded entry points by signature inspection and asserts
seed-sensitivity *in both directions*. Snakemake would replace the **dependency and staleness**
parts (its DAG knows which outputs are older than their inputs, and `--report` produces
provenance) but nothing else: no workflow engine checks that a different seed gives a
different answer, or that a path resolves inside the checkout. Adopt Snakemake for `scripts/`
ordering and staleness; keep both audits. Nextflow is a JVM/Groovy dependency for a pure-Python
project — no. DVC's repository has moved from `iterative/dvc` to `treeverse/dvc`; the project
is still active (push 2026-08-24, `3.67.1` 2026-03-31) but a change of stewarding organisation
is worth waiting out before making it load-bearing for `outputs/`.

### 2.6 Simulation-based inference and identifiability

| Need | Candidate | Maintenance | Verdict | Migration cost |
|---|---|---|---|---|
| SBC / TARP / coverage diagnostics (CIRCULARITY.md) | **`sbi`** | Apache-2.0 · `0.27.0` 2026-08-04 · 857★ · 64 issues | **Adopt the diagnostics** | 2–3 days |
| Amortised posterior over reporter kinetics | `sbi` (NPE/NLE/NRE/FMPE/NPSE) | as above | **Do not adopt yet** | — |
| Structural identifiability of the dynamic route (IDENTIFIABILITY §4b) | **`StructuralIdentifiability.jl`** | MIT · GH push 2026-08-23 · 129★ · 18 issues | **Adopt for one question** | 2–3 days |
| ″ | `SIAN-Julia` | MIT · 2026-02-16 · 30★ | **Do not adopt** (SI.jl subsumes) | — |
| ″ | `strike-goldd` (MATLAB) | GPL-3.0 · 2026-05-22 · 33★ | **Do not adopt** (MATLAB) | — |
| ″ | **`StrikePy`** | 2022-04-21 · 14★ · **dormant 4 years** | **Do not adopt** | — |
| Practical identifiability / profile likelihood | `pyPESTO` | active · 280★ | **Evaluate** | — |

**`sbi` — adopt the diagnostics, not the inference.** DeepWiki-verified: it ships
`sbi.diagnostics.run_sbc`, `sbi.diagnostics.tarp.run_tarp`, and expected-coverage via
`run_sbc(..., reduce_fns=...)`. `docs/research/CIRCULARITY.md` §5 already cites Talts et al.
(arXiv:1804.06788) for Tier 1, and `docs/research/PREDICTION.md` records that the model is
overconfident and that retuning the error bars would hide rather than fix it. SBC is the test
that distinguishes those two diagnoses formally. Adopting `run_sbc` and `run_tarp` against
the existing generator costs a dev dependency and a script; it does not require adopting NPE
or a neural density estimator, and it does not touch `estimator.py`.

Adopting `sbi`'s *inference* would mean fitting a posterior to the generator — which is
precisely the circularity CIRCULARITY.md warns about, now with a neural network in the
middle. Do not, until the generator is a *distribution* over panels (CIRCULARITY Fix 1)
rather than one panel.

**Structural identifiability.** IDENTIFIABILITY §4b proposes the observability route and
marks it *"my analysis, not a verified literature result"*. `StructuralIdentifiability.jl`
answers exactly that question mechanically: `@ODEmodel` with states, inputs and outputs;
`assess_identifiability` returns `:globally` / `:locally` / `:nonidentifiable` per parameter;
`assess_local_identifiability` covers state observability; `find_identifiable_functions`
returns the identifiable combinations, which is the more useful answer when individual
regulon activities are not recoverable but a function of them is. The seven-regulon crosstalk
model with four outputs is small enough to be tractable. There is no Python interface and no
statement of scale limits in the docs — but a one-off Julia script that prints a verdict is a
2–3 day cost against a question the repo has flagged as its biggest upside.

**`StrikePy` is the wrong horse.** It is the Python port of STRIKE-GOLDD, and it has not been
touched since 2022-04-21 (14 stars). The maintained STRIKE-GOLDD is MATLAB. Between a
dormant Python port and a maintained MATLAB toolbox, the maintained *Julia* package is the
better answer.

---

## 3. The Flapjack question

**Repositories.** Flapjack 1 was three repos under `flapjacksynbio`, all now consolidated.
`flapjack2` (MIT, push **2026-08-19**) is the current line, and its README records that
*"Development is now maintained by the Genetic Logic Lab at the University of Colorado
Boulder"*. `pyFlapjack` (MIT, push 2026-07-28, 7★, 12 issues) is the Python client.
`flapjack_api` (MIT, last push 2025-01-07, 24 issues) is Flapjack 1's backend. **Flapjack is
alive and actively maintained.** Citation: Yáñez Feliú et al., *ACS Synth. Biol.* 10(1):183–191,
2021, doi:10.1021/acssynbio.0c00554.

None of these repos is indexed on DeepWiki, so everything below was read from source.

### 3.1 The overlap is real, and larger than expected

`flapjack_api/analysis/analysis.py` (752 lines in FJ1, 781 in FJ2) implements fourteen
analyses: `Velocity`, `Expression Rate (indirect)`, `Expression Rate (direct)`,
`Expression Rate (inverse)`, `Mean/Max Expression`, `Mean/Max Velocity`, `Induction Curve`,
`Heatmap`, `Kymograph`, `Alpha`, `Rho`, `Background Correct`. Three of them compute promoter
activity:

- **indirect** — `dvaldt = savgol_filter(ival(time), pre_smoothing, 2, deriv=1)/dt;
  ksynth = dvaldt / sdensity(time)`. That is `(dF/dt)/OD` on total well fluorescence.
- **direct** — `wellfare.infer_synthesis_rate_onestep(cfp, cod, ttu, degr=self.degr,
  eps_L=..., positive=True)`: regularised linear inversion with a **degradation rate** and a
  **non-negativity constraint**.
- **inverse** — fits a 20-Gaussian synthesis profile through a forward model
  `dp/dt = OD·profile − γ·p` by Tikhonov-regularised `least_squares` with non-negative bounds.

It also has a background model with a media term (no cells) and a strain term (cells, no DNA)
computed per (assay, media, strain) — i.e. **the autofluorescence control this project's
`docs/DATA_INVENTORY.md` records as never measured** is a first-class concept in Flapjack's
schema. And `registry/upload.py` has Synergy, BMG and FluoPi parsers.

### 3.2 The uncomfortable part: their formulation is cleaner than ours ✅

`reporter.py::promoter_activity` computes `dR/dt + (μ + k_deg)·R` where `R = (RFU−bg)/(OD−blank)`.
Flapjack's indirect method computes `(dF/dt)/OD` on *total* fluorescence. These are the same
quantity. Writing `F = R·X`:

```
dF/dt = X·dR/dt + R·dX/dt = X(k_synth − (μ+k_deg)R) + R·μX = X·k_synth − k_deg·R·X
⟹  k_synth = (dF/dt)/X + k_deg·(F/X)
```

**The growth rate cancels.** μ does not appear on the right-hand side for *any* k_deg. It
does not appear in the maturation case either: `Fi = (dFm/dt + k_deg·Fm)/k_mat`, then
`k_synth = (dFi/dt + (k_mat + k_deg_i)·Fi)/X`.

✅ **Verified numerically.** On a simulated culture whose growth falls 0.5 → 0.05 /h under a
*constant* promoter activity of 3.0, with k_deg = 0.08:

| Route | recovered k_synth | noisy RMSE (40 reps, 2% OD / 3% RFU) |
|---|---|---|
| naive RFU/OD | 5.22 → 11.85 (a spurious 2.3× "induction") | — |
| `promoter_activity` (estimates μ) | 2.980 ± 0.020 | 0.240 |
| `(dF/dt)/X + k_deg·(F/X)` (no μ) | 2.977 ± 0.015 | **0.209** |

Both invert the confound correctly. The μ-free form is **13% lower RMSE** because it takes
one derivative (of F) instead of two (of R, and of log OD for μ), and because it never
incurs the R–μ correlation that `analysis/uncertainty.py::activity_uncertainty` needs 400
Monte-Carlo draws to handle.

**This does not mean the repo's estimator is wrong** — it is correct, and it produces μ,
which the gates, D2's dilution-confound regression, `estimator.py` and the reported tables
all need. But it means the honest claim is narrower than "we correct for dilution and the
field does not." The Alon convention `(dGFP/dt)/OD` — which GENERALIZATION.md §4 Rank 1
describes as *"the incumbent, widely used estimator omits the (μ + k_deg)·R dilution term
entirely"* — **is** dilution-corrected for a stable FP. What it omits is k_deg, not dilution.
That sentence should be corrected before it appears in a write-up, and the μ-free form should
be added to `reporter.py` as a second, cheaper estimator whose agreement with the existing one
is a free consistency check.

### 3.3 What Flapjack does not have

Read across the whole of `analysis/`, `registry/` and `plot/` in both FJ1 and FJ2:

- **No uncertainty of any kind.** `grep` for std / sem / confidence / bootstrap / interval /
  quantile across `analysis/` returns only background-subtraction thresholds and a
  `normalize_mean_std` helper. Every expression rate Flapjack returns is a bare point
  estimate. There is no cluster bootstrap, no interval, no `estimable` flag.
- **No refusal.** `bg_correct` prints `'No background data to subtract'` and returns `{}`,
  then subtracts zero — silently, exactly the failure `observation.py::correct_inner_filter`
  raises on.
- **No gates.** No optical-linearity check, no blank-well detection, no three-valued verdict.
- **No replicate discipline.** The plate is not the unit of replication anywhere in the
  schema's analysis path.
- **Known defects in the live path.** FJ2 ships `docs/analysis-smoothing-report.md` (dated
  2026-08-03), written by its own maintainers: in `velocity` and `expression_rate_indirect`,
  *"Both methods contain a block that smooths the measurement series and then never uses the
  result"*, and *"Choosing lowess changes nothing about the result."*
- **A GPL-3.0 dependency inside an MIT project.** The `direct` method requires
  `git+https://github.com/ibis-inria/wellFARE.git` — **GPL-3.0, PyPI `1.0.4` released
  2018-07-12, repository last pushed 2020-07-14.** One of Flapjack's three expression-rate
  methods rests on an unmaintained GPL package installed from a git URL, in an MIT repo.
- **Operationally heavy.** Django 3.0.5 (EOL), `channels` 2.4, Postgres, Redis, Docker
  Compose, WebSocket-based analysis dispatch. To compute one expression rate you stand up a
  server.

### 3.4 Verdict: interoperate, do not build on, do not stay ignorant

**Do not build on it.** The dependency is a Django/Postgres/Redis stack for a package whose
entire runtime dependency list is numpy/scipy/pandas/openpyxl/cobra. The analysis engine has
maintainer-documented dead code in exactly the two methods this project would use, and no
uncertainty machinery to attach the repo's contribution to. Inverting the dependency —
running Flapjack and calling out to `ystwin` for intervals — is worse: it makes every result
contingent on someone else's server.

**Do not stay separate either.** Flapjack's data model is the closest thing genetic-circuit
characterisation has to a standard, and this repo currently speaks only "one exporter, one
plate map, four construct names" (GENERALIZATION.md §3.4). Two concrete steps:

1. **Speak the schema.** `Study → Assay → Sample → Measurement`, with `Signal`, `Media`,
   `Strain`, `Vector`, `Dna`, `Chemical`, `Supplement` around it. `marpaia/flapjack-data`
   (MIT, push 2026-07-15) is exactly this: the Flapjack entities as plain Python dataclasses
   plus a `Storage` protocol, an `InMemoryStorage`, an optional SQLAlchemy/Postgres backend,
   an optional numpy/scipy analysis engine, and a `parity/` harness that diffs its outputs
   against a live Flapjack server. Zero stars and unproven — but it is a ~200-line
   dependency-free way to make `recover_activity`'s input type Flapjack-shaped, versus
   standing up Django. Read it before writing a schema by hand.
2. **Publish the comparison.** Run `ystwin` and Flapjack's three methods on the same four
   NewProtocol plates and put the four numbers side by side, scored with
   `analysis/nulls.skill_score`. That is a Tier 1 result available today with no bench work,
   and it is the honest way to state what this project adds: k_deg and maturation, a cluster
   bootstrap, and a refusal — on top of an inversion the field already performs.

**Chi.Bio** (`HarrisonSteel/ChiBio`, 27★, push 2026-02-10) is not a competitor: it is the
control software for turbidostat hardware with integrated optics. It is relevant for a
different reason — IDENTIFIABILITY §6 asks for 5–10 minute sampling over the first two hours,
which a plate reader on a 30-minute cadence cannot give and a Chi.Bio reactor can. File it as
a hardware route, not a software one. **DeCoDe** has no repository matching the
characterisation problem; it is protein-library design and is not relevant here.

---

## 4. What to adopt now, ranked

**1. `pytfa.ThermoModel` replaces `bridge/tmfa.py`.** *First step:* upgrade the pin to
`pytfa==0.9.4` (safe — `thermo/metabolite.py` unchanged since 2019-01-22) and write one
test that builds a `ThermoModel` on the glycolysis subsystem of Yeast9 and asserts the same
feasible growth rate as `ThermoModel.build(...).growth_rate` reports today. If they agree,
delete `tmfa.py`'s constraint construction and keep only `_apply` / `flux_range` /
`measurement_value` as a thin layer over pytfa. Cost 3–5 days; removes ~230 lines; gains TVA,
relaxation, transport ΔG, per-compartment ranges and ΔG° error bounds.

**2. Delete the factor-analysis path in `analysis/latent.py`.** *First step:* implement
IDENTIFIABILITY §4a's "cheap first test" as written — take `reporter_loadings` for the four
buildable channels, generate k-sparse non-negative module activities, and measure support
recovery with `scipy.optimize.nnls` and `Lasso(positive=True)` for k = 1, 2, 3. If support
recovery works at k ≤ 2, `fit_latent` and `select_dimension` are replaced by a known-operator
solve that is not Ledermann-bound, and 150 lines go away. Cost 1–2 days. **Highest
value-per-line in the document.**

**3. `equilibrator`'s `multicompartmental_standard_dg_prime` replaces `electrical_work`.**
*First step:* pick the three multi-compartment reactions `bridge/equilibrator.py` currently
refuses (`abs(net) > 1.0`), compute each both ways, and compare. Then delete the refusal.
Cost 1–2 days.

**4. `sbi`'s SBC and TARP diagnostics.** *First step:* `pip install sbi` as a dev
dependency; write scripts/run_sbc.py (proposed) that draws 200 panels from `generator/stress_panel.py`,
runs the existing fit on each, and passes the rank statistics to `sbi.diagnostics.run_sbc`.
This tests CIRCULARITY.md Tier 1 formally and gives PREDICTION.md's overconfidence finding a
citable diagnostic. Cost 2–3 days; no change to `src/`.

**5. `particles` under `estimator.py`.** *First step:* express the twin as a
`StateSpaceModel` subclass with `PX` built from `LogD(Normal(...))` for activity and
`IndepProd(Normal, MixMissing(Normal))` for the two channels, and check the posterior against
the current filter on `tests/test_estimator.py`'s fixtures. Then wire FFBS. Cost 2–3 days;
gains the smoother the contract names as missing. Keep `Innovation` exactly as it is.

**6. `StructuralIdentifiability.jl` on the dynamic route.** *First step:* write the
seven-regulon crosstalk model with four outputs as an `@ODEmodel` and run
`find_identifiable_functions`. Whatever it returns — even "nothing beyond one combination" —
converts IDENTIFIABILITY §4b from "my analysis, not a verified literature result" into a
mechanically checked verdict. Cost 2–3 days.

**7. sklearn group splitters under `analysis/splits.py`.** *First step:* replace
`_components`/`_candidates`/`_choose` with `StratifiedGroupKFold` / `LeaveOneGroupOut` and keep
`SplitNotPossible`, `partition_hash` and `split_manifest` on top. Cost 2–3 days; removes
300–400 lines.

**8. Snakemake over `scripts/`.** *First step:* one `Snakefile` with a rule per existing
script, `input:`/`output:` naming the CSVs in `outputs/`. Keep all three audits. Cost 3–5 days.

**9. Merge `generator/synergy.py` into `plate/synergy.py`** and take the authoritative
`Blank ` block-name signal for `blank_subtracted` from the former. Not a library adoption —
but it is the internal duplication most likely to produce two different answers from one file.
Cost 1 day.

---

## 5. What to keep hand-rolled, and why it earns its keep

**`Innovation` and the NIS diagnostic (~60 lines).** No maintained Python library ships NIS.
`filterpy.stats.NEES` is four lines in a package whose last release was 2018; Stone Soup has
OSPA and SIAP but not NIS; `particles` has neither. More than that, the *discipline* is
original: computing the predictive moments before reweighting, carrying ESS alongside every
residual, and documenting in the docstring that scaling σ by the observed value biases the
test *toward passing*. A library would give you the number and not the caveat.

**`analysis/uncertainty.py::fold_change` (523 lines).** Verified against scipy, arch and
pingouin: none does a two-stage cluster bootstrap, and none has a concept of refusing to
report. `MIN_PLATES_FOR_INTERVAL = 3` with a `method` string explaining the refusal, and
`excludes_unity` returning `False` when the interval is absent, are the load-bearing lines of
the whole repository. Adopting a library here would mean adopting a library that returns a
narrow interval from two clusters.

**`analysis/power.py` (513 lines).** `statsmodels.stats.power` has no clustered design, and
no closed-form power function can be told that the unit of analysis is a per-plate slope or
that `DEFAULT_WELL_CV = 0.052` was measured from zero-dose wells. `PowerEstimate` subclassing
`float` so a power carries its own Wilson interval is a better answer than any library gives.

**`plate/synergy.py` (298 lines).** Nothing parses kinetic Synergy exports —
not Plateo, not Opentrons, not `pyphotometrics` (which does not exist). Flapjack's parser
handles one export shape and refuses the rest. Content-located headers, repeat-read
disambiguation by gain, and a `raw_channel` that refuses a double blank subtraction are all
original and all correct.

**The gates and diagnostics (`gates/*`, `diagnostics/*`, `calib/*`, 1,038 lines).** Three-valued
verdicts carrying their own thresholds have no library equivalent in any ecosystem I checked.
This is the part of the codebase most worth other people copying.

**`analysis/nulls.py`'s non-time-series surrogates (321 lines).** `rotated_subspace`,
`matched_marginals` and `shuffled_loadings` are nulls for *this repository's specific claims*
— they encode which structure a given claim is supposed to rely on. `phase_randomised` alone
has a library equivalent, and even that is 12 lines against a Cython/igraph dependency.

**`qpcr.py` (506 lines), `viz/figures.py` (1,179 lines), `analysis/splits.py`'s refusal
vocabulary.** The qPCR ecosystem in Python is dead (nothing above 2 stars, nothing since 2021).
`figures.py`'s rule — *no error bar without an interval*, stated in words on the figure's face
— is an editorial position, not a plotting problem.

**`scripts/audit_determinism.py` (614 lines).** Discovering seeded entry points by signature
inspection and asserting seed-sensitivity in *both* directions, with `seed_sensitive=False`
declarations that fail if the function turns out to be sensitive after all, is a testing idea
I did not find in any library. Snakemake and DVC track provenance; neither catches a seed
threaded halfway.

---

## 6. Interoperability: what this should speak

Nothing here currently speaks any external format. Ranked by cost-to-value:

**1. Flapjack's measurement schema (highest value, lowest cost).** `Study → Assay → Sample →
Measurement` + `Signal` + the registry entities. `marpaia/flapjack-data` (MIT) supplies it as
plain dataclasses with no server. Making `recover_activity`'s input type Flapjack-shaped costs
a day and makes every result in this repo readable by anyone in the genetic-circuit
characterisation community. It also forces the `Media`/`Strain`/`Vector` distinction that
makes the missing autofluorescence control (`vector=null, strain=set`) a *named absence*
rather than a footnote.

**2. iGEM InterLab units — MEFL and microspheres.** GENERALIZATION.md §4 Rank 2 already makes
this case and it is repeated here because it is the cheapest item in either document: two
calibration plates (fluorescein series, silica microsphere series) per reader and gain, per
Beal et al. 2018 (PLoS ONE 13(6):e0199432) and Beal et al. 2020 (Commun Biol 3:512). Every
number this repo reports is currently `RFU · OD⁻¹ · h⁻¹` on one PMT at one gain — a quantity
no other lab can reproduce. `calib/od.py` and `calib/od_linearity.py` are already the right
shape to hold the calibration; only the plates are missing.

**3. SBML — already spoken, worth stating.** `cobra.io.read_sbml_model` is in the dependency
list and Yeast9 is read as SBML. The gap is the other direction: the *reporter* model
(`dR/dt = k_synth − (μ + k_deg)R`, with maturation) is expressible as a four-species SBML
model, and doing so would let it be simulated by anything in the SBML ecosystem, checked by
`StructuralIdentifiability.jl` (§2.6), and fitted by pyPESTO via PEtab. `python-libsbml`
(sbmlteam, active 2026-03-09, 6 open issues) is the writer. Medium cost, high leverage —
it is the same artefact that unblocks items 6 and 7 of §4.

**4. SBOL — worth speaking, not worth generating.** `pySBOL3` (MIT, 46★, push 2025-03-26, 47
open issues) and `SBOL-utilities` (MIT, 22★, push 2025-10-26, **122 open issues**) are
maintained but not thriving. SBOL describes the *constructs* — `UPRE1`, `UPRE2`,
`NativeYap1`, `AlteredYap1` — which this repo currently identifies by bare strings in
`qpcr.STRESSOR_FOR_CONSTRUCT` and `plate/layout.py`. The right level of ambition is an
**optional SBOL URI field on the reporter profile**, resolvable against SynBioHub
(BSD-2-Clause, 83★, active 2026-08-19), not an SBOL document generator. Low cost, and it is
what makes a result citable against a part rather than a nickname.

**5. The iGEM Registry — cite, do not integrate.** The Registry has no maintained programmatic
API worth depending on; SynBioHub is the machine-readable route to the same parts. Record
BBa_ identifiers as metadata on the reporter profile alongside the SBOL URI.

**What not to speak.** PEtab is tempting once an SBML reporter model exists, but it is a
parameter-estimation exchange format and this project's estimation problem is a filter, not a
fit — revisit only if pyPESTO's profile likelihood is adopted for practical identifiability.
HDF5/NetCDF for `outputs/` would break the "a stranger can read the tables" property that
`audit_reproducibility.py` exists to protect; CSV is the right choice and should stay.

---

## 7. Summary table

| Verdict | Modules / needs |
|---|---|
| **Adopt** | `pytfa.ThermoModel` (→ `bridge/tmfa.py`) · `equilibrator.multicompartmental_standard_dg_prime` (→ `electrical_work`) · `scipy.optimize.nnls` + `Lasso(positive=True)` (→ `analysis/latent.py`) · `sbi.diagnostics` (SBC/TARP) · `StructuralIdentifiability.jl` (one question) |
| **Partially adopt** | `particles` (filter core + FFBS, keep `Innovation`) · sklearn group splitters (→ `splits.py` mechanics) · `snakemake` (→ `scripts/` DAG, keep audits) · `pyunicorn` (cross-check only) · Flapjack schema via `flapjack-data` |
| **Do not adopt** | `filterpy` (2018) · `fancyimpute` (broken vs sklearn ≥1.6) · `dynamax` / `numpyro` / `pyro` / `stonesoup` · `arch` / `pingouin` / `scipy.stats.bootstrap` for clusters · `statsmodels.stats.power` · `nolds` (no surrogates) · `neurokit2` · `Plateo` · `pyphotometrics` (does not exist) · `StrikePy` (dormant) · `nextflow` · `dvc` (wait out the org move) · Flapjack as a platform |
| **Keep** | `Innovation`/NIS · `fold_change` cluster bootstrap and its refusal · `power.py` · `plate/synergy.py` · all gates and diagnostics · `qpcr.py` · `figures.py` · `audit_determinism.py` |

---

## Appendix: maintenance facts, checked 2026-08-26

| Package | Licence | Latest release | Last repo activity | Open issues |
|---|---|---|---|---|
| `particles` | MIT | `0.4` 2023-11-06 | 2026-02-19 | 13 |
| `filterpy` | MIT | `1.4.5` **2018-10-10** | 2024-02-07 | 85 |
| `dynamax` | MIT | `1.0.2` 2026-06-25 | 2026-08-01 | 79 |
| `numpyro` | Apache-2.0 | `0.21.0` 2026-05-02 | 2026-08-25 | 71 |
| `pyro` | Apache-2.0 | `1.9.1` 2024-06-02 | 2026-08-04 | 284 |
| `stonesoup` | MIT | `1.9.1` 2026-06-24 | 2026-08-24 | 141 |
| `scikit-learn` | BSD-3 | `1.9.0` 2026-06-02 | 2026-08-25 | 2124 |
| `fancyimpute` | Apache-2.0 | `0.7.0` **2021-10-25** | 2023-10-25 | 3 (one is "broken vs new sklearn") |
| `mofapy2` | LGPL-3.0 | `0.7.4` 2026-04-23 | 2026-08-25 | 6 |
| `statsmodels` | BSD-3 | `0.14.6` 2025-12-05 | 2026-08-25 | 2815 |
| `scipy` | BSD-3 | `1.18.1` 2026-08-21 | 2026-08-25 | 1841 |
| `arch` | NCSA | `8.0.0` 2025-10-21 | 2026-08-10 | 51 |
| `pingouin` | GPL-3.0 | `0.6.1` 2026-03-28 | 2026-04-05 | 14 |
| `pyunicorn` | BSD-3 | `1.0.0` 2026-07-17 | 2026-08-24 | 38 |
| `neurokit2` | MIT | `0.2.13` 2026-03-02 | 2026-08-07 | 13 |
| `nolds` | MIT | `0.6.3` 2025-11-30 | 2025-11-30 | 20 |
| `pytfa` | Apache-2.0 | `0.9.4` **2022-06-27** | 2026-07-27 | 12 |
| `equilibrator-api` | MIT | `0.7.0` 2026-04-16 | GitLab 2026-08-15 | — |
| `cobra` | LGPL-2/GPL-2 | `0.32.1` 2026-08-11 | active | — |
| `Plateo` | MIT | `0.3.1` 2025-05-29 | 2026-05-13 | 3 |
| `sbi` | Apache-2.0 | `0.27.0` 2026-08-04 | 2026-08-20 | 64 |
| `snakemake` | MIT | `9.25.2` 2026-08-18 | 2026-08-25 | 1095 |
| `nextflow` | Apache-2.0 | `26.04.6` 2026-07-09 | 2026-08-25 | 441 |
| `dvc` | Apache-2.0 | `3.67.1` 2026-03-31 | 2026-08-24 (**org moved to `treeverse`**) | 201 |
| `StructuralIdentifiability.jl` | MIT | — | 2026-08-23 | 18 |
| `SIAN-Julia` | MIT | — | 2026-02-16 | 1 |
| `strike-goldd` (MATLAB) | GPL-3.0 | — | 2026-05-22 | 0 |
| `StrikePy` | — | — | **2022-04-21** | — |
| `flapjack2` | MIT | — | 2026-08-19 | 6 |
| `pyFlapjack` | MIT | — | 2026-07-28 | 12 |
| `flapjack_api` (FJ1) | MIT | — | 2025-01-07 | 24 |
| `wellFARE` (Flapjack dep) | **GPL-3.0** | `1.0.4` **2018-07-12** | **2020-07-14** | 0 |
| `flapjack-data` | MIT | — | 2026-07-15 | 0 |
| `pySBOL3` | MIT | — | 2025-03-26 | 47 |
| `SBOL-utilities` | MIT | — | 2025-10-26 | 122 |
| `SynBioHub` | BSD-2 | — | 2026-08-19 | 49 |
| `python-libsbml` | LGPL | — | 2026-03-09 | 6 |
| `ChiBio` | — | — | 2026-02-10 | — |
| `pyphotometrics` | — | **not on PyPI** | — | — |
