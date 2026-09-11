# Training on your own generator, and what it licenses

## Verdict

**Synthetic teacher recovery checks an inference procedure under supplied assumptions; it does not validate those assumptions against biology.** The useful distinction is between three different problems:

- **Direct label leakage:** held-out states, controls, future observations, or product targets enter fitting or prediction when the protocol forbids them. The teacher/student tests explicitly challenge these routes.
- **Shared-model inverse crime:** teacher and student share state definitions, observation physics, a compatible model class, and the same metabolic backend and priors. Accurate recovery can be genuine while shared model errors remain invisible.
- **Selection bias:** a configuration is chosen to maximize a finite-sample score and that optimized score is reported on the same sample. **Simulation Optimization Bias** concerns this last problem, in expectation under an assumed distribution.

Jensen's inequality does **not** prove that every held-out transfer R², power figure, or identifiability count overestimates performance on real plates. An independent held-out evaluation can be valid for the declared simulator distribution; the sign and size of its biological reality gap remain unmeasured. This corrects the earlier pointwise sim-versus-real claim on this page.

Sampling more configurations can reduce expected optimization bias under the assumptions below. A **UCBOG** is an approximate upper confidence bound on a candidate's *expected optimality gap under that distribution*, not on its loss on every fresh draw and not on transfer to unknown biology. The authors explicitly caution that "in the current state-of-the-art there is no way to estimate a policy's transferability to a domain from an unknown distribution." The repository already has family and UCBOG checks; their workflow-specific limits are in §4.

---

## 1. The formal object

Let `ξ` be a generator configuration — a crosstalk matrix, a set of EC50s, a noise model — drawn from a fixed `p(ξ)`, and `J(θ, ξ)` a larger-is-better objective. For an iid sample `S_n = (ξ_1, ..., ξ_n)`, write `Ĵ_n(θ) = mean_i J(θ, ξ_i)`. Assume a common feasible set, finite relevant expectations and well-defined global optima (or the corresponding suprema). Each fixed candidate's sample objective must be unbiased for its population objective.

**Simulation Optimization Bias** (Muratore, Gienger & Peters, *IEEE TPAMI* 43(4):1172–1183, 2021, Table 1 and §2.2):

```
b_n = E_S[ max_θ Ĵ_n(θ) ] − max_θ E_ξ[ J(θ, ξ) ]  ≥ 0
      └─ expected sample optimum ─┘   └─ population optimum ─┘
```

The inequality is in **expectation over samples**, not for each fitted model or reported metric, and equality is possible. Unbiased sampling error does not remove optimization bias, but error alone does not imply strict positivity. In particular this is not a theorem comparing an arbitrary held-out simulator score with a real-plate score from a different distribution.

For a fixed candidate `θᶜ` independent of the reference sample, the **optimality gap** and its sample estimate are:

```
G(θᶜ)   = max_θ E_ξ[J(θ,ξ)] − E_ξ[J(θᶜ,ξ)]  ≥ 0
Ĝ_n(θᶜ) = max_θ Ĵ_n(θ)      − Ĵ_n(θᶜ)        ≥ 0
E_S[Ĝ_n(θᶜ)] = G(θᶜ) + b_n                  ≥ G(θᶜ)
```

There is **no pointwise guarantee `Ĝ_n ≥ G`**. The connection theorem (§2.4) says the expected sample gap exceeds the population gap by `b_n`; a realized gap mixes candidate suboptimality and sample variability. It is not the bias itself. A candidate fitted on separate data can be treated conditionally on that fit; using the reference sample to select it invalidates that independence argument. Local fits and clipped negative gaps also need separate qualifications.

This descends from stochastic programming. Mak, Morton & Wood, *Operations Research Letters* 24(1–2):47–56, 1999, prove `E z*_n ≤ z*` (Theorem 1) and `E z*_{n+1} ≥ E z*_n` (Theorem 2) for **minimization**. For maximization the signs reverse: expected sample-optimum bias is nonincreasing with iid sample size under their assumptions, without requiring convexity, unimodality or smoothness. This says neither that adding arbitrary hand-picked structures must help nor that a fitted model's held-out score or bootstrap bound must improve monotonically. Their intuition:

> "In solving the original problem SP, we must find a decision that hedges against all possible realizations of ξ̃. When calculating these lower bounds, we optimize over a subset of ξ̃'s support. **Because of this 'inside information', we over-optimize and, on average, obtain an optimistic objective value.**"

Muratore's CoRL 2018 paper gives the historical framing: the concept was "already formulated under the name 'optimality gap' by the optimization community in the 1990s… [but] has neither been transferred to robotics nor Reinforcement Learning yet." That is the paper's historical context, not evidence that every other field lacks such methods.

## 2. Why "it looked fine in simulation" is not biological validation

OpenAI's automatic domain randomization study (arXiv:1910.07113 — **preprint, no peer-reviewed venue**) illustrates the distinction. The two reported simulation means are relatively close, while the real-world means differ by about twelvefold; the table alone does not establish statistical equivalence in simulation:

| Training distribution | entropy | simulation | reality |
|---|---:|---:|---:|
| Manual randomization | −0.348 | 42.5 ± 0.7 | **2.7 ± 1.1** |
| Automatic, largest | 0.393 | 46.7 ± 0.5 | **32.0 ± 6.4** |

*"The policy trained with manual domain randomization achieves high performance in simulation. However, when deployed on the robot, it fails."* Their hypothesis names the mechanism: *"training on a maximally diverse distribution over environments leads to transfer via emergent meta-learning… **if the training distribution is so large that the model cannot memorize a special-purpose solution per environment due to its finite capacity.**"*

The counterpart result: Peng, Andrychowicz, Zaremba & Abbeel (ICRA 2018, arXiv:1710.06537) randomized 95 dynamics parameters and took real-world success from **0.0 ± 0.0** (single nominal model) to **0.89 ± 0.06**. Their framing is the one that applies here directly: *"policies are prone to exploiting idiosyncrasies of the simulator… Therefore, instead of training a policy under one particular dynamics model, we train a policy that can perform a task under a variety of different dynamics models."* Note their objective is an **expectation over the ensemble**, not a worst case — a risk-neutral average, which is the right default and worth stating explicitly rather than implying robustness guarantees the method does not give.

## 3. The trap on the other side: more noise is not safer

The term "reality gap" comes from Jakobi, Husbands & Harvey, ECAL 1995, LNCS 929:704–720 — and the paper says something more useful than it is usually cited for. Their §7.2 is titled **"The Noise Level Has to Be Right"**, and the failure is **non-monotone in both directions**:

- **Too little noise:** *"Evolution has taken advantage of the fact that, in a zero noise environment, the simulated robot will react identically in similar situations."*
- **Too much noise:** *"It is perhaps less obvious that evolution can also take advantage of **too much** noise in the simulation to produce networks that rely on the extra noise and are thus incapable of reproducing their behaviours in reality."*

Their conclusion: *"Simulation to situation correspondence seems to be maximised when the noise levels of the simulation have similar amplitudes to those observed in reality."*

**This is directly load-bearing for this repository.** The audit recorded in the README found the generator's noise **three to seven times too low** and corrected it against measured plates. That was the right fix and it was made in the right direction. The 1995 result says the target is *calibration to measured noise*, not maximisation — so the correction should stop where it stopped, and any future "widen it for robustness" instinct should be resisted unless the widened level is itself measured. SimOpt (Chebotar et al., ICRA 2019, arXiv:1810.05687) makes the same point from the other end: *"it is often disadvantageous to use overly wide distributions of simulation parameters as they can include scenarios with infeasible solutions that hinder successful policy learning, or lead to exceedingly conservative policies."*

## 4. What the existing generator gates already cover, and what they do not

`generator/panel_gates.py` implements physical, posterior-predictive and diversity gates, and the README's audit history shows them working — they caught a generator with no growth in it and a noise model off by 3–7×. That is real and it is more than most projects have.

These gates check selected physical, noise and diversity properties of generated panels; when a measured reference is supplied, they compare those summaries too. They do not establish that the downstream fitted model transfers to biology. The earlier statement that family checks were "missing entirely" is stale:

- `generator/families.py` supplies structural variants and a loading-noise control. `scripts/run_cross_family.py` scores a **fixed three-state configuration** against the rotated-subspace null. A separate dimension/design search on named training families is evaluated on different named families. Their mean difference changes difficulty as well as selection exposure; it is descriptive, not selection optimism. This is not `run_transfer.py`'s current nested-width pipeline.
- `analysis/optimism.py` supplies the SPOTA-style bootstrap estimator for a **separately sampled candidate**, held fixed independently of fresh reference family sets. Its approximate bound targets that candidate's expected configuration optimality gap under the uniform synthetic sampler, not the named-split solution's optimism, all fitted weights or all pipeline outputs. Full selected-pipeline optimism is pending independent evaluation. The historical unsuccessful null tests remain in [`../CROSS_FAMILY.md`](../CROSS_FAMILY.md); non-rejection establishes neither equivalence nor absence.
- `scripts/run_in_silico_loop.py` is a separate supervised teacher/student workflow. It splits whole episodes, scrambles labels for a negative control, varies observation growth independently, and scores an untrained saturating-equation teacher separately. It does **not** run SPOTA or train the student across the older panel-family registry. `teacher_recovery_passed` and `cross_equation_generalization_passed` are separate from `biological_validation`, which remains false.

Independent test data address test-set leakage; structural challenges expose selected shared assumptions. Neither makes the assumed generator distribution a measured distribution over real yeast.

## 5. A tiered evidence vocabulary

The repo needs a per-claim label, because the tables read as results and will be read as results. Each tier below has a citable precedent, so this is not invented nomenclature.

| Tier | Name | What it means | Precedent |
|---|---|---|---|
| **0** | **Asserted** | A supplied assumption, such as the fixed synthetic control coefficients in `generator/in_silico.py::TeacherParameters`; not a fitted biological law. | — |
| **1** | **Self-consistent in simulation** | Recovery on held-out synthetic data under stated assumptions. A point-estimation RMSE is not formal posterior SBC. | Related but distinct: Simulation-Based Calibration — Talts, Betancourt, Simpson, Vehtari & Gelman, arXiv:1804.06788 |
| **2** | **Robust across declared generator families** | Performance holds across specified structural perturbations. Report the family split and, if computed, the UCBOG with its objective and distribution. | Domain randomization (Peng, ICRA 2018); SOB/UCBOG (Muratore, TPAMI 2021) |
| **3** | **Predicted then measured** | A prediction was registered before the measurement and scored against it. | Sim-vs-Real Correlation Coefficient — Kadian et al., *IEEE RA-L* 2020, arXiv:1912.06321 |

Formal SBC is a posterior-calibration procedure involving repeated prior-predictive draws and rank diagnostics, not a synonym for teacher recovery. Its limitations are nevertheless instructive. Talts et al. state that SBC *"offers no guarantee that the posterior will cover the ground truth for any single observation or that the model will be rich enough to capture the truth at all."* Modrák et al. (*Bayesian Analysis* 20(2), 2025, doi:10.1214/23-BA1404), **Theorem 7**, give a degenerate case: an approximator that ignores the data and returns the prior passes classical rank-based SBC with uniform ranks. Their data-dependent test quantities, including the log-likelihood, address that SBC failure. The current student's RMSE, label-shuffle and causal-access tests address different questions and do not claim posterior calibration.

Their scope disclaimer is the sentence to put next to any Tier-1 claim: *"SBC as a simulation method has no way to inform us about a discrepancy between the process that generated real data and the assumptions of our statistical model."*

## 6. Why a clean simulation report is not reassurance

Frazier, Robert & Rousseau (*JRSS-B* 82(2):421–444, 2020, doi:10.1111/rssb.12356) prove what happens under misspecification, and it is worse than imprecision:

1. The posterior **concentrates** — tightly and confidently — on a **pseudo-true value** `θ*` that minimises summary-statistic distance, not on any truth.
2. That pseudo-true value **depends on the choice of metric**: *"ABC based on two different metrics… will produce two different pseudo-true values… This lies in stark contrast to the posterior concentration result… under correct model specification [which] concentrates on the same true value regardless of the choice of d(·,·)."*
3. Credible sets *"can yield credible sets with arbitrary levels of coverage"*, and the posterior is not even asymptotically Gaussian.
4. Regression adjustment makes it **worse while looking better**: adjusted posteriors have *"much smaller posterior variability and much shorter credible sets… this behavior **gives researchers a false sense of precision**, and leads to poor coverage rates."*

Cranmer, Brehmer & Louppe (*PNAS* 117(48):30055–30062, 2020) state the field-level version: *"**None of these diagnostics address the issues encountered if the model is misspecified** and the simulator is not an accurate description of the system being studied."*

## 7. Follow-ups and the existing protocol

### Fix 1 — State and challenge the training distribution

`generator/stress_panel.py` declares the baseline topology and response laws; `generator/families.py` already varies edges, falling limbs, adaptation, feedback and loadings. Broader perturbations need explicit biological motivation and a declared distribution, not arbitrary extra noise. Use the existing family results as a structural challenge, not as evidence that `run_training.py`, `run_transfer.py`, or the newer teacher/student runner all train over that registry.

An ensemble-trained student would be a separate experiment. Report its held-out families and objective, and distinguish changes in equation structure from changes in assay noise or input waveform. Mak/Morton/Wood's result applies to iid sample-optimum bias under a fixed distribution; it does not guarantee monotone biological transfer or improved scores merely by adding families.

### Fix 2 — Report a UCBOG only with the objective it bounds

The SPOTA procedure (Muratore, Treede, Gienger & Peters, CoRL 2018, PMLR 87:700–713) needs no real data to study an **assumed simulation distribution**. The existing `analysis/optimism.py` implementation:

1. Fits a **candidate** on `n_c` sampled panels.
2. Fits `n_G` **reference** solutions, each on a fresh set of `n_r` panels, initialised from the candidate.
3. Scores both on each reference's panels **with synchronised random seeds**, producing gaps `Ĵ(reference) − Ĵ(candidate)`. It retains raw gaps and clips negative ones for the bootstrap. Clipping prevents a negative estimate but does not turn a local reference fit into a global optimum.
4. Bootstraps the gap samples using the basic one-sided form `Ḡᵘ = 2Ḡ − Q_α[bootstrap means]`.

Published settings include `α = 0.05`, `B = 1000` bootstrap replications and `n_G = 20` references. Interpreting `P(G ≤ Ḡᵘ) ≈ 1 − α` needs the sampling, optimization and bootstrap assumptions; clipping and an atom at zero do not themselves establish coverage. The expected-bias theorem is not a monotonicity guarantee for this finite bootstrap bound. One reference family set supplies one gap; individual wells or overlapping stressor folds are not extra bootstrap units. The existing bound concerns the separately sampled candidate's dimension/design optimality gap, not the named-split choice's optimism or the full selected pipeline. Transporting fitted weights and repeating the original co-dose search would require a different evaluation.

### Fix 3 — Learn which claims the generator gets right

Koos, Mouret & Doncieux (*IEEE Trans. Evol. Comput.* 17(1):122–145, 2013; formal definition in arXiv:1307.1870) define a **transferability function**: *"a function that maps, for the whole search space, descriptors of solutions… to a transferability score that represents how well the simulation matches the reality."* Their argument for why this is the cheap move is exactly right here:

> "learning the fitness function is likely to be harder than learning the transferability function. Indeed, using a machine learning technique to learn the fitness function of a robot **is equivalent to automatically design a simulator**… On the contrary, predicting that a solution will not be transferable can often be done using a few simple criteria… **learning the transferability complements a state-of-the-art simulator instead of reinventing or improving it.**"

This maps onto the repo's existing philosophy almost exactly: it is a **gate**, learned from a handful of real comparisons, that refuses claims in the regions where the generator has been measured to be wrong. Their working budget was about **10 real evaluations**.

### Fix 4 — When real data exists, score the correlation, not the fit

Kadian et al.'s **SRCC** is the Pearson correlation between simulated and real performance across candidates. It answers "does my simulator *rank* things the way reality does?" — which is the question a design tool actually needs. It requires real rollouts, unlike the UCBOG. The reported jump from 0.18 to 0.844 after tuning simulation parameters is the empirical form of §1's theoretical warning: a simulator can look fine and rank almost at chance.

## 8. The cross-generator-family protocol, concretely

**Families.** The original proposal is now represented by the registry in `generator/families.py`; loading noise is deliberately a parameter-perturbation control rather than a structural variant:

1. **Baseline** — the literature topology as it stands.
2. **Edge-dropped** — remove one crosstalk edge (e.g. Rpn4 ← Hsf1).
3. **Edge-added** — add a plausible edge the literature does not record.
4. **Biphasic** — module-specific falling limbs instead of the baseline's shared viability factor.
5. **Adapting** — an endpoint attenuation of ESR, not an integrated adaptation ODE.
6. **Feedback** — a negative feedback edge handled by the family's linear solve.
7. **Loading noise** — loadings redrawn within a declared, asserted bracket, not measured uncertainty.

**Scoring.** `scripts/run_cross_family.py` uses fixed three-state transfer/null comparisons, separately selects a dimension/design on families 1–4 for evaluation on 5–7, and fits a separate SPOTA candidate on sampled families. It refits latent loadings within `leave_one_stressor_out`; it does not transport one fitted observer unchanged between families. Its channel R² comparator is the outer-training mean, and oracle scores are in-sample diagnostics. The named-family split is a structural holdout inspired by the M-open distinction (Bernardo & Smith; Vehtari & Ojanen, *Statist. Surv.* 6:142–228, 2012), not a demonstration that these families span biological truth.

**Reporting.** Use the authoritative JSON fields `transfer_configuration`, `selected_configuration` and `candidate_configuration` in the optimism table, plus `configuration` in the per-family table. The corrected default external candidate records three states for the first two and two for SPOTA; final retained-artifact adoption remains separate. Keep `split_selection_scope`, `bound_selection_scope`, the score target, reference-set resampling unit and coverage warning alongside any numerical result. `selected_pipeline_optimism_status` remains pending. The historical tests in `docs/CROSS_FAMILY.md` detected no superiority over the rotated-subspace null; this non-rejection remains visible without being promoted to equivalence, absence of information or a robust negative biological result.

**The newer loop asks a different question.** Its frozen-run review reports same-teacher recovery passing and the changed-equation challenge failing under the same RMSE criterion. Teacher and student still share the GEM, enzyme priors and observation/control assumptions. Neither this engineering recovery nor the conceptual stress atlas identifies missing biological pathways.

---

## Bibliography

**Simulation optimization bias**
- Muratore, Gienger & Peters. "Assessing Transferability From Simulation to Reality for Reinforcement Learning." *IEEE TPAMI* 43(4):1172–1183, 2021. doi:10.1109/TPAMI.2019.2952353. arXiv:1907.04685
- Muratore, Treede, Gienger & Peters. "Domain Randomization for Simulation-Based Policy Optimization with Transferability Assessment." CoRL 2018, PMLR 87:700–713
- Mak, Morton & Wood. "Monte Carlo bounding techniques for determining solution quality in stochastic programs." *Oper. Res. Lett.* 24(1–2):47–56, 1999. doi:10.1016/S0167-6377(98)00054-6 — note: **R. Kevin** Wood; the arXiv version of the TPAMI paper has a typo
- Muratore, Ramos, Turk, Yu, Gienger & Peters. "Robot Learning From Randomized Simulations: A Review." *Front. Robot. AI* 9:799893, 2022. doi:10.3389/frobt.2022.799893
- Hobbs & Hepenstal. "Is optimization optimistically biased?" *Water Resour. Res.* 25(2):152–160, 1989 — ⚠️ verified only via citing works
- Bayraksan & Morton. "Assessing solution quality in stochastic programs." *Math. Program.* 108(2–3):495–514, 2006 — ⚠️ verified only via citing works

**Domain randomization and the reality gap**
- Peng, Andrychowicz, Zaremba & Abbeel. "Sim-to-Real Transfer of Robotic Control with Dynamics Randomization." ICRA 2018:3803–3810. arXiv:1710.06537
- Tobin, Fong, Ray, Schneider, Zaremba & Abbeel. "Domain randomization for transferring deep neural networks from simulation to the real world." IROS 2017:23–30. arXiv:1703.06907
- OpenAI et al. "Solving Rubik's Cube with a Robot Hand." arXiv:1910.07113 — **preprint, no peer-reviewed venue**
- Mehta, Diaz, Golemo, Pal & Paull. "Active Domain Randomization." CoRL 2019 / PMLR 100:1162–1176, 2020. arXiv:1904.04762
- Jakobi, Husbands & Harvey. "Noise and the Reality Gap: The Use of Simulation in Evolutionary Robotics." ECAL 1995, LNCS 929:704–720. doi:10.1007/3-540-59496-5_337
- Chebotar, Handa, Makoviychuk, Macklin, Issac, Ratliff & Fox. "Closing the Sim-to-Real Loop." ICRA 2019:8973–8979. arXiv:1810.05687
- Muratore, Eilers, Gienger & Peters. "Data-Efficient Domain Randomization With Bayesian Optimization." *IEEE RA-L* 6(2):911–918, 2021. arXiv:2003.02471

**Measuring the gap**
- Koos, Mouret & Doncieux. "The Transferability Approach." *IEEE Trans. Evol. Comput.* 17(1):122–145, 2013. doi:10.1109/TEVC.2012.2185849. Companion: arXiv:1307.1870
- Kadian et al. "Sim2Real Predictivity: Does Evaluation in Simulation Predict Real-World Performance?" *IEEE RA-L* 2020. arXiv:1912.06321 — ⚠️ volume/pages unverified

**Simulation-based inference and misspecification**
- Talts, Betancourt, Simpson, Vehtari & Gelman. "Validating Bayesian Inference Algorithms with Simulation-Based Calibration." arXiv:1804.06788 — **arXiv only, never published**
- Modrák et al. "Simulation-Based Calibration Checking for Bayesian Computation: The Choice of Test Quantities Shapes Sensitivity." *Bayesian Analysis* 20(2), 2025. doi:10.1214/23-BA1404. arXiv:2211.02383
- Säilynoja, Bürkner & Vehtari. "Graphical test for discrete uniformity and its applications in goodness-of-fit evaluation and multiple sample comparison." *Stat. Comput.* 32(2):32, 2022. doi:10.1007/s11222-022-10090-6
- Cranmer, Brehmer & Louppe. "The frontier of simulation-based inference." *PNAS* 117(48):30055–30062, 2020. doi:10.1073/pnas.1912789117
- Frazier, Robert & Rousseau. "Model Misspecification in Approximate Bayesian Computation: Consequences and Diagnostics." *JRSS-B* 82(2):421–444, 2020. doi:10.1111/rssb.12356 — the arXiv record's DOI field is corrupt; use this one
- Hermans, Delaunoy, Rozet, Wehenkel, Begy & Louppe. "A Crisis In Simulation-Based Inference? Beware, Your Posterior Approximations Can Be Unfaithful." *TMLR* 2022 — note the **published** title differs from the arXiv "A Trust Crisis…" (arXiv:2110.06581)
- Lemos, Coogan, Hezaveh & Perreault-Levasseur. "Sampling-Based Accuracy Testing of Posterior Estimators for General Inference." ICML 2023, PMLR 202:19256–19273. arXiv:2302.03026
- Linhart, Gramfort & Rodrigues. "L-C2ST: Local Diagnostics for Posterior Approximations in Simulation-Based Inference." NeurIPS 2023. arXiv:2306.03580
- Ward, Cannon, Beaumont, Fasiolo & Schmon. "Robust Neural Posterior Estimation and Statistical Model Criticism." NeurIPS 2022. arXiv:2210.06564
- Cannon, Ward & Schmon. "Investigating the Impact of Model Misspecification in Neural Simulation-based Inference." arXiv:2209.01845 — **preprint, no venue**
- Vehtari & Ojanen. "A survey of Bayesian predictive methods for model assessment, selection and comparison." *Statist. Surv.* 6:142–228, 2012
