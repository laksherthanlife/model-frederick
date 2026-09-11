# Can four channels reach more than four modules?

## Verdict

**The ceiling is real, and the intuition that something is being missed is also right — but not where it was being looked for.**

No amount of clever design breaks the dimension bound. Every identifiability theorem in nonlinear ICA and in interventional causal representation learning independently assumes the mixing map is injective, which *forces* observed ≥ latent. Auxiliary richness is not the binding constraint: the 8 stressors × 7 doses design satisfies every counting condition in that literature comfortably. Adding doses, adding stressors, and adding auxiliary structure cannot buy latent dimensions.

What *is* being missed is that this repository asks **two different identifiability questions with two different ceilings and reports them as one.** One path knows the loadings and solves a rank problem. The other fits the loadings and is bound by a far harsher classical limit that nobody has applied to it: **with four channels, a fitted factor model supports at most one factor, not four.** And the path that knows the loadings is throwing that knowledge away when it hands over to the path that does not — which is where the recoverable ground actually is.

---

## 1. The conflation

Two modules in this package both answer to the word "identifiable". They are not the same problem.

**Path A — loadings known.** `analysis/design.py::fisher_information` builds `L' S^-1 L` from `generator/stress_panel.py::reporter_loadings`, and `module_standard_errors` inverts it. `analysis/sensor_selection.py` scores candidate reporter sets the same way. Here `L` is **given** — transcribed from the literature with PMIDs — so this is a linear inverse problem with a known operator, and identifiability is the rank of `L`. *k* channels give rank ≤ *k*, hence ≤ *k* modules. **The README's table (2 → 2 of 7, 4 → 4 of 7, 7 → 7 of 7) is correct within this framing.**

**Path B — loadings fitted.** `analysis/latent.py::_factorise` calls `np.linalg.svd` on the activity matrix and *estimates* the loadings. This is exploratory factor analysis, and it is bound by the Ledermann degrees-of-freedom limit, which nothing in the repository references.

`run_training.py` and `run_transfer.py` — every transfer and "general stress state" result — go through **Path B**. The identifiability tables come from **Path A**. The tables are therefore not a warrant for the transfer results, and no line of code connects them.

## 2. The Ledermann bound, and why it is much worse than four

Ledermann (1937). For *p* observed variables and *m* common factors, the free parameters (loadings less rotational indeterminacy, plus unique variances) must not exceed the distinct entries of the covariance matrix:

```
pm − m(m−1)/2 + p  ≤  p(p+1)/2      ⟺      (p − m)² ≥ p + m      ⟺      m ≤ (2p + 1 − √(8p+1))/2
```

I verified the algebra symbolically: `p(p+1)/2 − [pm − m(m−1)/2 + p]` expands to exactly `½[(p−m)² − (p+m)]`.

| observed channels *p* | 2 | 3 | **4** | **5** | 6 | 7 | 8 | 10 | 15 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bound | 0.44 | 1.00 | 1.63 | 2.30 | 3.00 | 3.73 | 4.47 | 6.00 | 10.00 |
| **max integer factors** | **0** | **1** | **1** | **2** | **3** | **3** | **4** | **6** | **10** |

**Five fluorophores, one spent on the growth reference, leaves four stress channels — which supports one fitted factor.** Not four. The build the README recommends is, on Path B, a one-dimensional stress state.

Two honest qualifications. The bound is a *counting* condition: non-negative degrees of freedom does not imply identifiability, and structural zeros in the loadings change the count — the generic-identifiability refinements are Shapiro (1985) and Bekker & ten Berge (1997), whose metadata I have but whose full text I could not read, so do not attribute specific theorems to them. And `m > p` is unreachable for a Gaussian factor model under *any* of these refinements.

### A minor related finding

`select_dimension`'s docstring says cross-validation avoids the perfect-fit degeneracy that defeats BIC. Mostly it does — across 30 seeds on known-rank data it picked correctly 22–29 times out of 30. But at *p* = 4–5 it occasionally (1–2 in 30) picks *k = p*, where a rank-*p* model of a *p*-column matrix reproduces any masked entry and the "held-out" error is not held out at all. Never observed at *p* = 7. Worth a guard — cap `max_states` at `p − 1`, or refuse `k = n_channels` outright — rather than a rewrite.

## 3. What the literature settles, and it is not encouraging

All conditions below were read from primary sources.

**Nonlinear ICA with auxiliary variables does not help.** iVAE (Khemakhem, Kingma, Monti, Hyvärinen, AISTATS 2020, arXiv:1907.04809) writes the model as `x ∈ R^d`, `z ∈ R^n` with **`n ≤ d` stated explicitly**, and its Theorem 1 assumption (ii) is that the mixing `f` is **injective**. Its counting condition — `nk+1` distinct auxiliary values making an `nk × nk` matrix invertible — is satisfied trivially here: 8 × 7 = 56 stressor-dose cells, plus time and nutrient level, comfortably exceed what is needed for any plausible module count. **Auxiliary richness is not the constraint; the dimension inequality is.** TCL (arXiv:1605.06336) pins feature dimension = data dimension. PCL (PMLR 54:460) assumes `f` bijective from `Rⁿ` onto `Rⁿ`. GCL (arXiv:1805.08651) sets observed = source dimension by construction. Lachapelle et al. (CLeaR 2022, arXiv:2107.10098), whose mechanism-sparsity route is the most interesting relaxation, still states "we assume `d_z ≤ d_x`".

One trap worth knowing: if the dose response is essentially a **Gaussian mean shift** (k = 1, T(z) = z), GCL's Theorem 2 says the variability assumption *cannot* hold and iVAE's Proposition 1 says the linear indeterminacy is *irreducible*. Variance or shape modulation (k ≥ 2) is needed, not just a shifting mean.

**Interventions do not help, for two independent reasons.** Ahuja, Mahajan, Wang & Bengio (ICML 2023, arXiv:2209.11924) require an injective decoder and, for polynomial mixing of degree *p*, `n ≥ Σ_{r=0}^{p} C(r+d−1, d−1)` — observed must be *much larger* than latent. Squires, Seigal, Bhate & Uhler (ICML 2023) assume `p ≥ d` literally, and prove **d interventions are sufficient and in the worst case necessary**. Buchholz et al. (NeurIPS 2023, arXiv:2306.02235) state "at least *d* interventions… **this cannot be weakened even for linear mixing functions**". Varıcı et al. (AISTATS 2024, arXiv:2310.15450) say that without `d ≥ n` plus a diffeomorphism, "identifiability is ill-posed".

And the counting is **per latent node**. A stressor is an exogenous perturbation hitting many modules at once — a multi-node *soft* intervention with unknown, overlapping targets. Seven doses of one stressor is **one target explored at seven strengths, not seven interventions**. No verified theorem covers the multi-node-unknown-target case at all.

**Multi-view sharing does not help.** Anandkumar et al. (JMLR 15:2773, arXiv:1210.7559) Theorem 3.6 requires each view's loading matrix to have **full column rank**, i.e. `d_t ≥ k` for *every* view. Three views of four features each buy `k ≤ 4`, not 12. And MOFA / MOFA+ (Argelaguet et al., MSB 2018; Genome Biology 2020) and Group Factor Analysis (Klami et al., IEEE TNNLS 2015) prove nothing here — the string "identifiab" appears **zero times** in either MOFA paper, and GFA's single occurrence is an admission: "Since U is not identifiable, before averaging we projected each U to a common coordinate system." Both rely on ARD/spike-and-slab priors for interpretability, not identification.

## 4. The escapes that are actually open

### 4a. Stop doing factor analysis — you already know the loadings (do this first)

This is the recoverable ground, and it costs no bench time.

Path B throws away exactly what Path A has. `reporter_loadings` is a **known, literature-grounded operator**. So the problem is not blind source separation; it is

> solve `y = L x + ε` for module activities `x`, with `L` known, `x ≥ 0`, and `x` sparse

which is non-negative least squares with a sparsity prior — not exploratory factor analysis. Non-negativity is free and physically true (a regulon cannot have negative activity), and it is a genuine identifiability constraint, not a regulariser. The one verified result permitting latent ≫ observed is **overcomplete/underdetermined ICA under non-Gaussianity** — Goyal, Vempala & Xiao, "Fourier PCA and Robust Tensor Decomposition", STOC 2014 (arXiv:1306.5825), whose abstract states plainly that "the number of component distributions *m* can be arbitrarily higher than the dimension *n*". Its honest limit: it identifies the **mixing matrix**; recovering the **per-sample activity vector** from `y = Lx` with `m > n` stays underdetermined at every sample and needs sparsity or non-negativity as the extra prior. Here `L` is already known, so only the second half of that problem remains.

**The catch, and it is specific and serious.** Sparsity fails for the module the project most wants. The README records that ESR is read by six of the seven reporters as crosstalk and drives five of them. A module that is on almost always is not sparse, and no sparsity-based recovery will separate it. Sparsity buys the *episodic* modules — UPR, oxidative, heat, iron — and not the ubiquitous one. The existing recommendation to build a dedicated `STRE-general` channel at 0.4% of the design score is therefore not a nice-to-have; it is the only route to ESR under either framing.

**Cheap first test.** Take `reporter_loadings` for the four buildable channels, generate module activities that are *k*-sparse and non-negative, and ask how often NNLS with an L1 penalty recovers the correct *support* among seven modules for k = 1, 2, 3. If support recovery works at k ≤ 2, four channels reach more than four modules in the only sense that matters — knowing which are on. Pure re-analysis; no wet lab.

### 4b. Use the time axis (untested, and the biggest upside)

`analysis/latent.py` collapses a time course into an observations-by-channels matrix and runs an SVD. **It throws time away, and the Ledermann bound is a bound on exactly that discarded-time object.**

The underlying system is dynamic, and the regulons have genuinely different time constants — UPR fast and sustained, ESR a fast transient that adapts within tens of minutes, Hsf1 fast-transient, Rpn4 slow because it is protein-stability controlled, iron slow. In linear system identification an **observable** system with a *single* output can recover *n* internal states: the binding constraint is the rank of the observability matrix `[C; CA; CA²; …]`, not the number of outputs. A dense time course through four channels is a much larger object than a four-column covariance matrix.

Honesty about status: **this is my analysis, not a verified literature result** — the sub-agent tasked with the dynamic/observability and delay-embedding literature did not complete. Two caveats I can state confidently anyway. Observability requires knowing `A` (the crosstalk *and* the kinetics); if `A` is also unknown the counting gets much harder, though a Hankel matrix still carries more structure than a covariance matrix. And **structural identifiability is routinely far better than practical identifiability** — the profile-likelihood literature (Raue et al., Bioinformatics 2009) exists precisely because systems-biology models are formally identifiable and practically not.

Worth pursuing because it is testable for free: the repository already has a mechanistic generator with per-module kinetics. Give the modules distinct time constants, simulate four channels densely, and ask whether a state-space fit recovers more than one module. If it does, the ceiling was an artefact of the static representation. If it does not, that is a real and citable negative.

### 4c. More channels, physically

Fluorescence **lifetime** (FLIM/phasor) is an axis orthogonal to emission spectrum, so it adds channels without adding spectral slots; hyperspectral unmixing with more detection bins than fluorophores is the cheaper cousin. Whether the "five fluorophores" cap is a hard spectral-overlap limit or an artefact of the reader's filter set is worth ten minutes with the instrument's specification, because every route above improves monotonically in *p* and nothing else does.

## 5. What would falsify each route

| Route | Falsified if |
|---|---|
| Known-L sparse recovery (4a) | Support recovery fails at k = 2 on simulated non-negative sparse activities with the real loadings and realistic noise |
| Dynamics / observability (4b) | A state-space fit to densely-sampled 4-channel simulated data recovers no more modules than the static SVD, with modules given genuinely distinct time constants |
| More channels (4c) | The reader cannot resolve lifetime or additional spectral bins |
| Anything auxiliary- or intervention-based | Already falsified — see §3. Do not spend time here. |

## 6. The minimal decisive experiment

If one thing can be run to find out whether more than four modules are reachable with four channels, it is **not** a bigger dose grid. It is a **densely-sampled time course under a stressor known to excite modules with different kinetics**, plus the dedicated ESR channel.

Concretely: two stressors chosen to differ in which fast/slow modules they excite (DTT for UPR, and a heat or osmotic shock for the fast-transient pair), a single mid-range dose each — the responsive window is already known to be ~0.5–1 mM — read at **5–10 minute intervals for the first two hours** and then at the current cadence. The dose axis buys nothing here; the early time axis is where the distinct time constants live and where the current protocol samples nothing. Build `STRE-general` alongside, since ESR is unreachable by any other route.

Then re-run the analysis two ways on the same data: the current static SVD, and a state-space fit. The comparison is the result.

---

## Bibliography

Verified from primary sources unless marked.

**The ceiling**
- Ledermann, W. "On the Rank of the Reduced Correlational Matrix in Multiple-Factor Analysis." *Psychometrika* 2(2):85–93, 1937. doi:10.1007/BF02288062
- Shapiro, A. "Identifiability of factor analysis: some results and open problems." *Linear Algebra Appl.* 70:1–7, 1985. doi:10.1016/0024-3795(85)90038-2 — ⚠️ metadata verified, full text unread
- Bekker, P.A. & ten Berge, J.M.F. "Generic global identification in factor analysis." *Linear Algebra Appl.* 264:255–263, 1997. doi:10.1016/S0024-3795(96)00363-1 — ⚠️ metadata verified, full text unread

**Nonlinear ICA**
- Khemakhem, Kingma, Monti & Hyvärinen. "Variational Autoencoders and Nonlinear ICA: A Unifying Framework." PMLR 108:2207–2217, 2020. arXiv:1907.04809
- Hyvärinen & Morioka. "Unsupervised Feature Extraction by Time-Contrastive Learning and Nonlinear ICA." NeurIPS 29, 2016. arXiv:1605.06336
- Hyvärinen & Morioka. "Nonlinear ICA of Temporally Dependent Stationary Sources." PMLR 54:460–469, 2017
- Hyvärinen, Sasaki & Turner. "Nonlinear ICA Using Auxiliary Variables and Generalized Contrastive Learning." PMLR 89:859–868, 2019. arXiv:1805.08651
- Lachapelle et al. "Disentanglement via Mechanism Sparsity Regularization." PMLR 177:428–484, 2022. arXiv:2107.10098
- Hälvä et al. "Disentangling Identifiable Features from Noisy Data with Structured Nonlinear ICA." NeurIPS 2021. arXiv:2106.09620

**Interventional causal representation learning**
- Ahuja, Mahajan, Wang & Bengio. "Interventional Causal Representation Learning." PMLR 202:372–407, 2023. arXiv:2209.11924
- Squires, Seigal, Bhate & Uhler. "Linear Causal Disentanglement via Interventions." PMLR 202:32540–32560, 2023. arXiv:2211.16467
- Buchholz et al. "Learning Linear Causal Representations from Interventions under General Nonlinear Mixing." NeurIPS 2023. arXiv:2306.02235
- Varıcı, Acartürk, Shanmugam & Tajer. "General Identifiability and Achievability for Causal Representation Learning." PMLR 238:2314–2322, 2024. arXiv:2310.15450
- Brehmer, de Haan, Lippe & Cohen. "Weakly supervised causal representation learning." NeurIPS 2022. arXiv:2203.16437
- Schölkopf et al. "Toward Causal Representation Learning." *Proc. IEEE* 109(5):612–634, 2021. doi:10.1109/JPROC.2021.3058954 — review, contains no identifiability theorem; cite for sparse mechanism shift only

**Multi-view factor analysis**
- Anandkumar, Ge, Hsu, Kakade & Telgarsky. "Tensor Decompositions for Learning Latent Variable Models." *JMLR* 15:2773–2832, 2014. arXiv:1210.7559
- Argelaguet et al. "Multi-Omics Factor Analysis." *Mol. Syst. Biol.* 14(6):e8124, 2018. doi:10.15252/msb.20178124
- Argelaguet et al. "MOFA+." *Genome Biology* 21:111, 2020. doi:10.1186/s13059-020-02015-1
- Klami, Virtanen, Leppäaho & Kaski. "Group Factor Analysis." *IEEE TNNLS* 26(9):2136–2147, 2015. arXiv:1411.5799

**The escape**
- Goyal, Vempala & Xiao. "Fourier PCA and Robust Tensor Decomposition." STOC 2014. doi:10.1145/2591796.2591875. arXiv:1306.5825
- Eriksson & Koivunen. "Identifiability, Separability, and Uniqueness of Linear ICA Models." *IEEE Signal Process. Lett.* 11(7):601–604, 2004. doi:10.1109/LSP.2004.830118 — ⚠️ metadata verified, full text unread
- Raue et al. "Structural and practical identifiability analysis of partially observed dynamical systems by exploiting the profile likelihood." *Bioinformatics* 25(15):1923–1929, 2009

**Unverified** — Anderson & Rubin (1956), "Statistical Inference in Factor Analysis", *Proc. Third Berkeley Symp.* Vol. 5:111–150. Not indexed in Crossref; not retrieved. Verify before citing.
