# Cross-family transfer, named-split selection and a separate candidate gap bound

## Verdict and current scope

**No detected superiority is not equivalence.** The historical comparisons below did not
reject any family's `rotated_subspace` null. That retains every unsuccessful comparison,
but does not establish that fitted and arbitrary axes contain "as much information", that
latent structure is absent, or that transfer is impossible. These are synthetic comparisons,
not biological validation. The former stronger negative/equivalence verdict is withdrawn.

`scripts/run_cross_family.py` now records three different configurations rather than
letting one dimension label stand for all of them:

| calculation | authoritative configuration field | scope |
| --- | --- | --- |
| per-family transfer and null | `cross_family_transfer.csv:configuration`; `cross_family_optimism.csv:transfer_configuration` | fixed three states and recommended co-doses |
| named-family dimension/design selection | `cross_family_optimism.csv:selected_configuration` | selected on baseline, edge-dropped, edge-added and biphasic; evaluated on adapting, feedback and loading-noise |
| SPOTA candidate | `cross_family_optimism.csv:candidate_configuration` | a separately selected candidate on sampled families, fixed independently of fresh reference family sets |

The corrected default-budget external candidate uses **three states** for fixed transfer
and the named-split selection, but **two states** for SPOTA. The configuration JSON carries
the actual co-doses, dose ladder, replicate count, observed channels, readers, stressors and
noise settings. Equal labels must not be inferred from `candidate_states`, nor from shared
panel/readers with `run_transfer.py`, whose current width selection is nested.

The transfer target is median stressor-held-out channel R² relative to **outer-training
channel means**; folds exclude the held-out stressor and co-doses containing it, and share
training data. Oracle R² is an in-sample diagnostic. The named-family mean difference is
descriptive and changes family difficulty as well as selection exposure. It is not an
estimate of selection optimism. The approximate SPOTA UCBOG concerns its separate fixed
candidate's expected configuration optimality gap under the declared uniform synthetic
sampler, not the named-split solution, each future loss or unknown biology.

**Selected-pipeline optimism remains pending independent evaluation.** Latent weights are
refitted within each stressor fold, not transported, and this script does not repeat the
original search that produced the recommended co-doses. The target, split and bound scopes,
resampling unit and coverage warning are now explicit output fields.

**Artifact status.** Final governed regeneration/adoption is still pending. The numerical
tables in §§2–5a below are retained as historical results, including the poor scores, all
null non-rejections, sparse gap samples and unstable selections. They are not updated
current results simply because the corrected external candidate has been inspected. No
larger-budget number is a substitute for selected-pipeline validation.

---

## 1. What a family is, and why not a re-draw

`generator/families.py`. A family is a structurally perturbed panel, not the same panel with
different noise. The distinction is the whole point. Perturbing the constants asks whether the
fitter tolerates noise in numbers it was handed; the question that decides whether a transfer
number is about yeast or about `stress_panel.py` is whether it tolerates the *topology* being
wrong. And a distribution to re-draw from largely does not exist:
`docs/research/PARAMETER_SOURCES.md` counts 113 free numbers in the panel, of which **103 have
no value in the literature in any form a reader could cite**.

Seven families, one per structural axis the panel commits to. Each states what it moves and why
that alternative is live, because a robustness claim is worth no more than the plausibility of
the configurations it survived.

| family | what it perturbs | why that is plausible | sourced / asserted |
| --- | --- | --- | --- |
| `baseline` | nothing | the origin the others are displacements from | — |
| `edge_dropped` | removes `heat → proteasome`, weight 0.40 (Hsf1 transcribing RPN4) | the edge is inferred from promoter content. Hahn 2006 (PMID 16556235) confirms an HSE in the RPN4 promoter by mutagenesis, which establishes the site is functional under heat, not that Hsf1 carries a fixed share of proteasome induction under every agent the panel routes through it. Boy-Marcotte 1999 (PMID 10411744) finds the two regulons overlap only slightly; Ciccarelli 2023 (PMID 37467033) calls the coupling compensatory rather than additive | the doubt is sourced; the choice of *this* edge rather than the weaker `iron ← oxidative` is deliberate — that one the panel itself flags as single-source and never replicated, and removing a 0.10 edge would test almost nothing |
| `edge_added` | adds `ESR → oxidative` at 0.25 (Msn2/4 transcribing YAP1) | the panel's own annotation records "no documented Msn2/4 transcription of YAP1" — an absence of evidence, not a measured zero. Msn2/4 sits upstream of most of the network, so this is the most consequential single addition available: every agent drives the ESR, and the edge carries that into the proteasome and iron arms through weights already recorded | the absence is sourced; the weight is asserted, set between the panel's weakest recorded cascade edge (0.10) and its strongest (0.45) so an unmeasured edge is not given more confidence than a measured one |
| `biphasic` | each module gets its own falling limb instead of the single viability factor they share: ESR 1.5× later, proteasome 1.3× later, metabolite pools 0.7× earlier | the baseline is *already* biphasic, but through a factor shared by every module, which freezes the direction of the activity vector above the peak — past the turnover every regulon falls in lockstep and the high-dose wells add no direction. A module-specific limb does not, and that is a difference a subspace method is supposed to notice. Uses `panel_calibration._model`, the closed form that module fits to real ladders | direction sourced (Gasch 2000: the ESR rises under exactly the conditions where global transcription is cut back; Rpn4 is degraded by the proteasome it induces, so impairment stabilises the factor); magnitude asserted — the measured evidence that the limb is not one number is a 14% spread in the growth-halving dose between UPRE1 and UPRE2 under the same agent (1.45 vs 1.65 mM DTT), which shows it varies without saying how far |
| `adapting` | the ESR read at 35% of its sustained level | Gasch 2000 measured the environmental stress response as a transient: it peaks in ten to twenty minutes and decays back toward baseline inside the hour even while the insult persists, while stressor-specific regulons stay on. A plate read hours later reads an adapted ESR; the panel encodes a sustained one | the transient is sourced; the surviving fraction is asserted — a promoter fusion averages over about 1/μ = 2.5 h, so it reports a time-average of a decaying signal rather than either endpoint |
| `feedback` | closes a loop: adds `proteasome → heat` at −0.25 against the recorded `heat → proteasome`, so the heat arm damps itself | Hsf1 is held off by the chaperones it induces — the canonical account of heat-shock attenuation (Krakowiak et al., eLife 2018; Zheng et al., eLife 2016). The panel records the forward half and *cannot* represent the return: `_propagate` resolves modules in one fixed order, so a loop cannot be written down there at all. That is a representational limit, not a claim that no loop exists | the loop's existence and sign are sourced; the return weight is asserted |
| `loading_noise` | redraws every crosstalk loading, each within a factor of 2 at 95%; topology untouched | **the control.** The crosstalk loadings are the largest block of the unsourced parameters, and the statement behind each is qualitative — "KAR2 carries both a UPRE and an HSE" — which fixes the sign and not the coefficient. If transfer moves as far here as under a rewiring, the structural families are measuring parameter sensitivity and the comparison says nothing about topology | the bracket is asserted and is close to the least such a sentence can mean. The reporter's own element stays at 1.0: it is a normalisation of the channel's units, not a measurement |

One axis is deliberately absent: a family perturbing the cascade *weights* continuously. It is
the same axis as `edge_dropped` — dropping an edge is the limit of shrinking it — so it would
add a draw rather than a structure.

Propagation is a linear solve, `(I − DW)a = Dd`, rather than the single ordered pass in
`stress_panel._propagate`. On the literature topology the two agree exactly (the recorded
cascade is acyclic and `_CASCADE_ORDER` is a topological order of it — pinned in
`tests/test_families.py`), and the solve is what lets `feedback` exist at all. Together with
the per-module limb reducing to the shared one at scale 1.0, that makes the baseline family
reproduce `panel_dataset` **reading for reading** — verified to floating point, because
otherwise every displacement below would be measured from a different origin.

## 2. Historical transfer per family, with the null for each

This historical default-budget table used `analysis.transfer.leave_one_stressor_out` at
three latent states over the twelve-stressor panel and eight readers of `run_transfer.py`,
four channels revealed, the recommended eight co-dosed pairs, measured noise and growth-rate
error. `rotated_subspace` null, 20 draws, p floor 0.048. The current outer-training-mean
score and explicit configuration metadata are described above; this old table is not a
reproduction of the current nested-width `run_transfer.py` pipeline.

| family | transfer R² | Δ vs baseline | oracle R² | alignment | null median | p | beats null | skill over null |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| `baseline` | **+0.0499** | — | 0.0921 | 0.540 | +0.0358 | 0.333 | no | +0.015 |
| `edge_dropped` | +0.0434 | −0.0065 | 0.0904 | 0.508 | +0.0310 | 0.381 | no | +0.013 |
| `edge_added` | +0.0511 | +0.0012 | 0.0979 | 0.545 | +0.0486 | 0.476 | no | +0.003 |
| `biphasic` | +0.0512 | +0.0013 | 0.0965 | 0.548 | +0.0404 | 0.381 | no | +0.011 |
| `adapting` | **+0.0247** | **−0.0251** | 0.0665 | 0.447 | +0.0326 | 0.619 | no | **−0.008** |
| `feedback` | +0.0428 | −0.0071 | 0.0831 | 0.529 | +0.0308 | 0.333 | no | +0.012 |
| `loading_noise` | +0.0502 | +0.0004 | 0.0917 | 0.544 | +0.0361 | 0.333 | no | +0.015 |

**No family showed detected superiority in this historical table.** The smallest p is
0.333 and the best skill over a null is +0.015. Those non-rejections neither establish
same-information equivalence nor prove absence of a useful effect.

**Historical comparability note.** `docs/NULL_RESULTS.md` once reported +0.0298 before
commit `7cbdf34` corrected H2O2 induction EC50 from an unsourced 0.5 mM to 0.15 mM
(Goulev 2017, PMID 28418333). The recorded one-constant reversal gave +0.0327, attributing
much of that shift to the EC50 correction. That document later recorded +0.0499 as well.
Both histories remain; neither old value stands in for a new run after changes to the
score baseline and selection protocol. An increased score that still fails to reject its
null is not automatically stronger evidence of absence.

**The oracle column is an in-sample diagnostic, not an independent ceiling estimate.**
The historical range 0.067–0.098 comes from a model allowed to see the nominally held-out
stressor. It can motivate checking the four-channel design, but does not prove how well
all estimators could perform or add held-out evidence to the transfer column.

## 3. Historical family sensitivities, without an equivalence claim

The per-family deltas describe sensitivity of a synthetic score under the declared
perturbations. Null non-rejection does not make those deltas proof about latent structure.

**One loading-noise draw moved less than the structural variants at the default budget.**
Ranked by absolute displacement: `adapting` 0.0251, `feedback` 0.0071, `edge_dropped` 0.0065,
`biphasic` 0.0013, `edge_added` 0.0012, `loading_noise` 0.0004. The control was the smallest
in that run. This does not establish insensitivity to constants in general or an exclusive
role for topology.

**The number may be an ESR axis.** `adapting` touches one module and halves the score. It also
has by far the largest effect on alignment — the fraction of the held-out stressor's variation
lying in the fitted subspace at all — dropping it 0.540 → 0.447, three times the next largest
move (`edge_dropped`, 0.540 → 0.508). That would say the shared direction the model finds is
largely the general stress response, which five of the eight readers load on, and that it is
contingent on the panel's decision to model the ESR as a level rather than a transient. Stated
as a hypothesis rather than a result because it does not replicate at `--scale`; see the caveat
below.

**Removing crosstalk lowers transfer, adding it does not raise it.** `edge_dropped` −0.0065 and
`feedback` −0.0071 both remove effective cascade gain and both lower the score; `edge_added`
+0.0012 adds gain and barely moves it. Consistent with a model whose transfer comes from shared
variance rather than from the specific wiring: taking shared variance away costs something,
adding more of it is already saturated.

**The per-family ranking is not stable between budgets.** At `--scale` (§5a), `adapting`'s
displacement is −0.0070 rather than −0.0251; `biphasic` moves +0.0065. Earlier prose wrongly
called `biphasic` the largest absolute mover and quoted the scale control as −0.0002: the
retained historical table gives `loading_noise` +0.0012. Those prose claims are withdrawn,
not repaired by choosing a different run. Every family still failed to reject its null in
both historical tables. Neither those non-rejections nor a small control displacement
establishes equivalence, biological robustness or a stable sensitivity ranking.

## 4. A configuration chosen on some families, applied to families it never saw

The M-open split of `docs/research/CIRCULARITY.md` §8 searches latent dimension (2, 3 or 4)
and two declared co-dosing designs (blocked or the recommended eight pairs). It does **not**
repeat the pair search that produced those recommended pairs or evaluate the current
nested-width `run_transfer.py` pipeline. Loadings are refitted per stressor fold.

In the historical default run, selection on `baseline`, `edge_dropped`, `edge_added` and
`biphasic` chose **3 states, recommended pairs**. That configuration was applied to
`adapting`, `feedback` and `loading_noise`:

| | historical mean transfer R² |
| --- | ---: |
| families the configuration was chosen on | **+0.0489** |
| families it never saw | **+0.0392** |
| descriptive chosen-on minus held-out difference | **+0.0096** |

The difference was roughly twenty percent of the chosen-on score. At the historical
`--scale` budget the three figures were +0.0462, +0.0355 and +0.0107, with four rather than
three states selected. Both poorer held-out-family scores remain on record.

**This is not selection optimism or the UCBOG.** Fixed named family sets change difficulty
as well as selection exposure; this comparison does not decompose them or estimate the
expected bias from choosing a configuration. Its units are named synthetic families, not
independent wells. In the current schema `selected_configuration` and
`split_selection_scope` identify this calculation, while
`selected_pipeline_optimism_status` explicitly remains pending. §5's candidate is a
different solution selected on sampled families and held fixed for its references.

## 5. The UCBOG, and what it does and does not cover

`analysis/optimism.py`, the SPOTA procedure of Muratore, Treede, Gienger & Peters (CoRL 2018,
PMLR 87:700–713): fit a separate candidate on n_c sampled families and hold it fixed
independently of the reference sets; fit n_G references, each on fresh n_r families and
initialised from that candidate; score both on each reference's own families under
synchronised seeds. Retain raw gaps, clip negative gaps for the bootstrap and use the basic
one-sided form `Ḡᵘ = 2Ḡ − Q_α`. One independent reference family set supplies one gap; wells
and overlapping stressor folds are not additional independent bootstrap units.

The following values are **historical**, not the corrected current candidate's bound:

| | |
| --- | --- |
| **UCBOG** | **+0.0088** |
| mean gap | +0.0048 |
| gap samples | 0, 0, 0, +0.0113, 0, +0.0206, +0.0066, 0 |
| clipped | 0% — no reference scored below the candidate |
| candidate | 2 states, recommended pairs |
| α | 0.05 |
| bootstrap replications *B* | 1000 |
| reference solutions *n_G* | 8 (published setting is 20; `--scale` uses it) |
| candidate families *n_c* | 4 |
| reference families *n_r* | 2 |
| distribution over families | uniform over the seven-family registry |

+0.0088 was about 18% of the historical fixed-configuration baseline score, for scale
only. The numerator and denominator describe different targets, not a measured fraction
of that score attributable to optimism.

**It is not the same quantity as §4's +0.0096, and the near-agreement is a coincidence.**
The UCBOG targets the population optimality gap `max_θ E[J(θ)] − E[J(θᶜ)]` under the
assumed family distribution. It is estimated from same-sample reference/candidate gaps,
whose expectation includes reference-sample optimization bias under the independence and
global-optimum assumptions. §4 instead compares two different fixed sets of families and
also includes their difficulty difference. At the published budget they diverge sixfold
(+0.0018 against +0.0107; §5a). The earlier use of that coincidence as a consistency check
was wrong; neither quantity is a measured sim-to-biology gap.

Three limits, and the first is the one that matters.

**The bound targets an expected gap under the assumed family distribution, not loss on each
new family or real plate.** The reported +0.0088 is an approximate confidence statement for
that narrow objective under the sampling/optimization assumptions, not proof that the score
is free of sampling artifacts. It says nothing about transfer to an unknown biological
distribution. The sampler is uniform over seven hand-built families; those weights are a
declared choice, not measured frequencies of biological mechanisms.

**It targets the separate fixed candidate's configuration optimality gap, not the pipeline
or the named-split selection's optimism.** The solution space is (dimension, design), with
two declared designs rather than the original co-dose search. `leave_one_stressor_out`
refits latent loadings per held-out stressor and does not carry them across families. A
small bound for this restricted objective would not be a small bound on all tables or
proof of a small sim-to-biology gap. Full selected-pipeline optimism remains pending an
independent evaluation of that pipeline.

**The expected-bias theorem is not a monotonicity theorem for this bootstrap bound.**
`tests/test_optimism.py` demonstrates a falling average over twelve seeds on a toy problem,
with some individual increases. That regression example does not establish the same behavior
for every objective, local optimizer or family distribution.

And one property of the gap samples that qualifies the value rather than the method: **five of
the eight are exactly zero**, because the reference selected the same (dimension, design) as the
candidate and the common random numbers then make the two scores identical. With six discrete
options a tie is common, so the gap distribution here is mostly an atom at zero with three
non-zero draws carrying the whole bound. The basic bootstrap does not automatically attain
nominal coverage in this sparse setting. `--scale` supplies more references (20 rather than
8) but 18 of its 20 gaps are zero. A richer candidate space would change the objective being
bounded; it is not a guaranteed statistical repair for sparse gap samples.

One historical symptom worth recording: the candidate chosen on four sampled families was
**2 states**, the named-split configuration of §4 was **3 states**, and at `--scale` those
became 2 and **4**. Selection was not stable across family samples or budgets. Instability
in an argmax is not itself a measurement of expected selection optimism; it is another
reason to retain the distinct configuration identities.

## 5a. Historical run at the published reference-count setting

`--scale`: six doses, six replicates, 40 null draws (p floor 0.024), n_G = 20, n_c = 8. About
18 minutes.

| family | transfer R² | Δ vs baseline | null median | p | beats null | skill over null |
| --- | ---: | ---: | ---: | ---: | :---: | ---: |
| `baseline` | +0.0420 | — | +0.0515 | 0.634 | no | **−0.010** |
| `edge_dropped` | +0.0409 | −0.0012 | +0.0550 | 0.683 | no | **−0.015** |
| `edge_added` | +0.0444 | +0.0024 | +0.0604 | 0.707 | no | **−0.017** |
| `biphasic` | +0.0485 | +0.0065 | +0.0631 | 0.634 | no | **−0.016** |
| `adapting` | +0.0351 | −0.0070 | +0.0500 | 0.634 | no | **−0.016** |
| `feedback` | +0.0392 | −0.0028 | +0.0427 | 0.537 | no | **−0.004** |
| `loading_noise` | +0.0432 | +0.0012 | +0.0491 | 0.634 | no | **−0.006** |

**Every historical skill-over-null value is negative.** The null median exceeded the fitted
score on all seven structures in this run. That observed ordering and the non-rejections
are retained; neither establishes equivalence, reverse-direction significance or a universal
absence of transfer information. The default and larger budgets are different runs, not a
proof that one is simply the small-sample version of the other.

| | default | `--scale` |
| --- | ---: | ---: |
| UCBOG | +0.0088 | **+0.0018** |
| mean gap | +0.0048 | +0.0009 |
| gaps exactly zero | 5 of 8 | 18 of 20 |
| clipped | 0% | 0% |
| candidate | 2 states, recommended pairs | 2 states, recommended pairs |
| n_G / n_c | 8 / 4 | 20 / 8 |
| descriptive named-split difference (§4), not selection optimism | +0.0096 | +0.0107 |

Two readings.

**The reported bound fell from +0.0088 to +0.0018 at the larger budget.** This is a pair of
observed bootstrap estimates, with both candidate and reference budgets changed, not a
demonstration of Mak, Morton & Wood's expected sample-optimum theorem. That theorem does
not require this direction on each run or for this finite bootstrap statistic.

**It fell partly for a reason that is not about bias.** Eighteen of the twenty references
selected the same (dimension, design) as the candidate, so their gap is exactly zero and the
bound rests on two non-zero draws out of twenty. With six discrete options, agreement is the
common case and the gap distribution is nearly an atom at zero. The +0.0018 is an
approximate historical bound estimate for a narrow objective; sparse samples do not
establish nominal coverage. It is not evidence that the *tables* carry only 0.2% optimism.

## 6. Tier

`docs/research/CIRCULARITY.md` §5 asks every claim for a tier. These are simulation-only
method and sensitivity checks. The unsuccessful superiority tests span seven declared
families at two historical budgets; that is the scope of the observed non-rejection, not
a Tier-2 equivalence or absence result. The old statement "fitted axes do no better" and
the promotion of that statement to a robust negative finding are withdrawn. Attaching a
bound for a different fixed candidate does not strengthen a null test into equivalence.
No biological or prospective-validation tier is earned here.

## 7. Reproducing

Use the registered `cross-family` recipe and the governed external-staging/review process
in [`REPRODUCING.md`](REPRODUCING.md) for retained results. `--scale` is a distinct budget,
not a replacement chosen because its historical UCBOG was smaller. Do not overwrite the
retained default with an exploratory larger-budget run or describe either as a bound on
the selected pipeline.

The targeted checks need no scientific artifact regeneration:

```bash
CI=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD/src" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 -B -m pytest -p no:cacheprovider tests/test_run_cross_family.py \
  tests/test_families.py tests/test_optimism.py -q
```

The reporting tests use temporary synthetic fixtures; they do not reproduce or adopt the
historical numerical tables.

### Which parameters these numbers depend on

Recorded because §2 had to correct `NULL_RESULTS.md` for going silently stale, and this file
is exposed to the same failure. Every table above is a function of the eight readers, and those
readers reach only these modules, cascade closure included:

`ESR`, `UPR`, `oxidative`, `heat`, `osmotic`, `proteasome`, `iron`, `dna_damage`

So a change to any **induction EC50, lethal dose, target weight, crosstalk loading or basal
floor touching those eight** moves these numbers and this file must be re-run. That is what
happened to `NULL_RESULTS.md`: the H2O2 induction EC50 reaches `oxidative`.

**No metabolite pool is reader-visible** — not `redox`, `peroxide`, `atp`, `ph` or `nadh`, and
not through the cascade either. Pool parameters therefore cannot move these tables, which is
worth stating rather than rediscovering: the pool EC50s were reworked while this file was being
written (`module_ec50` now refuses to guess a pool arm and raises instead), and the baseline
transfer R² was re-measured at **+0.0499**, unchanged to four decimals. That invariance is
structural, not luck.

## 8. What this leaves undone

- **Final governed artifact regeneration/review is pending.** The old transfer tables here
  and in `docs/NULL_RESULTS.md` remain explicitly historical; current configuration/target
  metadata and all unsuccessful rows must accompany any adopted replacement.
- **The historical five-specimen determinism gap is closed in source.**
  `scripts/audit_determinism.py` now declares specimens for `families.family_dataset`,
  `families.build_family`, `families.baseline_family`, `optimism.basic_bootstrap_upper` and
  `optimism.estimate_optimism`, including explicit seed-insensitive cases. Their presence
  is not a claim that the whole current audit passed; run the read-only checks rather than
  inheriting the old failure or declaring a new pass from this note.
- **The bound does not cover the fitted model.** Widening θ needs a scorer that transports
  loadings across families, which means a transfer function that does not refit per held-out
  stressor. That is a change to `analysis/transfer.py`, not to this module.
- **Power and identifiability are still single-configuration.** The families are built and the
  bound is generic; pointing them at `run_power.py` costs a scorer.
- **No family perturbs the culture context or the read schedule.** Those axes exist in
  `panel_dataset` and are orthogonal to structure, so mixing them in would confound "the
  topology was wrong" with "the plate was read early". A joint sweep is a separate job.
- **Nothing here is Tier 3.** Prospective predicted-then-measured evidence needs new
  measurements. A simulation-distribution gap bound does not supply biological validation.
