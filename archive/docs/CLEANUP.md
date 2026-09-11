# Cleanup pass — post-hardening review

> **ARCHIVED — a snapshot of 2026-08, not a description of the repository today.**
> §0's headline, "THE SUITE IS RED", was true at `2b5dcde` and has not been true since;
> run `python3 -m pytest tests/` for the current state. Two of this document's own
> resolution entries are themselves wrong — `docs/CITATION_AUDIT.md` rows 2 and 3 say
> which, and `archive/README.md` says why it is kept anyway.

Three passes, merged. The first covered ~21 fast commits on `hardening`. The second covered
the six after it. This third revision reconciles the report against `2b5dcde`, because the
owner acted on most of it in four commits while it was being written — so a document that
still read as a list of live defects would be exactly the stale prose §6 is about.

Findings are `N`-numbered; `S` numbers refer to the first revision and are summarised in
§0b. **Read §0 and §0a first: eight of the thirteen `N` findings are fixed.**

## 0. THE SUITE IS RED, and has been for four commits

**`python3 -m pytest tests/` at `2b5dcde`: 10 failed, 1,548 passed, exit code 1.**

Bisected to the commit, with the import path forced so each checkout tests its own `src`:

| commit | the four affected files |
| --- | --- |
| `50d0529` | **100 passed** |
| `7cbdf34` | **10 failed**, 90 passed |
| `abcad97` | 10 failed, 90 passed |
| `16c689b` | 10 failed, 90 passed |
| `2b5dcde` | 10 failed, 90 passed |
| `38e76ed` (HEAD as this closed) | **7 failed**, 93 passed |

**`7cbdf34` introduced it** — the commit that fixed three findings from this report also
changed `STRESSORS["H2O2"].ec50` from **0.5 to 0.15** (verified directly: 0.5 under a
`50d0529` import path, 0.15 under HEAD's) and left ten tests asserting the old panel. It has
been red through three subsequent commits.

`abcad97` — *"Tests for the parameter corrections, and the downstream effect they had"* — was
the catch-up, and it reached `tests/test_design_criterion.py` and `tests/test_stress_panel.py`
only. Four other files assert the same constants. Three of those ten have since been fixed;
seven remain at `38e76ed`:

```
tests/test_per_module_potency.py   test_a_pool_responds_below_the_regulon_it_shares_an_agent_with
                                   test_the_pool_leads_at_low_dose
tests/test_stress_model.py         test_a_module_every_training_stressor_drives_is_among_the_best_recovered
                                   test_that_module_carries_real_signal_not_just_rank
tests/test_three_sensor_build.py   test_pooling_the_same_wells_more_than_doubles_what_is_recovered
                                   test_a_stressor_off_its_axes_does_not[sorbitol]
                                   test_a_stressor_off_its_axes_does_not[caffeine]
```

These are not cosmetic. Two of them are *falsification* tests — `module_ec50("H2O2","peroxide")
< module_ec50("H2O2","oxidative")` now reads `0.175 < 0.15`, and
`test_a_stressor_off_its_axes_does_not` asserts a stressor off the build's axes recovers at
most 3 modules where it now recovers **7**. That second one is the three-sensor build's
honesty check: it exists to prove the build cannot see what it never measured, and it is
failing in the permissive direction. Whether the new EC50 is right and the test is stale, or
the recovery is now spurious, is a scientific judgement and not mine — but it cannot stay
unresolved, because a red suite is how S0 shipped last time.

**A second-order defect made this hard to find and will do so again.** `tests/conftest.py`
never puts `src` on `sys.path`, so the suite imports whatever `pip install -e` resolved —
always the main checkout. My first bisect ran five different commits in five worktrees and
got *identical* results at all five, including one I had measured green, because every run
was silently testing the working tree instead of its own checkout. Any CI job, worktree, or
`git bisect run` for this repo has the same hole: it tests the installed package, not the
revision. One `sys.path.insert` in `conftest.py` — which every script under `scripts/`
already does for itself — closes it.

## 0a. Status at `2b5dcde`

Four commits landed on this report: `7cbdf34` (three findings), `abcad97` (citation
follow-through), `16c689b` (the mu-free route into production), `2b5dcde` (partition
replication). Every row below was re-verified against the tree, not read off a commit
message.

| # | Finding | State |
| --- | --- | --- |
| **N0** | committed training tables were `--quick` presented as full | **FIXED** — `7cbdf34`. Tables regenerated at full settings and `run_training.py` now stamps `run_mode`/`n_doses`/`replicates` into all four. Verified twice: every row reads `full, 6, 6`, **and a fresh default-settings run reproduces all four byte for byte** — the condition that failed at the second revision. `NULL_RESULTS.md` rewritten, and it says plainly that an earlier version reported `--quick` output without saying so |
| **N0b** | seed spread printed, never persisted | **partly** — the tables now carry the run mode, so a reader can tell what produced them; the per-seed spread is still stdout only |
| **N1** | PMID substitution hit the comment that named the wrong identifier | **FIXED** (by me, `50d0529`) |
| **N1b** | citation scanner skipped `tests/` | **FIXED** — `abcad97`/`7cbdf34`. Both `tests/test_citations.py::BASES` and `scripts/refresh_citations.py::cited` now walk `("src", "scripts", "tests", "docs")`; the table grew 40 → **53** rows. (The helper was called `_sources` when this entry was written and is now the `BASES` constant, which also picked up `docs`.) Verified: no PMID cited anywhere is absent from it |
| **N1c** | `alkaline_ph` cites a paper that does not support its claim | **FIXED** — the cleavage claim now cites Xu & Mitchell 2001, PMID 11698381 ("Rim101p is activated by cleavage of its C-terminal region", in *S. cerevisiae*), and the pH half cites Lamb 2001, PMID 11050096. Both abstracts read and both added to the vendored table |
| **N1d** | `source` says "Takaine 2021", table says 2019 | **FIXED** — PubMed dates PMID 33654827 to 2019 Aug 5; the source now says 2019 and records that the PMID looks 2021 because Bio-protocol indexed it late |
| **N2** | two `d2_` tables predated the `not-estimable` guard; 15 wells positive off a NaN fit | **FIXED** — `7cbdf34`. Verified: **all eight** gate tables now reproduce byte for byte, and the bogus `promoter-driven` block is gone (2 remain on 20260803, both past the estimability guard) |
| **N3** | mu-free inversion had no production caller | **FIXED** — `16c689b`. `run_sensor_characterisation.py:274` uses it, falls back to the per-cell route when OD is non-positive, and keeps `per_cell_activity` beside it so the shift is measured every run rather than asserted once |
| **N3b** | scorer ignored `DoseFit.identifiable` | **FIXED, and past what I asked for** — `_split_by_identifiability`, separate `skill_biphasic_identifiable`, `n_identifiable_fits`/`n_refused_fits` columns, and a "no identifiable fit; nothing to score" refusal |
| **N4** | `--kind` typo zeroed a tracked table | **FIXED** (by me, `50d0529`) |
| **N5** | the scorer cannot say *why* a split is unscoreable | **FIXED** — both refusals now name the cause. `heldout_replicate` reports which held-out plate is absent from the readings table and why (no reporter-channel blank); `heldout_construct` reports that the fits are per construct and lists the test and fitted sets to show there is no overlap. It also confirms your "two blanks away" reading: **20260804 is missing as well as 20260728** |
| **N6** | `_score` collapsed the manifest with no `seed` filter | **FIXED** — `2b5dcde`, and closed before it could bite: `make_splits.py` now emits three seeds, and `_score` filters on `rows.seed`. The same commit found a confound I missed — seeded splits withheld two, three and four doses at different seeds, so partitions differed in training size as well as in which dose was held out — and replaced seed-sampling with scoring every interior rung |
| **N7** | `skill_score` can raise mid-report | **FIXED** — `skill_score_if_defined` returns `None` where the comparison is undefined; all six reporting call sites use it, and the strict version stays for code that must not proceed |

**And the headline result moved twice.** `16c689b` then `2b5dcde` took the interpolation win
from +0.189 → +0.305 → a median of **+0.236 over three partitions, range [−0.073, +0.301],
one of three losing outright.** The owner's own commit message calls that "the
single-partition error this exercise exists to catch, committed while catching it". My §3
N3b text claiming "+0.189, the first positive held-out skill score" is superseded; the
durable version is in `docs/NULL_RESULTS.md`.

**What did not move: the documentation.** Nine of the eleven §6 rows are still live, and two
are now wrong in a new direction because the code they describe got *fixed* — see §6.

## 0b. Resolved since the first revision

| Was | Now |
| --- | --- |
| **S0 (critical)** — `not-estimable` placed above `model-inadequate`, shadowing the 53-of-53 negative-activity finding; suite red | **Ordering fixed** in `f19ca4d`. `verdict` checks `negative_activity_fraction` first, with a comment stating why the order is load-bearing, and `tests/test_dilution_confound.py:269` pins it. Verified by re-running `run_gates.py`: both July `d2_g1passed_*` tables are byte-identical to the committed ones again. **But S0's consequence 1 was not closed** — the two NewProtocol tables were never regenerated. See N2 |
| The dropped `k_deg` term in `analysis/power._one_trial` | **Fixed** in `7a04ae6` (`signal * (growth + well.k_deg)`). No committed number moves — every construct in `generator/literature.py` still has `k_deg = 0.0`, so the change is latent by construction |
| `audit_claims.py` — 4 failures, all `check_docs_module_references` naming absent `.py` files | **Fixed** in `f19ca4d`. **37/37 pass.** `COMPARISON.md`, `research/{EXTERNAL_DATA,GENERALIZATION,TOOLING}.md` reworded, and `audit_claims.py` gained a tolerance for third-party and proposed paths |
| 11 truncated comments, 31 dead imports, `_BLOCK_WIDTH`, nested `solve_pool`, the duplicate wiring test, the `mixed` fixture | **All landed** in `f19ca4d`. Re-derived from the AST: no truncated ` ...` comment survives anywhere, and **one** dead import does — see the row below |
| The last pass's claim that `kinetic/carotenoid.py`'s *"now-unused `brentq` import was checked and is still used elsewhere in the file"* | **False.** `brentq` appears exactly once in the file, on its own import line. `solve_pool` was its only caller and the last pass deleted it. Removed here — the one thing the last pass got factually wrong about its own edits |
| §S7b — the invisible coupling between a comment's digits and a CSV row | **Broken exactly as predicted.** See N1. The prediction was right and the guard was never added |
| §4 — `analysis/uncertainty.py` and `estimator.py` no longer dead | Still true. Dead-module list re-derived from the AST graph: **the same ten modules, 1,558 LOC**, unchanged. Recommendations in §4 stand verbatim |

Still open from the last pass, unchanged and not repeated below: **S1** (finite R² written
into an unmeasurable regression), **S2** (`g1_optical` cannot tell an empty well from a
failed one), **S3** (`implied_min_k_deg` returns 0.0), **S4** (a wrong blank reports as a
plate with no analysable wells), **S5** (two `_robust_fold_change`), **S6** (three
definitions of "blank"), **S7** (citation scanning implemented twice), **S9**
(`matplotlib` is a runtime import declared as dev), **S10** (two public
`fit_od_linearity`). Every one was re-checked against the current tree and every one is
still present. **S1 is upgraded**: it is no longer a trap reasoned about a priori. On the
two NewProtocol plates it has already put 15 wells' worth of positive verdicts into a
tracked table off a NaN regression — see N2.

## 1. Verdict on the new commits

Everything below is stated against **`50d0529`**. The owner had uncommitted work in
`analysis/design.py`, `analysis/sensor_selection.py` and `generator/stress_panel.py` when
this closed — a derived `integration_hours` replacing the literal 3.0 h — which is outside
this review.

**The suite is green** — see §8 for the observed summary. Run on a machine with the plate
exports, both GSMMs and the `thermo` extras, so nothing skipped out. S0's ordering bug is
genuinely fixed and pinned by a test; its stale-artefact consequence is not (N2).

The science added is good and two of the six commits are careful work: the mu-free
inversion is correct algebra with the maturation case *refused* rather than approximated,
and `run_heldout_score.py` scores every prediction against two trivial baselines and prints
the sign. That is the right instinct twice over — and the second one has already paid off
once, because the script was rewritten mid-pass to score the repo's own biphasic form and
produced **the first positive held-out skill number in this project** (N3b).

What went wrong went wrong in the same two ways each time — **a number was recorded from a
run nobody can reproduce, and a bulk edit hit one occurrence it should not have.** Three of
the four highest findings below are a tracked artefact disagreeing with the code that
writes it.

- **N0 (critical).** All four committed `training_*.csv` tables and every number in
  `docs/NULL_RESULTS.md` §"Seed replication" come from `run_training.py --quick`, not from
  a full run. Proven: `--quick` reproduces all four byte for byte, the full run reproduces
  none of them. At full size one of that section's three conclusions is **false** and the
  script's own `!!` warning fires.
- **N1 (high).** The 14-PMID correction was applied by substitution across the whole file
  and hit the comment that existed to name the *wrong* identifier. The result asserted that
  a correctly-cited paper is a dermoscopy study, and the deliberate table row recording the
  original error was pruned in the same commit. Fixed here.
- **N2 (high).** S0's *ordering* was fixed but the tables were never regenerated. Two
  committed `d2_g1passed_*` NewProtocol tables still carry pre-guard verdicts, including
  **15 wells labelled `promoter-driven` — a positive finding — off a regression whose R² is
  NaN.** The code today correctly says `not-estimable`.
- **N3 (high).** `promoter_activity_from_total` has **no production caller.** The one place
  it applies exactly still uses the route the repo's own `FINDINGS.md` documents as 40%
  noisier.
- **N3b (high).** The rewritten `run_heldout_score.py` predicts from biphasic fits that
  `panel_calibration` marks `identifiable=False` — **7 of 8 of them** — and never reads the
  flag. The flag tracks the outcome exactly: skill **+0.189** where the fits are
  identifiable, **−1.913** where none are.
- **N4 (medium).** `run_heldout_score.py --kind <typo>` wrote a zero-byte
  `outputs/heldout_scores.csv` over the tracked one and exited 0. Fixed here.

Documentation drift is now the largest standing category and it has grown, not shrunk:
`docs/research/PARAMETER_SOURCES.md` asserts as of today that a defect fixed in `7adae00`
is "*not* yet fixed", and `docs/REPRODUCING.md` §4 states the exact opposite of what N0
found.

## 2. What I fixed this pass

All seven landed in `50d0529`, committed by the owner rather than by me.

**`src/ystwin/generator/stress_panel.py:55-63` — a comment asserting a correct citation is
a dermoscopy paper.** See N1 for the mechanism. Rewritten to seven true lines naming both
real identifiers (Takaine `30858198`, Yaginuma `25283467`) and saying where the record of
the retired one now lives. Both identifiers were already in the vendored table, so the
citation tests are unaffected.

**`tests/test_atp_sensor.py:245` — a live wrong PMID.** The class docstring still credited
`24815987` (Vasefi 2014, dermoscopy) with Yaginuma's Kd of 4.5 mM. Corrected to `25283467`.
The commit titled "all 33 now resolve on-topic" missed it because
`tests/test_citations.py::_sources()` walks only `src/` and `scripts/` — see N1b.

**Comment concision.** `reporter.promoter_activity_from_total`'s docstring was 40 lines
above a 12-line body, and the same derivation appears verbatim in three places
(`docs/FINDINGS.md`, the function, and `tests/test_reporter.py::TestTotalSignalForm`).
Cut to the algebra, the one-line reason it matters, and a pointer to the measurement.
`Raises:` and the maturation rationale kept intact — that is the load-bearing part.

**`scripts/run_heldout_score.py` — `--kind` is now validated.** An unknown value scored
nothing, wrote an empty frame over `outputs/heldout_scores.csv`, and exited 0. It now
prints the available kinds and returns 2, matching the missing-input refusal ten lines
above it. This is a control-flow change, not a mechanical one; called out here for that
reason. It protects a tracked artefact from a typo and moves no number. Also dropped an
`f` prefix on a string with no placeholders.

> This file was rewritten by its owner partway through the pass — two models instead of
> one, see N3b. The guard survived the rewrite unchanged, and I re-reviewed and re-ran the
> new version: it reproduces the on-disk `outputs/heldout_scores.csv` exactly. Everything
> below about this script describes the **new** version.

**`scripts/run_training.py` — two dead bindings.** `built_replicated` and
`minimal_replicated` were assigned and never read. Calls kept, bindings dropped, with one
line saying the spread is printed and not written. See N0b.

**`tests/test_reporter.py`** — the new import was inserted above `ReporterKinetics`,
breaking the file's existing order. Restored.

**`src/ystwin/kinetic/carotenoid.py` — one dead import the last pass left behind.**
`brentq` was the nested `solve_pool`'s only consumer; deleting the function orphaned it and
the last pass recorded the opposite. Removed. A full AST re-derivation over `src`, `scripts`
and `tests` finds no other unused import (every remaining hit is `from __future__ import
annotations`, which a name-occurrence scan cannot see).

## 3. Science bugs found and NOT fixed

Ranked. `N` numbers are new this pass; `S` numbers from the last pass are listed in §0 and
not restated.

### N0 (critical) — FIXED in `7cbdf34` — every committed training table and the whole seed-replication result is `--quick` output presented as a full run

`875807d` added `DATASET_SEEDS = (0, 1, 2)` and `replicated_transfer`, and recorded the
result in `docs/NULL_RESULTS.md` §"Seed replication: one build comparison does not survive
it". Those numbers, and all four tracked `training_*.csv` tables, are `--quick` output.

`--quick` is 4 doses × 3 replicates. The default is 6 × 6.

Proven by running both into `YSTWIN_OUTPUTS` and diffing, twice:

| | `--quick` | full run | committed |
| --- | --- | --- | --- |
| `training_transfer_wide.csv` | — | — | **byte-identical to `--quick`** |
| `training_transfer_build.csv` | — | — | **byte-identical to `--quick`** |
| `training_transfer_three_sensor.csv` | — | — | **byte-identical to `--quick`** |
| `training_recovery.csv` | — | — | **byte-identical to `--quick`** |

Two consecutive full runs agree byte for byte with each other, so this is not
nondeterminism and not an environment drift — installed `numpy 2.4.1 / scipy 1.17.0 /
pandas 2.3.3` match the `==` pins exactly. It is the flag.

What the full run says instead:

| | `NULL_RESULTS.md` (= `--quick`) | full run |
| --- | --- | --- |
| five-channel latent width | 5, 5, 5 | **5, 3, 5** |
| five-channel modules recovered | 5, 6, 5 → 5 [5, 6] | 4, 4, 5 → 4 [4, 5] |
| five-channel score on driven modules | 0.216, 0.268, 0.248 → 0.248 | 0.223, 0.222, 0.223 → 0.223 |
| three-sensor latent width | 4, 4, 4 | 4, 4, 4 |
| three-sensor modules recovered | 5, 6, 6 → 6 [5, 6] | 2, 1, 4 → **2 [1, 4]** |
| three-sensor score | 0.169, 0.185, 0.206 → 0.185 | 0.103, 0.091, 0.13 → **0.103** |

Conclusion by conclusion:

1. *"Latent width was stable at both builds… the width selection is not seed-sensitive"* —
   **false at full size.** The five-channel build picks 5, 3, 5, and
   `replicated_transfer`'s own `!! latent width is not stable across seeds` line fires. The
   document reports the absence of a warning the code emits.
2. *"Module count does not separate the builds. Both span 5 to 6 across seeds."* — **false
   at full size.** 4–5 against 1–4. They separate, in the direction the section says they
   do not.
3. *"The score does separate them… by roughly a third."* — **survives, and strengthens.**
   0.222–0.223 against 0.091–0.130 is a factor of two, not a third. The claim is right; the
   numbers quoted for it are not.

So the one section written to demonstrate that a single seed cannot be trusted rests on a
single undisclosed configuration, and the finding it reports as surviving replication is
the one it understates.

Two mechanisms let this through and both are cheap to close:

- **`run_training.py` prints nothing about which mode it ran in**, and writes nothing into
  the CSVs recording it. A `--quick` table and a full table are indistinguishable after the
  fact. One `print` at the top of `main`, or a `quick` column, ends this class of error.
- **`docs/REPRODUCING.md` §4 asserts the opposite**: *"the four `training_*` tables are
  unchecked… the second set was written by a full `run_training.py` rather than the
  `--quick` run used here, so a diff against `--quick` output would mean nothing."* A diff
  against `--quick` output is in fact exact, and it is the check that was needed.

Not fixed: regenerating the tables changes every number in four tracked artefacts and
rewrites a `NULL_RESULTS.md` section, which is squarely out of scope. This is the single
highest-value change available in this repo right now.

### N0b (medium) — partly fixed in `7cbdf34` — the seed-replication result exists only in terminal scrollback

`replicated_transfer` returns a combined DataFrame with a `seed` column. `main` assigned it
to `built_replicated` / `minimal_replicated` and never used it, so nothing reaches
`outputs/`. Every other claim in this repo is backed by a tracked table — `REPRODUCING.md`
opens with *"Every number quoted anywhere in this repo comes from a script in `scripts/`
writing a table into `outputs/`"* — and this one is not. I removed the dead bindings; the
gap is the missing `to_csv`, which would add an output artefact and so is not mine.

Related: seed 0 is trained twice per build, once by `train_and_score` and once inside
`replicated_transfer`. Harmless, and it roughly doubles the cost of the slowest script in
the repo.

### N1 (high) — FIXED in `50d0529` — the PMID correction was a whole-file substitution, and it hit the one comment where the old digits were load-bearing

`7adae00` corrected 14 identifiers. **All 13 replacements are right.** Verified against
E-utilities `esummary` — title, journal and year of every new row in
`data/citations/pmid_titles.csv` match the record exactly, and every author-year attribution
in the `source` strings is correct:

`10594829` Jung & Levin 1999 · `10604478` Beck & Hall 1999 · `12509465` Lamb & Mitchell 2003
· `12676948` Young 2003 · `21982714` Hung 2011 · `3043194` Thiele 1988 · `32130885` Pak 2020
· `8422683` Liao & Butow 1993 · `9271382` Zhao & Eide 1997 · `9407035` Stathopoulos & Cyert
1997 · `9409150` Thomas & Surdin-Kerjan 1997 · `9510529` Kwast 1998 · `9741624` Huang 1998.

The vendored table also matches the source exactly: 40 rows, 40 PMIDs cited in
`src/` + `scripts/`, no cited-but-absent, no absent-but-cited, and every `cited_in` cell
correct. Dropping `12588995` (War1p, weak-acid stress) was right — it was cited for
Rim101 and alkaline pH.

The defect is the 14th. The edit replaced `24815987` → `25283467` everywhere in
`stress_panel.py`, including inside this comment:

```
# QUEEN-2m's yeast demonstration is Takaine, not Yaginuma 2014 -- that one is real and
# correctly transcribed, but it is an E. coli paper, reached by PMID 25283467, which
# belongs to a dermoscopy study.
```

`25283467` **is** Yaginuma 2014 (*Sci Rep* 4:6522). The dermoscopy paper is `24815987`
(Vasefi 2014). So the comment named the correct citation and called it a dermoscopy study,
one sentence after calling it real and correctly transcribed. Its second paragraph —
*"That identifier stays in this comment and never inside a `source` string. Naming it here
keeps its row in `data/citations/pmid_titles.csv` non-stale"* — became false twice over:
the identifier it now names **is** in two `source` strings, and the row it claimed to keep
alive was deleted from the table in the same commit.

This is §S7b of the last pass, verbatim: *"Shortening that comment without preserving the
digits fails the suite with a message about re-running `refresh_citations.py` — which would
prune the row, deleting the record of the error the row exists to preserve."* That is what
happened, except the substitution kept the suite green, so nothing objected.

Fixed here as a comment correction. **Not** fixed: the row itself. Restoring `24815987` to
the table needs the digits back in `src/` or `scripts/`, and the only honest place for them
is a comment — which is the coupling that just failed. The record now lives in
`tests/test_citations.py:132`, `docs/SOURCE_AUDIT.md` §3.2 and
`docs/research/PARAMETER_SOURCES.md`. If that is acceptable, say so in
`test_citations.py`'s docstring; if not, the guard §S7b asked for — a test asserting the
digits *are* present in the module text — is the thing to add.

### N1b (high) — FIXED in `abcad97` — the citation guarantee stops at the `tests/` boundary, and a wrong PMID was sitting on the other side of it

`tests/test_citations.py::_sources()` globs `src/**/*.py` and `scripts/**/*.py`. `tests/`
is not scanned, by either the test or `scripts/refresh_citations.py`. Five PMIDs are cited
only from `tests/` and are outside every check in this repo:

| PMID | Cited in | Status |
| --- | --- | --- |
| `24815987` | `tests/test_atp_sensor.py:245` | **Wrong.** Vasefi 2014, dermoscopy, credited with Yaginuma's Kd. Fixed here |
| `9427394` | `tests/test_atp_sensor.py:30` | Kratzer & Schuller 1997 — unverified |
| `21335394` | `tests/test_thermodynamic_bridge.py` | unverified |
| `2566299`, `27317316` | `tests/test_physiology_validation.py` | unverified |

The one that was checkable was wrong, which is the argument for scanning `tests/` too. Not
fixed: adding `tests/` to `_sources()` requires four more table rows, which is a data
change. It is a two-line change plus one `refresh_citations.py` run, and it closes a hole
that has already produced one error.

### N1c (medium) — OPEN — one of the 14 corrections did half of what its own audit asked

`docs/SOURCE_AUDIT.md` §3.1 flags `alkaline_ph` twice. The identifier row says
`12588995` → **`12509465`**, and a `⚠` note directly under the table says:

> `alkaline_ph` is doubly weak: even the *intended* Lamb & Mitchell 2003 (MCB 23:677) is
> about Rim101 repressing `NRG1`/`SMP1`, not about its proteolytic activation at alkaline
> pH. The proteolysis result is Li et al. 2004 / Xu & Mitchell 2001. **Fix the identifier
> and the claim.**

`7adae00` fixed the identifier. `stress_panel.py:240` still reads *"Lamb & Mitchell 2003,
PMID 12509465 -- Rim101 is cleaved and activated at alkaline pH"*, and E-utilities confirms
`12509465` is *"The transcription factor Rim101p governs ion tolerance and cell
differentiation by direct repression of the regulatory genes NRG1 and SMP1"*. So the module
now cites a real, on-topic paper for a mechanism that paper does not report — which is
exactly the residue `test_citations.py`'s own docstring admits it cannot catch: *"Whether a
paper supports the sentence next to it is a judgement, and no test makes it."*

Not fixed: changing either the claim or the citation is a scientific edit, and the audit
already names the two candidate papers. This is the cheapest outstanding item in
`SOURCE_AUDIT.md` §3.1 and it is a two-field edit.

### N1d (low) — OPEN — the vendored table and a `source` string disagree on a year

`stress_panel.py:433` cites *"Takaine 2021, PMID 33654827"*. The table row for `33654827`
says **2019**, and E-utilities agrees. `docs/research/PARAMETER_SOURCES.md` already
flagged this. The refresh pins digits, not years, so nothing catches it.

### N2 (high) — FIXED in `7cbdf34` — two committed `d2_` tables predate the `not-estimable` guard, and 15 wells in them still report a positive finding from a regression that does not exist

The last pass listed "the committed tables no longer match what the code produces" as
consequence 1 of S0. The ordering fix closed it for the July plates — both
`d2_g1passed_2026070*` and both `d2_*__mCitrine.csv` now reproduce byte for byte — and left
it open for the other two. Checked by running `scripts/run_gates.py` into `YSTWIN_OUTPUTS`:

| Table | Result |
| --- | --- |
| `g1_*` × 4 | **identical** |
| `d2_g1passed_20260701_ER_preliminary_(RAW).csv` | **identical** — the S0 fix restored these |
| `d2_g1passed_20260709_ER_stress_1st_(RAW).csv` | **identical** |
| `d2_g1passed_20260722_…_NewProtocol_ANALYSED.csv` | **47 of 67 verdicts differ** |
| `d2_g1passed_20260803_…_Replicate3.csv` | **35 of 52 verdicts differ** |

Only the `verdict` column moves; every numeric column is byte-identical. The direction:

| | committed | code today |
| --- | --- | --- |
| 20260722 | `mixed` × 44, `promoter-driven` × 3 | `not-estimable` × 47 |
| 20260803 | `mixed` × 23, `promoter-driven` × 12 | `not-estimable` × 35 |

**The new verdict is the right one.** `dilution_r2` is **NaN in every one of those 82
rows** (`qss_fraction` runs 0.00–0.16, and the `usable.sum() < _MIN_REGRESSION_POINTS`
branch correctly writes `nan, nan`). So the guard is firing exactly where it was designed
to, and the committed tables are the pre-guard state.

That makes the interesting number the 15 `promoter-driven` rows. `promoter-driven` is a
*positive* label — "the induction is real and not a dilution artefact" — and those 15 wells
were carrying it off a regression whose R² is NaN. This is S1 from the last pass, no longer
hypothetical and no longer confined to the `np.std(x) < 1e-6` branch I could only reason
about: on the two NewProtocol plates the shortfall is *points*, not spread, and 15 wells
reported a finding from a fit that was never run. The `not-estimable` guard closes it, and
the tracked artefacts still show the open version.

Two things follow:

1. **Regenerate the two `d2_g1passed_*` NewProtocol tables.** One `run_gates.py` run. It
   changes a verdict column in a tracked artefact, so it is not mine, but unlike N0 there
   is no ambiguity about which output is correct.
2. `docs/REPRODUCING.md:150` claims *"all six `g1_`/`d2_` tables"* reproduce byte for byte.
   There are **eight**, and two do not. Add to §6.

### N3 (high) — FIXED in `16c689b` — the mu-free inversion has no production caller, and the one place it applies exactly still uses the noisier route

The algebra is right and the tests are the right tests. `k_synth = (dF/dt)/X + k_deg (F/X)`
follows from `R = F/X` by substitution, mu cancels for any `k_deg`, and
`test_agrees_with_the_per_cell_route_when_growth_is_known_exactly` pins the two routes
together at `rtol=0.02` — which is the check that matters, because it proves the new route
is a reformulation and not a different estimator. Refusing a maturing reporter instead of
returning the mature-only answer is the correct call and is the best-behaved refusal added
this pass.

But **nothing outside `tests/test_reporter.py` calls it.** Grep of the whole tree: five
production sites compute activity, and all five use `promoter_activity`.

| Caller | Kinetics | Has the total and the biomass in hand? | Route used |
| --- | --- | --- | --- |
| `scripts/run_sensor_characterisation.py:262` | `k_deg=0.0`, no maturation | **Yes** — `normalised = (fl - rb) / (od - ob)` at `:124`, both channels present | per-cell |
| `src/ystwin/analysis/uncertainty.py:445, :475` | passed in | yes, via the same caller | per-cell |
| `src/ystwin/diagnostics/dilution_confound.py:154` | `k_deg` from the report | yes | per-cell |
| `src/ystwin/generator/calibrate.py:155` | `k_deg` argument | yes | per-cell |
| `scripts/run_scenarios.py:62` | `k_deg=0.0` | simulated, mu exact — no benefit | per-cell |

`run_sensor_characterisation.py` is the one that matters. Its `KINETICS =
ReporterKinetics(k_deg=0.0)` means no maturation and the second term vanishes, so the total
form applies with no caveat at all — and its output, `outputs/sensor_characterisation.csv`,
is the table behind `DEFAULT_WELL_CV`, `MEASURED_ACTIVITY_CV`, every dose-response fold,
and now `run_heldout_score.py`. Every one of those is computed through the route
`docs/FINDINGS.md` reports as 40% rougher on the same wells.

`scripts/run_calibration_nis.py` does **not** compute activity from a trace despite
appearances; its only hit is a prior tuple at `:123`.

Nothing "silently prefers" one — the split is unambiguous and one-sided. `FINDINGS.md`
says *"`promoter_activity` is kept, not replaced"*, which is right for the maturation case
but reads as if the choice is made per call site. It is not made at all.

Two things to decide, neither of them mine:

1. Whether `run_sensor_characterisation.py:262` should switch. It would move
   `activity_late` in a tracked table and therefore `DEFAULT_WELL_CV` and every fold.
2. Where the switch belongs. Threading a `route=` flag through `activity_uncertainty` and
   `dilution_confound_report` is the wrong shape; if the total form is the better estimator
   whenever there is no maturation, `promoter_activity` should dispatch to it when handed
   the total, and the per-cell entry point should be the special case.

**A caveat on the 40% claim.** `FINDINGS.md` checks that the negative-activity fraction is
unchanged (11% on both routes) and correctly says why that check matters — a reformulation
that smoothed the refutation away would be hiding the finding. But that check was run on
the **NewProtocol** wells, where the fraction is 11%. The finding it protects lives on the
**July** plates, where it is 57–100% (§S0 of the last pass), and those plates were not
re-run. The right check is on the plate that carries the claim.

**Reproducibility.** The FINDINGS table (1637.3 / 1630.6 / 594.9 / 358.6 / 11%) is produced
by no script in `scripts/`. Same class as N0b: a number in a document with no reproducer.

### N3b (high) — FIXED in `7cbdf34`, result since superseded — the held-out scorer predicts from fits its own module marks unidentifiable, and that flag predicts exactly when the prediction fails

**`scripts/run_heldout_score.py` was rewritten on disk during this pass**, after my first
review of it. It now scores the repo's own biphasic form —
`(basal + amplitude·d/(ec50+d)) / (1 + (d/lethal)^hill)` from
`generator/panel_calibration.fit_dose_response` — with the log-linear fit demoted to a
third baseline. I re-reviewed and re-ran the new version; it reproduces the on-disk
`outputs/heldout_scores.csv` exactly. My `--kind` guard (N4) survived the rewrite.

The rewrite is the right move — the calibration module exists precisely because a plain
saturating Hill cannot describe the fall at high dose, so the biphasic form is the fair
model to score. And the result is the most interesting number in the repo:

| split | biphasic | log-linear | train-mean | nearest-dose | skill biphasic | skill log-linear |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| extrapolation | 941.1 | 504.2 | 807.4 | **323.1** | **−1.913** | −0.561 |
| interpolation | **294.5** | 374.0 | 435.3 | 363.3 | **+0.189** | −0.029 |

**That `+0.189` was the first positive held-out skill score anywhere in this project.** On
interpolation the biphasic form beat both trivial baselines; on extrapolation it was nearly
twice as bad as carrying the nearest dose forward.

> **Superseded.** `16c689b` moved the identifiable-only figure to +0.305, then `2b5dcde`
> replicated over three partitions and got a median **+0.236, range [−0.073, +0.301], with
> one partition losing outright.** A real signal, not an established win. The numbers in
> this section are the ones that were current when the defect below was found; the live
> ones are in `docs/NULL_RESULTS.md`. The defect itself was fixed — see §0a.

The defect is what separates those two rows, and the script has the answer in hand and
throws it away. `fit_dose_response` returns a `DoseFit` carrying `identifiable: bool` and
`unidentified_reason: str`, and `panel_calibration` exports `induction_is_identifiable` as
public API. Refitting the eight training ladders:

| split | construct | `identifiable` | `unidentified_reason` |
| --- | --- | --- | --- |
| extrapolation | AlteredYap1 | **False** | *"the fitted lethal dose 1.12 is below the inducing EC50"* |
| extrapolation | NativeYap1 | **False** | lethal 1.12 below EC50 1.58 |
| extrapolation | UPRE1 | **False** | lethal 1.83 below EC50 9.01 |
| extrapolation | UPRE2 | **False** | lethal 1.87 below EC50 7.30 |
| interpolation | AlteredYap1 | **False** | lethal 1.09 below EC50 1.11 |
| interpolation | NativeYap1 | **False** | lethal 1.07 below EC50 1.40 |
| interpolation | UPRE1 | **True** | — |
| interpolation | UPRE2 | **True** | — |

Seven of eight fits are degenerate in a specific, named way: the curve says the culture
dies before it induces. `_fit_biphasic` catches `ValueError` and `RuntimeError` and then
predicts from `fit.amplitude`, `fit.ec50` and `fit.lethal_dose` without ever reading
`fit.identifiable`. So on `extrapolation` all four predictions come from unidentifiable
fits and the model scores −1.913; on `interpolation` the two identifiable constructs are
present and it scores +0.189.

**And the commit message already contains the right explanation, derived by hand.** It says
extrapolation loses because of *"a lethal-dose term in the denominator estimated from data
that never observed lethality."* That is precisely what `unidentified_reason` says, in a
string, on the object: *"the fitted lethal dose 1.12 is below the inducing EC50."* The
reasoning is correct and the machine-readable form of it was already in hand. So this is not
a wrong conclusion — it is a conclusion argued in prose that the code could have printed,
and that a reader of `heldout_scores.csv` cannot see at all.

This is §5's argument in one place: the repo built a typed refusal with a reason string,
made it public, and the newest consumer reads neither. Three ways to use it, in increasing
order of ambition — report `n_identifiable` alongside `n_test`; score identifiable and
unidentifiable constructs separately; or refuse the biphasic row when no training ladder is
identifiable. Any of the three moves the explanation out of the commit message and into the
table. The middle one is also the experiment §7 asks for.

Not fixed: all three change what the table reports.

Two smaller notes on the same rewrite:

- `from ystwin.generator.panel_calibration import _model as _biphasic_model` — a **private**
  module-level function imported across a package boundary. `_model` is deliberately absent
  from that module's `__all__`, and `DoseFit` has no `predict`. The missing public method is
  the reason the script had to reach for the underscore; `DoseFit.predict(dose)` is the
  two-line fix and it belongs in `panel_calibration`, not here.
- `rmse` changed from `np.mean` to **`np.nanmean`**, so a NaN prediction is silently dropped
  from one model's error and not the other's, while `n_test` reports the loop count for
  both. Verified latent, not live: zero non-finite biphasic predictions on either split
  today. It is still an RMSE over an unstated denominator.

### N4 (medium) — FIXED in `50d0529` — `run_heldout_score.py` destroyed its own tracked output on a typo

`main` built `kinds = [args.kind] if args.kind else …` with no membership check, so
`--kind nonsense` scored nothing, reached `pd.DataFrame([]).to_csv(...)`, wrote a **zero-byte**
`outputs/heldout_scores.csv` over the tracked table, and returned 0. Reproduced, then
fixed, then restored from git.

This is precisely the failure `docs/REPRODUCING.md` §5 says the repo does not have: *"It
does not fall back, and it does not write an empty table — an empty gate report reads like
a plate that passed nothing rather than like a plate that was never there."* And
`run_gates.py:50` carries a comment making the same argument. The script inherited the
refusal *shape* from `run_gates` (§5) and not the rule about writing.

Worth a wider look than the one site: any script that writes a tracked CSV
unconditionally at the end of `main` has the same shape.

### N5 (medium) — OPEN — two of four held-out splits cannot be scored, and the script cannot say why

`outputs/heldout_scores.csv` reports two results and two refusals. Both refusals are
structural, not incidental, and the reasons are real:

**`heldout_construct`** withholds whole constructs — test `{AlteredYap1, UPRE2}`, train
`{NativeYap1, UPRE1}` — and both `_fit_linear` and `_fit_biphasic` fit **per construct**, so
`row.construct not in fits` is true for every test row, always. The model class cannot
address this split at all. The script says *"no test row had a fitted construct"*, which
reads as a data gap.

**`heldout_replicate`** holds out plate `20260728`, which is **absent from
`outputs/sensor_characterisation.csv`**: it is one of the two plates with no
reporter-channel blank, so the characterisation pipeline never produces rows for it. Zero
test rows, permanently. The script says *"not scoreable: 42 train, 0 test rows"*.

The join is doing the right thing in both cases — this is a genuine mismatch between the
manifest (four plates: `20260722`, `20260728`, `20260803`, `20260804`) and the readings
table (two: `20260722`, `20260803`). `docs/NULL_RESULTS.md` explains both correctly and at
the right length. The gap is that the *script* cannot: it reports a row count where it
could report "the held-out plate is not in the characterisation table", which it has the
data to say.

One consequence worth stating even so. `docs/SPLITS.md` requires the training side of
`heldout_replicate` to keep at least `uncertainty.MIN_PLATES_FOR_INTERVAL` = 3 plates, and
concludes *"Four plates is therefore the minimum."* The manifest satisfies that — train
`{20260722, 20260803, 20260804}`. The **scorer** does not: only two of those three plates
exist in `sensor_characterisation.csv`, so even with the held-out plate present the training
side would carry two plates, below this package's own threshold for stating a
batch-generalisation interval at all. The replicate claim is two blanks away from being
scoreable, not one.

### N6 (low) — FIXED in `2b5dcde` — the plate-key join holds, but for a reason nothing enforces

`_plate_key` takes the first eight digits of whatever it is handed. Checked against every
value both sides actually produce:

| Side | Producer | Value | `_plate_key` |
| --- | --- | --- | --- |
| manifest | `make_splits.py:100`, `path.name[:8]` via `RECORDED_PLATES` | `20260722.0` (pandas floats the column because `simulated_panel` rows have no plate) | `20260722` ✓ |
| manifest | same | `20260728.0`, `20260803.0`, `20260804.0` | ✓ all four |
| readings | `run_sensor_characterisation.py:251`, `path.stem[:22]` | `20260722_ER&OxidativeS`, `20260803_ER&oxidatives` | ✓ both |

So it holds for all four. It holds because both sides happen to start with the date, and
only one side enforces that: `make_splits` takes `path.name[:8]` and looks it up in
`RECORDED_PLATES`, so a non-date-leading filename drops out. `path.stem[:22]` enforces
nothing. A stem like `Plate1_20260728_…` yields `12026072` on the readings side and
`20260728` on the manifest side, and the docstring's own failure mode — *"a naive join
silently produces no rows, which reads as 'no overlap' rather than 'wrong key'"* — returns,
now one layer deeper. The authoritative signal is on the manifest side already:
`RECORDED_PLATES` is keyed by the eight-digit date, and resolving the readings plate
through it rather than through digit-slicing would make the invariant enforced instead of
observed.

Two smaller notes on the same function:

- `_plate_key` on a missing value gives `""`. `main` filters to `real_biosensor` before
  mapping, so no NaN reaches it today; nothing states that dependency.
- `keyed = {(plate, construct, dose): assignment}` collapses the manifest with **no `seed`
  filter**, and the manifest has a `seed` column. One seed exists per `split_kind` today
  and no key collides (checked: 0 conflicting assignments over 238 keys), so it is correct
  now. The day `make_splits.py` sweeps seeds, the last row silently wins and
  `rows.partition_hash.iloc[0]` labels the result with a partition it did not use.

### N7 (low) — OPEN, now four call sites — `skill_score` can raise from inside the reporting loop

`skill_score` raises `ValueError` when `baseline == perfect`, which is the right refusal.
`_score` calls it with `perfect=0.0` and `best_baseline = min(rmse)`, so a baseline that
predicts a test set exactly — one test row whose nearest training dose has the same
activity — crashes the script mid-report. Latent, not live: both current baselines are
hundreds of RFU/OD/h from zero.

## 4. Dead code

Re-derived from the AST import graph with relative imports resolved. **Unchanged from the
last pass: the same ten modules, 1,558 LOC, each imported by nothing but its own test.**
The table and recommendations in the previous revision stand; nothing was wired and nothing
was added. `PARKED.md` still covers three of the ten, and the same five sit outside it:
`generator/calibrate.py`, `calib/kdeg.py`, `calib/od_linearity.py`,
`analysis/deconvolve.py`, `bridge/tmfa.py`.

**`analysis/latent.py` is not dead and is not removable.** Asked to re-derive it: four
non-test importers — `analysis/recovery.py:18`, `analysis/transfer.py:20`,
`analysis/stress_model.py:22`, `scripts/run_transfer.py:24` — plus
`scripts/audit_determinism.py`, which pins `fit_latent` and `select_dimension` by
signature. It is the centre of the training path, not a leaf. The note that it is
"replaceable by NNLS" is about the algorithm, not about deadness, and acting on it would
move every number in `training_*.csv`, `transfer_by_stressor.csv` and
`module_recovery.csv`. Its `_factorise` is an iterated `np.linalg.svd` imputation, so the
recovered subspace is LAPACK-dependent where singular values are close — worth knowing
before anyone attributes a moved digit to the swap. (It is *not* the cause of N0: two full
runs here agree byte for byte.)

**New at function scope:** `reporter.promoter_activity_from_total` is exported, documented,
tested — and called by nothing outside its test. See N3.

## 5. Refusal convention

The catalogue from the last pass stands. Against its recommendation — raise for caller
errors, return a flagged result for per-unit refusals, NaN only for genuinely undefined
scalars, retire `None` for values — the new code scores as follows:

| New refusal | Convention | Verdict |
| --- | --- | --- |
| `promoter_activity_from_total` raises on maturation | Tier 1, caller error | **Correct.** "You called the wrong function", with the reason and the alternative in the message. The best refusal added this pass |
| `promoter_activity_from_total` raises on non-positive biomass | Tier 1 | **Correct**, and consistent with `naive_specific_fluorescence` |
| `DilutionConfoundReport.verdict → "not-estimable"` | Tier 2, but as a magic string | **Half.** It survives into a DataFrame column, which is the property that matters, but it carries no reason and callers must string-compare. `uncertainty.FoldChange.estimable` + `method` is the model this should have followed, and would have made S0 hard to write |
| `DoseFit.identifiable` + `unidentified_reason` (existing tier-2 refusal) | Tier 2, and the best-designed one in the package | **Ignored by its newest consumer.** `run_heldout_score._fit_biphasic` predicts from `fit.amplitude`/`ec50`/`lethal_dose` without reading either field, on fits that are `identifiable=False` 7 times out of 8, and the flag turns out to separate a +0.189 skill from a −1.913 one (N3b). A flagged result nobody checks is a NaN with extra steps — which is the whole argument for tier 2, running in reverse |
| `run_heldout_score._score → {"kind": …, "note": …}` | Tier 2, script-level | **Consistent, and deliberately so.** It is the same dict-with-a-`note` that `run_gates.run()` uses at `:50-52`, under the comment *"Never return None here. An absent plate and a plate that passed nothing look identical in a report that says nothing about either."* The new script inherited the right pattern |
| `run_heldout_score._score → None` when a kind has no rows | Tier: pure lookup | **The inconsistency, and it is inside one function.** The same function refuses one way with a `note` dict and another way with `None`, ten lines apart, and `main` dropped the `None`s with a truthiness filter. The only way to reach it was a `--kind` that is not in the manifest — a caller error, tier 1 — and that is the path that produced N4. My fix puts the tier-1 refusal in `main` where it belongs, which leaves this branch unreachable rather than corrected; the honest follow-up is to delete it |

So the *raising* tier now has a coherent convention and the new code follows it exactly.
The flagged-result tier is where it still breaks, and this pass found a sharper diagnosis
than last pass's "five conventions with no rule": the convention is fine, **the consumers
do not read it.** `DoseFit` has exactly the `estimable`/`reason` pair last pass recommended
building everywhere, and the newest caller ignores it (N3b). Adding the flag to
`DilutionConfoundReport` and `WellQualityResult` is still the right move, and it is now
clearly not sufficient on its own — whatever is added has to be *checked*, and the cheapest
way to make that stick is a test per consumer asserting the refusal path is honoured.

Secondary, and still worth doing: `"not-estimable"` is a magic string where a flag and a
reason belong, and `_score` uses two refusal conventions ten lines apart.

## 6. Stale prose

**This is now the largest standing category, and the one nothing acted on.** Eleven rows
were listed at the second revision; **nine are still live at `2b5dcde`**, and two of those
are wrong in a *new* direction because the code they describe was fixed underneath them.
The code moved in four commits over one morning; not one document number moved with it.

`scripts/audit_claims.py` passes 37/37 and none of this is visible to it — it checks module
references, output-table existence, documented commands, dependency pinning, home paths and
the README's test count, and no other number in any document. Everything here was found by
hand and re-verified against the artefact at `2b5dcde`.

**Re-verified at `2b5dcde`:**

| Where | Says | Reality |
| --- | --- | --- |
| ~~`docs/research/PARAMETER_SOURCES.md:240`~~ | ~~"…and **not** yet fixed: the 14 broken module PMIDs are still present"~~ | **RESOLVED** — `abcad97`. The assertion is gone from the file |
| `docs/REPRODUCING.md:150-152` | the four `training_*` tables *"are unchecked… written by a full `run_training.py` rather than the `--quick` run used here, so a diff against `--quick` output would mean nothing"* | **Wrong twice over now.** It was wrong at the second revision because they were `--quick`; it is wrong today because `7cbdf34` regenerated them at full settings *and* stamped `run_mode` into every row, so they are neither unchecked nor ambiguous. The sentence describes a state that has not existed in either direction |
| `docs/REPRODUCING.md:139` | *"…reproduces … all six `g1_`/`d2_` tables byte for byte"* | There are **eight**. At the second revision two did not reproduce; after `7cbdf34` **all eight do** — verified this pass. So the claim is now true of more tables than it names, and the count is still wrong |
| `docs/REPRODUCING.md` §4 | the script/output table, *"built by reading the `to_csv` / `save` calls"* | Missing `run_heldout_score.py` → `heldout_scores.csv`, a new script writing a **tracked** table. The audit checks docs→module existence, never module→docs coverage |
| `docs/FINDINGS.md:1283` | the negative-activity fraction is unchanged between activity routes, 11% either way | True on the NewProtocol wells; the finding it protects is on the July plates at 57–100%, which were not re-run — see N3 |
| ~~`docs/NULL_RESULTS.md` §"The first held-out prediction … and it loses"~~ | ~~one "model RMSE" column, log-linear only~~ | **Resolved during the pass.** `50d0529` rewrote the section as "wins inside the range, loses outside it", with both forms, both skills, and an explicit note on what the earlier version got wrong. Accurate as it stands |
| `docs/CLAIM_BOUNDARY.md:18` | *"One forward prediction has now been scored on held-out real data — **and it lost**."* | **Still stale after three corrections passed it by.** The interpolation result went +0.189 → +0.305 → a median **+0.236 over three partitions**, and `NULL_RESULTS.md` was rewritten twice while this paragraph was not touched once. The Tier-3 argument after it — scored after the fact rather than registered in advance — is untouched by any of it and still holds; it is only the "and it lost" that is wrong |

**Carried over, all re-verified as still wrong:**

| Where | Says | Committed artefact says |
| --- | --- | --- |
| `docs/CLAIM_BOUNDARY.md:75` | attribution power "0% / 4% / 46% / 98%" | `power_analysis.csv`: **0.8 / 5.0 / 44.2 / 95.8** → 1% / 5% / 44% / 96% |
| `docs/CLAIM_BOUNDARY.md:76` | "Collinearity 0.998 dose-only vs 0.570 crossed" | **0.9932 / 0.5874**. `FINDINGS.md:140` says outright that `0.998 / 0.570` came "from a parameter set no longer in the repository" — so this line now contradicts another document in the same directory, not just the CSV |
| `docs/FINDINGS.md:97` **and** `:113` | dose × 2 nutrients collinearity **0.694**, in prose and again in the table | **0.5874**. Two occurrences, not one, and `:140` twelve lines below corrects the *other* collinearity pair while leaving this one |
| `docs/FINDINGS.md:147` | "At 1.5x the current design has **4%** power" | its own table at `:112` says **5%**, matching the CSV |
| `docs/SOURCE_AUDIT.md:43, 145-171, 452` | §3.3 "HIGH", arguing `DEFAULT_WELL_CV = 0.10` is dressed up and probably inverted | the constant is **0.052** and measured. Reads as a live finding against code that no longer exists |
| `docs/ARCHITECTURE.md:8-9` | its whole scope line: "53 code modules + 10 package `__init__.py`, ~10,250 LOC), every script under `scripts/` (12), and the contracts the 1,231-test suite pins down" | measured at `2b5dcde`: **58** code modules, **11** `__init__.py`, **~14,120** LOC, **21** scripts (+1 parked), and the collected count in §8. **All five numbers wrong, and drifting further** — `analysis/optimism.py`, `generator/families.py` and `scripts/run_cross_family.py` all arrived this morning. The sentence is the document's own statement of what it covers |
| `docs/ARCHITECTURE.md:1120` | quotes the README as claiming "149 tests" | the README makes no test-count claim |
| `docs/ARCHITECTURE.md` §2.6, §6.5, §8, §2's graph | pre-cleanup state | listed in the previous revision's §6; still unaddressed. §2 still shows the `analysis.power → reporter` edge that `f19ca4d` removed |

**The pattern behind N0, and it is not one commit's problem.** Three scripts take a flag
that changes the size of the experiment, none of them records which setting produced a
table, and headline numbers exist for both settings:

- `run_training.py --quick` — 4×3 vs 6×6. **All four tracked tables and a whole
  `NULL_RESULTS.md` section are the `--quick` setting** (N0).
- `run_transfer.py --scale` — 4 doses × 3 replicates × 2 seeds vs 6 × 6 × 5. The tracked
  `transfer_by_stressor.csv` is the default setting: 12 rows, mean R² **0.043**. The
  headline transfer triple in `CLAIM_BOUNDARY.md:78` and `FINDINGS.md:876-878` (0.064 /
  0.353 / 0.539 at 144 / 936 / 216 wells) corresponds to no cell in it. I could not tie
  those three figures to any tracked artefact — which is the finding, not that they are
  wrong.
- The activity-route comparison in `FINDINGS.md` (N3) has no producing script at all.

`REPRODUCING.md` opens by promising the opposite: *"Every number quoted anywhere in this
repo comes from a script in `scripts/` writing a table into `outputs/`."* One `print` of
the resolved configuration at the top of each `main`, or a settings column in each table,
would make all three self-describing. It is the cheapest structural fix in this document
and it closes the highest-severity finding in it.

`outputs/audit_claims.csv` records the repo's own commit count (`18 commits`, now 22), so
running the audit as documented dirties a tracked artefact every few commits. Restored
rather than committed. A count that changes on every commit does not belong in a tracked
table.

Cheapest structural fixes, in order of value per line: point `check_test_count` at
`docs/*.md` as well as the README (catches `ARCHITECTURE.md:9`); have `run_training.py`
announce its mode (catches N0); add `tests/` to `test_citations.py::_sources()` (catches
N1b). The power figures need a check that resolves `outputs/*.csv` cells named in prose —
real work, not a tweak, and it would catch six of the eight rows above.

## 7. What I could not verify

- **Whether the July-plate negative-activity fraction survives the total-signal route.**
  Needs `run_sensor_characterisation`-style handling of the two 2026-07 exports, which is
  a different code path. This is the check `FINDINGS.md` should carry and does not (N3).
- ~~**Whether the full-run training numbers are the right ones.**~~ **RESOLVED.** The
  tables are now full-mode, stamped, and reproduce byte for byte at HEAD (§8).
- **`9427394`, `21335394`, `2566299`, `27317316`** — the four formerly `tests/`-only PMIDs.
  They are now *in* the vendored table (N1b), so they are pinned to a resolved title; I did
  not check those four titles against the claims they sit beside. The mechanism is fixed;
  these four rows are unaudited content inside it.
- **S1's `np.std(x) < 1e-6` branch on real plates.** Still reasoned about, not exercised —
  but S1's *other* branch is now confirmed live on real data (N2), so the a priori argument
  that "it almost certainly never fires" has lost its supporting intuition.
- ~~**Whether the biphasic model's `identifiable=True` advantage is real or luck**~~
  **ANSWERED, and the answer is "partly luck".** `2b5dcde` replicated it: median **+0.236,
  range [−0.073, +0.301]**, one of three partitions losing outright. `38e76ed` narrows it
  further — the win holds below 0.5 mM and not above. The experiment §7 asked for was run
  the same morning.
- **`analysis/experiment_design.py:4`** — *"Eight single-agent rays already span all seven
  modules"*, at 24 modules. Still needs the analysis re-run, not a word changed.
- **`docs/FINDINGS.md:97` and `:113`'s 0.694.** Which parameter set produced it is still
  unknown.
- **Whether the seven remaining test failures mean the tests are stale or the recovery is
  spurious** (§0). `module_ec50("H2O2","peroxide")` moved from below the regulon it shares an
  agent with to above it, and a build that should reach at most 3 modules off its axes now
  reaches 7. Which side is wrong is a scientific call on the new EC50, and I did not make
  it.

## 8. Verification

Everything here was observed, not inferred. Where a number is stamped with a commit it is
because HEAD moved during the pass — three times.

**`python3 -m pytest tests/ -p no:cacheprovider`, run once at the end, whole directory, at
`2b5dcde`:**

```
10 failed, 1548 passed in 457.72s (0:07:37)
EXIT=1
```

**RED.** See §0 for the bisect and the seven that remain at `38e76ed`. For the record, the
second revision's run at `50d0529` was **1,464 passed, 0 failed, exit 0** — so this is a
regression introduced during the morning, not a pre-existing condition. This run used a
single `-q` (from `addopts`) rather than the `-qq` the documented command produces, which is
why a summary line exists to quote this time; see the note below.

**`python3 scripts/audit_claims.py`** — **37/37 checks pass.** Unchanged.

**Reproduction runs**, all through `YSTWIN_OUTPUTS` into scratch directories:

| Script | Result |
| --- | --- |
| `run_gates.py` | **all eight** `g1_*`/`d2_g1passed_*` tables byte for byte — the two that failed at the second revision now pass (N2) |
| `run_d2.py` | both `d2_*__mCitrine.csv` byte for byte |
| `run_training.py` (default) | **all four `training_*.csv` byte for byte** at HEAD — the condition that failed at the second revision (N0) |
| `run_heldout_score.py` | reproduces `heldout_scores.csv` byte for byte |

**The `[5, 3, 5]` disagreement is settled.** `docs/NULL_RESULTS.md` records that a review
pass reported latent-width instability at full size "which I could not reproduce". Both
readings are correct at their own commits — run in a detached worktree with the import path
forced:

| commit | five-channel latent width | modules recovered | score on driven modules |
| --- | --- | --- | --- |
| `50d0529` | **5, 3, 5** — and the `!!` warning fires | 4, 4, 5 | 0.223, 0.222, 0.223 |
| `2b5dcde` | **5, 5, 5** | 5, 5, 5 | 0.255, 0.244, 0.249 |

The instability was a property of the **pre-correction panel**. `7cbdf34` moved the H2O2
EC50 and `16c689b` replaced the posterior-predictive reference table, and between them the
width selection became stable. So the observation was real, the non-reproduction was real,
and the correction is what separates them — worth recording in `NULL_RESULTS.md` in place of
"unresolved", because it means the width selection is now robust for a reason rather than by
luck. The `2b5dcde` column matches that document's table exactly.

**Citations.** 53 table rows against every PMID cited in `src/`, `scripts/` **and** `tests/`
— no strays in either direction (N1b). The 13 corrected identifiers were checked against
NCBI E-utilities `esummary` at the second revision; titles, journals and years match.

**A sharp edge in the documented command.** `pyproject.toml` sets
`addopts = "-q --strict-markers"`, so `REPRODUCING.md` §3's `pytest tests/ -q` runs at `-qq`
and pytest prints **no summary line at all** — no "N passed", no "N failed", and no
`short test summary info` block listing which tests failed. A red suite run the documented
way shows a wall of dots with a few `F`s in it and nothing else. That is very close to how
the previous session came to report "green" without evidence, and it is one character to fix.

**Working tree.** I modified `docs/CLEANUP.md` only. Nothing was committed by me. The owner
committed `7cbdf34`, `abcad97`, `16c689b`, `2b5dcde` and `38e76ed` during the pass, three of
them carrying this report's fixes. Five detached worktrees were created under the scratchpad
for the bisect and removed afterwards.
