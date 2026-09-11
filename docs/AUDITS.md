# Audits

Four scripts that interrogate the repository rather than the biology. Each writes a
pass/fail table into `outputs/` and exits non-zero if anything fails, so any of them can
be a CI step or a pre-release gate.

They exist because the failures they look for are silent. A stale table, a path only one
machine has, a seed that is accepted and never used — none of these break a test, and all
of them make a number wrong in a way that is invisible until someone else tries to
reproduce it.

| Script | Question | Table |
| --- | --- | --- |
| `audit_claims.py` | Does the README describe the code that exists? | `outputs/audit_claims.csv` |
| `audit_reproducibility.py` | Could a stranger get these numbers? | `outputs/audit_reproducibility.csv` |
| `audit_determinism.py` | Does a seeded result mean what it appears to? | `outputs/audit_determinism.csv` |
| `audit_output_tables.py` | Does a committed table still come out of its generator? | `outputs/audit_output_tables.csv` |

```bash
python3 scripts/audit_claims.py
python3 scripts/audit_reproducibility.py            # add --quiet for failures only
python3 scripts/audit_determinism.py                # about 10 s
python3 scripts/audit_output_tables.py              # about 20 s; add --quiet
python3 scripts/audit_reproducibility.py --offline-suite   # minutes; see below
```

Every table has the same five columns — `check`, `expected`, `actual`, `status`, `detail`
— and `status` is `PASS`, `FAIL`, or `SKIP`. A `SKIP` is a check that could not run here
(cobra absent, no git, an optional dependency missing). It does not fail the audit and it
is not a pass: nothing was demonstrated.

---

## `audit_output_tables.py`

**Added 2026-09-04, and it closes the unguarded end of the one verification chain this
repository has.** `audit_claims.py` binds a number in prose to a **cell** in a committed
table. If the cell itself has drifted from the code that produced it, marking the prose
against it verifies only that two wrong numbers agree.

That was live. `outputs/proteomics_gain.csv` did not reproduce from
`scripts/simulate_proteomics_gain.py` at the script's own default `--seed 0`, even though
every draw comes off `np.random.default_rng(seed)` — so the committed table was made by a
code state that no longer exists. Every cell differed, and one column changed **sign**
(−0.7438 against +0.5773 at `protein_assay_cv_log = 0.2`), moving the script's own printed
conclusion from "it stops being worth it past CV 0.20" to 0.30. Nothing could catch it: no
prose cited the table, so no `audit:value` marker pointed at it, and `audit_determinism.py`
asks whether **seeded code** is reproducible rather than whether a **committed artefact**
still matches its generator.

Each registered script is re-run with `YSTWIN_OUTPUTS` pointed at a temporary directory —
the override `ystwin.paths.outputs_dir` documents for exactly this — and every table it
writes is diffed cell by cell against the committed one. Its first run found a second
defect: seven scripts printed their destination with `out.relative_to(paths.REPO_ROOT)`,
which raises for a path outside the repository, so they could not run under that override at
all. `paths.display_path` replaced it.

Scripts needing the genome-scale model, the wet-lab plate exports or the network are exempt
**by name with a reason** in `EXEMPT`, which is `audit_determinism.py`'s convention reused
rather than reinvented. `tests/test_audit_output_tables.py` runs the fastest registered pairs
so the gate cannot silently stop working.

It does **not** check that a table is correct — only that it is the table its generator
produces today. The three compose: this one says the cell came from the code, `audit_claims`
says the prose came from the cell, and the suite says the code does what it claims.

---

## `audit_reproducibility.py`

**Asset resolvers.** Calls every resolver in `ystwin.paths`, recording which environment
variable it consults by watching `_resolve` get called rather than by reading the source.
An asset that resolves to `None` **passes** — this repository is meant to work on a machine
that has never seen a plate reader. An asset that resolves *outside the checkout with its
variable unset* fails: that path was found by a convention (a sibling directory, a folder
on someone's Desktop) that nobody else's machine satisfies. Fix by setting the variable, or
by vendoring the file, not by relaxing the check.

**Override documented.** Every variable the resolvers consult must be named in
`docs/REPRODUCING.md`. A variable nobody documents is a variable nobody sets.

**Dependency pins.** Each requirement must be pinned with `==` *and* the pinned version
must be what `importlib.metadata` reports as installed. A drift between the two is exactly
the silent divergence pinning exists to prevent, and it is invisible from the file alone.
Optional groups that are simply not installed are skipped.

**Solver.** Reports what `cobra` resolves to and fails if it is not GLPK. Not a
correctness question — every solver returns *an* optimum — but a degenerate LP has a face
of optima and each solver picks a different vertex. The FVA widths in
`outputs/d1_capacity_sweep_ec.csv` are GLPK's.

**Writer map.** Reads each script's actual `to_csv` / `save` calls from the syntax tree and
compares them against the table in §4 of `docs/REPRODUCING.md`. Catches a script that
started writing a new table, and a documented table nothing writes. The guide is allowed to
be more specific than the source (`d2_<plate>__<reporter>.csv` against `f"d2_{tag}.csv"`).

**Outputs tracked.** Every file in `outputs/` must be tracked by git or deliberately
ignored. A result nobody else can see is not a result.

**Tables parse.** Every `outputs/*.csv` must be non-empty and hold at least one data row.
A header-only file is what a crashed script leaves behind and reads exactly like a result.

**Freshness.** Current verification uses explicit producer, input, model, parameter and
runtime identities, semantic comparison, and reviewed reproduction receipts. See
[Current verification and governed regeneration](REPRODUCING.md#current-verification-and-governed-regeneration)
and `data/artifact_registry.json` for the maintained procedure.

*Historical caveat:* the former audit used modification times and was inert on fresh
clones. That limitation is retained as history, not as an exemption from current
content-identity and reproduction checks.

**Offline suite** (`--offline-suite`, off by default, minutes). Runs the suite twice: once
as configured, once with every data-locating variable pointed at a directory that does not
exist. The second run must have zero failures and at least as many skips as the first.
This is the check that actually proves the "works without the wet-lab data" claim. Pointing
the variables at nothing — rather than unsetting them — is what makes it meaningful, since
`paths._resolve` treats a variable that is set as authoritative and will not fall back to
the sibling-directory default. Narrow it with `--suite-target tests/test_x.py` to smoke-test
the check itself.

---

## `audit_determinism.py`

Runs each seeded entry point **twice with one seed** and once with each of three others,
comparing a bit-exact digest of the result (`float.hex`, `ndarray.tobytes`) rather than
printed digits.

**Same seed, same answer.** Otherwise the table cannot be regenerated and a diff against a
re-run reports noise as change. Usual causes: a call into the global `np.random`, a set
iterated in hash order, an `id()` leaking into an ordering.

**Different seed, different answer.** The check people leave out, and the one that catches
a seed threaded halfway — accepted, handed to a generator, never drawn from. Such a
function reproduces perfectly and reports zero spread across seeds by construction.

That second property is not universal, so a specimen may declare `seed_sensitive=False`
**with a reason** — appropriate for a function returning a discrete choice, where the seed
moves the scores and not the winner. The expectation is asserted in *both* directions: a
specimen declared insensitive that turns out to be sensitive fails too, and one declared
insensitive with no reason fails immediately. Declaring it cannot silence a finding.

**Coverage.** Entry points are discovered by inspecting signatures across `ystwin`, not
listed. A seeded function with no specimen fails the coverage check the day it appears.
That nag is the point: a hardcoded list would leave new code permanently outside the audit.

---

## Adding a check

**To `audit_reproducibility.py`:** write a `check_*` function returning a list of rows from
`_row(check, expected, actual, status, detail)`, take whatever it inspects as an argument
rather than reaching for a module-level constant — that is what makes it testable against a
fixture — and call it from `collect()`.

**To `audit_determinism.py`:** add a `Specimen` to `SPECIMENS`. `target` is the
`module::name` discovery will report; `call` receives the resolved object and a seed. Keep
the input small — the whole audit runs in about ten seconds and should stay there — and
make sure the statistic is not saturated. A power estimate pinned at 1.0 cannot change with
the seed, and the sensitivity check would then report a fact about your specimen rather
than about the function.

**Either way, add the pair of cases to `tests/test_audits.py`**: a fixture in `tmp_path`
that satisfies the invariant and one that violates it, asserting `PASS` and `FAIL`
respectively. A check with only the passing half is a check nobody has verified detects
anything. Those tests never run the real audits over the real repository — that would be
slow, and it would turn a genuine finding about `outputs/` into a broken test.
