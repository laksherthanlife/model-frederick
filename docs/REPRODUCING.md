# Reproducing this repository

Current retained artifacts are registered by full relative path in
`data/artifact_registry.json`, with producers, inputs, parameters, model identity and
semantic contracts. `data/current_claims.json` binds scientific conclusions to those
records. Historical frozen runs use their original source/runtime; local exploratory
outputs do not acquire a current scientific role merely by existing on disk.

`data/gem/MANIFEST.md` is also a byte-frozen dependency of the historical native audit.
Its recorded `ecYeastGEM_yeast902.xml.gz` checksum predates the documented uptake-bound
change and SBML reserialization. The current compressed and decompressed model identities,
the earlier documented checksum, and the exact Git history are retained separately in
`data/public_inputs.json`. The verified current decompressed SHA-256 is
`ec743590f952836126215643757044975510ce95d1c6a164cd4426d03a847f22`;
the historical manifest is not rewritten to make its old metadata current.

Nothing here guesses where data is. Locations resolve through `src/ystwin/paths.py`: an
environment variable if you set one, otherwise a vendored copy, otherwise the
sibling-directory convention. A variable you set is used *instead* of the defaults, not
before them — if it points at a directory that does not exist, the asset counts as absent
and you get told which variable to fix, rather than silently getting someone else's copy.

**If you have just cloned this and have no wet-lab data, start at
[§6](#6-what-a-fresh-clone-can-reproduce).** It is a per-claim map of what your clone can
and cannot re-derive, checked by cloning this repository and running every script with
every data variable pointed at a directory that does not exist — not by reading the code
and reasoning about it. Several claims are honestly out of reach without files that are not
ours to publish, and §6 says which, rather than leaving you to discover it one refusal at a
time.

---

## Current verification and governed regeneration

The quality runtime is Python 3.14.2 with the numerical/solver pins and hashed dependency
lock. Use `requirements-quality.txt` for the complete verification environment; the
workflow additionally provisions the declared Ubuntu source-model tools. Public input
identities are checked by `scripts/provision_public_inputs.py`. Private original exports
are separate from mandatory public-data coverage and are never loaded from `.ystwin.env`
in CI. Use an activated virtual environment containing the pinned dependencies, rather
than a user-site installation: historical replay launches isolated Python (`-I`), which
does not see user-site packages. A successful import in ordinary system Python therefore
does not establish that its replay child has the required runtime.

```bash
python3 scripts/verify_quality.py --install \
  --cache-dir /tmp/ystwin-public-inputs --output-dir /tmp/ystwin-quality-new
```

Use a fresh external report directory. This command expects a clean tracked tree, retains
every child-process exit code, runs all blocking checks, and checks the tree again after
verification. Local work-in-progress can run the checks individually:

```bash
export CI=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -B -m pytest tests -ra --require-public-inputs -p no:cacheprovider
python3 -B -m ruff check --no-cache .
python3 -B scripts/audit_claims.py
python3 -B scripts/audit_reproducibility.py
python3 -B scripts/replay_portable_evidence.py --historical
```

Audits are read-only unless an explicit diagnostic destination is supplied. Historical
replay executes verified original Git source with its recorded interpreter/package
versions; it does not fit again or turn frozen results into current validation. A changed
supplied historical data/metadata input is refused; use the verified original source
snapshot, not a new hash on the old result.

The final native-reconciliation producer and its original test were recovered from Git
objects and retained as inert `.py.source` evidence. Their exact identities and the
original-runtime attempt are recorded in `data/native_reconciliation_recovery.json`.
That attempt remains blocked before numerical work: the original checkpoint binds the
absolute checkout root, even when every code/input record matches. Source recovery is
not successful relocated replay. Do not patch the frozen validator, rewrite the original
checkpoint, substitute the current consolidated runner, or overwrite the live worktree
to bypass this boundary.

```bash
python3 scripts/regenerate_current_artifacts.py --plan
python3 scripts/regenerate_current_artifacts.py --run product-prediction \
  --output-dir /tmp/ystwin-product-staging-new
```

Staging does not replace retained results. A nonzero comparison exit may mean the producer
succeeded but values, missingness, units, row identities or interpretation changed; inspect
`execution_complete`, the producer exit code, contracts and differences in the receipt.
Only unchanged, verified reproductions are eligible for metadata-only adoption. Changed
results require an exact-hash review with reasons and explicit acknowledgment of negative
case changes, followed by explicit `--adopt-receipt`, `--accept-scientific-changes`,
`--review` and `--install`. The tool rejects stale code/input/runtime identities. Neither CI
nor the audits perform this promotion, and `--restamp` is disabled.

The detailed sections below include historical experiment notes. For a current invocation,
the registered recipe and its evidence scope take precedence over an old default command.

### Corrected comparison contracts during integration

The following corrected producers have been inspected through external candidates; final
retained-artifact regeneration/review is still a separate step. A source/runtime binding
that changed after staging makes a receipt stale even if its numbers still agree. Do not
rewrite its identity, install a candidate merely to make tests green, or treat an
unapproved review as adoption.

| workflow | complete population and required interpretation |
| --- | --- |
| `env-to-product` | all eleven PHB measurements and all six β-carotene requests; four answers and two `SetpointUnreachable` ethanol refusals at requested μ = 0.15 and 0.2543 /h; missing returned growth/content and layer flag on refusals, environment layer `not-run` rather than `reported` |
| `cross-family` | all seven families; JSON `configuration` / `transfer_configuration` for fixed three-state transfer, `selected_configuration` for the named-family search, and `candidate_configuration` for the separately sampled SPOTA solution; the corrected default candidate is three/three/two states, not one shared configuration |
| `flux-null-baseline` | five parameter-count rows in `flux_null_baseline.csv` as a random-design demonstration, plus four fixed-candidate rows in `flux_candidate_scores.csv`; both model and comparator use each fold's training data, with `baseline=fold_training_mean_log_flux` and `inference_status=pending_exchangeability` for every real-data candidate |

PHB content holdouts are retrospective within one Kocharin cohort; the shipped law factor
is reported separately, not scored as a held-out refit. Equal answered β-carotene values
are the fixed-gene μ-only construction, not evidence of biological invariance. Only
μ = 0.101 /h has a complete answered carbon pair; one answered and one refused at the other
rates is an unavailable comparison, not an invariant pair.

Cross-family scores use outer-training channel means, and the oracle is in-sample. The
named-family score difference is descriptive, while the approximate SPOTA bound concerns
only its independently fixed candidate's expected gap under the declared synthetic
sampler. `selected_pipeline_optimism_status` remains pending. Non-rejection of the rotated
null is not equivalence or absence of information. The historical default/scale tables in
[`CROSS_FAMILY.md`](CROSS_FAMILY.md) are not new current estimates.

For flux, the constant has zero skill to roundoff and all three nonconstant candidates
remain negative in the inspected candidate. No candidate significance threshold or
p-value is issued: condition summaries and feed labels do not establish exchangeable
biological units. Preserve failed fits/draws rather than filtering them, and never reuse
the random-design percentiles as a universal null bar. The old full-data-mean comparison
and roughly 0.08 identifiability rule are explicitly withdrawn in
[`EXTERNAL_PRODUCT_VALIDATION.md` §3](EXTERNAL_PRODUCT_VALIDATION.md).

`tests/test_env_to_product.py` deliberately reads the retained root artifact. Until its
governed replacement is installed, old six-answer/all-`reported` data must fail the new
refusal and scope assertions. A separate read-only check against an explicitly named
external candidate can test those expectations without replacing the root file, but that
pass is candidate validation, **not a passing retained-artifact gate**. Report both
outcomes. No missing row or layer state should be converted to a skip or an invariance
claim to avoid that block.

### COSMIC source supplements

The supplied version-of-record paper and supplements are now inspected and identified in
`data/cosmic_sources.json`. Keep their originals outside the repository. To inspect the
pinned workbook without modifying it, select the input location explicitly and use a new
external report filename:

```bash
python3 -B scripts/inspect_cosmic_data.py \
  --input-dir /absolute/path/to/cosmic-dfba-data \
  --report /tmp/cosmic-source-inspection-new.json
```

The reader verifies the original workbook SHA-256 and byte count, preserves source cells,
blank slots, signed values, units and vessel IDs, and does not infer missing mappings.
Renamed duplicate workbooks are verified as the same source, not additional replicates.
The source covers ten CHO vessels, eight factor combinations and 130 sampling records;
day 8 is absent from the supplied observations and is not interpolated.

ST2's process quantities are normalized measurements, while its phase fractions are
inferred by nonlinear regression. ST3 rates, ST4 task priorities/efficiencies and ST5 fluxes
are computed source outputs, not independent measured validation targets. ST5 supplies
reaction/time axes but no vessel identifier or stoichiometric model. Coded oxygen/feed
levels and normalized readouts cannot be substituted for physical yeast inputs.

The paper's predictive classifier uses **LDA followed by logistic fitting**; its PCA is
an exploratory analysis of computed fluxes. The repository's optional PCA/soft-label
allocation model remains a separately declared implementation, not an exact copy.
An exact COSMIC replay still requires the specific CHO metabolic/secretory model,
complete fitted kinetics/classifier and `dFBA_data`, physical scales and feed/retention
configuration. The kernel's `0.35` conversion is not transferred without its missing basis.
The inspector does not claim model reproduction or independent biological validation.

## 1. Install

```bash
python3 -m pip install --require-hashes -r requirements-quality.txt
python3 -m pip install --no-deps --no-build-isolation -e .
```

Python 3.11 or newer. Dependencies are pinned with `==`, deliberately: this project
reports FBA fluxes to four significant figures and picks latent dimensions by comparing
scores that differ in the third, so a cobrapy or numpy minor release moves digits that
appear in `outputs/`. See the comment in `pyproject.toml`.

`cobra` solves through **GLPK** here (via `optlang`/`swiglpk`). No commercial solver is
installed, and that is the reproducible choice — the FVA widths in
`outputs/d1_capacity_sweep_ec.csv` are the ones GLPK reports, and CPLEX or Gurobi can
return a different vertex of the same degenerate optimum.

The thermodynamic bridge needs two more packages, which are not part of the default
install because nothing on the main pipeline imports them:

```bash
python3 -m pip install -e ".[dev,thermo]"
```

Without them, `tests/test_thermodynamic_*.py`, `tests/test_tmfa.py` and
`tests/test_equilibrator_backend.py` skip.

---

## 2. Data you need, and where to put it

The first two rows are what the sensor-characterisation pipeline needs; everything below
them extends coverage. Nothing here is required to run the test suite. A missing asset
always produces a skip or a named refusal, never a wrong number.

**Both wet-lab plate sets are now committed as text**, and that is what makes the row below
the important one. `data/plates/` holds the four NewProtocol replicates *and* the two
2026-07 ER exports, every value the reader wrote, verified against the workbooks
value-for-value by `plate_readings.py --export` and carrying each workbook's sha256 so a
holder of the original can prove provenance without the original being here. Neither set of
workbooks can be redistributed — both carry a named private individual in
`docProps/core.xml` — and neither needs to be.

Until 2026-08-29 only the NewProtocol set had that treatment, so `paths.igem_results()` fell
back to `../igem-results`, a sibling of the checkout. `audit_reproducibility.py` failed on it
for months and was right to: results computed through a path found by convention are
reproducible on one machine. That fallback is gone, the variable is the only way to reach the
workbooks, and every script reads the committed text instead.

| Asset | Default location | Variable | Needed by |
| --- | --- | --- | --- |
| **Biosensor plate readings, as committed text** | **`data/plates/` — tracked, ships with the clone** | **—** | **`plate_readings.py`, `run_gates.py`, `run_d2.py`** |
| NewProtocol biosensor plates (4 Synergy exports) | `../plates/Biosensor Testing`, else `~/Desktop/Result & Analysis/Plate Reader Result/Biosensor Testing` | `YSTWIN_PLATES` | `run_sensor_characterisation.py`, `run_g4.py` |
| The two 2026-07 ER **workbooks** | none — set the variable | `YSTWIN_IGEM_RESULTS` | re-exporting or re-verifying `data/plates/`, and the parser tests |
| Crosstalk workbook `ER&OX-summary.xlsx` | none — set the variable | `YSTWIN_CROSSTALK_WORKBOOK` | re-deriving `data/crosstalk/` only; the derived tables are tracked |
| Gen5 `.xpt` instrument files | none — set the variable | `YSTWIN_GEN5_XPT` | `read_gen5_xpt.py` and `tests/test_gen5_xpt.py`; the raw source of every plate |
| Bio-Rad Cq exports (11 and 13 Aug) | `../qpcr/ER and Oxidative Stress Biosensors`, else `~/Desktop/Result & Analysis/rt-qPCR results/ER and Oxidative Stress Biosensors` | `YSTWIN_QPCR_RAW` | `run_g4.py` |
| Repaired 2026-07-24 Cq export | `data/qpcr_repaired/` | `YSTWIN_QPCR` | `run_g4.py` |
| ATP-sensor plates | `data/atp_sensor/` | `YSTWIN_ATP_SENSOR` | `analyse_atp_sensor.py`, `tests/test_atp_sensor.py` |
| Yeast9 `yeast-GEM.xml` | `data/gem/`, else `../dente/data/` | `YSTWIN_YEAST_GEM` | `run_scenarios.py`, `parked/run_d1.py`, the FBA and thermodynamics tests |
| GECKO `ecYeastGEM_batch.xml` | `data/gem/`, else `../dente/data/` | `YSTWIN_EC_YEAST_GEM` | `parked/run_d1.py`, the ec-model tests |
| **ModelSEED + pytfa thermo tables** | **`data/thermo/` — tracked, ships with the clone** | `YSTWIN_THERMO` redirects | `bridge/thermodynamic.py` and its tests |
| **BioNumbers export** | **`data/kaggle/` — tracked, ships with the clone** | **—** | `tests/test_constants_against_bionumbers.py` |
| **GEO + SGD external validation data** | **`data/external/` — tracked, ships with the clone** | **—** | `run_external_validation.py` and the external-validation tests |
| **Jalihal 2021 nutrient-signalling SBML** | **`data/native_reference_models/jalihal2021/` — tracked** | **—** | `tests/test_native_reference_models.py` |

`../` means a sibling of this checkout, i.e. inside the directory that contains
`model-v2`. Four gitignored directories remain under `data/` -- `atp_sensor/`,
`qpcr_repaired/`, `native_law_v2/granados/sealed/` and `.../development/private/` -- plus
the single file `holdout_transfer/source_index.json`. Each carries a comment in
`.gitignore` saying what it is and why it stays out.

Every public third-party download is now tracked and ships with the clone. Until
2026-09-11 the thermo, BioNumbers, GEO/SGD and Jalihal tables were ignored as
"re-fetchable", but the tests that read them guard with a module-level
`pytest.mark.skipif(not PATH.exists())`, so a clone without them ran fewer checks than the
author's machine without ever saying so. They are committed for that reason, not for
convenience.

What is still withheld is withheld deliberately: wet-lab exports are not redistributable
from here (get them from the team), and the sealed holdout must stay sealed or the
transfer claims that rest on it mean nothing. The genome-scale models are public (Yeast9
from the `SysBioChalmers/yeast-GEM` releases, the GECKO batch model from
`SysBioChalmers/ecModels`), as are the thermodynamic tables.

`YSTWIN_OUTPUTS` redirects where scripts write, which is the safe way to regenerate a
table without overwriting the tracked one:

```bash
YSTWIN_OUTPUTS=/tmp/check python3 scripts/run_gates.py
diff outputs/g1_20260701_ER_preliminary_\(RAW\).csv /tmp/check/
```

`outputs_dir()` is where some legacy scripts *read* as well as write. Redirecting it can
hide the committed inputs from `run_calibration.py`, `run_heldout_score.py`,
`run_training.py` and `make_figures.py`. Supply their required inputs explicitly; do not
use in-place regeneration as a verification shortcut.

The claim, reproducibility, determinism and output-table audits are read-only by default.
Use `--output-dir /tmp/ystwin-diagnostics` when a report is needed. That destination is
for diagnostics, not retained scientific results. CI never restamps provenance and fails
if verification modifies a tracked file.

---

## 3. Run the tests

```bash
python3 -m pytest tests/ -q
```

Several minutes, most of it SBML parsing and mixed-integer solves. With no wet-lab data
and no GSMM the suite still passes — the real-data tests skip. To see that for yourself,
point the variables at nothing:

```bash
YSTWIN_PLATES=/nowhere YSTWIN_IGEM_RESULTS=/nowhere YSTWIN_QPCR=/nowhere \
YSTWIN_QPCR_RAW=/nowhere YSTWIN_YEAST_GEM=/nowhere YSTWIN_EC_YEAST_GEM=/nowhere \
YSTWIN_ATP_SENSOR=/nowhere YSTWIN_THERMO=/nowhere \
  python3 -m pytest tests/ -q
```

Two markers are registered: `integration` (touches real data or the full GSMM) and
`slow`. `-m "not integration"` skips the first.

---

## 4. Which script writes which file

Built by reading the `to_csv` / `save` calls, not from memory. `<plate>` is the export's
filename stem with spaces and `&` substituted, so the set of files depends on which plates
you have.

| Script | Writes into `outputs/` | Needs | Measured |
| --- | --- | --- | --- |
| `run_gates.py` | `g1_<plate>.csv`, `d2_g1passed_<plate>.csv`, `gates_manifest.json` | `YSTWIN_IGEM_RESULTS` | seconds per export; both survey every `.xlsx` |
| `run_d2.py` | `d2_<plate>__<reporter>.csv` | `YSTWIN_IGEM_RESULTS` | seconds per export, and most are skipped |
| `run_sensor_characterisation.py` | `sensor_characterisation.csv`, `autofluorescence_sensitivity.csv`, `late_window_sensitivity.csv` | `YSTWIN_PLATES` | ~20 s |
| `plate_readings.py` | `sensor_characterisation.csv`, `autofluorescence_sensitivity.csv`, `late_window_sensitivity.csv` — all three only with `--write`, and the last only with `--window` | `data/plates/`, or `YSTWIN_PLATES` if you have it | ~15 s, or ~4 min with `--window`; `--export` needs the plates |
| `run_calibration.py` | `panel_calibration.csv` | `outputs/sensor_characterisation.csv` | seconds |
| `run_g4.py` | `g4_anchor_fold_change.csv`, `g4_reporter_response.csv`, `g4_rt_minus_qc.csv`, `g4_verdicts.csv` | `YSTWIN_QPCR`, `YSTWIN_QPCR_RAW`, `YSTWIN_PLATES` | ~20 s |
| `run_power.py` | `power_analysis.csv` | nothing — simulation | 81 s |
| `run_transfer.py` | `transfer_by_stressor.csv`, `module_recovery.csv` | nothing — simulation | 10 s; `--scale` is far longer |
| `run_training.py` | `training_recovery.csv`, `training_transfer_wide.csv`, `training_transfer_build.csv`, `training_transfer_three_sensor.csv`, `stress_model_wide.npz`, `stress_model_build.npz`, `stress_model_three_sensor.npz` | reads `sensor_characterisation.csv` if present, as the posterior-predictive reference | 39 s with `--quick`; minutes without |
| `run_scenarios.py` | `scenario_predictions.csv` | **GSMM** (`YSTWIN_YEAST_GEM`) | 18 s — most of it the SBML parse |
| `parked/run_d1.py` | `d1_capacity_sweep_ec.csv` | **both GSMMs** | slowest: parses Yeast9 *and* the GECKO model, then FVA |
| `register_prediction.py` | `registered_prediction_D<rate>.csv` | `data/carotenoid/elizondo2025_steady_states.tsv` (tracked) | instant |
| `analyse_atp_sensor.py` | `atp_sensor_<block>.png` — three, one per plate block, and only with `--figures`; prints otherwise | `YSTWIN_ATP_SENSOR` | seconds |
| `repair_qpcr_export.py` | its second argument | an unopenable Cq export | seconds |

The two GSMM scripts are the ones that need a genome-scale model, and the SBML parse
dominates both: `run_scenarios.py` reads Yeast9 once, `parked/run_d1.py` reads Yeast9 and
the 9.8 MB GECKO model.

Run order matters twice: `run_calibration.py` reads the table
`run_sensor_characterisation.py` writes and exits if it is missing, and `run_training.py`
uses the same table as its reference when it exists.

`gates_manifest.json` is the one non-CSV in that table, and it is the run record rather
than a result: one entry per declared gate export, each carrying `status` (`processed`,
`refused`, or an exclusion reason), `source_set` and the `source_sha256` the export was
matched on, above a `complete` flag and the processed/refused counts. A declared current
source that never reached the reader is written in as `refused` rather than left out, so
the file says which plates the committed `g1_*.csv` and `d2_g1passed_*.csv` tables do and
do not cover. `run_gates.py` exits **2** when any entry is refused — the tables it wrote
on that run are partial, and this manifest is where you find out which ones. An
unavailable source set is not a refusal on its own: the committed plate text is replayed
first, and only what survives that is judged.

Naming it in the table above is what the table owed a reader, and it is **not** what closes
`writer map: run_gates.py` in `scripts/audit_reproducibility.py`. That check reads the
documented side with `writers_from_doc`, whose filename pattern accepts `.csv`, `.npz`,
`.png` and `.xlsx` and nothing else, so the one `.json` in this table is invisible to it
however carefully it is written. The source side already finds `outputs/gates_manifest.json`,
and the check reports it as written but not documented. Adding `json` to that alternation
closes the failure and introduces no other, measured by running the comparison both ways
against this file; the fix belongs in the audit, not here. Until it lands, do not add a
`run_calibration_nis.py` row here either — it writes `nis_manifest.json`, so a row for it
would fail the same way for the same reason.

The three `atp_sensor_*.png` panels used to be the entry here: made out-of-band, in
`outputs/`, produced by nothing. They were moved to `archive/outputs/`, and
`analyse_atp_sensor.py --figures` now draws the part of them that comes from the plates
(§6 says which part, and why the rest cannot be). Its output is gitignored, like every PNG
here that is not part of the story set.

`qpcr_fold_change.csv` used to be the second entry here, tracked and irreproducible, and
`audit_reproducibility.py` failed on it for three months while this file documented the
failure in three places. **It is now deleted rather than documented**, because it was not
merely irreproducible — it was superseded. It carried the same nine columns as
`outputs/g4_anchor_fold_change.csv` over two of the three qPCR replicates; `run_g4.py`
computes the same quantity over all three, nothing in the repository read the old file, and
no claim in any document rested on it. A gap that a script can close is a gap; a duplicate
of a table that a script already writes is just an old copy.

The `.npz` stress models and the PNGs are gitignored on purpose — binary, regenerable, and
they diff to nothing readable. Everything else written into `outputs/` is tracked.

**Not every tracked table is current.** Re-running against the code as it stands today
reproduces `sensor_characterisation.csv`, `scenario_predictions.csv` and all six `g1_`/`d2_`
tables byte for byte, and `g4_anchor_fold_change.csv`, `g4_reporter_response.csv` and
`g4_verdicts.csv` as well. Four disagree, all because the library moved after the table was
written, not because the inputs did:

- `g4_rt_minus_qc.csv` — `qpcr.rt_minus_margin` now also emits `gdna_share` and
  `correctable`, which the committed file predates.
- `power_analysis.csv`, `module_recovery.csv`, `transfer_by_stressor.csv` — the `analysis/`
  modules behind them have changed; `module_recovery` and `transfer_by_stressor` differ
  substantially, not in the last digit.

`d1_capacity_sweep_ec.csv` and the four `training_*` tables are unchecked: the first needs
both GSMMs, and the second set was written by a full `run_training.py` rather than the
`--quick` run used here, so a diff against `--quick` output would mean nothing.

One correction to "byte for byte" above: `sensor_characterisation.csv` reproduces to 15
significant digits, not to the last bit. The tracked file was written with 15, so on 15 of
its 251 rows the largest `naive_late` values lose their final bit — `9053.58689279822` for a
value that is `9053.586892798221`. One part in 10^15, below anything this repository
reports, and measured rather than assumed: see §7.

Regenerate them once `analysis/` settles. Being able to see that they are stale is the
point of tracking them.

---

## 5. When something is missing

A script that cannot find its data prints one line naming the variable to set and exits 1.
It does not fall back, and it does not write an empty table — an empty gate report reads
like a plate that passed nothing rather than like a plate that was never there.

A library function refuses the same way, with `FileNotFoundError` from
`paths.require`. Callers that should skip instead check the resolver for `None`, which is
how the test fixtures work.

---

## 6. What a fresh clone can reproduce

Checked, not reasoned about. The table below comes from cloning this repository into a
temporary directory, copying in nothing but `data/plates/`, and running each script with
every data-locating variable pointed at a directory that does not exist:

```bash
git clone <this repo> /tmp/clone && cd /tmp/clone
export YSTWIN_PLATES=/nowhere YSTWIN_IGEM_RESULTS=/nowhere YSTWIN_QPCR=/nowhere \
       YSTWIN_QPCR_RAW=/nowhere YSTWIN_YEAST_GEM=/nowhere YSTWIN_EC_YEAST_GEM=/nowhere \
       YSTWIN_ATP_SENSOR=/nowhere YSTWIN_THERMO=/nowhere
```

A clone gets **five** data files under `data/` — four citation tables and
`data/physiology/chemostatData_VanHoek1998.tsv` — plus `data/plates/`. There is not one
tracked `.xlsx` and not one tracked `.xml`, and there should not be: see §7 for why the
plate workbooks are committed as text instead of verbatim.

### The real-data claims

| Claim | Command | Needs | Clone? |
| --- | --- | --- | :---: |
| README §1 — dose response is mostly dilution (the naive-vs-corrected folds) | `python3 scripts/plate_readings.py` | `data/plates/` | **yes** |
| README §3 — the autofluorescence sweep, and which calls it flips | `python3 scripts/plate_readings.py` | `data/plates/` | **yes** |
| `outputs/sensor_characterisation.csv` — 251 well rows | `python3 scripts/plate_readings.py --write` | `data/plates/` | **yes** |
| `outputs/panel_calibration.csv` — panel EC50s against the real ladders | `python3 scripts/run_calibration.py` | committed `sensor_characterisation.csv` | **yes** |
| `outputs/heldout_scores.csv`, `heldout_interpolation_by_dose.csv` | `python3 scripts/run_heldout_score.py` | committed `split_manifest.csv` + `sensor_characterisation.csv` | **yes** |
| `outputs/late_window_sensitivity.csv` — the window sweep | `python3 scripts/plate_readings.py --window` | `data/plates/` | **yes** |
| README §2 — per-plate fluorescence coverage (`docs/DATA_INVENTORY.md`) | read `data/plates/manifest.csv` | `data/plates/` | **yes** |
| `outputs/g1_*.csv`, `d2_g1passed_*.csv` — the optical gate | `python3 scripts/run_gates.py` | **raw plates** + `YSTWIN_IGEM_RESULTS` | no |
| README §4 — qPCR is 91.4 % genomic DNA (`outputs/g4_*.csv`) | `python3 scripts/run_g4.py` | **qPCR exports**, not redistributable | no |
| `outputs/nis_summary.csv`, `nis_innovations.csv` — is the filter honest? | `python3 scripts/run_calibration_nis.py` | `data/plates/` | **yes** |
| `outputs/scenario_predictions.csv` — FBA scenarios | `python3 scripts/run_scenarios.py` | **Yeast9 GSMM** | no |
| `outputs/d1_capacity_sweep_ec.csv` — FVA widths | `python3 scripts/parked/run_d1.py` | **both GSMMs** | no |

That second row needs a word of care, because it is the one place where the committed text
settles an argument the docs are still having. `manifest.csv` gives `n_wells` per block, and
it reads **87 fluorescence wells on all four plates** — not the 21 that
`docs/DATA_INVENTORY.md` tabulates twice. 21 is what the reader returned before the
down-the-page block layout was handled; 87 is what the plates actually contain. So the
coverage claim is auditable from a clone without running anything, and what it shows is the
corrected count. `docs/DATA_INVENTORY.md` says as much in its own prose ("Both carry all 87
wells, blanks included") while leaving the superseded tables in place, and only `20260804`
is genuinely unusable — it is already blank-subtracted, so there is no background left to
recover.

### The simulation claims

Every one of these is reproducible on a clone, because none of them touches wet-lab data.
That was already true before `data/plates/` existed and is not the gap this page is about.

| Claim | Command | Clone? |
| --- | --- | :---: |
| README §5 / the Ledermann bound — four channels reach one factor | `python3 scripts/run_power.py` | **yes** |
| `outputs/transfer_by_stressor.csv`, `module_recovery.csv` | `python3 scripts/run_transfer.py` | **yes** |
| `outputs/cross_family_*.csv` | `python3 scripts/run_cross_family.py` | **yes** |
| `outputs/training_*.csv` — the trained model | `python3 scripts/run_training.py --quick` | **yes** |

### Where a clone was quietly worse off — fixed at the source

`scripts/make_splits.py` **used to exit 0 on a clone and write a smaller table than the one
it replaced.** The manifest has two halves; without `YSTWIN_PLATES` only the simulated half
was built and the tracked file went from 1650 rows to 936 — the whole `real_biosensor`
dataset gone. It said so on stdout, which was the honest half. The dishonest half was the
exit status and the overwrite: the file afterwards looked like a manifest rather than half
of one, and `scripts/run_heldout_score.py` reads it. It died with
`'DataFrame' object has no attribute 'plate'`, three scripts away from the missing variable.

It now refuses, naming `YSTWIN_PLATES`, and leaves the manifest alone. `--panel-only` says
the smaller manifest is intended. The refusal is deliberately narrow: a clone with no
manifest at all builds the simulated half freely, because there is nothing to destroy and
refusing would strand the first run.

This page previously called that "the one place in the repository where absent data produces
a smaller result instead of a refusal". It was not the one place. `scripts/run_g4.py` did
the same thing with qPCR replicates — one of three instead of three, turning a measured
anchor effect into `NaN` and a replicate requirement into `0`, which in the table is
indistinguishable from "there is no anchor at all". Both now go through
`ystwin/gates/artefact.py::refuse_partial_rebuild`. See
[`AUDIT_2026_08.md`](AUDIT_2026_08.md) §1.

### What is still not reproducible, stated plainly

- **Nothing downstream of the qPCR exports.** The G4 anchor — the 91.4 % contamination
  figure, the strongest negative result here — cannot be re-derived from this repository by
  anyone. The Cq exports are wet-lab files and one of them only opens after
  `scripts/repair_qpcr_export.py` rebuilds it.
- **Nothing needing a genome-scale model.** Both are public and re-fetchable (§2), so this
  is a download rather than a dead end.
- ~~**The optical gate and the NIS filter check**~~ — **done.** `run_gates.py`,
  `run_calibration_nis.py` and `run_g4.py` now read through `ystwin/plate/replay.py`, which
  is the loader that used to live inside `plate_readings.py` as a private shim. Each prints
  which source it used. `run_gates.py` reads two plate sets and only one is committed, so
  the swap is per export, on whether the manifest covers the name. The NIS row in the table
  above said **raw plates** and `no` for months after that became false; it now says
  `data/plates/` and `yes`, checked the way this section says it checks — cloning this
  repository, exporting every data variable to a directory that does not exist, and running
  the script, which reads the committed text and reports the same 319 well-channels and the
  same 45.3 / 19.3 aggregates as the tracked `nis_channel_summary.csv`. `run_gates.py` stays
  `no`: it needs the two July exports, which have no committed replay.
- **`run_d2.py`**, which reads the two July Synergy exports. Those are not in
  `data/plates/` — the committed text covers the NewProtocol replicates and the BY4741
  control — so there is no replay to give it.
- ~~**The three `outputs/atp_sensor_*.png` figures**~~ — **half closed, and the half that
  is not is stated rather than left open.** The gap was invisible until August 2026 because
  `viz/figures.py` writes `out_dir / f"{stem}.svg"` and the resolved pattern, `*.png`,
  matched every PNG in the repository; attribution now comes from `figures_manifest.csv`
  and the three ATP panels had no row in it. `scripts/analyse_atp_sensor.py --figures` is
  the writer they never had, and it needs `YSTWIN_ATP_SENSOR`, so a clone still cannot run
  it. What it draws is **every panel of those renders whose numbers come from the plates**
  — the traces with their peak hours, the clock-time-against-matched-density collapse, and
  each peak against the hour its culture stopped growing — reproducing the archived values
  to the digit (3.66× at 10 h, 1.13× at OD 0.5, peaks at h12.6/h11.5/h9.0). What it does
  not draw is the three panels that were never measurements on these plates: a literature
  bar, a model-against-data comparison whose model constant has moved since, and an
  analytic curve with a free parameter. `docs/research/XPT_INVENTORY.md` records the split
  panel by panel, and the archived renders stay in `archive/outputs/` because a panel this
  repository cannot regenerate should not sit in `outputs/` looking like one it can.

The raw plates, the qPCR exports and the ATP-sensor plates are obtainable from the project
owner rather than from this repository. The genome-scale models and the thermodynamic
tables are public and §2 says where.

---

## 7. The committed plate readings — `data/plates/`

### Why the workbooks are not tracked, and it is not size

The four NewProtocol Synergy exports are 989 kB together, which is small enough to commit.
They are not committed because **every one of them names a private individual** in
`docProps/core.xml` (`cp:lastModifiedBy`). Committing the workbook publishes that name into
git history permanently, where removing it costs a history rewrite. Stripping the metadata
first would leave a binary that no longer matches the lab's own copy while still diffing to
nothing legible, which is worse than either alternative.

So what is tracked is the measurement matrix and nothing else: every value the reader
wrote, as text. No author, no revision history, no embedded workbook metadata — optical
densities, relative fluorescence, elapsed times, temperatures, and the per-construct sheets
that carry the dose ladder.

The arithmetic, for the record:

| | bytes | |
| --- | ---: | --- |
| the four biosensor `.xlsx` | 1 012 796 | binary, each one carries the PII |
| the BY4741 OD600 plate beside them | 11 467 | same directory, read by the same glob |
| **all five workbooks** | **1 024 263** | diffs to nothing legible |
| **`data/plates/`, 17 files** | **248 386** | text, no PII, diffs line by line |

24 % of the binary, and it is the honest artefact rather than a lossy summary of it. The
BY4741 plate is included because `run_sensor_characterisation.py` globs the whole directory
and the run is not the same without it — it is the reporter-free control that was never read
in the fluorescence channel, and `collect()` has to see it in order to drop it.

### Fidelity is checked, not asserted

`python3 scripts/plate_readings.py --verify` rebuilds every block from the committed text
and compares it to the workbook with `==` on each float, not `np.isclose`. All 13 blocks
across the 5 exports round-trip exactly: every value, every elapsed time, every well
ordering, and the block metadata beside them.

Two details make that possible and both are load-bearing:

- **Elapsed times are stored as `H:MM:SS`** and rebuilt with the same
  `h + m / 60 + s / 3600` expression `plate/synergy.py::_to_hours` uses.
  `round(t * 3600) / 3600` is the same real number and a different float.
- **One file per block, not per export.** `20260804` reads its two channels in two
  different well orders, and that ordering reaches the committed output — it sets the order
  of `shared` in `prepare()` and therefore the row order of `sensor_characterisation.csv`.
  The CSV header *is* the well ordering, so nothing has to record it separately.

Running the pipeline from the committed text gives a frame **bit-for-bit identical** to
running it from the workbooks — worst absolute difference exactly 0 across all 16 columns
and 251 rows.

Against the tracked `outputs/sensor_characterisation.csv` both routes land 1.8e-12 away, in
`naive_late`, on 15 of 251 rows. That gap is not a difference between the two sources: the
tracked file was written with 15 significant digits (`9053.58689279822` for a value that is
`9053.586892798221`), so it loses the last bit on the largest values. It is one part in
10^15 and below any precision this repository reports, but it means the committed table is
*not* byte-reproducible from either source, and §4's claim that it reproduces byte for byte
is true only to 15 digits.

### Schema

`manifest.csv` — one row per kinetic block, in the order `read_synergy_kinetic` produced
them. That order matters: `raw_channel` breaks a tie between two equally wide raw blocks by
taking the first, and `aligned` builds its column index by walking the blocks in order.

| column | meaning |
| --- | --- |
| `export` | the original workbook's filename, exactly. The logbook lookups in `plate/layout.py` are keyed on it, and `plate` in `sensor_characterisation.csv` is its first 22 characters |
| `channel` | unique label — `mCitrine[1]`/`[2]` when a fluorophore was read twice, since repeat reads differ in gain |
| `fluorophore` | base name without the repeat suffix |
| `optics` | excitation,emission or absorbance wavelength as exported, blank when the sheet did not carry it |
| `sheet` | the worksheet it came from, kept for provenance |
| `derived` | the block sits on a "for plotting" sheet with summary columns among the wells |
| `blank_subtracted` | the block holds negative values, which raw optics cannot |
| `n_times`, `n_wells` | shape, so a truncated file is visible without parsing it |
| `readings_file`, `dose_response_file` | the files below; `dose_response_file` is empty for an export with no per-construct sheets |
| `source_sha256`, `source_bytes` | the workbook's digest. A holder of the original can prove these came from it without the original being here |

`<export>__<channel>.csv` — one block. `elapsed_hms`, then `temperature_c` when the reader
recorded it, then one column per well **in the block's own order**. Optical densities are
exact at three decimals and relative fluorescence is integer, as the instrument reports
them.

`<export>__dose_response.csv` — the per-construct sheets: `construct`, `dose_mM`,
`elapsed_min`, `signal`. Kept because dropping them changes the answer rather than only the
provenance: without them `collect()` takes the dose ladder from the logbook registry
instead of from the plate, and the two are not guaranteed to agree.

### Using it

```bash
python3 scripts/plate_readings.py                   # reproduce; raw plates if present, ~15 s
python3 scripts/plate_readings.py --source committed  # force the committed text
python3 scripts/plate_readings.py --window          # also sweep the late window, ~4 min
python3 scripts/plate_readings.py --write           # also write the outputs/ tables
python3 scripts/plate_readings.py --export          # rebuild data/plates (needs YSTWIN_PLATES)
python3 scripts/plate_readings.py --verify          # check the text against the workbooks
```

`--window` is opt-in because it re-runs the whole pipeline at six late-window fractions,
and the cluster bootstrap behind each fold — not the reading — is what makes it minutes
rather than seconds. It is worth running once: where the late window starts
is the largest unforced choice in the dose response, and the sweep is what separates "the
folds are near 1.0" from "the folds are near 1.0 for one arbitrary window".

The script does not reimplement the analysis. It imports `run_sensor_characterisation.py`
and swaps the two module-level readers `collect()` reaches the workbooks through, so the
plate handling, the blank logic, the layout recovery and the ODE inversion are the same code
in both modes. Which source was used is printed above the table and again below it, because
a fallback that reads like a fresh derivation is worse than no fallback at all.

`--export` refuses to leave a mismatched export in place: if any block fails to round-trip
it exits non-zero and tells you not to commit it.

## 8. Calibrated sensor-to-GEM product challenge

After the editable development install, run:

```bash
python3 scripts/run_hybrid_benchmark.py --output-dir outputs/hybrid_benchmark/my_run --workers 4
```

The command fits only the selected sensor-training plates, withholds one complete plate,
generates training features through the reporter ODE, and couples inferred growth retention
to the frozen enzyme-constrained yeast model. The product panel is squalene, glycogen,
trehalose, glutathione, and glycerol. The same native enzyme-prior preparation is used for
all products. Product quantities are not fitting inputs.

The default scenarios vary carbon source, oxygen availability, temperature, DTT, and H2O2.
Their initial conditions and allocation policy are explicit in `protocol.json`. These are
allocation-conditioned model scenarios, not validated forecasts of realized production;
the policy is not adjusted to improve the external scores. Intracellular outputs are extra
net accumulation beyond the model's existing biomass quotas, and glycogen uses the frozen
model's glucose-equivalent molecular basis.

The output directory must be new or empty. The runner refuses to overwrite an existing run
and does not use `YSTWIN_OUTPUTS` to find its sensor inputs. Use `--sensor-table` explicitly
to choose another sensor-characterisation table. Worker processes have independent models;
only the parent writes artifacts, in deterministic request order. `--workers 1` runs serially.

Artifacts include:

- `protocol.json`, `sensor_model.json`, and `enzyme_inputs.csv`: inputs, hashes, calibration,
  held-out identity, numerical method, and modelling assumptions.
- `scenarios.csv` and `case_provenance.json`: rates, titres, consumed-substrate yields, and
  the constraints used to obtain them.
- `ablations.csv` and `sensitivities.csv`: no-learned-stress and allocation-policy checks.
  These are sensitivity experiments, not confidence intervals.
- `external_predictions.csv`, `external_scores.csv`, and `external_summary.csv`: independent
  targets scored without refitting, including unsupported, infeasible, invalid, and censored
  observations rather than dropping them.
- `external_baseline_scores.csv`, `external_baseline_summary.csv`, and
  `external_baseline_comparison.csv`: declared null/model-prior baselines, with comparisons
  made on the same scored rows.
- `validation_coverage.csv` and `product_scenarios.png`: measured coverage and the scenario plot.

External scores condition on measured physiology where available. They do not validate a
complete reporter-to-product chain: no benchmark source measures both on the same cultures.
Unknown detection limits are not zeros, unmeasured oxygen is not silently supplied, and
unresolved measurement bases are not declared comparable. Primary-source qualifications
are recorded in `data/hybrid_benchmark/sources.json`.

A numerical refinement can be reproduced without changing biological inputs:

```bash
python3 scripts/run_hybrid_benchmark.py --output-dir outputs/hybrid_benchmark/my_refinement \
  --workers 4 --environments glucose_aerobic glucose_DTT_1mM \
  --grid-points 13 --steps 480 --skip-sensitivity --skip-external --no-plots \
  --compare-to outputs/hybrid_benchmark/my_run
```

`numerical_refinement.csv` compares matching successful scenarios. The comparison refuses
changed data, sensor seeds, model implementations, environments, or allocation parameters;
it is not a way to present a changed model as a convergence check.

Focused verification:

```bash
python3 -m pytest -q tests/test_hybrid_stress.py tests/test_hybrid_benchmark.py \
  tests/test_hybrid_validation.py tests/test_product_panel.py tests/test_dynamic_batch.py \
  tests/test_hybrid_kernel_adversarial.py tests/test_uptake_cap.py
```

The full-model tests retain the physical flux acceptance tolerance. GLPK uses a stricter
local working precision to account for its internal scaling, records it in rate provenance,
and restores solver configuration afterwards; invalid physical solutions still refuse.

Assess a completed run without retraining, changing predictions, or overwriting its protocol:

```bash
python3 scripts/run_hybrid_benchmark.py --assess-run outputs/hybrid_benchmark/my_run \
  --numerical-run outputs/hybrid_benchmark/my_refinement \
  --assessment-output outputs/hybrid_benchmark/my_run/assessment.json
```

The assessment output must not already exist. Its execution census, numerical compatibility,
external coverage, and matched-row baseline comparisons are separate fields. `run_complete`
and a successful command exit mean the computation or report generation completed, **not**
that the model was scientifically validated. The assessment explicitly records
`validated_model: false`; neither successful simulations nor agreement with a simulator
establishes a measured reporter-to-product relationship.

## 10. In-silico teacher/student loop

This workflow learns from generated state, derivative, and control labels, not from
wet-lab sensor calibration or product amounts. It is separate from the retrospective
benchmark above. The virtual teacher's equations, coefficients, reporter response, and
noise model are recorded rather than presented as fitted biological facts.

The student has three supervised synthetic state coordinates (UPR, oxidative stress and
burden), two imposed normalized inputs, and four virtual reporter channels plus cell density.
The fixed teacher control coefficients are declared assumptions; the student fits separate
weights from their generated labels. The broader catalogue and optional mechanistic panel
routes are separate workflows, not additional states learned by this loop.

The state/control laboratory runs without a genome-scale model file:

```bash
python3 scripts/run_in_silico_loop.py --output-dir outputs/in_silico/learning_example \
  --skip-products
```

The full product panel uses the same frozen EC model and native protein-abundance prior as
[section 8](#8-calibrated-sensor-to-gem-product-challenge):

```bash
python3 scripts/run_in_silico_loop.py --output-dir outputs/in_silico/production_example \
  --workers 4 --production-seed 65021 \
  --modes student_prefix_forecast zero_state_baseline
```

Defaults train on complete bilinear/random episodes, evaluate separate trajectories and
an untrained changed-equation family, and audit the same hidden trajectories under
independently varied observation growth. The production cases cover squalene, glycogen,
trehalose, glutathione, and glycerol in aerobic glucose, oxygen-limited glucose, and aerobic
ethanol. The production seed and waveform are separate from the training seed. The student
receives only the first two hours of simulated culture observations and the known imposed
input schedule, then forecasts the remaining six hours. Controller updates are explicitly
sample-and-hold, not an unstated continuous-control approximation.

Use `--products` and `--environments` for a subset, `--workers 1` for serial execution,
`--no-plots` for tables only, and `--skip-nuisance-audit` only when that separate audit is
intentionally omitted. Adding `oracle_state_learned_control` to `--modes` isolates the
control decoder from state estimation. Each output directory must be new or empty;
existing runs are not overwritten. Generated run directories remain local and gitignored.

Artifacts:

- `teacher.json`, `student.json`, and `shuffled_student.json`: declared teacher parameters,
  fitted matrices/scaling, training identities, and the independently scrambled-label control.
- `split_manifest.csv`: complete-trajectory train/validation/test assignments.
- `state_control_scores.csv`: state, control, and prefix-only forecast errors, including raw
  state bounds, constrained-initialization adjustments, and rejected forecasts.
- `observation_invariance.csv`: paired same-state/different-growth checks, with noiseless and
  noisy observations. Rejections remain explicit; forecast means are conditional on acceptance.
- `protocol.json` and `model_source.zip`: seed/configuration, units, model and implementation
  hashes, base commit, and a snapshot of the model implementation used by that run.
- `product_scores.csv`, `product_summary.csv`, `product_trajectories.csv`, and `states_*.csv`:
  teacher/student process quantities and their errors. Specific rates distinguish structural
  biomass from total dry biomass; yield uses actually consumed substrate.
- `product_provenance.json`: reference and inferred control schedules, structural-biomass
  conversion, metabolic constraints, and the origin of the prefix observations.
- `assessment.json`: predeclared recovery checks, complete product coverage, explicit
  forecast-rejection counts, and the changed-equation challenge kept separate.
- `state_forecast.png` and `product_forecasts.png`: the corresponding trajectory plots.

`run_complete` describes execution. `assess_teacher_recovery` computes
`teacher_recovery_passed` from the learning, observation-invariance and selected-product
checks. It is `null` when required checks are omitted and no known check fails, and false
when any known recovery check fails. `cross_equation_generalization_passed` is a separate
challenge using the **same** forecast-RMSE criterion as reference recovery; it does not turn
same-teacher agreement into independent validation. `biological_validation` remains false.
A new run returns a nonzero exit code for a known recovery failure or failed product job;
the separate cross-equation failure is reported but does not set that execution status.

Review a completed run without refitting or changing its predictions:

```bash
python3 scripts/run_in_silico_loop.py --assess-run outputs/in_silico/production_example \
  --assessment-output outputs/in_silico/production_example/review.json
```

The review file must not exist. Source artifact hashes, `refitted: false` and
`predictions_changed: false` are recorded; existing protocol, scores and assessment remain
unchanged. Review-mode exit zero means report generation succeeded, not that either scientific
check passed. The local `closed_loop_04` review at
`outputs/stress_map_audit_02/teacher_recovery_review.json` passed teacher recovery and failed
cross-equation generalization under the same forecast-RMSE criterion. The review records both
errors and the failed check; its source run retains the old assessment schema and is not rewritten.

The process is synthesis-only. Its structural biomass excludes the two reserve quotas,
and explicit intracellular tracked product contributes to total dry mass once. It neither
restores a hidden basal pool nor opens an untracked reserve-withdrawal source. The molecular
basis remains the GEM's own species definition. See architecture section 10 for the
reference-equivalence and measurement-invariance contracts.

Focused verification:

```bash
python3 -m pytest -q tests/test_in_silico_teacher.py tests/test_in_silico_learning.py \
  tests/test_storage_basis.py tests/test_in_silico_loop.py tests/test_dynamic_batch.py
```

## 11. Stress-map and teacher-assumption audit

```bash
python3 scripts/audit_stress_map.py --output-dir outputs/stress_map_audit_example
```

This reads both GEMs resolved by `YSTWIN_EC_YEAST_GEM` and `YSTWIN_YEAST_GEM`, the public
`data/gem/stress_pathway_coverage.tsv` catalogue and the conceptual
`data/gem/stress_response_map.json`. The destination must be new or empty. `--no-render`
writes only the tables/report; `--map` and `--coverage` select explicit alternative inputs.
No model optimization, fitting or regulatory parameterization is performed.

Outputs are `scope.json`, `gene_associations.csv`, `pathway_scope.csv`,
`teacher_assumptions.csv`, and the `yeast_stress_expanded` / `yeast_stress_architecture`
figures in SVG and PNG. The report records input hashes, actual model counts, queried ORFs
and their current reaction associations, rather than inferring presence from gene symbols
or titles. The assumption inventory classifies every `TeacherParameters` field, including
the fixed synthetic control law; moving those constants into a file would not make them
learned biology.

The older panel has 24 catalogue entries (19 transcriptional programmes and five pools),
25 stressors, and only two optional mechanistic reporter routes. The new learner's three
synthetic aggregates are a separate scope. The atlas's branch inventory is conceptual, with
explicit unparameterized extensions. Historical diagram layers and Yeast9 reaction IDs
remain reference metadata, not ontological or EC-model coverage claims. A queried GPR
association does not establish signaling or regulation, and the count of catalogue items
with associations is not a percentage of biological pathway coverage.

The existing local run `outputs/stress_map_audit_02` contains the rendered audit. Its EC
and Yeast9 counts and associations are specific to the hashed model files, not biological
coverage scores. `tests/test_stress_map_audit.py` checks the audit on toy models without a
full GEM production benchmark.

## 12. Published HOG model and measured-data fitting

Install the normal pinned project dependencies. `python-libsbml` and `xlrd` are direct
requirements, so the legacy source workbooks and the kinetic model are reproducible on a
clean installation. No network is needed for the default fit when `data/hog2013/` is present.

```bash
python3 -m pip install -e ".[dev]"
python3 scripts/fetch_hog_reference.py
python3 scripts/run_hog_learning.py --output-dir /tmp/ystwin-hog-example --gem-check --plot
```

The fetcher verifies existing artifacts and refuses changed hashes or identities instead of
overwriting them. Its default assets are the WT SBML, S2 measurements, S4 preliminary
Westerns and Text S1 methods. `--assets` can additionally acquire the corrected mutant models,
raw metabolite workbook and enzyme-assay workbook registered in the script. Sources are
Petelenz-Kurdziel et al. 2013 (PMID 23762021) and its 2014 mutant-caption correction.
The article's primary XML supplies the Creative Commons Attribution license; its exact
version, DOI/download links, attribution and original artifact hashes are in `sources.json`.

The destination must be new or empty. `--source-dir` selects an explicit source bundle;
`--max-nfev` bounds numerical optimization. `--gem-check` additionally requires the EC GEM
resolved by `paths.ec_yeast_gem()`; without it the measured-data fit and waveform holdouts
still run. `--plot` uses the project's matplotlib development dependency.

### Frozen inputs and outputs

`data/hog2013/fit_protocol.json` records parameterization, units, source cells, splits and
normalization before scoring. The loader retains 407 observations: 13 fitting rows, 38
stronger-dose holdout rows, 52 other-assay holdout rows and 304 unused/diagnostic rows. The
S4 workbook explicitly says its preliminary data were not used for the original parameter
fit; this is not independence from the original model-development process.

The default estimates a relative activation/deactivation balance with the fast deactivation
rate held fixed. A two-rate training-only diagnostic could move both rates together almost
fivefold with little error change; those minute-sampled data do not establish two separate
fast kinetic constants. The protein measurement scale is estimated only on training data.
The published 90%-Hog1-peak assumption is retained in measurement provenance, not presented
as a direct absolute measurement. S4's time column is unlabelled; its interpretation as
relative minutes is explicit and recorded.

The run writes `protocol.json`, `observations.csv`, `fit.json`, `heldout_scores.csv`,
`heldout_predictions.csv`, `model_native_trajectories.csv`, `metabolic_signals.csv`, and
`validation.json`. Optional artifacts are `gem_capacity_check.csv` and `hog_learning.png`.
The original SBML variables and source normalization remain distinguishable from measured
values and inferred controls. A failed run does not write a false completion summary; use a
fresh destination for the next attempt.

The completed local run is `outputs/hog_learning_run_03`. Its `fit.json` gives the fitted
activation-balance multiplier; `validation.json` gives the measurement-block mean waveform
errors, and `heldout_scores.csv` preserves the per-trace scores. The refit does not improve
the stronger-dose prediction over the published parameters. These machine-generated artifacts
are the numerical record; results are not manually duplicated here.

The GEM check increases the GPD1 enzyme bound while preserving its independent hard cap and
shared protein pool, but glycerol capacity remains unchanged. `gem_capacity_check.csv` records
the enzyme bounds, growth, uptake, product capacity, units and solver. This is an exact endpoint
capacity probe under an explicitly chosen growth/carbon-budget scenario, not realized osmolyte
flux, measured YPD growth or validated titre prediction. An intracellular glycerol sink is
distinct from the default extracellular secretion task. A calibrated retention-demand,
permeability and allocation law is still missing.

### Verification

```bash
python3 -m pytest -q tests/test_fetch_hog_reference.py tests/test_kinetic_sbml.py \
  tests/test_hog_data.py tests/test_hog_learning.py tests/test_hog_mechanism.py \
  tests/test_hog_learning_runner.py tests/test_product_panel.py
```

The BDF/Radau full-source comparison runs with the normal dependencies. An optional independent
reaction-rate check uses RoadRunner, pinned in the `sbml-audit` extra; it is not a runtime
requirement and otherwise skips. Install the extra in your verification environment. This
command also runs the BDF/Radau comparison:

```bash
python3 -m pip install -e ".[dev,sbml-audit]"
python3 -m pytest -q tests/test_kinetic_sbml.py -k 'independent_roadrunner or implicit_solvers'
```

All 58 reaction rates are compared at seven source-trajectory states with the reference
engine's clock explicitly set. This is separate from biological validation. Strict-tolerance
RoadRunner CVODE integration failed during the source trajectory, so no successful full
RoadRunner trajectory comparison is claimed. The tested implicit-solver trajectory comparison
is BDF versus Radau on the unchanged source equations.

## 13. Native training before an unseen-product test

`data/hog2013/native_training_protocol.json` fixes the native selections, parameter groups,
obligation summaries and prediction interpretation. Training uses the native HOG observations
and the verified native chemostat table, not product measurements or product-specific kinetic
calibrations. The HOG source and host GEM remain explicitly identified priors.

Choose a fresh or empty directory **inside the repository** and set `NEW_NATIVE_RUN` to it:

```bash
python3 scripts/run_native_training.py --output-dir "$NEW_NATIVE_RUN"
```

The runner writes the source/code contract, input manifest, selected native observations,
native fits, native-condition diagnostics, derived obligations, project-read audit,
`checkpoint.json` and `checkpoint_receipt.json`. Keep the printed expected checkpoint digest
separately before introducing any task. The earlier provisional native run was superseded
following a neutral pre-target review; the corrected local checkpoint is in
`outputs/native_training_run_02`.

The guard permits only the explicitly registered project input/code files and the new run's
artifacts. Interpreter libraries are identified from the interpreter configuration. This is
scoped audit evidence, not an OS sandbox. The final read audit records which native files and
model modules entered the actual fit. Native control data used in fitting are not subsequently
advertised as an independent holdout.

### Task preparation, prediction, then scoring

Only after a verified freeze may a curator introduce independently sourced chemistry and
experimental conditions. The chemistry schema is `ChemicalTask`; it has no titre or fitted
allocation field. The condition bundle has explicit experiment IDs, chemical-file digests,
reported conditions and source metadata. Product outcomes belong in a separate label file.
Do not substitute a gene-copy count for measured catalytic capacity. A measured intermediate
also requires different interpretation from a terminal product.

The generic coordinator validates the checkpoint, its readout declaration and all loaded
model/evaluator identities before opening the condition bundle. Set the variables below to
the retained artifact paths/digests and to new output directories inside the repository:

```bash
python3 scripts/evaluate_frozen_transfer.py predict \
  --contract "$NATIVE_CONTRACT" --checkpoint "$NATIVE_CHECKPOINT" \
  --checkpoint-sha256 "$NATIVE_CHECKPOINT_SHA256" \
  --conditions "$TASK_CONDITIONS" --output-dir "$NEW_PREDICTION_RUN" --workers 2
```

This writes the task manifest, sealed interval predictions, scenario details and an evaluation
receipt. Retain the printed prediction digest. No outcome labels are arguments to prediction.
The frozen local hypothesis runs sealed their numerical predictions before an auxiliary
`prediction_details.json` export failed on an unbounded native flux limit. Those partial detail
files are not complete records. Receipts were recovered only after validating the sealed
artifacts against their retained digests; model, inputs, code and prediction values were not
changed. This reporting defect remains in that frozen revision rather than being silently
patched during the experiment.

The native physiological curves may report censored or unsupported support; fixed source
hard bounds may also make a condition infeasible. Neither failure authorizes an automatic
relaxation or a product-specific recalibration.

Only after predictions are sealed should the separately approved label files be prepared or
opened for scoring:

```bash
python3 scripts/evaluate_frozen_transfer.py score \
  --contract "$NATIVE_CONTRACT" --checkpoint "$NATIVE_CHECKPOINT" \
  --checkpoint-sha256 "$NATIVE_CHECKPOINT_SHA256" \
  --receipt "$EVALUATION_RECEIPT" --predictions-sha256 "$PREDICTIONS_SHA256" \
  --label-manifest "$HELDOUT_LABEL_MANIFEST" --output-dir "$NEW_SCORE_RUN"
```

Scoring verifies disjoint source identities, matching quantity/unit, exact experiment IDs,
unchanged model/code/task hashes and the sealed prediction bytes. It retains failures and
reports both total and conditional coverage. Envelope consistency is not point-prediction
accuracy; do not quote its upper boundary as a predicted titre. A failed native or downstream
transfer is a result, not a reason to change the checkpoint after viewing outcomes.

The completed strict-versus-hypothesis assessment is summarized in
`outputs/heldout_transfer_summary.json`. Strict chemistry remains blocked; the three explicitly
unverified redox variants retain both successful envelopes and native-feasibility failures.
Their separate score artifacts are `outputs/heldout_score_nad/score.json`,
`outputs/heldout_score_oxygen_peroxide/score.json`, and
`outputs/heldout_score_oxygen_water/score.json`. The verified measurement projection and its
primary-source conversions are in `data/holdout_transfer/labels.json` and
`data/holdout_transfer/label_provenance.json`. No model or hypothesis was selected using those
scores. The curator-exposure qualification remains part of the result.

Focused checks, using only neutral fixtures/native data:

```bash
python3 -m pytest -q tests/test_native_hog_data.py tests/test_hog_learning.py \
  tests/test_native_physiology.py tests/test_hog_mechanism.py tests/test_native_obligations.py \
  tests/test_chemical_task.py tests/test_blind_transfer.py tests/test_training_freeze.py \
  tests/test_native_training_runner.py tests/test_frozen_transfer_runner.py
```

## 14. Native reconciliation and biology-learning readiness

These are new diagnostic/audit entry points, not modifications of the frozen native model.
They do not fit product outcomes or silently adopt new host parameters. Set the output variables
to new or empty directories inside the repository:

```bash
python3 scripts/run_native_reconciliation.py \
  --output-dir "$NEW_NATIVE_RECONCILIATION_RUN" --timeout-s 15

python3 scripts/audit_biology_learning.py \
  --output-dir "$NEW_BIOLOGY_AUDIT_RUN" \
  --require-claim conditional_parameter_estimation \
  --require-claim independent_mechanistic_evidence

python3 scripts/check_biology_readiness.py \
  --output-dir "$NEW_BIOLOGY_READINESS_RUN" \
  --native-mode printed_rounding --timeout-s 15
```

The audit/readiness commands deliberately return exit code **2** when a requested biological
claim is unsupported. Their saved reports distinguish that expected refusal from solver or
program failure. Do not remove the required claim or treat a unit-test pass as biological
validation. The readiness command re-solves the native conditions and revalidates authoritative
source evidence instead of trusting saved booleans.

The audit's `--reuse-native-diagnostics` option can reuse the pinned retrospective numerical
report without re-running fits. Its receipt records the source digest, and reuse grants no new
validation credit. Source/code validation and release authorization are still recomputed.

Current results:

- `outputs/biology_learning_review.json` combines the findings and remaining evidence needs.
- `outputs/native_reconciliation_20260907_01/findings.json` describes the native discrepancies;
  its `final` subdirectory holds the authoritative scenario run and full compressed audits.
- `outputs/biology_learning_audit_02_gate_hardening/ledger.json` records the verified source
  inventory and evidence grades. Use `require_claim(..., root=root)` for authorization, not
  serialized claim flags.
- `outputs/biology_readiness_01/readiness.json` records the combined readiness refusal.
- `outputs/biology_observability_gauge.json` records the original scale-symmetry probe;
  `tests/test_native_observability.py` also contains the subsequent pool-counterflow test.
- `data/biochemical_evidence/enzyme_activity_sources.json` records primary enzyme evidence,
  assay units and explicit non-transferable or unidentified quantities.
- `data/validation_candidates/` contains the initial source-lineage and external-data metadata.
  Only the later approved development projections in section 15 have been released; other
  groups remain reserved and no independent signaling validation has been performed.

Focused checks:

```bash
python3 -m pytest -q tests/test_native_reconciliation.py \
  tests/test_biology_learning_audit.py tests/test_native_observability.py \
  tests/test_biology_readiness.py tests/test_native_obligations.py \
  tests/test_native_physiology.py
```

Rounding-consistent native cases are not an independent biological validation. Source zeros
remain censored, native strain/medium and dry-mass contexts remain different, and shared source
fitting/selection history must be included in the next validation design. The already-scored
lycopene results cannot be reused as a fresh blind test for changes motivated by those results.

## 15. Source-cohort development, not a new validation claim

The original strong native-law protocol remains ineligible. The separate
`data/native_law_v2/development_protocol.json` fixes the published-cohort estimand, source fields,
prefix/response windows, whole-group assignments, empirical families and development selection
rule. Only its approved projections are released; mixed raw data and private custody records
remain outside the learner allowlist and outside the commit.

The source/prediction sequence is retained in:

- `data/native_law_v2/granados/development/release_manifest.json`
- `data/native_law_v2/granados/development/executor_refresh_01/executor_refresh_admission.json`
- `outputs/native_population_development/native_training_01/development_prediction_freeze.json`
- `data/native_law_v2/granados/development/phase_c_01/release_manifest.json`

Metadata admission can be checked without opening response labels:

```bash
python3 scripts/run_native_population_development.py status \
  --allowlist data/native_law_v2/granados/development/release_manifest.json \
  --executor-admission data/native_law_v2/granados/development/executor_refresh_01/executor_refresh_admission.json
```

Do not rerun fitting to improve a score already viewed. The published development comparison
and every group/frame error are in
`outputs/native_population_development/native_training_01_phase_c_01_scoring/development_report.json`.
The selected family retains its training-bound warning. The independent data review and score
recomputation preserve missingness, cohort definitions, all families and the predeclared metric.

The original scorer's reusable cached-admission defect was fixed separately. Preserved code
snapshots and unchanged numerical verification are recorded in the local original
`outputs/native_population_development/native_training_01_phase_c_01_scoring_fix_v2/reusable_interface_fix_verification.json`.
Its path-neutral public envelope is included in the section 17 bundle; the original contains a
machine-specific authority root and is deliberately not committed. The current scorer can verify
the published result without fitting or new selection:

```bash
python3 scripts/score_native_population_development.py verify-published
python3 -m pytest -q tests/test_native_population_development.py \
  tests/test_native_population_adversarial.py tests/test_native_population_scoring.py \
  tests/test_native_population_scoring_adversarial.py
```

This remains retrospective processed-cohort development. It is not a certified prospective
forecast, a new independent biological-law test, or permission to expose the remaining groups.

## 16. Partial-order models and consolidated utilities

A partial order carries supported comparisons, not invented magnitudes or importance. The
new module keeps context, calibration, quantity, relation kind and source assumptions explicit.
Choose an existing output directory with no conflicting result files for each demonstration:

```bash
python3 scripts/run_order_robustness.py --output-dir "$NEW_ORDER_RUN"
python3 scripts/run_native_order_robustness.py --output-dir "$NEW_NATIVE_ORDER_RUN"
python3 -m pytest -q tests/test_partial_orders.py tests/test_partial_orders_adversarial.py \
  tests/test_native_order_robustness.py
```

The first command is synthetic engineering evidence only. The native demonstration uses the
verified native chemostat adapter; results are in `outputs/native_order_summary.csv`, with full
conditional certificates and sampling assumptions in `outputs/native_order_robustness.json`.
Order-only absolute magnitudes can remain unbounded. The printed-precision bands are not
statistical confidence intervals, and nonlinear samples are not global certificates.

The reviewed branch utility exports encoded turnover priors and a declared capacity sweep:

```bash
python3 scripts/ecmodel_isoprenoid_prior.py --output-dir "$NEW_PRIOR_EXPORT"
python3 -m pytest -q tests/test_ecmodel_isoprenoid_prior.py tests/test_published_content_query.py
```

No turnover prior is relabeled as a measured enzyme property. Every protein coefficient remains
visible; reverse arms stay in the inventory but do not set the forward sweep. Optional content
comparators are source-labelled comparisons, never fitting targets. The current published-cassette
reader's `rows_above_content` query uses the structured measurand field rather than caption
keywords or another product's volumetric titre.

## 17. Path-neutral public evidence and replay

Some exact historical checkpoint, transfer and audit records contain machine-specific storage
paths. Their original bytes remain local and unchanged rather than being committed or rewritten.
The public counterpart is `data/frozen_evidence/native_v1/manifest.json`. The replay entrypoint
pins that release independently and verifies its identity before using it.

The bundle has distinct original and public identities. Its declared storage-path transforms do
not change scientific values. Unknown active dependencies cannot fall through to lineage-only
status; phase schemas and the existing role/source contracts determine their treatment. Boundary
references do not imply that private source bytes or historical execution were verified.

For a fresh, non-existing output file in an existing directory:

```bash
python3 scripts/replay_portable_evidence.py --output "$NEW_PUBLIC_REPLAY_JSON"
python3 -m pytest -q tests/test_portable_evidence.py tests/test_portable_evidence_adversarial.py \
  tests/test_portable_replay.py
```

The default replay requires no private-original directory. It checks the bound native inputs,
code and runtime; reproduces fixed HOG predictions and native obligations; and recomputes recorded
exchange, transfer and population metrics, including failure rows and original denominators.
`outputs/consolidation_checks/public_replay_final.json` retains the verified result.

Original-byte integrity is explicitly `not_checked` in public mode. Only an explicit
`--private-original-root` with the actual retained files permits that additional verification.
No new fitting, selection, transfer LP evaluation, chronology proof or biological-readiness claim
is made. Refusal on a changed code/runtime binding is not permission to restamp the historical
record. The original audit mode in earlier sections still needs its local original evidence;
public native reconciliation can instead use `--artifact-source portable`.

The public exporter tests use complete standalone synthetic source graphs with explicit fixture
pins; they do not hide original-preservation or exporter checks behind missing-private-data skips.
The shipped release is write-once. Re-exporting original records is a separate operation requiring
the retained originals and a new destination, not an overwrite of the public release.

Full regression runs can use the standard serial command in the README. With the separately
available `pytest-xdist` test runner, the verified parallel invocation is:

```bash
python3 -m pytest tests -q -n 4 --dist=loadfile --no-loadscope-reorder
```

GEM test fixtures cache private parsed templates and provide independent consumers. Expensive
derived builders retain private caches; mutable models and views are copied for each test.
Higher-scoped fixtures use source factories rather than borrowing a function-scoped consumer.
This preserves original LP bounds and solver configuration instead of letting earlier tests'
mutations or solved bases define later inputs. `tests/test_model_fixture_isolation.py` checks
that isolation and the cache lifetimes; it does not assert that every possible solver stall has
been explained.

## 18. Where the NIS excess comes from

`outputs/nis_channel_summary.csv` reports aggregate mean NIS 45.3 in OD and 19.3 in RFU
against a target of 1, alongside `median_min_ess` 1.029 and `median_min_unique_ancestors`
1.0 over 319 well-channels. Those two facts together are the problem: an ensemble down to
one effective particle can inflate NIS on its own, so until the sweep below existed nothing
separated a badly declared measurement scale, a wrong process model, and an ensemble that
had simply collapsed.

```bash
python3 scripts/run_calibration_nis.py --decompose --output-dir /tmp/ystwin-nis-decomposition
```

Roughly 13 minutes on one core, all of it the arms — the baseline alone is 31 s. It reads
the committed `data/plates/` text, so it runs on a clone. The extra table is
`nis_error_decomposition.csv`; the ordinary five outputs are unchanged, and without
`--decompose` nothing extra is written. `--decompose-particles` and `--decompose-sigma`
choose the arms; the defaults are the sweep below, and the 500 in the particle list is the
baseline itself and is dropped rather than run twice.

Each arm moves exactly one declared quantity and holds the seed, wells, priors and
readings. A particle arm changes the ensemble size. A measurement arm scales `od_rel_sigma`
and `rfu_rel_sigma` together, which leaves the dynamics alone — `_priors_for` reads only
`k_deg` out of the process overrides, so initialization starts from the same state. Process
error is not an arm; it is what neither family removes.

### What the default sweep measured

| arm | OD mean NIS | OD excess removed | RFU mean NIS | RFU excess removed | median min ESS |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline, 500 particles | 45.31 | — | 19.32 | — | 1.03 |
| 2000 particles | 32.19 | 29.6 % | 8.30 | 60.2 % | 1.07 (×1.04) |
| 8000 particles | 22.90 | 50.6 % | 4.44 | 81.2 % | 1.31 (×1.27) |
| σ × 2 | 11.74 | 75.8 % | 5.57 | 75.1 % | 4.07 (×3.96) |
| σ × 4 | 1.91 | 98.0 % | 1.30 | 98.4 % | 30.49 (×29.6) |
| σ × 8 | 0.42 | 101.3 % | 0.33 | 103.7 % | 101.81 (×98.9) |

ESS is a property of the joint filter, so one column serves both channels. Shares are of
`baseline − 1`, not of the NIS, and they are unclipped: past 100 % the arm has overshot the
target rather than merely reached it.

### What that rules in, and what it does not

**Particle approximation is in, and more particles is not the repair.** Sixteen times the
ensemble removes half of OD's excess and leaves NIS at 22.9 — still an order of magnitude
above target — and the surviving ensemble is *still* collapsed: median minimum ESS 1.31 out
of 8000 draws, median minimum unique ancestors 1. The suspect this sweep existed to
adjudicate is real and is not an artifact of running at 500.

**A measurement-only account is arithmetically available and physically refused.**
Measurement variance is 68.8 % of OD's predictive variance, so closing the whole excess by
declaring more reader noise needs `od_rel_sigma` × 8.09 — 2 % becoming 16.2 %, against the
1.1 % `PLATE_READER_CV` this repository measured from the residual scatter of well A10.
The empirical arms agree with that arithmetic (σ × 8 lands at 0.42, just past 1), which
makes it a falsified explanation rather than an untested one.

**The two families are not separated by this design, and the table says so.** Both knobs
act on the same weight degeneracy: a looser likelihood keeps particles alive exactly as a
larger ensemble does. `min_ess_ratio` records it — the σ × 8 arm multiplies the surviving
ensemble by 98.9, which is why its 101.3 % is not the measurement model's alone. So the
residual left to process error is an elimination among three families of which two overlap,
and the honest reading is a bound, not a split.

None of this changes the verdict. Every well-channel stays `INCONCLUSIVE` at the declared
configuration, coverage stays 22.7 % against a nominal 95 %, and the decomposition says
where the excess is, never that the filter's variance is sound.
