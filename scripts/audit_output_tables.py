"""A committed table in ``outputs/`` must still come out of the script that wrote it.

COVERAGE IS PARTIAL AND IS ITSELF RATCHETED. Every committed table in ``outputs/`` now lands
in exactly one of four named places, and whatever is left over is COUNTED on every run rather
than asserted here, because a hand-written total in a docstring is the thing this gate exists
to distrust:

* :data:`REGISTRY` -- re-run and diffed on every invocation. Cheap, deterministic, no argv.
* :data:`DEFERRED` -- VERIFIED to reproduce, but kept out of the default sweep by a measured
  cost or by argv the no-argv contract cannot pass. ``--slow`` runs them. They carry their own
  ratchet, :data:`_DEFERRED_TABLE_CEILING`, so nothing can be parked here and forgotten.
* :data:`EXEMPT` -- cannot be regenerated on a bare checkout at all: an input that is not
  tracked, a write-once frozen run directory, or a table assembled out of other scripts'
  outputs through the very directory the override redirects.
* :data:`DRIFT` -- measured NOT to reproduce. Reported as a failing row, every run, until
  somebody regenerates the table or fixes the generator.

    python3 scripts/audit_output_tables.py [--quiet] [--only NAME] [--slow]

EXEMPT USED TO BE KEYED BY SCRIPT STEM, AND THAT SILENTLY DISABLED IT. The coverage check
compares TABLE names, so ``enrolled |= set(EXEMPT)`` matched nothing: six recorded exemptions
subtracted zero tables and every one of them was still counted as unenrolled. EXEMPT is keyed
by table path now -- the same currency the ratchet counts in -- which is why the measured
count moved when the dictionary was re-keyed rather than because coverage changed. Three of
the old reasons were also false and are deleted rather than softened: ``product_environment_
sweep`` and ``design_flux_experiment`` were exempted as "needs the genome-scale model" and
both reproduce in seconds with no model at all, and ``gem_environment_predictions`` claimed a
table was "impossible on a checkout without the model" when ``data/gem/yeast-GEM.xml.gz`` is
git-tracked. The LP cost is that one's real and sufficient blocker.

THE GAP THIS CLOSES, and it is the one that makes the rest of the audit mean something.
`scripts/audit_claims.py` binds a number in prose to a CELL in a committed table. That is the
whole verification chain, and it has an unguarded end: if the CELL itself has drifted from the
code that produced it, marking the prose against it verifies that two wrong numbers agree.

It is not hypothetical. ``outputs/proteomics_gain.csv`` did not reproduce from
`scripts/simulate_proteomics_gain.py` at its own documented default seed. The script is fully
seeded -- ``--seed`` defaults to 0 and every draw comes off ``np.random.default_rng(seed)`` --
so the committed table was produced by a code state that no longer exists. Every cell
differed, and at ``protein_assay_cv_log = 0.2`` the committed file said -0.7438 where the
generator says +0.5773: A SIGN FLIP, in the column that decides the script's own printed
conclusion ("it stops being worth it past CV 0.20" against 0.30). Nothing could catch it. No
prose cited the table, so no `audit:value` marker pointed at it; and
`scripts/audit_determinism.py` discovers seeded entry points inside ``ystwin`` rather than
checking committed artefacts, so a script that is perfectly deterministic today still passes
while its committed output is from last week.

HOW IT WORKS. For each registered ``(script, tables)`` pair the script is re-run with
``YSTWIN_OUTPUTS`` pointed at a temporary directory -- the override `ystwin.paths.outputs_dir`
documents for exactly this purpose -- and each regenerated table is compared numerically,
cell by cell, against the committed one. Shape, column names and non-numeric cells must match
exactly; numeric cells must agree to :data:`TOLERANCE`.

THE TIME BUDGET IS DECIDED, NOT ACCRETED. The default sweep is capped at roughly three
minutes of wall clock, which is what the enrolled set costs today; a generator whose measured
run exceeds ~30 s goes to DEFERRED instead, however cleanly it reproduces. That is why
``run_cross_family`` (144 s), ``run_training`` (80 s), ``run_power`` (71 s) and
``gem_environment_predictions`` (229 s) are not in REGISTRY even though all four reproduce.
``--slow`` is the periodic non-CI sweep for them and costs about 45 minutes, nearly all of it
``run_transfer``.

WHAT IT DOES NOT CHECK. That a table is CORRECT -- only that it is the table its generator
produces today. A wrong computation reproduces perfectly. The two audits compose: this one
says the cell came from the code, `audit_claims.py` says the prose came from the cell, and
the tests say the code does what it claims.

Emits ``outputs/audit_output_tables.csv`` and exits non-zero if any check fails.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import subprocess
import sys
import tempfile
from typing import NamedTuple

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

__all__ = [
    "DEFERRED",
    "DRIFT",
    "EXEMPT",
    "REGISTRY",
    "TOLERANCE",
    "check_pair",
    "compare_tables",
]

REPO = pathlib.Path(__file__).resolve().parents[1]

#: Relative tolerance on a numeric cell. Not zero: a table round-trips through CSV text and
#: an LP or an optimiser can land a last digit differently on a different BLAS. Tight enough
#: that a real change in a reported figure fails -- the proteomics sign flip this gate was
#: written for is a 1.3-unit move in a column whose values are tens.
TOLERANCE = 1e-9

#: ``{script stem: (table names it writes)}``. Every entry must be deterministic, runnable
#: from a bare checkout with NO argv -- no untracked input, no network -- and cheap: the whole
#: dictionary is a ~3-minute sweep and a generator measured above ~30 s belongs in DEFERRED.
#: Grow it: a generator that is in none of these four dictionaries is a table nothing
#: regenerates.
REGISTRY: dict[str, tuple[str, ...]] = {
    "batch_sufficiency": ("batch_sufficiency.csv",),
    "calibrate_glycogen": ("glycogen_calibration.csv",),
    "content_band_coverage": ("content_band_coverage.csv", "content_band_coverage_summary.csv"),
    "design_fedbatch_run": ("fedbatch_design.csv",),
    # Was EXEMPT as "needs the genome-scale model", which was simply false: it imports numpy,
    # pandas and `ystwin.pathway` only, and is a seeded permutation-power simulation.
    "design_flux_experiment": ("flux_experiment_design.csv",),
    "env_to_product": ("env_to_product.csv",),
    "environment_channel": ("environment_channel.csv", "environment_law_scores.csv"),
    "error_budget": ("error_budget.csv",),
    "estimator_accuracy": ("estimator_accuracy.csv",),
    # The four carotenoid_entry_* tables were byte-identical twins of fit_pathway_flux's, because
    # both scripts made the same score_product_validation call. Removed at source, 2026-09-11.
    "fit_carotenoid_kinetics": ("carotenoid_fit.csv", "carotenoid_heldout.csv",
                                "carotenoid_parsimony.csv", "carotenoid_strain_capacity.csv"),
    # Enrolled by fixing the generator rather than by exempting five tables: `fit_pathway_flux`
    # wrote nothing without --output-dir, so it now defaults to paths.outputs_dir() like its
    # siblings `predict_product` and `fit_carotenoid_kinetics`.
    "fit_pathway_flux": ("pathway_flux_ranking.csv", "pathway_flux_nested.csv",
                         "pathway_flux_scores.csv", "pathway_flux_comparisons.csv"),
    "flux_null_baseline": ("flux_null_baseline.csv", "flux_candidate_scores.csv"),
    "parameter_provenance": ("parameter_provenance.csv",),
    "predict_product": ("product_prediction.csv", "product_prediction_scores.csv",
                        "product_prediction_selection.csv",
                        "product_prediction_comparisons.csv"),
    # Was EXEMPT as "needs the genome-scale model": there is no cobra import and no
    # `paths.yeast_gem()` call in it, and it runs in two seconds against an absent GEM.
    "product_environment_sweep": ("product_environment_sweep.csv",
                                  "product_environment_coverage.csv"),
    "register_prediction": ("registered_prediction_D018.csv",),
    "run_calibration_nis": ("nis_channel_summary.csv", "nis_summary.csv",
                            "nis_innovations.csv", "nis_exclusions.csv"),
    # `run_d2`, `run_gates` and `run_calibration_nis` all run from a bare checkout because
    # `ystwin.plate.replay` swaps the workbook readers for the committed plate text, which is
    # tracked -- so they are not in EXEMPT beside `run_sensor_characterisation`, which has no
    # such fallback. One caveat, stated rather than discovered: on a machine that HAS the raw
    # Synergy workbooks the scripts prefer them, so these pairs check the replay path there
    # against a committed table and would fail if the two ever disagreed. That is a finding,
    # not a false alarm -- `replay`'s whole claim is that its text round-trips the workbooks.
    "run_d2": ("d2_20260701_ER_preliminary_(RAW)__mCitrine.csv",
               "d2_20260709_ER_stress_1st_(RAW)__mCitrine.csv"),
    "run_gates": ("g1_20260701_ER_preliminary_(RAW).csv",
                  "g1_20260709_ER_stress_1st_(RAW).csv",
                  "g1_20260722_ERandOxidativeStress_NewProtocol_ANALYSED.csv",
                  "g1_20260728_ERandOxidativeStress_NewProtocol_Replicate2.csv",
                  "g1_20260803_ERandoxidativestress_Replicate3.csv",
                  "g1_20260804_ER_Oxidative_Replicate4.csv",
                  "d2_g1passed_20260701_ER_preliminary_(RAW).csv",
                  "d2_g1passed_20260709_ER_stress_1st_(RAW).csv",
                  "d2_g1passed_20260722_ERandOxidativeStress_NewProtocol_ANALYSED.csv",
                  "d2_g1passed_20260728_ERandOxidativeStress_NewProtocol_Replicate2.csv",
                  "d2_g1passed_20260803_ERandoxidativestress_Replicate3.csv",
                  "d2_g1passed_20260804_ER_Oxidative_Replicate4.csv"),
    "run_scenarios": ("scenario_predictions.csv",),
    "score_phb_environment": ("phb_environment_score.csv",),
    "score_sensor_crosstalk": ("sensor_crosstalk.csv",),
    "score_targets": ("target_scoreboard.csv",),
    "simulate_proteomics_gain": ("proteomics_gain.csv",),
    "sweep_conditions": ("condition_sweep.csv",),
    "thiolase_threshold": ("thiolase_threshold.csv",),
}

#: The fastest registered scripts, for the smoke check in `tests/test_audit_output_tables.py`.
#: The suite runs a couple of pairs so the gate cannot silently stop working; the full sweep
#: belongs on the command line, where a fourteen-second fit is acceptable and in a unit test
#: it is not. Nothing from DEFERRED may ever appear here.
FAST = ("environment_channel", "score_phb_environment", "error_budget")


class Deferred(NamedTuple):
    """One generator that reproduces its tables but is not run by the default sweep.

    ``cost`` records the MEASURED reason it is deferred -- seconds on the clock, or the argv
    the no-argv contract cannot pass -- because "slow" and "needs a flag" are the two excuses
    that turn into permanent silence if they are never written down with a number.
    """

    script: str
    argv: tuple[str, ...]
    tables: tuple[str, ...]
    cost: str


#: ``{label: Deferred}``, run by ``--slow`` and by nothing else. EVERY ENTRY HERE HAS BEEN
#: SEEN TO REPRODUCE; that is what separates this dictionary from EXEMPT, which is for tables
#: this harness cannot regenerate at all. The label carries the argv so one script can appear
#: twice.
DEFERRED: dict[str, Deferred] = {
    "run_cross_family": Deferred(
        "run_cross_family", (), ("cross_family_transfer.csv", "cross_family_optimism.csv"),
        "measured 101-144 s; pure seeded simulation, no data files"),
    "run_training": Deferred(
        "run_training", (), ("training_recovery.csv", "training_transfer_build.csv",
                             "training_transfer_three_sensor.csv", "training_transfer_wide.csv"),
        "measured 55-81 s at the full default mode the committed run_mode column records"),
    "run_power": Deferred(
        "run_power", (), ("power_analysis.csv",), "measured 71 s"),
    "env_to_product_multi": Deferred(
        "env_to_product_multi", (), ("env_to_product_multi.csv", "stress_module_ledger.csv"),
        "measured 30-52 s over the tracked Yeast9 GSMM"),
    "stress_diversion_scan": Deferred(
        "stress_diversion_scan", (), ("stress_diversion.csv",),
        "measured 48 s; the GSMM it reads, data/gem/yeast-GEM.xml.gz, IS tracked"),
    "maintenance_scale": Deferred(
        "maintenance_scale", (), ("maintenance_scale.csv",),
        "measured 61 s of bisection over NGAM; the GSMM it reads IS tracked"),
    "gem_environment_predictions": Deferred(
        "gem_environment_predictions", (),
        ("gem_environment_predictions.csv", "gem_environment_scores.csv",
         "gem_environment_formulation.csv"),
        "measured 229 s over ~200 LPs. Cost is the whole blocker: the model is tracked"),
    "run_transfer": Deferred(
        "run_transfer", (), ("transfer_by_stressor.csv", "module_recovery.csv"),
        "measured 1232-1459 s -- the default run is a full design comparison and null sweep"),
    "run_heldout_score --target uncertainty": Deferred(
        "run_heldout_score", ("--target", "uncertainty"), ("fold_bootstrap_coverage.csv",),
        "argv-gated write (default --target reporter) plus a measured 69 s"),
    "run_heldout_score --target product": Deferred(
        "run_heldout_score", ("--target", "product"),
        ("product_matched_entry_scores.csv", "product_matched_entry_comparisons.csv"),
        "argv-gated write (default --target reporter); measured 21 s"),
    "gem_version_diff <candidate>": Deferred(
        "gem_version_diff", ("data/gem/yeast-GEM-9.1.1.xml.gz",),
        ("gem_version_diff.csv",),
        "argv-gated: the candidate model is a positional argument, and with none the parser "
        "exits 2. Measured 8 s and REPRODUCES -- both GEMs are tracked, so cost is not it"),
    "plate_readings --write --window": Deferred(
        "plate_readings", ("--write", "--window"),
        ("sensor_characterisation.csv", "autofluorescence_sensitivity.csv",
         "late_window_sensitivity.csv"),
        "argv-gated write, and the window sweep re-reads every plate six times: 204 s"),
}

#: Shared by the six tables whose duplicate writes were removed at source on 2026-09-11.
_DEDUPED = (
    "no longer produced: its writer emitted a byte-identical twin of a table another script "
    "owns, and that duplicate write was removed at source 2026-09-11. The committed copy is "
    "retained because data/artifact_registry.json receipts record the executions that DID write "
    "it, and editing execution.written_paths to erase that would falsify a real run")

#: ``{table path under outputs/: reason}``. NOT script stems -- see the docstring; keying this
#: by stem is how it came to subtract nothing. Reserved for tables this harness cannot
#: regenerate at all. Cost and argv are NOT reasons to be here; they are what DEFERRED is for.
#: An exemption is a recorded decision; silence is not.
EXEMPT: dict[str, str] = {
    "afl_circuit.csv":
        "needs the raw Gen5 .xpt instrument files, which are not redistributable and have no "
        "default path; measured exit 1, 'set YSTWIN_GEN5_XPT'",
    "autofluorescence_measured.csv":
        "needs the raw Gen5 .xpt instrument files; measured exit 1 on the missing "
        "20260807_BY4741_ER_Oxidative.xpt",
    "gain_linearity.csv":
        "needs the raw Gen5 .xpt instrument files; measured exit 1 in 0.3 s",
    "d1_capacity_sweep_ec.csv":
        "its writer is scripts/parked/run_d1.py, which the scripts/<stem>.py lookup cannot "
        "reach, and it is the ecModel FVA sweep audit_reproducibility.py names as the "
        "repository's own example of solver-dependent degenerate-LP vertices",
    "external_channel_response.csv":
        "reads the Gasch 2000 GEO series matrices under data/external/, which .gitignore "
        "excludes and no clone has: git ls-files data/external returns nothing",
    "external_conditions.csv":
        "reads the untracked, network-fetched data/external/ inputs (gitignored)",
    "external_module_map.csv":
        "reads the untracked, network-fetched data/external/ inputs (gitignored)",
    "external_module_response.csv":
        "reads the untracked, network-fetched data/external/ inputs (gitignored)",
    "external_nulls.csv":
        "reads the untracked, network-fetched data/external/ inputs (gitignored)",
    "external_redundancy.csv":
        "reads the untracked, network-fetched data/external/ inputs (gitignored)",
    "external_transfer_pairs.csv":
        "reads the untracked, network-fetched data/external/ inputs (gitignored)",
    "external_channel_response_independent.csv":
        "untracked network-fetched inputs, and the _independent suffix exists only under "
        "--independent",
    "external_conditions_independent.csv":
        "untracked network-fetched inputs, and the suffix exists only under --independent",
    "external_module_map_independent.csv":
        "untracked network-fetched inputs, and the suffix exists only under --independent",
    "external_module_response_independent.csv":
        "untracked network-fetched inputs, and the suffix exists only under --independent",
    "external_nulls_independent.csv":
        "untracked network-fetched inputs, and the suffix exists only under --independent",
    "external_redundancy_independent.csv":
        "untracked network-fetched inputs, and the suffix exists only under --independent",
    "external_transfer_pairs_independent.csv":
        "untracked network-fetched inputs, and the suffix exists only under --independent",
    "g4_verdicts.csv":
        "needs the Bio-Rad qPCR exports that live outside the repository; measured here as a "
        "DEGRADED table rather than a clean skip -- anchor_effect row 2 regenerated as nan",
    "g4_anchor_fold_change.csv":
        "two of three qPCR biological replicates live outside the repository and replicate 1 "
        "is under the untracked data/qpcr_repaired; a bare checkout raises SystemExit",
    "g4_reporter_response.csv":
        "the anchor half needs the untracked qPCR exports; with them absent the script exits "
        "1, and it passes on this machine only from an untracked data/qpcr_repaired",
    "g4_rt_minus_qc.csv":
        "needs the qPCR exports: measured 48 rows against the committed 144, because "
        "qpcr_files() silently drops the two replicates whose directories are absent",
    "heldout_scores.csv":
        "reads split_manifest.csv and sensor_characterisation.csv from the same outputs "
        "directory the override redirects; measured exit 2 on the missing manifest",
    "heldout_interpolation_by_dose.csv":
        "same chained input: --input-dir defaults to the redirected outputs directory, so the "
        "manifest it needs is never there; measured exit 2",
    "split_manifest.csv":
        "needs the wet-lab plate exports; the script REFUSES rather than shrink the manifest, "
        "measured exit 1, and the committed 2616-row manifest is a full-scope run",
    "split_manifest_summary.csv":
        "written by the same refused make_splits run: full scope requires the unvendored "
        "biosensor plate exports",
    "panel_calibration.csv":
        "reads outputs/sensor_characterisation.csv through the redirected outputs directory; "
        "measured exit 1 naming that missing input",
    "fold_multiplicity.csv":
        "reads outputs/sensor_characterisation.csv through the redirected outputs directory; "
        "measured FileNotFoundError",
    "product_scoreboard.csv":
        "assembles product_prediction.csv out of the redirected outputs directory rather than "
        "computing a table; measured FileNotFoundError",
    "figures_manifest.csv":
        "a manifest of figures assembled from other scripts' outputs, and it READS those from "
        "the redirected directory: measured exit 1, 'no G1 tables under <tmp>'",
    "audit_claims.csv":
        "the marker-to-cell audit over other scripts' committed tables; read-only by design "
        "without --output-dir, and enrolling it would make this gate depend on its verdict",
    "audit_determinism.csv":
        "an audit that assembles other code's behaviour rather than computing a table, and it "
        "writes nothing without --output-dir",
    "audit_reproducibility.csv":
        "the audit that assembles every other script's artifacts; write-once, --output-dir "
        "only, and it rejects any output directory inside the repository",
    "ecmodel_prior_port_01/ecmodel_isoprenoid_kcats.csv":
        "path mismatch with a frozen bundle: the generator writes this name FLAT into the "
        "outputs root, while the committed copy is sha256-pinned inside the port directory",
    "ecmodel_prior_port_01/enzyme_capacity_sweep.csv":
        "same frozen-bundle path mismatch; the flat regeneration is cell-identical to the "
        "committed port copy, so this is a keying gap and not drift",
    "native_order_observations.csv":
        "a write-once frozen run: --output-dir ignores YSTWIN_OUTPUTS, is rejected unless it "
        "sits inside outputs/, and the script refuses to overwrite an existing result",
    "native_order_summary.csv":
        "same write-once frozen run; measured exit 2 with no argv, and the paths are opened "
        "with mode 'x'",
    "native_reconciliation_20260907_01/final/case_summary.csv":
        "a historical frozen run directory shipped with completion.json and verification.json; "
        "--output-dir is required and case_summary.csv is opened with mode 'x'",
    "native_training_run_02/native_training_observations.csv":
        "a frozen, hash-sealed training run: --output-dir is required, and the directory's "
        "checkpoint_receipt.json pins a sha256 that regenerating would break",
    "partial_order_v2/partial_order_synthetic_summary.csv":
        "a frozen re-run whose own run_context.json records the exact command that made it; "
        "the script also refuses to write into a directory that already holds a result",
    # The six below shared one git blob each with a table that is still produced: four scripts made
    # the same score_product_validation call and dumped it under their own names. Fixed at source.
    "pathway_flux_selection.csv": _DEDUPED,
    "content_band_selection.csv": _DEDUPED,
    "carotenoid_entry_selection.csv": _DEDUPED,
    "carotenoid_entry_comparisons.csv": _DEDUPED,
    "carotenoid_entry_heldout.csv": _DEDUPED,
    "carotenoid_entry_scores.csv": _DEDUPED,
    "regulation_counterfactual.csv":
        "--output-dir is required and must be empty, and the run loads the GECKO and plain "
        "GEMs and solves FVA over both; never verified, minutes of LP",
    "regulation_depth_sweep.csv":
        "--output-dir is required and must be empty; the sweep is repeated FVA per scale over "
        "the GECKO model. Never verified, minutes of LP",
    "regulation_module_reach.csv":
        "--output-dir is required and must be empty; measured exit 2 with no argv, before any "
        "of the FVA over the regulation layer",
    "regulation_product_range.csv":
        "--output-dir is required and must be empty; flux_variability_analysis over 28 "
        "conditions on two genome-scale models. Never verified, minutes of LP",
    "state_allocation/synthetic_demo/state_allocation.csv":
        "a frozen demo run in its own subdirectory: the script requires --output-dir AND one "
        "of --config/--synthetic-demo, and has no route to that subpath",
    "native_reconciliation_20260907_01/case_summary.csv":
        "a dated experiment directory the producer will not write into twice: "
        "run_native_reconciliation.py opens case_summary.csv with mode 'x' and the run itself "
        "refuses first; measured ValueError, 'refusing to overwrite an existing experiment "
        "directory'. The run is investigation_pending and its isolated replay stopped before "
        "numerical execution (data/native_reconciliation_recovery.json)",
    "native_training_run_01/native_training_observations.csv":
        "run_native_training.py requires --output-dir and refuses a destination that is not "
        "empty; measured FileExistsError, 'use a new or empty native-training destination'. "
        "Reproducing it means a full training execution into a fresh directory, not a re-write "
        "of this one",
}

#: ``{table path: what was measured}``. A table that does NOT reproduce. This is the state
#: the gate exists to make loud, so each entry is reported as a FAILING row on every run --
#: not filed as an exemption, which would hide the one confirmed instance of the defect.
#: Empty since 2026-09-10: the one entry was the superseded root copy of
#: partial_order_synthetic_summary.csv, removed because partial_order_v2 is a strict superset.
DRIFT: dict[str, str] = {}


def compare_tables(committed: pathlib.Path, regenerated: pathlib.Path) -> str:
    """``""`` when the two tables agree, otherwise the first disagreement, in words.

    Returns the reason rather than raising, because a mismatch has to reach the report as a
    failing row -- dropping it would leave the table unguarded while the audit said every
    check passed, which is the failure mode this whole script exists to remove.
    """
    if not regenerated.exists():
        return f"the script did not write {regenerated.name}"
    if not committed.exists():
        return f"outputs/{committed.name} is not committed; the script writes it"
    old = pd.read_csv(committed)
    new = pd.read_csv(regenerated)
    if list(old.columns) != list(new.columns):
        return (f"columns changed: committed {list(old.columns)} against "
                f"regenerated {list(new.columns)}")
    if len(old) != len(new):
        return f"row count changed: committed {len(old)} against regenerated {len(new)}"
    for column in old.columns:
        left, right = old[column], new[column]
        # Booleans are compared as TEXT, not as numbers. `True - True` raises in numpy, and a
        # boolean column is a verdict rather than a measurement: it has no tolerance.
        numeric = (pd.api.types.is_numeric_dtype(left)
                   and pd.api.types.is_numeric_dtype(right)
                   and not pd.api.types.is_bool_dtype(left)
                   and not pd.api.types.is_bool_dtype(right))
        if numeric:
            close = ((left - right).abs() <= TOLERANCE * right.abs().clip(lower=1.0))
            close |= left.isna() & right.isna()
            if not bool(close.all()):
                row = int((~close).idxmax())
                return (f"column {column!r} row {row}: committed {left.iloc[row]!r}, "
                        f"regenerated {right.iloc[row]!r}")
        elif not left.astype(str).equals(right.astype(str)):
            row = int((left.astype(str) != right.astype(str)).idxmax())
            return (f"column {column!r} row {row}: committed {left.iloc[row]!r}, "
                    f"regenerated {right.iloc[row]!r}")
    return ""


def check_pair(script: str, tables: tuple[str, ...],
               argv: tuple[str, ...] = ()) -> list[dict]:
    """Re-run one script into a temporary directory and diff every table it owns.

    ``argv`` is empty for everything in REGISTRY -- a table whose generator needs a flag to
    write it is a table the no-argv contract does not cover, and that is a DEFERRED entry
    whose flags are recorded beside it rather than an undeclared default.
    """
    path = REPO / "scripts" / f"{script}.py"
    label = script + (" " + " ".join(argv) if argv else "")
    if not path.exists():
        return [{"script": label, "table": table, "ok": False,
                 "detail": f"scripts/{script}.py does not exist"} for table in tables]
    with tempfile.TemporaryDirectory() as directory:
        environment = dict(os.environ, YSTWIN_OUTPUTS=directory)
        completed = subprocess.run([sys.executable, str(path), *argv], env=environment,
                                   capture_output=True, text=True, cwd=REPO)
        if completed.returncode != 0:
            tail = " ".join((completed.stderr or "").split())[-300:]
            return [{"script": label, "table": table, "ok": False,
                     "detail": f"the script exited {completed.returncode}: {tail}"}
                    for table in tables]
        rows = []
        for table in tables:
            reason = compare_tables(REPO / "outputs" / table,
                                    pathlib.Path(directory) / table)
            rows.append({"script": label, "table": table, "ok": not reason,
                         "detail": reason or "reproduces"})
        return rows


#: Committed tables in ``outputs/`` that are in none of the four dictionaries.
#:
#: SET TO THE MEASURED COUNT, and it may only ever be lowered -- the same ratchet discipline
#: as `scripts/audit_claims.py`'s unmarked-claim ceilings, and for the same reason: a gate
#: whose coverage is a hand-written list is a gate whose coverage silently rots.
#:
#: It was 89 and it FAILED at 115, which is the behaviour to preserve: it was first written as
#: 84 -- a guess, made in the same edit that added the ratchet to stop guessed ceilings -- and
#: the ratchet failed on its own first run and said so. Both numbers here are read off a run
#: of this script, never typed from a count in somebody's head.
_UNENROLLED_TABLE_CEILING = 0

#: Tables parked in DEFERRED: verified, but checked only under ``--slow``. Ratcheted for the
#: same reason -- otherwise the unenrolled counter reads zero while reproducible tables go
#: unverified in every default run. It falls when a generator gets cheap enough to promote to
#: REGISTRY, or when an argv-gated write learns a default the way `fit_pathway_flux` just did.
#:
#: 23 is this new category's INITIAL value, read off a run of this script rather than counted
#: by hand -- the mistake the unenrolled ratchet already made once. Setting an initial value
#: is not a raise; every change after this one is, and the no-raise rule applies from here.
_DEFERRED_TABLE_CEILING = 23


def enrolled_tables() -> set[str]:
    """Every committed table this gate has an explicit decision about."""
    tables = {table for tables in REGISTRY.values() for table in tables}
    tables |= {table for entry in DEFERRED.values() for table in entry.tables}
    return tables | set(EXEMPT) | set(DRIFT)


def unenrolled_tables() -> list[str]:
    """Committed ``outputs/*.csv`` in neither REGISTRY, DEFERRED, EXEMPT nor DRIFT."""
    enrolled = enrolled_tables()
    tracked = subprocess.run(["git", "ls-files", "-z", "--", "outputs/*.csv"],
                             cwd=REPO, capture_output=True, text=True, check=True)
    return sorted(relative.removeprefix("outputs/") for relative in tracked.stdout.split("\0")
                  if relative and relative.removeprefix("outputs/") not in enrolled
                  and relative != "outputs/audit_output_tables.csv")


def check_coverage() -> list[dict]:
    """The two ratchets on the gate's own reach."""
    missing = unenrolled_tables()
    ok = len(missing) <= _UNENROLLED_TABLE_CEILING
    deferred = sum(len(entry.tables) for entry in DEFERRED.values())
    return [{
        "script": "(coverage)", "table": "", "ok": ok,
        "detail": (f"{len(missing)} committed tables enrolled in neither REGISTRY, DEFERRED, "
                   f"EXEMPT nor DRIFT (ceiling {_UNENROLLED_TABLE_CEILING})"
                   + ("" if ok else
                      f" -- ROSE. Enrol it or exempt it with a reason; never raise the "
                      f"ceiling. Newest unenrolled: {missing[:5]}")),
    }, {
        "script": "(deferred)", "table": "",
        "ok": deferred <= _DEFERRED_TABLE_CEILING,
        "detail": (f"{deferred} verified tables are checked only under --slow "
                   f"(ceiling {_DEFERRED_TABLE_CEILING})"
                   + ("" if deferred <= _DEFERRED_TABLE_CEILING else
                      " -- ROSE. Make the generator cheap or give it a default, then promote "
                      "it to REGISTRY; never raise the ceiling.")),
    }]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quiet", action="store_true", help="only print failures")
    parser.add_argument("--only", action="append", default=None,
                        help="script stem to check; repeatable. Default: all registered")
    parser.add_argument("--slow", action="store_true",
                        help="also re-run DEFERRED; about 45 minutes, mostly run_transfer")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        help="write diagnostics here; default verification is read-only")
    args = parser.parse_args(argv)

    wanted = REGISTRY if args.only is None else {
        name: tables for name, tables in REGISTRY.items() if name in set(args.only)}
    deferred = DEFERRED if args.only is None else {
        label: entry for label, entry in DEFERRED.items() if entry.script in set(args.only)}
    if args.only and not wanted and not deferred:
        print(f"no registered script matches {args.only}; have {sorted(REGISTRY)}",
              file=sys.stderr)
        return 2

    rows = []
    for script, tables in wanted.items():
        rows += check_pair(script, tables)
    for label, entry in deferred.items():
        if args.slow:
            rows += check_pair(entry.script, entry.tables, entry.argv)
        else:
            rows.append({"script": label, "table": "", "ok": True,
                         "detail": f"DEFERRED, not run without --slow: {entry.cost}"})
    for table, reason in sorted(EXEMPT.items()):
        rows.append({"script": table, "table": "", "ok": True,
                     "detail": f"EXEMPT: {reason}"})
    # A drift is a failing row every run, with no way to acknowledge it into silence: the only
    # exits are regenerating the table or fixing the generator.
    for table, finding in sorted(DRIFT.items()):
        rows.append({"script": "(drift)", "table": table, "ok": False,
                     "detail": f"DOES NOT REPRODUCE: {finding}"})
    if args.only is None:
        rows += check_coverage()

    frame = pd.DataFrame(rows)
    report = None
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        report = args.output_dir / "audit_output_tables.csv"
        frame.to_csv(report, index=False)

    failed = frame[~frame.ok]
    for row in frame.itertuples():
        if args.quiet and row.ok:
            continue
        mark = "ok  " if row.ok else "FAIL"
        label = f"{row.script}" + (f" -> {row.table}" if row.table else "")
        print(f"  [{mark}] {label}: {row.detail}")
    checked = int((frame.table != "").sum())
    destination = f" -> {report}" if report else " (read-only; no report written)"
    print(f"\n  {checked - len(failed)}/{checked} committed tables reproduce from their "
          f"generator{destination}")
    if len(failed):
        print("  A committed table is not what its script produces. Regenerate it and "
              "re-read every sentence that quotes it -- the numbers moved.")
    return 1 if len(failed) else 0


if __name__ == "__main__":
    raise SystemExit(main())
