from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import pathlib
import subprocess
import zipfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.analysis.in_silico import fit_student, score_student
from ystwin.fba.product_panel import ProductTask, install_product, prepare_panel_model, product_names
from ystwin.fba.solver import load_model
from ystwin.fba.storage import separate_storage_biomass
from ystwin.generator.in_silico import (
    CONTROL_NAMES,
    LATENT_NAMES,
    generate_episodes,
    observe_reporters,
    simulate_teacher,
)
from ystwin.in_silico import InSilicoEnvironment, simulate_controlled_product

_MODEL = None
_STUDENT = None
_OPTIONS = None


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _write_json(path, value):
    path.write_text(json.dumps(_plain(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def _worker_init(student, model_path, options):
    global _MODEL, _STUDENT, _OPTIONS
    model, _ = load_model(model_path)
    _MODEL, _ = prepare_panel_model(model)
    _STUDENT, _OPTIONS = student, options


def _interpolate(times, source_times, values):
    return np.column_stack([np.interp(times, source_times, values[:, column])
                            for column in range(values.shape[1])])


def _forecast_control_schedule(student, times, inputs, prefix_times, observations, control_times):
    forecast = student.forecast(times, inputs, prefix_times, observations)
    inferred = student.infer(prefix_times, observations)
    initial = student.initialize(inferred)
    future = times > prefix_times[-1]
    states = _interpolate(control_times, np.r_[prefix_times[-1], times[future]],
                          np.vstack([initial.state, forecast[future]]))
    initializations = []
    for index, hour in enumerate(control_times):
        if hour > prefix_times[-1]:
            continue
        available = int(np.searchsorted(prefix_times, hour, side="right"))
        if available == 0:
            raise ValueError("a control update cannot precede the observation history")
        initialized = student.initialize(inferred[:available])
        if not initialized.accepted:
            raise ValueError(f"incompatible state estimate at control update {hour:g} h")
        states[index] = initialized.state
        initializations.append({"control_time_h": float(hour),
                                "last_observation_h": float(prefix_times[available - 1]),
                                **initialized.diagnostics()})
    return forecast, student.predict_controls(states), initializations


def _environment(name, hours):
    base = InSilicoEnvironment(name, hours=hours)
    if name == "glucose_low_oxygen":
        return replace(base, oxygen_uptake=2.0)
    if name == "ethanol_aerobic":
        return replace(base, carbon_source="ethanol", substrate_g_per_l=15.343)
    if name != "glucose_aerobic":
        raise ValueError(f"unsupported environment {name}")
    return base


def _production_job(job):
    product_name, environment_name = job
    args, student = _OPTIONS, _STUDENT
    environment = _environment(environment_name, args.hours)
    installed, original_task = install_product(_MODEL, product_name)
    structural, basis = separate_storage_biomass(installed)
    task = ProductTask(
        original_task.name, original_task.reaction_id, original_task.metabolite_id,
        original_task.molar_mass_g_per_mol, original_task.carbon_atoms,
        "explicit tracked pool with storage-free structural biomass and separate total dry mass",
        ("Synthetic teacher/student comparison; no product-amount fitting",
         "Model molecular bases, including glucose-equivalent glycogen",
         "No remobilization in this synthesis-only control protocol"))
    teacher = generate_episodes(
        1, seed=args.production_seed, protocol=args.production_protocol,
        hours=args.hours, dt=args.dt, noise_cv=0.0)[0]
    times, inputs = teacher.times_h, teacher.inputs
    knots = np.linspace(0, args.hours, args.control_intervals + 1)
    true_controls = _interpolate(knots[:-1], times, teacher.controls)
    simulation_options = {
        "structural_mass_fraction": basis.structural_mass_fraction,
        "grid_points": args.grid_points, "output_dt": args.dt,
    }
    truth = simulate_controlled_product(structural, task, knots, true_controls,
                                         environment, **simulation_options)
    observation_latent = _interpolate(truth.times_h, times, teacher.latent)
    observations = observe_reporters(
        truth.times_h, observation_latent, truth.growth_rate_per_h,
        truth.structural_biomass_g_per_l, seed=args.seed + 9100, noise_cv=args.noise_cv)
    prefix = truth.times_h <= args.prefix_hours + 1e-12
    predicted_latent, control_schedule, initializations = _forecast_control_schedule(
        student, times, inputs, truth.times_h[prefix], observations[prefix], knots[:-1])
    predicted_controls = student.predict_controls(predicted_latent)
    schedules = {
        "student_prefix_forecast": control_schedule,
        "oracle_state_learned_control": _interpolate(
            knots[:-1], times, student.predict_controls(teacher.latent)),
        "zero_state_baseline": np.repeat(student.predict_controls(np.zeros((1, len(LATENT_NAMES)))),
                                          len(knots) - 1, axis=0),
    }
    teacher_summary = truth.summary()
    rows = [{"product": product_name, "environment": environment_name,
             "mode": "teacher", "status": "predicted", **teacher_summary}]
    trajectories = [{
        "product": product_name, "environment": environment_name, "mode": "teacher",
        "time_h": hour, "structural_biomass_g_l": truth.structural_biomass_g_per_l[index],
        "total_dry_biomass_g_l": truth.total_dry_biomass_g_per_l[index],
        "substrate_mmol_l": truth.substrate_mmol_per_l[index],
        "product_mmol_l": truth.product_mmol_per_l[index],
        "product_molar_mass_g_per_mol": task.molar_mass_g_per_mol,
    } for index, hour in enumerate(truth.times_h)]
    for mode, schedule in schedules.items():
        if mode not in args.modes:
            continue
        prediction = simulate_controlled_product(structural, task, knots, schedule,
                                                  environment, **simulation_options)
        summary = prediction.summary()
        row = {"product": product_name, "environment": environment_name,
               "mode": mode, "status": "predicted", **summary}
        for metric in ("titre_mg_per_l", "productivity_mg_per_l_h", "yield_g_per_g",
                       "mean_specific_rate_mmol_per_g_struct_h", "mean_specific_rate_mmol_per_g_total_dw_h",
                       "content_mg_per_g_total_dw"):
            actual, said = teacher_summary[metric], summary[metric]
            row[f"{metric}_teacher"] = actual
            row[f"{metric}_absolute_error"] = abs(said - actual)
            row[f"{metric}_absolute_percent_error"] = abs(said - actual) / abs(actual) * 100 if actual else np.nan
        rows.append(row)
        for index, hour in enumerate(prediction.times_h):
            trajectories.append({
                "product": product_name, "environment": environment_name, "mode": mode,
                "time_h": hour, "structural_biomass_g_l": prediction.structural_biomass_g_per_l[index],
                "total_dry_biomass_g_l": prediction.total_dry_biomass_g_per_l[index],
                "substrate_mmol_l": prediction.substrate_mmol_per_l[index],
                "product_mmol_l": prediction.product_mmol_per_l[index],
                "product_molar_mass_g_per_mol": task.molar_mass_g_per_mol,
            })
    inference = pd.DataFrame({"time_h": times})
    for index, name in enumerate(LATENT_NAMES):
        inference[f"{name}_teacher"] = teacher.latent[:, index]
        inference[f"{name}_student"] = predicted_latent[:, index]
    for index, name in enumerate(CONTROL_NAMES):
        inference[f"{name}_teacher"] = teacher.controls[:, index]
        inference[f"{name}_student"] = predicted_controls[:, index]
    return {
        "rows": rows, "trajectories": trajectories, "inference": inference.to_dict("records"),
        "provenance": {
            "product": product_name, "environment": asdict(environment),
            "teacher_episode_id": teacher.episode_id,
            "stress_input_times_h": times.tolist(), "stress_inputs": inputs.tolist(),
            "structural_mass_fraction": basis.structural_mass_fraction,
            "storage_basis": basis.metadata,
            "controller_update_times_h": knots.tolist(), "true_controls": true_controls.tolist(),
            "student_controls": schedules["student_prefix_forecast"].tolist(),
            "control_state_initializations": initializations,
            "state_table_basis": "raw causal prefix estimates followed by bounded forecasts; executed controls use bounded causal state estimates",
            "teacher_stages": truth.stage_provenance,
            "forecast_prefix_h": args.prefix_hours,
            "observations": "reporter ODE driven by the same reference culture's sampled structural growth/density",
            "product_training_labels": False,
        },
    }


def _safe_production_job(job):
    try:
        return _production_job(job)
    except (ValueError, RuntimeError) as exc:
        return {"rows": [{"product": job[0], "environment": job[1], "mode": "failed",
                           "status": "failed", "reason": f"{type(exc).__name__}: {exc}"}],
                "trajectories": [], "inference": [], "provenance": {"error": str(exc)}}


def _nuisance_audit(student, args):
    episodes = generate_episodes(min(4, args.test_episodes), seed=args.seed + 67000,
                                  protocol="pulse", hours=args.hours, dt=args.dt, noise_cv=0.0)
    rows = []
    for episode_index, episode in enumerate(episodes):
        times = episode.times_h
        profiles = {
            "assay": episode.growth_rate_per_h,
            **{f"constant_{value:g}": np.full(len(times), value) for value in (0.05, 0.15, 0.30, 0.45)},
            "ramp_up": np.linspace(0.05, 0.45, len(times)),
            "ramp_down": np.linspace(0.45, 0.05, len(times)),
        }
        for profile, growth in profiles.items():
            exposure = np.r_[0, np.cumsum(np.diff(times) * (growth[:-1] + growth[1:]) / 2)]
            density = 0.03 * np.exp(exposure)
            for noise in sorted({0.0, args.noise_cv}):
                observed = simulate_teacher(
                    times, episode.inputs, episode_id=episode.episode_id,
                    seed=args.seed + 68000 + episode_index, noise_cv=noise,
                    growth_rate_per_h=growth, cell_density=density)
                if not np.array_equal(observed.latent, episode.latent):
                    raise RuntimeError("the nuisance audit must preserve the underlying latent trajectory")
                row = {"episode_id": episode.episode_id, "growth_profile": profile,
                       "noise_cv": noise, "status": "predicted", "forecast_rmse": float("nan")}
                try:
                    row.update(score_student(student, [observed], prefix_h=args.prefix_hours).iloc[0].to_dict())
                    if bool(row.get("forecast_rejected", False)):
                        row["status"] = "rejected"
                except ValueError as exc:
                    estimate = student.infer(times, observed.observations)
                    row.update(status="rejected", reason=str(exc),
                               state_rmse=float(np.sqrt(np.mean((estimate - episode.latent) ** 2))))
                rows.append(row)
    return pd.DataFrame(rows)


def _plots(out, student, episode, prefix_h, trajectories=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prefix = episode.times_h <= prefix_h
    forecast = student.forecast(episode.times_h, episode.inputs,
                                episode.times_h[prefix], episode.observations[prefix])
    figure, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    for index, (name, axis) in enumerate(zip(LATENT_NAMES, axes)):
        axis.plot(episode.times_h, episode.latent[:, index], label="teacher", color="#171717")
        axis.plot(episode.times_h, forecast[:, index], label="learned forecast", color="#2369b4")
        axis.axvline(prefix_h, color="#777777", linestyle="--", label="observation cutoff")
        axis.set_ylabel(name)
        axis.grid(alpha=0.15)
        axis.legend(fontsize=8)
    axes[-1].set_xlabel("Time (h)")
    figure.suptitle("In-silico state recovery and prefix-only forecast")
    figure.savefig(out / "state_forecast.png", dpi=150)
    plt.close(figure)
    if trajectories is None or trajectories.empty:
        return
    selected = trajectories[trajectories.environment == trajectories.environment.iloc[0]]
    products = list(selected["product"].unique())
    figure, axes = plt.subplots(len(products), 1, figsize=(11, 3 * len(products)),
                                squeeze=False, constrained_layout=True)
    for product, axis in zip(products, axes[:, 0]):
        for mode, group in selected[selected["product"] == product].groupby("mode", sort=False):
            axis.plot(group.time_h, group.product_mmol_l * group.product_molar_mass_g_per_mol,
                      label=mode.replace("_", " "))
        axis.set_title(product)
        axis.set_ylabel("Tracked product (mg/L)")
        axis.set_xlabel("Time (h)")
        axis.grid(alpha=0.15)
        axis.legend(fontsize=7)
    figure.suptitle("Product-label holdouts: teacher vs learned controls in the same in-silico world")
    figure.savefig(out / "product_forecasts.png", dpi=150)
    plt.close(figure)


def assess_teacher_recovery(state_scores, product_scores, nuisance_scores, expected_product_cases):
    reference = state_scores[state_scores.split == "test_protocol"]
    shifted = state_scores[state_scores.split == "test_family_shift"]
    checks = {
        "no_rejected_reference_forecasts": bool(len(reference) and not reference.forecast_rejected.any()),
        "reference_mean_forecast_rmse_below_0_01": bool(reference.forecast_rmse.mean() < 0.01),
        "reference_forecast_beats_persistence": bool(reference.forecast_rmse.mean() < reference.persistence_rmse.mean()),
        "growth_invariance_noiseless_rmse_below_0_03": None,
        "no_rejected_nuisance_prefixes": None,
        "complete_product_coverage": None,
        "median_product_errors_below_5pct": None,
        "worst_product_errors_below_10pct": None,
    }
    report = {
        "schema_version": 2,
        "scope": "same_teacher_surrogate_recovery",
        "biological_validation": False,
        "shared_teacher_assumptions": ["latent-state definitions", "virtual reporter law",
                                       "metabolic control law", "GEM and enzyme priors", "process integration"],
        "interpretation": "Small product errors quantify student approximation error while shared biological/model errors cancel; they do not validate the shared assumptions.",
        "reference_forecast_rmse": float(reference.forecast_rmse.mean()),
        "reference_forecasts_rejected": int(reference.forecast_rejected.sum()),
        "shifted_equation_forecast_rmse": float(shifted.forecast_rmse.mean()),
        "shifted_equation_forecasts_rejected": int(shifted.forecast_rejected.sum()),
        "shifted_equation_challenge": "reported separately; deliberately outside the trained dynamics family",
    }
    if nuisance_scores is not None:
        noiseless = nuisance_scores[nuisance_scores.noise_cv == 0]
        maximum = float(noiseless.groupby("growth_profile").state_rmse.mean().max())
        checks["growth_invariance_noiseless_rmse_below_0_03"] = maximum < 0.03
        checks["no_rejected_nuisance_prefixes"] = bool((nuisance_scores.status == "predicted").all())
        report["worst_noiseless_growth_profile_state_rmse"] = maximum
    if product_scores is not None:
        predicted = product_scores[product_scores["mode"] == "student_prefix_forecast"]
        teachers = product_scores[product_scores["mode"] == "teacher"]
        case_keys = ["product", "environment"]
        checks["complete_product_coverage"] = bool(
            len(predicted) == len(teachers) == expected_product_cases
            and len(predicted[case_keys].drop_duplicates()) == expected_product_cases
            and set(map(tuple, predicted[case_keys].to_numpy())) == set(map(tuple, teachers[case_keys].to_numpy()))
            and (product_scores.status == "predicted").all())
        error_columns = [f"{name}_absolute_percent_error" for name in
                         ("titre_mg_per_l", "mean_specific_rate_mmol_per_g_total_dw_h", "yield_g_per_g")]
        if len(predicted) and all(name in predicted for name in error_columns):
            errors = predicted[error_columns].to_numpy(dtype=float)
            complete = bool(np.isfinite(errors).all())
            medians = predicted[error_columns].median().to_numpy(dtype=float)
            worst = predicted[error_columns].max().to_numpy(dtype=float)
            checks["median_product_errors_below_5pct"] = complete and bool((medians < 5).all())
            checks["worst_product_errors_below_10pct"] = complete and bool((worst < 10).all())
            report["median_product_error_pct"] = dict(zip(error_columns, medians))
            report["worst_product_error_pct"] = dict(zip(error_columns, worst))
        else:
            checks["median_product_errors_below_5pct"] = False
            checks["worst_product_errors_below_10pct"] = False
    cross_equation = {
        "held_out_equation_rows_present": bool(len(shifted)),
        "no_rejected_forecasts": bool(len(shifted) and not shifted.forecast_rejected.any()),
        "mean_forecast_rmse_below_same_0_01_criterion": bool(shifted.forecast_rmse.mean() < 0.01),
        "forecast_beats_persistence": bool(shifted.forecast_rmse.mean() < shifted.persistence_rmse.mean()),
    }
    report["checks"] = checks
    known_failure = any(value is False for value in checks.values())
    report["teacher_recovery_passed"] = False if known_failure else (
        True if all(value is not None for value in checks.values()) else None)
    report["cross_equation_checks"] = cross_equation
    report["cross_equation_generalization_passed"] = all(cross_equation.values())
    return report


def reassess_run(run_dir, assessment_output):
    run = pathlib.Path(run_dir).resolve()
    output = pathlib.Path(assessment_output).resolve()
    if output.exists():
        raise ValueError(f"refusing to overwrite an existing assessment: {output}")
    protocol_path = run / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if protocol.get("mode") != "in_silico_teacher_student" or protocol.get("run_complete") is not True:
        raise ValueError("review requires a completed in-silico teacher/student run")
    scores_path = run / "state_control_scores.csv"
    products_path = run / "product_scores.csv"
    nuisance_path = run / "observation_invariance.csv"
    scores = pd.read_csv(scores_path)
    products = pd.read_csv(products_path) if products_path.exists() else None
    nuisance = pd.read_csv(nuisance_path) if nuisance_path.exists() else None
    report = assess_teacher_recovery(
        scores, products, nuisance, len(protocol["products"]) * len(protocol["environments"]))
    inputs = [path for path in (protocol_path, scores_path, products_path, nuisance_path) if path.exists()]
    report["source_run"] = paths.display_path(run)
    report["source_hashes"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    report["refitted"] = False
    report["predictions_changed"] = False
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_json(output, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description="Train and challenge an in-silico stress-state/controller student")
    parser.add_argument("--output-dir", type=pathlib.Path)
    parser.add_argument("--assess-run", type=pathlib.Path)
    parser.add_argument("--assessment-output", type=pathlib.Path)
    parser.add_argument("--train-episodes", type=int, default=64)
    parser.add_argument("--test-episodes", type=int, default=16)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--production-seed", type=int, default=9000)
    parser.add_argument("--production-protocol", choices=("random", "pulse", "ramp", "periodic"), default="pulse")
    parser.add_argument("--modes", nargs="+", choices=("student_prefix_forecast", "oracle_state_learned_control", "zero_state_baseline"),
                        default=["student_prefix_forecast", "oracle_state_learned_control", "zero_state_baseline"])
    parser.add_argument("--hours", type=float, default=8.0)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--noise-cv", type=float, default=0.001)
    parser.add_argument("--prefix-hours", type=float, default=2.0)
    parser.add_argument("--control-intervals", type=int, default=4)
    parser.add_argument("--grid-points", type=int, default=5)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--products", nargs="+", choices=product_names())
    parser.add_argument("--environments", nargs="+",
                        choices=("glucose_aerobic", "glucose_low_oxygen", "ethanol_aerobic"),
                        default=["glucose_aerobic", "glucose_low_oxygen", "ethanol_aerobic"])
    parser.add_argument("--skip-products", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--skip-nuisance-audit", action="store_true")
    args = parser.parse_args(argv)
    if args.assess_run is not None:
        if args.assessment_output is None or args.output_dir is not None:
            parser.error("--assess-run requires --assessment-output and no --output-dir")
        report = reassess_run(args.assess_run, args.assessment_output)
        print(json.dumps(_plain(report), indent=2, sort_keys=True, allow_nan=False))
        return 0
    if args.output_dir is None or args.assessment_output is not None:
        parser.error("a new run requires --output-dir; --assessment-output is only for --assess-run")
    if min(args.train_episodes, args.test_episodes, args.control_intervals, args.workers) < 1:
        parser.error("episode counts, control intervals and workers must be positive")
    if (not np.isfinite((args.hours, args.dt, args.prefix_hours, args.noise_cv)).all()
            or not 0 < args.prefix_hours < args.hours or args.dt <= 0
            or args.noise_cv < 0 or args.grid_points < 3 or min(args.seed, args.production_seed) < 0):
        parser.error("require finite 0 < prefix-hours < hours, dt > 0, nonnegative noise/seeds, grid-points >= 3")
    out = args.output_dir.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SystemExit(f"Refusing to overwrite an existing run: {out}")
    out.mkdir(parents=True, exist_ok=True)
    training = generate_episodes(args.train_episodes, seed=args.seed, hours=args.hours,
                                 dt=args.dt, noise_cv=args.noise_cv)
    validation = generate_episodes(args.test_episodes, seed=args.seed + 1000, hours=args.hours,
                                   dt=args.dt, noise_cv=args.noise_cv)
    tests = generate_episodes(args.test_episodes, seed=args.seed + 2000, protocol="pulse",
                              hours=args.hours, dt=args.dt, noise_cv=args.noise_cv)
    family_shift = generate_episodes(args.test_episodes, seed=args.seed + 3000, family="saturating",
                                     protocol="periodic", hours=args.hours, dt=args.dt,
                                     noise_cv=args.noise_cv)
    student = fit_student(training, seed=args.seed)
    _write_json(out / "student.json", student.metadata)
    _write_json(out / "teacher.json", training[0].metadata)
    manifest = []
    scores = []
    for split, episodes in (("train", training), ("validation", validation),
                             ("test_protocol", tests), ("test_family_shift", family_shift)):
        manifest.extend({"episode_id": episode.episode_id, "split": split} for episode in episodes)
        if split != "train":
            scored = score_student(student, episodes, prefix_h=args.prefix_hours)
            scored["split"] = split
            scores.append(scored)
    if len(training) > 1:
        rng = np.random.default_rng(args.seed + 17000)
        permutations = [rng.permutation(len(training)) for _ in range(3)]
        corrupted = tuple(replace(
            episode, latent=training[permutations[0][index]].latent,
            latent_derivative=training[permutations[1][index]].latent_derivative,
            controls=training[permutations[2][index]].controls)
            for index, episode in enumerate(training))
        negative = fit_student(corrupted, seed=args.seed)
        _write_json(out / "shuffled_student.json", negative.metadata)
        negative_scores = score_student(negative, validation, prefix_h=args.prefix_hours)
        negative_scores["split"] = "negative_control_shuffled_labels"
        scores.append(negative_scores)
    pd.DataFrame(manifest).to_csv(out / "split_manifest.csv", index=False)
    score_table = pd.concat(scores, ignore_index=True)
    score_table.to_csv(out / "state_control_scores.csv", index=False)
    print("State/controller evaluation on whole held-out episodes:", flush=True)
    print(score_table.groupby("split").agg(
        episodes=("episode_id", "size"), state_rmse=("state_rmse", "mean"),
        control_nrmse=("control_normalized_rmse", "mean"), forecast_rmse=("forecast_rmse", "mean"),
        persistence_rmse=("persistence_rmse", "mean"), rejected=("forecast_rejected", "sum"),
    ).to_string(), flush=True)
    protocol = {
        "schema_version": 2,
        "mode": "in_silico_teacher_student", "biological_validation": False,
        "evidence_scope": "same_teacher_surrogate_recovery",
        "teacher_model_scope": training[0].metadata["model_scope"],
        "base_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=paths.REPO_ROOT,
                                      capture_output=True, text=True, check=True).stdout.strip(),
        "run_complete": False, "wet_lab_calibration_used": False,
        "learning_targets": ["causal observations to latent state", "latent transition law",
                             "latent state to metabolic controls"],
        "not_learned": ["GEM stoichiometry", "product amount labels", "real biological validity"],
        "training_episodes": args.train_episodes, "seed": args.seed,
        "production_seed": args.production_seed, "production_protocol": args.production_protocol,
        "production_modes": args.modes,
        "nuisance_audit": not args.skip_nuisance_audit,
        "hours": args.hours, "dt": args.dt, "noise_cv": args.noise_cv,
        "prefix_hours": args.prefix_hours, "control_intervals": args.control_intervals,
        "controller_execution": "sample-and-hold at explicitly declared update times",
        "forecast_access": "prefix observations and known input schedule only; no future observations or labels",
        "model_source": "declared synthetic dynamics and control laws; independently supplied GEM/proteome priors",
        "products": list(args.products or product_names()), "environments": args.environments,
        "workers": args.workers, "grid_points": args.grid_points,
        "latent_names": list(LATENT_NAMES), "control_names": list(CONTROL_NAMES),
        "units": {
            "latent": "dimensionless bounded state coordinates",
            "stress_inputs": "normalized imposed inputs, not a fitted mM-to-state map",
            "carbon_uptake": "mmol carbon per reference gDCW per hour",
            "oxygen_uptake": "mmol oxygen per reference gDCW per hour",
            "ngam_control": "mmol ATP per reference gDCW per hour; converted to structural basis before solving",
            "initial_biomass": "g structural biomass/L; zero initial tracked product pool",
        },
        "code_hashes": {name: hashlib.sha256((paths.REPO_ROOT / name).read_bytes()).hexdigest()
                        for name in ("src/ystwin/generator/in_silico.py", "src/ystwin/analysis/in_silico.py",
                                     "src/ystwin/fba/storage.py", "src/ystwin/in_silico.py",
                                     "src/ystwin/fba/dynamic.py", "src/ystwin/fba/dynamic_rates.py",
                                     "src/ystwin/fba/product_panel.py", "src/ystwin/fba/physiology.py",
                                     "src/ystwin/fba/solver.py", "src/ystwin/reporter.py",
                                     "src/ystwin/growth.py", "src/ystwin/readings.py",
                                     "scripts/run_in_silico_loop.py")},
    }
    with zipfile.ZipFile(out / "model_source.zip", "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in protocol["code_hashes"]:
            archive.write(paths.REPO_ROOT / name, arcname=name)
    _write_json(out / "protocol.json", protocol)
    nuisance = product_table = None
    if not args.skip_nuisance_audit:
        nuisance = _nuisance_audit(student, args)
        nuisance.to_csv(out / "observation_invariance.csv", index=False)
        protocol["nuisance_rejected_prefixes"] = int((nuisance.status == "rejected").sum())
        print("Same-state observation invariance across growth conditions:", flush=True)
        print(nuisance.groupby(["growth_profile", "noise_cv"]).agg(
            state_rmse=("state_rmse", "mean"), forecast_rmse=("forecast_rmse", "mean"),
            rejected=("status", lambda values: int((values == "rejected").sum())),
        ).to_string(), flush=True)
    if not args.no_plots:
        _plots(out, student, tests[0], args.prefix_hours)
    if not args.skip_products:
        model_path = paths.require(paths.ec_yeast_gem(), "frozen EC model", "YSTWIN_EC_YEAST_GEM")
        protocol["ec_model_sha256"] = hashlib.sha256(model_path.read_bytes()).hexdigest()
        jobs = [(product, environment) for product in protocol["products"] for environment in args.environments]
        rows, trajectories, provenance = [], [], []
        if args.workers == 1:
            _worker_init(student, model_path, args)
            results = map(_safe_production_job, jobs)
            pool = None
        else:
            pool = ProcessPoolExecutor(max_workers=args.workers,
                                       mp_context=multiprocessing.get_context("spawn"),
                                       initializer=_worker_init, initargs=(student, model_path, args))
            results = pool.map(_safe_production_job, jobs, chunksize=1)
        try:
            for job, result in zip(jobs, results):
                rows.extend(result["rows"])
                trajectories.extend(result["trajectories"])
                provenance.append(result["provenance"])
                pd.DataFrame(rows).to_csv(out / "product_scores.csv", index=False)
                pd.DataFrame(result["inference"]).to_csv(out / f"states_{job[0]}_{job[1]}.csv", index=False)
                display = pd.DataFrame(result["rows"])
                columns = [name for name in ("product", "environment", "mode", "status", "titre_mg_per_l",
                           "productivity_mg_per_l_h", "yield_g_per_g", "titre_mg_per_l_absolute_percent_error",
                           "reason") if name in display]
                print(display[columns].to_string(index=False), flush=True)
        finally:
            if pool is not None:
                pool.shutdown()
        trajectory_table = pd.DataFrame(trajectories)
        trajectory_table.to_csv(out / "product_trajectories.csv", index=False)
        _write_json(out / "product_provenance.json", provenance)
        if not args.no_plots:
            _plots(out, student, tests[0], args.prefix_hours, trajectory_table)
        product_table = pd.DataFrame(rows)
        if "titre_mg_per_l_absolute_percent_error" in product_table:
            product_table.groupby("mode").agg(
                cases=("product", "size"),
                median_titre_error_pct=("titre_mg_per_l_absolute_percent_error", "median"),
                worst_titre_error_pct=("titre_mg_per_l_absolute_percent_error", "max"),
                median_yield_error_pct=("yield_g_per_g_absolute_percent_error", "median"),
            ).to_csv(out / "product_summary.csv")
        protocol["failed_product_jobs"] = sum(row["status"] == "failed" for row in rows)
    assessment = assess_teacher_recovery(
        score_table, product_table, nuisance, len(protocol["products"]) * len(args.environments))
    _write_json(out / "assessment.json", assessment)
    protocol["run_complete"] = True
    protocol["teacher_recovery_passed"] = assessment["teacher_recovery_passed"]
    protocol["cross_equation_generalization_passed"] = assessment["cross_equation_generalization_passed"]
    _write_json(out / "protocol.json", protocol)
    print(f"Teacher recovery: {assessment['teacher_recovery_passed']}; cross-equation generalization: {assessment['cross_equation_generalization_passed']}; biological validation: not tested", flush=True)
    print(f"In-silico results: {out}", flush=True)
    return int(protocol.get("failed_product_jobs", 0) > 0 or assessment["teacher_recovery_passed"] is False)


if __name__ == "__main__":
    raise SystemExit(main())
