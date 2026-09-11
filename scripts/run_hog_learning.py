from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.analysis.hog_data import load_hog_data
from ystwin.analysis.hog_learning import evaluate_hog_fit, fit_hog_parameters
from ystwin.mech.hog import HogModel, HogProtocol, hog_metabolic_signals


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _json(path, value):
    path.write_text(json.dumps(_plain(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def _summary(scores):
    metrics = ["fitted_shape_rmse", "published_shape_rmse", "persistent_shape_rmse"]
    blocks = scores.groupby(["split", "experiment_id"], sort=True)[metrics].mean()
    result = {}
    for split, group in blocks.groupby(level="split", sort=True):
        result[split] = {**group.mean().to_dict(), "measurement_blocks": len(group),
                         "trace_count": int(scores["split"].eq(split).sum()),
                         "observations": int(scores.loc[scores["split"].eq(split), "observations"].sum()),
                         "aggregation": "equal measurement-block weights after averaging trace scores within block"}
    return result


def _gem_capacity_check(signals):
    from ystwin.fba.dynamic_rates import build_network_rates
    from ystwin.fba.product_panel import apply_relative_abundances, install_product, prepare_panel_model
    from ystwin.fba.solver import load_model

    model_path = paths.ec_yeast_gem()
    if model_path is None:
        raise FileNotFoundError("the frozen EC GEM is required for --gem-check")
    source_model, solver_settings = load_model(model_path)
    prepared, reference = prepare_panel_model(source_model)
    if reference.loc[reference.orf == "YDL022W", "gene"].tolist() != ["GPD1"]:
        raise ValueError("the structured proteomics mapping does not identify GPD1 as YDL022W")
    peak = signals.loc[signals.nacl_molar == 0.4].sort_values("gpd1_capacity_multiplier").iloc[-1]
    records = []
    for label, ratio in (("reference", 1.0), ("modeled_peak_induction", float(peak.gpd1_capacity_multiplier))):
        controlled, applied = apply_relative_abundances(prepared, reference, {"YDL022W": ratio})
        controlled, task = install_product(controlled, "glycerol", output_compartment="c")
        rates, closure = build_network_rates(
            controlled, "r_1714", "r_1992", "r_4046", product_reaction_id=task.reaction_id,
            uptake_vmax_mmol_per_gdcw_h=20.0 / 6.0, uptake_km_mmol_per_l=0.5,
            oxygen_lower_bound=-1000.0, allocation_fraction=1.0,
            fixed_growth_per_h=0.1, grid_points=3)
        growth, uptake, retained = rates.at_uptake(float(rates.uptake_grid[-1]))
        records.append({
            "condition": label, "modeled_time_s": float(peak.time_s) if ratio != 1 else 0.0,
            "gpd1_relative_abundance": ratio,
            "applied_enzyme_cap_mmol_per_gdcw": float(applied.iloc[0].applied_cap_mmol_per_gdcw),
            "growth_per_h": growth, "glucose_uptake_mmol_per_gdcw_h": uptake,
            "intracellular_glycerol_capacity_mmol_per_gdcw_h": retained,
            "product_reaction_id": task.reaction_id, "output_basis": task.output_basis,
            "ngam_mmol_per_gdcw_h": closure.ngam, "solver": str(solver_settings),
            "interpretation": "glucose-only, aerobic capacity probe at an explicitly chosen 0.1/h growth; not a condition-matched YPD prediction, realized flux or validated titre",
        })
    return pd.DataFrame(records)


def _plot(output, training, fitted_training, predictions, signals):
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    times = training.time_relative_s.to_numpy() / 60
    axes[0].errorbar(times, training.value.to_numpy() * 1e6,
                     yerr=training.sd.to_numpy() * 1e6, fmt="o", color="#203b55", markersize=4,
                     label="Measured, literature-scaled")
    axes[0].plot(times, fitted_training * 1e6, color="#bd6434", label="Fitted source model")
    axes[0].set(title="Fit: WT, 0.4 M NaCl", xlabel="Minutes after shock", ylabel="Processed Hog1PP (micromolar)")
    axes[0].legend(fontsize=8)
    held = predictions.loc[predictions["split"] == "heldout_dose"]
    for index, (_, group) in enumerate(held.groupby("experiment_id", sort=True)):
        color = ("#203b55", "#bd6434", "#62876e")[index % 3]
        x = group.time_relative_s.to_numpy() / 60
        axes[1].plot(x, group.value, "o", color=color, markersize=3, alpha=0.65,
                     label=f"Measured block {index + 1}")
        axes[1].plot(x, group.fitted_prediction, "-", color=color, linewidth=1.5)
    axes[1].set(title="Untuned holdout: WT, 0.8 M", xlabel="Minutes after shock", ylabel="Peak-scaled phospho-Hog1")
    axes[1].legend(fontsize=8)
    for dose, group in signals.groupby("nacl_molar", sort=True):
        axes[2].plot((group.time_s - 3600) / 60, group.gpd1_capacity_multiplier, label=f"{dose:g} M NaCl")
    axes[2].set(title="Source-model metabolic consequence", xlabel="Minutes after shock",
                ylabel="Gpd1 amount / time-matched control")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.axvline(0, color="#899099", linewidth=0.6)
    figure.savefig(output / "hog_learning.png", dpi=160)
    plt.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Fit a source-backed HOG response to measurements and score frozen assay/dose holdouts")
    parser.add_argument("--source-dir", type=Path, default=paths.data_dir() / "hog2013")
    parser.add_argument("--output-dir", type=Path, default=paths.outputs_dir() / "hog_learning")
    parser.add_argument("--max-nfev", type=int, default=30)
    parser.add_argument("--gem-check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args(argv)
    output = args.output_dir
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("use an empty or fresh --output-dir to preserve previous artifacts")
    protocol = json.loads((args.source_dir / "fit_protocol.json").read_text())
    dataset = load_hog_data(args.source_dir, preliminary_time_unit="minutes")
    splits = dataset.split()
    expected = {"fit": protocol["fit"]["expected_observations"],
                **{name: protocol["heldouts"][name]["expected_observations"]
                   for name in ("heldout_assay", "heldout_dose")}}
    if any(len(splits[name]) != count for name, count in expected.items()):
        raise ValueError("observations do not match the frozen protocol")
    model = HogModel.from_source(args.source_dir)
    output.mkdir(parents=True, exist_ok=True)
    _json(output / "protocol.json", {**protocol, "data_metadata": dataset.metadata,
                                     "max_nfev": args.max_nfev,
                                     "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    dataset.observations.to_csv(output / "observations.csv", index=False)
    print("Fitting measured 0.4 M WT Hog1 data; held-out values are not supplied to the optimizer.", flush=True)
    fit = fit_hog_parameters(model, splits["fit"],
                             parameter_names=tuple(protocol["fit"]["estimated_parameters"]),
                             fold_bounds=tuple(protocol["fit"]["search_fold_bounds"]), max_nfev=args.max_nfev)
    _json(output / "fit.json", {**fit.diagnostics, "parameter_values": fit.parameter_values,
                                "learned_quantity": protocol["fit"]["learned_quantity"]})
    heldout = pd.concat([splits["heldout_dose"], splits["heldout_assay"]], ignore_index=True)
    scores, predictions = evaluate_hog_fit(model, fit, heldout)
    scores.to_csv(output / "heldout_scores.csv", index=False)
    predictions.to_csv(output / "heldout_predictions.csv", index=False)
    times = np.unique(np.r_[np.arange(0.0, 14400.1, 60.0), 3605.0])
    control = model.simulate(times, protocol=HogProtocol(0.0), parameters=fit.parameter_values)
    signals, trajectories = [], []
    observable_names = ("Hog1PP_measured", "Gpd1_measured", "gpd1mRNA_measured", "glycerol_measured",
                        "glycerol_i", "glycerol_e", "relVM", "Fps1r", "OD")
    for dose in (0.4, 0.8):
        trajectory = model.simulate(times, protocol=HogProtocol(dose), parameters=fit.parameter_values)
        signals.append(pd.DataFrame({**hog_metabolic_signals(trajectory, control), "nacl_molar": dose}))
        trajectories.append(pd.DataFrame({"time_s": times, "nacl_molar": dose,
                                           **{name: trajectory.variables[name] for name in observable_names}}))
    signals = pd.concat(signals, ignore_index=True)
    signals.to_csv(output / "metabolic_signals.csv", index=False)
    pd.concat(trajectories, ignore_index=True).to_csv(output / "model_native_trajectories.csv", index=False)
    gem_results = None
    if args.gem_check:
        print("Checking the Gpd1-to-GEM capacity connection; intracellular retention is not secretion.", flush=True)
        gem_results = _gem_capacity_check(signals)
        gem_results.to_csv(output / "gem_capacity_check.csv", index=False)
    summary = {
        "protocol_id": protocol["protocol_id"], "fit": fit.diagnostics,
        "heldout_metrics": _summary(scores),
        "validation_scope": protocol["heldouts"]["scope"],
        "heldout_parameters_refitted": False, "mechanistic_structure_learned": False,
        "metabolic_bridge": protocol["metabolic_bridge"],
        "gem_capacity_check": None if gem_results is None else gem_results.to_dict("records"),
        "all_stress_biology_covered": False,
        "gaps": ["Other stress pathways and cross-pathway dynamics remain uncovered by this HOG fit",
                 "Fast phosphorylation/dephosphorylation rates are not separately resolved by minute-sampled data",
                 "The published OD driver, initial conditions and non-fitted reaction parameters remain assumptions/priors",
                 "The source model omits ATP/ADP and NAD(H) dynamics; the GEM check is not a validated cofactor-feedback model",
                 "No condition-matched product-titre or prospective new-laboratory validation has been performed"],
    }
    _json(output / "validation.json", summary)
    if args.plot:
        training = splits["fit"]
        grid = np.unique(np.r_[0.0, training.time_model_s.to_numpy()])
        fitted = model.simulate(grid, protocol=HogProtocol(0.4), parameters=fit.parameter_values)
        predicted = fit.observation_gain * np.interp(training.time_model_s, grid, fitted.variables["Hog1PP_measured"])
        _plot(output, training, predicted, predictions, signals)
    print(json.dumps(_plain(summary["heldout_metrics"]), indent=2))
    print(f"Artifacts: {paths.display_path(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
