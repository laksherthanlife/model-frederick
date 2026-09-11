"""Score a prediction on a held-out split of the real plates.

The first scored forward prediction in this project on data the model did not see. Every
other real-data result here is diagnostic -- G1 pass rates, variance shares, gate
verdicts, corrected dose responses -- which describe the plates rather than predict
anything about them. `docs/CLAIM_BOUNDARY.md` records the consequence: nothing is Tier 3,
predicted-then-measured.

This closes that on the smallest honest target. `outputs/split_manifest.csv` already
defines the partitions and refuses the leaky ones; until now nothing read it. Here the
manifest picks the held-out groups, a dose-response is fitted on the rest, and the
held-out activity is scored against baselines.

The baselines are the point. An RMSE alone says nothing -- the question is whether the
model beats predicting the training mean, or carrying the nearest measured dose forward.
Skill is reported against the better of the two, so a negative number means the model is
worse than the trivial answer and says so.

Usage: python scripts/run_heldout_score.py [--kind extrapolation]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths
from ystwin.analysis.nulls import skill_score_if_defined
from ystwin.generator.panel_calibration import _model as _biphasic_model
from ystwin.generator.panel_calibration import fit_dose_response

TARGET = "activity_late"


def _plate_key(value) -> str:
    """Normalise both sources onto the eight-digit date.

    The manifest carries the plate as a number and the characterisation table as a
    filename stem, so a naive join silently produces no rows -- which reads as "no
    overlap" rather than "wrong key".
    """
    text = str(value)
    digits = "".join(c for c in text if c.isdigit())
    return digits[:8]


def _fit_linear(train: pd.DataFrame) -> dict:
    """Log-linear activity against dose, per construct. Two parameters."""
    fits = {}
    for construct, group in train.groupby("construct"):
        x = np.log1p(group.dose_mM.to_numpy(dtype=float))
        y = group[TARGET].to_numpy(dtype=float)
        usable = np.isfinite(x) & np.isfinite(y)
        if usable.sum() < 3 or np.std(x[usable]) < 1e-12:
            continue
        fits[construct] = np.polyfit(x[usable], y[usable], 1)
    return fits


def _fit_biphasic(train: pd.DataFrame) -> dict:
    """The repo's own biphasic form: basal plus saturating induction, scaled by viability.

    `(basal + amplitude d/(ec50+d)) / (1 + (d/lethal)^hill)`, four parameters, fitted by
    generator/panel_calibration.fit_dose_response. This is the shape the measured ladders
    actually have -- the calibration module exists because a saturating Hill could not
    describe the fall at high dose -- so it is the fair model to score, and the linear fit
    is kept only as a second baseline.
    """
    fits = {}
    for construct, group in train.groupby("construct"):
        doses = group.dose_mM.to_numpy(dtype=float)
        activity = group[TARGET].to_numpy(dtype=float)
        usable = np.isfinite(doses) & np.isfinite(activity)
        if usable.sum() < 5:
            continue
        try:
            fits[construct] = fit_dose_response(doses[usable], activity[usable])
        except (ValueError, RuntimeError):
            continue
    return fits


def _why_unscoreable(kind, rows, readings, train, test):
    """Name the cause, not the row count.

    "0 test rows" reads as a data gap when both real causes are structural: a held-out
    plate that the characterisation pipeline never produced rows for, or a training side
    too small to fit. Both are visible here, so report them.
    """
    assigned = rows[rows.assignment.notna()]
    held = sorted({str(r.plate_key) for r in assigned.itertuples()
                   if r.assignment == "test"})
    present = {str(k) for k in readings.plate_key.unique()}
    missing = [p for p in held if p not in present]
    if missing:
        return (f"unscoreable: the manifest holds out plate(s) {missing}, which are absent "
                f"from the readings table (it carries {sorted(present)}). Those plates have "
                "no reporter-channel blank, so characterisation never produced rows for "
                "them -- zero test rows permanently, not a join fault")
    if test.empty:
        return (f"unscoreable: {len(train)} train rows and no test row survived the join on "
                "(plate, construct, dose); the manifest and the readings table disagree on "
                "at least one of those keys")
    return (f"unscoreable: {len(train)} training rows, below the 6 needed to fit both the "
            f"biphasic and log-linear forms ({len(test)} test rows were available)")


def _split_by_identifiability(fits: dict) -> tuple[set, set]:
    """Separate fits the calibration module trusts from those it does not.

    ``DoseFit.identifiable`` is the module's own verdict, with the reason attached, and
    the scorer was ignoring it. It should not: on these splits the flag predicts the
    outcome exactly. Where any fit is unidentifiable the reason is that the estimated
    lethal dose sits below the estimated EC50 -- the fit is saying the culture dies
    before it induces, which is not a dose-response and cannot be extrapolated from.
    """
    good = {c for c, f in fits.items() if getattr(f, "identifiable", False)}
    return good, set(fits) - good


def _baseline_prediction(train, doses, name):
    doses = np.asarray(doses, dtype=float)
    if name == "train_mean":
        return np.full(doses.shape, float(train[TARGET].mean()))
    if name != "nearest_dose":
        raise ValueError(f"unknown baseline {name!r}")
    means = train.groupby("dose_mM", sort=True)[TARGET].mean()
    distance = np.abs(doses[:, None] - means.index.to_numpy(dtype=float)[None, :])
    closest = np.isclose(distance, distance.min(axis=1, keepdims=True), rtol=0, atol=1e-12)
    return (closest @ means.to_numpy(dtype=float)) / closest.sum(axis=1)


def _model_prediction(fit, doses, name):
    doses = np.asarray(doses, dtype=float)
    if name == "loglinear":
        return fit[0] * np.log1p(doses) + fit[1]
    return np.asarray(_biphasic_model(doses, fit.ec50, fit.lethal_dose,
                                      fit.amplitude, fit.basal), dtype=float)


def _select_training_choices(train, kind="interpolation"):
    models, baselines, diagnostics = {}, {}, {}
    plate_column = "plate_key" if "plate_key" in train else "plate"
    group_column = plate_column if kind == "heldout_replicate" else "dose_mM"
    for construct, group in train.groupby("construct", sort=True):
        splits = []
        if group_column in group:
            for held in sorted(group[group_column].unique()):
                test_mask = group[group_column] == held
                train_mask = ~test_mask
                if kind == "extrapolation":
                    train_mask = group.dose_mM < held
                if kind == "interpolation" and held in (group.dose_mM.min(), group.dose_mM.max()):
                    continue
                if train_mask.sum() >= 5 and group.loc[train_mask, "dose_mM"].nunique() >= 4:
                    splits.append((group[train_mask], group[test_mask]))
        errors = {name: [] for name in ("train_mean", "nearest_dose", "loglinear", "biphasic")}
        for inner, validation in splits:
            doses = validation.dose_mM.to_numpy(dtype=float)
            truth = validation[TARGET].to_numpy(dtype=float)
            for name in ("train_mean", "nearest_dose"):
                predicted = _baseline_prediction(inner, doses, name)
                errors[name].append(float(np.mean((predicted - truth) ** 2)))
            for name, fitter in (("loglinear", _fit_linear), ("biphasic", _fit_biphasic)):
                fitted = fitter(inner).get(construct)
                if fitted is None or (name == "biphasic" and not fitted.identifiable):
                    errors[name].append(float("inf"))
                else:
                    predicted = _model_prediction(fitted, doses, name)
                    errors[name].append(float(np.mean((predicted - truth) ** 2)))
        scores = {name: float(np.mean(values)) if values else float("inf")
                  for name, values in errors.items()}
        models[construct] = min(("loglinear", "biphasic"), key=scores.get)
        baselines[construct] = min(("train_mean", "nearest_dose"), key=scores.get)
        diagnostics[construct] = {"inner_folds": len(splits), "losses": scores,
                                  "model_selection_estimable": np.isfinite(scores[models[construct]])}
    method = ("inner leave-plate-out on outer training only" if kind == "heldout_replicate"
              else "inner forward-dose splits on outer training only" if kind == "extrapolation"
              else "inner leave-dose-out on outer training only")
    return models, baselines, diagnostics, method


def _score(kind: str, manifest: pd.DataFrame, readings: pd.DataFrame,
           seed: int | None = None) -> dict | None:
    rows = manifest[manifest.split_kind == kind]
    if seed is not None:
        rows = rows[rows.seed == seed]
    if rows.empty:
        return None
    if rows.seed.nunique() != 1 or rows.partition_hash.nunique() != 1:
        raise ValueError("score one seed and partition at a time; merged partitions can leak")
    keyed = {}
    for row in rows.itertuples():
        if pd.isna(row.dose_mM):
            continue
        key = (row.plate_key, row.construct, round(float(row.dose_mM), 6))
        if key in keyed and keyed[key] != row.assignment:
            raise ValueError(f"conflicting train/test assignments for {key}")
        keyed[key] = row.assignment

    readings = readings.copy()
    readings["assignment"] = [
        keyed.get((k, c, round(float(d), 6)))
        for k, c, d in zip(readings.plate_key, readings.construct, readings.dose_mM)
    ]
    train = readings[readings.assignment == "train"].dropna(subset=[TARGET])
    test = readings[readings.assignment == "test"].dropna(subset=[TARGET])
    if len(train) < 6 or test.empty:
        return {"kind": kind, "note": _why_unscoreable(kind, rows, readings, train, test)}

    linear, biphasic = _fit_linear(train), _fit_biphasic(train)
    selected_models, selected_baselines, inner_diagnostics, selection_method = _select_training_choices(train, kind)
    identifiable, refused = _split_by_identifiability(biphasic)
    bi, lin, actual, mean_base, near_base, trusted = [], [], [], [], [], []
    selected_prediction, selected_baseline = [], []
    for row in test.itertuples():
        if row.construct not in biphasic or row.construct not in linear:
            continue
        fit = biphasic[row.construct]
        trusted.append(row.construct in identifiable)
        bi.append(_biphasic_model(row.dose_mM, fit.ec50, fit.lethal_dose,
                                  fit.amplitude, fit.basal))
        slope, intercept = linear[row.construct]
        lin.append(slope * np.log1p(row.dose_mM) + intercept)
        actual.append(getattr(row, TARGET))

        same = train[train.construct == row.construct]
        mean_base.append(float(_baseline_prediction(same, [row.dose_mM], "train_mean")[0]))
        near_base.append(float(_baseline_prediction(same, [row.dose_mM], "nearest_dose")[0]))
        selected_baseline.append(mean_base[-1] if selected_baselines[row.construct] == "train_mean"
                                 else near_base[-1])
        selected_prediction.append(bi[-1] if selected_models[row.construct] == "biphasic" else lin[-1])

    if not actual:
        # The model fits per construct, so a split withholding whole constructs leaves
        # every test row without a fit -- structural, not a data gap. Say which.
        return {"kind": kind, "note": (
            f"unscoreable by this model class: the fits are per construct, and this split "
            f"withholds whole constructs. Test carries {sorted(set(test.construct))}, "
            f"fitted {sorted(biphasic)} -- no overlap, so no test row can be predicted")}

    truth = np.asarray(actual, dtype=float)
    keep = np.asarray(trusted, dtype=bool)
    rmse = lambda p, m=None: float(np.sqrt(np.nanmean(  # noqa: E731
        ((np.asarray(p, dtype=float) - truth)[m if m is not None else slice(None)]) ** 2)))
    model, linear_rmse = rmse(bi), rmse(lin)
    mean_only, nearest_only = rmse(mean_base), rmse(near_base)
    best_baseline = rmse(selected_baseline)
    result = {
        "kind": kind, "n_train": len(train), "n_test": len(actual),
        "n_test_available": len(test), "n_unscoreable_test": len(test) - len(actual),
        "n_identifiable_fits": len(identifiable), "n_refused_fits": len(refused),
        "rmse_biphasic": model, "rmse_loglinear": linear_rmse,
        "rmse_train_mean": mean_only, "rmse_nearest_dose": nearest_only,
        "rmse_selected_baseline": best_baseline,
        "rmse_selected_model": rmse(selected_prediction),
        "skill_selected_model": skill_score_if_defined(rmse(selected_prediction), best_baseline, perfect=0.0),
        "skill_selected_vs_train_mean": skill_score_if_defined(rmse(selected_prediction), mean_only, perfect=0.0),
        "skill_selected_vs_nearest": skill_score_if_defined(rmse(selected_prediction), nearest_only, perfect=0.0),
        "selected_baselines": json.dumps(selected_baselines, sort_keys=True),
        "selected_models": json.dumps(selected_models, sort_keys=True),
        "baseline_selection": selection_method,
        "model_selection_estimable": all(bool(v["model_selection_estimable"]) for v in inner_diagnostics.values()),
        "skill_method": "1 - RMSE / RMSE of training-selected baseline",
        "skill_biphasic": skill_score_if_defined(model, best_baseline, perfect=0.0),
        "skill_loglinear": skill_score_if_defined(linear_rmse, best_baseline, perfect=0.0),
        "seed": int(rows.seed.iloc[0]),
        "partition_hash": rows.partition_hash.iloc[0][:12],
    }
    # The headline number should come only from fits the module vouches for. Scoring the
    # rest alongside is what turned a refusal into an apparent failure of the model.
    if keep.any():
        result["rmse_biphasic_identifiable"] = rmse(bi, keep)
        result["skill_biphasic_identifiable"] = skill_score_if_defined(
            rmse(bi, keep), rmse(selected_baseline, keep), perfect=0.0)
        result["n_test_identifiable"] = int(keep.sum())
    return result


def _exhaustive_interpolation(readings: pd.DataFrame) -> pd.DataFrame:
    """Score every interior rung, one at a time, rather than sampling seeds.

    Sampling was confounded: the seeded splits withheld two, three and four doses at
    different seeds, so partitions differed in training size as well as in which dose was
    held out, and a spread over them mixed the two. There are only five interior rungs, so
    the space is small enough to exhaust -- which makes the partitions comparable and the
    result complete rather than sampled.
    """
    from ystwin.analysis.splits import SplitNotPossible, make_split

    inventory = readings[["plate", "construct", "stressor", "dose_mM",
                          TARGET, "well"]].dropna()
    interior = sorted(inventory.dose_mM.unique())[1:-1]
    rows = []
    for dose in interior:
        try:
            split = make_split(inventory, "interpolation",
                               group_key=["plate", "construct", "dose_mM"],
                               within=["construct"], hold_out=dose, dataset="real")
        except SplitNotPossible as exc:
            rows.append({"held_out_dose": dose, "note": str(exc)[:70]})
            continue
        tag = split.assign_rows(inventory).to_numpy()
        train, test = inventory[tag == "train"], inventory[tag == "test"]
        if test.empty or len(train) < 20:
            rows.append({"held_out_dose": dose,
                         "note": f"{len(train)} train, {len(test)} test"})
            continue

        fits = _fit_biphasic(train)
        linear = _fit_linear(train)
        models, baselines, _, method = _select_training_choices(train)
        vouched, _ = _split_by_identifiability(fits)
        predicted, actual, nearest, baseline, selected = [], [], [], [], []
        for construct, group in test.groupby("construct"):
            if construct not in vouched:
                continue
            fit, same = fits[construct], train[train.construct == construct]
            for row in group.itertuples():
                predicted.append(_biphasic_model(row.dose_mM, fit.ec50, fit.lethal_dose,
                                                 fit.amplitude, fit.basal))
                actual.append(getattr(row, TARGET))
                nearest.append(float(_baseline_prediction(same, [row.dose_mM], "nearest_dose")[0]))
                baseline.append(float(_baseline_prediction(same, [row.dose_mM], baselines[construct])[0]))
                chosen_fit = fit if models[construct] == "biphasic" else linear[construct]
                selected.append(float(_model_prediction(chosen_fit, [row.dose_mM], models[construct])[0]))
        if not actual:
            rows.append({"held_out_dose": dose, "note": "no identifiable fit"})
            continue
        truth = np.asarray(actual, dtype=float)
        rmse = lambda p: float(np.sqrt(np.nanmean(  # noqa: E731
            (np.asarray(p, dtype=float) - truth) ** 2)))
        rows.append({
            "held_out_dose": dose, "n_test": len(actual),
            "n_fits_vouched": len(vouched),
            "rmse_biphasic": rmse(predicted), "rmse_nearest_dose": rmse(nearest),
            "rmse_selected_baseline": rmse(baseline), "rmse_selected_model": rmse(selected),
            "skill": skill_score_if_defined(rmse(predicted), rmse(baseline), perfect=0.0),
            "skill_selected_model": skill_score_if_defined(rmse(selected), rmse(baseline), perfect=0.0),
            "selected_baselines": json.dumps(baselines, sort_keys=True),
            "selected_models": json.dumps(models, sort_keys=True),
            "baseline_selection": method,
            "skill_method": "fixed biphasic vs training-selected baseline; RMSE skill",
            "comparison_cohort": "outer-training identifiable biphasic fits; excluded test rows counted separately",
            "n_test_available": len(test), "n_unscoreable_test": len(test) - len(actual),
            # Split.hash is a property, so the old `callable(split.hash)` guard was
            # always False and this column was always blank -- the one table in the
            # pipeline that cites the published claim recorded no provenance for it.
            "partition_hash": split.hash[:12],
        })
    return pd.DataFrame(rows)


def _score_products(out, n_draws=1000, confidence=0.95):
    from ystwin.pathway.flux import (
        carotenoid_measurements, score_product_validation, summarize_product_validation,
    )
    from ystwin.pathway.spec import load_pathway
    from ystwin.analysis.validation import matched_product_branches

    states = carotenoid_measurements()
    spec = load_pathway("beta_carotene")
    result = score_product_validation(spec, states, draws=n_draws, nominal_coverage=confidence, seed=0)
    summary = summarize_product_validation(result)
    matched = matched_product_branches(spec, states, result)
    matched_scores = pd.DataFrame([
        {"branch": branch, "n_independent_strains": group.strain.nunique(),
         "n_condition_predictions": len(group), "n_failed_predictions": int(group.status.ne("ok").sum()),
         "joint_rmse_log": float(np.sqrt(np.mean([
             np.mean(values.joint_log_mse.to_numpy()) for _, values in group.groupby("strain")])))}
        for branch, group in matched.groupby("branch", sort=True)
    ])
    matched.to_csv(out / "product_matched_entry_comparisons.csv", index=False)
    matched_scores.to_csv(out / "product_matched_entry_scores.csv", index=False)
    result.to_csv(out / "product_prediction.csv", index=False)
    summary.to_csv(out / "product_prediction_scores.csv", index=False)
    pd.DataFrame(result.attrs["selection_scores"]).to_csv(
        out / "product_prediction_selection.csv", index=False)
    pd.DataFrame(result.attrs["comparison_predictions"]).to_csv(
        out / "product_prediction_comparisons.csv", index=False)
    print(f"Nested product validation: {states.strain.nunique()} held-out strains, {len(states)} condition predictions, two paired products.")
    print("Entry-gene and model selection, entry and branch refits, and baseline fitting use only training-strain folds.")
    print(summary.to_string(index=False))
    print("Matched-entry ablation: same selected gene and entry law, only the branch changes.")
    print(matched_scores.to_string(index=False))
    print(f"Coverage is descriptive, not a validated {confidence:.0%} guarantee: {states.strain.nunique()} independent strains.")
    print("The shared validator refits entry and kinetics together while resampling training strains and observation uncertainty.")
    print(f"-> {out / 'product_prediction_scores.csv'}")
    return result


def _score_uncertainty(out, n_trials, n_resamples, confidence):
    from ystwin.analysis.uncertainty import fold_change_coverage

    table = pd.DataFrame([
        fold_change_coverage(n_plates, fold_cv=fold_cv, n_trials=n_trials,
                             n_resamples=n_resamples, confidence=confidence)
        for n_plates in (3, 6, 12) for fold_cv in (0.0, 0.15)
    ])
    table.to_csv(out / "fold_bootstrap_coverage.csv", index=False)
    print(table.to_string(index=False))
    print(f"-> {out / 'fold_bootstrap_coverage.csv'}")
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", default=None, help="score one split kind only")
    parser.add_argument("--target", choices=("reporter", "product", "uncertainty", "all"), default="reporter")
    parser.add_argument("--output-dir", type=pathlib.Path, default=None)
    parser.add_argument("--input-dir", type=pathlib.Path, default=None)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--coverage-trials", type=int, default=100)
    args = parser.parse_args()

    out = args.output_dir if args.output_dir is not None else paths.outputs_dir()
    out.mkdir(parents=True, exist_ok=True)
    if args.target in ("uncertainty", "all"):
        _score_uncertainty(out, args.coverage_trials, args.draws, args.confidence)
        if args.target == "uncertainty":
            return 0
    if args.target in ("product", "all"):
        _score_products(out, args.draws, args.confidence)
        if args.target == "product":
            return 0
    source = args.input_dir if args.input_dir is not None else paths.outputs_dir()
    manifest_path, readings_path = source / "split_manifest.csv", source / "sensor_characterisation.csv"
    for path in (manifest_path, readings_path):
        if not path.exists():
            print(f"missing {path.name}; run make_splits.py and "
                  "run_sensor_characterisation.py first", file=sys.stderr)
            return 2

    manifest = pd.read_csv(manifest_path)
    manifest = manifest[manifest.dataset == "real_biosensor"].copy()
    manifest["plate_key"] = manifest.plate.map(_plate_key)
    readings = pd.read_csv(readings_path)
    readings["plate_key"] = readings.plate.map(_plate_key)

    available = sorted(manifest.split_kind.unique())
    if args.kind and args.kind not in available:
        # Without this the run scores nothing, writes an empty heldout_scores.csv over
        # the tracked one, and exits 0 -- a typo in the flag reads as a result.
        print(f"no split kind {args.kind!r} in the manifest; it has {available}",
              file=sys.stderr)
        return 2
    kinds = [args.kind] if args.kind else available
    seeds = sorted(manifest.seed.unique())
    results = [r for r in (_score(k, manifest, readings, seed=sd)
                           for k in kinds for sd in seeds) if r]

    width = 84
    print("=" * width)
    print("Held-out prediction on the real plates, scored against trivial baselines")
    print("=" * width)
    # One line per (kind, seed), then the spread -- a single partition cannot tell a
    # real margin from a lucky held-out dose, and interpolation withholds a different
    # rung at every seed.
    for r in results:
        print()
        if "note" in r:
            print(f"  {r['kind']:18s} seed {r.get('seed', '?')}  {r['note']}")
            continue
        print(f"  {r['kind']:18s} seed {r['seed']}  {r['n_train']} train / "
              f"{r['n_test']} scored of {r['n_test_available']} test wells   partition {r['partition_hash']}")
        print(f"    RMSE  biphasic {r['rmse_biphasic']:9.1f}   log-linear {r['rmse_loglinear']:9.1f}"
              f"   train-mean {r['rmse_train_mean']:9.1f}   nearest-dose {r['rmse_nearest_dose']:9.1f}")
        print(f"    fits the calibration module vouches for: "
              f"{r['n_identifiable_fits']} of {r['n_identifiable_fits'] + r['n_refused_fits']}")
        print(f"    training-selected baseline: {r['selected_baselines']}")
        print(f"    training-selected model: {r['selected_models']}")
        for baseline_name, baseline_key in (("train mean", "skill_selected_vs_train_mean"),
                                            ("nearest dose", "skill_selected_vs_nearest")):
            value = r[baseline_key]
            print(f"    nested model vs {baseline_name}: "
                  + ("undefined" if value is None else f"{value:+.3f}"))
        for name, key in (("nested selected", "skill_selected_model"),
                          ("all fits       ", "skill_biphasic"),
                          ("log-linear     ", "skill_loglinear"),
                          ("identifiable   ", "skill_biphasic_identifiable")):
            if key not in r:
                print(f"    skill {name} --       -> no identifiable fit; nothing to score")
                continue
            if r[key] is None:
                print(f"    skill {name} --       -> selected baseline is perfect; skill undefined")
                continue
            verdict = "BEATS the training-selected baseline" if r[key] > 0 else "worse than trivial"
            print(f"    skill {name} {r[key]:+.3f}   -> {verdict}")

    frame = pd.DataFrame(results)
    frame.to_csv(out / "heldout_scores.csv", index=False)

    exhaustive = _exhaustive_interpolation(readings)
    exhaustive.to_csv(out / "heldout_interpolation_by_dose.csv", index=False)
    scored_rows = exhaustive[exhaustive.get("skill_selected_model").notna()] \
        if "skill_selected_model" in exhaustive else exhaustive.iloc[:0]
    if not scored_rows.empty:
        print()
        print("=" * width)
        print("Interpolation, every interior rung held out in turn (comparable partitions)")
        print("=" * width)
        print("  Common cohort: outer-training identifiable biphasic fits; exclusions are counted in the table.")
        print(f"  {'held out':>9} {'wells':>6} {'biphasic':>10} {'nested':>10} {'chosen base':>11} {'skill':>8}")
        for r in exhaustive.itertuples():
            if pd.isna(getattr(r, "skill_selected_model", np.nan)):
                print(f"  {r.held_out_dose:>9} {getattr(r, 'note', '')}")
                continue
            print(f"  {r.held_out_dose:>9.1f} {r.n_test:>6} {r.rmse_biphasic:>10.1f}"
                  f" {r.rmse_selected_model:>10.1f} {r.rmse_selected_baseline:>11.1f} {r.skill_selected_model:>+8.3f}")
        values = scored_rows.skill_selected_model
        wins = int((values > 0).sum())
        print(f"    nested-model median {values.median():+.3f}   "
              f"range [{values.min():+.3f}, {values.max():+.3f}]   "
              f"wins {wins}/{len(scored_rows)}")
        print(f"  -> {out / 'heldout_interpolation_by_dose.csv'}")

    scored = frame[frame.get("skill_selected_model").notna()] \
        if "skill_selected_model" in frame else frame.iloc[:0]
    if not scored.empty:
        print()
        print("  Across seeds, on models selected inside training folds, on the common fittable cohort:")
        for kind, group in scored.groupby("kind"):
            vals = group.skill_selected_model
            distinct = group.partition_hash.nunique()
            note = "" if distinct > 1 else "  (deterministic split: one partition)"
            print(f"    {kind:18s} skill {vals.median():+.3f} "
                  f"[{vals.min():+.3f}, {vals.max():+.3f}] over {len(vals)} seed(s)"
                  f"{note}")
    print(f"\n  -> {out / 'heldout_scores.csv'}")
    print("  Skill is 1.0 for a perfect prediction, 0.0 for matching the baseline, and")
    print("  negative for doing worse than it. A negative number is the honest outcome")
    print("  when a dose-response does not extrapolate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
