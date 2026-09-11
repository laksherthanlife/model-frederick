from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy.optimize import isotonic_regression, least_squares

from .. import reporter
from ..qpcr import STRESSOR_FOR_CONSTRUCT
from ..readings import SpecificFluorescence
from .latent import fit_latent
from .stress_model import StressModel

__all__ = ["StressEstimate", "CalibratedSensorModel", "fit_sensor_model"]

_AXES = {"UPR": ("UPRE1", "UPRE2"), "oxidative": ("NativeYap1", "AlteredYap1")}
_CONSTRUCTS = tuple(name for pair in _AXES.values() for name in pair)
_FEATURES = tuple(f"{name}_activity_fold" for name in _CONSTRUCTS)
_KEYS = ["plate", "construct", "stressor", "dose_mM"]
_COLUMNS = (*_KEYS, "well", "mu_late", "activity_late", "naive_late", "n_points")
_HEALTHY_GROWTH_RETENTION = 0.5
_SYNTHETIC_SIZE = 96
_SYNTHETIC_VALIDATION_SIZE = 24
_TIMES_H = np.linspace(0.0, 4.0, 25)
_KINETICS = reporter.ReporterKinetics(k_deg=0.0, k_mat=None)


@dataclass(frozen=True)
class StressEstimate:
    upr: float
    oxidative: float
    growth_retention: float
    source: str
    extrapolated: bool

    @property
    def modules(self) -> dict[str, float]:
        return {"UPR": self.upr, "oxidative": self.oxidative}


@dataclass(frozen=True)
class _AxisCalibration:
    stressor: str
    constructs: tuple[str, str]
    half_dose_mM: float
    max_dose_mM: float

    def induction(self, dose: float) -> float:
        return dose / (self.half_dose_mM + dose) * (
            self.half_dose_mM + self.max_dose_mM) / self.max_dose_mM

    def dose_at(self, induction: float) -> float:
        fraction = float(np.clip(induction, 0.0, 1.0)) * self.max_dose_mM / (
            self.half_dose_mM + self.max_dose_mM)
        return self.half_dose_mM * fraction / (1.0 - fraction)


@dataclass(frozen=True)
class _ConstructCalibration:
    basal_activity: float
    basal_growth: float
    induction_gain: float
    growth_doses: tuple[float, ...]
    growth_folds: tuple[float, ...]
    feature_bounds: tuple[float, float]

    def growth_at(self, dose: float) -> float:
        return float(np.interp(dose, self.growth_doses, self.growth_folds))


def _checked_seed(seed: int) -> int:
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    return int(seed)


def _sensor_frame(readings: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(readings, pd.DataFrame):
        raise TypeError("readings must be a sensor-characterisation DataFrame")
    if not readings.columns.is_unique:
        raise ValueError("sensor column names must be unique")
    missing = sorted(set(_COLUMNS) - set(readings.columns))
    if missing:
        raise ValueError(f"missing sensor columns: {missing}")
    return readings.loc[:, list(_COLUMNS)].copy()


def _conditions(readings: pd.DataFrame) -> pd.DataFrame:
    frame = _sensor_frame(readings)
    if frame.empty:
        raise ValueError("no sensor conditions supplied")
    for column in ("plate", "construct", "stressor", "well"):
        if not frame[column].map(lambda value: isinstance(value, str) and bool(value.strip())).all():
            raise ValueError(f"{column} identifiers must be nonempty strings")
    unknown = set(frame.construct) - set(_CONSTRUCTS)
    if unknown:
        raise ValueError(f"unsupported sensor constructs: {sorted(unknown)}")
    if not frame.stressor.eq(frame.construct.map(STRESSOR_FOR_CONSTRUCT)).all():
        raise ValueError("stressor/construct pairing must be UPRE1/UPRE2:DTT or Yap1:H2O2")
    if frame.duplicated(["plate", "well"]).any():
        raise ValueError("duplicate plate/well measurements are not independent replicates")
    for column in ("dose_mM", "mu_late", "activity_late", "naive_late", "n_points"):
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(float)
    if not (np.isfinite(frame.dose_mM) & (frame.dose_mM >= 0)).all():
        raise ValueError("dose_mM must be finite and nonnegative")
    rows = []
    for key, group in frame.groupby(_KEYS, sort=True):
        row = dict(zip(_KEYS, key))
        row["n_wells"] = len(group)
        for column in ("mu_late", "activity_late", "naive_late", "n_points"):
            values = group[column].to_numpy()
            row[column] = float(values.mean()) if np.isfinite(values).all() else np.nan
        points = group.n_points.to_numpy()
        valid_points = bool(np.isfinite(points).all() and np.all((points >= 5) & (points % 1 == 0)))
        row["reasons"] = [] if valid_points else ["insufficient_or_missing_timepoints"]
        for column, label in (("mu_late", "growth"), ("activity_late", "activity"),
                              ("naive_late", "specific_signal")):
            values = group[column].to_numpy()
            valid = bool(np.all(np.isfinite(values) & (values > 0)))
            row[f"{label}_valid"] = valid and valid_points
            if not valid:
                row["reasons"].append(f"{label}_nonpositive_or_missing")
        row["activity_valid"] &= row["specific_signal_valid"]
        rows.append(row)
    controls = {tuple(row[key] for key in _KEYS[:3]): row for row in rows if row["dose_mM"] == 0}
    for row in rows:
        control = controls.get(tuple(row[key] for key in _KEYS[:3]))
        for label, column, fold in (("activity", "activity_late", "activity_fold"),
                                     ("growth", "mu_late", "growth_retention")):
            paired = control is not None and control[f"{label}_valid"]
            row[f"{label}_paired"] = bool(paired and row[f"{label}_valid"])
            row[fold] = row[column] / control[column] if row[f"{label}_paired"] else np.nan
            if not paired:
                row["reasons"].append(f"missing_positive_within_plate_{label}_control")
    return pd.DataFrame(rows)


def _support(conditions: pd.DataFrame, pair: tuple[str, str]):
    windows, healthy_indices = {}, set()
    for plate, group in conditions[conditions.construct.isin(pair)].groupby("plate", sort=True):
        at_zero = group[group.dose_mM == 0]
        if len(at_zero) != 2 or not (at_zero.activity_paired & at_zero.growth_paired).all():
            windows[plate] = None
            continue
        ceiling = 0.0
        for dose, at_dose in group.groupby("dose_mM", sort=True):
            valid = at_dose.activity_paired & at_dose.growth_paired & (
                at_dose.growth_retention >= _HEALTHY_GROWTH_RETENTION)
            if len(at_dose) != 2 or not valid.all():
                break
            ceiling = float(dose)
            healthy_indices.update(at_dose.index)
        windows[plate] = [0.0, ceiling]
    ceilings = [window[1] for window in windows.values() if window is not None]
    if not ceilings or min(ceilings) <= 0:
        raise ValueError(f"no positive common healthy dose window with paired controls for {pair}")
    ceiling = min(ceilings)
    used = conditions.loc[sorted(healthy_indices)]
    used = used[used.dose_mM <= ceiling]
    if used.dose_mM.nunique() < 3:
        raise ValueError(f"induction for {pair} needs a control and at least two healthy positive doses")
    return ceiling, windows, used


def _simulate_features(panel, axes, doses) -> np.ndarray:
    features = {}
    for axis in axes.values():
        dose = float(doses.get(axis.stressor, 0.0))
        induction = axis.induction(dose)
        for construct in axis.constructs:
            calibration = panel[construct]
            mu = np.full_like(_TIMES_H, calibration.basal_growth * calibration.growth_at(dose))
            synthesis = np.full_like(
                _TIMES_H, calibration.basal_activity * (1 + calibration.induction_gain * induction))
            trace = reporter.simulate_reporter(
                _TIMES_H, synthesis, mu, _KINETICS,
                r0=calibration.basal_activity / calibration.basal_growth)
            activity = reporter.promoter_activity(
                _TIMES_H, SpecificFluorescence(trace), mu, _KINETICS)
            fold = float(np.mean(activity[int(0.75 * len(activity)):]) / calibration.basal_activity)
            if not np.isfinite(fold) or fold <= 0:
                raise ValueError(f"ODE-derived activity must be finite and positive for {construct}")
            features[construct] = fold
    return np.array([features[name] for name in _CONSTRUCTS])


def _synthetic_features(panel, axes, truth):
    return np.array([
        _simulate_features(panel, axes, {axis.stressor: axis.dose_at(state[index])
                                        for index, axis in enumerate(axes.values())})
        for state in truth
    ])


@dataclass(frozen=True)
class CalibratedSensorModel:
    feature_names: tuple[str, ...]
    metadata: dict = field(repr=False)
    _panel: dict[str, _ConstructCalibration] = field(repr=False)
    _axes: dict[str, _AxisCalibration] = field(repr=False)
    _readout: StressModel = field(repr=False)

    def infer(self, features: Mapping[str, float], *, source: str = "observed") -> StressEstimate:
        if not isinstance(features, Mapping):
            raise TypeError("features must be a mapping of named activity-fold channels")
        if set(features) != set(self.feature_names):
            raise ValueError(f"features require exactly these named channels: {self.feature_names}")
        values = []
        for name in self.feature_names:
            value = np.asarray(features[name])
            if value.ndim != 0 or value.dtype.kind not in "iuf":
                raise ValueError(f"feature {name} must be a numeric scalar")
            values.append(float(value))
        values = np.array(values)
        if not np.all(np.isfinite(values) & (values > 0)):
            raise ValueError("activity-fold features must be finite and strictly positive")
        if not isinstance(source, str) or not source.strip():
            raise ValueError("source must be a nonempty string")
        modules = np.clip(self._readout.predict_modules(values).to_numpy()[0], 0.0, 1.0)
        retention = 1.0
        for state, axis in zip(modules, self._axes.values()):
            dose = axis.dose_at(float(state))
            retention *= float(np.mean([self._panel[name].growth_at(dose) for name in axis.constructs]))
        bounds = np.array([self._panel[name].feature_bounds for name in _CONSTRUCTS])
        extrapolated = bool(np.any(values < bounds[:, 0] - 1e-8)
                            or np.any(values > bounds[:, 1] + 1e-8))
        return StressEstimate(float(modules[0]), float(modules[1]), float(retention), source, extrapolated)

    def predict_condition(
        self, stressor: str | None, dose_mM: float = 0.0, *, seed: int = 0,
    ) -> StressEstimate:
        _checked_seed(seed)
        supported = {axis.stressor: axis for axis in self._axes.values()}
        if stressor is not None and stressor not in supported:
            raise ValueError(f"unsupported stressor {stressor!r}; use DTT, H2O2, or None")
        dose = np.asarray(dose_mM)
        if dose.ndim != 0 or dose.dtype.kind not in "iuf":
            raise ValueError("dose_mM must be a finite nonnegative scalar")
        dose = float(dose)
        if not np.isfinite(dose) or dose < 0 or (stressor is None and dose != 0):
            raise ValueError("dose_mM must be finite, nonnegative, and zero when stressor is None")
        readings = _simulate_features(self._panel, self._axes, {stressor: dose})
        result = self.infer(
            dict(zip(self.feature_names, readings)),
            source=f"SIMULATED sensors: calibrated reporter ODE; {stressor or 'unstressed'} {dose:g} mM")
        return replace(result, extrapolated=result.extrapolated or bool(
            stressor is not None and dose > supported[stressor].max_dose_mM))

    def validation_table(self, readings: pd.DataFrame) -> pd.DataFrame:
        conditions = _conditions(readings)
        by_stressor = {axis.stressor: axis for axis in self._axes.values()}
        simulated, rows = {}, []
        for row in conditions.to_dict("records"):
            construct, stressor, dose = row["construct"], row["stressor"], row["dose_mM"]
            calibration = self._panel[construct]
            if (stressor, dose) not in simulated:
                simulated[stressor, dose] = _simulate_features(self._panel, self._axes, {stressor: dose})
            activity_fold = float(simulated[stressor, dose][_CONSTRUCTS.index(construct)])
            growth_fold = calibration.growth_at(dose)
            result = {key: row[key] for key in _KEYS}
            result.update({
                "n_wells": row["n_wells"],
                "held_out": row["plate"] not in self.metadata["training_plates"],
                "extrapolated": dose > by_stressor[stressor].max_dose_mM,
                "activity_valid": row["activity_paired"],
                "growth_valid": row["growth_paired"],
                "exclusion_reason": "; ".join(row["reasons"]),
                "activity_late_measured": row["activity_late"],
                "activity_late_predicted": activity_fold * calibration.basal_activity,
                "activity_fold_measured": row["activity_fold"],
                "activity_fold_predicted": activity_fold,
                "mu_late_measured": row["mu_late"],
                "mu_late_predicted": growth_fold * calibration.basal_growth,
                "growth_retention_measured": row["growth_retention"],
                "growth_retention_predicted": growth_fold,
            })
            for quantity, valid in (("activity_late", row["activity_valid"]),
                                     ("activity_fold", row["activity_paired"]),
                                     ("mu_late", row["growth_valid"]),
                                     ("growth_retention", row["growth_paired"])):
                error = result[f"{quantity}_predicted"] - result[f"{quantity}_measured"] if valid else np.nan
                result[f"{quantity}_error"] = error
                result[f"{quantity}_absolute_error"] = abs(error)
            rows.append(result)
        return pd.DataFrame(rows)


def fit_sensor_model(
    readings: pd.DataFrame, *, training_plates: tuple[str, ...] | None = None, seed: int = 0,
) -> CalibratedSensorModel:
    seed = _checked_seed(seed)
    frame = _sensor_frame(readings)
    if training_plates is not None:
        if (not isinstance(training_plates, tuple) or not training_plates
                or any(not isinstance(plate, str) or not plate.strip() for plate in training_plates)
                or len(set(training_plates)) != len(training_plates)):
            raise ValueError("training_plates must be a nonempty tuple of unique plate IDs")
        absent = set(training_plates) - set(frame.plate)
        if absent:
            raise ValueError(f"training plates absent from readings: {sorted(absent)}")
        frame = frame[frame.plate.isin(training_plates)].copy()
    conditions = _conditions(frame)
    plates = sorted(conditions.plate.unique())
    axes, panel, ranges, quality, activity_indices = {}, {}, {}, {}, set()
    for name, pair in _AXES.items():
        ceiling, windows, activity = _support(conditions, pair)
        stressor = STRESSOR_FOR_CONSTRUCT[pair[0]]
        dose = activity.dose_mM.to_numpy()
        channel = np.array([pair.index(construct) for construct in activity.construct])
        observed = activity.activity_fold.to_numpy()
        minimum = float(np.min(dose[dose > 0]))
        upper_gains = np.array([max(1e-10, float(observed[channel == c].max() - 1)) for c in (0, 1)])

        def residual(parameters):
            half, gains = parameters[0], parameters[1:]
            fraction = dose / (half + dose) * (half + ceiling) / ceiling
            return 1 + gains[channel] * fraction - observed

        fitted = least_squares(
            residual, np.r_[np.median(dose[dose > 0]), upper_gains / 2],
            bounds=(np.r_[minimum, 0.0, 0.0], np.r_[ceiling, upper_gains]),
            max_nfev=1000, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        if not fitted.success:
            raise ValueError(f"bounded induction calibration failed for {name}: {fitted.message}")
        if np.linalg.norm(fitted.x[1:]) < 1e-8:
            raise ValueError(f"{name} induction is not identifiable from healthy training conditions")
        axes[name] = _AxisCalibration(stressor, pair, float(fitted.x[0]), ceiling)
        activity_indices.update(activity.index)
        all_doses = conditions.loc[conditions.construct.isin(pair), "dose_mM"]
        ranges[stressor] = {
            "dose_mM": [0.0, ceiling],
            "measured_dose_mM": [float(all_doses.min()), float(all_doses.max())],
            "per_plate_supported_dose_mM": windows,
        }
        quality[name] = {
            "activity_fold_rmse": float(np.sqrt(np.mean(fitted.fun**2))),
            "half_dose_mM": float(fitted.x[0]),
            "half_dose_bounds_mM": [minimum, ceiling],
            "bounded_parameter_indices": [int(index) for index in np.flatnonzero(fitted.active_mask)],
            "n_paired_conditions": len(activity),
        }
        for index, construct in enumerate(pair):
            selected = activity[activity.construct == construct]
            basal = float(selected.loc[selected.dose_mM == 0, "activity_late"].mean())
            growth = conditions[(conditions.construct == construct) & conditions.growth_paired]
            controls = growth[growth.dose_mM == 0]
            knots = growth.groupby("dose_mM").growth_retention.agg(["mean", "count"])
            folds = np.clip(isotonic_regression(
                knots["mean"].to_numpy(), weights=knots["count"].to_numpy(), increasing=False).x, 0.0, 1.0)
            folds[0] = 1.0
            panel[construct] = _ConstructCalibration(
                basal, float(controls.mu_late.mean()), float(fitted.x[index + 1]),
                tuple(float(value) for value in knots.index), tuple(float(value) for value in folds),
                (min(1.0, float(selected.activity_fold.min())),
                 max(1.0, float(selected.activity_fold.max()))))
    rng = np.random.default_rng(seed)
    truth = rng.uniform(0.0, 1.0, (_SYNTHETIC_SIZE, 2))
    truth[:4] = [[0, 0], [1, 0], [0, 1], [1, 1]]
    synthetic = _synthetic_features(panel, axes, truth)
    latent = fit_latent(synthetic, n_states=2, seed=seed)
    basis = StressModel(
        list(_FEATURES), list(_AXES), synthetic.mean(axis=0), latent.loadings,
        np.var(latent.states, axis=0), max(latent.heldout_error, 1e-12),
        np.zeros((3, 2)), {})
    zero = basis.infer(np.ones(4))[0]
    weights, *_ = np.linalg.lstsq(basis.infer(synthetic) - zero, truth, rcond=None)
    basis = replace(basis, readout=np.vstack([weights, -zero @ weights]))
    probe_truth = rng.uniform(0.0, 1.0, (_SYNTHETIC_VALIDATION_SIZE, 2))
    probe = _synthetic_features(panel, axes, probe_truth)
    errors = basis.predict_modules(probe).to_numpy() - probe_truth
    rmse = {name: float(np.sqrt(np.mean(errors[:, index]**2))) for index, name in enumerate(_AXES)}
    for index, construct in enumerate(_CONSTRUCTS):
        low, high = panel[construct].feature_bounds
        panel[construct] = replace(panel[construct], feature_bounds=(
            min(low, float(synthetic[:, index].min())), max(high, float(synthetic[:, index].max()))))
    exclusions = []
    for index, row in conditions.iterrows():
        if index in activity_indices and row.growth_paired:
            continue
        reasons = list(row.reasons)
        stressor = row.stressor
        window = ranges[stressor]["per_plate_supported_dose_mM"].get(row.plate)
        if window is None:
            reasons.append("unpaired_construct_control")
        elif row.dose_mM > window[1]:
            reasons.append("outside_within_plate_healthy_window")
        if row.growth_paired and row.growth_retention < _HEALTHY_GROWTH_RETENTION:
            reasons.append("growth_below_half_within_plate_control")
        if row.dose_mM > ranges[stressor]["dose_mM"][1]:
            reasons.append("outside_common_training_dose_range")
        exclusions.append({
            **{key: row[key] for key in _KEYS}, "n_wells": int(row.n_wells),
            "activity_used": index in activity_indices, "growth_used": bool(row.growth_paired),
            "reasons": reasons or ["unpaired_construct_condition"],
        })
    metadata = {
        "training_plates": plates,
        "measurement_columns_used": list(_COLUMNS),
        "synthetic_training_size": _SYNTHETIC_SIZE,
        "synthetic_training_seed": seed,
        "synthetic_validation_size": _SYNTHETIC_VALIDATION_SIZE,
        "seed": seed,
        "synthetic_module_rmse": rmse,
        "observation_schema": {
            feature: {"construct": construct, "units": "dimensionless", "baseline": 1.0,
                      "definition": "activity_late (RFU/OD/h) divided by mean positive 0 mM control activity from the same plate and construct"}
            for feature, construct in zip(_FEATURES, _CONSTRUCTS)
        },
        "generator": "ystwin.reporter.simulate_reporter LSODA ODE, followed by promoter_activity inversion",
        "latent_model": "two-factor fit_latent basis and baseline-anchored supervised StressModel readout; labels are simulated induction fractions",
        "healthy_window": {
            "minimum_growth_retention": _HEALTHY_GROWTH_RETENTION,
            "definition": "contiguous doses with positive activity/signal and mu_late >= 0.5 times the measured same-plate 0 mM control in both constructs; intersect supported plate windows",
            "status": "declared operating criterion, not a measured viability threshold or score-tuned exclusion",
        },
        "calibrated_ranges": ranges,
        "calibration_quality": quality,
        "construct_calibration": {
            construct: {
                "basal_activity_RFU_per_OD_h": calibration.basal_activity,
                "basal_growth_per_h": calibration.basal_growth,
                "induction_increment_at_supported_limit": calibration.induction_gain,
                "growth_doses_mM": list(calibration.growth_doses),
                "growth_retention_knots": list(calibration.growth_folds),
                "activity_fold_bounds": list(calibration.feature_bounds),
            }
            for construct, calibration in panel.items()
        },
        "normalization_controls": [
            {**{key: row[key] for key in _KEYS[:3]},
             "activity_late": float(row.activity_late) if row.activity_valid else None,
             "mu_late": float(row.mu_late) if row.growth_valid else None}
            for _, row in conditions[conditions.dose_mM == 0].iterrows()
        ],
        "exclusions": exclusions,
        "assumptions": {
            "reporter_kinetics": {"k_deg_per_h": 0.0, "k_mat_per_h": None,
                                  "status": "fixed priors, not measured by the supplied summary table; instantaneous maturation and no active degradation"},
            "time_grid": {"duration_h": 4.0, "n_points": 25, "late_fraction": 0.75,
                          "status": "fixed simulation grid; duration and initial conditions are not identified by summary measurements"},
            "dynamics": "constant effective synthesis and mu_late over the simulated read, starting at basal activity / basal growth; transient reporter dilution is integrated, not an inferred culture time course",
            "induction": "shared Hill exponent 1 and shared bounded dose-shape per named construct pair; nonnegative per-construct induction gains; unity is the supported-dose endpoint, not full biological activation",
            "observation": "effective reader-unit activity scale; gain, dry weight per OD and autofluorescence are not separately identified or subtracted beyond the upstream table correction",
            "activity_estimator": "the table's total-signal activity and the simulated per-OD inversion express the same reporter balance; different numerical/noise biases are not identified from summary data",
            "growth": "monotone bounded isotonic dose-growth fits to within-plate mu_late ratios; measured stimulation above 1 remains in residuals but a growth cap cannot reward it",
            "growth_retention": "mean fitted construct-pair growth retention at the dose implied by the inferred axis; paired-axis retentions multiply for mixed states, an unmeasured independence assumption; empirical growth-cap effect, not a mechanistic ATP cost",
            "synthetic_design": "independent uniform induction states; no cross-axis activation; mixtures and module truth are simulation assumptions, not measured pathway labels",
            "prediction_seed": "condition predictions are deterministic noise-free ODE means; seed is validated but no observation-noise sample is drawn",
            "extrapolation": "dose outside the common healthy range or activity outside fitted/training channel bounds is flagged; module fractions are bounded to [0,1], growth interpolation holds endpoint values outside measured knots",
            "validation": "no refitting: raw predictions use training control scales; held-out controls define measured target folds only, never prediction normalization",
        },
        "evidence_status": {
            "sensor_calibration": "wet-lab sensor summary calibration; biological transfer must be assessed with validation_table on held-out plates",
            "latent_truth": "synthetic only; neither two-axis pathway ground truth nor mechanistic rate identifiability is established",
            "supported_modules": list(_AXES),
            "product_measurements_used": False,
            "product_predictions": "prospective SIMULATED sensors only; no measured product-associated sensor states",
            "growth_effect": "sensor-side phenotype association, not causal attribution or ATP-per-RFU conversion",
            "fitted_growth_retention_at_supported_endpoint": {
                name: float(np.mean([panel[construct].growth_at(axis.max_dose_mM)
                                     for construct in axis.constructs]))
                for name, axis in axes.items()
            },
        },
    }
    return CalibratedSensorModel(_FEATURES, metadata, panel, axes, basis)
