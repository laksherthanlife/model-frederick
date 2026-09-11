"""How well the activity inversion recovers an activity it was given, at a stated geometry.

`reporter.py::default_activity_window_h` chooses the smoothing window every dilution-
corrected number in this repository depends on, and it justified that choice with a table
in its own docstring:

    0.5 h -> 14.1%    2 h -> 2.6%    4 h -> 1.1%    8 h -> 0.6%    12 h -> 1.6%

**No script produced that table and no test checked it.** It was a bare assertion in a
comment, which is the failure this repository spends `audit_claims.py` on everywhere else.
This module is the producer. `scripts/estimator_accuracy.py` runs it and writes
`outputs/estimator_accuracy.csv`, and the docstring now quotes cells from that table
through an `audit:value` marker rather than remembering them.

**The table was also measured in a regime the plates are not in, and that is the finding.**
It was taken over a 24-hour run at 10-minute sampling -- 145 points -- where "one sixth of
the run" is 4 hours and sits in the flat part of the curve. The committed exports are
**25 points over 4.00 hours**, so ``duration / 6`` collapses onto the ``4 * dt`` floor and
returns 0.667 h: below the smallest entry in the table it is justified by. The only test of
the rule checked it at 12 h and 48 h, at six and twelve times the plates' density -- the
rule was tested only where it works.

**Everything here is grounded in the committed plates rather than invented**, because a
simulation that answers "how accurate is this estimator" with made-up noise answers a
different question. Read off `data/plates/` and `outputs/sensor_characterisation.csv`:

    geometry          25 points, dt = 1/6 h, duration 4.00 h, 1.54 doublings
    growth            mu falls 0.563 -> 0.217 /h across the run (medians over 251 wells)
    reader noise      1.1% multiplicative on OD and on RFU, as the residual of each
                      channel about a smooth fit in log space
    activity scale    ~1630 RFU/OD/h late, induced ~1.5x over basal -- UPRE2's real fold

Each is named at its constant below with where it came from.

**Maturation is simulated and then not inverted, deliberately.** mCitrine is a YFP and its
chromophore takes time to form; `promoter_activity_from_total` cannot invert maturation at
all -- it raises rather than returning the mature-only answer, because cancelling mu loses
the immature pool. The pipeline therefore runs with `ReporterKinetics(k_deg=0.0)` and no
maturation term, which is an assumption nobody had costed. Over a 24-hour run a 20-minute
maturation is 1.4% of the trace; over a 4-hour run it is 8% of it. :func:`maturation_bias`
measures what that assumption costs at each geometry, by generating with maturation and
inverting without it -- which is exactly what the analysis does to real wells.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from statistics import NormalDist

import numpy as np

from ..reporter import (
    ReporterKinetics,
    default_activity_window_h,
    promoter_activity_from_total,
    simulate_reporter,
)
from ..readings import CorrectedOD, CorrectedRFU

__all__ = [
    "PLATE_DOUBLINGS",
    "PLATE_DURATION_H",
    "PLATE_MU_EARLY",
    "PLATE_MU_LATE",
    "PLATE_POINTS",
    "PLATE_READER_CV",
    "AccuracyPoint",
    "Geometry",
    "NisArm",
    "accuracy_curve",
    "fold_recovery",
    "maturation_bias",
    "measurement_inflation_for_nis",
    "nis_error_decomposition",
    "prediction_metrics",
    "recovery_error",
    "true_activity",
]

# --------------------------------------------------------------------------- measured

PLATE_POINTS = 25
"""Timepoints in every committed export. `data/plates/*__mCitrine*.csv`, all four plates."""

PLATE_DURATION_H = 4.0
"""Hours from first read to last. 0:09:08 to 4:09:08 on 20260722, and the same elsewhere."""

PLATE_DOUBLINGS = 1.54
"""Biomass doublings across the run, well A10 of 20260722: OD 0.195 -> 0.568.

Recorded because it bounds what the correction can be measuring. The whole argument is
that growth dilutes the reporter; over one and a half doublings it has diluted it by a
factor of about three, and no smoothing choice makes a four-hour run a long one."""

PLATE_MU_EARLY = 0.563
"""Median `mu_max` over the 251 wells in `outputs/sensor_characterisation.csv`."""

PLATE_MU_LATE = 0.217
"""Median `mu_late` over the same wells. Growth more than halves across the run, which is
why the simulation ramps mu rather than holding it: a constant-mu trace would make the
inversion's job easier than it is."""

PLATE_READER_CV = 0.011
"""Multiplicative reader noise, per channel, per read.

Measured rather than assumed: the residual standard deviation of well A10's log OD about a
quadratic in time is 0.0107, and of its log RFU about a cubic is 0.0115. Both channels are
given the same 1.1%, which is the honest reading of two numbers that agree to a decimal
place. `reporter_cv` in `sensor_characterisation.csv` has a median of 0.0070, but that is a
within-plate dispersion across wells rather than a read-to-read one, so it is not the same
quantity and is not used here."""

MATURATION_HALF_LIFE_H = 1.0 / 3.0
"""Chromophore maturation half-time for a YFP, hours. 20 minutes.

Venus and mCitrine are the fast-maturing YFPs and are reported in the tens of minutes; the
value is order-of-magnitude and is used only to size a bias, never to correct one. What
matters for :func:`maturation_bias` is the ratio of this to the run length, and at 4 hours
20 minutes is 8% of the trace."""


@dataclass(frozen=True)
class Geometry:
    """One sampling design, as a plate reader would be configured.

    Args:
        label: Name for the report.
        n_points: Reads in the run.
        duration_h: First read to last.
    """

    label: str
    n_points: int
    duration_h: float

    @property
    def dt_h(self) -> float:
        return self.duration_h / (self.n_points - 1)

    def times(self) -> np.ndarray:
        return np.linspace(0.0, self.duration_h, self.n_points)

    @property
    def auto_window_h(self) -> float:
        """What `default_activity_window_h` picks here."""
        return default_activity_window_h(self.times())

    def summary(self) -> str:
        return (f"{self.label}: {self.n_points} points over {self.duration_h:g} h "
                f"(dt {self.dt_h * 60:.0f} min), auto window {self.auto_window_h:.3f} h")


PLATE = Geometry("committed plates", PLATE_POINTS, PLATE_DURATION_H)
"""The geometry every real result in this repository was computed at."""

DOCSTRING_REGIME = Geometry("the regime the table was measured in", 145, 24.0)
"""24 hours at 10-minute sampling, which is what `reporter.py`'s table was taken over."""


# --------------------------------------------------------------------------- the truth


def true_activity(times_h: np.ndarray, basal: float = 1100.0, induced_fold: float = 1.5,
                  lag_h: float = 0.5, rise_h: float = 0.75) -> np.ndarray:
    """Promoter activity under a stressor: basal, a lag, then a rise to a plateau.

    A logistic in time rather than a step, because a transcriptional response has a lag
    while the signal propagates and a finite rise while message and protein accumulate.
    The shape is what makes this a test of the estimator rather than of a constant: a
    smoother recovers a flat line at any window, and the error the table reports is
    entirely about how much of a real transient it flattens.

    Defaults are sized from the plates rather than chosen: `activity_late` has a median of
    1630 RFU/OD/h over 251 wells, and the one induction this repository will defend is
    UPRE2 at 1.0 mM DTT with a fold of 1.50, so a basal of 1100 rising 1.5x lands there.

    Args:
        times_h: Time grid.
        basal: Unstressed activity, RFU/OD/h.
        induced_fold: Plateau over basal.
        lag_h: Delay before the rise begins.
        rise_h: Time constant of the logistic rise.
    """
    t = np.asarray(times_h, dtype=float)
    return basal * (1.0 + (induced_fold - 1.0) / (1.0 + np.exp(-(t - lag_h) / rise_h)))


def _growth(times_h: np.ndarray, early: float = PLATE_MU_EARLY,
            late: float = PLATE_MU_LATE) -> np.ndarray:
    """Specific growth rate falling from ``early`` to ``late`` across the run.

    Linear in time, which is not a claim about the biology -- it is the simplest curve that
    reproduces the two numbers the plates actually pin, `mu_max` at the start and `mu_late`
    at the end. What matters to the estimator is that mu MOVES, because a stressor that
    slows growth is the entire confound being corrected for.
    """
    t = np.asarray(times_h, dtype=float)
    if t[-1] == t[0]:
        return np.full_like(t, early)
    return early + (late - early) * (t - t[0]) / (t[-1] - t[0])


def _observe(times_h: np.ndarray, kinetics: ReporterKinetics, rng: np.random.Generator,
             cv: float, activity: np.ndarray, mu: np.ndarray,
             initial_biomass: float = 0.195) -> tuple[np.ndarray, np.ndarray]:
    """Total fluorescence and biomass a reader would report, with multiplicative noise.

    Biomass starts at 0.195, well A10's first OD reading. The absolute scale cancels out of
    the inversion -- `k_synth = (dF/dt)/X` is homogeneous in it -- so it is here for
    realism rather than for the arithmetic.
    """
    per_cell = simulate_reporter(times_h, activity, mu, kinetics)
    biomass = initial_biomass * np.exp(
        np.concatenate([[0.0], np.cumsum(np.diff(times_h) * (mu[:-1] + mu[1:]) / 2.0)]))
    total = per_cell * biomass
    return (total * rng.lognormal(0.0, cv, total.shape),
            biomass * rng.lognormal(0.0, cv, biomass.shape))


# --------------------------------------------------------------------------- measurement


@dataclass(frozen=True)
class AccuracyPoint:
    """The inversion's error at one geometry and one window."""

    geometry: str
    n_points: int
    duration_h: float
    window_h: float
    window_points: int
    is_auto_window: bool
    median_relative_error: float
    late_window_relative_error: float
    n_replicates: int

    def summary(self) -> str:
        mark = "  <- auto" if self.is_auto_window else ""
        return (f"{self.window_h:6.3f} h ({self.window_points:3d} pts)  "
                f"whole trace {self.median_relative_error:6.1%}   "
                f"late window {self.late_window_relative_error:6.1%}{mark}")


def prediction_metrics(observed, predicted, total_variance, *, nominal_coverage: float = 0.95,
                       step_indices=None) -> dict[str, float | int]:
    """Coverage, NIS and signed lag-1 autocorrelation of *prior* predictions.

    Intervals use a Gaussian moment approximation, including measurement variance.
    They are prediction intervals, not posterior-state credible intervals. No normality
    or whiteness is established just by computing these statistics. Do not concatenate
    different wells or horizons into one series: summarize each separately first.

    Nonfinite values and nonpositive variances are counted as excluded, never as
    covered. Autocorrelation uses the full finite-series sum of squares, but its
    numerator only uses adjacent original steps. A gap must not become a lag-1 pair.
    Supply original integer ``step_indices`` when rows have already been omitted.
    A constant or fewer-than-three-valid-point series has undefined autocorrelation.
    """
    actual, mean, variance = (np.asarray(x, dtype=float)
                              for x in (observed, predicted, total_variance))
    if actual.ndim != 1 or actual.shape != mean.shape or actual.shape != variance.shape:
        raise ValueError("observed, predicted and total_variance must share a 1D shape")
    if not np.isfinite(nominal_coverage) or not 0 < nominal_coverage < 1:
        raise ValueError("nominal_coverage must be finite and in (0, 1)")
    indices = np.arange(actual.size) if step_indices is None else np.asarray(step_indices)
    if (indices.shape != actual.shape or not np.all(np.isfinite(indices))
            or np.any(indices != np.floor(indices)) or np.any(np.diff(indices) <= 0)):
        raise ValueError("step_indices must match the shape and be strictly increasing integers")
    usable = np.isfinite(actual) & np.isfinite(mean) & np.isfinite(variance) & (variance > 0)
    z = np.full(actual.shape, np.nan)
    z[usable] = (actual[usable] - mean[usable]) / np.sqrt(variance[usable])
    count = int(usable.sum())
    adjacent = usable[:-1] & usable[1:] & (np.diff(indices) == 1)
    autocorrelation = float("nan")
    if count >= 3:
        centred = z - np.mean(z[usable])
        denominator = float(centred[usable] @ centred[usable])
        if denominator > 0 and adjacent.any():
            autocorrelation = float(np.sum((centred[:-1] * centred[1:])[adjacent]) / denominator)
    critical = NormalDist().inv_cdf((1.0 + nominal_coverage) / 2.0)
    return {
        "n_predictions": count,
        "n_excluded": int(actual.size - count),
        "nominal_coverage": float(nominal_coverage),
        "coverage": float(np.mean(np.abs(z[usable]) <= critical)) if count else float("nan"),
        "mean_nis": float(np.mean(z[usable] ** 2)) if count else float("nan"),
        "mean_standardised": float(np.mean(z[usable])) if count else float("nan"),
        "lag1_autocorr": autocorrelation,
        "lag1_pairs": int(adjacent.sum()),
    }


# ------------------------------------------------------------- NIS error decomposition


@dataclass(frozen=True)
class NisArm:
    """One NIS run, differing from the baseline in exactly one declared quantity.

    Args:
        label: Name for the report.
        family: ``baseline``, ``particles`` or ``measurement``. Exactly one baseline.
        n_particles: Ensemble size the arm ran at.
        sigma_multiplier: Factor applied to both channels' declared relative sigma.
        mean_nis: The arm's aggregate mean NIS in this channel.
        measurement_variance_share: Median measurement share of predictive variance.
        min_ess: Median across wells of each well's minimum posterior ESS.
        min_unique_ancestors: Median across wells of surviving initial ancestors.
        coverage: Median empirical coverage of the nominal interval.
        well_channels: Well-channels the arm aggregated.
    """

    label: str
    family: str
    n_particles: int
    sigma_multiplier: float
    mean_nis: float
    measurement_variance_share: float
    min_ess: float
    min_unique_ancestors: float
    coverage: float
    well_channels: int


def measurement_inflation_for_nis(mean_nis: float, measurement_variance_share: float,
                                  *, target: float = 1.0) -> float:
    """The relative-sigma multiplier that would bring NIS to *target*, if measurement
    error were the whole story.

    Exact under one stated assumption: that inflating the declared measurement sigma
    leaves the prediction and the state variance where they are. It does not -- a looser
    likelihood reweights the ensemble -- which is why this is the prediction that the
    empirical measurement arm then falsifies or confirms, not the answer itself.

    With ``S`` the squared innovation and ``V = V_state + V_meas`` the predictive
    variance, ``NIS = S / V`` and ``f = V_meas / V``. Scaling sigma by ``c`` sends
    ``V -> (1 - f) V + c**2 f V``, so ``c**2 = (NIS / target - (1 - f)) / f``. At
    ``f = 1`` this is ``sqrt(NIS / target)``, and any smaller share demands more: an
    under-declared measurement scale is the *cheapest* explanation of an inflated NIS,
    so the multiplier this returns is a lower bound on what any measurement-only story
    would have to claim.

    Returns:
        The multiplier, which is below one for an over-dispersed arm and ``nan`` if the
        inputs are unusable. It saturates at ``0.0`` where no positive multiplier reaches
        *target* -- the filter is over-dispersed by more than its whole declared
        measurement variance, so removing all of it still lands below *target* -- and at
        ``inf`` where an above-target arm declares no measurement variance to scale.
    """
    if not np.isfinite(mean_nis) or not np.isfinite(measurement_variance_share):
        return float("nan")
    if not np.isfinite(target) or target <= 0:
        raise ValueError("target must be a positive finite NIS")
    if not 0.0 <= measurement_variance_share <= 1.0:
        raise ValueError("measurement_variance_share must be a fraction of predictive variance")
    if measurement_variance_share == 0.0:
        # No declared measurement variance to inflate; scaling it changes nothing.
        return float("inf") if mean_nis > target else 0.0
    squared = (mean_nis / target - (1.0 - measurement_variance_share)) / measurement_variance_share
    return float(np.sqrt(squared)) if squared > 0.0 else 0.0


def _ess_ratio(arm_min_ess: float, baseline_min_ess: float) -> float:
    """How far an arm moved the ensemble it was not supposed to move.

    ``nan`` rather than an infinity where the baseline reports no ESS at all: there is
    nothing to have moved away from, and a ratio would invent one.
    """
    if not np.isfinite(arm_min_ess) or not np.isfinite(baseline_min_ess):
        return float("nan")
    return arm_min_ess / baseline_min_ess if baseline_min_ess > 0.0 else float("nan")


def nis_error_decomposition(arms, *, target: float = 1.0) -> list[dict]:
    """What each arm removes of the baseline's NIS excess, arm by arm.

    NIS above *target* is unexplained squared innovation per unit declared predictive
    variance. The question this answers is which declared quantity it is unexplained
    *by*: raise the ensemble size and the particle approximation is what moved; raise
    the declared measurement sigma with the dynamics held and the measurement model is
    what moved; whatever neither buys is process-model error, by elimination among the
    three and only among the three.

    ``excess_removed`` is a share of ``baseline - target``, so it is comparable across
    channels whose baselines differ. It is signed and unclipped: an arm that makes NIS
    worse reports a negative share, which is a result and not an error. Where the
    baseline is already at or below *target* there is no excess to apportion and every
    share is ``nan`` rather than a fraction of a negative quantity.

    That elimination is only as clean as the arms are separate, and ``min_ess_ratio`` is
    the column that says whether they were. Both knobs act on the same weight degeneracy:
    a looser likelihood keeps particles alive exactly as a larger ensemble does, so a
    measurement arm whose ratio has left one has moved the particle approximation too and
    its ``excess_removed`` is not attributable to the measurement model alone. The ratio
    is reported rather than corrected for, because nothing here identifies the split.

    Args:
        arms: :class:`NisArm` values for one channel, exactly one of family ``baseline``.
        target: The NIS a correctly declared filter would average. One, by construction.

    Returns:
        One row per arm, baseline first, then in the order supplied.
    """
    arms = list(arms)
    baselines = [arm for arm in arms if arm.family == "baseline"]
    if len(baselines) != 1:
        raise ValueError("a decomposition needs exactly one baseline arm")
    if len({arm.label for arm in arms}) != len(arms):
        raise ValueError("arm labels must be unique")
    unknown = {arm.family for arm in arms} - {"baseline", "particles", "measurement"}
    if unknown:
        raise ValueError(f"unknown decomposition families: {sorted(unknown)}")
    baseline = baselines[0]
    excess = baseline.mean_nis - target
    rows = []
    for arm in [baseline] + [arm for arm in arms if arm is not baseline]:
        removed = float("nan")
        if np.isfinite(excess) and excess > 0.0 and np.isfinite(arm.mean_nis):
            removed = (baseline.mean_nis - arm.mean_nis) / excess
        rows.append({
            "label": arm.label, "family": arm.family, "n_particles": arm.n_particles,
            "sigma_multiplier": arm.sigma_multiplier, "well_channels": arm.well_channels,
            "mean_nis": arm.mean_nis, "nis_excess": arm.mean_nis - target,
            "excess_removed": removed,
            "measurement_variance_share": arm.measurement_variance_share,
            "implied_sigma_multiplier": measurement_inflation_for_nis(
                arm.mean_nis, arm.measurement_variance_share, target=target),
            "min_ess": arm.min_ess, "min_unique_ancestors": arm.min_unique_ancestors,
            "min_ess_ratio": _ess_ratio(arm.min_ess, baseline.min_ess),
            "coverage": arm.coverage, "baseline_has_excess": bool(excess > 0.0),
            "baseline_mean_nis": baseline.mean_nis, "target_nis": float(target),
        })
    return rows


def recovery_error(geometry: Geometry, window_h: float, *, cv: float = PLATE_READER_CV,
                   n_replicates: int = 200, seed: int = 0,
                   kinetics: ReporterKinetics | None = None,
                   generate_with: ReporterKinetics | None = None,
                   late_fraction: float = 0.75) -> tuple[float, float]:
    """Median relative error of the recovered activity, whole trace and late window.

    The late-window figure is the one the results rest on: `run_sensor_characterisation.py`
    reports `activity_late`, the mean of the recovered activity over the last quarter of
    the trace, and every fold change in `outputs/` is a ratio of two of those.

    Args:
        geometry: Sampling design.
        window_h: Smoothing window to invert at.
        cv: Multiplicative reader noise per channel per read.
        n_replicates: Independent noise draws.
        seed: Base seed.
        kinetics: What the INVERSION assumes. Defaults to the pipeline's own
            ``ReporterKinetics(k_deg=0.0)``.
        generate_with: What the SIMULATION uses. Defaults to ``kinetics``; pass a maturing
            reporter here and a non-maturing one above to measure what ignoring maturation
            costs, which is what :func:`maturation_bias` does.
        late_fraction: Where the late window starts, as a fraction of the points.
    """
    kinetics = ReporterKinetics(k_deg=0.0) if kinetics is None else kinetics
    generate_with = kinetics if generate_with is None else generate_with
    t = geometry.times()
    truth = true_activity(t)
    mu = _growth(t)
    late = slice(int(len(t) * late_fraction), None)

    whole, tail = [], []
    for replicate in range(n_replicates):
        rng = np.random.default_rng(seed + replicate)
        total, biomass = _observe(t, generate_with, rng, cv, truth, mu)
        recovered = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass), kinetics,
                                                 window_h=window_h)
        whole.append(np.median(np.abs(recovered - truth) / truth))
        # The late window enters as a MEAN and only then as a ratio, so the quantity that
        # matters is the error of that mean rather than the mean of the pointwise errors.
        # They are different numbers and the second is the larger one; reporting it would
        # overstate the damage to a fold change.
        tail.append(abs(np.mean(recovered[late]) - np.mean(truth[late])) / np.mean(truth[late]))
    return float(np.median(whole)), float(np.median(tail))


def accuracy_curve(geometry: Geometry, windows_h: tuple[float, ...] | None = None,
                   **kwargs) -> list[AccuracyPoint]:
    """The error against window at one geometry, including whatever the rule picks.

    The auto window is always included and flagged, because the point of the table is to
    say where on the curve the code actually sits rather than where it could sit.
    """
    auto = geometry.auto_window_h
    if windows_h is None:
        windows_h = (0.5, 1.0, 2.0, 4.0, 8.0, 12.0)
    # A window wider than the run is not a window; savgol clips it to the trace length and
    # the report would show several identical rows under different labels.
    candidates = sorted({w for w in windows_h if w <= geometry.duration_h} | {auto})
    points = []
    for window in candidates:
        whole, tail = recovery_error(geometry, window, **kwargs)
        n = int(round(window / geometry.dt_h))
        points.append(AccuracyPoint(
            geometry=geometry.label, n_points=geometry.n_points,
            duration_h=geometry.duration_h, window_h=window,
            window_points=n + (n % 2 == 0), is_auto_window=bool(np.isclose(window, auto)),
            median_relative_error=whole, late_window_relative_error=tail,
            n_replicates=kwargs.get("n_replicates", 200)))
    return points


def maturation_bias(geometry: Geometry, *, half_life_h: float = MATURATION_HALF_LIFE_H,
                    window_h: float | None = None, **kwargs) -> dict[str, float]:
    """What the pipeline's no-maturation assumption costs, at one geometry.

    Generated with a maturing reporter and inverted without one, which is what the analysis
    does to every real well: `run_sensor_characterisation.py` sets
    ``ReporterKinetics(k_deg=0.0)`` and `promoter_activity_from_total` refuses maturation
    outright, so there is no path by which a maturation term could reach a committed number.

    Returned as a bias rather than folded into :func:`accuracy_curve`, because it is not
    noise: it is the same sign on every well and every replicate, so it does not shrink
    with more plates and does not appear in a cluster bootstrap. A fold change is a ratio
    of two such biases and cancels much of it -- which is the reason to measure it rather
    than to assume either that it matters or that it does not.
    """
    if not np.isfinite(half_life_h) or half_life_h <= 0:
        raise ValueError("half_life_h must be finite and positive")
    if "generate_with" in kwargs:
        raise ValueError("maturation_bias sets generate_with; supply kinetics for the baseline loss")
    kinetics = kwargs.get("kinetics") or ReporterKinetics(k_deg=0.0)
    if kinetics.has_maturation:
        raise ValueError("maturation_bias requires a non-maturing inversion baseline")
    window_h = geometry.auto_window_h if window_h is None else window_h
    maturing = replace(kinetics, k_mat=float(np.log(2.0) / half_life_h))
    clean_whole, clean_late = recovery_error(geometry, window_h, **kwargs)
    biased_whole, biased_late = recovery_error(
        geometry, window_h, generate_with=maturing, **kwargs)
    return {
        "geometry": geometry.label,
        "window_h": window_h,
        "maturation_half_life_h": half_life_h,
        "half_life_as_fraction_of_run": half_life_h / geometry.duration_h,
        "late_error_no_maturation": clean_late,
        "late_error_with_maturation": biased_late,
        "added_by_maturation": biased_late - clean_late,
        "whole_trace_no_maturation": clean_whole,
        "whole_trace_with_maturation": biased_whole,
    }


def fold_recovery(geometry: Geometry, *, window_h: float | None = None,
                  true_fold: float = 1.5, dosed_mu_retained: float = 0.75,
                  half_life_h: float | None = None, cv: float = PLATE_READER_CV,
                  n_replicates: int = 200, seed: int = 0,
                  late_fraction: float = 0.75) -> dict[str, float]:
    """What the estimator returns for a fold change whose true value is known.

    **This is the quantity every committed result is, and the one worth measuring.**
    Nothing in `outputs/` is a single well's activity: `late_window_sensitivity.csv` holds
    ratios of a dosed well's late-window activity to its own plate's zero-dose control, and
    the calls this repository will defend are such ratios clearing 1.0 -- seven of them at
    four biological replicates, five at three, six at two. The count is deliberately not
    repeated here; `outputs/fold_multiplicity.csv` carries it and moves when the panel does.

    A ratio cancels anything common to both wells, which is why the ~10% maturation bias on
    a single well in :func:`maturation_bias` is not the number to quote. What it cannot
    cancel is the part that DIFFERS between the two, and two things differ by construction:
    the dosed well is induced, so its activity is changing while the control's is flat and a
    lag bites harder on a moving signal; and the dose slows growth, which is the entire
    confound the inversion exists to remove.

    Args:
        geometry: Sampling design.
        window_h: Smoothing window. Defaults to what the rule picks here.
        true_fold: The induction to recover. 1.5 is UPRE2 at 1.0 mM DTT, the one fold in
            this repository whose interval clears 1.0 at every window tried.
        dosed_mu_retained: Fraction of the control's growth rate the dosed well keeps. 0.75
            is the middle of what the plates show at a resolvable dose, and it is the whole
            reason a naive RFU/OD ratio overstates the fold.
        half_life_h: Maturation half-time to GENERATE with, or ``None`` to generate with
            the same instantaneous-maturation reporter the inversion assumes.
        cv: Reader noise per channel per read.
        n_replicates: Independent noise draws.
        seed: Base seed.
        late_fraction: Where the late window starts.
    """
    window_h = geometry.auto_window_h if window_h is None else window_h
    inverts_with = ReporterKinetics(k_deg=0.0)
    generates_with = (inverts_with if half_life_h is None else
                      ReporterKinetics(k_deg=0.0, k_mat=float(np.log(2.0) / half_life_h)))

    t = geometry.times()
    late = slice(int(len(t) * late_fraction), None)
    control_activity = true_activity(t, induced_fold=1.0)
    dosed_activity = true_activity(t, induced_fold=true_fold)
    control_mu = _growth(t)
    dosed_mu = control_mu * dosed_mu_retained

    # The fold the estimator is being asked for is the ratio of the two LATE-WINDOW MEAN
    # activities, not `true_fold` -- the logistic has not fully plateaued by the end of a
    # four-hour run, so the achievable answer is slightly below the asymptote. Scoring
    # against the asymptote would charge the estimator for the shape of the truth.
    target = float(np.mean(dosed_activity[late]) / np.mean(control_activity[late]))

    corrected, naive = [], []
    for replicate in range(n_replicates):
        rng = np.random.default_rng(seed + replicate)
        folds = []
        for activity, mu in ((dosed_activity, dosed_mu), (control_activity, control_mu)):
            total, biomass = _observe(t, generates_with, rng, cv, activity, mu)
            recovered = promoter_activity_from_total(t, CorrectedRFU(total), CorrectedOD(biomass), inverts_with,
                                                     window_h=window_h)
            folds.append((float(np.mean(recovered[late])), float(np.mean((total / biomass)[late]))))
        corrected.append(folds[0][0] / folds[1][0])
        naive.append(folds[0][1] / folds[1][1])

    corrected = np.asarray(corrected)
    return {
        "geometry": geometry.label,
        "window_h": float(window_h),
        "maturation_half_life_h": float("nan") if half_life_h is None else half_life_h,
        "target_fold": target,
        "recovered_fold": float(np.median(corrected)),
        "relative_error": float(np.median(corrected) / target - 1.0),
        "spread_of_recovered_fold": float(np.std(corrected)),
        "naive_fold": float(np.median(naive)),
        "n_replicates": n_replicates,
    }
