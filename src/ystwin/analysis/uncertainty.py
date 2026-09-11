"""Intervals for the dilution-corrected fold changes, and refusal when there are none.

The reported sensor-characterisation result is a ratio: recovered promoter activity at a
dose over the same quantity at zero dose. It is currently reported as a bare number --
``0.96``, ``1.01``, ``0.97`` -- and read as though those are distinguishable from 1.0.
Whether they are is a question about spread, and spread was never computed.

Two things decide the answer and neither is optional.

**The unit of replication is the plate, not the well.** Wells on one plate share an
inoculum, a medium batch, a reader, and a position in an incubator. Resampling them as if
independent inflates the effective sample size by the number of wells per plate and
produces an interval that is too narrow by roughly the square root of that factor.
``analysis/power.py`` already refuses that mistake for the design side; this module is the
matching refusal for the analysis side.

**With two plates there is no interval to compute.** A cluster bootstrap resamples
clusters, and two clusters admit three distinct resamples, of which two are degenerate.
Whatever comes out is a statement about which of two plates was drawn, not about the
population of plates. The honest output is a refusal, so
:func:`fold_change` returns an unestimable result below ``MIN_PLATES_FOR_INTERVAL`` rather
than a tight-looking number. That is the whole point of the module: a missing interval has
to be visible, because a fold change without one silently reads as certain.

The corresponding positive statement -- "the induction is *absent*, not merely
undetected" -- is an equivalence claim, not a failure to reject, and :func:`equivalence`
makes it with two one-sided tests against a margin the caller has to name in advance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from ..growth import growth_rate_uncertainty
from ..reporter import ReporterKinetics, default_activity_window_h, promoter_activity
from ..readings import CorrectedOD, SpecificFluorescence, require

__all__ = [
    "ActivityUncertainty",
    "EquivalenceVerdict",
    "FoldChange",
    "MIN_PLATES_FOR_INTERVAL",
    "activity_uncertainty",
    "equivalence",
    "fold_change",
    "fold_change_coverage",
    "window_sensitivity",
]

MIN_PLATES_FOR_INTERVAL = 3
"""Fewest biological replicates for which a cluster bootstrap is reported at all.

Three is already poor -- the resample distribution has 10 distinct outcomes and the
coverage of a percentile interval at that cluster count is well below nominal. It is set
here as the floor at which a *widened* interval is more informative than silence, not as a
count anyone should design to. Two is refused outright.
"""

_LOG_EPS = 1e-12


# ---------------------------------------------------------------------------
# fold change with a cluster bootstrap
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FoldChange:
    """A dose's activity relative to its zero-dose control, with what is known about it.

    Args:
        construct: Reporter the fold belongs to.
        dose: Dose in the stressor's own units.
        point: Ratio of the aggregated dosed activity to the aggregated control.
        low: Lower interval bound, or NaN when not estimable.
        high: Upper interval bound, or NaN when not estimable.
        method: How the interval was produced, or why it was not.
        n_plates: Biological replicates carrying both the dose and its control.
        n_plates_contributing: How many of those actually entered the estimate -- the
            per-plate folds that are finite and positive, which is exactly the set
            :func:`_combine` averages. **This is the number that governs the interval**,
            and it is not always ``n_plates``: a plate whose recovered activity is negative
            drops out of the geometric mean silently, and reporting the interval as though
            it were there inflates the evidence. See :func:`fold_change`.
        n_wells: Total wells aggregated, for the record.
        n_resamples: Bootstrap resamples drawn.
        estimable: False when the interval was refused. Callers that report a fold
            change must check this rather than formatting NaN bounds.
    """

    construct: str
    dose: float
    point: float
    low: float
    high: float
    method: str
    n_plates: int
    n_plates_contributing: int
    n_wells: int
    n_resamples: int
    estimable: bool
    confidence: float = 0.95
    is_reference: bool = False
    target: str = "geometric mean of positive within-plate dose/control activity ratios"
    resampling: str = "two_stage"
    coverage_status: str = "nominal, approximate; actual coverage requires repeated independent experiments"
    assumptions: tuple[str, ...] = (
        "independent exchangeable biological plates",
        "exchangeable wells within each plate and dose; dose and control wells unpaired",
        "positive per-plate ratios; exclusions condition the estimand on usable plates",
        "t-inflated percentile bootstrap is approximate and can double-count within-plate sampling variation",
    )

    @property
    def excludes_unity(self) -> bool:
        """Whether the interval separates this fold from 1.0.

        False whenever the interval is not estimable, which is the conservative
        reading: an absent interval is not evidence of an effect.
        """
        if not self.estimable or not np.isfinite([self.point, self.low, self.high]).all():
            return False
        if self.point <= 0 or self.low <= 0 or self.low > self.high:
            return False
        return bool(self.low > 1.0 or self.high < 1.0)

    def summary(self) -> str:
        if not self.estimable:
            return f"{self.construct} @ {self.dose:g}: {self.point:.2f} (no interval: {self.method})"
        verdict = "excludes 1.0" if self.excludes_unity else "includes 1.0"
        return (
            f"{self.construct} @ {self.dose:g}: {self.point:.2f} "
            f"[{self.low:.2f}, {self.high:.2f}] ({verdict}, n={self.n_plates_contributing}"
            + (f" of {self.n_plates}" if self.n_plates_contributing != self.n_plates else "")
            + " plates)"
        )


def _plate_fold(dosed: np.ndarray, control: np.ndarray) -> float:
    """One plate's fold: its own dosed mean over its own control mean."""
    denominator = float(np.mean(control))
    if not np.isfinite(denominator) or denominator < _LOG_EPS:
        return float("nan")
    return float(np.mean(dosed)) / denominator


def _combine(folds: np.ndarray) -> float:
    """Geometric mean of the per-plate folds.

    The fold is computed **within** each plate and only then combined, rather than as
    a ratio of pooled means. The difference is not cosmetic: a plate-level shift that
    multiplies a plate's dosed and control wells alike cancels exactly in a per-plate
    ratio, and does *not* cancel in a ratio of means across plates. Since inoculum
    density, medium batch and reader gain all act as exactly such a shift, pooling
    first would import batch variation into a quantity that is supposed to be free of
    it -- and would do so in a way that changes the point estimate, not merely its
    spread.

    This also matches the unit of analysis used on the design side: ``analysis/power.py``
    fits one slope per plate and tests the slopes, for the same reason.

    Geometric rather than arithmetic because folds are multiplicative: the mean of
    2x and 0.5x is 1x, not 1.25x.
    """
    usable = folds[np.isfinite(folds) & (folds > 0)]
    if usable.size == 0:
        return float("nan")
    return float(np.exp(np.mean(np.log(usable))))


def fold_change(
    readings: pd.DataFrame,
    construct: str,
    dose: float,
    value: str = "activity_late",
    control_dose: float = 0.0,
    n_resamples: int = 4000,
    confidence: float = 0.95,
    seed: int = 0,
) -> FoldChange:
    """Fold change over the zero-dose control, with a cluster bootstrap interval.

    Plates are resampled with replacement; wells are then resampled with replacement
    *within* each drawn plate. That is the two-stage nonparametric bootstrap for
    clustered data, and it reproduces both variance components -- between-plate and
    within-plate -- instead of only the second.

    The interval is taken on the log ratio and exponentiated back. A ratio is bounded
    below by zero and unbounded above, so a symmetric interval on the raw scale can
    reach negative folds; on the log scale it cannot, and the sampling distribution is
    much closer to symmetric there.

    Args:
        readings: Long-form table with at least ``plate``, ``dose_mM``, ``construct``
            and the column named by ``value``. A ``well`` column is used when present
            and its absence is recorded in the result's ``method``.
        construct: Which reporter.
        dose: The dosed condition.
        value: Column holding the quantity to ratio. ``activity_late`` is the
            dilution-corrected activity; ``naive_late`` is the uncorrected readout,
            which is worth computing alongside so the correction's size is visible.
        control_dose: The reference condition.
        n_resamples: Bootstrap draws.
        confidence: Nominal coverage.
        seed: Random seed.

    Raises:
        ValueError: if required columns are missing, or the construct or either dose
            is absent from the table.
    """
    if not 0 < confidence < 1:
        raise ValueError("confidence must be in (0, 1)")
    if not isinstance(n_resamples, (int, np.integer)) or n_resamples < 2:
        raise ValueError("n_resamples must be an integer of at least two")
    required = {"plate", "construct", "dose_mM", value}
    missing = required - set(readings.columns)
    if missing:
        raise ValueError(f"readings is missing required columns: {sorted(missing)}")

    subset = readings[readings.construct == construct]
    if subset.empty:
        raise ValueError(f"no rows for construct {construct!r}")
    dosed = subset[np.isclose(subset.dose_mM, dose)]
    control = subset[np.isclose(subset.dose_mM, control_dose)]
    if dosed.empty:
        raise ValueError(f"construct {construct!r} has no rows at dose {dose}")
    if control.empty:
        raise ValueError(f"construct {construct!r} has no rows at control dose {control_dose}")

    has_wells = "well" in readings.columns
    if dosed.plate.isna().any() or control.plate.isna().any():
        raise ValueError("plate identifiers must not be missing")
    if has_wells:
        if dosed.well.isna().any() or control.well.isna().any():
            raise ValueError("well identifiers must not be missing when wells are resampled")
        dosed = dosed.groupby(["plate", "well"], as_index=False)[value].mean()
        control = control.groupby(["plate", "well"], as_index=False)[value].mean()
    metadata = {"confidence": float(confidence),
                "resampling": "two_stage" if has_wells else "cluster_only"}
    shared = sorted(set(dosed.plate) & set(control.plate))
    if not shared:
        raise ValueError(
            f"construct {construct!r} has no plate carrying both dose {dose} and its "
            f"control; a fold change across different plates confounds the batch effect "
            f"with the dose"
        )
    dosed = dosed[dosed.plate.isin(shared)]
    control = control[control.plate.isin(shared)]

    dosed_by_plate = {p: g[value].to_numpy(dtype=float) for p, g in dosed.groupby("plate", sort=True)}
    control_by_plate = {p: g[value].to_numpy(dtype=float) for p, g in control.groupby("plate", sort=True)}
    per_plate = np.array([
        _plate_fold(dosed_by_plate[p], control_by_plate[p]) for p in shared
    ])
    point = _combine(per_plate)
    n_plates = len(shared)
    # THE REPLICATE COUNT THAT GOVERNS THE INTERVAL IS THE ONE THAT ENTERED THE ESTIMATE,
    # and until this line it was the one that entered the *experiment*.
    #
    # `_combine` takes the geometric mean over the per-plate folds that are finite and
    # positive; a plate whose recovered activity came out negative -- a dying culture
    # losing signal -- is dropped from it without a word. Everything downstream then used
    # `len(shared)`: the cluster bootstrap's t inflation, `MIN_PLATES_FOR_INTERVAL`, and
    # the `n_plates` a reader sees. So a fold resting on one plate was widened as though
    # it rested on four, which is backwards -- fewer clusters must widen it, not narrow it.
    #
    # It was invisible until the panel went from three plates to four, because the point
    # estimate does not move when a plate that contributes nothing is added while the t
    # quantile does: `t.ppf(0.975, 3)` is 3.18 against 4.30 at df = 2. AlteredYap1 at
    # 2.0 mM H2O2 in the widest windows went from [0.010, 1.349] to [0.013, 0.473] on that
    # alone -- a culture being killed by peroxide, turning into a resolvable call because
    # a fourth plate agreed that it was dying.
    contributing = int(np.count_nonzero(np.isfinite(per_plate) & (per_plate > 0)))
    n_wells = int(len(dosed) + len(control))
    if np.isclose(dose, control_dose) and np.isfinite(point) and point > 0:
        return FoldChange(
            construct, float(dose), 1.0, 1.0, 1.0,
            "normalization identity, not independent evidence of no response", n_plates,
            contributing, len(control), 0, True, confidence=float(confidence),
            is_reference=True, resampling="none", coverage_status="algebraic identity")

    if not (np.isfinite(point) and point > 0):
        # A FOLD WITH NO POINT ESTIMATE CANNOT HAVE AN INTERVAL AROUND IT, and without this
        # it got one. `_combine` takes a geometric mean over the per-plate folds that are
        # finite and positive, so a condition where every plate's recovered activity is
        # negative returns NaN -- and the bootstrap below then sailed past its own guard,
        # which only requires half the RESAMPLES to be positive, and returned
        # `estimable=True` with `point=nan` and two real-looking bounds.
        #
        # That is worse than a wrong number. `excludes_unity` reads `low > 1.0 or
        # high < 1.0` and knows nothing of the point, so AlteredYap1 at 2.0 mM H2O2 --
        # where the naive fold is -0.22 because the peroxide is killing the culture and
        # total fluorescence is falling -- came back as an interval of [0.003, 0.585]
        # entirely below 1.0, which reads as a confident finding of strong repression.
        #
        # The condition is real and worth reporting; what it is not is a fold change.
        #
        # Checked before the replicate count below, not after, because zero contributing
        # plates would trip that check too and would report the vaguer of the two facts.
        # "The geometric mean does not exist" is a stronger statement than "fewer than
        # three plates produced one", and it is the one a reader needs.
        return FoldChange(
            construct=construct, dose=float(dose), point=point,
            low=float("nan"), high=float("nan"),
            method=(
                "refused: no positive per-plate fold, so the geometric mean does not "
                "exist. A finite positive dose mean and a finite positive control mean "
                "are required; this is not evidence of induction or repression"
            ),
            n_plates=n_plates, n_plates_contributing=contributing,
            n_wells=n_wells, n_resamples=0, estimable=False, **metadata,
        )

    if contributing < MIN_PLATES_FOR_INTERVAL:
        detail = (f"{contributing} of {n_plates} biological replicate(s) produced a usable "
                  f"fold" if contributing != n_plates
                  else f"{n_plates} biological replicate(s)")
        return FoldChange(
            construct=construct, dose=float(dose), point=point,
            low=float("nan"), high=float("nan"),
            method=(
                f"refused: {detail}, "
                f"{MIN_PLATES_FOR_INTERVAL} needed for a cluster bootstrap"
            ),
            n_plates=n_plates, n_plates_contributing=contributing,
            n_wells=n_wells, n_resamples=0, estimable=False, **metadata,
        )

    rng = np.random.default_rng(seed)
    plates = np.asarray([plate for plate, fold in zip(shared, per_plate)
                         if np.isfinite(fold) and fold > 0], dtype=object)

    draws = np.empty(n_resamples, dtype=float)
    for i in range(n_resamples):
        picked = rng.choice(plates, size=contributing, replace=True)
        folds = np.empty(contributing, dtype=float)
        for j, label in enumerate(picked):
            pair = []
            for pool in (dosed_by_plate[label], control_by_plate[label]):
                pair.append(rng.choice(pool, size=pool.size, replace=True) if has_wells else pool)
            folds[j] = _plate_fold(pair[0], pair[1])
        draws[i] = _combine(folds) if np.isfinite(folds).all() and np.all(folds > 0) else np.nan

    finite = draws[np.isfinite(draws) & (draws > 0)]
    if finite.size < (1.0 + confidence) / 2.0 * n_resamples:
        return FoldChange(
            construct=construct, dose=float(dose), point=point,
            low=float("nan"), high=float("nan"),
            method=(f"refused: {n_resamples - finite.size} of {n_resamples} bootstrap ratios "
                    "are non-positive or undefined; discarding them exceeds the interval tail budget"),
            n_plates=n_plates, n_plates_contributing=contributing,
            n_wells=n_wells, n_resamples=int(finite.size), estimable=False, **metadata,
        )

    alpha = 1.0 - confidence
    log_draws = np.log(finite)
    lo, hi = np.quantile(log_draws, [alpha / 2.0, 1.0 - alpha / 2.0])

    # Percentile coverage falls off badly at small cluster counts. Widen by the ratio of
    # the t quantile to the normal one at this many clusters, which is the standard
    # small-sample correction and is honest about being approximate.
    from scipy import stats

    inflation = float(
        stats.t.ppf(1.0 - alpha / 2.0, df=max(contributing - 1, 1))
        / stats.norm.ppf(1.0 - alpha / 2.0)
    )
    centre = float(np.log(point)) if point > 0 else float(np.mean(log_draws))
    lo = centre + (lo - centre) * inflation
    hi = centre + (hi - centre) * inflation

    well_note = "wells resampled within plate" if has_wells else "plate means only, no well column"
    return FoldChange(
        construct=construct, dose=float(dose), point=point,
        low=float(np.exp(lo)), high=float(np.exp(hi)),
        method=(f"{metadata['resampling']} cluster bootstrap on log ratio, approximately t-inflated x{inflation:.2f} at "
                f"{contributing} contributing plate(s); {well_note}; nominal {confidence:.1%}, coverage not guaranteed"),
        n_plates=n_plates, n_plates_contributing=contributing,
        n_wells=n_wells, n_resamples=int(finite.size), estimable=True, **metadata,
    )


# ---------------------------------------------------------------------------
# equivalence
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EquivalenceVerdict:
    """Whether a fold is positively small, positively real, or simply unresolved.

    The three-way outcome mirrors the gates elsewhere in this package. "No induction
    at all" is a claim that needs evidence of its own, and a wide interval straddling
    the margin supports neither that claim nor its opposite.
    """

    construct: str
    dose: float
    verdict: Literal["EQUIVALENT", "DIFFERENT", "INCONCLUSIVE"]
    margin: tuple[float, float]
    fold: FoldChange
    method: str = "interval TOST against a prespecified practical-equivalence margin"

    @property
    def one_sided_alpha(self) -> float:
        return (1.0 - self.fold.confidence) / 2.0

    def summary(self) -> str:
        low, high = self.margin
        return (
            f"{self.construct} @ {self.dose:g}: {self.verdict} "
            f"against a no-induction margin of [{low:g}, {high:g}] "
            f"(interval TOST, one-sided alpha={self.one_sided_alpha:g}) -- {self.fold.summary()}"
        )


def equivalence(
    fold: FoldChange,
    margin: tuple[float, float] = (0.9, 1.1),
) -> EquivalenceVerdict:
    """Classify a fold change against a pre-named no-induction margin.

    Two one-sided tests, in interval form: the fold is *equivalent* to no induction
    when its whole confidence interval lies inside the margin, *different* when the
    interval excludes 1.0 entirely, and inconclusive otherwise.

    Naming the margin in advance is the load-bearing part. "Corrected fold 0.96, no
    induction at all" is only a result if 0.96 was going to count as no induction
    before the number was seen.

    Args:
        fold: Result of :func:`fold_change`.
        margin: Multiplicative bounds inside which a fold counts as no induction.

    Raises:
        ValueError: if the margin does not bracket 1.0, which would make the
            equivalence region unreachable by a null effect.
    """
    low, high = float(margin[0]), float(margin[1])
    if not np.isfinite([low, high]).all() or not (0 < low < 1.0 < high):
        raise ValueError(f"equivalence margin {margin} must bracket 1.0 with finite positive bounds")

    if (not fold.estimable or fold.is_reference
            or not np.isfinite([fold.point, fold.low, fold.high]).all()
            or min(fold.point, fold.low) <= 0 or fold.low > fold.high):
        return EquivalenceVerdict(fold.construct, fold.dose, "INCONCLUSIVE", (low, high), fold)
    if low <= fold.low and fold.high <= high:
        return EquivalenceVerdict(fold.construct, fold.dose, "EQUIVALENT", (low, high), fold)
    if fold.excludes_unity:
        return EquivalenceVerdict(fold.construct, fold.dose, "DIFFERENT", (low, high), fold)
    return EquivalenceVerdict(fold.construct, fold.dose, "INCONCLUSIVE", (low, high), fold)


# ---------------------------------------------------------------------------
# propagating measurement error into the recovered activity
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ActivityUncertainty:
    """Recovered promoter activity and the spread the inputs put on it."""

    times_h: np.ndarray
    activity: np.ndarray
    sigma: np.ndarray
    growth_rate_se: float
    reporter_cv: float
    method: str

    def band(self, z: float = 1.96) -> tuple[np.ndarray, np.ndarray]:
        """Pointwise interval. Pointwise, not simultaneous -- do not read it as a
        band that the whole curve stays inside."""
        return self.activity - z * self.sigma, self.activity + z * self.sigma


def _reporter_noise_cv(times_h: np.ndarray, reporter: np.ndarray, window_h: float) -> float:
    """Reader noise on the per-cell trace, from residuals about its own smooth.

    Estimated rather than asked for. A caller supplying a CV would be guessing, and the
    trace already carries the information: whatever the smooth does not explain at a
    scale shorter than the window is noise.
    """
    from scipy.signal import savgol_filter

    r = np.asarray(reporter, dtype=float)
    dt = float(np.median(np.diff(times_h)))
    window = max(int(round(window_h / dt)) | 1, 5)
    window = min(window, r.size if r.size % 2 == 1 else r.size - 1)
    smooth = savgol_filter(r, window_length=window, polyorder=min(3, window - 1))
    scale = float(np.median(np.abs(r)))
    if scale <= 0:
        return 0.0
    # Median absolute deviation: one bad timepoint should not set the noise level.
    residual = r - smooth
    mad = float(np.median(np.abs(residual - np.median(residual))))
    return float(1.4826 * mad / scale)


def activity_uncertainty(
    times_h: np.ndarray,
    reporter: SpecificFluorescence,
    optical_density: CorrectedOD,
    growth_rate: np.ndarray,
    kinetics: ReporterKinetics = ReporterKinetics(),
    window_h: float | None = None,
    n_draws: int = 400,
    method: Literal["delta", "monte-carlo"] = "monte-carlo",
    seed: int = 0,
) -> ActivityUncertainty:
    """Promoter activity with the error its two noisy inputs imply.

    ``promoter_activity`` computes ``dR/dt + (mu + k_deg) R``. Both ``R`` and ``mu``
    are estimates, and -- this is the part a naive propagation gets wrong -- they are
    **correlated**, because both are built from the same optical density trace:
    ``R = (RFU - background) / (OD - blank)`` and ``mu = d(ln OD)/dt``. An OD error
    pushes ``R`` down while pushing ``mu`` up, so treating them as independent
    misstates the variance in a direction that depends on which term dominates.

    The Monte-Carlo path handles that correctly by perturbing the *optical density*
    and letting both quantities move together, which is what physically happens. The
    delta path is a linearisation that assumes independence; it is much faster and is
    right to within a few percent when the growth-rate error is small, and it is
    documented here as the sweep-time approximation rather than the reported one.

    Args:
        times_h: Ascending time grid, hours.
        reporter: Per-cell reporter signal, ``(RFU - background) / (OD - blank)``.
        optical_density: The blank-corrected OD those values were built from. Required
            because the correlation cannot be reconstructed without it.
        growth_rate: Specific growth rate on the same grid.
        kinetics: Reporter constants.
        window_h: Smoothing window; defaults to the reporter module's own scaling.
        n_draws: Monte-Carlo draws. Ignored by the delta path.
        method: ``"monte-carlo"`` for reported numbers, ``"delta"`` for sweeps.
        seed: Random seed.

    Raises:
        ValueError: on shape mismatch or non-positive optical density.
    """
    t = np.asarray(times_h, dtype=float)
    r = require(reporter, SpecificFluorescence, name="reporter").array
    density = require(optical_density, CorrectedOD, name="optical_density")
    od = density.array
    mu = np.asarray(growth_rate, dtype=float)
    if not (t.shape == r.shape == od.shape == mu.shape):
        raise ValueError("times, reporter, optical density and growth rate must share a shape")
    if not np.all(od > 0):
        raise ValueError("optical density must be strictly positive")
    if t.ndim != 1 or t.size < 5 or not np.all(np.diff(t) > 0):
        raise ValueError("times must contain at least five increasing time points")
    if not np.isfinite(np.stack([t, r, od, mu])).all():
        raise ValueError("times, reporter, optical density and growth must be finite")
    if method == "monte-carlo" and (not isinstance(n_draws, (int, np.integer)) or n_draws < 2):
        raise ValueError("Monte Carlo uncertainty needs at least two integer draws")

    if window_h is None:
        window_h = default_activity_window_h(t)
    activity = promoter_activity(t, reporter, mu, kinetics, window_h=window_h)
    mu_se = growth_rate_uncertainty(t, density)
    cv = _reporter_noise_cv(t, r, window_h)

    if method == "delta":
        # d(activity)/d(mu) = R, and d(activity)/dR carries the derivative term; the
        # linearisation keeps the algebraic part and drops the smoothing operator's
        # effect on the noise, which is why it is the approximation.
        sigma = np.sqrt((r * mu_se) ** 2 + ((mu + kinetics.k_deg) * cv * np.abs(r)) ** 2)
        return ActivityUncertainty(t, activity, sigma, mu_se, cv, "delta method, inputs independent")

    if method != "monte-carlo":
        raise ValueError(f"unknown method {method!r}; use 'delta' or 'monte-carlo'")

    from ..growth import specific_growth_rate

    rng = np.random.default_rng(seed)
    background_free = r * od  # recover the numerator so OD can be perturbed coherently
    draws = np.empty((n_draws, t.size), dtype=float)
    od_cv = _reporter_noise_cv(t, od, window_h)
    numerator_cv = _reporter_noise_cv(t, background_free, window_h)
    for i in range(n_draws):
        od_i = od * (1.0 + rng.normal(0.0, od_cv, od.size)) if od_cv > 0 else od
        od_i = np.clip(od_i, 1e-9, None)
        numerator = background_free * (1.0 + rng.normal(0.0, numerator_cv, r.size))
        r_i = numerator / od_i
        # The growth rate must be re-estimated the way the real pipeline estimates it.
        # A raw finite difference on noisy log-OD is a far worse estimator than the
        # smoothed one, and using it here would charge the activity for noise that
        # ``specific_growth_rate`` never lets through -- inflating sigma several-fold.
        mu_i = specific_growth_rate(t, CorrectedOD(od_i), window_h=min(window_h, float(t[-1] - t[0]) / 2.0))
        draws[i] = promoter_activity(t, SpecificFluorescence(r_i), mu_i, kinetics, window_h=window_h)

    sigma = np.nanstd(draws, axis=0, ddof=1)
    return ActivityUncertainty(
        t, activity, sigma, mu_se, cv,
        f"monte-carlo, {n_draws} draws, OD perturbed so R and mu move together; "
        f"numerator CV {numerator_cv:.4g} estimated before OD division; "
        "independent numerator/OD measurement errors assumed",
    )


# ---------------------------------------------------------------------------
# how much the answer depends on the analysis window
# ---------------------------------------------------------------------------


def window_sensitivity(
    readings: pd.DataFrame,
    construct: str,
    dose: float,
    value_columns: tuple[str, ...] = ("activity_late", "naive_late"),
    control_dose: float = 0.0,
) -> pd.DataFrame:
    """Fold change under each available summary column, side by side.

    The reported numbers summarise the trace over ``t > 0.75 * duration``, a choice
    with no justification behind it. A fold that moves when the window moves is a
    property of the window. This reports the corrected and uncorrected folds together
    so the size of the dilution correction is visible next to the fold itself, which
    is the comparison the claim actually rests on.

    Args:
        readings: Long-form table as for :func:`fold_change`.
        construct: Which reporter.
        dose: The dosed condition.
        value_columns: Summary columns to compare.
        control_dose: The reference condition.
    """
    rows = []
    for column in value_columns:
        if column not in readings.columns:
            continue
        fold = fold_change(readings, construct, dose, value=column, control_dose=control_dose)
        rows.append({
            "construct": construct, "dose_mM": dose, "column": column,
            "fold": fold.point, "low": fold.low, "high": fold.high,
            "estimable": fold.estimable, "n_plates": fold.n_plates,
            "n_plates_contributing": fold.n_plates_contributing, "confidence": fold.confidence,
            "interval_target": fold.target, "interval_method": fold.method,
            "coverage_status": fold.coverage_status,
        })
    if not rows:
        raise ValueError(f"none of {value_columns} present in readings")
    return pd.DataFrame(rows)


def fold_change_coverage(n_plates, *, n_wells=6, true_fold=1.35, plate_cv=0.2,
                         fold_cv=0.15, well_cv=0.08, n_trials=100,
                         n_resamples=600, confidence=0.95, seed=0):
    from .power import PowerEstimate

    if any(not isinstance(n, (int, np.integer)) or n < 1 for n in (n_plates, n_wells, n_trials)):
        raise ValueError("plate, well and trial counts must be positive integers")
    if not np.isfinite([true_fold, plate_cv, fold_cv, well_cv]).all() or true_fold <= 0 or min(plate_cv, fold_cv, well_cv) < 0:
        raise ValueError("true_fold must be positive and noise scales finite and non-negative")
    rng = np.random.default_rng(seed)
    hits = issued = 0
    for _ in range(n_trials):
        level = np.exp(rng.normal(0.0, plate_cv, n_plates))
        fold = true_fold * np.exp(rng.normal(0.0, fold_cv, n_plates))
        means = level[:, None] * np.column_stack([np.ones(n_plates), fold])
        values = means[:, :, None] * np.exp(rng.normal(0.0, well_cv, (n_plates, 2, n_wells)))
        frame = pd.DataFrame({
            "plate": np.repeat(np.arange(n_plates), 2 * n_wells),
            "well": np.tile(np.arange(2 * n_wells), n_plates),
            "construct": "simulated", "dose_mM": np.tile(np.repeat([0.0, 1.0], n_wells), n_plates),
            "activity_late": values.ravel(),
        })
        result = fold_change(frame, "simulated", 1.0, n_resamples=n_resamples,
                             confidence=confidence, seed=int(rng.integers(0, 2**31 - 1)))
        issued += int(result.estimable)
        hits += int(result.estimable and result.low <= true_fold <= result.high)
    precision = PowerEstimate(hits, issued) if issued else None
    return {
        "n_plates": n_plates, "n_wells_per_dose": n_wells, "true_fold": true_fold,
        "plate_log_sd": plate_cv, "plate_fold_log_sd": fold_cv, "well_log_sd": well_cv,
        "n_trials": n_trials, "n_issued": issued, "n_covered": hits,
        "n_refused": n_trials - issued, "nominal_coverage": confidence,
        "actual_coverage": float(precision) if precision is not None else np.nan,
        "coverage_mc_low": precision.low if precision is not None else np.nan,
        "coverage_mc_high": precision.high if precision is not None else np.nan,
        "issued_and_covered_fraction": hits / n_trials, "n_resamples": n_resamples,
        "interval_target": "population geometric mean of within-plate dose/control ratios",
        "method": "two-stage cluster/well percentile bootstrap with approximate t inflation",
        "coverage_scope": "simulation under independent lognormal plate, plate-by-dose and well effects; not real-data calibration",
    }
