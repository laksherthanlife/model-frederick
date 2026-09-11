"""How many biological replicates to detect a real induction after growth correction.

Detection means the *recovered* promoter activity varies with dose, not that the raw
signal does -- growth inhibition alone guarantees the latter. Simulated at fixed growth
so the reporter reaches steady state, which is what the design requires anyway.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
from scipy import stats

from ..generator.culture import CultureParameters
from ..generator.design import Condition, collinearity
from ..growth import specific_growth_rate
from ..readings import CorrectedOD

__all__ = [
    "PowerEstimate",
    "detection_power",
    "discrimination_power",
    "power_curve",
    "ratio_discrimination_power",
    "ratio_growth_invariance",
    "replicates_needed",
]

_ALPHA = 0.05
_DURATION_H = 60.0
_N_TIMEPOINTS = 240

DEFAULT_SIMULATIONS = 200
"""Campaigns simulated per estimate when the caller does not ask for a precision.

Two hundred draws put a Monte Carlo standard error of about 3.5 percentage points on
an estimate near 0.5. That is fine for "is this design roughly usable" and far too
coarse to separate 95% from 99%, which is a comparison this module's output has been
used to make. Callers who need that separation should pass ``tolerance``.
"""


class PowerEstimate(float):
    """A power, carrying the precision of the simulation that produced it.

    Subclasses ``float`` so it compares and formats like the bare number it replaces,
    while making the Monte Carlo error reachable. A power quoted without one is a
    number whose last digit is noise: at the default simulation count two estimates
    can differ by seven percentage points and be the same power.

    Attributes:
        hits: Campaigns in which the effect was detected.
        n_simulations: Campaigns simulated.
        low: Lower Wilson score bound.
        high: Upper Wilson score bound.
    """

    __slots__ = ("hits", "n_simulations", "low", "high", "confidence", "method", "assumptions")

    def __new__(cls, hits: int, n_simulations: int, confidence: float = 0.95,
                method: str = "two-sided plate-slope t-test with positive direction required"):
        if not isinstance(n_simulations, (int, np.integer)) or n_simulations <= 0:
            raise ValueError("n_simulations must be positive and integer")
        if not isinstance(hits, (int, np.integer)) or not 0 <= hits <= n_simulations:
            raise ValueError("hits must be an integer between zero and n_simulations")
        if not 0 < confidence < 1:
            raise ValueError("confidence must be in (0, 1)")
        point = hits / n_simulations
        self = super().__new__(cls, point)
        self.hits = int(hits)
        self.n_simulations = int(n_simulations)
        self.confidence = float(confidence)
        self.method = method
        self.assumptions = ("independent simulated campaigns and biological plates",
                            "approximately Gaussian plate coefficients for the t-test",
                            "Wilson bounds quantify Monte Carlo error, not biological model validity")
        # Wilson score, not Wald. These estimates live at 0.0 and 1.0, where a Wald
        # interval has zero width and covers nothing.
        z = float(stats.norm.ppf(0.5 + confidence / 2.0))
        n = float(n_simulations)
        centre = (point + z * z / (2 * n)) / (1 + z * z / n)
        half = (z / (1 + z * z / n)) * np.sqrt(point * (1 - point) / n + z * z / (4 * n * n))
        self.low = float(max(0.0, centre - half))
        self.high = float(min(1.0, centre + half))
        return self

    @property
    def standard_error(self) -> float:
        """Monte Carlo standard error of the point estimate."""
        p = float(self)
        return float(np.sqrt(max(p * (1 - p), 0.0) / self.n_simulations))

    def summary(self) -> str:
        return (
            f"{float(self):.0%} [{self.low:.0%}, {self.high:.0%}] "
            f"({self.n_simulations} simulated campaigns)"
        )

    def __repr__(self) -> str:  # pragma: no cover - display only
        return f"PowerEstimate({float(self):.3f}, ci=[{self.low:.3f}, {self.high:.3f}], n={self.n_simulations})"


def _reporter_at(activity: float, growth: float, k_deg: float, t: float,
                 initial: float | None = None) -> float:
    """Closed form of dR/dt = k_synth - (mu+k_deg)R for constant inputs.

    Exact, and avoids an ODE solve per well; a power sweep runs tens of thousands.
    """
    loss = max(growth + k_deg, 1e-12)
    steady = activity / loss
    start = steady if initial is None else initial
    return float(steady + (start - steady) * np.exp(-loss * t))


def _scaled(params: CultureParameters, induction_fold: float) -> CultureParameters:
    """Reset the promoter peak to a chosen fold over basal, leaving growth untouched."""
    return replace(params, promoter_peak=params.promoter_basal * induction_fold)


def _expression_noise(params, log_noise):
    scale = float(np.exp(log_noise))
    return replace(params, promoter_basal=params.promoter_basal * scale,
                   promoter_peak=params.promoter_peak * scale)


def _positive_plate_test(coefficients):
    values = np.asarray(coefficients, dtype=float)
    if values.size < 2 or not np.isfinite(values).all():
        return False
    if np.std(values, ddof=1) < 1e-18:
        return bool(np.mean(values) > 1e-18)
    statistic, pvalue = stats.ttest_1samp(values, 0.0)
    return bool(pvalue < _ALPHA and statistic > 0)


def _observe(condition: Condition, params: CultureParameters, rng, reader_cv: float,
             growth_known: bool) -> tuple[float, float]:
    """One replicate: observed per-cell reporter, and the growth rate as it is known."""
    activity = params.promoter_activity_at(condition.dose_mM)
    growth = max(condition.growth_rate, 1e-6)
    reporter = _reporter_at(activity, growth, params.k_deg, _DURATION_H)
    signal = float(reporter * (1.0 + rng.normal(0.0, reader_cv)))
    if growth_known:
        return signal, growth
    # On a plate mu is estimated from a noisy OD trace, and that error propagates.
    window = np.linspace(0.0, 8.0, 50)
    od = 0.05 * np.exp(growth * window) * (1.0 + rng.normal(0.0, reader_cv, window.size))
    od = np.clip(od, 1e-6, None)
    estimated = float(np.nanmedian(specific_growth_rate(window, CorrectedOD(od), window_h=2.0)))
    return signal, max(estimated, 1e-6)


DEFAULT_WELL_CV = 0.052
"""Between-well biological spread within one plate, as a lognormal sigma.

A separate variance component from ``biological_cv``, and it must stay separate. The
plate-level term is shared by every condition on a plate, so it largely divides out of
a within-plate dose slope; the well-level term does not, and it lands directly on that
slope. Conflating them would make the model look adequate while the quantity that
actually limits power went unrepresented.

**Measured, not assumed.** From the zero-dose wells of the two NewProtocol plates that
carry a usable blank: three wells per construct per plate, differing only in which well
they are, which is exactly the quantity this models. The lognormal sigma of recovered
activity across those wells:

    AlteredYap1  0.038 / 0.083     NativeYap1  0.037 / 0.074
    UPRE1        0.006 / 0.052     UPRE2       0.052 / 0.057

Median 0.052 over eight construct-by-plate groups, range 0.006 to 0.083, and the plate
of 2026-08-03 is consistently the noisier of the two. Reproduce with the zero-dose rows
of ``outputs/sensor_characterisation.csv``.

Three wells per group makes each estimate noisy on its own; the median across eight is
what the constant takes. It is a floor in one specific sense worth stating: these are
wells of a shared overnight pipetted at one sitting, so anything that varies between
*sittings* is absent from it.

Independently corroborated: Kensy et al. 2009, the standard microtiter-plate validation
for parallel biomass and fluorescence measurement, reports under 5% standard deviation
between replicate wells of one clone.
"""


def _one_trial(
    params: CultureParameters, conditions, induction_fold: float,
    n_replicates: int, rng, reader_cv: float, biological_cv: float,
    growth_known: bool, well_cv: float = DEFAULT_WELL_CV,
) -> bool:
    """True when the per-plate dose slopes differ from zero across replicates.

    The unit of analysis is the plate, not the well. A biological replicate shares one
    batch effect across every condition on it, so treating wells as independent inflates
    the effective sample size by the number of conditions and overstates power severalfold.

    Two biological variance components, not one. The plate effect multiplies every
    condition on a plate alike and so mostly cancels out of that plate's dose slope --
    which is exactly why the plate is the right unit. What limits power is the *well*
    effect: each well is its own culture and its deviation lands on the slope directly.
    Model only the plate term and the dominant one is missing, so the replicate counts
    come back too small.
    """
    if n_replicates < 2:
        return False
    scaled = _scaled(params, induction_fold)
    slopes = []
    for _ in range(n_replicates):
        batch = _expression_noise(scaled, rng.normal(0.0, biological_cv))
        doses, recovered = [], []
        for condition in conditions:
            well = _expression_noise(batch, rng.normal(0.0, well_cv)) if well_cv > 0 else batch
            signal, growth = _observe(condition, well, rng, reader_cv, growth_known)
            # Steady-state inversion of dR/dt = k_synth - (mu + k_deg) R. The loss term
            # is the sum, so k_deg has to be carried: dropping it agrees with
            # reporter.promoter_activity only for a perfectly stable reporter, and every
            # construct in generator/literature.py currently has k_deg = 0.0.
            recovered.append(signal * (growth + well.k_deg))
            doses.append(condition.dose_mM)
        doses = np.asarray(doses, dtype=float)
        recovered = np.asarray(recovered, dtype=float)
        if np.std(doses) < 1e-15 or np.std(recovered) < 1e-18:
            return False
        slopes.append(stats.linregress(np.log1p(doses), recovered).slope)
    return _positive_plate_test(slopes)


def detection_power(
    params: CultureParameters,
    conditions,
    induction_fold: float,
    n_replicates: int,
    n_simulations: int = DEFAULT_SIMULATIONS,
    reader_cv: float = 0.02,
    biological_cv: float = 0.12,
    growth_known: bool = False,
    seed: int | None = None,
    well_cv: float = DEFAULT_WELL_CV,
    tolerance: float | None = None,
    max_simulations: int = 20_000,
    confidence: float = 0.95,
) -> PowerEstimate:
    """Fraction of simulated campaigns in which a true induction is detected.

    Returns a :class:`PowerEstimate`, which behaves as the float it always was and
    additionally carries the Monte Carlo interval. That interval is not decoration: at
    the default simulation count the estimate is good to about three percentage points,
    and this function's output has been used to compare designs that differ by less.

    Args:
        params: Strain kinetics; the promoter peak is overridden by ``induction_fold``.
        conditions: Planned design from :mod:`ystwin.generator.design`.
        induction_fold: True promoter induction at saturating dose. 1.0 means none,
            so the returned value is the false-positive rate.
        n_replicates: Biological replicates per condition.
        n_simulations: Simulated campaigns, when ``tolerance`` is not given.
        reader_cv: Per-reading noise.
        biological_cv: Between-plate spread in basal promoter activity.
        growth_known: True in a chemostat, where the dilution rate sets mu exactly. On a
            plate mu is estimated from a noisy OD trace and that error propagates.
        seed: Random seed.
        well_cv: Between-well spread within a plate. See :data:`DEFAULT_WELL_CV`.
        tolerance: If given, keep simulating until the Wilson interval is narrower
            than this, or ``max_simulations`` is reached. Ask for the precision the
            comparison needs rather than discovering afterwards that it was absent.
        max_simulations: Cap on the adaptive loop.
    """
    conditions = tuple(conditions)
    if not conditions or any(not np.isfinite([c.dose_mM, c.growth_rate]).all()
                             or c.dose_mM < 0 or c.growth_rate <= 0 for c in conditions):
        raise ValueError("conditions must have finite non-negative doses and positive growth")
    if not isinstance(n_replicates, (int, np.integer)) or n_replicates < 1:
        raise ValueError("n_replicates must be a positive integer")
    if not np.isfinite([reader_cv, biological_cv, well_cv, induction_fold]).all() or min(reader_cv, biological_cv, well_cv) < 0 or induction_fold <= 0:
        raise ValueError("noise must be non-negative and induction_fold positive, all finite")
    PowerEstimate(0, n_simulations, confidence)
    if tolerance is not None:
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("tolerance must be positive and finite")
        if not isinstance(max_simulations, (int, np.integer)) or max_simulations < 1:
            raise ValueError("max_simulations must be a positive integer")
        n_simulations = min(n_simulations, max_simulations)
    rng = np.random.default_rng(seed)

    def run(n: int) -> int:
        return sum(
            _one_trial(params, conditions, induction_fold, n_replicates, rng,
                       reader_cv, biological_cv, growth_known, well_cv)
            for _ in range(n)
        )

    if tolerance is None:
        return PowerEstimate(run(n_simulations), n_simulations, confidence)

    looks = 1 + int(np.ceil(np.log2(max_simulations / n_simulations)))
    look_confidence = 1 - (1 - confidence) / looks
    hits, drawn = run(n_simulations), n_simulations
    method = "two-sided plate-slope t-test; Wilson Monte Carlo intervals Bonferroni-adjusted over adaptive looks"
    estimate = PowerEstimate(hits, drawn, look_confidence, method)
    while estimate.high - estimate.low > tolerance and drawn < max_simulations:
        batch = min(drawn, max_simulations - drawn)
        hits += run(batch)
        drawn += batch
        estimate = PowerEstimate(hits, drawn, look_confidence, method)
    return estimate


def replicates_needed(
    params: CultureParameters,
    conditions,
    induction_fold: float,
    target_power: float = 0.8,
    max_replicates: int = 24,
    n_simulations: int = DEFAULT_SIMULATIONS,
    seed: int | None = None,
    require_confidence: bool = True,
    **kwargs,
) -> int | None:
    """Smallest replicate count reaching ``target_power``, or None within the cap.

    Walking upward and stopping at the first point estimate to clear the target is a
    biased-low estimator: with Monte Carlo noise on every step, the walk stops on the
    first *lucky* draw rather than the first adequate design, and the bias grows as the
    simulation count falls. Since the answer is the number of cultures somebody then
    goes and grows, erring low is the expensive direction.

    By default the *lower* confidence bound must clear the target, so the returned
    count is one the simulation actually supports.

    Args:
        params: Strain kinetics.
        conditions: Planned design.
        induction_fold: True induction at saturating dose.
        target_power: Power to reach.
        max_replicates: Cap on the search.
        n_simulations: Campaigns per candidate count.
        seed: Random seed.
        require_confidence: Require the Wilson lower bound to clear ``target_power``
            rather than the point estimate. Set False for the older, optimistic
            behaviour.
        **kwargs: Passed through to :func:`detection_power`.
    """
    if not 0 < target_power < 1:
        raise ValueError("target_power must be in (0, 1)")
    if not isinstance(max_replicates, (int, np.integer)) or max_replicates < 1:
        raise ValueError("max_replicates must be a positive integer")
    conditions = tuple(conditions)
    confidence = kwargs.pop("confidence", 0.95)
    if not 0 < confidence < 1:
        raise ValueError("confidence must be in (0, 1)")
    if require_confidence:
        confidence = 1 - (1 - confidence) / max_replicates
    for n in range(1, max_replicates + 1):
        power = detection_power(params, conditions, induction_fold, n,
                                n_simulations=n_simulations, seed=seed, confidence=confidence, **kwargs)
        reached = power.low if require_confidence else float(power)
        if reached >= target_power:
            return n
    return None


def power_curve(
    params: CultureParameters,
    conditions,
    induction_folds,
    replicate_counts,
    n_simulations: int = DEFAULT_SIMULATIONS,
    seed: int | None = None,
    **kwargs,
) -> pd.DataFrame:
    """Power across effect sizes and replicate counts for one design."""
    coll = collinearity(conditions) if len(conditions) >= 3 else float("nan")
    rows = []
    for fold in induction_folds:
        for n in replicate_counts:
            estimate = detection_power(params, conditions, fold, n,
                                       n_simulations=n_simulations, seed=seed, **kwargs)
            rows.append({
                "induction_fold": fold, "n_replicates": n,
                "power": float(estimate),
                "power_low": estimate.low, "power_high": estimate.high,
                "n_simulations": estimate.n_simulations,
                "mc_confidence": estimate.confidence, "test_method": estimate.method,
                "test_alpha": _ALPHA, "independent_unit": "biological plate",
                "inference_target": "probability of positive rejection; non-rejection is not equivalence",
                "n_conditions": len(conditions), "collinearity": coll,
            })
    return pd.DataFrame(rows)


def _plate_coefficients(dose, growth, signal, plate, include_growth=True):
    coefficients = []
    for label in np.unique(plate):
        mask = plate == label
        columns = [np.ones(mask.sum()), np.log1p(dose[mask])]
        if include_growth:
            columns.append(np.log(np.clip(growth[mask], 1e-9, None)))
        design = np.column_stack(columns)
        target = np.log(np.clip(signal[mask], 1e-30, None))
        fitted, _, rank, _ = np.linalg.lstsq(design, target, rcond=None)
        if rank < design.shape[1] or len(target) <= design.shape[1]:
            return []
        coefficients.append(float(fitted[1]))
    return coefficients


def _discrimination_trial(
    params: CultureParameters, conditions, induction_fold: float,
    n_replicates: int, rng, reader_cv: float, biological_cv: float,
    growth_known: bool, well_cv: float = DEFAULT_WELL_CV,
) -> bool:
    """True when dose keeps a significant coefficient with growth rate also in the model."""
    if n_replicates < 2:
        return False
    scaled = _scaled(params, induction_fold)
    rows = []
    for plate in range(n_replicates):
        batch = _expression_noise(scaled, rng.normal(0.0, biological_cv))
        for condition in conditions:
            well = _expression_noise(batch, rng.normal(0.0, well_cv)) if well_cv > 0 else batch
            signal, growth = _observe(condition, well, rng, reader_cv, growth_known)
            rows.append((condition.dose_mM, growth, signal, plate))
    dose, growth, signal, plate = (np.asarray(x, dtype=float) for x in zip(*rows))
    return _positive_plate_test(_plate_coefficients(dose, growth + params.k_deg, signal, plate))


def discrimination_power(
    params: CultureParameters,
    conditions,
    induction_fold: float,
    n_replicates: int,
    n_simulations: int = DEFAULT_SIMULATIONS,
    reader_cv: float = 0.02,
    biological_cv: float = 0.12,
    growth_known: bool = False,
    seed: int | None = None,
    well_cv: float = DEFAULT_WELL_CV,
) -> PowerEstimate:
    """Fraction of campaigns where dose survives with growth rate also in the model.

    This is what a decoupled design buys. Detecting a dose slope is easy on any design;
    showing the signal is driven by dose *rather than* by the growth rate the dose caused
    requires the two to vary independently, which collinearity forbids.
    """
    rng = np.random.default_rng(seed)
    hits = sum(
        _discrimination_trial(params, conditions, induction_fold, n_replicates, rng,
                              reader_cv, biological_cv, growth_known, well_cv)
        for _ in range(n_simulations)
    )
    return PowerEstimate(hits, n_simulations)


def _ratio_trial(
    params: CultureParameters, conditions, induction_fold: float, n_replicates: int,
    rng, reader_cv: float, biological_cv: float, stress_spec, reference_spec,
    unmixing_cv: float, reference_responds: bool, include_growth: bool,
) -> bool:
    """True when dose keeps a significant coefficient in the co-expressed ratio."""
    from ..generator.unmixing import ratio_series

    if n_replicates < 2:
        return False
    scaled = _scaled(params, induction_fold)
    rows = []
    for plate in range(n_replicates):
        batch_effect = float(np.exp(rng.normal(0.0, biological_cv)))
        for condition in conditions:
            activity = scaled.promoter_activity_at(condition.dose_mM) * batch_effect
            reference = scaled.promoter_activity_at(
                condition.dose_mM if reference_responds else 0.0
            ) * batch_effect
            growth = max(condition.growth_rate, 1e-6)
            ratio = ratio_series(activity, reference, growth, stress_spec, reference_spec)
            noise = rng.normal(0.0, np.hypot(reader_cv, unmixing_cv))
            rows.append((condition.dose_mM, growth, ratio * (1.0 + noise), plate))
    dose, growth, ratio, plate = (np.asarray(x, dtype=float) for x in zip(*rows))
    return _positive_plate_test(_plate_coefficients(dose, growth, ratio, plate, include_growth))


def ratio_discrimination_power(
    params: CultureParameters,
    conditions,
    induction_fold: float,
    n_replicates: int,
    stress_spec,
    reference_spec,
    n_simulations: int = DEFAULT_SIMULATIONS,
    reader_cv: float = 0.02,
    biological_cv: float = 0.12,
    unmixing_cv: float = 0.05,
    reference_responds: bool = False,
    include_growth: bool = False,
    seed: int | None = None,
) -> PowerEstimate:
    """Attribution power when a constitutive reference is co-expressed in the same cell.

    Growth cancels in the ratio, so collinearity between dose and growth stops mattering.
    Mismatched maturation between the two fluorophores and unmixing error are what remain.

    Args:
        stress_spec, reference_spec: Fluorophores carrying the sensor and the reference.
        unmixing_cv: Extra multiplicative error from spectral unmixing.
        reference_responds: If True the reference is itself stress-responsive, which
            divides the signal out. Present so that failure mode is testable.
        include_growth: Put growth back in the model. Off by default because mu cancels
            in the ratio by construction, so the term has no true effect -- and on a
            collinear design it still inflates the dose standard error, costing the very
            power the ratio was meant to buy. Verify the cancellation once with
            :func:`ratio_growth_invariance` rather than carrying the term every time.
    """
    rng = np.random.default_rng(seed)
    hits = sum(
        _ratio_trial(params, conditions, induction_fold, n_replicates, rng, reader_cv,
                     biological_cv, stress_spec, reference_spec, unmixing_cv,
                     reference_responds, include_growth)
        for _ in range(n_simulations)
    )
    return PowerEstimate(hits, n_simulations)


def ratio_growth_invariance(
    stress_activity: float,
    reference_activity: float,
    stress_spec,
    reference_spec,
    growth_rates=(0.40, 0.25, 0.10, 0.05),
) -> float:
    """Largest relative departure of the ratio across growth rates; 0 is exact.

    The check that licenses dropping growth from the model. It needs growth to vary
    independently of dose, which is what a nutrient axis provides -- once, as a
    validation, rather than in every experiment.
    """
    from ..generator.unmixing import ratio_series

    values = np.array([
        ratio_series(stress_activity, reference_activity, mu, stress_spec, reference_spec)
        for mu in growth_rates
    ])
    return float(np.max(np.abs(values - values.mean())) / values.mean())
