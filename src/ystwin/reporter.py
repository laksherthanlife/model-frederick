"""Reporter as a dynamic state: dR/dt = k_synth - (mu + k_deg) R, with optional maturation.

Dilution sets the gain. At quasi-steady state R -> k_synth/(mu + k_deg), so for a stable FP a
slowing culture shows rising RFU/OD with no change in promoter activity at all. Since stress
slows growth, a latent state fitted to raw RFU/OD is partly a measurement of 1/mu.
``promoter_activity`` inverts both dilution and maturation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .readings import CorrectedOD, CorrectedRFU, SpecificFluorescence, require
from scipy.integrate import solve_ivp
from scipy.signal import savgol_filter

__all__ = [
    "ReporterKinetics",
    "TARGET_ACTIVITY_WINDOW_H",
    "default_activity_window_h",
    "naive_specific_fluorescence",
    "promoter_activity",
    "promoter_activity_from_total",
    "simulate_reporter",
]


@dataclass(frozen=True)
class ReporterKinetics:
    """Per-reporter constants, in 1/h.

    Args:
        k_deg: Active degradation rate of the mature protein. ~0 for a stable FP
            such as GFP or mCitrine; large for a degron-tagged reporter.
        k_mat: Chromophore maturation rate. ``None`` means instantaneous.
        k_deg_immature: Degradation of the immature species; defaults to ``k_deg``.
    """

    k_deg: float = 0.0
    k_mat: float | None = None
    k_deg_immature: float | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.k_deg) or self.k_deg < 0:
            raise ValueError("k_deg must be finite and non-negative")
        if self.k_mat is not None and (not np.isfinite(self.k_mat) or self.k_mat <= 0):
            raise ValueError("maturation rate k_mat must be finite and positive, or None for instantaneous")
        if self.k_deg_immature is not None:
            if not np.isfinite(self.k_deg_immature) or self.k_deg_immature < 0:
                raise ValueError("k_deg_immature must be finite and non-negative")
            if self.k_mat is None:
                raise ValueError("k_deg_immature requires a maturation rate k_mat")

    @property
    def has_maturation(self) -> bool:
        return self.k_mat is not None

    @property
    def immature_loss(self) -> float:
        return self.k_deg if self.k_deg_immature is None else self.k_deg_immature


def naive_specific_fluorescence(rfu: CorrectedRFU,
                                optical_density: CorrectedOD) -> SpecificFluorescence:
    """RFU/OD -- the conventional readout, retained so the confound stays visible.

    Returns a :class:`~ystwin.readings.SpecificFluorescence` rather than a bare array, so
    that the one legitimate consumer of this quantity -- :func:`promoter_activity`, which
    divides the dilution back out -- receives it already labelled. The confound this
    function exists to demonstrate is that the ratio carries a ``1/mu`` gain; a type that
    said only "float array" would let it be read as an activity, which is the mistake.

    Raises:
        TypeError: if either channel is uncorrected or bare.
        ValueError: if any density is non-positive.
    """
    return require(rfu, CorrectedRFU, name="rfu").per(optical_density)


def _check_grid(t: np.ndarray, *series: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    if t.ndim != 1:
        raise ValueError("times must be a one-dimensional grid")
    if t.size < 5:
        raise ValueError("need at least 5 timepoints")
    if not np.all(np.isfinite(t)):
        raise ValueError("times must be finite")
    if np.any(np.diff(t) <= 0):
        raise ValueError("times must be strictly increasing")
    for s in series:
        values = np.asarray(s, dtype=float)
        if values.shape != t.shape:
            raise ValueError("all series must share the time grid's shape")
        if not np.all(np.isfinite(values)):
            raise ValueError("all series must be finite")
    return t


def simulate_reporter(
    times_h: np.ndarray,
    k_synth: np.ndarray,
    growth_rate: np.ndarray,
    kinetics: ReporterKinetics = ReporterKinetics(),
    r0: float | None = None,
    i0: float | None = None,
) -> np.ndarray:
    """Forward-simulate mature per-cell reporter concentration.

    Args:
        times_h: Ascending time grid, hours.
        k_synth: Promoter activity (reporter units per cell per hour).
        growth_rate: Specific growth rate ``mu`` on the same grid.
        kinetics: Reporter constants.
        r0: Initial mature concentration. Defaults to the quasi-steady state.
        i0: Initial immature concentration. Defaults to its quasi-steady state.
    """
    t = _check_grid(times_h, k_synth, growth_rate)
    k_syn = np.asarray(k_synth, dtype=float)
    mu = np.asarray(growth_rate, dtype=float)
    if np.any(k_syn < 0):
        raise ValueError("k_synth must be non-negative")
    for name, value in (("r0", r0), ("i0", i0)):
        if value is not None and (not np.isfinite(value) or value < 0):
            raise ValueError(f"{name} must be finite and non-negative")
    if not kinetics.has_maturation and i0 is not None:
        raise ValueError("i0 requires a maturing reporter; it cannot be used with k_mat=None")
    if r0 is None and mu[0] + kinetics.k_deg <= 0:
        raise ValueError("r0 is required when mature reporter has no positive initial loss rate")
    if (kinetics.has_maturation and i0 is None
            and mu[0] + float(kinetics.k_mat) + kinetics.immature_loss <= 0):
        raise ValueError("i0 is required when immature reporter has no positive initial loss rate")

    def synth(x): return float(np.interp(x, t, k_syn))
    def mu_at(x): return float(np.interp(x, t, mu))

    if kinetics.has_maturation:
        k_mat = float(kinetics.k_mat)
        if i0 is None:
            i0 = k_syn[0] / (mu[0] + k_mat + kinetics.immature_loss)
        if r0 is None:
            r0 = k_mat * i0 / (mu[0] + kinetics.k_deg)

        def rhs(x, y):
            imm, mat = y
            m = mu_at(x)
            return [
                synth(x) - (m + k_mat + kinetics.immature_loss) * imm,
                k_mat * imm - (m + kinetics.k_deg) * mat,
            ]

        y0 = [i0, r0]
    else:
        if r0 is None:
            r0 = k_syn[0] / (mu[0] + kinetics.k_deg)

        def rhs(x, y):
            return [synth(x) - (mu_at(x) + kinetics.k_deg) * y[0]]

        y0 = [r0]

    sol = solve_ivp(
        rhs, (t[0], t[-1]), y0, t_eval=t, method="LSODA", rtol=1e-9, atol=1e-12, max_step=float(np.diff(t).max()),
    )
    if not sol.success:
        raise RuntimeError(f"reporter integration failed: {sol.message}")
    return sol.y[-1]


def promoter_activity_from_total(
    times_h: np.ndarray,
    total_signal: CorrectedRFU,
    biomass: CorrectedOD,
    kinetics: ReporterKinetics = ReporterKinetics(),
    window_h: float | None = None,
    polyorder: int = 3,
) -> np.ndarray:
    """Promoter activity from total fluorescence, without needing a growth rate.

    The same quantity :func:`promoter_activity` returns. Writing the per-cell
    concentration as ``R = F / X`` gives ``dR/dt = (dF/dt)/X - R mu``, and substituting
    into ``k_synth = dR/dt + (mu + k_deg) R`` cancels mu exactly, for any ``k_deg``:

        k_synth = (dF/dt)/X + k_deg (F/X)

    So mu -- itself a smoothed derivative of a noisy optical trace -- is removed rather
    than estimated, and one numerical derivative is taken instead of two. Measured cost
    of the difference is in ``docs/FINDINGS.md``.

    Args:
        times_h: Ascending time grid, hours.
        total_signal: Background-corrected TOTAL fluorescence, not per cell.
        biomass: Blank-corrected optical density, same grid, strictly positive. Typed as
            a density because every caller in this repository passes one; a true dry
            weight in g/L would divide to different units and has no type here, because
            ``gdcw_per_od`` is the factor between them and nothing has measured it -- see
            :class:`~ystwin.readings.SpecificFluorescence`.
        kinetics: Reporter constants. Maturation is not handled here -- see below.
        window_h: Smoothing window; defaults to :func:`default_activity_window_h`.
        polyorder: Local polynomial order for the derivative.

    Raises:
        ValueError: on shape mismatch, non-positive biomass, or maturation kinetics.
            A maturation step needs the immature pool, which the cancellation above
            does not reach; use :func:`promoter_activity` for that case rather than
            silently returning the mature-only answer.
    """
    signal = require(total_signal, CorrectedRFU, name="total_signal").array
    x = require(biomass, CorrectedOD, name="biomass").array
    t = _check_grid(times_h, signal, x)
    if not np.all(x > 0):
        raise ValueError("biomass must be strictly positive; blank-correct first")
    if kinetics.has_maturation:
        raise ValueError(
            "maturation is not recoverable from the total-signal form: cancelling mu "
            "loses the immature pool. Use promoter_activity for a maturing reporter."
        )
    if window_h is None:
        window_h = default_activity_window_h(t)
    return _smooth_derivative(t, signal, window_h, polyorder) / x + kinetics.k_deg * signal / x


TARGET_ACTIVITY_WINDOW_H = 2.0
"""Window the inversion aims for, in hours, before the run's own floor and cap apply.

An ABSOLUTE duration, and that is the correction. This rule used to be purely relative --
one sixth of the run -- and `analysis/estimator_accuracy.py` shows that is the wrong
scaling: the recovered activity's error is governed by the window's duration in hours and
not by what fraction of the run it is. Two geometries six times apart in density agree
whenever their windows match in hours and disagree when they match in points:

    window 0.33 h   9.8% at 145 pts / 24 h    8.2% at 25 pts / 4 h
    window 1.17 h   5.4%                      4.5%
    window 2.17 h   2.3%                      2.0%
    window 3 POINTS 9.8% at 145 pts / 24 h    2.1% at 25 pts / 24 h   -- points do not

So a fraction of the run is right only at the run length it was tuned on. On the committed
plates -- 25 points over 4.00 h -- one sixth is 0.667 h, which collapses onto the floor and
lands at 8% where 1.33 h gives 3.2%. See `outputs/estimator_accuracy.csv`.

Two hours rather than four, though four scores better on a long run: at 4.00 h of plate the
cap below binds at 1.33 h either way, and on a 24 h run ``duration / 6`` is already 4 h and
outranks this. So the constant only ever decides SHORT runs, which is exactly where the old
rule failed, and it leaves every longer geometry where it was.
"""


def default_activity_window_h(times_h: np.ndarray) -> float:
    """Smoothing window for the activity inversion, scaled to the run and to the assay.

    The window trades noise rejection against bias. Measured by
    `analysis/estimator_accuracy.py` on induction pulses simulated at the committed plates'
    own geometry, growth and reader noise -- 25 points over 4.00 h, mu falling 0.563 to
    0.217 /h, 1.1% multiplicative per read -- the median relative error of the recovered
    activity, and the standard deviation of a recovered 1.5x FOLD, were:

        window    0.667 h   1.000 h   1.333 h   2.000 h
        activity     8.0%      4.6%      3.2%      2.0%
        fold sd    0.1157    0.0948    0.0786    0.0670

    **The fold's point estimate is accurate at every one of them** -- under 0.6% low
    throughout -- because a fold is a ratio of two wells on one plate and whatever the
    smoother does to the numerator it does to the denominator. What the window buys is
    precision, not accuracy, and that is what reaches an interval.

    The rule is therefore the larger of a sixth of the run and
    :data:`TARGET_ACTIVITY_WINDOW_H`, capped at a third of the run to limit bias against a
    real transient, with a sampling-interval floor that outranks the cap: an over-wide
    window is merely biased, while a window spanning fewer points than the polynomial order
    is not a fit at all.

    The table this docstring used to carry -- ``0.5 h -> 14.1% ... 12 h -> 1.6%`` -- was
    produced by nothing and checked by nothing, and was measured over a 24 h run at a
    density the plates are not at. The numbers above come from a module and a committed
    table, and `tests/test_estimator_accuracy.py` fails if they drift.
    """
    t = _check_grid(times_h)
    duration = float(t[-1] - t[0])
    dt = float(np.median(np.diff(t)))
    floor = 4.0 * dt
    return float(max(floor, min(max(duration / 6.0, TARGET_ACTIVITY_WINDOW_H, floor),
                                duration / 3.0)))


def _smooth_derivative(t: np.ndarray, y: np.ndarray, window_h: float, polyorder: int) -> np.ndarray:
    if not np.isfinite(window_h) or window_h <= 0:
        raise ValueError("window_h must be finite and positive")
    if (isinstance(polyorder, (bool, np.bool_))
            or not isinstance(polyorder, (int, np.integer)) or not 1 <= polyorder < t.size):
        raise ValueError("polyorder must be a positive integer smaller than the number of timepoints")
    intervals = np.diff(t)
    dt = float(np.median(intervals))
    uniform = np.allclose(intervals, dt, rtol=1e-7, atol=1e-12)
    if uniform:
        window = min(int(round(min(window_h, float(t[-1] - t[0])) / dt)), y.size)
        if window % 2 == 0:
            window += 1
        window = max(window, polyorder + 2)
        if window % 2 == 0:
            window += 1
        window = min(window, y.size if y.size % 2 == 1 else y.size - 1)
        if polyorder >= window:
            raise ValueError("polyorder requires a wider window than the time grid supports")
        return savgol_filter(y, window_length=window, polyorder=polyorder, deriv=1, delta=dt)

    minimum = polyorder + 2
    if minimum % 2 == 0:
        minimum += 1
    minimum = min(minimum, y.size)
    width = min(float(window_h), float(t[-1] - t[0]))
    derivative = np.empty_like(y, dtype=float)
    for index, centre in enumerate(t):
        lower = max(float(t[0]), min(float(centre) - width / 2.0, float(t[-1]) - width))
        left = int(np.searchsorted(t, lower, side="left"))
        right = int(np.searchsorted(t, lower + width, side="right"))
        while right - left < minimum:
            if left > 0 and (right == t.size or centre - t[left - 1] <= t[right] - centre):
                left -= 1
            else:
                right += 1
        offsets = t[left:right] - centre
        scale = float(np.max(np.abs(offsets)))
        design = np.polynomial.polynomial.polyvander(offsets / scale, polyorder)
        coefficients, _, rank, _ = np.linalg.lstsq(design, y[left:right] - y[index], rcond=None)
        if rank < polyorder + 1:
            raise ValueError("time grid cannot support the requested local polynomial derivative")
        derivative[index] = coefficients[1] / scale
    return derivative


def promoter_activity(
    times_h: np.ndarray,
    reporter: SpecificFluorescence,
    growth_rate: np.ndarray,
    kinetics: ReporterKinetics = ReporterKinetics(),
    window_h: float | None = None,
    polyorder: int = 3,
) -> np.ndarray:
    """Invert dilution (and maturation) to recover promoter activity from signal.

    This is the quantity a latent stress state may legitimately be conditioned on.
    Raw ``RFU/OD`` is not: it carries the ``1/mu`` gain.

    Args:
        times_h: Ascending time grid, hours.
        reporter: Background-corrected fluorescence divided by blank-corrected density,
            on the same grid. This argument's docstring read "per-cell mature reporter
            concentration" until 2026-08-31 while every caller passed RFU/OD, so what
            comes back is RFU/OD/h and not a molar rate. The type now says which of the
            two it is; ``docs/ARCHITECTURE.md`` §6.3 is where the discrepancy was found.
        growth_rate: Specific growth rate ``mu``, on the same grid.
        kinetics: Reporter constants.
        window_h: Smoothing window for the numerical derivative. Defaults to
            :func:`default_activity_window_h`.
        polyorder: Local polynomial order for the derivative.
    """
    r = require(reporter, SpecificFluorescence, name="reporter").array
    mu = np.asarray(growth_rate, dtype=float)
    t = _check_grid(times_h, r, mu)
    if window_h is None:
        window_h = default_activity_window_h(t)

    dr = _smooth_derivative(t, r, window_h, polyorder)
    balance = dr + (mu + kinetics.k_deg) * r
    if not kinetics.has_maturation:
        return balance

    # Mature balance gives the immature pool; differentiate once more for synthesis.
    immature = balance / float(kinetics.k_mat)
    di = _smooth_derivative(t, immature, window_h, polyorder)
    return di + (mu + float(kinetics.k_mat) + kinetics.immature_loss) * immature
