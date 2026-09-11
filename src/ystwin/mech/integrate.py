"""The stiff driver, the window-aware reduction set, and the tau/T audit that licenses it.

WHY A DRIVER AT ALL. :func:`ystwin.mech.state.catalogue_span_hours` measures the spread of this
architecture's time constants at **1.11e8** -- 6.25 ms for Ypd1 -> Ssk1 phosphotransfer (MEASURED
160 /s, Janiak-Spens 2005) against 193 h for plasmid loss at the fed-batch setpoint. An explicit
Runge-Kutta chooses its step from the fastest of those and the observable only ever sees the
slowest: ``solve_ivp(method="RK45")`` has a real-axis stability limit near |lambda|h <~ 3, so
lambda = 160 /s caps the step at 0.019 s and a 4.14 h read costs 7.7e5 steps for an answer that
changes on a timescale of hours. That is the textbook definition of stiffness and it has exactly
two legitimate answers -- an implicit method, and eliminating the fast states with a stated
justification. This module does both, because the residual spread after reduction is still
~1e2-1e3 and any sweep will push a state faster than expected.

THE REDUCTION CRITERION, AND WHAT IT IS NOT. A quasi-steady-state elimination replaces
``dx/dt = (u - x)/tau`` with ``x = u``. Over a window T the bias it puts on the TIME-INTEGRAL of
x -- which is exactly what a stable-FP plate read records -- is ``(tau/T)(1 - exp(-T/tau))``,
bounded above by ``tau/T``. This repository has measured its own floor for that observable:
``panel_experiment.OBSERVED_ACTIVITY_CV = 0.146``, across 28 matched conditions on two real
plates. Hence the rule, and it is imported rather than restated here:

    reduce a state to algebra iff tau/T < 0.146   (equivalently T/tau > 6.85)

The rule is symmetric and the audit applies both halves. A state far SLOWER than the window
cannot move enough to be seen either, so ``T/tau < 0.146`` licenses freezing it at its initial
condition, with the same bias bound the other way up. The irreducible band is therefore
``0.146 <= tau/T <= 6.85``, and everything outside it reduces -- in one of two different ways.

THREE THINGS THE RATIO TEST DOES NOT SETTLE, all of them recorded on every row:

1. **It is necessary and not sufficient.** A 2-minute Crz1 burst has tau/T = 0.008 on a plate
   read and passes easily. Reducing it is still wrong: Cai, Dalal & Elowitz 2008 (Nature
   455:485-90, PMID 18818649) MEASURED that calcium sets the burst FREQUENCY and not its
   duration, so the pulse train is a self-generated carrier and not a fast state tracking a slow
   input. Replacing it by its mean is a mean-field approximation, and Cai 2008 is a direct
   demonstration that the mean does not determine the output. So a state marked
   ``Encoding.FREQUENCY`` gets :attr:`Verdict.FREQUENCY_ENCODED` whatever its ratio, and the
   only route to reducing it is :attr:`Justification.MEAN_FIELD` -- which a caller must DECLARE
   with a reason and which this module never infers.

2. **tau/T is a systematic bias and 0.146 is a random CV.** A dose-correlated 10% bias in every
   well does not average out the way a 15% well-to-well spread does, so a row that passes by a
   small margin is passing against the wrong kind of yardstick. Every row therefore reports the
   exact bias alongside the ratio, so a marginal pass is visible rather than implied.

3. **The bias formula assumes a persistent driver.** It is derived for a step input. When a
   state's driver vanishes inside the window -- Pincus 2010 (PMID 20625545) MEASURED HAC1
   splicing deactivating within 2 h at 1.5 mM DTT, well inside a 4.14 h read -- the window that
   matters is the driver's, not the run's. :func:`audit_timescales` accepts
   ``driver_persists_h`` per state and scores against ``T_eff = min(T, driver)``. It is never
   assumed: with nothing declared, T_eff is the whole window.

THE WINDOW IS A DECLARED INPUT AND THE TWO REDUCTION SETS ARE NOT NESTED. Run the audit on
:data:`~ystwin.mech.state.TIMESCALE_CATALOGUE` and the plate read keeps the transcript pool, FP
maturation, the UPR and cross-protection while freezing plasmid loss; the 5-day fed-batch
eliminates every one of those and keeps plasmid loss. Neither set contains the other, so a model
tuned to one is mis-specified for the other. That is arithmetic on measured time constants, not
a modelling preference, which is why :class:`~ystwin.mech.state.Window` is an argument here
rather than a default.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.integrate import solve_ivp

from ..generator.panel_experiment import MEASURED_GROWTH_RATE_SE, OBSERVED_ACTIVITY_CV
from .state import Encoding, MechState, StateVar, Window

__all__ = [
    "ATOL",
    "AuditRow",
    "FALLBACK_METHOD",
    "GROWTH_FLOOR_PER_H",
    "Justification",
    "PINNED_METHOD",
    "REDUCTION_RATIO_CEILING",
    "RETAINED_RATIO_CEILING",
    "RTOL",
    "ReducedSystem",
    "ReductionRefused",
    "Trajectory",
    "Verdict",
    "audit_timescales",
    "freeze_bias",
    "integrate_window",
    "qss_bias",
    "reduce_for_window",
    "require_reducible",
]

PINNED_METHOD = "BDF"
"""The implicit method this package integrates with, pinned exactly as `fba/solver.py` pins the
LP solver and for the same reason: an explicit method on this system does not give a different
answer slowly, it gives a step size chosen by chemistry the observable cannot see."""

FALLBACK_METHOD = "LSODA"
"""Used only when BDF reports failure. LSODA switches between a non-stiff and a stiff method on
its own, so it recovers from a right-hand side whose stiffness turns on part way through a run."""

RTOL = 1e-6
"""Relative tolerance. ASSERTED as a numerical setting, not fitted: it is not a model parameter
and has no target to be counted against. `tests/test_mech_integrate.py` pins it by showing that
tightening it tenfold moves the reported observable by far less than the assay floor."""

ATOL = 1e-9
"""Default absolute tolerance. ASSERTED, with the same standing as RTOL. Pass a per-state array
or a name-keyed mapping to :func:`integrate_window` when the states differ in magnitude by more
than a few orders -- a single scalar is a claim that they do not."""

REDUCTION_RATIO_CEILING = OBSERVED_ACTIVITY_CV
"""tau/T below which a quasi-steady-state elimination is licensed: 0.146.

IMPORTED, never restated. It is `generator/panel_experiment.OBSERVED_ACTIVITY_CV`, the
plate-to-plate CV of dilution-corrected activity measured across 28 matched conditions on the
2026-07-22 and 2026-08-03 plates. Restating the digits here would let the two drift apart, and a
reduction criterion that no longer matches the instrument it was derived from is worse than none.
"""

RETAINED_RATIO_CEILING = 1.0 / REDUCTION_RATIO_CEILING
"""tau/T above which the state is too slow to move inside the window: 6.849. The other half of
the same inequality, so it carries no independent number."""

GROWTH_FLOOR_PER_H = MEASURED_GROWTH_RATE_SE
"""0.0117 /h, the standard error of an estimated specific growth rate across 210 real wells.
The floor for any ablation scored on GROWTH rather than on reporter activity. It is here so a
caller scoring a growth observable does not reach for the reporter CV, which is the substitution
this project's own review calls out by name."""


class ReductionRefused(RuntimeError):
    """A caller asked to reduce a state the audit does not license reducing.

    Deliberately not a ``ValueError``: the request is well-formed and the refusal is a finding
    about the model, not about the argument. `pathway/capacity.py::CapacityUnmeasured` is the
    precedent.
    """


class Verdict(enum.Enum):
    """What the audit says about one state in one window."""

    ELIMINATE = "eliminate"
    """tau/T < 0.146 over the whole interval, and the state relaxes toward a value its driver
    sets. Replace the differential equation by its algebraic quasi-steady solution."""

    FREEZE = "freeze"
    """T/tau < 0.146 over the whole interval. The state cannot move far enough inside the window
    to be seen; hold it at its initial condition. This is the OTHER reduction, and it is what
    makes the two windows' irreducible sets non-nested."""

    KEEP = "keep"
    """Inside the irreducible band 0.146 <= tau/T <= 6.849. It needs a differential equation."""

    STRADDLES = "straddles"
    """The measured time constant is an INTERVAL that crosses a threshold, so the sources
    disagree about whether the state is reducible. Refused rather than resolved: picking the
    convenient end of a 9x methodological disagreement is how a model acquires a bias nobody
    can find later. Report which measurement would settle it."""

    FREQUENCY_ENCODED = "frequency-encoded"
    """The ratio test passes and reduction is still wrong, because the state is a self-generated
    carrier whose frequency is the signal (Cai 2008, PMID 18818649). Reducible only under a
    declared :attr:`Justification.MEAN_FIELD`."""

    REFUSED = "refused"
    """No source measures this state's time constant, so it cannot be audited at all."""


class Justification(enum.Enum):
    """The ground a reduction rests on. Every reduction records exactly one."""

    RATIO = "ratio"
    """The measured time constant against the measured assay floor, for a state that relaxes
    toward a value its driver sets. The only justification this module grants on its own."""

    MEAN_FIELD = "mean-field"
    """Replacing a pulse train by its mean. NOT a quasi-steady-state reduction and not licensed
    by the ratio: Cai 2008 (PMID 18818649) measured that the mean does not determine the output,
    and Hansen & O'Shea 2013 measured that two promoters driven by the same train respond to
    different functions of it. A caller must declare it with a reason; it is never inferred."""


def qss_bias(tau_h: float, window_h: float) -> float:
    """Fractional bias a quasi-steady elimination puts on the time-integral of the state.

    Exact for a step input on ``dx/dt = (u - x)/tau`` from ``x(0) = 0``: the reduced answer is
    ``u*T`` and the true one is ``u*(T - tau(1 - exp(-T/tau)))``, so the reduced integral is high
    by ``(tau/T)(1 - exp(-T/tau))`` of itself. Bounded above by ``tau/T``, which is the ratio the
    criterion actually tests -- so the criterion is conservative, and by how much is visible.
    """
    if tau_h <= 0 or window_h <= 0:
        raise ValueError(f"tau and window must be positive, got tau={tau_h}, T={window_h}")
    return (tau_h / window_h) * (-math.expm1(-window_h / tau_h))


def freeze_bias(tau_h: float, window_h: float) -> float:
    """Fractional bias from holding a slow state at its initial condition instead of decaying it.

    The complement of :func:`qss_bias` on the same integral: the time-average of
    ``exp(-t/tau)`` over ``[0, T]`` is ``(tau/T)(1 - exp(-T/tau))``, so freezing at 1 is high by
    one minus that. To first order it is ``T/(2*tau)``, half the ``T/tau`` the criterion tests --
    conservative in the same direction as the fast half.
    """
    return 1.0 - qss_bias(tau_h, window_h)


@dataclass(frozen=True)
class AuditRow:
    """One state's tau/T ratio, verdict and justification, in one window.

    Args:
        state: The variable audited.
        window: The declared window it was audited against.
        effective_window_h: ``min(window.duration_h, driver_persists_h)``. Differs from the
            window only when a source MEASURES how long the state's driver is present.
        tau_low_h, tau_high_h: The time-constant interval in this window's growth band.
        ratio_low, ratio_high: ``tau/T_eff`` at each end of that interval.
        bias_low, bias_high: The exact time-integral bias the reduction would introduce, from
            :func:`qss_bias` for a fast state and :func:`freeze_bias` for a slow one.
        verdict: See :class:`Verdict`.
        justification: Set only when the verdict is a reduction.
        detail: Why, in one sentence, when the verdict is not a plain reduction.
    """

    state: StateVar
    window: Window
    effective_window_h: float
    tau_low_h: float | None
    tau_high_h: float | None
    ratio_low: float | None
    ratio_high: float | None
    bias_low: float | None
    bias_high: float | None
    verdict: Verdict
    justification: Justification | None = None
    detail: str = ""

    @property
    def reduces(self) -> bool:
        return self.verdict in (Verdict.ELIMINATE, Verdict.FREEZE)

    def line(self) -> str:
        """One fixed-width row of :meth:`ReducedSystem.audit_table`."""
        if self.ratio_low is None:
            ratio = "      REFUSED      "
            bias = "      --      "
        else:
            ratio = f"{self.ratio_low:>8.3g} - {self.ratio_high:<8.3g}"
            bias = f"{self.bias_low:>5.1%} -{self.bias_high:>6.1%}"
        just = self.justification.value if self.justification else ""
        return (f"  {self.state.name:<20s} {self.state.encoding.value:<12s} {ratio}  "
                f"{bias}  {self.verdict.value:<17s} {just}")


@dataclass(frozen=True)
class ReducedSystem:
    """The reduced system for one window, plus the audit table that licensed every cut.

    Both halves are returned together on purpose. A reduction set with no audit beside it is an
    assertion about which states matter, and this repository has been burned by exactly that
    shape of claim before.

    Args:
        window: The declared window.
        integrated: States that still need a differential equation, in catalogue order. This is
            the system :func:`integrate_window` should be handed.
        eliminated: States replaced by their algebraic quasi-steady solution.
        frozen: States held at their initial condition.
        audit: One :class:`AuditRow` per state, including the ones that were refused.
    """

    window: Window
    integrated: MechState
    eliminated: MechState
    frozen: MechState
    audit: tuple[AuditRow, ...]

    def row(self, name: str) -> AuditRow:
        for r in self.audit:
            if r.state.name == name:
                return r
        raise KeyError(f"no state named {name!r} in this audit")

    def verdict(self, name: str) -> Verdict:
        return self.row(name).verdict

    def audit_table(self) -> str:
        """The audit table, with the caveats that make the ratio readable, as printable text."""
        head = [
            f"tau/T audit against {self.window.name}: T = {self.window.duration_h:g} h, "
            f"mu = {self.window.growth_rate_low_per_h:g}-{self.window.growth_rate_high_per_h:g} /h",
            f"  criterion: eliminate iff tau/T < {REDUCTION_RATIO_CEILING:g}, freeze iff "
            f"T/tau < {REDUCTION_RATIO_CEILING:g} (panel_experiment.OBSERVED_ACTIVITY_CV, "
            f"28 matched conditions on two real plates)",
            "  THE RATIO TEST IS NECESSARY AND NOT SUFFICIENT. A 2 min Crz1 burst passes it and",
            "  reducing it is still wrong: calcium sets the FREQUENCY, not the duration",
            "  (Cai 2008, PMID 18818649). And tau/T is a SYSTEMATIC bias while 0.146 is a RANDOM",
            "  CV, so a marginal pass is not the same kind of quantity -- the bias column is the",
            "  exact figure the ratio bounds.",
            "",
            f"  {'state':<20s} {'encoding':<12s} {'tau/T':^19s}  {'bias':^14s}  "
            f"{'verdict':<17s} justification",
        ]
        body = [r.line() for r in self.audit]
        held = {v: sum(1 for r in self.audit if r.verdict is v)
                for v in (Verdict.KEEP, Verdict.STRADDLES,
                          Verdict.FREQUENCY_ENCODED, Verdict.REFUSED)}
        tail = [
            "",
            f"  eliminated {len(self.eliminated)}, frozen {len(self.frozen)}, "
            f"integrated {len(self.integrated)} -- of which {held[Verdict.KEEP]} in the "
            f"irreducible band, {held[Verdict.STRADDLES]} straddling, "
            f"{held[Verdict.FREQUENCY_ENCODED]} frequency-encoded, "
            f"{held[Verdict.REFUSED]} with no measured time constant",
        ]
        return "\n".join(head + body + tail)


def _effective_window(window: Window, name: str,
                      driver_persists_h: Mapping[str, float] | None) -> float:
    if not driver_persists_h or name not in driver_persists_h:
        return window.duration_h
    declared = float(driver_persists_h[name])
    if declared <= 0:
        raise ValueError(f"driver_persists_h[{name!r}] must be positive, got {declared}")
    return min(window.duration_h, declared)


def _audit_one(state: StateVar, window: Window,
               driver_persists_h: Mapping[str, float] | None) -> AuditRow:
    t_eff = _effective_window(window, state.name, driver_persists_h)
    if state.timescale_is_refused:
        return AuditRow(
            state=state, window=window, effective_window_h=t_eff,
            tau_low_h=None, tau_high_h=None, ratio_low=None, ratio_high=None,
            bias_low=None, bias_high=None, verdict=Verdict.REFUSED,
            detail="no measured time constant, so the ratio cannot be formed",
        )

    tau_low, tau_high = state.taus_in(window)
    r_low, r_high = tau_low / t_eff, tau_high / t_eff
    fast_all, fast_any = r_high < REDUCTION_RATIO_CEILING, r_low < REDUCTION_RATIO_CEILING
    slow_all, slow_any = r_low > RETAINED_RATIO_CEILING, r_high > RETAINED_RATIO_CEILING

    def audited(verdict, justification=None, detail="", biases=None):
        b_low, b_high = biases if biases else (None, None)
        return AuditRow(state=state, window=window, effective_window_h=t_eff,
                        tau_low_h=tau_low, tau_high_h=tau_high,
                        ratio_low=r_low, ratio_high=r_high,
                        bias_low=b_low, bias_high=b_high,
                        verdict=verdict, justification=justification, detail=detail)

    qss = (qss_bias(tau_low, t_eff), qss_bias(tau_high, t_eff))
    frozen = (freeze_bias(tau_high, t_eff), freeze_bias(tau_low, t_eff))

    if state.encoding is Encoding.FREQUENCY:
        return audited(Verdict.FREQUENCY_ENCODED, detail=(
            "the frequency is the signal, so the ratio does not license a reduction; declare "
            "Justification.MEAN_FIELD with a reason if the mean is to be used anyway"),
            biases=qss)

    if state.encoding is Encoding.INTEGRATING:
        if slow_all:
            return audited(Verdict.FREEZE, Justification.RATIO, biases=frozen)
        if slow_any:
            return audited(Verdict.STRADDLES, detail=(
                "the interval crosses T/tau = 0.146; the sources disagree about whether it moves "
                "inside the window"), biases=frozen)
        return audited(Verdict.KEEP, detail=(
            "accumulates, so it has no quasi-steady value to be eliminated to; its 1/mu is a "
            "doubling timescale, not a relaxation time"), biases=qss)

    if fast_all:
        return audited(Verdict.ELIMINATE, Justification.RATIO, biases=qss)
    if slow_all:
        return audited(Verdict.FREEZE, Justification.RATIO, biases=frozen)
    if fast_any:
        return audited(Verdict.STRADDLES, detail=(
            f"tau/T runs {r_low:.3g} to {r_high:.3g} across the sources, crossing "
            f"{REDUCTION_RATIO_CEILING:g}; the reduction is licensed by one measurement and "
            f"refused by the other"), biases=qss)
    if slow_any:
        return audited(Verdict.STRADDLES, detail=(
            f"tau/T runs {r_low:.3g} to {r_high:.3g}, crossing {RETAINED_RATIO_CEILING:.3g}; "
            f"the freeze is licensed by one measurement and refused by the other"), biases=frozen)
    return audited(Verdict.KEEP, biases=qss)


def audit_timescales(states: MechState, window: Window, *,
                     driver_persists_h: Mapping[str, float] | None = None) -> tuple[AuditRow, ...]:
    """Every state's tau/T ratio, exact bias, verdict and justification against one window.

    Args:
        states: The system to audit.
        window: The declared operating window. Not optional and not defaulted -- the irreducible
            set differs between this project's two vessels and the sets are not nested.
        driver_persists_h: Per-state hours over which the state's DRIVER is present, where a
            source measures it. Scoring uses ``min(window, driver)``, because the tau/T bias
            formula is derived for a step and a driver that vanishes inside the window makes the
            run length the wrong denominator. Omit it and the whole window is used.

    Returns:
        One :class:`AuditRow` per state, in the system's order.
    """
    return tuple(_audit_one(v, window, driver_persists_h) for v in states)


def reduce_for_window(states: MechState, window: Window, *,
                      driver_persists_h: Mapping[str, float] | None = None) -> ReducedSystem:
    """The reduced system for a window, and the audit table that licensed every cut.

    A state is eliminated to algebra, frozen at its initial condition, or integrated. Anything
    the audit refuses -- a straddling interval, a frequency-encoded carrier, an unmeasured time
    constant -- is INTEGRATED, because the safe direction when the criterion cannot decide is to
    keep the differential equation and pay for it.

    Args:
        states: The system to reduce.
        window: The declared operating window.
        driver_persists_h: See :func:`audit_timescales`.
    """
    rows = audit_timescales(states, window, driver_persists_h=driver_persists_h)
    eliminated = tuple(r.state for r in rows if r.verdict is Verdict.ELIMINATE)
    frozen = tuple(r.state for r in rows if r.verdict is Verdict.FREEZE)
    reduced = {v.name for v in eliminated} | {v.name for v in frozen}
    integrated = tuple(v for v in states if v.name not in reduced)
    return ReducedSystem(
        window=window,
        integrated=MechState(integrated),
        eliminated=MechState(eliminated),
        frozen=MechState(frozen),
        audit=rows,
    )


def require_reducible(state: StateVar, window: Window,
                      justification: Justification = Justification.RATIO, *,
                      reason: str = "",
                      driver_persists_h: Mapping[str, float] | None = None) -> AuditRow:
    """Return the audit row, or refuse -- the gate a caller goes through to reduce a state.

    :func:`reduce_for_window` decides for a whole system; this is for the caller who has already
    decided and wants the decision checked. It exists so that a reduction taken by hand somewhere
    in ``mech/`` fails loudly rather than silently disagreeing with the audit table.

    Args:
        state: The variable the caller wants to reduce.
        window: The declared operating window.
        justification: :attr:`Justification.RATIO` is checked against the audit.
            :attr:`Justification.MEAN_FIELD` is the only route to reducing a frequency-encoded
            state, requires ``reason``, and is recorded rather than verified -- this module
            cannot check a claim about a carrier it does not simulate.
        reason: Required with ``MEAN_FIELD``. What makes the mean adequate here, given that
            Cai 2008 measured that it is not adequate in general.
        driver_persists_h: See :func:`audit_timescales`.

    Raises:
        ReductionRefused: with the ratio, the floor and what would settle it.
    """
    row = _audit_one(state, window, driver_persists_h)

    if justification is Justification.MEAN_FIELD:
        if state.encoding is not Encoding.FREQUENCY:
            raise ReductionRefused(
                f"{state.name!r} is {state.encoding.value}-encoded, so a mean-field argument is "
                f"not what licenses reducing it; use Justification.RATIO and let the audit decide")
        if not reason.strip():
            raise ReductionRefused(
                f"{state.name!r}: Justification.MEAN_FIELD requires a reason. Cai 2008 "
                f"(PMID 18818649) measured that the mean of a frequency-modulated train does not "
                f"determine the output, so 'the ratio passes' is not one")
        if row.ratio_high is not None and row.ratio_high >= REDUCTION_RATIO_CEILING:
            raise ReductionRefused(
                f"{state.name!r}: tau/T = {row.ratio_low:.3g}-{row.ratio_high:.3g} does not pass "
                f"{REDUCTION_RATIO_CEILING:g}, so a mean-field reduction is not even fast enough "
                f"to be worth arguing about")
        return row

    if row.verdict is Verdict.REFUSED:
        refusal = state.tau_h
        raise ReductionRefused(
            f"{state.name!r}: {row.detail}. {refusal.reason} "
            f"(nearest source: {refusal.source})")
    if row.verdict is Verdict.FREQUENCY_ENCODED:
        raise ReductionRefused(
            f"{state.name!r} passes the ratio (tau/T = {row.ratio_low:.3g}) and is "
            f"frequency-encoded, so the ratio alone does not license the reduction. {row.detail}. "
            f"Source: {state.source}")
    if row.verdict is Verdict.STRADDLES:
        raise ReductionRefused(
            f"{state.name!r}: {row.detail}. Source: {state.source}")
    if row.verdict is Verdict.KEEP:
        raise ReductionRefused(
            f"{state.name!r}: tau/T = {row.ratio_low:.3g}-{row.ratio_high:.3g} against "
            f"{window.name} (T = {row.effective_window_h:g} h) sits inside the irreducible band "
            f"[{REDUCTION_RATIO_CEILING:g}, {RETAINED_RATIO_CEILING:.4g}]. Reducing it would bias "
            f"the time-integral by {row.bias_low:.1%}-{row.bias_high:.1%} against a "
            f"{REDUCTION_RATIO_CEILING:.1%} assay floor. {row.detail}")
    return row


@dataclass(frozen=True)
class Trajectory:
    """The result of one integration: the states over time, and what it cost.

    Args:
        t: Time grid, hours.
        y: ``(n_states, n_times)`` solution in the system's order.
        system: The states, so a trace can be pulled by name and labelled with its units.
        method: Which integrator actually ran -- BDF, or LSODA if BDF reported failure.
        n_rhs_evals: Right-hand-side evaluations. The honest cost measure, and what the
            stiff-versus-explicit comparison in the tests is scored on.
    """

    t: np.ndarray
    y: np.ndarray
    system: MechState
    method: str
    n_rhs_evals: int

    def of(self, name: str) -> np.ndarray:
        """One state's trace, by name."""
        return self.y[self.system.index(name)]

    def time_integral(self, name: str) -> float:
        """Trapezoidal integral of one state over the run, in state-units x hours.

        This is the observable, not a diagnostic: a stable fluorescent protein read at the end of
        a plate run records the accumulated integral of what produced it, which is why the
        reduction criterion is stated as a bias on this quantity and not on the endpoint.
        """
        return float(np.trapezoid(self.of(name), self.t))

    def time_average(self, name: str) -> float:
        """The time-integral divided by the run length."""
        return self.time_integral(name) / (self.t[-1] - self.t[0])


def _atol_array(atol: float | Sequence[float] | Mapping[str, float],
                system: MechState) -> float | np.ndarray:
    if isinstance(atol, Mapping):
        missing = sorted(set(system.names) - set(atol))
        if missing:
            raise KeyError(f"no absolute tolerance for {missing}")
        return np.array([float(atol[n]) for n in system.names], dtype=float)
    if isinstance(atol, (int, float)):
        return float(atol)
    arr = np.asarray(atol, dtype=float)
    if arr.shape != (len(system),):
        raise ValueError(f"atol has shape {arr.shape}, expected ({len(system)},)")
    return arr


def integrate_window(rhs: Callable[[float, np.ndarray], Sequence[float]],
                     y0: Mapping[str, float] | Sequence[float],
                     window: Window,
                     system: MechState, *,
                     n_points: int = 201,
                     method: str = PINNED_METHOD,
                     rtol: float = RTOL,
                     atol: float | Sequence[float] | Mapping[str, float] = ATOL,
                     nonnegative: bool = True,
                     t0_h: float = 0.0) -> Trajectory:
    """Integrate ``system`` across ``window`` with the pinned stiff method.

    Args:
        rhs: ``f(t, y) -> dy/dt`` over the system's states in order, t in hours.
        y0: Initial values, name-keyed or in the system's order.
        window: The declared window. Its ``duration_h`` sets the end of the run.
        system: The states being integrated, which supplies the order and the names.
        n_points: Output grid size. Output resolution only; the solver picks its own steps.
        method: Integrator. Defaults to :data:`PINNED_METHOD`. Passing an explicit method is
            allowed and is how the tests measure what an explicit one costs on this system.
        rtol, atol: Tolerances. ``atol`` may be a scalar, a per-state sequence, or a name-keyed
            mapping.
        nonnegative: Clip the state to zero INSIDE the right-hand side rather than clipping the
            solution afterwards. Clipping the state would make the trajectory discontinuous and
            defeat the implicit method's Jacobian; clipping the argument keeps the flow smooth
            and is what `kinetic/carotenoid.py` and `estimator._propagate` already do.
        t0_h: Start of the run. Defaults to 0.

    Returns:
        A :class:`Trajectory`, carrying which method actually ran.

    Raises:
        RuntimeError: if both the requested method and the fallback fail.
    """
    y_start = system.pack(y0) if isinstance(y0, Mapping) else np.asarray(y0, dtype=float)
    if y_start.shape != (len(system),):
        raise ValueError(f"initial state has shape {y_start.shape}, expected ({len(system)},)")

    def clipped(t, y):
        return rhs(t, np.maximum(y, 0.0)) if nonnegative else rhs(t, y)

    t_end = t0_h + window.duration_h
    grid = np.linspace(t0_h, t_end, int(n_points))
    tolerance = _atol_array(atol, system)

    attempts = [method] if method != PINNED_METHOD else [PINNED_METHOD, FALLBACK_METHOD]
    last = None
    for attempt in attempts:
        sol = solve_ivp(clipped, (t0_h, t_end), y_start, t_eval=grid, method=attempt,
                        rtol=rtol, atol=tolerance, dense_output=True)
        if sol.success:
            return Trajectory(t=sol.t, y=sol.y, system=system, method=attempt,
                              n_rhs_evals=int(sol.nfev))
        last = sol
    raise RuntimeError(
        f"integration failed under {attempts} across {window.name}: {last.message}")
