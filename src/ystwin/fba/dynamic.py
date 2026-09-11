"""Dynamic FBA: a batch or fed-batch run integrated over time, so titre means something.

Everything else in this package solves a STEADY STATE. `pathway/solve.py` closes
``d[X]/dt = 0`` and returns a content; `predict.py` reports that content and a rate. That is
the right arithmetic for a chemostat, where the pump fixes ``mu`` and nothing accumulates.

It is the wrong arithmetic for the vessel anyone actually manufactures in. In a batch or
fed-batch culture nothing is at steady state: substrate depletes, biomass accumulates,
``mu`` falls as the carbon runs out, and TITRE IS AN INTEGRAL rather than a concentration
held constant by a pump. Those are the quantities a process is judged on, and this module
computes them the only way they can be computed -- by integrating.

    dX/dt  =  mu(t) * X                       biomass, gDCW/L
    dS/dt  = -q_S(t) * X  +  feed(t)          substrate, mmol/L
    dP/dt  =  q_P(t) * X                      product, mmol/L   <- the titre

At each step the GEM is solved for ``mu`` and the exchange rates under the substrate
available at that instant, and the product rate comes from whichever layer the caller says
sets it. So this module supplies the DYNAMICS and never the biology: it integrates rates it
is given, and refuses when it is not given them.

WHY THIS IS NOT A REPLACEMENT FOR THE STEADY-STATE CHAIN. The two answer different
questions and neither subsumes the other. The chemostat chain answers "what content at this
growth rate", which is a mechanistic question with a calibrated answer. This answers "what
titre after 72 hours", which is a process question. What makes the second worth building is
that the DATA is different too: batch endpoint titres are far more numerous in the
literature than chemostat steady states, and they span strains and laboratories that no
chemostat series covers.

WHAT IT CANNOT DO, stated before anyone reads a number out of it. It has never been scored.
Every dataset vendored in this repository is a chemostat series except
`data/carotenoid_batch/`, which is endpoint titres with no time course, so this module can be
compared against a FINAL number and not against a trajectory. A trajectory that arrives at
the right endpoint by the wrong route is not detectable here, and saying so is cheaper than
discovering it later.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ..mech.adapters import simulate_protocol

__all__ = [
    "BatchResult",
    "FeedProfile",
    "simulate_batch",
    "simulate_fixed_volume_comparison",
    "simulate_protocol",
]


@dataclass(frozen=True)
class FeedProfile:
    """Substrate fed into the vessel, mmol/L/h, as a function of time.

    ``rate = 0.0`` is a BATCH culture: nothing is added, and the run ends when the substrate
    is gone. A positive rate is fed-batch. The distinction is one number rather than two
    code paths, because it is one number in the physics.
    """

    rate_mmol_per_l_h: float = 0.0
    start_h: float = 0.0
    end_h: float = float("inf")

    def __post_init__(self) -> None:
        if not np.isfinite(self.rate_mmol_per_l_h) or self.rate_mmol_per_l_h < 0:
            raise ValueError("feed rate_mmol_per_l_h must be finite and non-negative")
        if not np.isfinite(self.start_h) or self.start_h < 0:
            raise ValueError("feed start_h must be finite and non-negative")
        if np.isnan(self.end_h) or self.end_h < self.start_h:
            raise ValueError("feed end_h must not precede start_h")

    def integral(self, start_h: float, end_h: float) -> float:
        if self.rate_mmol_per_l_h == 0.0:
            return 0.0
        amount = float(self.rate_mmol_per_l_h) * max(
            0.0, min(end_h, self.end_h) - max(start_h, self.start_h))
        if not np.isfinite(amount):
            raise ValueError("integrated feed must be finite")
        return amount

    def at(self, hour: float) -> float:
        return (self.rate_mmol_per_l_h
                if self.start_h <= hour < self.end_h else 0.0)


@dataclass(frozen=True)
class BatchResult:
    """One simulated run, with the trajectory kept so the endpoint can be checked."""

    hours: np.ndarray
    biomass_g_per_l: np.ndarray
    substrate_mmol_per_l: np.ndarray
    product_mmol_per_l: np.ndarray
    growth_rate_per_h: np.ndarray
    product_molar_mass_g_per_mol: float
    substrate_molar_mass_g_per_mol: float | None = None
    total_feed_mmol_per_l: float = 0.0
    _integrated_uptake_mmol_per_l: float | None = field(default=None, repr=False)
    _biomass_exposure_g_h_per_l: float | None = field(default=None, repr=False)

    @property
    def simulation_route(self) -> str:
        return "fixed_volume_comparison"

    @property
    def substrate_consumed_mmol_per_l(self) -> float:
        if self._integrated_uptake_mmol_per_l is not None:
            return self._integrated_uptake_mmol_per_l
        return float(self.substrate_mmol_per_l[0] + self.total_feed_mmol_per_l
                     - self.substrate_mmol_per_l[-1])

    @property
    def average_productivity_mg_per_l_h(self) -> float:
        return self.titre_mg_per_l / float(self.hours[-1] - self.hours[0])

    @property
    def biomass_exposure_g_h_per_l(self) -> float:
        if self._biomass_exposure_g_h_per_l is not None:
            return float(self._biomass_exposure_g_h_per_l)
        return float(np.sum(self.biomass_g_per_l[:-1] * np.diff(self.hours)))

    @property
    def mean_specific_rate_mmol_per_gdcw_h(self) -> float:
        exposure = self.biomass_exposure_g_h_per_l
        return float(self.product_mmol_per_l[-1] / exposure) if exposure > 0 else float("nan")

    @property
    def yield_g_per_g(self) -> float:
        if self.substrate_molar_mass_g_per_mol is None:
            raise ValueError("substrate_molar_mass_g_per_mol is required to report a mass yield")
        consumed = self.substrate_consumed_mmol_per_l * self.substrate_molar_mass_g_per_mol
        return self.titre_mg_per_l / consumed if consumed > 0 else float("nan")

    @property
    def titre_mg_per_l(self) -> float:
        """The number a process is judged on."""
        return float(self.product_mmol_per_l[-1] * self.product_molar_mass_g_per_mol)

    @property
    def final_biomass_g_per_l(self) -> float:
        return float(self.biomass_g_per_l[-1])

    @property
    def content_mg_per_gdcw(self) -> float:
        """Titre divided by biomass -- the unit the carotenoid literature reports.

        Not the same quantity as the steady-state chain's content, and the difference is
        not cosmetic: this is an accumulated total over a run whose growth rate changed,
        while the chain's is a balance at one fixed growth rate.
        """
        return self.titre_mg_per_l / self.final_biomass_g_per_l

    @property
    def substrate_exhausted(self) -> bool:
        return bool(self.substrate_mmol_per_l[-1] <= 1e-9)

    def summary(self) -> str:
        return (f"{self.hours[-1]:.0f} h: {self.titre_mg_per_l:.1f} mg/L, "
                f"{self.final_biomass_g_per_l:.2f} gDCW/L, "
                f"{self.content_mg_per_gdcw:.2f} mg/gDCW"
                f"{'; substrate exhausted' if self.substrate_exhausted else ''}")


def simulate_fixed_volume_comparison(
    growth_and_rates,
    initial_biomass_g_per_l: float,
    initial_substrate_mmol_per_l: float,
    hours: float,
    product_molar_mass_g_per_mol: float,
    feed: FeedProfile | None = None,
    steps: int = 400,
    substrate_molar_mass_g_per_mol: float | None = None,
    *,
    rtol: float = 1e-7,
    atol: float = 1e-10,
) -> BatchResult:
    """Integrate a batch or fed-batch run.

    Args:
        growth_and_rates: ``f(substrate_mmol_per_l) -> (mu_per_h, q_substrate, q_product)``,
            with ``q_substrate`` positive for uptake and ``q_product`` positive for
            formation, both mmol/gDCW/h. THIS is where the GEM and the pathway layers enter;
            this module never decides biology. A caller that wants FBA to supply ``mu`` and
            ``q_substrate`` passes a closure that solves the model.
        initial_biomass_g_per_l: Inoculum.
        initial_substrate_mmol_per_l: Starting carbon.
        hours: Run length.
        product_molar_mass_g_per_mol: To report a titre in mg/L.
        feed: Fed-batch profile, or ``None`` for batch.
        steps: Output intervals and an upper bound on the adaptive integration step size.
            Local error is controlled by ``rtol`` and ``atol`` independently of this grid;
            feed boundaries and substrate depletion are resolved within each interval.

    Raises:
        ValueError: on a non-positive duration, biomass or step count.
    """
    if not np.isfinite(hours) or hours <= 0:
        raise ValueError(f"hours must be positive and finite, got {hours}")
    if not np.isfinite(initial_biomass_g_per_l) or initial_biomass_g_per_l <= 0:
        raise ValueError(
            f"initial_biomass_g_per_l must be positive, got {initial_biomass_g_per_l}. A "
            f"culture inoculated with nothing stays at nothing -- dX/dt = mu*X is zero for "
            f"all time at X = 0, so this is a silent no-growth run rather than an error")
    if not np.isfinite(initial_substrate_mmol_per_l) or initial_substrate_mmol_per_l < 0:
        raise ValueError("initial_substrate_mmol_per_l must be finite and non-negative")
    for name, value in (("product_molar_mass_g_per_mol", product_molar_mass_g_per_mol),
                        ("substrate_molar_mass_g_per_mol", substrate_molar_mass_g_per_mol)):
        if value is None and name == "substrate_molar_mass_g_per_mol":
            continue
        if value is None or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if isinstance(steps, bool) or not isinstance(steps, (int, np.integer)) or steps < 2:
        raise ValueError(f"steps must be an integer of at least 2, got {steps}")
    for name, value in (("rtol", rtol), ("atol", atol)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if rtol < 100 * np.finfo(float).eps:
        raise ValueError("rtol must be at least 100 times floating-point precision")
    if feed is None:
        feed = FeedProfile()
    if not isinstance(feed, FeedProfile):
        raise ValueError("feed must be a FeedProfile or None")
    hours = float(hours)
    initial_biomass_g_per_l = float(initial_biomass_g_per_l)
    initial_substrate_mmol_per_l = float(initial_substrate_mmol_per_l)
    total_feed = feed.integral(0.0, hours)
    supply = initial_substrate_mmol_per_l + total_feed
    if not np.isfinite(supply):
        raise ValueError("initial substrate plus integrated feed must be finite")
    if not np.isfinite(initial_biomass_g_per_l * hours):
        raise ValueError("integrated biomass exposure must be finite")
    if not callable(growth_and_rates):
        raise ValueError("growth_and_rates must be callable")

    def checked_rates(substrate, hour):
        observed = growth_and_rates(float(substrate))
        try:
            rates = np.asarray(observed, dtype=float)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"growth and substrate/product rates must be three finite numbers at {hour:g} h") from exc
        if rates.shape != (3,) or not np.all(np.isfinite(rates)) or np.any(rates < 0):
            raise ValueError(f"growth and substrate/product rates must be finite and non-negative, "
                             f"got {rates} at {hour:g} h")
        mu, q_s, q_p = rates
        if q_s == 0.0 and (mu > 0.0 or q_p > 0.0):
            raise ValueError(f"growth and product rates require substrate uptake at {hour:g} h")
        if substrate == 0.0 and np.any(rates != 0.0):
            raise ValueError(f"rates must be zero at depleted substrate at {hour:g} h")
        return mu, q_s, q_p

    t = np.linspace(0.0, hours, steps + 1)
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("hours and steps must define distinct output times")
    trajectory = np.zeros((5, steps + 1))
    state = np.array([initial_biomass_g_per_l, initial_substrate_mmol_per_l,
                      0.0, 0.0, 0.0])
    trajectory[:, 0] = state
    checked_rates(state[1], 0.0)
    feed_changes = (feed.start_h, feed.end_h) if total_feed > 0.0 else ()
    boundaries = sorted({0.0, hours, *(float(boundary) for boundary in feed_changes
                                     if 0.0 < boundary < hours)})

    def depleted(hour, values):
        return values[1]

    depleted.terminal = True
    depleted.direction = -1
    for start, end in zip(boundaries, boundaries[1:]):
        feed_rate = float(feed.at(start))
        cursor = start
        evaluations = 0

        def derivative(hour, values):
            nonlocal evaluations
            evaluations += 1
            if evaluations > max(200000, 8 * steps):
                raise RuntimeError("batch integration exceeded its evaluation budget; "
                                   "stiff or discontinuous starvation rates could not be resolved")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"integration state must remain finite at {hour:g} h")
            biomass, substrate = values[:2]
            # Substrate cannot go negative, and clamping the RATE rather than the STATE is the
            # difference between a run that stops growing when carbon runs out and one that
            # quietly integrates a negative concentration and reports a titre from it.
            mu, q_s, q_p = checked_rates(substrate, hour) if substrate >= 0.0 else (0.0, 0.0, 0.0)
            return np.array([mu * biomass, feed_rate - q_s * biomass,
                             q_p * biomass, q_s * biomass, biomass])

        while cursor < end:
            if state[1] == 0.0 and feed_rate == 0.0:
                checked_rates(0.0, cursor)
                samples = (t > cursor) & (t <= end)
                trajectory[:, samples] = state[:, None]
                trajectory[4, samples] += state[0] * (t[samples] - cursor)
                state[4] += state[0] * (end - cursor)
                break
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    solution = solve_ivp(
                        derivative, (cursor, end), state, method="RK23",
                        max_step=hours / steps, rtol=rtol, atol=atol,
                        events=depleted, dense_output=True)
            except FloatingPointError as exc:
                raise ValueError("rates or integration state exceeded finite numerical range") from exc
            if not solution.success:
                raise RuntimeError(f"batch integration failed at {solution.t[-1]:g} h: {solution.message}")
            reached = float(solution.t[-1])
            samples = (t > cursor) & (t <= reached)
            if np.any(samples):
                trajectory[:, samples] = solution.sol(t[samples])
            state = solution.y[:, -1].copy()
            if solution.status == 1:
                state[1] = 0.0
                checked_rates(0.0, reached)
                if reached <= cursor:
                    raise RuntimeError("substrate depletion prevents forward integration; rates are discontinuous at starvation")
                at_event = t == reached
                trajectory[:, at_event] = state[:, None]
            if np.any(state < 0.0) or np.any(trajectory[:, samples] < 0.0):
                raise RuntimeError("batch integration could not resolve non-negative inventories")
            cursor = reached

    X, S, P, consumed, exposure = trajectory
    supplied = np.array([initial_substrate_mmol_per_l + feed.integral(0.0, hour) for hour in t])
    ledger_error = S + consumed - supplied
    if not np.all(np.isfinite(trajectory)) or np.any(np.abs(ledger_error) > atol + rtol * supplied):
        raise RuntimeError("batch integration failed its substrate ledger")
    MU = np.array([checked_rates(substrate, hour)[0] for substrate, hour in zip(S, t)])
    return BatchResult(t, X, S, P, MU, float(product_molar_mass_g_per_mol),
                       substrate_molar_mass_g_per_mol, total_feed,
                       float(consumed[-1]), float(exposure[-1]))


simulate_batch = simulate_fixed_volume_comparison
