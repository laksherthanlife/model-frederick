"""Exponential fed-batch: the vessel that replaces the chemostat, with mu measured not assumed.

WHY THIS EXISTS. Every steady-state result in this package divides by a growth rate, and
until now that rate came from a chemostat pump -- ``mu = D``. The chemostat is the one piece
of equipment in the whole design the team does not have and would have to build. This module
supplies the same three things without it, and the case for each is checked in
`tests/test_fedbatch.py` rather than asserted here.

THE CONTROL LAW. Feed an exponentially increasing rate into a single vessel with no effluent
line::

    F(t) = F0 * exp(mu_set * t),    F0 = mu_set * X0 * V0 / (Yxs * Sf)

THE PART THAT IS NOT OBVIOUS, and the reason this works at all. The textbook objection is
that F0 needs the yield Yxs known in advance, and Yxs for these strains is both low and
uncertain. **The objection is wrong about the setpoint.** Write B = X*V for total biomass.
Under carbon limitation every fed molecule is consumed, so q_S = F*Sf/B, and

    d(ln q_S)/dt = mu_set - mu(q_S)

which is a negative feedback loop whose fixed point is mu = mu_set for ANY monotone
q_S -> mu map. Yxs and X0 enter only through F0 -- through the loop's INITIAL CONDITION, not
its fixed point. Get the yield wrong by a factor of six and the culture still converges on
the setpoint; it just starts further from it. That is measured in
:func:`tests.test_fedbatch` against a deliberately wrong yield.

WHAT IT COSTS, and these are the honest limits, all of them checked:

  * THE GENERATION BUDGET IS FINITE. g = log2(1 + (Vmax-V0)*Y*Sf/(X0*V0)), independent of
    mu_set -- the setpoint only sets how fast you spend it. Settling eats part of it,
    logarithmically in how badly Yxs is known. One mu per vessel-run, where a chemostat
    gives many by stepping D. **This is the real loss and it is not recoverable**; the
    compensation is that a run takes days rather than weeks.
  * AND THE YIELD ERROR IS ASYMMETRIC, for a reason that is NOT the settling algebra. The
    asymptotic settling law ``(1 - 1/k) * exp(-mu_set * t)``, k = Y_true/Y_assumed, is
    symmetric enough to be misleading: it says an over-estimate costs MORE generations than
    an under-estimate, and simulation says the opposite happens. The mechanism is the uptake
    ceiling. F0 sets the opening uptake exactly, ``q_S(0) = mu_set / Y_assumed``, with no
    dependence on the vessel at all -- so under-estimating the yield raises the opening
    uptake, and once ``q_S(0) > q_max`` the culture starts ALREADY past its ceiling, grows
    at mu_max rather than at the setpoint, and piles up carbon it cannot take up. It never
    enters the controlled regime, so no amount of settling time helps.
    Over-estimating the yield only lowers the opening uptake, which is always inside the
    ceiling and merely slow. **When the yield is uncertain, guess it HIGH** -- a feed that
    starts too lean converges; one that starts too rich never begins.
    :func:`design_fedbatch` refuses the second case at design time, which it can do because
    ``mu_set / Y_assumed`` is known before anything is inoculated.
  * THE SETPOINT MUST SIT BELOW mu_max WITH MARGIN. The convergence rate is
    lambda = q_S* * (dmu/dq_S) at the setpoint, which collapses as the uptake curve
    flattens. Push mu_set to mu_max and the feed stops setting mu; push past it and glucose
    accumulates without limit. :func:`design_fedbatch` REFUSES rather than returning a
    design that would do this, and it refuses on this repository's own upper calibration
    state -- see `REFUSES_THE_UPPER_ELIZONDO_STATE`.
  * SAMPLING IS A DISTURBANCE A CHEMOSTAT DOES NOT HAVE. A chemostat is sampled from its
    overflow for free; here every sample removes biomass and kicks mu above setpoint.
    Ignoring it costs ~3% RMSE on mu. :func:`reanchor_feed` is the fix and it is a line of
    arithmetic, not a piece of equipment.

AND THE THING THE CHEMOSTAT NEVER GAVE. A chemostat's D = F/V is ASSERTED from pump and
level calibration and is never independently measured. Here mu is measured after the fact
from the biomass series that is already being taken -- :func:`estimate_growth_rate` does OLS
on ln(X*V) and returns a standard error. **Do not assert mu from the pump. Take it from the
biomass.** The published RMSE for feedback-controlled fed-batch in S. cerevisiae is 15 +/- 3%
(Kottelat, Freeland & Dabros 2021, Processes 9:723), but that is the error of an ONLINE
heat-signal ESTIMATE of mu, not the physiological control error, and the two must not be
conflated -- that paper says the derivative "introduc[es] noise to the resulting estimate".
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "REFUSES_THE_UPPER_ELIZONDO_STATE",
    "FedBatchDesign",
    "GrowthRateEstimate",
    "SetpointAboveCapacity",
    "design_fedbatch",
    "estimate_growth_rate",
    "reanchor_feed",
    "simulate_fedbatch",
]

#: Recorded because it is a refusal this repository's own calibration set walks into. The
#: upper Elizondo state is mu = 0.254 /h at q_glucose = 2.59 g glucose/gDCW/h, and 0.254 /h
#: IS that strain's measured uptake ceiling -- the state sits ON mu_max with no margin. An
#: exponential fed-batch reproduces the lower state (mu = 0.101 /h) comfortably and CANNOT be
#: trusted to reproduce the upper one: at the ceiling the feed no longer sets mu.
REFUSES_THE_UPPER_ELIZONDO_STATE = (
    "mu = 0.254 /h is the upper Elizondo chemostat state and also that strain's measured "
    "mu_max, so a fed-batch designed for it has zero margin and its feed does not set mu"
)


class SetpointAboveCapacity(ValueError):
    """The setpoint is at or above what the strain can grow at, so the feed cannot set mu.

    Its own exception type because this is the ONE failure mode that turns an exponential
    fed-batch from a controller into an uncontrolled accumulation, and a caller should be
    able to catch it specifically rather than by matching on a message.
    """


@dataclass(frozen=True)
class GrowthRateEstimate:
    """A growth rate MEASURED from the biomass series, with the uncertainty it was measured to.

    The point of the whole module. A chemostat reports D with no error bar because D is an
    input; this is an output, so it has one.
    """

    mu_per_h: float
    standard_error_per_h: float
    points: int
    span_h: float

    @property
    def relative_standard_error(self) -> float:
        return self.standard_error_per_h / self.mu_per_h if self.mu_per_h else float("inf")

    def summary(self) -> str:
        return (f"mu = {self.mu_per_h:.4f} +/- {self.standard_error_per_h:.4f} /h "
                f"({100 * self.relative_standard_error:.2f}%) from {self.points} points "
                f"over {self.span_h:.1f} h")


@dataclass(frozen=True)
class FedBatchDesign:
    """A vessel plan: what to feed, for how long, and how much mu it can actually deliver."""

    mu_setpoint_per_h: float
    initial_volume_l: float
    max_volume_l: float
    initial_biomass_g_per_l: float
    feed_substrate_g_per_l: float
    assumed_yield_g_per_g: float
    mu_max_per_h: float

    @property
    def initial_feed_rate_l_per_h(self) -> float:
        """F0 = mu_set * X0 * V0 / (Yxs * Sf), in litres of feed per hour."""
        return (self.mu_setpoint_per_h * self.initial_biomass_g_per_l * self.initial_volume_l
                / (self.assumed_yield_g_per_g * self.feed_substrate_g_per_l))

    def feed_rate_l_per_h(self, hours: float) -> float:
        """The open-loop profile. Exponential, and that is the whole controller."""
        return self.initial_feed_rate_l_per_h * math.exp(self.mu_setpoint_per_h * hours)

    def generation_budget_at(self, yield_g_per_g: float) -> float:
        """Doublings the vessel holds at a given TRUE yield, independent of the setpoint.

        g = log2(1 + (Vmax - V0) * Y * Sf / (X0 * V0)). The setpoint sets how fast the budget
        is spent, never how large it is, which is why a low inoculum is the lever: X0 is in
        the denominator and nothing else on the right is free.

        Takes the yield as an ARGUMENT because the budget is a physical quantity set by what
        the strain actually does, while :attr:`assumed_yield_g_per_g` is what the planner
        guessed when setting F0. Conflating the two makes the budget appear to shrink when
        you under-estimate the yield, which is backwards -- under-estimating gives you MORE
        biomass per litre of feed than planned, not less.
        """
        if yield_g_per_g <= 0:
            raise ValueError(f"yield must be positive, got {yield_g_per_g}")
        addable = self.max_volume_l - self.initial_volume_l
        biomass_makeable = addable * yield_g_per_g * self.feed_substrate_g_per_l
        return math.log2(1.0 + biomass_makeable
                         / (self.initial_biomass_g_per_l * self.initial_volume_l))

    @property
    def generation_budget(self) -> float:
        """The PLANNED budget, at the assumed yield. See :meth:`generation_budget_at`."""
        return self.generation_budget_at(self.assumed_yield_g_per_g)

    @property
    def final_biomass_g_per_l(self) -> float:
        """X at the moment the vessel is full. Independent of the setpoint AND of X0."""
        return ((self.max_volume_l - self.initial_volume_l)
                * self.assumed_yield_g_per_g * self.feed_substrate_g_per_l
                / self.max_volume_l)

    def settling_generations(self, yield_error_factor: float,
                             tolerance: float = 0.02) -> float:
        """Generations spent converging, given a yield wrong by ``yield_error_factor``.

        The relative error in mu decays as ``(1 - 1/k) * exp(-mu_set * t)`` with
        ``k = Yxs_true / Yxs_assumed``, so the time to a tolerance is
        ``ln(|1 - 1/k| / tol) / mu_set`` and the GENERATIONS are that over ln(2) -- which
        removes mu_set entirely. **The cost of not knowing the yield is logarithmic**, which
        is why a factor-of-two error is affordable and a factor of ten still is.
        """
        if yield_error_factor <= 0:
            raise ValueError(f"yield_error_factor must be positive, got {yield_error_factor}")
        offset = abs(1.0 - 1.0 / yield_error_factor)
        if offset <= tolerance:
            return 0.0
        return math.log(offset / tolerance) / math.log(2.0)

    def usable_generations(self, yield_error_factor: float,
                           tolerance: float = 0.02) -> float:
        """Budget minus settling. What is left to actually measure in."""
        return self.generation_budget - self.settling_generations(
            yield_error_factor, tolerance)

    @property
    def initial_uptake_g_per_gdcw_h(self) -> float:
        """q_S(0) = mu_set / Y_assumed. The opening uptake, and the whole safety question.

        Note what is ABSENT: the vessel. F0 scales with X0*V0 and q_S(0) divides by the same
        X0*V0, so the opening uptake depends only on the setpoint and the assumed yield.
        A bigger vessel or a heavier inoculum cannot rescue a feed that starts too rich.
        """
        return self.mu_setpoint_per_h / self.assumed_yield_g_per_g

    def worst_tolerable_yield_overestimate(self, max_uptake_g_per_gdcw_h: float) -> float:
        """The LOWEST assumed yield whose opening uptake still fits under the ceiling.

        Solving ``mu_set / Y_assumed <= q_max`` for Y_assumed. Anything below this starts the
        culture past its uptake ceiling, which is the one error the controller cannot
        recover from. Reported as a yield rather than a ratio because the yield is what the
        planner actually types in.
        """
        if max_uptake_g_per_gdcw_h <= 0:
            raise ValueError(f"max_uptake must be positive, got {max_uptake_g_per_gdcw_h}")
        return self.mu_setpoint_per_h / max_uptake_g_per_gdcw_h

    @property
    def capacity_margin(self) -> float:
        """How far the setpoint sits below mu_max, as a fraction of mu_max."""
        return (self.mu_max_per_h - self.mu_setpoint_per_h) / self.mu_max_per_h

    def summary(self) -> str:
        return (f"mu_set {self.mu_setpoint_per_h:.3f} /h ({100 * self.capacity_margin:.0f}% "
                f"below mu_max {self.mu_max_per_h:.3f}); F0 "
                f"{1000 * self.initial_feed_rate_l_per_h:.2f} mL/h; "
                f"{self.generation_budget:.2f} generations in the vessel; "
                f"X_final {self.final_biomass_g_per_l:.1f} g/L")


def design_fedbatch(mu_setpoint_per_h: float,
                    mu_max_per_h: float,
                    *,
                    initial_volume_l: float,
                    max_volume_l: float,
                    initial_biomass_g_per_l: float,
                    feed_substrate_g_per_l: float,
                    assumed_yield_g_per_g: float,
                    max_uptake_g_per_gdcw_h: float,
                    minimum_margin: float = 0.10) -> FedBatchDesign:
    """Plan a run, or refuse to.

    Args:
        mu_setpoint_per_h: The growth rate to hold. This is the chemostat's D.
        mu_max_per_h: What the strain can actually do, MEASURED on this strain in this
            medium. There is no default: the failure mode this argument exists to prevent is
            precisely a setpoint chosen without reference to the strain.
        max_uptake_g_per_gdcw_h: The strain's uptake ceiling, measured. Required for the
            same reason and against a DIFFERENT failure: a setpoint can sit safely below
            mu_max and still open with a feed too rich to take up, if the assumed yield is
            low. Simulation found this one -- it is not the settling algebra, and the
            settling algebra points the other way.
        minimum_margin: How far below mu_max the setpoint must sit, as a fraction. The
            default 10% is a judgement, not a measurement, and is stated as one -- the
            convergence rate lambda = q_S * dmu/dq_S falls off smoothly, so there is no
            sharp threshold to read off. What IS measured is that at zero margin the feed
            fails to set mu at all.

    Raises:
        SetpointAboveCapacity: when the setpoint has less than ``minimum_margin`` headroom.
        ValueError: on a non-positive or self-contradictory vessel.
    """
    if mu_setpoint_per_h <= 0:
        raise ValueError(f"mu_setpoint must be positive, got {mu_setpoint_per_h}")
    if mu_max_per_h <= 0:
        raise ValueError(f"mu_max must be positive, got {mu_max_per_h}")
    if not 0.0 <= minimum_margin < 1.0:
        raise ValueError(f"minimum_margin is a fraction in [0, 1), got {minimum_margin}")
    if max_volume_l <= initial_volume_l:
        raise ValueError(
            f"max_volume_l ({max_volume_l}) must exceed initial_volume_l "
            f"({initial_volume_l}); the difference is the entire generation budget")
    for name, value in (("initial_volume_l", initial_volume_l),
                        ("initial_biomass_g_per_l", initial_biomass_g_per_l),
                        ("feed_substrate_g_per_l", feed_substrate_g_per_l),
                        ("assumed_yield_g_per_g", assumed_yield_g_per_g)):
        if value <= 0:
            raise ValueError(f"{name} must be positive, got {value}")

    if max_uptake_g_per_gdcw_h <= 0:
        raise ValueError(f"max_uptake must be positive, got {max_uptake_g_per_gdcw_h}")

    opening_uptake = mu_setpoint_per_h / assumed_yield_g_per_g
    if opening_uptake > max_uptake_g_per_gdcw_h:
        raise SetpointAboveCapacity(
            f"the feed opens at q_S(0) = mu_set/Y = {opening_uptake:.3g} g/gDCW/h, above "
            f"the strain's uptake ceiling {max_uptake_g_per_gdcw_h:g}. The culture would "
            f"start ALREADY past its ceiling, grow at mu_max instead of the setpoint, and "
            f"accumulate the carbon it cannot take up -- it never enters the regime the "
            f"exponential feed controls, so no run length fixes it. Note the vessel cannot "
            f"help: F0 scales with X0*V0 and q_S(0) divides it out again. Raise the assumed "
            f"yield to at least {mu_setpoint_per_h / max_uptake_g_per_gdcw_h:.4g} g/g, or "
            f"lower the setpoint. When the yield is uncertain, guess it HIGH: a lean feed "
            f"converges, a rich one never begins")

    margin = (mu_max_per_h - mu_setpoint_per_h) / mu_max_per_h
    if margin < minimum_margin:
        raise SetpointAboveCapacity(
            f"mu_setpoint {mu_setpoint_per_h:g} /h sits {100 * margin:.1f}% below mu_max "
            f"{mu_max_per_h:g} /h, under the {100 * minimum_margin:.0f}% margin this design "
            f"requires. An exponential feed sets mu through the negative feedback "
            f"d(ln q_S)/dt = mu_set - mu(q_S), whose convergence rate is q_S * dmu/dq_S -- "
            f"that goes to zero as the uptake curve flattens, so near mu_max the feed stops "
            f"setting mu, and past it glucose accumulates without bound. "
            f"Lower the setpoint, or measure a higher mu_max on this strain. "
            f"For the record: {REFUSES_THE_UPPER_ELIZONDO_STATE}")

    return FedBatchDesign(
        mu_setpoint_per_h=mu_setpoint_per_h,
        initial_volume_l=initial_volume_l,
        max_volume_l=max_volume_l,
        initial_biomass_g_per_l=initial_biomass_g_per_l,
        feed_substrate_g_per_l=feed_substrate_g_per_l,
        assumed_yield_g_per_g=assumed_yield_g_per_g,
        mu_max_per_h=mu_max_per_h)


def estimate_growth_rate(times_h, biomass_g_per_l, volumes_l=None) -> GrowthRateEstimate:
    """mu and its standard error, by OLS on ln(X*V). The measurement a chemostat never makes.

    TOTAL biomass X*V, not concentration: the vessel is filling, so X alone understates mu by
    exactly the volume's own growth. Passing no volumes asserts a constant one, which is a
    batch and not a fed-batch -- allowed, because the same estimator is the right one for a
    batch growth curve, but it is an assertion the caller is making.

    Args:
        times_h: Sample times. At least three, because two points have no residual and
            therefore no standard error -- and the standard error is the point.

    Raises:
        ValueError: fewer than three points, mismatched lengths, non-positive biomass, or
            times that do not advance.
    """
    times = [float(t) for t in times_h]
    biomass = [float(x) for x in biomass_g_per_l]
    if volumes_l is None:
        volumes = [1.0] * len(times)
    else:
        volumes = [float(v) for v in volumes_l]

    if not len(times) == len(biomass) == len(volumes):
        raise ValueError(
            f"times, biomass and volumes must be the same length, got "
            f"{len(times)}, {len(biomass)}, {len(volumes)}")
    if len(times) < 3:
        raise ValueError(
            f"need at least 3 points to estimate mu with a standard error, got "
            f"{len(times)}. Two points fit a line exactly and report no uncertainty, which "
            f"is the one thing this estimator exists to provide")
    if any(x <= 0 for x in biomass) or any(v <= 0 for v in volumes):
        raise ValueError("biomass and volume must be positive to take a logarithm")
    if len(set(times)) != len(times):
        raise ValueError("sample times must be distinct")

    total = [math.log(x * v) for x, v in zip(biomass, volumes)]
    n = len(times)
    mean_t = sum(times) / n
    mean_y = sum(total) / n
    sxx = sum((t - mean_t) ** 2 for t in times)
    if sxx == 0:
        raise ValueError("sample times must advance; every point is at the same time")
    slope = sum((t - mean_t) * (y - mean_y) for t, y in zip(times, total)) / sxx
    intercept = mean_y - slope * mean_t
    residuals = [y - (intercept + slope * t) for t, y in zip(times, total)]
    # n - 2 for the two fitted parameters. With exactly 3 points this is 1 degree of freedom,
    # which is thin but real -- and thin is the honest report, not a reason to hide it.
    dof = n - 2
    variance = sum(r * r for r in residuals) / dof
    return GrowthRateEstimate(mu_per_h=slope,
                              standard_error_per_h=math.sqrt(variance / sxx),
                              points=n,
                              span_h=max(times) - min(times))


def reanchor_feed(design: FedBatchDesign,
                  measured_biomass_g_per_l: float,
                  volume_l: float) -> float:
    """The feed rate that puts a SAMPLED vessel back on setpoint, L/h.

    A chemostat is sampled from its overflow for free. Here a sample removes biomass, which
    raises q_S = F*Sf/B for the biomass that remains and pushes mu above setpoint. Ignoring
    it costs a few percent RMSE on mu; re-anchoring the exponential on the biomass actually
    measured at each sample costs nothing, because the biomass is being measured anyway.

    This is F = mu_set * B / (Y * Sf) evaluated at the CURRENT total biomass rather than at
    the extrapolated one -- the same law, with the feedback closed by hand once per sample.
    """
    if measured_biomass_g_per_l <= 0 or volume_l <= 0:
        raise ValueError(
            f"biomass and volume must be positive, got {measured_biomass_g_per_l} and "
            f"{volume_l}")
    total_biomass_g = measured_biomass_g_per_l * volume_l
    return (design.mu_setpoint_per_h * total_biomass_g
            / (design.assumed_yield_g_per_g * design.feed_substrate_g_per_l))


def simulate_fedbatch(design: FedBatchDesign,
                      uptake_to_growth,
                      *,
                      true_yield_g_per_g: float,
                      max_uptake_g_per_gdcw_h: float,
                      hours: float,
                      steps: int = 4000):
    """Integrate the vessel, INCLUDING volume, and report whether mu found the setpoint.

    `fba/dynamic.py` integrates a batch at constant volume and cannot represent this: the
    whole mechanism is that V grows, so total biomass X*V grows and q_S = F*Sf/(X*V) closes
    the loop. That is why this is a separate integrator rather than a flag on that one.

    Args:
        uptake_to_growth: ``q_S (g glucose/gDCW/h) -> mu (/h)``. The strain's physiology, and
            deliberately arbitrary: the fixed-point argument holds for ANY monotone map, so
            the tests feed it a Crabtree kink to prove exactly that.
        true_yield_g_per_g: What the strain's yield REALLY is, as against
            ``design.assumed_yield_g_per_g`` used to set F0. Separating the two is the point
            of the simulator -- it is how "the yield barely matters" gets tested rather than
            repeated.
        max_uptake_g_per_gdcw_h: The strain's uptake ceiling. Required, because it is what
            makes the failure mode VISIBLE: offered carbon above this accumulates in the
            vessel instead of vanishing into an unbounded uptake, and that accumulation is
            the signature of a setpoint above capacity.

    Returns:
        A dict with the trajectory and the settled growth rate.
    """
    if hours <= 0:
        raise ValueError(f"hours must be positive, got {hours}")
    if steps < 2:
        raise ValueError(f"need at least 2 steps, got {steps}")
    if true_yield_g_per_g <= 0:
        raise ValueError(f"true_yield_g_per_g must be positive, got {true_yield_g_per_g}")
    if max_uptake_g_per_gdcw_h <= 0:
        raise ValueError(
            f"max_uptake must be positive, got {max_uptake_g_per_gdcw_h}")

    dt = hours / steps
    filled = False
    volume = design.initial_volume_l
    biomass = design.initial_biomass_g_per_l
    substrate_g = 0.0          # unconsumed carbon in the vessel, as a MASS
    times, mus, volumes, biomasses, substrates, uptakes = [], [], [], [], [], []

    for step in range(steps + 1):
        t = step * dt
        feed = design.feed_rate_l_per_h(t)
        total_biomass = biomass * volume
        supplied_g_per_h = feed * design.feed_substrate_g_per_l

        # Everything the cells could reach this step: what arrives, plus what is already
        # standing in the vessel. Capped by what they can actually take up -- the cap is
        # what lets carbon accumulate rather than disappear.
        offered_g_per_h = supplied_g_per_h + substrate_g / dt
        q_s = min(offered_g_per_h / total_biomass, max_uptake_g_per_gdcw_h)
        mu = uptake_to_growth(q_s)

        times.append(t)
        mus.append(mu)
        volumes.append(volume)
        biomasses.append(biomass)
        substrates.append(substrate_g / volume)
        uptakes.append(q_s)

        if step == steps:
            break

        consumed_g = q_s * total_biomass * dt
        substrate_g = max(0.0, substrate_g + supplied_g_per_h * dt - consumed_g)
        total_biomass = total_biomass + mu * total_biomass * dt
        volume = volume + feed * dt
        biomass = total_biomass / volume
        if volume >= design.max_volume_l:
            # Record the state that ENDED the run. Reading `vessel_full` off the last
            # appended volume without this reports False on every run that filled, because
            # the break lands after the append rather than before it.
            filled = True
            times.append(t + dt)
            mus.append(mu)
            volumes.append(volume)
            biomasses.append(biomass)
            substrates.append(substrate_g / volume)
            uptakes.append(q_s)
            break
    else:
        filled = volume >= design.max_volume_l

    settled = mus[-1]
    return {
        "times_h": times,
        "growth_rate_per_h": mus,
        "volume_l": volumes,
        "biomass_g_per_l": biomasses,
        "substrate_g_per_l": substrates,
        "uptake_g_per_gdcw_h": uptakes,
        "settled_growth_rate_per_h": settled,
        "setpoint_error_relative": (settled - design.mu_setpoint_per_h)
                                   / design.mu_setpoint_per_h,
        "vessel_full": filled,
        "residual_substrate_g_per_l": substrates[-1],
    }
