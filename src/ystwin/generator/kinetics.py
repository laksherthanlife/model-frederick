"""When a reporter reaches the level a steady-state model assumes it already has.

Every simulated reading was a steady state, and a real plate is read on a timetable.
Reporter accumulates as ``dR/dt = k - (mu + k_deg) R``, so it approaches ``k/(mu + k_deg)``
with a time constant of ``1/(mu + k_deg)``. Reading before that is reading a fraction of
the response.

The bias runs the wrong way, which is why it has to be modelled rather than noted. The
relaxation rate is the growth rate, so a stressed culture -- slow-growing, the one the
experiment is about -- takes longest to arrive. An early read understates stress most in
exactly the wells where stress is greatest.

A ratiometric sensor carries none of this. It reports an equilibrium between two forms of
one molecule and settles in seconds, so it is at its final value whenever the plate is read.
"""

from __future__ import annotations

import numpy as np

__all__ = ["reporter_at_time", "time_to_fraction"]


def reporter_at_time(
    time_h: float, steady_state: float, growth_rate: float,
    initial: float = 0.0, k_deg: float = 0.0,
) -> float:
    """Reporter level at one time, from the closed form of its own balance.

    Args:
        time_h: Hours since the step in promoter activity.
        steady_state: Level the reporter would reach given unlimited time.
        growth_rate: Specific growth rate, which is also the dilution rate.
        initial: Level when the step happened.
        k_deg: Loss beyond dilution; zero for a stable fluorescent protein.
    """
    if time_h < 0:
        raise ValueError(f"cannot read at a negative time, got {time_h}")
    relaxation = float(growth_rate) + float(k_deg)
    if relaxation <= 0:
        return float(initial)
    approach = 1.0 - np.exp(-relaxation * float(time_h))
    return float(initial + (steady_state - initial) * approach)


def time_to_fraction(fraction: float, growth_rate: float, k_deg: float = 0.0) -> float:
    """Hours to cover a given fraction of the distance to steady state.

    The practical form of the same fact: it answers when the plate can be read, and it
    answers longer for every culture the stressor has slowed.
    """
    if not 0.0 < fraction < 1.0:
        raise ValueError(f"fraction must be between 0 and 1, got {fraction}")
    relaxation = float(growth_rate) + float(k_deg)
    if relaxation <= 0:
        return float("inf")
    return float(-np.log(1.0 - fraction) / relaxation)
