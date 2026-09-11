"""Pin the LP solver, so a flux number means the same thing on two machines.

cobrapy picks its solver from whatever optlang finds installed. This package pins
``swiglpk``, so on a clean install that resolves to GLPK -- but a machine that also has
Gurobi or CPLEX silently gets a different simplex, and a different simplex breaks ties
among alternate optima differently. Nothing in the code would notice.

That matters here more than it looks. ``fba/fva.py`` exists because the product flux
range is *wide*: D1 measured relative width 1.000, meaning the feasible interval runs
from zero to the ceiling at every growth level. A wide interval is precisely the case
where the particular optimum a solver returns is arbitrary, so any single flux read off
one solve is a property of the solver as much as of the model.

The tolerance is pinned at cobrapy's own default rather than tightened. The goal is that
two machines agree, not that the numbers move; changing it would silently restate every
FBA result in the repository.
"""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass

import cobra

__all__ = [
    "PINNED_SOLVER",
    "PINNED_TOLERANCE",
    "SolverSettings",
    "configure",
    "load_model",
    "growth_or_none",
]

# GLPK because swiglpk is the pinned dependency; see pyproject.toml.
PINNED_SOLVER = "glpk"
PINNED_TOLERANCE = 1e-7

# FVA forks a worker pool by default. Worker count changes nothing about the answer in
# principle and everything about reproducing a run, so it is fixed at one.
FVA_PROCESSES = 1


@dataclass(frozen=True)
class SolverSettings:
    """What was actually configured, for stamping into a results table.

    Every table this repository publishes records the conditions it was produced under.
    A flux table without its solver is missing one of them.
    """

    solver: str
    tolerance: float
    cobra_version: str

    def __str__(self) -> str:
        return (f"{self.solver} tol={self.tolerance:g} cobra={self.cobra_version}")


def growth_or_none(model: cobra.Model) -> float | None:
    """Maximise the current objective, returning ``None`` when the LP does not solve.

    Exists because ``slim_optimize()`` signals failure with **nan**, not ``None``, and both
    of the obvious guards are silently wrong about nan:

        value is None   ->  False, so the "infeasible" branch is dead code
        value <= 0      ->  False, so a "cannot grow" check lets it straight through

    Three call sites in this package had one or the other. The audit layer's
    "refuted by stoichiometry rather than merely expensive" note could never fire, and a nan
    propagated into a reported growth cost instead. ``error_value=None`` is not the fix
    either -- it makes cobra raise ``Infeasible`` rather than return anything -- so the
    check has to be an explicit ``isnan``.

    Returns:
        The objective value, or ``None`` if the model is infeasible or unbounded.
    """
    value = model.slim_optimize()
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    return float(value)


def configure(model: cobra.Model, tolerance: float = PINNED_TOLERANCE) -> SolverSettings:
    """Pin ``model``'s solver and tolerance in place, and report what was set.

    Args:
        model: Model to configure.
        tolerance: Feasibility/optimality tolerance. Defaults to the pinned value.

    Raises:
        RuntimeError: if the pinned solver is unavailable. Falling back to whatever is
            installed is what this function exists to prevent -- a silent fallback would
            reintroduce the machine dependence while looking configured.
    """
    try:
        model.solver = PINNED_SOLVER
    except Exception as exc:  # optlang raises its own types for a missing interface
        raise RuntimeError(
            f"the pinned solver {PINNED_SOLVER!r} is unavailable ({exc}). Install it "
            "rather than letting cobrapy choose: a different simplex breaks ties among "
            "alternate optima differently, and the product flux range here is wide "
            "enough that the tie is the answer."
        ) from exc
    model.tolerance = tolerance
    return SolverSettings(PINNED_SOLVER, tolerance, cobra.__version__)


def load_model(path: str | pathlib.Path,
               tolerance: float = PINNED_TOLERANCE) -> tuple[cobra.Model, SolverSettings]:
    """Read an SBML model and pin its solver before anything solves with it.

    Returns the model and the settings, so a caller can record them beside its numbers.
    """
    model = cobra.io.read_sbml_model(str(path))
    return model, configure(model, tolerance)
