"""Emulator for the constraint-to-flux map, so state estimation stays affordable.

A particle filter needs ~1e5 LP solves per culture. The map from uptake bounds to optimal fluxes
is piecewise linear, so grid interpolation is exact inside each piece. It refuses to extrapolate:
outside the sampled box the active basis is unknown, and a confidently wrong flux is worse than
none.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import cobra
import numpy as np
from cobra.exceptions import OptimizationError
from scipy.interpolate import RegularGridInterpolator

__all__ = ["FluxSurrogate", "build_surrogate"]

_MAX_GRID_POINTS = 200_000


@dataclass(frozen=True)
class FluxSurrogate:
    """Interpolated stand-in for repeated LP solves over a box of uptake bounds."""

    inputs: dict[str, tuple[float, float]]
    outputs: tuple[str, ...]
    axes: tuple[np.ndarray, ...]
    interpolators: dict[str, RegularGridInterpolator] = field(repr=False)
    n_samples: int = 0
    n_infeasible: int = 0

    def _as_array(self, point: dict[str, float]) -> np.ndarray:
        missing = set(self.inputs) - set(point)
        if missing:
            raise KeyError(f"missing input bounds for {sorted(missing)}")
        values = []
        for name in self.inputs:
            lo, hi = self.inputs[name]
            v = np.asarray(point[name], dtype=float)
            if np.any(v < min(lo, hi) - 1e-9) or np.any(v > max(lo, hi) + 1e-9):
                raise ValueError(
                    f"{name}={point[name]} is outside the sampled range "
                    f"[{min(lo, hi)}, {max(lo, hi)}]; the surrogate does not extrapolate"
                )
            values.append(v)
        return np.stack(np.broadcast_arrays(*values), axis=-1)

    def predict(self, point: dict[str, float]) -> dict[str, float]:
        """Predicted fluxes at one set of uptake bounds."""
        query = self._as_array(point).reshape(1, -1)
        return {name: float(self.interpolators[name](query)[0]) for name in self.outputs}

    def predict_many(self, points: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Predicted fluxes at many sets of bounds at once."""
        query = self._as_array(points)
        return {name: np.asarray(self.interpolators[name](query)) for name in self.outputs}

    def summary(self) -> str:
        shape = "x".join(str(len(a)) for a in self.axes)
        return (
            f"surrogate over {list(self.inputs)} on a {shape} grid "
            f"({self.n_samples} solves, {self.n_infeasible} infeasible) "
            f"-> {list(self.outputs)}"
        )


def build_surrogate(
    model: cobra.Model,
    inputs: dict[str, tuple[float, float]],
    outputs: list[str],
    resolution: int = 7,
    objective: str | None = None,
    infeasible_fill: float = 0.0,
) -> FluxSurrogate:
    """Sample the LP on a grid of uptake bounds and build an interpolator.

    Args:
        model: Model to sample. Not modified.
        inputs: Reaction id -> (low, high) range for its lower bound.
        outputs: Reaction ids whose optimal flux should be emulated.
        resolution: Grid points per input axis. Cost grows as ``resolution ** len(inputs)``,
            so this is practical for a handful of inputs, not dozens.
        objective: Objective reaction; defaults to the model's own.
        infeasible_fill: Value recorded where the LP has no solution. Counted and
            reported rather than hidden -- an interpolation that straddles an
            infeasible node is not trustworthy.
    """
    for rid in list(inputs) + list(outputs):
        if rid not in model.reactions:
            raise KeyError(f"reaction {rid!r} is not in this model")
    if resolution < 2:
        raise ValueError("resolution must be at least 2")
    n_points = resolution ** len(inputs)
    if n_points > _MAX_GRID_POINTS:
        raise ValueError(
            f"{len(inputs)} inputs at resolution {resolution} needs {n_points} LP solves; "
            f"reduce resolution or the number of inputs (cap {_MAX_GRID_POINTS})"
        )

    axes = tuple(np.linspace(lo, hi, resolution) for lo, hi in inputs.values())
    mesh = np.meshgrid(*axes, indexing="ij")
    flat = np.stack([m.ravel() for m in mesh], axis=-1)

    results = {name: np.empty(flat.shape[0], dtype=float) for name in outputs}
    n_infeasible = 0
    with model as m:
        if objective is not None:
            m.objective = objective
        for i, row in enumerate(flat):
            for rid, value in zip(inputs, row):
                m.reactions.get_by_id(rid).lower_bound = float(value)
            # raise_error=True keeps infeasibility on the exception path.
            try:
                solution = m.optimize(raise_error=True)
            except OptimizationError:
                n_infeasible += 1
                for name in outputs:
                    results[name][i] = infeasible_fill
                continue
            for name in outputs:
                results[name][i] = solution.fluxes[name]

    grid_shape = tuple(len(a) for a in axes)
    interpolators = {
        name: RegularGridInterpolator(
            axes, values.reshape(grid_shape), method="linear", bounds_error=False, fill_value=None
        )
        for name, values in results.items()
    }
    return FluxSurrogate(
        inputs=dict(inputs),
        outputs=tuple(outputs),
        axes=axes,
        interpolators=interpolators,
        n_samples=int(flat.shape[0]),
        n_infeasible=n_infeasible,
    )
