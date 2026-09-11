"""Where a well sits on a plate, and why the generator has to know.

Edge wells evaporate faster than interior ones. The medium concentrates, the stressor
concentrates with it, and the culture behaves as though it had been dosed higher than the
pipette says. It is the most reproducible artifact in plate work, and the reason plate maps
keep controls out of the outer ring.

It cannot be folded into the noise term. Measurement noise is independent per well, so
replicates average it away; a position effect belongs to the position, so replicates in the
same ring inherit the same bias and averaging does nothing at all. A generator carrying
noise but no geometry overstates what replication buys.
"""

from __future__ import annotations

import numpy as np

__all__ = ["edge_multiplier", "edge_wells", "well_grid", "well_positions"]

_ROWS = 8
_COLUMNS = 12


def well_grid(rows: int = _ROWS, columns: int = _COLUMNS) -> list[str]:
    """Well names in fill order, A1 through the last row and column."""
    letters = [chr(ord("A") + r) for r in range(rows)]
    return [f"{letter}{c + 1}" for letter in letters for c in range(columns)]


def edge_wells(rows: int = _ROWS, columns: int = _COLUMNS) -> frozenset[str]:
    """The outer ring, where evaporation concentrates the medium."""
    letters = [chr(ord("A") + r) for r in range(rows)]
    ring = set()
    for r, letter in enumerate(letters):
        for c in range(columns):
            if r in (0, rows - 1) or c in (0, columns - 1):
                ring.add(f"{letter}{c + 1}")
    return frozenset(ring)


def well_positions(n_samples: int, rows: int = _ROWS, columns: int = _COLUMNS) -> list[str]:
    """Assign samples to wells in fill order.

    Deterministic on purpose. A design that puts every replicate of one treatment down a
    single column has confounded treatment with position, and that is only visible if the
    assignment is inspectable rather than drawn at random.
    """
    grid = well_grid(rows, columns)
    if n_samples > len(grid):
        raise ValueError(f"{n_samples} samples but a plate holds only {len(grid)}")
    return grid[:n_samples]


def edge_multiplier(wells, effect: float, rows: int = _ROWS, columns: int = _COLUMNS) -> np.ndarray:
    """Per-well scaling from evaporation, one for interior wells."""
    ring = edge_wells(rows, columns)
    return np.array([1.0 + effect if w in ring else 1.0 for w in wells])
