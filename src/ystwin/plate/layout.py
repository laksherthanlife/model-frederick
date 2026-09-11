"""Which well held which construct at which dose.

No plate map ships with the exports, so it is recovered: the per-construct sheets hold
(RFU-blank)/(OD-blank) per dose, and the construct-to-column assignment reproducing them
identifies the layout. The margin over the runner-up says how firmly. RECORDED_PLATES carries
the layouts transcribed from the logbook, which confirmed the recovery.
"""

from __future__ import annotations

import itertools
import pathlib
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

__all__ = [
    "NEWPROTOCOL_LAYOUT",
    "RECORDED_PLATES",
    "RecordedPlate",
    "blanks_for_export",
    "plate_key",
    "recorded_for_export",
    "recorded_well_roles",
    "LayoutRecovery",
    "PlateLayout",
    "recover_layout",
    "wells_for",
]


@dataclass(frozen=True)
class PlateLayout:
    """Construct-to-column-block and dose-to-row assignment for one plate design.

    Args:
        construct_columns: Construct -> the plate columns holding its replicates.
        dose_rows: Plate rows in ascending dose order.
    """

    construct_columns: dict[str, tuple[int, ...]]
    dose_rows: tuple[str, ...]

    def wells(self, construct: str, dose_index: int) -> list[str]:
        if construct not in self.construct_columns:
            raise KeyError(
                f"no columns recorded for {construct!r}; "
                f"have {sorted(self.construct_columns)}"
            )
        if not 0 <= dose_index < len(self.dose_rows):
            raise IndexError(
                f"dose index {dose_index} is outside this plate's "
                f"{len(self.dose_rows)} dose rows"
            )
        row = self.dose_rows[dose_index]
        return [f"{row}{c}" for c in self.construct_columns[construct]]


# Recovered from the 2026-07-22 plate; see module docstring for the fit quality.
NEWPROTOCOL_LAYOUT = PlateLayout(
    construct_columns={
        "UPRE1": (1, 2, 3),
        "UPRE2": (4, 5, 6),
        "NativeYap1": (7, 8, 9),
        "AlteredYap1": (10, 11, 12),
    },
    dose_rows=tuple("ABCDEFG"),
)


def wells_for(construct: str, dose_index: int, layout: PlateLayout = NEWPROTOCOL_LAYOUT) -> list[str]:
    """Wells holding one construct at one dose index."""
    return layout.wells(construct, dose_index)


@dataclass(frozen=True)
class LayoutRecovery:
    """A recovered layout with the evidence for it."""

    layout: PlateLayout
    median_relative_error: float
    runner_up_error: float
    n_conditions: int
    ranking: list[tuple[float, dict[str, int]]] = field(repr=False, default_factory=list)

    @property
    def margin_over_runner_up(self) -> float:
        """How many times better the winner fits than the next assignment."""
        if self.median_relative_error <= 0:
            return float("inf")
        return self.runner_up_error / self.median_relative_error

    def summary(self) -> str:
        cols = {c: v[0] for c, v in self.layout.construct_columns.items()}
        return (
            f"first columns {cols}, rows {''.join(self.layout.dose_rows)} ascending; "
            f"fit {self.median_relative_error:.2%} over {self.n_conditions} conditions, "
            f"{self.margin_over_runner_up:.1f}x better than the runner-up"
        )


def recover_layout(
    normalised: pd.DataFrame,
    derived: pd.DataFrame,
    block_width: int = 3,
) -> LayoutRecovery:
    """Find the construct-to-column assignment that reproduces the derived sheets.

    Args:
        normalised: Time-indexed ``(RFU - blank)/(OD - blank)`` per well.
        derived: Long frame with ``construct``, ``dose_mM``, ``time_h``, ``signal``
            from the per-construct sheets.
        block_width: Replicate wells per condition.
    """
    constructs = sorted(derived.construct.unique())
    rows = sorted({w[0] for w in normalised.columns})
    n_dose_rows = max(len(derived[derived.construct == c].dose_mM.unique()) for c in constructs)
    dose_rows = tuple(rows[:n_dose_rows])
    blocks = [
        tuple(range(start, start + block_width))
        for start in range(1, block_width * len(constructs) + 1, block_width)
    ]
    grid = normalised.index.to_numpy(dtype=float)

    def fit(assignment: dict[str, tuple[int, ...]]) -> tuple[float, int]:
        errors = []
        for construct, cols in assignment.items():
            doses = sorted(derived[derived.construct == construct].dose_mM.unique())
            for row, dose in zip(dose_rows, doses):
                wells = [f"{row}{c}" for c in cols if f"{row}{c}" in normalised.columns]
                if not wells:
                    continue
                observed = normalised[wells].mean(axis=1).to_numpy()
                g = derived[
                    (derived.construct == construct) & (derived.dose_mM == dose)
                ].sort_values("time_h")
                expected = np.interp(grid, g.time_h.to_numpy(), g.signal.to_numpy())
                usable = np.isfinite(observed) & np.isfinite(expected) & (np.abs(expected) > 1e-9)
                if usable.any():
                    errors.append(
                        float(np.median(
                            np.abs(observed[usable] - expected[usable]) / np.abs(expected[usable])
                        ))
                    )
        return (float(np.median(errors)) if errors else np.inf), len(errors)

    ranked = []
    for permutation in itertools.permutations(blocks):
        assignment = dict(zip(constructs, permutation))
        error, n = fit(assignment)
        ranked.append((error, assignment, n))
    ranked.sort(key=lambda r: r[0])

    best_error, best_assignment, n_conditions = ranked[0]
    runner_up = ranked[1][0] if len(ranked) > 1 else float("inf")
    return LayoutRecovery(
        layout=PlateLayout(construct_columns=best_assignment, dose_rows=dose_rows),
        median_relative_error=best_error,
        runner_up_error=runner_up,
        n_conditions=n_conditions,
        ranking=[(e, {c: v[0] for c, v in a.items()}) for e, a, _ in ranked[:5]],
    )


# --- transcribed from the wet-lab logbook -------------------------------------


@dataclass(frozen=True)
class RecordedPlate:
    """A plate exactly as the logbook records it.

    The numerical recovery above agrees with these on construct columns and dose
    rows. Blank positions it could not have recovered at all -- they are not in the
    derived sheets -- and they moved between plates, so they are transcribed.

    Args:
        layout: Construct columns and dose rows.
        blank_wells: Wells filled with medium only.
        doses_mM: Dose ladder per stressor, ascending, aligned to ``dose_rows``.
        note: Anything about the plate worth carrying with the numbers.
    """

    layout: PlateLayout
    blank_wells: tuple[str, ...]
    doses_mM: dict[str, tuple[float, ...]]
    note: str = ""

    def __post_init__(self) -> None:
        overlap = set(self.culture_wells) & set(self.blank_wells)
        if overlap:
            raise ValueError(f"recorded wells cannot be both culture and blank: {sorted(overlap)}")
        if len(self.culture_wells) != len(set(self.culture_wells)):
            raise ValueError("recorded culture wells must have one construct and dose assignment")
        if len(self.blank_wells) != len(set(self.blank_wells)):
            raise ValueError("recorded blank wells must be unique")

    @property
    def culture_wells(self) -> tuple[str, ...]:
        return tuple(well for dose in range(len(self.layout.dose_rows))
                     for construct in self.layout.construct_columns
                     for well in self.layout.wells(construct, dose))

    @property
    def well_roles(self) -> dict[str, str]:
        return {**dict.fromkeys(self.culture_wells, "culture"),
                **dict.fromkeys(self.blank_wells, "blank")}


_STANDARD = PlateLayout(
    construct_columns={
        "UPRE1": (1, 2, 3),
        "UPRE2": (4, 5, 6),
        "NativeYap1": (7, 8, 9),
        "AlteredYap1": (10, 11, 12),
    },
    dose_rows=tuple("ABCDEFG"),
)
# Logbook maps 5 mM in row G but stocks top out at 4 mM; exports say 4.
_DOSES = {
    "DTT": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
    "H2O2": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0),
}

RECORDED_PLATES = {
    "20260722": RecordedPlate(
        _STANDARD, ("H1", "H2", "H3"), _DOSES,
        "RNA extracted from rows A, C and D (0, 0.2, 0.5 mM) after the read. "
        "Native Yap1 culture was at OD 0.9 against 3.7 for others, so it was made "
        "up in 10 mL rather than 3 mL -- a different dilution from every other well.",
    ),
    "20260728": RecordedPlate(_STANDARD, ("H4", "H5", "H6"), _DOSES),
    # Logbook says H4-H6; the export has none and its H1-H3 read as medium.
    "20260803": RecordedPlate(
        _STANDARD, ("H1", "H2", "H3"), _DOSES,
        "Logbook records H4-H6; corrected to H1-H3 against the export, where H1-H3 "
        "are flat in the reporter channel and no H4-H6 wells exist.",
    ),
    "20260804": RecordedPlate(
        _STANDARD, ("H1", "H2", "H3"), _DOSES,
        "Logbook records H4-H6, as 20260803's entry did and for the same reason: the "
        "template was copied from 20260728, where H4-H6 really are the blanks. The "
        "instrument file settles it. In 20260804's .xpt every well of row H is present "
        "and flat, and the reporter channel separates them by a factor of six -- H1-H3 "
        "read 341 RFU, the medium's own autofluorescence, and H4-H12 read 56, the "
        "instrument's dark count, so H4-H12 hold no medium. Gen5 agrees: the "
        "blank-subtracted .xlsx export is the archive minus the per-timepoint mean of "
        "H1:H3, reproduced to the last bit in fluorescence.",
    ),
    "20260821": RecordedPlate(_STANDARD, ("H1", "H2", "H3"), _DOSES),
}


def recorded_well_roles(*frames: pd.DataFrame | None) -> dict[str, str] | None:
    roles = None
    for frame in frames:
        if frame is None:
            continue
        if not isinstance(frame, pd.DataFrame):
            raise TypeError("recorded well identity requires a DataFrame with well labels")
        declared = frame.attrs.get("well_roles")
        if declared is None:
            continue
        if not isinstance(declared, Mapping):
            raise ValueError("well_roles must be a mapping of well labels to recorded roles")
        if roles is None:
            roles = {}
        for well, role in declared.items():
            if not isinstance(well, str) or not well:
                raise ValueError("recorded well labels must be nonempty strings")
            if not isinstance(role, str) or role not in ("culture", "blank"):
                raise ValueError(f"unknown recorded role {role!r} for well {well!r}")
            if well in roles and roles[well] != role:
                raise ValueError(f"conflicting recorded identity for well {well!r}")
            roles[well] = role
    return roles


def plate_key(filename: str) -> str | None:
    """Which logbook plate an export is a measurement of, or ``None`` if it is not one.

    The identity that matters to the analysis is the *plate*, not the file: a biological
    replicate is a culture read once, and two files can describe one read. ``20260804`` is
    both -- the instrument's own ``.xpt`` and the blank-subtracted ``.xlsx`` Gen5 exported
    from it -- and counting them separately would report four replicates where three were
    grown.

    This is the same lookup :func:`recorded_for_export` and :func:`blanks_for_export`
    already do, returning the key instead of the entry, so "the same plate" means one thing
    across the pipeline rather than being re-decided by each caller.
    """
    stem = pathlib.Path(str(filename).strip()).name
    for date in RECORDED_PLATES:
        if stem.startswith(date):
            return date
    return None


def recorded_for_export(filename: str) -> "RecordedPlate | None":
    """The logbook entry for the plate an export came from, by leading date.

    Two of the four biosensor exports carry no per-construct derived sheets, so the
    layout and the dose ladder cannot be read off the plate. Both are recorded here,
    which is enough to place every well without inventing anything -- and the recorded
    blanks for 20260728 (H4-H6) are exactly the wells detection finds in that export,
    so the entry is corroborated rather than merely asserted.
    """
    key = plate_key(filename)
    return None if key is None else RECORDED_PLATES[key]


def blanks_for_export(filename: str) -> tuple[str, ...] | None:
    """Blank wells recorded for the plate an export came from, by leading date.

    Returns ``None`` when the export is not one of the recorded plates, so a caller
    can fall back to detection rather than silently using someone else's blanks.
    """
    plate = recorded_for_export(filename)
    return None if plate is None else plate.blank_wells
