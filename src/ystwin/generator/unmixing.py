"""Linear unmixing of co-expressed reporters, and what it does to the growth confound.

Co-expression changes the problem. Reporters in one cell share mu exactly, so a ratio
between two of them cancels the dilution term outright:

    R_i / R_j -> k_i / k_j    when their loss rates match

That is an exact cancellation, not a correction, and it needs no growth estimate. What
breaks it is mismatched maturation or stability between fluorophores, plus whatever error
the unmixing itself introduces.

Emission is modelled as a skewed Gaussian, which is the shape that makes adjacent
fluorophores overlap and unmixing necessary in the first place.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["PANEL", "FluorophoreSpec", "mixing_matrix", "ratio_series", "unmix"]


@dataclass(frozen=True)
class FluorophoreSpec:
    """One fluorophore's optics and reporter kinetics."""

    name: str
    excitation_nm: float
    emission_nm: float
    maturation_h: float = 0.4
    k_deg: float = 0.0
    emission_width_nm: float = 22.0
    brightness: float = 1.0

    def emission_at(self, wavelength: np.ndarray) -> np.ndarray:
        """Skewed-Gaussian emission: a narrow blue edge and a long red tail."""
        w = np.asarray(wavelength, dtype=float)
        offset = w - self.emission_nm
        width = np.where(offset < 0, self.emission_width_nm * 0.75, self.emission_width_nm * 1.6)
        return self.brightness * np.exp(-0.5 * (offset / width) ** 2)


# Maturation half-times are order-of-magnitude literature figures, not measured here.
PANEL = {
    "CFP": FluorophoreSpec("CFP", 433, 475, maturation_h=0.4, brightness=0.8),
    "eGFP": FluorophoreSpec("eGFP", 488, 507, maturation_h=0.4, brightness=1.0),
    "YFP": FluorophoreSpec("YFP", 514, 527, maturation_h=0.3, brightness=1.1),
    "mOrange2": FluorophoreSpec("mOrange2", 549, 565, maturation_h=2.3, brightness=0.7),
    "mApple": FluorophoreSpec("mApple", 568, 592, maturation_h=0.6, brightness=0.6),
}


def mixing_matrix(names, bands, panel=None) -> np.ndarray:
    """Fraction of each fluorophore's emission falling in each detection band.

    Args:
        names: Fluorophores present, in the order the unmixed output should take.
        bands: ``(low_nm, high_nm)`` detection windows.
        panel: Spec lookup; defaults to :data:`PANEL`.
    """
    panel = panel or PANEL
    if len(bands) < len(names):
        raise ValueError(
            f"unmixing {len(names)} fluorophores needs at least as many bands, got {len(bands)}"
        )
    matrix = np.zeros((len(bands), len(names)))
    for row, (low, high) in enumerate(bands):
        grid = np.linspace(low, high, 64)
        for col, name in enumerate(names):
            spec = panel[name] if isinstance(name, str) else name
            matrix[row, col] = float(np.trapezoid(spec.emission_at(grid), grid) / (high - low))
    return matrix


def unmix(measured: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Least-squares recovery of fluorophore abundances from banded readings."""
    measured = np.asarray(measured, dtype=float)
    solution, *_ = np.linalg.lstsq(np.asarray(matrix, dtype=float), measured, rcond=None)
    return solution


def _steady_state(activity: float, growth: float, spec: FluorophoreSpec) -> float:
    """Mature reporter at steady state, including the maturation step."""
    loss = max(growth + spec.k_deg, 1e-12)
    maturation = 1.0 / max(spec.maturation_h, 1e-12)
    immature = activity / (loss + maturation)
    return maturation * immature / loss


def ratio_series(
    stress_activity: float,
    reference_activity: float,
    growth_rate: float,
    stress_spec: FluorophoreSpec,
    reference_spec: FluorophoreSpec,
) -> float:
    """Steady-state ratio of two co-expressed reporters at one growth rate.

    Growth cancels exactly when the two share loss and maturation rates. It does not when
    they differ, which is why a reference reporter should be the same fluorophore class as
    the sensor wherever possible.
    """
    return _steady_state(stress_activity, growth_rate, stress_spec) / _steady_state(
        reference_activity, growth_rate, reference_spec
    )
