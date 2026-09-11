"""OD calibration from a dilution series: A_meas = A_true / (1 + k*A_true).

A dilution series is the only design in which true relative concentrations are known
independently of the reader. The deliverable is ``linear_range_max`` -- the measured OD where
deviation first exceeds a stated tolerance, which reduces exactly to tolerance/k.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from ..gates.g1_optical import OpticalQualityGate
from ..readings import CorrectedOD, RawOD, require

__all__ = ["ODCalibration", "fit_od_calibration"]

_MIN_DILUTION_SPREAD = 2.0  # highest/lowest dilution factor


def _saturating(true_od: np.ndarray, k: float) -> np.ndarray:
    return true_od / (1.0 + k * true_od)


@dataclass(frozen=True)
class ODCalibration:
    """Reader- and plate-specific absorbance response.

    Args:
        saturation_k: Curvature of the detector response, 1/OD. Zero is a perfectly
            linear reader.
        top_true_od: Fitted true OD of the undiluted stock.
        blank: Media-only absorbance that was subtracted before fitting.
        gdcw_per_od: Grams dry cell weight per linearised OD unit, if measured.
        residual_std: Standard deviation of the fit residuals, in OD units.
    """

    saturation_k: float
    top_true_od: float
    blank: float
    gdcw_per_od: float | None = None
    residual_std: float = 0.0

    def linearise(self, measured_above_blank: CorrectedOD) -> np.ndarray:
        """Blank-corrected measured OD -> true, dilution-linear OD.

        Raises:
            ValueError: if a reading sits at or beyond the detector's asymptote,
                where the inverse does not exist. Such a well carries no density
                information and must be excluded, not extrapolated.
        """
        a = require(measured_above_blank, CorrectedOD, name="measured_above_blank").array
        if self.saturation_k <= 0:
            return a
        asymptote = 1.0 / self.saturation_k
        if np.any(a >= asymptote * (1 - 1e-9)):
            raise ValueError(
                f"reading at or beyond the detector asymptote ({asymptote:.3f} OD); "
                "the inverse does not exist -- exclude the well rather than extrapolate"
            )
        return a / (1.0 - self.saturation_k * a)

    def linear_range_max(self, tolerance: float = 0.05) -> float:
        """Measured OD at which the response first departs from linearity by ``tolerance``.

        Args:
            tolerance: Fractional shortfall of measured against true OD to permit.
        """
        if not 0 < tolerance < 1:
            raise ValueError("tolerance must be in (0, 1)")
        if self.saturation_k <= 0:
            return float("inf")
        true_at_limit = tolerance / (self.saturation_k * (1.0 - tolerance))
        return float(_saturating(true_at_limit, self.saturation_k))

    def optical_gate(self, tolerance: float = 0.05, **overrides) -> OpticalQualityGate:
        """A G1 gate whose linear-range threshold comes from this calibration."""
        return OpticalQualityGate(od_linear_max=self.linear_range_max(tolerance), **overrides)

    def to_dcw(self, measured_above_blank: CorrectedOD) -> np.ndarray:
        """Blank-corrected measured OD -> dry cell weight, g/L.

        Raises:
            ValueError: if no dry-weight factor was measured. The conversion is
                strain- and medium-specific and is not guessed.
        """
        if self.gdcw_per_od is None:
            raise ValueError(
                "no gdcw_per_od was measured for this strain and medium; "
                "supply it to fit_od_calibration rather than assuming a literature value"
            )
        return self.linearise(measured_above_blank) * self.gdcw_per_od

    def summary(self) -> str:
        limit = self.linear_range_max()
        return (
            f"saturation k={self.saturation_k:.4f}/OD  "
            f"linear to {limit:.3f} measured OD (5% tolerance)  "
            f"stock true OD {self.top_true_od:.3f}  residual sd {self.residual_std:.4f}"
        )


def fit_od_calibration(
    dilution_factors: np.ndarray,
    measured_od: RawOD,
    blank: float,
    gdcw_per_od: float | None = None,
) -> ODCalibration:
    """Fit the detector response to a dilution series.

    Args:
        dilution_factors: Relative concentration of each well, stock = 1.0.
        measured_od: Raw absorbance of each well, blank included.
        blank: Media-only absorbance.
        gdcw_per_od: Optional measured dry-weight factor, g/L per linearised OD.
    """
    d = np.asarray(dilution_factors, dtype=float)
    a = require(measured_od, RawOD, name="measured_od").minus_blank(blank).array
    if d.shape != a.shape:
        raise ValueError("dilution factors and readings must have the same shape")
    if d.size < 3:
        raise ValueError("need at least 3 dilution steps")
    if not np.all(d > 0):
        raise ValueError("dilution factors must be positive")
    if d.max() / d.min() < _MIN_DILUTION_SPREAD:
        raise ValueError(
            f"dilution series spans only {d.max() / d.min():.2f}x; "
            f"need at least {_MIN_DILUTION_SPREAD}x for the curvature to be identifiable"
        )

    def model(dil, top, k):
        return _saturating(top * dil, k)

    top0 = float(a.max() / max(d.max(), 1e-9))
    fit, _ = curve_fit(
        model, d, a, p0=[top0, 0.1], bounds=([1e-9, 0.0], [np.inf, 100.0]), maxfev=20000
    )
    top, k = float(fit[0]), float(fit[1])
    residual = a - model(d, top, k)
    return ODCalibration(
        saturation_k=k,
        top_true_od=top,
        blank=float(blank),
        gdcw_per_od=gdcw_per_od,
        residual_std=float(np.std(residual)),
    )
