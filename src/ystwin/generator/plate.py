"""S0: assemble whole synthetic plates in the shape the real exports take.

Defaults reproduce the 2026-07 and 2026-08 plates, including constructs inoculated at different
densities because on the real plate they were. Every record carries source="synthetic".
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..observation import ReporterOptics
from ..plate.layout import NEWPROTOCOL_LAYOUT, PlateLayout
from ..qpcr import STRESSOR_FOR_CONSTRUCT
from .culture import CultureParameters, simulate_culture

__all__ = ["DEFAULT_DOSES", "DEFAULT_PANEL", "PlateConditions", "SyntheticPlate", "generate_plate"]

# Defaults set so an uncalibrated run reproduces the real plates' OD range.
DEFAULT_PANEL = {
    "UPRE1": CultureParameters(
        mu_max=0.36, growth_ic50=3.0, growth_hill=1.2,
        promoter_basal=1.0e-3, promoter_peak=1.55e-3, promoter_ec50=1.1, promoter_hill=1.8,
        carrying_capacity=3.0,
    ),
    "UPRE2": CultureParameters(
        mu_max=0.40, growth_ic50=0.75, growth_hill=1.1,
        promoter_basal=1.0e-3, promoter_peak=1.60e-3, promoter_ec50=0.9, promoter_hill=2.0,
        carrying_capacity=3.0,
    ),
    "NativeYap1": CultureParameters(
        mu_max=0.35, growth_ic50=1.4, growth_hill=1.6,
        promoter_basal=1.0e-3, promoter_peak=1.45e-3, promoter_ec50=0.7, promoter_hill=2.2,
        carrying_capacity=2.6,
    ),
    "AlteredYap1": CultureParameters(
        mu_max=0.31, growth_ic50=1.3, growth_hill=1.6,
        promoter_basal=1.0e-3, promoter_peak=1.50e-3, promoter_ec50=0.75, promoter_hill=2.2,
        carrying_capacity=2.6,
    ),
}

DEFAULT_DOSES = {
    "DTT": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
    "H2O2": (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0),
}
DEFAULT_OPTICS = ReporterOptics(
    gain=2.6e6, background=400.0, autofluorescence=900.0, inner_filter_coeff=None
)


@dataclass(frozen=True)
class PlateConditions:
    """How the plate was set up and read.

    Args:
        layout: Construct columns and dose rows.
        blank_wells: Wells filled with medium only.
        duration_h: Length of the read.
        n_timepoints: Readings taken.
        doses: Dose ladder per stressor.
        optics: Reporter measurement parameters.
        inoculum_od: Starting blank-corrected OD, per construct. A bare float
            applies to all; the real plates differed between constructs by 3.7x.
        od_blank: Medium-only absorbance.
        reader_cv: Multiplicative reader noise.
        well_cv: Well-to-well variation in starting biomass.
        gdcw_per_od: Dry weight per OD unit, for the biomass-to-OD conversion.
        seed: Random seed.
    """

    layout: PlateLayout = NEWPROTOCOL_LAYOUT
    blank_wells: tuple[str, ...] = ("H1", "H2", "H3")
    duration_h: float = 4.14
    n_timepoints: int = 25
    doses: dict[str, tuple[float, ...]] = field(default_factory=lambda: dict(DEFAULT_DOSES))
    optics: ReporterOptics = DEFAULT_OPTICS
    inoculum_od: dict[str, float] | float = 0.16
    od_blank: float = 0.098
    reader_cv: float = 0.01
    well_cv: float = 0.04
    gdcw_per_od: float = 0.42
    seed: int | None = None

    def inoculum_for(self, construct: str) -> float:
        if isinstance(self.inoculum_od, dict):
            return float(self.inoculum_od.get(construct, 0.16))
        return float(self.inoculum_od)


@dataclass(frozen=True)
class SyntheticPlate:
    """A generated plate plus the ground truth that produced it."""

    od: pd.DataFrame
    rfu: pd.DataFrame
    truth: pd.DataFrame
    blank_wells: tuple[str, ...]
    conditions: PlateConditions = field(repr=False)

    @property
    def times_h(self) -> np.ndarray:
        return self.od.index.to_numpy(dtype=float)


def generate_plate(
    panel: dict[str, CultureParameters],
    conditions: PlateConditions = PlateConditions(),
) -> SyntheticPlate:
    """Generate one plate from the mechanistic model.

    Args:
        panel: Construct -> kinetics.
        conditions: Plate setup and reader settings.
    """
    rng = np.random.default_rng(conditions.seed)
    t = np.linspace(0.0, conditions.duration_h, conditions.n_timepoints)
    optics = conditions.optics

    od_cols: dict[str, np.ndarray] = {}
    rfu_cols: dict[str, np.ndarray] = {}
    records = []

    for construct, params in panel.items():
        if construct not in conditions.layout.construct_columns:
            continue
        stressor = STRESSOR_FOR_CONSTRUCT.get(construct, "DTT")
        ladder = conditions.doses.get(stressor, DEFAULT_DOSES["DTT"])
        for dose_index, dose in enumerate(ladder):
            try:
                wells = conditions.layout.wells(construct, dose_index)
            except IndexError:
                break
            for well in wells:
                start_od = conditions.inoculum_for(construct) * float(
                    np.exp(rng.normal(0.0, conditions.well_cv))
                )
                biomass0 = start_od * conditions.gdcw_per_od
                out = simulate_culture(t, dose, params, biomass0)

                true_od = out["biomass"] / conditions.gdcw_per_od + conditions.od_blank
                signal = (optics.autofluorescence + optics.gain * out["reporter"]) * out["biomass"]
                true_rfu = optics.background + signal

                od_cols[well] = true_od * (1 + rng.normal(0, conditions.reader_cv, t.size))
                rfu_cols[well] = true_rfu * (1 + rng.normal(0, conditions.reader_cv, t.size))
                records.append({
                    "well": well, "construct": construct, "stressor": stressor,
                    "dose_mM": float(dose), "dose_index": dose_index,
                    "promoter_activity": float(out["promoter_activity"][0]),
                    # Unrestricted rate for this strain and dose, so it is seed-independent.
                    "growth_rate_max": float(params.growth_rate_at(dose)),
                    "initial_biomass": float(biomass0),
                    "source": "synthetic",
                })

    for well in conditions.blank_wells:
        od_cols[well] = conditions.od_blank * (1 + rng.normal(0, conditions.reader_cv, t.size))
        rfu_cols[well] = optics.background * (1 + rng.normal(0, conditions.reader_cv, t.size))

    index = pd.Index(t, name="time_h")
    order = sorted(od_cols, key=lambda w: (w[0], int(w[1:])))
    return SyntheticPlate(
        od=pd.DataFrame({w: od_cols[w] for w in order}, index=index),
        rfu=pd.DataFrame({w: rfu_cols[w] for w in order}, index=index),
        truth=pd.DataFrame(records),
        blank_wells=tuple(conditions.blank_wells),
        conditions=conditions,
    )
