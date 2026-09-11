"""Turning a glutathione probe's oxidation degree into a potential, with pH where it belongs.

Two published ranges exist for the same compartment -- around -290 to -306 mV assuming
pH 7.0, and -350 measured alongside a pHluorin that read 7.5. They are the same measurement.
The couple is GSSG + 2H+ + 2e- -> 2GSH, one proton per electron, so its potential moves
about 60 mV per pH unit and half a unit accounts for most of the gap.

The reason this matters is not tidiness. Roughly one pH unit of cytosolic acidification
under peroxide produces 60 mV of apparent shift, against a genuine response of 40 to 50, so
the artifact is larger than the signal and has the same sign. A glutathione sensor read
without a pH measurement cannot be interpreted under any stressor that acidifies.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "NERNST_MV_PER_PH",
    "ROGFP2_MIDPOINT_MV",
    "apparent_shift_from_ph",
    "e_gsh_from_oxd",
]

_GAS_CONSTANT = 8.314
_FARADAY = 96485.0
_TEMPERATURE_K = 303.15
_REFERENCE_PH = 7.0

NERNST_MV_PER_PH = 2.303 * _GAS_CONSTANT * _TEMPERATURE_K / _FARADAY * 1000.0

ROGFP2_MIDPOINT_MV = -280.0
"""roGFP2 half-oxidation potential at pH 7.0. Hanson 2004 puts it near -280 mV; the exact
figure shifts a few mV between reports and by 60 per pH unit, which is the larger effect."""


def e_gsh_from_oxd(oxidation_degree: float, ph: float = _REFERENCE_PH,
                   midpoint_mv: float = ROGFP2_MIDPOINT_MV) -> float:
    """Glutathione redox potential in mV, from the probe's oxidation degree and the pH.

    Both the probe midpoint and the glutathione couple carry two protons per two electrons,
    so both move together with pH and the correction is a single term.

    Args:
        oxidation_degree: Fraction of probe oxidised, strictly between 0 and 1.
        ph: Compartment pH. Assuming 7.0 where the cytosol is 7.5 mis-states the answer by
            30 mV, which is most of the gap between the two published literatures.
        midpoint_mv: Probe half-oxidation potential at pH 7.0.
    """
    if not 0.0 < oxidation_degree < 1.0:
        raise ValueError(
            f"oxidation degree must be between 0 and 1 exclusive, got {oxidation_degree}")
    nernst = _GAS_CONSTANT * _TEMPERATURE_K / (2 * _FARADAY) * 1000.0
    from_probe = midpoint_mv + nernst * np.log(oxidation_degree / (1.0 - oxidation_degree))
    return float(from_probe - NERNST_MV_PER_PH * (ph - _REFERENCE_PH))


def apparent_shift_from_ph(from_ph: float, to_ph: float) -> float:
    """Redox shift a pH change alone would appear to produce, in mV.

    Positive means it looks like an oxidation, which is what acidification does -- the same
    direction a real oxidant moves it, and the reason the two cannot be told apart without
    measuring pH.
    """
    return float(NERNST_MV_PER_PH * (from_ph - to_ph))
