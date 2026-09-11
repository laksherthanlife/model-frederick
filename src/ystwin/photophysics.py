"""How the fluorophore itself responds to its chemical environment. A leaf, on purpose.

The reporter is a protein, and a protein's brightness is not a property of the promoter
driving it. This module holds the part of the signal that belongs to the fluorophore --
today one titration curve -- so that both ends of the experiment can reach it.

Why it is its own module and not a section of something larger
--------------------------------------------------------------
:func:`citrine_ph_response` lived in `generator/context.py`, next to the culture the cells
grow in. That put a *measurement-channel* constant inside the *generator*, so the only way
for `observation.py` -- the module that turns twin state into an instrument reading -- to
reach it was to import the generator. Both sit at the same depth, so that import would have
worked and it would have been wrong: the measurement model would then depend on the
simulator, and every consumer of a predicted RFU would drag `stress_panel`, `culture` and
`chemostat_physiology` in behind it.

So this file imports nothing from `ystwin`. It is depth 0. `generator/context.py` and
`observation.py` both import *it*, which is the shape that lets cytosolic pH reach the
instrument channel without a cycle and without either module owning the other's constants.

What the seam is worth, in the numbers this project measures
------------------------------------------------------------
A full unit of cytosolic acidification, pH 7.0 to 6.0, removes **30%** of Citrine's
fluorescence with the promoter untouched. The induction folds this project reports are
around **1.5**. So an unmodelled pH shift is a third of the signal, and it is not noise: the
stressors that acidify the cytosol are the ones being dosed, so it has dose structure and it
enters every dose-response fit as a slope. That is the whole reason this is a seam rather
than a footnote.
"""

from __future__ import annotations

__all__ = [
    "CITRINE_PKA",
    "citrine_ph_response",
]

CITRINE_PKA = 5.7
"""Chromophore pKa of Citrine, the fluorophore these constructs carry.

**Measured.** Griesbeck et al. 2001, PMID 11387331, "Reducing the environmental sensitivity of
yellow fluorescent protein. Mechanism and applications": the Q69M substitution "confers a much
lower pK(a) (5.7) than for previous YFPs", along with indifference to chloride.

pH is a property of the culture and this is a property of the protein, which is why the two
now live apart. Its largest effect on this experiment is not on growth but on the reporter --
see the module docstring for what that costs."""


def citrine_ph_response(ph: float) -> float:
    """Fraction of Citrine's maximum fluorescence at a given pH.

    Single-protonation Henderson-Hasselbalch, ``1 / (1 + 10 ** (pKa - pH))``, at the measured
    :data:`CITRINE_PKA`. The pKa is measured; the single-site form is the standard model for a
    fluorescent-protein titration and is **asserted** here -- Griesbeck reports a pKa, not a
    Hill coefficient, so the steepness of the curve is not pinned by that paper.

    Takes the **intracellular** pH, not the medium pH. Those are different numbers and yeast
    holds them apart: the medium can be at 4 while the cytosol is near 7. Passing a medium pH
    here would predict a dark cell that is not dark.
    """
    if not 0.0 < ph < 14.0:
        raise ValueError(f"pH must be between 0 and 14, got {ph}")
    return float(1.0 / (1.0 + 10.0 ** (CITRINE_PKA - float(ph))))
