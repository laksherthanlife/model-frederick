"""Figures, kept out of the analysis so that a plot can never invent a number.

Nothing in this package computes a result. Every figure is handed a frame that came
either from a committed table in ``outputs/`` or from a generator whose output is
deterministic, and its only job is to render what it was given. That separation is the
reason the figures can be trusted: a rendering step that also fitted, interpolated or
defaulted would be a second analysis with no tests on it.

The one thing the figures *do* decide is what they refuse to draw. A fold change with no
estimable interval gets a point and an explicit note that no interval exists, never an
error bar; a dose whose recovered activity went negative gets a marked gap, never a
plotted zero. See :mod:`ystwin.viz.figures`.
"""

from __future__ import annotations

from .figures import (
    ARTEFACT_ZONES,
    CORRECTABLE_GDNA_CEILING,
    PALETTE,
    attribution_power,
    attribution_power_table,
    decoupling_designs,
    decoupling_grid_scatter,
    dose_response_correction,
    dose_response_table,
    g1_pass_rates,
    g1_summary,
    g4_contamination,
    g4_contamination_table,
    read_table,
    save_figure,
    series_colours,
)

__all__ = [
    "ARTEFACT_ZONES",
    "CORRECTABLE_GDNA_CEILING",
    "PALETTE",
    "attribution_power",
    "attribution_power_table",
    "decoupling_designs",
    "decoupling_grid_scatter",
    "dose_response_correction",
    "dose_response_table",
    "g1_pass_rates",
    "g1_summary",
    "g4_contamination",
    "g4_contamination_table",
    "read_table",
    "save_figure",
    "series_colours",
]
