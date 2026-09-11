"""Recombinant burden and proteome allocation, built as a DECLARED GENOTYPE AXIS.

**THIS BLOCK DOES NOT CLAIM CRITERION 1, AND THE REASON IS STRUCTURAL RATHER THAN
PROVISIONAL.** Criterion 1 is an unbroken chain from an environmental input to a
bioproduction observable. This block's only input is ``n_copies``, which is genotype.
Nothing in L0-L5 writes into it. The regulator that would make allocation
environment-responsive is Dot6/Tod6 -- the RP/RiBi repression arm of TORC1 -- and TORC1 is
REJECT in `REVISED_BUILD_LIST.md` §2 on three measurements, so that arm cannot be built
behind this one. `ARCHITECTURE_GAPS.md` 0.4 offered exactly two ways to close the loop and
one of them is now gone; this module takes the other one explicitly. The loop
burden -> stress -> adaptation -> burden is **open at both ends** here, and this docstring
is where that is said out loud rather than left implied by an absent import. The other end
is worth stating precisely, because the obvious version of it is wrong: the reverse arm,
stress multiplying the burden slope, HAS been measured, and against heat, which is stressor
6 of this panel's 25. What refuses it is the observable Farkas measured it on, not its
absence -- see :data:`BURDEN_STRESS_MULTIPLIER`, which names the one experiment that would
switch this block's environmental input on.

WHAT IS BUILT, AND WHAT IS ONLY REPORTED.

* **The growth arm is IN_CHAIN.** ``mu = mu_env * (1 - loss_per_copy * n)``. It is
  parameter-free: :data:`KAFRI_GROWTH_LOSS_PER_COPY` is one MEASURED number and nothing
  else enters. Its ablation is scored against this repository's own measured growth floor
  on 183 real wells -- see :func:`score_growth_ablation`. The verdict is CONDITIONAL on
  the declared genotype and the condition is the point of the re-scoping: the effect is
  EXACTLY ZERO at the panel strains' own single copy, first clears the 0.0117 /h floor at
  5 copies, and is 1.84x it at 8. Parameter-free is not condition-free: 0.01 is Kafri's
  RICH-MEDIUM slope, the same paper measured that the extent varies with medium without
  publishing a per-condition value, and that residual is :data:`BURDEN_STRESS_MULTIPLIER`.
  It is also worth saying what the 183 wells do and do not do in that ablation. They fix
  mu_bar and nothing else: the copy column is constant across them, so the two fitted models
  are observationally IDENTICAL in sample -- same residual RMS to machine precision, same
  predictions -- and the effect is the closed form mu_bar*loss*(query - reference)/(1 - loss*
  reference). That is the honest reading of a declared axis: the wells set the scale the
  floor is compared against, the axis supplies the increment, and no well discriminates the
  two models. :func:`score_growth_ablation` says so and a test pins it.
* **The allocation arm (phi_H per copy, b, r0, lag) is REPORTED.** It carries 8 free
  scalars against 0 countable targets and :data:`ALLOCATION_PARAMS` fails
  ``require_gate()`` on purpose. That is worse than the pH slice the gate was calibrated
  on (`ARCHITECTURE_GAPS.md` 0.3: 8 against 2). The 8 is an upper bound and the verdict
  does not need it: four of the eight are REFUSED with no value at all, and ``b`` and
  ``c_r`` are functions of two of the others, so the irreducible count is **two** swept
  axes -- ``phi_per_copy`` and ``phi_max_library`` -- against still zero targets. Nothing
  here feeds a prediction; the numbers exist so the next reader can see what would have to
  be measured.

THE TWO CORRECTIONS THIS MODULE CARRIES, both found by verification against the sources.

1. ``METZL_RAZ_RIBOSOME_SLOPE`` is 0.35 per GENERATION/h, not per mu. Against mu in 1/h it
   is 0.35/ln(2) = 0.504943. It is IMPORTED from `fba/stress_proteostasis.py`, which
   already carries the corrected value and its regression test, and is never retyped here.
2. ``b`` is SWEPT, not fixed. b = 0.72 is an arithmetic category error -- it divides
   Kafri's 18% growth loss AT 18 COPIES by Metzl-Raz's 25% proteome fraction AT THE
   HIGHEST-BURDEN STRAIN, which is ~13 copies at 1.9%/copy, an 18/13.2 = 1.37x inflation.
   :func:`refuted_burden_slopes` records that arithmetic so it cannot come back.

ONE MISATTRIBUTION WAS FOUND AND FIXED WHILE SOURCING THIS BLOCK, and it is resolved
against the primary record rather than against the brief. The harvest and
`REVISED_BUILD_LIST.md` attribute the >40% non-toxic ceiling to 'Kintaka/Moriya 2025'. The
paper is eLife 13:RP99572 (PMID 40960085, PMC12443478) and its author list, read off the
NCBI E-utilities esummary and efetch records for that PMID, is Fujita Y, Namba S, Kamada Y
and Moriya H -- Kintaka R is not on it. The quoted sentence is verbatim from the abstract
those records return, and the separate Kintaka 2016 citation (PMID 27538565) is correct and
untouched. :func:`misattributed_quotes` carries it beside what the record really says.
A SECOND MISATTRIBUTION was found in adversarial verification and is this module's own:
:data:`PHI_PER_COPY` asserted a MASS basis that neither Metzl-Raz nor Kafri states, and both
of its calibrations are PaxDb-anchored and molar-proportional. It is recorded, not rescaled.
The growth arm is untouched by it -- ``b * phi_per_copy`` is Kafri's growth-axis measurement
and never passes through a proteome fraction -- and the exposure is the REPORTED arm's
comparison against Eguchi's MASS ceiling, up to 0.54x at a 27 kDa product.
Farkas 2018's citation needed no correction and is reproduced as the harvest carries it:
eLife 7:e29845 (PMID 29377792). `tests/test_mech_burden.py` fails if either drifts.

THE THING THAT FELL OUT OF DOING IT PROPERLY, and it is why the growth arm passes its gate.
``b`` was DERIVED as (Kafri's 1% growth loss per copy)/(a proteome fraction per copy), and
the band's two ends are the same measured growth loss over two calibrations of the same
fraction: 1.9% by mass spectrometry (Metzl-Raz) and 2.3% by fluorescence (Kafri Fig 1E).
So ``b * phi_per_copy == 0.01 per copy`` identically across the band, the growth channel is
**invariant to where in the sweep you sit**, and the sweep is a reparameterisation rather
than a degree of freedom. b = 0.72 is dangerous precisely because it breaks that identity:
0.72 x 1.9% = 1.37% per copy, a 37% inflation of a number Kafri measured directly.

E8 IS RE-SCORED AS A FITTED CHECK, and its arithmetic does not reproduce. The architecture
bills "8 copies -> +14% lag predicted against +15% measured" as held out. It is not held
out: ``c_r`` is derived from the r0 endpoints of the same burden library whose lag was
measured in the same figure of the same paper. And +14% is not obtainable anywhere in
the two swept bands -- the exact algebra runs +17.4% (c_r = 0.079, phi = 1.9%/copy) to
+30.9% (c_r = 0.104, phi = 2.3%/copy), and the paper's own "~15% fall in r0 at 8
copies" gives 1/(1-0.15) = +17.6%. +15% is the first-order identity 1/(1-x) ~= 1+x
applied to that same 15%. :func:`metzl_raz_lag_check` prints every route side by side.

Time base is HOURS throughout, matching `generator/culture.py` and `mech/state.py`.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .. import paths
from ..fba.stress_proteostasis import (
    EGUCHI_BURDEN_LIMIT,
    METZL_RAZ_RIBOSOME_SLOPE,
    ribosome_allocation_penalty,
)
from ..pathway import proteome
from .ablation import (
    GROWTH_RATE_FLOOR,
    AblationResult,
    Fit,
    FittableModel,
    Observable,
    ablate,
)
from .integrate import Trajectory, integrate_window, reduce_for_window
from .params import Param, ParamRegistry, provenance_table
from .state import PLATE_READ_4H, Encoding, MechState, StateVar, Window

__all__ = [
    "ALLOCATION_PARAMS",
    "BURDEN_STATES",
    "BurdenAboveCeiling",
    "BurdenedGrowth",
    "EGUCHI_CEILING",
    "GROWTH_PARAMS",
    "KAFRI_GROWTH_LOSS_PER_COPY",
    "FUJITA_CEILING",
    "LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES",
    "METZL_RAZ_COPIES_AT_8",
    "PRODUCT_FRACTION",
    "ProductClass",
    "UnburdenedGrowth",
    "allocation_gate",
    "burden_slope_band",
    "burdened_growth_rate",
    "consistent_phi_per_copy",
    "copies_to_clear_the_floor",
    "destabilisation_invariance",
    "displacement_over_measured",
    "gem_burden_comparison",
    "growth_gate",
    "growth_loss",
    "growth_rate_observable",
    "integrate_product_fraction",
    "lag_ratio",
    "metzl_raz_lag_check",
    "misattributed_quotes",
    "paxdb_route_refused",
    "plate_growth_rates",
    "product_fraction_rhs",
    "proteome_fraction",
    "provenance",
    "reduction_audit",
    "refuted_burden_slopes",
    "require_admissible",
    "ribosome_reserve",
    "score_growth_ablation",
    "steady_state_fraction",
]


# --------------------------------------------------------------------------------------
# The growth arm: one measured number, and the refusal that rides with it
# --------------------------------------------------------------------------------------

GROWTH_PARAMS = ParamRegistry("mech/burden.py::growth (IN_CHAIN)")
"""The arm that enters a prediction. One MEASURED constant and one REFUSED multiplier."""

KAFRI_GROWTH_LOSS_PER_COPY = GROWTH_PARAMS.add(Param.measured(
    "growth_loss_per_copy", 0.01, "relative growth-rate loss per genomic copy",
    "Kafri, Metzl-Raz, Jona & Barkai 2016, Cell Reports 14:22-31 (PMID 26725116): ~1% "
    "growth-rate loss per genomic copy measured by competition (Fig 1G/H, S1D), on a "
    "pTDH3-mCherry library of 1 to ~20 copies verified by qPCR of genomic DNA. Cross-checks "
    "against the same paper's Fig 1F, where the 18-copy strain grows ~18% slower than the "
    "1-copy strain in YPD: 0.18/18 = 0.01. Both routes are in RICH MEDIUM -- Fig 1H's own "
    "legend says 'Cells were grown in YPD'. THE SHAPE generalises and THE SLOPE does not, "
    "and the paper says both in one sentence whose second half is the half that matters "
    "here: 'In all conditions tested, growth rate decreased linearly with increasing "
    "mCherry expression, although to different extents' (Figs 2A, 2C, S2), "
    "with 'cells growing in low-phosphate conditions appeared less sensitive to the "
    "introduced burden compared to cells growing in standard media, possibly reflecting "
    "their somewhat lower division rate'. THAT LAST CLAUSE IS THE ONE THAT BITES: it "
    "hypothesises the cost scales with division rate, which is exactly the invariance a "
    "constant RELATIVE loss assumes away. So no threshold is "
    "installed anywhere -- linearity holds in every condition -- but 0.01 is the rich-medium "
    "slope carried unchanged into every environment. NO PER-CONDITION SLOPE IS STATED IN "
    "THE TEXT; Fig 4D plots the 'absolute effect per copy number' for low phosphate, SC and "
    "low nitrogen, so the three numbers exist in a figure this module cannot read off. That "
    "residual is m_stress, and it is REFUSED rather than fitted",
    ci95=None))
"""0.01 per copy, MEASURED IN RICH MEDIUM. The ONLY number the growth prediction uses, and it
is on the growth axis directly rather than through a proteome fraction -- which is why the
b/phi decomposition cancels out of it. See :func:`consistent_phi_per_copy`. What it does NOT
carry is a condition, and Kafri measured that the extent varies with medium; the environment
term that would is :data:`BURDEN_STRESS_MULTIPLIER`, refused."""

BURDEN_STRESS_MULTIPLIER = GROWTH_PARAMS.add(Param.refused(
    "m_stress", "dimensionless multiplier on the burden slope",
    "Farkas, Kalapis, Bodi et al. 2018, eLife 7:e29845 (PMID 29377792), "
    "'Hsp70-associated chaperones have a critical role in buffering protein production "
    "costs', reports FOUR multipliers on the cost of yEVenus, not the two usually quoted: "
    "cycloheximide over 0.0018-0.18 ug/mL reaching 3.7-fold (Fig 2A, p<0.001); MPA "
    "2.0-fold (Fig 2B, p<0.001; the paper's own dose is INTERNALLY INCONSISTENT -- "
    "'0.30 ug/ml' in Materials and methods against '30 ug/ml' in the Fig 2B legend, a 100x "
    "spread, so neither is usable); AZC over 1-2.5 mM reaching >4-fold (Fig 3B, p<0.001); "
    "and HEAT over 30/37/40 degC reaching 2.8-fold at 40 (Fig 3A, p<0.001, >=12 technical "
    "measurements of 15 biological replicates per point)",
    reason="the observable, not the absence of a dose-response and not the absence of the "
           "stressor. Heat IS on this panel (stress_panel.STRESSORS['heat'], 'degC above "
           "30', ec50 6.0, lethal 18.0), so Farkas's +7 and +10 degC sit inside its range, "
           "and cycloheximide and AZC are real dose ladders. What every one of the four "
           "measures is a ratio of COLONY SIZE AFTER 48 h ON SOLID MEDIUM -- Farkas's "
           "Materials and methods, verbatim -- which is the same observable "
           "FARKAS_COLONY_SIZE_SLOPE is excluded for, and it convolves lag, growth rate and "
           "final density. Multiplying Kafri's liquid-competition SPECIFIC GROWTH RATE "
           "slope by a solid-medium colony-size ratio is the observable substitution this "
           "module refuses everywhere else. AND THE REFUSAL IS NOT THAT m = 1: Kafri 2016 "
           "measured the same multiplier ON THE RIGHT OBSERVABLE and found it is not, "
           "reporting that the linear burden holds 'in all conditions tested ... although to "
           "different extents', with low-phosphate cells 'less sensitive to the introduced "
           "burden compared to cells growing in standard media'. That is a MEDIUM x burden "
           "interaction in liquid on relative growth rate, and it is qualitative -- the paper "
           "states no per-condition slope in its text -- Fig 4D plots the per-copy effect for "
           "its three media, so the numbers exist only in a figure. So the block carries a "
           "rich-medium slope into every environment and m != 1 is measured, unread and "
           "refused.",
    missing="Farkas's own 30/37/40 degC ladder rerun with specific growth rate in liquid at "
            "each point, on a burden series. That is one experiment, on a stressor this "
            "panel already applies, and it is the single measurement that would turn this "
            "block's environmental input on. Until it exists m = 1 for all 25 stressors, "
            "and that is a refusal rather than a measurement."))
"""The burden x stress interaction. This is the loop this repository cannot close, and it is
the second reason criterion 1 is not claimed. Note what the refusal is NOT: it is not that
nobody measured the interaction, and it is not that the stressor is off the panel. Farkas
measured it against heat, which is stressor 6 of 25 here. It is that he measured it on a
different observable, and this module does not convert between observables. And it is NOT a
claim that the multiplier is 1: Kafri's own condition series measured that it is not, on the
right observable, without publishing a slope for any condition. So this refusal is the one
place where the growth arm's single free scalar actually sits, and it is a known non-unity."""

GROWTH_PARAMS.add_target(GROWTH_RATE_FLOOR)


def growth_loss(n_copies: float) -> float:
    """Relative growth-rate loss for ``n_copies`` of the cassette. Parameter-free.

    ``0.01 * n``, straight off Kafri's competition measurement. Linear with no threshold
    because Kafri measured linearity with no threshold anywhere inside the range; the
    ceiling is enforced separately by :func:`require_admissible`, which REFUSES rather than
    bending the line.

    Raises:
        ValueError: on a negative copy number, or on a load that would drive growth to zero
            or below -- at 100 copies this law predicts mu = 0, which is outside every
            measurement behind it.
    """
    n = float(n_copies)
    if n < 0:
        raise ValueError(f"copy number must be non-negative, got {n}")
    loss = float(KAFRI_GROWTH_LOSS_PER_COPY) * n
    if loss >= 1.0:
        raise ValueError(
            f"{n:g} copies drives the predicted growth loss to {loss:.3g}, at or past total "
            f"arrest. Kafri 2016's library runs 1 to ~20 copies; this law is a measured "
            f"slope inside that range and not an extrapolation to lethality")
    return loss


def burdened_growth_rate(mu_env: float, n_copies: float,
                         reference_copies: float = 1.0) -> float:
    """The environment's growth rate, taxed by the DECLARED genotype.

    Args:
        mu_env: Growth rate the environment sets, 1/h. Measured, not asserted -- on this
            project's plates it is the log-OD slope; in the vessel it is fitted off the
            biomass series.
        n_copies: Cassette copies. **Declared genotype. Nothing measures it here.**
        reference_copies: The genotype ``mu_env`` was measured on. The burden that strain
            already carries is absorbed into ``mu_env`` and is unidentifiable, so only the
            INCREMENT is predicted. Default 1, which is the panel strains' declared
            genotype, not a measurement of them.

    Raises:
        ValueError: if ``mu_env`` is not positive.
    """
    if mu_env <= 0:
        raise ValueError(f"mu_env must be positive, got {mu_env}")
    return float(mu_env) * (1.0 - growth_loss(n_copies)) / (1.0 - growth_loss(reference_copies))


# --------------------------------------------------------------------------------------
# The allocation arm: REPORTED, and it fails its own gate
# --------------------------------------------------------------------------------------

ALLOCATION_PARAMS = ParamRegistry("mech/burden.py::allocation (REPORTED, gate FAILS)")
"""The proteome-fraction arm. Kept visible rather than deleted; feeds nothing."""

PHI_PER_COPY = ALLOCATION_PARAMS.add(Param.swept(
    "phi_per_copy", "g product protein / g total cell protein, per genomic copy",
    "Two calibrations of the same quantity, and the sources agree they agree. Metzl-Raz "
    "et al. 2017, eLife 6:e28034 (PMID 28857745), Materials and Methods: one genomic "
    "pTDH3-mCherry copy is ~1.9% of the proteome by mass spectrometry, and the paper adds "
    "that this 'is in very good agreement with our previous fluorescence-based calibration "
    "(Kafri et al., 2016b)'. Kafri 2016 (PMID 26725116) Fig 1E: a single pTDH3-GFP "
    "construct is ~2% of the proteome with the dashed reference line at 2.3%. BASIS NOT "
    "STATED BY EITHER PAPER: Metzl-Raz's 'Protein abundance was obtained by summing the "
    "three most intense, unique peptides per protein', a top-3 intensity metric that scales "
    "with molar amount and not with mass; and Kafri Fig 1E's abundance axis is ppm from "
    "PaxDb/CYCLOPS (Table S3), which PHI_FROM_PAXDB itself calls MOLAR. See "
    "misattributed_quotes()",
    bounds=(0.019, 0.023),
    missing="SDS-PAGE densitometry against a BSA dilution series on THIS project's own "
            "construct, which is how Farkas 2018 got 3.7% for yEVenus. One gel -- and it is "
            "the same gel that would settle the mass-versus-molar basis above."))
"""1.9-2.3% of the proteome per copy. Compared against Eguchi's explicitly MASS ceiling and
fed to a MASS conversion, on a basis neither source states -- see :func:`paxdb_route_refused`
and the second row of :func:`misattributed_quotes`."""

BURDEN_SLOPE = ALLOCATION_PARAMS.add(Param.swept(
    "b", "relative growth loss per unit heterologous proteome mass fraction",
    "DERIVED then SWEPT: b = (Kafri's 1% growth loss per copy) / (phi_per_copy), evaluated "
    "at both ends of that band. 0.01/0.023 = 0.4348 and 0.01/0.019 = 0.5263, which "
    "REVISED_BUILD_LIST.md writes as 0.43-0.53. The endpoints are COMPUTED from the two "
    "rows above rather than typed, so b*phi_per_copy = 0.01/copy holds exactly at both. "
    "REFUSED: 0.72 (arithmetic category error) and 0.68 (a different observable) -- see "
    "refuted_burden_slopes()",
    bounds=(float(KAFRI_GROWTH_LOSS_PER_COPY) / 0.023,
            float(KAFRI_GROWTH_LOSS_PER_COPY) / 0.019),
    missing="nothing measurable would collapse this band, because it is not a measurement "
            "band: it is one measured growth loss expressed against two calibrations of "
            "one proteome fraction. Collapsing phi_per_copy collapses this."))
"""``float()`` raises. That is deliberate: a midpoint taken here is a fitted parameter
nobody counted, and this particular midpoint has already been mis-taken twice (0.53 as a
point, 0.72 as an inflation)."""

R0_UNBURDENED = ALLOCATION_PARAMS.add(Param.measured(
    "r0_0", 0.081, "fraction of proteome, non-translating ribosome reserve",
    "Metzl-Raz et al. 2017, eLife 6:e28034 (PMID 28857745): r0 = 8.1% in the wild type, "
    "verbatim, the intercept of phi_R = 0.35*(generations/h) + r0. The abstract's '~8%' is "
    "the same quantity rounded and is what fba/stress_proteostasis.py carries as "
    "METZL_RAZ_RIBOSOME_RESERVE; the two differ by 1.2% and are cross-checked in "
    "tests/test_mech_burden.py so they cannot drift apart"))

R0_HIGHEST_BURDEN = ALLOCATION_PARAMS.add(Param.measured(
    "r0_min", 0.055, "fraction of proteome",
    "Metzl-Raz et al. 2017 (PMID 28857745): r0 falls 8.1% -> 5.5% at the highest burden, "
    "verbatim. The other endpoint of the only measured r0 excursion in the literature"))

PHI_MAX_LIBRARY = ALLOCATION_PARAMS.add(Param.swept(
    "phi_max_library", "fraction of proteome at the highest-burden strain",
    "The paper's own three statements are not jointly consistent and this band is the "
    "inconsistency, not a measurement error. Metzl-Raz 2017 states (i) mCherry is '~25% of "
    "the total proteome in the highest burden cells', (ii) one copy is 1.9%, and (iii) "
    "8-copy cells show a '~15% decrease in the measured r0'. (i) with (ii) puts the top "
    "strain at ~13.2 copies, which forces a 19.5% r0 fall at 8 copies and contradicts "
    "(iii); "
    "(iii) with (ii) forces a top strain of ~17-18 copies, i.e. phi ~ 0.33. Two of the "
    "three can hold. Carried as the band 0.25-0.33",
    bounds=(0.25, 0.33),
    missing="the copy number of Metzl-Raz's highest-burden strain, stated as a number. "
            "That one integer settles which of the paper's three statements is loose."))

RIBOSOME_RESERVE_SLOPE = ALLOCATION_PARAMS.add(Param.swept(
    "c_r", "fall in r0 per unit heterologous proteome mass fraction",
    "DERIVED then SWEPT: c_r = (r0_0 - r0_min)/phi_max_library = 0.026/[0.33, 0.25] = "
    "[0.0788, 0.1040]. The architecture's 0.079 and 0.104 are these two ends. Computed "
    "from the rows above rather than typed",
    bounds=((0.081 - 0.055) / 0.33, (0.081 - 0.055) / 0.25),
    missing="same integer as phi_max_library. Nothing else moves it."))

RIBOSOME_SLOPE = ALLOCATION_PARAMS.add(Param.derived(
    "metzl_raz_ribosome_slope", METZL_RAZ_RIBOSOME_SLOPE, "per (1/h)",
    "IMPORTED from fba/stress_proteostasis.py::METZL_RAZ_RIBOSOME_SLOPE, never retyped. "
    "Metzl-Raz 2017 (PMID 28857745) Fig 2A measures phi_R = 0.35*X + 0.08 with X in "
    "GENERATIONS per hour; the paper converts it itself ('Delta_r/Delta_mu = 21/ln(2) "
    "[min]'), so against mu in 1/h the slope is 0.35/ln(2) = 0.504943"))

PROTEIN_MOLECULES_PER_CELL = ALLOCATION_PARAMS.add(Param.measured(
    "protein_molecules_per_cell", 5.0e7, "molecules/cell",
    "Ghaemmaghami et al. 2003, Nature 425:737-41 (PMID 14562106), quoted verbatim by "
    "Geiler-Samerotte et al. 2011 (PNAS 108:680-5, PMID 21187411): 'the estimated total "
    "protein content of a haploid yeast cell is ~50 "
    "million molecules per cell'. Carried so a mass fraction can be turned into a "
    "synthesis flux in molecules/s without inventing the conversion"))

FARKAS_COLONY_SIZE_SLOPE = ALLOCATION_PARAMS.add(Param.measured(
    "b_on_colony_size_48h", 0.025 / 0.037, "relative colony-size loss per proteome fraction",
    "Farkas et al. 2018, eLife 7:e29845 (PMID 29377792): yEVenus at 3.7% of the total proteome "
    "by SDS-PAGE densitometry against a 100-800 ng BSA series, costing 2.5% from a "
    "high-copy 2-micron plasmid and nothing detectable from a single copy. 0.025/0.037 = "
    "0.676. NAMED FOR ITS OBSERVABLE ON PURPOSE: colony size after 48 h on solid medium "
    "convolves lag, growth rate and final density, so it is an UPPER BOUND on b and not a "
    "replicate of Kafri's specific growth rate. It is deliberately outside the swept band"))

PHI_FROM_PAXDB = ALLOCATION_PARAMS.add(Param.refused(
    "phi_per_copy_from_paxdb", "g product / g total protein per copy",
    "pathway/proteome.py vendors PaxDb's integrated S. cerevisiae table and "
    "PROMOTER_ANCHOR_PPM carries TDH3 at 23256 molar ppm",
    reason="PaxDb's ppm counts MOLECULES and proteome.proteome_fraction() says so in its "
           "own docstring, while Eguchi's ceiling is a MASS fraction and says so explicitly "
           "('not \"15% molecules of total protein molecules\"'). Metzl-Raz's phi_H is NOT "
           "the mass counterexample this refusal used to claim -- it is PaxDb-anchored too, "
           "which is the second row of misattributed_quotes(). The conversion is the "
           "product's MW over "
           "AVERAGE_PROTEIN_MW_G_PER_MOL, which is 0.84x to 1.49x on this repository's own "
           "three carotenoid enzymes, and proteome.py's docstring records this exact "
           "molar-for-mass substitution having been made and fixed once already.",
    missing="the product's molar mass. With it the route is arithmetic; without it the "
            "route is a factor of up to 1.5 in either direction. See paxdb_route_refused()."))

TAI_GATE = ALLOCATION_PARAMS.add(Param.refused(
    "tai_burden_gate", "tRNA adaptation index",
    "The 0.37 threshold reaches this project through a narrative review that attributes it "
    "to Eguchi 2018 and Xu 2021; Eguchi's own finding is qualitative -- poorly "
    "codon-optimised isozymes never reached burden levels even from pTDH3 on a multicopy "
    "plasmid -- and 0.37 is the S. cerevisiae genome AVERAGE, not a measured threshold",
    reason="no tAI is computed anywhere in this package, so the gate has no input; and the "
           "number offered for it is a genome statistic quoted secondarily rather than a "
           "threshold anyone measured.",
    missing="a tAI for this project's cassettes and a measured burden-versus-tAI curve. "
            "Until both exist, applying the burden law is a declaration about the "
            "construct, not a derivation from it."))

B_SECRETED = ALLOCATION_PARAMS.add(Param.refused(
    "b_secreted", "relative growth loss per proteome fraction, secreted product",
    "Both measurements of the burden slope -- Kafri 2016 and Farkas 2018 -- used soluble "
    "CYTOSOLIC fluorescent proteins. Kintaka et al. 2016 (PMID 27538565) showed protein "
    "LOCALISATION processes have their own, far lower limits, and Eguchi's hexose "
    "transporters never became detectable at all",
    reason="there is no published burden slope for a secreted or membrane-bound product in "
           "S. cerevisiae, and the two that exist are for a compartment this one is not in.",
    missing="a copy-number or expression series on a secreted product with growth measured "
            "at each level. This is the constant fba/insulin.py would need."))

LAG_ASSAY_FLOOR = ALLOCATION_PARAMS.add(Param.refused(
    "lag_time_noise_floor", "relative CV of time-to-50%-OD-rise",
    "This repository has measured exactly two noise floors, both in "
    "generator/panel_experiment.py: OBSERVED_ACTIVITY_CV = 0.146 for dilution-corrected "
    "reporter activity and MEASURED_GROWTH_RATE_SE = 0.0117 /h for growth rate",
    reason="lag time after dilution from stationary is a third observable and neither floor "
           "is its noise. The plate reads behind both floors start in growth and carry no "
           "lag phase to replicate, so nothing here bounds it. Borrowing 0.146 for it is "
           "the substitution REVISED_BUILD_LIST.md records as the error that sank Route C2.",
    missing="replicate wells re-inoculated from stationary phase, with the spread of "
            "time-to-50%-OD-rise across them. That measurement also turns E8 from a fitted "
            "check into a scoreable ablation."))

METZL_RAZ_COPIES_AT_8 = 8
"""The copy number Metzl-Raz's recovery-time statement is about. An integer from the paper."""

LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES = 0.15
"""Metzl-Raz 2017 (PMID 28857745), verbatim: 'the recovery time of eight-copy burden cells
was prolonged by ~15%, consistent with the corresponding ~15% decrease in the measured r0'.
RECOVERY TIME, not lag from a fresh inoculum -- naming the observable matters here, because
it is the observable :data:`LAG_ASSAY_FLOOR` says nothing on this project's plates measures.
MEASURED, and the target E8 is scored against -- as a FITTED check, because c_r comes off the
same library in the same figure."""


def growth_gate():
    """Criterion (e) for the arm that enters a prediction."""
    return GROWTH_PARAMS.gate()


def allocation_gate():
    """Criterion (e) for the reported arm. It FAILS, and that is the finding."""
    return ALLOCATION_PARAMS.gate()


def consistent_phi_per_copy(b: float) -> float:
    """The proteome fraction per copy that a given ``b`` was computed against.

    ``b`` is not an independent constant: it is Kafri's measured 1%-per-copy growth loss
    divided by a calibration of the proteome fraction per copy. So a value of ``b`` names
    its own calibration, and the product ``b * phi`` is the measurement both share. This is
    the function that makes the invariance checkable rather than asserted.

    Raises:
        ValueError: if ``b`` is outside the swept band. Outside it the identity is exactly
            what breaks -- b = 0.72 implies 1.39% per copy, and Kafri measured 1%.
    """
    pinned = BURDEN_SLOPE.at(b)
    return float(KAFRI_GROWTH_LOSS_PER_COPY) / pinned


def burden_slope_band() -> tuple[float, float]:
    """``(low, high)`` for b, computed rather than typed. Rounds to the mandated 0.43-0.53."""
    return BURDEN_SLOPE.bounds


def proteome_fraction(n_copies: float, phi_per_copy: float) -> float:
    """``f_H``: the fraction of total protein SYNTHESIS the product takes.

    Args:
        n_copies: Declared genotype.
        phi_per_copy: A point pinned inside :data:`PHI_PER_COPY`'s band, on the record.

    The name is deliberate. Kafri 2016 Fig 3B measured that a ~20-fold DESTABILISED GFP
    costs the same fitness as a stable one at the same copy number, so the cost is on the
    synthesis process and not on the accumulated pool. ``f_H`` is a flux fraction; the pool
    is :data:`PRODUCT_FRACTION`, and it is not what growth is taxed on.
    """
    if n_copies < 0:
        raise ValueError(f"copy number must be non-negative, got {n_copies}")
    return float(n_copies) * PHI_PER_COPY.at(phi_per_copy)


def ribosome_reserve(f_H: float, c_r: float) -> float:
    """``r0 = r0_0 - c_r * f_H``. REPORTED: nothing downstream reads it.

    Raises:
        ValueError: if the load drives the reserve to zero or below, which is outside the
            8.1% -> 5.5% excursion the two endpoints were measured across.
    """
    r0 = float(R0_UNBURDENED) - RIBOSOME_RESERVE_SLOPE.at(c_r) * float(f_H)
    if r0 <= 0:
        raise ValueError(
            f"f_H = {f_H:g} drives the ribosome reserve to {r0:.4g}. Metzl-Raz measured r0 "
            f"from 8.1% down to 5.5%; a reserve at or below zero is an extrapolation past "
            f"every point behind the slope")
    return r0


def lag_ratio(f_H: float, c_r: float) -> float:
    """``t_lag / t_lag,0 = r0_0 / r0``. REPORTED, and its check is FITTED -- see the module
    docstring and :func:`metzl_raz_lag_check`."""
    return float(R0_UNBURDENED) / ribosome_reserve(f_H, c_r)


def metzl_raz_lag_check(n_copies: int = METZL_RAZ_COPIES_AT_8) -> pd.DataFrame:
    """E8, re-scored as a FITTED check, with every route to the number side by side.

    `ARCHITECTURE_GAPS.md` 5.5 makes two charges and both are reproduced here by running
    the arithmetic: the claimed +14% is not obtainable from the architecture's own swept
    values, and the check is not held out because ``c_r`` is derived from the r0 endpoints
    of the same burden library whose lag was measured in the same figure of the same paper.

    Returns:
        One row per route, with ``fractional_increase`` and what each route assumes.
    """
    rows = []
    for phi_max in PHI_MAX_LIBRARY.bounds:
        c_r = (float(R0_UNBURDENED) - float(R0_HIGHEST_BURDEN)) / phi_max
        for phi in PHI_PER_COPY.bounds:
            f_H = n_copies * phi
            rows.append({
                "route": f"algebra, c_r={c_r:.4f} (phi_max={phi_max:g}), phi={phi:g}/copy",
                "fractional_increase": lag_ratio(f_H, c_r) - 1.0,
                "held_out": False,
                "note": "c_r and the lag come from the same library and the same figure",
            })
    rows.append({
        "route": "direct from the paper's own '~15% decrease in the measured r0'",
        "fractional_increase": 1.0 / (1.0 - LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES) - 1.0,
        "held_out": False,
        "note": "no c_r at all; r0_0/r0 with the measured fall substituted directly",
    })
    rows.append({
        "route": "first-order identity 1/(1-x) ~= 1+x applied to that same 15%",
        "fractional_increase": LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES,
        "held_out": False,
        "note": "this is where the paper's 'consistent with' comes from, and it is the "
                "measured number restated rather than a prediction",
    })
    rows.append({
        "route": "MEASURED (Metzl-Raz 2017: eight-copy recovery time prolonged by ~15%)",
        "fractional_increase": LAG_MEASURED_FRACTIONAL_INCREASE_AT_8_COPIES,
        "held_out": True,
        "note": "the target. Its assay has no measured noise floor here -- see "
                "LAG_ASSAY_FLOOR -- so it cannot be scored under criterion (a)",
    })
    frame = pd.DataFrame(rows)
    frame.insert(0, "n_copies", n_copies)
    frame["architecture_claimed"] = 0.14
    frame["reproduces_the_claim"] = np.isclose(
        frame["fractional_increase"], 0.14, atol=5e-3)
    return frame


def refuted_burden_slopes() -> pd.DataFrame:
    """The two values of ``b`` this module refuses, with the arithmetic that refutes each.

    Kept as code rather than a comment because rule 3 of this build -- a number attributed
    to a paper that does not contain it -- has burned this project repeatedly, and 0.72 is
    exactly that shape: it is arithmetic on two papers' numbers that were never about the
    same strain.
    """
    phi_at_25_percent_copies = 0.25 / 0.019
    inflation = 18.0 / phi_at_25_percent_copies
    consistent = float(KAFRI_GROWTH_LOSS_PER_COPY) / 0.019
    return pd.DataFrame([
        {
            "value": 0.72,
            "verdict": "REFUTED -- arithmetic category error",
            "arithmetic": (
                f"0.18 (Kafri's growth loss AT 18 COPIES) / 0.25 (Metzl-Raz's proteome "
                f"fraction AT THE HIGHEST-BURDEN STRAIN) = 0.72. At 1.9%/copy, 25% is "
                f"{phi_at_25_percent_copies:.1f} copies, not 18. The inflation is "
                f"18/{phi_at_25_percent_copies:.1f} = {inflation:.3f}, and "
                f"{consistent:.4f} x {inflation:.3f} = {consistent * inflation:.3f}, which "
                f"reproduces 0.72 to rounding"),
            "in_the_band": False,
        },
        {
            "value": 0.68,
            "verdict": "EXCLUDED -- a different observable, so an upper bound",
            "arithmetic": (
                f"Farkas 2018: 0.025 / 0.037 = {float(FARKAS_COLONY_SIZE_SLOPE):.4f}, on "
                f"COLONY SIZE AFTER 48 h ON SOLID MEDIUM, which convolves lag, growth rate "
                f"and final density. Kafri's b is on specific growth rate in liquid "
                f"competition. Two incommensurable measurements agreeing is a coincidence, "
                f"not a replication"),
            "in_the_band": False,
        },
        {
            "value": consistent,
            "verdict": "IN THE BAND, at its 1.9%/copy end",
            "arithmetic": (
                f"{float(KAFRI_GROWTH_LOSS_PER_COPY):g} (measured per copy) / 0.019 "
                f"(Metzl-Raz mass spectrometry) = {consistent:.4f}. The other end is "
                f"{float(KAFRI_GROWTH_LOSS_PER_COPY):g}/0.023 = "
                f"{float(KAFRI_GROWTH_LOSS_PER_COPY) / 0.023:.4f}, Kafri's own fluorescence "
                f"calibration"),
            "in_the_band": True,
        },
    ])


def misattributed_quotes() -> pd.DataFrame:
    """Attributions this block inherited from the harvest that the primary record contradicts.

    Rule 3 of this build calls a number attributed to a paper that does not contain it the
    worst outcome available, so the one that was found is carried as data rather than as
    prose, with the record that settles it named.

    The second row is this module's own, found in adversarial verification rather than
    inherited: the mass-versus-molar basis of ``phi_per_copy`` was asserted here and is
    stated by neither source. It is recorded rather than acted on, because rescaling a band
    on an inference is the same act as asserting it was.

    The ROW THAT IS NOT HERE is the more useful half. Metzl-Raz's '~25% of the total
    proteome in the highest burden cells' was challenged during this build as possibly being
    the paper's INACTIVE-RIBOSOME ~25% (r0 ~8% over a ribosomal fraction r ~30%, which the
    paper also states, four times). It is not: the sentence 'a library of strains expressing
    increasing amounts of mCherry proteins was generated, with mCherry levels reaching ~25%
    of the total proteome in the highest burden cells' is verbatim in the Results. Both
    quantities are ~25% and they are different quantities; :data:`PHI_MAX_LIBRARY` uses the
    first, which is the correct one, and the challenge is recorded here so the next reader
    does not have to re-run it.
    """
    return pd.DataFrame([
        {
            "claim": "the >40% non-toxic ceiling is 'Kintaka/Moriya 2025'",
            "asserted_by": "harvest__recombinant-protein-burden...json; "
                           "REVISED_BUILD_LIST.md",
            "source_really_says": (
                "eLife 13:RP99572 (PMID 40960085) is by Fujita Y, Namba S, Kamada Y and "
                "Moriya H. Kintaka R is not an author. The finding is unaffected: 'These "
                "proteins can be expressed at levels exceeding 40% of total protein while "
                "maintaining yeast growth' is verbatim from its abstract"),
            "checked_against": "NCBI E-utilities esummary and efetch for PMID 40960085",
            "consequence": "citation corrected; the number stands (FUJITA_CEILING)",
        },
        {
            "claim": "phi_per_copy is a MASS fraction, not a molecule fraction",
            "asserted_by": "mech/burden.py::PHI_PER_COPY, this module, until it was checked",
            "source_really_says": (
                "neither paper states a basis and both are molar-proportional. Metzl-Raz "
                "2017 Methods: 'Protein abundance was obtained by summing the three most "
                "intense, unique peptides per protein', normalised to 'the median of 13 "
                "rich conditions from PaxBD'; Kafri 2016 Table S3, the source of Fig 1E's "
                "abundance axis: 'Normalized Protein Abundance Compiled from 16 Different "
                "Datasets Found in Two Databases in Parts per Million ... Data from PaxDB "
                "and CYCLOPS'. PaxDb ppm counts MOLECULES, which is the ground "
                "PHI_FROM_PAXDB is refused on"),
            "checked_against": "EuropePMC full text for PMID 28857745 and PMID 26725116",
            "consequence": (
                "NOT rescaled -- the band stays 0.019-0.023 and the growth arm is untouched, "
                "because b*phi_per_copy = 0.01/copy is on the growth axis and never passes "
                "through a proteome fraction. The exposure is on the REPORTED arm only: f_H "
                "is compared against Eguchi's explicitly MASS ceiling and fed to "
                "grams_per_gdcw_at_protein_fraction, and at Farkas's own 27 kDa fluorescent "
                "protein over a 50 kDa proteome average the factor is 0.54x. 8 copies would "
                "be 8.2% rather than 15.2% and require_admissible() would not raise"),
        },
    ])


def paxdb_route_refused() -> pd.DataFrame:
    """Why phi_H cannot be read off PaxDb today, with the size of the error.

    The architecture calls phi_H "estimable today with no wet-lab work" from
    ``n_copies x PaxDb(the promoter's own native gene)``. PaxDb is MOLAR and the ceiling this
    module scores against -- Eguchi's -- is by MASS, so the route needs the product's molar
    mass. This tabulates the conversion factor on the three heterologous enzymes this
    repository does carry masses for, which is the honest bound on how wrong the
    substitution is. Note what this refusal does NOT rest on any more:
    :data:`PHI_PER_COPY` is not the mass counterexample, and the second row of
    :func:`misattributed_quotes` says why.
    """
    from ..fba.stress_proteostasis import CRT_ENZYME_MASSES_KDA

    average_kda = proteome.AVERAGE_PROTEIN_MW_G_PER_MOL / 1000.0
    rows = [
        {
            "protein": name,
            "mw_kda": mw,
            "mass_over_molar_fraction": mw / average_kda,
            "source": "UniProt sequence mass, fba/stress_proteostasis.py::"
                      "CRT_ENZYME_MASSES_KDA",
        }
        for name, mw in sorted(CRT_ENZYME_MASSES_KDA.items())
    ]
    frame = pd.DataFrame(rows)
    frame["average_protein_kda"] = average_kda
    frame["refusal"] = PHI_FROM_PAXDB.reason
    return frame


# --------------------------------------------------------------------------------------
# Admissibility: refuse above the ceiling, never a threshold inside the range
# --------------------------------------------------------------------------------------

class BurdenAboveCeiling(ValueError):
    """The declared load is past the highest expression anyone has measured for its class.

    A refusal and not a saturation. Kafri 2016 measured linearity with no threshold anywhere
    inside the range, so bending the law near the ceiling would be inventing a shape the
    data contradicts; the honest move at the edge is to stop.
    """


class ProductClass:
    """Which measured ceiling applies. The gate is folding quality, not amount.

    There is no '% of proteome at which the stress response fires', in either direction:
    Fujita/Moriya 2025 (PMID 40960085) reached >40% with no proteotoxic response at all, and
    Geiler-Samerotte 2011 (PMID 21187411) paid 3.2% of growth rate for a misfolding variant
    held below 0.1% of total protein. NOT 'the full Hsf1 regulon', which is what this
    docstring said until it was checked: the paper reports 25 differentially regulated
    proteins of which 20 have Hsf1p-bound promoters, and its own abstract calls this
    'coordinated synthesis of interacting cytosolic chaperone proteins in the absence of a
    wider stress response'. So the class is a declared property of the product, not a dose.
    """

    NON_TOXIC = "non-toxic-cysteine-free"
    """Fujita, Namba, Kamada & Moriya 2025 (PMID 40960085): mox-YG and Gpm1-CCmut 'can be
    expressed at levels exceeding 40% of total protein while maintaining yeast growth'."""

    EGFP_LIKE = "egfp-like-or-glycolytic"
    """Eguchi et al. 2018, eLife 7:e34595 (PMID 30095406): 'the limits of some glycolytic
    proteins were up to 15% of the total cellular protein'. EGFP itself triggers the
    heat-shock response through cysteine-driven aggregates (Namba 2022, G3 12:jkac106,
    PMID 35485947), so it sits here."""

    ALL = (NON_TOXIC, EGFP_LIKE)


FUJITA_CEILING = ALLOCATION_PARAMS.add(Param.measured(
    "phi_ceiling_non_toxic", 0.40, "fraction of total protein",
    "Fujita, Namba, Kamada & Moriya 2025, eLife 13:RP99572 (PMID 40960085), abstract "
    "verbatim: mox-YG and Gpm1-CCmut 'can be expressed at levels "
    "exceeding 40% of total protein while maintaining yeast growth'. The burden state there "
    "is amino-acid/nitrogen depletion with ribosome biogenesis down and respiration up -- "
    "and NO proteotoxic response"))

EGUCHI_CEILING = ALLOCATION_PARAMS.add(Param.measured(
    "phi_ceiling_egfp_like", EGUCHI_BURDEN_LIMIT, "fraction of total protein",
    "Eguchi et al. 2018, eLife 7:e34595 (PMID 30095406): 'the limits of some glycolytic "
    "proteins were up to 15% of the total cellular protein'. IMPORTED from "
    "fba/stress_proteostasis.py::EGUCHI_BURDEN_LIMIT rather than retyped"))

_CEILINGS = {ProductClass.NON_TOXIC: FUJITA_CEILING, ProductClass.EGFP_LIKE: EGUCHI_CEILING}


def require_admissible(f_H: float, product_class: str = ProductClass.EGFP_LIKE) -> float:
    """Refuse a load past its class's measured ceiling. Returns ``f_H`` unchanged.

    Raises:
        BurdenAboveCeiling: above the ceiling.
        ValueError: on an unknown class or a negative load.
    """
    if product_class not in _CEILINGS:
        raise ValueError(
            f"product class {product_class!r} is not one of {ProductClass.ALL}. The ceiling "
            f"is a property of folding quality and cysteine content, not of amount, so it "
            f"has to be declared per product")
    if f_H < 0:
        raise ValueError(f"f_H must be non-negative, got {f_H}")
    ceiling = _CEILINGS[product_class]
    if f_H > float(ceiling):
        raise BurdenAboveCeiling(
            f"f_H = {f_H:.4g} is above the measured ceiling {float(ceiling):g} for "
            f"{product_class!r}. {ceiling.source}. Refused rather than saturated: Kafri 2016 "
            f"measured linearity with no threshold anywhere inside the range, so a Hill term "
            f"here would be a shape the data contradicts")
    return float(f_H)


# --------------------------------------------------------------------------------------
# The one state, and the window that decides whether it is one
# --------------------------------------------------------------------------------------

PRODUCT_FRACTION = StateVar(
    name="product_fraction",
    units="g product protein / g total cell protein",
    tau_growth_multiple=1.0,
    encoding=Encoding.LEVEL,
    source="DEFINITION plus one measured input. dphi/dt = f_H*mu - (mu + k_dH)*phi relaxes "
           "to f_H*mu/(mu+k_dH) with tau = 1/(mu+k_dH); at k_dH = 0, which is the stable "
           "product Kafri 2016 and Metzl-Raz 2017 both used, that is 1/mu. Audited at "
           "k_dH = 0 because that is the SLOWEST case: any product turnover only shortens "
           "tau and pushes the verdict further toward elimination",
    sweep_axis=True,
    note="Declared genotype axis. n_copies and the product's turnover are inputs nothing in "
         "L0-L5 writes into, which is why this block does not claim criterion 1. LEVEL and "
         "not INTEGRATING: it is a fraction relaxing toward f_H, not a pool accumulating.",
)
"""The product pool, as a fraction of total protein. **Growth is NOT taxed on it** -- that
is Kafri Fig 3B, and :func:`destabilisation_invariance` runs the test."""

BURDEN_STATES = MechState((PRODUCT_FRACTION,))
"""One state. The growth arm is algebra over a genotype constant and has none."""


def steady_state_fraction(f_H: float, mu: float, k_dH: float = 0.0) -> float:
    """``f_H * mu / (mu + k_dH)``. The algebraic value the state relaxes to."""
    if mu <= 0:
        raise ValueError(f"mu must be positive, got {mu}")
    if k_dH < 0:
        raise ValueError(f"k_dH must be non-negative, got {k_dH}")
    return float(f_H) * float(mu) / (float(mu) + float(k_dH))


def product_fraction_rhs(f_H: float, mu: float, k_dH: float = 0.0):
    """The right-hand side for :func:`integrate_product_fraction`, in 1/h.

    Args:
        f_H: Fraction of protein SYNTHESIS the product takes. Constant over a window,
            because its only input is genotype.
        mu: Specific growth rate, 1/h. Measured off the biomass series.
        k_dH: Product turnover, 1/h. A declared product property, 0 for a stable one.
    """
    if mu <= 0:
        raise ValueError(f"mu must be positive, got {mu}")
    if k_dH < 0:
        raise ValueError(f"k_dH must be non-negative, got {k_dH}")

    def rhs(t: float, y: Sequence[float]):
        return [f_H * mu - (mu + k_dH) * y[0]]

    return rhs


def integrate_product_fraction(window: Window, f_H: float, mu: float, *,
                               k_dH: float = 0.0, phi0: float = 0.0,
                               n_points: int = 201) -> Trajectory:
    """Integrate the pool across a declared window with the package's pinned stiff driver."""
    return integrate_window(
        product_fraction_rhs(f_H, mu, k_dH), {"product_fraction": float(phi0)},
        window, BURDEN_STATES, n_points=n_points)


def reduction_audit(window: Window):
    """Criterion (d) for this block's one state, against a declared window.

    Returns a :class:`~ystwin.mech.integrate.ReducedSystem` whose ``audit_table()`` carries
    the ratio, the exact time-integral bias and the verdict. The two verdicts differ between
    this project's vessels and that is the point of declaring the window.
    """
    return reduce_for_window(BURDEN_STATES, window)


def destabilisation_invariance(n_copies: float, mu: float, phi_per_copy: float,
                               k_dH_multiples: Sequence[float] = (0.0, 20.0)) -> pd.DataFrame:
    """E7: growth is invariant to product stability; the pool is not. Zero free parameters.

    Kafri 2016 Fig 3B measured that a ~20-fold less stable GFP 'had the same effect on
    fitness' at the same copy number. A model that taxes the accumulated pool gets that
    wrong by ``mu/(mu+k_dH)``; this one gets it right, because the tax is on the synthesis
    fraction. Both columns are computed here so the claim is checkable rather than asserted.

    IT IS NOT HELD OUT, and saying so costs nothing while implying otherwise would be the
    error this module spends a section on for E8. The synthesis-versus-pool structure was
    CHOSEN from Fig 3B -- :func:`proteome_fraction`'s own docstring says so -- and a check
    that a structure reproduces the measurement it was chosen from is a consistency check,
    not a prediction. What it does buy is that the choice is now legible in code and that the
    rejected alternative's error is a number rather than an assertion: 21x at Kafri's own
    20-fold variant, and ``1 + k_dH/mu`` in general.

    The counterfactual uses the ``b`` PAIRED with the pinned ``phi_per_copy``, not an endpoint
    of ``b``'s band. ``b`` and ``phi_per_copy`` are two spellings of one measurement and only
    the paired corners satisfy ``b*phi = 0.01``; taking an unpaired corner made the stable row
    disagree with ``growth_loss`` by up to 17% for a reason that has nothing to do with
    stability, which is the comparison this function exists to make.

    Args:
        n_copies: Declared genotype.
        mu: Growth rate, 1/h.
        phi_per_copy: A point pinned inside :data:`PHI_PER_COPY`.
        k_dH_multiples: Turnover rates as multiples of ``mu``. 20 is Kafri's variant.
    """
    f_H = proteome_fraction(n_copies, phi_per_copy)
    b = float(KAFRI_GROWTH_LOSS_PER_COPY) / PHI_PER_COPY.at(phi_per_copy)
    rows = []
    for multiple in k_dH_multiples:
        k_dH = float(multiple) * mu
        rows.append({
            "k_dH_over_mu": float(multiple),
            "steady_state_pool": steady_state_fraction(f_H, mu, k_dH),
            "growth_loss": growth_loss(n_copies),
            "pool_taxing_model_would_predict": b * steady_state_fraction(f_H, mu, k_dH),
        })
    frame = pd.DataFrame(rows)
    frame.insert(0, "f_H", f_H)
    frame.insert(0, "n_copies", n_copies)
    return frame


# --------------------------------------------------------------------------------------
# Criterion (a): the growth ablation, on real wells, against the measured growth floor
# --------------------------------------------------------------------------------------

_GATE1_TABLES = ("NewProtocol", "Replicate3")
_PLATE_READ_H = PLATE_READ_4H.duration_h
"""4.14 h, imported rather than retyped: `mech/state.py` already carries it and its
source (generator/plate.py::PlateConfig.duration_h, the real Synergy read length)."""


def plate_growth_rates() -> pd.DataFrame:
    """Per-well specific growth rates off the three committed gate-1 tables.

    ``mu = ln(od_fold_change) / 4.14 h`` on every well that passed gate 1, which is the same
    chord estimate `mech/state.py::PLATE_GROWTH_BAND_PER_H` is derived from and the same
    instrument the 0.0117 /h floor was measured on. It is a chord over the whole read rather
    than the slope of a fitted window, so it is the MEAN specific rate across the read.

    Returns:
        ``plate``, ``well``, ``growth_rate`` (1/h) and ``copies``, which is filled in by
        :func:`growth_rate_observable` from the declared genotype.

    Raises:
        FileNotFoundError: if the three tables are not in this checkout.
    """
    tables = sorted((paths.REPO_ROOT / "outputs").glob("g1_2026*.csv"))
    wanted = [t for t in tables if any(k in t.name for k in _GATE1_TABLES)]
    if len(wanted) != 3:
        raise FileNotFoundError(
            f"the three NewProtocol gate-1 tables are not in this checkout; found "
            f"{[t.name for t in wanted]}. They are what the growth floor and the growth "
            f"band were both measured on, and this ablation is not scoreable without them")
    frames = []
    for table in wanted:
        rows = pd.read_csv(table).query("passed")
        frames.append(pd.DataFrame({
            "plate": table.name,
            "well": rows["well"].astype(str),
            "growth_rate": np.log(rows["od_fold_change"].astype(float)) / _PLATE_READ_H,
        }))
    return pd.concat(frames, ignore_index=True)


@dataclass
class _LinearGrowth(FittableModel):
    """Shared fitting machinery: one scalar, solved in closed form rather than searched."""

    name: str = "growth"
    parameter_names: tuple[str, ...] = ("mu_env",)
    response_column: str = "growth_rate"

    def _weights(self, data: pd.DataFrame) -> np.ndarray:
        raise NotImplementedError

    def refit(self, data: pd.DataFrame) -> Fit:
        y = np.asarray(data[self.response_column], dtype=float)
        w = self._weights(data)
        mu_env = float(np.dot(w, y) / np.dot(w, w))
        residual = w * mu_env - y
        return Fit(
            model=self.name,
            parameters={"mu_env": mu_env},
            n_free=1,
            n_rows=len(data),
            residual_rms=float(np.sqrt(np.mean(residual ** 2))),
            converged=True,
        )

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        return float(parameters["mu_env"]) * self._weights(data)


@dataclass
class BurdenedGrowth(_LinearGrowth):
    """``mu = mu_env * (1 - 0.01 * copies)``. The full model, and parameter-free in the
    burden term: the only scalar it fits is the environment's own growth rate."""

    name: str = "burdened (Kafri 1%/copy)"

    def _weights(self, data: pd.DataFrame) -> np.ndarray:
        return np.array([1.0 - growth_loss(n) for n in data["copies"]], dtype=float)


@dataclass
class UnburdenedGrowth(_LinearGrowth):
    """``mu = mu_env``. The reduced model: the burden term deleted, the same one free
    scalar, so the comparison cannot be won by spending a parameter."""

    name: str = "unburdened"

    def _weights(self, data: pd.DataFrame) -> np.ndarray:
        return np.ones(len(data), dtype=float)


def growth_rate_observable(query_copies: float, *,
                           reference_copies: float = 1.0,
                           data: pd.DataFrame | None = None) -> Observable:
    """Growth rate at a DECLARED copy number, fitted on real wells.

    Both models are fitted to the same measured per-well growth rates, at the panel strains'
    declared genotype; the summary then evaluates each model at ``query_copies``. That is
    what ``Observable.summarise`` taking the model as well as the fit is for -- the copy
    number axis does not exist in this project's data and cannot be fitted, so the axis is
    declared and the models are evaluated along it.

    The effect this produces is
    ``mu_bar * loss * (query - reference) / (1 - loss * reference)``, so it depends on the
    INCREMENT in copy number and only weakly on the reference. At ``query == reference`` it
    is exactly zero, which is correct: the burden a fitted strain already carries is
    absorbed into ``mu_env`` and is not identifiable from these wells.

    Args:
        query_copies: The genotype to score. Declared.
        reference_copies: The genotype the fitted wells carry. Declared, default 1.
        data: Rows to fit. Defaults to :func:`plate_growth_rates`.
    """
    frame = (plate_growth_rates() if data is None else data.copy())
    frame["copies"] = float(reference_copies)
    query = pd.DataFrame({"copies": [float(query_copies)], "growth_rate": [np.nan]})

    def summarise(model: FittableModel, fit: Fit) -> float:
        return float(model.predict(query, fit.parameters)[0])

    return Observable(
        name=GROWTH_RATE_FLOOR.name,
        units="1/h",
        assay=GROWTH_RATE_FLOOR.assay,
        scored="absolute",
        data=frame,
        summarise=summarise,
        description=(
            f"specific growth rate predicted at {float(query_copies):g} declared cassette "
            f"copies, from wells carrying {float(reference_copies):g}"),
    )


def score_growth_ablation(query_copies: float, *, reference_copies: float = 1.0,
                          data: pd.DataFrame | None = None) -> AblationResult:
    """Criterion (a) for the growth arm: delete the burden term and see what moves.

    Scored against :data:`~ystwin.mech.ablation.GROWTH_RATE_FLOOR` -- 0.0117 /h, the
    standard error of a log-OD slope across 210 real wells -- and ABSOLUTELY, because that
    floor is an absolute rate error and turning it into a CV at the point of use is the
    substitution `mech/ablation.py` names.

    The answer is conditional and the condition is the whole re-scoping: at the panel's own
    single-copy genotype the effect is exactly zero, and it first clears the floor at 5
    copies. Below that this block is a correctly-sourced parameter with no measurable
    consequence in this design.

    WHAT THE 183 WELLS DO, stated because "scored on real wells" reads stronger than it is.
    The copy column is constant across them, so the two models' weight vectors differ by the
    scalar ``1 - loss*reference`` and their fits are the SAME fit reparameterised: identical
    residual RMS to machine precision and identical in-sample predictions. The wells supply
    ``mu_bar`` and nothing else, and the effect is exactly
    ``mu_bar * loss * (query - reference) / (1 - loss*reference)``. Nothing here is frozen and
    the parameter counts are equal, so the protocol in `mech/ablation.py` is satisfied -- but
    no well discriminates the two models, and a reader should not take this ablation for a
    measurement of the burden slope. It is a measurement of whether a DECLARED increment,
    priced at Kafri's slope, exceeds the growth assay's own error at this plate's growth rate.
    """
    observable = growth_rate_observable(
        query_copies, reference_copies=reference_copies, data=data)
    return ablate(BurdenedGrowth(), UnburdenedGrowth(), observable, GROWTH_RATE_FLOOR)


def copies_to_clear_the_floor(*, reference_copies: float = 1.0,
                              data: pd.DataFrame | None = None,
                              max_copies: int = 40) -> int:
    """The smallest declared copy number whose ablation beats the measured growth floor.

    Raises:
        ValueError: if no copy number up to ``max_copies`` clears it.
    """
    frame = plate_growth_rates() if data is None else data
    for n in range(int(math.ceil(reference_copies)) + 1, max_copies + 1):
        if score_growth_ablation(n, reference_copies=reference_copies,
                                 data=frame).clears_floor:
            return n
    raise ValueError(
        f"no copy number up to {max_copies} clears {GROWTH_RATE_FLOOR.noise_floor:g} /h "
        f"from a reference of {reference_copies:g}")


# --------------------------------------------------------------------------------------
# The GEM coupling, through fba/stress_proteostasis.py and not around it
# --------------------------------------------------------------------------------------

def displacement_over_measured(mu: float, phi_per_copy: float) -> float:
    """How much steeper the GEM-side displacement law is than the measured burden.

    `fba/stress_proteostasis.py::ribosome_allocation_penalty` assumes a heterologous
    fraction displaces every native sector including the ribosomes, giving a loss of
    ``f * (r0 + slope*mu)/(slope*mu)``. Kafri measured ``0.01 * n`` directly. The ratio is
    ``(1 + r0/(slope*mu)) * phi_per_copy / 0.01``, independent of copy number because both
    laws are linear -- so one number answers it at a given mu.

    It is 2.8-3.4x at the plate's MEDIAN growth rate and 4.9-5.9x at the fed-batch
    setpoint, both taken across the phi band; over the plate's full IQR growth band it
    widens to 2.7-3.8x. The displacement assumption is not the burden measurement, and
    this is by how much.
    """
    per_copy = PHI_PER_COPY.at(phi_per_copy)
    return ribosome_allocation_penalty(per_copy, mu) / float(KAFRI_GROWTH_LOSS_PER_COPY)


def gem_burden_comparison(model, n_copies: float, mu: float, phi_per_copy: float,
                          sink_reaction: str = "PROTSYN_hetprot",
                          constrain=None,
                          product_class: str = ProductClass.EGFP_LIKE) -> pd.DataFrame:
    """What the GEM charges for this genotype's protein load, beside what was measured.

    Calls into `fba/stress_proteostasis.py` for every stoichiometric quantity -- the amino
    acid composition off ``r_4047``, the translation ATP, the sink and the growth coupling
    -- and adds the one row that module leaves as NaN: Kafri's measured growth loss. The
    module already reports that its chemistry-only price is roughly half the measured cost;
    this puts the measured cost on the same axis as a copy number.

    Args:
        model: A cobra model ALREADY carrying the heterologous sink
            (``stress_proteostasis.add_heterologous_protein_sink``).
        n_copies: Declared genotype.
        mu: Growth rate the comparison is reported at, 1/h.
        phi_per_copy: A point pinned inside :data:`PHI_PER_COPY`.
        sink_reaction: The sink to load.
        constrain: Optional callable applied to the scoped model to set the uptake regime.
        product_class: Which measured ceiling applies. Declared, not inferred -- 8 copies at
            1.9%/copy is 15.2%, already past the EGFP-like ceiling and admissible only for a
            product declared non-toxic.
    """
    from ..fba.stress_proteostasis import (
        burden_curve,
        grams_per_gdcw_at_protein_fraction,
        read_protein_composition,
    )

    f_H = require_admissible(proteome_fraction(n_copies, phi_per_copy), product_class)
    composition = read_protein_composition(model)
    load = grams_per_gdcw_at_protein_fraction(composition, f_H)
    curve = burden_curve(model, [0.0, load], sink_reaction=sink_reaction,
                         constrain=constrain)
    gem_loss = 1.0 - float(curve.iloc[1]["relative_growth"])
    measured = growth_loss(n_copies)
    rows = [
        {"source": "this GEM: amino acids + translation ATP",
         "basis": "measured here, fba/stress_proteostasis.py",
         "relative_growth_loss": gem_loss},
        {"source": "Kafri 2016 (PMID 26725116), 1% per genomic copy",
         "basis": "MEASURED on the growth axis directly",
         "relative_growth_loss": measured},
        {"source": "proteome displacement, Metzl-Raz allocation law",
         "basis": "measured phi_R law + the displacement assumption",
         "relative_growth_loss": ribosome_allocation_penalty(f_H, mu)},
    ]
    frame = pd.DataFrame(rows)
    frame.insert(0, "f_H", f_H)
    frame.insert(0, "n_copies", n_copies)
    frame["over_measured"] = frame["relative_growth_loss"] / measured
    return frame


def provenance(heading_level: int = 3) -> str:
    """Both registries as markdown, with each arm's gate line at the top of its table."""
    return provenance_table(GROWTH_PARAMS, ALLOCATION_PARAMS, heading_level=heading_level)
