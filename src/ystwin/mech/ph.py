"""Medium pH and a weak acid, reaching growth, the reporter, the redox probe and
eleven enzyme reactions.

WHAT THIS IS FOR. `REVISED_BUILD_LIST.md` downgraded this block from BUILD FIRST to
**BUILD_REDUCED and gated it on a count**: `ARCHITECTURE_GAPS.md` 0.3 tallies the full pH
slice at **eight free scalars against two usable targets**, and criterion (e) refuses to
build a piece where free > targets. What follows is the shape that passes -- **two free
scalars against three independent targets** -- and every one of the six scalars that was
dropped is recorded in :data:`NOT_BUILT` with the measurement that would bring it back.

THE CHAIN, and the depth of each branch counted from the medium::

    pH_ex, A_e (total weak acid)                                        DEPTH 0
        |  Henderson-Hasselbalch, exact, MEASURED pKa
    [AH]_o = A_e / (1 + 10^(pH_ex - pKa))                               DEPTH 1
        |  only the NEUTRAL species crosses; P_AH MEASURED in yeast
    J_in = (S/V) P_AH ([AH]_o - [AH]_i)                                 DEPTH 2
        |
    A_i  (state, tau 0.27-20 min)  ->  [A-]_i trapped                   DEPTH 3
        |
    pH_c (ALGEBRA, not a state: tau 0.0017-0.0029 h at the working dose)  DEPTH 4
         COMPUTED at 40 mM acetate / pH_ex 4.0, and the same band in both windows
         because it is set by entry rather than by growth. Charge balance --
         beta (pH_0 - pH_c) = [A-]_i, which holds at every instant because only the
         NEUTRAL species crosses -- makes pH_c an exact function of A_i, and
         :func:`weak_acid_rhs` now conserves it. The band it used to carry
         (2.9-6.4 h on the plate, 10.3-15.4 h in fed-batch, tau ~ 1/mu) was the
         signature of a MISSING DILUTION TERM in that equation and is RETRACTED
        |
        +-- B1  Pma1 proton bill -> s_ATP(regime) -> mu                 DEPTH 5-7
        +-- B3  Citrine quench   -> RFU        THE MEASUREMENT CHANNEL  DEPTH 5-6
        +-- B4  apparent E_GSH   -> mV         THE SECOND ONE           DEPTH 5-6
        +-- B5  enzyme V_max     -> 11 reactions, MEASURED, unequal signs DEPTH 5-6
        +-- B2/B6/B7/B8  V-ATPase, Msn2, Rim101, Haa1                   NOT BUILT

B5 is where "many reactions at different depths" stops being a slogan: one ``pH_c`` sets
eleven glycolytic ``V_max`` values through Luzia 2022's MEASURED table (PMID 35429225, in
CEN.PK113-7D), and they do not share a sign -- the two directions of GAPDH move 83% down and
40% up over the SAME excursion. It costs no free scalar and no state. What it does not do is
score. At every dose this block works at, the branch is null at one bound (pH_c does not
move) and outside Luzia's assayed 6.19-7.9 at the other, which :func:`enzyme_branch_verdict`
computes and :func:`enzyme_branch_dose_ceiling_mM` puts a number on -- and its assay's
missing noise floor is named in :data:`NOT_BUILT` rather than replaced by the plate CV.

Two of the branches are on the instrument rather than on the cell: pH corrupts the two
channels that would measure the stress, and both correcting functions -- `photophysics.py`'s
``citrine_ph_response`` and `generator/redox.py`'s ``apparent_shift_from_ph`` -- already
existed in this repository and were **called by neither data-generating path**. Computing
``pH_c`` is what makes them reachable, and this module calls both.

THE TWO BOUNDS, WHICH ARE THIS MODULE'S ANSWER TO A REFUSED PARAMETER. A steady state at
fixed pH_c is impossible: at Abbott's measured half-yield acetate (105 mM, pH_ex 5.0) it
demands ``[A-]_i = 8013 mM``. So the cell must export anion or let pH_c fall, and *which*
is set by the Tpo2/Tpo3/Pdr12 capacity -- for which no measurement exists in any source
this pass could reach. The architecture's instruction was to sweep it. Sweeping an invented
range would be fitting with extra steps, so this module **brackets** it instead, with two
computable ends:

* :func:`no_pump_ph` -- nothing is exported; the cytosolic buffer absorbs every proton. One
  free scalar (``beta``). It is charge balance rather than a bound, so the repaired ODE
  settles ON it (5.67657 against 5.67594 at beta 200) instead of crossing it. Its transient
  is an unfitted match to Ullah 2012's 2 min READ -- and only that: it has no recovery arm,
  so it cannot match what Ullah measured after.
* :func:`proton_bill` -- pH_c is held at Orij's MEASURED resting value and the pump pays;
  report the bill. One free scalar (``V_cyt``).

Neither bound is the cell. Both are computed, both are reported, and the interval between
them is the honest size of the missing measurement. :data:`NOT_BUILT` names it.

WHAT IS MEASURED HERE AND WHAT IS NOT, because this block has a specific history of the
second being mistaken for the first. ``P_AH`` is MEASURED **in S. cerevisiae** (Gabba 2020,
PMID 31843263) for exactly three acids -- acetic, formic, lactic -- and this module
**REFUSES propionic, benzoic and sorbic by name** rather than substituting the vesicle
values from the same table, which are 192-707x larger. ``s_ATP`` is **not** measured: it is
the gradient of a linear program, it is a two-valued step at the onset of overflow, and
this module takes it from `bridge/maintenance_calibration.py` per regime rather than as the
scalar 0.004642, which is **2.32x too small at the operating point this repository
validates its GEM against**.

WHAT THE ABLATIONS ACTUALLY ABLATE, SAID BEFORE THEY ARE READ. All four score a CONSTANT
and not a state, and only the last integrates anything:
:func:`lactate_ablation` removes ``P_AH``, :func:`growth_separability_ablation` removes the
medium-pH partitioning and :func:`reporter_quench_ablation` removes the quench law. The
growth branch in particular is algebraically *identical* to ``mu_ref / (1 + dose/K)`` -- an
instantaneous Hill of n = 1 in total dose, verified to machine precision in
``tests/test_mech_ph.py``. What the mechanism buys there is not the shape, it is that
``K = 1 / (s_ATP V_cyt 10^(pH_c - pKa) f_undiss(pH_ex))`` is COMPUTED from a measured pKa
rather than fitted, so the EC50 becomes a measured function of medium pH -- 521 mM total
acetate at pH_ex 4.0 against 2884 mM at pH 5.5, a factor of 5.53 that a single EC50 on
total dose cannot produce at any value. That is a real and falsifiable claim and it is
still not a claim about a state. AND THIS BLOCK HAS NO ABLATION THAT REMOVES A STATE.
:func:`reporter_state_ablation` was that answer and it is RETRACTED: repaired, it scores
0.0051x the reporter floor and FAILS criterion (a), because what it used to score at 2.34x
was :func:`weak_acid_rhs`'s own charge-balance error. Two of the other three carry declared
conditions of their own -- see :func:`reporter_quench_ablation`, which a fitted incumbent
absorbs at any single medium pH, and :func:`growth_separability_ablation`, whose pass holds
on the fermentative branch and fails on the respiratory one. What survives unqualified is
:func:`lactate_ablation`, at 98.97x a floor measured inside the same experiment.

THE HEADLINE RESULT IS A MEASURED NEGATIVE. ``P_lactic = 0``: Gabba 2020 concludes "the
yeast plasma membrane is impermeable to lactic acid on timescales up to ~2.5 h". Abbott 2008
(PMID 18676708), independently and in chemostats, found that transferring lactate's
half-yield concentration from pH 5.0 to pH 3.5 by Henderson-Hasselbalch under-predicts by
**almost 9-fold (85 mM predicted, 750 mM observed)** while the same transfer works for
benzoate to within **9.4%**. This module gets both from ``P_AH`` with **no fitted parameter
on either side**, and :func:`lactate_ablation` scores it at **99x** the floor that
benzoate row measures. See :func:`abbott_half_yield_table`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.optimize import brentq, curve_fit

from ..bridge.maintenance_calibration import MaintenanceRegime, REFERENCE_REGIME
from ..generator.redox import NERNST_MV_PER_PH, apparent_shift_from_ph
from ..pathway.thermo_gate import CYTOSOLIC_VOLUMES_ML_PER_GDCW
from ..photophysics import citrine_ph_response
from .ablation import (
    GROWTH_RATE_FLOOR,
    REPORTER_ACTIVITY_FLOOR,
    AblationResult,
    Fit,
    FittableModel,
    Observable,
    ablate,
)
from .integrate import ReducedSystem, Trajectory, integrate_window, reduce_for_window
from .params import Param, ParamRegistry, Target, provenance_table
from .state import PLATE_READ_4H, Encoding, MechState, StateVar, Window

__all__ = [
    "ABBOTT_HALF_YIELD_COLUMNS",
    "ACIDS",
    "ENZYME_SIGN_ABOVE_REFERENCE",
    "ENZYME_VMAX_CHANGES",
    "LUZIA_ENZYMES",
    "CYTOSOLIC_BUFFER_CAPACITY",
    "CYTOSOLIC_VOLUME_ML_PER_GDCW",
    "LACTATE_TRANSFER_FLOOR",
    "NOT_BUILT",
    "PH_GATE",
    "PH_PARAMS",
    "REFUSED_ACIDS",
    "RESTING_CYTOSOLIC_PH",
    "SURFACE_TO_VOLUME_PER_CM",
    "AcidUnmeasured",
    "EnzymeBranchVerdict",
    "EnzymeIntervalUnmeasured",
    "NoPumpBound",
    "ProtonBill",
    "WeakAcid",
    "abbott_half_yield_table",
    "acid",
    "apparent_redox_shift_mv",
    "atp_borne_fraction_of_orij_slope",
    "atp_route_slope_decades_per_ph",
    "buffer_capacity",
    "charge_imbalance_mM",
    "citrine_quench",
    "entry_rate_constant_per_h",
    "enzyme_branch_dose_ceiling_mM",
    "enzyme_branch_in_range",
    "enzyme_branch_verdict",
    "enzyme_vmax_change",
    "growth_separability_ablation",
    "growth_separability_physicality",
    "lactate_ablation",
    "no_pump_ph",
    "orij_growth_per_h",
    "permeation_tau_s",
    "ph_states",
    "proton_bill",
    "proton_pumping_atp",
    "quench_identifiability",
    "provenance_summary",
    "reduction_audit",
    "reporter_quench_ablation",
    "reporter_state_ablation",
    "simulate_weak_acid",
    "trapped_anion",
    "trapped_total",
    "undissociated_fraction",
    "undissociated_outside",
    "weak_acid_rhs",
]

_LN10 = math.log(10.0)


class AcidUnmeasured(NotImplementedError):
    """An acid was asked for whose yeast permeability nobody has measured.

    ``NotImplementedError`` to match `params.NoSingleValue` and
    `kinetic/carotenoid.py::calibrated_kinetics`: the request is reasonable and the
    measurement is absent, which is a different fault from a malformed call.
    """


class EnzymeIntervalUnmeasured(NotImplementedError):
    """A pH interval, or a pH, at which nobody has measured a yeast enzyme's ``V_max``.

    Raised instead of interpolating :data:`ENZYME_VMAX_CHANGES` or extrapolating past
    Luzia 2022's assayed 6.19-7.9. Same grade as :class:`AcidUnmeasured`: the request is
    reasonable and the measurement is absent.
    """


# --------------------------------------------------------------------------------------
# The constants, graded. Two are free; the rest are pinned or live outside the gate.
# --------------------------------------------------------------------------------------

PH_PARAMS = ParamRegistry("mech/ph.py::weak_acid")
"""The built block's constants. Two free scalars against three independent targets --
:data:`PH_GATE` computes it at import and :meth:`ParamRegistry.require_gate` refuses if it
ever stops being true."""

NOT_BUILT = ParamRegistry("mech/ph.py::not-built")
"""Every constant the sources name that this module does not implement, and why.

Deliberately outside the gate, on `mech/population.py::NOT_BUILT`'s precedent: a refusal
that no equation consumes is a record rather than a degree of freedom. Six of these are the
scalars `ARCHITECTURE_GAPS.md` 0.3 counted, and dropping them is what makes the count pass.
"""

#: Every ci95 row here is the source's printed +/-, which is one sigma and not a 95% band.
_ONE_SIGMA = (
    " THE PRINTED INTERVAL IS THE SOURCE'S OWN +/- AND IS ONE SIGMA, NOT A 95% BAND: "
    "`params.Param.ci95` documents itself as a 95% interval, so this row is narrower than "
    "the field name claims by about 1.96x. Widening it would be inventing a distribution "
    "the paper does not print.")

#: No table is named here on purpose: every row below names its own locus, and the two that
#: are not Table 2 (the pKas, the geometry) previously read "Table 2 Table 1".
_GABBA = ("Gabba, Frallicciardi, van 't Klooster, Henderson, Syga, Mans, van Maris & "
          "Poolman 2020, Biophys J 118:422-434, PMID 31843263 (PMC6976801)")

PKA_ACETIC = PH_PARAMS.add(Param.measured(
    "pka_acetic", 4.76, "dimensionless",
    "MEASURED, standard thermodynamic value at 25 C, carried in " + _GABBA + " Table 1 and "
    "used by Abbott 2008 (PMID 18676708) for the same conversion this module runs"))

PKA_FORMIC = PH_PARAMS.add(Param.measured(
    "pka_formic", 3.75, "dimensionless", "MEASURED, standard; " + _GABBA + " Table 1"))

PKA_LACTIC = PH_PARAMS.add(Param.measured(
    "pka_lactic", 3.86, "dimensionless",
    "MEASURED, standard. Abbott 2008 (PMID 18676708) states 'a pKa value of 3.86 for lactic "
    "acid' in the very sentence that makes the Henderson-Hasselbalch prediction this module "
    "tests, so the test and the paper use the same number"))

PKA_BENZOIC = PH_PARAMS.add(Param.measured(
    "pka_benzoic", 4.20, "dimensionless",
    "MEASURED, standard. Checked against Abbott 2008 (PMID 18676708), which reports 2 mM "
    "benzoate at pH 5 and 0.3 mM at pH 3.5 as '0.27 mM and 0.25 mM undissociated acid': "
    "pKa 4.20 reproduces those as 0.2736 and 0.2501, so the paper's own arithmetic pins it"))

P_AH_ACETIC = PH_PARAMS.add(Param.measured(
    "permeability_acetic", 1.4e-5, "cm/s",
    "MEASURED in S. cerevisiae plasma membrane (strains RA380 and MG10), " + _GABBA + ", "
    "Table 2. The paper's prose transposes acetic and formic; its Table 2 is consistent with "
    "its own P_vesicle/P_yeast ratios (990/1.4 = 707, 210/1.1 = 191), so the table is right. "
    "DO NOT substitute the vesicle value 990e-5 cm/s from the same table. Table 2 prints "
    "1.4 +/- 0.2." + _ONE_SIGMA,
    ci95=(1.2e-5, 1.6e-5)))

P_AH_FORMIC = PH_PARAMS.add(Param.measured(
    "permeability_formic", 1.1e-5, "cm/s",
    "MEASURED in S. cerevisiae, " + _GABBA + ", Table 2, which prints 1.1 +/- 0.1."
    + _ONE_SIGMA,
    ci95=(1.0e-5, 1.2e-5)))

P_AH_LACTIC = PH_PARAMS.add(Param.measured(
    "permeability_lactic", 0.0, "cm/s",
    "MEASURED NEGATIVE, and it is this block's headline. " + _GABBA + ", whose Table 2 "
    "carries NO yeast entry for lactic acid at all, only the 6e-5 cm/s vesicle value, and "
    "whose abstract states: "
    "'we conclude that the yeast plasma membrane is impermeable to lactic acid on timescales "
    "up to ~2.5 h', across wild type and a manifold knockout lacking all putative lactic-acid "
    "transporters. It is an upper bound below that assay's detection rather than a proven "
    "zero; that vesicle value is already the smallest of the five acids measured there and "
    "165x below acetic's"))

CELL_VOLUME_FL = PH_PARAMS.add(Param.measured(
    "cell_volume", 81.9, "fL",
    "MEASURED WHOLE-CELL volume, " + _GABBA + ", its cell-geometry parameters: 'we consider a "
    "spherical cell with volume V0 = 81.9 fL, a volume Vr at zero-turgor pressure of 66.6 "
    "fL, and nonosmotic volume b = 0.65 Vr = 43.3 fL'. Cross-checked against the same "
    "paper's unstirred-layer analysis, which sets the cell's characteristic length (its "
    "diameter) to 5.4 um: a sphere of 81.9 fL is 5.388 um across. THIS is the volume that "
    "fixes the membrane AREA, and taking it from the same paper as P_AH is deliberate -- "
    "P_AH and the geometry it multiplies must come from one measurement or their product "
    "is a splice"))

CELL_WATER_VOLUME_FL = PH_PARAMS.add(Param.measured(
    "cell_water_volume", 38.6, "fL",
    "MEASURED OSMOTIC volume, " + _GABBA + ": 'the osmotic volume V0 - b = 38.6 fL is "
    "filled with an aqueous solution at pH ~6.5'. This is the compartment a concentration "
    "in this module is expressed in, NOT the cell. THE ADVERSARIAL PASS FOUND IT USED AS "
    "ONE: 38.6 fL was tagged `reference_cell_volume` and the membrane area was computed "
    "from a sphere of THAT volume, i.e. a 4.19 um cell, giving S/V = 1.431e4 /cm. The "
    "paper's cell is 5.4 um across, so the area belongs to V0 and the volume to V0 - b, "
    "and the correct ratio is 1.65x larger. Nothing scored moved -- no ablation here uses "
    "k_entry -- but the reported entry time constants did"))

_R_CM = (3.0 * float(CELL_VOLUME_FL) * 1e-12 / (4.0 * math.pi)) ** (1.0 / 3.0)
_SV_PER_CM = 4.0 * math.pi * _R_CM ** 2 / (float(CELL_WATER_VOLUME_FL) * 1e-12)

SURFACE_TO_VOLUME_PER_CM = PH_PARAMS.add(Param.derived(
    "surface_to_volume", _SV_PER_CM, "1/cm",
    f"DERIVED from Gabba 2020 (PMID 31843263) as AREA-OVER-CELL-WATER, which is what the "
    f"flux law needs: the neutral species crosses the membrane of the whole cell and "
    f"dissolves in the cell water. r0 = (3 V0/4pi)^(1/3) = {_R_CM:.6g} cm from V0 = "
    f"{float(CELL_VOLUME_FL):g} fL, A = 4 pi r0^2 = {4.0 * math.pi * _R_CM ** 2:.6g} cm2, "
    f"V_water = {float(CELL_WATER_VOLUME_FL):g} fL, A/V_water = {_SV_PER_CM:.6g} /cm. The "
    "sphere is the MINIMUM-AREA shape at fixed volume, so this is still a LOWER bound on "
    "S/V for a budding cell and therefore an UPPER bound on every entry time constant "
    "computed from it. ARCHITECTURE_TARGET.md quotes 1.39e4, the pH verification pass "
    "1.316e4 and this module previously 1.431e4; all three are 3/r on a volume that is the "
    "cell water rather than the cell"))

RESTING_CYTOSOLIC_PH = PH_PARAMS.add(Param.measured(
    "resting_cytosolic_ph", 7.08, "dimensionless",
    "MEASURED by ratiometric pHluorin in live growing cultures, Orij, Urbanus, Vizeacoumar, "
    "Giaever, Boone, Nislow, Brul & Smits 2012, Genome Biol 13:R80, PMID 23021432 "
    "(PMC3506951), read from the paper itself: 7.08 +/- 0.05 at pH_ex 5.0 over 2,088 "
    "biological replicates, 7.08 +/- 0.03 at pH_ex 3.0, 7.20 +/- 0.15 at pH_ex 7.5. In its "
    "own words, 'titration of pH_ex leads to only minute changes in pH_c'. That invariance "
    "is the answer to 'why does the medium-pH knob do nothing on its own': it acts only "
    "through weak-acid partitioning. The +/- 0.05 is the standard deviation across those "
    "2,088 replicates." + _ONE_SIGMA,
    ci95=(7.03, 7.13)))

ORIJ_INTERCEPT = PH_PARAMS.add(Param.measured(
    "orij_intercept", 7.24, "pH units",
    "MEASURED. Orij 2012 (PMID 23021432), read from the paper: '190 independent 12-hour "
    "pH_c-growth curves of the wild-type strain were used to fit the variables in "
    "mu = 10^((pH_c + a)/b) where a was determined to be -7.24 +/- 0.06, and b to 0.73 +/- "
    "0.08'. The sign is theirs; this module carries |a| and subtracts. Separately, and NOT "
    "the same analysis, the paper reports R^2 = 0.99 for the linear log(mu)-vs-pH_c relation "
    "over t = 2-16 h from 24 biological replicates. HELD OUT: nothing here is fitted to it."
    + _ONE_SIGMA,
    ci95=(7.18, 7.30)))

ORIJ_SLOPE = PH_PARAMS.add(Param.measured(
    "orij_slope", 0.73, "pH units per decade of mu",
    "MEASURED, Orij 2012 (PMID 23021432): b = 0.73 +/- 0.08, i.e. 1/0.73 = 1.370 decades of "
    "mu per unit pH_c, a factor of 23.4. Established in the CAUSAL direction in the same "
    "paper: titrating a pma1-007 hypomorph from pH_ex 5.0 to 2.0 with HCl moved pH_c 7.0 -> "
    "6.2 and slowed growth 'quantitatively matching the relationship... observed in wild "
    "type. This proves that pH_c directly controls growth rate.' The +/- 0.08 means a +/-0.1 "
    "pH_c uncertainty is a 1.4x band on mu, and this module reports the band." + _ONE_SIGMA,
    ci95=(0.65, 0.81)))

ABBOTT_REFERENCE_GLUCOSE = PH_PARAMS.add(Param.measured(
    "abbott_reference_q_glucose", 6.03, "mmol glucose/gDCW/h",
    "MEASURED, Abbott, Suir, van Maris & Pronk 2008, Appl Environ Microbiol 74:5759-5768, "
    "PMID 18676708 (PMC2547041), Table 1 reference column: anaerobic glucose-limited "
    "chemostat, CEN.PK113-7D, D = 0.10 /h, pH 5.0, 25 g/L glucose, no acid added. Table 1 "
    "prints -6.03 +/- 0.10." + _ONE_SIGMA,
    ci95=(5.93, 6.13)))

ANAEROBIC_CATABOLIC_ATP = PH_PARAMS.add(Param.derived(
    "anaerobic_catabolic_atp", 2.0 * float(ABBOTT_REFERENCE_GLUCOSE), "mmol ATP/gDCW/h",
    f"DERIVED from Abbott 2008 (PMID 18676708) Table 1's reference column: "
    f"{float(ABBOTT_REFERENCE_GLUCOSE):g} mmol glucose/gDCW/h x 2 ATP/glucose "
    f"(alcoholic fermentation, established stoichiometry) = "
    f"{2.0 * float(ABBOTT_REFERENCE_GLUCOSE):g} mmol ATP/gDCW/h. The whole catabolic ATP "
    "supply at Abbott's own reference condition, used here ONLY to express the pump bill as "
    "a percentage. It is NOT the 12 mmol/gDCW/h the architecture anchors Jmax_pump on -- that "
    "number came from Abbott 2007 stressed fluxes the verification pass could not retrieve"))

_LUZIA = ("Luzia, Lao-Martil, Savakis, van Heerden, van Riel & Teusink 2022, FEBS J "
          "289:7133-7160, PMID 35429225 (PMC9790636), doi 10.1111/febs.16459, in "
          "S. cerevisiae CEN.PK113-7D")

LUZIA_PH_LOW = PH_PARAMS.add(Param.measured(
    "luzia_assay_ph_low", 6.19, "dimensionless",
    "MEASURED, the low end of the pH range over which V_max was assayed for all ten "
    "enzymes, " + _LUZIA + ". Below it this module REFUSES to evaluate the enzyme branch "
    "rather than extrapolating a measured curve past its own data"))

LUZIA_PH_HIGH = PH_PARAMS.add(Param.measured(
    "luzia_assay_ph_high", 7.9, "dimensionless",
    "MEASURED, the high end of the same range, " + _LUZIA))

_LUZIA_CONFOUND = (" The authors' own caveat travels with every row: adjusting the buffer pH "
                   "also changed K+ by 11% and Na+ by 48% across the series, so pH and "
                   "monovalent cation composition are partially confounded in these assays.")

GAPDH_FWD_FAMINE = PH_PARAMS.add(Param.measured(
    "vmax_gapdh_fwd_7p1_to_6p4", -0.83, "fraction of V_max",
    "MEASURED, " + _LUZIA + ": GAPDH forward loses 83% of V_max over pH 7.1 -> 6.4, which is "
    "the feast/famine cytosolic excursion the same paper measured by pHluorin in a "
    "microfluidic device. The largest pH sensitivity in the pathway, in the enzyme "
    "previously identified as its flux bottleneck." + _LUZIA_CONFOUND))

GAPDH_REV_FAMINE = PH_PARAMS.add(Param.measured(
    "vmax_gapdh_rev_7p1_to_6p4", 0.40, "fraction of V_max",
    "MEASURED, " + _LUZIA + ": GAPDH REVERSE GAINS 40% over the same pH 7.1 -> 6.4 "
    "excursion that costs the forward direction 83%. The sign inverts between the two "
    "directions of one enzyme, which is why no single 'acid damages enzymes' scalar can "
    "stand in for this table." + _LUZIA_CONFOUND))

ALD_STEP = PH_PARAMS.add(Param.measured(
    "vmax_ald_6p8_to_6p2", -0.88, "fraction of V_max",
    "MEASURED, " + _LUZIA + ", pH 6.8 -> 6.2." + _LUZIA_CONFOUND))

HXK_STEP = PH_PARAMS.add(Param.measured(
    "vmax_hxk_6p8_to_6p2", -0.82, "fraction of V_max",
    "MEASURED, " + _LUZIA + ", pH 6.8 -> 6.2." + _LUZIA_CONFOUND))

PGM_STEP = PH_PARAMS.add(Param.measured(
    "vmax_pgm_6p8_to_6p2", -0.81, "fraction of V_max",
    "MEASURED, " + _LUZIA + ", pH 6.8 -> 6.2." + _LUZIA_CONFOUND))

GAPDH_FWD_STEP = PH_PARAMS.add(Param.measured(
    "vmax_gapdh_fwd_6p8_to_6p2", -0.70, "fraction of V_max",
    "MEASURED, " + _LUZIA + ", pH 6.8 -> 6.2." + _LUZIA_CONFOUND))

GAPDH_FWD_ALKALINE = PH_PARAMS.add(Param.measured(
    "vmax_gapdh_fwd_7p8_to_6p8", -1.0 + 1.0 / 6.0, "fraction of V_max",
    "MEASURED, " + _LUZIA + ": a 6-fold INCREASE in GAPDH forward V_max from pH 6.8 to 7.8, "
    "so the reverse traverse 7.8 -> 6.8 is -83.3%. THE INTERVAL MATTERS AND WAS GOT WRONG "
    "ONCE: the pH harvest recorded this 6-fold as spanning 'the full pH range', and the "
    "verification pass reading PMC9790636 found the paper says 6.8 to 7.8, which is not the "
    "full 6.19-7.9 assayed range." + _LUZIA_CONFOUND))

CYTOSOLIC_BUFFER_CAPACITY = PH_PARAMS.add(Param.swept(
    "beta_cytosolic", "mM per pH unit",
    "SWEPT over 90-200. Upper end: Kahm 2012 (PMID 22737060) Table S2 asserts 200 mM/pH, "
    "traceable only to Blatt & Slayman 1987 (Neurospora) and Grabe & Oster 2001 (a mammalian "
    "organelle model); Gerber 2016 independently asserts 200 mM/pH labelled 'experimental "
    "observation' IN YEAST. Lower end: Gabba 2020 (PMID 31843263) fitted 160-320 mM cytosolic "
    "phosphate, which converts to a buffering CAPACITY at pH 7 through "
    "beta = 2.303 C Ka[H]/(Ka+[H])^2 = 0.55 C for pKa2 = 7.2, giving 88-176 mM/pH. THE pKa "
    "IN THAT CONVERSION IS A SPLICE AND THE ADVERSARIAL PASS PRICED IT: the same paper sets "
    "the cytosolic buffer's pKa to a MEASURED 6.59 ('the pKa of the yeast cytosol as "
    "measured in [30]'), and at 6.59 the same 160-320 mM gives 68-136 mM/pH at pH 7.08, so "
    "the swept lower end of 90 is above the bottom of its own derivation. The bound is left "
    "at 90 because lowering it only deepens every acidification here and enlarges every "
    "effect scored, and this module scores at the conservative end. The often "
    "quoted 4-8x spread is a category error between a concentration and a capacity; the real "
    "spread is ~2x. Poznanski 2013 (PMID 23206695) is said to give 41 mM/pH on a yeast "
    "cytosol EXTRACT and that number could not be verified from the paper, so it is NOT used",
    bounds=(90.0, 200.0),
    missing="a pH-clamp titration of living S. cerevisiae -- a measured acid pulse of known "
            "size against a simultaneous ratiometric pHluorin read, which pins beta directly. "
            "The instrument does not exist here today: the manifest is single-channel "
            "mCitrine 480/530 and pHluorin needs dual excitation near 390/470"))

CYTOSOLIC_VOLUME_ML_PER_GDCW = PH_PARAMS.add(Param.swept(
    "cytosolic_volume", "mL cell water per gDCW",
    "SWEPT over the range this repository already declares for the same conversion: "
    f"`pathway/thermo_gate.py::CYTOSOLIC_VOLUMES_ML_PER_GDCW` = "
    f"{CYTOSOLIC_VOLUMES_ML_PER_GDCW}, whose own docstring says 'nothing in this repository "
    "measures it, so no single value is offered as a default'. It enters here only where a "
    "concentration becomes a per-gram flux, i.e. the Pma1 ATP bill; the partitioning results "
    "and the whole lactate test are independent of it",
    bounds=(min(CYTOSOLIC_VOLUMES_ML_PER_GDCW), max(CYTOSOLIC_VOLUMES_ML_PER_GDCW)),
    missing="a wet-weight/dry-weight determination with an inulin or dextran space to "
            "exclude the wall and the periplasm, on this project's own strain and medium. "
            "`calib/od.py` already refuses the neighbouring gdcw_per_od conversion for the "
            "same reason"))


# --------------------------------------------------------------------------------------
# Targets: the three independent measurements this block can be refuted by
# --------------------------------------------------------------------------------------

_ABBOTT_ASSAY = ("50%-biomass-yield titration in anaerobic glucose-limited chemostat, "
                 "CEN.PK113-7D, D = 0.10 /h, 25 g/L glucose")

LACTATE_TRANSFER_FLOOR = PH_PARAMS.add_target(Target(
    name="half_yield_transfer",
    observable=("half-yield weak-acid concentration transferred across a pH change by "
                "Henderson-Hasselbalch"),
    assay=_ABBOTT_ASSAY,
    noise_floor=0.09402301019348747,
    units="fraction",
    source=(
        "MEASURED ON THE POSITIVE CONTROL IN THE SAME EXPERIMENT, which is why it is not a "
        "borrowed floor. Abbott 2008 (PMID 18676708) reports benzoic acid at total "
        "concentrations of 2 mM (pH 5) and 0.3 mM (pH 3.5) giving 'the same degree of "
        "reduction of the biomass yield, thus confirming that for benzoic acid, toxicity is "
        "mediated predominantly by the undissociated species'. Transferring the pH-5.0 point "
        "by Henderson-Hasselbalch at pKa 4.20 predicts 0.328207 mM total at pH 3.5 against "
        "0.3 mM observed: a relative error of 0.094023. That is this assay's own demonstrated "
        "accuracy at transferring a half-yield concentration across 1.5 pH units when the "
        "mechanism IS the undissociated species, so it is the floor an effect on the same "
        "quantity must beat"),
))
"""9.4%, measured on benzoate, the acid the same paper shows obeys Henderson-Hasselbalch.

Criterion (a) forbids borrowing the plate-reader CV for an observable it did not measure. A
chemostat half-yield concentration has no floor on this platform, so one is named here from
the positive control inside the very experiment the effect is scored on."""

GROWTH_TARGET = PH_PARAMS.add_target(GROWTH_RATE_FLOOR)
"""0.0117 /h, ABSOLUTE. `generator/panel_experiment.py::MEASURED_GROWTH_RATE_SE`, the
standard error of the log-OD slope across 210 real wells. Scores the medium-pH separability
prediction of :func:`growth_separability_ablation`."""

REPORTER_TARGET = PH_PARAMS.add_target(REPORTER_ACTIVITY_FLOOR)
"""0.146 CV. Scores the Citrine-quench branch, which is the one place this block acts on the
instrument rather than on the cell."""

PH_GATE = PH_PARAMS.require_gate()
"""Criterion (e), computed at import: 2 free scalars against 3 independent targets.

`ARCHITECTURE_GAPS.md` 0.3 counted the unreduced slice at 8 against 2 and the gate would
have refused it. The six scalars that went are in :data:`NOT_BUILT`.

WHAT THE 3 IS AND IS NOT, because the gate counts NAMED ASSAYS and not DATASETS, and for
this block those differ. Exactly one of the three targets has data in hand:
``half_yield_transfer`` is Abbott 2008's published chemostats. The other two name real
instruments this project owns and rows it has never run -- there is no acid well and no
medium-pH well anywhere in this repository -- so today they can score a prediction and
cannot refute one. Worse for the count than that: **neither free scalar enters the one
target that has data.** ``beta`` and ``V_cyt`` are both absent from the lactate test, which
is parameter-free on both sides. So on a strict reading of criterion (e) -- targets that
could currently refute -- this block is 2 free against 1 and the gate would REFUSE it. It
passes on the reading the gate implements, which is the reading `ARCHITECTURE_GAPS.md` 0.3
used, and the difference is one 96-well plate. Recorded here rather than left for the next
reader to rediscover; `mech/params.py::Target` has no field for it."""


# --------------------------------------------------------------------------------------
# Refusals. Each one is a constant a source names, that no equation here consumes.
# --------------------------------------------------------------------------------------

NOT_BUILT.add(Param.refused(
    "jmax_pump", "mmol ATP/gDCW/h",
    "ARCHITECTURE_TARGET.md sweeps 12-24, anchored at the low end on 'the extra anaerobic "
    "catabolic ATP at Abbott's half-yield point, ~2 x delta-q_glucose ~ 12' and at the high "
    "end on Kahm 2012's I_Pma1max = 16 uA/cm2",
    reason="both anchors fail. The 12 rests on Abbott 2007 stressed fluxes (12.98, 21.45 "
           "mmol/gDCW/h) that the pH verification pass established appear in NO retrievable "
           "source -- Abbott 2007 (PMID 17484738) is Cloudflare-blocked and was never opened, "
           "and the baseline numbers quoted with them are in fact Abbott 2008's reference "
           "column. The 16 uA/cm2 is Neurospora crassa (Gradmann 1978, PMID 25343; Sanders, "
           "Hansen & Slayman 1981, PMID 6458045) and no S. cerevisiae electrophysiological "
           "I_max was found. Rather than sweep a band with no defensible ends, this module "
           "computes the pump flux the mass balance REQUIRES and reports it",
    missing="an S. cerevisiae plasma-membrane patch clamp giving I_max for Pma1, or "
            "institutional access to Abbott 2007's stressed chemostat flux table"))

NOT_BUILT.add(Param.refused(
    "pka_pma1", "dimensionless",
    "Zhao et al. 2021, Nat Commun 12:6439, PMID 34750373: the purified native S. cerevisiae "
    "Pma1 hexamer is autoinhibited at pH 7.4 by a C-terminal regulatory helix that becomes "
    "disordered at pH 6.0. ARCHITECTURE_TARGET.md sweeps pKa in [6.0, 7.4] and n in [1, 4]",
    reason="TWO pH points are a mechanism, not a titration: no pKa and no Hill coefficient "
           "exists. A midpoint swept over 1.4 pH units and a steepness swept 4x is the "
           "instantaneous Hill this architecture exists to replace, moved from dose-space to "
           "pH-space (ARCHITECTURE_GAPS.md 2.1). The claim that without it 'the loop has no "
           "gain' is also wrong: Kahm's thermodynamic term already carries "
           f"2.303 RT/F = {NERNST_MV_PER_PH:.1f} mV per pH unit acting in the correcting "
           "direction, which is a physical constant rather than a fitted one",
    missing="a pH titration of purified Pma1 ATPase activity between 5.5 and 7.5, which is "
            "one enzyme assay on protein Zhao 2021 already purified"))

NOT_BUILT.add(Param.refused(
    "n_pma1", "dimensionless",
    "The Hill coefficient of the same autoinhibition release, swept over [1, 4]",
    reason="see pka_pma1: Zhao 2021 reports no cooperativity coefficient",
    missing="the same titration"))

NOT_BUILT.add(Param.refused(
    "vmax_anion_export", "mmol/L cell water/h",
    "Tpo2/Tpo3/Pdr12, the ATP-driven anion exporters. ARCHITECTURE_TARGET.md DEPTH 3 has "
    "J_exp = Vmax_ex [A-]_i/(K_ex + [A-]_i) and tags both constants REFUSED while requiring "
    "the equation to be built, which ARCHITECTURE_GAPS.md 1.9 correctly notes makes the "
    "block raise on construction",
    reason="no Vmax for Tpo2, Tpo3 or Pdr12 exists in units that close the anion balance. "
           "Abbott 2008 (PMID 18676708) Table 7 measured that a tpo2 tpo3 double deletion has "
           "NO growth defect at 500 or 750 mM lactic acid at pH 3, while haa1 does not grow "
           "at all at 750 mM -- so even the sign of these transporters' contribution is not "
           "settled in the one organism. This module brackets the term between "
           "no_pump_ph and proton_bill instead of sweeping an invented range",
    missing="a transport assay giving Vmax and Km for Pdr12 or Tpo2 on their actual anion "
            "substrates. Failing that, the bracket collapses if intracellular acetate is "
            "measured directly -- the Haa1-based acetate biosensor of Mormino, Siewers & "
            "Nygard 2021 (FEMS Yeast Res, PMID 34477863) is buildable on this plate reader"))

NOT_BUILT.add(Param.refused(
    "k_anion_export", "mM",
    "The half-saturation of the same export term",
    reason="see vmax_anion_export", missing="the same transport assay"))

NOT_BUILT.add(Param.refused(
    "g_proton_leak", "uS/cm2",
    "Kahm 2012 Table S2 carries g_H,leak = 25 uS/cm2 and ARCHITECTURE_TARGET.md writes the "
    "leak as J_leak = g_H (pH_c - pH_c,ref)",
    reason="two independent faults. The value is Neurospora crassa (Sanders, Hansen & "
           "Slayman 1981, PMID 6458045) and no yeast value was found; and the equation is "
           "dimensionally broken -- a conductance times a pH is not a flux. The physical form "
           "g_H (V_m - E_H) (A/V)/F needs a membrane potential, which the target architecture "
           "never writes down (ARCHITECTURE_GAPS.md 3.2). Building the leak would require "
           "membrane_potential as well, so both are refused together",
    missing="an S. cerevisiae plasma-membrane proton conductance, measured with V_m"))

NOT_BUILT.add(Param.refused(
    "membrane_potential", "V",
    "V_m enters the Pma1 driving force as (V_m - E_H + dG_ATP/F)",
    reason="the S. cerevisiae literature spans ~130 mV and the two methods disagree "
           "systematically: microelectrode recordings give -45 to -70 mV, dye partitioning "
           "gives far more negative values, and Kahm 2012 initialises at -200 mV. Against "
           "dG_ATP/F = 414 mV that is a third of the whole driving force, so a swept V_m is a "
           "swept pump. The two measurements disagree; this module does not pick one",
    missing="a patch-clamp resting potential on this strain, or a calibrated "
            "voltage-sensitive dye cross-checked against one"))

NOT_BUILT.add(Param.refused(
    "dg_atp_cytosolic", "kJ/mol",
    "The free energy of ATP hydrolysis in the Pma1 driving force. Kahm 2012 asserts -40 "
    "kJ/mol from Gradmann 1978 (Neurospora); Ke 2013 uses -36.03 plus fixed adenylates",
    reason="no measured dG_ATP trajectory under weak acid exists, and holding it constant is "
           "wrong exactly where it matters: Ullah, Chandrasekaran, Brul & Smits 2013, Front "
           "Microbiol, PMID 23781215 MEASURED that ATP falls NON-monotonically -- 'a stronger "
           "reduction of ATP with growth reducing than with growth inhibitory concentrations' "
           "-- the opposite of the naive expectation",
    missing="a time-resolved cytosolic ATP/ADP/Pi determination under an acetic-acid dose "
            "ladder. The repository's QUEEN-2m reporter is the instrument for it"))

NOT_BUILT.add(Param.refused(
    "vacuolar_proton_flux", "mmol/L/h",
    "The V-ATPase sequestration branch (B2)",
    reason="TWO independent grounds, and the second is a contradiction in this organism. "
           "(i) The Grabe & Oster 2001 pump surface everyone imports is organism-unverified: "
           "the paper does not state which organism the whole-vacuole patch clamp came from. "
           "(ii) Martinez-Munoz & Kane 2008 (PMID 18502746) find vma mutants have a much more "
           "acidic cytosol and call the V-ATPase's role 'unexpectedly prominent'; Dechant "
           "2010 (PMID 20581803) find that deleting V-ATPase subunits or inhibiting with "
           "concanamycin A 'did not affect the dynamics of cytosolic pH' and conclude the "
           "causality runs the other way. Building the term commits you to Kane, omitting it "
           "commits you to Dechant. THIS MODULE COMMITS TO DECHANT, for a reason from Kane's "
           "own data: vma mutants have 65-75% less Pma1 activity in plasma-membrane fractions "
           "and mislocalised Pma1, so Kane's cytosolic-pH defect is explicable at the Pma1 "
           "level with no vacuolar term",
    missing="an S. cerevisiae whole-vacuole J_H(dpH, dPsi) surface, plus a vacuolar volume "
            "fraction and buffering capacity in this organism"))

NOT_BUILT.add(Param.refused(
    "rim101_activation", "dimensionless",
    "The alkaline branch: pH_ex -> Rim13 cleavage -> processed Rim101 -> V-ATPase subunit "
    "expression (B7)",
    reason="no kinetic model of the Rim101 pathway exists in ANY organism. Missing: the pH_ex "
           "threshold and steepness, the Rim13 cleavage rate, the half-life of processed "
           "Rim101, and the gain onto V-ATPase expression. Two structural facts for whoever "
           "builds it: it is a PROTEOLYSIS cascade, so its reset is degradation and "
           "resynthesis rather than a phosphatase -- a structurally different ODE from every "
           "other branch here; and the processed form is already detectable in log-phase YPD "
           "at pH ~5.5, so the input function is not a switch off zero",
    missing="any time-resolved measurement of processed Rim101 against a step in pH_ex"))

NOT_BUILT.add(Param.refused(
    "msn2_gain_from_ph", "nuclear fraction per pH unit",
    "The crosstalk edge pH_c -> V-ATPase disassembly -> Ras/PKA -> Msn2 -> STRE (B6). "
    "`generator/stress_panel.py`'s crosstalk graph has four driven_by weights across three "
    "modules and 'ph' appears in none of them",
    reason="Dechant 2010 (PMID 20581803) is entirely qualitative -- images and timings, no "
           "dose-response -- and no pH_c -> Msn2 -> STRE transfer function exists in any "
           "organism. Drawing the edge as zero would assert independence; drawing it with a "
           "guess is the failure mode this project exists to avoid",
    missing="titrate pH_c with the pma1-007 hypomorph exactly as Orij 2012 did, and read "
            "pHluorin and an STRE reporter IN THE SAME WELLS. Since STRE/HSP12 is this "
            "repository's most-used reporter, this is the highest-value missing measurement "
            "in the whole pathway"))

NOT_BUILT.add(Param.refused(
    "haa1_gain_from_anion", "dimensionless",
    "Trapped anion -> Haa1 nuclear relocation -> HRK1 (a Pma1 activator) and TPO2/TPO3 (B8). "
    "`generator/stress_panel.py` routes acetic acid to {ph, ESR, cell_wall} and has no Haa1 "
    "module at all, although Haa1 is the acetic-acid transcription factor in this organism",
    reason="no measured transfer function from [A-]_i to any Haa1 target transcript exists. "
           "This is the branch that would give pH_c a time constant of its own on the hours "
           "scale -- Ullah 2012 (PMID 23001666) MEASURED that pH_i RECOVERY, not the initial "
           "drop, determines growth inhibition, and Kahm 2012 measured Pma1 ACTIVITY rising "
           "over ~10 min and decaying over hours -- so its absence is why this module's pH_c "
           "carries no slow recovery arm",
    missing="RT-qPCR of TPO2/TPO3/HRK1 against an acetate dose ladder with pH_c read in the "
            "same wells; or the Haa1-based acetate biosensor of PMID 34477863"))

for _name, _pka, _logp in (("propionic", 4.87, 0.33), ("benzoic", 4.20, 1.87),
                           ("sorbic", 4.76, 1.33)):
    NOT_BUILT.add(Param.refused(
        f"permeability_{_name}", "cm/s",
        f"Abbott 2008 (PMID 18676708) quotes a measured half-yield concentration for this "
        f"acid at pH 5 and its log P = {_logp:g}, and ARCHITECTURE_TARGET.md's test E3 asks "
        f"for all four acids",
        reason=("no S. cerevisiae permeability exists for it. Gabba 2020 (PMID 31843263) "
                "measured the yeast membrane for acetic, formic and lactic ONLY; its benzoic "
                "value (1.0e-1 cm/s) is a POPE:POPG:POPC VESICLE, and the yeast membrane is "
                "192-707x less permeable than that bilayer in the same table. Substituting it "
                "is the exact error the same paper forbids, so E3 is executed here for "
                "ACETATE ONLY and the other three are refused by name, per "
                "REVISED_BUILD_LIST.md"),
        missing=(f"the stopped-flow pyranine assay of Gabba 2020 run on {_name} acid across "
                 f"the yeast plasma membrane. Note the paper's abstract states the ORDER "
                 f"benzoic > acetic > formic > lactic holds in yeast too, so a yeast benzoic "
                 f"number may exist in its Table 2 and was not extracted by this pass")))

NOT_BUILT.add(Param.refused(
    "enzyme_vmax_curve", "fraction of V_max per pH unit",
    "The CONTINUOUS V_max(pH) surface for all ten enzymes of " + _LUZIA + " -- HXK, "
    "PGI, PFK, ALD, TPI, GAPDH forward and reverse, PGM, ENO, PYK, PDC -- "
    "which is what a rate law needs to be multiplied by at an arbitrary pH_c",
    reason="the paper's per-enzyme curves were not extracted by this pass, so interpolating "
           "between the intervals it does report would be inventing the shape of a measured "
           "curve. :func:`enzyme_vmax_change` therefore serves ONLY the intervals the paper "
           "states and raises on every other pair. Note the branch has no time constant of "
           "its own -- a V_max modulation tracks pH_c instantaneously -- so building it adds "
           "no state, which is why it costs nothing to leave at interval resolution",
    missing="the MATLAB dataset at github.com/DavidLaoM/pH_kinetics, which is the only "
            "machine-readable artefact in this whole pathway sweep: data/, models/ and "
            "parameter_estimation/, last pushed 2022-06-02"))

NOT_BUILT.add(Param.refused(
    "enzyme_vmax_assay_floor", "fraction",
    "The noise floor criterion (a) would need to score the enzyme branch: the replicate "
    "spread of V_max in " + _LUZIA + "'s in-vivo-like-buffer assay",
    reason="not extracted by this pass. THE ASSAY IS NAMED AND THE FLOOR IS NOT INVENTED: "
           "the stopping rule forbids borrowing OBSERVED_ACTIVITY_CV, which is a plate "
           "reader's CV on dilution-corrected reporter activity and measures nothing about "
           "a maximum-likelihood V_max from a progression curve. So this block reports the "
           "enzyme branch and does NOT claim an ablation on it",
    missing="the replicate spread of V_max at one pH from PMC9790636's File S1, which is "
            "open access"))

NOT_BUILT.add(Param.refused(
    "permeability_co2", "cm/s",
    "Dissolved CO2 is the dominant acidifier of a dense culture: Orij 2012 (PMID 23021432) "
    "drove wild-type pH_c from 7.0 to 6.7 in mid-to-late exponential phase BEFORE glucose ran "
    "out, and observed it 'even in well shaken Erlenmeyer flask cultures'. In a 96-well plate "
    "that is worse, not better, so this is a density-dependent acidification in every well of "
    "every plate this project has run",
    reason="Gabba 2020 (PMID 31843263) states plainly that CO2 permeability 'has, to the best "
           "of our knowledge, not been measured in yeast' and substitutes 1e-3 cm/s from "
           "GUINEA PIG COLON apical epithelium. This module therefore CANNOT explain "
           "density-dependent pH_c, and says so rather than omitting the term silently",
    missing="a CO2 permeability for the S. cerevisiae plasma membrane"))


# --------------------------------------------------------------------------------------
# The acids
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class WeakAcid:
    """One acid, with the two constants that decide everything downstream.

    Args:
        name: Identifier used in tables and messages.
        pka: The dissociation constant, MEASURED.
        permeability: ``P_AH`` across the *yeast* plasma membrane, MEASURED. Zero is a
            legitimate value and is the point of :data:`ACIDS`'s lactic row.
        log_p: Octanol-water partition coefficient, carried because Abbott 2008 names it as
            the physical reason for the permeability ordering and it is the only thing that
            makes the lactate result predictable rather than merely observed.
    """

    name: str
    pka: Param
    permeability: Param
    log_p: float

    @property
    def permeates(self) -> bool:
        """Whether the neutral species measurably crosses the yeast plasma membrane."""
        return float(self.permeability) > 0.0


ACIDS: Mapping[str, WeakAcid] = {
    "acetic": WeakAcid("acetic", PKA_ACETIC, P_AH_ACETIC, -0.31),
    "formic": WeakAcid("formic", PKA_FORMIC, P_AH_FORMIC, -0.54),
    "lactic": WeakAcid("lactic", PKA_LACTIC, P_AH_LACTIC, -0.60),
}
"""The three acids with a MEASURED yeast permeability, and only those three.

log P for acetic (-0.31) and lactic (-0.60) is quoted by Abbott 2008 (PMID 18676708) in the
sentence that names lipid solubility as the reason lactate is an outlier; formic's -0.54 is
the standard value and is carried for ordering only."""

REFUSED_ACIDS: Mapping[str, str] = {
    "propionic": "permeability_propionic",
    "benzoic": "permeability_benzoic",
    "sorbic": "permeability_sorbic",
}
"""Acids on the repository's stress panel or in Abbott's series with no yeast permeability.

:func:`acid` raises :class:`AcidUnmeasured` for these, carrying the refusal from
:data:`NOT_BUILT`. Benzoic is here despite having a measured *pKa* -- its pKa is used for
the floor calculation, but its *permeability* is a vesicle number and this module will not
run the chain on it."""


def acid(name: str) -> WeakAcid:
    """The acid by name, or a refusal naming the measurement that is missing."""
    key = str(name).strip().lower()
    if key in ACIDS:
        return ACIDS[key]
    if key in REFUSED_ACIDS:
        refusal = NOT_BUILT[REFUSED_ACIDS[key]]
        raise AcidUnmeasured(
            f"{key!r} has no measured S. cerevisiae permeability. {refusal.reason} "
            f"What would close it: {refusal.missing}")
    raise AcidUnmeasured(
        f"no acid named {key!r}; this module carries {', '.join(sorted(ACIDS))} with a "
        f"measured yeast permeability and refuses {', '.join(sorted(REFUSED_ACIDS))} by name")


# --------------------------------------------------------------------------------------
# DEPTH 0-3: partitioning, entry and the trapped anion. No free scalars anywhere here.
# --------------------------------------------------------------------------------------

def undissociated_fraction(ph: float, pka: float) -> float:
    """Fraction of a monoprotic weak acid that is undissociated at this pH.

    ``1/(1 + 10^(pH - pKa))``. Exact, and the only place medium pH and dose meet -- which is
    why they are not separable knobs and why one Hill EC50 on dose cannot express them at
    any value.
    """
    return float(1.0 / (1.0 + 10.0 ** (float(ph) - float(pka))))


def undissociated_outside(total_mM: float, ph_ex: float, pka: float) -> float:
    """``[AH]_o``, the only species that crosses, in mM."""
    if total_mM < 0:
        raise ValueError(f"total acid must be non-negative, got {total_mM}")
    return float(total_mM) * undissociated_fraction(ph_ex, pka)


def trapped_total(ah_out_mM: float, ph_c: float, pka: float) -> float:
    """Total intracellular acid at partition equilibrium, mM.

    At equilibrium the NEUTRAL species is equal on both sides, so
    ``A_i = [AH]_o (1 + 10^(pH_c - pKa))``. At Abbott's measured half-yield acetate this
    returns 8051 mM, which is impossible -- the argument that there is no steady state at
    fixed pH_c, and the reason this module carries two bounds rather than one answer.
    """
    return float(ah_out_mM) * (1.0 + 10.0 ** (float(ph_c) - float(pka)))


def trapped_anion(ah_out_mM: float, ph_c: float, pka: float) -> float:
    """``[A-]_i``, the impermeant anion at partition equilibrium, mM.

    One proton was released for each one of these, so this is also the cumulative proton
    load the pump has to have handled.
    """
    return float(ah_out_mM) * 10.0 ** (float(ph_c) - float(pka))


def entry_rate_constant_per_h(acid_name: str) -> float:
    """``(S/V) P_AH`` in 1/h: the first-order rate at which the neutral species equilibrates."""
    a = acid(acid_name)
    return float(SURFACE_TO_VOLUME_PER_CM) * float(a.permeability) * 3600.0


def permeation_tau_s(acid_name: str) -> float:
    """Bare permeation time constant of the NEUTRAL species, seconds.

    3.0 s for acetic and 3.8 s for formic on this geometry. ``ARCHITECTURE_TARGET.md``
    quotes 10 ms and 9.8 ms for the same two, which `ARCHITECTURE_GAPS.md` 1.1 identified as
    the POPE:POPG:POPC vesicle numbers. The geometry-independent statement of that error is
    the paper's own P_vesicle/P_yeast: 707x for acetic and 191x for formic. Raises for an
    impermeant acid rather than returning inf.
    """
    k = entry_rate_constant_per_h(acid_name)
    if k <= 0.0:
        raise AcidUnmeasured(
            f"{acid_name!r} has P_AH = 0 (MEASURED negative, {_GABBA}), so the neutral "
            f"species has no permeation time constant. That is the result, not a gap")
    return float(3600.0 / k)


def buffer_capacity(beta_mM_per_ph: float, a_i_mM: float, ph_c: float, pka: float) -> float:
    """Total buffering capacity, mM/pH: the cytosol plus the trapped acid itself.

    The imported acid is a buffer in its own right, and at pH_c near its pKa a large trapped
    pool is comparable to the whole cytosolic capacity. The second term is the standard
    monoprotic capacity ``2.303 C f (1-f)``, the SAME formula that converts Gabba's fitted
    phosphate concentration into a capacity in :data:`CYTOSOLIC_BUFFER_CAPACITY`'s source, so
    it costs no new constant.
    """
    f = 1.0 - undissociated_fraction(ph_c, pka)
    return float(beta_mM_per_ph) + _LN10 * float(a_i_mM) * f * (1.0 - f)


# --------------------------------------------------------------------------------------
# DEPTH 4: the two bounds on cytosolic pH
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class NoPumpBound:
    """Where pH_c settles if nothing is exported and the cytosolic buffer takes every proton.

    Args:
        ph_c: The buffer-limited floor.
        drop: ``resting_ph - ph_c``.
        beta_mM_per_ph: The swept buffering capacity this was computed at.
        ah_out_mM: The undissociated concentration outside that drove it.
        trapped_anion_mM: ``[A-]_i`` at that floor, which equals the protons delivered.
    """

    ph_c: float
    drop: float
    beta_mM_per_ph: float
    ah_out_mM: float
    trapped_anion_mM: float

    def report(self) -> str:
        return (f"NO-PUMP BOUND  [AH]_o = {self.ah_out_mM:.4g} mM, beta = "
                f"{self.beta_mM_per_ph:g} mM/pH -> pH_c {self.ph_c:.4f} "
                f"(a drop of {self.drop:.3f}), [A-]_i = {self.trapped_anion_mM:.4g} mM")


def no_pump_ph(ah_out_mM: float, beta_mM_per_ph: float, pka: float, *,
               resting_ph: float | None = None) -> NoPumpBound:
    """The buffer-limited cytosolic pH with no export of any kind.

    Solves ``beta (pH_0 - pH_c) = [AH]_o 10^(pH_c - pKa)``: the protons the buffer absorbs
    equal the anion formed once the neutral species has equilibrated across the membrane.
    Both sides are exact; the only input that is not measured is ``beta``. THAT IDENTITY IS
    CHARGE BALANCE and it holds at every instant, not only at equilibrium, because nothing
    but the neutral species crosses -- so this is not really a bound but the constraint the
    dynamics must satisfy, and :func:`weak_acid_rhs` now satisfies it identically. The
    integrated trace settles ON this line rather than crossing it: 5.67657 against 5.67594
    at beta 200, 5.40601 against 5.40562 at beta 90, on a 4.14 h plate at 40 mM acetate /
    pH_ex 4.0. The remaining 4-6e-4 pH units are physical and not numerical -- growth holds
    ``A_i`` a hair BELOW partition equilibrium, since ``J_in = mu A_i > 0`` requires
    ``[AH]_i < [AH]_o`` -- so this is now a bound approached from the correct side. Before
    the repair the same trace crossed it by 0.30-0.32 pH units and kept falling.

    This is the FAST arm, and it is an EMERGENCE check rather than a fit --
    :func:`simulate_weak_acid` reaches 95% of this drop in 0.36-0.63 min across the swept
    beta at 40 mM acetate / pH_ex 4.0 -- the 0.37-0.75 min this said before was the default
    2001-point OUTPUT GRID (7.45 s a step), not the model, and the true band is 5-19%
    shorter. Against Ullah 2012 (PMID 23001666), whose operational
    read is "the minimum pHi reached within 2 min" at 1 s sampling. STATE THE MATCH EXACTLY,
    because the adversarial pass found this overclaimed as "reaches its minimum by 2 min":
    the integrated pH_c has NO interior minimum, it falls monotonically for the whole 4.14 h
    window, and the quantity that arrives by 2 min is 95% of the partition-equilibrium drop.
    The model therefore matches Ullah's 2 min READ and cannot reproduce the recovery that
    follows it -- Pma1 is what recovers pH_i in that paper, and this bound has no pump. See
    ``NOT_BUILT['haa1_gain_from_anion']``.

    Args:
        ah_out_mM: Undissociated acid outside, from :func:`undissociated_outside`.
        beta_mM_per_ph: Cytosolic buffering capacity. :data:`CYTOSOLIC_BUFFER_CAPACITY` is
            SWEPT, so pass a point from its ``grid()`` or ``at()`` and carry the band.
        pka: The acid's pKa.
        resting_ph: Starting pH_c. Defaults to Orij's MEASURED 7.08.
    """
    p0 = float(RESTING_CYTOSOLIC_PH) if resting_ph is None else float(resting_ph)
    beta = float(beta_mM_per_ph)
    ah = float(ah_out_mM)
    if beta <= 0:
        raise ValueError(f"beta must be positive, got {beta}")
    if ah <= 0:
        return NoPumpBound(p0, 0.0, beta, ah, 0.0)

    def residual(p: float) -> float:
        return beta * (p0 - p) - ah * 10.0 ** (p - float(pka))

    low = 1.0
    if residual(low) < 0.0:
        raise ValueError(
            f"[AH]_o = {ah:g} mM at pKa {pka:g} delivers more protons than beta = {beta:g} "
            f"mM/pH can absorb above pH 1; the no-pump bound is off the physical scale, "
            f"which is a finding about the dose rather than a solver failure")
    p = float(brentq(residual, low, p0 - 1e-12, xtol=1e-12))
    return NoPumpBound(p, p0 - p, beta, ah, ah * 10.0 ** (p - float(pka)))


@dataclass(frozen=True)
class ProtonBill:
    """What holding pH_c at its measured resting value costs, at 1 H+ per ATP.

    Args:
        mmol_atp_per_gdcw_h: The bill.
        trapped_anion_mM: ``[A-]_i`` the bound requires.
        fraction_of_catabolic_atp: Against :data:`ANAEROBIC_CATABOLIC_ATP`.
        growth_rate_per_h: The measured mu the bill was computed at; the bill is
            proportional to it, because dilution is the only route by which the anion leaves.
        cytosolic_volume_ml_per_gdcw: The swept axis this was computed at.
        physical: Whether ``trapped_anion_mM`` is below the osmotic scale a cell runs at.
    """

    mmol_atp_per_gdcw_h: float
    trapped_anion_mM: float
    fraction_of_catabolic_atp: float
    growth_rate_per_h: float
    cytosolic_volume_ml_per_gdcw: float
    physical: bool

    def report(self) -> str:
        flag = "" if self.physical else "  -- NOT PHYSICAL, see .physical"
        return (f"PUMP-HOLDS BOUND  [A-]_i = {self.trapped_anion_mM:.5g} mM at mu = "
                f"{self.growth_rate_per_h:g} /h, V_cyt = "
                f"{self.cytosolic_volume_ml_per_gdcw:g} mL/gDCW -> "
                f"{self.mmol_atp_per_gdcw_h:.4g} mmol ATP/gDCW/h = "
                f"{self.fraction_of_catabolic_atp:.1%} of the anaerobic catabolic supply"
                f"{flag}")


#: Above this an intracellular anion pool is not a cell. ASSERTED order-of-magnitude osmotic
#: scale, outside PH_PARAMS on purpose: it only flags, enters no prediction, costs no freedom.
_OSMOTIC_SCALE_MM = 500.0


def proton_pumping_atp(proton_flux: float, *, atp_per_proton: float) -> float:
    if (isinstance(proton_flux, (bool, np.bool_)) or not math.isfinite(proton_flux)
            or proton_flux < 0):
        raise ValueError("proton_flux must be finite and nonnegative")
    if (isinstance(atp_per_proton, (bool, np.bool_)) or not math.isfinite(atp_per_proton)
            or atp_per_proton <= 0):
        raise ValueError("atp_per_proton requires an explicit positive finite coupling ratio")
    return float(proton_flux) * float(atp_per_proton)


def proton_bill(ah_out_mM: float, growth_rate_per_h: float,
                cytosolic_volume_ml_per_gdcw: float, pka: float, *,
                ph_c: float | None = None) -> ProtonBill:
    """The Pma1 ATP bill required to hold pH_c, in mmol ATP/gDCW/h.

    The mass balance, not a capacity. In steady state every imported AH dissociates and its
    anion leaves only by growth dilution, so protons must be exported at exactly
    ``mu [A-]_i``; Pma1 is 1 H+ per ATP, so the proton flux IS the ATP flux and no extra
    constant is needed. Multiplying by the cell-water volume puts it in the units the GEM's
    NGAM reaction uses.

    Args:
        ah_out_mM: Undissociated acid outside.
        growth_rate_per_h: MEASURED mu. In this project mu comes off the biomass series; do
            not assert it from a pump.
        cytosolic_volume_ml_per_gdcw: A point on :data:`CYTOSOLIC_VOLUME_ML_PER_GDCW`.
        pka: The acid's pKa.
        ph_c: The pH being held. Defaults to Orij's MEASURED resting 7.08.
    """
    p = float(RESTING_CYTOSOLIC_PH) if ph_c is None else float(ph_c)
    anion = trapped_anion(ah_out_mM, p, pka)
    litres = float(cytosolic_volume_ml_per_gdcw) / 1000.0
    bill = proton_pumping_atp(float(growth_rate_per_h) * anion * litres, atp_per_proton=1.0)
    return ProtonBill(
        mmol_atp_per_gdcw_h=bill,
        trapped_anion_mM=anion,
        fraction_of_catabolic_atp=bill / float(ANAEROBIC_CATABOLIC_ATP),
        growth_rate_per_h=float(growth_rate_per_h),
        cytosolic_volume_ml_per_gdcw=float(cytosolic_volume_ml_per_gdcw),
        physical=anion <= _OSMOTIC_SCALE_MM,
    )


# --------------------------------------------------------------------------------------
# B1, B3, B4, B5, B9: the branches, each with its depth counted from the medium
# --------------------------------------------------------------------------------------

def growth_under_bill(mu_reference_per_h: float, bill: ProtonBill, *,
                      regime: MaintenanceRegime = REFERENCE_REGIME) -> float:
    """B1. Growth rate after the pump's ATP debit, 1/h.

    ``mu = mu_ref / (1 + s_ATP [A-]_i V_cyt)``, the self-consistent solution of
    ``mu = mu_ref - s_ATP mu [A-]_i V_cyt`` -- the bill is proportional to mu, so the naive
    subtraction double-counts.

    ``s_ATP`` is taken PER REGIME from `bridge/maintenance_calibration.py`, never as the
    scalar. It is a two-valued step at the onset of overflow, and at
    ``fba/physiology.REFERENCE_AEROBIC_BATCH`` -- the operating point this repository
    validates its GEM against -- the fermentative branch is 0.010778 against the respiratory
    0.004642, so the old scalar priced every ATP debit at 43% of its value there. It is also
    NOT a measurement: ``regime.provenance`` says MODEL-DERIVED / IN_SILICO, the gradient of
    a linear program, and no culture was read to obtain it.

    HOW MUCH OF THE MEASURED EFFECT THIS ACCOUNTS FOR, which the block reports nowhere else
    and the second adversarial pass computed. This law's own half-effect dose is **164 mM
    (V_cyt 2.7) to 444 mM (V_cyt 1.0) UNDISSOCIATED acetic acid**, against Abbott 2008's
    MEASURED half-yield of 38.35 mM and this repository's own panel EC50 of 38.72 mM -- so
    the ATP route is **4.3x to 11.6x too weak**, and at 60 mM total acetate at pH_ex 4.5 it
    predicts a growth reduction of 8.0-19.1% where the panel puts the half-effect. On the
    observable Abbott actually measured, the same statement is the bill being 6.6-17.9% of
    the catabolic ATP supply at his half-yield point, i.e. a 6.2-15.2% yield loss by
    Herbert-Pirt at 2 ATP/glucose against his measured 50%. The honest reading is the same
    one :func:`atp_borne_fraction_of_orij_slope` reaches on the other channel: this is a
    PARTIAL mechanism, worth 6-30% of what is measured, not the mechanism.
    """
    if mu_reference_per_h <= 0:
        raise ValueError(f"reference growth rate must be positive, got {mu_reference_per_h}")
    slope = regime.ngam_growth_slope_per_h_per_mmol_atp
    litres = bill.cytosolic_volume_ml_per_gdcw / 1000.0
    return float(mu_reference_per_h / (1.0 + slope * bill.trapped_anion_mM * litres))


def orij_growth_per_h(ph_c: float) -> float:
    """B9. Growth rate from cytosolic pH by Orij's MEASURED law, 1/h.

    ``mu = 10^((pH_c - 7.24)/0.73)``, Orij 2012 (PMID 23021432), fitted to 190 independent
    12 h pH_c-growth curves with R^2 = 0.99 and established in the CAUSAL direction by
    titrating pH_c with the pma1-007 hypomorph. 1.370 decades of mu per pH unit, a factor
    of 23.

    HELD OUT. Nothing in this module is fitted to it, and it is not used inside any
    prediction here -- :func:`growth_under_bill` is the mechanism and this is the check.
    Two caveats travel with it and both are load-bearing. It is a BIJECTION, so a mu
    inferred from OD and a pH_c derived from the law test nothing; use it only where every
    parameter upstream of pH_c is pinned independently. And its own
    ``b = 0.73 +/- 0.08`` makes a +/-0.1 pH_c uncertainty a 1.4x band on mu, so report the
    band -- :func:`orij_growth_band_per_h` does.
    """
    return float(10.0 ** ((float(ph_c) - float(ORIJ_INTERCEPT)) / float(ORIJ_SLOPE)))


def orij_growth_band_per_h(ph_c: float, ph_uncertainty: float = 0.1) -> tuple[float, float]:
    """Orij's law over its own reported parameter interval and a stated pH_c uncertainty."""
    a_lo, a_hi = ORIJ_INTERCEPT.ci95
    b_lo, b_hi = ORIJ_SLOPE.ci95
    values = [10.0 ** ((float(ph_c) + d - a) / b)
              for d in (-abs(ph_uncertainty), abs(ph_uncertainty))
              for a in (a_lo, a_hi) for b in (b_lo, b_hi)]
    return (float(min(values)), float(max(values)))


def atp_route_slope_decades_per_ph(ah_out_mM: float, ph_c: float,
                                   cytosolic_volume_ml_per_gdcw: float, pka: float, *,
                                   regime: MaintenanceRegime = REFERENCE_REGIME) -> float:
    """How many decades of mu the ATP route moves per unit pH_c, AND ITS SIGN IS NEGATIVE.

    ``d log10(mu)/d pH_c = -s_ATP V_cyt [A-]_i`` at fixed ``[AH]_o``, returned here as the
    signed slope. `ARCHITECTURE_TARGET.md` B9 states that the ATP route and Orij's law
    "double-count" and turns that into test E6. **In the presence of a weak acid they do not
    double-count, they oppose**: a HIGHER cytosolic pH deprotonates more of the imported
    acid, traps more anion and costs MORE pumping, while Orij's law says a higher pH_c grows
    faster. This module reports the opposition rather than adding the two.

    E6 as the architecture frames it -- reproduce Orij's 1.37 decades/pH from the ATP route
    alone -- is NOT EXECUTABLE here and that is a refusal, not an omission: with no weak acid
    present this slope is identically zero, because weak-acid influx is the only proton load
    this module carries. Reproducing Orij's slope would need the basal proton leak, whose
    conductance is Neurospora and whose driving force needs a membrane potential; both are
    in :data:`NOT_BUILT`. What IS computable is the magnitude in the presence of acid: at
    Abbott's half-yield acetate it is 0.086-0.233 decades/pH across the ``V_cyt`` sweep,
    i.e. 6-17% of Orij's slope, opposite in sign.
    """
    litres = float(cytosolic_volume_ml_per_gdcw) / 1000.0
    anion = trapped_anion(ah_out_mM, ph_c, pka)
    return float(-regime.ngam_growth_slope_per_h_per_mmol_atp * litres * anion)


def atp_borne_fraction_of_orij_slope(
        ah_out_mM: float, ph_c: float | None = None, *,
        regime: MaintenanceRegime = REFERENCE_REGIME) -> tuple[float, float]:
    """E6, RE-SCOPED: how much of the pH-growth coupling this module's ATP route accounts for.

    `REVISED_BUILD_LIST.md` Phase 1 re-scopes E6 from "held-out check the block passes" to
    "held-out measurement of how much of the pH-growth coupling is ATP-borne". This returns
    that fraction as a band over the swept ``V_cyt``: ``|d log10(mu)/d pH_c| / (1/b)``, with
    Orij's MEASURED ``b = 0.73``.

    THE ANSWER IS 6.3-17.0% AT ABBOTT'S HALF-YIELD ACETATE, and 0.0% with no acid present.
    That is not the 34% (0.473 decades) the build list quotes from the Central Carbon pass,
    and the difference is structural rather than arithmetic: **the only proton load this
    module carries is weak-acid influx**, while Orij titrated pH_c with HCl in a pma1-007
    hypomorph and no acid at all. Reproducing 0.473 decades from these equations needs
    ``V_cyt`` = 5.48 mL/gDCW, 2.03x the top of the band `pathway/thermo_gate.py` declares, so
    this pass does not adopt it. What would close the gap is the basal proton leak -- and its
    conductance is Neurospora and its driving force needs a membrane potential, both in
    :data:`NOT_BUILT`.
    """
    p = float(RESTING_CYTOSOLIC_PH) if ph_c is None else float(ph_c)
    fractions = [abs(atp_route_slope_decades_per_ph(ah_out_mM, p, v, float(PKA_ACETIC),
                                                    regime=regime)) * float(ORIJ_SLOPE)
                 for v in CYTOSOLIC_VOLUME_ML_PER_GDCW.bounds]
    return (min(fractions), max(fractions))


def citrine_quench(ph_c: float, ph_reference: float | None = None) -> float:
    """B3. Fraction of the Citrine signal that survives acidification to ``ph_c``.

    Calls `photophysics.py`'s ``citrine_ph_response`` at the MEASURED chromophore pKa of 5.7
    (Griesbeck 2001, PMID 11387331) and normalises to the resting pH. That function existed
    in this repository and was reachable from a single script and from NEITHER
    data-generating path; computing pH_c is what makes it live.

    The stake, in this project's own numbers: measured induction folds are around 1.5, and
    the no-pump bound at 40 mM acetate / pH_ex 4.0 returns 0.35-0.51 here across the beta
    sweep. So the artefact is one to two thirds of the signal, it has DOSE STRUCTURE because
    the stressors that acidify are the ones being titrated, and it enters every dose-response
    fit as a slope rather than as noise. Caveat at the point of use: the single-protonation
    form is ASSERTED -- Griesbeck reports a pKa and not a Hill coefficient.
    """
    ref = float(RESTING_CYTOSOLIC_PH) if ph_reference is None else float(ph_reference)
    return float(citrine_ph_response(float(ph_c)) / citrine_ph_response(ref))


def apparent_redox_shift_mv(ph_c: float, ph_reference: float | None = None) -> float:
    """B4. Redox shift that acidification alone would appear to produce, mV.

    Calls `generator/redox.py`'s ``apparent_shift_from_ph`` at its DERIVED Nernst slope of
    60.2 mV per pH unit at 303.15 K. That function, too, was called by no ``src`` module.

    Positive means it looks like an oxidation -- the same direction a real oxidant moves it.
    At the no-pump bound for 40 mM acetate this returns +52 to +101 mV across the beta and
    medium-pH range, against a genuine peroxide response of 40-50 mV: the artefact EXCEEDS
    the signal and carries the same sign, so the two cannot be told apart without measuring
    pH. Note the asymmetry that makes this sharper than "pH affects redox": GSH's thiolate
    at pKa ~8.9 moves ~10x over 6.4-7.4 while Tsa1's peroxidatic Cys47 at pKa 5.4 moves
    1.09x, which is inside its own rate constant's error bar -- so pH does real work on the
    glutathione arm and none that can be measured on the peroxiredoxin arm.
    """
    ref = float(RESTING_CYTOSOLIC_PH) if ph_reference is None else float(ph_reference)
    return float(apparent_shift_from_ph(ref, float(ph_c)))


LUZIA_ENZYMES = ("HXK", "PGI", "PFK", "ALD", "TPI", "GAPDH_fwd", "GAPDH_rev",
                 "PGM", "ENO", "PYK", "PDC")
"""The ELEVEN reaction directions of the TEN enzymes Luzia 2022 assayed against pH.

Read from PMC9790636 in the adversarial pass, and it corrects a count this module had
wrong. The paper's own denominator is ten -- "four out of ten enzymes exhibiting a decrease
in activity above 60%" -- and its Methods list PGK (22.5 U/mL), ADH (88 U/mL) and G3PDH
(0.6 U/mL) as COUPLING enzymes added in excess to read other reactions, not as enzymes whose
own V_max was measured against pH. Those three were previously carried here and are removed.

The branch that makes this block "one input reaching MANY reactions": one ``pH_c`` multiplies
eleven ``V_max`` values by eleven different measured factors, and they do not share a
sign."""

ENZYME_VMAX_CHANGES: Mapping[tuple[str, float, float], Param] = {
    ("GAPDH_fwd", 7.1, 6.4): GAPDH_FWD_FAMINE,
    ("GAPDH_rev", 7.1, 6.4): GAPDH_REV_FAMINE,
    ("ALD", 6.8, 6.2): ALD_STEP,
    ("HXK", 6.8, 6.2): HXK_STEP,
    ("PGM", 6.8, 6.2): PGM_STEP,
    ("GAPDH_fwd", 6.8, 6.2): GAPDH_FWD_STEP,
    ("GAPDH_fwd", 7.8, 6.8): GAPDH_FWD_ALKALINE,
}
"""``(enzyme, pH_from, pH_to) -> signed fractional change in V_max``, MEASURED.

Exactly the intervals Luzia 2022 states, and no others. Every key is a pair of pH values the
paper reports a change BETWEEN; :func:`enzyme_vmax_change` will not interpolate to a third."""

ENZYME_SIGN_ABOVE_REFERENCE: Mapping[str, int] = {
    "ENO": +1, "GAPDH_fwd": +1, "PGI": +1, "TPI": +1,
    "PFK": -1, "ALD": -1, "GAPDH_rev": -1, "PDC": -1,
}
"""Which way each enzyme's ``V_max`` moves as pH rises ABOVE the assay's reference pH 6.8.

MEASURED, Luzia 2022 (PMID 35429225). Eight of the eleven are resolved in sign there and
the other three are not carried. THE POINT IS THE DISAGREEMENT: the enzymes split, no clean
upper/lower-glycolysis boundary appears, and the two directions of GAPDH oppose each other --
so a scalar "acid inhibits metabolism" multiplier is refuted by the measurement it would
claim to summarise."""


def enzyme_branch_in_range(ph_c: float) -> bool:
    """Whether ``ph_c`` lies inside the 6.19-7.9 window Luzia 2022 actually assayed."""
    return float(LUZIA_PH_LOW) <= float(ph_c) <= float(LUZIA_PH_HIGH)


def enzyme_vmax_change(enzyme: str, ph_from: float, ph_to: float) -> float:
    """B5. Signed fractional change in one enzyme's ``V_max`` between two pH values.

    ``-0.83`` means V_max falls to 17% of its value at ``ph_from``. Serves only the intervals
    :data:`ENZYME_VMAX_CHANGES` carries and raises otherwise, because the shape of the curve
    between them was not extracted -- see ``NOT_BUILT['enzyme_vmax_curve']``.

    DEPTH 5 from the medium, DEPTH 6 counting the flux it multiplies, on the same footing as
    B3 and B4. It carries NO state and no time constant: a V_max modulation tracks pH_c
    instantaneously, which is exactly why the pH_c time constant this module computes becomes
    the rate-limiting step for the whole downstream response.
    """
    key = (str(enzyme), float(ph_from), float(ph_to))
    if key in ENZYME_VMAX_CHANGES:
        return float(ENZYME_VMAX_CHANGES[key])
    for ph in (ph_from, ph_to):
        if not enzyme_branch_in_range(ph):
            raise EnzymeIntervalUnmeasured(
                f"pH {float(ph):g} is outside Luzia 2022's assayed "
                f"{float(LUZIA_PH_LOW):g}-{float(LUZIA_PH_HIGH):g} (PMID 35429225), so no "
                f"measured V_max exists there for {enzyme}. Extrapolating a measured curve "
                f"past its own data is the failure this module refuses")
    available = sorted(f"{e} {a:g}->{b:g}" for e, a, b in ENZYME_VMAX_CHANGES)
    raise EnzymeIntervalUnmeasured(
        f"no measured V_max change for {enzyme} between pH {float(ph_from):g} and "
        f"{float(ph_to):g}. Measured intervals: {', '.join(available)}. "
        f"{NOT_BUILT['enzyme_vmax_curve'].missing} would close it")


@dataclass(frozen=True)
class EnzymeBranchVerdict:
    """Whether the enzyme branch is evaluable at each end of this module's bracket.

    Args:
        no_pump_ph_c: pH_c at the buffer-limited bound, over the swept beta.
        pump_holds_ph_c: pH_c at the pump-holds bound, which is Orij's measured resting value.
        no_pump_in_range: Whether the no-pump bound lies inside Luzia's assayed range.
        pump_holds_in_range: The same for the pumped bound.
    """

    no_pump_ph_c: tuple[float, float]
    pump_holds_ph_c: float
    no_pump_in_range: bool
    pump_holds_in_range: bool

    def report(self) -> str:
        return (f"ENZYME BRANCH  no-pump bound pH_c {self.no_pump_ph_c[0]:.3f}-"
                f"{self.no_pump_ph_c[1]:.3f} -> "
                f"{'inside' if self.no_pump_in_range else 'OUTSIDE'} Luzia's assayed range; "
                f"pump-holds bound pH_c {self.pump_holds_ph_c:.3f} -> "
                f"{'inside' if self.pump_holds_in_range else 'OUTSIDE'}, where nothing moves "
                f"and the branch is identically zero")


def enzyme_branch_dose_ceiling_mM(acid_name: str = "acetic") -> tuple[float, float]:
    """Largest ``[AH]_o`` whose no-pump bound still lands inside Luzia's assayed range, mM.

    COMPUTED by solving ``no_pump_ph([AH]_o, beta) = 6.19`` at each end of the swept
    buffering capacity: **2.98-6.61 mM undissociated acetic acid**, i.e. 4.6-10.2 mM total at
    pH_ex 4.5. Abbott's measured half-yield is 38.35 mM undissociated and this module's
    working dose (40 mM at pH_ex 4.0) is 34.08, so both sit 5-13x above the ceiling.

    Infinite for an acid with ``P_AH = 0``: nothing crosses, pH_c never leaves Orij's resting
    value, and the branch is in range at every dose -- which is the lactate result again.
    """
    a = acid(acid_name)
    if not a.permeates:
        return (math.inf, math.inf)
    out = []
    for beta in CYTOSOLIC_BUFFER_CAPACITY.bounds:
        out.append(float(brentq(
            lambda x: no_pump_ph(x, beta, float(a.pka)).ph_c - float(LUZIA_PH_LOW),
            1e-9, 1e4)))
    return (min(out), max(out))


def enzyme_branch_verdict(total_acid_mM: float, ph_ex: float, *,
                          acid_name: str = "acetic") -> EnzymeBranchVerdict:
    """Where this module's own two bounds sit relative to Luzia's assayed pH range.

    THE RESULT IS NEGATIVE AT BOTH ENDS AT EVERY DOSE THIS BLOCK WORKS AT, and it is
    computed rather than asserted. At the pump-holds bound pH_c does not move from Orij's
    resting 7.08, so every enzyme factor is exactly 1 and the branch contributes nothing. At
    the no-pump bound a dose above :func:`enzyme_branch_dose_ceiling_mM` -- 2.98-6.61 mM
    undissociated acetate, against Abbott's 38.35 mM half-yield -- drives pH_c below 6.19,
    where no yeast V_max was measured, and the branch refuses to extrapolate. Below that
    ceiling it IS evaluable, and this function reports which case a given dose is in rather
    than asserting either.

    So the enzyme arm is reported and not scored: unusable at one end because it is null and
    at the other because it would be extrapolation. Closing that gap is the same missing
    measurement -- which end the cell is at -- that :data:`NOT_BUILT` names for
    ``vmax_anion_export``.
    """
    a = acid(acid_name)
    pka = float(a.pka)
    resting = float(RESTING_CYTOSOLIC_PH)
    # P_AH = 0 means the neutral species never equilibrates, so there is no partition bound.
    if not a.permeates:
        bounds = (resting, resting)
    else:
        ah_o = undissociated_outside(total_acid_mM, ph_ex, pka)
        lo, hi = CYTOSOLIC_BUFFER_CAPACITY.bounds
        bounds = tuple(sorted(no_pump_ph(ah_o, b, pka).ph_c for b in (lo, hi)))
    return EnzymeBranchVerdict(
        no_pump_ph_c=bounds,
        pump_holds_ph_c=resting,
        no_pump_in_range=all(enzyme_branch_in_range(p) for p in bounds),
        pump_holds_in_range=enzyme_branch_in_range(resting),
    )


# --------------------------------------------------------------------------------------
# The states, the right-hand side, and the tau/T audit that decides whether they are states
# --------------------------------------------------------------------------------------

def _manifold_mode_tau_h(ah_out_mM: float, beta_mM_per_ph: float, pka: float, mu_per_h: float,
                         ph_c: float, k_entry_per_h: float) -> float:
    """Relaxation time of the ONE physical mode of the (A_i, pH_c) system, hours.

    Charge balance makes the system one-dimensional. :func:`weak_acid_rhs` conserves
    ``E = beta (pH_0 - pH_c) - [A-]_i`` up to dilution, the physical initial condition has
    ``E = 0``, so the trajectory never leaves the manifold ``pH_c = g(A_i)`` that
    :func:`no_pump_ph` solves. Differentiating that constraint gives ``dpH_c/dA_i = -f_d/D``
    with ``D`` the total buffering capacity, and substituting it into ``dA_i/dt`` leaves

        1/tau = k_entry (1 - f_d) (1 + 2.303 A_i f_d^2 / D) + mu

    which is the FAST eigenvalue of the full Jacobian evaluated on the manifold -- the other
    eigenvalue is the ``-mu`` decay of an imbalance the physical trajectory never carries.
    ``tests/test_mech_ph.py`` checks this hand-derived rate against :func:`_jacobian`.

    THIS IS AN UPPER BOUND ON pH_c's OWN TIME CONSTANT, and deliberately so: pH_c does not lag
    A_i at all, so its lag is exactly zero and what is returned instead is the rate at which
    its driver moves it. The audit is therefore conservative -- it predicts a reduction bias
    of 0.04-0.07% on the plate, where the measured bias of slaving pH_c instead of
    integrating it is 7.2e-6 at the pinned tolerances and 3.2e-11 at rtol 1e-12.

    Returns ``inf`` when the mode is STATIONARY. At ``P_AH = 0`` the influx is identically
    zero, A_i never leaves zero, pH_c never leaves ``pH_0``, and neither state moves at all;
    ``population.py`` returns ``math.inf`` on the same footing for a plasmid that never
    segregates.
    """
    k = float(k_entry_per_h)
    if k == 0.0:
        return math.inf
    f_u = undissociated_fraction(ph_c, pka)
    a_i = trapped_total(ah_out_mM, ph_c, pka)
    d = buffer_capacity(beta_mM_per_ph, a_i, ph_c, pka)
    rate = k * f_u * (1.0 + _LN10 * a_i * (1.0 - f_u) ** 2 / d) + float(mu_per_h)
    return math.inf if rate == 0.0 else 1.0 / rate


def _jacobian(a_i: float, ph_c: float, ah_out_mM: float, beta: float, pka: float,
              mu: float, k_entry: float) -> np.ndarray:
    """Numerical Jacobian of :func:`weak_acid_rhs` at one point.

    Kept as the CHECK on :func:`_manifold_mode_tau_h`'s hand algebra rather than as the
    source of a time constant: on the charge-balance manifold its fast eigenvalue must equal
    that rate and its other one must be -mu, and ``tests/test_mech_ph.py`` asserts both.
    """
    y0 = np.array([a_i, ph_c], dtype=float)
    base = np.array(weak_acid_rhs(0.0, y0, ah_out_mM=ah_out_mM, k_entry_per_h=k_entry,
                                  pka=pka, beta_mM_per_ph=beta, growth_rate_per_h=mu))
    out = np.zeros((2, 2))
    for j in range(2):
        step = 1e-6 * max(1.0, abs(y0[j]))
        y = y0.copy()
        y[j] += step
        moved = np.array(weak_acid_rhs(0.0, y, ah_out_mM=ah_out_mM, k_entry_per_h=k_entry,
                                       pka=pka, beta_mM_per_ph=beta, growth_rate_per_h=mu))
        out[:, j] = (moved - base) / step
    return out


def charge_imbalance_mM(a_i_mM: float, ph_c: float, beta_mM_per_ph: float, pka: float, *,
                        resting_ph: float | None = None) -> float:
    """``beta (pH_0 - pH_c) - [A-]_i``, the electroneutrality residual, mM.

    Zero is the physical state. Nothing but the NEUTRAL species crosses the membrane, so every
    anion made inside is balanced by one proton on the cytosolic buffer, and that identity holds
    at every instant rather than only at equilibrium -- it is charge balance, not a steady state.
    :func:`no_pump_ph` solves it as an algebraic equation and :func:`weak_acid_rhs` preserves it
    as an invariant: with the dilution term in place ``dE/dt = -mu E`` exactly, so a trajectory
    started at ``E = 0`` never leaves it and one perturbed off it is only diluted back.

    Args:
        a_i_mM: Total weak acid in the cell water.
        ph_c: Cytosolic pH.
        beta_mM_per_ph: Cytosolic buffering capacity.
        pka: The acid's pKa.
        resting_ph: pH_c before the acid arrived. Defaults to Orij's MEASURED 7.08.
    """
    p0 = float(RESTING_CYTOSOLIC_PH) if resting_ph is None else float(resting_ph)
    f_d = 1.0 - undissociated_fraction(float(ph_c), pka)
    return float(float(beta_mM_per_ph) * (p0 - float(ph_c)) - float(a_i_mM) * f_d)


def weak_acid_rhs(t: float, y, *, ah_out_mM: float, k_entry_per_h: float, pka: float,
                  beta_mM_per_ph: float, growth_rate_per_h: float,
                  resting_ph: float | None = None) -> list[float]:
    """``d(A_i, pH_c)/dt`` for the no-pump bound, in mM/h and pH/h.

    ::

        f_d      = 1 / (1 + 10^(pKa - pH_c))                dissociated fraction inside
        J_in     = k_entry ([AH]_o - A_i (1 - f_d))         only the neutral species crosses
        dA_i/dt  = J_in - mu A_i
        dpH_c/dt = (mu beta (pH_0 - pH_c) - f_d J_in) / (beta + 2.303 A_i f_d (1 - f_d))

    Both terms are derived rather than asserted. Protons appear only when imported AH
    dissociates: ``R = d[A-]/dt + mu [A-] = f_d J_in + A_i (df_d/dpH) dpH/dt``, and setting
    ``-beta dpH/dt = R`` gives the denominator, whose second term is the trapped acid's own
    buffering capacity. The ``mu beta (pH_0 - pH_c)`` term is the protons leaving with the
    anion: the buffer that took them is diluted at the same rate the pool it balances is.

    THE DILUTION TERM WAS MISSING AND THE SECOND ADVERSARIAL PASS MEASURED WHAT THAT COST.
    Without it the anion was diluted and the protons it left behind were not, so
    :func:`charge_imbalance_mM` grew in proportion to mu: at 40 mM acetate / pH_ex 4.0
    ``beta (pH_0 - pH_c)`` exceeded ``[A-]_i`` by 1.000x at mu = 0, 1.49x at mu = 0.101 and
    2.45x at mu = 0.245, and by 60-65x by the end of ``FEDBATCH_5D``. With the term the ratio
    is 1.000000 at every mu and in both windows, because ``dE/dt = -mu E`` identically --
    substitute the two lines above into ``E = beta (pH_0 - pH_c) - A_i f_d`` and every J
    cancels.

    WHAT THE REPAIR RETRACTS, stated here because this block's headline rested on the error.
    pH_c is no longer a dynamic state: the invariant makes it an exact algebraic function of
    ``A_i``, the integrated trace now settles ON :func:`no_pump_ph` (5.67657 against 5.67594
    at beta 200, 5.40601 against 5.40562 at beta 90) where before it crossed it by 0.31 pH,
    the pH_c mode's ``tau ~ 1/mu`` was the error's own signature and is gone, and
    :func:`reporter_state_ablation` -- which was scoring the charge-balance error itself --
    falls from 2.34x the reporter floor to 0.0051x and FAILS criterion (a).

    Args:
        t: Unused; the system is autonomous. Present for ``solve_ivp``.
        y: ``(A_i, pH_c)``.
        ah_out_mM: Undissociated acid outside.
        k_entry_per_h: ``(S/V) P_AH``, from :func:`entry_rate_constant_per_h`.
        pka: The acid's pKa.
        beta_mM_per_ph: Cytosolic buffering capacity.
        growth_rate_per_h: MEASURED mu; it dilutes the anion and the buffered protons alike.
        resting_ph: pH_c the buffer is referenced to. Defaults to Orij's MEASURED 7.08, and
            must be the same value the trajectory starts at or the invariant starts non-zero.
    """
    p0 = float(RESTING_CYTOSOLIC_PH) if resting_ph is None else float(resting_ph)
    a_i = max(float(y[0]), 0.0)
    ph = float(y[1])
    mu = float(growth_rate_per_h)
    beta = float(beta_mM_per_ph)
    f_d = 1.0 - undissociated_fraction(ph, pka)
    j_in = float(k_entry_per_h) * (float(ah_out_mM) - a_i * (1.0 - f_d))
    d_a = j_in - mu * a_i
    d_ph = (mu * beta * (p0 - ph) - f_d * j_in) / (beta + _LN10 * a_i * f_d * (1.0 - f_d))
    return [float(d_a), float(d_ph)]


def ph_states(window: Window, *, ah_out_mM: float, acid_name: str = "acetic",
              beta_bounds: tuple[float, float] | None = None) -> MechState:
    """``(A_i, pH_c)`` with time constants COMPUTED for this window, not tabulated.

    Both time constants are properties of the operating point rather than of the molecule,
    so they are derived here from the exact linearisation at the window's MEASURED growth
    band and across the swept buffering capacity -- the same shape as
    `mech/population.py::population_states`.

    * ``A_i`` -- ``tau = 1/(k_entry (1 - f_d(pH_c)) + mu)``, evaluated over pH_c in
      [5.4, 7.4]. The pH dependence is the whole story: the NEUTRAL species equilibrates in
      3 s, but filling the pool means filling the anion too, which is 175x larger at pH 7.
    * ``pH_c`` -- :func:`_manifold_mode_tau_h`, the ONE physical mode of the coupled system
      once charge balance has removed the other. At the working dose (40 mM acetate,
      pH_ex 4.0, ``[AH]_o`` = 34.08 mM) it is 0.0017-0.0029 h across the beta sweep, the
      SAME band in both windows because it is set by entry and not by growth. IT IS AN
      UPPER BOUND: pH_c does not lag A_i at all, and this is the rate at which its driver
      moves it. THE BAND IS DOSE-DEPENDENT and the returned ``source`` string computes it
      for whatever dose was passed rather than repeating this one.

    THIS RETRACTS THE BLOCK'S CRITERION-4 HEADLINE, which was that the pH_c verdict FLIPS
    between the two windows -- irreducible on the plate, algebra in the vessel. That
    flip was an artefact of the missing dilution term in :func:`weak_acid_rhs`: the old
    ``tau ~ 1/mu`` was the charge-balance error's own relaxation and nothing else. Repaired,
    pH_c ELIMINATES in BOTH windows, and it is not even a marginal elimination -- the
    reduction to the slaved pH is EXACT, not merely licensed by a ratio. Slaving pH_c to
    ``A_i`` reproduces the integrated Citrine endpoint to 1.0e-6 and its time-integral to
    7.2e-6 at the package's pinned tolerances, and to 4.4e-12 and 3.2e-11 at rtol 1e-12:
    what is left is the BDF tolerance and not model content.

    Args:
        window: The declared operating window; supplies the measured growth band.
        ah_out_mM: Undissociated acid outside. The time constants depend on the dose.
        acid_name: Which acid; sets ``k_entry``.
        beta_bounds: Buffering-capacity band. Defaults to
            :data:`CYTOSOLIC_BUFFER_CAPACITY`'s swept bounds.
    """
    beta_lo, beta_hi = beta_bounds or CYTOSOLIC_BUFFER_CAPACITY.bounds
    k = entry_rate_constant_per_h(acid_name)
    pka = float(acid(acid_name).pka)
    mu_lo = window.growth_rate_low_per_h
    mu_hi = window.growth_rate_high_per_h

    fast_taus = [1.0 / (k * undissociated_fraction(p, pka) + mu)
                 for p in (5.4, 7.4) for mu in (mu_lo, mu_hi)]
    manifold_taus = [_manifold_mode_tau_h(ah_out_mM, beta, pka, mu, p, k)
                     for beta in (beta_lo, beta_hi) for mu in (mu_lo, mu_hi)
                     for p in (5.6, 6.5)]

    if k > 0.0:
        # The band is computed at this dose, not quoted: it moves with [AH]_o, and a
        # hardcoded one goes stale silently at the next dose.
        ph_source = (
            f"DERIVED as the one physical mode of the coupled (A_i, pH_c) system at "
            f"[AH]_o = {ah_out_mM:g} mM, over beta {beta_lo:g}-{beta_hi:g} mM/pH and mu "
            f"{mu_lo:g}-{mu_hi:g} /h: {min(manifold_taus):.4g}-{max(manifold_taus):.4g} h. "
            f"Charge balance leaves the system one-dimensional, so this is ENTRY-driven and "
            f"not growth-driven -- the other eigenvalue, tau ~ 1/mu, belongs to a charge "
            f"imbalance the physical trajectory never carries. It is an UPPER BOUND on "
            f"pH_c's own time constant: pH_c does not lag A_i at all")
        ph_note = (
            "pH_c ELIMINATES in BOTH windows and the reduction is EXACT rather than "
            "licensed by a ratio: slaving it to the anion reproduces the integrated Citrine "
            "endpoint to 4.4e-12 at rtol 1e-12, which is integrator tolerance. The earlier "
            "'verdict FLIPS between the windows' criterion-4 headline is RETRACTED -- its "
            "tau ~ 1/mu was the missing dilution term in weak_acid_rhs and nothing else")
    else:
        ph_source = (
            f"DERIVED and INFINITE. {acid_name} has a MEASURED P_AH of zero ({_GABBA}), so "
            f"k_entry = 0, the proton load is identically zero, A_i never leaves zero and "
            f"pH_c never leaves its resting value. The mode is STATIONARY, not slow: "
            f"freezing pH_c at Orij's measured resting value is exact, and the true freeze "
            f"bias is 0 at any window length")
        ph_note = (
            "This is the headline stated in reduction language: 900 mM lactic acid at "
            "pH_ex 5.0 leaves pH_c exactly where it started in BOTH windows, with no fitted "
            "parameter anywhere. NOTE the audit row's bias column reads nan rather than 0 "
            "-- mech/integrate.py::qss_bias evaluates (tau/T)(1-exp(-T/tau)) as inf*0 at "
            "tau = inf instead of taking its limit of 1")

    return MechState((
        StateVar(
            name="A_i", units="mM total weak acid in cell water",
            tau_h=(min(fast_taus), max(fast_taus)),
            encoding=Encoding.LEVEL,
            source=(
                f"DERIVED at this window's measured growth band {mu_lo:g}-{mu_hi:g} /h from "
                f"the exact linearisation tau = 1/(k_entry (1 - f_d) + mu), with k_entry = "
                f"(S/V) P_AH = {k:.6g} /h from Gabba 2020's MEASURED P_AH (PMID 31843263) "
                f"and the cell geometry in the same paper, evaluated over pH_c 5.4-7.4"),
            constrained_by=(
                "NOT measured on this platform today. The named route is the Haa1-based "
                "acetate biosensor of Mormino, Siewers & Nygard 2021 (PMID 34477863), which "
                "is buildable on the 96-well reader this project already owns"),
            sweep_axis=True,
            note=("The pH dependence is the point: at pH_c 7.4 the pool is 437x the neutral "
                  "species and takes ~20 min to fill, at 5.4 only 5.4x and under a minute")),
        StateVar(
            name="pH_c", units="dimensionless",
            tau_h=(min(manifold_taus), max(manifold_taus)),
            encoding=Encoding.LEVEL,
            source=ph_source,
            constrained_by=(
                "ratiometric pHluorin, which is in `generator/stress_panel.py`'s REPORTERS "
                "table but is NOT an instrument here: it needs dual excitation near 390/470 "
                "and the manifest is single-channel mCitrine 480/530. New strain AND new "
                "optics, not a free readout"),
            sweep_axis=True,
            note=ph_note),
    ))


def simulate_weak_acid(window: Window, *, total_acid_mM: float, ph_ex: float,
                       beta_mM_per_ph: float, acid_name: str = "acetic",
                       growth_rate_per_h: float | None = None,
                       resting_ph: float | None = None,
                       n_points: int = 2001) -> Trajectory:
    """Integrate the no-pump bound across ``window`` with the package's pinned stiff driver.

    The system spans 3 to 5 orders of magnitude in time constant, so it is integrated with
    `mech/integrate.py`'s BDF rather than an explicit method, and the driver is not
    hand-rolled here.
    """
    a = acid(acid_name)
    pka = float(a.pka)
    ah_o = undissociated_outside(total_acid_mM, ph_ex, pka)
    mu = window.growth_rate_low_per_h if growth_rate_per_h is None else float(growth_rate_per_h)
    k = entry_rate_constant_per_h(acid_name)
    system = ph_states(window, ah_out_mM=max(ah_o, 1e-9), acid_name=acid_name)
    p0 = float(RESTING_CYTOSOLIC_PH) if resting_ph is None else float(resting_ph)

    def rhs(t, y):
        return weak_acid_rhs(t, y, ah_out_mM=ah_o, k_entry_per_h=k, pka=pka,
                             beta_mM_per_ph=beta_mM_per_ph, growth_rate_per_h=mu,
                             resting_ph=p0)

    return integrate_window(rhs, {"A_i": 0.0, "pH_c": p0}, window, system,
                            n_points=n_points, nonnegative=False,
                            atol={"A_i": 1e-6, "pH_c": 1e-9})


def reduction_audit(window: Window, *, ah_out_mM: float,
                    acid_name: str = "acetic") -> ReducedSystem:
    """Criterion (d) for this block, in this window, computed rather than claimed.

    Reduce a state to algebra iff ``tau/T < 0.146``; freeze it iff ``T/tau < 0.146``. The
    ratio and the exact bias are in :meth:`ReducedSystem.audit_table`. BOTH states ELIMINATE
    in BOTH windows, and the pH_c row is the repaired one: it read KEEP on the plate only
    while :func:`weak_acid_rhs` was violating charge balance. For an impermeant acid the
    pH_c row reads FREEZE instead, which is the lactate headline in reduction language.
    """
    return reduce_for_window(ph_states(window, ah_out_mM=ah_out_mM, acid_name=acid_name),
                             window)


# --------------------------------------------------------------------------------------
# Criterion (a): three ablations, three floors, three assays
# --------------------------------------------------------------------------------------

ABBOTT_HALF_YIELD_COLUMNS = ("acid", "ph_ex", "total_mM", "undissociated_mM",
                             "half_yield_total_mM")
"""Column contract of :func:`abbott_half_yield_table`."""


def abbott_half_yield_table() -> pd.DataFrame:
    """Every weak-acid half-yield concentration this pass read out of Abbott 2008 itself.

    Read from the paper (PMID 18676708, PMC2547041), not from a secondary summary, and each
    row reproduces the paper's own arithmetic. The lactate rows are its Results; the benzoate
    rows are the sentence that makes benzoate the positive control; the pH-5 row for acetate
    is quoted there from Abbott 2007 (PMID 17484738), whose own text this pass could NOT
    retrieve -- so that row carries the acid's dose and nothing else.

    Verbatim, for the record: "For cultures grown at pH 5, 900 mM of lactic acid (~61 mM
    undissociated acid) was required to decrease the biomass yield on glucose to 50% ... This
    led to the prediction that a lactic acid concentration of 85 mM (~60 mM undissociated
    acid at pH 3.5) should cause a 50% reduction ... However, experiments showed that the
    required concentration was almost 9-fold higher (750 mM) ... In contrast, benzoic acid at
    total concentrations of 2 mM (pH 5) and 0.3 mM (pH 3.5), corresponding to 0.27 mM and
    0.25 mM undissociated acid, respectively, showed the same degree of reduction of the
    biomass yield."
    """
    rows = [("lactic", 5.0, 900.0), ("lactic", 3.5, 750.0), ("lactic", 3.0, 500.0),
            ("benzoic", 5.0, 2.0), ("benzoic", 3.5, 0.3),
            ("acetic", 5.0, 105.0)]
    pkas = {"lactic": float(PKA_LACTIC), "benzoic": float(PKA_BENZOIC),
            "acetic": float(PKA_ACETIC)}
    data = []
    for name, ph, total in rows:
        undissociated = undissociated_outside(total, ph, pkas[name])
        data.append({"acid": name, "ph_ex": ph, "total_mM": total,
                     "undissociated_mM": undissociated, "half_yield_total_mM": total})
    return pd.DataFrame(data, columns=list(ABBOTT_HALF_YIELD_COLUMNS))


@dataclass(frozen=True)
class _HalfYieldTransfer(FittableModel):
    """Predict a half-yield concentration at one pH from the same acid's value at another.

    Both models below are PARAMETER-FREE: ``parameter_names`` is empty, ``refit`` searches
    nothing, and :func:`ablate`'s frozen-incumbent probe passes them trivially because there
    is nothing to freeze. That is the point of this ablation -- neither side spends a scalar.
    """

    name: str = "half-yield transfer"
    use_permeability: bool = True
    response_column: str = "half_yield_total_mM"
    parameter_names: tuple[str, ...] = ()

    def refit(self, data: pd.DataFrame) -> Fit:
        predicted = self.predict(data, {})
        observed = np.asarray(data[self.response_column], dtype=float)
        residual = float(np.sqrt(np.mean((predicted - observed) ** 2)))
        return Fit(model=self.name, parameters={}, n_free=0, n_rows=len(data),
                   residual_rms=residual, converged=True)

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        out = []
        for _, row in data.iterrows():
            out.append(self.transfer(str(row["acid"]), float(row["anchor_ph"]),
                                     float(row["anchor_total_mM"]), float(row["ph_ex"])))
        return np.asarray(out, dtype=float)

    def transfer(self, acid_name: str, anchor_ph: float, anchor_total_mM: float,
                 target_ph: float) -> float:
        """Half-yield TOTAL concentration at ``target_ph``, mM."""
        pka = {"lactic": float(PKA_LACTIC), "benzoic": float(PKA_BENZOIC),
               "acetic": float(PKA_ACETIC)}[acid_name]
        anchor_undissociated = undissociated_outside(anchor_total_mM, anchor_ph, pka)
        if self.use_permeability and acid_name in ACIDS and not ACIDS[acid_name].permeates:
            # P_AH = 0: the neutral species delivers nothing, so it is not the agent and
            # the only statement left is that the total transfers unchanged.
            return float(anchor_total_mM)
        return float(anchor_undissociated / undissociated_fraction(target_ph, pka))


def _lactate_frame() -> pd.DataFrame:
    """The one row both models are scored on: lactate, pH 5.0 anchor, pH 3.5 target."""
    return pd.DataFrame([{
        "acid": "lactic", "anchor_ph": 5.0, "anchor_total_mM": 900.0, "ph_ex": 3.5,
        "half_yield_total_mM": 750.0}])


def lactate_ablation() -> AblationResult:
    """THE HEADLINE. Removing the MEASURED permeability, scored on Abbott's own chemostats.

    Two parameter-free models predict the total lactic acid needed for a 50% biomass-yield
    reduction at pH 3.5, anchored on the measured pH-5.0 point:

    * FULL keeps ``P_AH``. Lactic acid's is a MEASURED ZERO (Gabba 2020, PMID 31843263), so
      the neutral species is not the agent, Henderson-Hasselbalch has nothing to transfer,
      and the model's remaining statement is that the total carries over: **900 mM**.
    * REDUCED drops it and assumes the undissociated species is the agent -- which is what a
      single EC50 on total dose implies, and what Abbott themselves wrote down as the null:
      **87.3 mM** (the paper rounds to 85).

    Observed: **750 mM**. Effect 9.31 relative, against a floor of 0.094 measured on
    benzoate IN THE SAME EXPERIMENT: **99x**.

    THE FULL ARM CARRIES AN UNSTATED POSITIVE, and the adversarial pass puts a number on
    it. ``P_AH = 0`` establishes only the NEGATIVE -- the neutral species is not the agent.
    It does not establish that the TOTAL is what transfers; that is a third hypothesis
    sitting beside "the anion is the agent", which predicts 2762 mM at pH 3.5 and 6919 mM at
    pH 3.0, residuals of 268% and 1284%. So the choice made here is the best of the three on
    Abbott's own rows by a wide margin -- but it is a choice, not a deduction from ``P_AH``,
    and :func:`lactate_ablation`'s 99x is a score on that choice.

    What must travel with it, because the full model is not exonerated. Its residual against
    the observation is 20% at pH 3.5 -- 2.1x that floor -- and 80% at the pH-3.0 point, where
    the reduced model's is 86%. So the defensible claim is the NEGATIVE: Henderson-Hasselbalch
    must fail for lactate and hold for benzoate, and ``P_AH`` predicts which is which with no
    fitted parameter. Abbott independently reached the same conclusion from the transcriptome
    -- "at pH 5, the yield decrease was caused mostly by osmotically induced glycerol
    production and not by the classic weak-acid action" -- and named the physical cause this
    module's constant encodes: "probably related to the low lipid solubility of lactic acid
    (log P = -0.60). Lipid solubility is strongly correlated with weak-organic-acid toxicity."
    """
    data = _lactate_frame()
    full = _HalfYieldTransfer(name="partition x measured P_AH", use_permeability=True)
    reduced = _HalfYieldTransfer(name="Henderson-Hasselbalch only", use_permeability=False)
    observable = Observable(
        name="half_yield_transfer", units="mM total lactic acid",
        assay=_ABBOTT_ASSAY, scored="relative", data=data,
        summarise=lambda model, fit: float(model.predict(data, fit.parameters)[0]),
        description=("total lactic acid for a 50% biomass-yield reduction at pH 3.5, "
                     "transferred from the measured pH-5.0 point"))
    return ablate(full, reduced, observable, floor=LACTATE_TRANSFER_FLOOR)


@dataclass(frozen=True)
class _SeparabilityGrowth(FittableModel):
    """Growth rate under an acetate dose, with or without the medium-pH partitioning."""

    name: str = "separability growth"
    partition: bool = True
    cytosolic_volume_ml_per_gdcw: float = 1.0
    reference_growth_per_h: float = 0.245
    regime: MaintenanceRegime = REFERENCE_REGIME
    response_column: str = "growth_rate_per_h"
    parameter_names: tuple[str, ...] = ()

    def refit(self, data: pd.DataFrame) -> Fit:
        predicted = self.predict(data, {})
        observed = np.asarray(data[self.response_column], dtype=float)
        return Fit(model=self.name, parameters={}, n_free=0, n_rows=len(data),
                   residual_rms=float(np.sqrt(np.mean((predicted - observed) ** 2))),
                   converged=True)

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        return np.asarray([self.mu(float(r["total_mM"]), float(r["ph_ex"]))
                           for _, r in data.iterrows()], dtype=float)

    def mu(self, total_mM: float, ph_ex: float) -> float:
        """Growth under the bill. Exactly ``mu_ref / (1 + dose/K)``, an n = 1 Hill in dose.

        Written down rather than left to be discovered: the bill is linear in ``[AH]_o`` and
        ``[AH]_o`` is linear in dose, so :func:`growth_under_bill`'s divisor form collapses
        to a Michaelis-Menten in total dose with
        ``K = 1/(s_ATP V_cyt 10^(pH_c - pKa) f_undiss(pH_ex))``. The mechanism sets K; it
        does not change the shape. No state is consulted here and pH_c is held at Orij's
        resting value, which is the pump-holds end of the bracket.
        """
        pka = float(PKA_ACETIC)
        # The reduced model ignores medium pH, which is literally what the package does
        # today: predict.py records CultureContext.ph_medium and it reaches nothing.
        ah_o = (undissociated_outside(total_mM, ph_ex, pka) if self.partition
                else undissociated_outside(total_mM, 5.0, pka))
        bill = proton_bill(ah_o, self.reference_growth_per_h,
                           self.cytosolic_volume_ml_per_gdcw, pka)
        return growth_under_bill(self.reference_growth_per_h, bill, regime=self.regime)


def growth_separability_physicality(
        total_mM: float = 40.0,
        ph_range: tuple[float, float] = (4.0, 5.5)) -> tuple[ProtonBill, ProtonBill]:
    """The two :class:`ProtonBill` ends that :func:`growth_separability_ablation` scores.

    Returned rather than discarded so the ``physical`` flag travels with the effect size.
    At the default 40 mM acetate both ends are ``physical = False`` -- 7120 mM and 1287 mM
    of trapped anion against the 500 mM osmotic scale -- so the growth ablation is a
    statement about the pump-holds bound, not about a cell that could exist.
    """
    pka = float(PKA_ACETIC)
    return tuple(
        proton_bill(undissociated_outside(total_mM, p, pka),
                    PLATE_READ_4H.growth_rate_low_per_h,
                    min(CYTOSOLIC_VOLUMES_ML_PER_GDCW), pka)
        for p in ph_range)


def growth_separability_ablation(
        cytosolic_volume_ml_per_gdcw: float | None = None,
        reference_growth_per_h: float | None = None,
        window: Window = PLATE_READ_4H,
        regime: MaintenanceRegime = REFERENCE_REGIME) -> AblationResult:
    """Medium pH and dose are not separable knobs, scored on growth rate.

    The design is the one the pH verification pass named: acetate crossed with medium pH.
    The observable is the SPREAD of predicted mu across medium pH at a fixed total dose --
    the full model's is 0.014-0.049 /h, the reduced model's is exactly zero, because
    ``CultureContext.ph_medium`` currently reaches nothing.

    Scored against ``MEASURED_GROWTH_RATE_SE`` = 0.0117 /h ABSOLUTE, the standard error of
    the log-OD slope across 210 real wells. Defaults are the CONSERVATIVE corner of both
    swept axes -- ``V_cyt`` at its 1.0 mL/gDCW low end and mu at the low end of
    ``PLATE_READ_4H``'s own measured band, imported from the window rather than retyped --
    so the reported ratio is the smallest the sweep produces, 1.21x. At the other
    corner it is 4.21x, and it clears at every point in between.

    THE RESULT IS REGIME-CONDITIONAL, AND THE REGIME IS NOW A DECLARED ARGUMENT rather than
    an undeclared axis every ratio was quietly evaluated on. Every number above is at
    ``REFERENCE_REGIME`` (fermentative, s_ATP 0.010778). On the respiratory branch
    (0.004642) the same conservative corner gives **0.00638 /h = 0.55x the floor and FAILS**,
    reaching 1.38x only at V_cyt = 2.7. `REVISED_BUILD_LIST.md` E2 attaches exactly this
    caveat -- "for a 96-well plate the oxygen regime is UNKNOWN, so s_ATP(regime) is a
    REFUSAL there until a kLa is measured, not a lookup" -- and this ablation is scored on a
    96-well plate. So the honest statement of the 1.21x is "1.21x IF the wells are
    fermentative", and the oxygen regime of a 96-well plate is not measured here.

    SAY THIS TOO: the piece removed is the partitioning CONSTANT, not a state. Both models
    here are instantaneous algebra -- see :meth:`_SeparabilityGrowth.mu`, which is an n = 1
    Hill in dose at a fixed pH_c -- and neither integrates ``A_i`` or ``pH_c``. The claim
    being scored is that K is a measured function of medium pH, not that a state is needed.

    AND SAY THIS, WHICH THE ADVERSARIAL PASS FOUND UNSTATED: the two ends of the scored
    spread sit at ``[A-]_i`` = 1287 mM (pH_ex 5.5) and 7120 mM (pH_ex 4.0), both flagged
    ``physical = False`` by :class:`ProtonBill` against the 500 mM osmotic scale. That is not
    a bug, it is the module's own headline restated -- holding pH_c at Orij's resting value
    under this much acid is impossible, which is exactly why the bracket exists -- but it
    means this ablation measures the PUMP-HOLDS BOUND and not the cell. The result carries
    the flag: :func:`growth_separability_physicality` returns both ends so a reader cannot
    take the 1.21x without it.

    SAY THIS PLAINLY: there is no acid well and no medium-pH well in this repository, so
    against the data in hand today the ablation moves nothing at all. This is a prediction
    about one 96-well plate that does not exist yet, and the plate is cheap -- acetic acid
    0/20/40/80/160 mM crossed with pH 4.0/4.5/5.0/5.5 in the existing background, OD600 only.
    The block predicts that wells sharing an UNDISSOCIATED concentration collapse onto one
    growth curve while wells sharing a TOTAL concentration do not; a single EC50 on dose
    predicts the opposite. Nothing else in the repository can distinguish those two.
    """
    v = (min(CYTOSOLIC_VOLUMES_ML_PER_GDCW) if cytosolic_volume_ml_per_gdcw is None
         else CYTOSOLIC_VOLUME_ML_PER_GDCW.at(cytosolic_volume_ml_per_gdcw))
    # Imported from the window rather than retyped: a hardcoded band goes stale silently.
    mu_ref = (window.growth_rate_low_per_h if reference_growth_per_h is None
              else float(reference_growth_per_h))
    doses = (0.0, 20.0, 40.0, 80.0, 160.0)
    phs = (4.0, 4.5, 5.0, 5.5)
    data = pd.DataFrame([{"total_mM": d, "ph_ex": p, "growth_rate_per_h": mu_ref}
                         for d in doses for p in phs])

    full = _SeparabilityGrowth(name="partitioned dose", partition=True,
                               cytosolic_volume_ml_per_gdcw=v,
                               reference_growth_per_h=mu_ref, regime=regime)
    reduced = _SeparabilityGrowth(name="total dose only", partition=False,
                                  cytosolic_volume_ml_per_gdcw=v,
                                  reference_growth_per_h=mu_ref, regime=regime)

    def spread(model: FittableModel, fit: Fit) -> float:
        return abs(model.mu(40.0, 5.5) - model.mu(40.0, 4.0))

    observable = Observable(
        name="growth_rate", units="1/h",
        assay=GROWTH_RATE_FLOOR.assay, scored="absolute", data=data, summarise=spread,
        description=("spread of predicted specific growth rate across medium pH 4.0-5.5 at "
                     "a fixed 40 mM total acetate"))
    return ablate(full, reduced, observable, floor=GROWTH_RATE_FLOOR)


@dataclass(frozen=True)
class _QuenchedReporter(FittableModel):
    """Dilution-corrected reporter activity, with or without the pH quench correction."""

    name: str = "quenched reporter"
    quench: bool = True
    beta_mM_per_ph: float = 200.0
    response_column: str = "activity"
    parameter_names: tuple[str, ...] = ()

    def refit(self, data: pd.DataFrame) -> Fit:
        predicted = self.predict(data, {})
        observed = np.asarray(data[self.response_column], dtype=float)
        return Fit(model=self.name, parameters={}, n_free=0, n_rows=len(data),
                   residual_rms=float(np.sqrt(np.mean((predicted - observed) ** 2))),
                   converged=True)

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        return np.asarray([self.activity(float(r["total_mM"]), float(r["ph_ex"]))
                           for _, r in data.iterrows()], dtype=float)

    def activity(self, total_mM: float, ph_ex: float) -> float:
        if not self.quench:
            return 1.0
        ah_o = undissociated_outside(total_mM, ph_ex, float(PKA_ACETIC))
        bound = no_pump_ph(ah_o, self.beta_mM_per_ph, float(PKA_ACETIC))
        return citrine_quench(bound.ph_c)


def reporter_quench_ablation(beta_mM_per_ph: float | None = None) -> AblationResult:
    """B3, scored on the instrument: pH corrupts the channel that would measure the stress.

    Full model applies ``citrine_ph_response`` at the pH_c the block computes; reduced model
    leaves the reporter alone, which is what both data-generating paths do today. Scored
    relatively against ``OBSERVED_ACTIVITY_CV`` = 0.146, the plate-to-plate CV of
    dilution-corrected activity across 28 matched conditions on two real plates.

    Evaluated at the no-pump bound, and the default ``beta`` is the end of the sweep that
    gives the SMALLEST effect: 3.38x the floor at 200 mM/pH, rising to 4.45x at 90.

    IT SCORES THE ALGEBRA, NOT THE STATE, AND THE TWO NO LONGER DIFFER. This used to read
    that they differ by more than the floor, because :func:`reduction_audit` KEPT ``pH_c`` on
    the 4.14 h plate at a bias of 53-74% and ``no_pump_ph`` was the reduction that audit
    refused. Both of those were the charge-balance defect: repaired, the audit ELIMINATES
    ``pH_c`` at a bias of 0.04-0.07%, and :func:`reporter_state_ablation` measures the
    difference at 0.0051x the floor. So this ablation scores a constant, like the other two.

    AND AGAINST THE ACTUAL INCUMBENT IT FAILS. The reduced model here is frozen at 1.0 at
    every dose, which is not what this repository would fit; the incumbent for a dose
    response is a fitted dose response. :func:`quench_identifiability` refits a
    two-parameter Hill in TOTAL dose to the quench curve itself and measures the residual:
    **0.0104-0.0350x the floor at any FIXED medium pH** over a 300-fold dose range -- which
    is every plate this repository owns -- so at equal parameter count the artefact is
    absorbed by the EC50 and is not identifiable there. It separates only ACROSS medium pH,
    where the same fit leaves **1.33x** rather than 3.38x. The 3.38x below is therefore the
    effect against a null nobody would fit, and the number that travels is 1.33x, on a
    medium-pH-crossed plate that does not exist yet. The module docstring already says the
    quench "enters every dose-response fit as a slope rather than as noise"; this is that
    sentence with a number on it.

    The bracket matters here more than anywhere else in the module. At the pump-holds bound
    pH_c does not move, the quench factor is 1, and this effect is exactly zero. So the
    reporter artefact is not a prediction of the block -- it is a prediction of one END of
    the block, and which end the cell is at is the measurement :data:`NOT_BUILT` names.
    """
    beta = (max(CYTOSOLIC_BUFFER_CAPACITY.bounds) if beta_mM_per_ph is None
            else CYTOSOLIC_BUFFER_CAPACITY.at(beta_mM_per_ph))
    data = pd.DataFrame([{"total_mM": d, "ph_ex": p, "activity": 1.0}
                         for d in (0.0, 20.0, 40.0) for p in (4.0, 4.5, 5.0, 5.5)])
    full = _QuenchedReporter(name="pH-quenched Citrine", quench=True, beta_mM_per_ph=beta)
    reduced = _QuenchedReporter(name="uncorrected Citrine", quench=False,
                                beta_mM_per_ph=beta)
    observable = Observable(
        name="reporter_activity", units="fraction of unquenched signal",
        assay=REPORTER_ACTIVITY_FLOOR.assay, scored="relative", data=data,
        summarise=lambda model, fit: model.activity(40.0, 4.0),
        description=("Citrine signal surviving the cytosolic acidification caused by 40 mM "
                     "acetate at medium pH 4.0, at the no-pump bound"))
    return ablate(full, reduced, observable, floor=REPORTER_ACTIVITY_FLOOR)


@dataclass(frozen=True)
class _StatefulReporter(FittableModel):
    """Citrine signal at the end of a plate read, with ``pH_c`` integrated or reduced."""

    name: str = "stateful reporter"
    integrate: bool = True
    beta_mM_per_ph: float = 200.0
    window: Window = PLATE_READ_4H
    response_column: str = "activity"
    parameter_names: tuple[str, ...] = ()

    def refit(self, data: pd.DataFrame) -> Fit:
        predicted = self.predict(data, {})
        observed = np.asarray(data[self.response_column], dtype=float)
        return Fit(model=self.name, parameters={}, n_free=0, n_rows=len(data),
                   residual_rms=float(np.sqrt(np.mean((predicted - observed) ** 2))),
                   converged=True)

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        return np.asarray([self.activity(float(r["total_mM"]), float(r["ph_ex"]))
                           for _, r in data.iterrows()], dtype=float)

    def activity(self, total_mM: float, ph_ex: float) -> float:
        """Fraction of signal surviving at the end of the window."""
        pka = float(PKA_ACETIC)
        if total_mM <= 0.0:
            return 1.0
        if not self.integrate:
            ah_o = undissociated_outside(total_mM, ph_ex, pka)
            return citrine_quench(no_pump_ph(ah_o, self.beta_mM_per_ph, pka).ph_c)
        trace = simulate_weak_acid(self.window, total_acid_mM=total_mM, ph_ex=ph_ex,
                                   beta_mM_per_ph=self.beta_mM_per_ph, n_points=201)
        return citrine_quench(float(trace.of("pH_c")[-1]))


def reporter_state_ablation(beta_mM_per_ph: float | None = None) -> AblationResult:
    """RETRACTED, AND THE RETRACTION IS THE RESULT: pH_c is not a state, so removing it costs
    0.0051x the reporter floor and this ablation FAILS criterion (a).

    It is kept, failing, rather than deleted, because what it measures is the repair. The
    design is unchanged: the full model integrates the two-state system across
    ``PLATE_READ_4H`` and reads the Citrine signal at 4.14 h; the reduced model replaces
    ``pH_c`` by :func:`no_pump_ph`. Both sides are parameter-free, so no fit can be blamed
    for the collapse.

    WHAT IT USED TO SAY AND WHY THAT WAS WRONG. It reported 0.342 = **2.34x the floor at
    beta = 200 mM/pH and 2.86x at 90**, and called that criterion (a)'s answer for the ODE.
    Every bit of that gap was :func:`weak_acid_rhs`'s missing dilution term: the integrated
    pH_c was drifting away from charge balance in proportion to mu, so the two arms were
    being compared at 5.37 against 5.68 pH units. With the term in place they agree, and the
    measured effect is **0.000746 = 0.0051x the floor at beta = 200 and 0.000589 = 0.0040x
    at 90** -- two to three orders of magnitude below the plate's own 14.6% CV.

    WHAT IS LEFT IS NOT THE STATE EITHER, AND IT IS NOT NUMERICAL NOISE. The reduced arm uses
    :func:`no_pump_ph`, which puts ``A_i`` at partition equilibrium with the medium, while
    growth holds it a hair below -- ``J_in = mu A_i > 0`` requires ``[AH]_i < [AH]_o``. That
    offset is 3.78e-4 of the signal at beta 200 and 2.06e-4 at 90, it does not move with the
    integrator tolerance, and it is the whole remaining effect. Slaving pH_c to ``A_i``
    exactly instead reproduces the integrated endpoint to 4.4e-12 at rtol 1e-12. So the
    residual measures the anion pool's own dilution offset, not a pH state.

    So this block has NO state whose removal is measurable on an instrument it owns, and
    :data:`NOT_BUILT` already names what would change that -- ``vmax_anion_export``. Give the
    cell a pump and pH_c stops being slaved to ``A_i``, because export breaks the one-to-one
    between anion made and proton delivered that makes charge balance an algebraic constraint
    here.

    The bracket caveat still applies and now cuts the same way: at the pump-holds end
    ``pH_c`` does not move at all, both models return 1, and the ablation is undefined rather
    than small.
    """
    beta = (max(CYTOSOLIC_BUFFER_CAPACITY.bounds) if beta_mM_per_ph is None
            else CYTOSOLIC_BUFFER_CAPACITY.at(beta_mM_per_ph))
    data = pd.DataFrame([{"total_mM": d, "ph_ex": p, "activity": 1.0}
                         for d in (0.0, 40.0) for p in (4.0, 4.5)])
    full = _StatefulReporter(name="integrated pH_c", integrate=True, beta_mM_per_ph=beta)
    reduced = _StatefulReporter(name="pH_c reduced to algebra", integrate=False,
                                beta_mM_per_ph=beta)
    observable = Observable(
        name="reporter_activity", units="fraction of unquenched signal",
        assay=REPORTER_ACTIVITY_FLOOR.assay, scored="relative", data=data,
        summarise=lambda model, fit: model.activity(40.0, 4.0),
        description=("Citrine signal at the end of a 4.14 h plate read under 40 mM acetate "
                     "at medium pH 4.0, with pH_c integrated versus reduced to algebra"))
    return ablate(full, reduced, observable, floor=REPORTER_ACTIVITY_FLOOR)


def quench_identifiability(beta_mM_per_ph: float | None = None,
                           ph_range: tuple[float, ...] = (4.0, 4.5, 5.0, 5.5),
                           n_doses: int = 40) -> tuple[float, float, float]:
    """What a FITTED incumbent does to :func:`reporter_quench_ablation`'s effect, computed.

    ``reporter_quench_ablation`` scores against a reduced model frozen at 1.0 at every dose,
    which is not the incumbent: the incumbent for a dose response is a fitted dose response.
    This refits a two-parameter Hill in TOTAL dose to the quench curve itself, over a 40-point
    log grid spanning 0.1-316 mM, and returns ``(worst per-pH, best per-pH, pooled)`` maximum
    residuals as multiples of ``OBSERVED_ACTIVITY_CV``.

    THE ANSWER, at beta = 200 mM/pH: **0.0104-0.0350 at any single medium pH, and 1.33 when
    the four medium pHs are pooled**. Every plate in this repository is at ONE medium pH, so
    at equal parameter count the artefact is absorbed by the fitted EC50 and is not
    identifiable there. It separates only across medium pH, and that plate does not exist yet.

    Args:
        beta_mM_per_ph: A point on :data:`CYTOSOLIC_BUFFER_CAPACITY`. Defaults to the end of
            the sweep giving the smallest effect, as the ablation does.
        ph_range: Medium pHs to fit at, and to pool across.
        n_doses: Points on the log dose grid.
    """
    beta = (max(CYTOSOLIC_BUFFER_CAPACITY.bounds) if beta_mM_per_ph is None
            else CYTOSOLIC_BUFFER_CAPACITY.at(beta_mM_per_ph))
    pka = float(PKA_ACETIC)
    doses = np.logspace(-1.0, 2.5, int(n_doses))

    def quench(ph_ex: float) -> np.ndarray:
        return np.asarray([citrine_quench(no_pump_ph(
            undissociated_outside(d, ph_ex, pka), beta, pka).ph_c) for d in doses])

    def hill(dose, ec50, n):
        return 1.0 / (1.0 + (dose / ec50) ** n)

    def worst_residual(x, y):
        fitted, _ = curve_fit(hill, x, y, p0=[40.0, 1.0], maxfev=400000)
        return float(np.max(np.abs(hill(x, *fitted) - y)))

    per_ph = [worst_residual(doses, quench(p)) / float(REPORTER_ACTIVITY_FLOOR.noise_floor)
              for p in ph_range]
    pooled = worst_residual(np.concatenate([doses] * len(ph_range)),
                            np.concatenate([quench(p) for p in ph_range]))
    return (min(per_ph), max(per_ph),
            pooled / float(REPORTER_ACTIVITY_FLOOR.noise_floor))


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------

def provenance_summary() -> str:
    """The gate, the provenance table and the refusals, as printable markdown."""
    lines = [PH_GATE.summary(), "",
             f"s_ATP is taken per regime from bridge/maintenance_calibration.py: "
             f"{REFERENCE_REGIME.summary()}", "",
             provenance_table(PH_PARAMS, NOT_BUILT, heading_level=3)]
    return "\n".join(lines)
