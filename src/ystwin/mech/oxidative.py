"""Oxidative / Yap1, written in the ONE variable whose entry rate has actually been measured.

WHAT MAKES THIS BLOCK DIFFERENT FROM THE OTHER NINE. Every other block in
`ARCHITECTURE_TARGET.md` L4 needs a sensor input function ``u_j(dose)`` that nobody has
measured, and writes one anyway. This block does not have to: H2O2 entry is physical, it is
first-order, and Branco, Marinho, Cyrne & Antunes 2004 (J Biol Chem 279:6501-6, PMID
14645222) measured the whole-culture consumption rate directly in *S. cerevisiae* with a
bolus -- the delivery mode this project actually uses. Selvaggio, Coelho & Salvador 2018
(Redox Biol 15:297-315, PMID 29304480) SI eq.919 refits that figure to
``(1.33 +/- 0.08)e-3 /s`` at ``7e6 cells/mL``. One number, right organism, right delivery.

WRITING THE MODEL IN ``E`` IS WHAT BUYS THAT. The alternative -- a cytosolic H2O2 state fed
through a membrane permeability -- requires a constant whose three published values span
**570x** for the same physical quantity: Goulev 2017 ``alpha = 1 /min`` (ASSERTED, from a
signalling latency), Tomalin 2016 ``0.331 /s`` (FITTED, and in *S. pombe*), Selvaggio 2018
``9.5 /s`` (DERIVED from Branco). A model written in the extracellular pool never forms that
constant at all. :data:`NOT_BUILT` records it as REFUSED with all three values attached.

THE THREE STATES, AND WHY THEY ARE THESE THREE. The build list says three states and the
harvest verifier says which three:

* ``h2o2_extracellular`` -- the well is not a reservoir. This is the piece the repository is
  missing everywhere today: a dose of 0.15 mM is currently 0.15 mM at hour eight.
* ``oxidative_mrna`` -- the transcriptional delay, MEASURED at ln2/40 min.
* ``oxidative_effector`` -- the regulon protein. It is the reporter on the plate and it is
  the carrier of cross-protection in the cell, and those are the same pool.

The two states an earlier draft of the architecture named -- cytosolic H2O2 and the Yap1
nuclear fraction -- are here as ALGEBRA, not because they are unimportant but because
criterion (d) forbids them a differential equation. Yap1 nuclear entry is 120 s (Goulev
2017, PMID 28418333, 30 s sampling) and membrane equilibration is 0.105-60 s -- the SLOWEST
of the three disputed permeabilities, so the reduction does not depend on which is right.
Against a 4.14 h read that is ``tau/T`` = 8.05e-3 and at most 4.03e-3, both far under the
0.146 ceiling. :func:`yap1_nuclear_fraction` is that algebra, named, with its measured time
constant on the record; :func:`cytosolic_gradient` REFUSES, because the GAIN it would need
-- as opposed to the speed -- is the 570x constant.

THE FIVE-STATE PEROXIREDOXIN CHEMISTRY IS REFUSED, NOT DEFERRED. `REVISED_BUILD_LIST.md` is
explicit and the arithmetic is in :data:`NOT_BUILT`: Selvaggio's Model 1 collapses to a
linear gain above ~10 uM external, and swapping the whole Tsa1 abundance 3.6x (Ghaemmaghami
2003 against the Ho 2018 median) moves its only continuous output by 3%. Nine cited rate
constants that change no observable this project can measure.

WHAT THIS BLOCK DOES AND DOES NOT REPRODUCE. Every number below was produced by running the
functions named beside it, and the failures are here for the same reason the passes are.

* **Criterion (a) PASSES, at 4.73-7.04x the floor**, on the inoculum-density prediction.
  Both models are refit on the same committed rows with the same three observation
  nuisances and then asked for a well that was never run: the endpoint RFU at OD600 0.05
  over the endpoint at OD600 0.80, same pipetted dose. A dose held forever has no density
  to depend on, so the incumbent answers exactly 1.0; the depleting pool answers 1.98-2.45
  across the three committed replicates and both ends of the swept cells/OD600 band. The
  floor is 0.146 propagated through a ratio to **0.2065**, because a ratio of two readings
  carries both readings' noise and borrowing the single-reading CV would be optimistic by
  41%. :func:`inoculum_spread_ablation`. Scored on a single inoculum against the plain
  0.146 instead (:func:`inoculum_ablation`) it still clears, at 1.77-3.07x.
  **It is a prediction and not a validation**: all three plates ran at one inoculum, and
  the discriminating experiment -- the same ladder at four inocula -- has not been run.
* **Criterion (c) EMERGENCE, and this is the strongest result here.** Let the plates pick
  the decay rate themselves, told nothing about Branco, and the three committed replicates
  choose **1.966, 0.882 and 0.669 /h** against the band ``k_ref`` predicts at their own
  inoculum, **0.876-3.283 /h**. Two land inside it and the third is 1.31x below its floor,
  and the fit improves 1.66-2.16x in RMS over no decay at all on every plate. A 2004
  measurement, a different lab, a different assay, 2026 traces. :func:`fitted_decay_rates`.
* **Criterion (a) does NOT clear on curve shape.** Leave-one-dose-out on the same ladders
  has the depleting model winning on 3 of 3 plates and 11 of 12 folds, by +0.037 to +0.070
  relative error = **0.25-0.48x the 0.146 floor**. A consistent win, below the noise.
  :func:`shape_comparison`. Reporting only the passes would be the exact failure
  `mech/ablation.py` was written to prevent.
* **One truncation decides both of the last two, and it is MEASURED, not chosen.** This
  block takes growth as an input and has no growth-inhibition term, so a well whose growth
  is more than halved is a well whose dilution correction it does not describe. Truncating
  each ladder at the construct's own measured growth-halving dose (0.79 mM on NativeYap1)
  drops the 1.0 mM rung. Keep that rung and the same fits prefer 0.635, 0.037 and
  0.000 /h -- two of three plates choose no decay at all. :data:`GROWTH_HALVING_MM`.
* **The BOLUS/PERFUSION GAP does NOT fall out, and that is the sharpest negative result
  here.** The repository attributes the gap -- Goulev's 0.6 mM perfusion arrest against a
  bolus that halves growth at 0.79 and 1.24 mM -- to Tomalin 2016's thiol buffering, in
  *S. pombe*, with an assigned mammalian thiol pool. A depleting well predicts the gap on a
  time-integrated-exposure basis at **4.48-14.97x** across the swept conversion and this
  window's own measured growth band -- 4.70-14.82 at the plates' median mu -- against a
  measured 1.32-2.07x: right sign, **7.2x too large** at the plate's own density. Tomalin's
  sink would be an ADDITIONAL consumption term, so it makes the overshoot worse rather than
  better. Both published explanations over-predict, in the same direction. A peak-concentration
  criterion predicts exactly 1.0x and misses by 1.32x, 5.4x closer than exposure.
  :func:`bolus_perfusion_gap`.
* **Cross-protection is carried with NO memory parameter.** The effector's only loss is
  dilution plus a MEASURED turnover, so its departure from a pure generation counter is
  ``k_deg(Ctt1)/mu`` = **1.38-2.14%** across Guan's own doubling times, and it outlives the
  transcript that made it by **3.21-4.41x** in generations, threshold-free. Guan 2012
  reported persistence in GENERATIONS (4-5), transmission to daughters that never saw the
  pre-treatment, and no requirement for new protein synthesis. The absolute 4-5 is NOT
  predicted and is not claimed: it needs a detection threshold on the protection, and this
  block refuses to supply one. :func:`cross_protection`.

* **The block's ONE structural assumption is declared, and the plates cannot test it.**
  ``k`` proportional to ``X`` to the first power is anchored at Branco's single density.
  The plates grow 3.8x inside a read, so the exponent is in principle visible without a
  second inoculum -- and it is not: freeing it buys 1.000, 1.021 and 1.124x in RMS, all
  under the 0.146 floor, while the exponent lands on 1.052, 2.916 and 6.000, the last at
  the search bound. Three replicates of one design, disagreeing 6x about a quantity assumed
  to be 1. UNTESTED, not confirmed; :func:`density_exponent_identifiability`, and
  :data:`NOT_BUILT`'s ``density_exponent``. **It needs the same four-inoculum plate the
  headline needs**, which is the strongest argument in this file for running that plate.

ONE MEASURED CONTRADICTION THIS BLOCK CANNOT RESOLVE, carried rather than hidden. Guan 2012
attributes acquired tolerance to catalase Ctt1p, i.e. to removal capacity. Branco 2004's own
result (c) is that there is **no** correlation between the intracellular capacity to remove
H2O2 and resistance to it, and its result (d) is that resistance is acquired by halving the
plasma-membrane permeability -- and is acquired by cycloheximide too, so without new protein.
Two papers, both *S. cerevisiae*, both MEASURED, pointing at different mechanisms. The map
from the effector pool to protection is therefore :data:`NOT_BUILT`, and what this block
demonstrates is that the CARRIER persists, not how much it protects.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .ablation import (
    REPORTER_ACTIVITY_FLOOR,
    AblationResult,
    Fit,
    FittableModel,
    NestedComparison,
    Observable,
    ablate,
    nested_out_of_sample,
    reporter_block,
)
from .integrate import ReducedSystem, Trajectory, integrate_window, reduce_for_window
from .params import Param, ParamRegistry, Target, provenance_table
from .state import PLATE_READ_4H, Encoding, MechState, StateVar, Window

__all__ = [
    "ADAPTED_PERMEABILITY_FOLD",
    "CELLS_PER_ML_PER_OD600",
    "CROSS_PROTECTION_FLOOR",
    "DRIVER_FLOOR",
    "CrossProtection",
    "BolusPerfusionGap",
    "GOULEV_PERFUSION_ARREST_MM",
    "GUAN_GROWTH_BAND_PER_H",
    "K_EX_MM",
    "K_REF_PER_H",
    "LETHAL_RATIO_FLOOR",
    "MRNA_DECAY_PER_H",
    "NOT_BUILT",
    "OXIDATIVE_GATE",
    "OXIDATIVE_PARAMS",
    "PLATE_INOCULUM_OD600",
    "PLATE_LETHAL_BOLUS_MM",
    "PLATE_MEDIAN_GROWTH_PER_H",
    "PROTEIN_TURNOVER_CTT1_PER_H",
    "PROTEIN_TURNOVER_TSA1_PER_H",
    "REPORTER_TARGET",
    "SPREAD_RATIO_FLOOR",
    "X_REF_CELLS_PER_ML",
    "YapReporter",
    "bolus_perfusion_gap",
    "GROWTH_HALVING_ALTERED_YAP1_MM",
    "GROWTH_HALVING_MM",
    "GROWTH_HALVING_NATIVE_YAP1_MM",
    "committed_yap1_blocks",
    "cross_protection",
    "cytosolic_gradient",
    "density_exponent_identifiability",
    "dose_decay_band_per_h",
    "dose_decay_rate_per_h",
    "dose_decay_from_density_per_h",
    "dose_half_life_h",
    "driver_horizon_h",
    "effector_trace",
    "extracellular_h2o2",
    "fitted_decay_rates",
    "gem_peroxide_uptake_mmol_per_gdcw_h",
    "inoculum_ablation",
    "inoculum_endpoint_observable",
    "inoculum_spread_ablation",
    "inoculum_spread_observable",
    "memory_persistence_generations",
    "oxidative_rhs",
    "oxidative_states",
    "provenance_summary",
    "reduction_audit",
    "shape_comparison",
    "simulate_oxidative",
    "time_integrated_dose",
    "transfer_comparison",
    "yap1_frame",
    "yap1_nuclear_fraction",
]

_LN2 = math.log(2.0)


class GradientRefused(NotImplementedError):
    """The cytosolic/extracellular H2O2 gradient was asked for and has no defensible value.

    ``NotImplementedError`` to match `params.NoSingleValue` and
    `mech/ph.py::AcidUnmeasured`: the question is well-formed and the literature answers it
    three times, 570x apart.
    """


# --------------------------------------------------------------------------------------
# The constants, graded. Two are free; everything the block actually runs on is pinned.
# --------------------------------------------------------------------------------------

OXIDATIVE_PARAMS = ParamRegistry("mech/oxidative.py::yap1_three_state")
"""The built block's constants. Two free scalars against three independent targets --
:data:`OXIDATIVE_GATE` computes criterion (e) at import and refuses if it stops being true."""

NOT_BUILT = ParamRegistry("mech/oxidative.py::not-built")
"""Every constant the sources name that no equation here consumes, and why.

Outside the gate on `mech/ph.py::NOT_BUILT`'s precedent: a refusal that nothing reads is a
record, not a degree of freedom. Nine of these are the Selvaggio peroxireduction cycle that
`REVISED_BUILD_LIST.md` instructs this block not to build."""

_BRANCO = ("Branco, Marinho, Cyrne & Antunes 2004, J Biol Chem 279:6501-6, PMID 14645222")
_SELVAGGIO = ("Selvaggio, Coelho & Salvador 2018, Redox Biol 15:297-315, PMID 29304480")
_GOULEV = "Goulev et al. 2017, eLife 6:e23971, PMID 28418333"
_CHRISTIANO = "Christiano et al. 2014, Cell Rep 9:1959-65, PMID 25466257, via SGD"
_GUAN = "Guan, Haroon, Bravo, Will & Gasch 2012, Genetics 192:495-505, PMID 22851651"

K_REF_PER_H = OXIDATIVE_PARAMS.add(Param.measured(
    "k_ref", 1.33e-3 * 3600.0, "1/h at X_ref",
    "MEASURED in the host and in the right delivery mode. " + _BRANCO + " Figure 1 followed a "
    "100 uM bolus on a 7e6 cells/mL culture; " + _SELVAGGIO + " SI eq.919 refits it to "
    "k_cells = (1.33 +/- 0.08)e-3 /s, which is 4.788 /h. Its own verification pass reproduced "
    "every step of that arithmetic. TWO CORRECTIONS TO ARCHITECTURE_TARGET.md L1, both "
    "checked here: (i) it tabulates 1.33e-3 /s as '0.0798 h^-1', and 0.0798 is the per-MINUTE "
    "value -- a 60x error in the table cell. It did NOT propagate: the document's own time "
    "constants (36 min at OD 0.05, 2.5 min at OD 0.80) are computed from the correct 4.788 /h. "
    "(ii) Those same time constants imply 4.4-4.9e7 cells/mL per OD600, which is ABOVE the "
    "BNID band this module sweeps, so the document's taus are 1.46x and 1.62x faster than "
    "the exponential end of the bracket gives -- 2.5 against 3.7 min and 36 against 58.5. "
    "This module uses 4.788 /h",
    ci95=(1.25e-3 * 3600.0, 1.41e-3 * 3600.0)))
"""4.788 /h at 7e6 cells/mL. The one measured number the whole block turns on."""

X_REF_CELLS_PER_ML = OXIDATIVE_PARAMS.add(Param.measured(
    "X_ref", 7.0e6, "cells/mL",
    "MEASURED. The density of the culture " + _BRANCO + " consumed its bolus in, as reported "
    "in " + _SELVAGGIO + " SI eq.919. It is a denominator, not a rate: k_ref is defined "
    "AT this density and the model scales linearly away from it"))

CELLS_PER_ML_PER_OD600 = OXIDATIVE_PARAMS.add(Param.swept(
    "cells_per_ml_per_od600", "cells/mL per OD600 unit",
    "SWEPT over a bracket both of whose ends are in one curated external database, the "
    "vendored BioNumbers Nov-2024 export: BNID 100986 gives 3.0e7 cells/mL at OD600 = 1 for "
    "budding yeast, and BNID 106301 gives 8.0e6 for a stationary-phase culture at 30 C. A "
    "3.75x band, and it is the only conversion standing between k_ref's cells/mL and this "
    "project's OD600. It is not collapsible by argument: the plates start exponential and "
    "end approaching stationary, so both ends of the band are visited inside one read",
    bounds=(8.0e6, 3.0e7),
    missing="a haemocytometer or Coulter count against OD600 on this project's own strain, "
            "medium and reader. `calib/od.py` already refuses the neighbouring gdcw_per_od "
            "conversion for the same reason"))
"""FREE #1. Every reported result here is a band across it, never a point."""

K_EX_MM = OXIDATIVE_PARAMS.add(Param.bounded(
    "K_ex", 0.15, "mM extracellular H2O2",
    "BOUNDED. " + _GOULEV + " measured Yap1-GFP nuclear entry partial at 0.1 mM and saturated "
    "above 0.2 mM in S288C under SUSTAINED microfluidic delivery, which is the one delivery "
    "mode in which a measured half-max IS a half-saturation on the instantaneous external "
    "concentration -- exactly the quantity this model needs. 0.15 mM is the value "
    "generator/stress_panel.py::STRESSORS['H2O2'] already carries for the same reason. It is "
    "BOUNDED and not MEASURED because Goulev measured the SENSOR and this constant drives "
    "TRANSCRIPTION, one unmeasured transfer function further on",
    bounds=(0.02, 0.6),
    missing="a YRE-reporter dose series under sustained delivery, which would measure the "
            "transcriptional half-saturation directly. The bracket ends are both measured, "
            "and both on this constant's own EXTERNAL axis: Kritsiligkou 2021, PMID 34118234 "
            "reports verbatim that 'the lowest concentration triggering a detectable response "
            "was ... 20 uM for HyPer7' for EXOGENOUSLY ADDED H2O2 on BY4742, and says in the "
            "same paragraph that the matching CYTOSOLIC concentration is 'in the low nanomolar "
            "range'. So 20 uM is a floor on the EXTERNAL dose that does anything at all, which "
            "is the axis this bracket is on, and it is emphatically not a cytosolic 20 uM. "
            "0.6 mM is where Goulev sees complete growth arrest under perfusion, by which "
            "point the regulon is saturated"))
"""FREE #2, and it is the block's only shape parameter. Never fitted here."""

MRNA_DECAY_PER_H = OXIDATIVE_PARAMS.add(Param.measured(
    "k_deg_mrna", _LN2 / (40.0 / 60.0), "1/h",
    "MEASURED. " + _GOULEV + " Methods tabulates the antioxidant mRNA decay rate as "
    "ln2/40 min = 0.0173 /min and attributes it to Geisberg, Moqtaderi, Fan, Ozsolak & "
    "Struhl 2014, Cell 156:812-24, PMID 24529382, mRNA half-lives in yeast. Its own "
    "verification pass confirmed the Goulev table verbatim. It is the antioxidant-regulon "
    "value Goulev's own model runs on and NOT a TSA1-specific measurement -- Goulev's species "
    "is a generic antioxidant transcript and neither paper reports TSA1 by name -- and it is "
    "why this block does not inherit mech/state.py's global `mrna` row, whose 3.6-32 min "
    "methodological span STRADDLES the reduction threshold"))
"""1.0397 /h. With the plate's measured mu it puts the transcript pool just INSIDE the
irreducible band, which is the whole reason the mRNA state is here.

THE NARROWEST CALL IN THE BLOCK, AND ITS SENSITIVITY IS COMPUTED RATHER THAN HOPED. Solving
``1/(k_deg + mu) = 0.146 T`` on PLATE_READ_4H gives an mRNA half-life of 29.5 min at
mu = 0.245 /h and 32.2 min at 0.363 /h: at any antioxidant transcript half-life BELOW that,
this state eliminates and the block is a two-state block. 40 min clears the ceiling by
1.18-1.29x and nothing else in the repository's mRNA range would. Measuring the TSA1 or CTT1
transcript's own half-life on this strain is therefore the cheapest experiment that could
remove a state from here."""

PROTEIN_TURNOVER_TSA1_PER_H = OXIDATIVE_PARAMS.add(Param.measured(
    "k_deg_tsa1", _LN2 / 10.9, "1/h",
    "MEASURED. " + _CHRISTIANO + ": Tsa1 protein half-life 10.9 h. Already dilution-corrected "
    "at source, so `mu + k_deg` composes correctly rather than double-counting"))

PROTEIN_TURNOVER_CTT1_PER_H = OXIDATIVE_PARAMS.add(Param.measured(
    "k_deg_ctt1", _LN2 / 108.8, "1/h",
    "MEASURED. " + _CHRISTIANO + ": Ctt1 protein half-life 108.8 h. It is 9.98x Tsa1's and 26x "
    "the whole plate read, which is why " + _GUAN + " could find the memory of a stress "
    "carried by this protein and not by new synthesis"))

ADAPTED_PERMEABILITY_FOLD = OXIDATIVE_PARAMS.add(Param.measured(
    "adapted_permeability_fold", 2.0, "dimensionless",
    "MEASURED, verbatim from the " + _BRANCO + " abstract: 'the plasma membrane permeability "
    "to H2O2 decreases by a factor of two upon acquisition of resistance to this agent by "
    "pre-exposing cells either to nonlethal doses of H2O2 or to cycloheximide'. Carried "
    "because it is the measured SIZE of the adaptation, and because the cycloheximide arm "
    "says the adaptation does not need new protein -- which is in tension with the Ctt1p "
    "route " + _GUAN + " measured. Nothing in this module consumes it; see NOT_BUILT's "
    "`effector_to_protection_map`"))

PLATE_MEDIAN_GROWTH_PER_H = 0.3227
"""The plates' own median specific growth rate, 1/h. MEASURED, not asserted.

`mech/state.py::PLATE_GROWTH_BAND_PER_H` derives it as ``ln(od_fold)/4.14 h`` over the 183
gate-1 wells of the three committed NewProtocol reads: median 0.3227, IQR 0.2451-0.3625.
Every ``growth_rate_per_h`` default below is this number rather than a repeated literal, and
every band this module quotes is swept over the window's 0.245-0.363 /h instead."""

PLATE_INOCULUM_OD600 = OXIDATIVE_PARAMS.add(Param.measured(
    "plate_inoculum_od600", 0.16, "OD600 above blank",
    "MEASURED on this project's own plates. generator/plate.py::PlateConditions carries it as "
    "the blank-corrected starting OD of the NewProtocol reads, and generator/calibrate.py "
    "derives the same figure as the median first-row OD minus the blank on the committed "
    "exports. Cites a plate, not a PMID"))

GOULEV_PERFUSION_ARREST_MM = OXIDATIVE_PARAMS.add(Param.measured(
    "perfusion_arrest_mM", 0.6, "mM external H2O2",
    "MEASURED. " + _GOULEV + ": complete growth arrest and death at 0.6 mM under continuously "
    "replenished microfluidic delivery, i.e. at a HELD external concentration and at "
    "single-cell density"))

PLATE_LETHAL_BOLUS_MM = OXIDATIVE_PARAMS.add(Param.measured(
    "bolus_lethal_mM", 1.0, "mM external H2O2",
    "MEASURED on this project's own plates, and recorded inside "
    "generator/stress_panel.py::STRESSORS['H2O2']: 'Lethal dose 1.0 mM IS measured -- growth "
    "halves at 0.79 and 1.24 mM on NativeYap1 and AlteredYap1'. The 1.0 is the declared "
    "midpoint of those two constructs"))

_HALVING_SOURCE = (
    "MEASURED on this project's own plates and recorded verbatim inside "
    "generator/stress_panel.py::STRESSORS['H2O2']: 'growth halves at 0.79 and 1.24 mM on "
    "NativeYap1 and AlteredYap1'. Cites a plate, not a PMID")

GROWTH_HALVING_NATIVE_YAP1_MM = OXIDATIVE_PARAMS.add(Param.measured(
    "growth_halving_nativeyap1_mM", 0.79, "mM external H2O2", _HALVING_SOURCE))

GROWTH_HALVING_ALTERED_YAP1_MM = OXIDATIVE_PARAMS.add(Param.measured(
    "growth_halving_alteredyap1_mM", 1.24, "mM external H2O2", _HALVING_SOURCE))

GROWTH_HALVING_MM = {
    "NativeYap1": float(GROWTH_HALVING_NATIVE_YAP1_MM),
    "AlteredYap1": float(GROWTH_HALVING_ALTERED_YAP1_MM),
}
"""The dose above which a construct's own growth is more than halved.

Load-bearing, not bookkeeping. This block's growth rate is an INPUT measured from the OD
series, so above the halving dose the model is being scored on wells whose dilution
correction it does not describe -- and on this ladder that is the 1.0 mM rung, which is
above NativeYap1's measured 0.79 mM. Every ladder used here is truncated by this MEASURED
constant rather than by a chosen cutoff, and the difference it makes is the whole of
:func:`fitted_decay_rates`."""


# --------------------------------------------------------------------------------------
# Targets: the three independent measurements this block can be refuted by
# --------------------------------------------------------------------------------------

REPORTER_TARGET = OXIDATIVE_PARAMS.add_target(REPORTER_ACTIVITY_FLOOR)
"""0.146 CV, the plate reader. Scores :func:`inoculum_ablation` and :func:`shape_comparison`."""

SPREAD_RATIO_FLOOR = Target(
    name="reporter_activity_ratio",
    observable="ratio of two dilution-corrected endpoint activities",
    assay=REPORTER_ACTIVITY_FLOOR.assay,
    noise_floor=REPORTER_ACTIVITY_FLOOR.noise_floor * math.sqrt(2.0),
    units=REPORTER_ACTIVITY_FLOOR.units,
    source=("DERIVED from " + REPORTER_ACTIVITY_FLOOR.source + ". A RATIO of two endpoint "
            "activities carries both readings' noise, so its floor is sqrt(2) x 0.146 = "
            "0.2065 rather than 0.146. Using the single-reading CV on a ratio would be "
            "optimistic by 41%, which is exactly the kind of borrowing criterion (a) forbids "
            "-- so :func:`inoculum_spread_ablation` is scored against this and not against "
            "REPORTER_ACTIVITY_FLOOR"),
)
"""0.2065. NOT registered as a fourth target: it is the same measurement propagated through
a ratio, and one measurement counted twice is the arithmetic the free-scalar gate exists to
prevent."""

_LETHAL_ASSAY = ("H2O2 dose ladder to a growth-halving endpoint: 96-well bolus on this "
                 "project's plates against a microfluidic perfusion to full arrest")

LETHAL_RATIO_FLOOR = OXIDATIVE_PARAMS.add_target(Target(
    name="bolus_perfusion_lethal_ratio",
    observable="bolus growth-halving dose divided by perfusion full-arrest dose",
    assay=_LETHAL_ASSAY,
    noise_floor=(1.24 - 0.79) / 0.79,
    units="fraction",
    source=(
        "MEASURED INSIDE THE SAME EXPERIMENT, which is why it is not a borrowed floor. "
        "generator/stress_panel.py::STRESSORS['H2O2'] records growth halving at 0.79 mM on "
        "NativeYap1 and 1.24 mM on AlteredYap1 -- two constructs, one plate design, one "
        "endpoint, disagreeing by 0.5696 of the smaller. That is this assay's own "
        "demonstrated spread on exactly the quantity the ratio is built from, so it is the "
        "floor a predicted ratio has to beat. The perfusion arm is " + _GOULEV + "'s 0.6 mM"),
))
"""0.5696 relative. The block FAILS this target and the failure is the point of
:func:`bolus_perfusion_gap`."""

_GUAN_ASSAY = ("survival of a severe H2O2 challenge in cells released from a mild NaCl "
               "pre-treatment, read in generations of outgrowth")

CROSS_PROTECTION_FLOOR = OXIDATIVE_PARAMS.add_target(Target(
    name="cross_protection_generations",
    observable=("generations of outgrowth over which acquired H2O2 tolerance persists, "
                "relative to what the transcript that produced it could support"),
    assay=_GUAN_ASSAY,
    noise_floor=0.5 / 4.5,
    units="fraction",
    source=(
        _GUAN + " reports the persistence as 'four to five generations', so the measurement "
        "itself is quoted to +/- 0.5 of its 4.5 midpoint = 0.1111 relative. That interval is "
        "the source's own and is the only floor this observable has on any platform; this "
        "project has never run the assay"),
))
"""0.1111 relative. Guan's own quoted spread, used as the floor for criterion (c)."""

OXIDATIVE_GATE = OXIDATIVE_PARAMS.require_gate()
"""Criterion (e), computed at import: **2 free scalars against 3 independent targets**.

The unreduced block would have had many more. Everything dropped is in :data:`NOT_BUILT`,
and the count of what it holds is what the reduction bought."""


# --------------------------------------------------------------------------------------
# Refusals. Every one is a constant a named source carries that no equation here consumes.
# --------------------------------------------------------------------------------------

NOT_BUILT.add(Param.refused(
    "k_inf_membrane_permeability", "1/s",
    _GOULEV + " (alpha = 1 /min = 0.0167 /s); Tomalin et al. 2016, Free Radic Biol Med "
    "95:333-48, PMID 26944189 (k_perm = 0.331 /s); " + _SELVAGGIO + " SI eq.928 (9.5 /s)",
    reason="three published values for one physical quantity, spanning 570x. Goulev's is "
           "ASSERTED from a signalling latency ('the order of magnitude of timings (~1 min) "
           "of nuclear relocation'), Tomalin's is FITTED and in S. pombe, and Selvaggio's is "
           "DERIVED from Branco but does not reproduce the very experiment it came from. "
           "That last is this block's harvest verification pass reporting its own run of "
           "Selvaggio Model 1 at Branco's 100 uM / 7e6 cells/mL (harvest/verify__yap1-*.json, "
           "unit_problems[3]): a 718 s predicted extracellular half-life against Branco's "
           "measured 521 s = ln2/1.33e-3, because Tsa1 saturates, the gradient collapses to "
           "3.5-fold and back-flux eats 28% of the net efflux. Refitting inside the model "
           "moves k_inf to ~13 /s -- a 1.4x self-inconsistency on top of the 570x dispute",
    missing="a direct permeability determination in S. cerevisiae. Writing the model in the "
            "extracellular pool removes the need for it entirely, which is why this block "
            "has no cytosolic H2O2 state"))

NOT_BUILT.add(Param.refused(
    "k_alt_alternative_sinks", "1/s",
    _SELVAGGIO + " SI eq.969: kAlt = 21 (Ccp1) + 3.1 (Ctt1) + 9.5 (efflux) = 34 /s",
    reason="61% of that sum is the Ccp1 term, which assumes H2O2 has free access to the "
           "mitochondrial intermembrane space, and Zimmermann et al. 2025, Free Radic Biol "
           "Med 226:408-20, PMID 39515595 measured in vivo that Tsa1 alone accounts for "
           "almost all exogenous H2O2 scavenging while the heme peroxidases show a defect "
           "only when all three are deleted together. Selvaggio's Model 2 also double-counts "
           "the 9.5 /s efflux, which is inside kAlt and added again as an explicit "
           "permeation term",
    missing="a compartment-resolved scavenging budget. Branco's own culture-level rate "
            "constant already contains every sink, which is why this block uses it whole"))

for _name, _units, _src in (
    ("kOx", "1/M/s", "Tairum 2016 via " + _SELVAGGIO + " SI line 976, 4.7e7; Ogusucu et al. "
     "2007, Free Radic Biol Med 42:326-34, PMID 17210445 gives 2.2e7 by another method"),
    ("kCond", "1/s", _SELVAGGIO + " SI eq.997, 0.58 +/- 0.046, FITTED"),
    ("kSulf", "1/M/s", _SELVAGGIO + " SI eq.997, 7 +/- 3.4, FITTED at +/-49%"),
    ("kSrx", "1/s", _SELVAGGIO + " SI eq.1031-1037, 5.9e-5"),
    ("kRed", "1/M/s", _SELVAGGIO + " SI eq.1023, 1.2e6 as a Trx1/Trx2 weighted mean"),
    ("PrxT", "uM", _SELVAGGIO + " SI eq.973, 31 uM from Ghaemmaghami 2003"),
    ("TrxT", "uM", _SELVAGGIO + " SI eq.1045, 1.5 uM"),
    ("VMaxApp", "mM/s", _SELVAGGIO + " SI eq.1049-1054, 1.0 mM/s"),
    ("KM_trxr", "uM", _SELVAGGIO + " SI eq.1054, 0.8 uM"),
):
    NOT_BUILT.add(Param.refused(
        f"prx_cycle_{_name}", _units, _src,
        reason="the Selvaggio five-state peroxiredoxin cycle is REFUSED by "
               "REVISED_BUILD_LIST.md, and the arithmetic behind that instruction is that it "
               "collapses above ~10 uM external to a linear gain, and that a 3.6x swap of the "
               "whole Tsa1 abundance (Ghaemmaghami 2003 against the Ho 2018 median of 21 "
               "datasets, PMID 29361465) moves its only continuous output by 3%. Every "
               "constant is real, cited and in the right organism; not one changes an "
               "observable this project can measure",
        missing="HyPer7 for cytosolic H2O2, an anti-PrxSO3 western for the hyperoxidised "
                "pool and an AMS gel shift for the disulfide. This project owns none of "
                "them, and its two ratiometric reporters are a class its stated assay kit "
                "cannot read"))

NOT_BUILT.add(Param.refused(
    "orp1_yap1_relay_rates", "1/M/s",
    "Delaunay, Isnard & Toledano 2000, EMBO J 19:5157-66, PMID 11013218 and Delaunay, "
    "Pflieger, Barrault, Vinh & Toledano 2002, Cell 111:471-81, PMID 12437921 establish the "
    "mechanism and identify the Cys36 sulfenate; neither reports a rate constant",
    reason="not one rate constant of the actual Yap1 sensor is measured in S. cerevisiae: "
           "not Gpx3/Orp1 + H2O2, not the Orp1-SOH -> Yap1-Cys598 relay, not the competing "
           "intramolecular resolution whose RATIO with it sets the sensor gain, not "
           "thioredoxin + oxidised Yap1, and not Yap1 nuclear import or export separately. "
           "This block therefore has NO Yap1 state -- adding one would be inventing a "
           "transfer function, and the 120 s relocation is 21-23x faster than the transcript "
           "state it drives (tau = 1/(k_deg + mu), so the multiple carries the growth band), "
           "so it belongs in the algebraic driver anyway",
    missing="a stopped-flow determination of the Orp1 relay and its competing resolution, "
            "and separate import and export constants for reduced and oxidised Yap1"))

NOT_BUILT.add(Param.refused(
    "yap1_yre_transfer_function", "M and dimensionless",
    "no source. Morgan, Banks, Toone, Raitt, Kuge & Johnston 1997, EMBO J 16:1035-44, PMID "
    "9118942 shows TRX2 keeps roughly two thirds of its H2O2 induction in yap1, through Skn7",
    reason="the Yap1 Kd for the YRE, its Hill coefficient and its cooperativity with Skn7 at "
           "TRX2 are all unmeasured, so the promoter this project reads is a two-input gate "
           "with unmeasured weights. Skn7 is worse: no quantitative constant of any kind "
           "exists for its oxidative arm, and the well-parameterised Sln1-Ypd1 phosphorelay "
           "route is a dead end because Skn7's oxidative function is independent of the "
           "receiver aspartate D427. Both stay lumped inside K_ex, which is graded BOUNDED "
           "for exactly this reason",
    missing="a YRE occupancy measurement, and a yap1 / skn7 pair of dose series on this "
            "project's own TRX2 reporter to split the gate"))

NOT_BUILT.add(Param.refused(
    "effector_to_protection_map", "dimensionless per relative unit",
    _GUAN + " (tolerance carried by long-lived Ctt1p) against " + _BRANCO + " result (c)",
    reason="the two measurements disagree about the mechanism. Guan attributes acquired "
           "tolerance to catalase Ctt1p, i.e. to removal capacity; Branco measured that there "
           "is NO correlation between the intracellular capacity to remove H2O2 and "
           "resistance to it, and that resistance is acquired by halving the plasma-membrane "
           "permeability -- including after cycloheximide, i.e. with no new protein at all. "
           "Both are S. cerevisiae and both are MEASURED. A map from this block's effector "
           "pool to a protection factor would have to pick one, silently",
    missing="a dose-response for survival against a MEASURED Ctt1p abundance in the same "
            "cells, which is what would tell the two mechanisms apart"))

NOT_BUILT.add(Param.refused(
    "basal_endogenous_h2o2_production", "mM/min",
    _GOULEV + " fits eps = 0.06 mM/min; Tomalin 2016 (PMID 26944189) fits 5.28 uM/s against a "
    "MOCK 1 nM steady-state data point",
    reason="no measurement of the basal cytosolic H2O2 production rate exists for "
           "S. cerevisiae. Both published values are fitted, and one of them is anchored on "
           "an invented datum. A model written in the extracellular pool does not need it: "
           "the well starts at the pipetted dose and endogenous production does not reach it",
    missing="a HyPer7 or roGFP2-Tsa2dCR steady-state calibration in unstressed cells"))

NOT_BUILT.add(Param.refused(
    "density_exponent", "dimensionless",
    "no source. " + _BRANCO + " measured its rate constant at ONE density, 7e6 cells/mL, so "
    "the scaling of that rate with cell number is outside every experiment behind it",
    reason="the block assumes k proportional to X to the first power, and the committed "
           "plates cannot test it. Freeing the exponent alongside k0 buys 1.000, 1.021 and "
           "1.124x in residual RMS on the three replicates -- every one below the 0.146 "
           "reporter floor -- while the exponent itself lands on 1.052, 2.916 and 6.000, the "
           "last AT THE SEARCH BOUND. Three replicates of one design disagreeing by 6x on a "
           "quantity whose assumed value is 1 is an unidentifiability result, not a "
           "measurement, and density_exponent_identifiability() is it run rather than quoted. "
           "1.0 stays because it is the form Branco's own whole-culture rate constant is "
           "written in, and it is DECLARED here rather than left silent",
    missing="the dose ladder repeated at four inocula -- the same plate criterion (a)'s "
            "headline prediction needs. One design settles both"))

NOT_BUILT.add(Param.refused(
    "thiol_proteome_sink", "uM and 1/M/s",
    "Tomalin 2016 (PMID 26944189) Tables S6/S7: [Pr-SH] = 13000 uM at k = 5e-4 /uM/s",
    reason="the repository currently cites Tomalin inside STRESSORS['H2O2'] to explain the "
           "bolus/perfusion gap, with no organism flag. Tomalin is S. pombe, its "
           "peroxiredoxin is Tpx1 rather than Tsa1, its thiol pool size is ASSIGNED from a "
           "mammalian source, and Selvaggio reproduces the same biphasic threshold with no "
           "such sink at all. Two published explanations, and this module implements "
           "neither: it implements mass balance on a finite well, measured in the host",
    missing="a thiol-proteome titration in S. cerevisiae, and the five discriminating "
            "predictions Selvaggio section 4.4 lists, run in this organism rather than in "
            "HEK293"))


# --------------------------------------------------------------------------------------
# L2, as algebra: the two states criterion (d) will not let this block integrate
# --------------------------------------------------------------------------------------

YAP1_RELOCATION_TAU_H = 120.0 / 3600.0
"""120 s. MEASURED, """ + _GOULEV + """, 30 s sampling. Not a registry Param because no
equation multiplies by it: it is here so the reduction that absorbed the Yap1 state has a
number attached, and `mech/state.py`'s own `yap1_nuclear` row carries the same figure."""


def yap1_nuclear_fraction(e_mM, k_ex: float | None = None):
    """The sensor input function, algebraic. This is the ``u_j`` the other blocks lack.

    ``u(E) = E / (E + K_ex)``, a Michaelis form with no Hill exponent, because
    :data:`NOT_BUILT`'s ``yap1_yre_transfer_function`` records that no cooperativity has
    been measured for Yap1 at a YRE and inventing one would be a free scalar wearing a
    mechanism's name. It is the same shape as `Goulev 2017`'s ``g(H) = H/(H + K)``.

    Reduced rather than integrated: relocation is 120 s against a 4.14 h read, so
    ``tau/T = 0.00805`` against the 0.146 ceiling, and against the transcript STATE it drives
    it is 21.4-23.4x faster across the window's measured growth band (22.0x at the median mu).
    It is 28.9x against ``1/k_deg_mrna`` alone, and that is the wrong comparison: the block's
    own transcript tau is ``1/(k_deg + mu)`` and the ``+ mu`` is not optional.
    :func:`reduction_audit` prints the row.

    Args:
        e_mM: Extracellular H2O2, mM. Scalar or array.
        k_ex: Half-saturation, mM. Defaults to :data:`K_EX_MM`.
    """
    k = float(K_EX_MM) if k_ex is None else float(k_ex)
    if k <= 0:
        raise ValueError(f"K_ex must be positive, got {k}")
    e = np.asarray(e_mM, dtype=float)
    return e / (e + k)


def cytosolic_gradient() -> float:
    """The extracellular/cytosolic H2O2 ratio. REFUSED, with all three published values.

    The SPEED of the equilibration is not in dispute in any way that matters -- even the
    slowest of the three published permeabilities gives 60 s, so the reduction is licensed
    whichever is right. The GAIN is what is refused, and it is refused twice over: by the
    570x spread and by the peroxiredoxin cycle this block does not build.

    Raises:
        GradientRefused: always.
    """
    refusal = NOT_BUILT["k_inf_membrane_permeability"]
    raise GradientRefused(
        "the cytosolic/extracellular H2O2 gradient has no defensible value here. "
        f"{refusal.reason} What would close it: {refusal.missing}. Sources: {refusal.source}. "
        "This block is written in the extracellular pool so that it never forms this ratio; "
        "the gain is absorbed into K_ex, which is graded BOUNDED and carries the absorption "
        "in its own `missing` field")


# --------------------------------------------------------------------------------------
# L1: the well is not a reservoir
# --------------------------------------------------------------------------------------

def dose_decay_from_density_per_h(cells_per_ml: float) -> float:
    if isinstance(cells_per_ml, (bool, np.bool_)) or not math.isfinite(cells_per_ml) or cells_per_ml < 0:
        raise ValueError("cells_per_ml must be a finite nonnegative physical density")
    return float(K_REF_PER_H) * float(cells_per_ml) / float(X_REF_CELLS_PER_ML)


def dose_decay_rate_per_h(od600: float, cells_per_ml_per_od600: float) -> float:
    """First-order H2O2 consumption rate at one cell density, 1/h.

    ``k = k_ref * X / X_ref`` with ``X = od600 * cells_per_ml_per_od600``. The linearity in
    ``X`` is anchored at Branco's single density and is this block's one load-bearing
    untested structural claim; it is stated here rather than buried.

    Args:
        od600: Blank-corrected optical density.
        cells_per_ml_per_od600: A point of :data:`CELLS_PER_ML_PER_OD600`'s swept band. Pass
            it explicitly -- ``float()`` on a SWEPT parameter raises, and that is deliberate.
    """
    if od600 < 0:
        raise ValueError(f"od600 must be non-negative, got {od600}")
    point = CELLS_PER_ML_PER_OD600.at(float(cells_per_ml_per_od600))
    return dose_decay_from_density_per_h(float(od600) * point)


def dose_decay_band_per_h(od600: float) -> tuple[float, float]:
    """``(low, high)`` decay rate across the whole :data:`CELLS_PER_ML_PER_OD600` band."""
    low, high = CELLS_PER_ML_PER_OD600.bounds
    return (dose_decay_rate_per_h(od600, low), dose_decay_rate_per_h(od600, high))


def extracellular_h2o2(t_h, dose_mM: float, k0_per_h: float, growth_rate_per_h: float,
                       density_exponent: float = 1.0):
    """``E(t)`` in closed form for a well whose biomass grows exponentially.

    ``dE/dt = -k0 e^{n mu t} E`` integrates exactly to
    ``E(t) = E0 exp(-(k0/(n mu))(e^{n mu t} - 1))``. Closed form rather than a solver because
    it is exact, and because the fitting path evaluates it thousands of times.

    At ``k0 = 0`` this returns the pipetted dose forever, which is precisely the incumbent
    the block is scored against -- so the ablation's nesting is a point of the parameter
    space rather than a limit.

    Args:
        density_exponent: ``n`` in ``k = k_ref (X/X_ref)^n``. **1.0 is the block's one
            load-bearing untested structural claim**, and it is untested because Branco
            measured at a single density. The argument exists so
            :func:`density_exponent_identifiability` can ask the committed plates about it
            rather than leaving the assumption undeclared; nothing else here passes anything
            but 1.0.
    """
    t = np.asarray(t_h, dtype=float)
    mu = float(growth_rate_per_h)
    if mu <= 0:
        raise ValueError(f"growth rate must be positive, got {mu}")
    if k0_per_h < 0:
        raise ValueError(f"decay rate must be non-negative, got {k0_per_h}")
    nu = float(density_exponent) * mu
    if nu == 0.0:
        return float(dose_mM) * np.exp(-float(k0_per_h) * t)
    return float(dose_mM) * np.exp(-(float(k0_per_h) / nu) * np.expm1(nu * t))


def dose_half_life_h(k0_per_h: float, growth_rate_per_h: float) -> float:
    """Time for the pipetted dose to halve, hours. Infinite at ``k0 = 0``.

    Solved from the closed form, so it carries the growth of the consuming biomass: the
    half-life is shorter than ``ln2/k0`` because the culture that eats the dose is growing
    while it eats it.
    """
    if k0_per_h <= 0:
        return float("inf")
    mu = float(growth_rate_per_h)
    inner = 1.0 + mu * _LN2 / float(k0_per_h)
    return float(math.log(inner) / mu)


def time_integrated_dose(dose_mM: float, k0_per_h: float, growth_rate_per_h: float,
                         window_h: float, n_points: int = 20001) -> float:
    """``int_0^T E dt``, mM*h. The exposure a bolus delivers over a window."""
    t = np.linspace(0.0, float(window_h), int(n_points))
    return float(np.trapezoid(extracellular_h2o2(t, dose_mM, k0_per_h, growth_rate_per_h), t))


# --------------------------------------------------------------------------------------
# The three states, and the audit that licenses their shape
# --------------------------------------------------------------------------------------

def oxidative_states(window: Window, *, od600: float,
                     effector_turnover_per_h: float | None = None) -> MechState:
    """``(h2o2_extracellular, oxidative_mrna, oxidative_effector)`` with time constants
    COMPUTED for this window rather than tabulated.

    All three constants are properties of the operating point: the dose's is set by the cell
    density that eats it and by how fast that density grows, and both pools' are set by the
    window's MEASURED growth band. The same shape as `mech/ph.py::ph_states` and
    `mech/population.py::population_states`.

    ``h2o2_extracellular`` is declared :attr:`~ystwin.mech.state.Encoding.INTEGRATING`, and
    that is a substantive choice rather than bookkeeping. The ratio test in
    `mech/integrate.py` is derived for a state relaxing FROM zero TOWARD a value its driver
    sets; this state relaxes from its pipetted initial condition toward a driver of exactly
    zero. Its quasi-steady value is therefore 0, eliminating it deletes the entire dose, and
    the elimination bias is 100% however small ``tau/T`` is -- while ``tau/T`` remains the
    correct bound on the bias of FREEZING it, which is what the incumbent does. INTEGRATING
    is the encoding whose audit path offers exactly those two verdicts and never ELIMINATE,
    so it is the honest classification: a state with no quasi-steady value to be reduced to.

    Args:
        window: The declared operating window; supplies the measured growth band.
        od600: Starting optical density. The dose's time constant depends on it.
        effector_turnover_per_h: Degradation beyond dilution. Defaults to Ctt1's MEASURED
            ln2/108.8 h, the cross-protection carrier; pass Tsa1's for the peroxiredoxin.
    """
    delta = (float(PROTEIN_TURNOVER_CTT1_PER_H) if effector_turnover_per_h is None
             else float(effector_turnover_per_h))
    mu_lo, mu_hi = window.growth_rate_low_per_h, window.growth_rate_high_per_h
    k_lo, k_hi = dose_decay_band_per_h(od600)
    t_end = window.duration_h
    # The dose's rate rises with the biomass, so its tau is shortest at the end of the run.
    dose_taus = [1.0 / (k * math.exp(mu * t))
                 for k in (k_lo, k_hi) for mu in (mu_lo, mu_hi) for t in (0.0, t_end)]
    m_taus = [1.0 / (float(MRNA_DECAY_PER_H) + mu) for mu in (mu_lo, mu_hi)]
    p_taus = [1.0 / (mu + delta) for mu in (mu_lo, mu_hi)]

    return MechState((
        StateVar(
            name="h2o2_extracellular", units="mM",
            tau_h=(min(dose_taus), max(dose_taus)),
            encoding=Encoding.INTEGRATING,
            source=(
                f"DERIVED at this window's measured growth band {mu_lo:g}-{mu_hi:g} /h and at "
                f"OD600 {od600:g}, from tau = 1/(k_ref (X/X_ref) e^(mu t)) with k_ref = "
                f"{float(K_REF_PER_H):.4g} /h MEASURED ({_BRANCO}, refit by {_SELVAGGIO} SI "
                f"eq.919) and the cells/OD600 conversion swept over "
                f"{CELLS_PER_ML_PER_OD600.bounds[0]:.3g}-{CELLS_PER_ML_PER_OD600.bounds[1]:.3g}"),
            constrained_by=(
                "directly assayable in the spent medium by ferrous-oxidation-xylenol orange "
                "or Amplex Red, neither of which this project owns today. Named because "
                "criterion (b) asks for the assay, not for a plate reader it already has"),
            note=("It is INTEGRATING because its driver is zero and its initial condition is "
                  "the signal; see this function's docstring. At the plates' own inoculum "
                  "and their MEDIAN mu of 0.3227 /h it falls to 2.1-35.6% of the pipetted "
                  "dose within one hour and to 0.00-0.05% by the end of the read; across "
                  "this window's whole growth band those widen to 1.9-37.1% and 0.00-0.19%. "
                  "Either way the endpoint RFU is an area rather than a level.")),
        StateVar(
            name="oxidative_mrna", units="relative to the driver's own units",
            tau_h=(min(m_taus), max(m_taus)),
            encoding=Encoding.LEVEL,
            source=(
                f"DERIVED. tau = 1/(k_deg + mu) with k_deg = ln2/40 min MEASURED "
                f"({_GOULEV} Methods, attributed to Geisberg 2014, PMID 24529382) and mu the "
                f"window's own {mu_lo:g}-{mu_hi:g} /h. The + mu is not cosmetic: at the "
                f"plate's growth band dilution is 19-26% of the transcript's total loss"),
            constrained_by="RT-qPCR of TRX2 on the NativeYap1 strain (ystwin/qpcr.py)",
            note=("This is why the block does not inherit mech/state.py's global `mrna` row: "
                  "that row spans 3.6-32 min across methods and STRADDLES the reduction "
                  "threshold, while the antioxidant regulon's own 40 min clears it -- by "
                  "1.18-1.29x and no more. MRNA_DECAY_PER_H carries the 29.5-32.2 min "
                  "half-life at which this state would eliminate instead.")),
        StateVar(
            name="oxidative_effector", units="relative to the driver's own units",
            tau_h=(min(p_taus), max(p_taus)),
            encoding=Encoding.LEVEL,
            source=(
                f"DERIVED. tau = 1/(mu + k_deg) with k_deg = {delta:.5g} /h MEASURED "
                f"({_CHRISTIANO}) and mu {mu_lo:g}-{mu_hi:g} /h. At Ctt1's 108.8 h half-life "
                f"dilution is 97.4-98.3% of the loss, which is the whole content of the "
                f"cross-protection result: the memory is a generation counter"),
            constrained_by="mCitrine RFU on the committed NativeYap1 plates for the reporter "
                           "arm; a Ctt1p western for the protective arm, which is not an "
                           "instrument here",
            note=("One state, two readouts. The reporter and the protective protein share "
                  "the promoter, the transcript and the loss law, and differ only in the "
                  "MEASURED turnover added to dilution.")),
    ))


def oxidative_rhs(t: float, y, *, k0_per_h: float, growth_rate_per_h: float,
                  k_ex_mM: float | None = None,
                  effector_turnover_per_h: float | None = None):
    """``d(E, M, P)/dt``. Hours and mM throughout.

    ``dE/dt = -k0 e^{mu t} E``; ``dM/dt = u(E) - (k_deg_m + mu) M``;
    ``dP/dt = M - (mu + delta) P``. The Yap1 nuclear fraction is inside ``u`` and the
    cytosolic pool is inside ``K_ex``, both by the reduction argument in the module
    docstring.
    """
    mu = float(growth_rate_per_h)
    delta = (float(PROTEIN_TURNOVER_CTT1_PER_H) if effector_turnover_per_h is None
             else float(effector_turnover_per_h))
    e, m, p = float(y[0]), float(y[1]), float(y[2])
    u = float(yap1_nuclear_fraction(e, k_ex_mM))
    return [
        -float(k0_per_h) * math.exp(mu * float(t)) * e,
        u - (float(MRNA_DECAY_PER_H) + mu) * m,
        m - (mu + delta) * p,
    ]


def simulate_oxidative(window: Window, *, dose_mM: float, od600: float,
                       cells_per_ml_per_od600: float,
                       growth_rate_per_h: float | None = None,
                       k_ex_mM: float | None = None,
                       effector_turnover_per_h: float | None = None,
                       initial_effector: float = 0.0,
                       n_points: int = 601) -> Trajectory:
    """Integrate the three states across ``window`` with the package's pinned stiff driver.

    The system spans two to three orders of magnitude in time constant inside one plate read
    -- 3.7 min for the dose at OD600 0.80 and the top of the swept band, against 2.7-4.0 h
    for the effector -- so it goes through `mech/integrate.py`'s BDF rather than a
    hand-rolled loop. 3.7 min and not the 2.5 min `ARCHITECTURE_TARGET.md` tabulates: that
    figure implies 4.4e7 cells/mL per OD600, above the band :data:`CELLS_PER_ML_PER_OD600`
    sweeps, and :data:`K_REF_PER_H` records the correction.

    Args:
        initial_effector: ``P(0)``. Non-zero is how a pre-treated culture is started, which
            is what :func:`cross_protection` uses.
    """
    mu = window.growth_rate_low_per_h if growth_rate_per_h is None else float(growth_rate_per_h)
    k0 = dose_decay_rate_per_h(od600, cells_per_ml_per_od600)
    system = oxidative_states(window, od600=od600,
                              effector_turnover_per_h=effector_turnover_per_h)

    def rhs(t, y):
        return oxidative_rhs(t, y, k0_per_h=k0, growth_rate_per_h=mu, k_ex_mM=k_ex_mM,
                             effector_turnover_per_h=effector_turnover_per_h)

    return integrate_window(
        rhs,
        {"h2o2_extracellular": float(dose_mM), "oxidative_mrna": 0.0,
         "oxidative_effector": float(initial_effector)},
        window, system, n_points=n_points,
        atol={"h2o2_extracellular": 1e-9, "oxidative_mrna": 1e-9,
              "oxidative_effector": 1e-9})


def reduction_audit(window: Window, *, od600: float,
                    effector_turnover_per_h: float | None = None) -> ReducedSystem:
    """Criterion (d) for this block in this window, computed rather than claimed.

    On the 4.14 h plate read all three states are KEPT. On the 5 d fed-batch both pools
    ELIMINATE and only the dose survives -- and it survives on a technicality worth stating,
    that an INTEGRATING state is never eliminated, because by then the dose is long gone.
    The two irreducible sets are therefore not nested, which is the point of making the
    window an argument.
    """
    return reduce_for_window(
        oxidative_states(window, od600=od600,
                         effector_turnover_per_h=effector_turnover_per_h), window)


# --------------------------------------------------------------------------------------
# The chain in closed form, so the fit is linear in everything the plate does not measure
# --------------------------------------------------------------------------------------

DRIVER_FLOOR = 1e-12
"""Fraction of the pipetted dose below which the driver is treated as off, for quadrature only.

Not a modelling constant: :func:`extracellular_h2o2` is exact and never uses it. It is the
threshold that sets :func:`driver_horizon_h`, and at 1e-12 the neglected driver is at most
``1e-12 (dose + K_ex) / K_ex`` of the peak -- 2e-10 at the widest corner of this block's own
dose ladder and K_ex bracket, which is 8 orders below the quadrature error it replaces."""


def driver_horizon_h(k0_per_h: float, growth_rate_per_h: float,
                     floor: float = DRIVER_FLOOR,
                     density_exponent: float = 1.0) -> float:
    """When the well's H2O2 has fallen to ``floor`` of the pipetted dose. Infinite at ``k0 = 0``.

    Inverts the closed form of :func:`extracellular_h2o2` exactly: with ``nu = n mu``,
    ``t = ln(1 + nu ln(1/floor) / k0) / nu``. Beyond it the driver is off and the chain
    relaxes analytically, which is what lets :func:`effector_trace` spend its whole grid on
    the part of the window where anything is happening. Infinite when a sub-linear exponent
    leaves the well with a floor of its own above ``floor``.
    """
    if not 0.0 < float(floor) < 1.0:
        raise ValueError(f"floor must be in (0, 1), got {floor}")
    if k0_per_h < 0:
        raise ValueError(f"decay rate must be non-negative, got {k0_per_h}")
    if k0_per_h == 0:
        return float("inf")
    mu = float(growth_rate_per_h)
    if mu <= 0:
        raise ValueError(f"growth rate must be positive, got {mu}")
    decades = math.log(1.0 / float(floor))
    nu = float(density_exponent) * mu
    if nu == 0.0:
        return float(decades / float(k0_per_h))
    inner = 1.0 + nu * decades / float(k0_per_h)
    return float("inf") if inner <= 0.0 else float(math.log(inner) / nu)


def effector_trace(t_h, *, dose_mM: float, k0_per_h: float, growth_rate_per_h: float,
                   k_ex_mM: float | None = None,
                   effector_turnover_per_h: float | None = None,
                   density_exponent: float = 1.0,
                   n_grid: int = 1201):
    """``P(t)`` from a zero start, by exact convolution of the driver with the chain kernel.

    ``P(t) = int_0^t u(E(s)) g(t - s) ds`` with
    ``g(x) = (e^{-a x} - e^{-b x}) / (b - a)``, ``a = k_deg_m + mu``, ``b = mu + delta``.
    The KERNEL is exact; the convolution against it is trapezoidal, so what remains is a
    quadrature tolerance rather than a stiff-solver one.

    THE GRID TRACKS THE DRIVER AND NOT THE WINDOW, and that is the whole of the accuracy.
    A uniform grid over the requested window is wrong at exactly this block's operating
    point -- a dose consumed in minutes, read out over days. A 1201-point grid over a 120 h
    fed-batch has 0.1 h spacing against a dose half-life of 0.042 h at OD600 0.80 and the
    top of the swept band, and MISSES the spike that drives everything. Measured against
    Radau at rtol 1e-11 over each window's own operating box -- both ends of the swept
    cells/OD600 band, OD600 0.05-0.80, dose 0.1-0.5 mM, at that window's MEASURED growth
    band -- a window-tracking grid is 1.29e-4 worst-case on the 4.14 h read and **1.10e-1 on
    the 120 h fed-batch, which is only 1.3x BELOW the 0.146 reporter floor and so not a
    tolerance at all**: a quadrature rule whose error sits within 1.3x of the assay it is
    scored against has no margin, whichever side of the floor it lands on. Spending the same
    1201 points on ``[0, driver_horizon_h]`` and relaxing the two-exponential chain
    analytically past it gives 1.48e-5 and 1.80e-5 -- 9844x and 8124x below the floor, from
    1130x and 1.3x below it.

    The point of the closed form is that ``P`` does not depend on any fitted scalar, so
    :class:`YapReporter`'s fit is a linear least squares rather than a search.
    """
    mu = float(growth_rate_per_h)
    delta = (float(PROTEIN_TURNOVER_CTT1_PER_H) if effector_turnover_per_h is None
             else float(effector_turnover_per_h))
    a = float(MRNA_DECAY_PER_H) + mu
    b = mu + delta
    t = np.atleast_1d(np.asarray(t_h, dtype=float))
    t_max = float(np.max(t))
    if t_max <= 0:
        return np.zeros_like(t)
    t_conv = min(t_max, driver_horizon_h(float(k0_per_h), mu,
                                         density_exponent=float(density_exponent)))
    grid = np.linspace(0.0, t_conv, int(n_grid))
    h = grid[1] - grid[0]
    u = yap1_nuclear_fraction(
        extracellular_h2o2(grid, dose_mM, k0_per_h, mu, density_exponent), k_ex_mM)
    degenerate = abs(b - a) < 1e-12
    kernel = (grid * np.exp(-a * grid) if degenerate
              else (np.exp(-a * grid) - np.exp(-b * grid)) / (b - a))
    m_kernel = np.exp(-a * grid)
    # Trapezoidal convolution. kernel[0] = 0 so P needs only the u[0] half weight; M's
    # kernel starts at 1, so it needs both ends.
    p_grid = np.maximum(np.convolve(u, kernel)[:len(grid)] * h - 0.5 * h * u[0] * kernel, 0.0)
    m_grid = np.maximum(
        np.convolve(u, m_kernel)[:len(grid)] * h - 0.5 * h * (u[0] * m_kernel + u), 0.0)
    out = np.interp(t, grid, p_grid)
    tail = t > t_conv
    if np.any(tail):
        dt = t[tail] - t_conv
        m_end, p_end = float(m_grid[-1]), float(p_grid[-1])
        out[tail] = (p_end + m_end * dt) * np.exp(-b * dt) if degenerate else (
            p_end * np.exp(-b * dt)
            + m_end * (np.exp(-a * dt) - np.exp(-b * dt)) / (b - a))
    return out


# --------------------------------------------------------------------------------------
# Criterion (a): two nested reporter models that differ in ONE thing
# --------------------------------------------------------------------------------------

_REQUIRED_COLUMNS = ("time_h", "dose_mM", "signal", "od600")


@dataclass(frozen=True)
class YapReporter(FittableModel):
    """A YRE reporter driven by a dose that either depletes or does not.

    Three shapes, and they differ only in where ``k0`` comes from:

    * ``decay="none"`` -- the incumbent's DOSE TREATMENT, not the whole incumbent. The dose
      is the pipetted scalar forever, which is what
      `generator/stress_panel.py::module_response` does today: it takes a scalar dose and no
      time argument at all. What it does NOT carry over is that module's ``viability``
      term, which turns the panel's curve down past its peak -- and neither side here
      carries it, exactly as `mech/ablation.SaturatingReporter` does not. It is NOT omitted
      because it is small: on the truncated ladder the panel's own viability runs 1.000,
      0.997, 0.982 and 0.850, so at the 0.5 mM rung the headline is quoted at it is a 15%
      suppression, ABOVE the 14.6% reporter CV. It is omitted because putting it back makes
      the INCUMBENT worse rather than better -- residual RMS 472 to 488, 280 to 307 and 248
      to 273 on the three committed plates -- and moves the fitted decay rate by only 6-7%
      (1.966, 0.882, 0.669 to 2.084, 0.933, 0.715). A dose-dependent suppression that
      flattens the top of a dose-response is the one confound that could mimic a depleting
      pool, and on these plates it does not. :data:`GROWTH_HALVING_MM` bounds the omission;
      that measurement is what shows it is not carrying the result.
    * ``decay="declared"`` -- ``k0`` from the MEASURED ``k_ref``, the row's own OD600 and a
      declared point of the swept cells/OD600 band. Nothing about the decay is fitted.
    * ``decay="fitted"`` -- ``k0`` searched. Used only by :func:`shape_comparison`, where
      the question is what decay rate the plates themselves prefer, and where the strict
      nesting ``{R0, k_basal, gain} < {R0, k_basal, gain, k0_per_h}`` is what
      `mech/ablation.nested_out_of_sample` requires.

    Three of the four scalars are pure observation nuisances -- an initial pool, a basal
    transcription rate and an RFU gain -- and the model is LINEAR in all three, so the fit
    is an exact least-squares solve rather than a search. Nothing in the block's physics is
    fitted in the ``declared`` and ``none`` shapes.
    """

    name: str = "YapReporter"
    decay: str = "declared"
    cells_per_ml_per_od600: float = 3.0e7
    growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H
    k_ex_mM: float | None = None
    effector_turnover_per_h: float | None = None
    density_exponent: float = 1.0
    response_column: str = "signal"

    def __post_init__(self) -> None:
        if self.decay not in ("none", "declared", "fitted"):
            raise ValueError(
                f"decay must be 'none', 'declared' or 'fitted', got {self.decay!r}")
        if self.decay == "declared":
            CELLS_PER_ML_PER_OD600.at(float(self.cells_per_ml_per_od600))
        if self.growth_rate_per_h <= 0:
            raise ValueError(
                f"growth rate must be positive, got {self.growth_rate_per_h}")

    @property
    def parameter_names(self) -> tuple[str, ...]:
        base = ("R0", "k_basal", "gain")
        return base + ("k0_per_h",) if self.decay == "fitted" else base

    def _k0(self, od600: float, parameters: Mapping[str, float]) -> float:
        if self.decay == "none":
            return 0.0
        if self.decay == "fitted":
            return max(float(parameters["k0_per_h"]), 0.0)
        return dose_decay_rate_per_h(od600, self.cells_per_ml_per_od600)

    def _basis(self, data: pd.DataFrame, k0_fitted: float | None) -> np.ndarray:
        """``(n_rows, 3)`` design matrix for ``[R0, k_basal, gain]`` at a given ``k0``."""
        missing = [c for c in _REQUIRED_COLUMNS if c not in data.columns]
        if missing:
            raise KeyError(
                f"a YapReporter frame needs {', '.join(_REQUIRED_COLUMNS)}; missing "
                f"{', '.join(missing)}")
        t = np.asarray(data["time_h"], dtype=float)
        mu = float(self.growth_rate_per_h)
        decay = np.exp(-mu * t)
        columns = np.empty((len(data), 3), dtype=float)
        columns[:, 0] = decay
        columns[:, 1] = (1.0 - decay) / mu
        effector = np.zeros(len(data), dtype=float)
        keys = list(zip(np.asarray(data["dose_mM"], dtype=float),
                        np.asarray(data["od600"], dtype=float)))
        for dose, od in sorted(set(keys)):
            mask = np.array([k == (dose, od) for k in keys])
            k0 = (float(k0_fitted) if k0_fitted is not None
                  else self._k0(od, {"k0_per_h": 0.0}))
            effector[mask] = effector_trace(
                t[mask], dose_mM=dose, k0_per_h=k0, growth_rate_per_h=mu,
                k_ex_mM=self.k_ex_mM,
                effector_turnover_per_h=self.effector_turnover_per_h,
                density_exponent=self.density_exponent)
        columns[:, 2] = effector
        return columns

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        k0 = (max(float(parameters["k0_per_h"]), 0.0) if self.decay == "fitted" else None)
        basis = self._basis(data, k0)
        weights = np.array([float(parameters[n]) for n in ("R0", "k_basal", "gain")])
        return basis @ weights

    def _linear_fit(self, data: pd.DataFrame, k0: float | None):
        basis = self._basis(data, k0)
        y = np.asarray(data[self.response_column], dtype=float)
        weights, *_ = np.linalg.lstsq(basis, y, rcond=None)
        residual = basis @ weights - y
        return weights, residual

    def refit(self, data: pd.DataFrame) -> Fit:
        if self.decay == "fitted":
            def cost(x):
                return self._linear_fit(data, max(float(x[0]), 0.0))[1]

            solution = least_squares(cost, [1.0], bounds=([0.0], [200.0]), max_nfev=2000)
            k0 = float(max(solution.x[0], 0.0))
            weights, residual = self._linear_fit(data, k0)
            values = dict(zip(("R0", "k_basal", "gain"), map(float, weights)))
            values["k0_per_h"] = k0
            converged = bool(solution.status > 0)
        else:
            weights, residual = self._linear_fit(data, None)
            values = dict(zip(("R0", "k_basal", "gain"), map(float, weights)))
            converged = True
        return Fit(model=self.name, parameters=values,
                   n_free=len(self.parameter_names), n_rows=len(data),
                   residual_rms=float(np.sqrt(np.mean(residual ** 2))),
                   converged=converged)


def committed_yap1_blocks(construct: str = "NativeYap1") -> tuple[str, ...]:
    """Committed exports carrying a dose ladder for ``construct``, in manifest order."""
    from ..plate import replay

    found = []
    for export in replay.exports():
        try:
            frame = replay.load_doses(str(export))
        except (ValueError, FileNotFoundError):
            continue
        if construct in set(frame.construct):
            found.append(str(export))
    return tuple(found)


def yap1_frame(export: str, construct: str = "NativeYap1",
               od600: float | None = None,
               max_dose_mM: float | None = None) -> pd.DataFrame:
    """One committed ladder as a :class:`YapReporter` frame, truncated by a MEASURED dose.

    The truncation is not a chosen cutoff. This block takes the growth rate as an INPUT
    measured from the OD series and has no growth-inhibition term, so a well whose growth is
    more than halved is a well whose dilution correction the model does not describe.
    :data:`GROWTH_HALVING_MM` carries that dose per construct, MEASURED on these plates --
    0.79 mM for NativeYap1 -- which drops the 1.0 mM rung and every rung above it. The 2 and
    4 mM wells are dropped by the same rule and are separately the ones
    `generator/panel_experiment.py::PanelDataset.n_unusable` documents as going NEGATIVE
    after the dilution correction.

    Args:
        max_dose_mM: Override the measured truncation. ``None`` uses the construct's own
            growth-halving dose and refuses a construct that has none.
    """
    if max_dose_mM is None:
        try:
            max_dose_mM = GROWTH_HALVING_MM[construct]
        except KeyError:
            raise KeyError(
                f"no measured growth-halving dose for {construct!r}; this repository has one "
                f"for {', '.join(sorted(GROWTH_HALVING_MM))}. Pass max_dose_mM explicitly and "
                f"say where it came from -- a cutoff chosen to make a fit work is the thing "
                f"this argument exists to prevent") from None
    block = reporter_block(export, construct)
    block = block[block.dose_mM <= float(max_dose_mM)].reset_index(drop=True)
    block = block.copy()
    block["od600"] = (float(PLATE_INOCULUM_OD600) if od600 is None else float(od600))
    return block


def inoculum_endpoint_observable(data: pd.DataFrame, *, od600: float, dose_mM: float,
                                 read_h: float = PLATE_READ_4H.duration_h) -> Observable:
    """Endpoint RFU at a declared dose and a declared INOCULUM, from wells read at another.

    The observable the harvest verifier named, and the one a constant dose cannot produce:
    a well dosed at 0.5 mM reaches a different endpoint at OD 0.05 than at OD 0.80 because
    the denser culture eats the dose sooner. Both models are fitted to the same real rows at
    the plates' own inoculum and then asked for a well that was never run, which is the same
    shape as `mech/population.py::chord_growth_observable`'s declared read length.

    Scored RELATIVELY, so the comparison is dimensionless -- but note that the two models
    fit DIFFERENT gains to the same rows, so the plate's scale does not cancel here and the
    effect moves with the replicate. :func:`inoculum_spread_observable` is the version in
    which that SCALE cancels exactly, and it is the one the headline number comes from.
    """
    for column in _REQUIRED_COLUMNS:
        if column not in data.columns:
            raise KeyError(f"a reporter frame needs {', '.join(_REQUIRED_COLUMNS)}")
    query = pd.DataFrame({"time_h": [float(read_h)], "dose_mM": [float(dose_mM)],
                          "signal": [0.0], "od600": [float(od600)]})

    def summarise(model: FittableModel, fit: Fit) -> float:
        return float(model.predict(query, fit.parameters)[0])

    return Observable(
        name=REPORTER_ACTIVITY_FLOOR.name, units="RFU",
        assay=REPORTER_ACTIVITY_FLOOR.assay, scored="relative", data=data,
        summarise=summarise,
        description=(f"endpoint RFU at {float(dose_mM):g} mM and inoculum OD600 "
                     f"{float(od600):g}, read at {float(read_h):g} h"))


def inoculum_ablation(export: str, *, od600: float, dose_mM: float = 0.5,
                      cells_per_ml_per_od600: float = 3.0e7,
                      construct: str = "NativeYap1",
                      growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
                      read_h: float = PLATE_READ_4H.duration_h) -> AblationResult:
    """THE HEADLINE. Delete the depleting dose and the inoculum stops mattering.

    Both models carry the same three fitted observation nuisances and are refit on the same
    committed rows; they differ in one thing only, whether ``k0`` is the MEASURED
    ``k_ref X/X_ref`` or zero. The reduced model's prediction is identical at every
    inoculum by construction, because a constant dose has no density to depend on.

    Args:
        export: A committed export name.
        od600: The inoculum the endpoint is predicted at. The plates ran at 0.16.
        dose_mM: Which rung of the ladder to quote.
        cells_per_ml_per_od600: A declared point of the swept band.
        growth_rate_per_h: MEASURED, from the plate band. Not fitted and not asserted.
    """
    data = yap1_frame(export, construct)
    full = YapReporter(name="depleting dose (k_ref measured)", decay="declared",
                       cells_per_ml_per_od600=cells_per_ml_per_od600,
                       growth_rate_per_h=growth_rate_per_h)
    reduced = YapReporter(name="dose held forever (incumbent)", decay="none",
                          growth_rate_per_h=growth_rate_per_h)
    observable = inoculum_endpoint_observable(data, od600=od600, dose_mM=dose_mM,
                                              read_h=read_h)
    return ablate(full, reduced, observable, REPORTER_ACTIVITY_FLOOR)


def inoculum_spread_observable(data: pd.DataFrame, *, dose_mM: float = 0.5,
                               od_low: float = 0.05, od_high: float = 0.80,
                               read_h: float = PLATE_READ_4H.duration_h) -> Observable:
    """The endpoint RFU at a low inoculum over the endpoint at a high one, same pipetted dose.

    What divides out exactly is the plate-to-plate SCALE: multiply ``R0``, ``k_basal`` and
    ``gain`` by any constant and this ratio is unchanged, which is what makes it
    replicate-robust where the single-inoculum endpoint is not. The gain alone does NOT
    divide out -- the two basal terms are common to numerator and denominator, so the ratio
    moves with ``gain`` relative to them (1.39 to 5.46 as the fitted gain is scaled 0.25x to
    4x on the first committed plate). The incumbent still predicts exactly 1.0 at any
    parameters, because a dose held forever has no cell density to depend on; the effect is
    therefore the whole spread -- but it is a fitted spread and not a parameter-free one.
    """
    if not 0 < float(od_low) < float(od_high):
        raise ValueError(f"need 0 < od_low < od_high, got {od_low} and {od_high}")
    query = pd.DataFrame({
        "time_h": [float(read_h), float(read_h)],
        "dose_mM": [float(dose_mM), float(dose_mM)],
        "signal": [0.0, 0.0],
        "od600": [float(od_low), float(od_high)]})

    def summarise(model: FittableModel, fit: Fit) -> float:
        predicted = model.predict(query, fit.parameters)
        return float(predicted[0] / predicted[1])

    return Observable(
        name=REPORTER_ACTIVITY_FLOOR.name, units="dimensionless",
        assay=REPORTER_ACTIVITY_FLOOR.assay, scored="relative", data=data,
        summarise=summarise,
        description=(f"endpoint RFU at OD600 {float(od_low):g} over OD600 {float(od_high):g}, "
                     f"same {float(dose_mM):g} mM pipetted dose, read at {float(read_h):g} h"))


def inoculum_spread_ablation(export: str, *, dose_mM: float = 0.5,
                             od_low: float = 0.05, od_high: float = 0.80,
                             cells_per_ml_per_od600: float = 3.0e7,
                             construct: str = "NativeYap1",
                             growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
                             read_h: float = PLATE_READ_4H.duration_h) -> AblationResult:
    """Criterion (a), in the form the gain cannot influence: the inoculum SPREAD.

    Same two models and the same fitted rows as :func:`inoculum_ablation`, scored on a
    dimensionless ratio instead of on an RFU. The reduced model's answer is 1.0 by
    construction and the full one's is the density spread the depleting pool produces.

    Scored against :data:`SPREAD_RATIO_FLOOR`, which is 0.146 propagated through a ratio to
    0.2065, and NOT against the single-reading CV -- a ratio carries both readings' noise.
    """
    data = yap1_frame(export, construct)
    full = YapReporter(name="depleting dose (k_ref measured)", decay="declared",
                       cells_per_ml_per_od600=cells_per_ml_per_od600,
                       growth_rate_per_h=growth_rate_per_h)
    reduced = YapReporter(name="dose held forever (incumbent)", decay="none",
                          growth_rate_per_h=growth_rate_per_h)
    observable = inoculum_spread_observable(data, dose_mM=dose_mM, od_low=od_low,
                                            od_high=od_high, read_h=read_h)
    return ablate(full, reduced, observable, SPREAD_RATIO_FLOOR)


def fitted_decay_rates(construct: str = "NativeYap1",
                       growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
                       max_dose_mM: float | None = None) -> pd.DataFrame:
    """CRITERION (c), and this is the strongest unfitted result in the block.

    Let the plates choose the dose decay rate for themselves -- ``decay="fitted"``, nothing
    told about Branco -- and compare what they choose with the band ``k_ref`` predicts at the
    plates' own inoculum.

    On the ladder truncated at NativeYap1's MEASURED growth-halving dose, the three committed
    replicates prefer **1.966, 0.882 and 0.669 /h** against a measured band of
    **0.876-3.283 /h**; two land inside it and the third is 1.31x below its floor. Nothing
    about the block was fitted to Branco 2004, the measurement is from 2004 in a different
    lab by a different assay, and these are 2026 reporter traces.

    It reverses if the truncation is dropped: including the 1.0 mM rung the same fits give
    0.635, 0.037 and 0.000 /h, i.e. two of three plates prefer no decay at all. That rung is
    above this construct's own measured growth-halving dose, so it is a well whose dilution
    correction the model does not describe -- which is exactly why
    :data:`GROWTH_HALVING_MM` and not a chosen cutoff decides the ladder.

    HOW MANY LAND INSIDE DEPENDS ON WHICH mu OF THE WINDOW'S OWN BAND IS PASSED, and the
    quoted result is the MIDDLE of that dependence rather than the best of it. The dose's
    decay competes with the same mu that dilutes the reporter, so a slower assumed mu makes
    the plates ask for a faster k0: at 0.245 /h they choose 2.700, 1.241 and 0.948 /h and all
    THREE land inside the band; at the median 0.3227 /h, 1.966, 0.882 and 0.669 and two do;
    at 0.363 /h, 1.659, 0.732 and 0.557 and one does. What is band-free and mu-free is the
    part that matters: no decay at all is rejected on every plate at every mu, and the rate
    the plates ask for is of order 1 /h -- which is what a 2004 whole-culture measurement
    predicts and what a held dose does not.

    AND THE UNFITTED RATE OUT-PREDICTS THE FITTED ONE ACROSS PLATES, which is a sharper form
    of the same claim and does not use the band's width at all. That comparison is
    :func:`transfer_comparison`, which runs it rather than quoting it.

    Returns one row per committed export with the fitted rate, both models' residuals, and
    the measured band the plate's own inoculum implies.
    """
    low, high = dose_decay_band_per_h(float(PLATE_INOCULUM_OD600))
    rows = []
    for export in committed_yap1_blocks(construct):
        data = yap1_frame(export, construct, max_dose_mM=max_dose_mM)
        fitted = YapReporter(decay="fitted", growth_rate_per_h=growth_rate_per_h).refit(data)
        held = YapReporter(decay="none", growth_rate_per_h=growth_rate_per_h).refit(data)
        k0 = float(fitted.parameters["k0_per_h"])
        rows.append({
            "export": export,
            "n_rows": len(data),
            "max_dose_mM": float(data.dose_mM.max()),
            "fitted_k0_per_h": k0,
            "measured_k0_low_per_h": low,
            "measured_k0_high_per_h": high,
            "inside_measured_band": bool(low <= k0 <= high),
            "rms_fitted": fitted.residual_rms,
            "rms_no_decay": held.residual_rms,
            "rms_improvement": held.residual_rms / fitted.residual_rms,
        })
    return pd.DataFrame(rows)


def transfer_comparison(construct: str = "NativeYap1",
                        growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
                        od600: float | None = None) -> pd.DataFrame:
    """Criterion (c) in its sharpest form: a 2004 measurement beats a rate fitted on plates.

    Hold one committed plate out, fit ``k0`` on the other two, and score the held-out plate
    with only its three observation nuisances refit. Then score the SAME held-out plate with
    Branco's ``k_ref`` at a declared point of the swept band, which has learned nothing from
    any plate at all. The band's width plays no part here -- each column is one number.

    Measured on this checkout: the transferred fit chooses 0.765, 1.179 and 1.351 /h for
    held-out RMS **285.6, 175.4 and 184.5** against the incumbent's 472.2, 280.2 and 247.7;
    Branco's rate at 8e6 cells/mL/OD600 gives **269.4, 169.3 and 152.0** -- better than the
    transferred fit on all three plates, and 1.63-1.75x better than the incumbent. At the
    band's other end, 3e7, it beats the incumbent on two plates of three and loses on the
    third (233.4, 243.6, 265.0), which is the swept axis doing what a swept axis does.

    Args:
        od600: The inoculum Branco's rate is evaluated at. Defaults to the plates' own.
    """
    od = float(PLATE_INOCULUM_OD600) if od600 is None else float(od600)
    exports = committed_yap1_blocks(construct)
    frames = {export: yap1_frame(export, construct) for export in exports}
    scorer = YapReporter(decay="fitted", growth_rate_per_h=growth_rate_per_h)
    low, high = CELLS_PER_ML_PER_OD600.bounds
    rows = []
    for held in exports:
        others = pd.concat([frames[e] for e in exports if e != held], ignore_index=True)
        transferred = float(scorer.refit(others).parameters["k0_per_h"])

        def rms_at(k0: float) -> float:
            return float(np.sqrt(np.mean(scorer._linear_fit(frames[held], k0)[1] ** 2)))

        rows.append({
            "held_out": held,
            "k0_transferred_per_h": transferred,
            "rms_transferred": rms_at(transferred),
            "rms_branco_low_conversion": rms_at(dose_decay_rate_per_h(od, low)),
            "rms_branco_high_conversion": rms_at(dose_decay_rate_per_h(od, high)),
            "rms_incumbent": rms_at(0.0),
        })
    frame = pd.DataFrame(rows)
    frame["branco_beats_transferred"] = (
        frame.rms_branco_low_conversion < frame.rms_transferred)
    return frame


def density_exponent_identifiability(
        construct: str = "NativeYap1",
        growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
        bounds: tuple[float, float] = (-2.0, 6.0)) -> pd.DataFrame:
    """Can the committed plates test the block's one structural assumption? Measured: NO.

    ``k = k_ref (X/X_ref)^n`` with ``n = 1`` is anchored at Branco's SINGLE density, and the
    harvest verifier named it as this block's one load-bearing untested claim. The plates
    grow 3.8x inside a 4.14 h read at their own median mu, so the consuming density does
    move within a well and ``n`` is in principle visible without a second inoculum.

    It is not visible in practice, and that is the result. Freeing ``n`` alongside ``k0``
    buys **1.000, 1.021 and 1.124x** in residual RMS on the three committed replicates --
    all three below the 0.146 reporter floor -- while ``n`` itself lands on **1.052, 2.916
    and 6.000**, the last of those AT THE SEARCH BOUND. Three replicates of one experimental
    design, and they do not agree to within a factor of six on a quantity whose assumed
    value is 1. So the linearity is neither confirmed nor refuted here; it is UNTESTED, and
    :data:`NOT_BUILT`'s ``density_exponent`` records that with the experiment that would
    settle it -- the same ladder at four inocula that criterion (a)'s headline already
    needs. One plate design answers both.

    Returns one row per committed export. ``clears_floor`` is the honest column: it is False
    on every plate.
    """
    exports = committed_yap1_blocks(construct)
    rows = []
    for export in exports:
        data = yap1_frame(export, construct)
        linear = YapReporter(decay="fitted", growth_rate_per_h=growth_rate_per_h)
        rms_linear = linear.refit(data).residual_rms

        def cost(x, frame=data):
            model = YapReporter(decay="fitted", growth_rate_per_h=growth_rate_per_h,
                                density_exponent=float(x[1]))
            return model._linear_fit(frame, max(float(x[0]), 0.0))[1]

        solution = least_squares(cost, [1.0, 1.0],
                                 bounds=([0.0, bounds[0]], [200.0, bounds[1]]),
                                 max_nfev=600)
        free_rms = float(np.sqrt(np.mean(cost(solution.x) ** 2)))
        exponent = float(solution.x[1])
        rows.append({
            "export": export,
            "k0_free_exponent_per_h": float(max(solution.x[0], 0.0)),
            "fitted_exponent": exponent,
            "at_search_bound": bool(min(abs(exponent - bounds[0]),
                                        abs(exponent - bounds[1])) < 1e-6),
            "rms_free_exponent": free_rms,
            "rms_linear": rms_linear,
            "rms_improvement": rms_linear / free_rms,
            "clears_floor": bool(rms_linear / free_rms - 1.0
                                 > REPORTER_ACTIVITY_FLOOR.noise_floor),
        })
    return pd.DataFrame(rows)


def shape_comparison(export: str, *, construct: str = "NativeYap1",
                     cells_per_ml_per_od600: float = 3.0e7,
                     growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H) -> NestedComparison:
    """THE ONE THAT FAILS, and it is reported for that reason.

    Leave-one-dose-out on a committed ladder, with the decay rate FITTED in the full model
    and zero in the reduced one, which is the strict nesting Tier 0.3 asks for. Two things
    come out: whether the depleting dose predicts a held-out rung better than a held one --
    it does not, above the floor -- and what decay rate the plates themselves prefer, which
    is worth having next to the measured 0.88-3.28 /h whatever the verdict is.
    """
    data = yap1_frame(export, construct)
    full = YapReporter(name="depleting dose (k0 fitted)", decay="fitted",
                       cells_per_ml_per_od600=cells_per_ml_per_od600,
                       growth_rate_per_h=growth_rate_per_h)
    reduced = YapReporter(name="dose held forever (incumbent)", decay="none",
                          growth_rate_per_h=growth_rate_per_h)
    return nested_out_of_sample(full, reduced, data, REPORTER_ACTIVITY_FLOOR)


# --------------------------------------------------------------------------------------
# The bolus / perfusion gap the repository cites and cannot model
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BolusPerfusionGap:
    """What a depleting extracellular pool predicts about the gap, against what is measured.

    Args:
        predicted_ratio: Bolus dose over perfusion dose for EQUAL time-integrated exposure.
        measured_ratio_low, measured_ratio_high: The two constructs' growth-halving doses
            over Goulev's perfusion arrest dose.
        exposure_bolus_mM_h, exposure_perfusion_mM_h: What the two deliveries integrate to.
    """

    od600: float
    cells_per_ml_per_od600: float
    growth_rate_per_h: float
    window_h: float
    predicted_ratio: float
    measured_ratio_low: float
    measured_ratio_high: float
    exposure_bolus_mM_h: float
    exposure_perfusion_mM_h: float

    @property
    def overshoot(self) -> float:
        """How many times too large the exposure prediction is, against its nearest measured end."""
        return self.predicted_ratio / self.measured_ratio_high

    @property
    def peak_criterion_error(self) -> float:
        """The rival criterion's error factor: a peak-concentration rule predicts exactly 1.0."""
        return self.measured_ratio_low / 1.0

    def report(self) -> str:
        return "\n".join([
            f"BOLUS / PERFUSION GAP at OD600 {self.od600:g}, "
            f"{self.cells_per_ml_per_od600:.3g} cells/mL/OD, mu {self.growth_rate_per_h:g} /h, "
            f"T = {self.window_h:g} h",
            f"  a {float(GOULEV_PERFUSION_ARREST_MM):g} mM PERFUSION delivers "
            f"{self.exposure_perfusion_mM_h:.4g} mM.h",
            f"  the same nominal BOLUS delivers {self.exposure_bolus_mM_h:.4g} mM.h",
            f"  so equal exposure needs a bolus {self.predicted_ratio:.3g}x the perfusion dose",
            f"  MEASURED ratio is {self.measured_ratio_low:.3g}-{self.measured_ratio_high:.3g}x "
            f"(growth halving at 0.79 and 1.24 mM against Goulev's 0.6 mM arrest)",
            f"  VERDICT: the sign is right and the size is not -- exposure overshoots its "
            f"nearest measured end by {self.overshoot:.3g}x. A peak-concentration criterion "
            f"predicts exactly 1.0x and misses its nearest end by only "
            f"{self.peak_criterion_error:.3g}x, i.e. it is "
            f"{self.overshoot / self.peak_criterion_error:.3g}x closer than exposure",
        ])

    def __str__(self) -> str:
        return self.report()


def bolus_perfusion_gap(*, od600: float | None = None,
                        cells_per_ml_per_od600: float = 3.0e7,
                        growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
                        window_h: float = PLATE_READ_4H.duration_h) -> BolusPerfusionGap:
    """Does the gap fall out of a depleting extracellular pool? Measured here: NO.

    `generator/stress_panel.py::STRESSORS['H2O2']` records the gap and attributes it to
    Tomalin 2016's thiol buffering, in *S. pombe*, with an assigned mammalian thiol pool.
    This block offers the alternative the assignment asked to test -- mass balance on a
    finite well, measured in the host -- and the answer is that mass balance predicts the
    gap with the right sign and the wrong size.

    Both published explanations fail in the SAME direction. Tomalin's thiol sink would be an
    ADDITIONAL consumption term on top of this one, so adding it makes the overshoot worse,
    not better. What the overshoot says is that the lethal endpoint is not set by
    time-integrated exposure; a peak-concentration criterion predicts a ratio of exactly 1.0
    and misses by only the measured 1.32-2.07x.

    THE SIZE OF THE FAILURE IS PROPORTIONAL TO ``window_h`` AND THAT IS NOT A MEASURED
    QUANTITY. The bolus exposure saturates within the dose's own half-life while the
    perfusion arm integrates forever, so the ratio grows almost linearly in T: 2.09x at
    T = 0.5 h, 3.64x at 1 h, 14.8x at 4.14 h and 85.9x at 24 h, at OD600 0.16 and 3e7
    cells/mL/OD. At T below about half an hour this block would CLEAR
    :data:`LETHAL_RATIO_FLOOR` rather than fail it. The default is the plate read because
    that is the window this project measures in; Goulev's arrest was scored on a
    microfluidic device over an unstated exposure, so the honest statement is that the
    exposure criterion fails on the plate's own window and that the comparison has no
    window-free form. The SIGN does not depend on T.
    """
    od = float(PLATE_INOCULUM_OD600) if od600 is None else float(od600)
    k0 = dose_decay_rate_per_h(od, cells_per_ml_per_od600)
    perfusion = float(GOULEV_PERFUSION_ARREST_MM)
    exposure_perfusion = perfusion * float(window_h)
    exposure_bolus = time_integrated_dose(perfusion, k0, growth_rate_per_h, window_h)
    return BolusPerfusionGap(
        od600=od, cells_per_ml_per_od600=float(cells_per_ml_per_od600),
        growth_rate_per_h=float(growth_rate_per_h), window_h=float(window_h),
        predicted_ratio=exposure_perfusion / exposure_bolus,
        measured_ratio_low=0.79 / perfusion, measured_ratio_high=1.24 / perfusion,
        exposure_bolus_mM_h=exposure_bolus, exposure_perfusion_mM_h=exposure_perfusion)


# --------------------------------------------------------------------------------------
# Criterion (c): cross-protection, carried with no memory parameter
# --------------------------------------------------------------------------------------

GUAN_GROWTH_BAND_PER_H = (_LN2 / (140.0 / 60.0), _LN2 / (90.0 / 60.0))
"""(0.2971, 0.4621) /h. DERIVED from BNID 108255 in the vendored BioNumbers Nov-2024 export:
a normal laboratory haploid doubles in ~90 min in YPD and ~140 min in synthetic medium.
`Guan 2012` reports its memory in GENERATIONS and not in hours, so the band is what converts
its result into this block's units, and the conversion is the block's whole claim."""


def memory_persistence_generations(growth_rate_per_h: float,
                                   fraction_remaining: float = 0.05) -> dict[str, float]:
    """Generations of outgrowth over which a carrier falls to ``fraction_remaining``.

    Two carriers, one threshold. ``g = -log2(f) / (1 + k_deg/mu)``, so the THRESHOLD divides
    out of the ratio between them and the ratio is a prediction with no free scalar in it at
    all. That is the quantity :func:`cross_protection` scores against Guan 2012.

    Returns a mapping with the protein route, the transcript route, their ratio, and the
    fractional departure of the protein route from a pure generation counter.
    """
    mu = float(growth_rate_per_h)
    if mu <= 0:
        raise ValueError(f"growth rate must be positive, got {mu}")
    if not 0.0 < float(fraction_remaining) < 1.0:
        raise ValueError(f"fraction must be in (0, 1), got {fraction_remaining}")
    decades = -math.log2(float(fraction_remaining))
    protein = decades / (1.0 + float(PROTEIN_TURNOVER_CTT1_PER_H) / mu)
    transcript = decades / (1.0 + float(MRNA_DECAY_PER_H) / mu)
    return {
        "protein_generations": protein,
        "transcript_generations": transcript,
        "ratio": protein / transcript,
        "departure_from_generation_counter": float(PROTEIN_TURNOVER_CTT1_PER_H) / mu,
    }


@dataclass(frozen=True)
class CrossProtection:
    """A mild pre-dose, a washout, a severe challenge -- integrated, not asserted.

    Args:
        effector_at_challenge: ``P`` carried into the challenge by the pre-treated culture.
        effector_naive: ``P`` in a culture that never saw the pre-dose. Zero by construction.
        generations_of_outgrowth: How long the washout lasted, in generations.
        persistence_ratio: Protein route over transcript route, from
            :func:`memory_persistence_generations`. Threshold-free.
    """

    pre_dose_mM: float
    challenge_mM: float
    generations_of_outgrowth: float
    growth_rate_per_h: float
    effector_at_challenge: float
    effector_naive: float
    effector_at_pretreatment_end: float
    persistence_ratio: float
    guan_generations: tuple[float, float] = (4.0, 5.0)

    @property
    def fraction_retained(self) -> float:
        if self.effector_at_pretreatment_end <= 0:
            return 0.0
        return self.effector_at_challenge / self.effector_at_pretreatment_end

    @property
    def departure_from_generation_counter(self) -> float:
        """How far the carrier's loss departs from pure dilution, as a fraction."""
        return float(PROTEIN_TURNOVER_CTT1_PER_H) / self.growth_rate_per_h

    def report(self) -> str:
        low, high = self.guan_generations
        ratio_lo, ratio_hi = (memory_persistence_generations(GUAN_GROWTH_BAND_PER_H[1])["ratio"],
                              memory_persistence_generations(GUAN_GROWTH_BAND_PER_H[0])["ratio"])
        floor = CROSS_PROTECTION_FLOOR.noise_floor
        return "\n".join([
            f"CROSS-PROTECTION  {self.pre_dose_mM:g} mM pre-dose, "
            f"{self.generations_of_outgrowth:g} generations of outgrowth, then "
            f"{self.challenge_mM:g} mM",
            f"  effector at the end of the pre-treatment  {self.effector_at_pretreatment_end:.6g}",
            f"  effector entering the challenge           {self.effector_at_challenge:.6g} "
            f"({self.fraction_retained:.1%} retained)",
            f"  a naive culture enters at                 {self.effector_naive:.6g}",
            "  WHAT IS PREDICTED, with no memory parameter anywhere in the block:",
            f"    the carrier's loss departs from pure dilution by "
            f"{self.departure_from_generation_counter:.2%} at mu = "
            f"{self.growth_rate_per_h:g} /h, and by "
            f"{float(PROTEIN_TURNOVER_CTT1_PER_H) / GUAN_GROWTH_BAND_PER_H[1]:.2%}-"
            f"{float(PROTEIN_TURNOVER_CTT1_PER_H) / GUAN_GROWTH_BAND_PER_H[0]:.2%} across "
            f"Guan's own doubling times (BNID 108255), so the memory is a GENERATION counter "
            f"to well inside the {floor:.1%} spread Guan quotes",
            f"    and it outlives the transcript that made it by {ratio_lo:.3g}-{ratio_hi:.3g}x "
            f"in generations over that band, threshold-free -- the detection threshold divides "
            f"out of the ratio",
            f"  MEASURED, and not fitted to: Guan 2012 (PMID 22851651) reports persistence of "
            f"{low:g}-{high:g} GENERATIONS, transmission to daughters that never saw the "
            f"pre-treatment, and no requirement for nascent protein synthesis -- three "
            f"independent signatures of a dilution-limited carrier",
            "  NOT PREDICTED, and not claimed: the absolute 4-5. It needs a detection "
            "threshold on the protection, and NOT_BUILT['effector_to_protection_map'] records "
            "why this block refuses to supply one",
        ])

    def __str__(self) -> str:
        return self.report()


def cross_protection(*, pre_dose_mM: float = 0.1, challenge_mM: float = 1.0,
                     generations_of_outgrowth: float = 4.5,
                     od600: float | None = None,
                     cells_per_ml_per_od600: float = 3.0e7,
                     growth_rate_per_h: float = PLATE_MEDIAN_GROWTH_PER_H,
                     pre_treatment_h: float = 2.0) -> CrossProtection:
    """Criterion (c): a mild dose leaves a carrier that a later severe dose meets.

    Integrated with the same three states and the same constants used everywhere else, and
    with NO memory parameter -- the persistence is ``1/(mu + k_deg)`` with ``k_deg`` the
    MEASURED Ctt1 turnover, so 97-98% of the loss is dilution and the memory is a generation
    counter. Guan 2012 measured that memory at 4-5 generations, transmitted to daughters
    that never saw the pre-treatment, and independent of new protein synthesis.

    THE PRE-DOSE HERE IS H2O2 AND GUAN'S WAS SALT. Guan's protocol is a mild NaCl
    pre-treatment followed by a severe H2O2 challenge, and this block has no osmotic input,
    so the leg simulated is the H2O2-to-H2O2 self-protection arm rather than Guan's own
    cross-protection arm. What is scored against Guan is the part that does not depend on
    which stress made the carrier: its loss law. The absolute retained fraction below is a
    declared scenario, not a reproduction of Guan's experiment.

    What this does NOT do is convert the carried effector into a protection factor. That map
    is :data:`NOT_BUILT` because Guan and Branco measured two different mechanisms for the
    same phenomenon; see the module docstring.
    """
    od = float(PLATE_INOCULUM_OD600) if od600 is None else float(od600)
    mu = float(growth_rate_per_h)
    pre_window = Window(name="PRE_TREATMENT", duration_h=float(pre_treatment_h),
                        growth_rate_low_per_h=mu, growth_rate_high_per_h=mu,
                        source="a declared pre-treatment leg, not a measured vessel")
    pre = simulate_oxidative(pre_window, dose_mM=pre_dose_mM, od600=od,
                             cells_per_ml_per_od600=cells_per_ml_per_od600,
                             growth_rate_per_h=mu)
    p_end = float(pre.of("oxidative_effector")[-1])
    # Washout: the dose is gone, so the carrier is lost by dilution plus its own turnover.
    outgrowth_h = float(generations_of_outgrowth) * _LN2 / mu
    p_challenge = p_end * math.exp(-(mu + float(PROTEIN_TURNOVER_CTT1_PER_H)) * outgrowth_h)
    return CrossProtection(
        pre_dose_mM=float(pre_dose_mM), challenge_mM=float(challenge_mM),
        generations_of_outgrowth=float(generations_of_outgrowth), growth_rate_per_h=mu,
        effector_at_challenge=p_challenge, effector_naive=0.0,
        effector_at_pretreatment_end=p_end,
        persistence_ratio=memory_persistence_generations(mu)["ratio"])


# --------------------------------------------------------------------------------------
# The GEM link. stress_ros.py already prices a dose; this converts into its units.
# --------------------------------------------------------------------------------------

def gem_peroxide_uptake_mmol_per_gdcw_h(*, dose_mM: float, od600: float,
                                        cells_per_ml_per_od600: float,
                                        gdcw_per_od: float = 0.42) -> float:
    """The uptake flux this block's consumption rate represents, in the GEM's own units.

    ``fba/stress_ros.py`` makes peroxide dosable and prices it; what it needs is a flux in
    mmol/gDCW/h, and what this block produces is a well concentration falling at ``k E``.
    The conversion is ``k E / (gdcw_per_od * od600)``, with ``k E`` in mmol/L/h and the
    denominator in gDCW/L.

    THE FLUX IS LARGE AND SHORT. At the plates' 0.5 mM rung and their own inoculum it is
    **6.51-24.43 mmol/gDCW/h** across the swept conversion, i.e. up to 2.4x the glucose
    uptake the price list in ``stress_ros`` was measured at -- across the dose's own 12.3 to
    42.3 minute half-life. It does not depend on the inoculum, because ``k`` scales with the density
    and so does the biomass it is divided by.

    THE GEM PRICES THIS AT ZERO IN A WILD TYPE, which is worth knowing before calling it.
    Run through ``fba/stress_ros.peroxide_dose_response`` on yeast-GEM v9.0.2 at glucose 10,
    both ends of the swept conversion cost **0.000** of growth with catalase present, because
    catalase is cofactor-free; delete it (YGR088W, with or without the peroxisomal YDR256C)
    and the same two fluxes cost **5.45% and 20.43%** of growth. The cost is entirely in the
    thiol routes. Peroxide toxicity is Fenton chemistry and a stoichiometric model has none
    of it, so this is a floor on the burden and not an estimate of it.

    Args:
        gdcw_per_od: `generator/plate.py::PlateConditions.gdcw_per_od`. ASSERTED there and
            carried here as an argument rather than re-declared, since this module does not
            own it.
    """
    if od600 <= 0:
        raise ValueError(f"od600 must be positive to form a per-gDCW flux, got {od600}")
    k = dose_decay_rate_per_h(od600, cells_per_ml_per_od600)
    return k * float(dose_mM) / (float(gdcw_per_od) * float(od600))


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------

def provenance_summary(heading_level: int = 3) -> str:
    """The block's constants, its refusals and its gate, as a markdown table."""
    return provenance_table(OXIDATIVE_PARAMS, NOT_BUILT, heading_level=heading_level)
