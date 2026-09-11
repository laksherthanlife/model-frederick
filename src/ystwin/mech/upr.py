"""ER folding and the UPR: unfolded client -> Kar2 titration -> Ire1 -> HAC1 splicing ->
Hac1 protein -> UPRE occupancy -> chaperone expression, with the entry function REFUSED.

This is the one block in the sweep whose readout already exists on this project's plates.
UPRE1 and UPRE2 are built strains (`ystwin/qpcr.py` maps both to Hac1/DTT), dosed on a
seven-rung DTT ladder and read at 10-minute intervals on FOUR committed exports. So unlike
every other mechanistic block here, this one is scored against measured data rather than
against a published figure, and :func:`score_dtt_ladder` does exactly that -- one fit, seven
held-out (plate, construct) blocks, every score printed.

**THE VERDICT, FIRST, BECAUSE IT IS NOT THE ONE THIS BLOCK WAS COMMISSIONED TO REACH.** The
two conditions the build list attached are met: the entry function is REFUSED
(:data:`DTT_ENTRY`, :class:`DoseEntry`) and the Hac1-protein and UPRE-occupancy states are
restored on two MEASURED binding constants, which removes two of Pincus's fitted ones. But
criterion (a) FAILS, and it fails on this project's own committed data rather than on an
argument: removing the one surviving ODE state moves the growth-corrected activity fold by
0.0029x the plate floor. The state is therefore not claimed. Two further results came out of
building it -- a wrong observation model in the first draft of this module, corrected below,
and a REFUSAL (:data:`DEACTIVATION_MECHANISM`) naming the measurement that would tell this
project what actually terminates its UPR.

WHAT IS BUILT. A reduction of Pincus et al. 2010 (PLoS Biol 8(7):e1000415, PMID 20625545)
Text S1, from fifteen states to six, in normalised units. The published model is retyped
from a two-page supplementary PDF -- there is no SBML, no MATLAB and no BioModels entry --
and it is adapted rather than adopted, because four things about it do not transfer:

1. **It has no growth term anywhere.** Every decay constant in it silently absorbs one mu.
   Its BiP decay ``g[B] = 1.39e-4 /s`` ("2 h") is sourced to Axelsen & Sneppen 2004, a
   modelling paper, and is contradicted by a direct SILAC measurement of **27.0 h** for
   Kar2 (Christiano et al. 2014, Cell Rep 9:1959-65, PMID 25466257). Two hours is close to
   dilution in a 90-minute culture, so the constant was standing in for the growth term the
   model does not have. Here the chaperone pool is lost at ``mu + k_deg`` with mu MEASURED
   per DOSE off the OD600 series (the three replicate wells of a condition share one fitted
   rate, which is how ``outputs/sensor_characterisation.csv`` reports it), and Christiano's rates are already dilution-corrected
   (his Methods correct Kdeg for Kdil at 2.5 and 4 h doubling), so the sum composes.
2. **There is no Hac1 protein state and no occupancy variable.** The step this project's
   owner named as the worked example -- protein folding, binding to a site, expression of
   another protein -- is the step Pincus collapses into one fitted rational function
   ``f(H) = H^2/(296.5 + 5.26*H + H^2)``. It is restored here, and it *removes* two fitted
   scalars rather than adding any: see :data:`HAC1_UPRE2_KD_NM` and :data:`HAC1_MAX_OCCUPANCY`.
3. **The dose axis is unmeasured.** There is no published mapping from mM DTT in the medium
   to molecules of DTT in the yeast ER, and no Ki for tunicamycin against yeast Alg7. Both
   entry points are :data:`DTT_ENTRY` and :data:`TUNICAMYCIN_ENTRY`, both REFUSED, and
   `ARCHITECTURE_GAPS.md` 0.5 is right that writing a saturating function there would be
   the status quo with a new name. What a caller gets instead is :class:`DoseEntry`, which
   makes the caller state on the record that the load-vs-dose curve is a declared axis and
   not a measurement. **The block predicts SHAPE, not absolute dose.**
4. **Four of its states are unobservable and buy nothing.** The U / U.B / Ud / Ud.B
   disulfide bookkeeping costs four states, three unmeasured rate constants and five orders
   of magnitude of stiffness, and no assay in this repository touches any of them.

WHAT WAS REFUSED OUTRIGHT, so nobody re-derives it. Every entry in :data:`REFUSED_FROM_PINCUS`
is a number the source supplies that this module will not use. The two that matter most:
basal HAC1 mRNA = 200 copies/cell is cited in Text S1 to "Walter, P., Personal
Communication" and sets the splicing substrate pool; and the Ire1-BiP ratio R = 12.5 is
inferred from **mammalian** IRE1-alpha co-immunoprecipitation (Bertolotti 2000) in a paper
whose own words are "no real quantitative measurements have been recorded". Writing the
chain in units of basal total HAC1 makes the first one never arise, and eliminating the
Ire1-BiP complex makes the second one never arise. **Stroberg et al. 2018 is not used at
all**: the verification pass retrieved the published Table 1 and found that not one value in
the harvested transcription of it appears in the paper.

THE STOPPING RULE, ROW BY ROW.

* **(a) ABLATION -- THIS BLOCK FAILS IT, and that is its most useful result.**
  :func:`upr_ablation` removes the chaperone feedback -- the only arm that lets the promoter
  rate depend on time at all -- with both models refit on the same committed rows
  through `mech/ablation.py`. On the growth-corrected activity fold at 1 mM it moves the
  reading by **0.00042 relative, 0.0029x the 0.146 plate floor**: 344 times below it. Nor
  does the sign rescue it -- across the seven held-out blocks the FROZEN model scores better
  on six of the seven, by at most 0.0023 relative RMS, 1.6% of the floor. Under this build's
  own stopping rule that means the ODE state MUST NOT be claimed, and it is not. What the
  block ships is the algebra with MEASURED binding constants in it, the state retained only
  as the ablation's comparison arm and priced in public.
* **(b) IDENTIFIABILITY.** Every state names the assay that constrains it or is a declared
  sweep axis; :data:`UPR_STATES` will not construct otherwise. The ER interior is honestly
  unconstrained and says so.
* **(c) EMERGENCE.** The held-out construct and the held-out plates: fit UPRE1 on one plate,
  predict the other SEVEN (plate, construct) blocks with no refit at all. The reporters share
  one Hac1 pool and one UPRE gain by construction, which is a mechanistic claim and not a
  fitted transform. MEASURED relative RMS **0.1417 (UPRE2, same plate), 0.0791 and 0.1410
  (20260728), 0.0994 and 0.0747 (20260803), 0.1562 and 0.1846 (20260804)** against the 0.146
  floor -- **five of seven clear**, both misses are the same plate, and
  :func:`score_dtt_ladder` prints all of them. The 20260728 pair is new to this pass and was
  not won by tuning: that plate was being excluded by ``read_h``, a column of the scored
  frame that no equation ever read -- see :data:`_PLATE_COLUMNS`.
  **AND IT IS NOT EVIDENCE OF MECHANISM, which is measured rather than conceded.** A bare
  three-parameter Hill in dose, ``1 + A*d^n/(K^n + d^n)`` -- no state, no mu, no measured
  binding constant, ONE FEWER free scalar -- fitted to the same twelve wells reproduces the
  whole held-out column: it differs from this block by at most **0.0033 relative RMS, 0.022x
  the floor**, on any of the eight blocks, and it scores BETTER on four of the seven held
  out. Criterion (c) as scored here says the dose-response SHAPE transfers between plates and
  constructs, which is worth knowing; it does not discriminate this mechanism from the
  incumbent algebra. :func:`score_dtt_ladder` records the comparison and test_mech_upr.py
  re-runs it.
* **(d) REDUCTION.** :func:`upr_reduction` runs `mech/integrate.py`'s tau/T audit. On the
  4.14 h plate read exactly one of the six states survives it -- the chaperone pool, at
  tau/T = 0.67-0.99 -- and :func:`reduction_bias` MEASURES what eliminating the other five
  costs on the number the plate actually produces: **0.0007 to 0.049 across the sub-lethal
  ladder, at most 0.34x the floor**, so the reduction is free to the assay that would have
  to see it. Read that next to (a) and the two together say something the ratio test alone
  could not: the five states it eliminates cost less than the plate can measure, and the one
  state it keeps buys less still. On a 4.14 h plate read this whole chain is algebra, and
  that is now a measurement on the observable rather than an inference from tau/T.
* **(e) FREE-SCALAR GATE.** :data:`UPR_PARAMS` (the six-state reference chain) declares
  **13 free scalars against 4 independent targets** and ``require_gate()`` on it RAISES.
  The shipped plate shape, :data:`PLATE_PARAMS` and :class:`PlateUpr`, declares **4 against
  the same 4** and passes. :data:`DROPPED_TO_PASS_THE_GATE` is the list of nine and the
  reason each one goes.

THE OBSERVABLE, AND THE DEFECT THAT MOVED IT. This block is scored on ``activity_late`` from
``outputs/sensor_characterisation.csv`` -- this repository's own growth-corrected promoter
activity, computed by ``scripts/run_sensor_characterisation.py`` from the total-signal route
with a time-resolved mu and the real OD600 traces. An earlier version of this module scored
instead on the accumulated increment of ``plate/replay.py::load_doses``'s ``signal`` column
under an ``exp(mu*t)`` biomass weight, on the ground that a TOTAL well fluorescence
integrates ``X(s)k(s)``. Those sheets do not hold a total. ``plate/dose_response.py`` says
in its own first line that they hold NORMALISED fluorescence, and the numbers agree: the
late-window mean of that column reproduces ``naive_late``, the per-cell column of the
committed file, to 0.7% on the fitted block and to 1.3% at worst across all three plates,
both constructs and all seven doses. A per-cell signal obeys ``dS/dt = k - mu*S``, so its
increment is not a biomass integral. The correction is not cosmetic -- mu falls from 0.214
to 0.161 /h across this ladder, and the dilution artefact that produced was being read as
mechanism: it manufactured an apparent 1.19x ablation where the corrected observable gives
0.0029x, and it was the whole of the 4.7x miss on the 20260803 plate, which is now among
the best-predicted. What it did NOT change is the identifiability: two of the four scalars
still go to the box on either observable, because twelve wells at four doses do not pin four
scalars. :func:`upre_activity_fold` is the corrected
observable and the old one is deleted rather than kept alongside.

A SECOND DEFECT, FOUND BY ASKING WHETHER THE DECLARED BASAL STATE IS ACTUALLY STATIONARY,
and it was not. The chain normalises its Hac1 pool so that ``p = P_MAX`` at FULL splicing,
because :data:`HAC1_ABUNDANCE_MOLECULES` is an induced abundance -- but the divisor that
did it was ``_spliced_fraction(1.0, mu)``, which is a ratio to that same quantity and so is
identically 1.0. The production term was therefore short by the spliced share at full
splicing, ``k_spl/(mu + gamma + k_spl) = 0.629`` at mu = 0.214 /h, while
:meth:`UprChain.basal` computed the initial condition the other way. The two disagreed by
1.59x, so the chain's own zero-dose control DECAYED 37% across the read and every
occupancy it reported was a real induction on top of a spurious relaxation. The divisor is
now derived inside :meth:`UprChain.basal` from the same basal transcription rate the
equations use, which makes it exact for any :data:`HAC1_ISOFORM_DECAY_RATIO` as well, and
the basal state is a fixed point to 1e-15. What it changed: :func:`reduction_bias`, the
criterion (d) number, from a flat 0.024-0.066 to a hump peaking at 0.049 near 0.2 mM and
collapsing to 0.0007 at 1 mM -- **0.45x the floor to 0.34x**, and the shape is now the one
the mechanism predicts, since both models saturate at the same occupancy at the top of the
ladder. What it did NOT change: anything scored on the plate. :class:`PlateUpr` never
called the chain, so the fit, the six ladder scores and the ablation are identical to five
decimal places.

UNITS. Hours throughout, and normalised state variables throughout: the chaperone pool, the
KAR2 message and total HAC1 are each 1 at the unstressed steady state, the ER client pool is
in units of the chaperone-client Kd, and the Hac1 pool is in units of its own DNA-binding
Kd. That is not cosmetic. Pincus's model is in molecule counts with volume-divided rate
constants, so an uncited ``V_ER = 2.15 um^3`` silently scales every bimolecular term, and
its basal BiP of 430,000 molecules does not match SGD's record of the source it cites
(Ghaemmaghami 2003 = 337,000). In Ho, Baryshnikova & Brown 2018 (PMID 29361465), as SGD
curates it, IRE1 spans 60-3,100 molecules/cell across ten datasets and KAR2 spans
1,318-337,000 across seventeen -- a 52-fold and a 256-fold spread. Normalising is what keeps
those spreads out of the answer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares

from ..plate import replay
from ..reporter import ReporterKinetics, SpecificFluorescence, promoter_activity
from .ablation import (
    REPORTER_ACTIVITY_FLOOR,
    AblationResult,
    Fit,
    FittableModel,
    Observable,
    ablate,
    free_scalar_gate,
)
from .integrate import (
    PINNED_METHOD,
    ReducedSystem,
    Trajectory,
    integrate_window,
    reduce_for_window,
)
from .params import (
    FreeScalarGate,
    Param,
    ParamRegistry,
    RefusedValue,
    Target,
    provenance_table,
)
from .state import PLATE_READ_4H, MechState, StateVar, Window

__all__ = [
    "CHAIN_QUADRATURE_POINTS",
    "DEACTIVATION_FRACTION",
    "DEACTIVATION_MECHANISM",
    "DOSE_LADDER_MM",
    "DROPPED_TO_PASS_THE_GATE",
    "DTT_ENTRY",
    "DoseEntry",
    "EXCESS_FLOOR",
    "EntryRefused",
    "FrozenChaperoneUpr",
    "HAC1_MAX_OCCUPANCY",
    "INTERNAL_GRID_POINTS",
    "LATE_WINDOW_FRACTION",
    "LadderScore",
    "PLATE_PARAMS",
    "PLATE_STATE",
    "PlateUpr",
    "REFUSED_FROM_PINCUS",
    "RESPONSE_DEACTIVATION_H",
    "STARTING_POINTS",
    "SUBLETHAL_DOSE_MM",
    "TUNICAMYCIN_ENTRY",
    "UPR_PARAMS",
    "UPR_STATES",
    "UPR_TARGETS",
    "DEACTIVATION_VERSUS_RELAXATION",
    "UprChain",
    "chaperone_relaxation_rate",
    "deactivation_rate",
    "deactivation_time",
    "plate_gate",
    "provenance",
    "read_duration_h",
    "reduction_bias",
    "score_dtt_ladder",
    "upr_ablation",
    "upr_reduction",
    "upre_activity_fold",
    "upre_activity_trace",
]

_LN2 = math.log(2.0)
_AVOGADRO = 6.02214076e23


# --------------------------------------------------------------------------------------
# The constants, graded
# --------------------------------------------------------------------------------------

UPR_PARAMS = ParamRegistry("mech/upr.py -- six-state reference chain")
"""Every constant the reference chain contains. Its gate FAILS on purpose: 13 free against 4."""

PLATE_PARAMS = ParamRegistry("mech/upr.py -- plate-window shape")
"""The shipped shape. Same measured rows, nine fewer free ones. Its gate passes."""


def _both(param: Param) -> Param:
    UPR_PARAMS.add(param)
    PLATE_PARAMS.add(param)
    return param


HAC1_MRNA_DECAY_PER_H = _both(Param.measured(
    "hac1_mrna_decay", 60.0 / 20.0, "1/h",
    "MEASURED. Sidrauski, Cox & Walter 1996, Cell 87:405-413, PMID 8898194, Fig 5c: HAC1 "
    "mRNA lifetime 20 min. INHERITED SECOND-HAND through Pincus 2010 Text S1 (PMID 20625545) "
    "-- Cell 1996 is not open access and the verification pass could not read Fig 5c itself. "
    "Pincus applies the same constant to the spliced and unspliced isoforms; that half is "
    "ASSERTED and is carried separately as hac1_isoform_decay_ratio"))

HAC1_SPLICING_PER_H = _both(Param.measured(
    "hac1_splicing_rate", 60.0 / 11.0, "1/h",
    "MEASURED. Kawahara, Yanagi, Yura & Mori 1997, Mol Biol Cell 8:1845-1862, PMID 9348528, "
    "Fig 3: 11 min for HAC1 splicing. This is the rate at FULL Ire1 activation, so the "
    "splicing flux here is k_spl * a * Hu with a the active-Ire1 fraction -- which is what "
    "removes Pincus's separate Ire1 copy number from the chain. The chemistry is much "
    "faster than 11 min (Korennykh 2011 measures k2 = 0.25 /s on the stem-loop, ~4 s), so "
    "the 11 min is substrate delivery and ligation, not cleavage"))

HAC1_PROTEIN_DECAY_PER_H = _both(Param.measured(
    "hac1_protein_decay", _LN2 / (1.25 / 60.0), "1/h",
    "MEASURED. Pal, Chakravarty, Chakraborty, Chakraborty & Das 2007, PMID 17108329: "
    "'the half-life of Hac1p in the parental W303 cells was 1-1.5 min', by cycloheximide "
    "protein-synthesis shutoff. The midpoint 1.25 min is used and the 1-1.5 min bracket is "
    "the ci95. CAVEAT the source carries: measured at 37 C after a temperature shift, in "
    "tunicamycin-induced cells",
    ci95=(_LN2 / (1.5 / 60.0), _LN2 / (1.0 / 60.0))))

KAR2_PROTEIN_DECAY_PER_H = _both(Param.measured(
    "kar2_protein_decay", _LN2 / 27.0, "1/h",
    "MEASURED. Christiano, Nagaraj, Frohlich & Walther 2014, Cell Rep 9:1959-65, "
    "PMID 25466257, dynamic SILAC: Kar2 half-life 27.0 h (confirmed through SGD's curation "
    "of that reference; the numbers live in a supplementary table). Their Methods correct "
    "Kdeg for dilution at 2.5 and 4 h doubling times, so it is a degradation-only rate and "
    "(mu + k_deg) composes correctly. At mu = 0.28 /h this is 8.4% of the loss -- chaperone "
    "turnover is essentially pure dilution, which is precisely the job Pincus's fitted 2 h "
    "constant was doing in a model with no growth term. MIND THE CONVENTION when quoting the "
    "size of that contradiction: 27.0 h is a HALF-life and Pincus's g[B] = 1.39e-4 /s is a "
    "RATE whose reciprocal, 2.0 h, is a MEAN lifetime. Compared like for like the two differ "
    "by 0.5004 / 0.025672 = 19.5x, not by the 27/2 = 13.5x an earlier draft of this module "
    "quoted -- that number divided a half-life by a mean lifetime"))

HAC1_UPRE2_KD_NM = _both(Param.measured(
    "hac1_upre2_kd", 427.0, "nM",
    "MEASURED. Fordyce, Pincus, Kimmig, Nelson, El-Samad, Walter & DeRisi 2012, PNAS "
    "109:E3084-93, PMID 23054834: 'Fits to the UPRE-2 binding data yielded a Kd of "
    "427 +/- 37 nM', by MITOMI 2.0. CAVEAT: yeast Hac1i protein, but expressed by in vitro "
    "translation in RABBIT RETICULOCYTE LYSATE and assayed on isolated oligonucleotide. "
    "DO NOT quote the micromolar figures in the same paper ('UPRE-2 Kd = 497 +/- 60 uM') -- "
    "they come from the systematic-mutation experiment and are ~1000x this one",
    ci95=(390.0, 464.0)))

HAC1_ABUNDANCE_MOLECULES = _both(Param.measured(
    "hac1_abundance_max", 8970.0, "molecules/cell",
    "MEASURED. Ho, Baryshnikova & Brown 2018, Cell Syst 6:192-205, PMID 29361465, via SGD's "
    "curation: HAC1 spans 1,778-8,970 molecules/cell with a median of 2,069 across five "
    "proteome datasets. The TOP of that range is used, and WHY MATTERS. An earlier draft "
    "said it is 'the closest thing to an induced abundance in the set'; that is not "
    "supported -- this project's own verification pass warns that Ho's HAC1 values come from "
    "conditions in which Hac1p should be nearly absent, so none of the five is an induced "
    "abundance and the set is not a safe calibration. The real reason the choice is "
    "admissible is the same degeneracy that frees the compartment (see HAC1_MAX_OCCUPANCY): "
    "theta enters the observable only as (b + theta)/(b + theta0), so a rescaled ceiling is "
    "absorbed by the fitted upre_basal_share. MEASURED across Ho's whole range -- refitting "
    "at the MEDIAN 2,069 (a 4.3x lower abundance, a 2.82x lower ceiling) moves the worst of "
    "the eight ladder scores by 0.0044 relative RMS, 0.030x the assay floor, and at the "
    "MINIMUM 1,778 by 0.00024, 0.0017x. test_mech_upr.py measures that rather than asserting "
    "it. The block therefore does not rest on which end of Ho's range is taken"))

CELL_VOLUME_UM3 = _both(Param.measured(
    "haploid_cell_volume", 42.0, "um^3",
    "MEASURED. BNID 100427, the same haploid cell volume mech/state.py already uses for its "
    "membrane A/V ratio. Used ONCE, to turn Ho 2018's molecule counts into the nM the "
    "Fordyce Kd is quoted in. The compartment is DECLARED WHOLE-CELL: Hac1 acts in the "
    "nucleus (~3 um^3), where the same count is ~14x more concentrated. See "
    "HAC1_MAX_OCCUPANCY for why that 14x does not move the prediction here"))

_HAC1_MAX_NM = (float(HAC1_ABUNDANCE_MOLECULES) / _AVOGADRO
                / (float(CELL_VOLUME_UM3) * 1e-15) * 1e9)

HAC1_MAX_OCCUPANCY = _both(Param.derived(
    "hac1_max_occupancy", _HAC1_MAX_NM / (float(HAC1_UPRE2_KD_NM) + _HAC1_MAX_NM),
    "fraction of UPRE sites bound",
    f"DERIVED from two MEASURED rows and one declared compartment. "
    f"{float(HAC1_ABUNDANCE_MOLECULES):g} molecules (Ho 2018, PMID 29361465) in "
    f"{float(CELL_VOLUME_UM3):g} um^3 (BNID 100427) is {_HAC1_MAX_NM:.1f} nM; against the "
    f"MITOMI Kd of {float(HAC1_UPRE2_KD_NM):g} nM (Fordyce 2012, PMID 23054834) the "
    f"single-site Langmuir occupancy is {_HAC1_MAX_NM / (float(HAC1_UPRE2_KD_NM) + _HAC1_MAX_NM):.4f}. "
    "THIS IS THE ROW THAT REPLACES PINCUS'S FITTED a0 = 296.5 AND a1 = 5.26. The declared "
    "compartment does not move the prediction, and the reason is NOT that the knee is out "
    "of reach -- an earlier draft of this row said so and it is false. At the fitted "
    "parameters the occupancy is hard against this ceiling from 0.5 mM up. What makes the "
    "compartment free is that the occupancy enters the observable only as "
    "(b + theta)/(b + theta0), so rescaling theta rescales the fitted upre_basal_share with "
    "it: refitting at a 3 um^3 nucleus raises the ceiling to 0.9208, raises basal_share by "
    "2.030x against a 2.029x ceiling, and moves the worst of the eight ladder scores by "
    "5.1e-05 relative RMS -- 0.00035x the assay floor. test_mech_upr.py MEASURES that "
    "degeneracy rather than asserting it"))

DTT_LETHAL_MM = _both(Param.measured(
    "dtt_lethal_dose", 1.55, "mM",
    "MEASURED on this project's own plates and carried in generator/stress_panel.py: growth "
    "falls to half of control at 1.45 mM on UPRE1 and 1.65 mM on UPRE2. It is the block's "
    "ADMISSIBILITY BOUND, not a model parameter -- the 2 and 5 mM rungs of the committed "
    "ladder sit above it and this module refuses to score them"))

# -- free rows, reference chain only ------------------------------------------------------

FOLDING_RATE_PER_H = UPR_PARAMS.add(Param.swept(
    "folding_rate", "1/h",
    "SWEPT. Pincus Text S1 asserts a 20 min folding time with NO citation; Stroberg 2018 "
    "asserts the same 20 min from Braakman & Hebert 2013, a MAMMALIAN review. No ER folding "
    "rate constant has ever been reported for S. cerevisiae. The axis is Pincus's own "
    "sensitivity range, 12-20 min. It sets only how fast the ER client pool relaxes, and "
    "the plate shape eliminates that pool, which is why this row is not in PLATE_PARAMS",
    bounds=(60.0 / 20.0, 60.0 / 12.0),
    missing="a pulse-chase folding rate for a defined ER client in S. cerevisiae"))

KAR2_MRNA_DECAY_PER_H = UPR_PARAMS.add(Param.asserted(
    "kar2_mrna_decay", 60.0 / 25.0, "1/h",
    "ASSERTED. Pincus Text S1 gives g[Bm] = 6.67e-4 /s (25 min) with no citation. Kept in "
    "the reference chain and ELIMINATED in the plate shape, where the KAR2 message is at "
    "quasi-steady state and the constant cancels exactly against its own basal transcription",
    missing="a KAR2 mRNA half-life in S. cerevisiae; the repo's own RT-qPCR could bound it"))

HAC1_ISOFORM_DECAY_RATIO = UPR_PARAMS.add(Param.asserted(
    "hac1_isoform_decay_ratio", 1.0, "dimensionless",
    "ASSERTED. Pincus sets the spliced and unspliced HAC1 decay constants equal. Axelsen & "
    "Sneppen 2004's whole result is that the response is optimal when the ACTIVE message "
    "decays fast relative to the passive one, and Cherry et al. 2019 (PMID 30874502) shows "
    "distinct decay pathways acting on the HAC1 splicing species, so 1.0 is a convention "
    "and not a measurement. It CANCELS in the plate shape: the occupancy there depends on "
    "the spliced fraction relative to its own value at full splicing, and this ratio scales "
    "both",
    missing="isoform-resolved HAC1 mRNA half-lives in S. cerevisiae"))

IRE1_HILL_N = _both(Param.bounded(
    "ire1_hill_n", 4.5, "dimensionless",
    "BOUNDED by MEASUREMENT, which is the one place a fitted Pincus parameter becomes a "
    "bracketed one. Korennykh et al. 2009, Nature 457:687-693, PMID 19079236, measured Hill "
    "coefficients for RNase activation by self-association of the yeast Ire1 cytosolic "
    "fragment IN VITRO: n = 2 for Ire1-KR and Ire1-KR24, n = 3.5-8 for Ire1-KR32. Pincus "
    "fits 4.5 and it sits inside that. The cooperativity is REQUIRED by the data, not "
    "decoration: Pincus Fig S9c-d shows n = 1 destroys both the switch-like ire1bipless dose "
    "response and the washout delay",
    bounds=(2.0, 8.0),
    missing="an in-vivo Hill coefficient for Ire1 activation in S. cerevisiae"))

BASAL_FOLDING_LOAD = _both(Param.swept(
    "basal_folding_load", "fraction of folding capacity",
    "SWEPT. The unstressed ER client influx as a fraction of the folding capacity that "
    "clears it. Pincus fits S_U = 310 molecules/s against a basal folding pool of 372,149 "
    "and no yeast measurement of either exists; his own sensitivity range is 250-370. Here "
    "it is written as a dimensionless saturation so the molecule counts, the uncited "
    "V_ER = 2.15 um^3 and the 256-fold KAR2 abundance spread all drop out. THIS IS THE "
    "PARAMETER A HETEROLOGOUS SECRETED PRODUCT WOULD ENTER THROUGH",
    bounds=(1e-3, 0.97),
    missing="a measured secretory flux into the ER for these strains"))

UPRE_BASAL_SHARE = _both(Param.swept(
    "upre_basal_share", "occupancy-equivalent",
    "SWEPT. The UPRE-independent share of the promoter rate, in units of occupancy: the "
    "transcription rate is (b + theta) rather than (1 + N*theta), which is the same "
    "one-parameter family written so that b -> 0 is a point of the box instead of a limit "
    "at infinity. ONE value serves both the KAR2 regulon and the reporter, which is the "
    "block's claim that both are read through the same Hac1 pool and the same site class -- "
    "and it is what makes the UPRE2 prediction held out rather than refitted. Fordyce 2012 "
    "qualifies it: Hac1 is the first bZIP known to use two site geometries, so equal gain "
    "is the null this ladder tests",
    bounds=(0.0, 50.0),
    missing="a measured basal-to-induced ratio for these two UPRE constructs"))

# -- refusals -----------------------------------------------------------------------------

DTT_ENTRY = _both(Param.refused(
    "dtt_entry", "fraction of oxidative folding capacity lost per mM",
    "REFUSED. Searched: Pincus 2010 Text S1 parameterises D0 (DTT molecules inside the ER) "
    "as an unvalued input 'adjusted as a function of time to represent the dose', anchored "
    "only against a fitted basal Ero1 pool that is internally inconsistent by 2.16x in its "
    "own paper and 300-800x above the measured Ero1 abundance; Merksamer 2008 (PMID "
    "19026441) doses 2 mM and measures the consequence, not the lumenal concentration",
    reason="there is no published mapping from mM DTT in the medium to DTT in the yeast ER "
           "lumen, so a saturating function written here would be the incumbent Hill with a "
           "new name -- which ARCHITECTURE_GAPS.md 0.5 names as the failure mode of the "
           "first architecture",
    missing="a lumenal thiol measurement against medium DTT in S. cerevisiae -- an eroGFP "
            "or roGFP calibration curve would do it. Until then use DoseEntry, which puts "
            "the declared axis on the record"))

TUNICAMYCIN_ENTRY = UPR_PARAMS.add(Param.refused(
    "tunicamycin_entry", "fractional excess client influx per ug/mL",
    "REFUSED. Searched: no Ki for tunicamycin against yeast Alg7 is reported anywhere. "
    "Merksamer 2008 (PMID 19026441) measures the CONSEQUENCE in yeast at 1 ug/mL -- ER "
    "redox does not move for 60 min and then rises over 4 h, against 3 min for DTT -- which "
    "fixes the timescale contrast but not the dose axis",
    reason="the second entry point is the whole reason this block is worth building, and it "
           "cannot be given a number. generator/stress_panel.py currently gives DTT and "
           "tunicamycin the same {'UPR': 1.0} target, which is why they return identical "
           "answers; this module separates them structurally (u_tm enters the client influx, "
           "d_dtt the folding capacity) and still refuses the scale",
    missing="a Ki for tunicamycin against S. cerevisiae Alg7, or a yeast plate of the two "
            "stressors at matched UPR output"))

IRE1_DEACTIVATION = UPR_PARAMS.add(Param.refused(
    "ire1_deactivation_rate", "1/h",
    "REFUSED. Pincus Text S1, verbatim: 'There are no known measurements for g[I1A]'; the "
    "value he uses, 196 /s, is set equal to the BiP dissociation constant, itself fixed by "
    "requiring an Ire1-BiP ratio taken from MAMMALIAN IRE1-alpha co-IP. mech/state.py's "
    "TIMESCALE_CATALOGUE already carries the same refusal as ire1_cluster: van Anken 2014 "
    "images clusters at a fixed 45 min timepoint and reports no formation or dissolution rate",
    reason="no formation or dissolution rate for Ire1 clusters is reported in any yeast paper",
    missing="a time course of Ire1 focus formation and dissolution after a DTT step. "
            "Because 196 /s is a 5 ms complex lifetime -- a partition function, not an "
            "off-rate -- this module treats Ire1 activation as a rapid equilibrium and says "
            "so, rather than integrating a state whose rate constant does not exist"))

ERAD_RATE = UPR_PARAMS.add(Param.refused(
    "erad_rate", "1/h",
    "REFUSED. Searched: Sosa-Carrillo 2023 (PMID 37231013) fits an ERAD constant separately "
    "at every induction level, which is the signature of a missing mechanism rather than a "
    "measurement of one; pcSecYeast excludes the UPR and EGAD explicitly",
    reason="no ERAD rate constant for a defined ER client in S. cerevisiae",
    missing="a pulse-chase degradation rate for a misfolded ER client. The chain therefore "
            "sends all client loss through folding and dilution, which OVERSTATES the "
            "client pool and is the conservative direction for a block whose output rises "
            "with it"))

BASAL_HAC1_COPIES = UPR_PARAMS.add(Param.refused(
    "basal_hac1_copies", "molecules/cell",
    "REFUSED. Pincus Text S1 cites '[Walter] Walter, P., Personal Communication' for basal "
    "HAC1 mRNA = 200 copies/cell, and derives b[H1u_m] = 0.167 from it (0.167/8.33e-4 = "
    "200.5). Ho 2018 gives HAC1 PROTEIN 1,778-8,970 molecules/cell but no message count",
    reason="unpublished, and it sets the splicing substrate pool, so it is load-bearing",
    missing="total HAC1 mRNA by RT-qPCR against a standard curve. This chain is written in "
            "units of basal total HAC1, so the number never arises -- the refusal is "
            "recorded so the next reader does not reinstate it"))

IRE1_BIP_RATIO = UPR_PARAMS.add(Param.refused(
    "ire1_bip_ratio", "dimensionless",
    "REFUSED. Pincus pins his BiP dissociation constant by requiring R = c[I1.B]*[B]/g[I1.B] "
    "= 12.5, inferred from Bertolotti et al. 2000, Nat Cell Biol 2:326-32 -- MAMMALIAN "
    "IRE1-alpha co-immunoprecipitation -- and even there the inference (co-IP fraction equals "
    "equilibrium ratio) is an assumption. It is not self-consistent with his own figures "
    "either: 0.0350 * 57,851 / 196 = 10.3, not 12.5",
    reason="Pincus's own words are 'no real quantitative measurements have been recorded', "
           "and the inference chain behind the number runs through a different organism",
    missing="a Kar2-Ire1 dissociation constant in S. cerevisiae. The reduction eliminates "
            "the Ire1-BiP complex specifically so that this row never enters an equation"))

DEACTIVATION_MECHANISM = UPR_PARAMS.add(Param.refused(
    "deactivation_mechanism", "1/h",
    "REFUSED. The chain has no term that can terminate the response on the timescale the "
    "plates show, and this module will not manufacture one. What CAN be said without an "
    "approximation is the ablation: removing the chaperone pool moves the committed "
    "activity fold by 0.0029x the plate floor (upr_ablation()), so the one feedback the "
    "chain does carry is not doing the work. What CANNOT yet be said is how much faster the "
    "real termination is, and the honest table is in DEACTIVATION_VERSUS_RELAXATION: on ONE "
    "of the ten measurable (plate, construct, dose) blocks -- 20260722/UPRE1 at 1 mM -- the "
    "measured decay exceeds the fastest relaxation the chain permits even at the "
    "conservative corner, by 1.63x; on the other nine the mu bracket is wide enough that it "
    "does not, so the comparison is suggestive and not a refutation. Pincus 2010 reproduces the "
    "deactivation only by giving BiP a 2 h decay (0.5004 /h, sourced to Axelsen & Sneppen "
    "2004, a modelling paper) which the SILAC measurement this module uses contradicts by "
    "19.5x -- rate against rate, 0.5004 / 0.025672. See KAR2_PROTEIN_DECAY_PER_H: quoting "
    "27/2 = 13.5x compares a half-life with a mean lifetime",
    reason="the chain's own feedback is MEASURED not to be doing the work -- that is the "
           "ablation, and it needs no bracket -- while what IS doing it cannot be named yet, "
           "because mu(t) is not resolved on the per-construct sheets and nine of the ten "
           "blocks are therefore undecided. The leading candidate is that the DRIVER "
           "decays -- DTT is oxidised and consumed in aerobic medium -- which would make "
           "the termination a property of the medium rather than of the cell. Writing a "
           "rate here would be inventing the answer to the one question these plates can "
           "settle",
    missing="two measurements, either of which closes it. (1) Medium DTT against time in a "
            "dosed well -- an Ellman/DTNB read on spent medium at three timepoints -- which "
            "separates a decaying driver from a faster chaperone loss. (2) The OD600 traces "
            "alongside the per-construct sheets, which would replace the mu_max/mu_late "
            "bracket with the time-resolved rate the committed pipeline already computes "
            "and collapse the table below to single numbers"))

DEACTIVATION_VERSUS_RELAXATION = (
    ("20260722", "UPRE1", 0.5, (0.357, 0.471), (0.235, 0.528), 0.68),
    ("20260722", "UPRE1", 1.0, (0.777, 1.495), (0.186, 0.477), 1.63),
    ("20260722", "UPRE2", 0.5, (0.150, 0.206), (0.243, 0.430), 0.35),
    ("20260722", "UPRE2", 1.0, (0.415, 0.777), (0.178, 0.437), 0.95),
    ("20260804", "UPRE1", 0.5, (0.032, 0.130), (0.262, 0.879), 0.04),
    ("20260804", "UPRE1", 1.0, (0.089, 0.481), (0.207, 0.858), 0.10),
    ("20260804", "UPRE2", 0.5, (-0.006, 0.146), (0.269, 1.043), -0.01),
    ("20260804", "UPRE2", 1.0, (0.103, 0.359), (0.249, 0.894), 0.11),
    ("20260803", "UPRE1", 0.5, (0.175, 0.603), (0.410, 1.007), 0.17),
    ("20260803", "UPRE1", 1.0, (0.135, 0.578), (0.306, 0.850), 0.16),
)
"""MEASURED decay of the activity fold against the chain's fastest relaxation, per block.

Columns: plate, construct, dose, ``deactivation_rate()`` as ``(slow, fast)``,
``chaperone_relaxation_rate()`` as ``(slow, fast)``, and the CONSERVATIVE ratio
``deact_slow / relax_fast``. A value above 1 means the response terminates faster than the
chain can possibly terminate it, whatever the true mu is inside its bracket.

**One row of ten is above 1.** That is the whole honest content: the comparison is
suggestive at the top sub-lethal rung of the best-behaved plate and inconclusive everywhere
else, because the mu bracket is wide. It is tabulated in full rather than summarised to its
best row, and test_mech_upr.py re-derives every entry. The two 20260803/UPRE2 rows are
absent because that block never rises 2% above its control after 2 h, so there is no decay
to fit and :func:`deactivation_rate` raises instead of returning a number.
"""

REFUSED_FROM_PINCUS = (
    "basal HAC1 mRNA = 200 copies/cell (Walter, P., Personal Communication)",
    "Ire1-BiP ratio R = 12.5 (Bertolotti 2000, mammalian IRE1-alpha co-IP)",
    "a separate Ero1 state (2.91e6 molecules from the table against ~1.35e6 in the text, and "
    "300-800x above the 3,534-10,125 molecules/cell Ho 2018 measures)",
    "V_ER = 2.15 um^3 (uncited, and it divides every bimolecular constant)",
    "the U / U.B / Ud / Ud.B disulfide bookkeeping (4 states, 3 unmeasured constants, "
    "5 orders of stiffness, zero observables)",
    "the min([I1A],[H1u]) splicing term (non-differentiable; with Hu ~200 copies against a "
    "Korennykh Km of ~6,500 molecules the enzyme is far below saturation, so the correct "
    "limit is the bilinear one)",
    "the RS_m / RS splicing-reporter chain (generator/reporter.py already models immature "
    "and mature FP with its own maturation constant; this module emits a rate, not an RFU)",
    "every numeric value in the Stroberg 2018 parameter table",
)
"""What the source supplies and this module will not use, with the reason, in one place."""

DROPPED_TO_PASS_THE_GATE = (
    ("folding_rate", "the ER client pool is eliminated to quasi-steady state (tau/T = "
                     "0.048-0.080) and the ER dilution term is dropped, after which the "
                     "folding rate cancels from the algebra entirely"),
    ("kar2_mrna_decay", "the KAR2 message is eliminated to quasi-steady state (tau/T = "
                        "0.101) and the constant cancels against its own basal transcription"),
    ("hac1_isoform_decay_ratio", "after elimination the occupancy depends on the spliced "
                                 "fraction relative to its own full-splicing value, and the "
                                 "ratio scales both -- it is not identified and not needed"),
    ("ire1_deactivation_rate", "Ire1 is a declared rapid equilibrium, so the rate constant "
                               "that does not exist never enters an equation"),
    ("erad_rate", "no ERAD term is written at all, so the constant that does not exist is not needed"),
    ("tunicamycin_entry", "the committed ladder is DTT only; the second entry point stays "
                          "REFUSED and unfitted, and the block cannot be scored on it here"),
    ("basal_hac1_copies", "the chain is in units of basal total HAC1"),
    ("ire1_bip_ratio", "the Ire1-BiP complex is eliminated"),
    ("deactivation_mechanism", "the shipped shape predicts the late-window dose response "
                               "and makes no claim about the transient, so the rate it "
                               "would need is not one of its scalars -- it is what the "
                               "block asks for next"),
)
"""The nine rows the plate shape drops, and why each one goes. Criterion (e), shown."""


# --------------------------------------------------------------------------------------
# The targets
# --------------------------------------------------------------------------------------

_PLATE_ASSAY = REPORTER_ACTIVITY_FLOOR.assay

UPRE2_HELD_OUT = Target(
    name="upre2_held_out",
    observable="UPRE2's DTT ladder predicted from a UPRE1-only fit, with no refit",
    assay=_PLATE_ASSAY,
    noise_floor=float(REPORTER_ACTIVITY_FLOOR.noise_floor),
    units="relative CV",
    source="generator/panel_experiment.py::OBSERVED_ACTIVITY_CV. The two constructs share "
           "one Hac1 pool and one upre_basal_share by construction, so this is a prediction "
           "and not a transform with a free parameter in it",
)

ACTIVITY_TRANSIENT = Target(
    name="activity_transient",
    observable="the fall in the dose-to-control growth-corrected activity fold over the "
               "second half of the read, which an instantaneous dose term cannot produce "
               "at any setting",
    assay=_PLATE_ASSAY,
    noise_floor=float(REPORTER_ACTIVITY_FLOOR.noise_floor),
    units="relative CV",
    source="generator/panel_experiment.py::OBSERVED_ACTIVITY_CV. Pincus 2010 Fig 1E "
           "MEASURED the transience directly (dTR/dt at 1.5 and 2.2 mM DTT rises and "
           "returns toward baseline within 4 h and 2 h); on the committed 20260722 UPRE1 "
           "ladder the activity fold at 1 mM falls from 2.84 at 2.0 h to 1.09 at 4.0 h "
           "read at mu_late, 2.15 to 1.25 read at mu_max. THE BLOCK CURRENTLY MISSES THIS "
           "TARGET and says so: see DEACTIVATION_MECHANISM. It is counted in the gate "
           "because it constrains the piece, not because the piece satisfies it",
)

PLATE_REPLICATE_HELD_OUT = Target(
    name="plate_replicate_held_out",
    observable="a second plate's ladder at its own MEASURED growth rate, from the same fit",
    assay=_PLATE_ASSAY,
    noise_floor=float(REPORTER_ACTIVITY_FLOOR.noise_floor),
    units="relative CV",
    source="generator/panel_experiment.py::OBSERVED_ACTIVITY_CV. mu enters every dilution "
           "term and is MEASURED per well off the OD600 series, never fitted. What makes "
           "these blocks held out is that they are different wells, NOT that mu differs: "
           "MEASURED, sweeping mu across the full committed range (0.204 to 0.405 /h) moves "
           "the predicted fold by at most 0.0155 relative, 0.106x the assay floor, and by "
           "under 0.001x at the two top rungs, because the occupancy is saturated there. An "
           "earlier draft claimed the differing growth rate was what made the plate a "
           "prediction rather than a replicate; it is not, and test_mech_upr.py measures it",
)

HAC1_DEACTIVATION_TIME = Target(
    name="hac1_deactivation_time",
    observable="time for HAC1 splicing to fall back toward baseline after a sub-lethal DTT step",
    assay="HAC1 splicing time course by Northern blot and splicing-reporter flow cytometry, "
          "W303a in 2x SDC (Pincus et al. 2010, PMID 20625545, Fig 1E)",
    noise_floor=2.0 / 3.0,
    units="relative CV",
    source="mech/state.py::TIMESCALE_CATALOGUE carries this as the MEASURED interval 2 h at "
           "1.5 mM DTT to 4 h at 2.2 mM. The floor is the width of that interval over its "
           "midpoint, (4 - 2)/3 = 0.667, because the measurement itself does not resolve a "
           "deactivation time more finely than the two doses it spans. It is deliberately "
           "NOT the 0.146 plate CV: that floor belongs to a different instrument",
)

UPRE1_LADDER_FITTED = Target(
    name="upre1_dose_ladder",
    observable="the UPRE1 sub-lethal DTT ladder the four scalars are fitted to",
    assay=_PLATE_ASSAY,
    noise_floor=float(REPORTER_ACTIVITY_FLOOR.noise_floor),
    units="relative CV",
    source="generator/panel_experiment.py::OBSERVED_ACTIVITY_CV. Declared FITTED and "
           "therefore not counted in the gate's denominator",
    fitted=True,
)

UPR_TARGETS = (UPRE2_HELD_OUT, ACTIVITY_TRANSIENT, PLATE_REPLICATE_HELD_OUT,
               HAC1_DEACTIVATION_TIME, UPRE1_LADDER_FITTED)
"""Four independent targets and one declared-fitted one."""

for _t in UPR_TARGETS:
    UPR_PARAMS.add_target(_t)
    PLATE_PARAMS.add_target(_t)


# --------------------------------------------------------------------------------------
# The entry function, refused
# --------------------------------------------------------------------------------------

class EntryRefused(RefusedValue):
    """A caller asked for a load from a dose without declaring how the two are related."""


@dataclass(frozen=True)
class DoseEntry:
    """A DECLARED dose axis, which is the only way past :data:`DTT_ENTRY`'s refusal.

    The refusal is the point and this class does not soften it. What it does is force the
    relationship between a millimolar dose and an ER load to be named, tagged and carried
    with the prediction, so that a reader can see it is an axis rather than a measurement.
    Constructing one with ``declared=False`` raises.

    Args:
        top_load: The client load, as a fraction of folding capacity, at ``top_dose_mM``.
            The load at intermediate doses is linear in dose between the basal load and
            this, which is the same one-scalar family as
            ``d_dtt = 1 - 1/(1 + D*dose)`` and is written this way so that the admissibility
            bound ``load < 1`` is a box constraint instead of a discovered singularity.
            SAY PLAINLY WHAT IS AND IS NOT REFUSED HERE. :data:`DTT_ENTRY` refuses the
            SCALE; the FORM is ASSERTED and this parameterisation does not hide it. Through
            :meth:`UprChain.simulate` the capacity loss is
            ``1 - basal_load/load``, which for a load linear in dose is EXACTLY the
            Michaelis-Menten ``D*dose/(1 + D*dose)`` with ``D = (top_load - basal_load) /
            (basal_load * top_dose_mM)`` -- verified identical to machine precision in
            test_mech_upr.py, and at the fitted axis ``D = 3.08 /mM``, i.e. the incumbent
            ``generator/stress_panel.py::module_response``'s ``dose/(potency + dose)`` at a
            potency of 0.324 mM. No admissible setting makes it something other than
            saturating. That is the residue of ARCHITECTURE_GAPS.md 0.5 this block does not
            remove, and it is why the block claims SHAPE and refuses the scale rather than
            claiming to have replaced the entry function with a mechanism.
        top_dose_mM: The dose ``top_load`` refers to. Must be at or below
            :data:`DTT_LETHAL_MM`, which is MEASURED on this project's plates.
        declared: Must be True. It exists so that the sentence "this is a declared axis, not
            a measurement" has to be typed at every call site.
        note: What the caller is claiming, for the record.
    """

    top_load: float
    top_dose_mM: float
    declared: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        if not self.declared:
            raise EntryRefused(
                "the mM-to-ER-load mapping for DTT is REFUSED -- there is no published "
                "mapping from medium DTT to lumenal DTT in S. cerevisiae, and writing a "
                "saturating function here is the incumbent Hill with a new name. Pass "
                "declared=True to state on the record that this is a DECLARED SENSITIVITY "
                "AXIS and that the block predicts SHAPE, not absolute dose. "
                f"{DTT_ENTRY.provenance()}")
        if not 0.0 < self.top_load < 1.0:
            raise ValueError(
                f"top_load must lie in (0, 1) -- it is the fraction of folding capacity the "
                f"client influx demands, and at 1 the machinery is exactly saturated and the "
                f"client pool diverges. Got {self.top_load!r}")
        if not 0.0 < self.top_dose_mM <= float(DTT_LETHAL_MM):
            raise ValueError(
                f"top_dose_mM must lie in (0, {float(DTT_LETHAL_MM):g}] mM. Above the "
                f"MEASURED lethal dose the growth rate has halved and the block is "
                f"inadmissible; got {self.top_dose_mM!r}")

    def client_load(self, dose_mM: float, basal_load: float) -> float:
        """Client load as a fraction of folding capacity, at one dose.

        Named ``client_load`` and not ``load``: tests/test_call_signatures.py resolves calls
        by bare name across the package, so a method called ``load`` on any library class
        makes every ``np.load(..., allow_pickle=...)`` in the repository look like a stale
        keyword. One collision, three false positives, and the fix belongs here.
        """
        if dose_mM < 0:
            raise ValueError(f"dose must be non-negative, got {dose_mM}")
        slope = (self.top_load - basal_load) / self.top_dose_mM
        return basal_load + slope * float(dose_mM)


# --------------------------------------------------------------------------------------
# The states
# --------------------------------------------------------------------------------------

UPR_STATES = MechState((
    StateVar(
        name="er_client", units="multiples of the chaperone-client Kd",
        tau_h=(12.0 / 60.0, 20.0 / 60.0),
        source="ASSERTED. The relaxation time of the unfolded pool is the folding time, "
               "which Pincus asserts at 20 min uncited and Stroberg 2018 asserts at the "
               "same 20 min from a mammalian review. The interval is Pincus's own "
               "sensitivity range, 12-20 min",
        sweep_axis=True,
        note="Under DTT the effective relaxation is 1/((1-d)k_f + mu), so the elimination "
             "the ratio licenses at basal is licensed less strongly at dose; "
             "reduction_bias() measures the consequence instead of arguing about it.",
    ),
    StateVar(
        name="hac1_unspliced", units="fraction of basal total HAC1",
        tau_h=(1.0 / (60.0 / 20.0 + 60.0 / 11.0), 1.0 / (60.0 / 20.0)),
        source="MEASURED. Sidrauski 1996, PMID 8898194: 20 min HAC1 mRNA lifetime. The fast "
               "end adds the splicing flux at full Ire1 activation (11 min, Kawahara 1997, "
               "PMID 9348528), which is the only other loss this species has",
        constrained_by="HAC1i exon-junction RT-qPCR on the UPRE1 strain (ystwin/qpcr.py)",
        note="The splicing fraction it reports is near-binary in vivo, so this assay is a "
             "pathway control -- it says the switch flipped, not how far.",
    ),
    StateVar(
        name="hac1_spliced", units="fraction of basal total HAC1",
        tau_h=(1.0 / (60.0 / 20.0), 1.0 / (60.0 / 20.0)),
        source="MEASURED. Sidrauski 1996, PMID 8898194, 20 min, applied to the spliced "
               "isoform on Pincus's assertion that the two decay alike -- an assertion "
               "Cherry 2019, PMID 30874502, contradicts",
        constrained_by="HAC1i exon-junction RT-qPCR on the UPRE1 strain (ystwin/qpcr.py)",
    ),
    StateVar(
        name="hac1_protein", units="multiples of the UPRE-2 Kd",
        tau_h=(1.0 / (_LN2 / (1.0 / 60.0)), 1.0 / (_LN2 / (1.5 / 60.0))),
        source="MEASURED. Pal 2007, PMID 17108329: Hac1p half-life 1-1.5 min in W303 by "
               "cycloheximide shutoff, so the mean lifetime is 1.44-2.16 min",
        constrained_by="mCitrine RFU on the UPRE1/UPRE2 plates, and only through the "
                       "occupancy theta = P/(K_P + P) -- the pool itself is not measured "
                       "anywhere in this project",
        note="This is the state the architecture had deleted. Restoring it costs two "
             "MEASURED constants and removes Pincus's two fitted ones (a0 = 296.5, a1 = 5.26).",
    ),
    StateVar(
        name="kar2_mrna", units="fold over basal",
        tau_h=(1.0 / (60.0 / 25.0), 1.0 / (60.0 / 25.0)),
        source="ASSERTED. Pincus Text S1 g[Bm] = 6.67e-4 /s (25 min), no citation",
        constrained_by="KAR2 RT-qPCR (ystwin/qpcr.py) -- the measurement UPR_ANCHOR ranks "
                       "first, and it has not been run",
    ),
    StateVar(
        name="chaperone", units="fold over basal",
        tau_growth_multiple=1.0,
        source="MEASURED. The pool is lost at mu + k_deg with k_deg = ln2/27.0 h "
               "(Christiano 2014, PMID 25466257, dilution-corrected). 1/mu alone is used "
               "here because StateVar carries a growth multiple rather than a sum; it "
               "OVERSTATES tau by 9% at mu = 0.28 /h and 25% at the fed-batch setpoint, "
               "which makes the state look slower and so is the conservative direction for "
               "a reduction test",
        constrained_by="KAR2 RT-qPCR (ystwin/qpcr.py), and indirectly the mCitrine RFU "
                       "trajectory, whose fall over the second half of the read is this "
                       "pool's negative feedback",
        note="On a 4.14 h plate read this is the ONLY one of the six that survives the "
             "tau/T audit.",
    ),
))
"""The six integrated states of the reference chain, in integration order.

Ire1 is deliberately absent. `mech/state.py`'s catalogue already carries ``ire1_cluster``
with a :class:`~ystwin.mech.state.RefusedTimescale` -- van Anken 2014 images clusters at a
fixed 45 min timepoint and reports no rate -- and Pincus's own deactivation constant is
196 /s, a 5 ms complex lifetime, which is a partition function rather than an off-rate. So
Ire1 activation is a DECLARED rapid equilibrium here: algebraic, cooperative, with the Hill
coefficient bracketed by Korennykh's in-vitro measurement. That is a declaration, not
something the ratio test licensed, and it is on the record in :data:`IRE1_DEACTIVATION`.
"""


RESPONSE_DEACTIVATION_H = {"er_client": 2.0, "hac1_unspliced": 2.0, "hac1_spliced": 2.0}
"""Pincus's MEASURED 2 h deactivation at 1.5 mM DTT, as a ``driver_persists_h`` argument.

Passing it is a DIFFERENT claim from the one this block makes and the difference matters, so
it is offered rather than applied. `mech/integrate.py`'s hook shortens the window when a
state's DRIVER vanishes inside the run. Here the driver is the DTT-imposed folding deficit
and, so far as anything measured says, DTT does not go anywhere: what Pincus measured
deactivating within 2 h is the RESPONSE. What shuts it off is exactly what
:data:`DEACTIVATION_MECHANISM` refuses to name -- this module's own ablation says it is not
the chaperone feedback. Scoring the ratio against 2 h
therefore asks whether the states would be reducible in a world where the stressor was
washed out, which is Pincus's washout experiment and not this project's plate. Under it the
two HAC1 species and the ER client stop being reducible -- run
``upr_reduction(driver_persists_h=RESPONSE_DEACTIVATION_H)`` to see it -- and that is the
honest upper bound on how much the reduction is resting on the ratio.
"""


def upr_reduction(window: Window = PLATE_READ_4H, *,
                  driver_persists_h: Mapping[str, float] | None = None) -> ReducedSystem:
    """The tau/T audit of :data:`UPR_STATES` against one window, through `mech/integrate.py`.

    With nothing declared the whole window is used, which is right here: the driver is the
    DTT-imposed folding deficit and it persists for the entire read. See
    :data:`RESPONSE_DEACTIVATION_H` for the stricter alternative and why it answers a
    different question.
    """
    return reduce_for_window(UPR_STATES, window,
                             driver_persists_h=driver_persists_h or None)


# --------------------------------------------------------------------------------------
# The chain
# --------------------------------------------------------------------------------------

_GAMMA_H = float(HAC1_MRNA_DECAY_PER_H)
_K_SPL = float(HAC1_SPLICING_PER_H)
_GAMMA_P = float(HAC1_PROTEIN_DECAY_PER_H)
_K_DEG_C = float(KAR2_PROTEIN_DECAY_PER_H)
_GAMMA_CM = float(KAR2_MRNA_DECAY_PER_H)
_ISOFORM = float(HAC1_ISOFORM_DECAY_RATIO)
_P_MAX = float(HAC1_MAX_OCCUPANCY) / (1.0 - float(HAC1_MAX_OCCUPANCY))




def _ire1_active(client: float, chaperone: float, hill_n: float) -> float:
    """Active Ire1 fraction: the Kar2 titration, as a cooperative rapid equilibrium.

    Free chaperone is the total divided by ``1 + client`` -- the rapid-equilibrium partition
    of a pool that binds client with the Kd the client is measured in -- so the signal Ire1
    reads is ``client / free chaperone`` and rises both when client rises and when chaperone
    is consumed by it. That is the titration mechanism, with no Ire1-BiP complex and no
    :data:`IRE1_BIP_RATIO` in it.
    """
    free = max(chaperone, 1e-12) / (1.0 + client)
    load = client / free
    powered = load ** hill_n
    return powered / (1.0 + powered)


def _spliced_fraction(active: float, mu: float) -> float:
    """Spliced HAC1 at quasi-steady state, relative to its value at full Ire1 activation.

    Both HAC1 species eliminate on any window this project runs, and in the ratio the basal
    transcription rate, the total pool size and :data:`HAC1_ISOFORM_DECAY_RATIO` all cancel.
    That cancellation is what makes the unpublished 200 copies/cell never arise.
    """
    return ((active / (mu + _GAMMA_H + _K_SPL * active))
            * (mu + _GAMMA_H + _K_SPL))


def _occupancy(spliced_relative: float) -> float:
    """UPRE occupancy: the single-site Langmuir Fordyce 2012 measured, nothing else."""
    q = _P_MAX * spliced_relative
    return q / (1.0 + q)


@dataclass(frozen=True)
class UprChain:
    """The six-state reference chain. Integrated, not reduced, through `mech/integrate.py`.

    This is the object the plate shape is measured AGAINST, not the object shipped into a
    prediction: it carries :data:`FOLDING_RATE_PER_H`, :data:`KAR2_MRNA_DECAY_PER_H` and
    :data:`HAC1_ISOFORM_DECAY_RATIO`, none of which is pinned by a measurement, and
    :data:`UPR_PARAMS` fails ``require_gate()`` because of them.

    Args:
        basal_load: Unstressed client influx as a fraction of folding capacity.
        hill_n: Ire1 cooperativity.
        basal_share: UPRE-independent share of the promoter rate, in occupancy units.
        folding_rate_per_h: A point on :data:`FOLDING_RATE_PER_H`'s swept axis.
        mu_per_h: MEASURED specific growth rate for the well being simulated.
    """

    basal_load: float
    hill_n: float
    basal_share: float
    folding_rate_per_h: float
    mu_per_h: float

    def __post_init__(self) -> None:
        if not 0.0 < self.basal_load < 1.0:
            raise ValueError(f"basal_load must lie in (0, 1), got {self.basal_load}")
        if self.mu_per_h <= 0:
            raise ValueError(
                f"mu must be positive and MEASURED off the biomass series, got {self.mu_per_h}")
        FOLDING_RATE_PER_H.at(self.folding_rate_per_h)
        low, high = IRE1_HILL_N.bounds
        if not low <= self.hill_n <= high:
            raise ValueError(
                f"hill_n = {self.hill_n} is outside Korennykh 2009's measured bracket "
                f"[{low:g}, {high:g}] for yeast Ire1 RNase activation (PMID 19079236). "
                f"{IRE1_HILL_N.provenance()}")

    def _client_steady(self, capacity_loss: float, chaperone: float,
                       client_influx: float) -> float:
        """Client pool at steady state, WITH the ER dilution term the plate shape drops."""
        k_f = self.folding_rate_per_h
        def imbalance(u: float) -> float:
            return (k_f * client_influx
                    - k_f * (1.0 - capacity_loss) * chaperone * u / (1.0 + u)
                    - self.mu_per_h * u)
        if imbalance(1e9) > 0:
            return 1e9
        return brentq(imbalance, 1e-12, 1e9, xtol=1e-14, rtol=1e-13)

    def basal(self) -> dict[str, float]:
        """The unstressed steady state, and the basal transcription rates it implies."""
        mu = self.mu_per_h
        u0 = self._client_steady(0.0, 1.0, self.basal_load)
        a0 = _ire1_active(u0, 1.0, self.hill_n)
        hu0 = 1.0 / (1.0 + _K_SPL * a0 / (mu + _GAMMA_H * _ISOFORM))
        hs0 = 1.0 - hu0
        theta0 = _occupancy(_spliced_fraction(a0, mu))
        p0 = _P_MAX * _spliced_fraction(a0, mu)
        return {
            "er_client": u0, "hac1_unspliced": hu0, "hac1_spliced": hs0,
            "hac1_protein": p0, "kar2_mrna": 1.0, "chaperone": 1.0,
            "theta0": theta0, "active0": a0,
            "beta_hac1": (mu + _GAMMA_H) * hu0 + _K_SPL * a0 * hu0,
            "hac1_spliced_full": (_K_SPL * ((mu + _GAMMA_H) * hu0 + _K_SPL * a0 * hu0)
                                  / ((mu + _GAMMA_H + _K_SPL)
                                     * (mu + _GAMMA_H * _ISOFORM))),
            "beta_kar2_mrna": (mu + _GAMMA_CM) / (self.basal_share + theta0),
            "beta_chaperone": mu + _K_DEG_C,
        }

    def rhs(self, client_influx: float, capacity_loss: float, *,
            growth_rate_per_h: float | None = None,
            synthesis_scale: float = 1.0, folding_scale: float = 1.0):
        """``f(t, y)`` over :data:`UPR_STATES`, for one dose. Time in hours."""
        mu = self.mu_per_h if growth_rate_per_h is None else float(growth_rate_per_h)
        if not math.isfinite(mu) or mu < 0:
            raise ValueError("current growth_rate_per_h must be finite and nonnegative")
        if not math.isfinite(client_influx) or client_influx < 0:
            raise ValueError("client_influx must be finite and nonnegative")
        if not all(math.isfinite(x) and 0 <= x <= 1
                   for x in (capacity_loss, synthesis_scale, folding_scale)):
            raise ValueError("capacity_loss, synthesis_scale and folding_scale must lie in [0, 1]")
        k_f = self.folding_rate_per_h
        b = self.basal_share
        k = self.basal()

        def f(_t, y):
            u, hu, hs, p, cm, c = np.maximum(np.asarray(y, dtype=float), 0.0)
            a = _ire1_active(u, max(c, 1e-9), self.hill_n)
            theta = p / (1.0 + p)
            return [
                k_f * client_influx - k_f * folding_scale * (1.0 - capacity_loss) * c * u / (1.0 + u) - mu * u,
                synthesis_scale * k["beta_hac1"] - (mu + _GAMMA_H) * hu - _K_SPL * a * hu,
                _K_SPL * a * hu - (mu + _GAMMA_H * _ISOFORM) * hs,
                synthesis_scale * _P_MAX * (self.mu_per_h + _GAMMA_P) * hs / k["hac1_spliced_full"]
                - (mu + _GAMMA_P) * p,
                synthesis_scale * k["beta_kar2_mrna"] * (b + theta) - (mu + _GAMMA_CM) * cm,
                synthesis_scale * k["beta_chaperone"] * cm - (mu + _K_DEG_C) * c,
            ]

        return f

    def simulate(self, entry: DoseEntry, dose_mM: float,
                 window: Window = PLATE_READ_4H, n_points: int = 84) -> Trajectory:
        """Integrate one dose across one window with the pinned stiff driver."""
        if dose_mM > float(DTT_LETHAL_MM):
            raise ValueError(
                f"{dose_mM:g} mM is above the MEASURED lethal dose of "
                f"{float(DTT_LETHAL_MM):g} mM, where growth has fallen to half of control. "
                f"The block is inadmissible there and will not pretend otherwise")
        load = entry.client_load(dose_mM, self.basal_load)
        capacity_loss = 0.0 if load <= self.basal_load else 1.0 - self.basal_load / load
        start = self.basal()
        y0 = {name: start[name] for name in UPR_STATES.names}
        return integrate_window(self.rhs(self.basal_load, capacity_loss), y0, window,
                                UPR_STATES, n_points=n_points)

    @staticmethod
    def occupancy(trajectory: Trajectory) -> np.ndarray:
        """UPRE occupancy along a trajectory: the block's output into the reporter."""
        p = trajectory.of("hac1_protein")
        return p / (1.0 + p)


# --------------------------------------------------------------------------------------
# The plate observable: the repository's own growth-corrected activity
# --------------------------------------------------------------------------------------

DOSE_LADDER_MM = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
"""The committed DTT ladder on the four NewProtocol exports."""

SUBLETHAL_DOSE_MM = tuple(d for d in DOSE_LADDER_MM if d <= float(DTT_LETHAL_MM))
"""The rungs at or below the MEASURED lethal dose. The block refuses the other two."""

_PLATE_KINETICS = ReporterKinetics(k_deg=0.0)
"""The reporter kinetics scripts/run_sensor_characterisation.py uses, reused verbatim.

``k_deg = 0`` is that script's own D2 finding -- mCitrine loss is pure dilution -- and
``k_mat`` is left unset because it is UNMEASURED on this instrument. `REVISED_BUILD_LIST.md`
experiment 1.2 is the plate that would supply it, and it is listed as a precondition for
this phase. Setting it here would be inventing it, and it would differ between this module
and the committed file, which is worse than leaving it out of both.
"""

LATE_WINDOW_FRACTION = 0.75
"""Where the late window starts, as a fraction of the read. A READING CONVENTION.

Not a model parameter and not in either registry, for the same reason
:data:`INTERNAL_GRID_POINTS` is not: it is a definition the observable already carries.
``scripts/run_sensor_characterisation.py::LATE_WINDOW_FRACTION = 0.75`` is what produced the
``activity_late`` column this block is scored on, so the model has to average over the same
window or it is predicting a different number. One difference is stated rather than hidden:
that script takes the last quarter of the GROWTH WINDOW, after ``growth_window()`` trims the
low-density head, and this module takes the last quarter of the read, because the
per-construct sheets carry no OD. test_mech_upr.py MEASURES what that is worth by sweeping
the start over 0.60-0.85 and checking the predicted fold moves by less than a twentieth of
the assay floor -- MEASURED at 0.0058 relative, 4.0% of it, across 0.60 to 0.85.
"""

_PLATE_COLUMNS = ("dose_mM", "signal", "mu_per_h", "mu_control_per_h")
"""What a scored block must carry, and nothing else.

``read_h`` used to be here and is gone. No equation read it: the window the late mean is
taken over is :attr:`PlateUpr.window`, a declared constant, and the read duration entered
only through ``plate/replay.py``, which carries per-construct sheets for three of the four
committed plates. So a column the model never used was excluding the whole 20260728
replicate -- twelve UPRE1 and twelve UPRE2 sub-lethal wells whose ``activity_late`` is
committed in the same file as everybody else's. :func:`read_duration_h` still exists and
:func:`score_dtt_ladder`'s tests still check it where it is knowable: all three sheets end
at 4.0 h against the declared 4.14 h window, and that 3.5% moves the predicted fold by
0.0012 relative, 0.84% of the assay floor.
"""

PLATE_STATE = MechState((UPR_STATES["chaperone"],))
"""What the plate shape integrates: the one state the tau/T audit keeps, and no other.

Built from :data:`UPR_STATES` rather than declared separately, so the shape cannot drift
away from the audit that licensed it -- test_mech_upr.py asserts it equals
``upr_reduction(PLATE_READ_4H).integrated``.
"""

STARTING_POINTS = (
    (0.05, 0.25, 3.0, 0.0),
    (0.20, 0.90, 6.0, 0.05),
)
"""Two DECLARED, fixed starting points for the four-scalar fit, best residual wins.

Not a search over starting points that nobody counted -- the set is a literal, it never
changes with the data, and it is here rather than inside ``refit`` so a reader can count it.

It earns its keep on this ladder, MEASURED: the first point converges to
``(0.0060, 0.0108, 5.140, 0.000)`` at relative RMS 0.0620 and the second to
``(0.2424, 0.9900, 8.000, 1.809)`` at 0.0374. The gap between the two minima is 0.0246, 0.17x
the assay floor, so the basin matters; and the two parameter sets are nowhere near each
other, which is the same non-identifiability :class:`PlateUpr` records. What does NOT move is
the prediction: the two optima's held-out scores differ by at most 0.017 relative RMS across
the seven held-out blocks.
"""

CHAIN_QUADRATURE_POINTS = 801
"""Output grid the six-state chain's occupancy is reported on inside :func:`reduction_bias`.

ASSERTED as a numerical setting, and pinned by measurement in test_mech_upr.py: halving it
moves the reported bias by less than a thousandth of the assay floor. It is finer than
:data:`INTERNAL_GRID_POINTS` because the chain's occupancy is a pulse and the shape's is
close to flat over the window that is averaged."""

INTERNAL_GRID_POINTS = 97
"""Output grid the chaperone pool is integrated on.

ASSERTED as a numerical setting, not a model parameter: test_mech_upr.py pins it by
quadrupling it and measuring that the predicted fold does not move. It is a fixed grid
rather than one derived from the requested rows because the rows of the fold observable
carry no time at all -- the window is a property of the plate, not of the request."""


def _plate_columns(data: pd.DataFrame) -> None:
    missing = [c for c in _PLATE_COLUMNS if c not in data.columns]
    if missing:
        raise KeyError(
            f"a UPR plate block needs {', '.join(_PLATE_COLUMNS)}; missing {missing}. "
            f"mu is MEASURED per dose off the OD600 series and travels with the rows -- it "
            f"is never a fitted scalar here")


def read_duration_h(export: str) -> float:
    """How long the plate ran, from the committed per-construct trace.

    NOT an input to any prediction -- see :data:`_PLATE_COLUMNS` for why it stopped being
    one. It exists so the declared window can be checked against the plates that can check
    it, and it raises for a plate with no committed dose-response sheet rather than
    substituting the declared duration for a measured one.
    """
    return float(replay.load_doses(export).time_h.max())


def _characterisation(export: str, construct: str) -> pd.DataFrame:
    """The committed per-well DTT rows of outputs/sensor_characterisation.csv."""
    from ystwin import paths

    table = pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")
    stamp = export.split("_")[0]
    rows = table[(table.plate.str.startswith(stamp)) & (table.construct == construct)
                 & (table.stressor == "DTT")]
    if rows.empty:
        raise KeyError(
            f"outputs/sensor_characterisation.csv carries no DTT rows for {stamp}/"
            f"{construct}. Every number this block is scored on comes from that file, so an "
            f"absent block is refused rather than reconstructed")
    return rows


def _measured_growth_rates(export: str, construct: str) -> dict[float, float]:
    """Per-dose specific growth rate, MEASURED, from outputs/sensor_characterisation.csv."""
    rows = _characterisation(export, construct)
    return {float(k): float(v) for k, v in rows.groupby("dose_mM").mu_late.mean().items()}


def upre_activity_fold(export: str, construct: str) -> pd.DataFrame:
    """The SCORED observable: growth-corrected promoter activity over its matched control.

    One row per well per sub-lethal dose, straight out of the ``activity_late`` column of
    ``outputs/sensor_characterisation.csv`` -- this repository's own canonical quantity,
    produced by ``scripts/run_sensor_characterisation.py`` from the total-signal route with a
    time-resolved mu and the real OD600 traces. Nothing here re-derives it, so the block is
    scored against the committed number rather than against a reconstruction of it.

    ``signal`` is ``activity_late(well) / mean activity_late(matched zero-dose wells)``. The
    RFU-per-unit-activity conversion cancels in the fold, which is what lets the model carry
    four scalars and no amplitude.

    WHY NOT THE RAW DOSE-RESPONSE SHEET, and it matters because an earlier version of this
    module did exactly that. It fitted the accumulated increment
    ``[F(t,d) - F(0,d)] / [F(t,0) - F(0,0)]`` from ``plate/replay.py::load_doses`` under an
    ``exp(mu*t)`` biomass weight, on the stated ground that a TOTAL well fluorescence
    integrates ``X(s)k(s)``. Those sheets do not hold a total: ``plate/dose_response.py``
    says in its first line that they hold NORMALISED fluorescence, and the numbers agree --
    the late-window mean of that ``signal`` column reproduces ``naive_late``, the per-cell
    column of the committed file, to 1.3% at worst across every plate that commits a sheet,
    every construct and every dose. A per-cell signal obeys
    ``dS/dt = k - mu*S``, so its increment is not the biomass integral and the weight had the
    dilution term backwards. The consequence was not cosmetic: a dose-dependent dilution
    artefact (mu falls from 0.214 to 0.161 /h across this ladder) was being read as
    mechanism, and it inflated the ablation by 400x.

    Only the sub-lethal rungs are returned. The 2 and 5 mM wells sit above the MEASURED
    1.55 mM lethal dose, where growth has halved and the measured activity reverses; the
    block declares itself inadmissible there rather than fitting through it.
    """
    rows = _characterisation(export, construct)
    control = rows[rows.dose_mM == 0.0]
    if control.empty:
        raise KeyError(f"{export}/{construct} has no zero-dose control to normalise against")
    reference = float(control.activity_late.mean())
    dosed = rows[rows.dose_mM.isin([d for d in SUBLETHAL_DOSE_MM if d > 0.0])]
    if dosed.empty:
        raise ValueError(f"{export}/{construct} has no sub-lethal dosed wells")
    return pd.DataFrame({
        "dose_mM": dosed.dose_mM.astype(float).to_numpy(),
        "signal": (dosed.activity_late.astype(float) / reference).to_numpy(),
        "mu_per_h": dosed.mu_late.astype(float).to_numpy(),
        "mu_control_per_h": np.full(len(dosed), float(control.mu_late.mean())),
        "well": dosed.well.to_numpy(),
    })


def upre_activity_trace(export: str, construct: str, *,
                        rate: str = "mu_late") -> pd.DataFrame:
    """The DIAGNOSTIC observable: the same activity fold, resolved in time.

    Recovered from the committed per-construct sheets with this repository's own estimator --
    ``reporter.promoter_activity`` at :data:`_PLATE_KINETICS`, the same call and the same
    kinetics ``scripts/run_sensor_characterisation.py`` makes -- so what comes back is
    ``dS/dt + mu*S`` on the per-cell signal those sheets carry.

    DIAGNOSTIC and not scored, because of one declared approximation: the script corrects
    with a time-resolved mu off the OD trace, the per-construct sheets carry no OD, and this
    function therefore applies one MEASURED scalar rate per dose. In the late window the two
    agree -- the late-window mean of what comes back reproduces ``activity_late_percell`` to
    within 1.4% on every plate that commits a sheet, every construct and every sub-lethal
    dose (2.4% including the two
    supra-lethal rungs the block refuses), a tenth of the assay floor -- and
    early in the read, where mu is still near ``mu_max``, they do not. That is why ``rate`` is
    an argument and why :func:`deactivation_rate` returns a BRACKET rather than a number.

    Args:
        rate: ``"mu_late"`` or ``"mu_max"``, both MEASURED columns of
            outputs/sensor_characterisation.csv. Neither is the true mu(t); the pair bracket
            it, and every claim this module makes off this trace is stated across both.
    """
    if rate not in ("mu_late", "mu_max"):
        raise ValueError(
            f"rate must name a MEASURED column of outputs/sensor_characterisation.csv -- "
            f"'mu_late' or 'mu_max'; got {rate!r}. A growth rate that is not on the plate is "
            f"not available to this module and will not be invented")
    frame = replay.load_doses(export)
    block = frame[frame.construct == construct]
    if block.empty:
        raise KeyError(
            f"{export} carries no construct {construct!r}; it has "
            f"{', '.join(sorted(pd.unique(frame.construct)))}")
    wide = block.pivot_table(index="time_h", columns="dose_mM", values="signal")
    times = wide.index.to_numpy(dtype=float)
    rates = {float(k): float(v) for k, v in
             _characterisation(export, construct).groupby("dose_mM")[rate].mean().items()}
    have = [d for d in SUBLETHAL_DOSE_MM if d in wide.columns and d in rates]
    if 0.0 not in have:
        raise KeyError(f"{export}/{construct} has no zero-dose control to normalise against")

    def activity(dose: float) -> np.ndarray:
        return promoter_activity(times, SpecificFluorescence(wide[dose].to_numpy(float)),
                                 np.full_like(times, rates[dose]), _PLATE_KINETICS)

    control = activity(0.0)
    rows = []
    for dose in have:
        if dose == 0.0:
            continue
        fold = activity(dose) / control
        for t, value in zip(times, fold):
            rows.append({"time_h": float(t), "dose_mM": float(dose),
                         "signal": float(value), "mu_per_h": rates[dose],
                         "mu_control_per_h": rates[0.0]})
    return pd.DataFrame(rows)


EXCESS_FLOOR = 0.02
"""Rows whose fold is within 2% of its control carry no decay and are dropped from the
log-linear fit in :func:`deactivation_rate`. ASSERTED as a reading convention; the fitted
rate moves by less than the bracket's own width across 0.01-0.05."""


def deactivation_rate(export: str, construct: str, dose_mM: float, *,
                      from_h: float | None = None) -> tuple[float, float]:
    """MEASURED first-order decay rate of the activity fold, as a bracket over mu.

    The fold's excess over its control is fitted log-linearly from ``from_h`` -- default half
    the read, which is the window :func:`upre_activity_trace`'s scalar-mu approximation is
    verified in -- to the end of the read.

    Returns:
        ``(slow, fast)``: the rate at ``mu_max`` and at ``mu_late``, in 1/h. The true rate
        lies between them because the true mu does. Every claim made off this number is made
        at the SLOW end, so the bracket cannot be doing the work.

    THIS IS THE NUMBER THAT REFUTES THE CHAPERONE FEEDBACK -- see
    :data:`DEACTIVATION_MECHANISM`.
    """
    out = []
    for rate in ("mu_max", "mu_late"):
        trace = upre_activity_trace(export, construct, rate=rate)
        rows = trace[np.isclose(trace.dose_mM, dose_mM)].sort_values("time_h")
        if rows.empty:
            raise KeyError(
                f"{export}/{construct} has no {dose_mM:g} mM rung below the lethal dose")
        floor = float(rows.time_h.max()) / 2.0 if from_h is None else float(from_h)
        rows = rows[rows.time_h >= floor]
        excess = rows.signal.to_numpy(dtype=float) - 1.0
        usable = excess > EXCESS_FLOOR
        if usable.sum() < 3:
            raise ValueError(
                f"{export}/{construct} at {dose_mM:g} mM never rises {EXCESS_FLOOR:.0%} "
                f"above its control after {floor:g} h, so there is no decay to fit and one "
                f"will not be manufactured")
        slope = float(np.polyfit(rows.time_h.to_numpy(dtype=float)[usable],
                                 np.log(excess[usable]), 1)[0])
        out.append(-slope)
    return out[0], out[1]


def chaperone_relaxation_rate(export: str, construct: str, dose_mM: float) -> tuple[float, float]:
    """The fastest deactivation the chain's MEASURED constants allow, at one well.

    The chaperone pool is lost at ``mu + k_deg`` and nothing in the chain is slower, so this
    is an upper bound on how fast the feedback can shut the response off. Returned as the
    same ``(slow, fast)`` bracket over the two MEASURED growth rates, so it is compared
    against :func:`deactivation_rate` corner for corner.
    """
    rows = _characterisation(export, construct)
    at = rows[np.isclose(rows.dose_mM, dose_mM)]
    if at.empty:
        raise KeyError(f"{export}/{construct} has no {dose_mM:g} mM rung")
    return (float(at.mu_late.mean()) + _K_DEG_C, float(at.mu_max.mean()) + _K_DEG_C)


# --------------------------------------------------------------------------------------
# The plate-window shape: one state, four scalars
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class PlateUpr(FittableModel):
    """The shape that passes the gate: one ODE state, four fitted scalars.

    Five of the six states of :class:`UprChain` are eliminated to algebra with the tau/T
    audit behind each cut (:func:`upr_reduction`) and the cost MEASURED rather than assumed
    (:func:`reduction_bias`). What survives is the chaperone pool, whose relaxation time
    ``1/(mu + k_deg)`` is MEASURED at both ends.

    The four scalars, and why each one is irreducibly free:

    * ``basal_load`` -- the unstressed ER client influx as a fraction of folding capacity.
      No yeast measurement exists (:data:`BASAL_FOLDING_LOAD`).
    * ``top_load`` -- the same quantity at 1 mM DTT, which IS the declared dose axis. The
      mapping is REFUSED (:data:`DTT_ENTRY`) and this is its one declared scalar.
    * ``hill_n`` -- Ire1 cooperativity, bracketed 2-8 in vitro (:data:`IRE1_HILL_N`).
    * ``basal_share`` -- the UPRE-independent share of the promoter rate, ONE value for both
      the KAR2 regulon and the reporter (:data:`UPRE_BASAL_SHARE`).

    Everything else in the chain is either MEASURED or cancels. The parameterisation by
    ``top_load`` rather than by a dose slope is what makes the admissibility bound
    (load < 1, or the folding machinery is exactly saturated and the client pool diverges) a
    box constraint the optimiser cannot walk out of.

    WHAT THE OUTPUT IS. The late-window mean of ``basal_share + theta``, over its own value
    in the matched control well -- the same functional of the same window that produced the
    ``activity_late`` column it is scored against. There is no amplitude and no ``R0``: the
    RFU-per-unit-activity conversion cancels in the fold, exactly.

    TWO OF THE FOUR SIT ON THE BOX, and it is stated here rather than left in the fitted
    dict. On the committed 20260722 UPRE1 fold the fit is ``basal_load`` = 0.2424,
    ``top_load`` = 0.9900, ``hill_n`` = 8.0000 and ``basal_share`` = 1.8093: the middle two
    are at their upper bounds. Twelve wells at four distinct doses do not identify four
    scalars, and this is what that looks like. Criterion (b) is still met and not by an
    excuse -- ``top_load`` IS the declared dose axis :data:`DTT_ENTRY` refuses to supply, and
    ``hill_n`` is bracketed by Korennykh 2009's in-vitro measurement -- so both are held by a
    declaration rather than by a fit. Neither should be quoted as though these rows measured
    it. What would fix it is the dose axis, not the model:
    ``outputs/sensor_characterisation.csv`` commits four sub-lethal rungs for this ladder.
    """

    name: str = "PlateUpr"
    parameter_names: tuple[str, ...] = ("basal_load", "top_load", "hill_n", "basal_share")
    response_column: str = "signal"
    feedback: bool = True
    top_dose_mM: float = 1.0
    window: Window = PLATE_READ_4H
    method: str = PINNED_METHOD

    def _rate(self, load: float, hill_n: float, mu: float):
        """Occupancy as a function of the chaperone pool, at one dose."""
        def theta_at(chaperone: float) -> float:
            c = max(chaperone, 1e-6)
            client = min(load / c, 1.0 - 1e-9)
            return _occupancy(_spliced_fraction(
                _ire1_active(client / (1.0 - client), c, hill_n), mu))
        return theta_at

    def _promoter_trace(self, load: float, hill_n: float, basal_share: float, mu: float,
                        theta0: float) -> tuple[np.ndarray, np.ndarray]:
        """Grid and ``basal_share + theta`` along it, with the chaperone pool integrated.

        With the feedback frozen there is no state left and the rate is a constant, which is
        the incumbent's algebra and is returned as such rather than integrated.

        The driver is `mech/integrate.py`'s :func:`~ystwin.mech.integrate.integrate_window`,
        not a local ``solve_ivp``: it supplies the pinned stiff method, its own fallback to
        LSODA, the tolerances and the non-negativity clip, and it takes :data:`PLATE_STATE`
        so the states integrated here are literally the ones the tau/T audit kept.
        ``method`` stays at :data:`~ystwin.mech.integrate.PINNED_METHOD` -- after five states
        go to algebra this is one non-stiff equation and LSODA agrees far below the assay
        floor, MEASURED in test_mech_upr.py, but nothing here costs enough to trade away the
        fallback the pinned default brings with it.
        """
        theta_at = self._rate(load, hill_n, mu)
        if not self.feedback:
            grid = np.linspace(0.0, self.window.duration_h, INTERNAL_GRID_POINTS)
            return grid, np.full(grid.shape, basal_share + theta_at(1.0))

        def rhs(_t, y):
            c = max(float(y[0]), 1e-6)
            return [(mu + _K_DEG_C)
                    * ((basal_share + theta_at(c)) / (basal_share + theta0) - c)]

        trajectory = integrate_window(rhs, {"chaperone": 1.0}, self.window, PLATE_STATE,
                                      n_points=INTERNAL_GRID_POINTS, method=self.method)
        return trajectory.t, np.array(
            [basal_share + theta_at(c) for c in trajectory.of("chaperone")])

    def _late_promoter_rate(self, load: float, hill_n: float, basal_share: float,
                            mu: float, theta0: float) -> float:
        """Mean promoter rate over the late window: what ``activity_late`` averages."""
        grid, trace = self._promoter_trace(load, hill_n, basal_share, mu, theta0)
        late = grid >= LATE_WINDOW_FRACTION * self.window.duration_h
        return float(np.mean(trace[late]))

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        _plate_columns(data)
        basal_load = float(parameters["basal_load"])
        top_load = float(parameters["top_load"])
        hill_n = float(parameters["hill_n"])
        basal_share = float(parameters["basal_share"])
        if not 0.0 < basal_load < top_load < 1.0:
            return np.full(len(data), np.inf)
        entry = DoseEntry(top_load=top_load, top_dose_mM=self.top_dose_mM, declared=True,
                          note="fitted axis; the mM-to-lumen mapping stays REFUSED")

        mu_control = float(np.asarray(data.mu_control_per_h, dtype=float)[0])
        theta0 = _occupancy(_spliced_fraction(
            _ire1_active(basal_load / (1.0 - basal_load), 1.0, hill_n), mu_control))
        control = self._late_promoter_rate(basal_load, hill_n, basal_share, mu_control,
                                           theta0)
        cache: dict[tuple[float, float], float] = {}
        out = []
        for dose, mu in zip(np.asarray(data.dose_mM, dtype=float),
                            np.asarray(data.mu_per_h, dtype=float)):
            key = (round(float(dose), 6), round(float(mu), 9))
            if key not in cache:
                load = min(entry.client_load(float(dose), basal_load), 1.0 - 1e-9)
                cache[key] = self._late_promoter_rate(load, hill_n, basal_share,
                                                      float(mu), theta0)
            out.append(cache[key] / max(control, 1e-12))
        return np.asarray(out, dtype=float)

    def bounds(self) -> tuple[list[float], list[float]]:
        """The box. Every edge is a declared axis or an admissibility bound, not a guess."""
        low = [1e-3, 1.1e-3, IRE1_HILL_N.bounds[0], UPRE_BASAL_SHARE.bounds[0]]
        high = [0.97, 0.99, IRE1_HILL_N.bounds[1], UPRE_BASAL_SHARE.bounds[1]]
        return low, high

    def refit(self, data: pd.DataFrame) -> Fit:
        """Least squares on the RELATIVE residual, because the declared noise is a CV.

        `generator/panel_experiment.py::OBSERVED_ACTIVITY_CV` is a coefficient of variation,
        so the measurement model on this observable is multiplicative and the correctly
        weighted residual is ``(predicted - observed)/observed``. It is also the quantity
        every score in this module reports, so the fit optimises what it is graded on
        instead of a differently weighted proxy.
        """
        _plate_columns(data)
        observed = np.asarray(data.signal, dtype=float)
        if not np.all(observed > 0):
            raise ValueError(
                "a fold of a growth-corrected activity must be positive; a non-positive row "
                "means the control activity was estimated at or below zero and the row is "
                "not usable")
        low, high = self.bounds()

        def residual(x):
            guess = dict(zip(self.parameter_names, (float(v) for v in x)))
            predicted = self.predict(data, guess)
            return (np.where(np.isfinite(predicted), predicted, 1e3) - observed) / observed

        best = None
        for start in STARTING_POINTS:
            solution = least_squares(residual, list(start), bounds=(low, high),
                                     xtol=1e-12, ftol=1e-12, max_nfev=3000)
            score = float(np.sqrt(np.mean(residual(solution.x) ** 2)))
            if best is None or score < best[0]:
                best = (score, solution)
        rms, solution = best
        fitted = dict(zip(self.parameter_names, (float(v) for v in solution.x)))
        return Fit(self.name, fitted, len(self.parameter_names), len(data), rms,
                   bool(solution.status > 0))


@dataclass(frozen=True)
class FrozenChaperoneUpr(PlateUpr):
    """:class:`PlateUpr` with the chaperone pool frozen: the piece criterion (a) removes.

    Everything else is identical -- the same measured constants, the same occupancy, the
    same four scalars, refit on the same rows. What goes is the ONE state, so the UPRE
    output becomes an instantaneous function of dose again. That is precisely the incumbent
    `generator/stress_panel.py::module_response` in mechanistic clothing, and it is what the
    ablation prices. On this project's own plates it prices at 0.0029x the floor: see
    :func:`upr_ablation`.
    """

    name: str = "FrozenChaperoneUpr"
    feedback: bool = False


# --------------------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------------------

def _relative_rms(predicted: np.ndarray, observed: np.ndarray) -> float:
    return float(np.sqrt(np.mean(((predicted - observed) / observed) ** 2)))


@dataclass(frozen=True)
class LadderScore:
    """One (plate, construct) ladder scored against a fit, with the floor it is scored on."""

    export: str
    construct: str
    held_out: bool
    n_rows: int
    growth_rate_per_h: float
    relative_rms: float
    floor: float = float(REPORTER_ACTIVITY_FLOOR.noise_floor)

    @property
    def clears_floor(self) -> bool:
        return self.relative_rms < self.floor

    def line(self) -> str:
        tag = "HELD OUT" if self.held_out else "fitted  "
        verdict = "clears" if self.clears_floor else "MISSES"
        return (f"  {tag} {self.export.split('_')[0]} {self.construct:<6s} n={self.n_rows:<4d}"
                f" mu0={self.growth_rate_per_h:.3f} /h   rel RMS {self.relative_rms:.4f}"
                f"  {verdict} {self.floor:g}")


def score_dtt_ladder(fit_export: str, fit_construct: str = "UPRE1",
                     held_out: Sequence[tuple[str, str]] = (),
                     model: PlateUpr | None = None) -> tuple[Fit, tuple[LadderScore, ...]]:
    """Fit one ladder, then predict the others with NO refit, and report every score.

    This is the only place in this build where a mechanistic prediction meets this project's
    own measured data, so it is reported whole: the held-out rows are scored with the fitted
    plate's parameters unchanged, and a miss is printed rather than dropped.

    MEASURED, fitting 20260722/UPRE1 and holding out the other SEVEN blocks: in sample
    0.0374, then 0.1417 (UPRE2, same plate), 0.0791 and 0.1410 (20260728), 0.0994 and 0.0747
    (20260803), 0.1562 and 0.1846 (20260804), against a 0.146 floor. FIVE of the seven
    held-out blocks clear it and two miss, and both that miss are the same plate. The
    in-sample number is four scalars on twelve wells at four doses and should not be read as
    a fit quality; the held-out column is the block's claim.

    The 20260728 pair joined the scored set when ``read_h`` left the frame. That plate
    commits no per-construct dose-response sheet, so a column no equation read was excluding
    twenty-four wells whose ``activity_late`` is in the same committed file as the rest --
    :data:`_PLATE_COLUMNS` records it. Both of its blocks clear the floor, which is a
    genuine addition to criterion (c) and not a re-scoring of anything already counted.

    WHAT THIS COLUMN DOES NOT SHOW, MEASURED. Fit ``1 + A*d^n/(K^n + d^n)`` -- a Hill in
    DOSE ALONE, three free scalars against this model's four, with no chaperone state, no
    measured Kd and no mu -- to the same twelve wells, and it gives A = 0.2557,
    K = 0.2036 mM, n = 3.170 and tracks this block across every one of the eight blocks to
    within **0.0033 relative RMS (0.022x the assay floor)**, scoring better on four of the
    seven held out. The shipped shape is itself close to that curve: a three-parameter Hill
    approximates :meth:`PlateUpr.predict` over the whole sub-lethal ladder to 0.0055
    relative RMS, **0.038x the floor**. So the held-out column below is evidence that a dose
    shape transfers across plates and constructs, NOT evidence that the mechanism behind it
    is the right one -- consistent with :func:`upr_ablation`, which finds the one ODE state
    worth 0.0029x the floor. Reported here rather than left for a reader to discover,
    because ARCHITECTURE_GAPS.md 0.5 names exactly this -- an instantaneous saturating
    function of dose wearing mechanism's clothes -- as the failure mode of the first
    architecture, and the honest position is that the plate cannot tell the two apart. What
    the block keeps over the Hill is not fit: it is two MEASURED binding constants where
    Pincus had two fitted ones, and a refusal where the incumbent has a free dose scale.

    The 20260803 export, which the earlier accumulation observable missed by 4.7x, is now the
    BEST-predicted plate. That is not a rescue by tuning: that plate's anomaly was in its raw
    accumulated control signal -- 565 RFU across the read against 2,458 on 20260722 while
    growing almost twice as fast -- and a growth-corrected activity fold does not inherit it.
    It is direct evidence that the observable, not the mechanism, was carrying that miss.
    """
    model = PlateUpr() if model is None else model
    data = upre_activity_fold(fit_export, fit_construct)
    fit = model.refit(data)
    scores = [LadderScore(fit_export, fit_construct, False, len(data),
                          float(data.mu_control_per_h.iloc[0]),
                          _relative_rms(model.predict(data, fit.parameters),
                                        np.asarray(data.signal, dtype=float)))]
    for export, construct in held_out:
        rows = upre_activity_fold(export, construct)
        scores.append(LadderScore(
            export, construct, True, len(rows), float(rows.mu_control_per_h.iloc[0]),
            _relative_rms(model.predict(rows, fit.parameters),
                          np.asarray(rows.signal, dtype=float))))
    return fit, tuple(scores)


def upr_ablation(export: str = "20260722_ER&OxidativeStress_NewProtocol_ANALYSED.xlsx",
                 construct: str = "UPRE1", *, dose_mM: float = 1.0) -> AblationResult:
    """Criterion (a): remove the chaperone pool, both models refit, scored on the plate floor.

    **IT FAILS, and the number is the block's most useful result.** The observable is the
    growth-corrected activity fold at the top sub-lethal rung -- the number a reader of this
    plate would quote -- and removing the one ODE state moves it by 0.00042 relative,
    **0.0029x the 0.146 plate CV**. Three hundred and forty times below the floor. Nor does
    the sign rescue it: across the seven held-out blocks the FROZEN model scores better on
    six of the seven, by at most 0.0023 relative RMS -- 1.6% of the floor.

    The failure needs no approximation of any kind: both models are refit on the committed
    ``activity_late`` column, and the summarised number is a fold of two entries in it.

    IT IS ALSO NOT AN ARTEFACT OF THE RUNG IT IS SCORED ON, which is worth checking because
    the default dose is exactly where :func:`reduction_bias` collapses -- at 1 mM both models
    sit against the same :data:`HAC1_MAX_OCCUPANCY`. Run at every sub-lethal rung the effect
    is 0.0062x the floor at 0.1 mM, 0.0039x at 0.2, 0.0041x at 0.5 and 0.0029x at 1.0: the
    state fails to clear by at least 161x everywhere on the ladder, and the worst case for
    the block is the BOTTOM rung, not the top. The reason a dose where the frozen model's
    promoter rate is visibly wrong in time still scores flat is that both models are REFIT --
    four scalars are ample to absorb the missing relaxation into the dose axis, which is the
    whole reason the protocol refits the reduced arm instead of freezing it.

    A SUGGESTION about why, which is deliberately not promoted to a reason. The pool's only
    losses are dilution and Christiano 2014's MEASURED 27.0 h degradation, so it cannot
    relax faster than ``mu + k_deg`` (:func:`chaperone_relaxation_rate`). On one of the ten
    measurable (plate, construct, dose) blocks the response decays faster than that even at
    the conservative corner of both brackets -- 1.63x, 20260722/UPRE1 at 1 mM -- and on the
    other nine the mu bracket is too wide to say. :data:`DEACTIVATION_VERSUS_RELAXATION` is
    that table in full, and :data:`DEACTIVATION_MECHANISM` is the refusal and the two
    measurements that would settle it.
    """
    data = upre_activity_fold(export, construct)
    at = data[np.isclose(data.dose_mM, dose_mM)].iloc[[0]][list(_PLATE_COLUMNS)]
    at = at.reset_index(drop=True)

    def summarise(model: FittableModel, fit: Fit) -> float:
        return float(model.predict(at, fit.parameters)[0])

    observable = Observable(
        name=REPORTER_ACTIVITY_FLOOR.name,
        units="growth-corrected activity over the matched zero-dose wells",
        assay=REPORTER_ACTIVITY_FLOOR.assay,
        scored="relative",
        data=data,
        summarise=summarise,
        description=(f"late-window activity fold at {dose_mM:g} mM DTT, from "
                     f"outputs/sensor_characterisation.csv"),
    )
    return ablate(PlateUpr(), FrozenChaperoneUpr(), observable)


def reduction_bias(entry: DoseEntry, doses: Sequence[float], *, mu_per_h: float,
                   hill_n: float = 4.5, basal_load: float = 0.242,
                   basal_share: float = 0.0,
                   folding_rate_per_h: float | None = None,
                   window: Window = PLATE_READ_4H) -> dict[float, float]:
    """Criterion (d), MEASURED on the observable: what the plate shape's reduction costs.

    The ratio test says the five fast states may go; it does not say by how much the answer
    moves, and `mech/integrate.py`'s own docstring is explicit that ``tau/T`` is a bound on
    the bias of a TIME-INTEGRAL and that a systematic bias is not the same kind of quantity
    as a random CV. So this integrates the six-state chain and the one-state shape at the
    same parameters and compares them on **the number the plate actually produces** -- the
    late-window mean promoter rate at a dose over the matched control -- rather than on an
    internal variable no assay touches.

    **The answer is that the reduction is free, and that is the finding.** At the fitted
    parameters and mu = 0.214 /h it runs 0.021 at 0.1 mM, 0.049 at 0.2 mM, 0.007 at 0.5 mM
    and 0.001 at 1.0 mM against a 0.146 floor -- at most **0.34x** it, so eliminating five of
    the six states costs less than the assay that would have to see it. It is deliberately
    NOT monotone in dose: it is a HUMP, peaking at 0.049 near 0.2 mM (0.0049 at 0.05 mM,
    0.00004 at 0.75 mM) and collapsing above it. What it tracks is dropping
    :data:`FOLDING_RATE_PER_H` -- the shape has no ER dilution term, so its client pool sits
    high and the cooperative Ire1 step raises that difference to the power ``hill_n`` -- and
    the collapse at the top is the two models saturating at the same
    :data:`HAC1_MAX_OCCUPANCY`, where a difference in client pool no longer moves the
    occupancy.

    Read it beside :func:`upr_ablation` and the pair say more than either alone. The five
    states the ratio test eliminates cost 0.34x the floor to eliminate; the one state it
    keeps buys 0.0029x the floor to keep. On this window the whole chain is algebra as far as
    the plate can tell, and that is measured on the observable rather than inferred from a
    time-constant ratio.

    Returns:
        Relative difference in the late-window activity fold, per dose.
    """
    k_f = float(FOLDING_RATE_PER_H.bounds[0]) if folding_rate_per_h is None else \
        FOLDING_RATE_PER_H.at(folding_rate_per_h)
    chain = UprChain(basal_load=basal_load, hill_n=hill_n, basal_share=basal_share,
                     folding_rate_per_h=k_f, mu_per_h=mu_per_h)
    shape = PlateUpr(window=window)
    read_h = float(window.duration_h)

    def chain_late(dose: float) -> float:
        trajectory = chain.simulate(entry, dose, window=window,
                                    n_points=CHAIN_QUADRATURE_POINTS)
        rate = basal_share + UprChain.occupancy(trajectory)
        late = trajectory.t >= LATE_WINDOW_FRACTION * read_h
        return float(np.mean(rate[late]))

    theta0 = _occupancy(_spliced_fraction(
        _ire1_active(basal_load / (1.0 - basal_load), 1.0, hill_n), mu_per_h))
    control_chain = chain_late(0.0)
    control_shape = shape._late_promoter_rate(basal_load, hill_n, basal_share, mu_per_h,
                                              theta0)
    out: dict[float, float] = {}
    for dose in doses:
        load = min(entry.client_load(dose, basal_load), 1.0 - 1e-9)
        full = chain_late(dose) / max(control_chain, 1e-12)
        reduced = shape._late_promoter_rate(load, hill_n, basal_share, mu_per_h,
                                            theta0) / max(control_shape, 1e-12)
        out[float(dose)] = abs(full - reduced) / max(reduced, 1e-12)
    return out


DEACTIVATION_FRACTION = 0.1
"""What "the response has deactivated" means here: the excess over basal has fallen to a
tenth of its own peak. Pincus reports the fact ("HAC1 splicing deactivates within 2 h at
1.5 mM DTT") without a threshold, so one has to be chosen and said out loud. ASSERTED as a
reading convention rather than a measurement, and :func:`deactivation_time` is a property of
the six-state chain, not of the plate -- the plate's own deactivation is
:func:`deactivation_rate` and it is faster than the chain can go."""


def deactivation_time(entry: DoseEntry, dose_mM: float, *, mu_per_h: float,
                      hill_n: float, basal_load: float, basal_share: float = 0.0,
                      folding_rate_per_h: float | None = None,
                      fraction: float = DEACTIVATION_FRACTION,
                      window: Window = PLATE_READ_4H) -> float:
    """Hours from the DTT step until the splicing response has fallen back, per the chain.

    The UNFITTED prediction of :data:`HAC1_DEACTIVATION_TIME`. Nothing in the plate fit sees
    a deactivation time -- the four scalars are fitted to late-window folds -- so a number
    that lands inside Pincus's measured 2-4 h bracket is a prediction and not a recovery.

    Raises:
        ValueError: if the response has not fallen back inside the window, which is itself
            informative and is not rounded to the window's end.
    """
    k_f = float(FOLDING_RATE_PER_H.bounds[0]) if folding_rate_per_h is None else \
        FOLDING_RATE_PER_H.at(folding_rate_per_h)
    chain = UprChain(basal_load=basal_load, hill_n=hill_n, basal_share=basal_share,
                     folding_rate_per_h=k_f, mu_per_h=mu_per_h)
    trajectory = chain.simulate(entry, dose_mM, window=window, n_points=401)
    excess = UprChain.occupancy(trajectory) - UprChain.occupancy(trajectory)[0]
    peak = int(np.argmax(excess))
    if excess[peak] <= 0:
        raise ValueError(
            f"the chain shows no induction at {dose_mM:g} mM, so there is nothing to "
            f"deactivate -- check that the declared entry puts the load above basal")
    fallen = np.where(excess[peak:] <= fraction * excess[peak])[0]
    if not len(fallen):
        raise ValueError(
            f"the response is still at more than {fraction:.0%} of its peak at the end of "
            f"{window.name} ({window.duration_h:g} h). That is Pincus's own 5 mM behaviour "
            f"-- maximal and sustained -- and it is reported rather than clipped")
    return float(trajectory.t[peak + int(fallen[0])])


def plate_gate() -> FreeScalarGate:
    """Criterion (e) for the shipped shape: fitted scalars against independent targets."""
    return free_scalar_gate(PlateUpr(), UPR_TARGETS, "mech/upr.py -- plate-window shape")


def provenance() -> str:
    """Both registries as markdown, with their gates on top."""
    return provenance_table(UPR_PARAMS, PLATE_PARAMS)
