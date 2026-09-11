"""Plasmid-bearing and plasmid-free cells, and the ethanol the fed-batch layer never counted.

TWO STATES, TWO PARAMETERS, AND THEY ARE THE ONLY MECHANISM LEFT IN THE VESSEL THAT MAKES THE
TITRE. `ARCHITECTURE_GAPS.md` 0.1 tabulates every state in the target architecture against the
``tau/T < 0.146`` rule and finds that on a 5-day fed-batch only X, F and g survive -- every
signalling, transcript, effector, UPR and maturation state is algebra there, including the
marquee cross-protection result. So this module is not a stress block competing with the other
nine: it is the whole of the mechanistic layer in the one vessel that produces product, and
`REVISED_BUILD_LIST.md` moves it from Phase 4 to Phase 2 for exactly that reason.

THE SYSTEM. Imanaka & Aiba's two-population mechanism
(Ann NY Acad Sci 1981;369:1-14, PMID 7020540) as instantiated for *S. cerevisiae* by Cheng
1992 (Ohio State, thesis Ch. III; published as Cheng, Huang & Yang, Biotechnol Bioeng
1997;56:23-31, PMID 18636606)::

    X+ --mu+--> (2 - p) X+  +  p X-        mu+ = (1 - b) mu_free
    X- --mu---> 2 X-                        mu- = mu_free

Reduced to the bearing fraction ``F = X+/(X+ + X-)``, which is the form that answers the
fed-batch question and appears in neither source::

    dF/dt = F [ (1 - F) (mu+ - mu-) - p mu+ ]  =  -F mu_free [ (1 - F) b + p (1 - b) ]
    dg/dt = mu_pop / ln 2,   mu_pop = mu_free (1 - F b)

**mu_free CANCELS OUT of dF/dg.** Divide the two and the growth rate disappears: the collapse
is exactly generation-invariant, so a feeding strategy buys wall-clock and never buys
generations. That is a property of the composed system rather than an input to it, and it is
the single most useful thing this layer tells a fed-batch designer.

WHAT WAS DROPPED, AND WHY. Cheng's system also carries S, P, Y_XS, Y_PS, k_s, t_lag and a
plasmid-free death rate k_d. Every one of them is gone:

* **k_d, and with it the come-back phenomenon.** Criterion (e) is the reason and it is
  arithmetic: keeping it makes three free scalars against two independent targets and
  :meth:`~ystwin.mech.params.ParamRegistry.require_gate` refuses the piece. It is also
  unmeasurable here -- Cheng FITTED 0.05-0.1 /h and no independent *S. cerevisiae*
  measurement of it exists in any source this project could reach. Dropping it costs the
  fall-then-rise of F in stationary phase that Cheng and Mead, Gardner & Oliver 1986
  independently reported, and it costs the pulse-starvation recommendation that rests on it.
  :data:`NOT_BUILT` records the refusal so that nobody re-derives it.
* **The expressing loss band as a SECOND scalar.** The AMB yEGFP3 paper's two-to-fivefold
  multiplier on ``p`` is real and measured, and registering it beside the non-expressing
  band made this piece three free scalars against two targets, which ``require_gate()``
  refused at import --
  correctly, on a miscount. Which band ``p`` is swept over is a fact about which plasmid is in
  the flask, not a degree of freedom the model gets to spend, so it lives in
  :data:`BAND_VARIANTS`, outside the gate. The fix was the accounting, not a weaker gate.
* **Cheng's product ODE.** ``Y_PS = 0.135 g beta-galactosidase / g glucose`` is an ONPG
  activity reading wearing a mass unit: integrated, it puts 36-51% of yeast dry weight into
  one heterologous enzyme, against a total cell protein of ~45% of DCW. Use the repository's
  own productivity and multiply by F.
* **Cheng's burden.** ``mu+_m`` is read off his Figure 3.6 in SELECTIVE medium and ``mu-_m``
  off Figure 3.7 in NON-selective medium, and the two media differ by 0.005% tryptophan in a
  trp1 auxotroph. That is a between-medium contrast reported as a within-strain plasmid
  burden -- structurally the growth-rate confound this repository already retracted once.
* **S, Y_XS, k_s, t_lag.** `generator/culture.py` and `fba/fedbatch.py` already integrate
  biomass and substrate, and t_lag is a curve-fitting device Cheng gives six different values
  for. Four states and no identifiability.
* **Any copy-number state.** HARD REFUSAL: no kinetic model of 2-micron copy-number control
  exists in any organism. Futcher 1986 is topological, Som 1988 is a verbal autoregulation
  model whose only number is a >=100-fold FLP repression, and the Jayaram review (PMID
  25541598) states the predicted replication intermediates have never been observed. This
  layer supplies F, not gene dosage: `pathway/capacity.py`'s ``CapacityUnmeasured`` and
  `pathway/published_cassettes.py`'s ``EPISOMAL_2U`` non-number both stand unchanged.
* **The bacterial partition kernel.** ``P(plasmid-free daughter) = 2^-2n`` at Karim's
  MEASURED 14-34 copies per haploid cell gives 3.73e-9 down to 3.39e-21 per generation
  against a measured 1.0e-2 to 5.7e-2 -- six to nineteen orders of magnitude. Futcher & Cox
  said the same from their own data in 1984. Do not port it.

THE UNIT TRAP, WHICH IS WORTH 44% AND HAS ALREADY BITTEN THIS REPOSITORY. The continuous ``p``
of the ODE and the per-generation loss rate ``theta`` that every wet-lab paper reports are not
the same number: ``theta = 1 - 2^-p``, so ``p = -log2(1 - theta) ~ 1.4427 theta``. See
:func:`loss_per_generation` and :func:`divisions_per_generation`, which are the only sanctioned
conversion in this module. `mech/state.py`'s ``plasmid_bearing`` row currently carries
``tau_growth_multiple = 1/(-ln(1-0.05)) = 19.4957``; the correct multiple for a per-GENERATION
loss rate is ``1/p = 1/(-log2(1-0.05)) = 13.5134``, and the two differ by exactly ``1/ln 2``.
That file is not this module's to edit; :func:`population_states` computes its own.

CRITERION (a), MEASURED ON THIS PROJECT'S OWN PLATES by :func:`criterion_a_table`, and it is a
duration switch rather than a verdict. Both channels are scored through
:mod:`ystwin.mech.ablation` against the two floors this platform has measured, with both models
refit on the same real rows and the swept axes DECLARED rather than fitted::

    reporter endpoint, 12 committed blocks    4.14 h: 0.0002-0.141x the floor  clears 0/12
    (floor 0.146 relative CV, mu = 0.3227)      24 h: 0.0177-0.635x at p_low   clears 0/12
                                                24 h: 1.58-2.85x at p_high     CLEARS 12/12
    chord growth rate, 183 committed wells    4.14 h: 0, to 2.1e-12x         clears 0/2
    (floor 0.0117 /h, absolute)                 24 h: 0.204-0.692x at b = 0.10 clears 0/2
                                                24 h: 1.23-2.57x at b = 0.234  CLEARS 2/2

Both the harvest ("8% and close to dead weight on a single overnight plate") and its verifier
("1.08x to 15x at 8 generations") were partly wrong, and so was this module's own first draft,
which quoted 1.69-2.85x for the whole reporter band when only its top end clears anything. On
the 4.14 h read the layer is absorbed into the reporter's own decay constant and, in the growth
channel, into the fitted ``mu_free`` EXACTLY -- a read whose length equals the length the model
was fitted on cannot see this block at all, which is a structural statement and not a small
number. **MEASURED, so that the caveat is a number and not a sentence**: across all twelve
committed blocks at both ends of the loss band the two models' residual RMS on the real rows
differ by at most 3.06%, i.e. 0.21x the 0.146 floor, and on thirteen of the twenty-four fits
(four of them inside 4e-6 of a tie) the model WITH the bearing fraction is the marginally worse
one. The 24 h verdict is therefore an extrapolation of two fits that every committed row rates
equally, evaluated against a floor that was itself measured on 4.14 h reads. It is also an
extrapolation in the GROWTH RATE it drives the segregation with: scored at the 0.3227 /h chord
of the 4.14 h read it is 1.58-2.85x the floor, and re-scored at the 0.1869 /h chord
`generator/culture.py` itself predicts for a 24 h read it is 1.05-1.83x -- still 12/12, with the
weakest block clearing by 5%. That does not make it wrong; it makes the two
experiments :data:`LOSS_PER_DIVISION` and :data:`BURDEN` each name in their ``missing`` field
-- one plate of replica-plating agar, one OD600 curve against an isogenic empty vector -- the
things that would settle it.

Extended to a 24 h read it clears both floors, at the top of the loss band and at the
top of the burden band respectively; :func:`burden_to_clear_the_growth_floor` puts the growth
crossing at b = 0.225 for the slowest-losing plasmid and b = 0.131 for the fastest. The switch
is duration, and it is measured here rather than argued.

AND THE TWO CHANNELS ARE NOT EQUALLY GOOD, which is the last thing to say about criterion (a).
The reporter channel scores what survives refitting a decay constant that BOTH models carry.
The growth channel does not: `generator/culture.py`'s committed logistic already drops the
chord rate by 0.1286 /h between 4.14 h and 24 h at the plate's own inoculum -- 11.0x the floor
and 4.3x this block's whole effect, same direction, same observable -- so a 24 h OD trace
cannot attribute a falling chord rate to segregation. See :class:`ConstantGrowth`, which
carries the arithmetic. The growth row is a sensitivity; the reporter row carries the claim.

THE FED-BATCH ANSWER, which is what this block was promoted for. Takeover is UNCONDITIONAL in
any non-selective medium -- ``s = -b/(1-b) <= 0`` and ``p > 0`` make both terms of ``dF/dt``
negative for every ``F`` in ``(0, 1]`` -- so a design controls how far, never whether, and the
axis is GENERATIONS. :func:`fedbatch_collapse` run on the vessel `scripts/design_fedbatch_run.py`
publishes says the tension plainly: the 0.5 g/L inoculum spends 6.63 generations and ends with
9.9-43.9% of its biomass-hours carrying no plasmid, and the 0.05 g/L inoculum the fed-batch
layer prefers BECAUSE it buys 9.94 generations ends with 12.8-65.1%. The generation budget is
simultaneously the thing that makes the biomass and the thing that destroys the plasmid, and
the current generator cannot express either half of that.

THE ETHANOL DIAGNOSTIC, which lives in this layer and costs nothing. The string "ethanol" does
not occur in `fba/fedbatch.py` or `fba/dynamic.py` and there is no ``ethanol`` key in the
25-stressor panel -- yet the design tool recommends a run whose harvest this repository's own
six calibration states put at **32-72 g/L of an autotoxin**. :class:`EthanolDiagnostic` reports
``E = r (B - B0) / V`` with ``r`` the measured spread of those six states, and reports nothing
else: no growth penalty, no uptake, no stripping. Those four constants plus the functional form
are :data:`NOT_BUILT` refusals, and `REVISED_BUILD_LIST.md` §3 A2 records why -- the two
glucose-medium anchors bound the inhibition slope from ABOVE, and four independent audits each
put the ablation below the floor. The diagnostic is ``LayerState.REPORTED``: a bug report
against the chemostat-to-fed-batch swap, not a stress pathway.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .. import paths
from ..plate import replay
from .ablation import (
    GROWTH_RATE_FLOOR,
    REPORTER_ACTIVITY_FLOOR,
    AblationResult,
    Fit,
    FittableModel,
    Observable,
    SaturatingReporter,
    ablate,
    reporter_block,
)
from .integrate import Trajectory, integrate_window
from .params import Param, ParamRegistry, RefusedValue
from .state import FEDBATCH_GROWTH_PER_H, Encoding, MechState, StateVar, Window

__all__ = [
    "BAND_VARIANTS",
    "BURDEN",
    "BURDEN_TARGET",
    "ConstantGrowth",
    "DECLARED_READ_H",
    "ETHANOL_GATE",
    "ETHANOL_G_PER_MOL",
    "ETHANOL_PARAMS",
    "ETHANOL_YIELD",
    "EthanolDiagnostic",
    "FOCI_PER_HAPLOID_NUCLEUS",
    "INOCULUM_BEARING_FRACTION",
    "FedBatchCollapse",
    "GLUCOSE_G_PER_MOL",
    "LOSS_PER_DIVISION",
    "LOSS_PER_DIVISION_EXPRESSING",
    "MEASURED_LOSS_PER_GENERATION",
    "NOT_BUILT",
    "PLATE_GLUCOSE_G_PER_L",
    "PLATE_READ_H",
    "POPULATION_GATE",
    "POPULATION_PARAMS",
    "REPORTED",
    "REPORTER_TARGET",
    "SegregatingCulture",
    "SegregatingReporter",
    "Takeover",
    "TwoPopulation",
    "bearing_endpoint_observable",
    "bearing_fraction_at",
    "binomial_partition_loss_per_generation",
    "burden_to_clear_the_growth_floor",
    "cheng_hohnholz_convergence",
    "chord_growth_observable",
    "committed_reporter_blocks",
    "criterion_a_table",
    "divisions_per_generation",
    "ethanol_yield_band",
    "fedbatch_collapse",
    "foci_partition_loss_per_generation",
    "foci_prediction_overlap",
    "loss_per_generation",
    "mean_bearing_fraction",
    "plate_chord_growth_rates",
    "population_states",
    "provenance_summary",
    "score_growth_ablation",
    "score_reporter_ablation",
    "simulate_population",
    "standing_ethanol_band",
]

_LN2 = math.log(2.0)

REPORTED = "reported"
"""``predict.LayerState.REPORTED`` as a bare string, to keep this module's import graph flat.

`ystwin/predict.py` pulls in the GEM, the pathway solver and the stress panel; a mechanistic
state file that imported it would make every one of those a prerequisite of an ODE. The value
is pinned against the enum by ``tests/test_mech_population.py``.
"""


# --------------------------------------------------------------------------------------
# The constants, graded
# --------------------------------------------------------------------------------------

POPULATION_PARAMS = ParamRegistry("mech/population.py::two_population")
"""The built system's constants. Two free scalars, two independent targets -- see
:meth:`~ystwin.mech.params.ParamRegistry.require_gate`, which this module calls at import."""

ETHANOL_PARAMS = ParamRegistry("mech/population.py::fedbatch_ethanol")
"""The diagnostic's constants. One free scalar and NO target on this platform, so it does not
pass a predictive gate -- which is precisely why it ships :data:`REPORTED` and feeds nothing."""

BAND_VARIANTS = ParamRegistry("mech/population.py::band-variants")
"""The SAME free scalar under a different construct, and therefore outside the gate.

``p_loss_per_division_expressing`` is not a second degree of freedom: a run models one
construct, and which band ``p`` is swept over is a fact about which plasmid is in the flask,
not a number the model gets to choose. Registering it beside :data:`LOSS_PER_DIVISION` in the
gated registry counted one scalar twice and made the piece 3-against-2, which
``require_gate()`` refuses -- correctly, on a miscount. The fix is the accounting, not a
weaker gate: ``tests/test_mech_population.py`` pins this registry outside
:data:`POPULATION_PARAMS` and pins the free count at two.
"""

NOT_BUILT = ParamRegistry("mech/population.py::not-built")
"""Refusals: everything the sources name and this module does not implement.

Deliberately outside both gates. A refusal that no equation consumes is a record, not a degree
of freedom -- `mech/ablation.py`'s ``ABLATION_PARAMS`` is kept out of its gate for the mirror
reason. Registering them here is what stops the next pass re-deriving each one from scratch.
"""

MEASURED_LOSS_PER_GENERATION = (0.010, 0.011, 0.015, 0.015, 0.017, 0.017, 0.020, 0.021,
                                0.023, 0.027, 0.028, 0.033, 0.034, 0.036, 0.046, 0.057)
"""Sixteen mitotic loss rates per generation, MEASURED. Hohnholz, Pohlmann & Achstetter 2017,
Yeast 34:267-275, PMID 28207166, Table 1: eight isomeric 2-micron/STB/HIS3 YEp vectors, all
4346 bp and differing only in fragment order, across BY4742 and SY992, non-selective YPDAU, on
average 505 cfu replica-plated per point and never fewer than 227. The 5.7x spread is across
plasmids that are byte-for-byte the same DNA."""

LOSS_PER_DIVISION = POPULATION_PARAMS.add(Param.swept(
    "p_loss_per_division", "per cell division",
    "SWEPT over the Hohnholz 2017 (PMID 28207166) Table 1 band converted division-wise: "
    "theta in [0.010, 0.057] per generation -> p = -log2(1-theta) in [0.014500, 0.084670]. "
    "Sixteen replica-plated estimates on 2-micron/STB/HIS3 YEps in BY4742 and SY992, "
    "non-expressing. No functional form for p exists: Parker & DiBiasio 1987 (PMID 18576378) "
    "measure only that loss FALLS with growth rate and is zero at D >= 0.3 /h, and that it "
    "RISES with expression level; both are signs with no curve behind them",
    bounds=(0.0144995696951151, 0.0846703239869906),
    missing="replica plating on the team's own construct: ~200 cells on YPD, 48 h, replica "
            "onto selective agar at t0 and after N generations, theta = 1 - (F/I)^(1/N). One "
            "plate of agar, no new equipment. Until then this is a sensitivity axis"))

LOSS_PER_DIVISION_EXPRESSING = BAND_VARIANTS.add(Param.swept(
    "p_loss_per_division_expressing", "per cell division",
    "SWEPT. Hohnholz, Pohlmann & Achstetter, Appl Microbiol Biotechnol 2017;101:8455-8463, "
    "PMID 29052760, MEASURED that inserting a TEF1p-yEGFP3 expression block into the same "
    "backbone doubles the loss rate in the best isoform and raises it four- to fivefold in "
    "the other five. Applied to the band ends: theta in [0.020, 0.285] -> p in [0.029146, "
    "0.483985]. Cheng's 25 C expressed p = 0.08 against 32 C repressed p = 0.03 is an "
    "independent 2.7x in the same direction",
    bounds=(0.0291463456595165, 0.4839848529963354),
    missing="the same replica plating, on the EXPRESSING construct rather than on an empty "
            "vector. The multiplier is measured; which end of it this construct sits at is not"))

BURDEN = POPULATION_PARAMS.add(Param.swept(
    "b_burden", "dimensionless fraction of mu_free",
    "SWEPT, and MARKER-DOMINATED so it cannot be inherited from any paper. Karim, Curran & "
    "Alper 2013, FEMS Yeast Res, PMID 23107142, MEASURED 2-micron LEU2 plasmids growing at "
    "76.6% of a control 0.303 /h (b = 0.234) and 2-micron KanMX at essentially zero; "
    "Hohnholz 2017 (PMID 28207166) MEASURED equal doubling times for a non-expressing HIS3 "
    "YEp in both media (b = 0), and PMID 29052760 >= 10% retardation once "
    "the expression block is added",
    bounds=(0.0, 0.234),
    missing="one OD600 growth curve of the bearing strain against an isogenic EMPTY-VECTOR "
            "control in the SAME medium, on the plate reader this project already owns. Do "
            "NOT inherit Cheng's mu+ - mu-: his two rates come from two different media"))

ETHANOL_G_PER_MOL = ETHANOL_PARAMS.add(Param.derived(
    "ethanol_molar_mass", 46.069, "g/mol",
    "DERIVED from IUPAC/CIAAW conventional atomic weights: C2H6O = 2(12.011) + 6(1.008) + "
    "15.999 = 46.069. pathway/environment_flux.py carries 46.068 for the same species; the "
    "2e-5 relative difference is rounding of the atomic weights and moves nothing here"))

GLUCOSE_G_PER_MOL = ETHANOL_PARAMS.add(Param.derived(
    "glucose_molar_mass", 180.156, "g/mol",
    "DERIVED from IUPAC/CIAAW conventional atomic weights: C6H12O6 = 6(12.011) + 12(1.008) + "
    "6(15.999) = 180.156, which is exactly the value pathway/environment_flux.py uses"))

PLATE_GLUCOSE_G_PER_L = 20.0
"""The sugar charged into a well: `generator/context.py::CultureContext.glucose_g_per_L`'s own
default. Carried as a plain float because it is another module's configuration, not a
measurement of this one -- ``tests/test_mech_population.py`` reads it back from `context.py`
and fails if the default moves."""

INOCULUM_BEARING_FRACTION = 0.95
"""F at inoculation, ASSERTED. The default of every entry point here, so it is named once.

WHY IT IS NOT REGISTERED, which is the honest half. ``Tag.ASSERTED`` is in
:data:`~ystwin.mech.params.FREE_TAGS`, so adding it to :data:`POPULATION_PARAMS` makes three
free scalars against two targets and ``require_gate()`` refuses the module at import. It is
carried as a bare named float for the same reason :data:`PLATE_GLUCOSE_G_PER_L` and
:data:`DECLARED_READ_H` are: it is an INITIAL CONDITION of a declared sweep-axis state, not a
constant of the ODE. An integrator who wants it swept must decide that question, not silence it.

IT IS MEASURABLE AND WAS MEASURED, by the same table the loss band comes from. Hohnholz 2017
(PMID 28207166) Table 1 reports the t0 percentage of plasmid-carrying cells for all sixteen
transformant sets: 85.3-96.0%, mean 89.6%. 0.95 sits at the optimistic end of that range and
:func:`fedbatch_collapse` inherits the difference directly -- MEASURED here, the plasmid-free
share of biomass-hours over the published 6.63-generation run spans 5.1-37.9% at F0 = 1.0,
9.9-43.9% at 0.95 and 14.8-49.4% at the measured mean 0.898. At the least severe corner
(p_low, b = 0) roughly half the headline 9.9% is this number and not the run. The criterion-(a)
reporter channel is EXACTLY blind to it (F0 rescales R0, k_basal and k_max inside their own
box, verified at 2.75x the floor for F0 in [0.80, 1.00]) and the growth channel moves the
verdict by 9% and does not flip it (2.38-2.62x the floor across [0.853, 1.00])."""

FOCI_PER_HAPLOID_NUCLEUS = (3, 5)
"""Fluorescent STB-reporter foci per haploid nucleus, MEASURED by microscopy (Velmurugan et al.,
via Liu, Sau, Ma & Jayaram 2014, Microbiol Spectr 2(5), PMID 25541598). The review's position is
load-bearing here: "plasmid segregation as a single clustered unit is a misimpression... each
plasmid focus appears to function as an INDEPENDENT entity in segregation." See
:func:`foci_partition_loss_per_generation`."""


# Refusals. Each one is a constant a source names, that this module does not implement.
NOT_BUILT.add(Param.refused(
    "k_d_plasmid_free_death", "1/h",
    "Cheng 1992 thesis Table 3.2 FITTED 0.1 / 0.1 / 0.05 /h at 25 / 28 / 32 C in selective "
    "medium and 0.0 in non-selective; the same authors' 1997 follow-up (PMID 18636606) then "
    "attributes the come-back to a death-rate difference caused by glucose starvation in a "
    "NON-selective feed, i.e. abandons the 0.0",
    reason="only fitted values exist and they contradict each other across the same authors' "
           "own two papers; no independent S. cerevisiae measurement was found anywhere. "
           "Keeping it would also make three free scalars against two targets, which "
           "criterion (e) refuses",
    missing="a measured differential death rate of plasmid-free cells under carbon "
            "starvation. Dropping it costs the fall-then-rise of F in stationary phase "
            "(Cheng's 'come-back phenomenon'; Mead, Gardner & Oliver 1986 report the same "
            "sign independently) and every pulse-starvation feeding recommendation"))

NOT_BUILT.add(Param.refused(
    "i_ethanol_growth_inhibition", "1/h per g/L",
    "REVISED_BUILD_LIST.md §3 A2: the advocate's own two glucose-medium anchors bound this "
    "slope from ABOVE rather than below, excluding the value that produced its headline",
    reason="the law mu = mu_c - i*max(0, E - E_theta) needs five unmeasured constants, and "
           "four independent audits each put its ablation below the 0.146 floor",
    missing="a measured mu(E) curve for THIS strain in THIS medium. Experiment 1.4 of the "
            "revised list -- three offline enzymatic ethanol reads on a run already planned "
            "-- is what would start it"))

NOT_BUILT.add(Param.refused(
    "E_theta_inhibition_threshold", "g/L",
    "REVISED_BUILD_LIST.md §3 A2 audits it at 30 g/L; no threshold measured on this "
    "project's strain or medium exists",
    reason="a threshold is a property of a strain and a medium, and quoting one from another "
           "strain is how this repository has been burned before",
    missing="the same mu(E) curve. Note the consequence: whether the plate layer is "
            "PROVABLY inert turns on this number, so what the code can honestly report is "
            "the stoichiometric ceiling and not the verdict"))

NOT_BUILT.add(Param.refused(
    "q_ethanol_uptake", "mmol/gDCW/h",
    "Elizondo & Saa 2025's six chemostat states (PMID 40891387) measure a NET ethanol flux, "
    "and at their own steady state E = q M X / D that flux stands at only 0.28-1.26 g/L -- "
    "see mech.population.standing_ethanol_band, which recomputes it from the same table. "
    "The diagnostic's r is derived at that scale and applied out to tens of g/L, where the "
    "net flux is a strong function of E",
    reason="no co-consumption rate is measured for these strains at these concentrations",
    missing="paired ethanol and biomass reads over a fed-batch. Until then E as computed "
            "here is an UPPER BOUND: every sink is set structurally to zero"))

NOT_BUILT.add(Param.refused(
    "k_strip_ethanol_evaporation", "1/h",
    "No stripping coefficient is measured for the vessel scripts/design_fedbatch_run.py "
    "specifies (2 L, aerated), and it depends on the gassing rate that run has not fixed",
    reason="a stripping rate is a property of a vessel and a sparge rate, neither of which "
           "this design has committed to",
    missing="an off-gas ethanol balance, or an offline read at t = 0, mid-run and harvest"))

NOT_BUILT.add(Param.refused(
    "dn_dt_copy_number", "copies/cell/h",
    "Searched: Futcher 1986 (J Theor Biol 119:197-204) is a TOPOLOGICAL model with no rate "
    "constants; Som et al. 1988 (Cell 52:27-37) is a verbal autoregulation model whose only "
    "number is a >=100-fold FLP repression; Liu, Sau, Ma & Jayaram 2014 (PMID 25541598) "
    "states the Futcher mechanism's predicted replication intermediates have never been "
    "observed. BioModels returns 0 hits for 'plasmid segregation'",
    reason="there is NO published kinetic model of 2-micron copy-number control in any "
           "organism, so any dn/dt written here would be invention",
    missing="a measured copy-number relaxation rate. This layer supplies F, not dosage: "
            "pathway/capacity.py's CapacityUnmeasured refusal is untouched by it"))

NOT_BUILT.add(Param.refused(
    "k_deg_selection_marker", "1/h",
    "Srienc, Campbell & Bailey 1986 (PMID 18555421) MEASURED that plasmid-free cells keep "
    "growing in selective medium on inherited marker protein -- 'inconsistent with the "
    "hypothesis that plasmid-free cells are unable to grow in selective medium' -- but the "
    "decay rate that sets the memory length is behind a paywall and was not read",
    reason="the mechanism is published and validated and the constant is not",
    missing="the marker protein's degradation rate. Only needed if this twin is ever asked "
            "to predict SELECTIVE-medium plateaus; a non-selective fed-batch does not use it"))

NOT_BUILT.add(Param.refused(
    "p_of_growth_rate", "dimensionless",
    "Parker & DiBiasio 1987 (PMID 18576378) MEASURED that a 2-micron plasmid was 'completely "
    "stable' at D = 0.3 and 0.37 /h and less stable below. Two endpoints, no curve",
    reason="two endpoints do not define a function, and the repository has already retracted "
           "one channel for being a growth-rate confound",
    missing="a p(mu) curve. THIS IS THE BIGGEST HOLE FOR A FED-BATCH: the vessel runs at "
            "0.101 /h by design, a third of the rate at which the only measurement says loss "
            "vanishes, so a constant-p model is optimistic in the regime that matters"))

NOT_BUILT.add(Param.refused(
    "p_of_temperature", "dimensionless",
    "Cheng 1992 Table 3.1 gives p = 0.08 / 0.04 / 0.03 at 25 / 28 / 32 C",
    reason="that ladder is a property of his MATalpha2/HMLalpha2 temperature-REGULATED "
           "promoter -- expression is on at 25 C and repressed at 32 C -- so it measures the "
           "expression effect, not a physical temperature law for plasmids",
    missing="loss rates at several temperatures under constitutive expression"))

NOT_BUILT.add(Param.refused(
    "structural_instability_rate", "per generation",
    "Hohnholz & Achstetter 2019 recovered 586 ampR clones after cultivation and found ~11% "
    "altered forms, calling it 'likely a lower limit'",
    reason="an endpoint with no generation count attached is not a rate",
    missing="altered-form fraction against generations. It is a SECOND loss mode the "
            "two-population model does not represent at all"))

BURDEN_TARGET = POPULATION_PARAMS.add_target(GROWTH_RATE_FLOOR)
"""Target 1: the specific growth rate, floor 0.0117 /h ABSOLUTE, on the plate reader.

It pins ``b`` and is EXACTLY BLIND to ``p``: with b = 0 the two populations grow at the same
rate, total OD is identical whatever F does, and the measured ablation ratio on eight
committed wells is 0.000. That blindness is what makes it independent of :data:`REPORTER_TARGET`
rather than a second look at the same number."""

REPORTER_TARGET = POPULATION_PARAMS.add_target(REPORTER_ACTIVITY_FLOOR)
"""Target 2: dilution-corrected mCitrine activity, floor 0.146 relative CV.

Bulk RFU/OD over a mixed culture is ``F(t) x R_per_bearing_cell(t)``, so this channel carries
``p*mu`` once ``b`` is pinned by :data:`BURDEN_TARGET`. It is also the channel the layer
exists to protect: without F, `reporter.py` has nowhere to put a decay that is segregation and
``promoter_activity()`` reads it as promoter downregulation -- the same shape of error as the
growth-rate confound that already retracted this repository's promoter channel."""

POPULATION_GATE = POPULATION_PARAMS.require_gate()
"""Criterion (e), computed at import and refusing rather than reporting: 2 free scalars
(``p``, ``b``) against 2 independent targets. Adding ``k_d`` makes it 3 against 2 and this
line raises, which is why ``k_d`` is in :data:`NOT_BUILT`."""

# --------------------------------------------------------------------------------------
# The unit bridge. The only sanctioned conversion in this module.
# --------------------------------------------------------------------------------------

def loss_per_generation(p_per_division: float) -> float:
    """``theta = 1 - 2^-p``: the ODE's continuous ``p`` as the number a wet lab reports.

    Exact when the two populations grow at the same rate, which is the regime Hohnholz
    explicitly measured and the regime the estimator ``theta = 1 - (F/I)^(1/N)`` assumes.
    Conflating the two symbols is a 44% error and it hides a real convergence: Cheng's fitted
    p = 0.08 becomes 5.39%/generation, inside Hohnholz's independently replica-plated
    [1.0, 5.7]%, while the unconverted 8.0% falls outside it.
    """
    if p_per_division < 0:
        raise ValueError(f"p must be non-negative, got {p_per_division}")
    return -math.expm1(-p_per_division * _LN2)


def divisions_per_generation(theta_per_generation: float) -> float:
    """``p = -log2(1 - theta)``: the inverse of :func:`loss_per_generation`."""
    if not 0.0 <= theta_per_generation < 1.0:
        raise ValueError(
            f"theta is a probability per generation and must be in [0, 1), got "
            f"{theta_per_generation}")
    return -math.log2(1.0 - theta_per_generation)


def foci_partition_loss_per_generation(foci: int) -> float:
    """``2^-2m`` for ``m`` independently partitioning foci. The reconciliation, not a fit.

    IN HOHNHOLZ'S UNITS, WHICH IS THE WHOLE POINT OF RETURNING THIS NUMBER. ``m`` foci
    replicate to ``2m`` copies before division, each daughter draws them independently, and
    ``2^-2m`` is the fraction of DAUGHTERS that get none -- exactly what
    ``theta = 1 - (F/I)^(1/N)`` estimates (Hohnholz 2017 Table 1, footnote e). The
    probability that a division produces a plasmid-free cell is twice this, ``2^(1-2m)``,
    because either daughter can be the empty one: the same per-division-versus-per-daughter
    factor of 2 that :func:`cheng_hohnholz_convergence` prices in Cheng's discrete kernel,
    and it must not be applied twice.

    The bacterial kernel counts MOLECULES and fails for yeast by six to nineteen orders of
    magnitude (:func:`binomial_partition_loss_per_generation`). Counting FOCI instead, at the
    3-5 per haploid nucleus measured by fluorescence microscopy, predicts 1.56e-2 down to
    9.77e-4 per generation -- against 1.0e-2 to 5.7e-2 measured by replica plating in a
    different laboratory by a different method. A number from a microscope predicting a number
    from an agar plate, with nothing fitted in between. It overlaps the measured band only at
    ``m = 3``, on ``[1.0e-2, 1.56e-2]``.

    STATUS, stated because it is not a result: ``m`` is measured for STB-reporter plasmids in
    the NATIVE context, not for an engineered HIS3 YEp, so this is a hypothesis carrying one
    unmeasured structural parameter. DAPI/LacO-array imaging of the team's own construct is
    what would convert it.
    """
    if foci < 1:
        raise ValueError(f"a segregating unit count must be at least 1, got {foci}")
    return 2.0 ** (-2 * int(foci))


def binomial_partition_loss_per_generation(copies: int) -> float:
    """The same kernel applied to individual molecules -- the form that is REFUTED for yeast.

    Kept executable rather than described so the refutation is a number: at Karim 2013's
    measured 14-34 copies per haploid cell this returns 3.73e-9 to 3.39e-21 against a measured
    1.0e-2 to 5.7e-2. The premise it rests on -- random segregation, no active partitioning --
    is false for the 2-micron plasmid, which carries Rep1/Rep2/STB.
    """
    return foci_partition_loss_per_generation(copies)


# --------------------------------------------------------------------------------------
# The two-population system
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Takeover:
    """Whether, and how fast, the plasmid-free population wins.

    Args:
        unconditional: True when ``F = 0`` is the only stable state, which is the case for
            every non-selective medium: the selective advantage of carrying the plasmid is
            ``s = -b/(1-b) <= 0`` and the loss term is negative too, so both terms in
            ``dF/dt`` push the same way and the PER-CAPITA collapse rate
            ``-(1/F) dF/dt = mu[(1-F) b + p(1-b)]`` rises monotonically as F falls. The
            absolute ``|dF/dt|`` does not: it peaks at ``F = (b + p(1-b))/(2b)`` -- 0.639 at
            the top of both swept bands, MEASURED -- and falls after, because F itself is
            shrinking.
        selective_advantage: ``s = (mu+ - mu-)/mu+``. Negative or zero without selection.
        interior_fraction: The stable interior equilibrium ``F* = 1 - p/s`` where one exists,
            which needs ``s > p > 0`` and therefore needs selection. ``None`` otherwise.
        detail: The one-line statement, with the arithmetic that produced it.
    """

    unconditional: bool
    selective_advantage: float
    interior_fraction: float | None
    detail: str

    def summary(self) -> str:
        return self.detail


@dataclass(frozen=True)
class TwoPopulation:
    """The system: bearing fraction and a generation clock, on a supplied host growth rate.

    Args:
        p_per_division: Segregational loss probability per division. A point taken from
            :data:`LOSS_PER_DIVISION` or :data:`LOSS_PER_DIVISION_EXPRESSING` via
            ``Param.at``, or from :func:`divisions_per_generation` on a measured theta.
        burden: ``b``, the fractional growth-rate penalty of the bearing population. A point
            from :data:`BURDEN`.

    The host rate ``mu_free`` is an ARGUMENT to every method rather than a field: this layer
    consumes the growth rate the repository already computes and never asserts one, which is
    the same discipline `fba/fedbatch.py` applies to the fed-batch setpoint.
    """

    p_per_division: float
    burden: float

    def __post_init__(self) -> None:
        if not 0.0 <= self.p_per_division < 1.0:
            raise ValueError(
                f"p is a probability per division and must be in [0, 1), got "
                f"{self.p_per_division}")
        if not 0.0 <= self.burden < 1.0:
            raise ValueError(
                f"burden is a fraction of mu_free and must be in [0, 1), got {self.burden}. "
                f"b = 1 would stop the bearing population growing at all, which is a lethal "
                f"construct rather than a burden")

    @property
    def loss_per_generation(self) -> float:
        """``theta``, the number a replica-plating experiment would report."""
        return loss_per_generation(self.p_per_division)

    @property
    def selective_advantage(self) -> float:
        """``s = (mu+ - mu-)/mu+ = -b/(1-b)``. Zero or negative in any non-selective medium."""
        return -self.burden / (1.0 - self.burden)

    @property
    def peak_loss_fraction(self) -> float:
        """Where ``|dF/dt|`` is largest: ``F* = (b + p(1-b))/(2b)``, or 1 when ``b = 0``.

        The distinction the word "accelerates" hides. The PER-CAPITA rate
        ``-(1/F) dF/dt = mu[(1-F) b + p(1-b)]`` rises monotonically all the way to F = 0, but
        the ABSOLUTE rate is that times F and turns over. At the top of both swept bands
        (p = 0.084670, b = 0.234) the turnover is at F = 0.6386, so most of a fed-batch's
        collapse happens on the DECELERATING side of it.
        """
        if self.burden == 0.0:
            return 1.0
        interior = (self.burden + self.p_per_division * (1.0 - self.burden)) / (2.0 * self.burden)
        return min(interior, 1.0)

    def bearing_growth_rate(self, mu_free_per_h: float) -> float:
        """``mu+ = (1 - b) mu_free``."""
        return (1.0 - self.burden) * float(mu_free_per_h)

    def population_growth_rate(self, fraction: float, mu_free_per_h: float) -> float:
        """``mu_pop = mu_free (1 - F b)``: what an OD600 slope actually measures.

        It RISES as F falls, which is the signature the growth channel carries and the reason
        a culture can appear to accelerate with no change in medium.
        """
        return float(mu_free_per_h) * (1.0 - _unit(fraction) * self.burden)

    def d_fraction_d_time(self, fraction: float, mu_free_per_h: float) -> float:
        """``dF/dt = -F mu_free [ (1 - F) b + p (1 - b) ]``, 1/h."""
        f = _unit(fraction)
        return -f * float(mu_free_per_h) * (
            (1.0 - f) * self.burden + self.p_per_division * (1.0 - self.burden))

    def d_fraction_d_generation(self, fraction: float) -> float:
        """``dF/dg``, with ``mu_free`` cancelled out. The generation-invariance, exactly.

        ``dF/dg = -ln2 F [(1-F) b + p (1-b)] / (1 - F b)``. No growth rate, no substrate, no
        feed rate appears -- so two fed-batches at different setpoints spend the SAME
        generations reaching the same F and differ only in wall-clock.
        """
        f = _unit(fraction)
        return -_LN2 * f * ((1.0 - f) * self.burden
                            + self.p_per_division * (1.0 - self.burden)) / (1.0 - f * self.burden)

    def rhs(self, mu_free_per_h: float):
        """``f(t, y) -> dy/dt`` over ``(fraction, generations)`` for :func:`integrate_window`."""
        def _rhs(_t, y):
            f = _unit(y[0])
            return [self.d_fraction_d_time(f, mu_free_per_h),
                    self.population_growth_rate(f, mu_free_per_h) / _LN2]
        return _rhs

    def fraction_after_generations(self, initial_fraction: float, generations: float) -> float:
        """F after ``g`` generations. Closed form at ``b = 0``, integrated otherwise.

        At ``b = 0`` this is exactly ``F0 * 2^(-p g)``, which is the same statement as
        ``theta = 1 - 2^-p`` and is what makes the unit bridge checkable rather than asserted.
        """
        f0 = _unit(initial_fraction)
        g = float(generations)
        if g < 0:
            raise ValueError(f"generations must be non-negative, got {g}")
        if self.burden == 0.0:
            return f0 * 2.0 ** (-self.p_per_division * g)
        if g == 0.0 or f0 == 0.0:
            return f0
        solution = solve_ivp(lambda _g, y: [self.d_fraction_d_generation(y[0])],
                             (0.0, g), [f0], rtol=1e-10, atol=1e-12)
        if not solution.success:
            raise RuntimeError(f"F(g) integration failed: {solution.message}")
        return _unit(float(solution.y[0, -1]))

    def generations_to_fraction(self, initial_fraction: float, target_fraction: float) -> float:
        """Generations for F to fall from one value to another. The units the field reports in.

        Closed form at ``b = 0``: ``g = log2(F0/F) / p``. Bisected otherwise, on a function
        that is strictly decreasing so the root is unique.
        """
        f0, target = _unit(initial_fraction), _unit(target_fraction)
        if not 0.0 < target < f0:
            raise ValueError(
                f"the target fraction must lie strictly between 0 and F0; got F0 = {f0} and "
                f"target = {target}")
        if self.burden == 0.0:
            if self.p_per_division == 0.0:
                raise ValueError(
                    "with p = 0 and b = 0 nothing is lost and F never falls; that is the "
                    "ablated model, not a run")
            return math.log2(f0 / target) / self.p_per_division
        from scipy.optimize import brentq
        high = 1.0
        while self.fraction_after_generations(f0, high) > target:
            high *= 2.0
            if high > 1e6:
                raise RuntimeError(
                    f"F does not reach {target} within 1e6 generations at p = "
                    f"{self.p_per_division}, b = {self.burden}")
        return float(brentq(
            lambda g: self.fraction_after_generations(f0, g) - target, 0.0, high, xtol=1e-10))

    def takeover(self, differential_death_per_h: float = 0.0,
                 mu_free_per_h: float = 1.0) -> Takeover:
        """The conditions under which the plasmid-free population takes over.

        Args:
            differential_death_per_h: ``k_d``. NOT a parameter of this module -- it is
                :data:`NOT_BUILT` and defaults to 0. It is accepted here so the criterion can
                be STATED in full, because "takeover is unconditional" is only true at
                ``k_d <= b mu_free`` and a reader deserves to see the boundary.
            mu_free_per_h: The host rate, needed only when ``k_d`` is non-zero, because
                ``k_d`` does not scale with mu and ``b mu`` does. That asymmetry is the entire
                mechanism behind the come-back phenomenon this module does not implement.

        Returns:
            A :class:`Takeover`.
        """
        s = self.selective_advantage
        mu = float(mu_free_per_h)
        kd = float(differential_death_per_h)
        if kd <= self.burden * mu:
            peak = (f"|dF/dt| itself peaks at F = {self.peak_loss_fraction:.4g} and falls after"
                    if self.burden > 0 else "|dF/dt| itself is largest at F = 1")
            return Takeover(
                unconditional=True, selective_advantage=s, interior_fraction=None,
                detail=(
                    f"UNCONDITIONAL. s = {s:.6g} <= 0 and p = {self.p_per_division:.6g} > 0, "
                    f"so both terms of dF/dt are negative for every F in (0, 1] and F = 0 is "
                    f"the only stable state. The PER-CAPITA collapse rate ACCELERATES as F "
                    f"falls, because the (1 - F) burden term grows; {peak}. A stable interior "
                    f"fraction F* = 1 - p/s "
                    f"needs s > p > 0, i.e. it needs the plasmid to confer an advantage, "
                    f"which a non-selective medium does not"))
        one_minus = self.p_per_division * (1.0 - self.burden) * mu / (kd - self.burden * mu)
        interior = 1.0 - one_minus
        if interior <= 0.0:
            return Takeover(
                unconditional=True, selective_advantage=s, interior_fraction=None,
                detail=(
                    f"UNCONDITIONAL even at k_d = {kd:g} /h: the interior root sits at "
                    f"{interior:.6g}, outside (0, 1), so selection is too weak to hold any "
                    f"fraction against p = {self.p_per_division:.6g}"))
        return Takeover(
            unconditional=False, selective_advantage=s, interior_fraction=interior,
            detail=(
                f"ARRESTED at F* = {interior:.6g}, because k_d = {kd:g} /h exceeds "
                f"b*mu = {self.burden * mu:.6g} /h. NOTE that k_d is a NOT_BUILT refusal: "
                f"only fitted values exist for it and the authors who fitted them "
                f"contradicted themselves across two papers"))


def _unit(value: float) -> float:
    """Clip a fraction into [0, 1]. Applied inside the right-hand side, never to a solution."""
    return min(max(float(value), 0.0), 1.0)


# --------------------------------------------------------------------------------------
# The states, and the tau/T audit that decides whether they are states at all
# --------------------------------------------------------------------------------------

def population_states(window: Window, *, p_bounds: tuple[float, float],
                      burden: float = 0.0) -> MechState:
    """``(plasmid_bearing, generations)`` with time constants computed for this window.

    The time constant of F is not a constant of the plasmid: it is ``1/(p mu)`` at ``F -> 1``
    and ``1/(mu (b + p(1-b)))`` at ``F -> 0``, both exact linearisations of ``dF/dt`` at the
    ends of the state's own range, and both evaluated at the window's MEASURED growth band.
    That is why this is a function of the window rather than a table row.

    IT DELIBERATELY DOES NOT REUSE `mech/state.py`'s ``plasmid_bearing``. That row carries
    ``tau_growth_multiple = 1/(-ln(1-0.05)) = 19.4957``, which is the time constant for a loss
    rate of 0.05 per E-FOLD of biomass; the published quantity is per GENERATION and gives
    ``1/(-log2(1-0.05)) = 13.5134``. The two differ by exactly ``1/ln 2 = 1.4427`` -- the same
    unit trap :func:`loss_per_generation` exists to close. No verdict flips at 0.05, but the
    number is wrong and this module computes its own.

    Args:
        window: The declared operating window. Supplies the growth band.
        p_bounds: ``(low, high)`` loss probability per division; pass
            ``LOSS_PER_DIVISION.bounds`` or the expressing band.
        burden: ``b``, which enters the fast end of the F interval.
    """
    p_low, p_high = (float(p_bounds[0]), float(p_bounds[1]))
    if not 0.0 < p_low <= p_high < 1.0:
        raise ValueError(f"p_bounds must be an increasing band inside (0, 1), got {p_bounds}")
    mu_low = window.growth_rate_low_per_h
    mu_high = window.growth_rate_high_per_h
    tau_fast = 1.0 / (mu_high * (burden + p_high * (1.0 - burden)))
    tau_slow = 1.0 / (mu_low * p_low * (1.0 - burden))
    return MechState((
        StateVar(
            name="plasmid_bearing", units="fraction of cells bearing the plasmid",
            tau_h=(tau_fast, tau_slow),
            source=(
                f"DERIVED at this window's measured growth band {mu_low:g}-{mu_high:g} /h "
                f"from p in [{p_low:g}, {p_high:g}] per division and b = {burden:g}. "
                f"tau = 1/(p(1-b)mu) as F -> 1 and 1/((b + p(1-b))mu) as F -> 0, the exact "
                f"linearisations of dF/dt at the two ends of the state's range. p comes from "
                f"Hohnholz 2017 (PMID 28207166) Table 1 through theta = 1 - 2^-p"),
            constrained_by=(
                "replica plating onto selective agar, theta = 1 - (F/I)^(1/N) -- an estimator "
                "verified against all sixteen rows of Hohnholz 2017 Table 1; and jointly with "
                "the reporter's own decay, dilution-corrected mCitrine RFU/OD on the 96-well "
                "reader"),
            sweep_axis=True,
            note=("The ONE state whose verdict flips between the two windows, and it flips "
                  "differently for an expressing construct: MEASURED by reduce_for_window, the "
                  "4.14 h read FREEZES it at tau/T = 7.86-68.0 for the non-expressing band and "
                  "STRADDLES at 1.37-33.8 for the expressing one, while the 5 d fed-batch KEEPS "
                  "it in both (0.974-5.69 and 0.170-2.83)")),
        StateVar(
            name="generations", units="generations",
            tau_growth_multiple=_LN2,
            encoding=Encoding.INTEGRATING,
            source=("DEFINITION, not a literature number. dg/dt = mu_pop/ln2, so one "
                    "generation is ln2/mu_pop. It is a state because every published effect "
                    "in this field is per-generation and this repository has no other clock"),
            constrained_by="OD600 on the 96-well plate reader"),
    ))


def simulate_population(model: TwoPopulation, window: Window, *,
                        mu_free_per_h: float, initial_fraction: float = INOCULUM_BEARING_FRACTION,
                        initial_generations: float = 0.0,
                        n_points: int = 201) -> Trajectory:
    """Integrate ``(F, g)`` across a declared window with the package's pinned stiff driver.

    Args:
        model: The system.
        window: The declared window; its ``duration_h`` sets the run.
        mu_free_per_h: The UNBURDENED host growth rate. Read it off the biomass series or the
            GEM setpoint -- this module never asserts one.
        initial_fraction: F at inoculation. Defaults to
            :data:`INOCULUM_BEARING_FRACTION`, which is ASSERTED and carries its own
            measured range and sensitivity -- below 1.0 because a transformant preculture has
            already spent generations: Hohnholz's own industrial framing counts 50-60 of them
            from the primary seed lot to the final fermenter.
        initial_generations: g at inoculation, for continuing a seed train.
        n_points: Output grid size only; the solver picks its own steps.
    """
    if mu_free_per_h <= 0:
        raise ValueError(f"mu_free must be positive, got {mu_free_per_h}")
    system = population_states(
        window, p_bounds=(max(model.p_per_division, 1e-9), max(model.p_per_division, 1e-9)),
        burden=model.burden)
    return integrate_window(
        model.rhs(mu_free_per_h),
        {"plasmid_bearing": _unit(initial_fraction), "generations": float(initial_generations)},
        window, system, n_points=n_points)


# --------------------------------------------------------------------------------------
# The fed-batch ethanol diagnostic. REPORTED, and it feeds nothing.
# --------------------------------------------------------------------------------------

def ethanol_yield_band(path=None) -> tuple[float, float]:
    """``r`` in g ethanol per g DCW, DERIVED from this repository's own six calibration states.

    ``r = q_ethanol * M_ethanol / mu`` on each row of
    ``data/carotenoid/elizondo2025_steady_states.tsv`` -- three carotenoid strains at two
    dilution rates, the same six states `pathway/calibrations.py` fits its entry-flux scalar
    on. Computed here rather than transcribed, so the band cannot drift from the table.

    The band is 1.325 to 2.921 g/g. THE CAVEAT TRAVELS WITH IT: these are chemostat steady
    states whose own standing ethanol is only 0.28-1.26 g/L
    (:func:`standing_ethanol_band`, recomputed from the same six rows), and the diagnostic
    applies them out to tens of g/L where the net flux is a strong function of E. That is why
    the result is a band, why it is REPORTED, and why experiment 1.4 of the revised list --
    three offline enzymatic reads on a run already planned -- is what converts it into a
    measurement.
    """
    source = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv" if path is None \
        else path
    with open(source, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no calibration states in {source}")
    yields = [float(r["q_ethanol"]) * float(ETHANOL_G_PER_MOL) / 1000.0 / float(r["mu_per_h"])
              for r in rows]
    return (min(yields), max(yields))


def standing_ethanol_band(path=None) -> tuple[float, float]:
    """``E = q M X / D`` on the same six rows: the concentration ``r`` was actually measured at.

    The scale caveat, computed rather than transcribed, because it is the one number that says
    how far the diagnostic extrapolates. Each Elizondo state is a chemostat at steady state, so
    its standing ethanol is its production flux divided by its dilution rate: 0.279 to 1.263
    g/L across the six. :meth:`EthanolDiagnostic.at_harvest_g_per_l` then applies the same
    ``r`` at 32-72 g/L -- a 25- to 250-fold extrapolation in the variable the net flux depends
    on, which is why ``q_ethanol_uptake`` is a :data:`NOT_BUILT` refusal and E is an upper
    bound.
    """
    source = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv" if path is None \
        else path
    with open(source, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no calibration states in {source}")
    standing = [float(r["q_ethanol"]) * float(ETHANOL_G_PER_MOL) / 1000.0
                / float(r["mu_per_h"]) * float(r["biomass_g_per_l"]) for r in rows]
    return (min(standing), max(standing))


ETHANOL_YIELD = ETHANOL_PARAMS.add(Param.swept(
    "r_ethanol_per_biomass", "g ethanol / g DCW",
    "SWEPT over the MEASURED spread of this repository's own six calibration states "
    "(Elizondo & Saa 2025, PMID 40891387, S4; three strains at D = 0.101 and 0.254 /h), "
    "r = q_ethanol * 46.069 / 1000 / mu -- see mech.population.ethanol_yield_band, which "
    "recomputes it from data/carotenoid/elizondo2025_steady_states.tsv rather than "
    "transcribing it. Six values: 2.920, 2.921, 2.742, 2.870, 1.325, 2.636 g/g",
    bounds=(1.325042521108872, 2.9208101337521857),
    missing="offline enzymatic ethanol at t = 0, mid-run and harvest on the fed-batch "
            "scripts/design_fedbatch_run.py already specifies. Three samples, one kit. It "
            "also measures the net sink, which is the constant this band is weakest on"))

ETHANOL_GATE = ETHANOL_PARAMS.gate()
"""Criterion (e) for the diagnostic: 1 free scalar (``r``) against 0 independent targets on
this platform, so it does NOT pass. That is the whole reason it ships :data:`REPORTED` --
nothing it computes enters a prediction, so the free scalar buys no unearned agreement. The
predictive form, with a growth penalty, is refused in :data:`NOT_BUILT`.

Computed HERE and not beside the population gate: ``r`` is registered below the class that
uses it, and a gate taken before its own parameter was registered reported 0 free scalars and
passed. `tests/test_mech_population.py` pins it as failing."""


@dataclass(frozen=True)
class EthanolDiagnostic:
    """Accumulated ethanol on the fed-batch layer, as a REPORTED band and nothing else.

    ``E = r (B - B0) / V``. Exact rather than a reduction: the only source term is proportional
    to ``dB/dt`` and BOTH sinks are structurally zero -- ``q_ethanol_uptake`` and ``k_strip``
    are :data:`NOT_BUILT` refusals -- so the integral is closed-form and needs no state. The
    consequence is that E as reported here is an **upper bound**: any real sink lowers it.

    Args:
        yield_low_g_per_g, yield_high_g_per_g: The band from :func:`ethanol_yield_band`.

    The class has no method that returns a growth rate, a viability or a titre. That is the
    design: `REVISED_BUILD_LIST.md` §3 splits the ethanol arm into a growth-inhibition law
    (REJECTED, five unmeasured constants) and this accumulation (BUILD, zero), and the split
    is enforced by there being nothing here to call.
    """

    yield_low_g_per_g: float = field(default_factory=lambda: ethanol_yield_band()[0])
    yield_high_g_per_g: float = field(default_factory=lambda: ethanol_yield_band()[1])
    layer_state: str = REPORTED

    def __post_init__(self) -> None:
        if not 0 < self.yield_low_g_per_g <= self.yield_high_g_per_g:
            raise ValueError(
                f"the ethanol yield band must be positive and ordered, got "
                f"({self.yield_low_g_per_g}, {self.yield_high_g_per_g})")
        if self.layer_state != REPORTED:
            raise ValueError(
                f"this diagnostic is {REPORTED!r} by construction and cannot be promoted by "
                f"setting a field: promoting it would need the growth penalty, whose five "
                f"constants are NOT_BUILT refusals. Got {self.layer_state!r}")

    def accumulated_g_per_l(self, biomass_made_g: float, volume_l: float) -> tuple[float, float]:
        """``(low, high)`` g/L of ethanol from a mass of biomass made, in a final volume."""
        if biomass_made_g < 0:
            raise ValueError(f"biomass made must be non-negative, got {biomass_made_g}")
        if volume_l <= 0:
            raise ValueError(f"volume must be positive, got {volume_l}")
        return (self.yield_low_g_per_g * biomass_made_g / volume_l,
                self.yield_high_g_per_g * biomass_made_g / volume_l)

    def at_harvest_g_per_l(self, design) -> tuple[float, float]:
        """``(low, high)`` g/L at the moment the vessel is full.

        Args:
            design: A :class:`ystwin.fba.fedbatch.FedBatchDesign`. Duck-typed rather than
                imported so this module does not pull in cobra through the fba package.

        The number this exists to produce. On the design ``scripts/design_fedbatch_run.py``
        publishes today -- 2 L vessel, 0.5 g/L inoculum, X_final 24.5 g/L -- it is 32.1 to
        70.8 g/L, against a fed-batch layer that carries no ethanol term at all.
        """
        made = (design.final_biomass_g_per_l * design.max_volume_l
                - design.initial_biomass_g_per_l * design.initial_volume_l)
        return self.accumulated_g_per_l(made, design.max_volume_l)

    def plate_ceiling_g_per_l(self, glucose_g_per_l: float = PLATE_GLUCOSE_G_PER_L) -> float:
        """The most ethanol a well can physically contain, from its own sugar charge.

        Parameter-free: 2 mol ethanol per mol glucose is the Gay-Lussac stoichiometry, so at
        `context.py`'s default 20 g/L the ceiling is 111.01 mmol/L glucose -> 10.229 g/L
        ethanol. That is 3.1x to 7.0x BELOW the fed-batch band, which is the whole point --
        the vessel that makes the titre sits several-fold above the highest concentration the
        instrument this project owns could ever produce.

        WHAT THIS DOES NOT SAY. It does not say the plate is provably safe. That claim needs
        an inhibition threshold, and ``E_theta_inhibition_threshold`` is a :data:`NOT_BUILT`
        refusal. What the code can honestly report is the ceiling; the verdict needs a
        measurement nobody here has.
        """
        if glucose_g_per_l < 0:
            raise ValueError(f"glucose must be non-negative, got {glucose_g_per_l}")
        return (glucose_g_per_l / float(GLUCOSE_G_PER_MOL)) * 2.0 * float(ETHANOL_G_PER_MOL)

    def growth_penalty_per_h(self, ethanol_g_per_l: float) -> float:
        """The refusal, at the point a caller would reach for it.

        Raises:
            RefusedValue: always, naming the five constants and the experiment that would
                supply them.
        """
        del ethanol_g_per_l
        raise RefusedValue(
            "there is no ethanol growth penalty in this layer, and that is a verdict rather "
            "than an omission. The law mu = mu_c - i*max(0, E - E_theta) needs FIVE "
            "unmeasured constants -- the slope i, the threshold E_theta, the functional form, "
            "an uptake sink q_ethanol_uptake and a stripping rate k_strip -- and all five are "
            "recorded in mech.population.NOT_BUILT with what would close each. Four "
            "independent audits each put the resulting ablation below the 0.146 floor. This "
            "diagnostic is LayerState.REPORTED: it says how much ethanol is there and refuses "
            "to say what it does. Experiment 1.4 (three offline enzymatic reads on the "
            "already-planned fed-batch) is what would change that")


# --------------------------------------------------------------------------------------
# Criterion (a): the two channels, on real rows, against the two measured floors
# --------------------------------------------------------------------------------------

PLATE_READ_H = 4.14
"""Length of the committed Synergy reads, hours. `generator/plate.py::PlateConditions.duration_h`
and `mech/state.py::PLATE_READ_4H` both carry it; pinned here by test rather than imported so
this module's ablation does not depend on the plate package."""

DECLARED_READ_H = 24.0
"""The read length the ablation is EVALUATED at. A DECLARED axis, exactly as `mech/burden.py`
declares a copy number: no read this long is committed as text in this repository, so this is a
question put to a model fitted on 4.14 h rows and not a measurement.

It is not an arbitrary length. `docs/research/XPT_INVENTORY.md` records
``24h_continuous_plasmidloss.xpt`` -- 23.5 h, 141 reads, 78 occupied wells, OD600 with RFP --
as "the only measurement in the project of a term, plasmid loss, that the reporter ODE assumes
away", and it is not committed. At the measured plate band 24 h buys 8.5 to 12.6 generations
against the 1.5 to 2.2 of the read this project has committed, which is the whole of why the
verdict below switches."""

_GATE1_TABLES = ("NewProtocol", "Replicate3")


def plate_chord_growth_rates() -> pd.DataFrame:
    """Per-well specific growth rates off the three committed NewProtocol gate-1 tables.

    ``mu = ln(od_fold_change) / 4.14 h`` on every well that passed gate 1 -- the same chord
    estimate `mech/state.py::PLATE_GROWTH_BAND_PER_H` is derived from, on the same instrument
    the 0.0117 /h floor was measured on.

    Returns:
        ``plate``, ``well``, ``growth_rate`` (1/h) and ``read_duration_h``, which
        :func:`chord_growth_observable` moves to the declared read length.

    Raises:
        FileNotFoundError: if the three tables are not in this checkout.
    """
    tables = sorted((paths.REPO_ROOT / "outputs").glob("g1_2026*.csv"))
    wanted = [t for t in tables if any(k in t.name for k in _GATE1_TABLES)]
    if len(wanted) != 3:
        raise FileNotFoundError(
            f"the three NewProtocol gate-1 tables are not in this checkout; found "
            f"{[t.name for t in wanted]}. They are what the growth floor and the plate growth "
            f"band were both measured on, and this ablation is not scoreable without them")
    frames = []
    for table in wanted:
        rows = pd.read_csv(table).query("passed")
        frames.append(pd.DataFrame({
            "plate": table.name,
            "well": rows["well"].astype(str),
            "growth_rate": np.log(rows["od_fold_change"].astype(float)) / PLATE_READ_H,
            "read_duration_h": PLATE_READ_H,
        }))
    return pd.concat(frames, ignore_index=True)


@lru_cache(maxsize=512)
def _fraction_grid(p_per_division: float, burden: float, mu_free_per_h: float,
                   initial_fraction: float, t_max_h: float,
                   n_points: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """``(t, F)`` on a fixed grid, cached, for the ``b > 0`` case that has no closed form.

    NOT through :func:`~ystwin.mech.integrate.integrate_window`, and the reason is stated
    rather than left to be noticed: this is the inner loop of a ``least_squares`` fit,
    re-entered once per trial ``mu_free`` on a read length that is a fit argument and not
    a declared :class:`~ystwin.mech.state.Window`. The pinned driver runs the tau/T audit
    on every call, which is the right thing to do once for a trajectory and the wrong
    thing to do a thousand times inside a fit. :func:`simulate_population` -- the public
    trajectory, and the only one whose output is reported as a time course -- does go
    through it. The same applies to :meth:`TwoPopulation.fraction_after_generations`,
    which integrates in the GENERATION domain that ``integrate_window`` cannot express.
    """
    system = TwoPopulation(p_per_division, burden)
    grid = np.linspace(0.0, t_max_h, n_points)
    solution = solve_ivp(
        lambda _t, y: [system.d_fraction_d_time(y[0], mu_free_per_h)],
        (0.0, t_max_h), [initial_fraction], t_eval=grid, method="BDF", rtol=1e-9, atol=1e-12)
    if not solution.success:
        raise RuntimeError(f"F(t) integration failed: {solution.message}")
    return (tuple(grid.tolist()), tuple(solution.y[0].tolist()))


def bearing_fraction_at(model: TwoPopulation, times_h, *, mu_free_per_h: float,
                        initial_fraction: float = INOCULUM_BEARING_FRACTION) -> np.ndarray:
    """``F(t)`` at the given times. Closed form at ``b = 0``, integrated and cached otherwise.

    At ``b = 0`` this is exactly ``F0 exp(-p mu t)``, which is
    :meth:`TwoPopulation.fraction_after_generations` written in wall-clock -- the same
    statement as ``theta = 1 - 2^-p`` and the reason the unit bridge is checkable.
    """
    t = np.atleast_1d(np.asarray(times_h, dtype=float))
    if np.any(t < 0):
        raise ValueError("times must be non-negative")
    if mu_free_per_h <= 0:
        raise ValueError(f"mu_free must be positive, got {mu_free_per_h}")
    f0 = _unit(initial_fraction)
    if model.burden == 0.0:
        return f0 * np.exp(-model.p_per_division * float(mu_free_per_h) * t)
    span = float(max(t.max(), 1e-6))
    grid, values = _fraction_grid(model.p_per_division, model.burden, float(mu_free_per_h),
                                  f0, span, 513)
    return np.interp(t, np.asarray(grid), np.asarray(values))


def mean_bearing_fraction(model: TwoPopulation, *, mu_free_per_h: float, duration_h: float,
                          initial_fraction: float = INOCULUM_BEARING_FRACTION) -> float:
    """Time-average of ``F`` over a read of ``duration_h``. What a chord growth rate sees.

    The observed chord rate over a read is ``mu_free (1 - b Fbar)``, so this is the whole of
    the segregation signal in the growth channel -- and it is exactly ``F0`` when ``p = 0``
    and ``b = 0``, which is why the reduced model is a constant.
    """
    if duration_h <= 0:
        raise ValueError(f"duration must be positive, got {duration_h}")
    f0 = _unit(initial_fraction)
    if model.burden == 0.0:
        rate = model.p_per_division * float(mu_free_per_h)
        if rate * duration_h < 1e-12:
            return f0
        return f0 * (-math.expm1(-rate * duration_h)) / (rate * duration_h)
    grid, values = _fraction_grid(model.p_per_division, model.burden, float(mu_free_per_h),
                                  f0, float(duration_h), 513)
    return float(np.trapezoid(np.asarray(values), np.asarray(grid)) / duration_h)


@dataclass(frozen=True)
class SegregatingCulture(FittableModel):
    """``mu_obs(T) = mu_free (1 - b Fbar(T))``: the chord growth rate a read of length T sees.

    One fitted scalar, ``mu_free``, and the two population constants DECLARED rather than
    fitted -- they are swept axes, and fitting them here would be scoring a sweep against the
    rows that set it. The fit is one-dimensional and nonlinear only because ``Fbar`` depends
    on ``mu_free`` through the generations the read buys.
    """

    population: TwoPopulation | None = None
    initial_fraction: float = INOCULUM_BEARING_FRACTION
    name: str = "two-population (F falls over the read)"
    parameter_names: tuple[str, ...] = ("mu_free",)
    response_column: str = "growth_rate"

    def __post_init__(self) -> None:
        if self.population is None:
            raise ValueError(
                "SegregatingCulture needs a TwoPopulation: p and b are declared axes and "
                "this model does not invent them")

    def _weights(self, data: pd.DataFrame, mu_free: float) -> np.ndarray:
        if "read_duration_h" not in data.columns:
            raise KeyError("a growth frame needs a read_duration_h column")
        return np.array([
            1.0 - self.population.burden * mean_bearing_fraction(
                self.population, mu_free_per_h=mu_free, duration_h=float(t),
                initial_fraction=self.initial_fraction)
            for t in data["read_duration_h"]], dtype=float)

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        mu_free = float(parameters["mu_free"])
        return mu_free * self._weights(data, mu_free)

    def refit(self, data: pd.DataFrame) -> Fit:
        y = np.asarray(data[self.response_column], dtype=float)
        start = float(np.mean(y)) or 0.1
        solution = least_squares(
            lambda x: self.predict(data, {"mu_free": float(x[0])}) - y, [start],
            bounds=([1e-6], [5.0]), max_nfev=2000)
        mu_free = float(solution.x[0])
        residual = self.predict(data, {"mu_free": mu_free}) - y
        return Fit(self.name, {"mu_free": mu_free}, 1, len(data),
                   float(np.sqrt(np.mean(residual ** 2))), bool(solution.status > 0))


@dataclass(frozen=True)
class ConstantGrowth(FittableModel):
    """``mu_obs = mu_free``: this block deleted, which is the criterion-(a) reduced model.

    It is the two-population system at ``F`` held at its initial value -- a constant bearing
    fraction rescales ``mu_free`` and nothing else, which is exactly why a plate read cannot
    see the block and a long one can.

    IT IS THE ABLATED MODEL AND NOT A TRANSCRIPTION OF `generator/culture.py`, which since
    2026-09-05 integrates a logistic against a substrate-limited carrying capacity. That
    changes the incumbent's biomass curve for a SUBSTRATE reason and still leaves it with
    nowhere to put an ``F(t)``: what is deleted here is the segregation weight and only that,
    with both models refit on the same real 4.14 h wells.

    AND THAT IS EXACTLY WHY THE 24 h GROWTH PASS DOES NOT IDENTIFY THIS BLOCK, which is a
    limit on the interpretation rather than on the ablation. Deleting one piece is the right
    ablation; it does not show the piece is the only thing that could produce the signal.
    MEASURED at the plate's own committed inoculum (`generator/plate.py`: OD 0.16 x 0.42
    gDCW/OD = 0.0672 g/L), `generator/culture.py::simulate_culture` ALREADY predicts the chord
    rate falling from 0.3155 /h at 4.14 h to 0.1869 /h at 24 h -- a drop of 0.1286 /h, which
    is 11.0x the 0.0117 floor and 4.3x this block's entire 24 h effect, in the SAME direction
    on the SAME observable. So a 24 h OD trace could not attribute its falling chord rate to
    segregation: substrate exhaustion moves it several times further, and the repository
    already models that. The reporter channel does not have this problem -- there the
    competing decay is fitted as ``lambda_`` in BOTH models and what is scored is what
    survives that refit. Read the growth row as a sensitivity, and read
    :func:`score_reporter_ablation` as the channel that carries the claim.
    """

    name: str = "constant rate (incumbent)"
    parameter_names: tuple[str, ...] = ("mu_free",)
    response_column: str = "growth_rate"

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        return np.full(len(data), float(parameters["mu_free"]), dtype=float)

    def refit(self, data: pd.DataFrame) -> Fit:
        y = np.asarray(data[self.response_column], dtype=float)
        mu_free = float(np.mean(y))
        return Fit(self.name, {"mu_free": mu_free}, 1, len(data),
                   float(np.sqrt(np.mean((y - mu_free) ** 2))), True)


def chord_growth_observable(read_h: float = DECLARED_READ_H,
                            data: pd.DataFrame | None = None) -> Observable:
    """The chord growth rate a read of ``read_h`` would report, from wells read for 4.14 h.

    Both models are fitted to the same real per-well rates and then evaluated at the declared
    read length, which is the only axis on which the two differ: on the 4.14 h rows the full
    model's falling ``F`` is absorbed into its own ``mu_free`` and the two are indistinguishable
    by construction. Scored ABSOLUTELY, because 0.0117 /h is an absolute rate error.
    """
    frame = plate_chord_growth_rates() if data is None else data.copy()
    query = pd.DataFrame({"read_duration_h": [float(read_h)], "growth_rate": [np.nan]})

    def summarise(model: FittableModel, fit: Fit) -> float:
        return float(model.predict(query, fit.parameters)[0])

    return Observable(
        name=GROWTH_RATE_FLOOR.name, units="1/h", assay=GROWTH_RATE_FLOOR.assay,
        scored="absolute", data=frame, summarise=summarise,
        description=(f"chord specific growth rate over a {float(read_h):g} h read, from wells "
                     f"read for {PLATE_READ_H:g} h"))


def score_growth_ablation(p_per_division: float, burden: float, *,
                          read_h: float = DECLARED_READ_H,
                          initial_fraction: float = INOCULUM_BEARING_FRACTION,
                          data: pd.DataFrame | None = None) -> AblationResult:
    """Criterion (a) on the growth channel: delete the second population and see what moves.

    The verdict is a DURATION SWITCH and not a property of the block. At ``b = 0`` the effect
    is exactly zero however large ``p`` is -- two populations growing at the same rate make one
    OD trace -- so this channel is blind to segregation and sensitive only to the burden the
    segregation reveals.
    """
    observable = chord_growth_observable(read_h, data)
    full = SegregatingCulture(TwoPopulation(p_per_division, burden),
                              initial_fraction=initial_fraction)
    return ablate(full, ConstantGrowth(), observable, GROWTH_RATE_FLOOR)


@dataclass(frozen=True)
class SegregatingReporter(SaturatingReporter):
    """The incumbent's reporter algebra multiplied by ``F(t)``: what a MIXED culture reads.

    Bulk RFU/OD over a culture that is losing its plasmid is ``F(t) x R_per_bearing_cell(t)``,
    and `reporter.py` has nowhere to put the ``F(t)``. It is the same five fitted scalars, so
    the comparison cannot be won by spending a parameter, and it is exactly nested: at
    ``p = 0, b = 0`` the factor is the constant ``F0``, which rescales ``R0``, ``k_basal`` and
    ``k_max`` inside their own box.

    IT IS NOT A REPARAMETERISATION OF THE DECAY. ``F(t)`` multiplies the accumulation term as
    well as the initial condition, so ``R0 e^{-lt} + (k/l)(1 - e^{-lt})`` times ``e^{-p mu t}``
    is not any member of that family: the incumbent plateaus and this one turns over. The
    ablation measures exactly how much of that difference survives refitting.
    """

    population: TwoPopulation | None = None
    mu_free_per_h: float = 0.0
    initial_fraction: float = INOCULUM_BEARING_FRACTION
    name: str = "SaturatingReporter x F(t)"

    def __post_init__(self) -> None:
        if self.population is None or self.mu_free_per_h <= 0:
            raise ValueError(
                "SegregatingReporter needs a TwoPopulation and a positive measured mu_free; "
                "p, b and mu are declared here, never fitted")

    def _curve(self, t, d, R0, k_basal, k_max, K_dose, lambda_):
        base = SaturatingReporter._curve(self, t, d, R0, k_basal, k_max, K_dose, lambda_)
        return base * bearing_fraction_at(
            self.population, t, mu_free_per_h=self.mu_free_per_h,
            initial_fraction=self.initial_fraction)


def bearing_endpoint_observable(data: pd.DataFrame, *, read_h: float = DECLARED_READ_H,
                                dose_mM: float | None = None) -> Observable:
    """Dilution-corrected endpoint signal at the top dose, at a declared read length.

    NOT fold induction. ``F(t)`` divides out of a ratio taken at one time exactly, so fold
    induction is blind to this block by construction and scoring it would measure only how the
    fit moved. The endpoint signal is the number the reporter layer actually predicts and the
    one a segregating culture reads low.
    """
    for column in ("time_h", "dose_mM", "signal"):
        if column not in data.columns:
            raise KeyError(f"a reporter block needs time_h, dose_mM and signal; got "
                           f"{', '.join(map(str, data.columns))}")
    rungs = np.unique(np.asarray(data.dose_mM, dtype=float))
    top = float(rungs.max()) if dose_mM is None else float(dose_mM)
    if top not in set(rungs.tolist()):
        raise ValueError(f"dose {top:g} is not on this ladder")
    query = pd.DataFrame({"time_h": [float(read_h)], "dose_mM": [top], "signal": [0.0]})

    def summarise(model: FittableModel, fit: Fit) -> float:
        return float(model.predict(query, fit.parameters)[0])

    return Observable(
        name=REPORTER_ACTIVITY_FLOOR.name, units="RFU at the top dose",
        assay=REPORTER_ACTIVITY_FLOOR.assay, scored="relative", data=data,
        summarise=summarise,
        description=f"endpoint signal at {top:g} mM read at {float(read_h):g} h")


def score_reporter_ablation(export: str, construct: str, p_per_division: float,
                            burden: float, *, mu_free_per_h: float,
                            read_h: float = DECLARED_READ_H,
                            initial_fraction: float = INOCULUM_BEARING_FRACTION) -> AblationResult:
    """Criterion (a) on the reporter channel, on one committed (plate, construct) block.

    Args:
        export: A committed export name; :func:`ystwin.mech.ablation.reporter_block` reads it.
        construct: One of the four on those plates.
        p_per_division, burden: Declared points of the two swept axes.
        mu_free_per_h: The MEASURED plate growth rate. Not fitted and not asserted -- pass a
            point from :func:`plate_chord_growth_rates` or the plate band.
        read_h: The declared read length the endpoint is evaluated at.
        initial_fraction: F at inoculation.
    """
    data = reporter_block(export, construct)
    full = SegregatingReporter(population=TwoPopulation(p_per_division, burden),
                               mu_free_per_h=float(mu_free_per_h),
                               initial_fraction=initial_fraction)
    return ablate(full, SaturatingReporter(), bearing_endpoint_observable(data, read_h=read_h),
                  REPORTER_ACTIVITY_FLOOR)


def committed_reporter_blocks() -> tuple[tuple[str, str], ...]:
    """The (export, construct) blocks this repository can replay, in manifest order."""
    pairs = []
    for export in replay.exports():
        try:
            frame = replay.load_doses(str(export))
        except (ValueError, FileNotFoundError):
            continue
        for construct in sorted(pd.unique(frame.construct)):
            pairs.append((str(export), str(construct)))
    return tuple(pairs)


def burden_to_clear_the_growth_floor(p_per_division: float, *,
                                     read_h: float = DECLARED_READ_H,
                                     data: pd.DataFrame | None = None,
                                     steps: int = 25) -> float | None:
    """Smallest burden in the swept band whose growth ablation clears the floor, or ``None``.

    Scanned rather than solved: the effect is monotone in ``b`` but the fit is nonlinear, and
    a scan over the DECLARED band cannot report a value from outside it.
    """
    frame = plate_chord_growth_rates() if data is None else data
    low, high = BURDEN.bounds
    for step in range(steps + 1):
        b = low + (high - low) * step / steps
        if score_growth_ablation(p_per_division, b, read_h=read_h, data=frame).clears_floor:
            return b
    return None


def criterion_a_table(*, mu_free_per_h: float | None = None,
                      read_hours: tuple[float, ...] = (PLATE_READ_H, DECLARED_READ_H),
                      burdens: tuple[float, ...] = (0.0, 0.10, 0.234)) -> str:
    """Both channels scored at both read lengths, computed here rather than transcribed.

    Args:
        mu_free_per_h: The MEASURED plate rate the reporter channel converts ``p`` with.
            Defaults to the median of :func:`plate_chord_growth_rates`, which is the same
            0.3227 /h `mech/state.py::PLATE_GROWTH_BAND_PER_H` reports.
        read_hours: Read lengths to score at. The first is the one this project has run.
        burdens: Points of the swept burden band for the growth channel.
    """
    growth_rows = plate_chord_growth_rates()
    mu = float(np.median(growth_rows.growth_rate)) if mu_free_per_h is None \
        else float(mu_free_per_h)
    blocks = committed_reporter_blocks()
    lines = [
        f"CRITERION (a), measured: {len(growth_rows)} wells and {len(blocks)} reporter blocks, "
        f"mu = {mu:.4f} /h",
        f"  reporter floor {REPORTER_ACTIVITY_FLOOR.noise_floor:g} relative CV; growth floor "
        f"{GROWTH_RATE_FLOOR.noise_floor:g} /h absolute",
    ]
    for read in read_hours:
        lines.append(f"  read {read:g} h")
        for p in LOSS_PER_DIVISION.bounds:
            ratios = [score_reporter_ablation(export, construct, p, 0.0, mu_free_per_h=mu,
                                              read_h=read).ratio
                      for export, construct in blocks]
            lines.append(
                f"    reporter  p={p:.4f}  {min(ratios):.3g}-{max(ratios):.3g}x the floor, "
                f"clears {sum(r > 1.0 for r in ratios)}/{len(ratios)}")
        for b in burdens:
            ratios = [score_growth_ablation(p, b, read_h=read, data=growth_rows).ratio
                      for p in LOSS_PER_DIVISION.bounds]
            lines.append(
                f"    growth    b={b:.3f}      {min(ratios):.3g}-{max(ratios):.3g}x the floor, "
                f"clears {sum(r > 1.0 for r in ratios)}/{len(ratios)}")
    return "\n".join(lines)


# --------------------------------------------------------------------------------------
# The fed-batch answer: when the plasmid-free population takes over the vessel
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class FedBatchCollapse:
    """What a fed-batch loses to segregation, in the units a run is planned in.

    Args:
        generations: Doublings the run spends. The vessel's own budget, not its wall-clock:
            `fba/fedbatch.py::FedBatchDesign.generation_budget_at` computes it from the feed
            and the inoculum, and it is the quantity dF/dg is invariant in.
        final_fraction: F when the vessel is full.
        cell_hour_fraction: ``(1/T) int F dt``. The fraction of cell-hours that carried the
            plasmid, which is what multiplies a constant per-cell rate.
        biomass_weighted_fraction: ``int F X dt / int X dt`` over THIS run's own biomass.
            THE ONE THAT MULTIPLIES A PER-CELL PRODUCT RATE: biomass grows through the run, so
            the last generations dominate the product and they are the ones with the fewest
            plasmids left. Read it as the fraction of biomass-hours that carried the plasmid.
            It equals the ratio to a stable population's titre only at ``b = 0``, where the two
            cultures share the same ``X(t)``; at ``b > 0`` the plasmid-free cells also make the
            culture grow faster, and this number does not credit that extra biomass.
        generations_to_half: Doublings until half the population has lost the plasmid.
        hours_to_half: The same, at the growth rate supplied. Wall-clock, and a feeding
            strategy moves this and never moves ``generations_to_half``.
    """

    generations: float
    final_fraction: float
    cell_hour_fraction: float
    biomass_weighted_fraction: float
    generations_to_half: float
    hours_to_half: float

    @property
    def titre_retained(self) -> float:
        """The biomass-weighted fraction, named for what it multiplies in the outcome panel."""
        return self.biomass_weighted_fraction

    def summary(self) -> str:
        return (
            f"over {self.generations:.2f} generations F falls to {self.final_fraction:.3f}; "
            f"biomass-weighted mean F = {self.biomass_weighted_fraction:.3f}, so a per-cell "
            f"product rate is realised on {100 * self.biomass_weighted_fraction:.1f}% of the "
            f"biomass-hours; half the culture is plasmid-free after "
            f"{self.generations_to_half:.2f} generations ({self.hours_to_half:.1f} h here)")


def _wall_clock_over_generations(model: TwoPopulation, initial_fraction: float,
                                 mu_free_per_h: float, grid: np.ndarray):
    """``(F, hours)`` on a generation grid. ``dt = ln2 dg / mu_pop``, trapezoid in g.

    Factored out so ``hours_to_half`` is the SAME integral whether or not the run was long
    enough to reach it. The closed form it replaced, ``g ln2 / (mu_free (1 - b))``, holds F
    at 1 and so overstates the wall-clock at b = 0.234 by 7.6% at the bottom of the loss band
    and 8.5% at the top -- MEASURED against this integral at F0 = 0.95.
    """
    fractions = np.array([model.fraction_after_generations(initial_fraction, float(g))
                          for g in grid])
    rates = float(mu_free_per_h) * (1.0 - fractions * model.burden)
    hours = np.zeros_like(grid)
    hours[1:] = np.cumsum(_LN2 * np.diff(grid) / (0.5 * (rates[1:] + rates[:-1])))
    return fractions, hours


def fedbatch_collapse(model: TwoPopulation, *, generations: float,
                      mu_free_per_h: float = FEDBATCH_GROWTH_PER_H,
                      initial_fraction: float = INOCULUM_BEARING_FRACTION,
                      n_points: int = 2001) -> FedBatchCollapse:
    """Run the two-population system over a fed-batch and report what it costs.

    THE CONDITION FOR TAKEOVER IS NOT A THRESHOLD, which is the finding: in any non-selective
    medium ``s = -b/(1-b) <= 0`` and ``p > 0``, so both terms of ``dF/dt`` are negative for
    every ``F`` in ``(0, 1]`` and ``F = 0`` is the only stable state -- see
    :meth:`TwoPopulation.takeover`. What a design controls is not WHETHER but HOW FAR, and the
    axis is GENERATIONS: ``dF/dg`` has no ``mu`` in it, so a slower feed buys wall-clock and
    buys no plasmids. The lever is the vessel's generation budget, which
    `fba/fedbatch.py::FedBatchDesign.generation_budget_at` shows is set by the inoculum.

    Args:
        model: The system.
        generations: Doublings the run spends. Pass the design's own budget.
        mu_free_per_h: The MEASURED unburdened rate; sets wall-clock only.
        initial_fraction: F at inoculation, after the seed train.
        n_points: Output grid for the two integrals.
    """
    if generations <= 0:
        raise ValueError(f"generations must be positive, got {generations}")
    grid = np.linspace(0.0, float(generations), int(n_points))
    fractions, hours = _wall_clock_over_generations(
        model, initial_fraction, mu_free_per_h, grid)
    biomass = 2.0 ** grid
    cell_hours = float(np.trapezoid(fractions, hours) / hours[-1])
    weighted = float(np.trapezoid(fractions * biomass, hours)
                     / np.trapezoid(biomass, hours))
    if model.p_per_division == 0.0 and model.burden == 0.0:
        half = math.inf
    else:
        half = model.generations_to_fraction(initial_fraction, 0.5 * initial_fraction)
    if math.isinf(half):
        hours_to_half = math.inf
    elif half <= grid[-1]:
        hours_to_half = float(np.interp(half, grid, hours))
    else:
        far = np.linspace(0.0, float(half), int(n_points))
        hours_to_half = float(_wall_clock_over_generations(
            model, initial_fraction, mu_free_per_h, far)[1][-1])
    return FedBatchCollapse(
        generations=float(generations),
        final_fraction=float(fractions[-1]),
        cell_hour_fraction=cell_hours,
        biomass_weighted_fraction=weighted,
        generations_to_half=float(half),
        hours_to_half=hours_to_half)


# --------------------------------------------------------------------------------------
# The two unfitted convergences, computable so they are checkable
# --------------------------------------------------------------------------------------

def cheng_hohnholz_convergence() -> dict[str, float]:
    """Cheng's 1992 FITTED ``p`` against Hohnholz's 2017 REPLICA-PLATED ``theta``.

    Two laboratories, twenty-five years apart, different vector, different marker, different
    method -- one a curve fit to a batch fermentation, the other a colony count on agar. They
    agree only after the unit conversion, and that is the point of reporting it:
    ``theta(0.08) = 5.394e-2`` sits inside the measured ``[1.0e-2, 5.7e-2]``, while the
    unconverted 8.0e-2 sits outside it.

    AND IT SURVIVES THE SECOND AMBIGUITY IN THE SAME SYMBOL, which the harvest's verifier
    raised and priced at 39%. Cheng writes the Imanaka-Aiba kernel discretely, as
    ``X+ -> (2 - p) X+ + p X-``: read literally, one synchronous division retains
    ``1 - p/2``, so ``theta = p/2 = 0.500 p``, while the continuous ODE that is meant to
    represent it gives ``theta = 1 - 2^-p ~ 0.693 p``. The two readings of Cheng's fitted
    0.08 are 4.000e-2 and 5.394e-2, and **both land inside the replica-plated band** while
    the unconverted 8.0e-2 lands outside under either -- so the convergence is a fact about
    the number and not about which kernel is read into it.

    IT TOUCHES NO CONSTANT THIS MODULE BUILDS. ``p`` here is DEFINED by this module's own
    ODE, and :data:`LOSS_PER_DIVISION` is Hohnholz's measured ``theta`` pushed through
    :func:`divisions_per_generation`, the exact inverse of that same definition. The
    ambiguity can only bite on an imported point estimate, which is the one place it is
    reported and the one place it is not consumed.

    STATUS: a convergence of two independent parameter estimates, not a criterion-(c) pass of
    the assembled model. This repository holds no dynamic plasmid-loss dataset the model could
    be scored against -- bulk RFU/OD over a mixed culture is ``F(t) x R_per_bearing_cell(t)``,
    one number for two time-varying unknowns, so the team's own 24 h plate cannot separate
    them as collected.
    """
    p_cheng = 0.08
    theta = loss_per_generation(p_cheng)
    discrete = 0.5 * p_cheng
    low, high = min(MEASURED_LOSS_PER_GENERATION), max(MEASURED_LOSS_PER_GENERATION)
    return {
        "cheng_p_per_division": p_cheng,
        "cheng_theta_per_generation": theta,
        "cheng_theta_discrete_kernel": discrete,
        "hohnholz_low": low,
        "hohnholz_high": high,
        "converted_lands_inside": float(low <= theta <= high),
        "discrete_reading_lands_inside": float(low <= discrete <= high),
        "unconverted_lands_inside": float(low <= p_cheng <= high),
        "kernel_readings_differ_by": theta / discrete - 1.0,
    }


def foci_prediction_overlap() -> dict[str, float]:
    """The microscopy-to-agar prediction, and how much of the measured band it covers.

    ``m = 3-5`` foci per haploid nucleus, measured by fluorescence microscopy, put through an
    independent-partition kernel with nothing fitted, predict ``theta`` in
    ``[9.77e-4, 1.56e-2]`` per generation. Replica plating in a different laboratory measured
    ``[1.0e-2, 5.7e-2]``. The bands overlap on ``[1.0e-2, 1.56e-2]``, i.e. only the ``m = 3``
    end of the microscopy count reaches the measured band at all -- which is a weaker
    convergence than the per-DIVISION form ``2^(1-2m)`` would suggest, and that form is not
    the quantity replica plating reports.

    STATUS: a hypothesis with one unmeasured structural parameter -- ``m`` is measured for
    STB-reporter plasmids in the native context, not for an engineered HIS3 YEp. Reported so
    the arithmetic is checkable, not claimed as a result.
    """
    low_m, high_m = FOCI_PER_HAPLOID_NUCLEUS
    predicted = (foci_partition_loss_per_generation(high_m),
                 foci_partition_loss_per_generation(low_m))
    measured = (min(MEASURED_LOSS_PER_GENERATION), max(MEASURED_LOSS_PER_GENERATION))
    overlap_low = max(predicted[0], measured[0])
    overlap_high = min(predicted[1], measured[1])
    return {
        "predicted_low": predicted[0],
        "predicted_high": predicted[1],
        "measured_low": measured[0],
        "measured_high": measured[1],
        "overlaps": float(overlap_high > overlap_low),
        "overlap_low": overlap_low,
        "overlap_high": overlap_high,
    }


def provenance_summary() -> str:
    """One block of text: both gates, the takeover verdict, and the diagnostic's two numbers."""
    band = ethanol_yield_band()
    diagnostic = EthanolDiagnostic(band[0], band[1])
    lines = [
        f"POPULATION  {POPULATION_GATE.summary()}",
        f"ETHANOL     {ETHANOL_GATE.summary()}",
        f"            plate ceiling {diagnostic.plate_ceiling_g_per_l():.3f} g/L at "
        f"{PLATE_GLUCOSE_G_PER_L:g} g/L glucose; r = {band[0]:.4f}-{band[1]:.4f} g/g",
        f"REFUSED     {', '.join(sorted(p.name for p in NOT_BUILT.refusals()))}",
    ]
    return "\n".join(lines)
