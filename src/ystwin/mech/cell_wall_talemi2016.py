"""Talemi 2016's yeast osmo-stat, transcribed BY HAND IN THIS REPOSITORY for its Slt2 arm.

WHAT THIS FILE IS. Talemi SR, Tiger CF, Andersson M, Babazadeh R, Welkenhuysen N, Klipp E,
Hohmann S & Schaber J 2016, "Systems Level Analysis of the Yeast Osmo-Stat", Sci Rep 6:30950
(doi:10.1038/srep30950, PMID 27515486, PMC4981887) ships its selected model as BioModels
**MODEL1606100000**, an SBML Level 2 Version 4 document with 11 species, 12 reactions, 100
global parameters, 8 function definitions, 29 assignment rules and 3 rate rules. That file is
vendored and checksummed at ``data/native_reference_models/talemi2016/``, sha256
``c4937f14...835e5dc``, under **CC0 1.0** -- the deposit's own RDF carries the BioModels
dedication of "all copyright and related or neighbouring rights to this encoded model ... to the
public domain". The article and its supplementary files are CC BY 4.0. What follows is a hand
transcription of the deposit into this engine's idiom. **The executable is ours**, which is why
it carries ``AVAILABILITY = "transcribed_here"`` and never ``local_verified`` -- that level in
this repository means a checksummed publisher file loaded unchanged by the shared SBML loader,
and the shared loader refuses this one (:data:`SBML_LOADER_REFUSAL`).

Nothing here was typed from memory. Every constant, initial amount, compartment size, switch and
stoichiometric coefficient is re-read out of the vendored SBML by :func:`audit_deposit` at
construction and compared to the literal written here; a mismatch raises rather than running.

THE SIGN IS THE FIRST THING, BECAUSE IT READS BACKWARDS. This axis is **not** a cell-wall-damage
sensor. In this source the CWI pathway mediates **HYPO**-osmotic adaptation: Slt2 activation is
driven by cell volume through ``v9 = k9 * Vos``, so it rises when the cell **swells** and falls
when it shrinks, while HOG answers the opposite sign. The article states it directly -- a fixed
threshold gives "basal activation at the initial volume and quasi-linear activation upon volume
decrease below or increase above the initial volume, i.e., the activation threshold for the HOG
and the CWI pathway, respectively". :data:`SIGN_OF_THE_AXIS` carries that sentence, and
:func:`reproduction_report` computes both directions so the claim is arithmetic and not a
paraphrase: across ``Vos/Vos0`` from 0.90 to 1.10 the Slt2 rest point rises 10.7% -> 124.7% of
maximal phosphorylation while the HOG signal falls 0.1240 -> 0.0055. Anyone wiring this axis to
a cell-wall-damaging agent would be wiring it to the wrong sign of the wrong stimulus.

WHY THE AXIS IS NAMED ``cell_wall_slt2`` AND NOT ``cell_wall``. `generator/stress_panel.py`
lines 247-249 declare the cell_wall module as **Slt2/Mpk1 -> Rlm1** against an **RLM1 box**, with
``congo_red`` and ``caffeine`` as its stressors. The deposit contains none of that. ``grep -ci``
over the pinned SBML returns **0** for each of Rlm1, RLM1, congo, caffeine, calcofluor,
caspofungin, zymolyase, SDS, Pkc1, Bck1, Wsc1, Mid2, "damage" and "transcription". The model
terminates at phosphorylated Slt2 as a kinase, with no output arm at all. This follows the
``carbon_camp_input`` precedent exactly: an axis that transcribes one END of a declared module
gets its own name so it cannot be read as having closed the module. :data:`NOT_DRIVEN` refuses
the missing output and the missing input by name, and adopting this transcription retires
neither.

WHAT REPRODUCES, stated before the code because it is the result. Four independent checks, all
recomputed by :func:`reproduction_report` and asserted by
``tests/test_cell_wall_talemi2016_transcription.py``. Nothing here was tuned.

1. **Table S6's published closed forms reproduce the deposit's stored initial state to
   floating-point identity.** Table S6 prints the initial conditions as formulas, not numbers --
   ``Slt2PP|0 = 3230 * f_N2uM * [Slt2PP]0 * maxHog1nucf`` with ``[Slt2PP]0 = 0.246`` from data.
   Evaluated here that is 0.0357820196004481 uM against a deposited 1.05556957821322/29.5 =
   0.0357820196004481, a relative 1.7e-15. The same holds for Slt2, Hog1, Hog1PP, Fps1, Fps1P,
   Vos0, f_N2uM, ci0, cin0, VP0 and the Area -- twelve rows, worst relative 1.7e-15.
2. **Table S6's two published steady-state roots reproduce.** ``Slt2Signal|0`` and
   ``Hog1Signal|0`` are printed as the positive root of a quadratic; evaluated here they give
   0.24492386146460 and 0.01404533383442 against deposited 0.244923861463865 and
   0.0140453338344206, relative 3.0e-12 and 1.2e-16. The 3.0e-12 is the deposit's own root
   solve, not a disagreement about the formula.
3. **Table S8's five "Calculated" rows reproduce from their published steady-state definitions.**
   k4, km7, k6b, k13 and k15 are printed as formulas rather than values; recomputed here, four of
   them give the deposit's stored numbers to a worst relative **1.9e-15**. The fifth, k13, lands
   at 3.0e-12 and is audited against a looser tolerance for a stated reason and not a convenient
   one: k13 is the only calculated row that reads ``Slt2Signal|0``, and it inherits that root's
   own 3.0e-12 exactly. This check matters more than it looks -- it exercises the Fps1, glycerol,
   Hog1 and Slt2 rate laws simultaneously, because each formula is a different one set to zero.
4. **The deposited initial state is a fixed point of eleven of the twelve transcribed states**,
   at a worst relative rate of **2.4e-14 /s**. The twelfth is external glycerol and it is not a
   defect of the transcription -- see the fourth source defect below.

AND THE PUBLISHED EXPERIMENT REPRODUCES. Table S7 publishes the protocol itself: S1 = 800000 uM
sorbitol, S2 = 270000 uM, ``toff`` = 120 s, ``ts`` = 840 s, mixing time ``tm`` = 10 s. Run at
those five numbers -- the source's own defaults, none of them chosen here -- this transcription
gives the paper's central result: relative Slt2 phosphorylation rises from its 24.6% baseline to
**89.2%**, peaking 86 s after the dilution, while relative Hog1PP saturates at ~100% on the
hyper-shock. The post-dilution volume reaches 1.012 x V0, i.e. above the pre-stress volume,
which is precisely the condition Figure 1c says is required to prompt a hypo-osmotic shock.
Sweeping the dilution time reproduces the article's structure across the whole panel:

    dilution at     45 s     90 s      4 min    14 min    30 min    45 min
    V_max / V0     0.983    0.983      0.983     1.012     1.036     1.037
    Slt2P peak     25.1%    25.8%      38.9%     89.2%     84.7%     84.1%

The 45 s and 90 s columns are flat against a 24.6% baseline and the 14/30/45 min columns are
large, which is the article's own reading of Figure 1a-c. The **4 min column is the one the
paper is about**: it shows a real Slt2 peak (1.58x baseline) at a volume that never exceeds V0,
and Table 1's "4MiP" column records that only the models carrying the sensitizer could produce
it. This transcription is the sensitized model and it produces it.

THE SENSITIZER SWITCH IS NOT A MODEL COMPARISON, AND THIS SAYS SO RATHER THAN IMPLYING IT. The
deposit carries ``SW_Sensitizer``, and setting it to 0 does delete v14, v15 and v16 and collapse
v11 to Table S5's non-sensitized ``V11-a``. It does **not** reproduce the paper's non-sensitized
models, which were re-estimated. At the deposit's own constants the unsensitized signal has no
comparable rest point: Table S6's own non-dagger ``Slt2Signal|0`` row evaluates to **11913**,
five orders of magnitude above the sensitized 0.2449, because vmax11 = 0.713 alone cannot balance
v9 = 2.278 and the balance falls to the linear k10 = 1.31e-4 term instead. The signal then sits
~660000-fold above km11, is saturated, and cannot respond to anything: the 4', 14' and 30' peaks
come out at 24.2%, 23.7% and 23.6% against a 24.6% baseline, i.e. flat and slightly *negative*.
That is directionally the paper's conclusion and it is arrived at by a different mechanism, so
:func:`sensitizer_switch_comparison` returns it labelled as a switch experiment and never as a
reproduction of Table 1's ranking.

FIVE DEFECTS IN THE SOURCE, FOUND BY TRANSCRIBING IT. None is repaired here.

* **The accession the paper prints does not exist.** The Supplementary "Simulation Instructions"
  says the selected model is in BioModels as "MODEL1604100004".
  ``https://www.ebi.ac.uk/biomodels/MODEL1604100004`` answers **HTTP 404**. The model is
  MODEL1606100000, which no page of the article names. The deposit's identity was therefore
  established against the authors' own Supplementary S2 SBML instead, and it holds exactly: same
  compartments, same 11 species with identical ``initialAmount``, ``boundaryCondition`` and
  ``hasOnlySubstanceUnits``, same 100 parameters, and **13 strings differ with zero values
  differing** -- every one a scientific-notation reformat. :data:`DEPOSIT_IDENTITY` records it
  and :func:`audit_deposit` re-runs it against the second vendored file.
* **Table S8 and the deposit disagree about vmax16, and the main text takes Table S8's side.**
  Table S8 prints ``vmax16 = 1``, method "set", and the Results say "except parameter k16 which
  we set to one". The deposit writes ``v16_k = 3.26314``. Every other numeric Table S8 row --
  k0, k1, vmax2, km2, k3, k5, k6, vmax7, k8, k9, k10, vmax11, km11, k12, k14, ki16, n16 -- equals
  the deposit exactly. So it is 17 of 18. The disagreement is measured rather than deplored:
  under the published protocol the Slt2 peak moves 89.152% -> 89.801% at 14 min, 38.948% ->
  39.150% at 4 min, and not at all at 30 min. :data:`TABLE_S8_VMAX16` keeps the published reading
  runnable and :func:`reproduction_report` reports both.
* **v16 is a floor, which is why nobody noticed.** With n16 = 999.898 the term is a step: it is
  off above ki16 = 3.08897 and full on below it, so it does not set a rate, it pins the
  Sensitizer's minimum. Across the published protocols the Sensitizer bottoms out at 3.1109 with
  the deposited vmax16 and 3.1073 with the published 1 -- both just above the same floor. Delete
  the reaction entirely and the floor goes to 2.4856 and the 14' peak to 115.7%, so v16 is
  strongly load-bearing while *its own constant* is nearly inert. Those are different statements
  and this module keeps them apart.
* **External glycerol has no removal term, so the unstressed model has no fixed point.** Glyex is
  a dynamic species produced by v8 and consumed by nothing, and at the deposited state v8 is not
  zero: 4542.6 amount/s, because Glyin/Vos = 180000 uM stands against Glyex/Vex = 1800 uM. Table
  S6 says as much ("the external glycerol acts basically as a sink"). The cost is bounded and
  measured rather than assumed: over 24 unstressed hours Glyex/Vex goes 1800 -> 9646 uM against
  a ce0 of 260000, and the eleven other states settle and stay -- Vos moves 1.5e-5 relative and
  relative Slt2 phosphorylation 24.6000% -> 24.5993%. Nothing here is ever called a steady state
  of the whole system; the residual is reported per state.
* **The VP assignment rule is dead arithmetic, and Table S7 shows what it should be.** Table S7
  defines the non-turgid volume once, as ``VP0 = V0*exp(-P0/eps)``. The deposit adds a second
  rule ``VP = V0*(V0/VP0)^(-ModelValue_7/eps)`` in which ``ModelValue_7`` and ``eps`` are two
  copies of 14.3, so the exponent is -1 and the rule is the identity on VP0. It is carried here,
  and :func:`audit_deposit` checks that the two copies are still equal -- if either ever moved
  the rule would silently stop being the identity and the turgor reference would shift.

TWO SMALLER THINGS, RECORDED BECAUSE THEY COST READING TIME. The article calls the CWI driver
``CWISignal`` and every table and the deposit call it ``Slt2Signal``; the tables write ``Slt2PP``
for the species the deposit names ``Slt2P``. And Table S6's header says concentrations are in
"(umoL/fL)" while Table S8's says "umol/L"; the numbers settle it -- f_N2uM = 1e21/(mol*Vos0)
converts a molecule count into uM at Vos0 in fL, so uM is right and Table S6's header is a typo.

THE ABUNDANCE THIS ARM IS SCALED BY HAS NO SOURCE AT ALL. Table S7 gives Hog1t = 0.3821 uM
("6788 * f_N2uM") and Fps1t = 0.051 uM ("907 * f_N2uM"), both cited to
``http://yeastgfp.yeastgenome.org/`` -- a URL, not a paper, with no interval. The **Slt2 total
does not appear in Table S7 at all**: its 3230 molecules occur only inside Table S6's formula,
uncited. So the scale of the entire arm this module transcribes is the one abundance the source
does not source. All three are registered ASSERTED for that reason and none is graded MEASURED.

THE FREE-SCALAR GATE REFUSES THIS PIECE, AND SHOULD. Forty-three free scalars against **zero**
independent targets. Seventeen of Table S8's rows record their own method as "Estimated", fitted
in COPASI to the Hog1/Slt2/volume/glycerol time series; the published readouts this module
reproduces -- the 24.6% and 6.1% baselines, the 4'/14'/30' Slt2 peaks -- are the same data, so
they are registered ``fitted=True`` and cannot score the model. Table S7's biophysics (Lp, P0,
eps, fmin) is imported from Schaber 2010 and is genuinely external, but it is external to the
CWI arm's constants, not a target for them. :func:`gate` returns that verdict instead of hiding
it. A faithful transcription of a fitted model is evidence about the transcription.
"""
from __future__ import annotations

import hashlib
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from scipy.integrate import solve_ivp

from .. import paths
from .contracts import Compartment
from .integrate import FALLBACK_METHOD, PINNED_METHOD
from .params import Param, ParamRegistry, Target, provenance_table

__all__ = [
    "ARTIFACTS",
    "ARTIFACTS_NOT_VENDORED",
    "AVAILABILITY",
    "AXIS",
    "CALCULATED_ROWS",
    "DEPOSIT_IDENTITY",
    "DEPOSIT_SHIPS_UNSTRESSED",
    "DEPOSIT_SWITCHES",
    "LICENCE",
    "NATIVE_TIME_UNIT",
    "NOT_DRIVEN",
    "OsmoStat",
    "PANEL_CHANNEL",
    "PUBLISHED",
    "SBML_LOADER_REFUSAL",
    "SECONDS_PER_HOUR",
    "SIGN_OF_THE_AXIS",
    "SOURCE",
    "SOURCE_DEFECTS",
    "STATE_NAMES",
    "STATE_VARIABLES",
    "Shock",
    "TABLE_S8_VMAX16",
    "TALEMI_PARAMS",
    "audit_deposit",
    "calculated_constants",
    "gate",
    "initial_state",
    "provenance",
    "reproduction_report",
    "sensitizer_switch_comparison",
    "volume_response",
]

AVAILABILITY = "transcribed_here"
"""The level this module carries, and the reason it is not ``local_verified``.

``local_verified`` in this repository means a checksummed publisher file loaded unchanged by
``mech/kinetic_sbml.py``. That loader refuses this file (:data:`SBML_LOADER_REFUSAL`), so the
executable below is ours and inherits our bugs. Nothing here may be promoted past this level.
"""

AXIS = "cell_wall_slt2"
"""The axis name, deliberately not ``cell_wall``. See the module docstring and :data:`NOT_DRIVEN`.

The ``carbon_camp_input`` precedent: an axis that transcribes one end of a panel module is named
for the end it has, so that adopting it cannot be read as having closed the module.
"""

SOURCE = ("Talemi SR, Tiger CF, Andersson M, Babazadeh R, Welkenhuysen N, Klipp E, Hohmann S, "
          "Schaber J 2016, 'Systems Level Analysis of the Yeast Osmo-Stat', Sci Rep 6:30950, "
          "doi:10.1038/srep30950, PMID 27515486, PMC4981887; BioModels MODEL1606100000")

LICENCE = (
    "Two licences. The encoded model is CC0 1.0 "
    "(https://creativecommons.org/publicdomain/zero/1.0/), declared in the deposit's own RDF: "
    "'To the extent possible under law, all copyright and related or neighbouring rights to this "
    "encoded model have been dedicated to the public domain worldwide.' The article and its "
    "supplementary files are CC BY 4.0, declared in the article's own <permissions> block. "
    "Redistribution of both is therefore permitted, the second with attribution. This is a "
    "stronger position than jalihal2021/, which carries no licence and is custody-only")

ARTIFACTS = MappingProxyType({
    "sbml_biomodels":
        ("data/native_reference_models/talemi2016/MODEL1606100000_url.xml",
         "c4937f14c2d8e36b13a449b5ac9771d4e17426543298d7a2c9ca9c4df835e5dc"),
    "sbml_authors_supplementary_s2":
        ("data/native_reference_models/talemi2016/srep30950-s2-OsmoStat_Final.xml",
         "6d18c0ddaf42c9ea3d5d18d42e14157b833be5284ed639fd8d9d85e7a9edfde3"),
    "article_xml":
        ("data/native_reference_models/talemi2016/PMC4981887.xml",
         "ff98279a6e3d443c761b5395dbb746f793a57ae5d0c27171d62676a03f7272c6"),
    "tables_s4_s8_text":
        ("data/native_reference_models/talemi2016/srep30950-s1-p21-33-tables-S4-S8.txt",
         "15d4fc1d551ba8eadab9a53be86f06ecd27b557245c90f5765941fe0e96e2175"),
    "calcofluor_module_text":
        ("data/native_reference_models/talemi2016/srep30950-s1-p06-07-calcofluor-module.txt",
         "93e1b248685239cf6a8d2ca5351a29bdd55e22b2caef540dc9c5daf5fdfb88b7"),
})
"""The five vendored files this transcription was made from, by content.

The last two are DERIVED and ours: ``pdftotext -layout`` (poppler 26.04.0) over the page ranges
named in their filenames. They carry Tables S4-S8 and the calcofluor module so a reader can check
the rows quoted above without a 12.5 MB download; they are corroboration, not the anchor. The
anchor is :data:`ARTIFACTS_NOT_VENDORED`, plus the fact that every Table S6 and S8 row quoted is
independently confirmed by recomputing a 15-significant-figure deposit value from it.
"""

ARTIFACTS_NOT_VENDORED = MappingProxyType({
    "supplementary_pdf": (
        "srep30950-s1.pdf, 12491506 bytes, sha256 "
        "97e1c32b05f447767a40f8a51a00684b5f6885c03b7647ea049baadd85030a91; inside "
        "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4981887/supplementaryFiles. Not "
        "copied into the repository: at 12.5 MB it is larger than all of "
        "data/native_reference_models put together and would become the largest tracked file "
        "here. Pinned by hash so the two text extracts above can be regenerated and diffed"),
})
"""Upstream files pinned by hash and deliberately not vendored, with the reason for each."""

DEPOSIT_IDENTITY = (
    "The Supplementary Information's 'Simulation Instructions' cites BioModels accession "
    "MODEL1604100004, which answers HTTP 404; the model is MODEL1606100000, named nowhere in the "
    "article. Identity is therefore established by comparison against the authors' own "
    "Supplementary S2 SBML: same 4 compartments, same 11 species with identical initialAmount, "
    "boundaryCondition and hasOnlySubstanceUnits, same 100 global parameters, and 13 parameter "
    "STRINGS differing with ZERO values differing -- every one a scientific-notation reformat "
    "such as '1E-9' against '1e-09'. audit_deposit() re-runs this comparison")
"""Why a second SBML is vendored, and what the comparison found."""

SBML_LOADER_REFUSAL = (
    "mech/kinetic_sbml.py:332 -- UnsupportedSBMLError('rate and algebraic rules are unsupported; "
    "only parameter assignment rules are supported'). Three rateRules trigger it (Vos, HOGSignal, "
    "Slt2Signal) and, unlike williamson2009's seven, NONE of these is a literal zero: they are "
    "the model's three central ODEs. There is nothing to audit away here and nothing was "
    "stripped; the loader is right to refuse the file and this module transcribes it instead")
"""The shared loader refuses this deposit, for a stronger reason than the last one it refused."""

NATIVE_TIME_UNIT = "second"
"""The source's time base, kept internally.

Two independent confirmations: the SBML declares no time unit, so SBML Level 2's default second
applies; and Table S8's own header says "mass is in grams and time in seconds". Volumes are fL
and concentrations uM throughout.
"""

SECONDS_PER_HOUR = 3600.0
"""The conversion :meth:`OsmoStat.timecourse` applies once, at the boundary.

Three published time bases collide in this architecture (min^-1, s^-1, h^-1) and the resulting
60x has already cost errors. The published protocol row :class:`Shock` keeps the source's
seconds verbatim, because renaming Table S7's own toff/ts/tm would make that check unreadable.
"""

SIGN_OF_THE_AXIS = (
    "CWI answers a volume INCREASE and HOG answers a volume DECREASE. Talemi 2016, Results: a "
    "fixed threshold 'translates into basal activation at the initial volume and quasi-linear "
    "activation upon volume decrease below or increase above the initial volume, i.e., the "
    "activation threshold for the HOG and the CWI pathway, respectively.' Mechanically it is "
    "v9 = k9 * Vos: the Slt2 driver is proportional to osmotically active cell volume. This axis "
    "is therefore a HYPO-osmotic adaptation arm and MUST NOT be described, wired or dosed as a "
    "cell-wall-damage sensor")
"""The counterintuitive direction, quoted from the source and checked by :func:`volume_response`."""

_SBML_NS = "http://www.sbml.org/sbml/level2/version4"
_MATHML_NS = "http://www.w3.org/1998/Math/MathML"

_TABLE_S7 = "Talemi 2016 Table S7, 'Auxiliary variables, physical quantities and their Definition/value'"
_TABLE_S8 = "Talemi 2016 Table S8, 'Reaction rate constants and model parameters'"
_TABLE_S6 = "Talemi 2016 Table S6, 'State variables and their initial conditions'"
_ESTIMATED = ("method column reads 'Estimated': fitted in COPASI to the Hog1PP, Slt2PP, single-cell "
              "volume and intracellular glycerol time series of this paper, no interval reported")

TALEMI_PARAMS = ParamRegistry("mech/cell_wall_talemi2016.py -- osmo-stat, Slt2/CWI arm")


def _estimated(name: str, value: float, units: str, sbml_id: str, row: str) -> Param:
    """One of Table S8's seventeen fitted constants: registered ASSERTED, never MEASURED."""
    return TALEMI_PARAMS.add(Param.asserted(
        name, value, units, f"{_TABLE_S8} row '{row}', SBML id {sbml_id!r}; {_ESTIMATED}"))


def _auxiliary(name: str, value: float, units: str, sbml_id: str, remark: str) -> Param:
    """One of Table S7's declared quantities. ASSERTED unless it carries its own citation."""
    return TALEMI_PARAMS.add(Param.asserted(
        name, value, units, f"{_TABLE_S7} row '{remark}', SBML id {sbml_id!r}"))


# ------------------------------------------------------------------------------------------
# Table S7: the physical constants and the biophysics.
# ------------------------------------------------------------------------------------------

R = _auxiliary("R", 8.314, "J/(mol*K)", "R", "R [J/mol/K], Gas constant")
T_KELVIN = _auxiliary("T", 303.15, "K", "T", "T [K], corresponds to 30 degC")
AVOGADRO = _auxiliary("mol", 6.022e23, "1/mol", "mol", "Mol, Mole number")
C2P = _auxiliary("f_c2p", 1e-9, "MPa/uM", "c2p",
                 "f_c2p, factor converting concentrations to pressures in MPa")
LP = TALEMI_PARAMS.add(Param.asserted(
    "Lp", 0.013, "um/(MPa*s)",
    f"{_TABLE_S7} row 'Lp [um/Mpa/s], Hydraulic conductivity (estimate from data from [10])', "
    "SBML id 'Lp'; reference 10 is Eriksson et al. 2010 Lab Chip 10:617, a microfluidic device "
    "paper, and the row's own word is 'estimate'"))
P0 = TALEMI_PARAMS.add(Param.asserted(
    "P0", 0.61, "MPa",
    f"{_TABLE_S7} row 'P0 [MPa], Initial turgor pressure', cited to reference 7 = Schaber et al. "
    "2010 Eur Biophys J 39:1547 (doi 10.1007/s00249-010-0612-0), SBML id 'P0'; imported by "
    "Talemi from that paper without an interval"))
EPS = TALEMI_PARAMS.add(Param.asserted(
    "eps", 14.3, "MPa",
    f"{_TABLE_S7} row 'eps, Membrane rigidity', cited to Schaber et al. 2010 Eur Biophys J "
    "39:1547, SBML id 'eps'; imported without an interval"))
F_MIN = TALEMI_PARAMS.add(Param.asserted(
    "fmin", 0.41, "dimensionless",
    f"{_TABLE_S7} row 'fmin, Minimal cell volume (as fraction of total)', cited to Schaber et al. "
    "2010 Eur Biophys J 39:1547, SBML id 'minf'"))
F_NUC = _auxiliary("fn", 0.8, "dimensionless", "maxHog1nucf",
                   "fn, Fraction of activated Hog1 molecule in the nucleus upon maximal activation")
V0 = _auxiliary("V0", 50.0, "fL", "V0", "V0 [fL], Initial total cell volume")
MEDIUM_VOLUME_FACTOR = _auxiliary("medium_volume_factor", 1000.0, "dimensionless", "Vex",
                                  "Vmedium [fL] = 1000*V0, External volume")
CE0 = _auxiliary("ce0", 260000.0, "uM", "ce0",
                 "ce0 [uM], Initial osmolarity of the medium, cited to Schaber et al. 2010")
T_OFF = _auxiliary("toff", 120.0, "s", "toff", "toff [s], Time before first osmotic stress")
T_S = _auxiliary("ts_published", 840.0, "s", "ts",
                 "ts [s], Time between two consecutive osmotic stresses")
T_MIX = _auxiliary("tm", 10.0, "s", "tm", "tm [s], Mixing time of sorbitol in the medium")
S1_PUBLISHED = _auxiliary("S1", 800000.0, "uM", "s1",
                          "S1 [uM], Sorbitol concentration for first osmotic stress")
S2_PUBLISHED = _auxiliary("S2", 270000.0, "uM", "s2",
                          "S2 [uM], Sorbitol concentration for second osmotic stress")

# ------------------------------------------------------------------------------------------
# Table S6 and S7: the abundances and the two data-derived baselines.
# ------------------------------------------------------------------------------------------

HOG1_MOLECULES = TALEMI_PARAMS.add(Param.asserted(
    "hog1_molecules", 6788.0, "molecules/cell",
    f"{_TABLE_S7} row 'Hog1t [uM] = 0.3821', remark '6788 * f_N2uM: molecule numbers from "
    "http://yeastgfp.yeastgenome.org/'. A URL is not a citation this repository can resolve, "
    "there is no interval, and no paper is named, so this is ASSERTED and not MEASURED",
    missing="a named S. cerevisiae abundance measurement for Hog1 with an interval, in the "
            "working strain and growth condition"))
SLT2_MOLECULES = TALEMI_PARAMS.add(Param.asserted(
    "slt2_molecules", 3230.0, "molecules/cell",
    f"{_TABLE_S6} Slt2 and Slt2PP rows, '3230 * f_N2uM * ...'. UNLIKE Hog1t and Fps1t this "
    "number is NOT in Table S7 and carries no source at all -- not even the yeastGFP URL. It is "
    "the scale of the entire arm this module transcribes and it is the one abundance the source "
    "does not source",
    missing="a named S. cerevisiae Slt2/Mpk1 abundance with an interval. Its absence is why no "
            "absolute Slt2 concentration may be claimed from this transcription"))
FPS1_MOLECULES = TALEMI_PARAMS.add(Param.asserted(
    "fps1_molecules", 907.0, "molecules/cell",
    f"{_TABLE_S7} row 'Fps1t [uM] = 0.051', remark '907 f_N2uM: molecule numbers from "
    "http://yeastgfp.yeastgenome.org/'. Same grade and same reason as hog1_molecules",
    missing="a named S. cerevisiae Fps1 abundance with an interval"))
HOG1PP_REL_INI = TALEMI_PARAMS.add(Param.asserted(
    "hog1pp_rel_ini", 0.061, "dimensionless",
    f"{_TABLE_S6} Hog1PP row: 'It was derived from data that 6.1 % of the maximal "
    "phosphorylation is the steady state activation', SBML id 'Hog1PPrelIni'. Derived from THIS "
    "paper's own western blots, which is why it is also a fitted target and not an independent "
    "one"))
SLT2P_REL_INI = TALEMI_PARAMS.add(Param.asserted(
    "slt2p_rel_ini", 0.246, "dimensionless",
    f"{_TABLE_S6} Slt2PP row: 'It was derived from data that 24.6 % of the maximal "
    "phosphorylation is the steady state activation', SBML id 'Slt2PrelIni'. Same provenance and "
    "same limitation as hog1pp_rel_ini"))
GLYIN_INITIAL = TALEMI_PARAMS.add(Param.asserted(
    "glyin_initial", 180000.008668013, "uM",
    f"{_TABLE_S6} Glyin row, SBML id 'Metabolite_1'. The table's VALUE column prints 180000 and "
    "its own remark derives '1/18/29.5*1e8' from 0.1 mM/OD, 18e6 cells/mL and 29.5 fL, which is "
    "188324 -- a 4.6% disagreement between a table cell and the derivation beside it. The "
    "deposit uses the printed 180000, and so does this transcription"))
GLYEX_RATIO = TALEMI_PARAMS.add(Param.asserted(
    "glyex_over_glyin_initial", 0.01, "dimensionless",
    f"{_TABLE_S6} Glyex row states 'assumed to be 1000 times lower than Glyin'. The DEPOSIT is "
    "100 times lower: Glyex 9.00000043340065e7/50000 = 1800.00008668013 uM against Glyin "
    "180000.008668013 uM. The deposited ratio is used, because it is the file that runs, and the "
    "printed ratio is recorded here rather than silently reconciled"))

# ------------------------------------------------------------------------------------------
# Table S8: the seventeen estimated rate constants, plus the three switches.
# ------------------------------------------------------------------------------------------

K0 = _estimated("k0", 59.901, "uM/s", "k0", "k0, Hog1Signal production rate constant")
K1 = _estimated("k1", 25.4711, "1/s", "k1", "k1, Hog1Signal degradation rate constant")
VMAX2 = _estimated("vmax2", 2.15338, "1/(fL*s)", "k2",
                   "vmax2, Volume mediated Hog1Signal degradation vmax")
KM2 = _estimated("km2", 0.000939165, "uM", "Km2",
                 "km2, Volume mediated Hog1Signal degradation Michaelis constant")
K3 = _estimated("k3", 0.231769, "1/(uM*s)", "v3_k", "k3, Hog1 phosphorylation rate constant")
K5 = _estimated("k5", 0.12348, "1/(uM*s)", "v5_k", "k5, Hog1PP mediated Fps1 closure rate constant")
K6 = _estimated("k6", 1.1077e-06, "1/(uM*s)", "v6_k",
                "k6, Slt2PP mediated Fps1 dephosphorylation (opening)")
VMAX7 = _estimated("vmax7", 849.986, "uM/s", "v7_k", "vmax7, Hog1PP mediated glycerol production vmax")
K8 = _estimated("k8", 0.000776772, "fL/(um^2*s)", "v8_ks", "k8, Fps1 facilitated glycerol diffusion constant")
K9 = _estimated("k9", 0.0772144, "uM/(fL*s)", "k9",
                "k9, Volume mediated Slt2Signal production rate constant")
K10 = _estimated("k10", 0.000131323, "1/s", "k10", "k10, Slt2Signal degradation rate constant")
VMAX11 = _estimated("vmax11", 0.713377, "uM/s", "k11",
                    "vmax11, Sensitizer/Slt2PP mediated Slt2Signal degradation vmax")
KM11 = _estimated("km11", 0.0180575, "uM", "Km11",
                  "km11, Sensitizer/Slt2PP mediated Slt2Signal degradation km")
K12 = _estimated("k12", 0.00929813, "1/(uM*s)", "v12_k",
                 "k12, Slt2Signal mediated Slt2 phosphorylation rate constant")
K14 = _estimated("k14", 0.113591, "1/s", "v14_k",
                 "k14, Slt2PP mediated sensitizer production rate constant")
KI16 = _estimated("ki16", 3.08897, "uM", "v16_Ki", "Ki16, Sensitizer auto inhibitory feedback constant")
N16 = _estimated("n16", 999.898, "dimensionless", "v16_n", "n16, Sensitizer auto inhibitory feedback power")
SENSITIZER_INITIAL = _estimated("sensitizer_initial", 3.42838, "uM", "Sensitizer_0",
                                "Sensitizer initial concentration, estimated by the model (Table S6)")

VMAX16 = TALEMI_PARAMS.add(Param.asserted(
    "vmax16", 3.26314, "uM/s",
    "the DEPOSIT's SBML id 'v16_k'. It contradicts the paper twice: Table S8 prints 'vmax16 = 1, "
    "set' and the Results say 'except parameter k16 which we set to one'. Every other numeric "
    "Table S8 row equals the deposit exactly, so this is the one disagreement, and the deposited "
    "value is used because it is the file that runs. See TABLE_S8_VMAX16 for the published "
    "reading, which stays runnable",
    missing="nothing measurable -- this is a bookkeeping disagreement inside one paper, and only "
            "the authors can say which number made the figures"))

TABLE_S8_VMAX16 = 1.0
"""The published reading of vmax16, kept runnable so the disagreement can be priced.

Under the published protocol it moves the relative Slt2 phosphorylation peak from 89.152% to
89.801% at a 14 min dilution and from 38.948% to 39.150% at 4 min, and not at all at 30 min --
because v16 is a floor on the Sensitizer rather than a rate, and both readings hold the same
floor. :func:`reproduction_report` computes all six numbers.
"""

DEPOSIT_SWITCHES = MappingProxyType({
    "SW_Sensitizer": 1.0,
    "SW_HI": 0.0,
    "SW_SI": 0.0,
})
"""The deposit's three model-selection switches, which encode the paper's Table 1 result.

``SW_Sensitizer = 1`` selects Table S5's ``V11-b`` (sensitized) over ``V11-a``, and turns on
v14/v15/v16. ``SW_HI = 0`` and ``SW_SI = 0`` select ``V12-a`` and ``V3-a``, i.e. **no crosstalk
in either direction** between Hog1 and Slt2 -- which is the article's own conclusion: "The best
ranked model supports the absence of direct crosstalk between Hog1 and Slt2". They are exponents
on ``(1 + Ki*Inh^n)`` in the deposit's ``_2p_Mod_MA_Inh_KO``, so a zero makes the inhibition
factor exactly 1 rather than approximately so.
"""

CALCULATED_ROWS = MappingProxyType({
    "k4": "k3*Hog1Signal_0*Hog1_0 / (Hog1PP_0 * (1+ki3*Slt2PP_0^n3)^SW_SI)",
    "km7": "vmax7*Hog1PP_0*Vos_0 / (Fps1Open_0*k8*A*(Glyin_0-Glyex_0)) - Hog1PP_0",
    "k6b": "(k5*Hog1PP_0*Fps1_0 - k6*Slt2PP_0*Fps1P_0) / Fps1P_0",
    "k13": "k12*Slt2Signal_0*Slt2_0 / (Slt2PP_0 * (1+ki12*Hog1PP_0^n12)^SW_HI)",
    "k15": "(k14*Slt2PP_0 + v16(Sensitizer_0)) / Sensitizer_0",
})
"""Table S8's five rows whose Value column prints a formula rather than a number.

They are NOT registered as parameters, deliberately. Each is fully determined by the registered
constants plus the source's own steady-state assumption, so registering them would count five
degrees of freedom that do not exist; and grading them DERIVED would launder arithmetic on
seventeen ASSERTED rows into a PINNED tag, which is exactly the mistake ``mech/params.py``'s
docstring records. :func:`calculated_constants` computes them and :func:`audit_deposit` checks
all five against the deposit's stored values -- four at a worst relative 1.9e-15, and k13 at
3.0e-12 because it is the only one that reads Slt2Signal|0 and inherits that root's own error.

Table S8's Method column has k6 and k6b's cells crossed: k6 is marked "Calculated" and prints a
number, k6b is marked "Estimated" and prints a formula. The formula is right -- k6b recomputes
from the Fps1 steady state to 1.9e-15 -- so the two Method cells are transposed.
"""

DEPOSIT_SHIPS_UNSTRESSED = (
    "MODEL1606100000 ships the protocol SWITCHED OFF: s1 = s2 = 0 and ts = 800000 s, so no "
    "sorbitol is ever added and the dilution is pushed past any run. Table S7's own defaults are "
    "S1 = 800000 uM, S2 = 270000 uM and ts = 840 s -- i.e. the 14 minute experiment -- and toff = "
    "120 s and tm = 10 s agree with the deposit. This is not a defect; it is why the deposited "
    "initial state is a rest point at all, and it is why Shock.deposited() and Shock.published() "
    "are two different constructors rather than one default")
"""What the deposit ships, against what Table S7 publishes. The difference is the whole protocol."""

SOURCE_DEFECTS = MappingProxyType({
    "biomodels_accession_404": (
        "The Supplementary 'Simulation Instructions' cites BioModels MODEL1604100004, which "
        "answers HTTP 404. The model is MODEL1606100000 and the article names it nowhere. See "
        "DEPOSIT_IDENTITY for how identity was established instead"),
    "vmax16_published_as_one": (
        "Table S8 prints vmax16 = 1, method 'set', and the Results repeat it ('parameter k16 "
        "which we set to one'). The deposit writes 3.26314. 17 of Table S8's 18 numeric rows "
        "match the deposit exactly; this is the eighteenth"),
    "v16_is_a_floor_not_a_rate": (
        "n16 = 999.898 makes v16 a step at ki16 = 3.08897, so it pins the Sensitizer's minimum "
        "rather than setting a rate. Deleting the reaction moves the floor 3.1109 -> 2.4856 and "
        "the 14' Slt2 peak 89.2% -> 115.7%, while changing its constant from 3.26314 to the "
        "published 1 moves the peak by 0.65 percentage points. The reaction is load-bearing and "
        "its own constant is nearly inert, which is presumably why the disagreement survived"),
    "external_glycerol_has_no_sink": (
        "Glyex is produced by v8 and consumed by nothing, and v8 is 4542.6 amount/s at the "
        "deposited state, so the unstressed model has no fixed point. Over 24 h Glyex/Vex goes "
        "1800 -> 9646 uM against ce0 = 260000; the other eleven states settle, Vos moving 1.5e-5 "
        "relative. Table S6 says the external pool 'acts basically as a sink'"),
    "vp_rule_is_the_identity": (
        "Table S7 defines VP0 = V0*exp(-P0/eps) once. The deposit adds VP = "
        "V0*(V0/VP0)^(-ModelValue_7/eps) with ModelValue_7 and eps two copies of 14.3, so the "
        "exponent is -1 and the rule is the identity on VP0. audit_deposit checks the two copies "
        "are still equal"),
    "k6b_method_column_transposed": (
        "Table S8's Method column has k6 and k6b's cells crossed: k6 is marked 'Calculated' and "
        "prints a number, k6b is marked 'Estimated' and prints a formula. The formula is the "
        "right one -- k6b recomputes from the Fps1 steady state to 1.9e-15 -- so the two Method "
        "cells are transposed and k6 is the estimated row"),
    "glyin_table_cell_disagrees_with_its_own_derivation": (
        "Table S6's Glyin value column prints 180000 uM and its remark derives 1/18/29.5*1e8 = "
        "188324 uM from the same three numbers, a 4.6% disagreement. The deposit uses 180000"),
    "glyex_ratio_disagrees_with_the_deposit": (
        "Table S6 says external glycerol is 'assumed to be 1000 times lower than Glyin'. The "
        "deposit makes it 100 times lower (1800.00008668013 against 180000.008668013 uM)"),
    "naming_drift": (
        "The article calls the CWI driver CWISignal; every table and the deposit call it "
        "Slt2Signal. Tables S4-S8 write Slt2PP for the species the deposit names Slt2P. Table "
        "S6's header says concentrations are in 'umoL/fL' and Table S8's says 'umol/L'; f_N2uM = "
        "1e21/(mol*Vos0) settles it as uM"),
})
"""Eight things wrong in the source, each found by transcribing it. None is repaired here."""

# ------------------------------------------------------------------------------------------
# What this axis cannot reach, refused by name.
# ------------------------------------------------------------------------------------------

NOT_DRIVEN = ParamRegistry("mech/cell_wall_talemi2016.py::no-engine-channel")
"""The panel's declared input and output, the medium species, and the volume bridge.

Every row here is REFUSED, so ``float()`` on it raises with the reason and the missing
measurement rather than defaulting. A later agent wiring this axis cites these by name; adopting
the transcription retires none of them.
"""

NOT_DRIVEN.add(Param.refused(
    "promoter_response.rlm1_box", "dimensionless promoter occupancy per uM Slt2PP",
    "generator/stress_panel.py:247-249 declares the cell_wall module as 'Slt2/Mpk1 -> Rlm1' "
    "against an 'RLM1 box', citing Jung & Levin 1999 PMID 10594829; the transcribed source is "
    "Talemi 2016 MODEL1606100000",
    reason="The deposit ends at phosphorylated Slt2 as a KINASE and has no output arm of any "
           "kind. grep -ci over the pinned SBML returns 0 for Rlm1, RLM1 and 'transcription', "
           "and 0 over the article body and the supplementary text as well -- the strings do not "
           "occur anywhere in this source. Jung & Levin 1999 establishes that Rlm1 carries most "
           "CWI transcription; it reports no coefficient and is not a parameterised model. So "
           "the panel's declared reporter has no gain anywhere in this axis's evidence",
    missing="a published parameterised RLM1-box transcription model, or an Slt2 dose response on "
            "an RLM1-box reporter in the working strain, with an interval. Note this would also "
            "need the sign settled: the only Slt2 activation this source models is HYPO-osmotic"))

NOT_DRIVEN.add(Param.refused(
    "cell_wall_damage_dose", "ug/mL congo red or mM caffeine",
    "generator/stress_panel.py STRESSORS['congo_red'] (50 ug/mL) and STRESSORS['caffeine'] "
    "(10 mM), both routed to the cell_wall module; refused from the source side here",
    reason="Talemi's model has NO cell-wall-damage input. Its single stimulus is external "
           "osmolarity, entering through cen as a sorbitol concentration, and its single driver "
           "of Slt2 is cell volume. Congo red, caffeine, calcofluor white, caspofungin, "
           "zymolyase and SDS occur zero times in the deposit. Wiring a damaging agent to this "
           "axis would also invert its meaning: v9 = k9*Vos activates Slt2 when the cell SWELLS, "
           "and the source's CWI arm is a hypo-osmotic adaptation arm, not a damage sensor",
    missing="a published parameterised S. cerevisiae model in which a cell-wall-damaging agent "
            "activates Slt2 -- AND the transport law plus elemental counts a new EXTRACELLULAR "
            "row needs. Talemi's own calcofluor module is the nearest candidate and is refused "
            "separately below"))

NOT_DRIVEN.add(Param.refused(
    "medium.cell_wall_damaging_agent", "uM",
    "the contract side of cell_wall_damage_dose: mech/contracts.py EXTRACELLULAR is exactly "
    "eight species (glucose, ethanol, nitrogen, oxygen, glycerol, peroxide, acetate, osmolyte) "
    "and none is a cell-wall-damaging agent",
    reason="Adding a species to EXTRACELLULAR is one line of Python and a whole experiment: "
           "every row carries a carbon and a nitrogen count, joins the conservation audit, and "
           "needs a transport law before any flux may read it. Congo red and calcofluor white "
           "are not metabolised and not transported; they bind chitin and glucan at the wall, "
           "which is a compartment this contract does not have at all. The gap is on both sides "
           "at once, exactly as it is for calcium",
    missing="a cell-wall compartment or a wall-binding term, the elemental counts and transport "
            "law a new EXTRACELLULAR row needs, and a source whose Slt2 activation reads the "
            "dose. None of the three exists here"))

NOT_DRIVEN.add(Param.refused(
    "calcofluor_activating_module", "uM calcofluor white and 8 rate constants",
    "Talemi 2016 Supplementary Information, 'Calcofluor mediated Slt2 activating module' "
    "(vendored as srep30950-s1-p06-07-calcofluor-module.txt) and Figure 6 of the article",
    reason="This is the source's OWN cell-wall-damage input and the deposit does not contain it. "
           "The paper states plainly that 'No model inside the models ensemble was designed such "
           "that can respond to the presence of the calcofluor in the medium', so it bolted on a "
           "separate five-species module (Calcofluor, CALSignal, Degrader, Slt2, Slt2P) whose "
           "parameters were estimated from ONE experiment -- 0.11 uM calcofluor white two hours "
           "after a 0.8 M sorbitol shock. Those parameters are in no table: Table S8 ends at "
           "n16 and the strings k17, k18, v17, v18 occur nowhere. They appear only in the "
           "supplementary text, whose own Figure S6 caption numbers two different reactions "
           "'v4' ('Phosphorylated Slt2 activates the species Degarder, v4, which induces the "
           "CALSignal decay, v4'), and the deposit carries no Calcofluor, CALSignal or Degrader "
           "species to arbitrate against. Assigning its eight constants to its five rate laws "
           "would be inventing kinetics",
    missing="the calcofluor module as a deposited executable, or an unambiguous printed "
            "assignment of its eight constants to its five rate laws. Note that even with it, "
            "the panel's stressors are congo red and caffeine and not calcofluor white"))

NOT_DRIVEN.add(Param.refused(
    "osmotically_active_cell_volume_bridge", "fL per unit engine cell_volume_ratio",
    "mech/engine.py:203 carries 'cell_volume_ratio', Variable('dimensionless', "
    "Compartment.CELL_WATER, 'relative osmotically accessible volume'), driven at "
    "mech/engine.py:499 by mech/signalling.py's NativeHogResponse vVos reaction; Talemi's v9 = "
    "k9*Vos needs Vos in fL",
    reason="This is the one axis in this engine whose driver is ALREADY a dimensionless state, "
           "and that is exactly what makes the bridge tempting and wrong. The engine's ratio is "
           "produced by a DIFFERENT published osmotic model -- Petelenz-Kurdziel et al. 2013 "
           "doi:10.1371/journal.pcbi.1003084, deposit PetelenzKuehn_osmoadaptation_WT -- whose "
           "reference cellvol is 2.40e-4 in its own units against Talemi's 29.5 fL. Writing Vos "
           "= ratio * 29.5 fL asserts that the two models' 'osmotically active volume' is the "
           "same quantity at their respective references, which nothing here shows and which "
           "their different turgor laws (Talemi eps*ln(V/VP0); Petelenz-Kurdziel a piecewise "
           "linear -2.9 term) give reason to doubt. Talemi ALSO carries its own Vos ODE, its own "
           "Turgor, its own Fps1 and its own glycerol balance, so adopting it whole would run "
           "two published osmotic models side by side and adopting only its Slt2 arm would drive "
           "it from a volume its own equations did not produce. That choice is a wiring decision "
           "with a cost either way and it is NOT made here",
    missing="either a demonstration that the two deposits' osmotically active volumes agree at "
            "their references -- same definition, same solid fraction, same measurement -- or a "
            "decision, taken with the cost stated, to replace the engine's osmotic block with "
            "Talemi's entire one rather than to bridge them"))

NOT_DRIVEN.add(Param.refused(
    "external_osmolarity_basal", "uM",
    "Talemi Table S7 row 'ce0 [uM] = 260000, Initial osmolarity of the medium'; the engine's "
    "EXTRACELLULAR does carry an 'osmolyte' species, so this is the dose side and not the "
    "species side",
    reason="Talemi's cen is an ADDED sorbitol concentration on top of a declared 260000 uM total "
           "medium osmolarity. The engine's EXTRACELLULAR 'osmolyte' is one deliberately "
           "ion-agnostic species with no declared total osmolarity behind it and no ionic, "
           "sugar or amino-acid contributions enumerated, so there is no number to add the dose "
           "TO. Reading engine osmolyte as the whole medium osmolarity would assert that the "
           "rest of the medium contributes nothing osmotically, which is false for any real "
           "medium and would put the turgor reference in the wrong place -- ci0 = ce0 + "
           "P0/(f_c2p*R*T) is where ce0 enters, so an error here moves the resting turgor",
    missing="a declared total osmolarity for this repository's media, with the contribution of "
            "each component, so a dose keyed into Control.feed_mM['osmolyte'] can be added to a "
            "basal rather than replacing it"))

PANEL_CHANNEL = MappingProxyType({
    "module": "cell_wall",
    "axis": AXIS,
    "answered": MappingProxyType({}),
    "forfeited": MappingProxyType({"congo_red": 1.0, "caffeine": 1.0}),
    "note": (
        "generator/stress_panel.py routes congo_red and caffeine to the cell_wall module and "
        "declares its output as an RLM1 box. This build answers NEITHER stressor and supplies "
        "NO promoter output. What it supplies is twelve states and a published, reproduced "
        "volume-to-Slt2 arm whose stimulus is osmotic and whose sign is the opposite of a "
        "damage sensor's. It is registered under its own axis name for that reason"),
})
"""The honest mapping from this axis to `generator/stress_panel.py`: zero panel coverage."""

STATE_VARIABLES = (
    ("Vos", "fL", Compartment.CELL_WATER, "osmotically active cell volume, Table S4 rate rule"),
    ("HOGSignal", "uM", Compartment.REGULATORY, "Table S4 rate rule; the article calls it HOGSignal"),
    ("Slt2Signal", "uM", Compartment.REGULATORY,
     "Table S4 rate rule, d/dt = v9 - v10 - v11; the article calls it CWISignal"),
    ("Hog1", "umol", Compartment.CELL_WATER, "unphosphorylated Hog1, Table S4"),
    ("Hog1PP", "umol", Compartment.CELL_WATER, "double-phosphorylated Hog1, Table S4"),
    ("Slt2", "umol", Compartment.CELL_WATER, "unphosphorylated Slt2, Table S4"),
    ("Slt2P", "umol", Compartment.CELL_WATER,
     "phosphorylated Slt2; Tables S4-S8 write Slt2PP for this species"),
    ("Glyin", "umol", Compartment.CELL_WATER, "intracellular glycerol, Table S4"),
    ("Glyex", "umol", Compartment.MEDIUM, "extracellular glycerol, Table S4; a pure sink"),
    ("Fps1", "umol", Compartment.REGULATORY, "open Fps1 aquaglyceroporin, Table S4"),
    ("Fps1P", "umol", Compartment.REGULATORY, "closed (phosphorylated) Fps1, Table S4"),
    ("Sensitizer", "umol", Compartment.REGULATORY,
     "the hypothetical entity modulating Slt2Signal degradation, Table S4 dagger row"),
)
"""The twelve states this axis adds, in `mech/contracts.py`'s vocabulary.

Two mappings are declared rather than obvious. The nine species are AMOUNTS, not concentrations:
the deposit sets ``hasOnlySubstanceUnits="true"`` on every one and divides by a compartment
explicitly inside each rate law, so anything wiring this beside ``calcium_ke2013`` -- whose four
states are concentrations -- must not mix them. And the source's ``Membrane`` compartment has no
counterpart in :class:`~ystwin.mech.contracts.Compartment`; Fps1, Fps1P and Sensitizer are put
under ``REGULATORY`` because they gate or modulate a flux rather than being cell-water pools, and
the deposit gives Membrane size 1 so no volume scaling is lost by the choice.
"""

STATE_NAMES = tuple(row[0] for row in STATE_VARIABLES)
"""The integration order used by :meth:`OsmoStat.rhs` and :func:`initial_state`."""

PUBLISHED = MappingProxyType({
    "slt2p_relative_baseline_percent": 24.6,
    "hog1pp_relative_baseline_percent": 6.1,
    "hog1_total_uM": 0.3821,
    "fps1_total_uM": 0.051,
    "vos0_fL": 29.5,
    "turgor0_MPa": 0.61,
    "s1_uM": 800000.0,
    "s2_uM": 270000.0,
    "t_off_s": 120.0,
    "t_s_s": 840.0,
    "vmax16_table_s8": 1.0,
})
"""Every source number this transcription is scored against, with where it is printed.

The first two are Table S6's data-derived baselines, the next two Table S7's abundance rows
rounded as the table prints them, and the rest Table S7's geometry and protocol. All of them
were fitted to, or asserted by, this same paper, which is why every :class:`Target` registered
below carries ``fitted=True`` and the gate has nothing in its denominator.
"""

for _target in (
    Target("slt2p_baseline", "relative Slt2 phosphorylation at rest",
           "anti-phospho-p44/42 western blot, normalised to the maximal signal", 0.05,
           "fraction of maximal",
           f"{_TABLE_S6}: '24.6 % of the maximal phosphorylation is the steady state "
           "activation'. The model's Slt2/Slt2PP split was SET from this number, so it cannot "
           "also score the model", fitted=True),
    Target("hog1pp_baseline", "relative Hog1 phosphorylation at rest",
           "anti-phospho-p38 western blot, normalised to the maximal signal", 0.05,
           "fraction of maximal",
           f"{_TABLE_S6}: '6.1 % of the maximal phosphorylation is the steady state "
           "activation'. Set, not predicted", fitted=True),
    Target("slt2_peak_4min", "the 4 min hyper-hypo Slt2 activation peak (Table 1's '4MiP')",
           "anti-phospho western blot time series over a 0.8 M -> 0.27 M sorbitol dilution", 0.05,
           "fraction of maximal",
           "Talemi 2016 Table 1 and Figs 3-4; this is one of the series the model was estimated "
           "against, and the sensitizer was introduced because the ensemble could not fit it",
           fitted=True),
):
    TALEMI_PARAMS.add_target(_target)
del _target


# ------------------------------------------------------------------------------------------
# The transcription.
# ------------------------------------------------------------------------------------------

def _digest(relative: str) -> str:
    path = paths.data_dir().parent / relative
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _derived_geometry() -> dict:
    """Table S7's auxiliary quantities, every one evaluated from a registered constant."""
    v0 = float(V0)
    vos0 = v0 * (1.0 - float(F_MIN))
    return {
        "Vb": v0 * float(F_MIN),
        "Vos0": vos0,
        "Vex": float(MEDIUM_VOLUME_FACTOR) * v0,
        "f_N2uM": 1e21 / (float(AVOGADRO) * vos0),
        "VP0": v0 * math.exp(-float(P0) / float(EPS)),
        "A0": (36.0 * math.pi) ** (1.0 / 3.0) * v0 ** (2.0 / 3.0),
        "ci0": float(CE0) + float(P0) / (float(C2P) * float(R) * float(T_KELVIN)),
    }


def initial_state() -> dict:
    """Table S6's published closed forms, evaluated. Amounts in umol, volumes in fL.

    Every row is the formula the table prints, not a number read out of the deposit -- which is
    what makes the agreement in :func:`reproduction_report` a check rather than a tautology. The
    two signals are the positive root of the quadratic Table S6 prints for each.

    One reading had to be resolved and the deposit resolves it: Table S6 renders the Slt2 row as
    ``3230 * f_N2uM * 1 - [Slt2PP]0 * fn``, whose grouping the PDF's layout loses. Read as
    ``(1 - [Slt2PP]0) * fn`` it gives 0.10967 uM; read as ``1 - [Slt2PP]0 * fn`` it gives
    0.146037185686382 uM, and the deposit stores 4.30809697774826/29.5 = 0.146037185686382.
    """
    g = _derived_geometry()
    f, fn, vos0 = g["f_N2uM"], float(F_NUC), g["Vos0"]
    hog1_total, slt2_total = float(HOG1_MOLECULES) * f, float(SLT2_MOLECULES) * f
    hog1pp_c = hog1_total * float(HOG1PP_REL_INI) * fn
    slt2p_c = slt2_total * float(SLT2P_REL_INI) * fn
    glyin_c = float(GLYIN_INITIAL)
    return {
        "Vos": vos0,
        "HOGSignal": _signal_root(float(K1), float(K0), float(VMAX2) * vos0, float(KM2)),
        "Slt2Signal": _signal_root(float(K10), float(K9) * vos0,
                                   float(VMAX11) * float(SENSITIZER_INITIAL), float(KM11)),
        "Hog1": (hog1_total - hog1pp_c) * vos0,
        "Hog1PP": hog1pp_c * vos0,
        "Slt2": (slt2_total - slt2p_c) * vos0,
        "Slt2P": slt2p_c * vos0,
        "Glyin": glyin_c * vos0,
        "Glyex": glyin_c * float(GLYEX_RATIO) * g["Vex"],
        "Fps1": float(FPS1_MOLECULES) * f / 2.0,
        "Fps1P": float(FPS1_MOLECULES) * f / 2.0,
        "Sensitizer": float(SENSITIZER_INITIAL),
    }


def _signal_root(linear: float, production: float, vmax: float, km: float) -> float:
    """The positive root of ``production = linear*S + vmax*S/(km+S)``, as Table S6 prints it.

    Table S6 gives Hog1Signal|0 and Slt2Signal|0 as an explicit quadratic root, and both are the
    same shape: a constant production balanced by one first-order and one saturating loss. The
    PDF renders both denominators as ``2*k`` with the subscript lost; only ``2*linear`` closes
    the quadratic, and it reproduces both deposited values.
    """
    b = linear * km + vmax - production
    return (-b + math.sqrt(b * b + 4.0 * linear * production * km)) / (2.0 * linear)


def calculated_constants(*, sensitizer: bool = True, vmax16: float | None = None) -> dict:
    """Table S8's five "Calculated" rows, evaluated from the registered constants.

    Each is the source's own steady-state assumption on a different rate law, so recomputing all
    five exercises the Hog1, Fps1, glycerol and Slt2 blocks at once. :data:`CALCULATED_ROWS`
    carries the printed formulas and says why these are not registered as parameters.
    """
    g = _derived_geometry()
    y0 = initial_state()
    vos0, area = g["Vos0"], g["A0"]
    hog1_c, hog1pp_c = y0["Hog1"] / vos0, y0["Hog1PP"] / vos0
    slt2_c, slt2p_c = y0["Slt2"] / vos0, y0["Slt2P"] / vos0
    fps1, fps1p = y0["Fps1"], y0["Fps1P"]
    glyin_c, glyex_c = y0["Glyin"] / vos0, y0["Glyex"] / g["Vex"]
    fps1_open = fps1 / (fps1 + fps1p)
    sw_hi, sw_si = DEPOSIT_SWITCHES["SW_HI"], DEPOSIT_SWITCHES["SW_SI"]
    return {
        "k4": (float(K3) * y0["HOGSignal"] * hog1_c
               / (hog1pp_c * (1.0 + 0.1 * slt2p_c ** 0.1) ** sw_si)),
        "km7": (vos0 * float(VMAX7) * hog1pp_c
                / (fps1_open * float(K8) * area * (glyin_c - glyex_c)) - hog1pp_c),
        "k6b": (float(K5) * hog1pp_c * fps1 - float(K6) * slt2p_c * fps1p) / fps1p,
        "k13": (float(K12) * y0["Slt2Signal"] * slt2_c
                / (slt2p_c * (1.0 + 0.1 * hog1pp_c ** 0.1) ** sw_hi)),
        "k15": ((float(K14) * slt2p_c + _v16(y0["Sensitizer"], sensitizer, vmax16))
                / y0["Sensitizer"] if sensitizer else 0.0),
    }


def _v16(sensitizer_amount: float, on: bool, vmax16: float | None) -> float:
    """``SW * vmax16 / (1 + (Sensitizer/ki16)^n16)``, evaluated in log space.

    Identical arithmetic, not a changed law: with n16 = 999.898 the direct power overflows a
    Python float for Sensitizer above 6.28 and raises OverflowError, where the mathematical limit
    of the expression is exactly 0. Nothing is clipped and no tolerance is introduced.
    """
    if not on:
        return 0.0
    k = float(VMAX16) if vmax16 is None else float(vmax16)
    if sensitizer_amount <= 0.0:
        return k
    exponent = float(N16) * math.log(sensitizer_amount / float(KI16))
    return 0.0 if exponent > 709.0 else k / (1.0 + math.exp(exponent))


@dataclass(frozen=True)
class Shock:
    """Table S7's ``cen`` row: an osmotic sorbitol shock, in the source's own seconds.

    The deposit ships ``s1 = s2 = 0`` and ``ts = 800000`` -- the shock switched off and the
    dilution pushed past any run -- so the file as deposited is the unstressed state.
    :meth:`published` restores Table S7's own five numbers, which are the 14 min experiment.

    Times stay in seconds because they are transcribed published values; converting Table S7's
    ``toff = 120`` into hours would make the audit against that row unreadable.
    ``OsmoStat.timecourse`` takes hours and converts once, at the boundary.
    """

    s1_uM: float = 0.0
    s2_uM: float = 0.0
    t_off_s: float = float(T_OFF)
    t_s_s: float = 800000.0

    @classmethod
    def deposited(cls) -> Shock:
        """Exactly what MODEL1606100000 ships: no shock, no dilution."""
        return cls(s1_uM=0.0, s2_uM=0.0, t_off_s=float(T_OFF), t_s_s=800000.0)

    @classmethod
    def published(cls, *, dilute_after_s: float | None = None) -> Shock:
        """Table S7's S1, S2, toff and ts. ``dilute_after_s`` overrides ts only."""
        return cls(s1_uM=float(S1_PUBLISHED), s2_uM=float(S2_PUBLISHED), t_off_s=float(T_OFF),
                   t_s_s=float(T_S) if dilute_after_s is None else float(dilute_after_s))

    def added_osmolarity_uM(self, time_s: float) -> float:
        """``cen(t)``, transcribed from Table S7 and from the deposit's assignment rule."""
        if time_s < self.t_off_s:
            return 0.0
        if time_s < self.t_s_s + self.t_off_s:
            return self.s1_uM * (1.0 - math.exp(-(time_s - self.t_off_s) / float(T_MIX)))
        decay = math.exp(-((time_s - self.t_s_s) - self.t_off_s) / float(T_MIX))
        return (self.s1_uM - self.s2_uM) * decay + self.s2_uM


@dataclass(frozen=True)
class OsmoStat:
    """MODEL1606100000's twelve-state system, with the deposit's three switches exposed.

    Args:
        sensitizer: ``SW_Sensitizer``. True is the deposit and the paper's selected model.
            Setting it False deletes v14, v15 and v16 and collapses v11 to Table S5's ``V11-a``.
            It does NOT give the paper's non-sensitized models, which were re-estimated -- see
            :func:`sensitizer_switch_comparison`, which measures how far apart the two are.
        hog1_inhibits_slt2: ``SW_HI``. The deposit sets it to 0, which is the article's own
            model-selection result, not a simplification made here.
        slt2_inhibits_hog1: ``SW_SI``. The deposit sets it to 0, for the same reason.
        vmax16: overrides the deposited 3.26314. Pass :data:`TABLE_S8_VMAX16` for the published
            reading; see :data:`SOURCE_DEFECTS`.

    Construction runs :func:`audit_deposit`, so an instance cannot exist unless the vendored SBML
    still hashes to what was transcribed and every literal below still matches it.
    """

    sensitizer: bool = True
    hog1_inhibits_slt2: bool = False
    slt2_inhibits_hog1: bool = False
    vmax16: float | None = None

    def __post_init__(self) -> None:
        audit_deposit()

    @property
    def geometry(self) -> dict:
        return _derived_geometry()

    def initial_vector(self) -> np.ndarray:
        """:func:`initial_state` in :data:`STATE_NAMES` order.

        The signals are re-rooted when the sensitizer is switched off, using Table S6's OWN
        non-dagger ``Slt2Signal|0`` row -- the source prints both variants and this uses the one
        that matches the switch rather than carrying the sensitized baseline into a model that
        has no sensitizer.
        """
        state = dict(initial_state())
        if not self.sensitizer:
            state["Slt2Signal"] = _signal_root(
                float(K10), float(K9) * self.geometry["Vos0"], float(VMAX11), float(KM11))
        return np.array([state[name] for name in STATE_NAMES], dtype=float)

    def constants(self) -> dict:
        """The five calculated rows, at this instance's switches."""
        base = calculated_constants(sensitizer=self.sensitizer, vmax16=self.vmax16)
        if not self.sensitizer:
            y0 = self.initial_vector()
            idx = {name: i for i, name in enumerate(STATE_NAMES)}
            vos0 = self.geometry["Vos0"]
            base = dict(base)
            base["k13"] = (float(K12) * y0[idx["Slt2Signal"]] * (y0[idx["Slt2"]] / vos0)
                           / (y0[idx["Slt2P"]] / vos0))
        return base

    def fluxes(self, time_s: float, y: np.ndarray, shock: Shock, constants: dict) -> dict:
        """The twelve reaction rates plus dVos/dt, all in amount/s and fL/s.

        Transcribed one for one from the deposit's kinetic laws and its eight function
        definitions; the ``Membrane`` compartment has size 1 and is written explicitly wherever
        the deposit writes it, so the factor is visible rather than cancelled by hand.
        """
        idx = {name: i for i, name in enumerate(STATE_NAMES)}
        g = self.geometry
        vos, hog_signal, slt2_signal = y[idx["Vos"]], y[idx["HOGSignal"]], y[idx["Slt2Signal"]]
        hog1, hog1pp = y[idx["Hog1"]], y[idx["Hog1PP"]]
        slt2, slt2p = y[idx["Slt2"]], y[idx["Slt2P"]]
        glyin, glyex = y[idx["Glyin"]], y[idx["Glyex"]]
        fps1, fps1p, sens = y[idx["Fps1"]], y[idx["Fps1P"]], y[idx["Sensitizer"]]
        membrane, vex = 1.0, g["Vex"]

        total_volume = g["Vb"] + vos
        area = (36.0 * math.pi) ** (1.0 / 3.0) * total_volume ** (2.0 / 3.0)
        turgor = (float(EPS) * math.log(total_volume / g["VP0"])
                  if total_volume > g["VP0"] else 0.0)
        osmin_c = (glyin + (g["ci0"] - float(GLYIN_INITIAL)) * g["Vos0"]) / vos
        osmex_c = (float(CE0) + shock.added_osmolarity_uM(time_s) + glyex / vex
                   - float(GLYIN_INITIAL) * float(GLYEX_RATIO))
        fps1_open = (fps1 / membrane) / ((fps1 / membrane) + (fps1p / membrane))
        sw_hi = 1.0 if self.hog1_inhibits_slt2 else 0.0
        sw_si = 1.0 if self.slt2_inhibits_hog1 else 0.0
        sw_sens = 1.0 if self.sensitizer else 0.0

        v0 = float(K0)
        v1 = float(K1) * hog_signal
        v2 = float(VMAX2) * vos * hog_signal / (float(KM2) + hog_signal)
        v3 = (vos * float(K3) * hog_signal * (hog1 / vos)
              / (1.0 + 0.1 * (slt2p / vos) ** 0.1) ** sw_si)
        v4 = constants["k4"] * hog1pp
        v5 = membrane * float(K5) * (hog1pp / vos) * (fps1 / membrane)
        v6 = membrane * (slt2p / vos) * float(K6) * (fps1p / membrane)
        v6b = constants["k6b"] * fps1p
        v7 = vos * float(VMAX7) * (hog1pp / vos) / (constants["km7"] + (hog1pp / vos))
        v8 = fps1_open * float(K8) * area * ((glyin / vos) - (glyex / vex))
        v9 = float(K9) * vos
        v10 = float(K10) * slt2_signal
        v11 = (float(VMAX11) * ((sens / membrane) ** sw_sens) * slt2_signal
               / (float(KM11) + slt2_signal))
        v12 = (vos * float(K12) * slt2_signal * (slt2 / vos)
               / (1.0 + 0.1 * (hog1pp / vos) ** 0.1) ** sw_hi)
        v13 = constants["k13"] * slt2p
        v14 = membrane * float(K14) * (slt2p / vos) * sw_sens
        v15 = membrane * constants["k15"] * (sens / membrane) * sw_sens
        v16 = membrane * _v16(sens / membrane, self.sensitizer, self.vmax16)
        d_vos = -float(LP) * area * (
            turgor + float(C2P) * float(R) * float(T_KELVIN) * (osmex_c - osmin_c))
        return {"v0": v0, "v1": v1, "v2": v2, "v3": v3, "v4": v4, "v5": v5, "v6": v6,
                "v6b": v6b, "v7": v7, "v8": v8, "v9": v9, "v10": v10, "v11": v11, "v12": v12,
                "v13": v13, "v14": v14, "v15": v15, "v16": v16, "dVos": d_vos,
                "Turgor": turgor, "Area": area, "Fps1Open": fps1_open,
                "Osmin_uM": osmin_c, "Osmex_uM": osmex_c}

    def rhs(self, time_s: float, y: np.ndarray, shock: Shock, constants: dict) -> np.ndarray:
        """Table S4's ODE system, in :data:`STATE_NAMES` order. Seconds in, amount/s out."""
        f = self.fluxes(time_s, y, shock, constants)
        return np.array([
            f["dVos"],
            f["v0"] - f["v1"] - f["v2"],
            f["v9"] - f["v10"] - f["v11"],
            -f["v3"] + f["v4"],
            f["v3"] - f["v4"],
            -f["v12"] + f["v13"],
            f["v12"] - f["v13"],
            f["v7"] - f["v8"],
            f["v8"],
            -f["v5"] + f["v6"] + f["v6b"],
            f["v5"] - f["v6"] - f["v6b"],
            f["v14"] + f["v16"] - f["v15"],
        ], dtype=float)

    def steady_state_residual(self) -> dict:
        """Every state's rate at the deposited initial vector, absolute and relative.

        Reported per state rather than as one worst case, because one of the twelve is genuinely
        not at rest: see :data:`SOURCE_DEFECTS` ``external_glycerol_has_no_sink``.
        """
        y0 = self.initial_vector()
        rates = self.rhs(0.0, y0, Shock.deposited(), self.constants())
        return {name: {"rate_per_s": float(rate),
                       "relative_per_s": float(abs(rate) / abs(value)) if value else float(abs(rate))}
                for name, rate, value in zip(STATE_NAMES, rates, y0)}

    def timecourse(self, *, hours: float, shock: Shock | None = None,
                   points: int = 2001, method: str = PINNED_METHOD) -> dict:
        """Integrate for ``hours``, converting to the source's seconds once, here.

        Returns the state traces plus the three relative readouts the deposit's own Fit*
        assignment rules define, which are what the article's figures plot.
        """
        if hours <= 0:
            raise ValueError(f"hours must be positive, got {hours!r}")
        span_s = float(hours) * SECONDS_PER_HOUR
        protocol = Shock.deposited() if shock is None else shock
        constants = self.constants()
        times = np.linspace(0.0, span_s, int(points))
        last = None
        for attempt in (method, FALLBACK_METHOD):
            solution = solve_ivp(
                lambda t, y: self.rhs(t, y, protocol, constants),
                (0.0, span_s), self.initial_vector(), method=attempt, t_eval=times,
                rtol=1e-10, atol=1e-12, max_step=float(T_MIX) / 2.0)
            last = solution
            if solution.success:
                break
        if last is None or not last.success:
            raise RuntimeError(f"Talemi osmo-stat integration failed: {last.message if last else ''}")
        traces = {name: last.y[i] for i, name in enumerate(STATE_NAMES)}
        return {
            "time_s": last.t,
            "time_h": last.t / SECONDS_PER_HOUR,
            "states": traces,
            "slt2p_relative_percent": self.relative_slt2_phosphorylation(traces["Slt2P"]),
            "hog1pp_relative_percent": self.relative_hog1_phosphorylation(traces["Hog1PP"]),
            "volume_relative": (self.geometry["Vb"] + traces["Vos"]) / float(V0),
            "method": last.method if hasattr(last, "method") else method,
        }

    def relative_slt2_phosphorylation(self, slt2p_amount):
        """The deposit's ``FitSlt2Prel`` rule: percent of maximal phosphorylation.

        ``100 * Slt2P / (fn * (Slt2_0 + Slt2P_0))``. It reads 24.6 at the deposited state, which
        is Table S6's published baseline, and it is NOT bounded by 100 -- a large hypo-shock
        drives it past 120% in this model.
        """
        y0 = initial_state()
        return 100.0 * np.asarray(slt2p_amount) / (float(F_NUC) * (y0["Slt2"] + y0["Slt2P"]))

    def relative_hog1_phosphorylation(self, hog1pp_amount):
        """The deposit's ``FitHog1PPrel`` rule. Reads 6.1 at the deposited state."""
        y0 = initial_state()
        return 100.0 * np.asarray(hog1pp_amount) / (float(F_NUC) * (y0["Hog1"] + y0["Hog1PP"]))


def volume_response(ratios=(0.90, 0.95, 1.00, 1.05, 1.10)) -> dict:
    """The sign, computed. Slt2 rises with volume and the HOG signal falls with it.

    Each row holds ``Vos`` at ``ratio * Vos0`` and solves the two signals' rest points and the
    Slt2/Slt2P split analytically, so the direction is read off the equations rather than off an
    integration whose transient could be mistaken for the answer.
    """
    g = _derived_geometry()
    y0 = initial_state()
    slt2_total_amount = y0["Slt2"] + y0["Slt2P"]
    k13 = calculated_constants()["k13"]
    rows = {}
    for ratio in ratios:
        vos = g["Vos0"] * float(ratio)
        slt2_signal = _signal_root(float(K10), float(K9) * vos,
                                   float(VMAX11) * float(SENSITIZER_INITIAL), float(KM11))
        hog_signal = _signal_root(float(K1), float(K0), float(VMAX2) * vos, float(KM2))
        slt2p = slt2_total_amount * float(K12) * slt2_signal / (float(K12) * slt2_signal + k13)
        rows[float(ratio)] = {
            "v9_uM_per_s": float(K9) * vos,
            "slt2_signal_uM": slt2_signal,
            "hog_signal_uM": hog_signal,
            "slt2p_relative_percent": 100.0 * slt2p / (float(F_NUC) * slt2_total_amount),
        }
    values = [rows[float(r)] for r in ratios]
    return {
        "rows": rows,
        "slt2_monotone_increasing": all(
            a["slt2p_relative_percent"] < b["slt2p_relative_percent"]
            for a, b in zip(values, values[1:])),
        "hog_monotone_decreasing": all(
            a["hog_signal_uM"] > b["hog_signal_uM"] for a, b in zip(values, values[1:])),
        "statement": SIGN_OF_THE_AXIS,
    }


def sensitizer_switch_comparison(*, dilution_minutes=(4.0, 14.0, 30.0)) -> dict:
    """Run the deposit's ``SW_Sensitizer`` both ways, and say what the comparison is NOT.

    It is not the paper's model ranking. Talemi re-estimated every non-sensitized model; flipping
    the switch keeps the sensitized model's constants, and Table S6's own non-dagger
    ``Slt2Signal|0`` row shows how far that moves the baseline. The returned record carries the
    ratio so the difference cannot be read as a like-for-like comparison.
    """
    sensitized, plain = OsmoStat(), OsmoStat(sensitizer=False)
    idx = STATE_NAMES.index("Slt2Signal")
    baseline_on = float(sensitized.initial_vector()[idx])
    baseline_off = float(plain.initial_vector()[idx])
    rows = {}
    for minutes in dilution_minutes:
        shock = Shock.published(dilute_after_s=60.0 * float(minutes))
        hours = (shock.t_off_s + shock.t_s_s + 1800.0) / SECONDS_PER_HOUR
        after = {}
        for label, model in (("sensitized", sensitized), ("no_sensitizer", plain)):
            run = model.timecourse(hours=hours, shock=shock)
            window = run["time_s"] >= shock.t_off_s + shock.t_s_s
            after[label] = float(run["slt2p_relative_percent"][window].max())
        rows[float(minutes)] = after
    return {
        "baseline_slt2_signal_uM": {"sensitized": baseline_on, "no_sensitizer": baseline_off},
        "baseline_ratio": baseline_off / baseline_on,
        "km11_multiples_at_baseline": {"sensitized": baseline_on / float(KM11),
                                       "no_sensitizer": baseline_off / float(KM11)},
        "peak_percent_after_dilution": rows,
        "this_is_not_the_paper_ranking": (
            "Talemi re-estimated the non-sensitized models; this flips one switch at the "
            "sensitized model's own constants. At those constants vmax11 = 0.713 alone cannot "
            "balance v9 = k9*Vos0 = 2.278, so the balance falls to the linear k10 term and the "
            "signal sits far above km11, saturated and unable to respond. The flat peaks below "
            "are that saturation, not the paper's threshold argument"),
    }


def audit_deposit() -> dict:
    """Re-read the vendored SBML and check this transcription against it, or refuse to run.

    Seven checks, in the order a reader debugging a failure wants them: both SBML files' sha256; the species set with every
    initialAmount, boundaryCondition and hasOnlySubstanceUnits; the four compartment sizes; every
    registered parameter literal against the deposit's own value; the three switches; the five
    calculated rows against the deposit's stored values; and the two copies of ``eps`` that make
    the VP rule an identity. It also re-runs the identity comparison against the authors' own
    Supplementary S2 SBML, which is what establishes which deposit this is at all.
    """
    for key in ("sbml_biomodels", "sbml_authors_supplementary_s2", "article_xml",
                "tables_s4_s8_text", "calcofluor_module_text"):
        relative, expected = ARTIFACTS[key]
        actual = _digest(relative)
        if actual != expected:
            raise ValueError(
                f"{relative} has sha256 {actual}, expected {expected}. This transcription was "
                "made from the expected bytes; a changed artifact means it is stale")

    deposit = _read_sbml(ARTIFACTS["sbml_biomodels"][0])
    authors = _read_sbml(ARTIFACTS["sbml_authors_supplementary_s2"][0])

    expected_species = {name: (amount, boundary, True) for name, amount, boundary in (
        ("Glyin", 5.31000025570638e6, False), ("Osmin", 1.48097706012933e7, True),
        ("Hog1", 10.7219289272667, False), ("Hog1PP", 0.550073729657921, False),
        ("Slt2", 4.30809697774826, False), ("Slt2P", 1.05556957821322, False),
        ("Glyex", 9.00000043340065e7, False), ("Osmex", 1.3e10, True),
        ("Fps1", 0.0255278667484759, False), ("Fps1P", 0.0255278667484759, False),
        ("Sensitizer", 3.42838, False))}
    if deposit["species"] != expected_species:
        raise ValueError(
            f"the deposit's species block is not what was transcribed: deposit "
            f"{deposit['species']}, transcription {expected_species}")
    expected_compartments = {"Vos": 29.5, "Vex": 50000.0, "V": 50.0, "Membrane": 1.0}
    if deposit["compartments"] != expected_compartments:
        raise ValueError(f"compartment sizes changed: {deposit['compartments']}")

    literals = {
        "R": R, "T": T_KELVIN, "mol": AVOGADRO, "c2p": C2P, "Lp": LP, "P0": P0, "eps": EPS,
        "minf": F_MIN, "maxHog1nucf": F_NUC, "V0": V0, "ce0": CE0, "toff": T_OFF,
        "tm": T_MIX, "Hog1PPrelIni": HOG1PP_REL_INI,
        "Slt2PrelIni": SLT2P_REL_INI, "k0": K0, "k1": K1, "k2": VMAX2, "Km2": KM2,
        "v3_k": K3, "v5_k": K5, "v6_k": K6, "v7_k": VMAX7, "v8_ks": K8, "k9": K9, "k10": K10,
        "k11": VMAX11, "Km11": KM11, "v12_k": K12, "v14_k": K14, "v16_k": VMAX16,
        "v16_Ki": KI16, "v16_n": N16, "Sensitizer_0": SENSITIZER_INITIAL,
    }
    for sbml_id, param in literals.items():
        deposited = deposit["parameters"][sbml_id]
        expected = 0.0 if param is None else float(param)
        if deposited != expected:
            raise ValueError(
                f"parameter {sbml_id!r}: deposit {deposited!r}, transcription {expected!r}")
    for sbml_id, expected in (("SW_Sensitizer", DEPOSIT_SWITCHES["SW_Sensitizer"]),
                              ("SW_HI", DEPOSIT_SWITCHES["SW_HI"]),
                              ("SW_SI", DEPOSIT_SWITCHES["SW_SI"])):
        if deposit["parameters"][sbml_id] != expected:
            raise ValueError(f"switch {sbml_id!r}: deposit {deposit['parameters'][sbml_id]!r}, "
                             f"DEPOSIT_SWITCHES {expected!r}")
    if deposit["parameters"]["ModelValue_7"] != deposit["parameters"]["eps"]:
        raise ValueError(
            "the deposit's VP assignment rule is only the identity on VP0 because ModelValue_7 "
            f"and eps are two copies of one number; they now differ. {SOURCE_DEFECTS['vp_rule_is_the_identity']}")
    shipped = Shock(s1_uM=deposit["parameters"]["s1"], s2_uM=deposit["parameters"]["s2"],
                    t_off_s=deposit["parameters"]["toff"], t_s_s=deposit["parameters"]["ts"])
    if shipped != Shock.deposited():
        raise ValueError(
            f"the deposit's shipped protocol is {shipped}, not the unstressed one this "
            f"transcription reads. {DEPOSIT_SHIPS_UNSTRESSED}")

    calculated = calculated_constants()
    calculated_check = {}
    for name, sbml_id, tolerance in (("k4", "k4", 1e-14), ("km7", "Km7", 1e-14),
                                     ("k6b", "k6b", 1e-14), ("k13", "k13", 1e-11),
                                     ("k15", "k15", 1e-14)):
        deposited = deposit["parameters"][sbml_id]
        relative = abs(calculated[name] - deposited) / abs(deposited)
        if relative > tolerance:
            raise ValueError(
                f"Table S8's calculated row {name!r} recomputes to {calculated[name]!r} against a "
                f"deposited {deposited!r}, relative {relative:.3e} against a tolerance of "
                f"{tolerance:.0e}. {CALCULATED_ROWS[name]}")
        calculated_check[name] = {"recomputed": calculated[name], "deposited": deposited,
                                  "relative": relative, "tolerance": tolerance}

    rate_rule_targets = sorted(deposit["rate_rules"])
    if rate_rule_targets != ["HOGSignal", "Slt2Signal", "Vos"]:
        raise ValueError(f"the deposit's rate rules are {rate_rule_targets}, not the three "
                         "transcribed here; the ODE system has changed")

    if deposit["species"] != authors["species"] or deposit["compartments"] != authors["compartments"]:
        raise ValueError(
            "the BioModels deposit and the authors' Supplementary S2 SBML disagree about species "
            f"or compartments. {DEPOSIT_IDENTITY}")
    value_mismatch = {name for name in set(deposit["parameters"]) | set(authors["parameters"])
                      if deposit["parameters"].get(name) != authors["parameters"].get(name)}
    if value_mismatch:
        raise ValueError(
            f"the two SBML files disagree in VALUE on {sorted(value_mismatch)}; the identity "
            "claim in DEPOSIT_IDENTITY rested on zero value differences")
    reformatted = sorted(name for name in deposit["raw_parameters"]
                         if deposit["raw_parameters"][name] != authors["raw_parameters"].get(name))
    return {
        "sha256": {key: ARTIFACTS[key][1] for key in ARTIFACTS},
        "species": deposit["species"],
        "compartments": deposit["compartments"],
        "parameters_checked": len(literals) + 3,
        "calculated_rows": calculated_check,
        "rate_rules": rate_rule_targets,
        "assignment_rules": len(deposit["assignment_rules"]),
        "reactions": deposit["reactions"],
        "identity": {
            "statement": DEPOSIT_IDENTITY,
            "parameter_strings_reformatted": reformatted,
            "parameter_values_differing": 0,
        },
    }


def _read_sbml(relative: str) -> dict:
    root = ET.parse(paths.data_dir().parent / relative).getroot()
    model = root.find(f"{{{_SBML_NS}}}model")
    species, compartments, parameters, raw = {}, {}, {}, {}
    for element in model.find(f"{{{_SBML_NS}}}listOfSpecies"):
        species[element.get("id")] = (
            float(element.get("initialAmount")),
            element.get("boundaryCondition") == "true",
            element.get("hasOnlySubstanceUnits") == "true")
    for element in model.find(f"{{{_SBML_NS}}}listOfCompartments"):
        compartments[element.get("id")] = float(element.get("size"))
    for element in model.find(f"{{{_SBML_NS}}}listOfParameters"):
        parameters[element.get("id")] = float(element.get("value"))
        raw[element.get("id")] = element.get("value")
    rate_rules, assignment_rules = [], []
    for rule in model.find(f"{{{_SBML_NS}}}listOfRules"):
        target = rule.get("variable")
        if rule.tag == f"{{{_SBML_NS}}}rateRule":
            rate_rules.append(target)
        else:
            assignment_rules.append(target)
    reactions = [element.get("id") for element in model.find(f"{{{_SBML_NS}}}listOfReactions")]
    return {"species": species, "compartments": compartments, "parameters": parameters,
            "raw_parameters": raw, "rate_rules": rate_rules,
            "assignment_rules": assignment_rules, "reactions": reactions}


def reproduction_report() -> dict:
    """Recompute every number the module docstring claims. Nothing here is cached or stored."""
    audit = audit_deposit()
    model = OsmoStat()
    deposit = _read_sbml(ARTIFACTS["sbml_biomodels"][0])
    g, y0 = _derived_geometry(), initial_state()
    vos0, vex = g["Vos0"], g["Vex"]

    closed_form = {}
    for label, computed, deposited in (
            ("f_N2uM", g["f_N2uM"], deposit["parameters"]["N2uM"]),
            ("ci0", g["ci0"], deposit["parameters"]["ci0"]),
            ("cin0", g["ci0"] - float(GLYIN_INITIAL), deposit["parameters"]["cin0"]),
            ("VP0", g["VP0"], deposit["parameters"]["VP0"]),
            ("A0", g["A0"], deposit["parameters"]["A0"]),
            ("Vos0", vos0, deposit["compartments"]["Vos"]),
            ("Vex", vex, deposit["compartments"]["Vex"]),
            ("Hog1", y0["Hog1"], deposit["species"]["Hog1"][0]),
            ("Hog1PP", y0["Hog1PP"], deposit["species"]["Hog1PP"][0]),
            ("Slt2", y0["Slt2"], deposit["species"]["Slt2"][0]),
            ("Slt2P", y0["Slt2P"], deposit["species"]["Slt2P"][0]),
            ("Fps1", y0["Fps1"], deposit["species"]["Fps1"][0]),
            ("Slt2Signal", y0["Slt2Signal"], deposit["parameters"]["Slt2Signal"]),
            ("HOGSignal", y0["HOGSignal"], deposit["parameters"]["HOGSignal"])):
        closed_form[label] = {"published_formula": computed, "deposited": deposited,
                              "relative": abs(computed - deposited) / abs(deposited)}

    residual = model.steady_state_residual()
    without_glyex = {name: row for name, row in residual.items() if name != "Glyex"}

    protocol = {}
    for minutes in (0.75, 1.5, 4.0, 14.0, 30.0, 45.0):
        shock = Shock.published(dilute_after_s=60.0 * minutes)
        hours = (shock.t_off_s + shock.t_s_s + 1800.0) / SECONDS_PER_HOUR
        run = model.timecourse(hours=hours, shock=shock)
        window = run["time_s"] >= shock.t_off_s + shock.t_s_s
        peak_index = int(np.argmax(run["slt2p_relative_percent"][window]))
        protocol[minutes] = {
            "slt2p_peak_percent": float(run["slt2p_relative_percent"][window].max()),
            "seconds_after_dilution": float(run["time_s"][window][peak_index]
                                            - (shock.t_off_s + shock.t_s_s)),
            "hog1pp_max_percent": float(run["hog1pp_relative_percent"].max()),
            "volume_max_after_dilution": float(run["volume_relative"][window].max()),
            "sensitizer_min": float(run["states"]["Sensitizer"].min()),
        }

    vmax16_rows = {}
    published_reading = OsmoStat(vmax16=TABLE_S8_VMAX16)
    for minutes in (4.0, 14.0, 30.0):
        shock = Shock.published(dilute_after_s=60.0 * minutes)
        hours = (shock.t_off_s + shock.t_s_s + 1800.0) / SECONDS_PER_HOUR
        pair = {}
        for label, instance in (("deposit_3_26314", model), ("table_s8_one", published_reading)):
            run = instance.timecourse(hours=hours, shock=shock)
            window = run["time_s"] >= shock.t_off_s + shock.t_s_s
            pair[label] = {"peak_percent": float(run["slt2p_relative_percent"][window].max()),
                           "sensitizer_min": float(run["states"]["Sensitizer"].min())}
        pair["peak_difference_points"] = abs(pair["deposit_3_26314"]["peak_percent"]
                                             - pair["table_s8_one"]["peak_percent"])
        vmax16_rows[minutes] = pair

    unstressed = model.timecourse(hours=24.0, shock=Shock.deposited(), points=241)
    return {
        "audit": audit,
        "availability": AVAILABILITY,
        "axis": AXIS,
        "scope_limit": PANEL_CHANNEL["note"],
        "published_baselines": {
            "slt2p_relative_percent": float(model.relative_slt2_phosphorylation(y0["Slt2P"])),
            "hog1pp_relative_percent": float(model.relative_hog1_phosphorylation(y0["Hog1PP"])),
            "hog1_total_uM": float(HOG1_MOLECULES) * g["f_N2uM"],
            "fps1_total_uM": float(FPS1_MOLECULES) * g["f_N2uM"],
            "published": dict(PUBLISHED),
        },
        "table_s6_closed_forms": closed_form,
        "worst_closed_form_relative": max(row["relative"] for row in closed_form.values()),
        "worst_initial_state_relative": max(
            row["relative"] for name, row in closed_form.items()
            if name not in ("Slt2Signal", "HOGSignal")),
        "signal_roots": {name: closed_form[name] for name in ("Slt2Signal", "HOGSignal")},
        "table_s8_calculated_rows": audit["calculated_rows"],
        "worst_calculated_relative": max(row["relative"]
                                         for row in audit["calculated_rows"].values()),
        "steady_state_residual": residual,
        "worst_relative_rate_excluding_glyex": max(row["relative_per_s"]
                                                   for row in without_glyex.values()),
        "glyex_rate_per_s": residual["Glyex"]["rate_per_s"],
        "unstressed_24h_drift": {
            "glyex_uM": (float(unstressed["states"]["Glyex"][0] / vex),
                         float(unstressed["states"]["Glyex"][-1] / vex)),
            "vos_relative_change": float(abs(unstressed["states"]["Vos"][-1] - vos0) / vos0),
            "slt2p_relative_percent": (float(unstressed["slt2p_relative_percent"][0]),
                                       float(unstressed["slt2p_relative_percent"][-1])),
        },
        "published_protocol": protocol,
        "vmax16_disagreement": vmax16_rows,
        "sign": volume_response(),
        "gate": gate().summary(),
        "defects": dict(SOURCE_DEFECTS),
        "refusals": tuple(sorted(param.name for param in NOT_DRIVEN.refusals())),
    }


def gate():
    """Criterion (e), computed rather than claimed. It refuses, and it should."""
    return TALEMI_PARAMS.gate()


def provenance(heading_level: int = 3) -> str:
    """Both registries as one markdown table: the constants, and what is refused."""
    return provenance_table(TALEMI_PARAMS, NOT_DRIVEN, heading_level=heading_level)
