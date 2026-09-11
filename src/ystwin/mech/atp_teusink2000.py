"""Teusink 2000's yeast glycolysis, transcribed BY HAND IN THIS REPOSITORY.

WHAT THIS FILE IS. Teusink et al. 2000, "Can yeast glycolysis be understood in terms of in
vitro kinetics of the constituent enzymes? Testing biochemistry" (Eur J Biochem 267:5313-29,
DOI 10.1046/j.1432-1327.2000.01527.x, PMID 10951190) is deposited as BioModels
BIOMD0000000064. **The deposit is CC0 1.0** -- its own notes carry the public-domain
dedication verbatim -- so unlike ``data/native_reference_models/jalihal2021/`` it is
redistributable, and it is vendored and hash-pinned at :data:`SOURCE_PATH` /
:data:`SOURCE_DIGEST`. The article itself is Wiley and closed; it was NOT retrievable
(``pdfdirect`` answers HTTP 403), which bounds what the citations below can say -- see
CITATION PRECISION.

**THE EXECUTABLE PYTHON IS OURS.** It is a hand transcription in the idiom of
``mech/ph_ke2013.py`` and ``mech/heat_zheng2016.py``, carrying

    AVAILABILITY = "transcribed_here"

and never ``local_verified``. It is a transcription rather than a load because THIS
REPOSITORY'S OWN LOADER REFUSES THE FILE: ``mech/kinetic_sbml.py`` raises
``UnsupportedSBMLError('assignment target ADP must be a nonconstant global parameter')``,
because the deposit's three assignment rules target SPECIES (ATP, ADP, AMP) rather than
global parameters. That refusal is CORRECT and is left standing -- the loader was not
loosened to admit this file. The three rules are pure algebra (:func:`adenylate`) and are
transcribed instead.

WHAT REPRODUCES. :func:`steady_state` integrates the 14 dynamic species from the deposit's
own initial concentrations and lands on the paper's Table 4, to the precision Table 4 is
printed at. Recomputed by ``tests/test_atp_teusink2000_transcription.py``, never asserted
here (steady-state residual 4.1e-13 mM/min):

    quantity          Teusink Table 4      here        difference
    ATP               2.51 mM              2.5084      -0.06%
    ADP               1.29 mM              1.2921      +0.16%
    AMP               0.30 mM              0.2995      -0.16%
    energy charge     0.769                0.7694      +0.05%
    NAD               1.55 mM              1.5456      -0.29%
    NADH              0.04 mM              0.04444     rounds to the printed 0.04
    P2G               0.04 mM              0.04484     rounds to the printed 0.04
    PYR               8.52 mM              8.5232      +0.04%
    glucose flux      88 mM/min            88.15       +0.17%
    ethanol flux      129 mM/min           129.22      +0.17%
    glycerol flux     18.2 mM/min          18.202      +0.01%
    G6P               1.07 mM              1.0332      **known deposit divergence**

Two rows LOOK large in relative terms and are not: NADH and P2G are printed to two decimals,
so 0.04444 and 0.04484 ARE the published 0.04. Reporting those as 11% and 12% errors would
be reporting the article's print precision as a defect of this port.

The G6P row is the one real disagreement, and it is NOT a transcription error and is not
tuned away: the BioModels curator states it in the deposit notes ("differs slightly from the
one given in table 4") and gives the encoded value as 1.03. This transcription lands on
1.0332, i.e. it reproduces the CC0 ENCODING exactly and inherits the encoding's 3.4% gap to
the printed article. It is carried forward labelled, exactly as ``ph_ke2013.py`` carries
Ke's own Table S5 rest-point inconsistency.

THE ONE DOSABLE ARM, AND WHERE IT STOPS WORKING. Lowering GLCo -- the only panel-dosable
input -- moves the pool monotonically and substantially: ATP 2.508 -> 2.351 -> 1.692 -> 0.780
mM and energy charge 0.769 -> 0.742 -> 0.616 -> 0.396 at GLCo = 50, 20, 5 and 2 mM. **Below
about 1.5 mM the model does not settle**: the integration does not converge in any practical
horizon, which is the source's own well-known "turbo design" instability, not a bug in this
port. The validated domain of this transcription is therefore GLCo >= 2 mM, and deep glucose
starvation is OUTSIDE it. :func:`steady_state` raises rather than returning a number there.

WHAT THIS IS EVIDENCE ABOUT. The port, and nothing else. Table 4 is the SOURCE MODEL'S OWN
STEADY-STATE OUTPUT, not an independent measurement of this project's strain or regime, so
reproducing it scores the transcription and cannot score the biology. :func:`gate` therefore
registers ZERO independent targets against **91** free scalars and REFUSES, which is the same
verdict ``heat_zheng2016.py`` returns for the same reason. The 91 is the deposit's 85 valued
parameters (70 reaction-local + 15 global), plus its two constant species SUM_P and F26BP,
plus the four REFUSED rows below, which are free precisely because nothing pins them.

WHY THIS IS A PRIOR-REPLACEMENT AND NOT A NEW AXIS. ``mech/engine.py`` already carries atp,
adp and amp as dynamic intracellular species with a full stoichiometry (``adenylate_kinase``
at engine.py:274, plus upper/lower glycolysis, ppp, tca, respiration, maintenance). Every
constant driving them is an engine PRIOR in ``_PRIORS``. What this file supplies is the
published alternative, WITH ONE HONEST QUALIFICATION: the grade does not improve. Following
``ph_ke2013.py``, which grades all 78 of Ke's transcribed constants ASSERTED, every constant
here is ASSERTED too. The gain is provenance and an executable reproduction -- an ASSERTED
row citing a pinned CC0 deposit that reproduces its own published table beats an ASSERTED
row stamped "PRIOR, not a published constant" -- not a promotion to MEASURED.

THE ADENYLATE-CONSERVATION QUESTION, ANSWERED BEFORE TRANSCRIPTION RATHER THAN AFTER.
The brief located an "existing adenylate conservation check at engine.py:367". **There is
no such check.** engine.py:367 is the thermodynamic report's pH-match guard for
``adenylate_kinase``; the only conservation checks the engine enforces are carbon and
nitrogen (engine.py:235-238). The question was therefore answered by reading the
stoichiometry directly, and the answer is favourable but not free:

* ATP + ADP + AMP + cAMP is conserved by EVERY metabolic reaction in the engine
  (``adenylate_kinase``, ``maintenance``, ``regulatory_cost``, both glycolysis halves, ppp,
  tca, respiration, ``camp_synthesis``/``camp_hydrolysis`` all balance to zero).
* The single source term is ``growth``, which adds ``initial_atp + initial_adp +
  initial_amp`` as ADP per gDW of new biomass. So the engine's adenylate pool is an
  INTENSIVE conserved moiety -- conserved per gram dry weight, diluted by growth -- while
  Teusink's ``SUM_P`` is a strict conserved moiety in a non-growing preparation. The two are
  compatible in KIND. See :data:`TRANSFER_ASSUMPTIONS`.
* The two disagree in VALUE by 2.6-fold. The engine's priors give 0.0032 mmol/gDW, which at
  its ``cell_water_l_gdw`` = 0.002 L/gDW is **1.6 mM** against Teusink's SUM_P = **4.1 mM**.
* Teusink has NO cAMP. The engine's cAMP is a fourth adenylate reservoir with no counterpart
  in the source, so ``SUM_P`` cannot be mapped onto the engine's three-species sum alone.

FOUR THINGS THE SOURCE DOES NOT SUPPLY. All four are :meth:`Param.refused`, so they raise on
use instead of defaulting, and none of them is filled here:

1. **A stress -> ATP-demand coefficient.** ``KATPASE`` = 33.7 /min is a single fixed
   first-order constant with no stress, temperature, osmotic or nutrient dependence anywhere
   in the deposit. This file supplies the POOL, not the coupling, and
   ``bridge/latent_bridge.py`` already refuses the same constant.
2. **The mM-cytosolic to mmol/gDW bridge.** Teusink is per litre of cytosol; the engine is
   per gDW. The engine's ``cell_water_l_gdw`` = 0.002 L/gDW is an ASSERTED engine prior, not
   a Teusink number, and :func:`to_mmol_per_gdw` refuses rather than laundering it as
   measured.
3. **Any respiratory ATP yield.** The preparation is anaerobic and fermentative; there is no
   oxygen species, no respiratory chain and no P/O ratio in the deposit. Combining these
   glycolytic rate laws with the engine's ``respiratory_atp_yield`` is a transfer, declared.
4. **Any growth term.** No biomass, no dilution, no ``mu``. Teusink's steady state is a
   non-growing one.

CITATION PRECISION, stated because rule 3 of this build is a number attributed to a paper
that does not contain it. The article is paywalled and was not retrieved, so **no constant
below is attributed to a printed table row it was not verified against**. Every one cites
the pinned CC0 deposit by SBML reaction id and parameter id, which is exact and checkable
against :data:`SOURCE_DIGEST`. Only two article-level attributions are made, and both are
made by the deposit's own notes rather than by us: the Vmax set derives from the paper's
Table 1 in U/mg-protein through a stated conversion of "approx. 270", and the steady state
compared against is the paper's Table 4. The per-row split between Teusink's own in-vitro
assays and values he cites from earlier literature is therefore UNVERIFIED HERE, which is
the second reason every row is ASSERTED rather than MEASURED.

TWO MISPRINTS THE DEPOSIT RECORDS IN THE ARTICLE, transcribed in the CORRECTED form the
deposit uses, and named here so the divergence from the printed paper is not silent:
``vADH``'s last denominator term reads ``bpq/(Kib*Kiq*Kp)``, and ``vPFK``'s allosteric
``L`` is ``L0*(..)^2*(..)^2*(..)^2``, not ``L0*(..)^2*(..)^2*(..)``. The deposit also uses
the ADH equilibrium constant for the FORWARD direction, 1/1.45e4 = 6.9e-5, where the article
prints the reverse.

PANEL LIMIT, so nobody reads this as closing a module. Of the panel's two atp stressors only
``glucose_starvation`` (-0.60) is dosable, through ``GLCo`` -> ``external.glucose``. The
heavier ``antimycin_A`` arm (-0.85) is both unrepresentable in ``mech/contracts.py`` AND
outside this source's scope, since the preparation has no respiration to inhibit.

UNITS. Concentrations are mM (the deposit's cytosol compartment is 1.0 L, so its
``cytosol *`` factors are numerically inert and are kept for fidelity). Time is MINUTES, as
every ``mMpermin`` and ``permin`` unit in the deposit prints; no 60x conversion happens here,
because three time bases already collide in this architecture and the conversion belongs to
whoever wires this in, declared, not to the transcription.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .params import Param, ParamRegistry

__all__ = [
    "ADENYLATE_KINASE_EQUILIBRIUM_ONLY",
    "AVAILABILITY",
    "BOUNDARY_CONCENTRATIONS",
    "CYTOSOLIC_WATER_BRIDGE",
    "GROWTH_DILUTION",
    "INITIAL_CONCENTRATIONS",
    "NATIVE_TIME_UNIT",
    "PANEL_LIMIT",
    "PUBLISHED_FLUXES",
    "PUBLISHED_TABLE4",
    "RESPIRATORY_ATP_YIELD",
    "SOURCE_DIGEST",
    "SOURCE_PATH",
    "SPECIES_ORDER",
    "STRESS_ATP_DEMAND",
    "TEUSINK_PARAMS",
    "TRANSFER_ASSUMPTIONS",
    "SteadyState",
    "adenylate",
    "energy_charge",
    "gate",
    "provenance",
    "rates",
    "rhs",
    "steady_state",
    "to_mmol_per_gdw",
    "verify_source",
]

AVAILABILITY = "transcribed_here"

NATIVE_TIME_UNIT = "minute"
"""The deposit's time base, kept. Every rate here is per minute and nothing is converted."""

SOURCE_PATH = "data/native_reference_models/teusink2000/BIOMD0000000064.xml"

SOURCE_DIGEST = "d6a75df939d6707fdb904b44c1f530c8dc3600fd44267ff8e8f38dac560b4361"
"""sha256 of the vendored CC0 deposit. Redistributable outright: the file's own notes carry
the CC0 1.0 dedication, so unlike jalihal2021 this asset needs no custody-only handling."""

SOURCE_LICENCE = "CC0 1.0 Universal (public-domain dedication, stated in the deposit's notes)"

_CITE = ("Teusink et al. 2000, Eur J Biochem 267:5313-29, DOI 10.1046/j.1432-1327.2000.01527.x, "
         "PMID 10951190; encoded as BioModels BIOMD0000000064 (CC0 1.0)")


def _p(reaction: str, parameter: str) -> str:
    """The source location every constant carries: deposit, reaction id, parameter id."""
    return f"{_CITE}; {SOURCE_PATH} reaction <{reaction}> local parameter <{parameter}>"


def _g(parameter: str) -> str:
    return f"{_CITE}; {SOURCE_PATH} global parameter <{parameter}>"


def _s(species: str) -> str:
    return f"{_CITE}; {SOURCE_PATH} species <{species}> initialConcentration"


_VMAX_NOTE = (" Vmax rows derive from the article's Table 1 in U/mg-protein through the "
              "conversion of 'approx. 270' that the deposit's notes state; the pinned value "
              "is the encoded one, and the article was not retrievable to check the row.")

TEUSINK_PARAMS = ParamRegistry("mech/atp_teusink2000.py -- Teusink 2000 glycolysis and the AXP pool")

# ASSERTED throughout, on the ph_ke2013 precedent: a value read out of a published deposit
# whose per-row provenance is not verified here is ASSERTED with the citation attached.
SUM_P = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.sum_p", value=4.1, units="mM",
    source=_s("SUM_P") + ". The total adenylate moiety, held constant. The engine's own "
             "initial_atp+initial_adp+initial_amp is 0.0032 mmol/gDW = 1.6 mM at its "
             "cell_water_l_gdw prior, i.e. 2.6-fold smaller than this.",
    missing="a measured total adenylate pool in this project's strain and fed-batch regime"))

KEQ_AK = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.keq_ak", value=0.45, units="dimensionless",
    source=_g("KeqAK") + ". Adenylate kinase is held at EQUILIBRIUM here, not given a rate; "
                         "the engine instead carries a kinetic adenylate_kinase reaction. "
                         "DIRECTION, verified from the rules rather than assumed: this is "
                         "ATP*AMP/ADP^2 = 0.45, i.e. Keq for 2 ADP -> ATP + AMP. The engine "
                         "writes the SAME reaction the other way round (engine.py:274, "
                         "ATP + AMP -> 2 ADP), so anyone wiring these together must invert "
                         "it or they will invert the adenylate charge response.",
    missing="an adenylate kinase rate constant, which this source does not contain"))

KEQ_TPI = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.keq_tpi", value=0.045, units="dimensionless",
    source=_g("KeqTPI") + ". Triose phosphate isomerase is likewise at equilibrium: TRIO is "
                          "one lumped pool split into GAP and DHAP by this constant."))

F26BP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.f26bp", value=0.02, units="mM",
    source=_s("F26BP") + ". Held constant by the deposit; it is the PFK activator and the "
                         "only place a signalling input could enter this model.",
    missing="a published F2,6BP response to any stressor this panel can dose"))

# --- PFK allosteric block: global parameters ---
GR = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.gR", value=5.12, units="dimensionless", source=_g("gR")))
KM_PFK_F6P = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.KmF6P", value=0.1, units="mM", source=_g("KmPFKF6P")))
KM_PFK_ATP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.KmATP", value=0.71, units="mM", source=_g("KmPFKATP")))
L_ZERO = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.Lzero", value=0.66, units="dimensionless", source=_g("Lzero")))
CI_PFK_ATP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.CiATP", value=100.0, units="dimensionless", source=_g("CiPFKATP")))
KI_PFK_ATP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.KiATP", value=0.65, units="mM", source=_g("KiPFKATP")))
C_PFK_AMP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.CAMP", value=0.0845, units="dimensionless", source=_g("CPFKAMP")))
K_PFK_AMP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.KAMP", value=0.0995, units="mM", source=_g("KPFKAMP") +
    ". This is the ONLY place AMP feeds back on flux in the whole model."))
C_PFK_F26BP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.CF26BP", value=0.0174, units="dimensionless", source=_g("CPFKF26BP")))
K_PFK_F26BP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.KF26BP", value=0.000682, units="mM", source=_g("KPFKF26BP")))
C_PFK_F16BP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.CF16BP", value=0.397, units="dimensionless", source=_g("CPFKF16BP")))
K_PFK_F16BP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.KF16BP", value=0.111, units="mM", source=_g("KPFKF16BP")))
C_PFK_ATP = TEUSINK_PARAMS.add(Param.asserted(
    name="teusink.pfk.CATP", value=3.0, units="dimensionless", source=_g("CPFKATP")))

# --- per-reaction local parameters, reaction by reaction ---
_LOCAL: dict[str, tuple[tuple[str, float, str], ...]] = {
    "vGLT": (("VmGLT", 97.264, "mmol/min"), ("KmGLTGLCo", 1.1918, "mM"),
             ("KeqGLT", 1.0, "dimensionless"), ("KmGLTGLCi", 1.1918, "mM")),
    "vGLK": (("VmGLK", 226.452, "mM/min"), ("KmGLKGLCi", 0.08, "mM"),
             ("KmGLKATP", 0.15, "mM"), ("KeqGLK", 3800.0, "dimensionless"),
             ("KmGLKG6P", 30.0, "mM"), ("KmGLKADP", 0.23, "mM")),
    "vPGI": (("VmPGI_2", 339.677, "mM/min"), ("KmPGIG6P_2", 1.4, "mM"),
             ("KeqPGI_2", 0.314, "dimensionless"), ("KmPGIF6P_2", 0.3, "mM")),
    "vGLYCO": (("KGLYCOGEN_3", 6.0, "mM/min"),),
    "vTreha": (("KTREHALOSE", 2.4, "mM/min"),),
    "vPFK": (("VmPFK", 182.903, "mM/min"),),
    "vALD": (("VmALD", 322.258, "mM/min"), ("KmALDF16P", 0.3, "mM"),
             ("KeqALD", 0.069, "dimensionless"), ("KmALDGAP", 2.0, "mM"),
             ("KmALDDHAP", 2.4, "mM"), ("KmALDGAPi", 10.0, "mM")),
    "vGAPDH": (("VmGAPDHf", 1184.52, "mM/min"), ("KmGAPDHGAP", 0.21, "mM"),
               ("KmGAPDHNAD", 0.09, "mM"), ("VmGAPDHr", 6549.8, "mM/min"),
               ("KmGAPDHBPG", 0.0098, "mM"), ("KmGAPDHNADH", 0.06, "mM")),
    "vPGK": (("VmPGK", 1306.45, "mM/min"), ("KmPGKP3G", 0.53, "mM"),
             ("KmPGKATP", 0.3, "mM"), ("KeqPGK", 3200.0, "dimensionless"),
             ("KmPGKBPG", 0.003, "mM"), ("KmPGKADP", 0.2, "mM")),
    "vPGM": (("VmPGM", 2525.81, "mM/min"), ("KmPGMP3G", 1.2, "mM"),
             ("KeqPGM", 0.19, "dimensionless"), ("KmPGMP2G", 0.08, "mM")),
    "vENO": (("VmENO", 365.806, "mM/min"), ("KmENOP2G", 0.04, "mM"),
             ("KeqENO", 6.7, "dimensionless"), ("KmENOPEP", 0.5, "mM")),
    "vPYK": (("VmPYK", 1088.71, "mM/min"), ("KmPYKPEP", 0.14, "mM"),
             ("KmPYKADP", 0.53, "mM"), ("KeqPYK", 6500.0, "dimensionless"),
             ("KmPYKPYR", 21.0, "mM"), ("KmPYKATP", 1.5, "mM")),
    "vPDC": (("VmPDC", 174.194, "mM/min"), ("nPDC", 1.9, "dimensionless"),
             ("KmPDCPYR", 4.33, "mM")),
    "vSUC": (("KSUCC", 21.4, "1/min"),),
    "vADH": (("VmADH", 810.0, "mM/min"), ("KiADHNAD", 0.92, "mM"),
             ("KmADHETOH", 17.0, "mM"), ("KeqADH", 6.9e-05, "dimensionless"),
             ("KmADHNAD", 0.17, "mM"), ("KmADHNADH", 0.11, "mM"),
             ("KiADHNADH", 0.031, "mM"), ("KmADHACE", 1.11, "mM"),
             ("KiADHACE", 1.1, "mM"), ("KiADHETOH", 90.0, "mM")),
    "vG3PDH": (("VmG3PDH", 70.15, "mM/min"), ("KmG3PDHDHAP", 0.4, "mM"),
               ("KmG3PDHNADH", 0.023, "mM"), ("KeqG3PDH", 4300.0, "dimensionless"),
               ("KmG3PDHGLY", 1.0, "mM"), ("KmG3PDHNAD", 0.93, "mM")),
    "vATP": (("KATPASE", 33.7, "1/min"),),
}

K: dict[str, float] = {}
for _reaction, _rows in _LOCAL.items():
    for _pid, _value, _units in _rows:
        _note = _VMAX_NOTE if _pid.startswith("Vm") else ""
        if _pid == "KATPASE":
            _note = (" The whole ATP demand of the model, first-order in ATP and FIXED. It "
                     "carries no stress, temperature or nutrient dependence of any kind.")
        if _pid == "KeqADH":
            _note = (" The deposit's notes record that the article prints Keq for the REVERSE "
                     "direction (1.45e4); this is the forward 1/Keq.")
        TEUSINK_PARAMS.add(Param.asserted(
            name=f"teusink.{_reaction}.{_pid}", value=_value, units=_units,
            source=_p(_reaction, _pid) + _note))
        K[_pid] = _value

# --- the four quantities the source does not contain ---
STRESS_ATP_DEMAND = TEUSINK_PARAMS.add(Param.refused(
    name="teusink.stress_to_atp_demand", units="dimensionless per stress unit",
    source=f"{_CITE}. Searched: every one of the deposit's 85 valued parameters, all 17 "
           "kinetic laws and the model notes. KATPASE is a single fixed first-order constant "
           "and no stress, temperature, osmotic or nutrient term appears anywhere in the "
           "model. bridge/latent_bridge.py already refuses this same coefficient.",
    reason="the source supplies the adenylate POOL, not a coupling from any stress axis to "
           "ATP demand; inventing one is inventing kinetics",
    missing="a published, parameterised S. cerevisiae law relating a dosable stressor to "
            "ATP turnover, or a direct ATPase-flux measurement under that stressor"))

CYTOSOLIC_WATER_BRIDGE = TEUSINK_PARAMS.add(Param.refused(
    name="teusink.cytosol_l_per_gdw", units="L cytosol/gDW",
    source=f"{_CITE}. The deposit is per litre of cytosol (compartment 'cytosol', size 1.0) "
           "and contains no biomass, no dry weight and no cell count. The engine's "
           "cell_water_l_gdw = 0.002 L/gDW is an engine _PRIORS row, not a Teusink number. "
           "The notes' 'approx. 270' converts U/mg-protein to mM/min and is a protein-per-"
           "cytosol-volume factor, not a dry-weight bridge.",
    reason="converting mM cytosolic to mmol/gDW through an asserted engine prior would "
           "launder that prior as if the published source had supplied it",
    missing="a measured cytosolic water volume per gram dry weight for this strain and "
            "regime, which would also close the engine's own asserted cell_water_l_gdw"))

RESPIRATORY_ATP_YIELD = TEUSINK_PARAMS.add(Param.refused(
    name="teusink.respiratory_atp_yield", units="mol ATP/mol NADH",
    source=f"{_CITE}. The preparation is a non-respiring, anaerobic, fermentative one: the "
           "deposit has no oxygen species, no respiratory chain, no mitochondrion and no P/O "
           "ratio. Its NADH is reoxidised only by vADH and vG3PDH.",
    reason="combining these glycolytic rate laws with the engine's respiratory ATP yield is "
           "a transfer between regimes the source never occupies",
    missing="a parameterised S. cerevisiae respiratory model, or a measured P/O ratio in the "
            "aerobic fed-batch this project runs"))

GROWTH_DILUTION = TEUSINK_PARAMS.add(Param.refused(
    name="teusink.growth_dilution", units="1/min",
    source=f"{_CITE}. The steady state is a NON-GROWING one: no biomass state, no dilution "
           "term and no mu appears in any of the 17 rate laws. The engine by contrast adds "
           "initial_atp+initial_adp+initial_amp as ADP per gDW of new biomass in its growth "
           "reaction, making its adenylate pool intensive rather than strictly conserved.",
    reason="the source's conserved SUM_P holds only because nothing grows; asserting a "
           "dilution rate for it would be inventing the term the source omits",
    missing="a published glycolytic parameterisation measured in a growing culture"))

ADENYLATE_KINASE_EQUILIBRIUM_ONLY = (
    "Teusink holds adenylate kinase at equilibrium (KeqAK = 0.45, an algebraic rule). The "
    "engine runs it as a kinetic reaction with an asserted rate of 20 mmol/gDW/h. This "
    "source cannot supply that rate: an equilibrium constant is not a rate constant.")

TRANSFER_ASSUMPTIONS = {
    "non_growing_to_fed_batch": (
        "Teusink's steady state is non-growing; this project's regime is a growing fed-batch. "
        "SUM_P is strictly conserved there and only intensively conserved (per gDW, diluted by "
        "growth) in the engine. Declared, not repaired -- see GROWTH_DILUTION."),
    "fermentative_to_aerobic": (
        "The preparation is anaerobic and fermentative. Any use of these rate laws alongside "
        "the engine's respiration reaction is a transfer across regimes -- see "
        "RESPIRATORY_ATP_YIELD."),
    "mM_to_mmol_per_gdw": (
        "Per litre of cytosol here, per gram dry weight in the engine. The bridge is REFUSED, "
        "not defaulted -- see CYTOSOLIC_WATER_BRIDGE and to_mmol_per_gdw."),
    "no_camp_reservoir": (
        "The engine's cAMP is a fourth adenylate reservoir that Teusink does not have, so "
        "SUM_P = 4.1 mM does not map onto the engine's atp+adp+amp sum alone."),
    "adenylate_kinase_equilibrium": ADENYLATE_KINASE_EQUILIBRIUM_ONLY,
}

PANEL_LIMIT = (
    "Of the panel's two atp stressors only glucose_starvation (-0.60) is dosable, through "
    "GLCo -> external.glucose. antimycin_A (-0.85) is unrepresentable in mech/contracts.py AND "
    "outside this source's scope: the preparation has no respiration to inhibit. This file "
    "does not close the atp module.")

SPECIES_ORDER = ("GLCi", "G6P", "F6P", "F16P", "TRIO", "BPG", "P3G", "P2G",
                 "PEP", "PYR", "ACE", "P", "NAD", "NADH")
"""The 14 dynamic species. ATP, ADP and AMP are NOT states: they are algebraic in P."""

INITIAL_CONCENTRATIONS = {
    "GLCi": 0.087, "G6P": 2.45, "F6P": 0.62, "F16P": 5.51, "TRIO": 0.96, "BPG": 0.0,
    "P3G": 0.9, "P2G": 0.12, "PEP": 0.07, "PYR": 1.85, "ACE": 0.17, "P": 6.31,
    "NAD": 1.2, "NADH": 0.39,
}
"""Deposit initialConcentration rows, mM. P = 2*ATP + ADP is the high-energy phosphate pool."""

BOUNDARY_CONCENTRATIONS = {"GLCo": 50.0, "ETOH": 50.0, "GLY": 0.15}
"""Boundary species, held. GLCo is the one dosable input; Glyc/Trh/CO2/SUCC are pure sinks."""

VALIDATED_GLCO_FLOOR_mM = 2.0
"""The lowest glucose this transcription was VERIFIED to settle at. Below roughly 1.5 mM the
source's turbo-design instability stops it converging, so steady_state refuses rather than
returning an unconverged number. A numerical domain limit, not a fitted parameter."""

PUBLISHED_TABLE4 = {
    "ATP": 2.51, "ADP": 1.29, "AMP": 0.30, "NAD": 1.55, "NADH": 0.04,
    "G6P": 1.07, "F6P": 0.11, "F16P": 0.6, "P3G": 0.36, "P2G": 0.04,
    "PEP": 0.07, "PYR": 8.52, "ACE": 0.17,
}
"""Article Table 4 steady-state concentrations in mM, as the deposit's notes reprint them.
G6P is the row the deposit itself flags as diverging: it encodes 1.03 against the printed 1.07."""

DEPOSIT_G6P = 1.03
"""The curator's own value for the one row the CC0 encoding and the article disagree on."""

PUBLISHED_FLUXES = {"glucose": 88.0, "ethanol": 129.0, "glycogen": 6.0,
                    "trehalose": 4.8, "glycerol": 18.2, "succinate": 3.6}
"""Article Table 4 steady-state fluxes, mM/min, as the deposit's notes reprint them.
'trehalose' is the G6P flux through that branch, i.e. twice vTreha."""


def verify_source(root: Path | str | None = None) -> str:
    """sha256 the vendored CC0 deposit and refuse a mismatch. Raises if it is absent."""
    base = Path(root) if root is not None else Path(__file__).resolve().parents[3]
    path = base / SOURCE_PATH
    if not path.is_file():
        raise FileNotFoundError(
            f"{SOURCE_PATH} is not present. It is CC0 1.0 and redistributable; fetch it from "
            "https://www.ebi.ac.uk/biomodels/BIOMD0000000064 if this is a partial checkout")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != SOURCE_DIGEST:
        raise ValueError(
            f"{SOURCE_PATH} sha256 is {digest}, expected {SOURCE_DIGEST}. The transcription "
            "below was made from the pinned bytes and cannot claim a changed file")
    return digest


def to_mmol_per_gdw(concentration_mM: float) -> float:
    """REFUSED. There is no Teusink bridge from mM cytosolic to mmol/gDW; see the module
    docstring. Call this to get the refusal with its citation, never a silent default."""
    return concentration_mM * float(CYTOSOLIC_WATER_BRIDGE)


def adenylate(phosphate_mM: float | np.ndarray, *, sum_p: float = 4.1,
              keq_ak: float = 0.45) -> tuple:
    """The deposit's three assignment rules, verbatim, as pure algebra.

    ``kinetic_sbml.from_sbml`` refuses these rules because they target species rather than
    global parameters, which is why they are transcribed here rather than loaded. Adenylate
    kinase is held at equilibrium: ADP solves the quadratic, ATP and AMP follow.

    Args:
        phosphate_mM: the high-energy phosphate pool P = 2*ATP + ADP, mM.
        sum_p: total adenylate moiety, mM (deposit species SUM_P, constant).
        keq_ak: adenylate kinase equilibrium constant (deposit global KeqAK).

    Returns:
        ``(ATP, ADP, AMP)`` in mM.
    """
    p = phosphate_mM
    denominator = 1.0 - 4.0 * keq_ak
    root = p ** 2 * denominator + 2.0 * sum_p * p * (4.0 * keq_ak - 1.0) + sum_p ** 2
    adp = (sum_p - np.sqrt(np.maximum(root, 0.0))) / denominator
    atp = (p - adp) / 2.0
    amp = sum_p - atp - adp
    return atp, adp, amp


def energy_charge(atp: float, adp: float, *, sum_p: float = 4.1) -> float:
    """(ATP + ADP/2) / SUM_P. Atkinson's charge on the deposit's own pool."""
    return (atp + 0.5 * adp) / sum_p


def _pfk_R(f6p: float, atp: float) -> float:
    """Deposit function definition R_PFK, with the article misprint corrected as the deposit
    records: R = 1 + l1 + l2 + gR*l1*l2."""
    return (1.0 + f6p / float(KM_PFK_F6P) + atp / float(KM_PFK_ATP)
            + float(GR) * (f6p / float(KM_PFK_F6P)) * (atp / float(KM_PFK_ATP)))


def _pfk_T(atp: float) -> float:
    """Deposit function definition T_PFK."""
    return 1.0 + float(C_PFK_ATP) * (atp / float(KM_PFK_ATP))


def _pfk_L(atp: float, amp: float, f16p: float, f26bp: float) -> float:
    """Deposit function definition L_PFK. The third bracket is squared, which is the second
    misprint the deposit records against the article."""
    a = ((1.0 + float(CI_PFK_ATP) * (atp / float(KI_PFK_ATP)))
         / (1.0 + atp / float(KI_PFK_ATP))) ** 2
    b = ((1.0 + float(C_PFK_AMP) * (amp / float(K_PFK_AMP)))
         / (1.0 + amp / float(K_PFK_AMP))) ** 2
    c = ((1.0 + float(C_PFK_F26BP) * f26bp / float(K_PFK_F26BP)
          + float(C_PFK_F16BP) * f16p / float(K_PFK_F16BP))
         / (1.0 + f26bp / float(K_PFK_F26BP) + f16p / float(K_PFK_F16BP))) ** 2
    return float(L_ZERO) * a * b * c


def rates(state: dict, *, glco: float | None = None, etoh: float | None = None,
          gly: float | None = None, f26bp: float | None = None,
          sum_p: float | None = None, keq_ak: float | None = None) -> dict:
    """Every one of the deposit's 17 kinetic laws, transcribed verbatim, in mM/min.

    Args:
        state: the 14 dynamic species in mM, keyed by :data:`SPECIES_ORDER`.
        glco, etoh, gly: boundary concentrations, mM. Default to
            :data:`BOUNDARY_CONCENTRATIONS`. ``glco`` is the one dosable input.
        f26bp, sum_p, keq_ak: the deposit's constants, overridable for a declared sweep.

    Returns:
        A dict of reaction id -> rate in mM/min, plus ATP/ADP/AMP in mM.
    """
    glco = BOUNDARY_CONCENTRATIONS["GLCo"] if glco is None else glco
    etoh = BOUNDARY_CONCENTRATIONS["ETOH"] if etoh is None else etoh
    gly = BOUNDARY_CONCENTRATIONS["GLY"] if gly is None else gly
    f26bp = float(F26BP) if f26bp is None else f26bp
    sum_p = float(SUM_P) if sum_p is None else sum_p
    keq_ak = float(KEQ_AK) if keq_ak is None else keq_ak

    s = state
    glci, g6p, f6p, f16p = s["GLCi"], s["G6P"], s["F6P"], s["F16P"]
    trio, bpg, p3g, p2g = s["TRIO"], s["BPG"], s["P3G"], s["P2G"]
    pep, pyr, ace, phos = s["PEP"], s["PYR"], s["ACE"], s["P"]
    nad, nadh = s["NAD"], s["NADH"]
    atp, adp, amp = adenylate(phos, sum_p=sum_p, keq_ak=keq_ak)

    keq_tpi = float(KEQ_TPI)
    gap = (keq_tpi / (1.0 + keq_tpi)) * trio
    dhap = (1.0 / (1.0 + keq_tpi)) * trio

    v = {}
    v["vGLT"] = ((K["VmGLT"] / K["KmGLTGLCo"]) * (glco - glci / K["KeqGLT"])
                 / (1.0 + glco / K["KmGLTGLCo"] + glci / K["KmGLTGLCi"]
                    + 0.91 * glco * glci / (K["KmGLTGLCo"] * K["KmGLTGLCi"])))
    v["vGLK"] = ((K["VmGLK"] / (K["KmGLKGLCi"] * K["KmGLKATP"]))
                 * (glci * atp - g6p * adp / K["KeqGLK"])
                 / ((1.0 + glci / K["KmGLKGLCi"] + g6p / K["KmGLKG6P"])
                    * (1.0 + atp / K["KmGLKATP"] + adp / K["KmGLKADP"])))
    v["vPGI"] = ((K["VmPGI_2"] / K["KmPGIG6P_2"]) * (g6p - f6p / K["KeqPGI_2"])
                 / (1.0 + g6p / K["KmPGIG6P_2"] + f6p / K["KmPGIF6P_2"]))
    v["vGLYCO"] = K["KGLYCOGEN_3"]
    v["vTreha"] = K["KTREHALOSE"]
    r_pfk = _pfk_R(f6p, atp)
    v["vPFK"] = (K["VmPFK"] * float(GR) * (f6p / float(KM_PFK_F6P)) * (atp / float(KM_PFK_ATP))
                 * r_pfk / (r_pfk ** 2 + _pfk_L(atp, amp, f16p, f26bp) * _pfk_T(atp) ** 2))
    v["vALD"] = ((K["VmALD"] / K["KmALDF16P"])
                 * (f16p - gap * dhap / K["KeqALD"])
                 / (1.0 + f16p / K["KmALDF16P"] + gap / K["KmALDGAP"] + dhap / K["KmALDDHAP"]
                    + gap * dhap / (K["KmALDGAP"] * K["KmALDDHAP"])
                    + f16p * gap / (K["KmALDGAPi"] * K["KmALDF16P"])))
    v["vGAPDH"] = ((K["VmGAPDHf"] * gap * nad / (K["KmGAPDHGAP"] * K["KmGAPDHNAD"])
                    - K["VmGAPDHr"] * bpg * nadh / (K["KmGAPDHBPG"] * K["KmGAPDHNADH"]))
                   / ((1.0 + gap / K["KmGAPDHGAP"] + bpg / K["KmGAPDHBPG"])
                      * (1.0 + nad / K["KmGAPDHNAD"] + nadh / K["KmGAPDHNADH"])))
    v["vPGK"] = ((K["VmPGK"] / (K["KmPGKP3G"] * K["KmPGKATP"]))
                 * (K["KeqPGK"] * bpg * adp - p3g * atp)
                 / ((1.0 + bpg / K["KmPGKBPG"] + p3g / K["KmPGKP3G"])
                    * (1.0 + atp / K["KmPGKATP"] + adp / K["KmPGKADP"])))
    v["vPGM"] = ((K["VmPGM"] / K["KmPGMP3G"]) * (p3g - p2g / K["KeqPGM"])
                 / (1.0 + p3g / K["KmPGMP3G"] + p2g / K["KmPGMP2G"]))
    v["vENO"] = ((K["VmENO"] / K["KmENOP2G"]) * (p2g - pep / K["KeqENO"])
                 / (1.0 + p2g / K["KmENOP2G"] + pep / K["KmENOPEP"]))
    v["vPYK"] = ((K["VmPYK"] / (K["KmPYKPEP"] * K["KmPYKADP"]))
                 * (pep * adp - pyr * atp / K["KeqPYK"])
                 / ((1.0 + pep / K["KmPYKPEP"] + pyr / K["KmPYKPYR"])
                    * (1.0 + atp / K["KmPYKATP"] + adp / K["KmPYKADP"])))
    v["vPDC"] = (K["VmPDC"] * (pyr ** K["nPDC"] / K["KmPDCPYR"] ** K["nPDC"])
                 / (1.0 + pyr ** K["nPDC"] / K["KmPDCPYR"] ** K["nPDC"]))
    v["vSUC"] = K["KSUCC"] * ace
    v["vADH"] = -((K["VmADH"] / (K["KiADHNAD"] * K["KmADHETOH"]))
                  * (nad * etoh - nadh * ace / K["KeqADH"])
                  / (1.0 + nad / K["KiADHNAD"]
                     + K["KmADHNAD"] * etoh / (K["KiADHNAD"] * K["KmADHETOH"])
                     + K["KmADHNADH"] * ace / (K["KiADHNADH"] * K["KmADHACE"])
                     + nadh / K["KiADHNADH"]
                     + nad * etoh / (K["KiADHNAD"] * K["KmADHETOH"])
                     + K["KmADHNADH"] * nad * ace
                     / (K["KiADHNAD"] * K["KiADHNADH"] * K["KmADHACE"])
                     + K["KmADHNAD"] * etoh * nadh
                     / (K["KiADHNAD"] * K["KmADHETOH"] * K["KiADHNADH"])
                     + nadh * ace / (K["KiADHNADH"] * K["KmADHACE"])
                     + nad * etoh * ace / (K["KiADHNAD"] * K["KmADHETOH"] * K["KiADHACE"])
                     + etoh * nadh * ace
                     / (K["KiADHETOH"] * K["KiADHNADH"] * K["KmADHACE"])))
    v["vG3PDH"] = ((K["VmG3PDH"] / (K["KmG3PDHDHAP"] * K["KmG3PDHNADH"]))
                   * (dhap * nadh - gly * nad / K["KeqG3PDH"])
                   / ((1.0 + dhap / K["KmG3PDHDHAP"] + gly / K["KmG3PDHGLY"])
                      * (1.0 + nadh / K["KmG3PDHNADH"] + nad / K["KmG3PDHNAD"])))
    v["vATP"] = K["KATPASE"] * atp
    v["ATP"], v["ADP"], v["AMP"] = atp, adp, amp
    return v


def rhs(_t: float, y: np.ndarray, **kwargs) -> np.ndarray:
    """d[species]/dt in mM/min, from the deposit's reaction stoichiometry.

    vSUC consumes 2 ACE, 3 NAD and 4 P per succinate, which is the stoichiometry the deposit
    states it restored to match the publication; vTreha consumes 2 G6P and 1 P.
    """
    s = dict(zip(SPECIES_ORDER, y))
    v = rates(s, **kwargs)
    return np.array([
        v["vGLT"] - v["vGLK"],
        v["vGLK"] - v["vPGI"] - v["vGLYCO"] - 2.0 * v["vTreha"],
        v["vPGI"] - v["vPFK"],
        v["vPFK"] - v["vALD"],
        2.0 * v["vALD"] - v["vGAPDH"] - v["vG3PDH"],
        v["vGAPDH"] - v["vPGK"],
        v["vPGK"] - v["vPGM"],
        v["vPGM"] - v["vENO"],
        v["vENO"] - v["vPYK"],
        v["vPYK"] - v["vPDC"],
        v["vPDC"] - v["vADH"] - 2.0 * v["vSUC"],
        (-v["vGLK"] - v["vGLYCO"] - v["vTreha"] - v["vPFK"] + v["vPGK"] + v["vPYK"]
         - 4.0 * v["vSUC"] - v["vATP"]),
        -v["vGAPDH"] + v["vADH"] + v["vG3PDH"] - 3.0 * v["vSUC"],
        v["vGAPDH"] - v["vADH"] - v["vG3PDH"] + 3.0 * v["vSUC"],
    ])


@dataclass(frozen=True)
class SteadyState:
    """One integrated steady state, with the derived quantities this build exists for."""

    concentrations: dict
    fluxes: dict
    atp: float
    adp: float
    amp: float
    energy_charge: float
    residual: float
    glco: float

    def table4(self) -> dict:
        """The rows comparable with :data:`PUBLISHED_TABLE4`."""
        rows = {k: self.concentrations[k] for k in
                ("G6P", "F6P", "F16P", "P3G", "P2G", "PEP", "PYR", "ACE", "NAD", "NADH")}
        rows.update({"ATP": self.atp, "ADP": self.adp, "AMP": self.amp})
        return rows


def steady_state(*, glco: float | None = None, horizon_min: float = 2000.0,
                 rtol: float = 1e-10, atol: float = 1e-12,
                 residual_tol: float = 1e-6, **kwargs) -> SteadyState:
    """Integrate the deposit's own initial condition to steady state.

    Args:
        glco: extracellular glucose, mM. Defaults to the deposit's 50 mM. This is the ONLY
            panel-dosable input, and the deposit is known to stall near 1 mM.
        horizon_min: integration horizon in MINUTES, the deposit's native time base.
        residual_tol: the largest |d[species]/dt| that still counts as a steady state, mM/min.
            Exceeding it RAISES: a value read off an unconverged trajectory would be named
            for the horizon it was taken at, which is the defect this module's docstring
            criticises in Krakowiak's reporter sweep.
    """
    glco = BOUNDARY_CONCENTRATIONS["GLCo"] if glco is None else glco
    if glco < VALIDATED_GLCO_FLOOR_mM:
        raise ValueError(
            f"GLCo = {glco} mM is below this transcription's validated floor of "
            f"{VALIDATED_GLCO_FLOOR_mM} mM. The source model does not settle there -- its "
            "own 'turbo design' instability, not a defect of this port -- so no steady state "
            "is returned rather than one that did not converge. Deep glucose starvation is "
            "outside the source's validated domain")
    y0 = np.array([INITIAL_CONCENTRATIONS[name] for name in SPECIES_ORDER], dtype=float)
    solution = solve_ivp(lambda t, y: rhs(t, y, glco=glco, **kwargs), (0.0, horizon_min), y0,
                         method="LSODA", rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError(f"Teusink transcription did not integrate at GLCo = {glco} mM: "
                           f"{solution.message}")
    final = solution.y[:, -1]
    if not np.all(np.isfinite(final)):
        raise RuntimeError(f"Teusink transcription diverged at GLCo = {glco} mM")
    concentrations = dict(zip(SPECIES_ORDER, final))
    v = rates(concentrations, glco=glco, **kwargs)
    residual = float(np.max(np.abs(rhs(horizon_min, final, glco=glco, **kwargs))))
    if residual > residual_tol:
        raise RuntimeError(
            f"no steady state at GLCo = {glco} mM within {horizon_min} min: largest "
            f"|d[species]/dt| is {residual:.3g} mM/min against a tolerance of "
            f"{residual_tol:g}. Returning this would report a horizon, not a steady state")
    atp, adp, amp = v["ATP"], v["ADP"], v["AMP"]
    sum_p = kwargs.get("sum_p") or float(SUM_P)
    fluxes = {"glucose": v["vGLT"], "ethanol": v["vADH"], "glycogen": v["vGLYCO"],
              "trehalose": 2.0 * v["vTreha"], "glycerol": v["vG3PDH"],
              "succinate": v["vSUC"], "atpase": v["vATP"]}
    return SteadyState(concentrations=concentrations, fluxes=fluxes, atp=float(atp),
                       adp=float(adp), amp=float(amp),
                       energy_charge=energy_charge(float(atp), float(adp), sum_p=sum_p),
                       residual=residual, glco=glco)


def gate():
    """Criterion (e), computed rather than claimed.

    ZERO independent targets are registered on purpose. Table 4 is the source model's OWN
    steady-state output, so reproducing it scores this transcription and cannot score the
    biology; registering it as a target would let the port pay for itself.
    """
    return TEUSINK_PARAMS.gate()


def provenance() -> dict:
    """What this file may claim, in one record."""
    return {
        "availability": AVAILABILITY,
        "source": _CITE,
        "licence": SOURCE_LICENCE,
        "source_path": SOURCE_PATH,
        "source_sha256": SOURCE_DIGEST,
        "loader_refusal": ("mech/kinetic_sbml.py raises UnsupportedSBMLError: the three "
                           "assignment rules target species, not global parameters. The "
                           "loader was NOT loosened; the rules are transcribed as algebra."),
        "time_unit": NATIVE_TIME_UNIT,
        "valued_parameters": 85,
        "valued_parameter_split": "70 reaction-local + 15 global; SUM_P and F26BP are "
                                  "constant species, counted separately",
        "registered_rows": len(TEUSINK_PARAMS),
        "by_tag": TEUSINK_PARAMS.by_tag(),
        "validated_glucose_floor_mM": VALIDATED_GLCO_FLOOR_mM,
        "refusals": tuple(p.name for p in TEUSINK_PARAMS.refusals()),
        "transfer_assumptions": TRANSFER_ASSUMPTIONS,
        "panel_limit": PANEL_LIMIT,
        "gate": gate().summary(),
    }
