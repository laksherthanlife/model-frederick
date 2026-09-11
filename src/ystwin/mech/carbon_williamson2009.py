"""Williamson 2009's Complete cAMP Model, transcribed BY HAND IN THIS REPOSITORY.

WHAT THIS FILE IS. Williamson T, Schwartz J-M, Kell DB & Stateva L 2009, "Deterministic
mathematical models of the cAMP pathway in *Saccharomyces cerevisiae*", BMC Systems Biology
3:70 (doi:10.1186/1752-0509-3-70, PMID 19607691, PMC2719611, **CC BY 2.0**) ships its
Complete cAMP Model as Additional file 1, an SBML Level 2 Version 1 document with 24 species,
20 reactions and 41 global parameters -- 39 kinetic constants plus two assignment-rule outputs. That file is vendored and checksummed at
``data/native_reference_models/williamson2009/1752-0509-3-70-S1.xml``, sha256
``e608d665...bcd1b782``. What follows is a hand transcription of it into this engine's idiom.
**The executable is ours**, which is why it carries ``AVAILABILITY = "transcribed_here"`` and
never ``local_verified`` -- that level in this repository means a checksummed publisher file
loaded unchanged by the shared SBML loader, and the shared loader refuses this one (below).

Nothing here was typed from memory. Every constant, initial concentration and stoichiometric
coefficient in this module is re-read out of the vendored SBML by :func:`audit_deposit` at
construction and compared to the literal written here; a mismatch raises rather than running.

WHAT REPRODUCES. Williamson prints two figures of the Complete cAMP Model: Figure 8 (the fit,
cAMP against time under 5 mM glucose at 60 s then 100 mM at 240 s) and Figure 9 (four panels
of predicted species and rates over 600 s under the same protocol). Thirty-six values read off
those five panels are reproduced here, under the deposited parameterisation, to a **mean 0.30%
and a worst 1.61% of the panel's own y-axis full scale** -- which is inside the precision at
which a printed figure can be read at all. :func:`reproduction_report` recomputes all of them
and ``tests/test_carbon_williamson2009_transcription.py`` asserts them; none was tuned.

A second, independent check that needs no figure: **the deposited initial concentrations are a
steady state of the deposited rate laws**, and this transcription lands them at a maximum
relative rate of 3.5e-14 /s. That is a machine-precision residual across seventeen dynamic
species, so it exercises every rate law and nearly every stoichiometric coefficient at once.
It does NOT exercise the PKA stoichiometry -- see the next paragraph -- and this module says so
rather than letting one check appear to cover everything.

THE ONE TRANSCRIPTION ERROR THIS CAUGHT, recorded because it is the reason both checks exist.
``PKAact`` writes ``<speciesReference species="C" stoichiometry="2"/>``: one holoenzyme releases
TWO catalytic subunits, and ``PKAdeact`` consumes two. Read as 1:1 the steady-state residual is
unchanged -- at a fixed point the two fluxes balance either way -- and Figure 8's first peak is
only 16% high at 2.87e-3 mM, but the second runs to 5.41e-3 mM against a published 3.18e-3, a
1.70x overshoot that is still rising when the published window closes. Only the dynamics catch it. :data:`STOICHIOMETRY` is
audited against the deposit for exactly this reason.

WHAT THIS CLOSES, AND -- SAID FIRST -- WHAT IT DOES NOT. This closes the **cAMP/PKA input side
only**. The stress panel's declared carbon module is Snf1 -> Adr1/Cat8 -> CSRE, and Williamson
contains no Snf1, no Adr1, no Cat8, no promoter and no transcription of any kind: the model
terminates at free PKA catalytic subunit. Nothing found anywhere supplies that promoter-output
arm, so :data:`PROMOTER_COUPLING_NOT_BUILT` refuses it by name and adopting this module changes
how ``pka`` is COMPUTED without touching that standing refusal. Do not describe this as closing
the carbon module.

WHY IT IS STILL WORTH HAVING. ``mech/engine.py`` computes cAMP from ``camp_synthesis`` = 0.005
mmol/gDW/h times a glucose signal times an energy term (line 500) and destroys it at
``camp_hydrolysis`` = 10 /h times cAMP (line 501), with ``pka`` relaxing first-order at
``pka_tau`` = 0.03 h. **PKA appears in neither cAMP flux.** There is no negative feedback, so
that form is monotone in glucose and cannot produce a glucose-induced cAMP transient at all.
Williamson supplies the feedback with published constants: PKA phosphorylates Pde1, whose kcat
rises from 1.11 to 25.25 /s on phosphorylation, and phosphorylates Cdc25, throttling Ras2
activation. Under this transcription those two loops turn a 20-fold glucose step into a 1.29x
cAMP peak that decays -- the behaviour the engine's present form is structurally unable to show.

THE FOUR CONDITIONS THE BUILD WAS SET, AND WHAT EACH ONE TURNED OUT TO BE.

(a) **THE SHARED LOADER REFUSES THIS FILE, AND STILL DOES.** ``mech/kinetic_sbml.py`` line 332
    raises ``UnsupportedSBMLError("rate and algebraic rules are unsupported; only parameter
    assignment rules are supported")`` on the seven ``rateRule``s. This was verified live, not
    inferred. Nothing was stripped and the loader was not weakened. Instead :func:`audit_deposit`
    is an adapter that reads all seven, checks that each one's MathML is the literal ``0`` and
    that its target species is already flagged ``boundaryCondition="true"``, and only then
    interprets them as what they say: those seven species are constant. If a rule is ever
    anything other than a literal zero on a boundary species, the audit raises and this module
    does not run. See :data:`RATE_RULE_AUDIT_RULE`.

(b) **THE CITED SOURCE CANNOT ARBITRATE, BECAUSE IT IS A SOFTWARE PAPER.** Table 4 gives
    glucose transport KM 0.08 mM, V 1.7 mmol/s and Ki 0.91, all three attributed to reference
    **[37]**. The deposit gives ``GlucTrans_Km_19`` 1.7, ``GlucTrans_V_19`` 0.062 and
    ``GlucTrans_Ki_19`` 0.91. Reference [37] in this article is *Schmidt H & Jirstrand M 2006,
    "Systems Biology Toolbox for MATLAB", Bioinformatics 22:514* -- the estimation package the
    authors used. It contains no transport kinetics, no Km and no Vmax. **There is therefore no
    cited measurement to resolve the disagreement against**, and :data:`TRANSPORT_CITED_SOURCE`
    is REFUSED with that said out loud rather than a number being chosen by preference.
    What DID arbitrate is the source's own printed output, which requires inventing nothing: run
    both readings against all 36 figure landmarks. The deposit reading gives mean 0.30% and worst
    1.61% of full scale; the Table 4 reading gives mean 1.48% and worst 6.92%, and misses the
    pre-second-pulse Gpa2 state by 5-7% of full scale on four separate landmarks. The deposit
    reading is used, the Table 4 reading stays runnable as :data:`TABLE_4_TRANSPORT`, and
    :func:`reproduction_report` reports both. Note the disagreement is not a clean transposition:
    the file's Km, 1.7, is the paper's V, but the file's V, 0.062, is not the paper's Km, 0.08.

(c) **THREE NEGATIVE INITIAL CONCENTRATIONS, CLIPPED AND DECLARED.** ``Gpr1Glucout``
    (-2.393e-33), ``Gpa2a`` (-4.997e-33) and ``Gpa2aKrh`` (-1.547e-32) are negative in the
    deposit -- solver residue from the steady-state solve that produced the initial state, on
    species whose true starting value at zero glucose is zero. They are clipped to 0. Two further
    species, ``Glucout`` (8.446e-32) and ``Glucin`` (4.710e-32), are positive numerical zeros and
    are left exactly as deposited. The cost of the clip is measured, not assumed: over the
    Figure 8 protocol it moves cAMP by at most 1.3e-12 mM, a relative 7.4e-10, which is below
    the integrator's own tolerance. :data:`NEGATIVE_INITIAL_CONCENTRATIONS` records all three.

(d) **SECONDS.** The SBML declares no unit definitions, so SBML Level 2's default time unit,
    the second, applies; Table 4's rate units are ``l/s`` and ``mmol/s`` and both figure x-axes
    are labelled "time (s)". :data:`SECONDS_PER_HOUR` performs the 3600x that
    ``signalling.NativeHogResponse`` already performs for the HOG source, and every public entry
    point here takes hours and converts once, at the boundary.

FIVE DEFECTS IN THE SOURCE, FOUND BY TRANSCRIBING IT. None is repaired here.

* **Figure 9C's caption has its two colours the wrong way round.** It reads "Levels of active
  (blue trace) and inactive (green trace) PKA". The plotted green trace starts at 0.20e-4 mM,
  peaks at 1.47e-4 and ends at 1.04e-4; the blue starts at 0.77e-4, bottoms at 0.145e-4 and ends
  at 0.35e-4. The deposit's own initial state has C = 0.2047e-4 and PKAi = 0.7651e-4, and
  conservation (PKAi + C/2 is exactly 8.675e-5 here for the whole run) caps PKAi at 0.8675e-4 --
  so the green trace, which reaches 1.47e-4, cannot be inactive PKA. Green is the free catalytic
  subunit C and blue is PKAi, the opposite of the caption. This transcription reproduces both
  traces to better than 0.6% of full scale under that reading and cannot under the printed one.
* **Figure 8's caption cites the wrong reference for its own data.** The caption credits the
  dotted trace to "[35]", which is Mendes 1993, the GEPASI software paper. The Methods credit the
  cAMP time course to [38], Rolland et al. 2000 (Mol Microbiol 38:348), which is the paper that
  actually performed the 5 mM / 100 mM glucose pulse. [37] is misused the same way in Table 4.
* **The Complete cAMP Model's initial concentrations are published nowhere but the SBML.** No
  table in the article carries them; Tables 2 and 3 belong to the simplified models. The protein
  totals they imply -- Gpa2 0.1126 uM, Ras2 0.4826 uM, Cdc25 18 nM, Pde1 34.5 nM, Cyr1 10 nM --
  have no citation and no measurement behind them anywhere in the paper. They are ASSERTED here,
  every one, and :data:`ABSOLUTE_ABUNDANCE_SCALE` refuses the check that would ground them.
* **Table 1 does not describe the model Table 4 and the SBML describe.** Table 1's row for the
  Complete cAMP Model reads "27 parameters, 15 variables". Table 4 lists 39 parameters and the
  deposit declares 39 kinetic constants, 17 dynamic species and 7 held ones. Whatever Table 1
  counted, it was not this deposit, and every count here is taken from the file that runs.
* **``Sdc25`` is a dangling species.** It is declared, given 1e-5, flagged boundary and held by a
  rateRule, and appears in no reaction and no rule. Sdc25p is the second Ras GEF; the model
  declares it and never uses it. It is carried, held constant, and named as unused.

THE FREE-SCALAR GATE REFUSES THIS PIECE, AND SHOULD. Sixty-three declared constants -- 39 rate
parameters and 24 initial concentrations -- against **zero independent targets**. Thirty-five of
Table 4's 39 rows record their own source as "This work", fitted by simulated annealing to one
cAMP time course; three cite a software paper and one cites Sass 1986. The only
real measurement anywhere near this model is that same time course, so it is registered
``fitted=True`` and does not count. Figure 9's four panels are the fitted model's own simulation
output, not measurements: reproducing them is evidence about the transcription and about nothing
else. :func:`gate` returns that verdict instead of hiding it.
"""
from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from scipy.integrate import solve_ivp

from .. import paths
from .contracts import ScientificRefusal
from .params import Param, ParamRegistry, Target, provenance_table

__all__ = [
    "ABSOLUTE_ABUNDANCE_SCALE",
    "AVAILABILITY",
    "ARTIFACTS",
    "BOUNDARY_SPECIES",
    "CONSERVED_MOIETIES",
    "DEPOSIT_TRANSPORT",
    "DYNAMIC_SPECIES",
    "GLUCOSE_PULSE_PROTOCOL",
    "INITIAL_CONCENTRATIONS",
    "LICENCE",
    "NEGATIVE_INITIAL_CONCENTRATIONS",
    "PROMOTER_COUPLING_NOT_BUILT",
    "PUBLISHED_LANDMARKS",
    "RATE_LAWS",
    "RATE_RULE_AUDIT_RULE",
    "SBML_LOADER_REFUSAL",
    "SECONDS_PER_HOUR",
    "SOURCE",
    "STOICHIOMETRY",
    "TABLE_4_TRANSPORT",
    "TRANSPORT_CITED_SOURCE",
    "TransportReading",
    "UNUSED_SPECIES",
    "WILLIAMSON_PARAMS",
    "WilliamsonCAMP",
    "audit_deposit",
    "cell_volume_litres",
    "gate",
    "glucose_pulse_timecourse",
    "landmark_table",
    "provenance",
    "reproduction_report",
]

AVAILABILITY = "transcribed_here"
"""The level this module carries, and the reason it is not ``local_verified``.

``local_verified`` in this repository means a checksummed publisher file loaded unchanged by
``mech/kinetic_sbml.py``. That loader refuses this file (:data:`SBML_LOADER_REFUSAL`), so the
executable below is ours and inherits our bugs. Nothing here may be promoted past this level.
"""

SOURCE = ("Williamson T, Schwartz J-M, Kell DB, Stateva L 2009, 'Deterministic mathematical "
          "models of the cAMP pathway in Saccharomyces cerevisiae', BMC Syst Biol 3:70, "
          "doi:10.1186/1752-0509-3-70, PMID 19607691, PMC2719611")

LICENCE = ("CC BY 2.0 (https://creativecommons.org/licenses/by/2.0/), declared in the article's "
           "own <permissions> block: 'Copyright (c) 2009 Williamson et al; licensee BioMed "
           "Central Ltd. This is an Open Access article distributed under the terms of the "
           "Creative Commons Attribution License ... provided the original work is properly "
           "cited.' Redistribution of Additional file 1 is therefore permitted, with attribution")

ARTIFACTS = {
    "sbml_additional_file_1":
        ("data/native_reference_models/williamson2009/1752-0509-3-70-S1.xml",
         "e608d6651dd3332c0a2f146bdd9402c77ee07f7c6358284bef911da5bcd1b782"),
    "article_xml":
        ("data/native_reference_models/williamson2009/PMC2719611.xml",
         "73769eb3c3cb45b415de3ede58dbad3082191cf90788313cb3cc7ad28a13e58d"),
}
"""The two upstream files this transcription was made from, by content."""

SBML_LOADER_REFUSAL = (
    "mech/kinetic_sbml.py:332 -- UnsupportedSBMLError('rate and algebraic rules are "
    "unsupported; only parameter assignment rules are supported'). Verified live against this "
    "file. Seven rateRules trigger it; audit_deposit() reads them rather than stripping them")

NATIVE_TIME_UNIT = "second"
"""The source's time base, kept internally. SBML Level 2's default time unit, corroborated by
Table 4's l/s and mmol/s and by both figure x-axes reading 'time (s)'."""

SECONDS_PER_HOUR = 3600.0
"""The conversion every public entry point applies once, at the boundary. Three published time
bases collide in this architecture and the resulting 60x has already cost errors."""

_SBML_NS = "http://www.sbml.org/sbml/level2"
_MATHML_NS = "http://www.w3.org/1998/Math/MathML"

_TABLE_4 = "Williamson 2009 Table 4, 'Parameters of the complete cAMP pathway Model'"
_FITTED = ("fitted by simulated annealing (SBToolbox SBparameterEstimation) to ONE cAMP time "
           "course, Rolland et al. 2000 PMID 11069660; Table 4 records its source as 'This work'")

WILLIAMSON_PARAMS = ParamRegistry("mech/carbon_williamson2009.py -- Gpr1/Gpa2-Ras/cAMP/PKA")


def _fitted(name: str, value: float, units: str, sbml_id: str, table_4_row: str) -> Param:
    """One of Williamson's 37 estimated rate constants: registered ASSERTED, never MEASURED."""
    return WILLIAMSON_PARAMS.add(Param.asserted(
        name=f"williamson.{name}", value=value, units=units,
        source=f"{_TABLE_4}, '{table_4_row}' row; SBML global parameter {sbml_id}. {_FITTED}",
        missing="an independent kinetic determination of this step in S. cerevisiae"))


# ----------------------------------------------------------------------------------------
# The 41 global parameters of Additional file 1, in the order the file declares them.
# Table 4 rounds several of these; the deposit's digits are used because the deposit ran.
# ----------------------------------------------------------------------------------------

GPR1_GLUCOSE_ASSOCIATION = _fitted(
    "gpr1_glucose_association_k1", 0.00301659503091382, "L^2 mmol^-1 s^-1",
    "Gpr1Glucass_k_1", "Gpr1 Glucose association k1 = 0.003")
GPR1_GLUCOSE_DISSOCIATION = _fitted(
    "gpr1_glucose_dissociation_k2", 0.14284392829434, "s^-1",
    "Gpr1Glucdiss_k_2", "Gpr1 Glucose dissociation k1 = 0.14")
GPA2_ACTIVATION_KA = _fitted(
    "gpa2_activation_kA", 57682.6153976306, "L^2 mmol^-1 s^-1",
    "Gpa2act_kA_3", "Gpa2 activation kA = 57682.6")
GPA2_DEACTIVATION_KA = _fitted(
    "gpa2_deactivation_kA", 12989.4420113722, "L^2 mmol^-1 s^-1",
    "Gpa2deact_kA_4", "Gpa2 deactivation kA = 12989.4")
GPA2_DEACTIVATION_K1 = _fitted(
    "gpa2_deactivation_k1", 0.899003720389271, "s^-1",
    "Gpa2deact_k1_4", "Gpa2 deactivation kF = 0.899")
GPA2_KRH_ASSOCIATION = _fitted(
    "gpa2_krh_association_k5", 391089.569083089, "L^2 mmol^-1 s^-1",
    "Gpa2Krhass_k_5", "Gpa2-Krh association kF = 391089.6")
GPA2_KRH_DISSOCIATION = _fitted(
    "gpa2_krh_dissociation_k6", 6.12210079340681, "s^-1",
    "Gpa2Krhdiss_k_6", "Gpa2-Krh dissociation kF = 6.12")
RAS2_ACTIVATION_KCAT = _fitted(
    "ras2_activation_kcat", 0.740604088325018, "s^-1",
    "Ras2act_kcat_7", "Ras2 activation kcat = 0.74")
RAS2_ACTIVATION_KM = _fitted(
    "ras2_activation_Km", 0.00138648047637647, "mmol L^-1",
    "Ras2act_Km_7", "Ras2 activation KM = 1.38e-3")
RAS2_ACTIVATION_KD = _fitted(
    "ras2_activation_Kd", 0.0440833580413153, "mmol L^-1",
    "Ras2act_Kd_7", "Ras2 activation Kd = 0.044")
RAS2_ACTIVATION_A = _fitted(
    "ras2_activation_a", 32.9031890215006, "dimensionless",
    "Ras2act_a_7", "Ras2 activation a = 32.9")
RAS2_ACTIVATION_B = _fitted(
    "ras2_activation_b", 63.7501519114594, "dimensionless",
    "Ras2act_b_7", "Ras2 activation b = 63.8")
RAS2_DEACTIVATION_KA = _fitted(
    "ras2_deactivation_kA", 519.811493442033, "L^2 mmol^-1 s^-1",
    "Ras2deact_kA_8", "Ras2 deactivation kA = 519.8")
RAS2_DEACTIVATION_K1 = _fitted(
    "ras2_deactivation_k1", 0.042059375613429, "s^-1",
    "Ras2deact_k1_8", "Ras2 deactivation kF = 0.042")
CDC25_PHOSPHORYLATION_KCAT = _fitted(
    "cdc25_phosphorylation_kcat", 0.185477275623681, "s^-1",
    "Cdc25phos_kcat_9", "Cdc25 phosphorylation kcat = 0.18")
CDC25_PHOSPHORYLATION_KM = _fitted(
    "cdc25_phosphorylation_Km", 0.00516736255009312, "mmol L^-1",
    "Cdc25phos_Km_9", "Cdc25 phosphorylation KM = 5.2e-3")
CDC25_DEPHOSPHORYLATION_KCAT = _fitted(
    "cdc25_dephosphorylation_kcat", 2.51636177378834, "s^-1",
    "Cdc25dephos_kcat_10", "Cdc25 dephosphorylation kcat = 2.52")
CDC25_DEPHOSPHORYLATION_KM = _fitted(
    "cdc25_dephosphorylation_Km", 0.0162516819523974, "mmol L^-1",
    "Cdc25dephos_Km_10", "Cdc25 dephosphorylation KM = 1.6e-2")
CAMP_SYNTHESIS_KCAT_GPA2 = _fitted(
    "camp_synthesis_kcat_gpa2", 2933.87233320588, "dimensionless",
    "cAMPsynth_kcatGpa2_11", "cAMP synthesis kcatGpa2 = 2933.9")
CAMP_SYNTHESIS_KCAT_RAS2 = _fitted(
    "camp_synthesis_kcat_ras2", 650.0, "dimensionless",
    "cAMPsynth_kcatRas2_11", "cAMP synthesis kcatRas2 = 650")
CAMP_SYNTHESIS_KM = _fitted(
    "camp_synthesis_Km", 0.004, "mmol L^-1",
    "cAMPsynth_Km_11", "cAMP synthesis KM = 4e-3")

CAMP_HYDROLYSIS_PDE2_KCAT = _fitted(
    "camp_hydrolysis_pde2_kcat", 1.0, "s^-1",
    "cAMPhydroPde2_kcat_12", "cAMP hydrolysis (Pde2) kcat = 1")
CAMP_HYDROLYSIS_PDE2_KM = WILLIAMSON_PARAMS.add(Param.asserted(
    name="williamson.camp_hydrolysis_pde2_Km", value=0.002, units="mmol L^-1",
    source=f"{_TABLE_4}, 'cAMP hydrolysis (Pde2) KM = 0.002' row; SBML cAMPhydroPde2_Km_12. "
           "The ONLY row of Table 4 attributed to a kinetic measurement rather than to this "
           "work: reference [15], Sass P et al. 1986 PNAS 83:9303, PMID 3025832, the PDE2 "
           "high-affinity phosphodiesterase. THIS TRANSCRIPTION DID NOT OPEN SASS 1986; the "
           "attribution is Williamson's and is recorded, not verified, so the grade stays "
           "ASSERTED rather than being promoted on a citation nobody here checked",
    missing="reading Sass 1986 and confirming a 2 uM Km for Pde2p in matched conditions"))

CAMP_HYDROLYSIS_PDE1P_KCAT = _fitted(
    "camp_hydrolysis_pde1p_kcat", 25.2539854673553, "s^-1",
    "cAMPhydroPde1p_kcat_13", "cAMP hydrolysis (Pde1p) kcat = 25.25")
CAMP_HYDROLYSIS_PDE1P_KM = _fitted(
    "camp_hydrolysis_pde1p_Km", 6.08741428019625e-07, "mmol L^-1",
    "cAMPhydroPde1p_Km_13", "cAMP hydrolysis (Pde1p) KM = 6e-7")
CAMP_HYDROLYSIS_PDE1_KCAT = _fitted(
    "camp_hydrolysis_pde1_kcat", 1.11283526581817, "s^-1",
    "cAMPhydroPde1_kcat_14", "cAMP hydrolysis (Pde1) kcat = 1.1")
CAMP_HYDROLYSIS_PDE1_KM = _fitted(
    "camp_hydrolysis_pde1_Km", 0.845082028313041, "mmol L^-1",
    "cAMPhydroPde1_Km_14", "cAMP hydrolysis (Pde1) KM = 0.85")
PDE1_PHOSPHORYLATION_KCAT = _fitted(
    "pde1_phosphorylation_kcat", 6.8221366652867, "s^-1",
    "Pde1phos_kcat_15", "Pde1 phosphorylation kcat = 6.82")
PDE1_PHOSPHORYLATION_KM = _fitted(
    "pde1_phosphorylation_Km", 0.00864811688536481, "mmol L^-1",
    "Pde1phos_Km_15", "Pde1 phosphorylation KM = 8.6e-3")
PDE1_DEPHOSPHORYLATION_KCAT = _fitted(
    "pde1_dephosphorylation_kcat", 2.41696125657814, "s^-1",
    "Pde1dephos_kcat_16", "Pde1 dephosphorylation kcat = 2.4")
PDE1_DEPHOSPHORYLATION_KM = _fitted(
    "pde1_dephosphorylation_Km", 0.00107471629124917, "mmol L^-1",
    "Pde1dephos_Km_16", "Pde1 dephosphorylation KM = 1.07e-3")
PKA_ACTIVATION_KF = _fitted(
    "pka_activation_kF", 761676017.429656, "mmol^-4 L^4 s^-1",
    "PKAact_kF_17", "PKA activation kF = 7.6e8")
PKA_ACTIVATION_KI = _fitted(
    "pka_activation_kI", 100.0, "L mmol^-1",
    "PKAact_kI_17", "PKA activation kI = 100")
PKA_DEACTIVATION_KF = _fitted(
    "pka_deactivation_kF", 19.9369792588504, "mmol^-1 L s^-1",
    "PKAdeact_kF_18", "PKA deactivation kF = 19.9")
PKA_DEACTIVATION_KA = _fitted(
    "pka_deactivation_kA", 22139.0512490006, "L mmol^-1",
    "PKAdeact_kA_18", "PKA deactivation kA = 2.2e4")

_TRANSPORT_DISAGREEMENT = (
    "TABLE 4 AND THE DEPOSIT DISAGREE. Table 4 prints 'Glucose transport KM 0.08 mM' and "
    "'Glucose transport V 1.7 mmol/s'; the SBML has GlucTrans_Km_19 = 1.7 and GlucTrans_V_19 = "
    "0.062. Not a clean transposition: the file's Km IS the paper's V, but the file's V is not "
    "the paper's Km. Arbitrated against the source's own Figures 8 and 9, not by preference -- "
    "see TRANSPORT_CITED_SOURCE and reproduction_report()['transport_arbitration']")

GLUCOSE_TRANSPORT_V = WILLIAMSON_PARAMS.add(Param.asserted(
    name="williamson.glucose_transport_V", value=0.062, units="mmol s^-1",
    source=f"SBML GlucTrans_V_19 = 0.062; {_TABLE_4} prints 1.7. {_TRANSPORT_DISAGREEMENT}",
    missing="a glucose-transport Vmax measured in the strain and medium Williamson modelled"))
GLUCOSE_TRANSPORT_KM = WILLIAMSON_PARAMS.add(Param.asserted(
    name="williamson.glucose_transport_Km", value=1.7, units="mmol L^-1",
    source=f"SBML GlucTrans_Km_19 = 1.7; {_TABLE_4} prints 0.08. {_TRANSPORT_DISAGREEMENT}",
    missing="a glucose-transport KM measured in the strain and medium Williamson modelled"))
GLUCOSE_TRANSPORT_KI = WILLIAMSON_PARAMS.add(Param.asserted(
    name="williamson.glucose_transport_Ki", value=0.91, units="dimensionless",
    source=f"SBML GlucTrans_Ki_19 = 0.91, and {_TABLE_4} prints 0.91 -- the ONE transport "
           "constant the paper and the deposit agree on. Table 4 attributes it to [37], which "
           "supplies no kinetics; see TRANSPORT_CITED_SOURCE",
    missing="the transport measurement Table 4's reference [37] does not contain"))
GLUCOSE_UTILISATION_K = _fitted(
    "glucose_utilisation_k", 0.03, "s^-1",
    "GlucUtil_k_20", "Glucose Utilisation kF = 0.03")

DEPOSIT_TRANSPORT = (0.062, 1.7, 0.91)
"""``(V, Km, Ki)`` as deposited. The reading this module runs; see :class:`TransportReading`."""

TABLE_4_TRANSPORT = (1.7, 0.08, 0.91)
"""``(V, Km, Ki)`` as Table 4 prints them. Kept runnable so the alternative is not deleted."""


# ----------------------------------------------------------------------------------------
# The 24 initial concentrations. Published NOWHERE but the SBML -- no table in the article
# carries them, and the protein totals they imply have no citation anywhere in the paper.
# ----------------------------------------------------------------------------------------

_IC_SOURCE = ("SBML Additional file 1 listOfSpecies initialConcentration. NOT PRINTED IN THE "
              "ARTICLE: Tables 2 and 3 belong to the simplified models and no table gives the "
              "Complete cAMP Model's initial state. The deposit's initial vector IS a steady "
              "state of its own rate laws (verified here to 3.5e-14 /s relative), but the "
              "abundances themselves carry no citation and no measurement")

_INITIAL_RAW = {
    "Rgs2": 1e-05, "Cyr1": 1e-05, "Pde2": 0.00016, "Ira": 1.1e-05, "Mih1": 1e-05,
    "PPA": 0.000224, "Sdc25": 1e-05,
    "Gpr1": 9.99999999999999e-06, "Glucout": 8.44618033564067e-32,
    "Gpr1Glucout": -2.39301410052969e-33, "Gpa2i": 0.0001126,
    "Gpa2a": -4.99721839585425e-33, "Krh": 0.0001, "Gpa2aKrh": -1.54688941759641e-32,
    "Ras2i": 0.000437206156106628, "Ras2a": 4.53938438933718e-05,
    "Cdc25": 1.2215083012251e-05, "Cdc25P": 5.78491698774897e-06,
    "cAMP": 0.000825917614307686, "Pde1": 3.34314121580972e-05,
    "PKAi": 7.65148419569979e-05, "C": 2.04703160860042e-05,
    "Glucin": 4.71036809955077e-32, "Pde1P": 1.06858784190283e-06,
}
"""Exactly as deposited, negatives included. :data:`INITIAL_CONCENTRATIONS` is the clipped one."""

NEGATIVE_INITIAL_CONCENTRATIONS = MappingProxyType({
    "Gpr1Glucout": -2.39301410052969e-33,
    "Gpa2a": -4.99721839585425e-33,
    "Gpa2aKrh": -1.54688941759641e-32,
})
"""Condition (c). Solver residue on three species whose value at zero glucose is zero; clipped
to 0.0, and the clip costs at most 7.4e-10 relative in cAMP over the Figure 8 protocol."""

NUMERICAL_ZERO_INITIAL_CONCENTRATIONS = MappingProxyType({
    "Glucout": 8.44618033564067e-32, "Glucin": 4.71036809955077e-32,
})
"""Positive numerical zeros, left exactly as deposited. Recorded so the clip cannot creep."""

INITIAL_CONCENTRATIONS = MappingProxyType(
    {name: (0.0 if value < 0.0 else value) for name, value in _INITIAL_RAW.items()})
"""The deposited initial state with the three negatives clipped to zero. mmol/L."""

for _name, _value in _INITIAL_RAW.items():
    # One registry row per initial concentration, so the gate counts all 24 rather than
    # hiding them inside a dict literal. Clipping is noted on the three that need it.
    _clip = (" CLIPPED TO 0.0 HERE: deposited negative, see NEGATIVE_INITIAL_CONCENTRATIONS"
             if _value < 0.0 else "")
    WILLIAMSON_PARAMS.add(Param.asserted(
        name=f"williamson.initial.{_name}", value=INITIAL_CONCENTRATIONS[_name],
        units="mmol L^-1", source=f"{_IC_SOURCE}; species {_name} = {_value!r}.{_clip}",
        missing=f"an absolute {_name} quantification in the modelled strain and condition"))
del _name, _value, _clip


# ----------------------------------------------------------------------------------------
# What the source does not contain. These raise on use rather than defaulting.
# ----------------------------------------------------------------------------------------

NOT_IN_SOURCE = ParamRegistry("mech/carbon_williamson2009.py -- refused, absent from Williamson")

PROMOTER_COUPLING_NOT_BUILT = NOT_IN_SOURCE.add(Param.refused(
    name="williamson.snf1_adr1_cat8_csre_gain", units="CSRE output per unit PKA",
    source="The strings Snf1, Adr1, Cat8, Mig1, CSRE, Msn2 and even 'transcription' occur "
           "ZERO times in either pinned artifact -- verified by grep against the article "
           f"full text and the deposited SBML, both listed in ARTIFACTS. {SOURCE}",
    reason="The stress panel's declared carbon module is Snf1 -> Adr1/Cat8 -> CSRE and this "
           "model terminates at the free PKA catalytic subunit: no promoter, no transcript, no "
           "Snf1 at all. Adopting Williamson changes how pka is COMPUTED and closes the "
           "cAMP/PKA INPUT side only. The standing PROMOTER_COUPLING_NOT_BUILT refusal on the "
           "output arm is untouched, and this module must never be described as closing the "
           "carbon module",
    missing="a parameterised dynamic model of Snf1 -> Adr1/Cat8 -> CSRE promoter output in "
            "S. cerevisiae, which this programme's survey did not find in any repository"))

TRANSPORT_CITED_SOURCE = NOT_IN_SOURCE.add(Param.refused(
    name="williamson.glucose_transport_cited_measurement", units="mmol L^-1 and mmol s^-1",
    source="Table 4 attributes glucose transport KM, V and Ki to reference [37]. Reference [37] "
           "of this article is Schmidt H & Jirstrand M 2006, 'Systems Biology Toolbox for "
           "MATLAB: a computational platform for research in systems biology', Bioinformatics "
           "22:514, PMID 16317076 -- the parameter-estimation package the authors used",
    reason="The citation is broken. A software paper contains no transport kinetics, so the "
           "disagreement between Table 4 (KM 0.08, V 1.7) and the deposit (Km 1.7, V 0.062) has "
           "NO cited measurement to be resolved against, and no number may be chosen by "
           "preference or by averaging. What arbitrates here instead is the source's own "
           "printed Figures 8 and 9, which requires inventing nothing -- but a figure the model "
           "was fitted to is not the measurement this row is refusing",
    missing="a glucose-transport KM and Vmax measured in S. cerevisiae in the strain, carbon "
            "source and cell density Williamson's Complete cAMP Model represents"))

ABSOLUTE_ABUNDANCE_SCALE = NOT_IN_SOURCE.add(Param.refused(
    name="williamson.absolute_abundance_scale", units="molecules per cell per mmol L^-1",
    source="No table in Williamson 2009 carries the Complete cAMP Model's initial "
           "concentrations, and the protein totals they set -- Gpa2 0.1126 uM, Ras2 0.4826 uM, "
           "Cdc25 18 nM, Pde1 34.5 nM, Cyr1 10 nM, Krh 0.1 uM -- are cited to nothing",
    reason="Without an absolute anchor the deposited pools are a self-consistent scale, not "
           "abundances: PKA activation kF = 7.6e8 mmol^-4 L^4 s^-1 is a fitted lumped constant "
           "on that scale and is not a physical fourth-order rate. Reading any of these as "
           "measured protein levels is the mistake this row exists to prevent downstream",
    missing="matched absolute quantification of Gpa2, Ras2, Cdc25, Pde1, Pde2, Cyr1 and Krh in "
            "one sample, in the glucose-starved condition the deposited state stands for"))

CELL_VOLUME = NOT_IN_SOURCE.add(Param.refused(
    name="williamson.cell_volume_litres", units="L",
    source="SBML listOfCompartments: a single compartment of size 1 with no unit declared, and "
           "no volume anywhere in the article. Table 4 nonetheless gives V in mmol/s and three "
           "rate constants in L^2/mmol/s, which only balance against a volume the model omits",
    reason="Converting this model's mM to the engine's mmol/gDW needs a cell water volume per "
           "gram dry weight that Williamson never states. The repository's own "
           "cell_water_l_gdw is an asserted prior, so routing through it would launder one "
           "assertion through another. The dimensionless PKA activation fraction is the "
           "coupling this module offers precisely because it needs no volume",
    missing="the cell water volume per gDW of the culture Williamson's model represents, in "
            "the same condition as its deposited steady state"))

STRESS_TO_PKA_DEMAND = NOT_IN_SOURCE.add(Param.refused(
    name="williamson.pka_to_growth_or_burden", units="per-hour growth per unit active PKA",
    source="Williamson 2009 contains no growth rate, no biomass, no protein synthesis and no "
           "heterologous product. Its two PKA substrates are Pde1 and Cdc25, both internal to "
           "the cAMP loop it closes",
    reason="This model supplies the INPUT to pka and nothing downstream of it. A coefficient "
           "turning active PKA into growth, ribosome content, burden or product would be "
           "invented here, and inventing it is the one thing this build was forbidden",
    missing="a measurement in S. cerevisiae relating PKA catalytic-subunit activity to growth "
            "rate or to recombinant protein output at fixed carbon uptake"))

UNUSED_SPECIES = ("Sdc25",)
"""Declared, valued at 1e-5, flagged boundary, held by a rateRule -- and in no reaction and no
rule. Sdc25p is the second Ras GEF; the model declares it and never uses it."""


# ----------------------------------------------------------------------------------------
# The transcribed model: species, rate laws, stoichiometry. All audited against the deposit.
# ----------------------------------------------------------------------------------------

BOUNDARY_SPECIES = ("Rgs2", "Cyr1", "Pde2", "Ira", "Mih1", "PPA", "Sdc25")
"""The seven species held constant. Each is ``boundaryCondition="true"`` AND carries a
``rateRule`` of literal 0 -- the two together are the whole content of condition (a)."""

DYNAMIC_SPECIES = (
    "Gpr1", "Glucout", "Gpr1Glucout", "Gpa2i", "Gpa2a", "Krh", "Gpa2aKrh", "Ras2i", "Ras2a",
    "Cdc25", "Cdc25P", "cAMP", "Pde1", "PKAi", "C", "Glucin", "Pde1P",
)
"""The seventeen integrated states, in the deposit's own declaration order."""

RATE_RULE_AUDIT_RULE = (
    "every rateRule must (1) target a species, (2) have MathML that is the single literal 0, "
    "and (3) target a species already flagged boundaryCondition='true'. Only then is it what "
    "it says -- a held constant. Anything else raises and this module does not run")


@dataclass(frozen=True)
class TranscribedReaction:
    """One reaction of Table 5, its deposited rate law, and the stoichiometry it carries.

    ``printed`` is Table 5's rate law in plain ASCII where Table 5 typesets one, and the
    deposit's MathML rendered infix where Table 5 shows only a graphic. ``stoichiometry`` maps
    species to signed coefficients and is checked against the SBML by :func:`audit_deposit`.
    """

    sbml_id: str
    table_5_row: str
    printed: str
    stoichiometry: tuple[tuple[str, int], ...]

    def __post_init__(self) -> None:
        if not self.stoichiometry:
            raise ValueError(f"{self.sbml_id}: a reaction moves at least one species")


RATE_LAWS = (
    TranscribedReaction(
        "Gpr1Glucass", "Gpr1-Glucose Association", "kf*[Gpr1]*[Glucout]",
        (("Gpr1", -1), ("Glucout", -1), ("Gpr1Glucout", 1))),
    TranscribedReaction(
        "Gpr1Glucdiss", "Gpr1-Glucose dissociation", "kf*[Gpr1Glucout]",
        (("Gpr1Glucout", -1), ("Gpr1", 1), ("Glucout", 1))),
    TranscribedReaction(
        "Gpa2act", "Gpa2 activation", "kA*[Gpr1Glucout]*[Gpa2i]",
        (("Gpa2i", -1), ("Gpa2a", 1))),
    TranscribedReaction(
        "Gpa2deact", "Gpa2 deactivation", "(kA*[Rgs2] + k1)*[Gpa2a]",
        (("Gpa2a", -1), ("Gpa2i", 1))),
    TranscribedReaction(
        "Gpa2Krhass", "Gpa2-Krh association", "kf*[Gpa2a]*[Krh]",
        (("Gpa2a", -1), ("Krh", -1), ("Gpa2aKrh", 1))),
    TranscribedReaction(
        "Gpa2Krhdiss", "Gpa2-Krh dissociation", "kf*[Gpa2aKrh]",
        (("Gpa2aKrh", -1), ("Gpa2a", 1), ("Krh", 1))),
    TranscribedReaction(
        "Ras2act", "Ras2 activation",
        "(kcat*[Cdc25]*[Ras2i]/Km)*(1 + b*[Glucin]/(Kd*a)) / "
        "(1 + [Glucin]/Kd + ([Ras2i]/Km)*(1 + [Glucin]/(Kd*a)))",
        (("Ras2i", -1), ("Ras2a", 1))),
    TranscribedReaction(
        "Ras2deact", "Ras2 deactivation", "(kA*[Ira] + k1)*[Ras2a]",
        (("Ras2a", -1), ("Ras2i", 1))),
    TranscribedReaction(
        "Cdc25phos", "Cdc25 phosphorylation", "[C]*[Cdc25]*kcat/(Km + [Cdc25])",
        (("Cdc25", -1), ("Cdc25P", 1))),
    TranscribedReaction(
        "Cdc25dephos", "Cdc25 dephosphorylation", "[Mih1]*[Cdc25P]*kcat/(Km + [Cdc25P])",
        (("Cdc25P", -1), ("Cdc25", 1))),
    TranscribedReaction(
        "cAMPsynth", "cAMP synthesis",
        "[Cyr1]*(kcatGpa2*[Gpa2a] + kcatRas2*[Ras2a])/Km",
        (("cAMP", 1),)),
    TranscribedReaction(
        "cAMPhydroPde2", "cAMP hydrolysis (Pde2)", "[Pde2]*[cAMP]*kcat/(Km + [cAMP])",
        (("cAMP", -1),)),
    TranscribedReaction(
        "cAMPhydroPde1p", "cAMP hydrolysis (Pde1P)", "[Pde1P]*[cAMP]*kcat/(Km + [cAMP])",
        (("cAMP", -1),)),
    TranscribedReaction(
        "cAMPhydroPde1", "cAMP hydrolysis (Pde1)", "[Pde1]*[cAMP]*kcat/(Km + [cAMP])",
        (("cAMP", -1),)),
    TranscribedReaction(
        "Pde1phos", "Pde1 phosphorylation", "[C]*[Pde1]*kcat/(Km + [Pde1])",
        (("Pde1", -1), ("Pde1P", 1))),
    TranscribedReaction(
        "Pde1dephos", "Pde1 dephosphorylation", "[PPA]*[Pde1P]*kcat/(Km + [Pde1P])",
        (("Pde1P", -1), ("Pde1", 1))),
    TranscribedReaction(
        "PKAact", "PKA activation", "[PKAi]*[cAMP]^4*kF/(1 + kI*[Krh])",
        (("PKAi", -1), ("C", 2))),
    TranscribedReaction(
        "PKAdeact", "PKA deactivation", "kF*[C]^2*(1 + kA*[Krh])",
        (("C", -2), ("PKAi", 1))),
    TranscribedReaction(
        "GlucTrans", "Glucose transport (reversible)",
        "V*([Glucout]/Km - [Glucin]/Km) / "
        "(1 + [Glucout]/Km + [Glucin]/Km + Ki*[Glucout]*[Glucin]/Km^2)",
        (("Glucout", -1), ("Glucin", 1))),
    TranscribedReaction(
        "GlucUtil", "Glucose metabolism", "kf*[Glucin]",
        (("Glucin", -1),)),
)

STOICHIOMETRY = MappingProxyType(
    {reaction.sbml_id: MappingProxyType(dict(reaction.stoichiometry))
     for reaction in RATE_LAWS})
"""Signed coefficients by reaction. ``PKAact`` releasing TWO catalytic subunits is the entry the
dynamics catch and the steady-state residual does not; see the module docstring."""

CONSERVED_MOIETIES = MappingProxyType({
    "gpr1_total": (("Gpr1", 1), ("Gpr1Glucout", 1)),
    "gpa2_total": (("Gpa2i", 1), ("Gpa2a", 1), ("Gpa2aKrh", 1)),
    "krh_total": (("Krh", 1), ("Gpa2aKrh", 1)),
    "ras2_total": (("Ras2i", 1), ("Ras2a", 1)),
    "cdc25_total": (("Cdc25", 1), ("Cdc25P", 1)),
    "pde1_total": (("Pde1", 1), ("Pde1P", 1)),
    "pka_total": (("PKAi", 1), ("C", 0.5)),
})
"""Seven moieties the stoichiometry conserves. ``pka_total`` weights C by 1/2 because one
holoenzyme releases two catalytic subunits; a 1:1 reading breaks this invariant."""

GLUCOSE_PULSE_PROTOCOL = ((60.0, 5.0), (240.0, 100.0))
"""``(time_s, Glucout_mM)``, from Williamson's own Methods: "5 mM glucose was added to
glucose-starved cell suspension after 60 seconds, followed by the addition of 100 mM glucose
after 240 seconds". Applied as state resets; the deposited SBML carries no events."""


# ----------------------------------------------------------------------------------------
# Condition (a): the audited adapter. Nothing is stripped and the shared loader is untouched.
# ----------------------------------------------------------------------------------------

def _digest(relative: str) -> str:
    path = paths.data_dir().parent / relative
    if not path.exists():
        raise FileNotFoundError(
            f"the transcription source {relative} is missing; this module is a transcription OF "
            "that file and will not run without it")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mathml_is_literal_zero(math_element) -> bool:
    """True only when the MathML body is a single ``<cn>`` whose value is exactly zero."""
    children = [child for child in math_element]
    if len(children) != 1 or children[0].tag != f"{{{_MATHML_NS}}}cn":
        return False
    return float(children[0].text.strip()) == 0.0


def audit_deposit() -> dict:
    """Re-read the vendored SBML and check this transcription against it, or refuse to run.

    Six checks, all of which must pass: the file's sha256; the seven rateRules against
    :data:`RATE_RULE_AUDIT_RULE`; the species set and every initial concentration; every one of
    the 41 global parameter values against the literals registered above; and every
    stoichiometric coefficient of all 20 reactions. Returns the audit record.

    This is the adapter condition (a) asks for. The rateRules are READ, not stripped: a rule
    that is a literal zero on a species already flagged ``boundaryCondition`` means that species
    is constant, which is exactly how :class:`WilliamsonCAMP` treats those seven. A rule that is
    anything else raises here.
    """
    relative, expected = ARTIFACTS["sbml_additional_file_1"]
    actual = _digest(relative)
    if actual != expected:
        raise ValueError(
            f"{relative} has sha256 {actual}, expected {expected}. This transcription was made "
            "from the expected bytes; a changed deposit means the transcription is stale")
    root = ET.parse(paths.data_dir().parent / relative).getroot()
    model = root.find(f"{{{_SBML_NS}}}model")

    species = {}
    for element in model.find(f"{{{_SBML_NS}}}listOfSpecies"):
        species[element.get("id")] = {
            "initial": float(element.get("initialConcentration")),
            "boundary": element.get("boundaryCondition") == "true",
        }
    expected_species = set(BOUNDARY_SPECIES) | set(DYNAMIC_SPECIES)
    if set(species) != expected_species:
        raise ValueError(
            f"the deposit declares species {sorted(species)}, this transcription declares "
            f"{sorted(expected_species)}")
    for name, deposited in species.items():
        if deposited["initial"] != _INITIAL_RAW[name]:
            raise ValueError(
                f"initial concentration of {name}: deposit {deposited['initial']!r}, "
                f"transcription {_INITIAL_RAW[name]!r}")
        if deposited["boundary"] != (name in BOUNDARY_SPECIES):
            raise ValueError(f"{name}: boundaryCondition disagrees with BOUNDARY_SPECIES")

    rate_rules, assignment_rules = {}, {}
    for rule in model.find(f"{{{_SBML_NS}}}listOfRules"):
        target = rule.get("variable")
        math = rule.find(f"{{{_MATHML_NS}}}math")
        if rule.tag == f"{{{_SBML_NS}}}rateRule":
            if not _mathml_is_literal_zero(math):
                raise ValueError(
                    f"rateRule on {target!r} is not the literal 0. {RATE_RULE_AUDIT_RULE}. It "
                    "must not be stripped to make the file load")
            if target not in species or not species[target]["boundary"]:
                raise ValueError(
                    f"rateRule on {target!r} targets something that is not a boundaryCondition "
                    f"species. {RATE_RULE_AUDIT_RULE}")
            rate_rules[target] = 0.0
        elif rule.tag == f"{{{_SBML_NS}}}assignmentRule":
            assignment_rules[target] = True
        else:
            raise ValueError(f"unexpected rule {rule.tag} on {target!r}")
    if set(rate_rules) != set(BOUNDARY_SPECIES):
        raise ValueError(
            f"the deposit holds {sorted(rate_rules)} by rateRule; BOUNDARY_SPECIES is "
            f"{sorted(BOUNDARY_SPECIES)}")

    deposited_parameters = {
        element.get("id"): float(element.get("value"))
        for element in model.find(f"{{{_SBML_NS}}}listOfParameters")}
    transcribed = _sbml_parameter_literals()
    for sbml_id, value in transcribed.items():
        if deposited_parameters[sbml_id] != value:
            raise ValueError(
                f"parameter {sbml_id}: deposit {deposited_parameters[sbml_id]!r}, transcription "
                f"{value!r}")
    unchecked = set(deposited_parameters) - set(transcribed) - set(assignment_rules)
    if unchecked:
        raise ValueError(f"deposited parameters not transcribed: {sorted(unchecked)}")

    for reaction in model.find(f"{{{_SBML_NS}}}listOfReactions"):
        deposited_stoichiometry: dict[str, int] = {}
        for side, sign in (("Reactants", -1), ("Products", 1)):
            listing = reaction.find(f"{{{_SBML_NS}}}listOf{side}")
            for reference in (listing if listing is not None else []):
                coefficient = int(float(reference.get("stoichiometry") or 1))
                name = reference.get("species")
                deposited_stoichiometry[name] = (
                    deposited_stoichiometry.get(name, 0) + sign * coefficient)
        mine = dict(STOICHIOMETRY[reaction.get("id")])
        if deposited_stoichiometry != mine:
            raise ValueError(
                f"reaction {reaction.get('id')}: deposit stoichiometry "
                f"{deposited_stoichiometry}, transcription {mine}")

    return {
        "sbml": {"path": relative, "sha256": actual},
        "article": {"path": ARTIFACTS["article_xml"][0],
                    "sha256": _digest(ARTIFACTS["article_xml"][0])},
        "species_checked": len(species),
        "parameters_checked": len(transcribed),
        "reactions_checked": len(STOICHIOMETRY),
        "rate_rules": {"count": len(rate_rules), "all_literal_zero": True,
                       "all_on_boundary_species": True, "targets": sorted(rate_rules),
                       "rule": RATE_RULE_AUDIT_RULE, "stripped": False},
        "assignment_rules": sorted(assignment_rules),
        "shared_loader": SBML_LOADER_REFUSAL,
        "unused_species": list(UNUSED_SPECIES),
    }


def _sbml_parameter_literals() -> dict[str, float]:
    """Every global parameter this transcription types, keyed by its SBML id."""
    return {
        "Gpr1Glucass_k_1": float(GPR1_GLUCOSE_ASSOCIATION),
        "Gpr1Glucdiss_k_2": float(GPR1_GLUCOSE_DISSOCIATION),
        "Gpa2act_kA_3": float(GPA2_ACTIVATION_KA),
        "Gpa2deact_kA_4": float(GPA2_DEACTIVATION_KA),
        "Gpa2deact_k1_4": float(GPA2_DEACTIVATION_K1),
        "Gpa2Krhass_k_5": float(GPA2_KRH_ASSOCIATION),
        "Gpa2Krhdiss_k_6": float(GPA2_KRH_DISSOCIATION),
        "Ras2act_kcat_7": float(RAS2_ACTIVATION_KCAT),
        "Ras2act_Km_7": float(RAS2_ACTIVATION_KM),
        "Ras2act_Kd_7": float(RAS2_ACTIVATION_KD),
        "Ras2act_a_7": float(RAS2_ACTIVATION_A),
        "Ras2act_b_7": float(RAS2_ACTIVATION_B),
        "Ras2deact_kA_8": float(RAS2_DEACTIVATION_KA),
        "Ras2deact_k1_8": float(RAS2_DEACTIVATION_K1),
        "Cdc25phos_kcat_9": float(CDC25_PHOSPHORYLATION_KCAT),
        "Cdc25phos_Km_9": float(CDC25_PHOSPHORYLATION_KM),
        "Cdc25dephos_kcat_10": float(CDC25_DEPHOSPHORYLATION_KCAT),
        "Cdc25dephos_Km_10": float(CDC25_DEPHOSPHORYLATION_KM),
        "cAMPsynth_kcatGpa2_11": float(CAMP_SYNTHESIS_KCAT_GPA2),
        "cAMPsynth_kcatRas2_11": float(CAMP_SYNTHESIS_KCAT_RAS2),
        "cAMPsynth_Km_11": float(CAMP_SYNTHESIS_KM),
        "cAMPhydroPde2_kcat_12": float(CAMP_HYDROLYSIS_PDE2_KCAT),
        "cAMPhydroPde2_Km_12": float(CAMP_HYDROLYSIS_PDE2_KM),
        "cAMPhydroPde1p_kcat_13": float(CAMP_HYDROLYSIS_PDE1P_KCAT),
        "cAMPhydroPde1p_Km_13": float(CAMP_HYDROLYSIS_PDE1P_KM),
        "cAMPhydroPde1_kcat_14": float(CAMP_HYDROLYSIS_PDE1_KCAT),
        "cAMPhydroPde1_Km_14": float(CAMP_HYDROLYSIS_PDE1_KM),
        "Pde1phos_kcat_15": float(PDE1_PHOSPHORYLATION_KCAT),
        "Pde1phos_Km_15": float(PDE1_PHOSPHORYLATION_KM),
        "Pde1dephos_kcat_16": float(PDE1_DEPHOSPHORYLATION_KCAT),
        "Pde1dephos_Km_16": float(PDE1_DEPHOSPHORYLATION_KM),
        "PKAact_kF_17": float(PKA_ACTIVATION_KF),
        "PKAact_kI_17": float(PKA_ACTIVATION_KI),
        "PKAdeact_kF_18": float(PKA_DEACTIVATION_KF),
        "PKAdeact_kA_18": float(PKA_DEACTIVATION_KA),
        "GlucTrans_V_19": float(GLUCOSE_TRANSPORT_V),
        "GlucTrans_Km_19": float(GLUCOSE_TRANSPORT_KM),
        "GlucTrans_Ki_19": float(GLUCOSE_TRANSPORT_KI),
        "GlucUtil_k_20": float(GLUCOSE_UTILISATION_K),
    }


def cell_volume_litres() -> float:
    """Refuses. Williamson's compartment is size 1 with no unit; see :data:`CELL_VOLUME`."""
    return float(CELL_VOLUME)


# ----------------------------------------------------------------------------------------
# The model.
# ----------------------------------------------------------------------------------------

@dataclass(frozen=True)
class TransportReading:
    """One of the two published glucose-transport parameterisations, named by where it is printed.

    Condition (b) has no cited arbiter (:data:`TRANSPORT_CITED_SOURCE`), so both readings stay
    runnable and :func:`reproduction_report` scores both against the source's own figures.
    """

    label: str
    v: float
    km: float
    ki: float

    def __post_init__(self) -> None:
        if self.label not in ("deposit", "table_4"):
            raise ValueError(
                f"a transport reading names where it is printed, one of 'deposit' or 'table_4'; "
                f"got {self.label!r}. A third reading would be a fitted parameter")


DEPOSIT_READING = TransportReading("deposit", *DEPOSIT_TRANSPORT)
TABLE_4_READING = TransportReading("table_4", *TABLE_4_TRANSPORT)


@dataclass(frozen=True)
class WilliamsonCAMP:
    """The seventeen-state Complete cAMP Model of Additional file 1, in Python.

    Args:
        transport: Which published glucose-transport reading to run. Defaults to the deposit,
            which reproduces the source's own figures to a mean 0.30% of full scale against the
            Table 4 reading's 1.48%; both are kept runnable and neither was tuned.
        clip_negative_initials: True clips the three deposited negatives to zero, which is
            condition (c) and costs at most 7.4e-10 relative in cAMP. False runs the deposit
            byte-for-byte, and is how that cost is measured rather than asserted.
    """

    transport: TransportReading = DEPOSIT_READING
    clip_negative_initials: bool = True

    def __post_init__(self) -> None:
        audit_deposit()

    @property
    def state_names(self) -> tuple[str, ...]:
        return DYNAMIC_SPECIES

    def initial_state(self) -> np.ndarray:
        source = INITIAL_CONCENTRATIONS if self.clip_negative_initials else _INITIAL_RAW
        return np.array([source[name] for name in DYNAMIC_SPECIES], dtype=float)

    def reaction_rates(self, y) -> dict[str, float]:
        """All 20 rate laws of Table 5, in the deposit's units of mmol/L and seconds."""
        s = dict(zip(DYNAMIC_SPECIES, y))
        s.update({name: _INITIAL_RAW[name] for name in BOUNDARY_SPECIES})
        glucin = s["Glucin"]
        kd, a, b = (float(RAS2_ACTIVATION_KD), float(RAS2_ACTIVATION_A),
                    float(RAS2_ACTIVATION_B))
        ras_km = float(RAS2_ACTIVATION_KM)
        v, km, ki = self.transport.v, self.transport.km, self.transport.ki
        return {
            "Gpr1Glucass": float(GPR1_GLUCOSE_ASSOCIATION) * s["Gpr1"] * s["Glucout"],
            "Gpr1Glucdiss": float(GPR1_GLUCOSE_DISSOCIATION) * s["Gpr1Glucout"],
            "Gpa2act": float(GPA2_ACTIVATION_KA) * s["Gpr1Glucout"] * s["Gpa2i"],
            "Gpa2deact": (float(GPA2_DEACTIVATION_KA) * s["Rgs2"]
                          + float(GPA2_DEACTIVATION_K1)) * s["Gpa2a"],
            "Gpa2Krhass": float(GPA2_KRH_ASSOCIATION) * s["Gpa2a"] * s["Krh"],
            "Gpa2Krhdiss": float(GPA2_KRH_DISSOCIATION) * s["Gpa2aKrh"],
            "Ras2act": ((float(RAS2_ACTIVATION_KCAT) * s["Cdc25"] * s["Ras2i"] / ras_km)
                        * (1.0 + b * glucin / (kd * a))
                        / (1.0 + glucin / kd
                           + (s["Ras2i"] / ras_km) * (1.0 + glucin / (kd * a)))),
            "Ras2deact": (float(RAS2_DEACTIVATION_KA) * s["Ira"]
                          + float(RAS2_DEACTIVATION_K1)) * s["Ras2a"],
            "Cdc25phos": (s["C"] * s["Cdc25"] * float(CDC25_PHOSPHORYLATION_KCAT)
                          / (float(CDC25_PHOSPHORYLATION_KM) + s["Cdc25"])),
            "Cdc25dephos": (s["Mih1"] * s["Cdc25P"] * float(CDC25_DEPHOSPHORYLATION_KCAT)
                            / (float(CDC25_DEPHOSPHORYLATION_KM) + s["Cdc25P"])),
            "cAMPsynth": (s["Cyr1"] * (float(CAMP_SYNTHESIS_KCAT_GPA2) * s["Gpa2a"]
                                       + float(CAMP_SYNTHESIS_KCAT_RAS2) * s["Ras2a"])
                          / float(CAMP_SYNTHESIS_KM)),
            "cAMPhydroPde2": (s["Pde2"] * s["cAMP"] * float(CAMP_HYDROLYSIS_PDE2_KCAT)
                              / (float(CAMP_HYDROLYSIS_PDE2_KM) + s["cAMP"])),
            "cAMPhydroPde1p": (s["Pde1P"] * s["cAMP"] * float(CAMP_HYDROLYSIS_PDE1P_KCAT)
                               / (float(CAMP_HYDROLYSIS_PDE1P_KM) + s["cAMP"])),
            "cAMPhydroPde1": (s["Pde1"] * s["cAMP"] * float(CAMP_HYDROLYSIS_PDE1_KCAT)
                              / (float(CAMP_HYDROLYSIS_PDE1_KM) + s["cAMP"])),
            "Pde1phos": (s["C"] * s["Pde1"] * float(PDE1_PHOSPHORYLATION_KCAT)
                         / (float(PDE1_PHOSPHORYLATION_KM) + s["Pde1"])),
            "Pde1dephos": (s["PPA"] * s["Pde1P"] * float(PDE1_DEPHOSPHORYLATION_KCAT)
                           / (float(PDE1_DEPHOSPHORYLATION_KM) + s["Pde1P"])),
            "PKAact": (s["PKAi"] * s["cAMP"] ** 4 * float(PKA_ACTIVATION_KF)
                       / (1.0 + float(PKA_ACTIVATION_KI) * s["Krh"])),
            "PKAdeact": (float(PKA_DEACTIVATION_KF) * s["C"] ** 2
                         * (1.0 + float(PKA_DEACTIVATION_KA) * s["Krh"])),
            "GlucTrans": (v * (s["Glucout"] / km - glucin / km)
                          / (1.0 + s["Glucout"] / km + glucin / km
                             + ki * s["Glucout"] * glucin / km ** 2)),
            "GlucUtil": float(GLUCOSE_UTILISATION_K) * glucin,
        }

    def rhs(self, _t: float, y) -> np.ndarray:
        """d[species]/dt in mmol/L/s, from :data:`STOICHIOMETRY` times the rate laws."""
        rates = self.reaction_rates(y)
        index = {name: position for position, name in enumerate(DYNAMIC_SPECIES)}
        derivative = np.zeros(len(DYNAMIC_SPECIES))
        for reaction, coefficients in STOICHIOMETRY.items():
            rate = rates[reaction]
            for species, coefficient in coefficients.items():
                derivative[index[species]] += coefficient * rate
        return derivative

    def steady_state_residual(self, y=None) -> dict:
        """How far the deposited initial state is from resting under these rate laws.

        The check that needs no figure. The deposit's initial vector is a steady state its
        authors solved for, so a faithful transcription must land it at zero; this one lands it
        at 3.5e-14 /s relative, across seventeen states at once.
        """
        state = self.initial_state() if y is None else np.asarray(y, dtype=float)
        derivative = self.rhs(0.0, state)
        scale = np.maximum(np.abs(state), 1e-12)
        relative = np.abs(derivative) / scale
        worst = int(np.argmax(relative))
        return {
            "max_absolute_per_s": float(np.max(np.abs(derivative))),
            "max_relative_per_s": float(np.max(relative)),
            "worst_species": DYNAMIC_SPECIES[worst],
            "note": "does NOT test the PKA stoichiometry: at a fixed point PKAact and PKAdeact "
                    "balance under a 1:1 reading too. Only the dynamics catch that",
        }

    def conserved_moieties(self, y) -> dict[str, float]:
        s = dict(zip(DYNAMIC_SPECIES, np.asarray(y, dtype=float)))
        return {name: float(sum(weight * s[species] for species, weight in terms))
                for name, terms in CONSERVED_MOIETIES.items()}

    def simulate(self, *, duration_h: float, pulses=GLUCOSE_PULSE_PROTOCOL,
                 sample_step_s: float = 0.1) -> dict:
        """Integrate the glucose-pulse protocol. ``duration_h`` is HOURS; see condition (d).

        ``pulses`` is a sequence of ``(time_s, Glucout_mM)`` state resets, defaulting to
        Williamson's own 5 mM at 60 s and 100 mM at 240 s. The deposited SBML carries no events,
        so the pulses are applied here as resets between integration segments.
        """
        duration_s = float(duration_h) * SECONDS_PER_HOUR
        if not np.isfinite(duration_s) or duration_s <= 0.0:
            raise ValueError(f"duration_h must be finite and positive, got {duration_h!r}")
        schedule = sorted((float(t), float(level)) for t, level in pulses)
        if schedule and (schedule[0][0] < 0.0 or schedule[-1][0] >= duration_s):
            raise ValueError(
                f"pulse times must lie inside (0, {duration_s} s); got "
                f"{[t for t, _ in schedule]}")
        edges = [0.0] + [t for t, _ in schedule] + [duration_s]
        levels = [None] + [level for _, level in schedule]
        index = {name: position for position, name in enumerate(DYNAMIC_SPECIES)}

        state = self.initial_state()
        times, trajectory = [], []
        for segment, (start, stop) in enumerate(zip(edges[:-1], edges[1:])):
            if levels[segment] is not None:
                state = state.copy()
                state[index["Glucout"]] = levels[segment]
            solution = solve_ivp(self.rhs, (start, stop), state, method="LSODA",
                                 rtol=1e-10, atol=1e-16, dense_output=True, max_step=0.5)
            if not solution.success:
                raise ScientificRefusal(
                    f"integration failed on [{start}, {stop}] s: {solution.message}")
            grid = np.arange(start, stop + 0.5 * sample_step_s, sample_step_s)
            times.append(grid)
            trajectory.append(solution.sol(grid))
            state = solution.y[:, -1].copy()
        time_s = np.concatenate(times)
        values = np.concatenate(trajectory, axis=1)
        return {
            "time_s": time_s, "time_h": time_s / SECONDS_PER_HOUR,
            "species": {name: values[index[name]] for name in DYNAMIC_SPECIES},
            "transport_reading": self.transport.label,
        }

    def pka_active_fraction(self, y) -> float:
        """``percentFreeC``/100: the deposit's own assignmentRule, ``C/(C + 2*PKAi)``.

        The dimensionless quantity this module offers as the coupling to the engine's ``pka``
        signal. Dimensionless is the point -- it needs no cell volume, which Williamson does
        not supply (:data:`CELL_VOLUME`), the way ``NativeHogResponse`` couples Hog1 as a
        fraction rather than a molarity.
        """
        s = dict(zip(DYNAMIC_SPECIES, np.asarray(y, dtype=float)))
        total = s["C"] + 2.0 * s["PKAi"]
        if total <= 0.0:
            raise ScientificRefusal("no PKA in the state; the fraction is undefined")
        return float(s["C"] / total)


def glucose_pulse_timecourse(*, transport: TransportReading = DEPOSIT_READING,
                             duration_h: float = 600.0 / SECONDS_PER_HOUR,
                             clip_negative_initials: bool = True) -> dict:
    """Williamson's Figure 8 / Figure 9 protocol, 600 s by default. Convenience wrapper."""
    model = WilliamsonCAMP(transport=transport,
                           clip_negative_initials=clip_negative_initials)
    run = model.simulate(duration_h=duration_h)
    run["rates"] = _rate_traces(model, run)
    return run


def _rate_traces(model: WilliamsonCAMP, run: dict) -> dict:
    """cAMP synthesis and the two hydrolysis arms, which is what Figure 9A plots."""
    stacked = np.array([run["species"][name] for name in DYNAMIC_SPECIES])
    computed = [model.reaction_rates(stacked[:, column])
                for column in range(stacked.shape[1])]
    return {
        "cAMP_synthesis": np.array([r["cAMPsynth"] for r in computed]),
        "cAMP_hydrolysis_pde1": np.array(
            [r["cAMPhydroPde1"] + r["cAMPhydroPde1p"] for r in computed]),
        "cAMP_hydrolysis_pde2": np.array([r["cAMPhydroPde2"] for r in computed]),
    }


# ----------------------------------------------------------------------------------------
# The reproduction. Thirty-six values read off Williamson's Figures 8 and 9.
# ----------------------------------------------------------------------------------------

_FIG8 = ("Williamson 2009 Figure 8, solid blue trace (the model), read off the printed panel; "
         "y-axis 0 to 4e-3 mM")
_FIG9A = "Williamson 2009 Figure 9A, read off the printed panel; y-axis 0 to 2.5e-4 mmol/s"
_FIG9B = "Williamson 2009 Figure 9B, read off the printed panel; y-axis 0 to 1.2e-4 mM"
_FIG9C = ("Williamson 2009 Figure 9C, read off the printed panel; y-axis 0 to 1.5e-4 mM. "
          "COLOURS TAKEN AS PLOTTED, NOT AS CAPTIONED: green is C and blue is PKAi, the "
          "opposite of the caption -- see the module docstring")
_FIG9D_TOP = "Williamson 2009 Figure 9D top panel; y-axis 5e-5 to 9.5e-5 mM"
_FIG9D_BOTTOM = "Williamson 2009 Figure 9D bottom panel; y-axis 6e-6 to 12e-6 mM"


@dataclass(frozen=True)
class Landmark:
    """One value read off a printed figure, with the panel scale it was read against.

    ``full_scale`` is the panel's own y-axis span. A figure read-off cannot be more precise
    than a couple of percent of that, so it is the denominator the reproduction is scored in;
    scoring in percent of the value itself would flatter the large landmarks and punish the
    small ones for no reason connected to how well the transcription runs.
    """

    key: str
    value: float
    units: str
    full_scale: float
    where: str


PUBLISHED_LANDMARKS = (
    Landmark("cAMP_basal", 0.85e-3, "mM", 4e-3, f"{_FIG8}; the flat trace before t = 60 s"),
    Landmark("cAMP_peak1", 2.47e-3, "mM", 4e-3, f"{_FIG8}; first peak, near t = 105 s"),
    Landmark("cAMP_t240", 1.35e-3, "mM", 4e-3, f"{_FIG8}; trough at t = 240 s"),
    Landmark("cAMP_peak2", 3.18e-3, "mM", 4e-3, f"{_FIG8}; second peak, near t = 275 s"),
    Landmark("cAMP_t420", 2.18e-3, "mM", 4e-3, f"{_FIG8}; t = 420 s, the panel's last point"),
    Landmark("synth_basal", 0.74e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; blue at t = 0"),
    Landmark("pde1_basal", 0.27e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; red at t = 0"),
    Landmark("pde2_basal", 0.47e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; green at t = 0"),
    Landmark("synth_t240", 1.40e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; blue at t = 240 s"),
    Landmark("pde1_t240", 0.76e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; red at t = 240 s"),
    Landmark("pde2_t240", 0.68e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; green at t = 240 s"),
    Landmark("synth_t600", 2.03e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; blue at t = 600 s"),
    Landmark("pde1_t600", 1.22e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; red at t = 600 s"),
    Landmark("pde2_t600", 0.81e-4, "mmol/s", 2.5e-4, f"{_FIG9A}; green at t = 600 s"),
    Landmark("Gpa2i_t240", 1.015e-4, "mM", 1.2e-4, f"{_FIG9B}; blue at t = 240 s"),
    Landmark("Krh_t240", 0.905e-4, "mM", 1.2e-4, f"{_FIG9B}; red at t = 240 s"),
    Landmark("Gpa2aKrh_t240", 0.09e-4, "mM", 1.2e-4, f"{_FIG9B}; cyan at t = 240 s"),
    Landmark("Gpa2i_t600", 0.45e-4, "mM", 1.2e-4, f"{_FIG9B}; blue at t = 600 s"),
    Landmark("Krh_t600", 0.49e-4, "mM", 1.2e-4, f"{_FIG9B}; red at t = 600 s"),
    Landmark("Gpa2aKrh_t600", 0.51e-4, "mM", 1.2e-4, f"{_FIG9B}; cyan at t = 600 s"),
    Landmark("Gpa2a_t600", 0.163e-4, "mM", 1.2e-4, f"{_FIG9B}; green at t = 600 s"),
    Landmark("C_t0", 0.20e-4, "mM", 1.5e-4, f"{_FIG9C}; green at t = 0"),
    Landmark("PKAi_t0", 0.77e-4, "mM", 1.5e-4, f"{_FIG9C}; blue at t = 0"),
    Landmark("C_t240", 0.60e-4, "mM", 1.5e-4, f"{_FIG9C}; green at t = 240 s"),
    Landmark("PKAi_t240", 0.57e-4, "mM", 1.5e-4, f"{_FIG9C}; blue at t = 240 s"),
    Landmark("C_peak", 1.47e-4, "mM", 1.5e-4, f"{_FIG9C}; green maximum, near t = 295 s"),
    Landmark("PKAi_min", 0.145e-4, "mM", 1.5e-4, f"{_FIG9C}; blue minimum, near t = 295 s"),
    Landmark("C_t600", 1.04e-4, "mM", 1.5e-4, f"{_FIG9C}; green at t = 600 s"),
    Landmark("PKAi_t600", 0.35e-4, "mM", 1.5e-4, f"{_FIG9C}; blue at t = 600 s"),
    Landmark("Ras2a_basal", 4.55e-5, "mM", 9.5e-5, f"{_FIG9D_TOP}; t = 0"),
    Landmark("Ras2a_peak", 9.35e-5, "mM", 9.5e-5, f"{_FIG9D_TOP}; maximum, near t = 120 s"),
    Landmark("Ras2a_t600", 5.1e-5, "mM", 9.5e-5, f"{_FIG9D_TOP}; t = 600 s"),
    Landmark("Cdc25_t0", 12.2e-6, "mM", 12e-6, f"{_FIG9D_BOTTOM}; blue at t = 0"),
    Landmark("Cdc25P_t0", 5.8e-6, "mM", 12e-6, f"{_FIG9D_BOTTOM}; green at t = 0"),
    Landmark("Cdc25_t600", 5.6e-6, "mM", 12e-6, f"{_FIG9D_BOTTOM}; blue at t = 600 s"),
    Landmark("Cdc25P_t600", 12.4e-6, "mM", 12e-6, f"{_FIG9D_BOTTOM}; green at t = 600 s"),
)
"""Every published value this transcription is scored on. Read off the printed panels of
Figures 8 and 9 before the model was run, and none of them was adjusted afterwards."""

FIGURE_READ_TOLERANCE = 0.02
"""How close a landmark has to land, as a fraction of its panel's y-axis full scale. Two percent
of a rendered 300-pixel panel is about six pixels, which is the honest precision of reading a
printed figure -- the tolerance is set by the measurement, not by what the model achieved."""


def _computed_landmarks(transport: TransportReading) -> dict[str, float]:
    run = glucose_pulse_timecourse(transport=transport)
    t, s, r = run["time_s"], run["species"], run["rates"]

    def at(trace, when):
        return float(trace[int(np.argmin(np.abs(t - when)))])

    first = (t >= 60.0) & (t <= 240.0)
    second = t > 240.0
    camp = s["cAMP"]
    return {
        "cAMP_basal": at(camp, 59.9), "cAMP_peak1": float(camp[first].max()),
        "cAMP_t240": at(camp, 240.0), "cAMP_peak2": float(camp[second].max()),
        "cAMP_t420": at(camp, 420.0),
        "synth_basal": at(r["cAMP_synthesis"], 0.0),
        "pde1_basal": at(r["cAMP_hydrolysis_pde1"], 0.0),
        "pde2_basal": at(r["cAMP_hydrolysis_pde2"], 0.0),
        "synth_t240": at(r["cAMP_synthesis"], 240.0),
        "pde1_t240": at(r["cAMP_hydrolysis_pde1"], 240.0),
        "pde2_t240": at(r["cAMP_hydrolysis_pde2"], 240.0),
        "synth_t600": at(r["cAMP_synthesis"], 600.0),
        "pde1_t600": at(r["cAMP_hydrolysis_pde1"], 600.0),
        "pde2_t600": at(r["cAMP_hydrolysis_pde2"], 600.0),
        "Gpa2i_t240": at(s["Gpa2i"], 240.0), "Krh_t240": at(s["Krh"], 240.0),
        "Gpa2aKrh_t240": at(s["Gpa2aKrh"], 240.0),
        "Gpa2i_t600": at(s["Gpa2i"], 600.0), "Krh_t600": at(s["Krh"], 600.0),
        "Gpa2aKrh_t600": at(s["Gpa2aKrh"], 600.0), "Gpa2a_t600": at(s["Gpa2a"], 600.0),
        "C_t0": at(s["C"], 0.0), "PKAi_t0": at(s["PKAi"], 0.0),
        "C_t240": at(s["C"], 240.0), "PKAi_t240": at(s["PKAi"], 240.0),
        "C_peak": float(s["C"].max()), "PKAi_min": float(s["PKAi"].min()),
        "C_t600": at(s["C"], 600.0), "PKAi_t600": at(s["PKAi"], 600.0),
        "Ras2a_basal": at(s["Ras2a"], 0.0), "Ras2a_peak": float(s["Ras2a"].max()),
        "Ras2a_t600": at(s["Ras2a"], 600.0),
        "Cdc25_t0": at(s["Cdc25"], 0.0), "Cdc25P_t0": at(s["Cdc25P"], 0.0),
        "Cdc25_t600": at(s["Cdc25"], 600.0), "Cdc25P_t600": at(s["Cdc25P"], 600.0),
    }


def landmark_table(transport: TransportReading = DEPOSIT_READING) -> dict:
    """Published against computed for all 36 landmarks, with both error metrics."""
    computed = _computed_landmarks(transport)
    rows = {}
    for landmark in PUBLISHED_LANDMARKS:
        here = computed[landmark.key]
        rows[landmark.key] = {
            "published": landmark.value, "here": here, "units": landmark.units,
            "error_percent_of_value": 100.0 * (here - landmark.value) / landmark.value,
            "error_percent_of_full_scale":
                100.0 * abs(here - landmark.value) / landmark.full_scale,
            "where": landmark.where,
        }
    scale_errors = [row["error_percent_of_full_scale"] for row in rows.values()]
    return {
        "transport_reading": transport.label,
        "rows": rows,
        "mean_percent_of_full_scale": float(np.mean(scale_errors)),
        "max_percent_of_full_scale": float(np.max(scale_errors)),
        "worst_landmark": max(rows, key=lambda key: rows[key]["error_percent_of_full_scale"]),
        "landmarks_within_tolerance": sum(
            1 for value in scale_errors if value <= 100.0 * FIGURE_READ_TOLERANCE),
        "landmarks": len(scale_errors),
    }


def reproduction_report() -> dict:
    """Everything this transcription is scored on, computed rather than asserted.

    Four blocks: the audit of the deposit, the steady-state residual, the 36-landmark figure
    reproduction under both transport readings, and the clip cost. The test file asserts these
    numbers; none of them was tuned, and the tolerances were written after they were computed.
    """
    model = WilliamsonCAMP()
    deposit = landmark_table(DEPOSIT_READING)
    table_4 = landmark_table(TABLE_4_READING)

    clipped = glucose_pulse_timecourse(clip_negative_initials=True)["species"]["cAMP"]
    unclipped = glucose_pulse_timecourse(clip_negative_initials=False)["species"]["cAMP"]
    conserved_start = model.conserved_moieties(model.initial_state())
    run = model.simulate(duration_h=600.0 / SECONDS_PER_HOUR)
    stacked = np.array([run["species"][name] for name in DYNAMIC_SPECIES])
    conserved_drift = {
        name: float(np.max(np.abs(
            np.array([model.conserved_moieties(stacked[:, column])[name]
                      for column in range(0, stacked.shape[1], 50)]) - conserved_start[name])))
        for name in CONSERVED_MOIETIES
    }

    return {
        "availability": AVAILABILITY,
        "source": SOURCE,
        "licence": LICENCE,
        "audit": audit_deposit(),
        "steady_state_residual": model.steady_state_residual(),
        "figure_reproduction": {"deposit": deposit, "table_4": table_4},
        "transport_arbitration": {
            "question": "Table 4 says KM 0.08 / V 1.7; the deposit says Km 1.7 / V 0.062",
            "cited_source_can_arbitrate": False,
            "why": TRANSPORT_CITED_SOURCE.reason,
            "arbitrated_against": "the source's own printed Figures 8 and 9, all 36 landmarks",
            "deposit_mean_percent_of_full_scale": deposit["mean_percent_of_full_scale"],
            "deposit_max_percent_of_full_scale": deposit["max_percent_of_full_scale"],
            "table_4_mean_percent_of_full_scale": table_4["mean_percent_of_full_scale"],
            "table_4_max_percent_of_full_scale": table_4["max_percent_of_full_scale"],
            "chosen": DEPOSIT_READING.label,
            "not_a_transposition": "the file's Km IS the paper's V (1.7), but the file's V "
                                   "(0.062) is not the paper's Km (0.08)",
        },
        "negative_initial_clip": {
            "clipped": dict(NEGATIVE_INITIAL_CONCENTRATIONS),
            "left_as_deposited": dict(NUMERICAL_ZERO_INITIAL_CONCENTRATIONS),
            "max_absolute_camp_difference_mM": float(np.max(np.abs(clipped - unclipped))),
            "max_relative_camp_difference": float(
                np.max(np.abs((clipped - unclipped) / clipped))),
        },
        "conservation": {"initial": conserved_start, "max_drift_over_600_s": conserved_drift},
        "scope": {
            "closes": "the cAMP/PKA input side: how pka is COMPUTED",
            "does_not_close": "Snf1 -> Adr1/Cat8 -> CSRE, the panel's declared carbon module",
            "refusal": PROMOTER_COUPLING_NOT_BUILT.reason,
        },
        "gate": gate().summary(),
        "biological_validation": False,
        "source_parameters_refitted": False,
    }


def gate():
    """Criterion (e) for this piece. It REFUSES, and the docstring says why it should."""
    return WILLIAMSON_PARAMS.gate()


def provenance() -> str:
    """The graded parameter tables, in the shape every other mech module reports them."""
    return provenance_table(WILLIAMSON_PARAMS, NOT_IN_SOURCE)


WILLIAMSON_PARAMS.add_target(Target(
    name="rolland2000_camp_timecourse",
    observable="cytosolic cAMP after 5 mM glucose at 60 s and 100 mM at 240 s",
    assay="acetylated cAMP radioimmunoassay on glucose-starved cell suspension, converted from "
          "nmol/g wet weight to nM by Williamson's Eqn 5 (Cw 0.15, Vc 2.68e-6 L per 1e7 cells)",
    noise_floor=1e-4,
    units="mM",
    source="Rolland F et al. 2000, Mol Microbiol 38:348, PMID 11069660, cited as [38] by "
           "Williamson's Methods and misprinted as [35] in the Figure 8 caption",
    fitted=True))
"""The only real measurement anywhere near this model -- and the one its 35 fitted parameters were
estimated against, so it is registered ``fitted=True`` and does not count in the gate."""


if __name__ == "__main__":  # pragma: no cover
    import json

    print(json.dumps(reproduction_report(), indent=2, default=float))
