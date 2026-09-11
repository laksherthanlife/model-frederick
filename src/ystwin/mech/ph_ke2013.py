"""Ke 2013's ion-regulation ODE system, transcribed BY HAND IN THIS REPOSITORY.

WHAT THIS FILE IS, AND WHAT IT IS NOT. Ke, Ingram & Haynes 2013, "An integrative model of
ion regulation in yeast" (PLoS Comput Biol 9: e1002879, PMC3547829, CC-BY) is the only
complete published *S. cerevisiae* ODE system in which intracellular pH is an OUTPUT, and
it ships **no executable artifact of any kind** -- no SBML, no MATLAB, nothing in BioModels
or JWS Online (`data/parameter_evidence.json` records that search). What follows is a
hand transcription of its printed equations and tables into runnable Python. The equations
and every number are Ke's, read off the seven supplementary PDFs this repository has pinned
and checksummed. **The executable file is ours.** That is why it carries

    AVAILABILITY = "transcribed_here"

and not ``local_verified``, which in this repository means a checksummed file that the
upstream authors published. A transcription is not upstream provenance. It earns its place
by reproducing the source's own published numbers -- see the next section -- not by the fact
that it runs.

WHAT REPRODUCES AND WHAT DOES NOT, stated before the code because it is the result.
Ke's unstressed condition is external Na+ 5 mM, K+ 1 mM, pH 6.5, and the article reports
four numbers for it. This transcription, run to steady state with the reductions declared
below, gives:

    quantity                     Ke 2013            here          difference
    intracellular pH             7.14 (Fig. 7E)     7.267         +0.127 pH units
    membrane potential           ~ -94 mV (KCl)     -110.9 mV     -16.9 mV
    [Na+]cyt                     100 mM  (Table S5) 90.3 mM       -9.7%
    [K+]cyt                      200 mM  (Table S5) 209.7 mM      +4.8%
    pH_i after a step to pH 8.0  7.25 (Fig. 7E)     7.342         +0.092
    the SHIFT that step causes   +0.11              +0.075        68% of it

The direction, the mechanism and the order of magnitude carry; the levels do not match to
the printed precision. That is a RESULT about the transcription, not a defect to be tuned
away, and nothing here was fitted to close it. The Table S5 rows are a genuinely
independent check: the article states those initial concentrations ARE the unstressed
steady state, so a transcription landing within 5-10% of them without being told to is
evidence the flux laws were read correctly. The residual disagreement is bounded below by
the three ambiguities in the next section, only one of which we could resolve.

THREE PRINTED FORMS THAT DO NOT SURVIVE DIMENSIONAL ANALYSIS, and what was done with each.
Each is recorded verbatim in :data:`READING_CHOICES` alongside the reading actually used,
and each reading is switchable, so no correction is hidden inside a function body.

1. Eqn 2.38, H+ uptake: printed ``k * (R T ln([H+]ext/[H+]int) - Em F)``. The bracket is
   then J/mol, but Table S2 gives k_H_uptake in ``1e-18 mol/(s*V)``, so the printed product
   is 96485x too large. Taken as ``k * (R T / F ln(...) - Em)``, which is the only reading
   consistent with the paper's own unit column. The printed-law equilibrium has
   J_H_uptake = 6.07e3 amol/s and a matching J_Pma1 of 6.07e3 amol/s, at pH_i = 4.31 --
   2.83 pH units below Ke's printed 7.14. This is a mathematical equilibrium, not a
   physiological prediction; its existence does not repair the dimensional discrepancy.
2. The Tok1p gating rates: printed ``exp(+/- l * Em / (R T))``. Em is in volts and R T is in
   J/mol, so the exponent is ~1e-5 and the gating is voltage-INDEPENDENT, contradicting the
   same page's statement that "the rates of these transitions are governed by membrane
   potential". Taken as ``exp(+/- l F Em / (R T))``, matching Eqn 2.19's own GHK exponent
   two paragraphs later. This one costs 27 mV: the printed reading gives Em = -138.5 mV and
   pH_i 7.252, the corrected one -110.9 mV and pH_i 7.267. BOTH are reported.
3. Eqn 2.8: its second line expands J_Ena1,K with ``[Na+]``, while J0_Ena1,K one paragraph
   earlier is defined with ``[K+]``. Read as ``[K+]``; a K+ flux driven by the Na+ gradient
   is a typesetting slip, not a model.

THE TWO SCALARS KE 2013 STRUCTURALLY DOES NOT CONTAIN. This model supplies the missing
STRUCTURE for the 'ph' programme and not the missing numbers, and the shortfall is not an
oversight this file can repair:

* **No cytosolic buffering capacity, at all.** The string "buffer" appears ZERO times in the
  article XML and in all seven supplements. Eqn 2.40 writes ``dH+/dt`` on FREE protons with
  no buffering term, so every proton entering the cytosol moves the pH. A model that omits a
  quantity is not evidence the quantity is small. Its omission limits physiological
  interpretation, but cannot by itself explain a disagreement between two implementations
  that both omit it. :data:`NOT_IN_SOURCE` refuses it by name.
* **A fixed 300 mM anion, with no reference and no measurement.** Table S4's ``Anion`` row
  ("Intracellular negative charges that balances H+, K+ and Na+") is the ONLY row in that
  table with an empty reference column, and Eqn 2.1 makes the membrane potential the
  difference between the cation total and that number, divided by a capacitance -- so Em is
  0.192 V per attomole of net charge and 300 mM is doing all the work. It is ASSERTED here
  with that said out loud, never MEASURED.

Both stay bracketed. Neither is fitted here and neither should be fitted downstream.

WHAT COULD NOT BE TRANSCRIBED, and why -- recorded in :data:`UNTRANSCRIBED`, never guessed:

* the whole Hog1p/glycerol module (Text S1 3.1). Its rate constants are in Table S3, but
  ``Vcyt``, ``Vnuc`` and the delay ``tau`` appear in the equations and in NO table; they
  live in the upstream Zi 2010 model, which this repository has not pinned;
* the Nrg1p ODE (Eqn 3.5), which needs ``KmNrg1,pH``. Table S3 lists ``kNrg1,pH`` (5.25e-8
  mM/s) and no Km at all, and a rate constant is not a Michaelis constant;
* consequently ENA1 mRNA and Ena1p (Eqn 3.6), which are downstream of Nrg1p. Table S3 also
  prints ``KmENA1,Nrg1`` as "1.143.0*1e-4", which is not a number;
* the split of ``kNha1`` into the ``kNha1_low`` and ``kNha1_high`` that Eqns 2.16-2.17
  require. Table S2 gives one ``kNha1``. Both are set to it here and the cost is MEASURED
  rather than assumed: driving ``kNha1_low`` from 0 to 10x moves pH_i by 0.004 units.

TWO FURTHER INTERNAL INCONSISTENCIES THIS TRANSCRIPTION FOUND, both computed in
``tests/test_ke2013_transcription.py`` rather than asserted here:

* the volume module is not at rest at Table S5's initial volume. At Volume_cyt = 20.1062
  um^3 the printed Eqns 4.2-4.4 give D_Pressure = -1.51e6 J/m^3, not 0, because
  [Osmo_Cyt] = 1176 mM against [Osmo_Ext] = 250 mM and a turgor of only 0.875e6 J/m^3. The
  printed equations settle at 1.32x the published initial volume;
* the earlier claim of an acidic nonexistence cutoff is withdrawn. A mixed-unit residual
  test and ill-conditioned ion coordinates mistook solver failures for missing roots.
  Positive, charge-scaled coordinates find an equilibrium at external pH 5.5, with
  intracellular pH about 4.68, without refitting a parameter. This is not physiological
  validation: the fixed-volume, fixed-regulator reduction supplies no validated acidic
  viability threshold, and ``acidic_limit`` now explicitly refuses that quantity;
* Table S5's ``Init_CNon = 0`` is not the steady state of Eqn 3.4 at Table S5's own Ca2+ and
  the Ppz level Eqn 2.21 gives. The block relaxes to 1.20e-5 mM activated calcineurin, and
  Crz1p with it to 2.54x the Table S5 row that the same table's footnote calls the
  unstressed steady state.

THE REDUCTIONS THE REPRODUCTION RUN MAKES, all three declared and all three forced by the
list above rather than chosen for convenience:

1. the Hog1p/glycerol states are held at Table S5, which the article itself states are the
   unstressed steady state. Untranscribable otherwise (Vcyt, Vnuc, tau);
2. the cell volume is held at Table S5's 20.1062 um^3, because the printed volume module
   does not rest there. Legitimate ONLY for the iso-osmotic comparison used here -- a step
   in external pH at constant osmolarity. It is NOT legitimate for the NaCl, sorbitol or
   KCl figures, and this module therefore does not claim them;
3. ``Ena1(t)/Ena1_0`` is held at 1, since Eqn 3.6 is untranscribable. The article attributes
   the alkaline pH_i rise to "reductions in the activities of those secondary transporters
   that utilize the H+ gradient", i.e. to the block that IS transcribed, so the comparison
   tests the paper's own stated mechanism.

UNITS. Amounts are attomoles (1 amol = 1e-18 mol), which is the unit Ke's Table S2 writes as
"1e-18 mol"; volumes are um^3, and 1 amol/um^3 = 1 mM exactly, which is why the concentration
conversion below is a plain division. Time is seconds throughout, as Table S2 and S3 print.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares
from scipy.special import expit

from .contracts import Compartment, ScientificRefusal, finite
from .params import Param, ParamRegistry, provenance_table

__all__ = [
    "AVAILABILITY",
    "AS_PRINTED",
    "Environment",
    "EquilibriumNotConverged",
    "IonState",
    "KE2013_PARAMS",
    "NOT_IN_SOURCE",
    "ODE_SYSTEM",
    "ODE_TOTALS",
    "PUBLISHED",
    "READING_CHOICES",
    "Reading",
    "STATE_VARIABLES",
    "STEADY_STATE_TOLERANCE",
    "STEADY_STATE_RELATIVE_TOLERANCE",
    "UNIT_CONSISTENT",
    "UNTRANSCRIBED",
    "acidic_limit",
    "alkaline_step",
    "calcineurin_rhs",
    "fluxes",
    "ion_rhs",
    "membrane_potential_V",
    "pressure_difference",
    "reproduction_report",
    "steady_state",
    "transcription_record",
    "turgor_pressure",
    "write_transcription_record",
]

AVAILABILITY = "transcribed_here"
"""The level this run introduces, and the reason it is not ``local_verified``.

``local_verified`` means a byte-checksummed file the upstream authors published.
``equations_only`` -- what `data/parameter_evidence.json` records for Ke 2013 -- means the
equations are pinned but nothing runs. This level sits between them: the EQUATIONS are
upstream and checksummed, the EXECUTABLE is this repository's, and a reader must be able to
tell those apart at a glance. Nothing transcribed here may ever be promoted past it.
"""

SOURCE = ("Ke R, Ingram PJ, Haynes K 2013, 'An integrative model of ion regulation in "
          "yeast', PLoS Comput Biol 9(1): e1002879, doi:10.1371/journal.pcbi.1002879, "
          "PMC3547829, CC-BY-4.0")

ARTIFACTS = {
    "article_xml": ("data/native_reference_models/ke2013/PMC3547829.xml",
                    "9d9004232e30e6f12e8db8b8140aa9ac55c59069ded74f1db7a9ebf213967bb2"),
    "table_s1_odes": ("data/native_reference_models/ke2013/pcbi.1002879.s009.pdf",
                      "d3d77b430d6985816fb7fd2c779b8dbb9ca1c987f541eef0716b7309bbefc9b4"),
    "table_s2_transporters": ("data/native_reference_models/ke2013/pcbi.1002879.s010.pdf",
                              "53c0423bc25fe35b703529b4da8d446864e9d1f4f8e98a6d1a724729a47264ad"),
    "table_s3_signalling": ("data/native_reference_models/ke2013/pcbi.1002879.s011.pdf",
                            "f6de58f7b3f1f1fb700ca19027adfc94e12433d9ce312426e41f6411da235ad4"),
    "table_s4_physiology": ("data/native_reference_models/ke2013/pcbi.1002879.s012.pdf",
                            "cb01b6953f12f1e74ab7e140aef0dc69a20aebf80adf0ac3afeff9b5ed778af6"),
    "table_s5_initial": ("data/native_reference_models/ke2013/pcbi.1002879.s013.pdf",
                         "edc0da0939d90e90f622d852cd4053f83a2a15b56f73f6173805ee6a981f464c"),
    "table_s6_sensitivity": ("data/native_reference_models/ke2013/pcbi.1002879.s014.pdf",
                             "d638c12eb3ad983988ba2690fa0c9bbe90a58534d8a79cc4c88cbd2d2250e8ce"),
    "text_s1_derivations": ("data/native_reference_models/ke2013/pcbi.1002879.s015.pdf",
                            "b0b69fd25483b508e66ffcfdcb04decf8cddeacdeee651c15b5a19b5eab661aa"),
}
"""The eight pinned upstream files, by role, with the sha256 already in
`data/parameter_evidence.json`. The test re-hashes them: a transcription whose source has
moved underneath it is not a transcription of anything."""

RECORD_PATH = Path("data/transcriptions/ke2013.json")

PUBLISHED = {
    "unstressed_ph_i": 7.14,
    "alkaline_ph_i": 7.25,
    "unstressed_em_mV": -94.0,
    "kcl_em_mV": 20.0,
    "unstressed_na_mM": 100.0,
    "unstressed_k_mM": 200.0,
}
"""Every number this transcription is scored against, with where it is printed.

``unstressed_ph_i`` / ``alkaline_ph_i``: Results, alkaline section -- "intracellular pH
increased from 7.14 to 7.25 upon pH 8.0 stress (Fig. 7E)". ``unstressed_em_mV``: Results,
KCl section -- "restored the membrane potential to the unstressed level (around -94 mV)".
``kcl_em_mV``: same paragraph, "the membrane potential was up to +20 mV" under 0.8 M KCl;
recorded but NOT claimed here, because 0.8 M KCl is an osmotic stress and reduction (2)
holds the volume fixed. ``unstressed_na_mM`` / ``unstressed_k_mM``: Table S5, whose footnote
states these are the unstressed steady state.
"""


KE2013_PARAMS = ParamRegistry("mech/ph_ke2013.py::ke2013-transcription")
"""Every number Ke 2013 prints that this transcription uses.

EVERY ROW IS ASSERTED, AND THAT IS DELIBERATE. Ke's Tables S2 and S3 carry a provenance
column of their own -- "Taken from [n]" against "Estimated from [n]", with an asterisk on
the rows "adjusted during the model integration step". That column is carried through
VERBATIM in each ``source`` string, because it is the only provenance these numbers have and
losing it would be the whole failure this repository exists to avoid. It is not enough to
promote a row: the upstream references behind "Taken from" were not traced to their own
primary measurements here, most are not *S. cerevisiae*, and an asterisked row is a fit to
Ke's own integration step. A defined SI constant is ASSERTED here too, for the reason
`params.py` gives: grading a constant MEASURED was one of the flattering grades that failed.
"""

NOT_IN_SOURCE = ParamRegistry("mech/ph_ke2013.py::absent-from-ke2013")
"""The two quantities Ke 2013 structurally does not contain. See the module docstring."""

_S2 = "Ke 2013 Table S2 (transporter module, pinned s010.pdf), row"
_S3 = "Ke 2013 Table S3 (transcriptional module, pinned s011.pdf), row"
_S4 = "Ke 2013 Table S4 (basic physiological parameters, pinned s012.pdf), row"
_S5 = "Ke 2013 Table S5 (initial concentrations, pinned s013.pdf), row"
_STAR = (" ASTERISKED in the source: 'the values of these parameters were adjusted during "
         "the model integration step such that simulation results are consistent with "
         "experimental data', i.e. fitted by Ke, not measured")

# Table S4. Physical constants and cell geometry.
R_GAS = KE2013_PARAMS.add(Param.asserted(
    "R", 8.31447, "J/(mol*K)", f"{_S4} 'R, Gas constant, 8.31447' (no reference column)"))
TEMPERATURE_K = KE2013_PARAMS.add(Param.asserted(
    "T", 310.0, "K", f"{_S4} 'T, Absolute temperature, 310 K' (no reference column)"))
FARADAY = KE2013_PARAMS.add(Param.asserted(
    "F", 96485.0, "C/mol", f"{_S4} 'F, Faraday constant, 96485 C/mol' (no reference column)"))
G0_ATP = KE2013_PARAMS.add(Param.asserted(
    "G0_ATP", -36.03e3, "J/mol",
    f"{_S4} 'G0_ATP, Gibbs energy for ATP hydrolysis, -36.03*1e+3 J' -- the printed unit is "
    "J, not J/mol; it is used against R*T in every flux law, so J/mol is the only "
    "dimensionally possible reading"))
MEMBRANE_CAPACITANCE = KE2013_PARAMS.add(Param.asserted(
    "Cm", 1e-2, "F/m^2", f"{_S4} 'Cm, Membrane capacitance, 1e-2 farad/m2' [2]"))
CELL_SURFACE_UM2 = KE2013_PARAMS.add(Param.asserted(
    "Smem", 50.2655, "um^2", f"{_S4} 'Smem, Cell surface area, 50.2655 um2'"))
CELL_VOLUME_UM3 = KE2013_PARAMS.add(Param.asserted(
    "VolumeCell", 33.510, "um^3", f"{_S4} 'VolumeCell, Cell volume, 33.510 um3'"))
PERMEABILITY_K = KE2013_PARAMS.add(Param.asserted(
    "PK", 1e-6, "um/s", f"{_S4} 'PK, Membrane permeability to potassium ion, 1e-6 um/s' [67,68]"))
PERMEABILITY_NA = KE2013_PARAMS.add(Param.asserted(
    "PNa", 3.8e-8, "um/s", f"{_S4} 'PNa, Membrane permeability to sodium ion, 3.8*1e-8 um/s' [67,68]"))
ATP_MM = KE2013_PARAMS.add(Param.asserted(
    "ATP", 2.6, "mM",
    f"{_S4} '[ATP], ATP concentration, 2.6' [70,71]. The table prints the unit as 'mM-1' for "
    "this row and for ADP and Pi; mM is the only reading that makes the flux laws close"))
ADP_MM = KE2013_PARAMS.add(Param.asserted(
    "ADP", 1.0, "mM", f"{_S4} '[ADP], ADP concentration, 1.0' [70,71]; printed unit 'mM-1'"))
PI_MM = KE2013_PARAMS.add(Param.asserted(
    "Pi", 3.0, "mM", f"{_S4} '[Pi], Phosphate concentration, 3.0' [70,71]; printed unit 'mM-1'"))
HYDRAULIC_PERMEABILITY = KE2013_PARAMS.add(Param.asserted(
    "Lp", 1.19e6, "m^4/(J*s)", f"{_S4} 'Lp, Hydraulic membrane permeability, 1.19*1e+6 m4/(J*s)' [51]"))
GEOMETRIC_FACTOR = KE2013_PARAMS.add(Param.asserted(
    "GEK", 7.85e-11, "m^2", f"{_S4} 'GEK, Geometrical factor, 7.85*1e-11 m2' [51]"))
OSMO_CYT_0 = KE2013_PARAMS.add(Param.asserted(
    "Osmocyt0", 600.0, "mM", f"{_S4} 'Osmocyt0, Total concentration of initial intracellular osmolytes, 600 mM' [51]"))
OSMO_EXT_0 = KE2013_PARAMS.add(Param.asserted(
    "Osmoext0", 250.0, "mM", f"{_S4} 'Osmoext0, Total concentration of initial extracellular osmolytes, 250 mM' [51]"))
TURGOR_0 = KE2013_PARAMS.add(Param.asserted(
    "Turgor0", 0.875e6, "J/m^3", f"{_S4} 'Turgor0, Initial turgor pressure, 0.875*1e+6 J/m3' [51]"))
R_VOL = KE2013_PARAMS.add(Param.asserted(
    "rvol", 0.63, "dimensionless", f"{_S4} 'rvol, Volumetric elastic modulus, 0.63' [51]"))
ANION_MM = KE2013_PARAMS.add(Param.asserted(
    "Anion", 300.0, "mM",
    f"{_S4} 'Anion, Intracellular negative charges that balances H+, K+ and Na+, 300 mM'. "
    "THE ONLY ROW IN TABLE S4 WITH AN EMPTY REFERENCE COLUMN. No measurement of it appears "
    "anywhere in the article or in the seven supplements, and Eqn 2.1 makes the membrane "
    "potential the difference between the cation total and this number",
    missing="a determination of the non-permeant intracellular anion charge for the working "
            "strain and medium; Em here moves 0.192 V per attomole of net charge, so this "
            "single unmeasured scalar sets the whole membrane potential"))

# Table S2. Transporter module. The provenance column travels in the source string.
K_ENA1 = KE2013_PARAMS.add(Param.asserted(
    "kEna1", 1.98e-4, "amol/(s*mM^2)", f"{_S2} 'kEna1*, 1.98*1e-4, 1e-18 mol/(s*mM2), Estimated from [27]'.{_STAR}"))
KM_ENA1_NA = KE2013_PARAMS.add(Param.asserted(
    "KmEna1,Na", 200.0, "mM", f"{_S2} 'KmEna1,Na, 200, mM, Taken from [58]'"))
KM_ENA1_K = KE2013_PARAMS.add(Param.asserted(
    "KmEna1,K", 300.0, "mM", f"{_S2} 'KmEna1,K, 300, mM, Estimated from [58]'"))
K_NHA1 = KE2013_PARAMS.add(Param.asserted(
    "kNha1", 1.092e5, "amol/(s*mM^5)", f"{_S2} 'kNha1*, 1.092*1e+5, 1e-18 mol/(s*mM5), Estimated from [12]'.{_STAR}"))
KM_NHA1_NA = KE2013_PARAMS.add(Param.asserted(
    "KmNha1,Na", 12.7, "mM", f"{_S2} 'KmNha1,Na, 12.7, mM, Taken from [36]'"))
KM_NHA1_HIGH_K = KE2013_PARAMS.add(Param.asserted(
    "KmNha1_high,K", 12.4, "mM", f"{_S2} 'KmNha1_high,K, 12.4, mM, Taken from [36]'"))
KM_NHA1_LOW_K = KE2013_PARAMS.add(Param.asserted(
    "KmNha1_low,K", 1240.0, "mM", f"{_S2} 'KmNha1_low,K, 1240' -- units and reference columns both EMPTY"))
KM_NHA1_HOG1 = KE2013_PARAMS.add(Param.asserted(
    "KmNha1,Hog1", 1e-4, "mM", f"{_S2} 'KmNha1,Hog1, 1e-4, mM' -- reference column EMPTY"))
PS_TOK1 = KE2013_PARAMS.add(Param.asserted(
    "Ps,Tok1", 4.6642, "amol*m^3/(s*um^3)",
    f"{_S2} 'Ps,Tok1*, 4.6642, 1e-18 m/s, Estimated from [59]'.{_STAR} The printed unit m/s "
    "is a permeability per area, but Eqn 2.19 divided by F must give amol/s and the paper "
    "calls it 'the TOTAL permeability of Tok1p at the plasma membrane', i.e. area already "
    "included; read as 1e-18 m^3/s"))
ALPHA_TOK1 = KE2013_PARAMS.add(Param.asserted(
    "alphaTok1", 231.6, "dimensionless", f"{_S2} 'aTok1, 231.6, Taken from [38]'"))
BETA_TOK1 = KE2013_PARAMS.add(Param.asserted(
    "betaTok1", 8122.7, "dimensionless", f"{_S2} 'bTok1, 8122.7, Taken from [38]'"))
L_TOK1_EXT = KE2013_PARAMS.add(Param.asserted(
    "lTok1,ext", 0.0, "dimensionless", f"{_S2} 'lTok1,ext, 0.0, Taken from [38]'"))
L_TOK1_INT = KE2013_PARAMS.add(Param.asserted(
    "lTok1,int", 0.76875, "dimensionless", f"{_S2} 'lTok1,int, 0.76875, Taken from [38]'"))
K_TOK1_1 = KE2013_PARAMS.add(Param.asserted(
    "kTok1,1", 3.4e7, "1/s", f"{_S2} 'kTok1,1, 3.4*1e+7, s-1, Taken from [38]'"))
K_TOK1_2 = KE2013_PARAMS.add(Param.asserted(
    "kTok1,2", 3.4e7, "1/s", f"{_S2} 'kTok1,2, 3.4*1e+7, s-1, Taken from [38]'"))
K_TOK1_OR = KE2013_PARAMS.add(Param.asserted(
    "kTok1,OR", 1e4, "1/s", f"{_S2} 'kTok1,OR, 1e+4, s-1, Taken from [38]'"))
K_TOK1_RO = KE2013_PARAMS.add(Param.asserted(
    "kTok1,RO", 1e4, "1/s", f"{_S2} 'kTok1,RO, 1e+4, s-1, Taken from [38]'"))
KM_TOK1_HOG1 = KE2013_PARAMS.add(Param.asserted(
    "KmTok1,Hog1", 2e-5, "mM", f"{_S2} 'KmTok1,Hog1, 2*1e-5' -- units and reference columns both EMPTY"))
K_TRK = KE2013_PARAMS.add(Param.asserted(
    "kTrk", 54.0, "amol/(s*V^2)", f"{_S2} 'kTrk*, 54.0, 1e-18 mol/(s*V2), Estimated from [6,55]'.{_STAR}"))
KM_TRK_HIGH_K = KE2013_PARAMS.add(Param.asserted(
    "KmTrk_high,K", 0.01, "mM", f"{_S2} 'KmTrk_high,K, 0.01, mM, Taken from [39]'"))
KM_TRK_HIGH_NA = KE2013_PARAMS.add(Param.asserted(
    "KmTrk_high,Na", 100.0, "mM", f"{_S2} 'KmTrk_high,Na, 100.0, mM' -- reference column EMPTY"))
KM_TRK_MEDIUM_K = KE2013_PARAMS.add(Param.asserted(
    "KmTrk_medium,K", 0.4, "mM", f"{_S2} 'KmTrk_medium,K, 0.4, mM, Taken from [39]'"))
KM_TRK_MEDIUM_NA = KE2013_PARAMS.add(Param.asserted(
    "KmTrk_medium,Na", 40.0, "mM", f"{_S2} 'KmTrk_medium,Na, 40.0, mM' -- reference column EMPTY"))
KM_TRK_PPZ = KE2013_PARAMS.add(Param.asserted(
    "KmTrk,Ppz", 1e-5, "mM", f"{_S2} 'KmTrk,Ppz, 1e-5, mM, Estimated from [9,55]'"))
KM_TRK_CN = KE2013_PARAMS.add(Param.asserted(
    "KmTrk,Cn", 5e-4, "mM", f"{_S2} 'KmTrk,Cn, 5*1e-4, mM, Estimated from [9,55]'"))
PPZ_TOTAL = KE2013_PARAMS.add(Param.asserted(
    "Ppz0", 6.607e-5, "mM", f"{_S2} 'Ppz0, 6.607*1e-5, mM, [60]' -- the reference column carries the citation with no Taken/Estimated verb"))
KM_PPZ = KE2013_PARAMS.add(Param.asserted(
    "KmPpz", 7.94e-5, "mM", f"{_S2} 'KmPpz, 7.94*1e-5, mM' -- reference column EMPTY"))
K_NSC1 = KE2013_PARAMS.add(Param.asserted(
    "kNSC1", 6e-4, "amol/(s*mM)", f"{_S2} 'kNHS1*, 6*1e-4, 1e-18 mol/(s*mM)' -- printed 'kNHS1', used as kNSC1 in Eqns 2.31-2.32; reference column EMPTY.{_STAR}"))
KM_NSC1_K = KE2013_PARAMS.add(Param.asserted(
    "KmNSC1,K", 60.0, "mM", f"{_S2} 'KmNSC1,K, 60.0, mM, Taken from [61]'"))
KM_NSC1_NA = KE2013_PARAMS.add(Param.asserted(
    "KmNSC1,Na", 100.0, "mM", f"{_S2} 'KmNSC1,Na, 100.0, mM' -- reference column EMPTY"))
K_H_PROD = KE2013_PARAMS.add(Param.asserted(
    "KH_prod", 5.0, "amol/s", f"{_S2} 'KH_prod, 5, 1e-18 mol/s' -- reference column EMPTY"))
K_H_UPTAKE = KE2013_PARAMS.add(Param.asserted(
    "kH_uptake", 10.0, "amol/(s*V)", f"{_S2} 'kH_uptake, 10.0, 1e-18 mol/(s*V)' -- reference column EMPTY. The printed unit is what forces reading choice 1"))
K_PMA1 = KE2013_PARAMS.add(Param.asserted(
    "kPma1", 4.8e4, "amol/(s*mM^2)", f"{_S2} 'kPma1*, 4.8*1e+4, 1e-18 mol/(s*mM2), Estimated from [45]'.{_STAR}"))

# Table S3. Only the calcineurin rows are used; the Hog1p and Nrg1p rows are untranscribable.
C_CA = KE2013_PARAMS.add(Param.asserted(
    "CCa", 2.5e-6, "mM/s", f"{_S3} 'CCa, 2.5*1e-6, mM s-1' -- reference column EMPTY"))
D_CA = KE2013_PARAMS.add(Param.asserted(
    "dCa", 0.05, "1/s", f"{_S3} 'dCa, 0.05, s-1, Estimated from [25]'"))
H_NA_CYT = KE2013_PARAMS.add(Param.asserted(
    "h_Na_cyt", 12.0, "dimensionless", f"{_S3} 'h_Na_cyt*, 12, Estimated from [25,27]'.{_STAR}"))
H_NA_EXT = KE2013_PARAMS.add(Param.asserted(
    "h_Na_ext", 2.0, "dimensionless", f"{_S3} 'h_Na_ext*, 2, Estimated from [25,27]'.{_STAR}"))
K_CA_CYT = KE2013_PARAMS.add(Param.asserted(
    "kCa,cyt", 8.0e-6, "mM/s", f"{_S3} 'kCa,cyt*, 8.0*1e-6, mM s-1, Estimated from [25,27]'.{_STAR}"))
K_CA_EXT = KE2013_PARAMS.add(Param.asserted(
    "kCa,ext", 2.4e-5, "mM/s", f"{_S3} 'kCa,ext*, 2.4*1e-5, mM s-1, Estimated from [25,27]'.{_STAR}"))
KM_CA_CYT = KE2013_PARAMS.add(Param.asserted(
    "KmCa_cyt", 180.0, "mM", f"{_S3} 'KmCa_cyt, 180.0, mM, Estimated from [25,27]'"))
KM_CA_EXT = KE2013_PARAMS.add(Param.asserted(
    "KmCa_ext", 500.0, "mM", f"{_S3} 'KmCa_ext, 500.0, mM, Estimated from [25,27]'"))
K_CA_PH = KE2013_PARAMS.add(Param.asserted(
    "kCa,pH", 1.5e-5, "mM/s", f"{_S3} 'kCa,pH*, 1.5*1e-5, mM s-1' -- reference column EMPTY.{_STAR}"))
K_CN_A = KE2013_PARAMS.add(Param.asserted(
    "kCN,a", 1.0e9, "1/(mM^3*s)", f"{_S3} 'kCN,a, 1.0e+9, mM-3 s-1, Estimated from [62]'"))
K_CN_DA = KE2013_PARAMS.add(Param.asserted(
    "kCN,da", 5.0e-4, "1/s", f"{_S3} 'kCN,da, 5.0*1e-4, s-1, Estimated from [62]'"))
K_CN_PPZ_DA = KE2013_PARAMS.add(Param.asserted(
    "kCN_Ppz,da", 200.0, "1/(mM*s)", f"{_S3} 'kCN_Ppz,da*, 200.0, mM-1 s-1, Estimated from [63]'.{_STAR}"))
C_CRZ1 = KE2013_PARAMS.add(Param.asserted(
    "CCrz1", 7.866e-7, "mM/s", f"{_S3} 'CCrz1, 7.866*1e-7, mM s-1' -- reference column EMPTY"))
D_CRZ1 = KE2013_PARAMS.add(Param.asserted(
    "dCrz1", 4.1e-3, "1/s", f"{_S3} 'dCrz1, 4.1*1e-3, s-1' -- reference column EMPTY"))
K_CRZ1 = KE2013_PARAMS.add(Param.asserted(
    "kCrz1", 0.1, "1/s", f"{_S3} 'kCrz1*, 0.1, s-1' -- reference column EMPTY.{_STAR}"))

# Table S5. Initial concentrations, which the source states ARE the unstressed steady state.
INIT_H_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_H", 1e-4, "mM", f"{_S5} 'Init_H, Initial intracellular H+ concentration, 1e-4 (pH7), mM'"))
INIT_K_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_K", 200.0, "mM", f"{_S5} 'Init_K, Initial intracellular K+ concentration, 200, mM'"))
INIT_NA_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_Na", 100.0, "mM", f"{_S5} 'Init_Na, Initial intracellular Na+ concentration, 100, mM'"))
INIT_VOLUME_UM3 = KE2013_PARAMS.add(Param.asserted(
    "Init_CytVolume", 20.1062, "um^3", f"{_S5} 'Init_CytVolume, Initial volume for the osmotically changeable compartment, 20.1062, um3'"))
INIT_HOG1PPC_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_Hog1PPc", 4.443e-6, "mM", f"{_S5} 'Init_Hog1PPc, Initial phosphorylated cytosolic Hog1p concentration, 4.443*1e-6, mM'"))
INIT_GLYCEROL_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_Glycerol", 576.0, "mM", f"{_S5} 'Init_Glycerol, Initial intracellular glycerol concentration, 576, mM'"))
INIT_CA_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_Ca", 5e-5, "mM", f"{_S5} 'Init_Ca, Initial Ca2+ concentration, 5*1e-5, mM'"))
INIT_CNOFF_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_CNoff", 1.1628e-3, "mM", f"{_S5} 'Init_CNoff, Initial inactive calcineurin concentration, 1.1628*1e-3, mM'"))
INIT_CN_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_CNon", 0.0, "mM", f"{_S5} 'Init_CNon, Initial active calcineurin concentration, 0, mM'"))
INIT_CRZ1_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_Crz1", 1.916e-4, "mM", f"{_S5} 'Init_Crz1, Initial Crz1p concentration, 1.916*1e-4, mM'"))
INIT_ENA1_MM = KE2013_PARAMS.add(Param.asserted(
    "Init_Ena1", 1.1364e-4, "mM", f"{_S5} 'Init_Ena1, Initial Ena1p protein concentration, 1.1364*1e-4, mM'"))

NOT_IN_SOURCE.add(Param.refused(
    "cytosolic_buffer_capacity", "mM per pH unit",
    "The proton balance Eqn 2.40 would need it to convert a proton flux into a pH change",
    reason="Ke 2013 does not contain the quantity in any form. The string 'buffer' appears "
           "ZERO times in the article XML and in all seven supplementary PDFs, and Eqn 2.40 "
           "is written on FREE protons -- dH+/dt = -J_Pma1 + 3/2 J_Nha1,Na + 3/2 J_Nha1,K + "
           "J_H_uptake + J_H_production -- with no buffering term at all. An omitted "
           "quantity is not a small quantity. Its absence limits physiological validity, "
           "but does not explain the reproduction discrepancy with a source that also "
           "omits buffering",
    missing="a cytosolic buffering capacity titrated in this strain and medium, e.g. by the "
            "pHluorin nigericin-clamp route this repository already uses for pH_c"))
NOT_IN_SOURCE.add(Param.refused(
    "anion_export_capacity", "amol/s",
    "The counter-anion flux that would let the 300 mM Anion row move at all",
    reason="Table S4's Anion is a CONSTANT with an empty reference column. There is no "
           "anion transport equation anywhere in Text S1, so the source cannot supply a "
           "capacity even in principle; charge balance is closed by assumption",
    missing="the same transport assay `mech/ph.py`'s vmax_anion_export refuses for; the two "
            "refusals are the same missing measurement seen from two models"))


@dataclass(frozen=True)
class TranscribedEquation:
    """One printed equation, what it says, and whether it runs here.

    ``printed`` is the equation as typeset in the pinned PDF, in plain ASCII and with
    nothing rearranged -- it is the audit trail back to the page. ``status`` is one of
    ``executable``, ``recorded`` (transcribed but not wired into any run here) or
    ``untranscribed`` (something needed is not in the source), and ``note`` says why.
    """

    label: str
    where: str
    printed: str
    status: str
    note: str = ""
    ode_count: int = 1

    def __post_init__(self) -> None:
        if int(self.ode_count) < 1:
            raise ValueError(f"{self.label}: an entry stands for at least one ODE")
        if self.status not in ("executable", "recorded", "untranscribed"):
            raise ValueError(f"{self.label}: status {self.status!r} is not one of "
                             "executable, recorded, untranscribed")
        if self.status == "untranscribed" and not self.note.strip():
            raise ValueError(f"{self.label}: an untranscribed equation must say what is "
                             "missing; that is the entire content of the status")


ODE_SYSTEM = (
    TranscribedEquation(
        "dH/dt", "Table S1; Text S1 Eqn 2.40",
        "d[H+]/dt = -J_Pma1 + (3/2)*J_Nha1,Na + (3/2)*J_Nha1,K + J_H_uptake + J_H_production",
        "executable",
        "on the AMOUNT of H+ (attomoles), not the concentration -- Text S1 2.5 is explicit. "
        "The 3/2 is Nha1p's 3 H+ per 2 K+/Na+. No buffering term exists in the source"),
    TranscribedEquation(
        "dNa/dt", "Table S1; Text S1 Eqn 2.40",
        "d[Na+]/dt = -J_Ena1,Na - J_Nha1,Na + J_Trk,Na + J_NSC1,Na + J_diff,Na",
        "executable"),
    TranscribedEquation(
        "dK/dt", "Table S1; Text S1 Eqn 2.40",
        "d[K+]/dt = -J_Ena1,K - J_Nha1,K + J_Trk1,K + J_NSC1,K - J_Tok1 + J_diff,K",
        "executable"),
    TranscribedEquation(
        "dVolume_cyt/dt", "Table S1; Text S1 Eqn 4.1",
        "d Volume_cyt/dt = -G_EK * Lp * D_Pressure",
        "recorded",
        "transcribed as turgor_pressure() and pressure_difference(), but NOT integrated in "
        "the reproduction run: at Table S5's own initial volume D_Pressure is -1.51e6 J/m^3 "
        "rather than 0, so the printed module does not rest where the paper starts it"),
    TranscribedEquation(
        "dPbs2/dt, dPbs2PP/dt", "Table S1; Text S1 Eqn 3.1",
        "d[Pbs2]/dt = -K_pho^Pbs2*[Pbs2]/(1+(Ptur/alphaHog1)^8) + K_depho^Pbs2*[Pbs2PP] "
        "- [Pbs2]*Vratio  (and the mirror image for [Pbs2PP])",
        "untranscribed", ode_count=2,
        note="needs Ptur from the volume module, which does not rest at Table S5, and Vratio, "
        "which is d(Volume)/dt / Volume. Adopted by Ke from Zi 2010 [50]"),
    TranscribedEquation(
        "dHog1c/dt, dHog1PPc/dt, dHog1n/dt, dHog1PPn/dt", "Table S1; Text S1 Eqn 3.1",
        "d[Hog1c]/dt = -K_pho^Hog1*[Pbs2PP]*[Hog1c] + K_depho^Hog1PPc*[Hog1PPc] "
        "- K_imp^Hog1c*[Hog1c] + K_exp^Hog1n*[Hog1n]*Vnuc/Vcyt - [Hog1c]*Vratio  (and the "
        "three companions)",
        "untranscribed", ode_count=4,
        note="Vcyt and Vnuc appear in all four equations and in NO table: not in S2, S3, S4 or "
        "S5. They live in the upstream Zi 2010 model [50], which this repository has not "
        "pinned. Guessing a nuclear volume fraction would be inventing a parameter"),
    TranscribedEquation(
        "dGlyint/dt, dYt/dt", "Table S1; Text S1 Eqn 3.1",
        "d[Glyint]/dt = Ks0^Glyc + Ks1^Glyc*(totalHog1PP)^4/(betaHog1^4+(totalHog1PP)^4) "
        "+ Ks2^Glyc*[Yt] - (Kexp0^Glyc + Kexp1^Glyc*(Ptur)^12/((gammaHog1)^12+(Ptur)^12)) "
        "*[Glyint] - [Glyint]*Vratio",
        "untranscribed", ode_count=2,
        note="downstream of the Hog1p states and of Ptur; totalHog1PP also carries the printed "
        "conversion (([Hog1PPc]*Vcyt+[Hog1PPn]*Vnuc)*602)/6780, which needs the same two "
        "absent volumes"),
    TranscribedEquation(
        "dz1..dz4/dt", "Table S1; Text S1 Eqn 3.1",
        "d[z1]/dt = 4*([Hog1PPn]-[z1])/tau ;  d[z_i]/dt = 4*([z_{i-1}]-[z_i])/tau",
        "untranscribed", ode_count=4,
        note="tau appears in four equations and in no parameter table. Table S3's Hog1 block "
        "lists every rate constant Zi 2010 [50] supplies EXCEPT this delay"),
    TranscribedEquation(
        "dCa/dt", "Table S1; Text S1 Eqn 3.3",
        "d[Ca2+]/dt = CCa - dCa*[Ca2+] + kCa,cyt*[Na+]cyt^h/( [Na+]cyt^h + KmCa,cyt^h ) "
        "+ kCa,ext*[Na+]ext^h/( [Na+]ext^h + KmCa,ext^h ) + kCa,pH*(pH_ext - 6.5) "
        "- 3*kCN,a*[Ca2+]^3*[CNoff] + 3*kCN,da*[CN] - [Ca2+]*Vratio",
        "executable",
        "the Hill exponents differ between the two saturating terms: h_Na_cyt = 12 for the "
        "cytosolic one, h_Na_ext = 2 for the external one, both asterisked as fitted"),
    TranscribedEquation(
        "dCNoff/dt, dCN/dt", "Table S1; Text S1 Eqn 3.4",
        "d[CNoff]/dt = -kCN,a*[Ca2+]^3*[CNoff] + kCN,da*[CN] + kCN_Ppz,da*[Ppz]*[CN] "
        "- [CNoff]*Vratio ;  d[CN]/dt = the negative of the first three terms - [CN]*Vratio",
        "executable", ode_count=2,
        note="Table S5's Init_CNon = 0 is NOT the steady state of this pair at Table S5's own "
        "Ca2+ and the Ppz that Eqn 2.21 gives; the block relaxes away from it"),
    TranscribedEquation(
        "dCrz1/dt", "Table S1; Text S1 Eqn 3.4",
        "d[Crz1]/dt = CCrz1 - dCrz1*[Crz1] + kCrz1*[CN] - [Crz1]*Vratio",
        "executable"),
    TranscribedEquation(
        "dNrg1/dt", "Table S1; Text S1 Eqn 3.5",
        "d[Nrg1]/dt = CNrg1*KmNrg1,pH/(KmNrg1,pH + k_pHext) - dNrg1*[Nrg1] - [Nrg1]*Vratio, "
        "with k_pHext = max(pH_Ext - 6, 0)",
        "untranscribed",
        "KmNrg1,pH is in the equation and in no table. Table S3's Nrg1p block lists CNrg1, "
        "dNrg1 and kNrg1,pH (5.25e-8 mM/s) -- a rate constant, in rate-constant units, "
        "which cannot stand in for a Michaelis constant in a dimensionless ratio"),
    TranscribedEquation(
        "dENA1mRNA/dt, dEna1/dt", "Table S1; Text S1 Eqn 3.6",
        "d[ENA1mRNA]/dt = CENA1,Nrg1*KmENA1,Nrg1^h/(KmENA1,Nrg1^h + [Nrg1]^h) "
        "+ kENA1,Crz1*[Crz1] + kENA1,Hog1*[z2] - dENA1mRNA*[ENA1mRNA] ; "
        "d[Ena1]/dt = ktEna1*[ENA1mRNA] - dEna1*[Ena1] - [Ena1]*Vratio",
        "untranscribed", ode_count=2,
        note="downstream of Nrg1p and of the Hog1p chain variable z2, both untranscribable. "
        "Table S3 also prints KmENA1,Nrg1 as '1.143.0*1e-4', which is not a number. The "
        "reproduction run therefore holds Ena1(t)/Ena1_0 at 1"),
)
"""Table S1 in full: TWENTY-THREE ODEs, grouped into thirteen entries.

Seven ODEs are executable here (H+, Na+, K+, Ca2+, CNoff, CN, Crz1), one more -- the volume
equation -- is transcribed but deliberately not integrated, and fifteen refuse. The fifteen
are not a shortfall in the reading. Every one needs a scalar that appears in an equation and
in no table, and the module docstring names each. :data:`ODE_TOTALS` does the arithmetic so
the counts in these docstrings cannot drift away from the data."""

ODE_TOTALS = {
    status: sum(eq.ode_count for eq in ODE_SYSTEM if eq.status == status)
    for status in ("executable", "recorded", "untranscribed")
}
"""7 executable, 1 recorded, 15 untranscribed -- 23, which is Table S1's own row count."""

STATE_VARIABLES = (
    ("h_amol", "amol", Compartment.CELL_WATER, "H+ amount; pH_i = -log10(h_amol/V * 1e-3)"),
    ("na_amol", "amol", Compartment.CELL_WATER, "Na+ amount"),
    ("k_amol", "amol", Compartment.CELL_WATER, "K+ amount"),
    ("volume_um3", "um^3", Compartment.CELL_WATER, "the osmotically changeable compartment"),
    ("ca_mM", "mM", Compartment.CELL_WATER, "cytosolic Ca2+, a CONCENTRATION in the source"),
    ("cn_off_mM", "mM", Compartment.REGULATORY, "inactive calcineurin"),
    ("cn_mM", "mM", Compartment.REGULATORY, "activated calcineurin"),
    ("crz1_mM", "mM", Compartment.REGULATORY, "Crz1p"),
)
"""Name, units, compartment and role, in `mech/contracts.py`'s vocabulary.

The mixed basis is the source's, not a slip: Text S1 2.5 tracks the three cations as
AMOUNTS because the fluxes are amounts per second, while every signalling species in
Tables S3 and S5 is a concentration. Anything reading this module must not mix them."""

READING_CHOICES = (
    {"equation": "Eqn 2.38, H+ uptake",
     "printed": "J_H_uptake = kH_uptake * (R*T*ln([H+]ext/[H+]int) - Em*F)",
     "problem": "the bracket is J/mol but Table S2 gives kH_uptake in 1e-18 mol/(s*V); the "
                "printed product is larger than every other flux by the factor F = 96485",
     "used": "J_H_uptake = kH_uptake * (R*T/F*ln([H+]ext/[H+]int) - Em)",
     "authority": "the paper's own unit column for kH_uptake",
     "cost": "at the printed-law equilibrium J_H_uptake and J_Pma1 are both about "
             "6.07e3 amol/s, with J_H_production = 5. The equilibrium has pH_i = 4.31, "
             "2.83 units below the published 7.14; it is not a physiological prediction"},
    {"equation": "Tok1p gating rates, Text S1 2.2.4 (unnumbered, after Eqn 2.18)",
     "printed": "k'_Tok1,1 = kTok1,1*exp(-lTok1,ext*Em/(R*T)), and the three companions",
     "problem": "Em is in volts and R*T is in J/mol, so the exponent is ~1e-5 and the gating "
                "is voltage-independent -- contradicting the same page's 'the rates of these "
                "transitions are governed by membrane potential'",
     "used": "the same expressions with F*Em/(R*T), matching Eqn 2.19's own GHK exponent",
     "authority": "Eqn 2.19 three paragraphs later, and the Loukin & Saimi [38] model Ke "
                  "says it is modifying",
     "cost": "27 mV. As printed: Em = -138.5 mV, pH_i = 7.252. As used: -110.9 mV, 7.267. "
             "Both are computed and both are reported"},
    {"equation": "Eqn 2.8, K+ flux through Ena1p",
     "printed": "J_Ena1,K = P_Ena1,K*(Ena1(t)/Ena1_0)*kEna1*([Na+]cyt*[ATP] - "
                "[Na+]ext*[ADP]*[Pi]*exp((G0_ATP-F*Em)/(R*T)))",
     "problem": "the second line expands a K+ flux with the Na+ gradient",
     "used": "the same with [K+], as J0_Ena1,K is defined one paragraph earlier",
     "authority": "Eqn 2.8's own first line, J_Ena1,K = P_Ena1,K*(Ena1(t)/Ena1_0)*J0_Ena1,K",
     "cost": "not separately quantified; a K+ flux driven by the Na+ gradient is a "
             "typesetting slip rather than an alternative model"},
)
"""Three printed forms that do not survive dimensional analysis, verbatim, with the reading
used and what it costs. Choices 1 and 2 are switchable through :class:`Reading`."""

UNTRANSCRIBED = tuple(eq for eq in ODE_SYSTEM if eq.status == "untranscribed")
"""Six entries, fifteen ODEs, each naming the scalar that is absent from every table."""


@dataclass(frozen=True)
class Reading:
    """Which reading of the three ambiguous printed forms a run uses.

    Defaults are :data:`UNIT_CONSISTENT`. Nothing here silently corrects the paper: every
    field names a printed form, and :data:`AS_PRINTED` reproduces the page exactly so the
    cost of each correction is a number this module computes rather than a claim it makes.

    Args:
        h_uptake_per_faraday: True divides Eqn 2.38's bracket by F, giving volts to match
            Table S2's unit for kH_uptake. False is the page.
        tok1_gating_faraday: True puts F into the Tok1p gating exponents. False is the page,
            and makes the gating voltage-independent.
        nha1_low_rate / nha1_high_rate: Eqns 2.16-2.17 need two rate constants where Table
            S2 prints one. Both default to kNha1 and the sensitivity is measured, not assumed.
    """

    h_uptake_per_faraday: bool = True
    tok1_gating_faraday: bool = True
    nha1_low_rate: float = float(K_NHA1)
    nha1_high_rate: float = float(K_NHA1)

    def __post_init__(self) -> None:
        finite(self.nha1_low_rate, "nha1_low_rate", minimum=0.0)
        finite(self.nha1_high_rate, "nha1_high_rate", minimum=0.0)


UNIT_CONSISTENT = Reading()
AS_PRINTED = Reading(h_uptake_per_faraday=False, tok1_gating_faraday=False)


@dataclass(frozen=True)
class Environment:
    """The medium, in the three quantities Ke's flux laws read.

    Ke's unstressed condition, stated in Results: "the external Na+ and K+ concentrations
    and the external pH are assumed to be 5 mM, 1 mM, and pH 6.5, respectively".
    """

    ph_ext: float = 6.5
    na_ext_mM: float = 5.0
    k_ext_mM: float = 1.0

    def __post_init__(self) -> None:
        finite(self.ph_ext, "ph_ext", minimum=0.0, maximum=14.0)
        finite(self.na_ext_mM, "na_ext_mM", minimum=0.0)
        finite(self.k_ext_mM, "k_ext_mM", minimum=0.0)

    @property
    def h_ext_mM(self) -> float:
        return 10.0 ** (-self.ph_ext) * 1e3


@dataclass(frozen=True)
class IonState:
    """The three cation amounts and the compartment volume they live in."""

    h_amol: float
    na_amol: float
    k_amol: float
    volume_um3: float = float(INIT_VOLUME_UM3)

    def __post_init__(self) -> None:
        # Solver output is numpy-typed; coerce so downstream JSON and bools stay native.
        for field in ("h_amol", "na_amol", "k_amol"):
            object.__setattr__(self, field, finite(getattr(self, field), field, positive=True))
        object.__setattr__(self, "volume_um3",
                           finite(self.volume_um3, "volume_um3", positive=True))

    @property
    def h_mM(self) -> float:
        return self.h_amol / self.volume_um3

    @property
    def na_mM(self) -> float:
        return self.na_amol / self.volume_um3

    @property
    def k_mM(self) -> float:
        return self.k_amol / self.volume_um3

    @property
    def ph_i(self) -> float:
        return -math.log10(self.h_mM * 1e-3)

    @classmethod
    def from_table_s5(cls) -> IonState:
        """Table S5's initial concentrations, converted to amounts at its initial volume."""
        volume = float(INIT_VOLUME_UM3)
        return cls(float(INIT_H_MM) * volume, float(INIT_NA_MM) * volume,
                   float(INIT_K_MM) * volume, volume)


_R = float(R_GAS)
_T = float(TEMPERATURE_K)
_RT = _R * _T
_F = float(FARADAY)
_CAPACITANCE_F = float(MEMBRANE_CAPACITANCE) * float(CELL_SURFACE_UM2) * 1e-12
_AMOL = 1e-18


def membrane_potential_V(state: IonState, anion_mM: float = float(ANION_MM)) -> float:
    """Eqn 2.1. Em = (H+ + K+ + Na+ - Anion) * F / (Cm * S_Mem), on AMOUNTS.

    The paper writes ``Anion_cyt`` as a concentration in Table S4 and subtracts it from a sum
    of amounts in Eqn 2.1; it is converted here at the current volume. Under the fixed-volume
    reproduction the two readings coincide exactly, which is why this ambiguity costs nothing
    for the comparison actually made -- and it is recorded rather than resolved.
    """
    net_amol = math.fsum((state.h_amol, state.na_amol, state.k_amol,
                         -anion_mM * state.volume_um3))
    return net_amol * _AMOL * _F / _CAPACITANCE_F


def _boltzmann(em_V: float) -> float:
    """The ``exp((G0_ATP - F*Em)/(R*T))`` factor shared by Eqns 2.7, 2.8, 2.15 and 2.37."""
    return math.exp((float(G0_ATP) - _F * em_V) / _RT)


def fluxes(state: IonState, env: Environment, *, hog1ppc_mM: float = float(INIT_HOG1PPC_MM),
           cn_mM: float = float(INIT_CN_MM), ena1_ratio: float = 1.0,
           reading: Reading = UNIT_CONSISTENT, anion_mM: float = float(ANION_MM)) -> dict:
    """Every transporter flux of Text S1 sections 2.1-2.4, in amol/s, at one state.

    Returns a mapping whose keys are the flux names Table S1's ODEs use, plus ``Em`` (V),
    ``P_Tok1_O`` and ``Ppz`` for inspection. Positive means the direction Ke's ODE signs
    assume, so ``J_Pma1`` positive is protons LEAVING the cell.
    """
    h_c, na_c, k_c = state.h_mM, state.na_mM, state.k_mM
    em = membrane_potential_V(state, anion_mM)
    boltz = _boltzmann(em)
    z = em * _F / _RT

    # Ena1p: Eqns 2.5-2.8. The K+ line uses [K+], see READING_CHOICES entry 3.
    ena_den = 1.0 + k_c / float(KM_ENA1_K) + na_c / float(KM_ENA1_NA)
    drive = float(ADP_MM) * float(PI_MM) * boltz
    j_ena1_na = ((na_c / float(KM_ENA1_NA)) / ena_den * ena1_ratio * float(K_ENA1)
                 * (na_c * float(ATP_MM) - env.na_ext_mM * drive))
    j_ena1_k = ((k_c / float(KM_ENA1_K)) / ena_den * ena1_ratio * float(K_ENA1)
                * (k_c * float(ATP_MM) - env.k_ext_mM * drive))

    # Nha1p: Eqns 2.9-2.17. Hog1p phosphorylation moves it to the low-K+-affinity state.
    p_low = hog1ppc_mM / (float(KM_NHA1_HOG1) + hog1ppc_mM)
    p_high = float(KM_NHA1_HOG1) / (float(KM_NHA1_HOG1) + hog1ppc_mM)
    den_hi = 1.0 + k_c / float(KM_NHA1_HIGH_K) + na_c / float(KM_NHA1_NA)
    den_lo = 1.0 + k_c / float(KM_NHA1_LOW_K) + na_c / float(KM_NHA1_NA)
    v_k = (reading.nha1_low_rate * p_low * (k_c / float(KM_NHA1_LOW_K)) / den_lo
           + reading.nha1_high_rate * p_high * (k_c / float(KM_NHA1_HIGH_K)) / den_hi)
    v_na = (reading.nha1_low_rate * p_low * (na_c / float(KM_NHA1_NA)) / den_lo
            + reading.nha1_high_rate * p_high * (na_c / float(KM_NHA1_NA)) / den_hi)
    j_nha1_k = v_k * (env.h_ext_mM ** 3 * k_c ** 2 - h_c ** 3 * env.k_ext_mM ** 2 * boltz)
    j_nha1_na = v_na * (env.h_ext_mM ** 3 * na_c ** 2 - h_c ** 3 * env.na_ext_mM ** 2 * boltz)

    # Tok1p: Eqns 2.18-2.20. The four-state distribution solves exactly with R_I = 1.
    gate = _F if reading.tok1_gating_faraday else 1.0
    r_i = 1.0
    occ_o = float(K_TOK1_RO) / float(K_TOK1_OR) * r_i
    k1p = float(K_TOK1_1) * math.exp(-float(L_TOK1_EXT) * em * gate / _RT)
    k1m = float(ALPHA_TOK1) * float(K_TOK1_1) * math.exp(float(L_TOK1_EXT) * em * gate / _RT)
    k2p = float(BETA_TOK1) * float(K_TOK1_2) * math.exp(-float(L_TOK1_INT) * em * gate / _RT)
    k2m = float(K_TOK1_2) * math.exp(float(L_TOK1_INT) * em * gate / _RT)
    r_io = k1p * env.k_ext_mM / k1m
    r_o = k2p * r_io / (k2m * k_c)
    p_open = occ_o / (occ_o + r_i + r_io + r_o)
    if z > 0:
        ghk = z * (k_c - math.exp(-z) * env.k_ext_mM) / -math.expm1(-z)
    elif z < 0:
        ghk = z * (k_c * math.exp(z) - env.k_ext_mM) / math.expm1(z)
    else:
        ghk = k_c - env.k_ext_mM
    i_tok1 = float(PS_TOK1) * _F * ghk
    j_tok1 = (float(KM_TOK1_HOG1) / (hog1ppc_mM + float(KM_TOK1_HOG1))
              * p_open * i_tok1 / _F)

    # Ppz and the Trk system: Eqns 2.21-2.28. Ppz is the pH-sensing branch of the model.
    ppz = float(KM_PPZ) / (float(KM_PPZ) + h_c) * float(PPZ_TOTAL)
    trk_den = ppz / float(KM_TRK_PPZ) + cn_mM / float(KM_TRK_CN) + 1.0
    p_medium = (ppz / float(KM_TRK_PPZ)) / trk_den
    p_high_aff = (cn_mM / float(KM_TRK_CN) + 1.0) / trk_den
    den_h = 1.0 + env.k_ext_mM / float(KM_TRK_HIGH_K) + env.na_ext_mM / float(KM_TRK_HIGH_NA)
    den_m = 1.0 + env.k_ext_mM / float(KM_TRK_MEDIUM_K) + env.na_ext_mM / float(KM_TRK_MEDIUM_NA)
    j_trk_k = (((env.k_ext_mM / float(KM_TRK_MEDIUM_K)) / den_m * p_medium
                + (env.k_ext_mM / float(KM_TRK_HIGH_K)) / den_h * p_high_aff)
               * float(K_TRK) * em ** 2)
    j_trk_na = (((env.na_ext_mM / float(KM_TRK_MEDIUM_NA)) / den_m * p_medium
                 + (env.na_ext_mM / float(KM_TRK_HIGH_NA)) / den_h * p_high_aff)
                * float(K_TRK) * em ** 2)

    # NSC1: Eqns 2.29-2.32, and passive leakage: Eqn 2.39. H+ does not leak in this model.
    nsc_den = 1.0 + env.k_ext_mM / float(KM_NSC1_K) + env.na_ext_mM / float(KM_NSC1_NA)
    ez = math.exp(z)
    j_nsc1_k = ((env.k_ext_mM / float(KM_NSC1_K)) / nsc_den * float(K_NSC1)
                * (env.k_ext_mM - k_c * ez))
    j_nsc1_na = ((env.na_ext_mM / float(KM_NSC1_NA)) / nsc_den * float(K_NSC1)
                 * (env.na_ext_mM - na_c * ez))
    j_diff_k = float(PERMEABILITY_K) * float(CELL_SURFACE_UM2) * (env.k_ext_mM - k_c * ez)
    j_diff_na = float(PERMEABILITY_NA) * float(CELL_SURFACE_UM2) * (env.na_ext_mM - na_c * ez)

    # H+ production, extrusion and uptake: Eqns 2.33, 2.37, 2.38.
    j_pma1 = float(K_PMA1) * (h_c * float(ATP_MM) - env.h_ext_mM * drive)
    gradient = math.log(env.h_ext_mM / h_c)
    if reading.h_uptake_per_faraday:
        j_h_uptake = float(K_H_UPTAKE) * (_RT / _F * gradient - em)
    else:
        j_h_uptake = float(K_H_UPTAKE) * (_RT * gradient - em * _F)

    return {"Em": em, "J_Ena1_Na": j_ena1_na, "J_Ena1_K": j_ena1_k, "J_Nha1_Na": j_nha1_na,
            "J_Nha1_K": j_nha1_k, "J_Tok1": j_tok1, "J_Trk_Na": j_trk_na, "J_Trk_K": j_trk_k,
            "J_NSC1_Na": j_nsc1_na, "J_NSC1_K": j_nsc1_k, "J_diff_Na": j_diff_na,
            "J_diff_K": j_diff_k, "J_Pma1": j_pma1, "J_H_uptake": j_h_uptake,
            "J_H_production": float(K_H_PROD), "P_Tok1_O": p_open, "Ppz": ppz}


def _ion_balance_terms(f):
    return (
        (-f["J_Pma1"], 1.5 * f["J_Nha1_Na"], 1.5 * f["J_Nha1_K"],
         f["J_H_uptake"], f["J_H_production"]),
        (-f["J_Ena1_Na"], -f["J_Nha1_Na"], f["J_Trk_Na"], f["J_NSC1_Na"], f["J_diff_Na"]),
        (-f["J_Ena1_K"], -f["J_Nha1_K"], f["J_Trk_K"], f["J_NSC1_K"], -f["J_Tok1"], f["J_diff_K"]),
    )


def ion_rhs(state: IonState, env: Environment, **kwargs) -> tuple[float, float, float]:
    """Eqn 2.40 / Table S1's first three rows, in amol/s. No buffering term: see the docstring."""
    return tuple(math.fsum(terms) for terms in _ion_balance_terms(fluxes(state, env, **kwargs)))


def calcineurin_rhs(ca_mM: float, cn_off_mM: float, cn_mM: float, crz1_mM: float, *,
                    ppz_mM: float, na_cyt_mM: float, env: Environment,
                    v_ratio_per_s: float = 0.0) -> tuple[float, float, float, float]:
    """Eqns 3.3-3.4: cytosolic Ca2+, the two calcineurin pools and Crz1p, in mM/s.

    ``v_ratio_per_s`` is Ke's ``Vratio`` = (dVolume/dt)/Volume, the dilution term every
    signalling ODE in Table S1 carries. It is 0 under the fixed-volume reduction and that is
    the ONLY setting used here, because the printed volume module does not rest at Table S5.
    The initial Ca2+ spike Eqn 3.2 imposes at stress onset is a separate, discontinuous
    reset and is deliberately not applied inside a right-hand side.
    """
    h_cyt = float(H_NA_CYT)
    h_ext = float(H_NA_EXT)
    cyt_term = (float(K_CA_CYT) * na_cyt_mM ** h_cyt
                / (na_cyt_mM ** h_cyt + float(KM_CA_CYT) ** h_cyt))
    ext_term = (float(K_CA_EXT) * env.na_ext_mM ** h_ext
                / (env.na_ext_mM ** h_ext + float(KM_CA_EXT) ** h_ext))
    activation = float(K_CN_A) * ca_mM ** 3 * cn_off_mM
    deactivation = float(K_CN_DA) * cn_mM + float(K_CN_PPZ_DA) * ppz_mM * cn_mM
    d_ca = (float(C_CA) - float(D_CA) * ca_mM + cyt_term + ext_term
            + float(K_CA_PH) * (env.ph_ext - 6.5)
            - 3.0 * activation + 3.0 * float(K_CN_DA) * cn_mM - ca_mM * v_ratio_per_s)
    d_cn_off = -activation + deactivation - cn_off_mM * v_ratio_per_s
    d_cn = activation - deactivation - cn_mM * v_ratio_per_s
    d_crz1 = (float(C_CRZ1) - float(D_CRZ1) * crz1_mM + float(K_CRZ1) * cn_mM
              - crz1_mM * v_ratio_per_s)
    return d_ca, d_cn_off, d_cn, d_crz1


def turgor_pressure(volume_um3: float) -> float:
    """Eqn 4.2, in J/m^3. Turgor falls as the compartment shrinks and floors at zero.

    Read from the rendered page, not from a text extractor: pdftotext returns the numerator
    with its superscript transposed, as ``Volume_cyt - Volume_cyt^0``, which reverses the
    sign of the whole response and makes the volume module unconditionally unstable.
    """
    v0 = float(INIT_VOLUME_UM3)
    span = v0 - float(R_VOL) * v0
    return max(float(TURGOR_0) * (1.0 - (v0 - volume_um3) / span), 0.0)


def pressure_difference(state: IonState, env: Environment = Environment(), *,
                        glycerol_mM: float = float(INIT_GLYCEROL_MM),
                        added_osmolarity_mM: float = 0.0,
                        na_ext_0_mM: float = 5.0, k_ext_0_mM: float = 1.0) -> float:
    """Eqns 4.3-4.4, in J/m^3. ``d Volume_cyt/dt = -G_EK * Lp * D_Pressure`` (Eqn 4.1).

    Evaluated at Table S5's own initial state this returns -1.51e6 rather than 0: the
    printed osmotic balance is 1176 mM inside against 250 mM outside, and a turgor of
    0.875e6 J/m^3 does not close a 926 mM gap. That is why the reproduction run holds the
    volume fixed, and the number is asserted in the test rather than described here.
    ``added_osmolarity_mM`` is the non-ionic osmolyte a sorbitol experiment would add.
    """
    v0 = float(INIT_VOLUME_UM3)
    osmo_ext = (float(OSMO_EXT_0) + added_osmolarity_mM - (na_ext_0_mM + k_ext_0_mM)
                + (env.na_ext_mM + env.k_ext_mM))
    osmo_cyt = ((float(OSMO_CYT_0) - (float(INIT_NA_MM) + float(INIT_K_MM))) * v0
                / state.volume_um3 + (state.na_mM + state.k_mM) + glycerol_mM)
    return (osmo_ext - osmo_cyt) * _RT + turgor_pressure(state.volume_um3)


def refuse_untranscribed(label: str) -> None:
    """Raise :class:`ScientificRefusal` for any equation :data:`UNTRANSCRIBED` names.

    The house idiom: a caller reaching for the Hog1p module or the ENA1 chain should hit the
    absent scalar by name rather than a plausible default.
    """
    for equation in UNTRANSCRIBED:
        if label == equation.label or label in equation.label:
            raise ScientificRefusal(
                f"Ke 2013 {equation.label} ({equation.where}) is transcribed but not "
                f"executable here: {equation.note}")
    raise ValueError(f"{label!r} is not an untranscribed Ke 2013 equation; "
                     f"choose from {[e.label for e in UNTRANSCRIBED]}")


STEADY_STATE_TOLERANCE = 1e-8
"""Absolute balance tolerance, in amol/s for each of H+, Na+ and K+."""
STEADY_STATE_RELATIVE_TOLERANCE = 1e-10
"""Relative tolerance on the sum of absolute flux terms in each ion balance."""


class EquilibriumNotConverged(ScientificRefusal):
    status = "numerical_nonconvergence"
    nonexistence_proven = False

    def __init__(self, env, best_scaled_residual):
        self.best_scaled_residual = best_scaled_residual
        super().__init__(
            f"Ke ion equilibrium did not converge at external pH {env.ph_ext}; "
            f"best tolerance-scaled balance residual {best_scaled_residual:.3g} exceeds one. "
            "This is a numerical failure, not proof of root nonexistence or biological lethality.")


def steady_state(env: Environment, *, reading: Reading = UNIT_CONSISTENT,
                 hog1ppc_mM: float = float(INIT_HOG1PPC_MM),
                 cn_mM: float = float(INIT_CN_MM), ena1_ratio: float = 1.0,
                 volume_um3: float | None = None, guess: IonState | None = None) -> IonState:
    """Find an equilibrium of the explicitly reduced, fixed-volume ion equations.

    Log ion fractions and a voltage-scaled log total keep amounts positive without
    asking two large, nearly cancelling alkali amounts to determine the Newton step.
    All three acceptance residuals are amount fluxes, with absolute and relative
    tolerances; the optimizer's success flag does not establish an equilibrium.
    A returned root is not a validated biological state or a proof of uniqueness.
    """
    if not isinstance(env, Environment) or not isinstance(reading, Reading):
        raise ValueError("env and reading require Environment and Reading instances")
    if (env.na_ext_mM, env.k_ext_mM) != (5.0, 1.0):
        raise ScientificRefusal("the fixed-volume, fixed-regulator reproduction supports pH changes at the source's Na=5 mM and K=1 mM only; osmotic/ionic protocols require the missing coupled volume and signalling dynamics")
    if guess is not None and not isinstance(guess, IonState):
        raise ValueError("guess must be an IonState")
    volume = finite(float(INIT_VOLUME_UM3) if volume_um3 is None else volume_um3,
                    "volume_um3", positive=True)
    kwargs = dict(hog1ppc_mM=finite(hog1ppc_mM, "hog1ppc_mM", minimum=0.0),
                  cn_mM=finite(cn_mM, "cn_mM", minimum=0.0),
                  ena1_ratio=finite(ena1_ratio, "ena1_ratio", minimum=0.0), reading=reading)
    anion = float(ANION_MM) * volume
    charge_per_volt = _CAPACITANCE_F / (_AMOL * _F)

    def state_from_coordinates(x):
        total = anion * math.exp(charge_per_volt * x[2] / anion)
        alkali = total * expit(-x[0])
        return IonState(total * expit(x[0]), alkali * expit(x[1]),
                        alkali * expit(-x[1]), volume)

    def balances(state):
        terms = _ion_balance_terms(fluxes(state, env, **kwargs))
        return (np.array([math.fsum(row) for row in terms]),
                np.array([math.fsum(abs(value) for value in row) for row in terms]))

    starts = ([guess] if guess is not None else []) + [IonState.from_table_s5()]
    best = math.inf
    for seed in starts:
        amounts = np.array([seed.h_amol, seed.na_amol, seed.k_amol]) * volume / seed.volume_um3
        start = [math.log(amounts[0] / (amounts[1] + amounts[2])),
                 math.log(amounts[1] / amounts[2]),
                 anion / charge_per_volt * math.log(math.fsum(amounts) / anion)]
        try:
            _, scale = balances(state_from_coordinates(start))
            scale = np.maximum(scale, 1.0)
            solved = least_squares(
                lambda x: balances(state_from_coordinates(x))[0] / scale,
                start, x_scale=[1.0, 1.0, 0.1], max_nfev=2000,
                ftol=1e-13, xtol=1e-13, gtol=1e-13)
            state = state_from_coordinates(solved.x)
            residual, gross = balances(state)
            tolerance = STEADY_STATE_TOLERANCE + STEADY_STATE_RELATIVE_TOLERANCE * gross
            score = float(np.max(np.abs(residual) / tolerance))
        except (ValueError, OverflowError, FloatingPointError):
            continue
        if math.isfinite(score):
            best = min(best, score)
            if score <= 1.0:
                return state
    raise EquilibriumNotConverged(env, best)


def alkaline_step(*, ph_ext_from: float = 6.5, ph_ext_to: float = 8.0,
                  minutes: float = 200.0, reading: Reading = UNIT_CONSISTENT) -> dict:
    """Ke's Fig. 7E: sit at the unstressed steady state, step external pH, read pH_i.

    The step is iso-osmotic, which is the one perturbation the fixed-volume reduction is
    legitimate for. Returns the before and after states and the shift, in the units Ke
    prints them in.
    """
    before_env = Environment(ph_ext=ph_ext_from)
    after_env = Environment(ph_ext=ph_ext_to)
    before = steady_state(before_env, reading=reading)
    seconds = finite(minutes, "minutes", positive=True) * 60.0

    def rhs(_t, y):
        return list(ion_rhs(IonState(y[0], y[1], y[2], before.volume_um3), after_env,
                            reading=reading))

    solved = solve_ivp(rhs, (0.0, seconds),
                       [before.h_amol, before.na_amol, before.k_amol],
                       method="LSODA", rtol=1e-9, atol=1e-13)
    if not solved.success:
        raise ScientificRefusal(f"the alkaline step did not integrate: {solved.message}")
    after = IonState(*solved.y[:, -1], before.volume_um3)
    return {"before": before, "after": after, "minutes": minutes,
            "ph_i_before": before.ph_i, "ph_i_after": after.ph_i,
            "ph_i_shift": after.ph_i - before.ph_i,
            "em_mV_before": membrane_potential_V(before) * 1e3,
            "em_mV_after": membrane_potential_V(after) * 1e3}


def calcineurin_steady_state(env: Environment, ion: IonState) -> dict:
    """Where Eqns 3.3-3.4 settle, against Table S5's stated unstressed steady state.

    Table S5's footnote claims its initial concentrations ARE the unstressed steady state.
    For the calcineurin block they are not, and this returns both so the gap is a number.
    """
    ppz = fluxes(ion, env)["Ppz"]
    total_cn = float(INIT_CNOFF_MM) + float(INIT_CN_MM)

    def activated(ca: float) -> float:
        """Eqn 3.4 at rest, using the conservation d[CNoff]/dt + d[CN]/dt = 0 at Vratio = 0."""
        on = float(K_CN_A) * ca ** 3
        off = float(K_CN_DA) + float(K_CN_PPZ_DA) * ppz
        return total_cn * on / (on + off)

    def residual(ca: float) -> float:
        d_ca, _d_off, _d_cn, _d_crz = calcineurin_rhs(
            ca, total_cn - activated(ca), activated(ca), float(INIT_CRZ1_MM), ppz_mM=ppz,
            na_cyt_mM=ion.na_mM, env=env)
        return d_ca

    ca = brentq(residual, 1e-9, 1e-1, xtol=1e-18, rtol=1e-14)
    cn = activated(ca)
    crz1 = (float(C_CRZ1) + float(K_CRZ1) * cn) / float(D_CRZ1)
    return {"converged": bool(abs(residual(ca)) < 1e-12), "ppz_mM": ppz,
            "ca_mM": ca, "cn_mM": cn, "crz1_mM": crz1,
            "table_s5_ca_mM": float(INIT_CA_MM), "table_s5_cn_mM": float(INIT_CN_MM),
            "table_s5_crz1_mM": float(INIT_CRZ1_MM),
            "crz1_fold_vs_table_s5": crz1 / float(INIT_CRZ1_MM)}


_ACIDIC_LIMIT_REASON = (
    "A physiological acidic limit is not identified by these reduced ion equations. "
    "The former cutoff was a root-finder convergence boundary, not proof of root "
    "nonexistence. Buffering, coupled regulation and a validated viability criterion "
    "remain missing; changing numerical search resolution cannot supply them.")


def acidic_limit(reading: Reading = UNIT_CONSISTENT, *, resolution: float = 0.01) -> dict:
    """Refuse the former convergence-derived cutoff; no physiological limit is identified."""
    raise ScientificRefusal(_ACIDIC_LIMIT_REASON)


def reproduction_report() -> dict:
    """Every comparison this transcription is scored on, computed rather than asserted.

    This is the function that decides whether the transcription has earned its place. It
    reports the mismatches as well as the matches; ``tests/test_ke2013_transcription.py``
    asserts the numbers, and none of them was tuned.
    """
    unstressed = Environment()
    rows = {}
    for name, reading in (("unit_consistent", UNIT_CONSISTENT), ("as_printed_tok1",
                                                                Reading(tok1_gating_faraday=False))):
        state = steady_state(unstressed, reading=reading)
        step = alkaline_step(reading=reading)
        rows[name] = {
            "ph_i": state.ph_i, "em_mV": membrane_potential_V(state) * 1e3,
            "na_mM": state.na_mM, "k_mM": state.k_mM,
            "ph_i_alkaline_200min": step["ph_i_after"], "ph_i_shift": step["ph_i_shift"],
            "em_mV_alkaline_200min": step["em_mV_after"],
        }

    here = rows["unit_consistent"]
    comparison = {
        "unstressed_ph_i": (PUBLISHED["unstressed_ph_i"], here["ph_i"]),
        "unstressed_em_mV": (PUBLISHED["unstressed_em_mV"], here["em_mV"]),
        "unstressed_na_mM": (PUBLISHED["unstressed_na_mM"], here["na_mM"]),
        "unstressed_k_mM": (PUBLISHED["unstressed_k_mM"], here["k_mM"]),
        "alkaline_ph_i": (PUBLISHED["alkaline_ph_i"], here["ph_i_alkaline_200min"]),
        "alkaline_ph_i_shift": (PUBLISHED["alkaline_ph_i"] - PUBLISHED["unstressed_ph_i"],
                                here["ph_i_shift"]),
    }

    sensitivity = {}
    for label, reading in (("nha1_low_rate_zero", Reading(nha1_low_rate=0.0)),
                           ("nha1_low_rate_10x", Reading(nha1_low_rate=10.0 * float(K_NHA1)))):
        state = steady_state(unstressed, reading=reading)
        sensitivity[label] = {"ph_i": state.ph_i, "em_mV": membrane_potential_V(state) * 1e3}
    sensitivity["nha1_split_ph_span"] = abs(
        sensitivity["nha1_low_rate_zero"]["ph_i"] - sensitivity["nha1_low_rate_10x"]["ph_i"])

    as_printed = {"has_fixed_point": True}
    try:
        printed_state = steady_state(unstressed, reading=AS_PRINTED)
        printed_flux = fluxes(printed_state, unstressed, reading=AS_PRINTED)
        as_printed.update({
            "ph_i": printed_state.ph_i,
            "em_mV": membrane_potential_V(printed_state) * 1e3,
            "J_H_uptake_amol_per_s": printed_flux["J_H_uptake"],
            "J_Pma1_amol_per_s": printed_flux["J_Pma1"],
            "ph_i_below_published": PUBLISHED["unstressed_ph_i"] - printed_state.ph_i,
        })
    except (ScientificRefusal, OverflowError) as exc:
        as_printed = {"has_fixed_point": None, "status": "numerically_unresolved",
                      "why": str(exc)}

    table_s5 = IonState.from_table_s5()
    acidic = steady_state(Environment(ph_ext=5.5))
    return {
        "readings": rows,
        "acidic_counterexample": {
            "external_ph": 5.5, "intracellular_ph": acidic.ph_i,
            "kind": "reduced-model equilibrium, not a physiological prediction",
            "biological_validation": False,
        },
        "published_vs_here": comparison,
        "sensitivity": sensitivity,
        "both_forms_as_printed": as_printed,
        "volume_module": {
            "d_pressure_at_table_s5_J_per_m3": pressure_difference(table_s5),
            "turgor_at_table_s5_J_per_m3": turgor_pressure(table_s5.volume_um3),
            "note": "Eqn 4.1 would move the volume at Table S5's own initial state; the "
                    "printed module does not rest where the paper starts it",
        },
        "calcineurin_module": calcineurin_steady_state(
            unstressed, steady_state(unstressed)),
        "acidic_limit": {"status": "refused", "limit_ph_ext": None,
                         "reason": _ACIDIC_LIMIT_REASON},
        "numerical_contract": {
            "balance_units": "amol/s", "absolute_tolerance": STEADY_STATE_TOLERANCE,
            "relative_tolerance": STEADY_STATE_RELATIVE_TOLERANCE,
            "relative_scale": "sum of absolute ion-balance flux terms",
            "source_parameters_refitted": False, "biological_validation": False,
            "failure_semantics": "nonconvergence is unresolved, not proof of root nonexistence",
        },
    }


def transcription_record() -> dict:
    """The whole transcription as a JSON-shaped record, for `data/transcriptions/ke2013.json`.

    Includes the fragment a later merge drops into `data/parameter_evidence.json`. That
    fragment records ``availability = "transcribed_here"`` and never ``local_verified``.
    """
    report = reproduction_report()
    return {
        "schema_version": 1,
        "transcription_id": "ke2013",
        "availability": AVAILABILITY,
        "availability_meaning": (
            "the equations are upstream and byte-checksummed; the executable file is this "
            "repository's. Never promote to local_verified, which means a checksummed file "
            "the upstream authors published"),
        "source": SOURCE,
        "license": "CC-BY-4.0; attribution retained, equations and parameter values are Ke's",
        "transcribed_by": "this repository, src/ystwin/mech/ph_ke2013.py",
        "method": (
            "Tables S1-S5 read from the pinned PDFs. Text extraction alone was not "
            "sufficient and was not trusted: Eqns 2.1, 2.34 and 2.36 return empty from "
            "pdftotext, and Eqn 4.2 returns with its superscript transposed, which reverses "
            "the sign of the turgor response. Every equation used here was re-read from the "
            "page rendered at 190-400 dpi"),
        "artifacts": {role: {"path": path, "sha256": digest}
                      for role, (path, digest) in ARTIFACTS.items()},
        "ode_totals": dict(ODE_TOTALS),
        "equations": [
            {"label": eq.label, "where": eq.where, "printed": eq.printed,
             "status": eq.status, "ode_count": eq.ode_count, "note": eq.note}
            for eq in ODE_SYSTEM],
        "reading_choices": list(READING_CHOICES),
        "state_variables": [
            {"name": name, "units": units, "compartment": compartment.value, "role": role}
            for name, units, compartment, role in STATE_VARIABLES],
        "parameters": [
            {"name": p.name, "value": p.value, "units": p.units, "tag": p.tag,
             "provenance_verbatim": p.source} for p in KE2013_PARAMS.params],
        "absent_from_source": [
            {"name": p.name, "units": p.units, "tag": p.tag, "reason": p.reason,
             "missing": p.missing} for p in NOT_IN_SOURCE.params],
        "published_targets": dict(PUBLISHED),
        "reproduction": report,
        "parameter_evidence_fragment": {
            "models": {
                "ke2013_transcription": {
                    "source_id": "ke2013",
                    "programmes": ["ph"],
                    "availability": AVAILABILITY,
                    "equation_asset": ARTIFACTS["table_s1_odes"][0],
                    "executable_asset": "src/ystwin/mech/ph_ke2013.py",
                    "gap": (
                        "A hand transcription made in this repository, not an upstream "
                        "artifact. Seven of Table S1's twenty-three ODEs run and one more "
                        "(the volume equation) is transcribed but not integrated; fifteen "
                        "refuse, each because a scalar in the equation is in no table "
                        "(Vcyt, Vnuc and tau for the Hog1p module; KmNrg1,pH for Nrg1p, and "
                        "so the ENA1 chain downstream of it). Three printed forms fail "
                        "dimensional analysis and both readings of each are computed. It "
                        "reproduces the direction and 68% of the magnitude of Ke's alkaline "
                        "pH_i shift and sits 0.127 pH units and 17 mV off the printed "
                        "unstressed levels; nothing was fitted to close that. The two "
                        "scalars this programme needs are still absent: Ke has no buffering "
                        "term at all and closes charge balance with an unmeasured 300 mM "
                        "anion"),
                }
            }
        },
    }


def write_transcription_record(path: Path | str = RECORD_PATH) -> Path:
    """Recompute the record and write it. The test re-runs this and diffs, so the JSON on
    disk can never drift away from the module that produced it."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(transcription_record(), indent=2, sort_keys=False) + "\n")
    return target


def provenance() -> "object":
    """The graded parameter table, in the same shape every other mech module reports it."""
    return provenance_table(KE2013_PARAMS.params + NOT_IN_SOURCE.params)


if __name__ == "__main__":  # pragma: no cover
    print(json.dumps(reproduction_report(), indent=2, default=float))
