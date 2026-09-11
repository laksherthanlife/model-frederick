"""Ke 2013's calcineurin/Crz1p block, WIRED as a stress axis. Nothing here is re-transcribed.

WHAT THIS FILE IS. `mech/ph_ke2013.py` already carries a hand transcription of Ke, Ingram &
Haynes 2013 (PLoS Comput Biol 9: e1002879, PMC3547829, CC-BY 4.0), including Eqns 3.3-3.4 --
cytosolic Ca2+, the two calcineurin pools and Crz1p -- as :func:`~mech.ph_ke2013.calcineurin_rhs`,
with all fifteen Table S3 rate constants and four Table S5 initial conditions already registered
in ``KE2013_PARAMS``. Nothing in ``src/`` imported it. This module IMPORTS that block and gives
it the shape a stress axis needs: a state vector, per-hour derivatives, an engine-facing driver
and the back-coupling into the ion system. **Every rate law and every number is imported, not
retyped.** Two small pieces of Ke algebra that live inline inside
:func:`~mech.ph_ke2013.fluxes` and are not exported -- Eqn 2.21 (Ppz) and Eqn 2.24's Trk
affinity split -- are rewritten here, and ``tests/test_calcium_ke2013_transcription.py`` asserts
both against the audited implementation rather than trusting the retype. So: zero new vendoring,
zero new licence exposure, zero new transcription risk, one new pinned-source dependency
(the same eight files ``ph_ke2013.ARTIFACTS`` already hashes).

WHAT REPRODUCES, stated before the code because it is the result. Ke publishes the calcium
block's rate constants in Table S3 (pinned ``s011.pdf``) and its initial concentrations in
Table S5 (pinned ``s013.pdf``), two separate supplementary files, and Table S5's footnote states
those initials ARE the unstressed steady state. Two of the three rest points that Table S3 alone
implies land on Table S5's own rows:

    quantity              from Table S3            Table S5 row     difference
    [Ca2+]  = CCa/dCa     5.000000e-05 mM          5.0e-05 mM       0.00000%   (exact)
    [Crz1p] = CCrz1/dCrz1 1.918537e-04 mM          1.916e-04 mM     +0.13%
    [CN]                  1.20e-05 mM (relaxed)    0 mM             does NOT reproduce

That is a genuine independent check on the reading: nothing was fitted, the two tables were read
off different PDFs, and a mistyped constant would not land on a third file's number to five
figures. It is a check on the TRANSCRIPTION, not evidence that the model is right.

WHAT DOES NOT REPRODUCE, CARRIED FORWARD UNCHANGED AND NOT TUNED. Table S5's ``Init_CNon = 0``
is not the rest point of Eqn 3.4 at Table S5's own Ca2+ and the Ppz that Eqn 2.21 gives: at the
published initials ``d[CNoff]/dt = -1.45e-07 mM/s`` rather than zero, the block relaxes to
1.20e-05 mM activated calcineurin, and Crz1p follows to 2.54x its Table S5 row. Note the sharper
form of the contradiction that this module's cross-check exposes: Table S5's Crz1p row IS
exactly ``CCrz1/dCrz1``, i.e. the Crz1p you get at CN = 0, so Ke plainly intended CN = 0 as the
unstressed rest point -- and Eqn 3.4 as printed refuses it. This is the SOURCE's inconsistency.
It is already regression-tested at ``tests/test_ke2013_transcription.py:206``, it is re-asserted
here, and nothing in this file may be adjusted to close it.

THE DOMAIN FINDING, WHICH THE WIRING AGENT MUST READ BEFORE TOUCHING ``engine.py``. Eqn 3.3's
alkaline driver is LINEAR and unbounded below: ``+ kCa,pH*(pH_ext - 6.5)``, with
kCa,pH = 1.5e-05 mM/s. Below external pH 6.5 it is a SINK, and it is six times larger per pH
unit than CCa, the only zero-order source. Setting Ca2+ to zero and asking where total
production vanishes gives

    pH_floor = 6.5 - (CCa + cyt_Na_term + ext_Na_term)/kCa,pH = 6.3327

and below that pH the equation has no non-negative rest point at all -- it drives cytosolic
Ca2+ negative. This is arithmetic on Ke's own printed constants, not a solver artefact, and it
is not repairable without inventing a floor Ke does not print. :data:`PH_PRODUCTION_FLOOR`
records it and :meth:`CalciumAxis.rates_per_h` RAISES below it rather than clipping.

The consequence is uncomfortable and is the honest headline of this build: ``Control.ph`` is
gated [3.0, 8.0] at ``engine.py:724``, and every protocol in this repository runs at pH 4.0-5.0.
**The calcium axis is undefined at the repository's own operating pH.** It becomes defined only
when the panel's ``sodium_hydroxide`` stressor ("pH units above 5", full dose 3.0) is dosed above
+1.33 units, and enters Ke's own validated window [6.5, 8.0] only above +1.5 units. Across that
window the axis is worth having: Crz1p swings 29.6x from pH 6.5 to 8.0. Below it, the correct
behaviour is to refuse, and a wiring that silently clamps pH to 6.5 to keep the axis alive would
be inventing the one number Ke does not supply.

WHAT IS NOT DRIVEN BY THE ENGINE, said plainly rather than left to read as driven. Eqn 3.3 has
four Ca2+ sources and this engine can move exactly one of them.

* ``kCa,pH*(pH_ext - 6.5)`` IS driven, by ``Control.ph``. This is the whole axis.
* the two sodium arms are NOT. ``kCa,cyt`` (Hill h = 12 in cytosolic Na+) and ``kCa,ext``
  (h = 2 in external Na+) have no engine signal behind them: ``mech/contracts.py``'s
  EXTRACELLULAR is exactly eight species and none is sodium, and ``osmolyte`` is deliberately
  ion-agnostic, so it cannot stand in for Na+ without asserting an ionic fraction nobody
  measured. Both arms therefore run at KE'S OWN unstressed baseline -- Table S5's
  ``Init_Na`` = 100 mM inside, the article's 5 mM outside -- and are CONSTANTS in every engine
  run. :data:`NOT_DRIVEN` refuses the two engine channels by name. At that baseline the two
  arms contribute 6.9e-09 and 2.4e-09 mM/s against CCa's 2.5e-06, i.e. 0.37% of the
  unstressed Ca2+ source between them, so holding them costs little at rest -- but it also
  means this axis CANNOT answer the panel's NaCl -> calcium weight of 0.25.
* ``CCa`` is Ke's constant basal influx and has no input at all.

WHAT MAY NEVER BE CLAIMED FROM THIS BUILD.

1. **It does not answer ``calcium_chloride`` (panel weight 1.00, the axis's heaviest
   stressor).** Eqn 3.3 contains no extracellular-calcium term of any kind -- Ke's calcium is
   a closed cytosolic pool with a constant source. A CaCl2 dose cannot enter this model even in
   principle, and the only dynamic yeast model that does take extracellular Ca2+ as an input
   (Cui & Kaandorp 2006) needs a medium calcium species ``EXTRACELLULAR`` does not have. The
   honest panel channel for this axis is ``sodium_hydroxide -> calcium`` at +0.35, and nothing
   else. :data:`PANEL_CHANNEL` states that, with the two weights it forfeits.
2. **Crz1p here has no transcriptional output.** The panel declares this module as
   ``Calcineurin -> Crz1`` reading a CDRE, and Ke's own Crz1p -> ENA1 arm (Eqn 3.6) is already
   recorded UNTRANSCRIBED in ``ph_ke2013``: Table S3 prints ``KmENA1,Nrg1`` as ``1.143.0*1e-4``,
   which is not a number, and the arm is downstream of the untranscribable Nrg1p ODE anyway.
   crz1_mM is a regulator concentration and a CDRE reporter is NOT derivable from it here.
3. **Eqn 3.2's stress-onset Ca2+ spike is not applied.** Ke imposes a discontinuous Ca2+ reset
   at the instant stress begins. A continuously-controlled fed-batch has no such instant, so the
   trigger time is undefined; picking one would be inventing a protocol event.
4. **The alkaline pH module of the panel is a different axis and stays refused.** Rim101 has no
   kinetic model in any organism (``mech/ph.py:650``), and ``sodium_hydroxide -> alkaline_ph``
   at 1.00 is untouched by this build.

WHAT THIS AXIS DOES GIVE BACK, so it is not a dangling appendage. Activated calcineurin enters
the ion system through ``KmTrk,Cn`` (Table S2, 5e-4 mM): Eqn 2.24 splits the Trk system between
its medium- and high-affinity states as ``p_high = (CN/KmTrk,Cn + 1)/(Ppz/KmTrk,Ppz + CN/KmTrk,Cn
+ 1)``, so calcineurin pushes Trk toward high affinity and changes K+/Na+ uptake.
:meth:`CalciumAxis.trk_affinity_split` exposes exactly that term. Ppz itself (Eqn 2.21) IS
engine-drivable, from cytosolic pH, which the engine already computes.

Ke's OWN published test of that back-coupling is reproduced here and is reported rather than
scored, because the paper states it qualitatively: "The intracellular Na+ did not show notable
difference in cells treated with and without FK506 (Fig. 7A), suggesting that calcineurin
activation is not required for K+/Na+ homeostasis for alkaline pH adaptation." Under this
module's coupled fixed point at pH 8.0, blocking calcineurin moves intracellular Na+ by +10.0%
and K+ by -5.7%, against a 4.7-fold rise in Ca2+ and a 48-fold rise in activated calcineurin
over the same step. The direction and the smallness carry; "notable" is not a number in the
source, so this is REPORTED and is not claimed as a quantitative reproduction.

UNITS AND TIME BASE. Every state here is a concentration in mM, which is Ke's basis for the
signalling module (Text S1 3.3-3.4 and Tables S3/S5) -- unlike the three cations, which
``ph_ke2013`` carries as attomole AMOUNTS. Ke's time base is seconds; the engine's is hours, so
:meth:`CalciumAxis.rates_per_h` multiplies by :data:`SECONDS_PER_HOUR`, the same factor
``signalling.NativeHogResponse`` already applies to its SBML source.

THE FREE-SCALAR ACCOUNTING IS ``ph_ke2013``'s, NOT A SECOND ONE. Every constant this axis reads
is already registered in ``KE2013_PARAMS`` and already counted there. Registering them again
would double-count the same rows in the gate, so :data:`USED_PARAMETERS` holds references to the
existing Param objects and the only new registry here, :data:`NOT_DRIVEN`, contains refusals --
which carry no value and cannot inflate a free-scalar count.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

from .contracts import Compartment, ScientificRefusal, finite
from .params import Param, ParamRegistry
from .ph_ke2013 import (
    ARTIFACTS as KE2013_ARTIFACTS,
    AVAILABILITY as KE2013_AVAILABILITY,
    C_CA,
    C_CRZ1,
    D_CA,
    D_CRZ1,
    Environment,
    H_NA_CYT,
    H_NA_EXT,
    INIT_CA_MM,
    INIT_CN_MM,
    INIT_CNOFF_MM,
    INIT_CRZ1_MM,
    INIT_NA_MM,
    IonState,
    KM_CA_CYT,
    KM_CA_EXT,
    KM_PPZ,
    KM_TRK_CN,
    KM_TRK_PPZ,
    K_CA_CYT,
    K_CA_EXT,
    K_CA_PH,
    K_CN_A,
    K_CN_DA,
    K_CN_PPZ_DA,
    K_CRZ1,
    PPZ_TOTAL,
    SOURCE as KE2013_SOURCE,
    calcineurin_rhs,
    calcineurin_steady_state,
    fluxes,
    steady_state,
)

__all__ = [
    "ARTIFACTS",
    "AVAILABILITY",
    "CalciumAxis",
    "CalciumState",
    "KE_VALIDATED_PH",
    "NOT_APPLIED",
    "NOT_DRIVEN",
    "PANEL_CHANNEL",
    "PH_PRODUCTION_FLOOR",
    "PUBLISHED",
    "SECONDS_PER_HOUR",
    "SOURCE",
    "STATE_VARIABLES",
    "USED_PARAMETERS",
    "alkaline_response",
    "coupled_state",
    "fk506_comparison",
    "ppz_mM",
    "production_floor_ph",
    "reproduction_report",
    "table_s5_cross_check",
    "trk_affinity_split",
]

AVAILABILITY = KE2013_AVAILABILITY
"""Inherited verbatim from ``ph_ke2013``: ``transcribed_here``.

This module adds no availability of its own because it adds no transcription. The equations are
Ke's, the executable is ``ph_ke2013``'s, and this file is the wiring. Wiring an existing
transcription cannot promote it, and nothing here may ever be recorded above the level of the
block it imports.
"""

SOURCE = KE2013_SOURCE
ARTIFACTS = KE2013_ARTIFACTS
"""The same eight pinned upstream files ``ph_ke2013`` hashes. No new artifact is vendored."""

SECONDS_PER_HOUR = 3600.0
"""Ke prints s^-1 throughout Tables S2-S3; the engine's state vector is per hour."""

KE_VALIDATED_PH = (6.5, 8.0)
"""The only external pH values Ke simulates: unstressed 6.5 and the Fig. 7 alkaline stress 8.0.

Ke's Results state the unstressed condition is "external Na+ and K+ concentrations and the
external pH are assumed to be 5 mM, 1 mM, and pH 6.5", and Fig. 7 is the pH 8.0 step. There is
no acidic simulation anywhere in the article, which is why the sign problem below pH 6.5 never
surfaces in the source.
"""

PUBLISHED = {
    "table_s5_init_ca_mM": 5e-5,
    "table_s5_init_cnoff_mM": 1.1628e-3,
    "table_s5_init_cnon_mM": 0.0,
    "table_s5_init_crz1_mM": 1.916e-4,
    "unstressed_ph_ext": 6.5,
    "alkaline_ph_ext": 8.0,
}
"""Every source number this wiring is scored against, with where it is printed.

The four ``table_s5_*`` rows are Table S5 (pinned ``s013.pdf``), whose footnote states they are
the unstressed steady state -- which is what makes them an independent target for constants read
out of Table S3. The two pH values are the Results section's own unstressed and stressed
conditions.
"""

USED_PARAMETERS = MappingProxyType({
    param.name: param for param in (
        # Table S3, Eqn 3.3: the cytosolic Ca2+ balance.
        C_CA, D_CA, K_CA_CYT, KM_CA_CYT, H_NA_CYT, K_CA_EXT, KM_CA_EXT, H_NA_EXT, K_CA_PH,
        # Table S3, Eqn 3.4: calcineurin activation and Crz1p.
        K_CN_A, K_CN_DA, K_CN_PPZ_DA, C_CRZ1, D_CRZ1, K_CRZ1,
        # Table S2, Eqns 2.21 and 2.24: the Ppz sensor and the Trk back-coupling.
        KM_PPZ, PPZ_TOTAL, KM_TRK_PPZ, KM_TRK_CN,
        # Table S5: the four initial conditions, plus the held cytosolic Na+ baseline.
        INIT_CA_MM, INIT_CNOFF_MM, INIT_CN_MM, INIT_CRZ1_MM, INIT_NA_MM,
    )
})
"""The twenty-four Ke rows this axis reads, by reference to ``KE2013_PARAMS``.

Every one is already registered and already counted there, and every one is ASSERTED -- Ke's own
"Taken from" / "Estimated from" provenance column travels in each ``source`` string, and nine of
these rows are asterisked in the source as "adjusted during the model integration step", i.e.
fitted by Ke rather than measured. Reading them through this mapping keeps that visible.
"""

NOT_DRIVEN = ParamRegistry("mech/calcium_ke2013.py::no-engine-channel")
"""The engine inputs Eqn 3.3 asks for and ``mech/contracts.py`` cannot supply."""

NOT_DRIVEN.add(Param.refused(
    "extracellular_calcium_mM", "mM",
    "Ke 2013 Text S1 Eqn 3.3, read in full: the four Ca2+ source terms are CCa (a constant), "
    "kCa,cyt in cytosolic Na+, kCa,ext in external Na+, and kCa,pH in external pH. There is no "
    "extracellular-calcium term",
    reason="The source structurally cannot take a calcium dose: its cytosolic Ca2+ pool has a "
           "constant basal influx and no medium calcium at all, so CaCl2 has nowhere to enter. "
           "Independently, mech/contracts.py's EXTRACELLULAR is exactly eight species (glucose, "
           "ethanol, nitrogen, oxygen, glycerol, peroxide, acetate, osmolyte) and none is "
           "calcium, so even a source that took the dose could not be given it. This is why the "
           "panel's calcium_chloride stressor (weight 1.00) is NOT answered by this build",
    missing="both halves: a medium calcium species in the EXTRACELLULAR contract, AND a "
            "published parameterised S. cerevisiae model whose Ca2+ influx depends on it. Cui & "
            "Kaandorp 2006 is the only candidate for the second and needs the first"))

NOT_DRIVEN.add(Param.refused(
    "cytosolic_sodium_mM", "mM",
    "Ke 2013 Table S3 rows kCa,cyt (8.0e-6 mM/s) and KmCa_cyt (180 mM) with h_Na_cyt = 12; "
    "Eqn 3.3's first saturating source term",
    reason="The engine has no intracellular sodium state. mech/contracts.py's INTRACELLULAR "
           "carries no Na+, and 'osmolyte' is deliberately ion-agnostic -- reading it as sodium "
           "would assert an ionic fraction of the osmolyte pool that nothing in this repository "
           "measures. The arm therefore runs at Ke's OWN Table S5 baseline Init_Na = 100 mM and "
           "is a constant in every engine run, not a signal",
    missing="an intracellular sodium species in the shared state vector, fed by a medium sodium "
            "species and a transport law; ph_ke2013 has the transport laws but the engine has "
            "neither species"))

NOT_DRIVEN.add(Param.refused(
    "external_sodium_mM", "mM",
    "Ke 2013 Table S3 rows kCa,ext (2.4e-5 mM/s) and KmCa_ext (500 mM) with h_Na_ext = 2; "
    "Eqn 3.3's second saturating source term",
    reason="Same absence as the cytosolic arm, one compartment out: EXTRACELLULAR has no sodium "
           "species, so a NaCl dose cannot be keyed into Control.feed_mM. The arm runs at the "
           "article's unstressed external 5 mM and is a constant. This is why the panel's NaCl "
           "-> calcium weight of 0.25 is NOT answered by this build either",
    missing="a medium sodium species in the EXTRACELLULAR contract; the same contracts decision "
            "that blocks iron, copper and the CaCl2 arm"))

NOT_DRIVEN.add(Param.refused(
    "crz1_to_cdre_gain", "reporter units per mM Crz1p",
    "Ke 2013 Text S1 Eqn 3.6 and Table S3's ENA1 block, both already recorded UNTRANSCRIBED in "
    "mech/ph_ke2013.py's ODE_SYSTEM",
    reason="Crz1p is where Ke's calcium block ends as an executable object. Its only printed "
           "output arm is ENA1 transcription, which is untranscribable twice over: it is "
           "downstream of the Nrg1p ODE (which needs KmNrg1,pH, a Michaelis constant that "
           "appears in Eqn 3.5 and in no table) and Table S3 prints its own KmENA1,Nrg1 as "
           "'1.143.0*1e-4', which is not a number. The panel declares this module against a "
           "CDRE reporter and no CDRE gain exists in this source at all",
    missing="a published parameterised CDRE (PMC1/CMK2) transcription model, or a measured "
            "Crz1p-to-CDRE dose response in the working strain. Neither is what Ke 2013 is"))


NOT_APPLIED = MappingProxyType({
    "eqn_3_2_stress_onset_calcium_spike": (
        "Ke imposes a discontinuous cytosolic Ca2+ reset at the instant a stress begins. It is "
        "transcribed nowhere in this repository and is deliberately not applied: a "
        "continuously-controlled fed-batch has no stress-onset instant, so the trigger time is "
        "undefined, and choosing one would invent a protocol event rather than read a constant. "
        "The consequence is that this axis reproduces Ke's ADAPTED levels and not his transient "
        "peaks, and no transient calcium amplitude may be claimed from it."),
    "table_s5_init_cnon_zero": (
        "Table S5's Init_CNon = 0 is carried forward UNCHANGED even though it is not the rest "
        "point of Eqn 3.4. See the module docstring: the block relaxes to 1.20e-05 mM activated "
        "calcineurin and Crz1p to 2.54x its Table S5 row. This is the source's inconsistency and "
        "is regression-tested in both this module's test and "
        "tests/test_ke2013_transcription.py:206. It must not be tuned away."),
})
"""Two pieces of Ke's own model that run nowhere here, each with the reason it does not."""


PANEL_CHANNEL = MappingProxyType({
    "module": "calcium",
    "answered": MappingProxyType({"sodium_hydroxide": 0.35}),
    "forfeited": MappingProxyType({"calcium_chloride": 1.00, "NaCl": 0.25}),
    "note": (
        "generator/stress_panel.py routes three stressors to the calcium module. This build "
        "answers exactly one of them, and only across part of its dose range: sodium_hydroxide "
        "is dosed in 'pH units above 5' with a full dose of 3.0, so it reaches Control.ph 8.0, "
        "but doses below +1.33 units leave Control.ph under PH_PRODUCTION_FLOOR where Eqn 3.3 "
        "has no non-negative rest point, and doses below +1.5 units are outside Ke's own "
        "validated window. calcium_chloride (1.00) and NaCl (0.25) are forfeited for the "
        "reasons NOT_DRIVEN gives: no medium calcium and no sodium anywhere in the contract. "
        "The module's declared CDRE output is forfeited too -- see NOT_DRIVEN's "
        "crz1_to_cdre_gain. This build supplies four regulator states, not panel coverage."),
})
"""The honest mapping from this axis to `generator/stress_panel.py`, including what it loses."""


STATE_VARIABLES = (
    ("ca_mM", "mM", Compartment.CELL_WATER, "cytosolic free Ca2+, Eqn 3.3"),
    ("cn_off_mM", "mM", Compartment.REGULATORY, "inactive calcineurin, Eqn 3.4"),
    ("cn_mM", "mM", Compartment.REGULATORY, "activated calcineurin, Eqn 3.4"),
    ("crz1_mM", "mM", Compartment.REGULATORY, "Crz1p, Eqn 3.4"),
)
"""The four states this axis adds, in `mech/contracts.py`'s vocabulary.

All four are CONCENTRATIONS, which is Ke's basis for the signalling module. ``ph_ke2013``'s
three cations are AMOUNTS in the same source; anything wiring both must not mix them.
"""


@dataclass(frozen=True)
class CalciumState:
    """The four states of Eqns 3.3-3.4, in mM.

    Non-negativity is enforced at construction rather than clipped inside a right-hand side: a
    negative calcineurin pool is a solver telling you the axis was integrated outside its domain,
    and silencing it would hide exactly the pH-floor failure this module exists to declare.
    """

    ca_mM: float
    cn_off_mM: float
    cn_mM: float
    crz1_mM: float

    def __post_init__(self) -> None:
        for field in ("ca_mM", "cn_off_mM", "cn_mM", "crz1_mM"):
            object.__setattr__(self, field,
                               finite(getattr(self, field), field, minimum=0.0))

    @classmethod
    def from_table_s5(cls) -> CalciumState:
        """Ke's Table S5 rows, which its own footnote calls the unstressed steady state.

        For this block they are not; see :data:`NOT_APPLIED`. They are still the published
        starting point and are used unchanged.
        """
        return cls(float(INIT_CA_MM), float(INIT_CNOFF_MM), float(INIT_CN_MM),
                   float(INIT_CRZ1_MM))

    @property
    def calcineurin_total_mM(self) -> float:
        """Ke's conserved calcineurin moiety: Eqn 3.4 moves protein between the two pools only."""
        return self.cn_off_mM + self.cn_mM

    @property
    def activated_fraction(self) -> float:
        total = self.calcineurin_total_mM
        if total <= 0.0:
            raise ScientificRefusal("an activated fraction of an empty calcineurin pool is not "
                                    "defined; Ke's Table S5 total is 1.1628e-3 mM")
        return self.cn_mM / total

    def as_tuple(self) -> tuple[float, float, float, float]:
        return (self.ca_mM, self.cn_off_mM, self.cn_mM, self.crz1_mM)


def ppz_mM(cytosolic_ph: float) -> float:
    """Eqn 2.21: Ppz = KmPpz/(KmPpz + [H+]cyt) * Ppz0, in mM.

    Retyped from the copy inline in :func:`~mech.ph_ke2013.fluxes`, which does not export it;
    the test asserts this against that implementation rather than trusting the retype.
    """
    ph = finite(cytosolic_ph, "cytosolic_ph", minimum=0.0, maximum=14.0)
    h_cyt_mM = 10.0 ** (-ph) * 1e3
    return float(KM_PPZ) / (float(KM_PPZ) + h_cyt_mM) * float(PPZ_TOTAL)


def trk_affinity_split(cn_mM: float, ppz_mM_value: float) -> tuple[float, float]:
    """Eqn 2.24's Trk medium/high affinity split, the calcium axis's back-coupling.

    Returns ``(p_medium, p_high)``. Ppz pushes Trk to medium affinity and calcineurin pushes it
    back to high affinity, so this is the term through which Crz1's upstream reaches K+/Na+ uptake.
    """
    cn = finite(cn_mM, "cn_mM", minimum=0.0)
    ppz = finite(ppz_mM_value, "ppz_mM", minimum=0.0)
    denominator = ppz / float(KM_TRK_PPZ) + cn / float(KM_TRK_CN) + 1.0
    return ((ppz / float(KM_TRK_PPZ)) / denominator,
            (cn / float(KM_TRK_CN) + 1.0) / denominator)


def production_floor_ph(*, na_cyt_mM: float = float(INIT_NA_MM),
                        na_ext_mM: float = Environment().na_ext_mM) -> float:
    """The external pH below which Eqn 3.3 has no non-negative rest point.

    At Ca2+ = 0 every removal term in Eqn 3.3 vanishes, so total production is
    ``CCa + cyt_Na_term + ext_Na_term + kCa,pH*(pH_ext - 6.5)``. Where that is negative the
    equation drives cytosolic Ca2+ below zero, whatever the solver does. This is arithmetic on
    Ke's printed constants; the source never simulates below pH 6.5 and so never meets it.
    """
    na_c = finite(na_cyt_mM, "na_cyt_mM", minimum=0.0)
    na_e = finite(na_ext_mM, "na_ext_mM", minimum=0.0)
    h_c, h_e = float(H_NA_CYT), float(H_NA_EXT)
    cyt = float(K_CA_CYT) * na_c ** h_c / (na_c ** h_c + float(KM_CA_CYT) ** h_c)
    ext = float(K_CA_EXT) * na_e ** h_e / (na_e ** h_e + float(KM_CA_EXT) ** h_e)
    return 6.5 - (float(C_CA) + cyt + ext) / float(K_CA_PH)


PH_PRODUCTION_FLOOR = production_floor_ph()
"""6.3327 at Ke's own unstressed sodium. Below it this axis refuses; it does not clip.

``Control.ph`` is gated [3.0, 8.0] and every protocol in this repository runs at 4.0-5.0, so the
axis is undefined at this repository's operating pH. That is the finding, not a bug to route
around.
"""


class CalciumAxis:
    """Ke 2013's Eqns 3.3-3.4 as a stress axis: four states, one engine driver, one back-coupling.

    The right-hand side is :func:`~mech.ph_ke2013.calcineurin_rhs` unchanged -- this class adds
    the unit conversion, the domain guard, and the two held sodium baselines the engine cannot
    supply. Nothing it does is fitted, and it introduces no constant of its own.

    Args:
        na_cyt_mM: The HELD cytosolic Na+ Eqn 3.3's h = 12 arm reads. Defaults to Ke's Table S5
            ``Init_Na``. It is an argument only so the test can measure what holding it costs;
            an engine caller has nothing to pass, which is what :data:`NOT_DRIVEN` records.
        na_ext_mM: The HELD external Na+ the h = 2 arm reads, defaulting to the article's
            unstressed 5 mM. Same status.
        allow_outside_validated_ph: Ke simulates only pH 6.5 and 8.0. Between
            :data:`PH_PRODUCTION_FLOOR` and 6.5 the equations are defined but unvalidated;
            this flag decides whether that band returns numbers or raises. It never opens the
            band below the floor, which is not a matter of validation but of sign.
    """

    def __init__(self, *, na_cyt_mM: float = float(INIT_NA_MM),
                 na_ext_mM: float = Environment().na_ext_mM,
                 allow_outside_validated_ph: bool = True) -> None:
        self.na_cyt_mM = finite(na_cyt_mM, "na_cyt_mM", minimum=0.0)
        self.na_ext_mM = finite(na_ext_mM, "na_ext_mM", minimum=0.0)
        self.allow_outside_validated_ph = bool(allow_outside_validated_ph)
        self.ph_floor = production_floor_ph(na_cyt_mM=self.na_cyt_mM, na_ext_mM=self.na_ext_mM)
        self.provenance = MappingProxyType(dict(USED_PARAMETERS))
        self.metadata = MappingProxyType({
            "source": SOURCE,
            "availability": AVAILABILITY,
            "executable_block": "mech/ph_ke2013.py::calcineurin_rhs (imported, not retyped)",
            "equations": ("Text S1 Eqn 3.3 (cytosolic Ca2+), Eqn 3.4 (calcineurin pools and "
                          "Crz1p), Eqn 2.21 (Ppz), Eqn 2.24 (Trk affinity split)"),
            "engine_driver": "Control.ph, through Eqn 3.3's kCa,pH*(pH_ext - 6.5) term only",
            "held_constant": ("cytosolic Na+ at Table S5's 100 mM and external Na+ at the "
                              "article's 5 mM; neither has an engine channel"),
            "time_base": "source seconds; rates_per_h multiplies by 3600",
            "validated_ph": KE_VALIDATED_PH,
            "production_floor_ph": self.ph_floor,
        })

    def initial_state(self) -> CalciumState:
        return CalciumState.from_table_s5()

    def check_ph(self, medium_ph: float) -> float:
        """Refuse a medium pH this axis is not defined at, before any integration begins."""
        ph = finite(medium_ph, "medium_ph", minimum=0.0, maximum=14.0)
        if ph < self.ph_floor:
            raise ScientificRefusal(
                f"external pH {ph:.4f} is below {self.ph_floor:.4f}, where Ke 2013's Eqn 3.3 has "
                "no non-negative rest point: its kCa,pH*(pH_ext - 6.5) term is linear and "
                "unbounded below, and at Ca2+ = 0 total production is already negative. Ke never "
                "simulates an acidic condition, so the source supplies no floor and clamping to "
                "one would invent the constant it does not print. Every protocol in this "
                "repository runs at pH 4-5, so the calcium axis is undefined there")
        low, high = KE_VALIDATED_PH
        if not self.allow_outside_validated_ph and not low <= ph <= high:
            raise ScientificRefusal(
                f"external pH {ph:.4f} is outside Ke 2013's simulated window {KE_VALIDATED_PH}; "
                "the article reports only the unstressed pH 6.5 and the Fig. 7 pH 8.0 step")
        return ph

    def rates_per_h(self, state: CalciumState, *, medium_ph: float,
                    cytosolic_ph: float | None = None,
                    ppz_mM_value: float | None = None) -> tuple[float, float, float, float]:
        """Eqns 3.3-3.4 in mM/h, for ``(ca_mM, cn_off_mM, cn_mM, crz1_mM)``.

        Exactly one of ``cytosolic_ph`` (the engine-natural input, which Eqn 2.21 turns into Ppz)
        or ``ppz_mM_value`` (Ppz directly, for reproducing Ke against his own ion state) must be
        given. ``Vratio`` is 0 throughout: Ke's dilution term needs the volume module, which
        ``ph_ke2013`` records as not resting at Table S5's own initial volume.
        """
        if not isinstance(state, CalciumState):
            raise ValueError("state must be a CalciumState")
        if (cytosolic_ph is None) == (ppz_mM_value is None):
            raise ValueError("give exactly one of cytosolic_ph or ppz_mM_value; Ppz is Eqn 2.21 "
                             "of cytosolic pH and passing both invites two disagreeing values")
        ph_ext = self.check_ph(medium_ph)
        ppz = ppz_mM(cytosolic_ph) if ppz_mM_value is None else finite(
            ppz_mM_value, "ppz_mM_value", minimum=0.0)
        env = Environment(ph_ext=ph_ext, na_ext_mM=self.na_ext_mM, k_ext_mM=1.0)
        per_s = calcineurin_rhs(*state.as_tuple(), ppz_mM=ppz, na_cyt_mM=self.na_cyt_mM,
                                env=env, v_ratio_per_s=0.0)
        return tuple(rate * SECONDS_PER_HOUR for rate in per_s)

    def trk_affinity_split(self, state: CalciumState, *, cytosolic_ph: float | None = None,
                           ppz_mM_value: float | None = None) -> tuple[float, float]:
        """The back-coupling, as ``(p_medium, p_high)``. See :func:`trk_affinity_split`."""
        if (cytosolic_ph is None) == (ppz_mM_value is None):
            raise ValueError("give exactly one of cytosolic_ph or ppz_mM_value")
        ppz = ppz_mM(cytosolic_ph) if ppz_mM_value is None else finite(
            ppz_mM_value, "ppz_mM_value", minimum=0.0)
        return trk_affinity_split(state.cn_mM, ppz)


def table_s5_cross_check() -> dict:
    """Table S3's rate constants against Table S5's initial concentrations, two pinned PDFs.

    Table S5's footnote states its rows ARE the unstressed steady state, so each rest point that
    Table S3 alone implies is an independent target. Two of the three land on it and one does
    not; nothing here is fitted and none of the three may be adjusted to improve the agreement.
    """
    ca_from_rates = float(C_CA) / float(D_CA)
    crz1_from_rates = float(C_CRZ1) / float(D_CRZ1)
    env = Environment()
    ion = IonState.from_table_s5()
    ppz = fluxes(ion, env)["Ppz"]
    published = CalciumState.from_table_s5()
    per_s = calcineurin_rhs(*published.as_tuple(), ppz_mM=ppz, na_cyt_mM=ion.na_mM, env=env)
    return {
        "ca": {"from_table_s3": ca_from_rates, "table_s5": float(INIT_CA_MM),
               "relative_difference": ca_from_rates / float(INIT_CA_MM) - 1.0,
               "arithmetic": "CCa/dCa, the rest point of Eqn 3.3 with every other term zero"},
        "crz1": {"from_table_s3": crz1_from_rates, "table_s5": float(INIT_CRZ1_MM),
                 "relative_difference": crz1_from_rates / float(INIT_CRZ1_MM) - 1.0,
                 "arithmetic": "CCrz1/dCrz1, the rest point of Eqn 3.4 at CN = 0 -- which is "
                               "Table S5's own Init_CNon, so this is the rest point Ke intended"},
        "calcineurin": {"d_cn_off_per_s_at_table_s5": per_s[1],
                        "table_s5_cn_mM": float(INIT_CN_MM),
                        "reproduces": bool(abs(per_s[1]) < 1e-15),
                        "arithmetic": "Eqn 3.4 evaluated at Table S5's own rows; a rest point "
                                      "would give zero and this does not"},
        "ppz_mM_at_table_s5": ppz,
        "rhs_at_table_s5_per_s": per_s,
    }


def coupled_state(medium_ph: float, *, fk506: bool = False, tolerance: float = 1e-12,
                  max_iterations: int = 60) -> dict:
    """The calcium block and ``ph_ke2013``'s ion equations solved together at one external pH.

    The two systems are mutually coupled: activated calcineurin sets the Trk affinity split
    (Eqn 2.24) which sets cytosolic Na+, and cytosolic Na+ feeds Eqn 3.3's h = 12 arm while
    cytosolic H+ sets Ppz. This is a fixed-point iteration between them, under exactly the three
    reductions ``ph_ke2013`` declares -- fixed volume, held Hog1p/glycerol, Ena1 ratio 1.

    ``fk506`` reproduces Ke's Fig. 7 inhibitor arm by blocking calcineurin ACTIVITY: the ion
    system sees CN = 0 while Eqn 3.4 still reports what would have been activated. That is the
    reading Ke's Table 1 gives -- "Addition of FK506 blocks the activity of calcineurin".
    """
    ph = CalciumAxis().check_ph(medium_ph)
    env = Environment(ph_ext=ph)
    coupled_cn = 0.0
    ion = None
    block = None
    for _ in range(int(max_iterations)):
        ion = steady_state(env, cn_mM=coupled_cn, guess=ion)
        block = calcineurin_steady_state(env, ion)
        target = 0.0 if fk506 else block["cn_mM"]
        if abs(target - coupled_cn) <= tolerance * max(abs(target), 1e-18):
            coupled_cn = target
            break
        coupled_cn = target
    else:
        raise ScientificRefusal(
            f"the ion/calcineurin fixed point did not converge at external pH {ph}; a "
            "non-converged pair is not a steady state and is not reported as one")
    ppz = block["ppz_mM"]
    # Eqn 3.4 conserves the calcineurin moiety, so the inactive pool is the Table S5 total less
    # the activated one; calcineurin_steady_state solves under that same conservation.
    total_cn = float(INIT_CNOFF_MM) + float(INIT_CN_MM)
    state = CalciumState(block["ca_mM"], total_cn - block["cn_mM"], block["cn_mM"],
                         block["crz1_mM"])
    return {"ph_ext": ph, "fk506": bool(fk506), "state": state, "ion": ion, "ppz_mM": ppz,
            "coupling_cn_mM": coupled_cn, "na_cyt_mM": ion.na_mM, "k_cyt_mM": ion.k_mM,
            "ph_i": ion.ph_i,
            "trk_split": trk_affinity_split(coupled_cn, ppz),
            "crz1_fold_vs_table_s5": block["crz1_fold_vs_table_s5"]}


def alkaline_response(*, ph_from: float = 6.5, ph_to: float = 8.0) -> dict:
    """Ke's own unstressed and Fig. 7 alkaline conditions, and what the axis does between them.

    The fold changes are this module's OUTPUT, not a published number: Ke plots Fig. S7's enzyme
    concentrations but prints no calcineurin or Crz1p value anywhere in the article or the seven
    supplements. They are reported so a drift in the wiring fails a test, and they must not be
    cited as a reproduction of anything.
    """
    before = coupled_state(ph_from)
    after = coupled_state(ph_to)
    return {
        "ph_from": ph_from, "ph_to": ph_to,
        "ca_fold": after["state"].ca_mM / before["state"].ca_mM,
        "cn_fold": after["state"].cn_mM / before["state"].cn_mM,
        "crz1_fold": after["state"].crz1_mM / before["state"].crz1_mM,
        "activated_fraction_before": before["state"].activated_fraction,
        "activated_fraction_after": after["state"].activated_fraction,
        "before": before, "after": after,
    }


def fk506_comparison(*, ph_ext: float = 8.0) -> dict:
    """Ke's published Fig. 7A claim, recomputed: does blocking calcineurin move Na+ at pH 8.0?

    The article states it does not -- "The intracellular Na+ did not show notable difference in
    cells treated with and without FK506 (Fig. 7A), suggesting that calcineurin activation is not
    required for K+/Na+ homeostasis for alkaline pH adaptation". "Notable" is not a number in the
    source, so this is REPORTED, not scored: what the test pins is that the effect stays small
    against the 48-fold calcineurin activation driving it, and its sign.
    """
    wild_type = coupled_state(ph_ext, fk506=False)
    inhibited = coupled_state(ph_ext, fk506=True)
    return {
        "ph_ext": ph_ext,
        "na_cyt_mM": {"wild_type": wild_type["na_cyt_mM"], "fk506": inhibited["na_cyt_mM"]},
        "k_cyt_mM": {"wild_type": wild_type["k_cyt_mM"], "fk506": inhibited["k_cyt_mM"]},
        "relative_na_change": inhibited["na_cyt_mM"] / wild_type["na_cyt_mM"] - 1.0,
        "relative_k_change": inhibited["k_cyt_mM"] / wild_type["k_cyt_mM"] - 1.0,
        "ph_i": {"wild_type": wild_type["ph_i"], "fk506": inhibited["ph_i"]},
        "trk_high_affinity": {"wild_type": wild_type["trk_split"][1],
                              "fk506": inhibited["trk_split"][1]},
        "published_claim": ("Ke 2013 Results, alkaline section: intracellular Na+ shows no "
                            "notable difference with and without FK506 (Fig. 7A). Qualitative "
                            "in the source; reported here, not scored"),
    }


def reproduction_report() -> dict:
    """Everything this wiring is scored on, computed rather than asserted.

    ``tests/test_calcium_ke2013_transcription.py`` asserts every number below. The report
    deliberately carries the failure -- Table S5's calcineurin row -- alongside the two
    successes, because a report that only listed what worked would not be evidence.
    """
    return {
        "source": SOURCE,
        "availability": AVAILABILITY,
        "cross_check": table_s5_cross_check(),
        "alkaline": alkaline_response(),
        "fk506": fk506_comparison(),
        "ph_production_floor": PH_PRODUCTION_FLOOR,
        "engine_ph_gate": (3.0, 8.0),
        "validated_ph": KE_VALIDATED_PH,
        "refusals": tuple(sorted(param.name for param in NOT_DRIVEN.params)),
        "not_applied": tuple(sorted(NOT_APPLIED)),
        "parameters_used": len(USED_PARAMETERS),
    }
