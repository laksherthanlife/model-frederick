"""A Python transcription of the Zheng 2016 / Krakowiak 2018 Hsf1-Hsp70 feedback loop.

**THIS FILE IS OUR ARTIFACT, NOT AN UPSTREAM ONE.** It is a hand port of a published MATLAB
file into this engine's idiom. It is recorded at ``availability = "local_transcription"`` --
never ``local_verified``, which in this repository means a checksummed publisher file loaded
unchanged by the shared SBML loader. What is checksummed here is the MATLAB the port was made
FROM (:data:`SOURCE_DIGESTS`, verified at construction). The Python is ours and carries our
bugs, which is exactly why it is not allowed to inherit the upstream level.

WHY THE PORT IS LEGITIMATE AND THE LABEL IS THE WHOLE QUESTION. A prior survey of this
programme enumerated BioModels by TAXONOMY:4932, all 676 JWS Online models, GitHub SBML and
Europe PMC, and found nothing for the yeast heat-shock response that this engine's SBML loader
could run. Zheng 2016 comes closest: it ships seven real ``.m`` files, recorded at
``local_source_only`` because MATLAB is not the shared runtime. Porting them is transcription of
running code, not retyping equations out of a PDF, and the result is testable against the
source's own published numbers. The only thing that was ever wrong was calling such a thing a
pinned published model.

THREE THINGS THAT CHANGE WHAT THIS TRANSCRIPTION MAY CLAIM.

1. **ONE SOURCE, NOT TWO.** Krakowiak 2018 Source code 1 is byte-identical to Zheng 2016's
   ``titration_YFP_FB.m`` -- sha256 ``5998094f...c5e38e`` for both files. The two papers are
   one model. They must not be recorded as independent corroboration of the same constants,
   and ``mech/state.py``'s ``hsf1_free`` (tau = 1/k2, k2 = 2.783 /min) is downstream of exactly
   this: Krakowiak reprints k2 unchanged, and reprinting is not replication. What Krakowiak
   genuinely adds is one deleted term (the feedback-severed right-hand side), a five-fold lower
   transcription gain, and two new sweep drivers.

2. **THE GAIN IS PUBLISHED TWICE AT A FACTOR OF FIVE.** beta = 1.778 /min (Zheng Table 1) and
   0.3557 /min (Krakowiak's Methods table, printed in the column beside Zheng's) for the same
   model structure, with no interval and no stated reason; the source code writes the second as
   ``1.7783/5``. :data:`TRANSCRIPTION_GAIN` is therefore SWEPT, ``float()`` on it raises, and
   every entry point here takes ``gain`` as a required keyword. There is no default, because a
   default would silently pick a side of a five-fold published disagreement.

3. **TEMPERATURE STAYS REFUSED.** Both papers enter temperature only as an asserted initial
   unfolded-protein load -- 0.5183, 4.4492, 10.5141, 20.0397 a.u. at 25, 35, 39, 43 C, under the
   comment "From Temp vs UP relationship". Neither derives, fits nor measures that relationship;
   Zheng's Methods say the loads other than 39 C were chosen so as to "generally capture steady
   state YFP reporter outputs", which is choosing the input to fit the output it then predicts.
   :data:`TEMPERATURE_TO_UNFOLDED_PROTEIN` is REFUSED, this module has no temperature argument
   anywhere, and unfolded protein is an INPUT.

WHAT THE REPRODUCTION FOUND, INCLUDING THE PARTS THAT FAILED. Every number below is recomputed
by ``tests/test_heat_transcription.py`` and recorded in ``data/transcriptions/zheng2016.json``,
so a drift in the port fails a test rather than surviving in this paragraph.

* Zheng's three published parameter-screen acceptance criteria reproduce under his own gain:
  Hsf1 is 98.4% bound at 25 C; the complex falls to 4.7% of basal by 5 min (required <=10%);
  it recovers to 90% at 41.7 min and stands at 91.7% at 60 min (required >=90% within 60 min).
* **Krakowiak's gain FAILS Zheng's criterion 3 on the byte-identical model file**: 15.4% at
  60 min against a required 90%, a 5.8-fold shortfall, with no recovery inside the window. The
  five-fold split is therefore not a choice between two equally licensed numbers -- one of them
  contradicts the screen that selected its own co-parameters, and neither paper says so.
* Zheng Figure 1D's "temperature-dependent plateaus" hold under his gain (log2 fold 1.15, 2.08,
  2.89 at 35/39/43 C, every curve down to 0.1% of its peak transcription rate by 120 min) and
  fail at 43 C under the reduced gain, which is still running at **74% of peak rate** when the
  published window closes. Its height, 2.75 against 2.89, barely moves -- so the failure is
  invisible in the plotted curve and visible only in its slope.
* Krakowiak's stated reason for lowering the gain reproduces: the wild type / feedback-severed
  separation moves from 18 min to 84 min at a 25% threshold, a 4.7-fold delay.
* Krakowiak Source code 3's published axis identifies which gain made his figures: the basal
  sweep tops out at 9.24 under the reduced gain, inside the published ceiling of 12, and at
  22.67 under Zheng's, off it by 1.9-fold.
* Krakowiak Figure 6A's maximal-output claim reproduces monotonically in k2. **Its deactivation
  claim does not.** Deactivation time is non-monotone in affinity: three-fold higher affinity is
  SLOWER than wild type (24.8 min against 17.9) and five-fold lower is FASTER (15.6 min),
  reversing the published direction at both modest perturbations, across every shut-off
  threshold from 1% to 10% of peak rate. Only the 50-fold case that Source code 4 actually
  ships behaves as stated. Neither paper ships code for the 3x and 5x panels.

TWO DEFECTS IN THE SOURCE, FOUND BY PORTING IT.

* ``titration_YFP_FB.m`` line 30 subtracts ``k5*HSP_UP`` from FREE unfolded protein, which the
  published ``d[UP]/dt`` does not. The shipped code therefore destroys client at twice the rate
  its own paper's equations do. Both conventions run here (``client_loss=``); the difference is
  under 5e-5 in relative reporter output over 240 min, so this is reported as a divergence and
  not as a changed answer.
* ``kdil = 0`` makes the reporter a pure integrator with no steady state, yet Source code 3
  reads its sweep as "steady-state level of YFP reporter at 25 C". At k2 = 1e4 the basal fold
  change is 5.84, 9.24, 14.06 and 20.34 at horizons of 120, 240, 480 and 960 min. That figure's
  y-axis is set by the simulation horizon, not by the model, so nothing here is ever called a
  steady state; the readout is named for the horizon it was taken at.

THE FREE-SCALAR GATE REFUSES THIS PIECE, AND SHOULD. Fifteen free scalars against zero
independent targets: every figure reproduced above is one the parameters were selected against,
so none of them can score the model. :func:`gate` returns that verdict rather than hiding it. A
faithful port of a fitted model is evidence about the port, and about nothing else.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .. import paths
from .integrate import FALLBACK_METHOD, PINNED_METHOD
from .params import Param, ParamRegistry, Target

__all__ = [
    "AFFINITY_MULTIPLIERS",
    "BASAL_SWEEP_HORIZON_MIN",
    "HEAT_PARAMS",
    "HeatShockTranscription",
    "PUBLISHED_AXIS_MAXIMUM",
    "PUBLISHED_FOLD_AXIS",
    "PUBLISHED_GAINS",
    "PUBLISHED_UNFOLDED_LOADS",
    "SCREEN_UNFOLDED_LOAD_39C",
    "ScreenCriteria",
    "SOURCE_DIGESTS",
    "TEMPERATURE_TO_UNFOLDED_PROTEIN",
    "TRANSCRIPTION_GAIN",
    "TRANSCRIPTION_RECORD_PATH",
    "affinity_timecourse",
    "basal_versus_dissociation",
    "feedback_separation",
    "gain",
    "gate",
    "provenance",
    "reporter_fold_change",
    "screen_criteria",
    "transcription_record",
]

NATIVE_TIME_UNIT = "minute"
"""The source's time base, kept. Every rate here is per minute and no conversion happens: three
published time bases collide elsewhere in this architecture and the 60x has already cost errors."""

TRANSCRIPTION_RECORD_PATH = "data/transcriptions/zheng2016.json"

SOURCE_DIGESTS = {
    "data/native_reference_models/zheng2016/titration_YFP_FB.m":
        "5998094f1af311d73ecbff972f48d754a0c1048cf261df0bbe6f0307aee5c38e",
    "data/native_reference_models/krakowiak2018/elife-31668-code1-v3.m":
        "5998094f1af311d73ecbff972f48d754a0c1048cf261df0bbe6f0307aee5c38e",
    "data/native_reference_models/krakowiak2018/elife-31668-code2-v3.m":
        "8f2d4e977ffdc1315fb02f4a9c0ede6010480db4778d73eacbb64670fc7b0c93",
    "data/native_reference_models/zheng2016/Fig1_IPMS.m":
        "43914d010b1c6e9a76250c9110de4578b4c710a534a29d98b7d5c084d777381b",
    "data/native_reference_models/zheng2016/Fig2_YFP_Reporter.m":
        "6058697b10b9e376803ed12514383c3740a0cae5222da800e431e4032453ae84",
    "data/native_reference_models/krakowiak2018/elife-31668-code3-v3.m":
        "39564aad098ed62d9b5be773393ac5c55355ff558365dd357ead8c0a1a34b421",
    "data/native_reference_models/krakowiak2018/elife-31668-code4-v3.m":
        "9817d2cf18ab4bf0b9bb33cbd0745bacd8b89ecc45d52403b622b58d4a9db6df",
}
"""The upstream MATLAB this file was transcribed from, by content. The first two entries are the
byte-identity finding: Krakowiak Source code 1 IS Zheng's ``titration_YFP_FB.m``."""

BYTE_IDENTICAL_PAIR = (
    "data/native_reference_models/zheng2016/titration_YFP_FB.m",
    "data/native_reference_models/krakowiak2018/elife-31668-code1-v3.m",
)

_ZHENG = ("Zheng et al. 2016, eLife 5:e18638, PMID 27831465, PMC5127643; Table 1, "
          "DOI 10.7554/eLife.18638.021, and the shipped MATLAB, DOI 10.7554/eLife.18638.025")
_KRAKOWIAK = ("Krakowiak et al. 2018, eLife 7:e31668, PMID 29393852, PMC5809143; Materials and "
              "methods model parameter table and Source code 1-4, DOIs 10.7554/eLife.31668.014 "
              "through .017")

HEAT_PARAMS = ParamRegistry("mech/heat_zheng2016.py -- Hsf1-Hsp70 transcription")

# Every rate below is a point from a 7e5-set parameter screen, not a measurement. ASSERTED is
# the honest grade: the screen reports which values appear often, never a per-value interval.
CLIENT_ON_RATE = HEAT_PARAMS.add(Param.asserted(
    name="heat.hsp70_client_on_rate_k1_k3", value=166.8, units="min^-1 a.u.^-1",
    source=f"{_ZHENG}, k1,k3 row. Selected by a parameter screen over 7e5 combinations, in "
           "arbitrary abundance units; k3 is SET EQUAL to k1 by the papers' own assumption that "
           "Hsp70's on-rate is the same for Hsf1 and for unfolded protein.",
    missing="a binding assay in a declared absolute abundance scale"))

HSF1_OFF_RATE = HEAT_PARAMS.add(Param.asserted(
    name="heat.hsp70_hsf1_off_rate_k2", value=2.783, units="min^-1",
    source=f"{_ZHENG}, k2 row. Reprinted unchanged by Krakowiak from the byte-identical model "
           "file, which is not an independent second determination. mech/state.py's hsf1_free "
           "carries tau = 1/k2 = 21.6 s from this number.",
    missing="an independent dissociation-rate measurement for Hsp70-Hsf1"))

CLIENT_OFF_RATE = HEAT_PARAMS.add(Param.asserted(
    name="heat.hsp70_client_off_rate_k4", value=0.0464, units="min^-1",
    source=f"{_ZHENG}, k4 row.",
    missing="a client-specific dissociation measurement"))

CLIENT_REFOLDING_RATE = HEAT_PARAMS.add(Param.asserted(
    name="heat.client_refolding_k5", value=4.642e-7, units="min^-1",
    source=f"{_ZHENG}, k5 row. Table 1 rounds to 4.64e-7; every shipped driver writes 4.642e-7 "
           "and the code value is used here because the code produced the figures.",
    missing="a condition-matched client removal assay"))

HSF1_DNA_KD = HEAT_PARAMS.add(Param.asserted(
    name="heat.hsf1_dna_kd", value=0.0022, units="a.u.",
    source=f"{_ZHENG}, Kd row. An arbitrary-unit response constant, not a molar affinity: the "
           "whole state vector is anchored only to [HSP]o = 1 and a 500:1 Hsp70:Hsf1 ratio.",
    missing="absolute Hsf1 and DNA-binding calibration"))

HILL_N = HEAT_PARAMS.add(Param.asserted(
    name="heat.hill_n", value=3.0, units="dimensionless",
    source=f"{_ZHENG}, n (fixed) row. A fixed modelling convention motivated by Hsf1 "
           "trimerisation (Sorger & Nelson 1989), not a fitted or measured response exponent.",
    missing="a measured cooperativity for HSE-driven transcription"))

REPORTER_DILUTION = HEAT_PARAMS.add(Param.asserted(
    name="heat.reporter_dilution_kdil", value=0.0, units="min^-1",
    source=f"{_ZHENG}, kdil (fixed) row: 'Approximate to 0 because of short experiments'. This "
           "makes the reporter a pure integrator with NO steady state, which invalidates "
           "Krakowiak Source code 3's own comment; see basal_versus_dissociation.",
    missing="growth and reporter-turnover observations in the assay condition"))

HSP70_HSF1_RATIO = HEAT_PARAMS.add(Param.asserted(
    name="heat.hsp70_hsf1_abundance_ratio", value=500.0, units="dimensionless",
    source=f"{_ZHENG}, Table 2 and Methods: [HSP.Hsf1]o = 1/500, sourced to 'the range of "
           "measured values' in Chong 2015 and Kulak 2014 without a value being taken from "
           "either. Krakowiak's table prints the same number as 0.002.",
    missing="a matched absolute Hsp70 and Hsf1 quantification in one sample"))

TRANSCRIPTION_GAIN = HEAT_PARAMS.add(Param.swept(
    name="heat.transcription_gain_beta", units="min^-1",
    source=f"PUBLISHED TWICE AT A FACTOR OF FIVE FOR ONE MODEL STRUCTURE. {_ZHENG} prints "
           f"beta = 1.778 (code: 1.7783). {_KRAKOWIAK} prints 0.3557 in the column beside it "
           "(code: 1.7783/5), with no interval and no stated reason; the prose says only that "
           "the first value 'exaggerated' the feedback. Neither is an identified transcription "
           "rate, so this axis has no midpoint and no default.",
    bounds=(1.7783 / 5.0, 1.7783),
    missing="an absolute transcript and reporter calibration -- molecules of HSE-YFP per minute "
            "per active Hsf1 trimer -- able to decide between the two published gains"))

TEMPERATURE_TO_UNFOLDED_PROTEIN = HEAT_PARAMS.add(Param.refused(
    name="heat.temperature_to_unfolded_protein",
    units="a.u. unfolded protein per degree Celsius",
    source="Zheng Fig2_YFP_Reporter.m and Krakowiak Source code 3 and 4 carry the comment 'From "
           "Temp vs UP relationship'; Zheng's Methods state its form but never print its "
           f"coefficients, and neither paper measures unfolded protein. {_ZHENG}; {_KRAKOWIAK}.",
    reason="Temperature enters both papers ONLY as four asserted initial loads (0.5183, 4.4492, "
           "10.5141, 20.0397 a.u. at 25, 35, 39, 43 C). Zheng's Methods say the exponential was "
           "chosen to pass through the screen's 39 C value with the other temperatures picked "
           "to 'generally capture steady state YFP reporter outputs' -- the input chosen to fit "
           "the output it is then used to predict. There is no temperature-entry law to promote.",
    missing="a measured temperature-to-unfolded-protein relationship in S. cerevisiae, in an "
            "abundance scale shared with the model's Hsp70 and Hsf1 pools"))

ABSOLUTE_ABUNDANCE_SCALE = HEAT_PARAMS.add(Param.refused(
    name="heat.absolute_abundance_scale", units="molecules per cell per a.u.",
    source=f"{_ZHENG}, Tables 1 and 2: every state is in arbitrary units anchored only to "
           "[HSP]o = 1 and the 500:1 ratio above.",
    reason="With no absolute anchor, Kd = 0.0022 a.u. is not a molar affinity and "
           "k1 = 166.8 min^-1 a.u.^-1 is not a physical on-rate. Reading either as a physical "
           "constant is the one mistake this module exists to prevent downstream.",
    missing="absolute Hsp70, Hsf1 and unfolded-protein quantification in one matched sample"))

PUBLISHED_GAINS = {"zheng2016": 1.7783, "krakowiak2018": 1.7783 / 5.0}
"""The two published transcription gains, by the paper that printed each. Code values, not the
rounded table values 1.778 and 0.3557."""

PUBLISHED_UNFOLDED_LOADS = {25: 0.5183, 35: 4.4492, 39: 10.5141, 43: 20.0397}
"""Asserted initial unfolded-protein loads, keyed by the temperature the source ATTACHES to each.
The key is a label, not an input: :data:`TEMPERATURE_TO_UNFOLDED_PROTEIN` refuses the map."""

for _label, _load in PUBLISHED_UNFOLDED_LOADS.items():
    # Each load is a free scalar of this piece, not an observation. Registered one by one so
    # the gate counts all four rather than hiding them inside a dict literal.
    HEAT_PARAMS.add(Param.asserted(
        name=f"heat.unfolded_protein_load_{_label}c", value=_load, units="a.u.",
        source=f"Zheng Fig2_YFP_Reporter.m and Krakowiak Source code 3 and 4, under the comment "
               f"'From Temp vs UP relationship'; {_ZHENG}. Table 2 rounds the 39 C value to "
               "10.51. Only the 39 C load was used by the parameter screen; Zheng's Methods say "
               "the others were chosen to 'generally capture steady state YFP reporter outputs'.",
        missing="a measured unfolded-protein load at this temperature in the model's scale"))
del _label, _load

SCREEN_UNFOLDED_LOAD_39C = 10.0
"""One load, two temperatures, in the source's own code: ``Fig1_IPMS.m`` line 33 sets UPo = 10
under a heading reading "Temp = 35 C", and ``Hsf1_Phosph.m`` line 31 writes the same number as
``UPo_39``. Fig2_YFP_Reporter.m gives 35 C the load 4.4492 and 39 C the load 10.5141, and the
screen criteria are stated for a 25 -> 39 C shift, so Fig1_IPMS.m's comment is the stale one."""

BASAL_SWEEP_HORIZON_MIN = 240.0
PUBLISHED_AXIS_MAXIMUM = 12.0
"""Krakowiak Source code 3's horizon and its plotted y-limit. Both are reported because with
kdil = 0 the second is a consequence of the first rather than of the model."""

AFFINITY_MULTIPLIERS = {"k2_over_3": 1.0 / 3.0, "wild_type": 1.0, "k2_times_5": 5.0,
                        "k2_times_50": 50.0}
"""Multipliers on k2. Only ``k2_times_50`` is shipped as code (Source code 4, Figure 6A); the
3x and 5x panels of Figure 6-figure supplement 1 are described in prose with no script."""

_STATE_NAMES = ("HSP", "Hsf1", "UP", "HSP_Hsf1", "HSP_UP", "YFP")
_CLIENT_LOSS = ("source_code", "published_equations")


def gain(source_id: str) -> float:
    """One of the two published transcription gains, pinned on the record.

    The only sanctioned way past :data:`TRANSCRIPTION_GAIN`'s refusal to collapse. Naming the
    paper rather than typing a float is what keeps the five-fold split visible at every call
    site instead of becoming an anonymous literal.
    """
    if source_id not in PUBLISHED_GAINS:
        raise KeyError(
            f"gain must name the paper that published it, one of {sorted(PUBLISHED_GAINS)}; got "
            f"{source_id!r}. The gain is SWEPT over a five-fold published disagreement and has "
            "no default -- see TRANSCRIPTION_GAIN")
    return TRANSCRIPTION_GAIN.at(PUBLISHED_GAINS[source_id])


def _verify_sources() -> dict[str, str]:
    """Checksum every upstream MATLAB file this port was made from, or refuse to run."""
    root = paths.data_dir().parent
    for relative, expected in SOURCE_DIGESTS.items():
        path = root / relative
        if not path.exists():
            raise FileNotFoundError(
                f"the transcription source {relative} is missing; this module is a port OF that "
                "file and will not run without it")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(
                f"{relative} has sha256 {actual}, expected {expected}. The port was made from "
                "the expected bytes; a changed source means the transcription is stale")
    first, second = BYTE_IDENTICAL_PAIR
    if SOURCE_DIGESTS[first] != SOURCE_DIGESTS[second]:
        raise ValueError("the byte-identity finding is the reason these two papers are one "
                         "source; it cannot be edited away in the digest table")
    return dict(SOURCE_DIGESTS)


def transcription_record() -> dict:
    """The recorded transcription: what it was made from, what it refuses, what it reproduced."""
    path = paths.data_dir() / "transcriptions" / "zheng2016.json"
    record = json.loads(path.read_text())
    if record["availability"] == "local_verified":
        raise ValueError(
            "a transcription is our artifact and can never be availability=local_verified, "
            "which means a checksummed upstream file loaded unchanged by the shared loader")
    if record["availability"] != "local_transcription" or record["is_upstream_artifact"]:
        raise ValueError("the transcription record must declare itself a transcription")
    return record


@dataclass(frozen=True)
class HeatShockTranscription:
    """The six-state Hsf1-Hsp70 loop of ``titration_YFP_FB.m``, in Python.

    States, in the source's order: free Hsp70, free Hsf1, unfolded protein, the Hsp70-Hsf1
    complex, the Hsp70-client complex, and the HSE-YFP reporter. All in arbitrary units on the
    source's own scale, all in minutes.

    Args:
        gain: beta, in min^-1. Required, and there is no default -- see :func:`gain`.
        k2: The Hsp70-Hsf1 off-rate. Defaults to the published point; the affinity figures
            multiply it.
        feedback: False deletes the ``beta*hill`` term from ``d[HSP]/dt`` only, which is exactly
            what Krakowiak Source code 2 does and all that separates it from Source code 1.
        client_loss: ``"source_code"`` ports ``titration_YFP_FB.m`` line 30, which subtracts
            ``k5*[HSP.UP]`` from FREE client as well as from the complex. ``"published_equations"``
            ports Zheng's Methods, which do not. The shipped code loses client at twice the
            published rate; the reporter difference is under 5e-5 relative over 240 min.
    """

    gain: float
    k2: float = float(HSF1_OFF_RATE)
    feedback: bool = True
    client_loss: str = "source_code"

    def __post_init__(self) -> None:
        for name in ("gain", "k2"):
            value = getattr(self, name)
            if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite positive rate, got {value!r}")
        if self.client_loss not in _CLIENT_LOSS:
            raise ValueError(
                f"client_loss must be one of {_CLIENT_LOSS}; got {self.client_loss!r}. The two "
                "are the source's code and the source's paper, and they disagree")
        _verify_sources()

    @property
    def state_names(self) -> tuple[str, ...]:
        return _STATE_NAMES

    def rhs(self, _t: float, y) -> np.ndarray:
        """``titration_YFP_FB.m``, line for line."""
        hsp, hsf1, up, hsp_hsf1, hsp_up, _yfp = y
        k1 = float(CLIENT_ON_RATE)
        k3 = k1
        k4 = float(CLIENT_OFF_RATE)
        k5 = float(CLIENT_REFOLDING_RATE)
        kd = float(HSF1_DNA_KD)
        n = float(HILL_N)
        hill = self.gain * hsf1 ** n / (kd ** n + hsf1 ** n)
        free_client_loss = k5 * hsp_up if self.client_loss == "source_code" else 0.0
        return np.array([
            self.k2 * hsp_hsf1 - k1 * hsp * hsf1 + k4 * hsp_up - k3 * hsp * up + k5 * hsp_up
            + (hill if self.feedback else 0.0),
            self.k2 * hsp_hsf1 - k1 * hsp * hsf1,
            k4 * hsp_up - k3 * hsp * up - free_client_loss,
            k1 * hsp * hsf1 - self.k2 * hsp_hsf1,
            k3 * hsp * up - k4 * hsp_up - k5 * hsp_up,
            hill - float(REPORTER_DILUTION) * y[5],
        ])

    def initial_state(self, unfolded_protein: float, *, reporter: float) -> np.ndarray:
        """Zheng Table 2, with the client load supplied rather than derived from a temperature."""
        load = float(unfolded_protein)
        if isinstance(unfolded_protein, bool) or not np.isfinite(load) or load < 0:
            raise ValueError(
                f"unfolded_protein must be a finite non-negative load in a.u., got "
                f"{unfolded_protein!r}. It is an INPUT here: the temperature that would produce "
                "it is refused, see TEMPERATURE_TO_UNFOLDED_PROTEIN")
        return np.array([1.0, 0.0, load, 1.0 / float(HSP70_HSF1_RATIO), 0.0, float(reporter)])

    def simulate(self, times_min, initial_state, *, rtol: float = 1e-10,
                 atol: float = 1e-14) -> np.ndarray:
        """Integrate on the source's minute grid with this engine's pinned stiff pair.

        Returns an ``(len(times_min), 6)`` array in :data:`_STATE_NAMES` order. Tolerances are
        far tighter than MATLAB ode23s's defaults on purpose; every reported figure was checked
        to five significant digits against rtol 1e-3 and against Radau.
        """
        grid = np.asarray(times_min, dtype=float)
        if grid.ndim != 1 or grid.size < 2 or not np.all(np.diff(grid) > 0):
            raise ValueError("times_min must be a strictly increasing 1-D grid of minutes")
        y0 = np.asarray(initial_state, dtype=float)
        if y0.shape != (len(_STATE_NAMES),):
            raise ValueError(f"initial state has shape {y0.shape}, expected (6,)")
        last = None
        for method in (PINNED_METHOD, FALLBACK_METHOD):
            solution = solve_ivp(self.rhs, (grid[0], grid[-1]), y0, t_eval=grid, method=method,
                                 rtol=rtol, atol=atol)
            if solution.success:
                return solution.y.T
            last = solution
        raise RuntimeError(f"heat transcription failed to integrate: {last.message}")


def _complex_ratio(gain_value: float, load: float, *, minutes: float, points: int,
                   **kwargs) -> tuple[np.ndarray, np.ndarray]:
    """``Fig1_IPMS.m``: the shocked Hsp70-Hsf1 complex over its unshocked control."""
    model = HeatShockTranscription(gain=gain_value, **kwargs)
    times = np.linspace(0.0, minutes, points)
    control = model.simulate(times, model.initial_state(0.0, reporter=0.0))
    shocked = model.simulate(times, model.initial_state(load, reporter=0.0))
    return times, shocked[:, 3] / control[:, 3]


@dataclass(frozen=True)
class ScreenCriteria:
    """Zheng's three published parameter-screen acceptance criteria, scored.

    From the Methods, verbatim in substance: Hsf1 must be bound to Hsp70 in basal conditions at
    25 C; the complex must dissociate to <=10% of basal within five minutes of a 25 -> 39 C
    shift; and it must re-associate to >=90% of its initial level within 60 min.

    These are the criteria the screen used to SELECT k1-k5, beta and Kd. Passing them shows the
    port is faithful. It is not evidence that the model is right, and :func:`gate` says so.
    """

    gain_source: str
    unfolded_protein: float
    bound_fraction_basal: float
    ratio_at_5_min: float
    ratio_at_60_min: float
    recovery_time_min: float | None

    @property
    def criterion_1_bound_at_basal(self) -> bool:
        return self.bound_fraction_basal >= 0.9

    @property
    def criterion_2_dissociates(self) -> bool:
        return self.ratio_at_5_min <= 0.10

    @property
    def criterion_3_reassociates(self) -> bool:
        return self.recovery_time_min is not None and self.recovery_time_min <= 60.0

    @property
    def passes(self) -> bool:
        return (self.criterion_1_bound_at_basal and self.criterion_2_dissociates
                and self.criterion_3_reassociates)

    def summary(self) -> str:
        verdict = "PASSES" if self.passes else "FAILS"
        recovery = ("never inside 60 min" if self.recovery_time_min is None
                    else f"{self.recovery_time_min:.2f} min")
        return (f"Zheng screen criteria under beta from {self.gain_source} "
                f"(UP = {self.unfolded_protein:g} a.u.): {verdict} -- "
                f"basal bound fraction {self.bound_fraction_basal:.5f} (>=0.9), "
                f"ratio at 5 min {self.ratio_at_5_min:.5f} (<=0.10), "
                f"re-association {recovery} (<=60 min), "
                f"ratio at 60 min {self.ratio_at_60_min:.5f}")


def screen_criteria(*, gain_source: str,
                    unfolded_protein: float = PUBLISHED_UNFOLDED_LOADS[39],
                    points: int = 6001) -> ScreenCriteria:
    """Score Zheng's three published acceptance criteria under one of the two published gains.

    The result that matters: they pass under Zheng's own gain and criterion 3 FAILS under
    Krakowiak's, on the byte-identical model file, at 0.154 against a required 0.90.
    """
    gain_value = gain(gain_source)
    times, ratio = _complex_ratio(gain_value, float(unfolded_protein), minutes=60.0, points=points)
    model = HeatShockTranscription(gain=gain_value)
    basal = model.simulate(times, model.initial_state(0.0, reporter=0.0))
    total_hsf1 = basal[:, 1] + basal[:, 3]
    trough = int(np.argmin(ratio))
    above = np.flatnonzero(ratio[trough:] >= 0.9)
    recovery = None
    if above.size:
        index = trough + int(above[0])
        recovery = float(np.interp(0.9, [ratio[index - 1], ratio[index]],
                                   [times[index - 1], times[index]])) if index else float(times[0])
    return ScreenCriteria(
        gain_source=gain_source,
        unfolded_protein=float(unfolded_protein),
        bound_fraction_basal=float(np.min(basal[:, 3] / total_hsf1)),
        ratio_at_5_min=float(np.interp(5.0, times, ratio)),
        ratio_at_60_min=float(ratio[-1]),
        recovery_time_min=recovery,
    )


PUBLISHED_FOLD_AXIS = (-0.5, 4.0)
"""``Fig2_YFP_Reporter.m``'s plotted y-limits for log2 reporter fold change."""


def reporter_fold_change(*, gain_source: str, minutes: float = 120.0,
                         points: int = 12001) -> dict:
    """``Fig2_YFP_Reporter.m``, Zheng Figure 1D: reporter fold change at the four asserted loads.

    Keyed by the temperature label the source attaches to each load, which is a label and not an
    input. Published claim: "induction curves that reached temperature-dependent plateaus", where
    Zheng's own gloss on a plateau is that "no more YFP is being produced in the simulation" --
    so ``terminal_rate_per_min`` is the quantity that decides it, not the curve's height.
    """
    gain_value = gain(gain_source)
    model = HeatShockTranscription(gain=gain_value)
    times = np.linspace(0.0, float(minutes), int(points))
    reporter0 = 3.0
    out = {"time_min": times}
    for label, load in PUBLISHED_UNFOLDED_LOADS.items():
        trace = model.simulate(times, model.initial_state(load, reporter=reporter0))[:, 5]
        fold = np.log2(trace / reporter0)
        out[label] = {
            "reporter": trace,
            "log2_fold": fold,
            "log2_fold_at_end": float(fold[-1]),
            "terminal_rate_per_min": float((trace[-1] - trace[-2]) / (times[-1] - times[-2])),
            "within_published_axis": bool(PUBLISHED_FOLD_AXIS[0] <= fold.min()
                                          and fold.max() <= PUBLISHED_FOLD_AXIS[1]),
        }
    return out


def feedback_separation(*, gain_source: str, load: float = PUBLISHED_UNFOLDED_LOADS[39],
                        minutes: float = 240.0, points: int = 241,
                        thresholds=(0.10, 0.25, 0.50)) -> dict:
    """Krakowiak Figure 1B and Figure 1-figure supplement 1: wild type against severed feedback.

    Source code 2 differs from Source code 1 by one deleted term, so the same transcription runs
    both. Returns the first minute at which the severed run exceeds the wild type by each
    fractional threshold -- the quantity Krakowiak says the original gain got wrong.
    """
    gain_value = gain(gain_source)
    times = np.linspace(0.0, float(minutes), int(points))
    intact = HeatShockTranscription(gain=gain_value)
    severed = HeatShockTranscription(gain=gain_value, feedback=False)
    start = intact.initial_state(float(load), reporter=3.0)
    wild = intact.simulate(times, start)[:, 5]
    cut = severed.simulate(times, start)[:, 5]
    excess = (cut - wild) / np.maximum(wild, 1e-12)
    separation = {}
    for threshold in thresholds:
        reached = np.flatnonzero(excess >= threshold)
        separation[threshold] = float(times[reached[0]]) if reached.size else None
    return {"time_min": times, "wild_type": wild, "feedback_severed": cut,
            "separation_min": separation,
            "severed_over_wild_type_at_end": float(cut[-1] / wild[-1])}


def basal_versus_dissociation(*, gain_source: str,
                              load: float = PUBLISHED_UNFOLDED_LOADS[25],
                              horizon_min: float = BASAL_SWEEP_HORIZON_MIN,
                              points: int = 20) -> dict:
    """Krakowiak Source code 3 / Figure 2-figure supplement 1: basal reporter versus k2.

    The source calls its readout a steady state. It is not one: kdil = 0 makes the reporter an
    integrator, so the returned fold change is named for the horizon it was taken at, and the
    horizon is an argument precisely so the dependence is visible.
    """
    gain_value = gain(gain_source)
    times = np.linspace(0.0, float(horizon_min), int(horizon_min) + 1)
    rates = np.logspace(0.0, 4.0, int(points))
    reporter0 = 3.0
    folds = []
    for rate in rates:
        model = HeatShockTranscription(gain=gain_value, k2=float(rate))
        folds.append(model.simulate(times, model.initial_state(float(load), reporter=reporter0))
                     [-1, 5] / reporter0)
    folds = np.asarray(folds)
    return {"k2_per_min": rates, "fold_change_at_horizon": folds,
            "horizon_min": float(horizon_min),
            "monotone_increasing": bool(np.all(np.diff(folds) > 0)),
            "within_published_axis": bool(folds.max() <= PUBLISHED_AXIS_MAXIMUM)}


def affinity_timecourse(*, gain_source: str, multiplier: float,
                        basal_load: float = PUBLISHED_UNFOLDED_LOADS[25],
                        shock_load: float = PUBLISHED_UNFOLDED_LOADS[39],
                        minutes: float = 240.0, points: int = 24001,
                        shutoff_fraction: float = 0.05) -> dict:
    """Krakowiak Source code 4 / Figure 6A: two-stage run at one Hsp70-Hsf1 affinity.

    The source's protocol exactly: run 240 min at the basal load, take the reporter value at the
    end as the initial reporter for a 240-min run at the shock load. ``deactivation_min`` is the
    first minute at which transcription has fallen below ``shutoff_fraction`` of its own peak.
    """
    gain_value = gain(gain_source)
    k2 = float(HSF1_OFF_RATE) * float(multiplier)
    model = HeatShockTranscription(gain=gain_value, k2=k2)
    basal_times = np.linspace(0.0, float(minutes), int(minutes) + 1)
    basal = model.simulate(basal_times, model.initial_state(float(basal_load), reporter=3.0))[-1, 5]
    times = np.linspace(0.0, float(minutes), int(points))
    trace = model.simulate(times, model.initial_state(float(shock_load), reporter=basal))[:, 5]
    rate = np.gradient(trace, times)
    off = np.flatnonzero((rate < shutoff_fraction * rate.max()) & (times > 1.0))
    return {"time_min": times, "reporter": trace, "k2_per_min": k2,
            "basal_reporter": float(basal), "reporter_at_end": float(trace[-1]),
            "peak_rate_per_min": float(rate.max()),
            "deactivation_min": float(times[off[0]]) if off.size else None}


def gate():
    """Criterion (e), computed. It REFUSES, and that is the correct verdict.

    Fifteen free scalars against zero independent targets. Every published figure this module
    reproduces is one the parameters were selected against -- Zheng's screen chose k1-k5, beta
    and Kd to satisfy criteria 1-3, and Krakowiak re-fitted beta to Figure 1B -- so all four
    targets are registered ``fitted=True`` and none of them counts. A faithful port of a fitted
    model is evidence about the port.
    """
    return HEAT_PARAMS.gate()


def provenance() -> str:
    """One line naming what this is, so a reader never has to infer the level from the code."""
    record = transcription_record()
    return (f"TRANSCRIPTION ({record['availability']}), not an upstream model. Ported from "
            f"{len(SOURCE_DIGESTS)} checksummed MATLAB files across {_ZHENG.split(',')[0]} and "
            f"{_KRAKOWIAK.split(',')[0]}, which are ONE source: Krakowiak Source code 1 is "
            "byte-identical to Zheng's titration_YFP_FB.m. Transcription gain published twice "
            "at a factor of five and therefore SWEPT; temperature entry REFUSED; free-scalar "
            f"gate {gate().summary()}")


_SCREEN_ASSAY = ("serial 3xFLAG/V5 immunoprecipitation followed by mass spectrometry, triplicate "
                 "per time point (Zheng Figure 1C)")

HEAT_PARAMS.add_target(Target(
    name="zheng_screen_dissociation_and_recovery",
    observable="normalised [Hsp70-Hsf1] over a 25 -> 39 C shift, shocked over unshocked",
    assay=_SCREEN_ASSAY,
    noise_floor=0.10,
    units="normalised ratio",
    source=f"{_ZHENG}, Methods 'Parameter assignments and initial conditions'. THE FLOOR HERE IS "
           "THE PUBLISHED ACCEPTANCE TOLERANCE, NOT A MEASURED ASSAY CV: neither paper reports a "
           "replicate CV for the IP/MS ratio.",
    fitted=True))

HEAT_PARAMS.add_target(Target(
    name="zheng_figure_1D_reporter_plateaus",
    observable="log2 HSE-YFP fold change at four unfolded-protein loads over 120 min",
    assay="flow cytometry of an integrated 4xHSE-YFP reporter, population medians",
    noise_floor=0.10,
    units="log2 fold change",
    source=f"{_ZHENG}, Figure 1D and Fig2_YFP_Reporter.m. Fitted: Zheng's Methods state that the "
           "loads other than 39 C were chosen to 'generally capture steady state YFP reporter "
           "outputs', so the input was selected against this output.",
    fitted=True))

HEAT_PARAMS.add_target(Target(
    name="krakowiak_figure_1B_feedback_severed",
    observable="HSE-YFP in wild type against feedback-severed cells over a heat shock time course",
    assay="flow cytometry of an integrated 4xHSE-YFP reporter, three biological replicates",
    noise_floor=0.10,
    units="a.u. reporter",
    source=f"{_KRAKOWIAK}, Figure 1B and Figure 1-figure supplement 1. Fitted: beta was reduced "
           "five-fold specifically to match this separation.",
    fitted=True))

HEAT_PARAMS.add_target(Target(
    name="krakowiak_source_code_3_basal_sweep",
    observable="basal HSE-YFP fold change against the Hsp70-Hsf1 dissociation rate",
    assay="simulation only; Figure 2-figure supplement 1 plots no experimental points",
    noise_floor=0.10,
    units="fold change",
    source=f"{_KRAKOWIAK}, Source code 3. Fitted in the sense that matters here: it inherits the "
           "screened parameter set wholesale and adds no measurement of its own.",
    fitted=True))
