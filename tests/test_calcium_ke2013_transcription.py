"""Does the wired calcium axis reproduce Ke 2013, and does it refuse what Ke 2013 cannot say?

THE TEST THAT EARNS THIS FILE ITS PLACE is
:func:`test_table_s3_rate_constants_reproduce_table_s5_initial_concentrations`. Ke, Ingram &
Haynes publish the calcium block's rate constants in Table S3 (pinned ``s011.pdf``) and its
initial concentrations in Table S5 (pinned ``s013.pdf``), and Table S5's footnote states those
initials ARE the unstressed steady state. So each rest point implied by Table S3 alone is an
independent target sitting in a different file. Two of the three land on it:

    [Ca2+]  = CCa/dCa     = 5.000000e-05 mM  against Table S5's 5.0e-05 mM   (exact)
    [Crz1p] = CCrz1/dCrz1 = 1.918537e-04 mM  against Table S5's 1.916e-04 mM (+0.13%)

and the third does NOT: Table S5's ``Init_CNon = 0`` is not the rest point of Eqn 3.4, which is
asserted here as a finding with its number rather than tuned away. Nothing in
`src/ystwin/mech/calcium_ke2013.py` was fitted, and the tolerances below were written after the
numbers were computed.

THE SECOND THING THIS FILE PROTECTS is that the wiring did not quietly retype anything. The
module imports Ke's equations from `mech/ph_ke2013.py` instead of re-transcribing them, with two
exceptions -- Eqn 2.21 (Ppz) and Eqn 2.24 (the Trk affinity split) live inline inside
``ph_ke2013.fluxes`` and are not exported. Both retypes are checked against that audited
implementation here, the second by solving the audited Trk flux at two calcineurin levels, so a
mistyped constant fails rather than propagates.

THE THIRD is the domain. Eqn 3.3's alkaline driver is linear and unbounded below, so beneath
external pH 6.3327 the equation has no non-negative rest point at all -- and every protocol in
this repository runs at pH 4-5. The axis must REFUSE there, not clip, and that is asserted.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import pytest

from ystwin.mech.contracts import ScientificRefusal
from ystwin.mech.params import RefusedValue, Tag
from ystwin.mech import calcium_ke2013 as calcium
from ystwin.mech import ph_ke2013 as ke

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def report():
    return calcium.reproduction_report()


# --------------------------------------------------------------------------------------
# The reproduction: Table S3's constants against Table S5's concentrations.
# --------------------------------------------------------------------------------------

def test_table_s3_rate_constants_reproduce_table_s5_initial_concentrations(report):
    """Two pinned supplementary PDFs, read independently, agree to five figures and to 0.13%.

    This is the whole evidence that the fifteen Table S3 constants were read off the page
    correctly: a mistyped rate constant would not land on a third file's number. It says nothing
    about whether Ke's model is right, only that this repository is running Ke's model.
    """
    cross = report["cross_check"]
    assert cross["ca"]["from_table_s3"] == pytest.approx(cross["ca"]["table_s5"], rel=1e-12)
    assert abs(cross["ca"]["relative_difference"]) < 1e-12
    assert cross["crz1"]["from_table_s3"] == pytest.approx(1.918537e-4, rel=1e-6)
    assert cross["crz1"]["from_table_s3"] == pytest.approx(cross["crz1"]["table_s5"], rel=2e-3)
    assert 0.0 < cross["crz1"]["relative_difference"] < 2e-3


def test_the_crz1_agreement_is_the_rest_point_ke_intended(report):
    """CCrz1/dCrz1 is Eqn 3.4's rest point at CN = 0, which is Table S5's own Init_CNon.

    This is what sharpens the third row's failure below from an oddity into a contradiction: Ke
    plainly intended CN = 0 as the unstressed rest point -- his Crz1p row is the Crz1p that CN = 0
    produces -- and his own Eqn 3.4 refuses it.
    """
    cross = report["cross_check"]
    assert cross["calcineurin"]["table_s5_cn_mM"] == 0.0
    at_zero_calcineurin = (float(ke.C_CRZ1) + float(ke.K_CRZ1) * 0.0) / float(ke.D_CRZ1)
    assert at_zero_calcineurin == pytest.approx(cross["crz1"]["from_table_s3"], rel=1e-12)


def test_table_s5_is_not_the_rest_point_of_the_calcineurin_pair_and_is_not_tuned(report):
    """The third cross-check FAILS, and carrying that forward unchanged is the point.

    ``d[CNoff]/dt = -1.4535e-07 mM/s`` at Ke's own published initials. The same defect is already
    pinned at ``tests/test_ke2013_transcription.py:206`` through the 2.54x Crz1p relaxation; it is
    re-asserted from this side so a future edit cannot close it here and leave that one passing.
    """
    cross = report["cross_check"]
    assert cross["calcineurin"]["reproduces"] is False
    assert cross["calcineurin"]["d_cn_off_per_s_at_table_s5"] == pytest.approx(-1.4535e-7, rel=1e-6)
    assert report["alkaline"]["before"]["crz1_fold_vs_table_s5"] == pytest.approx(2.535, abs=0.02)
    assert "must not be tuned away" in calcium.NOT_APPLIED["table_s5_init_cnon_zero"]


def test_the_pinned_sources_still_hash_to_what_was_read():
    """Wiring a transcription of a byte stream that has moved is wiring nothing."""
    assert calcium.ARTIFACTS is ke.ARTIFACTS
    for role, (path, digest) in calcium.ARTIFACTS.items():
        blob = (ROOT / path).read_bytes()
        assert hashlib.sha256(blob).hexdigest() == digest, f"{role} ({path}) has changed"


# --------------------------------------------------------------------------------------
# The two retypes, checked against the audited implementation they were copied from.
# --------------------------------------------------------------------------------------

def test_ppz_matches_the_audited_copy_of_eqn_2_21():
    """Eqn 2.21 is retyped here because ``ph_ke2013.fluxes`` does not export it. Verify it."""
    ion = ke.IonState.from_table_s5()
    audited = ke.fluxes(ion, ke.Environment())["Ppz"]
    assert calcium.ppz_mM(ion.ph_i) == pytest.approx(audited, rel=1e-12)
    acidic, alkaline = calcium.ppz_mM(6.0), calcium.ppz_mM(8.0)
    assert acidic < alkaline < float(ke.PPZ_TOTAL)


def test_trk_split_matches_the_audited_copy_of_eqn_2_24():
    """Solve the audited Trk flux at two calcineurin levels and recover the split from it.

    ``J_Trk_K = (A*p_medium + B*p_high) * kTrk * Em^2`` with A and B fixed by the medium and the
    membrane potential, so two calls to the audited ``fluxes`` at one fixed ion state determine
    the retyped split exactly. This is a check on the copy, not on Ke.
    """
    ion = ke.IonState.from_table_s5()
    env = ke.Environment()
    ppz = ke.fluxes(ion, env)["Ppz"]
    cn_high = 5e-4
    at_zero = ke.fluxes(ion, env, cn_mM=0.0)
    at_high = ke.fluxes(ion, env, cn_mM=cn_high)
    medium_zero, high_zero = calcium.trk_affinity_split(0.0, ppz)
    medium_high, high_high = calcium.trk_affinity_split(cn_high, ppz)
    # Two linear equations in the two unknown medium/high affinity coefficients of Eqn 2.24.
    determinant = medium_zero * high_high - medium_high * high_zero
    assert abs(determinant) > 1e-12
    coefficient_medium = (at_zero["J_Trk_K"] * high_high - at_high["J_Trk_K"] * high_zero) / determinant
    coefficient_high = (medium_zero * at_high["J_Trk_K"] - medium_high * at_zero["J_Trk_K"]) / determinant
    for cn in (0.0, 1e-5, 5e-4, 5e-3):
        medium, high = calcium.trk_affinity_split(cn, ppz)
        predicted = coefficient_medium * medium + coefficient_high * high
        assert predicted == pytest.approx(ke.fluxes(ion, env, cn_mM=cn)["J_Trk_K"], rel=1e-9)


def test_the_split_is_a_partition_and_calcineurin_moves_it_toward_high_affinity():
    """Ke's Table 1: calcineurin "switches Trk system from a medium affinity state to a high
    affinity state". The retyped term must have that sign, and must sum to one."""
    ppz = calcium.ppz_mM(7.0)
    previous = -math.inf
    for cn in (0.0, 1e-5, 1e-4, 1e-3):
        medium, high = calcium.trk_affinity_split(cn, ppz)
        assert medium + high == pytest.approx(1.0, rel=1e-12)
        assert high > previous
        previous = high


# --------------------------------------------------------------------------------------
# The domain: where Eqn 3.3 stops being defined, and what the engine's pH window means.
# --------------------------------------------------------------------------------------

def test_the_axis_refuses_below_the_ph_where_eqn_3_3_has_no_nonnegative_rest_point(report):
    """Below external pH 6.3327 Ke's own printed constants drive cytosolic Ca2+ negative.

    ``kCa,pH`` = 1.5e-5 mM/s per pH unit against ``CCa`` = 2.5e-6 mM/s, so the linear alkaline
    driver overwhelms the only zero-order source 0.167 pH units below 6.5. Ke never simulates an
    acidic condition and prints no floor, so the axis refuses rather than clipping to one.
    """
    assert report["ph_production_floor"] == pytest.approx(6.3327, abs=1e-3)
    axis = calcium.CalciumAxis()
    state = axis.initial_state()
    for ph in (3.0, 4.0, 5.0, 6.0, 6.33):
        with pytest.raises(ScientificRefusal, match="no non-negative rest point"):
            axis.rates_per_h(state, medium_ph=ph, cytosolic_ph=7.0)
    for ph in (6.5, 7.0, 8.0):
        rates = axis.rates_per_h(state, medium_ph=ph, cytosolic_ph=7.0)
        assert len(rates) == 4 and all(math.isfinite(rate) for rate in rates)


def test_the_axis_is_undefined_at_this_repositorys_own_operating_ph(report):
    """The honest headline: ``Control.ph`` is gated [3, 8] and every protocol here runs at 4-5.

    That is below the floor, so wiring this axis into the engine gives an axis that refuses on
    every existing protocol and comes alive only under a large ``sodium_hydroxide`` dose. A
    wiring that clamped pH upward to keep it alive would be inventing the constant Ke omits.
    """
    low_gate, high_gate = report["engine_ph_gate"]
    assert (low_gate, high_gate) == (3.0, 8.0)
    assert low_gate < report["ph_production_floor"] < high_gate
    assert report["validated_ph"] == (6.5, 8.0)
    assert report["ph_production_floor"] < 6.5
    with pytest.raises(ScientificRefusal, match="pH 4-5"):
        calcium.CalciumAxis().check_ph(5.0)


def test_the_unvalidated_band_between_the_floor_and_ke_can_be_closed():
    """Between the floor and Ke's own pH 6.5 the equations are defined but never simulated.

    The default lets that band run and says so; ``allow_outside_validated_ph=False`` closes it.
    Neither setting opens anything below the floor, which is a matter of sign, not validation.
    """
    strict = calcium.CalciumAxis(allow_outside_validated_ph=False)
    assert strict.check_ph(7.0) == 7.0
    with pytest.raises(ScientificRefusal, match="outside Ke 2013's simulated window"):
        strict.check_ph(6.4)
    with pytest.raises(ScientificRefusal, match="no non-negative rest point"):
        strict.check_ph(6.0)
    assert calcium.CalciumAxis().check_ph(6.4) == 6.4


# --------------------------------------------------------------------------------------
# The refusals: what this axis may not be asked, and what it may not be said to answer.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["extracellular_calcium_mM", "cytosolic_sodium_mM",
                                  "external_sodium_mM", "crz1_to_cdre_gain"])
def test_every_missing_engine_channel_raises_rather_than_defaulting(name):
    """A refused parameter has no float, and says what is missing and what would close it."""
    param = calcium.NOT_DRIVEN[name]
    assert param.tag == Tag.REFUSED
    assert param.value is None
    assert param.reason.strip() and param.missing.strip() and param.source.strip()
    with pytest.raises(RefusedValue):
        float(param)


def test_the_calcium_chloride_stressor_is_forfeited_not_answered():
    """Eqn 3.3 has no extracellular-calcium term, so the axis's heaviest stressor is unanswerable.

    Asserted structurally rather than in prose: the refusal exists, the panel declaration lists
    calcium_chloride as forfeited at its full weight, and it is not in the answered set.
    """
    channel = calcium.PANEL_CHANNEL
    assert channel["module"] == "calcium"
    assert channel["answered"] == {"sodium_hydroxide": 0.35}
    assert channel["forfeited"]["calcium_chloride"] == 1.00
    assert channel["forfeited"]["NaCl"] == 0.25
    assert "calcium_chloride" not in channel["answered"]
    assert "no extracellular-calcium term" in calcium.NOT_DRIVEN["extracellular_calcium_mM"].source


def test_the_panel_declaration_matches_the_panel_itself():
    """The forfeited and answered weights must be the panel's actual numbers, not remembered ones."""
    from ystwin.generator.stress_panel import STRESSORS
    for stressor, weight in dict(calcium.PANEL_CHANNEL["answered"]).items():
        assert STRESSORS[stressor].targets["calcium"] == pytest.approx(weight)
    for stressor, weight in dict(calcium.PANEL_CHANNEL["forfeited"]).items():
        assert STRESSORS[stressor].targets["calcium"] == pytest.approx(weight)
    routed = {name for name, s in STRESSORS.items() if "calcium" in s.targets}
    declared = set(calcium.PANEL_CHANNEL["answered"]) | set(calcium.PANEL_CHANNEL["forfeited"])
    assert routed == declared, "a stressor routed to the calcium module is unaccounted for"


def test_the_two_sodium_arms_are_held_constants_and_not_engine_signals():
    """They run at Ke's own unstressed baseline, and holding them costs 0.37% of the Ca2+ source.

    Measured rather than assumed: the two arms together contribute 9.3e-09 mM/s against CCa's
    2.5e-06 mM/s at rest. Small is not the same as driven, which is why they are still refused.
    """
    axis = calcium.CalciumAxis()
    assert axis.na_cyt_mM == float(ke.INIT_NA_MM) == 100.0
    assert axis.na_ext_mM == 5.0
    h_c, h_e = float(ke.H_NA_CYT), float(ke.H_NA_EXT)
    cyt = float(ke.K_CA_CYT) * 100.0 ** h_c / (100.0 ** h_c + float(ke.KM_CA_CYT) ** h_c)
    ext = float(ke.K_CA_EXT) * 5.0 ** h_e / (5.0 ** h_e + float(ke.KM_CA_EXT) ** h_e)
    assert (cyt + ext) / float(ke.C_CA) == pytest.approx(0.0037, abs=5e-4)
    state = axis.initial_state()
    base = axis.rates_per_h(state, medium_ph=7.0, cytosolic_ph=7.0)[0]
    doubled = calcium.CalciumAxis(na_cyt_mM=200.0).rates_per_h(
        state, medium_ph=7.0, cytosolic_ph=7.0)[0]
    assert doubled != base


def test_eqn_3_2_stress_onset_spike_is_recorded_as_not_applied():
    """Ke's discontinuous Ca2+ reset has no trigger time in a continuously-controlled vessel."""
    note = calcium.NOT_APPLIED["eqn_3_2_stress_onset_calcium_spike"]
    assert "no stress-onset instant" in note
    assert "no transient calcium amplitude may be claimed" in note


def test_the_wiring_cannot_promote_the_transcriptions_availability():
    """Wiring an existing transcription is not new provenance, and the label must not drift."""
    assert calcium.AVAILABILITY == ke.AVAILABILITY == "transcribed_here"
    assert calcium.SOURCE is ke.SOURCE
    assert "CC-BY-4.0" in calcium.SOURCE


def test_the_axis_registers_no_new_free_scalar():
    """Every constant read here is already counted in ``KE2013_PARAMS``; refusals carry no value.

    A second registry over the same rows would double-count them in the free-scalar gate, so the
    only registry this module owns holds refusals.
    """
    for name, param in calcium.USED_PARAMETERS.items():
        assert ke.KE2013_PARAMS[name] is param
        assert param.tag == Tag.ASSERTED
    assert len(calcium.USED_PARAMETERS) == 24
    assert all(param.tag == Tag.REFUSED for param in calcium.NOT_DRIVEN.params)
    assert not set(calcium.USED_PARAMETERS) & {p.name for p in calcium.NOT_DRIVEN.params}


# --------------------------------------------------------------------------------------
# The axis as an axis: units, states, and the back-coupling that stops it being an appendage.
# --------------------------------------------------------------------------------------

def test_rates_are_per_hour_and_are_the_imported_right_hand_side_times_3600():
    """The engine's time base is hours and Ke's is seconds; nothing else changes in between."""
    axis = calcium.CalciumAxis()
    state = axis.initial_state()
    ppz = calcium.ppz_mM(7.2)
    per_h = axis.rates_per_h(state, medium_ph=7.0, ppz_mM_value=ppz)
    per_s = ke.calcineurin_rhs(*state.as_tuple(), ppz_mM=ppz, na_cyt_mM=100.0,
                               env=ke.Environment(ph_ext=7.0), v_ratio_per_s=0.0)
    assert calcium.SECONDS_PER_HOUR == 3600.0
    for hourly, secondly in zip(per_h, per_s):
        assert hourly == pytest.approx(secondly * 3600.0, rel=1e-12)


def test_the_four_states_are_declared_and_start_at_table_s5():
    names = tuple(name for name, _units, _compartment, _role in calcium.STATE_VARIABLES)
    assert names == ("ca_mM", "cn_off_mM", "cn_mM", "crz1_mM")
    assert all(units == "mM" for _n, units, _c, _r in calcium.STATE_VARIABLES)
    state = calcium.CalciumState.from_table_s5()
    assert state.as_tuple() == (float(ke.INIT_CA_MM), float(ke.INIT_CNOFF_MM),
                                float(ke.INIT_CN_MM), float(ke.INIT_CRZ1_MM))
    assert state.calcineurin_total_mM == pytest.approx(1.1628e-3, rel=1e-12)
    with pytest.raises(ValueError):
        calcium.CalciumState(-1e-9, 1e-3, 0.0, 1e-4)


def test_the_alkaline_response_is_this_modules_output_and_not_a_published_number(report):
    """Ke prints no calcineurin or Crz1p value anywhere, so these folds are pinned, not claimed.

    They are asserted only so a drift in the wiring fails a test. What they show is that the axis
    is not inert: half the calcineurin pool activates across Ke's own pH 6.5 -> 8.0 step.
    """
    alkaline = report["alkaline"]
    assert alkaline["ca_fold"] == pytest.approx(4.667, rel=5e-3)
    assert alkaline["cn_fold"] == pytest.approx(48.23, rel=5e-3)
    assert alkaline["crz1_fold"] == pytest.approx(29.58, rel=5e-3)
    assert alkaline["activated_fraction_before"] == pytest.approx(0.0104, abs=5e-4)
    assert alkaline["activated_fraction_after"] == pytest.approx(0.4997, abs=5e-3)


def test_calcineurin_reaches_the_ion_system_through_the_trk_affinity_split(report):
    """The back-coupling is real and is measured, so the axis is not a dangling appendage."""
    alkaline = report["alkaline"]
    before_high = alkaline["before"]["trk_split"][1]
    after_high = alkaline["after"]["trk_split"][1]
    assert before_high == pytest.approx(0.2067, abs=2e-3)
    assert after_high == pytest.approx(0.3396, abs=2e-3)
    assert after_high > before_high


def test_fk506_reproduces_ke_figure_7as_direction_and_smallness_and_is_not_scored(report):
    """Ke: intracellular Na+ shows "no notable difference" with and without FK506 at pH 8.0.

    Recomputed here as +10.0% in Na+ and -5.7% in K+, against a 48-fold calcineurin activation
    driving them. The published claim is qualitative -- "notable" is not a number in the source --
    so what is asserted is that the effect stays an order of magnitude below its driver and that
    blocking calcineurin moves Trk toward medium affinity, raising Na+. It is reported, not scored.
    """
    fk = report["fk506"]
    assert fk["ph_ext"] == 8.0
    assert fk["relative_na_change"] == pytest.approx(0.100, abs=5e-3)
    assert fk["relative_k_change"] == pytest.approx(-0.0575, abs=5e-3)
    assert abs(fk["relative_na_change"]) < 0.15
    assert fk["trk_high_affinity"]["fk506"] < fk["trk_high_affinity"]["wild_type"]
    assert fk["na_cyt_mM"]["fk506"] > fk["na_cyt_mM"]["wild_type"]
    assert "not scored" in fk["published_claim"]


def test_the_coupled_fixed_point_refuses_rather_than_reporting_a_non_converged_pair():
    """A fixed-point iteration that has not converged is not a steady state."""
    with pytest.raises(ScientificRefusal, match="did not converge"):
        calcium.coupled_state(8.0, tolerance=0.0, max_iterations=1)


def test_rates_require_exactly_one_of_cytosolic_ph_or_ppz():
    """Two disagreeing Ppz values in one call is how a silent inconsistency gets in."""
    axis = calcium.CalciumAxis()
    state = axis.initial_state()
    with pytest.raises(ValueError, match="exactly one"):
        axis.rates_per_h(state, medium_ph=7.0)
    with pytest.raises(ValueError, match="exactly one"):
        axis.rates_per_h(state, medium_ph=7.0, cytosolic_ph=7.0, ppz_mM_value=1e-5)
