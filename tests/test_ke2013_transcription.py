"""Does the Ke 2013 transcription reproduce Ke 2013?

THE TEST THAT EARNS THE FILE ITS PLACE is :func:`test_alkaline_step_reproduces_the_published_shift`
plus its three companions below it. A transcription that merely runs is worth nothing: this
one is scored against four numbers Ke, Ingram & Haynes actually printed -- the unstressed
intracellular pH (7.14, Fig. 7E), the pH after a step to external pH 8.0 (7.25, same
figure), the unstressed membrane potential (~ -94 mV, Results/KCl) and Table S5's own
unstressed cation levels.

TWO OF THOSE FOUR REPRODUCE AND TWO DO NOT, and the failing pair is asserted here as a
finding with both numbers rather than deleted or tuned. Nothing in
`src/ystwin/mech/ph_ke2013.py` was fitted; the tolerances below were written after the
numbers were computed, and they are deliberately tight enough that a later edit which
changes the transcription will break them.

The remaining tests protect the things that make the first four meaningful: that the pinned
PDFs still hash to what was read, that the equations which could not be transcribed refuse
instead of defaulting, that the record on disk is what the module computes, and that this
artifact can never be promoted to an availability level that would claim upstream
provenance for our own re-implementation.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ystwin.mech.contracts import ScientificRefusal
from ystwin.mech.params import RefusedValue, Tag
from ystwin.mech import ph_ke2013 as ke

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def report():
    return ke.reproduction_report()


@pytest.fixture(scope="module")
def record():
    return ke.transcription_record()


# --------------------------------------------------------------------------------------
# The reproduction. Four published numbers; two hold, two do not, and both are reported.
# --------------------------------------------------------------------------------------

def test_alkaline_step_reproduces_the_published_shift(report):
    """Ke Fig. 7E: "intracellular pH increased from 7.14 to 7.25 upon pH 8.0 stress".

    THIS IS THE ONE THAT PASSES. The step is iso-osmotic, which is the single perturbation
    the fixed-volume reduction is legitimate for, and the mechanism the paper attributes the
    rise to -- reduced H+ uptake through the collapsed proton gradient -- is inside the block
    that was transcribable. The direction and the order of magnitude carry; 68% of the
    magnitude does, and the shortfall is stated rather than closed.
    """
    published, here = report["published_vs_here"]["alkaline_ph_i_shift"]
    assert published == pytest.approx(0.11, abs=1e-9)
    assert here > 0.0, (
        f"Ke reports intracellular pH RISING by {published:+.3f} on a step to external pH "
        f"8.0; this transcription gives {here:+.4f}. A sign disagreement would mean the "
        "proton balance was mis-transcribed, not that the levels differ")
    assert 0.5 * published < here < 1.5 * published, (
        f"published shift {published:+.3f} pH units, transcription {here:+.4f} -- "
        f"{here / published:.0%} of it")
    assert here == pytest.approx(0.0747, abs=5e-4), (
        "the reproduced shift moved; recompute reproduction_report() and say why")


def test_table_s5_cation_levels_reproduce_without_being_told_them(report):
    """Table S5's footnote states its initial concentrations ARE the unstressed steady state.

    THIS ONE PASSES TOO, and it is the strongest evidence the flux laws were read correctly:
    the solver was started from those values but is free to leave them, and 22 transporter
    parameters had to be right for it to come back within 10%.
    """
    na_published, na_here = report["published_vs_here"]["unstressed_na_mM"]
    k_published, k_here = report["published_vs_here"]["unstressed_k_mM"]
    assert abs(na_here - na_published) / na_published < 0.11, (
        f"Table S5 Init_Na = {na_published} mM, transcription {na_here:.2f} mM")
    assert abs(k_here - k_published) / k_published < 0.06, (
        f"Table S5 Init_K = {k_published} mM, transcription {k_here:.2f} mM")
    assert na_here + k_here == pytest.approx(300.0, abs=0.5), (
        "the cation total should sit on the 300 mM Anion row, because Eqn 2.1 makes Em the "
        "difference between them and 0.192 V per attomole punishes any other total")


def test_absolute_ph_does_not_reproduce_and_this_is_by_how_much(report):
    """THIS ONE IS A REPORTED FAILURE TO REPRODUCE, kept as an assertion so it cannot be
    quietly lost.

    Ke prints an unstressed intracellular pH of 7.14 and an alkaline-stressed 7.25. This
    transcription gives 7.267 and 7.342: both 0.13 pH units alkaline of the source. The
    cause remains unresolved among the transcription's declared reductions and reading
    choices. Missing buffering limits physiological interpretation but cannot explain a
    disagreement with a source that also omits buffering. No coefficient is fitted to close
    this numerical reproduction gap.
    """
    published_rest, here_rest = report["published_vs_here"]["unstressed_ph_i"]
    published_alk, here_alk = report["published_vs_here"]["alkaline_ph_i"]
    assert here_rest > published_rest, (
        "the disagreement changed sign; that is a different finding and needs re-reading")
    assert here_rest - published_rest == pytest.approx(0.127, abs=5e-3), (
        f"Ke 2013 unstressed pH_i = {published_rest}, transcription {here_rest:.4f}: "
        f"a difference of {here_rest - published_rest:+.4f} pH units")
    assert here_alk - published_alk == pytest.approx(0.092, abs=5e-3), (
        f"Ke 2013 pH 8.0 pH_i = {published_alk}, transcription {here_alk:.4f}: "
        f"a difference of {here_alk - published_alk:+.4f} pH units")


def test_membrane_potential_does_not_reproduce_and_this_is_by_how_much(report):
    """Ke, Results/KCl: "restored the membrane potential to the unstressed level (around
    -94 mV)". The transcription gives -110.9 mV under the unit-consistent reading of the
    Tok1p gating and -138.5 mV as printed, so the published number is not bracketed by the
    two readings -- it lies outside both, on the depolarised side.

    ALSO A REPORTED FAILURE. Em here is the difference between the cation total and Table
    S4's unmeasured 300 mM anion, divided by a capacitance: 0.192 V per attomole. A 17 mV
    disagreement is 0.09 attomoles out of 6032, i.e. 0.0015% of the charge. The missing
    anion measurement prevents a physiological validation, not numerical reproduction of
    equations with a declared anion value. The reproduction discrepancy remains unresolved.
    """
    published, here = report["published_vs_here"]["unstressed_em_mV"]
    as_printed = report["readings"]["as_printed_tok1"]["em_mV"]
    assert here == pytest.approx(-110.9, abs=0.5), (
        f"Ke 2013 unstressed Em ~ {published} mV, transcription {here:.2f} mV: "
        f"{here - published:+.2f} mV")
    assert as_printed < here < published, (
        f"the two readings of the Tok1p gating give {as_printed:.1f} and {here:.1f} mV and "
        f"the published {published} mV sits outside both, not between them")
    assert abs(here - published) < 20.0, "the disagreement grew; re-read the flux laws"


# --------------------------------------------------------------------------------------
# The three ambiguous printed forms, and what each reading costs.
# --------------------------------------------------------------------------------------

def test_the_printed_h_uptake_law_gives_an_impossible_cytosol(report):
    """Eqn 2.38 as typeset multiplies a J/mol bracket by a mol/(s*V) constant.

    Taken literally the model still has a fixed point, which is why this had to be computed
    rather than assumed: it is at pH_i 4.31, nearly three units below Ke's own 7.14 and
    below anything a living yeast cytosol occupies, with J_H_uptake at 6.07e3 amol/s against
    a J_H_production of 5. That is the evidence for reading the bracket in volts.
    """
    printed = report["both_forms_as_printed"]
    assert printed["has_fixed_point"] is True
    assert printed["ph_i"] == pytest.approx(4.313, abs=5e-3)
    assert printed["ph_i_below_published"] > 2.5
    assert printed["J_H_uptake_amol_per_s"] > 1e3 * float(ke.K_H_PROD)


def test_the_tok1_gating_correction_costs_27_millivolts(report):
    """Reading 2 is the expensive one, so both readings are computed and both are reported."""
    unit_consistent = report["readings"]["unit_consistent"]
    as_printed = report["readings"]["as_printed_tok1"]
    assert unit_consistent["em_mV"] - as_printed["em_mV"] == pytest.approx(27.5, abs=1.0)
    assert abs(unit_consistent["ph_i"] - as_printed["ph_i"]) < 0.02, (
        "the Tok1p reading moves the membrane potential a great deal and the pH hardly at "
        "all; if that stops being true the coupling has changed")


def test_the_unresolved_nha1_rate_split_is_bounded_not_assumed(report):
    """Eqns 2.16-2.17 need kNha1_low and kNha1_high; Table S2 prints one kNha1.

    Both are set to it, and the cost of that is measured here rather than waved away:
    driving the low-state rate from 0 to 10x moves pH_i by under 0.005 units, which is 4% of
    the 0.127 disagreement with the published value. The unresolved split is therefore not
    the explanation for the gap.
    """
    assert report["sensitivity"]["nha1_split_ph_span"] < 0.005
    assert report["sensitivity"]["nha1_split_ph_span"] > 0.0, (
        "a split that changes nothing at all would mean the low-affinity state is not wired")


# --------------------------------------------------------------------------------------
# Two internal inconsistencies the transcription found in the source.
# --------------------------------------------------------------------------------------

def test_the_volume_module_does_not_rest_at_table_s5(report):
    """Eqns 4.2-4.4 at Table S5's own initial volume give D_Pressure = -1.51e6 J/m^3, not 0.

    [Osmo_Cyt] is 1176 mM (300 mM non-ionic + 300 mM Na+/K+ + 576 mM glycerol) against
    [Osmo_Ext] 250 mM, and Table S4's turgor of 0.875e6 J/m^3 does not close a 926 mM gap.
    This is why reduction (2) holds the volume fixed, and why this module refuses the NaCl,
    sorbitol and KCl figures rather than reporting bad numbers for them.
    """
    d_pressure = report["volume_module"]["d_pressure_at_table_s5_J_per_m3"]
    assert d_pressure == pytest.approx(-1.512e6, rel=5e-3)
    assert abs(d_pressure) > report["volume_module"]["turgor_at_table_s5_J_per_m3"]


def test_osmotic_stresses_are_refused_rather_than_reported(report):
    """The reduction that makes the alkaline comparison legal makes the osmotic ones illegal,
    and the module says so by raising instead of returning a number."""
    with pytest.raises(ScientificRefusal, match="fixed-volume"):
        ke.steady_state(ke.Environment(ph_ext=6.5, k_ext_mM=801.0))


def test_table_s5_is_not_the_steady_state_of_the_calcineurin_block(report):
    """Table S5's footnote claims its rows are the unstressed steady state. For Eqn 3.4 they
    are not: with Init_CNon = 0 the activation term is unopposed, and the block relaxes to
    1.2e-5 mM activated calcineurin, carrying Crz1p to 2.5x its Table S5 row."""
    block = report["calcineurin_module"]
    assert block["converged"] is True
    assert block["table_s5_cn_mM"] == 0.0
    assert block["cn_mM"] > 1e-6
    assert block["crz1_fold_vs_table_s5"] == pytest.approx(2.54, abs=0.05)


# --------------------------------------------------------------------------------------
# The transcription itself: provenance, refusals, and what may never be claimed.
# --------------------------------------------------------------------------------------

def test_the_pinned_sources_still_hash_to_what_was_read():
    """A transcription of a byte stream that has moved is a transcription of nothing."""
    for role, (path, digest) in ke.ARTIFACTS.items():
        blob = (ROOT / path).read_bytes()
        assert hashlib.sha256(blob).hexdigest() == digest, f"{role} ({path}) has changed"


def test_availability_is_never_upstream_provenance(record):
    """The whole point of the new level. ``local_verified`` means a checksummed file the
    upstream authors published; this is our own re-implementation and must not borrow it."""
    assert ke.AVAILABILITY == "transcribed_here"
    assert record["availability"] == "transcribed_here"
    fragment = record["parameter_evidence_fragment"]["models"]["ke2013_transcription"]
    assert fragment["availability"] == "transcribed_here"
    assert fragment["source_id"] == "ke2013"
    assert fragment["executable_asset"] == "src/ystwin/mech/ph_ke2013.py"
    for metadata in (record, fragment):
        assert metadata["availability"] not in {"local_verified", "local_source_only"}


def test_every_ode_in_table_s1_is_accounted_for():
    """Twenty-three rows, none quietly dropped: 7 run, 1 is recorded, 15 refuse."""
    assert sum(ke.ODE_TOTALS.values()) == 23
    assert ke.ODE_TOTALS == {"executable": 7, "recorded": 1, "untranscribed": 15}
    assert len(ke.UNTRANSCRIBED) == 6


@pytest.mark.parametrize("label", [eq.label for eq in ke.UNTRANSCRIBED])
def test_untranscribed_equations_refuse_by_name(label):
    """The house idiom: hit the absent scalar, not a plausible default."""
    with pytest.raises(ScientificRefusal) as raised:
        ke.refuse_untranscribed(label)
    assert "not executable here" in str(raised.value)


def test_the_two_absent_scalars_refuse_rather_than_defaulting():
    """Ke supplies the structure and not these two numbers, and neither is invented here."""
    names = {p.name for p in ke.NOT_IN_SOURCE.params}
    assert names == {"cytosolic_buffer_capacity", "anion_export_capacity"}
    for param in ke.NOT_IN_SOURCE.params:
        assert param.tag == Tag.REFUSED
        with pytest.raises(RefusedValue):
            float(param)
    buffering = next(p for p in ke.NOT_IN_SOURCE.params
                     if p.name == "cytosolic_buffer_capacity")
    assert "ZERO times" in buffering.reason


def test_the_anion_is_asserted_with_its_empty_reference_column_said_out_loud():
    """Table S4's Anion row is the only row in that table with no reference, and it sets the
    entire membrane potential. It is ASSERTED, never MEASURED, and it names what is missing."""
    assert ke.ANION_MM.tag == Tag.ASSERTED
    assert float(ke.ANION_MM) == 300.0
    assert "EMPTY REFERENCE COLUMN" in ke.ANION_MM.source
    assert ke.ANION_MM.missing


def test_every_parameter_carries_ke_s_own_provenance_column_verbatim():
    """"Taken from" against "Estimated from", the asterisk on the fitted rows, and an
    explicit note where Ke's reference column is blank. Losing that column would turn a
    fitted integration constant into an anonymous literal."""
    saved = json.loads((ROOT / ke.RECORD_PATH).read_text())["parameters"]
    provenance = {param["name"]: param["provenance_verbatim"] for param in saved}
    assert len(provenance) == len(saved)
    assert set(provenance) == {param.name for param in ke.KE2013_PARAMS.params}
    for param in ke.KE2013_PARAMS.params:
        assert param.tag == Tag.ASSERTED, (
            f"{param.name} is graded {param.tag}; every Ke row is ASSERTED because the "
            "upstream references behind 'Taken from' were not traced here")
        assert param.units.strip()
        assert "Ke 2013 Table" in param.source
        assert param.source == provenance[param.name], (
            f"{param.name} lost Ke's recorded provenance column")
    starred = [p for p in ke.KE2013_PARAMS.params if "ASTERISKED" in p.source]
    assert len(starred) >= 12, (
        "Tables S2 and S3 asterisk the rows adjusted during Ke's own model-integration step; "
        "those are fits and the record has to keep saying so")


def test_the_record_on_disk_is_what_the_module_computes(record, tmp_path):
    """The JSON is generated, never hand-edited, so it cannot drift from the code."""
    on_disk = json.loads((ROOT / ke.RECORD_PATH).read_text())
    assert on_disk["parameters"] == record["parameters"]
    assert on_disk["equations"] == record["equations"]
    assert on_disk["reading_choices"] == record["reading_choices"]
    assert on_disk["artifacts"] == record["artifacts"]
    for key in ("unstressed_ph_i", "unstressed_em_mV", "alkaline_ph_i_shift"):
        written = on_disk["reproduction"]["published_vs_here"][key]
        computed = record["reproduction"]["published_vs_here"][key]
        assert written[0] == pytest.approx(computed[0])
        assert written[1] == pytest.approx(computed[1], rel=1e-6)
    assert on_disk["reproduction"]["both_forms_as_printed"] == pytest.approx(
        record["reproduction"]["both_forms_as_printed"], rel=1e-8)
    assert on_disk["reproduction"]["acidic_limit"] == record["reproduction"]["acidic_limit"]
    assert on_disk["reproduction"]["numerical_contract"] == record["reproduction"]["numerical_contract"]
    written_path = ke.write_transcription_record(tmp_path / "ke2013.json")
    assert json.loads(written_path.read_text())["availability"] == "transcribed_here"


def test_the_steady_state_is_actually_a_steady_state():
    """The root find is asserted rather than trusted: all three balances at zero."""
    env = ke.Environment()
    state = ke.steady_state(env)
    d_h, d_na, d_k = ke.ion_rhs(state, env)
    np.testing.assert_allclose((d_h, d_na, d_k), 0.0, atol=1e-9, rtol=0.0)
    assert ke.membrane_potential_V(state) == pytest.approx(-0.1109, abs=5e-4)


def test_the_flux_signs_match_the_biology_ke_describes():
    """A transcription with a sign error can still land on a fixed point, so the directions
    are checked separately: Pma1p extrudes protons, the Trk system takes K+ in, Tok1p lets
    it out, and Ena1p and Nha1p both export Na+."""
    env = ke.Environment()
    state = ke.steady_state(env)
    flux = ke.fluxes(state, env)
    assert flux["J_Pma1"] > 0.0
    assert flux["J_Trk_K"] > 0.0
    assert flux["J_Tok1"] > 0.0
    assert flux["J_Ena1_Na"] > 0.0
    assert flux["J_Nha1_Na"] > 0.0
    assert flux["Em"] < 0.0
    assert 0.0 < flux["P_Tok1_O"] < 1.0
    assert flux["J_Pma1"] > flux["J_H_production"], (
        "at rest Pma1p must carry more than the constant production, because Eqn 2.38's "
        "uptake returns protons the pump has to move again")


def test_ph_is_an_output_not_an_input():
    """The reason this source was worth transcribing for the 'ph' programme at all: nothing
    sets pH_i, it falls out of the proton balance, and it moves when the medium moves."""
    acidic = ke.steady_state(ke.Environment(ph_ext=5.5))
    neutral = ke.steady_state(ke.Environment(ph_ext=6.5))
    alkaline = ke.steady_state(ke.Environment(ph_ext=7.5))
    assert acidic.ph_i < neutral.ph_i < alkaline.ph_i
    assert not math.isclose(acidic.ph_i, alkaline.ph_i, abs_tol=0.02)


@pytest.mark.parametrize("voltage", [-1e-8, 0.0, 1e-8])
def test_tok1_current_has_the_analytic_zero_voltage_limit(voltage):
    initial = ke.IonState.from_table_s5()
    charge = voltage * ke._CAPACITANCE_F / (ke._AMOL * ke._F)
    total = float(ke.ANION_MM) * initial.volume_um3 + charge
    state = ke.IonState(initial.h_amol, total - initial.h_amol - initial.k_amol,
                        initial.k_amol, initial.volume_um3)
    flow = ke.fluxes(state, ke.Environment())
    fraction = float(ke.KM_TOK1_HOG1) / (float(ke.INIT_HOG1PPC_MM) + float(ke.KM_TOK1_HOG1))
    expected = float(ke.PS_TOK1) * fraction * flow["P_Tok1_O"] * (state.k_mM - 1.0)
    assert flow["J_Tok1"] == pytest.approx(expected, rel=1e-6)


def test_an_optimizer_success_flag_cannot_claim_a_root_or_its_nonexistence(monkeypatch):
    monkeypatch.setattr(ke, "least_squares", lambda fun, x0, **kwargs: SimpleNamespace(
        x=np.array(x0), success=True), raising=False)
    with pytest.raises(ScientificRefusal) as failure:
        ke.steady_state(ke.Environment(ph_ext=5.5))
    assert getattr(failure.value, "status", None) == "numerical_nonconvergence"
    assert failure.value.nonexistence_proven is False
    assert failure.value.best_scaled_residual > 1.0


def test_osmotic_reproduction_is_refused_by_scope_before_solving(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("the fixed-volume reproduction must reject an osmotic protocol first")

    monkeypatch.setattr(ke, "root", forbidden, raising=False)
    monkeypatch.setattr(ke, "least_squares", forbidden, raising=False)
    with pytest.raises(ScientificRefusal, match="fixed-volume"):
        ke.steady_state(ke.Environment(k_ext_mM=801.0))


def test_a_solver_convergence_boundary_is_not_a_physiological_acid_limit():
    with pytest.raises(ScientificRefusal, match="not identified"):
        ke.acidic_limit()


@pytest.mark.parametrize("ph_ext,reading", [(6.5, ke.UNIT_CONSISTENT),
                                            (5.5, ke.UNIT_CONSISTENT),
                                            (6.5, ke.AS_PRINTED)])
def test_equilibria_remain_stationary_under_independent_time_integration(ph_ext, reading):
    step = ke.alkaline_step(ph_ext_from=ph_ext, ph_ext_to=ph_ext, reading=reading)
    assert step["ph_i_after"] == pytest.approx(step["ph_i_before"], abs=1e-7)
    assert step["em_mV_after"] == pytest.approx(step["em_mV_before"], abs=1e-5)
    for field in ("h_amol", "na_amol", "k_amol"):
        assert getattr(step["after"], field) == pytest.approx(
            getattr(step["before"], field), rel=1e-7)
