"""Does the Talemi 2016 transcription reproduce Talemi 2016, and does it refuse what it cannot?

THE TESTS THAT EARN THE FILE ITS PLACE are
:func:`test_table_s6_closed_forms_reproduce_the_deposited_initial_state`,
:func:`test_table_s8_calculated_rows_reproduce_from_their_published_formulas` and
:func:`test_the_published_protocol_reproduces_the_four_minute_slt2_peak`. A transcription that
merely runs is worth nothing. This one is scored against Table S6's published closed forms
(twelve initial-state rows plus two quadratic roots), Table S8's five "Calculated" rows, and the
article's own experiment run at Table S7's own five protocol numbers.

WHAT REPRODUCES: the twelve initial-state rows to a worst relative 1.7e-15, the two roots to
3.0e-12 and 1.2e-16, the five calculated rows to 1.9e-15 (four of them) and 3.0e-12 (k13, which
reads a root), the deposited state as a fixed point of eleven of twelve states at 2.4e-14 /s,
and the published protocol's whole shape -- flat at 45 s and 90 s, a real 1.58x peak at 4 min
below the pre-stress volume, and 89.2%/84.7%/84.1% at 14/30/45 min above it. Nothing in
`src/ystwin/mech/cell_wall_talemi2016.py` was fitted; the tolerances here were written after the
numbers were computed and are tight enough that an edit which changes the transcription breaks
them.

THE SIGN TEST IS A NEGATIVE CONTROL AS MUCH AS A REPRODUCTION.
:func:`test_slt2_answers_a_volume_increase_and_hog1_answers_a_decrease` fails if anyone ever
"fixes" this axis into a cell-wall-damage sensor, because a damage sensor does not have this
sign. :func:`test_the_sensitizer_switch_is_not_the_paper_ranking` pins the reason the deposit's
own switch cannot stand in for the paper's model comparison, so that a later reader cannot quote
the flat peaks as a reproduction of Table 1.

THE REFUSALS ARE TESTED THE SAME WAY THE REPRODUCTIONS ARE: every one raises on ``float()``, the
strings the source does not contain are checked to be genuinely absent from the pinned artifacts
rather than merely asserted to be, and the availability level cannot be promoted.

FOUR TAMPER TESTS stage a modified copy of the deposit, re-pin its hash so the audit actually
reaches the rules, and drive four different branches: a changed constant, a changed initial
amount, a changed switch and a changed compartment size all raise. A fifth leaves the hash
pinned and confirms that, without the re-pin, nothing gets past the checksum at all.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest

from ystwin.mech import cell_wall_talemi2016 as cw
from ystwin.mech.params import RefusedValue, Tag

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def report():
    return cw.reproduction_report()


@pytest.fixture(scope="module")
def model():
    return cw.OsmoStat()


# --------------------------------------------------------------------------------------
# The reproduction.
# --------------------------------------------------------------------------------------

def test_table_s6_closed_forms_reproduce_the_deposited_initial_state(report):
    """Twelve published formulas against twelve stored numbers, to floating point.

    Table S6 prints the initial conditions as formulas over Table S7's abundances and Table S6's
    two data-derived baselines. Evaluating them here has to land on the deposit's stored
    ``initialAmount`` values, and a misread formula cannot do that by accident at 1e-15.
    """
    rows = report["table_s6_closed_forms"]
    assert set(rows) >= {"Slt2", "Slt2P", "Hog1", "Hog1PP", "Fps1", "Vos0", "Vex",
                         "f_N2uM", "ci0", "cin0", "VP0", "A0"}
    assert report["worst_initial_state_relative"] < 5e-15
    for name in ("Slt2", "Slt2P", "Hog1", "Hog1PP", "Fps1", "f_N2uM", "ci0", "VP0", "A0"):
        assert rows[name]["relative"] < 5e-15, name


def test_the_slt2_row_grouping_is_the_one_the_deposit_resolves():
    """``1 - [Slt2PP]0*fn``, not ``(1 - [Slt2PP]0)*fn``. The PDF loses the grouping.

    The alternative reading is not absurd -- it is the obvious one -- and it is wrong by 25%,
    which the deposit settles. Pinned so a later reader cannot re-derive it the other way.
    """
    y0 = cw.initial_state()
    geometry = cw._derived_geometry()
    slt2_total = float(cw.SLT2_MOLECULES) * geometry["f_N2uM"]
    fn, rel = float(cw.F_NUC), float(cw.SLT2P_REL_INI)
    deposited = 4.30809697774826
    kept = slt2_total * (1.0 - rel * fn) * geometry["Vos0"]
    rejected = slt2_total * (1.0 - rel) * fn * geometry["Vos0"]
    assert y0["Slt2"] == pytest.approx(kept, rel=1e-15)
    assert kept == pytest.approx(deposited, rel=1e-14)
    assert abs(rejected - deposited) / deposited == pytest.approx(0.249, abs=0.01)


def test_table_s6_signal_roots_reproduce(report):
    """Both quadratic roots, and the reason one of them is only 3.0e-12."""
    roots = report["signal_roots"]
    assert roots["HOGSignal"]["relative"] < 1e-14
    assert roots["Slt2Signal"]["relative"] < 1e-11
    assert roots["Slt2Signal"]["published_formula"] == pytest.approx(0.244923861463865, rel=1e-11)


def test_table_s8_calculated_rows_reproduce_from_their_published_formulas(report):
    """k4, km7, k6b, k13 and k15: five rate laws set to zero, five stored numbers hit.

    This is the check that exercises the most of the transcription at once, because each row is
    a different block's steady-state condition. k13 inherits Slt2Signal|0's 3.0e-12 and is
    allowed exactly that, not a round number chosen to make it pass.
    """
    rows = report["table_s8_calculated_rows"]
    assert set(rows) == {"k4", "km7", "k6b", "k13", "k15"}
    for name in ("k4", "km7", "k6b", "k15"):
        assert rows[name]["relative"] < 5e-15, name
    assert rows["k13"]["relative"] < 1e-11
    assert rows["k13"]["relative"] > 1e-13, (
        "k13 landing at machine precision would mean it no longer reads Slt2Signal|0")


def test_k6b_is_calculated_and_table_s8_says_estimated(report):
    """The Method column's k6/k6b cells are transposed, shown rather than asserted."""
    rows = report["table_s8_calculated_rows"]
    assert rows["k6b"]["relative"] < 5e-15
    assert "transposed" in cw.SOURCE_DEFECTS["k6b_method_column_transposed"]


def test_the_deposited_initial_state_is_a_fixed_point_of_eleven_of_twelve_states(report, model):
    """Eleven states at rest to 2.4e-14 /s. The twelfth is a published defect, not ours."""
    assert report["worst_relative_rate_excluding_glyex"] < 1e-13
    residual = report["steady_state_residual"]
    assert residual["Glyex"]["rate_per_s"] == pytest.approx(4542.57, rel=1e-3)
    assert residual["Glyex"]["rate_per_s"] > 0.0, "external glycerol only accumulates"
    for name in cw.STATE_NAMES:
        if name != "Glyex":
            assert residual[name]["relative_per_s"] < 1e-13, name


def test_the_external_glycerol_drift_is_bounded_and_measured(report):
    """24 unstressed hours: Glyex 1800 -> 9646 uM, and the eleven others do not move."""
    drift = report["unstressed_24h_drift"]
    start, end = drift["glyex_uM"]
    assert start == pytest.approx(1800.0, rel=1e-6)
    assert end == pytest.approx(9646.0, rel=1e-3)
    assert end < float(cw.CE0) / 20.0, "the drift must stay small against ce0 = 260000 uM"
    assert drift["vos_relative_change"] < 1e-4
    assert drift["slt2p_relative_percent"][1] == pytest.approx(24.6, abs=1e-3)


def test_the_two_published_baselines_reproduce_exactly(report):
    """24.6% Slt2 and 6.1% Hog1, through the deposit's own Fit* readout rules."""
    baselines = report["published_baselines"]
    assert baselines["slt2p_relative_percent"] == pytest.approx(24.6, abs=1e-9)
    assert baselines["hog1pp_relative_percent"] == pytest.approx(6.1, abs=1e-9)
    assert baselines["hog1_total_uM"] == pytest.approx(cw.PUBLISHED["hog1_total_uM"], abs=5e-5)
    assert baselines["fps1_total_uM"] == pytest.approx(cw.PUBLISHED["fps1_total_uM"], abs=5e-4)


def test_the_published_protocol_reproduces_the_four_minute_slt2_peak(report):
    """Table S7's own S1, S2, toff, tm; only the dilution time is swept.

    The article's structure, all four claims at once: a real Slt2 peak at 4 min *below* the
    pre-stress volume, which is the observation the sensitizer was introduced for; large peaks at
    14/30/45 min *above* it; nothing at 45 s and 90 s; and Hog1 saturating on the hyper-shock.
    """
    rows = report["published_protocol"]
    assert rows[0.75]["slt2p_peak_percent"] == pytest.approx(25.13, abs=0.1)
    assert rows[1.5]["slt2p_peak_percent"] == pytest.approx(25.82, abs=0.1)
    assert rows[4.0]["slt2p_peak_percent"] == pytest.approx(38.95, abs=0.2)
    assert rows[14.0]["slt2p_peak_percent"] == pytest.approx(89.15, abs=0.3)
    assert rows[30.0]["slt2p_peak_percent"] == pytest.approx(84.72, abs=0.3)
    assert rows[45.0]["slt2p_peak_percent"] == pytest.approx(84.14, abs=0.3)
    for minutes in (0.75, 1.5, 4.0):
        assert rows[minutes]["volume_max_after_dilution"] < 1.0, minutes
    for minutes in (14.0, 30.0, 45.0):
        assert rows[minutes]["volume_max_after_dilution"] > 1.0, minutes
    assert rows[4.0]["slt2p_peak_percent"] / 24.6 > 1.5
    assert rows[1.5]["slt2p_peak_percent"] / 24.6 < 1.1
    for row in rows.values():
        assert row["hog1pp_max_percent"] > 95.0


def test_the_fourteen_minute_run_is_table_s7s_own_default(model):
    """``Shock.published()`` takes ts from Table S7, which prints 840 s = 14 min."""
    shock = cw.Shock.published()
    assert shock.t_s_s == 840.0
    assert shock.s1_uM == 800000.0 and shock.s2_uM == 270000.0
    assert shock.t_off_s == 120.0
    assert cw.Shock.deposited().s1_uM == 0.0
    assert cw.Shock.deposited().t_s_s == 800000.0


def test_the_shock_ramp_and_dilution_match_table_s7(model):
    """``cen(t)``: zero before toff, a tm-limited ramp to S1, then a decay to S2."""
    shock = cw.Shock.published()
    assert shock.added_osmolarity_uM(0.0) == 0.0
    assert shock.added_osmolarity_uM(119.9) == 0.0
    assert shock.added_osmolarity_uM(120.0) == 0.0
    assert shock.added_osmolarity_uM(120.0 + 10.0) == pytest.approx(800000.0 * (1 - np.exp(-1)))
    assert shock.added_osmolarity_uM(120.0 + 500.0) == pytest.approx(800000.0, rel=1e-15)
    assert shock.added_osmolarity_uM(120.0 + 840.0 + 500.0) == pytest.approx(270000.0, rel=1e-15)


# --------------------------------------------------------------------------------------
# The sign, which is the thing most likely to be got wrong later.
# --------------------------------------------------------------------------------------

def test_slt2_answers_a_volume_increase_and_hog1_answers_a_decrease(report):
    """The counterintuitive direction, as arithmetic on the rest points.

    This test exists to fail if anyone ever re-reads this axis as a cell-wall-damage sensor. A
    damage sensor is not monotone increasing in cell volume.
    """
    sign = report["sign"]
    assert sign["slt2_monotone_increasing"]
    assert sign["hog_monotone_decreasing"]
    rows = sign["rows"]
    assert rows[0.90]["slt2p_relative_percent"] == pytest.approx(10.70, abs=0.1)
    assert rows[1.00]["slt2p_relative_percent"] == pytest.approx(24.60, abs=1e-6)
    assert rows[1.10]["slt2p_relative_percent"] == pytest.approx(124.73, abs=0.5)
    assert rows[0.90]["hog_signal_uM"] > 20 * rows[1.10]["hog_signal_uM"]


def test_the_slt2_driver_is_proportional_to_volume_and_nothing_else(model):
    """``v9 = k9 * Vos``: the whole input side of this axis, in one flux."""
    y0 = model.initial_vector()
    constants = model.constants()
    base = model.fluxes(0.0, y0, cw.Shock.deposited(), constants)["v9"]
    swollen = y0.copy()
    swollen[cw.STATE_NAMES.index("Vos")] *= 1.2
    assert model.fluxes(0.0, swollen, cw.Shock.deposited(), constants)["v9"] == pytest.approx(
        1.2 * base, rel=1e-14)
    shrunk = y0.copy()
    shrunk[cw.STATE_NAMES.index("Slt2P")] *= 3.0
    assert model.fluxes(0.0, shrunk, cw.Shock.deposited(), constants)["v9"] == pytest.approx(base)


def test_the_relative_readout_is_not_bounded_by_one_hundred(report):
    """A large hypo-shock drives the Fit rule past 100%, so it is not a saturating fraction."""
    assert report["sign"]["rows"][1.10]["slt2p_relative_percent"] > 100.0


# --------------------------------------------------------------------------------------
# What the switches do, and what they do not.
# --------------------------------------------------------------------------------------

def test_the_deposit_selects_no_crosstalk_in_either_direction():
    """SW_HI = SW_SI = 0 is the article's model-selection result, not a simplification here."""
    assert cw.DEPOSIT_SWITCHES["SW_HI"] == 0.0
    assert cw.DEPOSIT_SWITCHES["SW_SI"] == 0.0
    assert cw.DEPOSIT_SWITCHES["SW_Sensitizer"] == 1.0
    default = cw.OsmoStat()
    assert default.sensitizer and not default.hog1_inhibits_slt2 and not default.slt2_inhibits_hog1


def test_the_sensitizer_switch_is_not_the_paper_ranking():
    """Flipping SW_Sensitizer moves the baseline 48000-fold, so it is not a like-for-like model.

    Pinned because the flat peaks it produces look like the paper's conclusion and are not: the
    unsensitized signal saturates 660000-fold above km11 and cannot respond to anything.
    """
    comparison = cw.sensitizer_switch_comparison()
    assert comparison["baseline_ratio"] > 1e4
    assert comparison["km11_multiples_at_baseline"]["no_sensitizer"] > 1e5
    assert comparison["km11_multiples_at_baseline"]["sensitized"] < 100
    peaks = comparison["peak_percent_after_dilution"]
    for minutes in (4.0, 14.0, 30.0):
        assert peaks[minutes]["no_sensitizer"] < 24.6, minutes
        assert peaks[minutes]["sensitized"] > peaks[minutes]["no_sensitizer"], minutes
    assert "not the paper" in comparison["this_is_not_the_paper_ranking"].lower()


def test_v16_is_a_floor_and_its_own_constant_is_nearly_inert(report):
    """The reaction is load-bearing; the number Table S8 and the deposit disagree about is not."""
    rows = report["vmax16_disagreement"]
    assert rows[14.0]["peak_difference_points"] < 1.0
    assert rows[30.0]["peak_difference_points"] < 1e-2
    for pair in rows.values():
        for reading in ("deposit_3_26314", "table_s8_one"):
            assert pair[reading]["sensitizer_min"] > float(cw.KI16)
    deleted = cw.OsmoStat(vmax16=0.0)
    shock = cw.Shock.published()
    hours = (shock.t_off_s + shock.t_s_s + 1800.0) / cw.SECONDS_PER_HOUR
    run = deleted.timecourse(hours=hours, shock=shock)
    window = run["time_s"] >= shock.t_off_s + shock.t_s_s
    assert float(run["states"]["Sensitizer"].min()) == pytest.approx(2.486, abs=0.02)
    assert float(run["slt2p_relative_percent"][window].max()) == pytest.approx(115.7, abs=1.0)


def test_the_published_vmax16_disagrees_with_the_deposit_and_both_stay_runnable():
    """Table S8 says 1 and the Results repeat it; the deposit says 3.26314."""
    assert float(cw.VMAX16) == 3.26314
    assert cw.TABLE_S8_VMAX16 == 1.0
    assert cw.PUBLISHED["vmax16_table_s8"] == 1.0
    assert "vmax16_published_as_one" in cw.SOURCE_DEFECTS
    text = (ROOT / cw.ARTIFACTS["tables_s4_s8_text"][0]).read_text(errors="replace")
    assert "vmax16" in text and "3.08897" in text and "999.898" in text
    assert "3.26314" not in text, "Table S8 must not contain the deposited vmax16"


# --------------------------------------------------------------------------------------
# The deposit, its identity, and the audit.
# --------------------------------------------------------------------------------------

def test_the_pinned_artifacts_still_hash_to_what_was_transcribed():
    for key, (relative, expected) in cw.ARTIFACTS.items():
        path = ROOT / relative
        assert path.exists(), relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, key


def test_the_two_sbml_files_agree_in_value_and_differ_only_in_formatting(report):
    """The deposit's identity, established by comparison because the printed accession 404s."""
    identity = report["audit"]["identity"]
    assert identity["parameter_values_differing"] == 0
    assert len(identity["parameter_strings_reformatted"]) == 13
    assert "MODEL1604100004" in identity["statement"]
    assert "404" in identity["statement"]


def test_the_deposit_ships_the_protocol_switched_off(report):
    """s1 = s2 = 0 and ts = 800000, which is why the initial state is a rest point."""
    assert cw.Shock.deposited().s1_uM == 0.0
    assert "800000" in cw.DEPOSIT_SHIPS_UNSTRESSED
    assert float(cw.T_S) == 840.0


def test_the_audit_covers_every_registered_literal(report):
    audit = report["audit"]
    assert audit["parameters_checked"] >= 37
    assert sorted(audit["rate_rules"]) == ["HOGSignal", "Slt2Signal", "Vos"]
    assert len(audit["reactions"]) == 12
    assert set(audit["species"]) == {
        "Glyin", "Osmin", "Hog1", "Hog1PP", "Slt2", "Slt2P", "Glyex", "Osmex",
        "Fps1", "Fps1P", "Sensitizer"}
    for _amount, _boundary, substance_units in audit["species"].values():
        assert substance_units, "every species is an amount, not a concentration"


def test_the_deposits_three_rate_rules_are_the_odes_not_held_constants():
    """Unlike williamson2009's seven literal zeros, none of these can be audited away."""
    assert "NONE of these is a literal zero" in cw.SBML_LOADER_REFUSAL
    text = (ROOT / cw.ARTIFACTS["sbml_biomodels"][0]).read_text()
    assert text.count("<rateRule") == 3


@pytest.fixture
def tampered(tmp_path, monkeypatch):
    """Stage a modified deposit and re-pin its hash, so the audit reaches the rules."""
    def stage(old: str, new: str):
        relative, _ = cw.ARTIFACTS["sbml_biomodels"]
        source = (ROOT / relative).read_text()
        assert source.count(old) == 1, f"{old!r} occurs {source.count(old)} times"
        staged = tmp_path / "data" / "native_reference_models" / "talemi2016"
        staged.mkdir(parents=True, exist_ok=True)
        for key, (rel, _digest) in cw.ARTIFACTS.items():
            target = tmp_path / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((ROOT / rel).read_bytes())
        body = source.replace(old, new)
        (tmp_path / relative).write_text(body)
        digest = hashlib.sha256(body.encode()).hexdigest()
        patched = dict(cw.ARTIFACTS)
        patched["sbml_biomodels"] = (relative, digest)
        monkeypatch.setattr(cw, "ARTIFACTS", patched)
        monkeypatch.setattr(cw.paths, "data_dir", lambda: tmp_path / "data")
    return stage


def test_a_changed_constant_is_caught_by_the_audit(tampered):
    tampered('id="k9" metaid="COPASI55" name="k9" value="0.0772144"',
             'id="k9" metaid="COPASI55" name="k9" value="0.0772145"')
    with pytest.raises(ValueError, match="k9"):
        cw.audit_deposit()


def test_a_changed_initial_amount_is_caught_by_the_audit(tampered):
    tampered('initialAmount="4.30809697774826"', 'initialAmount="4.30809697774827"')
    with pytest.raises(ValueError, match="species block"):
        cw.audit_deposit()


def test_a_changed_switch_is_caught_by_the_audit(tampered):
    tampered('id="SW_HI" metaid="COPASI73" name="SW_HI" value="0"',
             'id="SW_HI" metaid="COPASI73" name="SW_HI" value="1"')
    with pytest.raises(ValueError, match="SW_HI"):
        cw.audit_deposit()


def test_a_changed_compartment_size_is_caught_by_the_audit(tampered):
    tampered('id="Membrane" metaid="COPASI4" name="Membrane" size="1"',
             'id="Membrane" metaid="COPASI4" name="Membrane" size="2"')
    with pytest.raises(ValueError, match="compartment sizes"):
        cw.audit_deposit()


def test_a_changed_deposit_trips_the_checksum_first(tmp_path, monkeypatch):
    relative, _ = cw.ARTIFACTS["sbml_biomodels"]
    for key, (rel, _digest) in cw.ARTIFACTS.items():
        target = tmp_path / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT / rel).read_bytes())
    body = (ROOT / relative).read_text().replace('value="0.0772144"', 'value="0.0772145"')
    (tmp_path / relative).write_text(body)
    monkeypatch.setattr(cw.paths, "data_dir", lambda: tmp_path / "data")
    with pytest.raises(ValueError, match="sha256"):
        cw.audit_deposit()


# --------------------------------------------------------------------------------------
# What the source does not contain, refused rather than defaulted.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "promoter_response.rlm1_box",
    "cell_wall_damage_dose",
    "medium.cell_wall_damaging_agent",
    "calcofluor_activating_module",
    "osmotically_active_cell_volume_bridge",
    "external_osmolarity_basal",
])
def test_the_refusals_raise_rather_than_defaulting(name):
    param = cw.NOT_DRIVEN[name]
    assert param.tag == Tag.REFUSED
    assert param.reason and param.missing
    with pytest.raises(RefusedValue):
        float(param)


@pytest.mark.parametrize("absent", [
    "Rlm1", "RLM1", "congo", "Congo", "caffeine", "Calcofluor", "calcofluor",
    "caspofungin", "zymolyase", "Pkc1", "Bck1", "Wsc1", "Mid2", "transcription",
])
def test_the_strings_the_refusals_rest_on_are_genuinely_absent_from_the_deposit(absent):
    """The refusals claim the deposit has no Rlm1 and no damage input. Checked, not asserted."""
    text = (ROOT / cw.ARTIFACTS["sbml_biomodels"][0]).read_text()
    assert absent not in text


def test_the_calcofluor_module_is_in_the_paper_and_not_in_the_deposit():
    """The source's own cell-wall-damage input exists, is out of the deposit, and is refused."""
    module_text = (ROOT / cw.ARTIFACTS["calcofluor_module_text"][0]).read_text(errors="replace")
    assert "Calcofluor mediated Slt2 activating module" in module_text
    assert "CALSignal" in module_text and "Degrader" in module_text
    deposit = (ROOT / cw.ARTIFACTS["sbml_biomodels"][0]).read_text()
    for species in ("Calcofluor", "CALSignal", "Degrader"):
        assert species not in deposit
    tables = (ROOT / cw.ARTIFACTS["tables_s4_s8_text"][0]).read_text(errors="replace")
    for constant in ("k17", "k18", "v17", "v18"):
        assert constant not in tables
    with pytest.raises(RefusedValue):
        float(cw.NOT_DRIVEN["calcofluor_activating_module"])


def test_the_sign_statement_is_quoted_from_the_article_and_the_article_is_pinned():
    """`SIGN_OF_THE_AXIS` quotes the Results; the quote has to be in the vendored article."""
    article = (ROOT / cw.ARTIFACTS["article_xml"][0]).read_text(errors="replace")
    assert "activation upon volume decrease below or increase above the initial volume" in article
    assert "volume decrease below or increase above the initial volume" in cw.SIGN_OF_THE_AXIS
    assert "MUST NOT be described, wired or dosed as a cell-wall-damage sensor" in cw.SIGN_OF_THE_AXIS


def test_the_axis_is_not_called_cell_wall():
    """The ``carbon_camp_input`` precedent, enforced."""
    assert cw.AXIS == "cell_wall_slt2"
    assert cw.PANEL_CHANNEL["module"] == "cell_wall"
    assert cw.PANEL_CHANNEL["answered"] == {}
    assert set(cw.PANEL_CHANNEL["forfeited"]) == {"congo_red", "caffeine"}


def test_the_state_variables_are_amounts_and_declare_their_compartments():
    assert len(cw.STATE_VARIABLES) == 12
    assert cw.STATE_NAMES[:3] == ("Vos", "HOGSignal", "Slt2Signal")
    units = {row[0]: row[1] for row in cw.STATE_VARIABLES}
    assert units["Vos"] == "fL" and units["Slt2Signal"] == "uM"
    assert all(units[name] == "umol" for name in
               ("Hog1", "Hog1PP", "Slt2", "Slt2P", "Glyin", "Glyex", "Fps1", "Fps1P",
                "Sensitizer"))
    source = (ROOT / "src/ystwin/mech/cell_wall_talemi2016.py").read_text()
    assert "The nine species are AMOUNTS, not concentrations" in source


# --------------------------------------------------------------------------------------
# Provenance, grading and the gate.
# --------------------------------------------------------------------------------------

def test_every_registered_parameter_is_graded_and_none_claims_to_be_measured():
    """Not one constant here is MEASURED, and the three abundances say why in their source."""
    for param in cw.TALEMI_PARAMS.params:
        assert param.tag in Tag.ALL
        assert param.tag not in (Tag.MEASURED, Tag.DERIVED, Tag.BORROWED), param.name
        assert param.source.strip()
    assert "yeastgfp" in cw.TALEMI_PARAMS["hog1_molecules"].source
    assert "NOT in Table S7" in cw.TALEMI_PARAMS["slt2_molecules"].source


def test_the_gate_refuses_and_every_target_is_fitted():
    result = cw.gate()
    assert not result.passes
    assert len(result.free) == len(cw.TALEMI_PARAMS.params)
    assert result.targets == ()
    assert set(result.fitted_targets) == {"slt2p_baseline", "hog1pp_baseline", "slt2_peak_4min"}
    for target in cw.TALEMI_PARAMS.targets:
        assert target.fitted


def test_the_calculated_rows_are_not_registered_as_parameters():
    """Registering them would count five degrees of freedom that do not exist."""
    for name in cw.CALCULATED_ROWS:
        assert name not in cw.TALEMI_PARAMS


def test_this_transcription_can_never_claim_upstream_provenance(report):
    assert cw.AVAILABILITY == "transcribed_here"
    assert report["availability"] == "transcribed_here"
    assert "local_verified" not in report["availability"]


def test_the_licence_records_both_licences():
    assert "CC0" in cw.LICENCE and "CC BY 4.0" in cw.LICENCE
    assert "public domain" in cw.LICENCE


def test_time_enters_in_hours_and_is_converted_once(model):
    """The public entry point takes hours; Table S7's protocol row keeps the source's seconds."""
    assert cw.SECONDS_PER_HOUR == 3600.0
    run = model.timecourse(hours=0.5, points=11)
    assert run["time_s"][-1] == pytest.approx(1800.0)
    assert run["time_h"][-1] == pytest.approx(0.5)
    with pytest.raises(ValueError):
        model.timecourse(hours=0.0)


def test_the_provenance_table_renders_both_registries():
    table = cw.provenance()
    assert "rlm1_box" in table
    assert "slt2_molecules" in table
    assert "REFUSED" in table and "ASSERTED" in table


def test_the_scope_limit_travels_with_the_report(report):
    note = report["scope_limit"]
    assert "answers NEITHER stressor" in note
    assert "NO promoter output" in note
    assert set(report["refusals"]) == {param.name for param in cw.NOT_DRIVEN.refusals()}
