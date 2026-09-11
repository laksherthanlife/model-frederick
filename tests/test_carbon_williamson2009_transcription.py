"""Does the Williamson 2009 transcription reproduce Williamson 2009?

THE TESTS THAT EARN THE FILE ITS PLACE are
:func:`test_thirty_six_published_figure_landmarks_reproduce` and
:func:`test_the_deposited_initial_state_is_a_steady_state`. A transcription that merely runs is
worth nothing. This one is scored against thirty-six values read off Williamson's own printed
Figures 8 and 9 -- cAMP, the three cAMP fluxes, the four Gpa2-module species, both PKA pools,
Ras2a and both Cdc25 forms -- and against the deposit's own initial vector, which its authors
solved to a steady state and which a faithful transcription therefore has to leave at rest.

WHAT REPRODUCES: all 36 landmarks, to a mean 0.30% and a worst 1.61% of the panel's own y-axis
full scale, and the residual is 3.5e-14 /s relative. Nothing in
`src/ystwin/mech/carbon_williamson2009.py` was fitted; the tolerances below were written after
the numbers were computed and are tight enough that a later edit which changes the transcription
breaks them.

THE NEGATIVE CONTROLS matter as much as the reproductions, because a check that cannot fail is
not a check. :func:`test_the_pka_stoichiometry_is_load_bearing` shows the deposit's
``stoichiometry="2"`` moves the second cAMP peak by 1.7x AND that the steady-state residual is
blind to it, which is the limit of that check stated as a test rather than as a hope.
:func:`test_the_table_4_transport_reading_is_measurably_worse` shows the arbitration of
condition (b) is decisive rather than a preference. Four tamper tests stage a modified copy of
the deposit, re-pin its hash so the audit actually reaches the rules, and drive each branch of
the condition (a) adapter: a nonzero rateRule, a rateRule on a dynamic species, a deleted
rateRule and a changed stoichiometric coefficient all raise. A fifth leaves the hash pinned and
confirms that, without the re-pin, nothing gets past the checksum at all.

The rest protect what makes those meaningful: that the pinned SBML still hashes to what was
transcribed, that the shared loader still refuses the file for the reason claimed, that what the
source does not contain refuses instead of defaulting, and that this artifact can never be
promoted to an availability level claiming upstream provenance for our own re-implementation.
"""
from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from ystwin.mech import carbon_williamson2009 as wm
from ystwin.mech.kinetic_sbml import KineticModel, UnsupportedSBMLError
from ystwin.mech.params import RefusedValue, Tag

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def report():
    return wm.reproduction_report()


@pytest.fixture(scope="module")
def model():
    return wm.WilliamsonCAMP()


# --------------------------------------------------------------------------------------
# The reproduction.
# --------------------------------------------------------------------------------------

def test_thirty_six_published_figure_landmarks_reproduce(report):
    """Williamson Figures 8 and 9, all five panels, under the deposited parameterisation.

    Scored in percent of each panel's own y-axis full scale, because that is the precision a
    printed figure can be read at; 2% of a rendered 300-pixel panel is about six pixels. Every
    one of the 36 lands inside it, and the aggregate is an order of magnitude better than that.
    """
    deposit = report["figure_reproduction"]["deposit"]
    assert deposit["landmarks"] == 36
    assert deposit["landmarks_within_tolerance"] == 36
    assert deposit["mean_percent_of_full_scale"] == pytest.approx(0.296, abs=0.02)
    assert deposit["max_percent_of_full_scale"] == pytest.approx(1.61, abs=0.05)
    assert deposit["worst_landmark"] == "pde2_t240"
    for key, row in deposit["rows"].items():
        assert row["error_percent_of_full_scale"] <= 100.0 * wm.FIGURE_READ_TOLERANCE, key


def test_the_headline_camp_transient_matches_figure_8(report):
    """Figure 8's four landmark cAMP values, named individually so a regression says which.

    This is the behaviour the engine's present cAMP form structurally cannot produce: a peak
    that comes back down, then a larger peak that also comes back down.
    """
    rows = report["figure_reproduction"]["deposit"]["rows"]
    assert rows["cAMP_basal"]["here"] == pytest.approx(8.259e-4, rel=1e-3)
    assert rows["cAMP_peak1"]["here"] == pytest.approx(2.491e-3, rel=2e-3)
    assert rows["cAMP_t240"]["here"] == pytest.approx(1.332e-3, rel=2e-3)
    assert rows["cAMP_peak2"]["here"] == pytest.approx(3.223e-3, rel=2e-3)
    # The transient is the point: both peaks decay rather than running away.
    assert rows["cAMP_t240"]["here"] < 0.6 * rows["cAMP_peak1"]["here"]
    assert rows["cAMP_t420"]["here"] < 0.7 * rows["cAMP_peak2"]["here"]


def test_the_deposited_initial_state_is_a_steady_state(model):
    """The check that needs no figure, and exercises every rate law at once.

    Williamson's Methods describe the procedure -- find a steady state, then set all
    concentrations to it -- and the deposit's initial vector is the result. A transcription with
    one mistyped constant or one wrong coefficient would not land seventeen coupled states at a
    machine-precision residual.
    """
    residual = model.steady_state_residual()
    assert residual["max_relative_per_s"] < 1e-12
    assert residual["max_absolute_per_s"] < 1e-15


def test_the_stoichiometry_conserves_the_seven_moieties(report):
    """Seven totals the deposit's stoichiometry conserves, over the full 600 s protocol."""
    initial = report["conservation"]["initial"]
    assert initial["gpa2_total"] == pytest.approx(1.126e-4, rel=1e-9)
    assert initial["pka_total"] == pytest.approx(8.675e-5, rel=1e-9)
    for name, drift in report["conservation"]["max_drift_over_600_s"].items():
        assert drift < 1e-15, name


# --------------------------------------------------------------------------------------
# Negative controls. A check that cannot fail is not a check.
# --------------------------------------------------------------------------------------

def test_the_pka_stoichiometry_is_load_bearing(model, report):
    """``PKAact`` releases TWO catalytic subunits, and only the dynamics notice.

    Read 1:1 the steady-state residual is unchanged, because at a fixed point the two fluxes
    balance either way -- so this test also pins the limit of the residual check, which the
    module's own note states. The second cAMP peak is where the error shows: 5.4e-3 mM against
    a published 3.18e-3, and it never comes back down.
    """
    assert dict(wm.STOICHIOMETRY["PKAact"]) == {"PKAi": -1, "C": 2}
    assert dict(wm.STOICHIOMETRY["PKAdeact"]) == {"C": -2, "PKAi": 1}

    class OneToOne(wm.WilliamsonCAMP):
        def rhs(self, t, y):
            rates = self.reaction_rates(y)
            index = {name: i for i, name in enumerate(wm.DYNAMIC_SPECIES)}
            derivative = np.zeros(len(wm.DYNAMIC_SPECIES))
            for reaction, coefficients in wm.STOICHIOMETRY.items():
                scale = 0.5 if reaction in ("PKAact", "PKAdeact") else 1.0
                for species, coefficient in coefficients.items():
                    step = coefficient if species != "C" else coefficient * scale
                    derivative[index[species]] += step * rates[reaction]
            return derivative

    wrong = OneToOne()
    assert wrong.steady_state_residual()["max_relative_per_s"] < 1e-12
    run = wrong.simulate(duration_h=420.0 / wm.SECONDS_PER_HOUR)
    after = run["time_s"] > 240.0
    peak = float(run["species"]["cAMP"][after].max())
    published = 3.18e-3
    assert peak > 1.5 * published
    assert peak > 1.5 * report["figure_reproduction"]["deposit"]["rows"]["cAMP_peak2"]["here"]


def test_the_table_4_transport_reading_is_measurably_worse(report):
    """Condition (b), arbitrated against the source's own figures and not by preference.

    The cited source cannot decide it -- reference [37] is a MATLAB package -- so the deposit
    reading is chosen on the only evidence that requires inventing nothing: it reproduces
    Williamson's own printed output and the Table 4 reading does not, by a factor of five in
    mean error and four in the worst landmark.
    """
    arbitration = report["transport_arbitration"]
    assert arbitration["cited_source_can_arbitrate"] is False
    assert arbitration["chosen"] == "deposit"
    assert arbitration["deposit_mean_percent_of_full_scale"] < 0.5
    assert arbitration["table_4_mean_percent_of_full_scale"] > 1.0
    assert (arbitration["table_4_mean_percent_of_full_scale"]
            > 4.0 * arbitration["deposit_mean_percent_of_full_scale"])
    table_4 = report["figure_reproduction"]["table_4"]
    assert table_4["max_percent_of_full_scale"] > 100.0 * wm.FIGURE_READ_TOLERANCE
    assert table_4["landmarks_within_tolerance"] < 36
    # Not a transposition: the file's Km is the paper's V, but its V is not the paper's Km.
    assert wm.DEPOSIT_TRANSPORT[1] == wm.TABLE_4_TRANSPORT[0]
    assert wm.DEPOSIT_TRANSPORT[0] != wm.TABLE_4_TRANSPORT[1]


def test_clipping_the_three_negative_initials_changes_nothing_measurable(report):
    """Condition (c): the cost of the clip is measured, not assumed."""
    clip = report["negative_initial_clip"]
    assert set(clip["clipped"]) == {"Gpr1Glucout", "Gpa2a", "Gpa2aKrh"}
    assert all(value < 0.0 for value in clip["clipped"].values())
    assert all(abs(value) < 1e-30 for value in clip["clipped"].values())
    assert clip["max_relative_camp_difference"] < 1e-8
    assert set(clip["left_as_deposited"]) == {"Glucout", "Glucin"}
    assert all(value > 0.0 for value in clip["left_as_deposited"].values())


# --------------------------------------------------------------------------------------
# Condition (a): the seven rateRules are audited, not stripped.
# --------------------------------------------------------------------------------------

def test_the_shared_loader_still_refuses_this_file():
    """The premise of the whole hand transcription, checked rather than quoted."""
    relative, _ = wm.ARTIFACTS["sbml_additional_file_1"]
    with pytest.raises(UnsupportedSBMLError, match="rate and algebraic rules are unsupported"):
        KineticModel.from_sbml(str(ROOT / relative))


def test_the_seven_rate_rules_are_audited_and_kept(report):
    """All seven are literal zeros on species already flagged ``boundaryCondition``."""
    rules = report["audit"]["rate_rules"]
    assert rules["count"] == 7
    assert rules["all_literal_zero"] is True
    assert rules["all_on_boundary_species"] is True
    assert rules["stripped"] is False
    assert tuple(sorted(rules["targets"])) == tuple(sorted(wm.BOUNDARY_SPECIES))
    assert report["audit"]["assignment_rules"] == ["percentActiveGpa2", "percentFreeC"]


def _staged_deposit(tmp_path, monkeypatch, tampered: str):
    """Stage a modified copy of the deposit and re-pin its hash, so the audit reaches the rules.

    Re-pinning is the point: without it every tamper trips the checksum first and the rateRule
    branches are never reached, which is a test that looks green and checks nothing.
    """
    relative = wm.ARTIFACTS["sbml_additional_file_1"][0]
    staged = tmp_path / "data" / "native_reference_models" / "williamson2009"
    staged.mkdir(parents=True)
    target = staged / Path(relative).name
    target.write_text(tampered)
    (staged / "PMC2719611.xml").write_bytes(
        (ROOT / wm.ARTIFACTS["article_xml"][0]).read_bytes())
    monkeypatch.setattr(wm.paths, "data_dir", lambda: tmp_path / "data")
    monkeypatch.setitem(
        wm.ARTIFACTS, "sbml_additional_file_1",
        (relative, hashlib.sha256(target.read_bytes()).hexdigest()))


def _deposit_text() -> str:
    return (ROOT / wm.ARTIFACTS["sbml_additional_file_1"][0]).read_text()


def test_a_rate_rule_that_is_not_a_literal_zero_is_refused(tmp_path, monkeypatch):
    """Condition (a), branch one. A nonzero rule is a rule, and it raises rather than dropping."""
    text = _deposit_text()
    tampered = text.replace(
        '<rateRule variable="Rgs2">\n        <math xmlns="http://www.w3.org/1998/Math/MathML">\n'
        '          <cn type="integer"> 0 </cn>',
        '<rateRule variable="Rgs2">\n        <math xmlns="http://www.w3.org/1998/Math/MathML">\n'
        '          <cn type="integer"> 1 </cn>', 1)
    assert tampered != text
    _staged_deposit(tmp_path, monkeypatch, tampered)
    with pytest.raises(ValueError, match="is not the literal 0"):
        wm.audit_deposit()


def test_a_rate_rule_on_a_non_boundary_species_is_refused(tmp_path, monkeypatch):
    """Condition (a), branch two. Zero on a DYNAMIC species would silently freeze a state."""
    text = _deposit_text()
    tampered = text.replace('<rateRule variable="Rgs2">', '<rateRule variable="cAMP">', 1)
    assert tampered != text
    _staged_deposit(tmp_path, monkeypatch, tampered)
    with pytest.raises(ValueError, match="not a boundaryCondition"):
        wm.audit_deposit()


def test_a_deleted_rate_rule_is_refused(tmp_path, monkeypatch):
    """Condition (a), branch three: stripping a rule is exactly what must not go unnoticed."""
    text = _deposit_text()
    tampered = text.replace(
        '<rateRule variable="Rgs2">\n        <math xmlns="http://www.w3.org/1998/Math/MathML">\n'
        '          <cn type="integer"> 0 </cn>\n        </math>\n      </rateRule>\n', "", 1)
    assert tampered != text
    _staged_deposit(tmp_path, monkeypatch, tampered)
    with pytest.raises(ValueError, match="by rateRule"):
        wm.audit_deposit()


def test_a_changed_stoichiometry_in_the_deposit_is_refused(tmp_path, monkeypatch):
    """The PKA coefficient the dynamics depend on cannot drift away from the file unnoticed."""
    text = _deposit_text()
    tampered = text.replace('<speciesReference species="C" stoichiometry="2"/>',
                            '<speciesReference species="C" stoichiometry="3"/>', 1)
    assert tampered != text
    _staged_deposit(tmp_path, monkeypatch, tampered)
    with pytest.raises(ValueError, match="stoichiometry"):
        wm.audit_deposit()


def test_a_changed_deposit_trips_the_checksum_first(tmp_path, monkeypatch):
    """And with the hash left pinned, no tamper gets as far as the rules at all."""
    relative = wm.ARTIFACTS["sbml_additional_file_1"][0]
    staged = tmp_path / "data" / "native_reference_models" / "williamson2009"
    staged.mkdir(parents=True)
    (staged / Path(relative).name).write_text(_deposit_text().replace("0.062", "0.08", 1))
    monkeypatch.setattr(wm.paths, "data_dir", lambda: tmp_path / "data")
    with pytest.raises(ValueError, match="sha256"):
        wm.audit_deposit()


def test_the_audit_covers_every_constant_and_coefficient(report):
    """Nothing typed into the module is unchecked against the deposit."""
    audit = report["audit"]
    assert audit["species_checked"] == 24
    assert audit["parameters_checked"] == 39
    assert audit["reactions_checked"] == 20
    assert audit["unused_species"] == ["Sdc25"]


def test_the_pinned_sbml_still_hashes_to_what_was_transcribed():
    for label, (relative, expected) in wm.ARTIFACTS.items():
        path = ROOT / relative
        assert path.exists(), label
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, label


def test_a_changed_constant_is_caught_by_the_audit(monkeypatch):
    """The audit is the reason no literal in the module can silently drift from the deposit."""
    original = wm._sbml_parameter_literals

    def drifted():
        values = original()
        values["cAMPhydroPde1p_kcat_13"] = 25.25
        return values

    monkeypatch.setattr(wm, "_sbml_parameter_literals", drifted)
    with pytest.raises(ValueError, match="cAMPhydroPde1p_kcat_13"):
        wm.audit_deposit()


def test_the_deposit_carries_seven_rate_rules_and_two_assignment_rules():
    """Read straight out of the XML, so the audit's own claim has an independent witness."""
    relative, _ = wm.ARTIFACTS["sbml_additional_file_1"]
    namespace = "{http://www.sbml.org/sbml/level2}"
    model = ET.parse(ROOT / relative).getroot().find(f"{namespace}model")
    rules = list(model.find(f"{namespace}listOfRules"))
    kinds = [rule.tag.replace(namespace, "") for rule in rules]
    assert kinds.count("rateRule") == 7
    assert kinds.count("assignmentRule") == 2


# --------------------------------------------------------------------------------------
# What the source does not contain refuses, instead of defaulting.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "williamson.snf1_adr1_cat8_csre_gain",
    "williamson.glucose_transport_cited_measurement",
    "williamson.absolute_abundance_scale",
    "williamson.cell_volume_litres",
    "williamson.pka_to_growth_or_burden",
])
def test_the_refusals_raise_rather_than_defaulting(name):
    param = wm.NOT_IN_SOURCE[name]
    assert param.tag == Tag.REFUSED
    assert param.value is None
    assert param.reason.strip() and param.missing.strip()
    with pytest.raises(RefusedValue):
        float(param)


def test_the_promoter_arm_is_refused_by_name():
    """The scope limit, enforced rather than only written down.

    This module closes the cAMP/PKA input side. Snf1 -> Adr1/Cat8 -> CSRE is not in Williamson
    and is not built here, and calling ``float()`` on it says so with the citation attached.
    """
    with pytest.raises(RefusedValue, match="Snf1"):
        float(wm.PROMOTER_COUPLING_NOT_BUILT)


def test_the_cell_volume_helper_refuses():
    with pytest.raises(RefusedValue, match="compartment"):
        wm.cell_volume_litres()


def test_the_transport_reading_cannot_be_a_third_number():
    """Two readings are published. A third would be a fitted parameter."""
    with pytest.raises(ValueError, match="names where it is printed"):
        wm.TransportReading("split_the_difference", 0.88, 0.89, 0.91)


# --------------------------------------------------------------------------------------
# What this artifact is allowed to claim.
# --------------------------------------------------------------------------------------

def test_the_gate_refuses_and_the_one_target_is_fitted():
    """Sixty-three free scalars against zero independent targets, reported not hidden.

    The only measurement near this model is the cAMP time course its 35 estimated parameters
    were fitted to, so it cannot score them. Figure 9 is the fitted model's own output.
    """
    gate = wm.gate()
    assert not gate.passes
    assert len(gate.free) == 63
    assert gate.targets == ()
    assert gate.fitted_targets == ("rolland2000_camp_timecourse",)
    assert "REFUSED" in gate.summary()


def test_every_registered_parameter_is_graded_and_none_claims_to_be_measured():
    """Not one constant of this model is a measurement, and none is graded as one."""
    tags = {param.tag for param in wm.WILLIAMSON_PARAMS.params}
    assert tags == {Tag.ASSERTED}
    assert {param.tag for param in wm.NOT_IN_SOURCE.params} == {Tag.REFUSED}
    assert len(wm.WILLIAMSON_PARAMS) == 63


def test_this_transcription_can_never_claim_upstream_provenance(report):
    """``local_verified`` means a publisher file the shared loader runs unchanged. This is not."""
    assert wm.AVAILABILITY == "transcribed_here"
    assert report["availability"] == "transcribed_here"
    assert report["biological_validation"] is False
    assert report["source_parameters_refitted"] is False
    assert "CC BY 2.0" in wm.LICENCE


def test_the_scope_limit_travels_with_the_report(report):
    scope = report["scope"]
    assert "cAMP/PKA input side" in scope["closes"]
    assert "CSRE" in scope["does_not_close"]


def test_the_pka_coupling_is_dimensionless_and_matches_the_deposits_own_rule(model):
    """``percentFreeC`` is the deposit's assignmentRule; this is it, as a fraction.

    Dimensionless is deliberate: it needs no cell volume, which Williamson does not supply.
    """
    basal = model.pka_active_fraction(model.initial_state())
    assert basal == pytest.approx(0.1180, abs=5e-4)
    run = model.simulate(duration_h=600.0 / wm.SECONDS_PER_HOUR)
    stacked = np.array([run["species"][name] for name in wm.DYNAMIC_SPECIES])
    peak = max(model.pka_active_fraction(stacked[:, column])
               for column in range(0, stacked.shape[1], 20))
    assert 0.0 < basal < peak < 1.0
    assert peak > 5.0 * basal


def test_time_enters_in_hours_and_is_converted_once(model):
    """Condition (d). The source runs in seconds; every entry point here takes hours."""
    assert wm.NATIVE_TIME_UNIT == "second"
    assert wm.SECONDS_PER_HOUR == 3600.0
    run = model.simulate(duration_h=600.0 / wm.SECONDS_PER_HOUR)
    assert run["time_s"][-1] == pytest.approx(600.0)
    assert run["time_h"][-1] == pytest.approx(600.0 / 3600.0)
    with pytest.raises(ValueError, match="pulse times"):
        model.simulate(duration_h=30.0 / wm.SECONDS_PER_HOUR)


def test_the_provenance_table_renders_both_registries():
    table = wm.provenance()
    assert "REFUSED" in table
    assert "williamson.glucose_transport_Km" in table
    assert "williamson.snf1_adr1_cat8_csre_gain" in table
