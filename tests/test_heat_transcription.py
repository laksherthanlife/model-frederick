"""The Zheng 2016 / Krakowiak 2018 transcription, scored against the source's own figures.

A transcription earns its place by reproducing the paper it came from, so most of this file
re-derives published numbers. Several assertions deliberately encode a FAILURE, because the
whole point of running a reproduction is that a mismatch is a result:

* Krakowiak's transcription gain fails Zheng's own published parameter-screen criterion 3.
* Under that gain the 43 C reporter curve has not plateaued when Zheng's published window ends.
* Krakowiak's "increased affinity gives faster deactivation" does not hold at 3-fold.
* The free-scalar gate REFUSES this piece, and the test asserts the refusal rather than
  routing around it.

Every number asserted here is also written down in ``data/transcriptions/zheng2016.json``, and
the last test in the file checks the record against the live computation, so a drift in the port
falls out as a failure rather than surviving in a JSON file nobody re-runs.
"""
from __future__ import annotations

import hashlib
import json
import math
import pathlib
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ystwin import paths
from ystwin.mech import heat_zheng2016 as heat
from ystwin.mech.params import RefusedValue, SweptValue, Tag

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
GAINS = ("zheng2016", "krakowiak2018")


# --- what this artifact is allowed to claim ---------------------------------------------


def test_the_transcription_is_never_recorded_as_a_verified_upstream_model():
    record = heat.transcription_record()
    assert record["availability"] == "local_transcription"
    assert record["availability"] != "local_verified"
    assert record["is_upstream_artifact"] is False
    assert record["produced_by"] == "src/ystwin/mech/heat_zheng2016.py"
    # The upstream files keep their own, lower level; the port does not inherit it.
    assert record["upstream_availability"] == "local_source_only"


def test_the_record_refuses_to_load_if_it_ever_claims_the_verified_level(tmp_path, monkeypatch):
    record = json.loads((paths.data_dir() / "transcriptions" / "zheng2016.json").read_text())
    record["availability"] = "local_verified"
    fake = tmp_path / "transcriptions"
    fake.mkdir()
    (fake / "zheng2016.json").write_text(json.dumps(record))
    monkeypatch.setattr(paths, "data_dir", lambda: tmp_path)
    with pytest.raises(ValueError, match="never be availability=local_verified"):
        heat.transcription_record()


def test_every_transcription_source_is_checksummed_and_present():
    for relative, expected in heat.SOURCE_DIGESTS.items():
        path = REPO_ROOT / relative
        assert path.exists(), relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, relative


@pytest.mark.parametrize("source_id", GAINS)
def test_transcription_identifiers_match_the_pinned_primary_article(source_id):
    source = next(row for row in heat.transcription_record()["sources"]
                  if row["source_id"] == source_id)
    article = ET.parse(REPO_ROOT / source["local_xml"])
    identifiers = {node.get("pub-id-type"): node.text
                   for node in article.findall("./front/article-meta/article-id")}
    for kind in ("pmid", "pmcid", "doi"):
        assert source[kind].lower() == identifiers[kind].lower()
    assert any(f"PMID {identifiers['pmid']}" in param.source for param in heat.HEAT_PARAMS.params)


def test_the_two_papers_are_one_source():
    """Krakowiak Source code 1 IS Zheng's titration_YFP_FB.m, byte for byte."""
    zheng, krakowiak = heat.BYTE_IDENTICAL_PAIR
    assert (REPO_ROOT / zheng).read_bytes() == (REPO_ROOT / krakowiak).read_bytes()
    assert heat.SOURCE_DIGESTS[zheng] == heat.SOURCE_DIGESTS[krakowiak]
    assert heat.SOURCE_DIGESTS[zheng].startswith("5998094f1af311d73ecbff972f48d754")
    record = heat.transcription_record()
    assert "not be recorded as independent corroboration" in record["one_source_not_two"]["consequence"]


def test_a_changed_source_file_stops_the_transcription(tmp_path, monkeypatch):
    monkeypatch.setitem(heat.SOURCE_DIGESTS,
                        "data/native_reference_models/zheng2016/Fig1_IPMS.m", "0" * 64)
    with pytest.raises(ValueError, match="the transcription is stale"):
        heat.HeatShockTranscription(gain=heat.PUBLISHED_GAINS["zheng2016"])


# --- the two refusals -------------------------------------------------------------------


def test_temperature_entry_stays_refused():
    param = heat.TEMPERATURE_TO_UNFOLDED_PROTEIN
    assert param.tag == Tag.REFUSED
    assert param.value is None
    with pytest.raises(RefusedValue, match="temperature_to_unfolded_protein"):
        float(param)
    # No entry point may take a temperature: the load is an input, always.
    for name in ("screen_criteria", "reporter_fold_change", "feedback_separation",
                 "basal_versus_dissociation", "affinity_timecourse"):
        signature = getattr(heat, name).__doc__ or ""
        assert "celsius" not in signature.lower().replace(" c", "")
    assert heat.ABSOLUTE_ABUNDANCE_SCALE.tag == Tag.REFUSED


def test_the_gain_has_no_default_and_no_midpoint():
    assert heat.TRANSCRIPTION_GAIN.tag == Tag.SWEPT
    with pytest.raises(SweptValue, match="FACTOR OF FIVE"):
        float(heat.TRANSCRIPTION_GAIN)
    assert heat.TRANSCRIPTION_GAIN.bounds == pytest.approx((1.7783 / 5.0, 1.7783))
    assert heat.gain("zheng2016") == pytest.approx(1.7783)
    assert heat.gain("krakowiak2018") == pytest.approx(1.7783 / 5.0)
    assert heat.gain("zheng2016") / heat.gain("krakowiak2018") == pytest.approx(5.0)
    with pytest.raises(KeyError, match="name the paper"):
        heat.gain("1.7783")


@pytest.mark.parametrize("entry", ["screen_criteria", "reporter_fold_change",
                                   "feedback_separation", "basal_versus_dissociation"])
def test_no_entry_point_defaults_the_gain(entry):
    with pytest.raises(TypeError):
        getattr(heat, entry)()


# --- Zheng's three published parameter-screen criteria ----------------------------------


@pytest.fixture(scope="module")
def screen():
    return {source: heat.screen_criteria(gain_source=source) for source in GAINS}


def test_criterion_1_hsf1_is_bound_at_basal(screen):
    """'Hsf1 must be bound to Hsp70 in basal conditions at 25 C.'"""
    for source in GAINS:
        assert screen[source].criterion_1_bound_at_basal
        assert screen[source].bound_fraction_basal == pytest.approx(0.98359, rel=1e-3)


def test_criterion_2_the_complex_dissociates_within_five_minutes(screen):
    """'...must dissociate to <=10% of its basal level within five minutes.' Both gains pass."""
    assert screen["zheng2016"].ratio_at_5_min == pytest.approx(0.04729, rel=1e-3)
    assert screen["krakowiak2018"].ratio_at_5_min == pytest.approx(0.00918, rel=1e-3)
    for source in GAINS:
        assert screen[source].criterion_2_dissociates


def test_criterion_3_reassociation_holds_for_zheng(screen):
    """'...must re-associate to >=90% of its initial level within 60 min.'"""
    zheng = screen["zheng2016"]
    assert zheng.criterion_3_reassociates
    assert zheng.recovery_time_min == pytest.approx(41.71, rel=2e-3)
    assert zheng.ratio_at_60_min == pytest.approx(0.91719, rel=1e-3)
    assert zheng.passes


def test_criterion_3_FAILS_under_the_krakowiak_gain(screen):
    """THE SHARPEST RESULT. A five-fold lower gain on the byte-identical model file fails the
    acceptance criterion that selected its own co-parameters, and no paper says so."""
    krakowiak = screen["krakowiak2018"]
    assert not krakowiak.criterion_3_reassociates
    assert krakowiak.recovery_time_min is None
    assert krakowiak.ratio_at_60_min == pytest.approx(0.15385, rel=1e-3)
    assert 0.90 / krakowiak.ratio_at_60_min == pytest.approx(5.85, rel=1e-2)
    assert not krakowiak.passes
    assert "FAILS" in krakowiak.summary()


def test_the_source_attaches_one_load_to_two_temperatures():
    """The load 10 a.u. is 35 C in one shipped driver and 39 C in another, while a third gives
    35 C the load 4.4492. Read directly out of the files rather than asserted in prose."""
    root = REPO_ROOT / "data/native_reference_models/zheng2016"
    ipms = root.joinpath("Fig1_IPMS.m").read_text(errors="replace").splitlines()
    phospho = root.joinpath("Hsf1_Phosph.m").read_text(errors="replace").splitlines()
    reporter = root.joinpath("Fig2_YFP_Reporter.m").read_text(errors="replace").splitlines()
    heading = next(i for i, ln in enumerate(ipms) if "INITIAL CONDITIONS" in ln and "35" in ln)
    load = next(ln for ln in ipms[heading:] if ln.strip().startswith("UPo"))
    assert load.split("=")[1].strip().startswith("10;")
    assert any(ln.strip().startswith("UPo_39") and "10;" in ln for ln in phospho)
    assert any(ln.strip().startswith("UPo_35 = 4.4492") for ln in reporter)
    assert any(ln.strip().startswith("UPo_39 = 10.5141") for ln in reporter)


def test_the_screen_load_and_the_table_load_agree_on_the_criteria():
    """Whichever of the two the 39 C run is meant to use, the three verdicts are the same."""
    for source in GAINS:
        table = heat.screen_criteria(gain_source=source)
        screened = heat.screen_criteria(gain_source=source,
                                        unfolded_protein=heat.SCREEN_UNFOLDED_LOAD_39C)
        assert abs(table.ratio_at_5_min - screened.ratio_at_5_min) < 0.006
        assert table.criterion_2_dissociates == screened.criterion_2_dissociates
        assert table.criterion_3_reassociates == screened.criterion_3_reassociates
    assert heat.SCREEN_UNFOLDED_LOAD_39C != heat.PUBLISHED_UNFOLDED_LOADS[39]


# --- Zheng Figure 1D: temperature-dependent plateaus ------------------------------------


@pytest.fixture(scope="module")
def figure_1d():
    return {source: heat.reporter_fold_change(gain_source=source) for source in GAINS}


def test_figure_1d_is_ordered_by_load_and_fits_its_published_axis(figure_1d):
    for source in GAINS:
        folds = [figure_1d[source][label]["log2_fold_at_end"] for label in (25, 35, 39, 43)]
        assert folds == sorted(folds)
        assert all(figure_1d[source][label]["within_published_axis"] for label in (25, 35, 39, 43))
    expected = {35: 1.1490, 39: 2.0817, 43: 2.8870}
    for label, value in expected.items():
        assert figure_1d["zheng2016"][label]["log2_fold_at_end"] == pytest.approx(value, rel=1e-3)


def test_the_plateaus_hold_under_the_zheng_gain(figure_1d):
    """Zheng's own gloss: a plateau means 'no more YFP is being produced in the simulation'."""
    for label in (35, 39, 43):
        panel = figure_1d["zheng2016"][label]
        rate = np.gradient(panel["reporter"], figure_1d["zheng2016"]["time_min"])
        assert rate[-1] / rate.max() < 0.005, label


def test_the_43C_plateau_FAILS_under_the_krakowiak_gain(figure_1d):
    """A second, independent contradiction of the reduced gain: the height barely moves, so the
    failure is invisible in the plotted curve and lives entirely in its slope."""
    panel = figure_1d["krakowiak2018"][43]
    times = figure_1d["krakowiak2018"]["time_min"]
    rate = np.gradient(panel["reporter"], times)
    assert rate[-1] / rate.max() == pytest.approx(0.7408, rel=1e-2)
    assert panel["terminal_rate_per_min"] == pytest.approx(0.1122, rel=1e-2)
    zheng = figure_1d["zheng2016"][43]
    assert panel["terminal_rate_per_min"] / zheng["terminal_rate_per_min"] > 100
    assert abs(panel["log2_fold_at_end"] - zheng["log2_fold_at_end"]) < 0.15
    # The lower temperatures still plateau, so this is specific to the largest load.
    for label in (35, 39):
        other = figure_1d["krakowiak2018"][label]
        other_rate = np.gradient(other["reporter"], times)
        assert other_rate[-1] / other_rate.max() < 0.01, label


# --- Krakowiak Figure 1B / Figure 1-figure supplement 1 ---------------------------------


def test_severing_the_feedback_is_one_deleted_term():
    """Source code 2 differs from Source code 1 only in d[HSP]/dt, so one class runs both."""
    zheng = REPO_ROOT / "data/native_reference_models/krakowiak2018/elife-31668-code1-v3.m"
    severed = REPO_ROOT / "data/native_reference_models/krakowiak2018/elife-31668-code2-v3.m"
    intact_lines = [ln for ln in zheng.read_text().splitlines() if ln.strip()]
    severed_lines = [ln for ln in severed.read_text().splitlines() if ln.strip()]
    differing = [(a, b) for a, b in zip(intact_lines, severed_lines) if a != b]
    assert len(differing) == 2  # the function name, and d[HSP]/dt
    assert "beta*Hsf1^nH" in differing[1][0] and "beta*Hsf1^nH" not in differing[1][1]

    gain_value = heat.gain("zheng2016")
    times = np.linspace(0.0, 5.0, 51)
    intact = heat.HeatShockTranscription(gain=gain_value)
    cut = heat.HeatShockTranscription(gain=gain_value, feedback=False)
    start = intact.initial_state(heat.PUBLISHED_UNFOLDED_LOADS[39], reporter=3.0)
    # Only the chaperone pool differs; the reporter equation is untouched in the source.
    assert not np.allclose(intact.simulate(times, start)[:, 0], cut.simulate(times, start)[:, 0])


def test_the_reduced_gain_delays_the_wild_type_versus_severed_separation():
    """Krakowiak's stated reason for the change: the first gain 'exaggerated' the feedback and
    separated the two curves too early. That is the only published justification for it."""
    fast = heat.feedback_separation(gain_source="zheng2016")
    slow = heat.feedback_separation(gain_source="krakowiak2018")
    assert fast["separation_min"] == {0.10: 13.0, 0.25: 18.0, 0.50: 22.0}
    assert slow["separation_min"] == {0.10: 70.0, 0.25: 84.0, 0.50: 105.0}
    for threshold in (0.10, 0.25, 0.50):
        assert slow["separation_min"][threshold] > fast["separation_min"][threshold]
    assert slow["separation_min"][0.25] / fast["separation_min"][0.25] == pytest.approx(4.67, rel=1e-2)
    assert fast["severed_over_wild_type_at_end"] == pytest.approx(14.52, rel=1e-3)
    assert slow["severed_over_wild_type_at_end"] == pytest.approx(3.121, rel=1e-3)


# --- Krakowiak Source code 3 / Figure 2-figure supplement 1 -----------------------------


def test_the_basal_sweep_identifies_which_gain_made_the_published_figure():
    """The published y-axis stops at 12. Only the reduced gain fits under it."""
    zheng = heat.basal_versus_dissociation(gain_source="zheng2016")
    krakowiak = heat.basal_versus_dissociation(gain_source="krakowiak2018")
    assert zheng["monotone_increasing"] and krakowiak["monotone_increasing"]
    assert float(krakowiak["fold_change_at_horizon"].max()) == pytest.approx(9.2427, rel=1e-3)
    assert float(zheng["fold_change_at_horizon"].max()) == pytest.approx(22.6695, rel=1e-3)
    assert krakowiak["within_published_axis"]
    assert not zheng["within_published_axis"]
    assert float(zheng["fold_change_at_horizon"].max()) / heat.PUBLISHED_AXIS_MAXIMUM > 1.8


def test_the_published_basal_axis_is_set_by_the_horizon_not_by_the_model():
    """kdil = 0 makes the reporter a pure integrator, so Source code 3's 'steady-state level'
    does not exist. Doubling its 240-minute horizon takes the curve off its own axis."""
    assert float(heat.REPORTER_DILUTION) == 0.0
    folds = {}
    for horizon in (120.0, 240.0, 480.0, 960.0):
        sweep = heat.basal_versus_dissociation(gain_source="krakowiak2018", horizon_min=horizon)
        folds[horizon] = float(sweep["fold_change_at_horizon"][-1])
    assert folds[120.0] == pytest.approx(5.8414, rel=1e-3)
    assert folds[240.0] == pytest.approx(9.2427, rel=1e-3)
    assert folds[480.0] == pytest.approx(14.0637, rel=1e-3)
    assert folds[960.0] == pytest.approx(20.3394, rel=1e-3)
    assert folds[240.0] <= heat.PUBLISHED_AXIS_MAXIMUM < folds[480.0]
    assert sorted(folds.values()) == list(folds.values())


# --- Krakowiak Source code 4 / Figure 6A ------------------------------------------------


@pytest.fixture(scope="module")
def affinity():
    return {source: {name: heat.affinity_timecourse(gain_source=source, multiplier=multiplier)
                     for name, multiplier in heat.AFFINITY_MULTIPLIERS.items()}
            for source in GAINS}


def test_lower_affinity_raises_the_maximal_output_under_both_gains(affinity):
    """The half of the published claim that reproduces, and it reproduces monotonically."""
    order = ("k2_over_3", "wild_type", "k2_times_5", "k2_times_50")
    for source in GAINS:
        outputs = [affinity[source][name]["reporter_at_end"] for name in order]
        assert outputs == sorted(outputs), source
    assert affinity["zheng2016"]["wild_type"]["reporter_at_end"] == pytest.approx(12.7623, rel=1e-3)
    assert affinity["zheng2016"]["k2_times_50"]["reporter_at_end"] == pytest.approx(20.7189, rel=1e-3)
    assert affinity["krakowiak2018"]["wild_type"]["reporter_at_end"] == pytest.approx(12.6355, rel=1e-3)
    assert affinity["krakowiak2018"]["k2_times_50"]["reporter_at_end"] == pytest.approx(17.0594, rel=1e-3)


def test_the_shipped_fifty_fold_case_deactivates_more_slowly(affinity):
    """Figure 6A itself -- the only affinity panel either paper ships code for."""
    for source in GAINS:
        wild = affinity[source]["wild_type"]
        mutant = affinity[source]["k2_times_50"]
        assert mutant["deactivation_min"] > wild["deactivation_min"], source
        assert mutant["basal_reporter"] > wild["basal_reporter"], source


def test_deactivation_time_is_NOT_monotone_in_affinity():
    """The half of the published claim that does NOT reproduce. 'Increased affinity showed faster
    deactivation kinetics' reverses at 3-fold, and 'decreasing the affinity showed slower
    deactivation' reverses at 5-fold, under every shut-off threshold tested."""
    for fraction in (0.01, 0.02, 0.05, 0.10):
        times = {name: heat.affinity_timecourse(gain_source="zheng2016", multiplier=multiplier,
                                                shutoff_fraction=fraction)["deactivation_min"]
                 for name, multiplier in heat.AFFINITY_MULTIPLIERS.items()}
        assert times["k2_over_3"] > times["wild_type"], fraction
        assert times["k2_times_50"] > times["wild_type"], fraction
        ordered = [times[name] for name in ("k2_over_3", "wild_type", "k2_times_5", "k2_times_50")]
        assert ordered != sorted(ordered, reverse=True), fraction
    at_five_percent = {name: heat.affinity_timecourse(gain_source="zheng2016",
                                                      multiplier=multiplier)["deactivation_min"]
                       for name, multiplier in heat.AFFINITY_MULTIPLIERS.items()}
    assert at_five_percent["k2_over_3"] == pytest.approx(24.82, rel=1e-2)
    assert at_five_percent["wild_type"] == pytest.approx(17.94, rel=1e-2)
    assert at_five_percent["k2_times_5"] == pytest.approx(15.60, rel=1e-2)


# --- the source's own equations disagree with the source's own code ---------------------


def test_the_shipped_code_loses_client_at_twice_the_published_rate():
    """titration_YFP_FB.m line 30 subtracts k5*[HSP.UP] from FREE unfolded protein; Zheng's
    published d[UP]/dt does not. Structural, and numerically almost invisible."""
    source_text = (REPO_ROOT / "data/native_reference_models/zheng2016/titration_YFP_FB.m").read_text()
    up_line = next(ln for ln in source_text.splitlines() if ln.strip().startswith("ty(3)"))
    assert "- k5*HSP_UP" in up_line

    state = np.array([0.4, 1e-3, 5.0, 1e-3, 6.0, 3.0])
    gain_value = heat.gain("zheng2016")
    shipped = heat.HeatShockTranscription(gain=gain_value).rhs(0.0, state)
    published = heat.HeatShockTranscription(gain=gain_value,
                                            client_loss="published_equations").rhs(0.0, state)
    k5 = float(heat.CLIENT_REFOLDING_RATE)
    total_shipped = shipped[2] + shipped[4]
    total_published = published[2] + published[4]
    assert total_shipped == pytest.approx(-2.0 * k5 * state[4], rel=1e-6)
    assert total_published == pytest.approx(-1.0 * k5 * state[4], rel=1e-6)
    # Every other state is untouched by the choice.
    assert np.allclose(np.delete(shipped, 2), np.delete(published, 2))


def test_the_client_equation_choice_does_not_change_any_reported_answer():
    times = np.linspace(0.0, 240.0, 2401)
    for source in GAINS:
        gain_value = heat.gain(source)
        shipped = heat.HeatShockTranscription(gain=gain_value)
        published = heat.HeatShockTranscription(gain=gain_value,
                                                client_loss="published_equations")
        start = shipped.initial_state(heat.PUBLISHED_UNFOLDED_LOADS[39], reporter=3.0)
        a = shipped.simulate(times, start)[:, 5]
        b = published.simulate(times, start)[:, 5]
        assert np.max(np.abs(a - b) / b) < 5e-5, source


def test_an_unknown_client_loss_convention_is_refused():
    with pytest.raises(ValueError, match="the source's code and the source's paper"):
        heat.HeatShockTranscription(gain=heat.gain("zheng2016"), client_loss="mine")


# --- invariants the port must not break -------------------------------------------------


def test_total_hsf1_is_conserved():
    """d[Hsf1]/dt + d[HSP.Hsf1]/dt is identically zero in the source, so total Hsf1 is a
    conserved quantity and its drift measures the port's integration error, nothing else."""
    times = np.linspace(0.0, 240.0, 2401)
    for source in GAINS:
        model = heat.HeatShockTranscription(gain=heat.gain(source))
        start = model.initial_state(heat.PUBLISHED_UNFOLDED_LOADS[43], reporter=3.0)
        trace = model.simulate(times, start)
        total = trace[:, 1] + trace[:, 3]
        assert np.max(np.abs(total - total[0]) / total[0]) < 1e-6, source


def test_the_reporter_never_decreases():
    """kdil = 0 with a non-negative Hill term: the reporter is an integrator by construction."""
    times = np.linspace(0.0, 240.0, 2401)
    model = heat.HeatShockTranscription(gain=heat.gain("krakowiak2018"))
    trace = model.simulate(times, model.initial_state(heat.PUBLISHED_UNFOLDED_LOADS[39],
                                                      reporter=3.0))
    assert np.all(np.diff(trace[:, 5]) >= -1e-9)


def test_the_reported_criteria_survive_the_matlab_default_tolerance():
    """MATLAB ode23s's defaults are rtol 1e-3 / atol 1e-6, and the published figures were made
    with them. Nothing asserted here is an artefact of this port's much tighter settings."""
    times = np.linspace(0.0, 60.0, 6001)
    for source, expected, passes in (("zheng2016", 0.91719, True),
                                     ("krakowiak2018", 0.15385, False)):
        model = heat.HeatShockTranscription(gain=heat.gain(source))
        loose_control = model.simulate(times, model.initial_state(0.0, reporter=0.0),
                                       rtol=1e-3, atol=1e-6)
        loose_shock = model.simulate(times,
                                     model.initial_state(heat.PUBLISHED_UNFOLDED_LOADS[39],
                                                         reporter=0.0), rtol=1e-3, atol=1e-6)
        ratio = loose_shock[-1, 3] / loose_control[-1, 3]
        # Agreement to 0.3% at MATLAB's own defaults, and the verdict is nowhere near the edge.
        assert ratio == pytest.approx(expected, rel=5e-3), source
        assert bool(ratio >= 0.9) is passes, source


def test_an_unfolded_protein_load_must_be_a_finite_non_negative_input():
    model = heat.HeatShockTranscription(gain=heat.gain("zheng2016"))
    for bad in (-1.0, math.inf, math.nan, True):
        with pytest.raises(ValueError, match="finite non-negative load"):
            model.initial_state(bad, reporter=3.0)


# --- criterion (e), and what it is allowed to conclude ----------------------------------


def test_the_free_scalar_gate_refuses_this_piece():
    """The correct verdict, asserted rather than routed around: a port of a fitted model scored
    on the fitted model's own figures has no independent target and cannot be licensed by them."""
    result = heat.gate()
    assert not result.passes
    assert len(result.free) == 15
    assert result.targets == ()
    assert len(result.fitted_targets) == 4
    assert "REFUSED" in result.summary()
    assert not heat.HEAT_PARAMS.independent_targets()


def test_every_registered_target_is_marked_fitted():
    for target in heat.HEAT_PARAMS.targets:
        assert target.fitted, target.name
        assert target.source.strip()


def test_no_registered_parameter_is_graded_measured_or_derived():
    """Not one number in this model is measured in S. cerevisiae. Fifteen screened or fixed
    values and two refusals; a MEASURED grade appearing here would be a mistake, not a win."""
    tags = {param.tag for param in heat.HEAT_PARAMS.params}
    assert tags <= {Tag.ASSERTED, Tag.SWEPT, Tag.REFUSED}
    assert Tag.MEASURED not in tags and Tag.DERIVED not in tags
    assert len(heat.HEAT_PARAMS.refusals()) == 2


def test_provenance_says_what_this_is_in_one_line():
    line = heat.provenance()
    assert line.startswith("TRANSCRIPTION (local_transcription), not an upstream model.")
    assert "byte-identical" in line
    assert "REFUSED" in line


# --- the record and the code cannot drift apart -----------------------------------------


def test_the_recorded_reproduction_matches_the_live_computation():
    record = heat.transcription_record()
    targets = {entry["id"]: entry for entry in record["reproduction"]["targets"]}

    screen = {source: heat.screen_criteria(gain_source=source) for source in GAINS}
    c3 = targets["zheng_screen_criterion_3"]["measured"]
    assert screen["zheng2016"].ratio_at_60_min == pytest.approx(c3["beta_zheng"]["ratio_at_60_min"],
                                                                rel=1e-3)
    assert c3["beta_zheng"]["passes"] is True
    assert screen["krakowiak2018"].ratio_at_60_min == pytest.approx(
        c3["beta_krakowiak"]["ratio_at_60_min"], rel=1e-3)
    assert c3["beta_krakowiak"]["passes"] is False
    assert c3["beta_krakowiak"]["recovery_time_min"] is None

    c2 = targets["zheng_screen_criterion_2"]["measured"]
    assert screen["zheng2016"].ratio_at_5_min == pytest.approx(c2["ratio_at_5_min_beta_zheng"],
                                                               rel=1e-3)
    assert screen["krakowiak2018"].ratio_at_5_min == pytest.approx(
        c2["ratio_at_5_min_beta_krakowiak"], rel=1e-3)

    separation = targets["krakowiak_figure_1_supplement_1_separation"]["measured"]
    for source, key in (("zheng2016", "beta_zheng"), ("krakowiak2018", "beta_krakowiak")):
        live = heat.feedback_separation(gain_source=source)
        assert live["separation_min"][0.25] == separation["separation_25pct_min"][key]
        assert live["severed_over_wild_type_at_end"] == pytest.approx(
            separation["no_feedback_over_wild_type_at_240_min"][key], rel=1e-3)

    basal = targets["krakowiak_figure_2_supplement_1_basal_versus_dissociation"]["measured"]
    for source, key in (("zheng2016", "beta_zheng"), ("krakowiak2018", "beta_krakowiak")):
        live = heat.basal_versus_dissociation(gain_source=source)
        assert float(live["fold_change_at_horizon"].max()) == pytest.approx(
            basal["fold_change_at_k2_1e4"][key], rel=1e-3)
    assert basal["published_axis_maximum"] == heat.PUBLISHED_AXIS_MAXIMUM

    outputs = targets["krakowiak_figure_6A_affinity"]["measured"]["reporter_at_240_min"]
    for source, key in (("zheng2016", "beta_zheng"), ("krakowiak2018", "beta_krakowiak")):
        for name, multiplier in heat.AFFINITY_MULTIPLIERS.items():
            live = heat.affinity_timecourse(gain_source=source, multiplier=multiplier)
            assert live["reporter_at_end"] == pytest.approx(outputs[key][name], rel=1e-3), name

    gate = record["free_scalar_gate"]
    assert gate["verdict"] == "REFUSED"
    assert gate["free_scalars"] == len(heat.gate().free)
    assert gate["independent_targets"] == 0


def test_the_record_names_every_source_it_was_transcribed_from():
    record = heat.transcription_record()
    recorded = {entry["path"]: entry["sha256"] for entry in record["transcribed_from"]}
    assert recorded == heat.SOURCE_DIGESTS
    for entry in record["transcribed_from"]:
        assert (REPO_ROOT / entry["path"]).stat().st_size == entry["bytes"]
    for entry in record["not_transcribed"]:
        assert (REPO_ROOT / entry["path"]).exists()
        assert entry["path"] not in heat.SOURCE_DIGESTS


def test_the_record_keeps_both_published_gains_and_the_discrepancies():
    record = heat.transcription_record()
    beta = record["parameters"]["transcription_gain_beta"]
    assert beta["zheng2016"]["code_value"] == pytest.approx(heat.PUBLISHED_GAINS["zheng2016"])
    assert beta["krakowiak2018"]["code_value"] == pytest.approx(
        heat.PUBLISHED_GAINS["krakowiak2018"], rel=1e-4)
    assert beta["ratio"] == 5.0
    assert beta["status"] == "declared_choice_required"
    found = {entry["id"] for entry in record["discrepancies_found_while_transcribing"]}
    assert {"up_equation_paper_versus_code", "fig1_ipms_temperature_label",
            "no_steady_state_exists", "table_rounding"} <= found
    refused = {entry["id"] for entry in record["refusals"]}
    assert "heat.temperature_to_unfolded_protein" in refused
    assert all(entry["value"] is None for entry in record["refusals"])
