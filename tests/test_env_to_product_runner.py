"""Run the requested comparison afresh, retaining refusals without inventing content."""

from __future__ import annotations

import importlib.util
import pathlib
import warnings

import numpy as np
import pandas as pd
import pytest

from ystwin.generator.context import CultureContext
from ystwin.pathway import calibrations
from ystwin.pathway.solve import ImplausibleContent
from ystwin.predict import Environment, Genotype, LayerState, SetpointUnreachable


REQUESTED_CASES = [
    (carbon, rate)
    for carbon in ("glucose", "ethanol")
    for rate in (0.101, 0.15, 0.2543)
]
PREDICTION_COLUMNS = ["growth_rate_per_h", "content_mg_per_gdw"]


def _load_script():
    filename = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "env_to_product.py"
    spec = importlib.util.spec_from_file_location("env_to_product_runner", filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def script():
    with warnings.catch_warnings():
        return _load_script()


@pytest.fixture
def run_comparison(script, monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path / "outputs"))
    predictor = script.predict_product
    calls = []

    def record_prediction(spec, genotype, environment, flux, kinetics, **kwargs):
        calls.append((spec, genotype, environment, flux, kinetics, kwargs))
        return predictor(spec, genotype, environment, flux, kinetics, **kwargs)

    monkeypatch.setattr(script, "predict_product", record_prediction)

    def run():
        assert script.main() == 0
        frame = pd.read_csv(tmp_path / "outputs" / "env_to_product.csv")
        return frame, capsys.readouterr().out, calls

    return run


def test_real_runner_retains_all_requests_and_phb_after_unreachable_setpoint(run_comparison):
    frame, _, calls = run_comparison()
    phb = frame[frame["product"] == "phb"]
    carotene = frame[frame["product"] == "beta_carotene"]

    assert len(frame) == 17
    assert len(phb) == 11
    assert list(zip(carotene.carbon_source, carotene.mu_per_h)) == REQUESTED_CASES
    assert len(calls) == len(REQUESTED_CASES)
    for (carbon, rate), (spec, genotype, environment, flux, kinetics, options) in zip(
            REQUESTED_CASES, calls):
        assert spec.product == "beta_carotene"
        assert genotype == Genotype(entry_expression=1.0)
        assert environment == Environment(
            context=CultureContext(carbon_source=carbon), growth_rate_setpoint_per_h=rate)
        assert flux is calibrations.BETA_CAROTENE_FLUX
        assert kinetics is calibrations.BETA_CAROTENE_KINETICS
        assert options == {"mode": "empirical"}

    refused = carotene[carotene.refused.notna()]
    answered = carotene[carotene.refused.isna()]
    assert list(zip(refused.carbon_source, refused.mu_per_h)) == [
        ("ethanol", 0.15), ("ethanol", 0.2543)]
    assert len(answered) == 4
    assert refused.refused.eq("SetpointUnreachable").all()
    assert refused.refusal_reason.str.contains("maximum growth rate of 0.140 /h").all()
    assert refused.refusal_reason.str.contains("no held steady state").all()
    assert refused[PREDICTION_COLUMNS].isna().all().all()
    assert refused.environment_layer_state.eq(LayerState.NOT_RUN).all()
    assert refused.environment_sets_the_number.isna().all()
    assert answered.refusal_reason.isna().all()
    assert np.isfinite(answered[PREDICTION_COLUMNS]).all().all()
    assert answered.content_mg_per_gdw.gt(0).all()
    assert answered.growth_rate_per_h.equals(answered.mu_per_h)
    assert answered.environment_layer_state.eq(LayerState.REPORTED).all()
    assert answered.environment_sets_the_number.eq(False).all()
    assert carotene.prediction_mode.eq("empirical").all()
    assert carotene.prediction_scope.eq(
        "fixed-gene comparison, not environment validation").all()
    low, high = calibrations.BETA_CAROTENE_KINETICS["lycopene"].growth_rate_range
    assert answered.growth_rate_per_h.between(low, high).all()
    assert phb[["measured_mg_per_gdw", "predicted_mu_only",
                "predicted_with_environment"]].notna().all().all()


def test_summary_counts_refusals_and_never_infers_invariance_from_one_answer(run_comparison):
    _, output, _ = run_comparison()

    assert "4 answered, 2 refused out of 6 requested cases" in output
    assert "mu = 0.101 /h: 2 answered, 0 refused" in output
    for rate in (0.15, 0.2543):
        line = next(line for line in output.splitlines() if f"mu = {rate:g} /h:" in line)
        assert "1 answered, 1 refused" in line
        assert "matched-carbon comparison unavailable" in line
    assert "SetpointUnreachable" in output
    assert "not environment validation" in output
    assert "Missing refused cases cannot establish carbon-source invariance" in output
    assert "carbon source reaches nothing" not in output
    assert "1 means the carbon source changed nothing" not in output


def test_phb_is_retrospective_within_one_cohort_and_keeps_loo_values(run_comparison, script):
    frame, output, _ = run_comparison()
    phb = frame[frame["product"] == "phb"]
    states = pd.read_csv(script.STATES, sep="\t")

    assert phb.state_id.tolist() == states.state_id.tolist()
    assert np.array_equal(phb.measured_mg_per_gdw, states.phb_mg_per_gdw)
    assert phb.prediction_scope.str.contains("retrospective leave-one-state-out").all()
    assert phb.prediction_scope.str.contains("one Kocharin 2013 source cohort").all()
    assert phb.prediction_scope.str.contains("not independent external validation").all()
    assert "retrospective leave-one-state-out within one Kocharin 2013 source cohort" in output
    assert "not independent external validation" in output
    errors = {}
    for column in ("predicted_mu_only", "predicted_with_environment"):
        errors[column] = np.exp(np.abs(np.log(phb[column] / phb.measured_mg_per_gdw)))
    assert np.median(errors["predicted_mu_only"]) == pytest.approx(1.5052, abs=5e-5)
    assert np.median(errors["predicted_with_environment"]) == pytest.approx(1.3201, abs=5e-5)
    assert errors["predicted_with_environment"].max() > 2


@pytest.mark.parametrize("refusal_type", [SetpointUnreachable, ImplausibleContent])
def test_a_typed_refusal_retains_its_reason_and_does_not_stop_later_cases(
        script, monkeypatch, refusal_type):
    predictor = script.predict_product
    attempted = []
    reason = "scientific refusal retained verbatim"

    def refuse_first(spec, genotype, environment, flux, kinetics, **kwargs):
        attempted.append((environment.context.carbon_source, environment.growth_rate_setpoint_per_h))
        if len(attempted) == 1:
            raise refusal_type(reason)
        return predictor(spec, genotype, environment, flux, kinetics, **kwargs)

    monkeypatch.setattr(script, "predict_product", refuse_first)
    rows = script._carotene_rows()

    assert attempted == REQUESTED_CASES
    assert rows.iloc[0].refused == refusal_type.__name__
    assert rows.iloc[0].refusal_reason == reason
    assert rows.iloc[0][PREDICTION_COLUMNS].isna().all()
    assert pd.isna(rows.iloc[0].environment_sets_the_number)
    assert rows.iloc[1].refused == ""
    assert rows.iloc[1].content_mg_per_gdw > 0


def test_all_refused_still_writes_phb_and_six_missing_prediction_rows(
        script, monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path / "outputs"))

    def refuse_every_case(*args, **kwargs):
        raise SetpointUnreachable("no supported held state")

    monkeypatch.setattr(script, "predict_product", refuse_every_case)
    assert script.main() == 0
    output = capsys.readouterr().out
    frame = pd.read_csv(tmp_path / "outputs" / "env_to_product.csv")
    carotene = frame[frame["product"] == "beta_carotene"]

    assert len(frame) == 17
    assert list(zip(carotene.carbon_source, carotene.mu_per_h)) == REQUESTED_CASES
    assert carotene.refused.eq("SetpointUnreachable").all()
    assert carotene.refusal_reason.eq("no supported held state").all()
    assert carotene[PREDICTION_COLUMNS].isna().all().all()
    assert carotene.environment_layer_state.eq(LayerState.NOT_RUN).all()
    assert carotene.environment_sets_the_number.isna().all()
    assert "0 answered, 6 refused out of 6 requested cases" in output
    assert output.count("matched-carbon comparison unavailable") == 3


@pytest.mark.parametrize("error_type", [RuntimeError, ValueError, TypeError, KeyError])
def test_programming_errors_are_not_mislabeled_as_scientific_refusals(
        script, monkeypatch, error_type):
    error = error_type("unexpected predictor error")

    def broken_prediction(*args, **kwargs):
        raise error

    monkeypatch.setattr(script, "predict_product", broken_prediction)
    with pytest.raises(error_type) as raised:
        script._carotene_rows()
    assert raised.value is error


def test_import_does_not_change_the_warning_policy(script):
    with warnings.catch_warnings():
        before = list(warnings.filters)
        _load_script()
        assert warnings.filters == before


def test_prediction_warnings_are_not_suppressed(script, monkeypatch):
    def warning_prediction(*args, **kwargs):
        warnings.warn("unexpected scientific warning", RuntimeWarning)

    monkeypatch.setattr(script, "predict_product", warning_prediction)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(RuntimeWarning, match="unexpected scientific warning"):
            script._carotene_rows()
