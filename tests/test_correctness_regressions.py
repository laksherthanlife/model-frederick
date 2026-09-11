from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis.uncertainty import fold_change
from ystwin.generator.context import CultureContext
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.capacity import CapacityUnmeasured
from ystwin.pathway.enzyme_capacity import EnzymeCapacityUnmeasured
from ystwin.pathway.solve import NodeKinetics
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, SetpointUnreachable, predict_product


def _predict(*, genotype=None, environment=None, kinetics=None, **kwargs):
    return predict_product(
        load_pathway("beta_carotene"),
        Genotype(1.0) if genotype is None else genotype,
        Environment(growth_rate_setpoint_per_h=0.18) if environment is None else environment,
        BETA_CAROTENE_FLUX,
        BETA_CAROTENE_KINETICS if kinetics is None else kinetics,
        **kwargs,
    )


@pytest.mark.parametrize("stressor", [None, "DTT", "H2O2"])
@pytest.mark.parametrize("context,setpoint", [
    (CultureContext(carbon_source="ethanol"), 0.18),
    (CultureContext(carbon_source="ethanol"), 0.14),
    (CultureContext(oxygen=0.0), 0.18),
    (CultureContext(temperature_c=46.0), 0.18),
    (CultureContext(growth_phase="stationary"), 0.18),
])
def test_unreachable_growth_is_independent_of_zero_dose_stressor_label(context, setpoint, stressor):
    with pytest.raises(SetpointUnreachable):
        _predict(environment=Environment(context, stressor=stressor,
                                         growth_rate_setpoint_per_h=setpoint))


@pytest.mark.parametrize("stressor", [None, "DTT", "H2O2"])
def test_a_reachable_ethanol_setpoint_still_returns_the_held_rate(stressor):
    result = _predict(environment=Environment(CultureContext(carbon_source="ethanol"),
                                              stressor=stressor,
                                              growth_rate_setpoint_per_h=0.12))
    assert result.growth_rate_per_h == 0.12
    assert result.content_mg_per_gdcw > 0.0


@pytest.mark.parametrize("numerator,control", [
    (-2.0, -1.0), (-2.0, 1.0), (2.0, -1.0),
    (np.nan, 1.0), (2.0, np.nan), (np.inf, 1.0), (2.0, np.inf),
])
def test_a_fold_requires_finite_positive_numerator_and_control(numerator, control):
    readings = pd.DataFrame([
        {"plate": plate, "well": well, "construct": "UPRE1", "dose_mM": dose,
         "activity_late": value}
        for plate in ("p1", "p2", "p3")
        for well, dose, value in (("A1", 0.0, control), ("B1", 1.0, numerator))
    ])
    result = fold_change(readings, "UPRE1", 1.0, n_resamples=40)
    assert not result.estimable
    assert result.n_plates_contributing == 0
    assert np.isnan(result.point)
    assert np.isnan(result.low) and np.isnan(result.high)


def test_a_refused_cassette_cannot_return_the_calibration_strains_prediction():
    with pytest.raises(CapacityUnmeasured):
        _predict(genotype=Genotype(1.0, cassette={"crtYB": 4.0}))


def test_cassette_refusal_uses_the_exception_not_its_wording(monkeypatch):
    def refuse(*args, **kwargs):
        raise CapacityUnmeasured("measurement unavailable")

    monkeypatch.setattr("ystwin.predict.capacity_for", refuse)
    with pytest.raises(CapacityUnmeasured, match="measurement unavailable"):
        _predict(genotype=Genotype(1.0, cassette={"crtYB": 1.0}))


def test_missing_enzyme_measurements_are_not_swallowed():
    with pytest.raises(EnzymeCapacityUnmeasured, match="kcat_per_s"):
        _predict(enzyme_capacity={"enzyme": "crtYB", "kcat_per_s": None,
                                  "enzyme_mmol_per_gdcw": 1e-6})


@pytest.mark.parametrize("value", [np.nan, np.inf, -1.0])
def test_invalid_enzyme_measurements_are_not_swallowed(value):
    with pytest.raises(ValueError, match="kcat_per_s"):
        _predict(enzyme_capacity={"enzyme": "crtYB", "kcat_per_s": value,
                                  "enzyme_mmol_per_gdcw": 1e-6})


@pytest.mark.parametrize("kwargs,required", [
    ({"audit_reaction": "product"}, "audit_model"),
    ({"audit_glucose_uptake": 10.0}, "audit_model"),
    ({"thermo_model": object()}, "thermo"),
    ({"thermo_background_m": {"s_0189": 1e-4}}, "thermo"),
    ({"thermo_volumes": (1.0,)}, "thermo"),
    ({"refuse_thermodynamically_blocked": True}, "thermo"),
    ({"refuse_thermodynamically_blocked": False}, "thermo"),
])
def test_dependent_prediction_inputs_cannot_be_silently_discarded(kwargs, required):
    with pytest.raises(ValueError, match=required):
        _predict(**kwargs)


def test_kinetics_for_an_absent_node_are_rejected():
    with pytest.raises(ValueError, match="not_a_node"):
        _predict(kinetics={**BETA_CAROTENE_KINETICS,
                           "not_a_node": NodeKinetics(rate_constant=2.0)})


@pytest.mark.parametrize("field,value", [
    ("ph_medium", 6.0), ("ph_cytosolic", 6.5), ("glucose_g_per_L", 10.0),
])
def test_unmodelled_context_inputs_are_explicitly_rejected(field, value):
    context = replace(CultureContext(), **{field: value})
    with pytest.raises(ValueError, match=field):
        _predict(environment=Environment(context, growth_rate_setpoint_per_h=0.18))


def test_the_public_predictor_does_not_acknowledge_an_unvalidated_latent_input_for_the_caller():
    from ystwin.bridge.latent_bridge import LatentState

    with pytest.raises(ValueError, match="latent"):
        _predict(latent=LatentState(0.02, "UPRE2", 4.0))


def test_a_structurally_failed_audit_blocks_a_numeric_prediction_even_with_benign_prose(monkeypatch):
    verdict = SimpleNamespace(within_envelope=False, notes=(), regulation="unstressed",
                              headroom=0.5, summary=lambda: "audit completed")
    monkeypatch.setattr("ystwin.fba.audit.audit_predicted_flux", lambda *a, **k: verdict)
    with pytest.raises(ValueError, match="envelope"):
        _predict(audit_model=object(), audit_reaction="product", audit_glucose_uptake=10.0)


def test_a_failed_mechanistic_gate_does_not_return_an_ordinary_prediction():
    from ystwin.mech.chain import SweepPoint

    with pytest.raises(ValueError):
        _predict(mech=SweepPoint.midpoint())


def test_the_default_is_an_explicit_empirical_result_with_input_provenance():
    genotype = Genotype(1.0, label="calibration-relative strain")
    environment = Environment(growth_rate_setpoint_per_h=0.18)
    result = _predict(genotype=genotype, environment=environment)

    assert result.mode == "empirical"
    assert result.genotype is genotype
    assert result.environment is environment
    assert result.content_mg_per_gdcw > 0.0


def test_legacy_refusal_diagnostics_are_not_ordinary_numeric_prediction_objects():
    result = _predict(genotype=Genotype(1.0, cassette={"crtYB": 4.0}), mode="legacy")

    assert result.mode == "legacy"
    assert not result.supported
    assert not hasattr(result, "content_mg_per_gdcw")
    assert not hasattr(result, "rate_mmol_per_gdcw_h")
    assert not hasattr(result, "solution")
    assert result.calculation.mode == "legacy"
    assert result.calculation.content_mg_per_gdcw == _predict().content_mg_per_gdcw
    assert result.summary().startswith("LEGACY DIAGNOSTIC, not a supported prediction")


def test_unknown_prediction_modes_are_rejected():
    with pytest.raises(ValueError, match="mode"):
        _predict(mode="silently_ignore_refusals")


def test_a_gate_failure_is_authoritative_even_without_refusal_words(monkeypatch):
    from ystwin.mech.chain import SweepPoint

    gate = SimpleNamespace(passes=False, free=("a", "b"), targets=("t",),
                           summary=lambda: "assembly assessed")
    monkeypatch.setattr("ystwin.mech.chain.assembled_gate", lambda *a, **k: gate)
    with pytest.raises(ValueError, match="assembly assessed"):
        _predict(mech=SweepPoint.midpoint())


def test_an_unphysical_bill_is_refused_even_if_the_free_scalar_gate_passes(monkeypatch):
    from ystwin.mech.chain import SweepPoint

    monkeypatch.setattr("ystwin.mech.chain.require_assembled_gate", lambda: None)
    result = SimpleNamespace(bill_is_physical=False)
    monkeypatch.setattr("ystwin.predict._mech_layers", lambda *a, **k: result)
    with pytest.raises(ValueError, match="not physical"):
        _predict(mech=SweepPoint.midpoint())


def test_a_blocked_thermodynamic_step_is_not_just_a_note_by_default(monkeypatch):
    from ystwin.pathway.thermo_gate import (
        Feasibility, GateReport, StepEnergy, ThermodynamicallyBlocked,
    )

    step = StepEnergy("lycopene", "CRTI", Feasibility.CANNOT_RUN, 1.0, 2.0)
    report = GateReport("beta_carotene", 1.0, 303.15, (step,))
    reports = tuple(replace(report, cytosolic_volume_ml_per_gdcw=volume)
                    for volume in (1.0, 2.0, 2.7))
    monkeypatch.setattr("ystwin.predict.gate_across_volumes", lambda *a, **k: reports)
    monkeypatch.setattr(GateReport, "summary", lambda _: "chemistry assessed")

    with pytest.raises(ThermodynamicallyBlocked):
        _predict(thermo=object())
    with pytest.raises(ValueError, match="mode='legacy'"):
        _predict(thermo=object(), refuse_thermodynamically_blocked=False)
    legacy = _predict(thermo=object(), mode="legacy", refuse_thermodynamically_blocked=False)
    assert not legacy.supported
    assert legacy.calculation.thermodynamics[0].cannot_run == (step,)


def test_diagnostic_prose_cannot_overrule_a_structured_audit_success(monkeypatch):
    verdict = SimpleNamespace(within_envelope=True, notes=("not REFUSED",),
                              regulation="unstressed", headroom=10.0,
                              summary=lambda: "REFUSED is not this result")
    monkeypatch.setattr("ystwin.fba.audit.audit_predicted_flux", lambda *a, **k: verdict)
    result = _predict(audit_model=object(), audit_reaction="product", audit_glucose_uptake=10.0)

    assert result.flux_audit is verdict
    assert result.content_mg_per_gdcw > 0.0


@pytest.mark.parametrize("cassette", [{"crtE": 10.0}, {"not_a_pathway_gene": 2.0}])
def test_gene_copy_inputs_without_a_calibrated_route_are_rejected(cassette):
    with pytest.raises(ValueError, match="cassette dosage"):
        _predict(genotype=Genotype(1.0, cassette=cassette))


@pytest.mark.parametrize("node,parameters", [
    ("phytoene", NodeKinetics(rate_constant=1.0)),
    ("lycopene", replace(BETA_CAROTENE_KINETICS["lycopene"], degradation_rate_per_h=0.1)),
    ("beta_carotene", NodeKinetics(vmax_per_growth=1.0, km=1.0)),
])
def test_supplied_kinetics_must_match_the_nodes_declared_rate_law_and_fate(node, parameters):
    with pytest.raises(ValueError, match="not used"):
        _predict(kinetics={**BETA_CAROTENE_KINETICS, node: parameters})
