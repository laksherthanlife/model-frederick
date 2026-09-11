from dataclasses import replace

import numpy as np
import pytest

from ystwin import predict
from ystwin.mech import contracts, engine


MEDIUM = {"glucose": 50.0, "nitrogen": 20.0, "oxygen": 0.2, "osmolyte": 250.0}


@pytest.fixture(scope="module")
def physical_case():
    parameters = engine.EngineParameters.prior()
    genotype = contracts.Genotype.prior()
    protocol = contracts.Protocol(
        (0.0, 0.13, 0.3, 0.5),
        (contracts.Control(0.0, oxygen_transfer_per_h=30.0, oxygen_saturation_mM=0.2),),
        events=(contracts.Event("sample", 0.13, withdraw_l=0.01),),
        source="entrypoint regression protocol",
    )
    initial = engine.initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=0.2,
                                medium_mM=MEDIUM)
    return protocol, genotype, initial, parameters


@pytest.mark.parametrize("name", [
    "Control", "Event", "Protocol", "PhysicalState", "ObservationModel", "SimulationResult",
])
def test_predict_reexports_the_shared_contract_classes_without_shadow_schemas(name):
    assert getattr(predict, name) is getattr(contracts, name)


def test_the_physical_genotype_and_initializer_are_the_shared_objects():
    assert predict.ProtocolGenotype is contracts.Genotype
    assert predict.ProtocolGenotype is not predict.Genotype
    assert predict.EngineParameters is engine.EngineParameters
    assert predict.initialize is engine.initialize


def test_protocol_adapter_forwards_every_supplied_object_and_solver_option(monkeypatch, physical_case):
    protocol, genotype, initial, parameters = physical_case
    observation = contracts.ObservationModel.prior(seed=19)
    tolerances = {"biomass_gdw": 1e-12}
    recorded = {}
    sentinel = object()

    def simulate(*args, **kwargs):
        recorded["args"] = args
        recorded["kwargs"] = kwargs
        return sentinel

    monkeypatch.setattr(engine, "simulate", simulate)
    result = predict.predict_protocol(
        protocol, genotype, initial, parameters, observation=observation,
        rtol=2e-7, atol=tolerances, max_step_h=0.025, method="Radau",
    )

    assert result is sentinel
    assert all(actual is expected for actual, expected in zip(
        recorded["args"], physical_case, strict=True))
    assert recorded["kwargs"] == {"observation": observation, "rtol": 2e-7,
                                  "atol": tolerances, "max_step_h": 0.025, "method": "Radau"}
    assert recorded["kwargs"]["atol"] is tolerances


def test_protocol_adapter_matches_engine_truth_observations_and_validity(physical_case):
    observation = contracts.ObservationModel.prior(seed=5)
    expected = engine.simulate(*physical_case, observation=observation, rtol=2e-7)
    result = predict.predict_protocol(*physical_case, observation=observation, rtol=2e-7)

    assert isinstance(result, contracts.SimulationResult)
    assert not isinstance(result, predict.ProductPrediction)
    assert result.validity == expected.validity
    assert not result.validity.biological_validation
    assert result.validity.parameter_basis == "prior-conditional; not experimentally validated"
    assert result.prior_parameters == expected.prior_parameters
    assert result.variables == expected.variables
    assert result.final_state == expected.final_state
    assert result.diagnostics["protocol_source"] == physical_case[0].source
    assert result.diagnostics["events_applied"] == ("sample",)
    assert result.diagnostics["rtol"] == 2e-7
    np.testing.assert_array_equal(result.times_h, expected.times_h)
    for name in expected.truth:
        np.testing.assert_array_equal(result.trace(name), expected.trace(name))
    for channel in expected.observations.values:
        np.testing.assert_array_equal(result.observations.values[channel],
                                      expected.observations.values[channel])
        np.testing.assert_array_equal(result.observations.missing[channel],
                                      expected.observations.missing[channel])


def test_protocol_prediction_does_not_use_the_flat_empirical_ceiling(monkeypatch, physical_case):
    def forbidden(*args, **kwargs):
        raise AssertionError("the physical protocol cannot fall back to the empirical pools law")

    for name in ("predict_product", "solve_pathway", "predict_flux_from_expression"):
        monkeypatch.setattr(predict, name, forbidden)
    reference = predict.predict_protocol(*physical_case)
    protocol, genotype, _, parameters = physical_case
    knockout = genotype.with_copies(crtyb=0)
    initial = engine.initialize(parameters, knockout, volume_l=1.0, biomass_gdw_l=0.2,
                                medium_mM=MEDIUM)
    changed = predict.predict_protocol(protocol, knockout, initial, parameters)

    assert reference.trace("internal.beta_carotene")[-1] > 0.0
    np.testing.assert_array_equal(changed.trace("internal.beta_carotene"), 0.0)
    np.testing.assert_array_equal(changed.trace("gene.crtyb.protein"), 0.0)
    assert changed.trace("internal.ggpp")[-1] > reference.trace("internal.ggpp")[-1]
    for result in (reference, changed):
        for element in ("carbon", "nitrogen", "live_biomass"):
            np.testing.assert_allclose(result.trace(f"balance.{element}_residual"), 0.0, atol=3e-9)


def test_protocol_prediction_requires_explicit_state_genotype_and_parameters(physical_case):
    protocol, genotype, initial, parameters = physical_case
    with pytest.raises(TypeError):
        predict.predict_protocol(protocol, genotype, initial)
    with pytest.raises(ValueError, match="typed Protocol and Genotype"):
        predict.predict_protocol(protocol, predict.Genotype(1.0), initial, parameters)
    with pytest.raises(TypeError, match="kinetics"):
        predict.predict_protocol(*physical_case, kinetics={})


def test_protocol_prediction_propagates_engine_refusals_without_a_numeric_fallback(monkeypatch, physical_case):
    refusal = contracts.ScientificRefusal("unsupported declared intervention")

    def refuse(*args, **kwargs):
        raise refusal

    monkeypatch.setattr(engine, "simulate", refuse)
    with pytest.raises(contracts.ScientificRefusal) as raised:
        predict.predict_protocol(*physical_case)
    assert raised.value is refusal


def test_protocol_prediction_enforces_the_engines_temperature_domain(physical_case):
    protocol, genotype, initial, parameters = physical_case
    unsupported = replace(protocol, controls=(replace(protocol.controls[0], temperature_c=45.0),))
    with pytest.raises(contracts.ScientificRefusal):
        predict.predict_protocol(unsupported, genotype, initial, parameters)


def test_protocol_restart_uses_the_shared_final_state_and_event_receipts():
    parameters = engine.EngineParameters.prior()
    genotype = contracts.Genotype(())
    initial = contracts.PhysicalState(0.0, 1.0, 0.0, {"glucose": 1.0})
    event = contracts.Event("carbon bolus", 0.2, add_mmol={"glucose": 2.0})
    first = predict.predict_protocol(
        contracts.Protocol((0.0, 0.2), (contracts.Control(0.0),), events=(event,)),
        genotype, initial, parameters,
    )
    second = predict.predict_protocol(
        contracts.Protocol((0.2, 0.5), (contracts.Control(0.0),), events=(event,)),
        genotype, first.final_state, parameters,
    )

    np.testing.assert_array_equal(first.trace("external.glucose"), [1.0, 3.0])
    np.testing.assert_array_equal(second.trace("external.glucose"), [3.0, 3.0])
    assert second.final_state.event_receipts == (event,)


def test_generator_observations_reweight_physical_particles_and_resume_through_predict(monkeypatch):
    from ystwin import hybrid
    from ystwin.fba import allocation, dynamic
    from ystwin.generator import culture
    from ystwin.mech.inference import PhysicalParticleFilter

    def forbidden(*args, **kwargs):
        raise AssertionError("physical inference must not enter a comparison or a standalone allocation solve")

    for module, names in (
        (culture, ("simulate_empirical_comparison", "simulate_culture")),
        (hybrid, ("simulate_hybrid_comparison", "simulate_hybrid_case")),
        (dynamic, ("simulate_fixed_volume_comparison", "simulate_batch")),
        (allocation, ("solve_allocation",)),
    ):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)
    parameters = engine.EngineParameters.prior()
    genotype = contracts.Genotype.prior()
    observation = contracts.ObservationModel.prior(seed=31)
    control = contracts.Control(0.0, oxygen_transfer_per_h=30.0,
                                oxygen_saturation_mM=0.2, folding_inhibition=0.2)
    sample = contracts.Event("prefix-sample", 0.01, withdraw_l=0.01)
    prefix = contracts.Protocol((0.0, 0.01, 0.02), (control,), events=(sample,))
    # The learner's finite support is declared independently of the evaluator's
    # hidden state. Only OD/RFU observations cross the inference interface.
    ensemble = tuple(engine.initialize(parameters, genotype, volume_l=1.0,
                                       biomass_gdw_l=density, medium_mM=MEDIUM)
                     for density in (0.1, 0.3))
    generated = culture.simulate_protocol(
        prefix, genotype, engine.initialize(parameters, genotype, volume_l=1.0,
                                            biomass_gdw_l=0.25, medium_mM=MEDIUM),
        parameters, observation=observation,
    )
    readings = {name: float(values[-1]) for name, values in generated.observations.values.items()}
    learner = PhysicalParticleFilter(
        ensemble, parameters=parameters, genotype=genotype, observation_model=observation,
        measurement_noise={"od": 0.02, "rfu": 0.5}, weights=(0.5, 0.5), resample_threshold=0.0,
    )
    posterior = learner.step(prefix, readings)
    assert posterior.n_channels_used == 2
    assert posterior.weights[1] > posterior.weights[0]
    assert posterior.weights != (0.5, 0.5)
    assert generated.observations.metadata["truth_feedback"] is False
    assert generated.observations.metadata["units"] == {"od": "OD", "rfu": "RFU"}

    future = contracts.Protocol((0.02, 0.03, 0.04), (control,))
    forecast = learner.forecast(future)
    assert learner.time_h == 0.02  # Forecasts must not assimilate or advance the filter.
    assert forecast.weights == posterior.weights
    for particle, expected in zip(posterior.particles, forecast.particles, strict=True):
        result = predict.predict_protocol(
            future, genotype, particle.state, particle.parameters,
            observation=replace(observation, initial_mature_mmol=particle.mature_mmol),
        )
        assert result.final_state == expected.state
        assert result.final_state.event_receipts == (sample,)
        assert particle.parameters is parameters
        for name in result.truth:
            np.testing.assert_array_equal(result.trace(name), expected.simulation.trace(name))
        for channel in ("od", "rfu"):
            np.testing.assert_array_equal(result.observations.values[channel],
                                          expected.simulation.observations.values[channel])
        assert result.variables["internal.beta_carotene"].units == "mmol"
        assert result.variables["flux.lcy"].units == "mmol/h"
        assert result.variables["content.beta_carotene"].units == "mmol/gDW"
        assert result.variables["content.beta_carotene"].basis == "content per structural biomass"
        assert result.variables["biomass_gdw"].basis == "structural biomass"
        assert result.validity.biological_validation is False
        assert result.prior_parameters
    assert "does not identify" in forecast.conditional_validity


def test_synthetic_stress_coordinates_are_not_implicitly_physical_inputs(physical_case):
    from ystwin.mech.adapters import simulate_protocol
    from ystwin.mech.inference import PhysicalParticleFilter

    protocol, genotype, initial, parameters = physical_case
    coordinates = {"upr": 0.2, "oxidative": 0.1, "burden": 0.3}
    with pytest.raises(ValueError, match="shared initial_state requires PhysicalState"):
        simulate_protocol(protocol, genotype, coordinates, parameters)
    with pytest.raises(ValueError, match="typed PhysicalState"):
        predict.predict_protocol(protocol, genotype, coordinates, parameters)
    learner = PhysicalParticleFilter(
        (initial,), parameters=parameters, genotype=genotype,
        observation_model=contracts.ObservationModel.prior(),
        measurement_noise={"od": 0.02, "rfu": 0.5}, weights=(1.0,),
    )
    with pytest.raises(contracts.ScientificRefusal, match="unidentified measurement maps"):
        learner.step(protocol, coordinates)
    assert learner.time_h == initial.time_h
