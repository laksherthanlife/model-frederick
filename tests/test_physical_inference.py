from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest

from ystwin.estimator import Innovation, ParticleFilter
from ystwin.mech.contracts import (
    Control, Event, Genotype, ObservationModel, PhysicalState, Protocol, ScientificRefusal,
)
from ystwin.mech.engine import EngineParameters, _Kernel, initialize, simulate
from ystwin.mech.inference import PhysicalParticleFilter
from ystwin.mech.params import Param, RefusedValue


MEDIUM = {"glucose": 50.0, "nitrogen": 20.0, "oxygen": 0.2, "osmolyte": 250.0, "acetate": 2.0}
CONTROL = Control(0.0, oxygen_transfer_per_h=20.0, oxygen_saturation_mM=0.2)


def optical_values(model, **values):
    return replace(model, values={
        **model.values,
        **{name: Param.asserted(name, value, model.values[name].units, "explicit synthetic test intervention")
           for name, value in values.items()},
    })


@pytest.fixture(scope="module")
def case():
    parameters = EngineParameters.prior()
    genotype = Genotype.prior()
    states = []
    for density in (0.1, 0.3):
        state = initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=density,
                           medium_mM=MEDIUM)
        state = replace(
            state, dead_biomass_gdw=0.02 * density,
            dead_intracellular_mmol={"atp": 1e-7 * density},
            intracellular_mmol={**state.intracellular_mmol, "acetate": 0.05 * density,
                                "beta_carotene": 0.01 * density, "glycogen": 0.1 * density},
            expression={**state.expression, "reporter": replace(
                state.expression["reporter"], active_mmol=0.7 * state.expression["reporter"].protein_mmol)},
        )
        states.append(state)
    return dict(ensemble=tuple(states), parameters=parameters, genotype=genotype,
                observation_model=ObservationModel.prior(seed=27),
                measurement_noise={"od": 0.04, "rfu": 0.3}, weights=(0.3, 0.7), seed=19)


def make_filter(case, **overrides):
    return PhysicalParticleFilter(**{**case, **overrides})


def interval(start=0.0, stop=0.002):
    return Protocol((start, stop), (CONTROL,))


def assert_same_particles(first, second):
    assert first == second
    for a, b in zip(first, second, strict=True):
        assert a.history == b.history
        assert a.mature_history_mmol == b.mature_history_mmol
        assert a.parameters is b.parameters
        assert a.expected_readings == b.expected_readings
        if a.simulation is not None:
            assert set(a.simulation.truth) == set(b.simulation.truth)
            for name in a.simulation.truth:
                np.testing.assert_array_equal(a.simulation.truth[name], b.simulation.truth[name])


def snapshot(pf):
    return pf.particles, pf.weights, pf.last_result, pf.n_resamples, deepcopy(pf._rng.bit_generator.state)


def assert_unchanged(pf, before):
    assert pf.particles is before[0]
    assert pf.weights == before[1]
    assert pf.last_result is before[2]
    assert pf.n_resamples == before[3]
    assert pf._rng.bit_generator.state == before[4]


def test_complete_histories_and_optical_memory_match_actual_engine_exactly(case):
    parameters = case["parameters"]
    alternatives = (parameters, parameters.with_values(glucose_transport=Param.asserted(
        "glucose_transport", 8.0, "mmol/gDW/h", "explicit synthetic parameter ensemble")))
    mature = tuple(state.expression["reporter"].protein_mmol * fraction
                   for state, fraction in zip(case["ensemble"], (0.1, 0.4), strict=True))
    start_event = Event("initial-dose", 0.0, add_mmol={"osmolyte": 0.1})
    event = Event("sample-and-feed", 0.001, withdraw_l=0.1, add_volume_l=0.04,
                  add_mmol={"glucose": 0.2, "osmolyte": 0.3})
    flow = replace(CONTROL, feed_l_h=0.1, outflow_l_h=0.1, feed_mM={"glucose": 30.0, "nitrogen": 5.0})
    protocol = Protocol((0.0, 0.001, 0.002),
                        (flow, replace(flow, time_h=0.001, temperature_c=33.0, ph=4.0)),
                        mode="chemostat", events=(start_event, event))
    pf = make_filter(case, parameters=alternatives, initial_mature_mmol=mature)
    result = pf.step(protocol, {"od": None, "rfu": None})
    assert result.time_h == 0.002
    assert result.weights == case["weights"]
    assert result.n_channels_used == 0
    assert result.innovations == ()
    assert result.prior_ess == result.posterior_ess == result.post_resample_ess
    assert not result.resampled
    assert result.parent_indices == (0, 1)
    for i, particle in enumerate(result.particles):
        direct = simulate(protocol, case["genotype"], case["ensemble"][i], alternatives[i],
                          observation=replace(case["observation_model"], initial_mature_mmol=mature[i]))
        kernel = _Kernel(alternatives[i], case["genotype"])
        assert len(kernel.names) > 100
        assert particle.parameters is alternatives[i]
        assert particle.ancestor == i
        assert particle.state == direct.final_state
        assert tuple(state.time_h for state in particle.history) == protocol.times_h
        for j, state in enumerate(particle.history):
            assert isinstance(state, PhysicalState)
            np.testing.assert_array_equal(kernel.pack(state), [direct.truth[name][j] for name in kernel.names])
            assert state.event_receipts == ((start_event,) if j == 0 else (start_event, event))
        for name in direct.truth:
            np.testing.assert_array_equal(particle.simulation.truth[name], direct.truth[name])
        np.testing.assert_array_equal(particle.mature_history_mmol, direct.observations.metadata["mature_mmol"])
        for channel in ("od", "rfu"):
            assert particle.expected_readings[channel] == direct.observations.values[channel][-1]
        assert not particle.simulation.validity.biological_validation
    assert "conditional" in result.conditional_validity
    assert "independently validate" in result.conditional_validity
    assert PhysicalParticleFilter is not ParticleFilter


def test_optics_change_readings_not_truth_and_simulation_noise_is_never_likelihood_noise(case):
    clear = case["observation_model"]
    noisy = replace(optical_values(clear, od_noise_sd=10.0, rfu_noise_sd=1000.0,
                                  missing_probability=1.0), seed=99)
    changed = optical_values(noisy, rfu_gain=2e9, inner_filter=30.0, maturation=0.5,
                             gDW_per_OD_l=0.7)
    baseline = make_filter(case).step(interval(), {})
    without_draws = make_filter(case, observation_model=noisy).step(interval(), {})
    different_optics = make_filter(case, observation_model=changed).step(interval(), {})
    assert_same_particles(baseline.particles, without_draws.particles)
    for i, (base, altered) in enumerate(zip(baseline.particles, different_optics.particles, strict=True)):
        assert base.history == altered.history
        assert base.mature_mmol != altered.mature_mmol
        assert base.expected_readings["rfu"] != altered.expected_readings["rfu"]
        assert base.expected_readings["od"] != altered.expected_readings["od"]
        for name in base.simulation.truth:
            np.testing.assert_array_equal(base.simulation.truth[name], altered.simulation.truth[name])
        expected = simulate(
            interval(), case["genotype"], case["ensemble"][i], case["parameters"],
            observation=optical_values(changed, od_noise_sd=0.0, rfu_noise_sd=0.0, missing_probability=0.0),
        )
        for channel in ("od", "rfu"):
            assert altered.expected_readings[channel] == expected.observations.values[channel][-1]
            assert not altered.simulation.observations.missing[channel].any()
    assert float(noisy.values["missing_probability"]) == 1.0
    assert float(noisy.values["rfu_noise_sd"]) == 1000.0
    assert clear.initial_mature_mmol == changed.initial_mature_mmol == 0.0


def test_weighted_likelihood_and_innovations_are_prequential_for_all_channels(case):
    pf = make_filter(case, weights=(0.2, 0.8), resample_threshold=0.0)
    prior = pf.forecast(interval())
    readings = {name: prior.particles[0].expected_readings[name] + 0.01 for name in ("od", "rfu")}
    before = np.array(pf.weights)
    log_weights = np.log(before)
    for name, reading in readings.items():
        predicted = np.array([p.expected_readings[name] for p in prior.particles])
        log_weights -= 0.5 * ((reading - predicted) / case["measurement_noise"][name])**2
    expected_weights = np.exp(log_weights - log_weights.max())
    expected_weights /= expected_weights.sum()
    result = pf.step(interval(), readings)
    np.testing.assert_allclose(result.weights, expected_weights, rtol=1e-14, atol=0.0)
    assert result.prior_ess == pytest.approx(1.0 / (before @ before))
    assert result.posterior_ess == pytest.approx(1.0 / (expected_weights @ expected_weights))
    assert result.post_resample_ess == result.posterior_ess
    assert not result.resampled
    assert result.n_channels_used == 2
    for innovation in result.innovations:
        assert isinstance(innovation, Innovation)
        predicted = np.array([p.expected_readings[innovation.channel] for p in prior.particles])
        mean = before @ predicted
        assert innovation.predicted == pytest.approx(mean, rel=1e-14)
        assert innovation.state_variance == pytest.approx(before @ (predicted - mean)**2, rel=1e-14)
        assert innovation.measurement_variance == case["measurement_noise"][innovation.channel]**2
        assert innovation.predicted != pytest.approx(expected_weights @ predicted)
        assert innovation.prior_ess == result.prior_ess
        assert innovation.posterior_ess == result.posterior_ess
        assert innovation.nis == pytest.approx(
            (readings[innovation.channel] - mean)**2 / innovation.total_variance)
    missing = pf.step(interval(0.002, 0.003), {"od": None})
    assert missing.weights == result.weights
    assert missing.innovations == ()
    assert missing.n_channels_used == 0
    assert not missing.resampled
    assert missing.prior_ess == missing.posterior_ess == missing.post_resample_ess == result.post_resample_ess
    assert len(missing.particles[0].history) == 3
    next_protocol = interval(0.003, 0.004)
    next_prior = pf.forecast(next_protocol)
    predicted = np.array([p.expected_readings["rfu"] for p in next_prior.particles])
    rfu_only = pf.step(next_protocol, {"od": None, "rfu": float(predicted[1])})
    assert rfu_only.n_channels_used == 1
    assert len(rfu_only.innovations) == 1
    assert rfu_only.innovations[0].channel == "rfu"
    assert rfu_only.innovations[0].predicted == pytest.approx(np.array(missing.weights) @ predicted)
    assert rfu_only.prior_ess == missing.post_resample_ess


@pytest.mark.parametrize("offset", [None, 0.0, np.nextafter(1.0, 0.0)], ids=["seeded", "zero", "near_one"])
@pytest.mark.parametrize("weights, target", [
    pytest.param((0.1, 0.2, 0.3, 0.4), 2, id="posterior_zeros"),
    pytest.param((1.0, 0.0, 0.0, 0.0), 0, id="trailing_zeros"),
    pytest.param((0.0, 0.0, 1.0, 0.0), 2, id="leading_and_trailing_zeros"),
    pytest.param((0.0, 0.0, 0.0, 1.0), 3, id="leading_zeros"),
])
def test_resampling_moves_complete_histories_parameters_and_maturation_with_ancestry(
    case, monkeypatch, offset, weights, target,
):
    class ResamplingRng(np.random.Generator):
        calls = 0

        def random(self):
            self.calls += 1
            draw = super().random()
            return draw if offset is None else offset

    parameters = tuple(case["parameters"].with_values(mu_max=Param.asserted(
        "mu_max", value, "1/h", "explicit synthetic parameter support")) for value in (0.2, 0.3, 0.4, 0.5))
    states = tuple(initialize(p, case["genotype"], volume_l=1.0, biomass_gdw_l=density, medium_mM=MEDIUM)
                   for p, density in zip(parameters, (0.05, 0.1, 0.2, 0.4), strict=True))
    options = dict(ensemble=states, parameters=parameters, weights=weights,
                   measurement_noise={"od": 1e-7, "rfu": 0.3}, resample_threshold=0.9,
                   initial_mature_mmol=tuple(s.expression["reporter"].protein_mmol * 0.1 for s in states))
    instrumented, plain = make_filter(case, **options), make_filter(case, **options)
    for pf in (instrumented, plain):
        monkeypatch.setattr(pf, "_rng", ResamplingRng(pf._rng.bit_generator))
    prior = instrumented.forecast(interval())
    assert instrumented._rng.calls == plain._rng.calls == 0
    reading = {"od": prior.particles[target].expected_readings["od"], "rfu": None}
    diagnosed = instrumented.step(interval(), reading, collect_innovations=True)
    undiagnosed = plain.step(interval(), reading, collect_innovations=False)
    assert_same_particles(diagnosed.particles, undiagnosed.particles)
    assert_same_particles(diagnosed.particles, (prior.particles[target],) * 4)
    expected_rng = np.random.default_rng(case["seed"])
    expected_rng.random()
    assert instrumented._rng.calls == plain._rng.calls == 1
    assert instrumented._rng.bit_generator.state == plain._rng.bit_generator.state == expected_rng.bit_generator.state
    assert instrumented.n_resamples == plain.n_resamples == 1
    assert diagnosed.weights == undiagnosed.weights == (0.25,) * 4
    assert diagnosed.prior_ess == pytest.approx(1.0 / np.sum(np.array(weights)**2))
    assert diagnosed.posterior_ess == undiagnosed.posterior_ess == 1.0
    assert diagnosed.post_resample_ess == undiagnosed.post_resample_ess == 4.0
    assert diagnosed.resampled and undiagnosed.resampled
    assert diagnosed.n_particles == diagnosed.prior_unique_particles == 4
    assert diagnosed.unique_particles == diagnosed.unique_ancestors == 1
    assert diagnosed.parent_indices == undiagnosed.parent_indices == (target,) * 4
    assert diagnosed.n_channels_used == undiagnosed.n_channels_used == 1
    assert undiagnosed.innovations == ()
    innovation = diagnosed.innovations[0]
    assert innovation.posterior_ess == 1.0
    assert innovation.post_resample_ess == 4.0
    assert innovation.unique_ancestors == innovation.unique_particles == 1
    for particle in diagnosed.particles:
        assert particle.ancestor == target
        assert particle.parameters is parameters[target]
        assert particle.history == prior.particles[target].history
        assert particle.mature_history_mmol == prior.particles[target].mature_history_mmol
    before = snapshot(instrumented)
    for _ in range(3):
        assert instrumented.unique_particles == instrumented.unique_ancestors == 1
        assert instrumented.effective_sample_size == 4.0
    assert_unchanged(instrumented, before)
    next_protocol = interval(0.002, 0.003)
    continued = instrumented.step(next_protocol, {})
    uninstrumented = plain.step(next_protocol, {}, collect_innovations=False)
    assert_same_particles(continued.particles, uninstrumented.particles)
    assert instrumented._rng.calls == plain._rng.calls == 1
    assert instrumented._rng.bit_generator.state == plain._rng.bit_generator.state == expected_rng.bit_generator.state
    assert continued.unique_ancestors == continued.unique_particles == 1
    assert continued.prior_ess == continued.posterior_ess == continued.post_resample_ess == 4.0
    assert not continued.resampled
    assert all(p.ancestor == target for p in continued.particles)


@pytest.mark.parametrize("leading_zeros", [0, 2])
def test_resampling_preserves_zero_weight_support_when_cumulative_sum_rounds_down(case, monkeypatch, leading_zeros):
    class NearOneRng:
        def random(self):
            return np.nextafter(1.0, 0.0)

    weights = (0.0,) * leading_zeros + (0.1,) * 10 + (0.0,) * 2
    pf = make_filter(case, ensemble=(case["ensemble"][0],) * len(weights), weights=weights,
                     resample_threshold=1.0)
    monkeypatch.setattr(pf, "_rng", NearOneRng())
    result = pf.step(interval(), {"od": 0.0})
    assert result.resampled
    assert result.prior_ess == result.posterior_ess == pytest.approx(10.0)
    assert result.post_resample_ess == pytest.approx(len(weights))
    assert all(weights[index] > 0 for index in result.parent_indices)
    assert result.parent_indices[0] == leading_zeros
    assert result.parent_indices[-1] == leading_zeros + 9
    assert tuple(particle.ancestor for particle in result.particles) == result.parent_indices


def test_repeated_forecasts_are_pure_and_restart_keeps_one_clock_and_event_receipts(case):
    pf = make_filter(case, resample_threshold=0.0)
    sample = Event("boundary-sample", 0.001, withdraw_l=0.1)
    first_protocol = Protocol((0.0, 0.001), (CONTROL,), events=(sample,))
    first = pf.step(first_protocol, {})
    dose = Event("later-dose", 0.002, add_mmol={"osmolyte": 0.2})
    second_protocol = Protocol((0.001, 0.002, 0.003),
                               (CONTROL, replace(CONTROL, time_h=0.002, temperature_c=32.0)),
                               events=(sample, dose))
    before = snapshot(pf)
    a, b = pf.forecast(second_protocol), pf.forecast(second_protocol)
    assert a == b
    assert_same_particles(a.particles, b.particles)
    assert_unchanged(pf, before)
    zero = pf.forecast(second_protocol, horizon_h=0.0)
    assert zero.particles is pf.particles
    assert zero.time_h == pf.time_h
    assert zero.weights == pf.weights
    assert zero.prior_ess == zero.posterior_ess == zero.post_resample_ess
    prefix = pf.forecast(second_protocol, horizon_h=0.001)
    assert prefix.time_h == 0.002
    assert_unchanged(pf, before)
    advanced = pf.step(second_protocol, {})
    assert_same_particles(a.particles, advanced.particles)
    for i, particle in enumerate(advanced.particles):
        direct = simulate(second_protocol, case["genotype"], first.particles[i].state,
                          first.particles[i].parameters,
                          observation=replace(case["observation_model"],
                                              initial_mature_mmol=first.particles[i].mature_mmol))
        assert particle.state == direct.final_state
        assert particle.mature_mmol == direct.observations.metadata["final_mature_mmol"]
        assert particle.mature_mmol > first.particles[i].mature_mmol > 0.0
        assert tuple(state.time_h for state in particle.history) == (0.0, 0.001, 0.002, 0.003)
        assert particle.state.event_receipts == (sample, dose)
        assert particle.history[1].event_receipts == (sample,)
        assert particle.history[2].event_receipts == (sample, dose)
        assert particle.ancestor == i


@pytest.mark.parametrize("horizon", [-1.0, np.nan, np.inf, True])
def test_invalid_forecast_horizons_do_not_change_filter(case, horizon):
    pf = make_filter(case)
    before = snapshot(pf)
    with pytest.raises(ValueError, match="forecast horizon_h"):
        pf.forecast(interval(), horizon_h=horizon)
    assert_unchanged(pf, before)


def test_forecast_refuses_gaps_backwards_unsupported_horizons_and_zero_time_events(case):
    pf = make_filter(case)
    before = snapshot(pf)
    with pytest.raises(ValueError, match="absolute clock"):
        pf.forecast(interval(0.001, 0.002))
    with pytest.raises(ValueError, match="within the supplied protocol"):
        pf.forecast(interval(), horizon_h=0.003)
    pending = Protocol((0.0, 0.002), (CONTROL,), events=(Event("initial-sample", 0.0, withdraw_l=0.1),))
    with pytest.raises(ValueError, match="pending event"):
        pf.forecast(pending, horizon_h=0.0)
    with pytest.raises(ValueError):
        Protocol((-0.001, 0.0), (CONTROL,))
    assert_unchanged(pf, before)
    later = make_filter(case, ensemble=tuple(replace(state, time_h=0.1) for state in case["ensemble"]))
    before = snapshot(later)
    with pytest.raises(ValueError, match="no backward time"):
        later.step(interval(), {})
    with pytest.raises(ValueError, match="no backward time"):
        later.forecast(interval())
    assert_unchanged(later, before)


def test_zero_horizon_keeps_declared_initial_states_without_invented_readings(case):
    pf = make_filter(case)
    before = snapshot(pf)
    result = pf.forecast(interval(), horizon_h=0.0)
    assert result.particles is pf.particles
    assert tuple(p.state for p in result.particles) == case["ensemble"]
    assert result.time_h == 0.0
    with pytest.raises(ValueError, match="positive shared-engine protocol interval"):
        result.particles[0].expected_readings
    assert_unchanged(pf, before)


def test_particle_diversity_counts_hidden_physical_pools_not_a_four_state_projection(case):
    state = case["ensemble"][0]
    changed = replace(state, intracellular_mmol={**state.intracellular_mmol, "glycogen": 0.02})
    pf = make_filter(case, ensemble=(state, changed))
    assert pf.particles[0].state.biomass_gdw == pf.particles[1].state.biomass_gdw
    assert pf.particles[0].state.expression == pf.particles[1].state.expression
    assert pf.unique_particles == pf.unique_ancestors == 2
    assert make_filter(case, ensemble=(state, state)).unique_particles == 1


@pytest.mark.parametrize("missing", ["ensemble", "parameters", "genotype", "observation_model", "measurement_noise", "weights"])
def test_constructor_requires_explicit_state_parameter_and_measurement_support(case, missing):
    arguments = {name: value for name, value in case.items() if name != missing}
    with pytest.raises(TypeError):
        PhysicalParticleFilter(**arguments)


@pytest.mark.parametrize("noise", [{"od": 0.1}, {"od": 0.1, "rfu": 1.0, "a450": 0.1}, None])
def test_unknown_noise_maps_are_refused(case, noise):
    with pytest.raises(ScientificRefusal, match="explicit Gaussian SDs"):
        make_filter(case, measurement_noise=noise)


@pytest.mark.parametrize("value", [0.0, -1.0, np.nan, np.inf, True, 1e308, 1e-200])
def test_invalid_or_unresolvable_noise_variances_are_rejected(case, value):
    with pytest.raises(ValueError):
        make_filter(case, measurement_noise={"od": value, "rfu": 1.0})


@pytest.mark.parametrize("weights", [(0.0, 0.0), (-1.0, 2.0), (np.nan, 1.0), (np.inf, 1.0), (1.0,), (True, 1.0)])
def test_invalid_weights_never_fall_back_to_uniform(case, weights):
    with pytest.raises(ValueError):
        make_filter(case, weights=weights)


def test_unidentified_optics_and_reporter_identity_refuse_instead_of_producing_numbers(case):
    with pytest.raises(ScientificRefusal, match="measurement map is unidentified"):
        make_filter(case, observation_model=None)
    for reporter in ("absent", "hsp70", None, ["reporter"]):
        with pytest.raises(ScientificRefusal, match="reporter gene"):
            make_filter(case, observation_model=replace(case["observation_model"], reporter=reporter))
    with pytest.raises(ScientificRefusal, match="missing optical/maturation calibration"):
        ObservationModel("reporter", {})
    unknown = Param.refused("inner_filter", "gDW/mmol", "not identified in this synthetic fixture",
                            reason="no attenuation calibration", missing="paired pigment and fluorescence measurements")
    with pytest.raises(RefusedValue):
        replace(case["observation_model"], values={**case["observation_model"].values, "inner_filter": unknown})


def test_full_state_clock_receipt_expression_and_optical_initialization_contracts(case):
    a, b = case["ensemble"]
    for ensemble, message in (
        ((), "complete PhysicalState"),
        ((a, (0.1, 0.0, 0.0, 0.3)), "complete PhysicalState"),
        ((a, replace(b, time_h=0.1)), "one absolute clock"),
        ((a, replace(b, event_receipts=(Event("already-applied", 0.0),))), "same event receipts"),
        ((a, replace(b, expression={})), "complete expression states"),
    ):
        with pytest.raises(ValueError, match=message):
            make_filter(case, ensemble=ensemble)
    with pytest.raises(ValueError, match="one explicit set per particle"):
        make_filter(case, parameters=(case["parameters"],))
    for mature in ((-1.0, 0.0), (0.0,), (1.0, 1.0)):
        with pytest.raises(ValueError):
            make_filter(case, initial_mature_mmol=mature)


@pytest.mark.parametrize("measurements, error", [
    ({"od600": 0.2}, ScientificRefusal), ({"a450": None}, ScientificRefusal),
    ({"rfu": np.nan}, ValueError), ({"od": np.inf}, ValueError), ({"od": True}, ValueError),
    (None, ValueError),
])
def test_bad_measurements_are_rejected_without_advancing_the_filter(case, measurements, error):
    pf = make_filter(case)
    before = snapshot(pf)
    with pytest.raises(error):
        pf.step(interval(), measurements)
    assert_unchanged(pf, before)


def test_zero_prior_mass_is_not_resurrected_by_observations(case):
    pf = make_filter(case, weights=(1.0, 0.0), resample_threshold=0.0)
    predicted = pf.forecast(interval())
    result = pf.step(interval(), {"od": predicted.particles[1].expected_readings["od"]})
    assert result.weights == (1.0, 0.0)
    assert result.prior_ess == result.posterior_ess == result.post_resample_ess == 1.0
    assert result.unique_ancestors == 2
    assert not result.resampled


def test_likelihood_failure_and_shared_engine_refusal_are_atomic_without_fallback(case):
    pf = make_filter(case)
    before = snapshot(pf)
    with pytest.raises(ScientificRefusal, match="prior-weight fallback"):
        pf.step(interval(), {"od": 1e308})
    assert_unchanged(pf, before)
    unsupported = Protocol((0.0, 0.002), (replace(CONTROL, temperature_c=45.0),))
    with pytest.raises(ScientificRefusal, match="operating window"):
        pf.step(unsupported, {})
    assert_unchanged(pf, before)
    a, b = case["ensemble"]
    pf = make_filter(case, ensemble=(a, replace(b, volume_l=0.05)))
    before = snapshot(pf)
    protocol = Protocol((0.0, 0.002), (CONTROL,), events=(Event("too-large-sample", 0.001, withdraw_l=0.1),))
    with pytest.raises(ValueError, match="empties the reactor"):
        pf.step(protocol, {})
    assert_unchanged(pf, before)
