import importlib
from dataclasses import replace

import numpy as np
import pytest

from ystwin.mech import engine
from ystwin.mech.contracts import (
    Control, Event, Genotype, ObservationModel, PhysicalState, Protocol,
)
from ystwin.mech.params import Tag
from ystwin.mech.state import TIMESCALE_CATALOGUE


MODULES = ("ystwin.mech.adapters", "ystwin.generator.culture", "ystwin.fba.dynamic", "ystwin.hybrid")


@pytest.fixture
def request_values():
    parameters = engine.EngineParameters.prior()
    genotype = Genotype.prior()
    state = engine.initialize(parameters, genotype, volume_l=0.01, biomass_gdw_l=0.2,
                              medium_mM={"glucose": 30, "nitrogen": 10, "oxygen": 0.2,
                                         "osmolyte": 250})
    protocol = Protocol((0, 0.015, 0.03), (Control(0, oxygen_transfer_per_h=20,
                                                oxygen_saturation_mM=0.2),),
                        events=(Event("sample", 0.015, withdraw_l=0.001),))
    return protocol, genotype, state, parameters


@pytest.mark.parametrize("module_name", MODULES)
def test_protocol_entry_is_exactly_one_shared_engine_call_without_rebuilding_inputs(module_name, monkeypatch, request_values):
    module = importlib.import_module(module_name)
    original = engine.simulate
    calls, outputs = [], []
    def recording(*args, **kwargs):
        calls.append((args, kwargs))
        result = original(*args, **kwargs)
        outputs.append(result)
        return result
    monkeypatch.setattr(engine, "simulate", recording)
    observation = ObservationModel.prior()
    result = module.simulate_protocol(*request_values, observation=observation, rtol=2e-6)
    assert len(calls) == 1
    assert all(actual is supplied for actual, supplied in zip(calls[0][0], request_values, strict=True))
    assert calls[0][1]["observation"] is observation
    assert calls[0][1]["rtol"] == 2e-6
    assert result is outputs[0]
    np.testing.assert_array_equal(result.times_h, request_values[0].times_h)
    assert result.diagnostics["events_applied"] == ("sample",)


@pytest.mark.parametrize("module_name", MODULES)
def test_protocol_adapter_refuses_missing_or_legacy_shaped_inputs_without_choosing_priors(module_name, monkeypatch, request_values):
    module = importlib.import_module(module_name)
    def not_called(*args, **kwargs):
        pytest.fail("a rejected shared request must not construct priors or enter a solver")
    monkeypatch.setattr(engine.EngineParameters, "prior", not_called)
    monkeypatch.setattr(engine, "simulate", not_called)
    with pytest.raises(TypeError):
        module.simulate_protocol(request_values[0])
    for position, invalid in enumerate(({"times_h": (0, 1)}, {"crte": 1}, {"volume_l": 1}, {})):
        values = list(request_values)
        values[position] = invalid
        with pytest.raises(ValueError, match="shared"):
            module.simulate_protocol(*values)
    with pytest.raises(ValueError, match="ObservationModel"):
        module.simulate_protocol(*request_values, observation={})


def test_all_protocol_entries_return_identical_physical_truth_and_do_not_call_comparisons(monkeypatch, request_values):
    from ystwin import hybrid
    from ystwin.fba import dynamic
    from ystwin.generator import culture

    expected = engine.simulate(*request_values)
    def forbidden(*args, **kwargs):
        pytest.fail("shared protocol dispatch entered a comparison simulator")
    for module, name in ((hybrid, "simulate_hybrid_comparison"),
                         (dynamic, "simulate_fixed_volume_comparison"),
                         (culture, "simulate_empirical_comparison")):
        monkeypatch.setattr(module, name, forbidden)
        result = module.simulate_protocol(*request_values)
        for key in expected.truth:
            np.testing.assert_array_equal(result.truth[key], expected.truth[key])
        assert result.validity == expected.validity
        assert isinstance(result.final_state, PhysicalState)


def test_existing_fixed_volume_entry_is_an_explicitly_named_unchanged_comparison():
    from ystwin.fba import dynamic

    assert dynamic.simulate_batch is dynamic.simulate_fixed_volume_comparison
    def rates(substrate):
        factor = substrate / (1 + substrate)
        return 0.1 * factor, factor, 0.02 * factor
    args = (rates, 0.1, 10.0, 1.0, 100.0)
    comparison = dynamic.simulate_fixed_volume_comparison(*args, steps=20)
    compatibility = dynamic.simulate_batch(*args, steps=20)
    np.testing.assert_array_equal(comparison.biomass_g_per_l, compatibility.biomass_g_per_l)
    np.testing.assert_array_equal(comparison.product_mmol_per_l, compatibility.product_mmol_per_l)
    assert comparison.simulation_route == "fixed_volume_comparison"


def test_existing_empirical_culture_entry_is_an_explicitly_named_unchanged_comparison():
    from ystwin.generator import culture

    assert culture.simulate_culture is culture.simulate_empirical_comparison
    parameters = culture.CultureParameters(0.3, 1, 1, 1e-3, 2e-3, 1, 1, 2)
    args = (np.linspace(0, 1, 6), 0.2, parameters, 0.05)
    comparison = culture.simulate_empirical_comparison(*args)
    compatibility = culture.simulate_culture(*args)
    for name in ("biomass", "growth_rate", "promoter_activity", "reporter"):
        np.testing.assert_array_equal(comparison[name], compatibility[name])
    assert comparison["simulation_route"] == "empirical_culture_comparison"


def test_hybrid_comparison_name_preserves_the_old_entry_without_a_shadow_protocol_schema():
    from ystwin import hybrid
    from ystwin.mech import adapters

    assert hybrid.simulate_hybrid_case is hybrid.simulate_hybrid_comparison
    assert hybrid.simulate_protocol is adapters.simulate_protocol
    assert hybrid.simulate_hybrid_comparison.__name__ == "simulate_hybrid_comparison"


def test_hsf1_timescale_is_a_fitted_source_prior_not_a_measured_binding_rate():
    parameter = engine.EngineParameters.prior().values["hsf1_tau"]
    row = TIMESCALE_CATALOGUE["hsf1_free"]
    assert parameter.tag == Tag.ASSERTED
    assert "FITTED" in row.source
    assert "arbitrary" in row.source
    assert "10.7554/eLife.18638.021" in row.source
    assert "data/parameter_evidence.json" in parameter.source
    assert "MEASURED." not in row.source
    assert float(parameter) == pytest.approx(row.tau_h[0])


def test_shared_clock_accepts_a_history_checkpoint_via_a_different_adapter(request_values):
    from ystwin.fba import dynamic
    from ystwin.generator import culture

    protocol, genotype, state, parameters = request_values
    first = culture.simulate_protocol(protocol, genotype, state, parameters)
    resumed = dynamic.simulate_protocol(Protocol((0.03, 0.04), protocol.controls), genotype,
                                        first.final_state, parameters)
    direct = engine.simulate(Protocol((0.03, 0.04), protocol.controls), genotype,
                             first.final_state, parameters)
    assert resumed.diagnostics["events_applied"] == ("sample",)
    for name in direct.truth:
        np.testing.assert_array_equal(resumed.truth[name], direct.truth[name])
    with pytest.raises(ValueError, match="absolute time"):
        culture.simulate_protocol(replace(protocol, times_h=(0, 0.04)), genotype,
                                  first.final_state, parameters)
