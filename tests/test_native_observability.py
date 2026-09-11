from __future__ import annotations

import libsbml
import numpy as np
import pytest

from ystwin import paths
from ystwin.mech.hog import HogModel, HogProtocol


@pytest.fixture(scope="module")
def source_response():
    model = HogModel.from_source()
    times = np.array([0.0, 3600.0, 3720.0, 5400.0, 9000.0, 14400.0])
    trajectory = model.simulate(times, protocol=HogProtocol(0.4), rtol=1e-9, atol=1e-12)
    return model, trajectory


def gauge_parameters(model, scale):
    return {"kv18f_1": scale * model.parameter_values["kv18f_1"],
            "kv6_1": model.parameter_values["kv6_1"] / scale}


@pytest.mark.parametrize("scale", [0.5, 2.0, 10.0])
def test_absolute_enzyme_and_catalytic_scale_have_an_exact_source_rhs_gauge(source_response, scale):
    model, trajectory = source_response
    kinetic = model.kinetic_model
    enzyme = kinetic.species_ids.index("Gpd1")
    for time_s, state in zip(trajectory.times_s, trajectory.states):
        transformed = state.copy()
        transformed[enzyme] *= scale
        original_rhs = kinetic.rhs(time_s, state)
        changed_rhs = kinetic.rhs(time_s, transformed, gauge_parameters(model, scale))
        changed_rhs[enzyme] /= scale
        np.testing.assert_allclose(changed_rhs, original_rhs, rtol=1e-12, atol=1e-13)


def test_free_protein_observation_gain_hides_the_absolute_enzyme_scale(source_response):
    model, original = source_response
    kinetic = model.kinetic_model
    enzyme = kinetic.species_ids.index("Gpd1")
    initial = kinetic.initial_state()
    initial[enzyme] *= 2.0
    changed = model.simulate(original.times_s, protocol=HogProtocol(0.4),
                              parameters=gauge_parameters(model, 2.0), initial_state=initial,
                              rtol=1e-9, atol=1e-12)
    expected_states = original.states.copy()
    expected_states[:, enzyme] *= 2.0
    peak = np.maximum(np.max(np.abs(expected_states), axis=0), 1e-20)
    assert np.max(np.abs(changed.states - expected_states) / peak) < 1e-7
    np.testing.assert_allclose(changed.variables["Gpd1_measured"] / 2.0,
                               original.variables["Gpd1_measured"], rtol=1e-7, atol=1e-10)
    for observable in ("Hog1PP_measured", "glycerol_measured", "glycerol_e"):
        np.testing.assert_allclose(changed.variables[observable], original.variables[observable],
                                   rtol=1e-7, atol=1e-10)


def test_the_rate_parameter_is_not_merely_unused_by_the_model(source_response):
    model, _ = source_response
    kinetic = model.kinetic_model
    state = kinetic.initial_state()
    index = kinetic.species_ids.index("glycerol_i")
    original = kinetic.rhs(0.0, state)
    changed = kinetic.rhs(0.0, state, {"kv6_1": model.parameter_values["kv6_1"] / 2.0})
    assert abs(changed[index] - original[index]) > 1e-12


def test_pool_balances_cannot_resolve_common_inward_and_outward_counterflow():
    document = libsbml.readSBMLFromFile(str(paths.data_dir() / "hog2013" / "model_wt.xml"))
    source = document.getModel()
    species_ids = [species.getId() for species in source.getListOfSpecies()]

    def column(reaction_id):
        reaction = source.getReaction(reaction_id)
        result = np.zeros(len(species_ids))
        for references, sign in ((reaction.getListOfReactants(), -1), (reaction.getListOfProducts(), 1)):
            for reference in references:
                result[species_ids.index(reference.getSpecies())] += sign * reference.getStoichiometry()
        return result

    outward, inward = column("v13a"), column("v13b")
    np.testing.assert_array_equal(outward + inward, np.zeros(len(species_ids)))
    np.testing.assert_array_equal(column("v13aBatch") + column("v13bBatch"), np.zeros(len(species_ids)))
    original = 2.0 * outward + 1.0 * inward
    changed = 7.0 * outward + 6.0 * inward
    np.testing.assert_array_equal(changed, original)
    assert np.count_nonzero(outward) == 2
