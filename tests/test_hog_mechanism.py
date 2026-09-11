from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import libsbml
import numpy as np
import pytest

from ystwin.mech.hog import HogModel, HogProtocol, hog_glycerol_balance, hog_metabolic_signals


def test_hog_model_uses_published_state_network_and_explicit_osmotic_protocol():
    model = HogModel.from_source()
    times = np.array([0.0, 3600.0, 3602.5, 3605.0, 3720.0, 4200.0, 7200.0])
    result = model.simulate(times, protocol=HogProtocol(nacl_molar=0.4))
    assert result.states.shape == (len(times), 29)
    np.testing.assert_allclose(result.variables["input"][:4], [0, 0, 0.5, 1])
    assert result.variables["parameter_97"][0] == pytest.approx(0.8)
    assert np.isfinite(result.states).all()
    assert result.variables["Hog1PP_measured"][4] > result.variables["Hog1PP_measured"][1]
    assert result.variables["Hog1PP_measured"][-1] < result.variables["Hog1PP_measured"][4]
    assert result.metadata["hog_protocol"]["nacl_molar"] == 0.4
    assert result.metadata["biological_validation"] is False


def test_dose_changes_dynamics_without_mutating_published_parameters():
    model = HogModel.from_source()
    before = model.parameter_values
    times = np.array([0, 3600, 3720, 4200, 7200], dtype=float)
    low = model.simulate(times, protocol=HogProtocol(nacl_molar=0.4))
    high = model.simulate(times, protocol=HogProtocol(nacl_molar=0.8))
    assert high.variables["parameter_97"][0] == pytest.approx(1.6)
    assert not np.allclose(low.variables["Hog1PP_measured"], high.variables["Hog1PP_measured"])
    assert model.parameter_values == before
    assert not np.array_equal(low.variables["glycerol_i"], low.variables["glycerol_e"])


@pytest.mark.parametrize("dose", [-0.1, np.nan, np.inf])
def test_invalid_osmotic_inputs_are_refused(dose):
    with pytest.raises(ValueError):
        HogProtocol(nacl_molar=dose)


def test_enzyme_control_does_not_reinterpret_an_unmapped_regulator_as_capacity():
    times = np.array([0.0, 3600.0])
    variables = {"Gpd1_measured": np.array([1.0, 2.0]), "Fps1r": np.array([0.25, -1e-18]),
                 "Hog1": np.array([0.4, 0.2]), "Hog1PP": np.array([0.05, 0.25]),
                 "glycerol_measured": np.array([0.05, 0.8]), "glycerol_e": np.array([0.002, 0.003])}
    current = SimpleNamespace(times_s=times, variables=variables)
    control = SimpleNamespace(times_s=times, variables={**variables,
                              "Gpd1_measured": np.ones(2), "Fps1r": np.full(2, 0.25)})
    signals = hog_metabolic_signals(current, control)
    np.testing.assert_allclose(signals["gpd1_capacity_multiplier"], [1, 2])
    assert "fps1_transport_activity_multiplier" not in signals
    assert current.variables["Fps1r"][1] == -1e-18


def test_kinetic_overrides_cannot_silently_replace_the_input_protocol():
    model = HogModel.from_source()
    with pytest.raises(ValueError, match="protocol"):
        model.simulate([0, 3600], protocol=HogProtocol(nacl_molar=0.4),
                       parameters={"parameter_97": 0.1})


@pytest.fixture(scope="module")
def native_hog():
    model = HogModel.from_source()
    times = np.array([0., 3600., 3602.5, 3605., 3720., 4200., 5400., 7200., 9000., 14400.])
    return model, model.simulate(times, protocol=HogProtocol(0.4))


def test_wt_source_roles_and_units_are_audited_not_inferred_from_names():
    path = Path(__file__).resolve().parents[1] / "data" / "hog2013" / "model_wt.xml"
    document = libsbml.readSBML(str(path))
    source = document.getModel()
    assert source.getCompartment("intra").getConstant()
    assert source.getCompartment("intra").getSize() == 0.00025
    assert source.getCompartment("extra").getSize() == 0.5
    assert source.getSpecies("cellvol").getCompartment() == "intra"
    assert not source.getSpecies("cellvol").getConstant()
    assert source.getUnitDefinition("volume").getUnit(0).getScale() == -3
    assert source.getNumUnitDefinitions() == 1
    expected = {
        "v6": {"trioseP": -1., "glycerol_i": 1.},
        "v6b": {"trioseP": -1., "glycerol_i": 1.},
        "v13a": {"glycerol_i": -1., "glycerol_e": 1.},
        "v13b": {"glycerol_e": -1., "glycerol_i": 1.},
        "v13aBatch": {"glycerol_e": 1.},
        "v13bBatch": {"glycerol_e": -1.},
        "vVglyci": {"glycerol_i": -1.},
        "vVos": {"cellvol": 1.},
        "v1": {"glucose_e": -1., "glucose_i": 1.},
        "v1Batch": {"glucose_e": -1.},
    }
    for rid, stoichiometry in expected.items():
        reaction = source.getReaction(rid)
        actual = {
            ref.getSpecies(): sign * ref.getStoichiometry()
            for refs, sign in ((reaction.getListOfReactants(), -1),
                               (reaction.getListOfProducts(), 1))
            for ref in refs
        }
        assert actual == stoichiometry
        assert not reaction.getKineticLaw().isSetSubstanceUnits()
        assert not reaction.getKineticLaw().isSetTimeUnits()
    assert source.getReaction("v13a").getReversible()
    assert not source.getReaction("v13b").getReversible()
    assert [s.getSpecies() for s in source.getReaction("v13a").getListOfModifiers()] == ["Fps1r"]
    assert "intra *" in libsbml.formulaToL3String(source.getReaction("v6").getKineticLaw().getMath())
    assert "intra *" not in libsbml.formulaToL3String(source.getReaction("v1").getKineticLaw().getMath())
    assert libsbml.formulaToL3String(source.getRuleByVariable("glycerol_measured").getMath()) == (
        "glycerol_i * cellvol / cellvol_init")


def test_native_glycerol_chain_rule_is_exact_even_on_a_sparse_shock_grid(native_hog):
    model, trajectory = native_hog
    before = trajectory.states.copy()
    result = hog_glycerol_balance(model, trajectory)
    k = model.kinetic_model
    derivative = np.stack([k.rhs(t, state) for t, state in zip(trajectory.times_s, trajectory.states)])
    v = trajectory.variables
    chain = (derivative[:, k.species_ids.index("glycerol_i")] * v["cellvol"]
             + v["glycerol_i"] * derivative[:, k.species_ids.index("cellvol")]) / v["cellvol_init"]
    rate = {rid: trajectory.reaction_rates[:, k.reaction_ids.index(rid)] for rid in k.reaction_ids}
    factor = v["cellvol"] / (v["cellvol_init"] * v["intra"])
    rates = result.reference_rates
    np.testing.assert_allclose(rates["retention"], chain, rtol=1e-13, atol=1e-18)
    np.testing.assert_allclose(rates["synthesis"], factor * (rate["v6"] + rate["v6b"]), rtol=1e-13)
    np.testing.assert_allclose(rates["passive_transport"], factor * rate["v13a"], rtol=1e-13)
    np.testing.assert_allclose(rates["active_import"], factor * rate["v13b"], rtol=1e-13)
    np.testing.assert_allclose(rates["glucose_uptake"], factor * rate["v1"], rtol=1e-13)
    np.testing.assert_allclose(rates["volume_dilution"], -factor * rate["vVglyci"], rtol=1e-13)
    np.testing.assert_allclose(rates["volume_dilution"] + rates["volume_chain_rule"], 0., atol=1e-18)
    np.testing.assert_allclose(rates["retention"], rates["synthesis"] - rates["outward_transport"]
                               + rates["inward_transport"], rtol=1e-12, atol=1e-18)
    np.testing.assert_allclose(rates["balance_residual"], 0., atol=1e-18)
    assert not np.allclose(np.gradient(v["glycerol_measured"], trajectory.times_s), chain, atol=1e-6)
    assert abs(rates["volume_dilution"][3]) > rates["synthesis"][3]
    assert np.any(rates["retention"] < 0)
    np.testing.assert_array_equal(trajectory.states, before)
    assert result.metadata["absolute_flux_identified"] is False
    assert result.metadata["reference_rate_unit"] == "mol/L_reference/s (source metabolite convention)"
    assert result.metadata["native_training_scope"] == "glucose_only"
    assert "volume" in result.metadata["required_absolute_calibration"]
    assert result.metadata["source_model_sha256"] == trajectory.metadata["source_model_sha256"]


def test_batch_corrections_and_imposed_od_are_not_intracellular_fluxes(native_hog):
    model, trajectory = native_hog
    balance = hog_glycerol_balance(model, trajectory)
    k = model.kinetic_model
    v = trajectory.variables
    rates = {rid: trajectory.reaction_rates[:, k.reaction_ids.index(rid)] for rid in k.reaction_ids}
    external = balance.extracellular_rates
    np.testing.assert_allclose(external["net_transport"], (rates["v13a"] - rates["v13b"]) / v["extra"])
    np.testing.assert_allclose(external["batch_correction"],
                               (rates["v13aBatch"] - rates["v13bBatch"]) / v["extra"])
    np.testing.assert_allclose(external["net_accumulation"],
                               external["net_transport"] + external["batch_correction"])
    np.testing.assert_allclose(external["glucose_uptake"], (rates["v1"] + rates["v1Batch"]) / v["extra"])
    np.testing.assert_allclose(external["glucose_batch_correction"], rates["v1Batch"] / v["extra"])
    np.testing.assert_allclose(external["balance_residual"], 0., atol=1e-20)
    np.testing.assert_allclose(external["glucose_balance_residual"], 0., atol=1e-20)
    np.testing.assert_allclose(rates["v1Batch"], v["extra"] * v["v1speed"]
                               * (v["cellnum"] / v["initcellnum"] - 1))
    assert not np.allclose(external["glucose_uptake"], balance.reference_rates["glucose_uptake"])
    assert "imposed" in balance.metadata["growth_convention"]


def test_scale_free_commitment_uses_new_synthesis_not_retention_or_import(native_hog):
    model, trajectory = native_hog
    balance = hog_glycerol_balance(model, trajectory)
    fraction = balance.carbon_commitment()
    r = balance.reference_rates
    np.testing.assert_allclose(fraction, 3 * r["synthesis"] / (6 * r["glucose_uptake"]), rtol=1e-14)
    assert 0 < fraction[1] < fraction[5] < 1
    assert r["inward_transport"][5] > 0
    assert r["retention"][-1] < 0 < fraction[-1]
    assert not np.allclose(fraction, .5 * r["retention"] / r["glucose_uptake"])
    altered = replace(balance, reference_rates={**r, "active_import": 100 * r["active_import"]})
    np.testing.assert_array_equal(altered.carbon_commitment(), fraction)
    assert "medium" in balance.metadata["medium_dependence"]
    assert balance.metadata["carbon_origin_identified"] is False


@pytest.mark.parametrize("key,value", [("synthesis", -1.), ("glucose_uptake", 0.),
                                        ("glucose_uptake", -1.), ("synthesis", np.inf),
                                        ("synthesis", 100.)])
def test_carbon_credit_or_unsupported_commitment_is_refused_without_clipping(native_hog, key, value):
    model, trajectory = native_hog
    balance = hog_glycerol_balance(model, trajectory)
    altered = replace(balance, reference_rates={**balance.reference_rates,
                                                key: np.full(len(balance.times_s), value)})
    with pytest.raises(ValueError):
        altered.carbon_commitment()


def test_diagnostics_honor_effective_protocol_and_parameter_overrides_without_mutation():
    model = HogModel.from_source()
    parameters = {"kv6_1": model.parameter_values["kv6_1"] * 1.5}
    protocol = HogProtocol(0.8, shock_time_s=20.)
    trajectory = model.simulate([0, 20, 22.5, 25, 120], protocol=protocol, parameters=parameters)
    balance = hog_glycerol_balance(model, trajectory)
    k = model.kinetic_model
    overrides = {**parameters, "parameter_97": 1.6, "t_stress": 20.}
    derivative = np.stack([k.rhs(t, state, overrides) for t, state in zip(trajectory.times_s, trajectory.states)])
    v = trajectory.variables
    chain = (derivative[:, k.species_ids.index("glycerol_i")] * v["cellvol"]
             + v["glycerol_i"] * derivative[:, k.species_ids.index("cellvol")]) / v["cellvol_init"]
    np.testing.assert_allclose(balance.reference_rates["retention"], chain, rtol=1e-13, atol=1e-18)
    balance.metadata["parameter_overrides"]["kv6_1"] = -1
    assert trajectory.metadata["parameter_overrides"] == parameters
    assert model.parameter_values["kv6_1"] != parameters["kv6_1"]


def test_negative_irreversible_synthesis_cannot_hide_inside_a_positive_total():
    model = HogModel.from_source()
    trajectory = model.simulate([0.], protocol=HogProtocol(.4),
                                parameters={"kv6_1": -.1 * model.parameter_values["kv6_1"]})
    balance = hog_glycerol_balance(model, trajectory)
    assert balance.reference_rates["synthesis"][0] > 0
    with pytest.raises(ValueError, match="synthesis"):
        balance.carbon_commitment()


def test_negative_roundoff_regulator_is_preserved_and_signed_transport_is_conserved():
    model = HogModel.from_source()
    state = model.kinetic_model.initial_state({"Fps1r": -1e-18})
    trajectory = model.simulate([0], protocol=HogProtocol(.4), initial_state=state)
    balance = hog_glycerol_balance(model, trajectory)
    r = balance.reference_rates
    assert trajectory.variables["Fps1r"][0] == -1e-18
    assert r["passive_transport"][0] < 0
    assert r["outward_transport"][0] == 0
    np.testing.assert_allclose(r["outward_transport"] - r["inward_transport"],
                               r["passive_transport"] - r["active_import"], rtol=1e-14)
    assert np.isfinite(balance.carbon_commitment()).all()


def test_native_balance_refuses_wrong_source_and_nonphysical_reference_volume(native_hog):
    model, trajectory = native_hog
    wrong = replace(trajectory, metadata={**trajectory.metadata, "source_model_sha256": "other"})
    with pytest.raises(ValueError, match="source"):
        hog_glycerol_balance(model, wrong)
    invalid = model.simulate([0], protocol=HogProtocol(.4), parameters={"cellvol_init": -1.})
    with pytest.raises(ValueError, match="volume"):
        hog_glycerol_balance(model, invalid)
