from dataclasses import replace

import numpy as np
import pytest

from ystwin.mech.contracts import (
    Compartment, Control, Event, ExpressionState, Genotype, ObservationModel, PhysicalState,
    Promoter, Protocol, ScientificRefusal, ThermodynamicRequest, prior,
)
from ystwin.mech.engine import EngineParameters, initialize, simulate
from ystwin.mech.params import Param, RefusedValue, SweptValue, Tag


@pytest.mark.parametrize("value", [True, np.bool_(False), float("nan"), float("inf"), 1j, "1"])
def test_shared_physical_numbers_do_not_accept_boolean_complex_string_or_nonfinite(value):
    with pytest.raises(ValueError):
        Control(value)
    with pytest.raises(ValueError):
        PhysicalState(0.0, value, 0.0, {})


def test_concentration_constructor_is_an_explicit_millimolar_litre_amount_conversion():
    state = PhysicalState.from_concentrations(volume_l=0.2, biomass_gdw_l=1.5,
                                             medium_mM={"glucose": 12.0})
    assert state.biomass_gdw == pytest.approx(0.3)
    assert state.extracellular_mmol["glucose"] == pytest.approx(2.4)
    assert state.intracellular_mmol["glucose"] == 0.0
    with pytest.raises(ScientificRefusal, match="unsupported"):
        PhysicalState.from_concentrations(volume_l=1, biomass_gdw_l=1,
                                         medium_mM={"glucose_g_l": 12})


def test_contracts_snapshot_mutable_inputs_instead_of_aliasing_them():
    feed = {"glucose": 10.0}
    event_amounts = {"nitrogen": 2.0}
    control = Control(0, feed_mM=feed)
    event = Event("pulse", 0.5, add_mmol=event_amounts)
    state = PhysicalState(0, 1, 0, feed)
    feed["glucose"] = 999
    event_amounts["nitrogen"] = 999
    assert control.feed_mM["glucose"] == 10
    assert state.extracellular_mmol["glucose"] == 10
    assert event.add_mmol["nitrogen"] == 2
    with pytest.raises(TypeError):
        control.feed_mM["glucose"] = 7


@pytest.mark.parametrize("times", [(0,), (0, 0), (1, 0), (0, np.nan), (-1, 1)])
def test_protocol_has_one_strict_absolute_clock(times):
    with pytest.raises(ValueError):
        Protocol(times, (Control(0),))


def test_control_and_event_boundaries_are_independent_of_read_grid():
    protocol = Protocol((4, 8), (Control(0), Control(5.25)),
                        events=(Event("pulse", 6.125), Event("sample", 6.125)))
    assert protocol.boundaries_h == (4.0, 5.25, 6.125, 8.0)
    assert protocol.control_at(5.249).time_h == 0
    assert protocol.control_at(5.25).time_h == 5.25
    with pytest.raises(ValueError, match="ambiguous"):
        Protocol((0, 1), (Control(0), Control(0)))
    with pytest.raises(ValueError, match="unique"):
        Protocol((0, 1), (Control(0),), events=(Event("same", 0), Event("same", 1)))
    with pytest.raises(ValueError, match="inside"):
        Protocol((0, 1), (Control(0),), events=(Event("late", 2),))


@pytest.mark.parametrize("mode,control", [
    ("batch", Control(0, feed_l_h=0.1)),
    ("fedbatch", Control(0, outflow_l_h=0.1)),
    ("chemostat", Control(0, feed_l_h=0.1, outflow_l_h=0.2)),
])
def test_vessel_names_enforce_authoritative_flow_semantics(mode, control):
    with pytest.raises(ValueError):
        Protocol((0, 1), (control,), mode=mode)


def test_unimplemented_inputs_are_refused_by_typed_species_and_modes():
    with pytest.raises(ScientificRefusal, match="vessel"):
        Protocol((0, 1), (Control(0),), mode="perfusion_cell_retention")
    for name in ("DTT", "tunicamycin", "sorbate", "amino_acids", "galactose"):
        with pytest.raises(ScientificRefusal, match="unsupported"):
            Control(0, feed_mM={name: 1.0})
    with pytest.raises(ScientificRefusal, match="copy-number"):
        Genotype((), inheritance="episomal")


def test_required_parameters_never_pick_a_midpoint_or_silently_assume_units():
    parameters = EngineParameters.prior()
    values = dict(parameters.values)
    values.pop("snf1_tau")
    with pytest.raises(ScientificRefusal, match="snf1_tau"):
        EngineParameters(values, allow_prior=True)
    with pytest.raises(ValueError, match="units"):
        parameters.with_values(snf1_tau=prior("snf1_tau", 6, "min", "test"))
    with pytest.raises(ValueError, match="Param"):
        parameters.with_values(snf1_tau=0.1)
    with pytest.raises(SweptValue):
        parameters.with_values(snf1_tau=Param.swept("snf1_tau", "h", "test axis", bounds=(0.01, 1)))
    with pytest.raises(RefusedValue):
        parameters.with_values(snf1_tau=Param.refused("snf1_tau", "h", "not measured",
                                                    reason="no kinetics", missing="activation trace"))
    with pytest.raises(ScientificRefusal, match="allow_prior"):
        EngineParameters(parameters.values)


def test_parameter_provenance_survives_overrides_without_mutating_the_base():
    base = EngineParameters.prior()
    changed = base.with_values(glucose_transport=Param.measured("transport", 5, "mmol/gDW/h",
                                                               "synthetic fixture, not a publication"))
    assert float(base.values["glucose_transport"]) == 12
    assert changed.values["glucose_transport"].tag == Tag.MEASURED
    assert base.values["snf1_tau"].tag == Tag.ASSERTED
    assert "PRIOR" in base.values["snf1_tau"].source
    with pytest.raises(TypeError):
        base.values["snf1_tau"] = changed.values["glucose_transport"]


def test_promoters_require_known_signal_roles_and_explicit_signed_gain():
    basal = prior("basal", 0.5, "dimensionless", "test")
    promoter = Promoter(basal, {"pka": prior("pka", -0.5, "dimensionless", "test repression")})
    assert promoter.activity({"pka": 1}) == 0
    assert promoter.activity({"pka": 0}) == 0.5
    with pytest.raises(ScientificRefusal, match="transcription factor"):
        Promoter(basal, {"not_a_tf": basal})
    gene = Genotype.prior().genes[0]
    with pytest.raises(ValueError, match="integer"):
        replace(gene, copies=True)
    with pytest.raises(ValueError, match="coding sequence"):
        replace(gene, transcript_nt=1)
    with pytest.raises(ValueError, match="unique"):
        Genotype((gene, gene))


def test_expression_active_is_a_subset_and_state_receipts_are_historical():
    with pytest.raises(ValueError, match="subset"):
        ExpressionState(protein_mmol=1, active_mmol=2)
    with pytest.raises(ValueError, match="positive live"):
        PhysicalState(0, 1, 0, {}, intracellular_mmol={"atp": 1})
    with pytest.raises(ValueError, match="clock"):
        PhysicalState(0, 1, 0, {}, event_receipts=(Event("future", 1),))
    with pytest.raises(ValueError, match="fractions"):
        PhysicalState(0, 1, 0, {}, signals={"pka": 2})


def test_optical_calibration_is_not_defaulted_from_physical_biomass():
    with pytest.raises(ScientificRefusal, match="calibration"):
        ObservationModel("reporter", {})
    model = ObservationModel.prior()
    with pytest.raises(RefusedValue):
        replace(model, values={**model.values, "rfu_gain": Param.refused(
            "gain", "RFU L/mmol", "not measured", reason="reader-specific", missing="purified standard")})
    with pytest.raises(ValueError, match="at most"):
        replace(model, values={**model.values, "missing_probability": prior(
            "missing", 1.1, "dimensionless", "invalid")})


def test_scientific_refusals_do_not_silently_request_a_different_model():
    p, g = EngineParameters.prior(), Genotype(())
    state = PhysicalState(0, 1, 0, {})
    with pytest.raises(ScientificRefusal, match="temperature"):
        simulate(Protocol((0, 1), (Control(0, temperature_c=50),)), g, state, p)
    with pytest.raises(ScientificRefusal, match="pH"):
        simulate(Protocol((0, 1), (Control(0, ph=2),)), g, state, p)
    with pytest.raises(ValueError, match="absolute time"):
        simulate(Protocol((1, 2), (Control(0),)), g, state, p)
    with pytest.raises(ValueError, match="empties"):
        simulate(Protocol((0, 1), (Control(0),), events=(Event("empty", 0.5, withdraw_l=1),)),
                 g, state, p)


def test_complete_thermodynamic_coverage_is_an_explicit_refusal():
    from ystwin.bridge.thermodynamic import ThermodynamicData

    p = replace(EngineParameters.prior(), thermodynamics=ThermodynamicRequest(
        ThermodynamicData({}, 0.0), "empty thermodynamic fixture", require_complete=True))
    with pytest.raises(ScientificRefusal, match="complete thermodynamic"):
        simulate(Protocol((0, 1), (Control(0),)), Genotype(()), PhysicalState(0, 1, 0, {}), p)


def test_result_every_variable_has_units_compartment_and_separated_observation():
    p, g = EngineParameters.prior(), Genotype.prior()
    initial = initialize(p, g, volume_l=0.001, biomass_gdw_l=0.2,
                         medium_mM={"glucose": 50, "nitrogen": 10, "oxygen": 0.2, "osmolyte": 200})
    result = simulate(Protocol((0, 0.02), (Control(0),)), g, initial, p)
    assert set(result.truth) == set(result.variables)
    assert result.variables["internal.atp"].units == "mmol"
    assert result.variables["internal.atp"].compartment is Compartment.CELL_WATER
    assert result.variables["gene.crte.active"].basis == "subset of protein"
    assert result.variables["medium_mM.oxygen"].units == "mM"
    assert result.observations is None
    assert not result.validity.biological_validation
    assert result.validity.parameter_basis.startswith("prior-conditional")
    assert "snf1_tau" in result.prior_parameters
    with pytest.raises(ValueError):
        result.trace("internal.atp")[0] = 1
    with pytest.raises(TypeError):
        result.truth["unexpected"] = np.zeros(2)


def test_active_acid_export_requires_complete_explicit_parameter_provenance():
    from ystwin.mech.contracts import AcidExportParameters

    assert EngineParameters.prior().acid_export is None
    export = AcidExportParameters.prior()
    enabled = EngineParameters.prior(acid_export=export)
    assert enabled.acid_export is export
    for name in ("vmax", "km_anion", "atp_per_proton", "atp_per_anion"):
        assert getattr(export, name).tag == Tag.ASSERTED
    with pytest.raises(ValueError, match="units"):
        replace(export, vmax=prior("vmax", 1, "mmol/L/h", "test wrong basis"))
    with pytest.raises(RefusedValue):
        replace(export, km_anion=Param.refused("km", "mM", "unmeasured export affinity",
                                             reason="no transport assay", missing="anion uptake/export kinetics"))
    with pytest.raises(ValueError, match="AcidExportParameters"):
        replace(enabled, acid_export={"vmax": 1})


def test_constraints_do_not_guess_activities_or_molecular_identities_from_total_pools():
    from ystwin.bridge.thermodynamic import ThermodynamicData

    data = ThermodynamicData({"lycopene_c": 0, "betacarotene_c": -10}, 1)
    with pytest.raises(ScientificRefusal, match="activity_coefficients"):
        ThermodynamicRequest(data, "synthetic fixture", mode="constrain_supported")
    with pytest.raises(ScientificRefusal, match="identity mapping"):
        ThermodynamicRequest(data, "synthetic fixture", reactions=("adenylate_kinase",))
    for reaction in ("pdh", "acid_export", "respiration", "crti"):
        with pytest.raises(ScientificRefusal, match="molecular stoichiometry"):
            ThermodynamicRequest(data, "synthetic fixture", reactions=(reaction,))
    with pytest.raises(ValueError, match="one-to-one"):
        ThermodynamicRequest(data, "synthetic fixture", metabolite_ids={"lycopene": "same", "beta_carotene": "same"})
