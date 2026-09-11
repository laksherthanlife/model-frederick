import hashlib
import re
from dataclasses import replace

import numpy as np
import pytest

from ystwin.analysis.parameter_evidence import (
    JALIHAL_SBML_SHA256, load_parameter_evidence, repository_root, source_path,
)
from ystwin.mech import oxidative, ph, upr
from ystwin.mech.contracts import (
    PROMOTER_COUPLING_NOT_BUILT, PROMOTER_UNCOUPLED_SIGNALS, SIGNALS,
    AcidExportParameters, Control, Event, Genotype, ObservationModel, PhysicalState, Promoter,
    ScientificRefusal,
    Protocol, ThermodynamicRequest, prior, sourced_promoter_response,
)
from ystwin.mech.engine import EngineParameters, _Kernel, initialize, simulate
from ystwin.mech.params import RefusedValue, Tag
from ystwin.mech.signalling import CARBON_NITROGEN_FLUX_COUPLING, NativeHogResponse


MEDIUM = {"glucose": 50.0, "nitrogen": 20.0, "oxygen": 0.2, "osmolyte": 250.0}
AERATED = Control(0.0, oxygen_transfer_per_h=30.0, oxygen_saturation_mM=0.2)


def parameters_with(**values):
    parameters = EngineParameters.prior()
    return parameters.with_values(**{name: prior(name, value, parameters.values[name].units,
                                                   "test intervention") for name, value in values.items()})


def run_case(*, medium=None, controls=(AERATED,), events=(), times=(0.0, 0.25, 0.5, 1.0),
             parameters=None, genotype=None, mode="batch", state=None, **options):
    parameters = EngineParameters.prior() if parameters is None else parameters
    genotype = Genotype.prior() if genotype is None else genotype
    initial = initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=0.2,
                         medium_mM=MEDIUM if medium is None else medium,
                         time_h=times[0]) if state is None else state
    return simulate(Protocol(times, controls, mode=mode, events=events), genotype, initial,
                    parameters, **options)


@pytest.fixture(scope="module")
def reference():
    return run_case()


def test_shared_upr_current_growth_does_not_redefine_reference_synthesis():
    chain = upr.UprChain(0.2, 4.5, 0.1, 4.0, 0.3)
    state = np.array([chain.basal()[name] for name in upr.UPR_STATES.names])
    reference = np.array(chain.rhs(0.2, 0.0)(0.0, state))
    stopped = np.array(chain.rhs(0.2, 0.0, growth_rate_per_h=0.0)(0.0, state))
    np.testing.assert_allclose(reference, 0.0, atol=1e-13)
    np.testing.assert_allclose(stopped - reference, 0.3 * state, atol=1e-13)


def test_shared_peroxide_uses_current_physical_density_not_optical_calibration():
    density = float(oxidative.X_REF_CELLS_PER_ML)
    assert oxidative.dose_decay_from_density_per_h(density) == float(oxidative.K_REF_PER_H)
    assert oxidative.dose_decay_from_density_per_h(2 * density) == 2 * float(oxidative.K_REF_PER_H)
    for value in (-1.0, np.nan, np.inf, True):
        with pytest.raises(ValueError):
            oxidative.dose_decay_from_density_per_h(value)


def test_all_reaction_columns_and_full_trajectory_conserve_elements(reference):
    names = reference.diagnostics["stoichiometric_state_names"]
    matrix = reference.diagnostics["stoichiometric_matrix"]
    carbon = np.array([reference.variables[name].carbon for name in names])
    nitrogen = np.array([reference.variables[name].nitrogen for name in names])
    carbon[names.index("ledger.co2")] = 1.0
    np.testing.assert_allclose(carbon @ matrix, 0, atol=1e-11)
    np.testing.assert_allclose(nitrogen @ matrix, 0, atol=1e-11)
    for element in ("carbon", "nitrogen", "live_biomass"):
        np.testing.assert_allclose(reference.trace(f"balance.{element}_residual"), 0, atol=2e-9)
    assert reference.diagnostics["state_projection"] == "none"


@pytest.mark.parametrize("mode", ["batch", "fedbatch", "chemostat"])
def test_feed_and_flow_have_analytic_amount_not_just_concentration_solutions(mode):
    rate = 0.2 if mode != "batch" else 0.0
    control = Control(0, feed_l_h=rate, outflow_l_h=rate if mode == "chemostat" else 0.0,
                      feed_mM={"glucose": 20, "nitrogen": 8})
    times = np.linspace(0, 2, 9)
    result = simulate(Protocol(times, (control,), mode=mode), Genotype(()),
                      PhysicalState(0, 1, 0, {"glucose": 10, "nitrogen": 4}), EngineParameters.prior())
    if mode == "batch":
        concentration = np.full_like(times, 10.0)
        amount = concentration
        volume = np.ones_like(times)
    elif mode == "fedbatch":
        amount = 10 + 4 * times
        volume = 1 + rate * times
        concentration = amount / volume
    else:
        volume = np.ones_like(times)
        concentration = 20 - 10 * np.exp(-rate * times)
        amount = concentration
    np.testing.assert_allclose(result.trace("volume_l"), volume, rtol=2e-6)
    np.testing.assert_allclose(result.trace("external.glucose"), amount, rtol=2e-6)
    np.testing.assert_allclose(result.trace("medium_mM.glucose"), concentration, rtol=2e-6)
    np.testing.assert_allclose(result.trace("balance.carbon_residual"), 0, atol=1e-10)
    np.testing.assert_allclose(result.trace("balance.nitrogen_residual"), 0, atol=1e-10)


@pytest.mark.parametrize("mode", ["batch", "fedbatch", "chemostat"])
def test_live_reactor_carbon_nitrogen_and_birth_budgets_survive_feed_death_and_sampling(mode):
    flow = 0.05 if mode != "batch" else 0.0
    control = replace(AERATED, feed_l_h=flow, outflow_l_h=flow if mode == "chemostat" else 0.0,
                      feed_mM={"glucose": 100, "nitrogen": 50, "acetate": 2, "osmolyte": 250})
    result = run_case(controls=(control,), mode=mode,
                      events=(Event("dose", 0.13, add_volume_l=0.1, add_mmol={"peroxide": 0.1, "acetate": 3}),
                              Event("sample", 0.43, withdraw_l=0.2)),
                      times=(0, 0.13, 0.43, 0.7))
    for name in ("carbon", "nitrogen", "live_biomass"):
        np.testing.assert_allclose(result.trace(f"balance.{name}_residual"), 0, atol=3e-9)
    assert result.trace("ledger.biomass_died")[-1] > 0
    assert result.trace("ledger.biomass_washed_out")[-1] > 0
    assert result.trace("ledger.biomass_formed")[-1] > 0


def test_washout_is_not_growth_and_death_moves_inventories_to_dead_cells():
    p = parameters_with(mu_max=0, death_baseline=0.1, death_starvation=0, death_stress=0)
    control = replace(AERATED, feed_l_h=0.2, outflow_l_h=0.2)
    result = run_case(parameters=p, genotype=Genotype(()), controls=(control,), mode="chemostat",
                      times=np.linspace(0, 1, 6))
    t = result.times_h
    np.testing.assert_allclose(result.trace("biomass_gdw"), 0.2 * np.exp(-0.3 * t), rtol=3e-6)
    np.testing.assert_allclose(result.trace("dead_biomass_gdw"), 0.2 * np.exp(-0.2 * t) * (1 - np.exp(-0.1 * t)), rtol=3e-6, atol=1e-11)
    assert np.all(result.trace("growth_per_h") == 0)
    np.testing.assert_allclose(result.trace("death_per_h"), 0.1)
    np.testing.assert_allclose(result.trace("washout_per_h"), 0.2)
    assert result.trace("dead.atp")[-1] > 0
    assert result.trace("ledger.biomass_formed")[-1] == 0


def test_exact_event_time_ordering_final_time_events_and_restart_receipts():
    p, g = EngineParameters.prior(), Genotype(())
    state = PhysicalState(0, 1, 0, {"glucose": 10})
    events = (Event("sample", 0.37, withdraw_l=0.2),
              Event("feed", 0.37, add_volume_l=0.5, add_mmol={"glucose": 3}),
              Event("last", 1, add_mmol={"glucose": 2}))
    result = simulate(Protocol((0, 0.369, 0.37, 1), (Control(0),), events=events), g, state, p)
    np.testing.assert_allclose(result.trace("external.glucose"), (10, 10, 11, 13))
    np.testing.assert_allclose(result.trace("volume_l"), (1, 1, 1.3, 1.3))
    assert result.diagnostics["events_applied"] == ("sample", "feed", "last")
    restart = simulate(Protocol((1, 2), (Control(0),), events=(events[-1],)), g, result.final_state, p)
    np.testing.assert_allclose(restart.trace("external.glucose"), (13, 13))
    with pytest.raises(ValueError, match="receipt"):
        simulate(Protocol((1, 2), (Control(0),), events=(replace(events[-1], add_mmol={"glucose": 4}),)),
                 g, result.final_state, p)


def test_dense_readout_does_not_create_a_second_clock_or_lose_a_narrow_intervention():
    controls = (AERATED, replace(AERATED, time_h=0.217, folding_inhibition=0.9),
                replace(AERATED, time_h=0.229))
    sparse = run_case(controls=controls, times=(0, 0.5, 1))
    dense = run_case(controls=controls, times=np.linspace(0, 1, 21))
    for name in sparse.truth:
        np.testing.assert_array_equal(sparse.trace(name), dense.trace(name)[[0, 10, 20]])
    untreated = run_case(times=(0, 0.5, 1))
    assert sparse.trace("upr.hac1_protein")[1] > untreated.trace("upr.hac1_protein")[1]
    assert sparse.diagnostics["event_boundaries_h"] == (0, 0.217, 0.229, 1)


def test_history_restart_and_time_translation_agree_without_reinitializing_native_states():
    p, g = EngineParameters.prior(), Genotype.prior()
    full = run_case(parameters=p, genotype=g, controls=(AERATED, replace(AERATED, time_h=0.4, folding_inhibition=0.8)),
                    events=(Event("salt", 0.4, add_mmol={"osmolyte": 200}),), times=(0, 0.4, 1))
    first = run_case(parameters=p, genotype=g, events=(Event("salt", 0.4, add_mmol={"osmolyte": 200}),), times=(0, 0.4))
    second = run_case(parameters=p, genotype=g, state=first.final_state,
                      controls=(replace(AERATED, time_h=0.4, folding_inhibition=0.8),),
                      events=(Event("salt", 0.4, add_mmol={"osmolyte": 200}),), times=(0.4, 1))
    for name in full.diagnostics["stoichiometric_state_names"]:
        np.testing.assert_allclose(full.trace(name)[-1], second.trace(name)[-1], rtol=2e-5, atol=2e-10)
    shifted = run_case(times=(10, 10.25, 10.5, 11), controls=(replace(AERATED, time_h=10),))
    unshifted = run_case()
    for name in ("biomass_gdw", "signal.hog1", "upr.hac1_protein", "internal.beta_carotene"):
        np.testing.assert_allclose(shifted.trace(name), unshifted.trace(name), rtol=2e-6, atol=1e-10)


def test_prior_conditional_signalling_has_the_requested_carbon_nitrogen_and_stress_routes(reference):
    carbon = run_case(medium={**MEDIUM, "glucose": 0, "ethanol": 50})
    nitrogen = run_case(medium={**MEDIUM, "nitrogen": 0})
    oxidative_stress = run_case(medium={**MEDIUM, "peroxide": 0.2})
    heat = run_case(controls=(replace(AERATED, temperature_c=37),))
    osmotic = run_case(medium={**MEDIUM, "osmolyte": 650})
    er = run_case(controls=(replace(AERATED, folding_inhibition=0.8),))
    assert carbon.trace("signal.pka")[-1] < reference.trace("signal.pka")[-1]
    assert carbon.trace("signal.snf1")[-1] > reference.trace("signal.snf1")[-1]
    assert carbon.trace("external.ethanol")[-1] < 50
    assert nitrogen.trace("signal.torc1")[-1] < reference.trace("signal.torc1")[-1]
    assert oxidative_stress.trace("signal.yap1")[-1] > reference.trace("signal.yap1")[-1]
    assert heat.trace("signal.hsf1")[-1] > reference.trace("signal.hsf1")[-1]
    assert osmotic.trace("signal.hog1")[-1] > reference.trace("signal.hog1")[-1]
    assert er.trace("upr.hac1_protein")[-1] > reference.trace("upr.hac1_protein")[-1]
    for altered, gene in ((oxidative_stress, "antioxidant"), (heat, "hsp70"), (osmotic, "gpd1"), (er, "reporter")):
        assert altered.trace(f"promoter.{gene}")[-1] > reference.trace(f"promoter.{gene}")[-1]
        assert altered.trace(f"gene.{gene}.mrna")[-1] > 0
        assert altered.trace(f"flux.{gene}.protein_synthesis")[-1] > 0
        assert np.all(altered.trace(f"gene.{gene}.active") <= altered.trace(f"gene.{gene}.protein"))


_CARBON_FLUX_SIGN = {"glucose_transport": 1, "ethanol_oxidation": 1, "glycogen_mobilization": 1,
                     "growth": 1, "glycogen_synthesis": -1}
_RESPIRATORY = {**MEDIUM, "glucose": 1.0, "ethanol": 50.0}
_CARBON_ROWS = ("carbon.pka.w_pka_camp", "carbon.snf1.w_mig_snf", "nitrogen.tor.w_sch9_torc")


def test_the_declared_carbon_nitrogen_flux_map_is_complete_and_signed():
    declared = {flux for fluxes in CARBON_NITROGEN_FLUX_COUPLING.values() for flux in fluxes}
    assert declared == set(_CARBON_FLUX_SIGN)
    assert set(PROMOTER_UNCOUPLED_SIGNALS) <= set(CARBON_NITROGEN_FLUX_COUPLING)
    assert set(CARBON_NITROGEN_FLUX_COUPLING) <= set(SIGNALS)


@pytest.mark.parametrize("signal", sorted(CARBON_NITROGEN_FLUX_COUPLING))
def test_each_carbon_and_nitrogen_state_moves_every_flux_it_declares(signal):
    medium = _RESPIRATORY if signal == "snf1" else MEDIUM
    live = run_case(medium=medium)
    silenced = run_case(medium=medium,
                        genotype=replace(Genotype.prior(), signal_activity={signal: 0.0}))
    assert live.trace(f"signal.{signal}")[-1] > 0.5
    assert silenced.trace(f"signal.{signal}")[-1] == pytest.approx(0.0, abs=1e-9)
    for flux in CARBON_NITROGEN_FLUX_COUPLING[signal]:
        change = live.trace(f"flux.{flux}")[-1] - silenced.trace(f"flux.{flux}")[-1]
        assert _CARBON_FLUX_SIGN[flux] * change > 0


def test_no_promoter_in_the_shipped_genotype_reads_a_carbon_state(reference):
    genotype = Genotype.prior()
    assert genotype.promoter_drivers == ("torc1", "yap1", "hsf1", "hog1", "upre")
    assert set(genotype.promoter_drivers).isdisjoint(PROMOTER_UNCOUPLED_SIGNALS)
    drivers = {name: float(reference.trace(f"driver.{name}")[-1]) for name in SIGNALS}
    off = {**drivers, **dict.fromkeys(PROMOTER_UNCOUPLED_SIGNALS, 0.0)}
    on = {**drivers, **dict.fromkeys(PROMOTER_UNCOUPLED_SIGNALS, 1.0)}
    for gene in genotype.genes:
        assert gene.promoter.activity(off) == gene.promoter.activity(on)
        assert gene.promoter.activity(drivers) == gene.promoter.activity(on)
    coupled = Promoter(prior("basal", 0.0, "dimensionless", "test"),
                       {"snf1": prior("snf1", 1.0, "dimensionless", "test")})
    assert coupled.activity(on) == 1.0 and coupled.activity(off) == 0.0


def test_the_carbon_promoter_refusal_names_its_measurement_and_matches_the_inventory():
    inventory = load_parameter_evidence()
    rows = [row for row in inventory.to_dict()["parameters"]
            if row["programme"] in ("carbon_pka_snf1", "nitrogen_tor")]
    assert tuple(row["id"] for row in rows) == _CARBON_ROWS
    for row in rows:
        assert not sourced_promoter_response(row)
        assert row["units_status"] == "unspecified_in_source"
        assert row["uncertainty"]["kind"] == "not_reported"
        assert row["identifiability"] == "unresolved"
    for programme in ("carbon_pka_snf1", "nitrogen_tor"):
        assert inventory.programme(programme)["blockers"]
    assert sourced_promoter_response({"units_status": "source_declared",
                                      "uncertainty": {"kind": "confidence_interval"}})
    assert not sourced_promoter_response({"units_status": "source_declared",
                                          "uncertainty": {"kind": "not_reported"}})
    with pytest.raises(ValueError, match="mapping"):
        sourced_promoter_response("carbon.pka.w_pka_camp")


def test_the_carbon_promoter_refusal_quotes_the_inventory_rows_it_rejected():
    rows = [row for row in load_parameter_evidence().to_dict()["parameters"]
            if row["programme"] == "carbon_pka_snf1"]
    assert len(PROMOTER_COUPLING_NOT_BUILT) == len(PROMOTER_UNCOUPLED_SIGNALS)
    for signal in PROMOTER_UNCOUPLED_SIGNALS:
        refusal = PROMOTER_COUPLING_NOT_BUILT[f"promoter_response.{signal}"]
        assert refusal.tag == Tag.REFUSED
        with pytest.raises(RefusedValue, match="interval"):
            float(refusal)
        assert "occupancy" in refusal.units and "interval" in refusal.missing
        assert any(f"{row['id']} = {row['value']}" in refusal.source for row in rows)
        assert JALIHAL_SBML_SHA256 in refusal.source


def test_the_carbon_refusal_reads_the_pinned_source_model_it_cites():
    inventory = load_parameter_evidence().to_dict()
    digests = {asset["path"]: asset["sha256"]
               for asset in inventory["sources"]["jalihal2021"]["local_artifacts"]}
    source = {}
    for name in ("variables.txt", "parameters.txt"):
        relative = f"data/native_reference_models/jalihal2021/{name}"
        payload = source_path(repository_root(), relative).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == digests[relative]
        source[name] = payload.decode("utf-8")
    equations = dict(line.split("\t", 1) for line in source["variables.txt"].splitlines() if line.strip())
    constants = dict(line.split("\t", 1) for line in source["parameters.txt"].splitlines() if line.strip())
    assert len(equations) == 25 and {"PKA", "Snf1", "Mig1", "Dot6"} <= set(equations)
    assert not [name for name, rhs in equations.items()
                if name != "Mig1" and re.search(r"\bMig1\b", rhs)]
    assert re.search(r"\bSch9\s*\*\s*PKA\b", equations["Dot6"])
    audited = {row["id"].rsplit(".", 1)[-1] for row in inventory["parameters"]
               if row["source_id"] == "jalihal2021"}
    refusals = [PROMOTER_COUPLING_NOT_BUILT[f"promoter_response.{signal}"]
                for signal in PROMOTER_UNCOUPLED_SIGNALS]
    for name in ("k_transcription", "k_mRNA_degr", "w_gln1_gln3"):
        assert all(name in refusal.reason and constants[name] in refusal.reason
                   for refusal in refusals)
        assert name not in audited


def test_the_nitrogen_promoter_coupling_is_a_declared_prior_that_moves_the_product():
    genotype = Genotype.prior()
    driven = [gene for gene in genotype.genes if "torc1" in gene.promoter.responses]
    assert [gene.name for gene in driven] == ["crte", "crti", "crtyb"]
    for gene in driven:
        response = gene.promoter.responses["torc1"]
        assert response.tag == Tag.ASSERTED
        assert "PRIOR" in response.source and response.missing.strip()
    ladder = [run_case(genotype=replace(genotype, signal_activity={"torc1": gain}))
              for gain in (0.0, 0.25)] + [run_case(genotype=genotype)]
    for name in ("signal.torc1", "promoter.crte", "gene.crte.mrna", "flux.crte.mrna_synthesis"):
        values = [result.trace(name)[-1] for result in ladder]
        assert all(a < b for a, b in zip(values, values[1:]))
    assert ladder[0].trace("promoter.crte")[-1] == pytest.approx(float(driven[0].promoter.basal))


def test_the_carbon_and_nitrogen_states_are_graded_dose_responses_not_switches():
    glucose = [run_case(medium={**MEDIUM, "glucose": level}) for level in (1.0, 25.0, 100.0)]
    nitrogen = [run_case(medium={**MEDIUM, "nitrogen": level}) for level in (1.0, 5.0, 60.0)]
    snf1 = [result.trace("signal.snf1")[-1] for result in glucose]
    pka = [result.trace("signal.pka")[-1] for result in glucose]
    torc1 = [result.trace("signal.torc1")[-1] for result in nitrogen]
    assert all(a > b for a, b in zip(snf1, snf1[1:]))
    assert all(a < b for a, b in zip(pka, pka[1:]))
    assert all(a < b for a, b in zip(torc1, torc1[1:]))
    for result in (*glucose, *nitrogen):
        for name in ("pka", "snf1", "torc1"):
            trace = result.trace(f"signal.{name}")
            assert np.all(trace >= 0.0) and np.all(trace <= 1.0)


def test_the_carbon_and_nitrogen_states_relax_on_their_declared_time_constants():
    times = (0.0, 0.02, 0.05, 0.25, 1.0)
    fast = run_case(times=times)
    slow = run_case(times=times,
                    parameters=parameters_with(pka_tau=0.3, snf1_tau=1.0, torc1_tau=1.0))
    for name in ("pka", "snf1", "torc1"):
        quick, lagged = fast.trace(f"signal.{name}"), slow.trace(f"signal.{name}")
        assert quick[0] == pytest.approx(lagged[0], abs=1e-12)
        assert np.all(lagged[1:] < quick[1:])
        assert quick[-1] > quick[1]


def test_source_hog_adapter_reuses_native_rates_and_cancels_absolute_protein_units():
    adapter = NativeHogResponse()
    ref = adapter.reference
    fraction = ref["Hog1PP"] / adapter.total
    rates = adapter.rates(0, fraction=fraction, volume_ratio=1, osmolarity_M=ref["OsmoE"],
                          glycerol_external_M=ref["glycerol_e"], glycerol_internal_M=ref["glycerol_i"],
                          temperature_c=ref["vV_T"] - 273.15)
    native = adapter.model.reaction_rates(0, adapter.initial, {"parameter_97": 0})
    scale = 3600 / (ref["intra"] * adapter.total)
    assert rates[0] == pytest.approx(native[adapter.reactions["v16f"]] * scale)
    assert rates[1] == pytest.approx(native[adapter.reactions["v16r"]] * scale)
    derivative = adapter.model.rhs(0, adapter.initial, {"parameter_97": 0})
    assert rates[2] == pytest.approx(3600 * derivative[adapter.indices["cellvol"]] / ref["cellvol"])
    more_glycerol = adapter.rates(0, fraction=fraction, volume_ratio=1, osmolarity_M=ref["OsmoE"],
                                 glycerol_external_M=ref["glycerol_e"], glycerol_internal_M=ref["glycerol_i"] + 0.1,
                                 temperature_c=ref["vV_T"] - 273.15)
    assert more_glycerol[2] > rates[2]


def test_hog_gpd1_glycerol_feedback_and_carbon_cost_are_one_connected_system():
    genotype = Genotype.prior()
    extra_gpd1 = genotype.with_copies(gpd1=8)
    control = run_case(genotype=genotype, medium={**MEDIUM, "osmolyte": 650}, times=(0, 1, 2))
    changed = run_case(genotype=extra_gpd1, medium={**MEDIUM, "osmolyte": 650}, times=(0, 1, 2))
    assert changed.trace("internal.glycerol")[-1] > control.trace("internal.glycerol")[-1]
    assert changed.trace("signal.hog1")[-1] < control.trace("signal.hog1")[-1]
    assert changed.trace("cell_volume_ratio")[-1] > control.trace("cell_volume_ratio")[-1]
    assert changed.trace("flux.glycerol_synthesis")[-1] > control.trace("flux.glycerol_synthesis")[-1]
    np.testing.assert_allclose(changed.trace("balance.carbon_residual"), 0, atol=3e-9)


def test_ph_transport_preserves_acid_and_buffer_charge_and_changes_promoter_growth():
    g = Genotype.prior()
    reporter = next(x for x in g.genes if x.name == "reporter")
    responsive = replace(reporter, promoter=Promoter(prior("basal", 0, "dimensionless", "test"),
                                                    {"acid": prior("acid", 1, "dimensionless", "test")}))
    g = replace(g, genes=tuple(responsive if x.name == "reporter" else x for x in g.genes))
    low = run_case(genotype=g, medium={**MEDIUM, "acetate": 10}, controls=(replace(AERATED, ph=4),))
    high = run_case(genotype=g, medium={**MEDIUM, "acetate": 10}, controls=(replace(AERATED, ph=6),))
    assert low.trace("ph_c")[-1] < high.trace("ph_c")[-1]
    assert low.trace("growth_per_h")[-1] < high.trace("growth_per_h")[-1]
    assert low.trace("promoter.reporter")[-1] > high.trace("promoter.reporter")[-1]
    p = EngineParameters.prior()._resolved
    for result in (low, high):
        acid = result.trace("content.acetate") / p["cell_water_l_gdw"]
        residuals = [ph.charge_imbalance_mM(a, x, p["buffer_capacity"], float(ph.PKA_ACETIC),
                                           resting_ph=p["resting_ph"])
                     for a, x in zip(acid, result.trace("ph_c"), strict=True)]
        np.testing.assert_allclose(residuals, 0, atol=1e-8)
        np.testing.assert_allclose(result.trace("external.acetate") + result.trace("internal.acetate")
                                   + result.trace("dead.acetate"), 10, atol=1e-8)


def test_atp_and_pyridine_cofactors_are_conserved_pools_with_separate_growth_synthesis(reference):
    p = EngineParameters.prior()._resolved
    for pool, initial in ((("atp", "adp", "amp", "camp"), p["initial_atp"] + p["initial_adp"] + p["initial_amp"]),
                          (("nad", "nadh"), p["initial_nad"] + p["initial_nadh"]),
                          (("nadp", "nadph"), p["initial_nadp"] + p["initial_nadph"])):
        total = sum(reference.trace(f"{prefix}.{name}") for prefix in ("internal", "dead") for name in pool)
        np.testing.assert_allclose(total - total[0], initial * reference.trace("ledger.biomass_formed"), atol=1e-10)
    atp_change = reference.trace("internal.atp") + reference.trace("dead.atp") - reference.trace("internal.atp")[0]
    np.testing.assert_allclose(atp_change, reference.trace("ledger.atp_produced") - reference.trace("ledger.atp_consumed"), atol=1e-9)
    assert reference.trace("ledger.nadph_consumed")[-1] > 0
    assert reference.trace("flux.respiration")[-1] > 0
    assert reference.trace("flux.fermentation")[-1] > 0
    assert reference.trace("flux.ppp")[-1] > 0
    assert reference.trace("internal.glycogen")[-1] > 0


def test_anoxia_blocks_oxygen_costs_and_desaturation_not_all_carbon_flux(reference):
    result = run_case(medium={**MEDIUM, "oxygen": 0}, controls=(Control(0),))
    assert np.all(result.trace("flux.respiration") == 0)
    assert np.all(result.trace("flux.crti") == 0)
    assert result.trace("internal.beta_carotene")[-1] == 0
    assert result.trace("internal.phytoene")[-1] > 0
    assert result.trace("external.ethanol")[-1] > 0
    assert result.trace("flux.fermentation")[-1] > 0
    assert reference.trace("internal.beta_carotene")[-1] > 0


@pytest.mark.parametrize("gene,accumulates", [("crte", "ipp"), ("crti", "phytoene"), ("crtyb", "ggpp")])
def test_product_bottlenecks_are_genotype_specific_not_a_flat_ceiling(gene, accumulates, reference):
    knockout = run_case(genotype=Genotype.prior().with_copies(**{gene: 0}))
    assert knockout.trace("internal.beta_carotene")[-1] == 0
    assert knockout.trace(f"internal.{accumulates}")[-1] > reference.trace(f"internal.{accumulates}")[-1]
    assert knockout.trace(f"gene.{gene}.protein")[-1] == 0
    assert knockout.trace(f"gene.{gene}.mrna")[-1] == 0


def test_native_competition_and_domain_specific_activity_change_product(reference):
    less_native = run_case(parameters=parameters_with(native_isoprenoid=0))
    cyclase_dead = run_case(parameters=parameters_with(lcy_kcat=0))
    overexpression = run_case(genotype=Genotype.prior().with_copies(crte=4, crti=4, crtyb=4))
    assert less_native.trace("internal.beta_carotene")[-1] > reference.trace("internal.beta_carotene")[-1]
    assert less_native.trace("internal.native_isoprenoid")[-1] == 0
    assert cyclase_dead.trace("internal.beta_carotene")[-1] == 0
    assert cyclase_dead.trace("internal.lycopene")[-1] > reference.trace("internal.lycopene")[-1]
    assert overexpression.trace("internal.beta_carotene")[-1] > reference.trace("internal.beta_carotene")[-1]


def test_expression_costs_are_atom_and_atp_debits_and_copies_act_before_protein():
    p, g = EngineParameters.prior(), Genotype.prior()
    state = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM=MEDIUM)
    base = _Kernel(p, g)
    more = _Kernel(p, g.with_copies(crte=4))
    y = base.pack(state)
    _, baseline = base.evaluate(0, y, AERATED)
    _, increased = more.evaluate(0, y, AERATED)
    assert increased["flux.crte.mrna_synthesis"] == pytest.approx(4 * baseline["flux.crte.mrna_synthesis"])
    assert increased["flux.crte.protein_synthesis"] == baseline["flux.crte.protein_synthesis"]
    assert increased["flux.crte"] == baseline["flux.crte"]
    col = base.reaction_names.index("crte.protein_synthesis")
    assert base.matrix[base.index["internal.atp"], col] < 0
    assert base.matrix[base.index["internal.nadph"], col] < 0
    assert base.matrix[base.index["external.nitrogen"], col] < 0
    assert base.matrix[base.index["internal.acetyl"], col] < 0
    costs = base.matrix[base.index["internal.atp"]]
    assert costs[base.reaction_names.index("crte.folding")] < 0


def test_carotenoid_desaturation_has_existing_oxidase_ros_stoichiometry(reference):
    names = reference.diagnostics["stoichiometric_state_names"]
    reactions = reference.diagnostics["reaction_names"]
    matrix = reference.diagnostics["stoichiometric_matrix"]
    col = reactions.index("crti")
    assert matrix[names.index("external.oxygen"), col] == -4
    assert matrix[names.index("internal.ros"), col] == 4
    assert matrix[names.index("internal.nadph"), col] == 0


def test_observation_gain_maturation_and_missingness_cannot_feed_back_into_truth(reference):
    clear = ObservationModel.prior(seed=3)
    missing = replace(clear, values={**clear.values,
                                    "rfu_gain": prior("gain", 2e9, "RFU L/mmol", "test"),
                                    "missing_probability": prior("missing", 1, "dimensionless", "test")})
    measured = run_case(observation=clear)
    absent = run_case(observation=missing)
    for name in reference.truth:
        np.testing.assert_array_equal(reference.trace(name), measured.trace(name))
        np.testing.assert_array_equal(reference.trace(name), absent.trace(name))
    assert np.isfinite(measured.observations.values["rfu"]).all()
    assert np.isnan(absent.observations.values["rfu"]).all()
    assert absent.observations.missing["rfu"].all()
    mature = measured.observations.metadata["mature_mmol"]
    assert mature[0] == 0
    assert mature[-1] > 0
    assert np.all(mature <= measured.trace("gene.reporter.protein"))


def test_optical_gain_and_ph_quench_apply_to_reporter_not_blank_or_truth():
    obs = ObservationModel.prior()
    low = run_case(observation=obs, medium={**MEDIUM, "acetate": 5})
    double = replace(obs, values={**obs.values, "rfu_gain": prior("gain", 2e9, "RFU L/mmol", "test")})
    high = run_case(observation=double, medium={**MEDIUM, "acetate": 5})
    auto = float(obs.values["autofluorescence"]) * low.trace("biomass_gdw") / low.trace("volume_l")
    attenuation = np.exp(-float(obs.values["inner_filter"]) * low.trace("content.beta_carotene"))
    blank = float(obs.values["rfu_blank"])
    np.testing.assert_allclose(high.observations.values["rfu"] - blank - auto * attenuation,
                               2 * (low.observations.values["rfu"] - blank - auto * attenuation), atol=1e-9)
    np.testing.assert_array_equal(low.trace("ph_c"), high.trace("ph_c"))


def test_observation_maturation_respects_sampling_event_and_oxygen_dependency():
    obs = ObservationModel.prior()
    sampled = run_case(observation=obs, events=(Event("sample", 0.5, withdraw_l=0.5),), times=(0, 0.499999, 0.5, 1))
    mature = sampled.observations.metadata["mature_mmol"]
    assert mature[2] == pytest.approx(mature[1] / 2, rel=1e-5)
    anoxic = run_case(observation=obs, medium={**MEDIUM, "oxygen": 0}, controls=(Control(0),))
    assert np.all(anoxic.observations.metadata["mature_mmol"] == 0)


def test_thermodynamic_audit_uses_existing_gate_and_keeps_unknown_blocked_and_uncovered_distinct():
    from ystwin.bridge.thermodynamic import ThermodynamicData
    from ystwin.pathway.thermo_gate import Feasibility, Unresolved

    p, g = EngineParameters.prior(), Genotype.prior()
    state = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM=MEDIUM)
    state = replace(state, intracellular_mmol={**state.intracellular_mmol, "lycopene": 1e-5, "beta_carotene": 1e-5})
    blocked = replace(p, thermodynamics=ThermodynamicRequest(
        ThermodynamicData({"lycopene_c": 0, "betacarotene_c": 100}, 1.0), "synthetic uphill fixture"))
    result = run_case(parameters=blocked, state=state, times=(0, 0.05))
    report = result.diagnostics["thermodynamic_reports"][0]
    assert report.feasibility == Feasibility.CANNOT_RUN
    assert report.dg_kj_per_mol == pytest.approx(100)
    assert result.validity.thermodynamic_covered == ("lcy",)
    assert "crti" in result.validity.thermodynamic_uncovered
    unknown = replace(p, thermodynamics=ThermodynamicRequest(ThermodynamicData({}, 0), "missing fixture"))
    unknown_result = run_case(parameters=unknown, state=state, times=(0, 0.05))
    assert unknown_result.diagnostics["thermodynamic_reports"][0].feasibility == Feasibility.CANNOT_SAY
    assert unknown_result.diagnostics["thermodynamic_reports"][0].reason == Unresolved.NO_ENERGY
    assert not result.validity.biological_validation


def test_solver_convergence_on_product_and_fast_signals_with_event_boundaries():
    events = (Event("salt", 0.313, add_mmol={"osmolyte": 300}),)
    coarse = run_case(events=events)
    tight = run_case(events=events, rtol=1e-8, max_step_h=0.025)
    independent = run_case(events=events, rtol=1e-8, max_step_h=0.025, method="Radau")
    for name in ("biomass_gdw", "internal.atp", "internal.beta_carotene", "signal.hog1", "upr.hac1_protein"):
        np.testing.assert_allclose(coarse.trace(name), tight.trace(name), rtol=5e-4, atol=2e-9)
        np.testing.assert_allclose(tight.trace(name), independent.trace(name), rtol=2e-5, atol=2e-10)


def test_starved_native_regulatory_tail_records_roundoff_without_moving_any_inventory():
    result = run_case(genotype=Genotype(()), medium={**MEDIUM, "glucose": 5, "nitrogen": 0},
                      times=(0, 1, 24, 120), max_step_h=1.0)
    assert result.diagnostics["amount_projection"] == "none"
    assert all(value >= 0 for value in result.final_state.upr.values())
    for name in ("carbon", "nitrogen", "live_biomass"):
        np.testing.assert_allclose(result.trace(f"balance.{name}_residual"), 0, atol=1e-8)


def test_an_initially_unreachable_oxygen_ros_branch_wakes_at_the_declared_supply_event():
    result = run_case(medium={**MEDIUM, "oxygen": 0},
                      controls=(Control(0), replace(AERATED, time_h=0.4)),
                      times=(0, 0.2, 0.4, 0.6, 1))
    assert np.all(result.trace("internal.ros")[:3] == 0)
    assert np.all(result.trace("internal.beta_carotene")[:3] == 0)
    assert result.trace("internal.ros")[-1] > 0
    assert result.trace("internal.beta_carotene")[-1] > 0
    assert result.trace("flux.respiration")[-1] > 0


def test_storage_is_remobilized_when_the_shared_carbon_signal_falls():
    p, g = EngineParameters.prior(), Genotype.prior()
    state = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM={**MEDIUM, "glucose": 0})
    state = replace(state, intracellular_mmol={**state.intracellular_mmol, "glycogen": 0.1})
    result = run_case(parameters=p, genotype=g, state=state)
    assert result.trace("flux.glycogen_mobilization")[-1] > 0
    assert result.trace("internal.glycogen")[-1] < 0.1
    assert result.trace("flux.lower_glycolysis")[-1] > 0
    np.testing.assert_allclose(result.trace("balance.carbon_residual"), 0, atol=1e-9)


def test_hsp70_feedback_and_inherited_antioxidant_change_recovery_not_the_input_protocol(reference):
    hsp = run_case(genotype=Genotype.prior().with_copies(hsp70=8))
    assert hsp.trace("signal.hsf1")[-1] < reference.trace("signal.hsf1")[-1]
    assert hsp.trace("unfolded_fraction")[-1] < reference.trace("unfolded_fraction")[-1]
    p, g = EngineParameters.prior(), Genotype.prior()
    state = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM={**MEDIUM, "peroxide": 0.02})
    protein = state.expression["antioxidant"]
    memory = replace(state, expression={**state.expression, "antioxidant": replace(
        protein, protein_mmol=8 * protein.protein_mmol, active_mmol=8 * protein.active_mmol)})
    naive = run_case(parameters=p, genotype=g, state=state)
    protected = run_case(parameters=p, genotype=g, state=memory)
    assert protected.trace("internal.ros")[1] < naive.trace("internal.ros")[1]
    assert protected.trace("gene.antioxidant.active")[1] > naive.trace("gene.antioxidant.active")[1]
    assert protected.trace("external.peroxide")[0] == naive.trace("external.peroxide")[0]


def test_optical_noise_is_seeded_and_reads_do_not_change_maturation_integration():
    model = ObservationModel.prior(seed=27)
    model = replace(model, values={**model.values,
                                   "rfu_noise_sd": prior("noise", 1, "RFU", "test"),
                                   "missing_probability": prior("missing", 0.5, "dimensionless", "test")})
    first = run_case(observation=model)
    repeated = run_case(observation=model)
    np.testing.assert_equal(first.observations.values["rfu"], repeated.observations.values["rfu"])
    clear = ObservationModel.prior()
    sparse = run_case(observation=clear, times=(0, 0.5, 1))
    dense = run_case(observation=clear, times=np.linspace(0, 1, 11))
    np.testing.assert_allclose(sparse.observations.metadata["mature_mmol"],
                               dense.observations.metadata["mature_mmol"][[0, 5, 10]], rtol=1e-12, atol=1e-20)


def test_existing_proton_bill_and_dynamic_pump_share_the_stoichiometric_cost_not_a_growth_penalty():
    bill = ph.proton_bill(0.01, 0.2, 2.0, float(ph.PKA_ACETIC))
    flux = bill.trapped_anion_mM * 0.2 * 2.0 / 1000.0
    assert bill.mmol_atp_per_gdcw_h == ph.proton_pumping_atp(flux, atp_per_proton=1.0)
    assert ph.proton_pumping_atp(flux, atp_per_proton=2.0) == pytest.approx(2 * bill.mmol_atp_per_gdcw_h)
    for value in (float("nan"), -1.0):
        with pytest.raises(ValueError):
            ph.proton_pumping_atp(value, atp_per_proton=1.0)


def test_paired_acid_export_pumps_protons_recovers_ph_and_pays_once_from_the_atp_inventory():
    export = AcidExportParameters.prior()
    active_parameters = replace(EngineParameters.prior(), acid_export=export)
    inactive_parameters = replace(active_parameters, acid_export=replace(export, vmax=prior(
        "disabled", 0, "mmol/gDW/h", "test intervention")))
    active = run_case(parameters=active_parameters, medium={**MEDIUM, "acetate": 5},
                      controls=(replace(AERATED, ph=4),))
    inactive = run_case(parameters=inactive_parameters, medium={**MEDIUM, "acetate": 5},
                        controls=(replace(AERATED, ph=4),))
    assert active.trace("ph_c")[-1] > inactive.trace("ph_c")[-1]
    assert active.trace("internal.acetate")[-1] < inactive.trace("internal.acetate")[-1]
    assert active.trace("flux.acid_export")[-1] > 0
    assert active.trace("ph_pump_atp_mmol_h")[-1] > 0
    cost = float(export.atp_per_proton) + float(export.atp_per_anion)
    np.testing.assert_allclose(active.trace("ph_pump_atp_mmol_h") + active.trace("anion_export_atp_mmol_h"),
                               cost * active.trace("flux.acid_export"))
    names = active.diagnostics["stoichiometric_state_names"]
    column = active.diagnostics["reaction_names"].index("acid_export")
    assert active.diagnostics["stoichiometric_matrix"][names.index("internal.atp"), column] == -cost
    np.testing.assert_allclose(active.trace("ledger.protons_pumped"), active.trace("ledger.acetate_exported"))
    np.testing.assert_allclose(active.trace("external.acetate") + active.trace("internal.acetate")
                               + active.trace("dead.acetate"), 5, atol=1e-8)
    atp_change = active.trace("internal.atp") + active.trace("dead.atp") - active.trace("internal.atp")[0]
    np.testing.assert_allclose(atp_change, active.trace("ledger.atp_produced") - active.trace("ledger.atp_consumed"), atol=1e-9)
    for name in ("carbon", "nitrogen"):
        np.testing.assert_allclose(active.trace(f"balance.{name}_residual"), 0, atol=1e-9)
    assert active.diagnostics["ph_model"] == "paired_acid_export_prior"


def test_acid_export_has_no_flux_or_atp_bill_without_substrate_or_atp():
    p = replace(EngineParameters.prior(), acid_export=AcidExportParameters.prior())
    free = run_case(parameters=p)
    assert np.all(free.trace("flux.acid_export") == 0)
    g = Genotype.prior()
    initial = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM={**MEDIUM, "acetate": 5})
    pools = {**initial.intracellular_mmol, **dict.fromkeys(("atp", "adp", "amp", "camp"), 0.0)}
    powerless = run_case(parameters=p, genotype=g, state=replace(initial, intracellular_mmol=pools),
                         times=(0, 0.01, 0.03))
    assert np.all(powerless.trace("flux.acid_export") == 0)
    assert np.all(powerless.trace("ph_pump_atp_mmol_h") == 0)


def test_fast_partition_has_one_owner_and_cofactor_pools_stay_dynamic(reference):
    partition = reference.diagnostics["pool_partition"]
    assert "internal.atp" in partition.retained_fast_pool_candidates
    assert "internal.nadph" in partition.integrated
    assert "internal.glycogen" in partition.integrated
    assert "ph_c" in partition.fast_algebraic
    assert not set(partition.integrated).intersection(partition.fast_algebraic)
    acid = run_case(medium={**MEDIUM, "acetate": 5}, times=(0, 0.001, 0.01, 0.1))
    np.testing.assert_allclose(acid.trace("acid.intracellular_neutral_mmol")
                               + acid.trace("acid.intracellular_anion_mmol"), acid.trace("internal.acetate"), atol=1e-12)
    np.testing.assert_allclose(acid.trace("acid.buffered_protons_mmol"),
                               acid.trace("acid.intracellular_anion_mmol"), atol=1e-10)
    assert not np.allclose(acid.trace("internal.atp"), acid.trace("internal.atp")[0], rtol=1e-3)


def thermodynamic_parameters(mode, energies, *, reactions=("lcy",), metabolite_ids=None, ph_value=7.08):
    from ystwin.bridge.thermodynamic import ThermodynamicData

    options = {} if metabolite_ids is None else {"metabolite_ids": metabolite_ids}
    species = ("lycopene", "beta_carotene") if reactions == ("lcy",) else ("atp", "adp", "amp")
    coefficients = {name: prior(name, 1.0, "dimensionless", "explicit ideal accessible-pool activity test prior") for name in species}
    return replace(EngineParameters.prior(), thermodynamics=ThermodynamicRequest(
        ThermodynamicData(energies, 1.0, ph=ph_value), "synthetic energy fixture, not a published value",
        mode=mode, reactions=reactions, activity_coefficients=coefficients, **options))


def test_thermodynamic_constraint_blocks_uphill_lcy_during_integration_not_after_the_result():
    g = Genotype.prior()
    energies = {"lycopene_c": 0, "betacarotene_c": 100}
    constrained = thermodynamic_parameters("constrain_supported", energies)
    audit = thermodynamic_parameters("audit", energies)
    state = initialize(constrained, g, volume_l=1, biomass_gdw_l=0.2, medium_mM=MEDIUM)
    state = replace(state, intracellular_mmol={**state.intracellular_mmol, "lycopene": 1e-5, "beta_carotene": 1e-5})
    blocked = run_case(parameters=constrained, state=state, genotype=g, times=(0, 0.05, 0.1))
    unconstrained = run_case(parameters=audit, state=state, genotype=g, times=(0, 0.05, 0.1))
    assert np.all(blocked.trace("flux.lcy") == 0)
    assert unconstrained.trace("flux.lcy")[-1] > 0
    assert blocked.trace("internal.beta_carotene")[-1] < unconstrained.trace("internal.beta_carotene")[-1]
    assert blocked.diagnostics["thermodynamic_constrained_reactions"] == ("lcy",)
    assert unconstrained.diagnostics["thermodynamic_constrained_reactions"] == ()
    np.testing.assert_allclose(blocked.trace("balance.carbon_residual"), 0, atol=1e-9)


def test_adenylate_kinase_direction_constraint_uses_dynamic_atp_adp_amp_not_external_signals():
    ids = {"atp": "ATP", "amp": "AMP", "adp": "ADP"}
    p = thermodynamic_parameters("constrain_supported", {"ATP": 0, "AMP": 0, "ADP": 50},
                                  reactions=("adenylate_kinase",), metabolite_ids=ids)
    g = Genotype.prior()
    initial = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM=MEDIUM)
    initial = replace(initial, intracellular_mmol={**initial.intracellular_mmol, "amp": 0.001, "adp": 1e-8})
    result = run_case(parameters=p, genotype=g, state=initial, times=(0, 0.01, 0.03))
    assert result.trace("flux.adenylate_kinase")[0] == 0
    assert result.diagnostics["thermodynamic_constrained_reactions"] == ("adenylate_kinase",)
    for name in ("carbon", "nitrogen"):
        np.testing.assert_allclose(result.trace(f"balance.{name}_residual"), 0, atol=1e-9)
    assert "internal.atp" in result.diagnostics["pool_partition"].integrated


def test_thermodynamic_constraints_refuse_unknown_energies_and_mismatched_conditions():
    missing = thermodynamic_parameters("constrain_supported", {})
    with pytest.raises(ScientificRefusal, match="energy"):
        run_case(parameters=missing, times=(0, 0.01))
    ids = {"atp": "ATP", "amp": "AMP", "adp": "ADP"}
    wrong_ph = thermodynamic_parameters("constrain_supported", {"ATP": 0, "AMP": 0, "ADP": 0},
                                        reactions=("adenylate_kinase",), metabolite_ids=ids, ph_value=6.0)
    with pytest.raises(ScientificRefusal, match="pH"):
        run_case(parameters=wrong_ph, times=(0, 0.01))
    positive = thermodynamic_parameters("constrain_supported", {"lycopene_c": 0, "betacarotene_c": -10})
    with pytest.raises(ScientificRefusal, match="temperature"):
        run_case(parameters=positive, controls=(replace(AERATED, temperature_c=37),), times=(0, 0.01))


def test_zero_product_thermodynamic_limit_has_no_pseudocount_or_invented_inventory():
    p = thermodynamic_parameters("constrain_supported", {"lycopene_c": 0, "betacarotene_c": -10})
    result = run_case(parameters=p, times=(0, 0.05, 0.1))
    assert result.trace("internal.beta_carotene")[0] == 0
    assert result.trace("internal.beta_carotene")[-1] > 0
    assert result.diagnostics["amount_projection"] == "none"
    np.testing.assert_allclose(result.trace("balance.carbon_residual"), 0, atol=1e-9)


def test_refuted_thermodynamic_estimates_cannot_be_revived_by_constraint_mode_or_zero_products():
    from ystwin.bridge.thermodynamic import RefutedEnergy

    p = thermodynamic_parameters("constrain_supported", {"lycopene_c": 0, "betacarotene_c": -10})
    request = p.thermodynamics
    data = replace(request.data, refuted={"betacarotene_c": RefutedEnergy(
        "fixture", "beta carotene", "inconsistent energy estimates", "synthetic refusal fixture")})
    p = replace(p, thermodynamics=replace(request, data=data))
    with pytest.raises(ScientificRefusal, match="refuted_uncertainty"):
        run_case(parameters=p, times=(0, 0.01))


def test_reverse_adenylate_flux_is_also_blocked_if_its_direction_is_uphill():
    ids = {"atp": "ATP", "amp": "AMP", "adp": "ADP"}
    p = thermodynamic_parameters("constrain_supported", {"ATP": 0, "AMP": 0, "ADP": -50},
                                  reactions=("adenylate_kinase",), metabolite_ids=ids)
    g = Genotype.prior()
    state = initialize(p, g, volume_l=1, biomass_gdw_l=0.2, medium_mM=MEDIUM)
    state = replace(state, intracellular_mmol={**state.intracellular_mmol, "atp": 1e-8, "amp": 1e-8, "adp": 8e-4})
    constrained = _Kernel(p, g)
    unconstrained = _Kernel(replace(p, thermodynamics=replace(p.thermodynamics, mode="audit")), g)
    vector = constrained.pack(state)
    assert unconstrained.evaluate(0, vector, AERATED)[1]["flux.adenylate_kinase"] < 0
    assert constrained.evaluate(0, vector, AERATED)[1]["flux.adenylate_kinase"] == 0


def test_partial_thermodynamic_constraints_report_each_reaction_on_the_single_clock():
    p = thermodynamic_parameters("constrain_supported", {"lycopene_c": 0, "betacarotene_c": -10})
    result = run_case(parameters=p, times=(0, 0.01, 0.03))
    assert result.diagnostics["thermodynamic_report_coordinates"] == ((0.0, "lcy"), (0.01, "lcy"), (0.03, "lcy"))
    assert "pdh" in result.validity.thermodynamic_uncovered
    assert "thermodynamic_activity.lycopene" in result.provenance
    assert not result.validity.biological_validation
