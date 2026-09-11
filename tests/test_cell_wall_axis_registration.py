"""The cell_wall_slt2 registration, and the measurement that keeps it undriven.

``mech/signalling.py``'s registration for this axis does not merely assert that wiring it
would double-count the osmotic response; it quotes numbers. This module recomputes every one
of them, so the claim cannot rot into a story about a run nobody can reproduce. The engine
side is run through the ordinary public entrypoint at the ordinary prior genotype, which is
the only configuration the claim is about.
"""

from __future__ import annotations

import numpy as np
import pytest

from ystwin.mech import cell_wall_talemi2016 as talemi
from ystwin.mech import contracts, engine, signalling


STEP_MM = 750.0
BASAL_OSMOLYTE_MM = 250.0
SAMPLE_TIMES_H = (0.0, 0.05, 0.1, 0.5, 6.0)
MEDIUM = {"glucose": 50.0, "nitrogen": 20.0, "oxygen": 0.2}

# The five Talemi states this engine already integrates, and what integrates them.
ALREADY_INTEGRATED = {
    "Vos": "cell_volume_ratio",
    "Hog1": "signal.hog1",
    "Hog1PP": "signal.hog1",
    "Glyin": "internal.glycerol",
    "Glyex": "external.glycerol",
}
ONLY_HERE = ("HOGSignal", "Slt2Signal", "Slt2", "Slt2P", "Fps1", "Fps1P", "Sensitizer")


@pytest.fixture(scope="module")
def axis():
    return signalling.TRANSCRIBED_AXES["cell_wall_slt2"]


@pytest.fixture(scope="module")
def engine_step():
    """The engine's own osmotic block under a +750 mM step, via the public entrypoint."""
    parameters = engine.EngineParameters.prior()
    genotype = contracts.Genotype.prior()
    protocol = contracts.Protocol(
        SAMPLE_TIMES_H,
        (contracts.Control(0.0, oxygen_transfer_per_h=30.0, oxygen_saturation_mM=0.2),),
        source="cell_wall_slt2 double-count measurement")
    medium = dict(MEDIUM, osmolyte=BASAL_OSMOLYTE_MM + STEP_MM)
    initial = engine.initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=0.2,
                                medium_mM=medium)
    result = engine.simulate(protocol, genotype, initial, parameters)
    return {"hog1": np.asarray(result.truth["signal.hog1"]),
            "volume": np.asarray(result.truth["cell_volume_ratio"]),
            "truth": result.truth, "validity": result.validity}


@pytest.fixture(scope="module")
def talemi_step():
    """Talemi's own osmotic block under the same step, from its own deposited rest state."""
    model = talemi.OsmoStat()
    shock = talemi.Shock(s1_uM=STEP_MM * 1e3, s2_uM=STEP_MM * 1e3, t_off_s=0.0, t_s_s=1.0e6)
    run = model.timecourse(hours=max(SAMPLE_TIMES_H), shock=shock, points=601)
    index = [int(np.argmin(np.abs(run["time_h"] - t))) for t in SAMPLE_TIMES_H]
    # FitHog1PPrel divides by fn*(Hog1+Hog1PP); undoing fn puts it on the engine's pool basis.
    pool = np.asarray(run["hog1pp_relative_percent"]) / 100.0 * float(talemi.F_NUC)
    return {"hog1": pool[index], "volume": np.asarray(run["volume_relative"])[index],
            "slt2p": np.asarray(run["slt2p_relative_percent"])[index]}


class TestTheRegistrationIsShapedLikeTheOthers:
    def test_the_axis_is_registered_under_its_own_name(self, axis):
        assert axis.axis == "cell_wall_slt2"
        assert axis.module == "mech/cell_wall_talemi2016.py"
        assert axis.states == talemi.STATE_VARIABLES

    def test_every_state_is_namespaced_typed_and_sourced(self, axis):
        assert axis.state_names == tuple(f"axis.cell_wall_slt2.{n}" for n in talemi.STATE_NAMES)
        for name, units, compartment, meaning in axis.states:
            assert isinstance(compartment, contracts.Compartment)
            assert units.strip() and meaning.strip() and name.strip()

    def test_the_engine_exposes_them_on_a_basis_that_cannot_integrate(self, axis):
        exposed = {name: var for name, var in engine.REGISTERED_AXIS_VARIABLES.items()
                   if name.startswith("axis.cell_wall_slt2.")}
        assert set(exposed) == set(axis.state_names)
        assert {var.basis for var in exposed.values()} == {"registered, not integrated"}

    def test_the_registered_names_do_not_collide_with_any_integrated_state(self, engine_step):
        assert not set(engine.REGISTERED_AXIS_VARIABLES) & set(engine_step["truth"])

    def test_it_is_undriven_here_and_in_the_flux_map(self, axis):
        assert axis.readers == () and axis.driven is False
        assert signalling.TRANSCRIBED_AXIS_FLUX_COUPLING["cell_wall_slt2"] == ()

    @pytest.mark.parametrize("qualified", [
        "AXIS_PROMOTER_NOT_BUILT::promoter_response.rlm1_box",
        "MEDIUM_SPECIES_NOT_CARRIED::medium.cell_wall_damaging_agent",
    ])
    def test_each_cited_refusal_resolves_and_raises_on_use(self, qualified):
        row = contracts.refused_row(qualified)
        with pytest.raises(Exception):
            float(row)

    def test_the_declared_promoter_is_backed_by_a_refused_coefficient(self, axis):
        assert axis.panel_promoter == "RLM1 box"
        assert axis.promoter_refusals == ("AXIS_PROMOTER_NOT_BUILT::promoter_response.rlm1_box",)

    def test_the_claimed_missing_medium_species_really_is_absent(self, axis):
        assert axis.missing_medium == (
            "MEDIUM_SPECIES_NOT_CARRIED::medium.cell_wall_damaging_agent",)
        assert "cell_wall_damaging_agent" not in contracts.EXTRACELLULAR

    def test_the_engine_result_reports_it_as_registered_and_undriven(self, engine_step):
        line = [s for s in engine_step["validity"].unsupported if "TRANSCRIBED_AXES" in s]
        assert len(line) == 1
        assert "cell_wall_slt2" in line[0] and "REGISTERED AND UNDRIVEN" in line[0]


class TestTheDoubleCountIsMeasuredNotAsserted:
    """The input exists and the output consumer exists; this is why it stays unwired."""

    def test_the_engine_already_answers_an_osmolyte_dose(self):
        parameters = engine.EngineParameters.prior()
        genotype = contracts.Genotype.prior()
        protocol = contracts.Protocol(
            (0.0, 0.05, 0.1), (contracts.Control(0.0, oxygen_transfer_per_h=30.0,
                                                 oxygen_saturation_mM=0.2),), source="ladder")
        peaks = []
        for osmolyte in (BASAL_OSMOLYTE_MM, BASAL_OSMOLYTE_MM + STEP_MM):
            initial = engine.initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=0.2,
                                        medium_mM=dict(MEDIUM, osmolyte=osmolyte))
            truth = engine.simulate(protocol, genotype, initial, parameters).truth
            peaks.append((float(np.max(truth["signal.hog1"])),
                          float(np.min(truth["cell_volume_ratio"]))))
        unstressed, stepped = peaks
        assert unstressed[0] < 0.1 < stepped[0]
        assert stepped[1] < 0.6 < unstressed[1]

    @pytest.mark.parametrize("state,integrated_by", sorted(ALREADY_INTEGRATED.items()))
    def test_five_of_the_twelve_states_are_already_integrated(self, state, integrated_by,
                                                              engine_step):
        assert state in talemi.STATE_NAMES
        assert integrated_by in engine_step["truth"]

    @pytest.mark.parametrize("state", ONLY_HERE)
    def test_the_other_seven_have_no_engine_counterpart(self, state, engine_step):
        assert state in talemi.STATE_NAMES
        assert not [n for n in engine_step["truth"] if state.lower() in n.lower()]

    def test_both_blocks_agree_that_hog1_saturates_on_the_same_step(self, engine_step,
                                                                    talemi_step):
        assert engine_step["hog1"].max() == pytest.approx(0.890, abs=0.02)
        assert talemi_step["hog1"].max() == pytest.approx(0.781, abs=0.02)
        assert abs(engine_step["hog1"].max() - talemi_step["hog1"].max()) < 0.15

    def test_and_disagree_permanently_about_the_adaptation(self, engine_step, talemi_step):
        engine_final, talemi_final = float(engine_step["hog1"][-1]), float(talemi_step["hog1"][-1])
        assert engine_final == pytest.approx(0.760, abs=0.02)
        assert talemi_final == pytest.approx(0.134, abs=0.02)
        assert engine_final / talemi_final > 4.0

    def test_the_two_cell_volumes_disagree_by_more_than_an_order_of_magnitude(self, engine_step,
                                                                              talemi_step):
        engine_excursion = 1.0 - float(engine_step["volume"][-1])
        talemi_excursion = 1.0 - float(talemi_step["volume"][-1])
        assert float(engine_step["volume"][-1]) == pytest.approx(0.595, abs=0.02)
        assert float(talemi_step["volume"][-1]) == pytest.approx(0.968, abs=0.01)
        assert engine_excursion / talemi_excursion > 10.0

    def test_talemi_adapts_and_the_engine_does_not_within_the_run(self, engine_step,
                                                                  talemi_step):
        recovered = ((talemi_step["volume"][-1] - talemi_step["volume"].min())
                     / (1.0 - talemi_step["volume"].min()))
        stuck = ((engine_step["volume"][-1] - engine_step["volume"].min())
                 / (1.0 - engine_step["volume"].min()))
        assert recovered == pytest.approx(0.889, abs=0.02)
        assert stuck == pytest.approx(0.191, abs=0.02)
        assert recovered / stuck > 4.0


class TestTheVolumeBridgeIsPricedNotJustRefused:
    """Wiring only the Slt2 arm avoids the duplicated states and costs this much instead."""

    def test_the_arm_reads_a_different_number_off_each_models_volume(self, engine_step,
                                                                     talemi_step):
        geometry = talemi._derived_geometry()
        engine_ratio = float(engine_step["volume"][-1])
        talemi_vos = (float(talemi_step["volume"][-1]) * float(talemi.V0)
                      - geometry["Vb"]) / geometry["Vos0"]
        rows = talemi.volume_response(ratios=(engine_ratio, talemi_vos))["rows"]
        off_engine = rows[engine_ratio]["slt2p_relative_percent"]
        off_talemi = rows[talemi_vos]["slt2p_relative_percent"]
        assert off_engine == pytest.approx(2.75, abs=0.2)
        assert off_talemi == pytest.approx(14.80, abs=0.5)
        assert off_talemi / off_engine > 4.5

    def test_the_engines_volume_falls_outside_the_band_the_sign_was_checked_over(self,
                                                                                 engine_step):
        assert float(engine_step["volume"].min()) < 0.90

    def test_the_bridge_and_the_dose_are_both_refused_at_the_source(self):
        for name in ("osmotically_active_cell_volume_bridge", "external_osmolarity_basal"):
            with pytest.raises(Exception):
                float(talemi.NOT_DRIVEN[name])

    def test_the_axis_answers_neither_panel_stressor(self):
        assert talemi.PANEL_CHANNEL["answered"] == {}
        assert set(talemi.PANEL_CHANNEL["forfeited"]) == {"congo_red", "caffeine"}
