"""Williamson's cAMP block as an INSTRUMENT against the engine's, under one glucose step.

WHAT THIS FILE IS FOR. ``mech/carbon_williamson2009.py`` is a transcription nothing consumes,
and ``tests/test_carbon_williamson2009_transcription.py`` only ever asks whether it reproduces
its own source: it never imports the engine. This file is the other half, the one
``tests/test_cell_wall_axis_registration.py`` established for Talemi -- drive BOTH blocks under
the same published protocol and record what the comparison says about OUR engine. It says
something negative, and the negative is the deliverable.

THE PROTOCOL IS WILLIAMSON'S OWN, taken from :data:`carbon_williamson2009.GLUCOSE_PULSE_PROTOCOL`
rather than retyped: 5 mM glucose to a glucose-starved suspension at 60 s, 100 mM at 240 s, 600 s
of run. Williamson gets it as two ``Glucout`` resets; the engine gets it as two typed
``contracts.Event`` boluses into 1 L, which land it at 5.000 mM and 99.965 mM -- the 0.035 mM gap
is the glucose the engine consumed in between, and it is the whole difference between the two
deliveries. Time enters the engine in hours and Williamson in seconds, converted once at
``SECONDS_PER_HOUR``.

THE RESULT, MEASURED. Williamson answers the first pulse with a peak 3.016x its basal cAMP,
47.9 s after the step, which then falls 69.6% of the way back to basal by the time the second
pulse arrives; the second pulse gives a second peak 1.294x the first, 37.1 s after the step,
which has given up 50.2% of its excursion by 600 s. The engine's cAMP does none of it. It is
monotone non-decreasing at every one of the 601 samples, its maximum inside each pulse window is
the LAST sample of that window rather than an interior peak, and the fraction of any excursion
returned toward baseline is exactly 0.0. Run six times longer it is still monotone; run at five
times the cell density it is still monotone. The engine has no transient at all, so peak height
and time to peak are not measurable on it -- only censored at the end of the run, at 3.478e-4
mmol/gDW and a PKA fraction of 0.740.

WHY, ALSO MEASURED, NOT ASSERTED. The engine destroys cAMP first-order at ``camp_hydrolysis``,
and the effective clearance constant recomputed from the trajectory is 10.000/h at every sample,
in both the prior run and a PKA-null counterfactual. Williamson's effective clearance over the
same protocol runs from 213.9/h to 418.0/h, a 1.955x span, as the phosphorylated fraction of
Pde1 rises from 0.031 to 0.187 against a kcat that goes 1.113 -> 25.254 /s on phosphorylation.
The second PKA loop is running too: phosphorylated Cdc25 goes 0.321 -> 0.689 and active Ras2,
having peaked at 0.193, is back to 0.105 while glucose is still 100 mM. Twelve of Williamson's
seventeen states -- the whole G-protein input and both feedback loops -- match nothing anywhere
in the engine's trajectory, which carries the input (Glucout, Glucin), the pool (cAMP) and the
output (PKA, as a fraction) and nothing in between. Two consequences: the engine's cAMP pool turns over in 360 s, which is 7.5x longer than the
whole 47.9 s in which Williamson rises to its first peak, and the engine's only PKA-to-cAMP path
is worth 0.24%. Deleting the pka signal entirely cuts glucose transport 3.96x and internal
glucose 10.5x and moves cAMP at 600 s by a factor 1.0024, without changing the shape or the
clearance constant by anything.

WHAT THIS CANNOT SETTLE, SAID BEFORE THE NUMBERS ARE USED. The unit bridge is refused at both
ends -- ``williamson.cell_volume_litres`` and ``williamson.absolute_abundance_scale`` -- so
Williamson's mM and the engine's mmol/gDW are not convertible and no absolute level is compared
here. Every comparison in this file is either a SHAPE (transient vs monotone, adaptation vs
none) or a per-time quantity (a first-order clearance constant, a time to peak), both of which
survive the missing volume. The one pair that is directly comparable in kind, the dimensionless
PKA fractions, ends within 1.232x -- and that closeness is NOT evidence: the engine's ``camp_k``
is an asserted prior and every Williamson constant is fitted, which is why
:func:`carbon_williamson2009.gate` refuses. Nothing here licenses wiring the axis; it stays
registered and undriven, and the two promoter refusals it carries still raise.

No tolerance below was chosen before its number was computed, and every number is recomputed
from the two models on each run rather than being read back from a stored artifact.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ystwin.mech import carbon_williamson2009 as wm
from ystwin.mech import contracts, engine, signalling
from ystwin.mech.params import RefusedValue, Tag

PROTOCOL_END_S = 600.0
SAMPLE_STEP_S = 1.0
FIRST_PULSE_S, SECOND_PULSE_S = (time for time, _ in wm.GLUCOSE_PULSE_PROTOCOL)
FIRST_LEVEL_MM, SECOND_LEVEL_MM = (level for _, level in wm.GLUCOSE_PULSE_PROTOCOL)
STARVED_MEDIUM = {"glucose": 0.0, "nitrogen": 20.0, "oxygen": 0.2}
VOLUME_L = 1.0
BIOMASS_GDW_L = 0.2

# The five Williamson states this engine already carries, and what carries them.
ALREADY_CARRIED = {"Glucout": "external.glucose", "Glucin": "internal.glucose",
                   "cAMP": "internal.camp", "PKAi": "signal.pka", "C": "signal.pka"}
NO_ENGINE_COUNTERPART = ("Gpr1", "Gpr1Glucout", "Gpa2i", "Gpa2a", "Krh", "Gpa2aKrh", "Ras2i",
                         "Ras2a", "Cdc25", "Cdc25P", "Pde1", "Pde1P")


def _engine_run(*, genotype=None, duration_s=PROTOCOL_END_S, step_s=SAMPLE_STEP_S,
                biomass_gdw_l=BIOMASS_GDW_L) -> dict:
    """Williamson's pulse protocol on the engine, through the ordinary public entrypoint.

    The two pulses are boluses into ``VOLUME_L`` sized to Williamson's own levels, applied at
    his own times converted once from seconds to hours.
    """
    genotype = contracts.Genotype.prior() if genotype is None else genotype
    times_s = np.arange(0.0, duration_s + 0.5 * step_s, step_s)
    events = (
        contracts.Event("williamson_pulse_1", FIRST_PULSE_S / wm.SECONDS_PER_HOUR,
                        add_mmol={"glucose": FIRST_LEVEL_MM * VOLUME_L}),
        contracts.Event("williamson_pulse_2", SECOND_PULSE_S / wm.SECONDS_PER_HOUR,
                        add_mmol={"glucose": (SECOND_LEVEL_MM - FIRST_LEVEL_MM) * VOLUME_L}),
    )
    protocol = contracts.Protocol(
        tuple(times_s / wm.SECONDS_PER_HOUR),
        (contracts.Control(0.0, oxygen_transfer_per_h=30.0, oxygen_saturation_mM=0.2),),
        events=events,
        source="Williamson 2009 Figure 8 glucose pulse, driven on this engine for comparison")
    parameters = engine.EngineParameters.prior()
    initial = engine.initialize(parameters, genotype, volume_l=VOLUME_L,
                                biomass_gdw_l=biomass_gdw_l, medium_mM=STARVED_MEDIUM)
    result = engine.simulate(protocol, genotype, initial, parameters)
    return {"time_s": times_s, "truth": result.truth, "variables": result.variables,
            "validity": result.validity, "parameters": parameters}


def _shape(time_s, values) -> dict:
    """The shape metrics both blocks are scored on. No units enter, so no bridge is needed.

    ``recovered`` is the fraction of a window's excursion given back by its end; a block with
    no interior peak scores 0.0 and reports the window edge as its argmax.
    """
    values = np.asarray(values, dtype=float)
    first = int(np.argmin(np.abs(time_s - FIRST_PULSE_S)))
    second = int(np.argmin(np.abs(time_s - SECOND_PULSE_S)))
    peak_one = first + int(np.argmax(values[first:second + 1]))
    peak_two = second + int(np.argmax(values[second:]))
    baseline = float(values[first])
    def recovered(peak_index, end_index):
        span = float(values[peak_index]) - baseline
        return 0.0 if span <= 0.0 else (float(values[peak_index]) - float(values[end_index])) / span
    return {
        "baseline": baseline,
        "peak_one": float(values[peak_one]), "t_peak_one_s": float(time_s[peak_one]),
        "at_second_pulse": float(values[second]),
        "peak_two": float(values[peak_two]), "t_peak_two_s": float(time_s[peak_two]),
        "final": float(values[-1]),
        "recovered_one": recovered(peak_one, second),
        "recovered_two": recovered(peak_two, len(values) - 1),
        "argmax_s": float(time_s[int(np.argmax(values))]),
        "monotone": bool(np.all(np.diff(values) >= 0.0)),
        "still_rising_at_the_end": bool(np.diff(values)[-1] > 0.0),
    }


@pytest.fixture(scope="module")
def engine_run():
    return _engine_run()


@pytest.fixture(scope="module")
def engine_hour():
    """The same protocol, integrated six times longer, to ask whether it EVER turns over."""
    return _engine_run(duration_s=6.0 * PROTOCOL_END_S, step_s=5.0)


@pytest.fixture(scope="module")
def engine_dense():
    """Five times the cell density: the one way the engine could peak is substrate depletion."""
    return _engine_run(biomass_gdw_l=5.0 * BIOMASS_GDW_L)


@pytest.fixture(scope="module")
def engine_pka_null():
    """pka held at zero through the typed genotype, not by editing anything."""
    genotype = replace(contracts.Genotype.prior(), signal_activity={"pka": 0.0},
                       source="pka-null counterfactual for the cAMP feedback measurement")
    return _engine_run(genotype=genotype)


@pytest.fixture(scope="module")
def williamson():
    """Williamson's own block under its own protocol, with the rate laws recomputed per sample."""
    model = wm.WilliamsonCAMP()
    run = model.simulate(duration_h=PROTOCOL_END_S / wm.SECONDS_PER_HOUR)
    species = run["species"]
    stacked = np.array([species[name] for name in wm.DYNAMIC_SPECIES])
    sampled, clearance, phosphorylated = [], [], []
    for column in range(0, stacked.shape[1], 10):
        rates = model.reaction_rates(stacked[:, column])
        hydrolysis = (rates["cAMPhydroPde1"] + rates["cAMPhydroPde1p"] + rates["cAMPhydroPde2"])
        sampled.append(run["time_s"][column])
        clearance.append(hydrolysis / species["cAMP"][column] * wm.SECONDS_PER_HOUR)
        phosphorylated.append(species["Pde1P"][column]
                              / (species["Pde1"][column] + species["Pde1P"][column]))
    return {
        "time_s": run["time_s"], "camp_mM": species["cAMP"], "species": species,
        "pka_fraction": species["C"] / (species["C"] + 2.0 * species["PKAi"]),
        "rate_time_s": np.array(sampled),
        "clearance_per_h": np.array(clearance),
        "pde1p_fraction": np.array(phosphorylated),
    }


class TestBothBlocksSeeTheSameGlucoseStep:
    """Before any comparison: the two really are being driven by one protocol."""

    def test_the_pulses_are_read_off_the_module_not_retyped_here(self):
        assert wm.GLUCOSE_PULSE_PROTOCOL == ((60.0, 5.0), (240.0, 100.0))
        assert (FIRST_PULSE_S, SECOND_PULSE_S) == (60.0, 240.0)
        assert (FIRST_LEVEL_MM, SECOND_LEVEL_MM) == (5.0, 100.0)

    def test_the_engine_is_driven_in_hours_and_williamson_in_seconds(self, engine_run):
        assert wm.NATIVE_TIME_UNIT == "second" and wm.SECONDS_PER_HOUR == 3600.0
        assert engine_run["time_s"][-1] == pytest.approx(PROTOCOL_END_S)
        assert len(engine_run["time_s"]) == 601

    def test_the_engine_reaches_williamsons_two_glucose_levels(self, engine_run):
        glucose = np.asarray(engine_run["truth"]["medium_mM.glucose"])
        assert glucose[0] == 0.0 and glucose[59] == 0.0
        assert glucose[60] == pytest.approx(FIRST_LEVEL_MM, rel=1e-12)
        assert glucose[240] == pytest.approx(99.9655, rel=1e-4)
        # The only difference between a bolus and Williamson's reset: 0.035 mM consumed by 240 s.
        assert SECOND_LEVEL_MM - glucose[240] == pytest.approx(0.0345, abs=5e-3)


class TestTheEngineProducesNoTransient:
    """The recorded negative. Nothing here is a tolerance on an interior peak, because there is none."""

    def test_the_engines_camp_is_monotone_through_both_pulses(self, engine_run):
        camp = _shape(engine_run["time_s"], engine_run["truth"]["content.camp"])
        assert camp["monotone"] is True
        assert camp["still_rising_at_the_end"] is True
        assert camp["argmax_s"] == PROTOCOL_END_S

    def test_every_window_maximum_is_the_window_edge_rather_than_a_peak(self, engine_run):
        camp = _shape(engine_run["time_s"], engine_run["truth"]["content.camp"])
        assert camp["t_peak_one_s"] == SECOND_PULSE_S
        assert camp["t_peak_two_s"] == PROTOCOL_END_S
        assert camp["peak_one"] == camp["at_second_pulse"]
        assert camp["peak_two"] == camp["final"]

    def test_nothing_is_returned_toward_baseline_in_either_window(self, engine_run):
        camp = _shape(engine_run["time_s"], engine_run["truth"]["content.camp"])
        pka = _shape(engine_run["time_s"], engine_run["truth"]["signal.pka"])
        # Non-degenerate: both windows rose, so a zero here is no recovery rather than no excursion.
        assert camp["peak_one"] > camp["baseline"] and camp["peak_two"] > camp["baseline"]
        assert camp["recovered_one"] == 0.0 and camp["recovered_two"] == 0.0
        assert pka["recovered_one"] == 0.0 and pka["recovered_two"] == 0.0

    def test_the_censored_levels_are_recorded_since_no_peak_can_be(self, engine_run):
        camp = np.asarray(engine_run["truth"]["content.camp"])
        pka = np.asarray(engine_run["truth"]["signal.pka"])
        assert camp[240] == pytest.approx(1.5374e-4, rel=2e-3)
        assert camp[-1] == pytest.approx(3.4781e-4, rel=2e-3)
        assert pka[240] == pytest.approx(0.39184, rel=2e-3)
        assert pka[-1] == pytest.approx(0.74047, rel=2e-3)

    def test_the_starved_rest_state_carries_no_camp_and_no_pka_at_all(self, engine_run):
        camp = np.asarray(engine_run["truth"]["content.camp"])
        pka = np.asarray(engine_run["truth"]["signal.pka"])
        # Exactly zero, not small: cAMP synthesis is strictly proportional to the glucose signal.
        assert camp[0] == 0.0 and camp[60] == 0.0
        assert pka[0] == 0.0 and pka[60] == 0.0

    def test_a_fold_change_is_undefined_on_the_engine_and_defined_on_williamson(self, engine_run,
                                                                                williamson):
        camp = _shape(engine_run["time_s"], engine_run["truth"]["content.camp"])
        native = _shape(williamson["time_s"], williamson["camp_mM"])
        assert camp["baseline"] == 0.0
        assert native["baseline"] == pytest.approx(8.2592e-4, rel=2e-3)
        assert native["peak_one"] / native["baseline"] == pytest.approx(3.016, rel=5e-3)

    def test_it_is_still_monotone_six_times_longer(self, engine_hour, engine_run):
        camp = _shape(engine_hour["time_s"], engine_hour["truth"]["content.camp"])
        assert camp["monotone"] is True
        assert camp["argmax_s"] == 6.0 * PROTOCOL_END_S
        assert camp["final"] == pytest.approx(4.5619e-4, rel=2e-3)
        # Williamson's whole protocol ends with the engine at 76% of the plateau it is climbing to.
        at_600 = float(np.asarray(engine_run["truth"]["content.camp"])[-1])
        assert at_600 / camp["final"] == pytest.approx(0.7624, rel=5e-3)

    def test_it_is_still_monotone_at_five_times_the_cell_density(self, engine_dense):
        camp = _shape(engine_dense["time_s"], engine_dense["truth"]["content.camp"])
        glucose = np.asarray(engine_dense["truth"]["medium_mM.glucose"])
        assert camp["monotone"] is True and camp["argmax_s"] == PROTOCOL_END_S
        # Depletion is the only route to a peak this form has, and it does not open: 1% consumed.
        assert glucose[-1] == pytest.approx(99.005, rel=1e-3)
        assert camp["final"] == pytest.approx(3.4768e-4, rel=2e-3)


class TestWilliamsonProducesOneUnderTheSameStep:
    """The instrument side: what the engine is being scored against, in its own units."""

    def test_the_first_pulse_gives_a_peak_that_comes_back_down(self, williamson):
        camp = _shape(williamson["time_s"], williamson["camp_mM"])
        assert camp["peak_one"] == pytest.approx(2.4910e-3, rel=2e-3)
        assert camp["t_peak_one_s"] - FIRST_PULSE_S == pytest.approx(47.9, abs=1.0)
        assert camp["recovered_one"] == pytest.approx(0.696, rel=1e-2)
        assert camp["monotone"] is False

    def test_the_second_pulse_gives_a_second_peak_that_also_comes_back_down(self, williamson):
        camp = _shape(williamson["time_s"], williamson["camp_mM"])
        assert camp["peak_two"] == pytest.approx(3.2226e-3, rel=2e-3)
        assert camp["t_peak_two_s"] - SECOND_PULSE_S == pytest.approx(37.1, abs=1.0)
        assert camp["peak_two"] / camp["peak_one"] == pytest.approx(1.294, rel=5e-3)
        assert camp["recovered_two"] == pytest.approx(0.502, rel=1e-2)

    def test_the_dimensionless_pka_fraction_overshoots_twice(self, williamson):
        fraction = _shape(williamson["time_s"], williamson["pka_fraction"])
        assert fraction["baseline"] == pytest.approx(0.11799, rel=2e-3)
        assert fraction["peak_one"] == pytest.approx(0.56814, rel=2e-3)
        assert fraction["t_peak_one_s"] == pytest.approx(131.1, abs=1.0)
        assert fraction["peak_two"] == pytest.approx(0.84252, rel=2e-3)
        assert fraction["t_peak_two_s"] == pytest.approx(295.2, abs=1.0)
        assert fraction["final"] == pytest.approx(0.60095, rel=2e-3)


class TestWhyTheEngineCannot:
    """The structural difference, measured off both trajectories rather than read off the source."""

    def test_the_engines_camp_clearance_is_one_constant_and_williamsons_is_not(self, engine_run,
                                                                               williamson):
        internal = np.asarray(engine_run["truth"]["internal.camp"])
        hydrolysis = np.asarray(engine_run["truth"]["flux.camp_hydrolysis"])
        present = internal > 0.0
        clearance = hydrolysis[present] / internal[present]
        declared = float(engine_run["parameters"].values["camp_hydrolysis"])
        assert declared == 10.0
        assert clearance.min() == pytest.approx(declared, rel=1e-9)
        assert clearance.max() == pytest.approx(declared, rel=1e-9)
        native = williamson["clearance_per_h"]
        assert native.min() == pytest.approx(213.9, rel=5e-3)
        assert native.max() == pytest.approx(418.0, rel=5e-3)
        assert native.max() / native.min() == pytest.approx(1.955, rel=5e-3)

    def test_the_engines_camp_pool_turns_over_slower_than_williamsons_whole_rise(self, williamson,
                                                                                 engine_run):
        declared = float(engine_run["parameters"].values["camp_hydrolysis"])
        turnover_s = wm.SECONDS_PER_HOUR / declared
        camp = _shape(williamson["time_s"], williamson["camp_mM"])
        assert turnover_s == pytest.approx(360.0)
        assert turnover_s / (camp["t_peak_one_s"] - FIRST_PULSE_S) == pytest.approx(7.5, rel=2e-2)
        assert declared / williamson["clearance_per_h"].max() == pytest.approx(0.0239, rel=1e-2)

    def test_the_feedback_williamson_has_is_pde1_phosphorylation(self, williamson):
        phosphorylated = williamson["pde1p_fraction"]
        basal = int(np.argmin(np.abs(williamson["rate_time_s"] - FIRST_PULSE_S)))
        assert phosphorylated[basal] == pytest.approx(0.03097, rel=5e-3)
        assert phosphorylated.max() == pytest.approx(0.1866, rel=5e-3)
        kcat = float(wm.CAMP_HYDROLYSIS_PDE1P_KCAT) / float(wm.CAMP_HYDROLYSIS_PDE1_KCAT)
        assert kcat == pytest.approx(22.69, rel=5e-3)

    def test_the_second_feedback_loop_throttles_ras_through_cdc25(self, williamson):
        species, times = williamson["species"], williamson["time_s"]
        basal = int(np.argmin(np.abs(times - FIRST_PULSE_S)))
        cdc25p = species["Cdc25P"] / (species["Cdc25"] + species["Cdc25P"])
        active_ras = species["Ras2a"] / (species["Ras2a"] + species["Ras2i"])
        assert cdc25p[basal] == pytest.approx(0.3214, rel=5e-3)
        assert cdc25p[-1] == pytest.approx(0.6895, rel=5e-3)
        # Active Ras2 peaks at 0.193 and is back to 0.105 while glucose is still 100 mM.
        assert active_ras.max() == pytest.approx(0.1930, rel=5e-3)
        assert active_ras[-1] == pytest.approx(0.1051, rel=5e-3)

    @pytest.mark.parametrize("state", NO_ENGINE_COUNTERPART)
    def test_twelve_of_the_seventeen_states_have_no_engine_counterpart(self, state, engine_run):
        assert state in wm.DYNAMIC_SPECIES
        assert not [name for name in engine_run["truth"] if state.lower() in name.lower()]

    @pytest.mark.parametrize("state,carried_by", sorted(ALREADY_CARRIED.items()))
    def test_the_other_five_are_the_input_the_output_and_nothing_between(self, state, carried_by,
                                                                        engine_run):
        assert state in wm.DYNAMIC_SPECIES
        assert carried_by in engine_run["truth"]
        assert set(NO_ENGINE_COUNTERPART) | set(ALREADY_CARRIED) == set(wm.DYNAMIC_SPECIES)

    def test_deleting_pka_from_the_engine_moves_camp_by_a_quarter_of_a_percent(self, engine_run,
                                                                               engine_pka_null):
        prior, null = engine_run["truth"], engine_pka_null["truth"]
        assert float(np.asarray(null["signal.pka"])[-1]) == 0.0
        transport = (float(np.asarray(prior["flux.glucose_transport"])[-1])
                     / float(np.asarray(null["flux.glucose_transport"])[-1]))
        internal_glucose = (float(np.asarray(prior["content.glucose"])[-1])
                            / float(np.asarray(null["content.glucose"])[-1]))
        camp = (float(np.asarray(null["content.camp"])[-1])
                / float(np.asarray(prior["content.camp"])[-1]))
        assert transport == pytest.approx(3.958, rel=5e-3)
        assert internal_glucose == pytest.approx(10.48, rel=5e-3)
        # A 3.96x change upstream is worth 0.24% here, and PKA reaches cAMP by no other path.
        assert camp == pytest.approx(1.0024, rel=5e-4)

    def test_the_pka_null_run_has_the_same_shape_and_the_same_clearance(self, engine_pka_null):
        camp = _shape(engine_pka_null["time_s"], engine_pka_null["truth"]["content.camp"])
        internal = np.asarray(engine_pka_null["truth"]["internal.camp"])
        hydrolysis = np.asarray(engine_pka_null["truth"]["flux.camp_hydrolysis"])
        present = internal > 0.0
        assert camp["monotone"] is True and camp["argmax_s"] == PROTOCOL_END_S
        assert camp["recovered_one"] == 0.0 and camp["recovered_two"] == 0.0
        assert (hydrolysis[present] / internal[present]).max() == pytest.approx(10.0, rel=1e-9)

    def test_the_engines_answer_to_the_twenty_fold_step_is_one_michaelis_ratio(self, engine_run):
        truth = engine_run["truth"]
        glucose = np.asarray(truth["medium_mM.glucose"])
        synthesis = np.asarray(truth["flux.camp_synthesis"])
        km = float(engine_run["parameters"].values["glucose_k"])
        before, after = glucose[239], glucose[241]
        saturating = (after / (after + km)) / (before / (before + km))
        assert after / before == pytest.approx(20.13, rel=1e-3)
        assert synthesis[241] / synthesis[239] == pytest.approx(1.1894, rel=1e-3)
        assert saturating == pytest.approx(synthesis[241] / synthesis[239], rel=1e-3)


class TestWhatThisComparisonCannotSettle:
    """Every limit the numbers above have, as a test rather than as a caveat in prose."""

    @pytest.mark.parametrize("name", ["williamson.cell_volume_litres",
                                      "williamson.absolute_abundance_scale"])
    def test_the_unit_bridge_is_refused_at_both_ends(self, name):
        param = wm.NOT_IN_SOURCE[name]
        assert param.tag == Tag.REFUSED and param.value is None
        with pytest.raises(RefusedValue):
            float(param)

    def test_no_absolute_level_is_compared_because_the_bases_differ(self, engine_run):
        engine_units = engine_run["variables"]["content.camp"].units
        axis = signalling.TRANSCRIBED_AXES["carbon_camp_input"]
        native_units = {name: units for name, units, _, _ in axis.states}["cAMP"]
        assert engine_units == "mmol/gDW" and native_units == "mM"
        # The only thing that would convert one to the other is the volume Williamson omits.
        with pytest.raises(RefusedValue, match="compartment"):
            wm.cell_volume_litres()

    def test_the_two_endpoint_pka_fractions_are_close_and_that_is_not_evidence(self, engine_run,
                                                                               williamson):
        engine_fraction = float(np.asarray(engine_run["truth"]["signal.pka"])[-1])
        native_fraction = float(williamson["pka_fraction"][-1])
        assert engine_fraction / native_fraction == pytest.approx(1.232, rel=5e-3)
        assert engine_run["parameters"].values["camp_k"].tag == Tag.ASSERTED
        gate = wm.gate()
        assert not gate.passes and gate.targets == ()

    def test_the_paths_differ_even_where_the_endpoints_agree(self, engine_run, williamson):
        engine_fraction = np.asarray(engine_run["truth"]["signal.pka"])
        native = _shape(williamson["time_s"], williamson["pka_fraction"])
        assert engine_fraction[240] / native["at_second_pulse"] == pytest.approx(1.119, rel=5e-3)
        # The engine never visits either overshoot it would have had to pass through.
        assert engine_fraction.max() < native["peak_two"]
        assert float(engine_fraction[131]) < native["peak_one"]

    def test_the_negative_does_not_license_wiring_the_axis(self):
        axis = signalling.TRANSCRIBED_AXES["carbon_camp_input"]
        assert axis.module == "mech/carbon_williamson2009.py"
        assert axis.readers == () and axis.driven is False
        assert signalling.TRANSCRIBED_AXIS_FLUX_COUPLING["carbon_camp_input"] == ()
        assert axis.panel_promoter == "CSRE"
        for qualified in axis.promoter_refusals:
            with pytest.raises(Exception):
                float(contracts.refused_row(qualified))

    def test_the_engine_result_still_reports_the_axis_registered_and_undriven(self, engine_run):
        line = [s for s in engine_run["validity"].unsupported if "TRANSCRIBED_AXES" in s]
        assert len(line) == 1
        assert "carbon_camp_input" in line[0] and "REGISTERED AND UNDRIVEN" in line[0]
