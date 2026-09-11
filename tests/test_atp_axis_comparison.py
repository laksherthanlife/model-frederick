"""The atp_glycolysis registration, and the two measurements that keep it undriven.

``mech/signalling.py``'s registration for this axis makes two claims that are not opinions:
that wiring Teusink's pool in is blocked at the units, and that wiring his ``KeqAK`` in
without inverting it would invert the adenylate charge response. Both are recomputed here
against the engine's ordinary public entrypoint at the ordinary prior genotype, so neither
can rot into a story about a run nobody can reproduce. ``atp_teusink2000`` is offered as the
published alternative to 19 engine ``_PRIORS`` rows -- 11 rate constants found from the
engine's own stoichiometry plus 8 named scalars, every one stamped "PRIOR, not a published
constant" and COUNTED below rather than quoted -- so this file is where it has to earn that,
exactly as ``tests/test_cell_wall_axis_registration.py`` prices the Talemi arm.

WHAT IS COMPARED, AND WHAT IS REFUSED. ``teusink.cytosol_l_per_gdw`` is REFUSED, so the
deposit's mM and the engine's mmol/gDW are NOT interconvertible and no test below converts
them. Every comparison here is dimensionless -- energy charge, ATP:ADP:AMP ratios, and the
sign of a net rate. The absolute pools (SUM_P = 4.1 mM against the engine's 0.0032 mmol/gDW)
are left uncompared on purpose; running them together through the engine's asserted
``cell_water_l_gdw`` is precisely the laundering the module refuses.

WHAT THE COMPARISON FOUND, all recomputed below and none of it asserted here:

  glucose, mM    engine charge    Teusink charge    engine ATP:ADP    Teusink ATP:ADP
  49.60          0.9724           0.7692            26.68             1.940
  19.61          0.9734           0.7411            27.61             1.693
   4.67          0.9782           0.6051            33.03             0.974

The disagreement is not a calibration offset. Over the one arm both models share, the
engine's charge RISES 0.6% while Teusink's FALLS 21.3%: the dose responses have opposite
sign, which is a larger finding than the level gap and is measured in
:class:`TestThePoolDisagrees`. The engine's own prior INITIAL condition, by contrast, sits
within 1.6% of the published charge and is abandoned by its own dynamics inside fifteen
minutes.

THE SIGN PROBLEM is sharper still, and :class:`TestTheKeqDirectionIsATrap` records it as a
trap rather than a claim. Adenylate kinase is exactly charge-neutral in both models, so the
orientation error CANNOT show up in the energy charge; it shows up in the split and in the
direction of the net flux. At every state the engine occupies, the correctly inverted
constant drives the reaction forward (ATP consumed) and the un-inverted one drives it
backward (ATP regenerated) -- opposite signs, everywhere, with the un-inverted magnitude
2.21x the correct one.

AND THE TRAP HIDES FROM THE OBVIOUS CHECK. Done in concentrations, the two orientations
agree in sign at every one of the same states, because the engine's pool sits at a
mass-action quotient of 7.24, outside the 0.45-to-2.22 band where they can differ. The
engine's rate law is not written in concentrations: its ``adenylate_kinase`` flux law is
in occupancy coordinates, where the same pool quotes 1.002. A concentration-space sanity check would
therefore pass a wiring that is inverted where it actually runs, which is why this file
computes the driving force in the coordinates the engine really uses.

That 7.24-against-1.002 gap is also a THIRD refusal, beyond the two the registration names:
there is no single number that is "the engine's adenylate kinase Keq" for ``KeqAK`` to
replace, because ``saturation`` is not proportional to its argument and the quotient
depends on which coordinates it is taken in. Bridging them needs ``energy_k``, an engine
PRIOR, on top of the already-refused ``cytosol_l_per_gdw``.

No tolerance below was chosen before its number was computed, and every number is recomputed
from both models on each run rather than read back from a stored artifact.
"""

from __future__ import annotations

import numpy as np
import pytest

from ystwin.mech import atp_teusink2000 as teusink
from ystwin.mech import contracts, engine, signalling
from ystwin.mech.signalling import saturation


SAMPLE_TIMES_H = (0.0, 0.25, 1.0, 2.0, 4.0, 8.0)
READ_INDEX = 1
"""t = 0.25 h: the engine's adenylate split has relaxed and the medium has barely moved."""

GLUCOSE_LADDER_mM = (50.0, 20.0, 5.0, 2.0)
MEDIUM = {"nitrogen": 20.0, "oxygen": 0.2}
PAIRED_RUNGS = (50.0, 20.0, 5.0)
"""The rungs whose engine residual is still inside the transcription's validated domain."""

KEQ_AK = float(teusink.KEQ_AK)
ENGINE_DIRECTION = {"internal.atp": -1, "internal.amp": -1, "internal.adp": 2}


def _charge(atp, adp, amp):
    """Atkinson's charge on whichever three-species pool is handed in, dimensionless."""
    return (atp + 0.5 * adp) / (atp + adp + amp)


def _driving_forces(atp, adp, amp, energy_k):
    """The engine's adenylate_kinase net-rate numerator under three equilibrium choices.

    Its sign alone fixes the direction of the reaction at a state, which is the whole of
    the claim being checked. Positive is forward, ATP + AMP -> 2 ADP, i.e. ATP consumed.
    """
    a, d, m = (saturation(atp, energy_k), saturation(adp, energy_k),
               saturation(amp, energy_k))
    return {
        # As the engine's adenylate_kinase flux law writes it: no Keq, i.e. an implicit one.
        "as_written": a * m - d * d,
        # KeqAK inverted for the engine's direction: Keq_forward = 1/0.45.
        "inverted": a * m - KEQ_AK * d * d,
        # The trap: 0.45 adopted as though it were the engine's own forward Keq.
        "un_inverted": a * m - d * d / KEQ_AK,
        "quotient": a * m / (d * d),
    }


def _concentration_forces(atp, adp, amp):
    """The same three numerators in mass-action concentrations, which is NOT the engine's
    law. Kept only to show that the obvious sanity check fails to see the trap."""
    return {"as_written": atp * amp - adp * adp,
            "inverted": atp * amp - KEQ_AK * adp * adp,
            "un_inverted": atp * amp - adp * adp / KEQ_AK,
            "quotient": atp * amp / (adp * adp)}


def _engine_run(glucose_mM):
    """One engine trajectory through the ordinary public entrypoint at the prior genotype."""
    parameters = engine.EngineParameters.prior()
    genotype = contracts.Genotype.prior()
    protocol = contracts.Protocol(
        SAMPLE_TIMES_H,
        (contracts.Control(0.0, oxygen_transfer_per_h=30.0, oxygen_saturation_mM=0.2),),
        source="atp_glycolysis prior-replacement measurement")
    initial = engine.initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=0.2,
                                medium_mM=dict(MEDIUM, glucose=glucose_mM))
    result = engine.simulate(protocol, genotype, initial, parameters)
    truth, energy_k = result.truth, float(parameters.values["energy_k"])
    atp, adp, amp = (np.asarray(truth[f"content.{n}"]) for n in ("atp", "adp", "amp"))
    camp = np.asarray(truth["content.camp"])
    j = READ_INDEX
    return {
        "residual_glucose_mM": float(truth["medium_mM.glucose"][j]),
        "atp": float(atp[j]), "adp": float(adp[j]), "amp": float(amp[j]),
        "camp_fraction": float(camp[j] / (atp[j] + adp[j] + amp[j] + camp[j])),
        "charge": float(_charge(atp[j], adp[j], amp[j])),
        "charge_trajectory": _charge(atp, adp, amp),
        "mass_action": float(atp[j] * amp[j] / adp[j] ** 2),
        "forces": _driving_forces(float(atp[j]), float(adp[j]), float(amp[j]), energy_k),
        "concentration_forces": _concentration_forces(float(atp[j]), float(adp[j]),
                                                      float(amp[j])),
        "flux_ak": float(truth["flux.adenylate_kinase"][j]),
        "truth": truth, "validity": result.validity,
    }


@pytest.fixture(scope="module")
def axis():
    return signalling.TRANSCRIBED_AXES["atp_glycolysis"]


@pytest.fixture(scope="module")
def engine_ladder():
    return {glucose: _engine_run(glucose) for glucose in GLUCOSE_LADDER_mM}


@pytest.fixture(scope="module")
def teusink_ladder(engine_ladder):
    """Teusink run at the engine's OWN residual glucose, not at a nominal dose."""
    return {glucose: teusink.steady_state(glco=engine_ladder[glucose]["residual_glucose_mM"])
            for glucose in PAIRED_RUNGS}


@pytest.fixture(scope="module")
def adenylate_priors():
    """Every engine constant that sets the adenylate pool, found from the stoichiometry."""
    parameters = engine.EngineParameters.prior()
    kernel = engine._Kernel(parameters, contracts.Genotype.prior())
    pool = ("internal.atp", "internal.adp", "internal.amp")
    reactions = [name for name, stoichiometry in kernel.stoichiometry.items()
                 if any(species in stoichiometry for species in pool)]
    constants = [name for name in reactions if name in engine._PRIORS]
    scalars = ["initial_atp", "initial_adp", "initial_amp", "energy_k", "maintenance_atp",
               "growth_atp", "respiratory_atp_yield", "camp_k"]
    stamp = "PRIOR, not a published constant"
    everything = constants + scalars
    return {"reactions": reactions, "reaction_constants": constants, "scalars": scalars,
            "stamped": [n for n in everything if parameters.values[n].source.startswith(stamp)],
            "unstamped": [n for n in everything
                          if not parameters.values[n].source.startswith(stamp)]}


@pytest.fixture(scope="module")
def engine_priors():
    """The engine's asserted initial adenylate condition, before any dynamics touch it."""
    values = engine.EngineParameters.prior().values
    atp, adp, amp = (float(values[f"initial_{n}"]) for n in ("atp", "adp", "amp"))
    return {"atp": atp, "adp": adp, "amp": amp, "charge": float(_charge(atp, adp, amp)),
            "mass_action": atp * amp / adp ** 2}


class TestTheRegistrationIsShapedLikeTheOthers:
    def test_the_axis_is_registered_under_its_own_name(self, axis):
        assert axis.axis == "atp_glycolysis"
        assert axis.module == "mech/atp_teusink2000.py"
        assert teusink.SOURCE_DIGEST in axis.source

    def test_every_state_is_namespaced_typed_and_sourced(self, axis):
        assert axis.state_names == tuple(f"axis.atp_glycolysis.{n}"
                                         for n in teusink.SPECIES_ORDER)
        for name, units, compartment, meaning in axis.states:
            assert isinstance(compartment, contracts.Compartment)
            assert units == "mM" and meaning.strip() and name.strip()

    def test_the_engine_exposes_them_on_a_basis_that_cannot_integrate(self, axis):
        exposed = {name: var for name, var in engine.REGISTERED_AXIS_VARIABLES.items()
                   if name.startswith("axis.atp_glycolysis.")}
        assert set(exposed) == set(axis.state_names)
        assert {var.basis for var in exposed.values()} == {"registered, not integrated"}

    def test_the_pool_the_axis_is_named_for_is_not_even_one_of_its_states(self, axis):
        """ATP, ADP and AMP are algebraic in P, which is why nothing can read them directly."""
        registered = {name.rsplit(".", 1)[1] for name in axis.state_names}
        assert not registered & {"ATP", "ADP", "AMP"}
        assert "P" in registered

    def test_the_registered_names_do_not_collide_with_any_integrated_state(self, engine_ladder):
        assert not set(engine.REGISTERED_AXIS_VARIABLES) & set(engine_ladder[50.0]["truth"])

    def test_it_is_undriven_here_and_in_the_flux_map(self, axis):
        assert axis.readers == () and axis.driven is False
        assert signalling.TRANSCRIBED_AXIS_FLUX_COUPLING["atp_glycolysis"] == ()

    def test_the_claimed_missing_medium_species_really_is_absent(self, axis):
        assert axis.missing_medium == (
            "MEDIUM_SPECIES_NOT_CARRIED::medium.respiratory_inhibitor",)
        assert "respiratory_inhibitor" not in contracts.EXTRACELLULAR
        with pytest.raises(Exception):
            float(contracts.refused_row(axis.missing_medium[0]))

    def test_no_promoter_refusal_is_owed_because_atp_is_a_metabolite_pool(self, axis):
        """Unlike cell_wall_slt2's RLM1 box, this axis owes no transcriptional refusal."""
        assert axis.panel_promoter == "none"
        assert axis.promoter_refusals == ()

    def test_the_engine_result_reports_it_as_registered_and_undriven(self, engine_ladder):
        line = [s for s in engine_ladder[50.0]["validity"].unsupported
                if "TRANSCRIBED_AXES" in s]
        assert len(line) == 1
        assert "atp_glycolysis" in line[0] and "REGISTERED AND UNDRIVEN" in line[0]


class TestWhatTheReplacementWouldBeReplacing:
    """The engine's adenylate pool costed in priors, derived rather than quoted."""

    def test_the_reactions_that_move_the_pool_are_found_not_listed(self, adenylate_priors):
        """Read off the engine's own stoichiometry, so a new reaction cannot hide from it."""
        assert len(adenylate_priors["reactions"]) == 35
        assert "adenylate_kinase" in adenylate_priors["reactions"]
        assert len(adenylate_priors["reaction_constants"]) == 11

    def test_every_constant_behind_them_is_stamped_an_unpublished_prior(self,
                                                                       adenylate_priors):
        assert len(adenylate_priors["stamped"]) == 19
        assert not adenylate_priors["unstamped"]

    def test_the_one_adenylate_cost_that_is_not_a_bare_prior_shows_what_promotion_looks_like(
            self):
        """translation_atp already carries a citation, which is the grade this file's
        source could reach and, per the module's own docstring, does not."""
        translation = engine.EngineParameters.prior().values["translation_atp"]
        assert not translation.source.startswith("PRIOR, not a published constant")
        assert "Stouthamer" in translation.source
        assert teusink.TEUSINK_PARAMS.by_tag()["MEASURED"] == 0

    def test_the_growth_reaction_is_the_only_adenylate_source_in_the_engine(self):
        """Teusink's SUM_P is strictly conserved; the engine's is conserved per gDW only."""
        kernel = engine._Kernel(engine.EngineParameters.prior(), contracts.Genotype.prior())
        pool = ("internal.atp", "internal.adp", "internal.amp", "internal.camp")
        column = sum(kernel.matrix[kernel.index[name]] for name in pool)
        moving = {name: value for name, value in zip(kernel.reaction_names, column)
                  if abs(value) > 1e-12}
        assert set(moving) == {"growth"}
        with pytest.raises(Exception):
            float(teusink.GROWTH_DILUTION)


class TestThePoolDisagrees:
    """Comparison 1. Both pools run to their own resting state on glucose and compared."""

    def test_the_transcription_lands_on_table_4_before_anything_is_compared_to_it(self):
        reference = teusink.steady_state()
        assert reference.atp == pytest.approx(teusink.PUBLISHED_TABLE4["ATP"], abs=0.01)
        assert reference.adp == pytest.approx(teusink.PUBLISHED_TABLE4["ADP"], abs=0.01)
        assert reference.amp == pytest.approx(teusink.PUBLISHED_TABLE4["AMP"], abs=0.01)
        assert reference.energy_charge == pytest.approx(0.769, abs=0.001)

    def test_the_engines_own_priors_start_within_two_percent_of_the_published_charge(
            self, engine_priors):
        """The asserted initial condition very nearly IS Teusink's published steady state."""
        assert engine_priors["charge"] == pytest.approx(0.7812, abs=0.001)
        assert abs(engine_priors["charge"] / 0.769376 - 1.0) < 0.02
        assert engine_priors["mass_action"] == pytest.approx(0.400, abs=0.001)

    def test_and_the_engines_own_dynamics_abandon_it_inside_fifteen_minutes(
            self, engine_ladder, engine_priors):
        relaxed = engine_ladder[50.0]
        assert relaxed["charge"] == pytest.approx(0.9724, abs=0.01)
        assert relaxed["charge"] - engine_priors["charge"] > 0.18
        assert relaxed["mass_action"] == pytest.approx(7.24, abs=0.4)
        assert relaxed["mass_action"] / KEQ_AK > 14.0

    def test_the_engine_charge_is_stationary_over_the_run_so_the_read_time_is_not_a_choice(
            self, engine_ladder):
        trajectory = engine_ladder[50.0]["charge_trajectory"][1:]
        assert float(np.max(trajectory) - np.min(trajectory)) < 0.01
        assert float(np.min(trajectory)) > 0.96

    @pytest.mark.parametrize("glucose,engine_charge,teusink_charge", [
        (50.0, 0.9724, 0.7692), (20.0, 0.9734, 0.7411), (5.0, 0.9782, 0.6051)])
    def test_the_two_energy_charges_disagree_at_every_matched_glucose(
            self, engine_ladder, teusink_ladder, glucose, engine_charge, teusink_charge):
        measured, published = engine_ladder[glucose], teusink_ladder[glucose]
        assert measured["charge"] == pytest.approx(engine_charge, abs=0.01)
        assert published.energy_charge == pytest.approx(teusink_charge, abs=0.005)
        assert measured["charge"] - published.energy_charge > 0.19

    @pytest.mark.parametrize("glucose,engine_ratio,teusink_ratio", [
        (50.0, 26.68, 1.940), (20.0, 27.61, 1.693), (5.0, 33.03, 0.974)])
    def test_the_atp_to_adp_ratio_disagrees_by_more_than_an_order_of_magnitude(
            self, engine_ladder, teusink_ladder, glucose, engine_ratio, teusink_ratio):
        measured, published = engine_ladder[glucose], teusink_ladder[glucose]
        assert measured["atp"] / measured["adp"] == pytest.approx(engine_ratio, rel=0.05)
        assert published.atp / published.adp == pytest.approx(teusink_ratio, rel=0.01)
        assert (measured["atp"] / measured["adp"]) / (published.atp / published.adp) > 13.0

    def test_the_adp_to_amp_ratio_is_the_one_that_nearly_agrees_and_then_stops(
            self, engine_ladder, teusink_ladder):
        """At 50 mM the two are 15% apart; by 5 mM Teusink has halved and the engine has not."""
        top = engine_ladder[50.0]["adp"] / engine_ladder[50.0]["amp"]
        bottom = engine_ladder[5.0]["adp"] / engine_ladder[5.0]["amp"]
        assert top == pytest.approx(3.686, rel=0.05)
        assert bottom == pytest.approx(4.071, rel=0.05)
        assert teusink_ladder[50.0].adp / teusink_ladder[50.0].amp == pytest.approx(4.311, rel=0.01)
        assert teusink_ladder[5.0].adp / teusink_ladder[5.0].amp == pytest.approx(2.165, rel=0.01)
        assert abs(top / (teusink_ladder[50.0].adp / teusink_ladder[50.0].amp) - 1.0) < 0.2
        assert bottom / (teusink_ladder[5.0].adp / teusink_ladder[5.0].amp) > 1.8

    def test_the_charge_response_to_the_one_shared_dose_has_the_OPPOSITE_SIGN(
            self, engine_ladder, teusink_ladder):
        """The sharpest result in this file: the level gap is a disagreement, this is a
        contradiction. Glucose falls 49.6 -> 4.7 mM and the two charges move apart."""
        measured = np.array([engine_ladder[g]["charge"] for g in PAIRED_RUNGS])
        published = np.array([teusink_ladder[g].energy_charge for g in PAIRED_RUNGS])
        assert np.all(np.diff(measured) > 0.0)
        assert np.all(np.diff(published) < 0.0)
        assert measured[-1] / measured[0] - 1.0 == pytest.approx(0.0060, abs=0.004)
        assert published[-1] / published[0] - 1.0 == pytest.approx(-0.2133, abs=0.005)

    def test_an_eighth_of_the_engines_adenylate_sits_where_the_source_has_no_species(
            self, engine_ladder):
        """SUM_P cannot map onto atp+adp+amp alone, priced rather than merely declared."""
        assert engine_ladder[50.0]["camp_fraction"] == pytest.approx(0.130, abs=0.02)
        assert engine_ladder[5.0]["camp_fraction"] == pytest.approx(0.110, abs=0.02)
        assert teusink.TRANSFER_ASSUMPTIONS["no_camp_reservoir"].strip()


class TestOnlyDimensionlessQuantitiesWereCompared:
    """The unit caveat, enforced rather than written down."""

    def test_the_bridge_that_would_be_needed_is_refused_at_the_source(self):
        with pytest.raises(Exception):
            float(teusink.CYTOSOLIC_WATER_BRIDGE)
        with pytest.raises(Exception):
            teusink.to_mmol_per_gdw(float(teusink.SUM_P))

    def test_so_the_absolute_pools_are_left_uncompared(self):
        """SUM_P is mM and the engine's initial pool is mmol/gDW; nothing above converts."""
        assert teusink.SUM_P.units == "mM"
        assert engine.EngineParameters.prior().values["initial_atp"].units == "mmol/gDW"
        assert teusink.TRANSFER_ASSUMPTIONS["mM_to_mmol_per_gdw"].strip()

    def test_the_engines_deepest_rung_falls_outside_the_transcriptions_domain(
            self, engine_ladder):
        """At the 2 mM rung the engine has already drawn the medium below the validated
        floor, so no paired number exists there and none is invented."""
        residual = engine_ladder[2.0]["residual_glucose_mM"]
        assert residual == pytest.approx(1.754, abs=0.05)
        assert residual < teusink.VALIDATED_GLCO_FLOOR_mM
        with pytest.raises(ValueError, match="validated floor"):
            teusink.steady_state(glco=residual)


class TestTheKeqDirectionIsATrap:
    """Comparison 2, recorded as a trap rather than claimed as a result.

    The sign of a net rate at a state is fixed by its numerator alone, so the SIGN claim
    below is complete. The size of the charge excursion that would follow is NOT claimed:
    that needs a rewired integration inside ``mech/engine.py``, which this file does not
    own and does not perform.
    """

    def test_the_deposits_constant_is_the_equilibrium_of_two_adp_to_atp_plus_amp(
            self, teusink_ladder):
        """Verified from the source's own output at every dose, not read off the label."""
        for state in teusink_ladder.values():
            assert state.atp * state.amp / state.adp ** 2 == pytest.approx(KEQ_AK, rel=1e-6)
        assert KEQ_AK == pytest.approx(0.45, abs=1e-12)

    def test_the_engine_writes_that_same_reaction_reversed(self):
        kernel = engine._Kernel(engine.EngineParameters.prior(), contracts.Genotype.prior())
        assert kernel.stoichiometry["adenylate_kinase"] == ENGINE_DIRECTION

    def test_adenylate_kinase_is_charge_neutral_so_the_trap_cannot_show_in_the_charge(self):
        """Both models conserve ATP + ADP/2 across this reaction exactly, so the error is
        invisible in the charge and visible only in the split and the flux direction."""
        kernel = engine._Kernel(engine.EngineParameters.prior(), contracts.Genotype.prior())
        stoichiometry = kernel.stoichiometry["adenylate_kinase"]
        assert stoichiometry["internal.atp"] + 0.5 * stoichiometry["internal.adp"] == 0.0
        phosphate = 2.0 * teusink.PUBLISHED_TABLE4["ATP"] + teusink.PUBLISHED_TABLE4["ADP"]
        charges = [teusink.energy_charge(*teusink.adenylate(phosphate, keq_ak=keq)[:2])
                   for keq in (KEQ_AK, 1.0, 1.0 / KEQ_AK)]
        assert charges[0] == pytest.approx(charges[1]) == pytest.approx(charges[2])

    def test_un_inverted_the_split_moves_adp_and_amp_in_opposite_directions(self):
        """At the deposit's own phosphate pool, the un-inverted constant deflates ADP to
        0.63x and inflates AMP to 1.80x -- and AMP is the only feedback in the model."""
        reference = teusink.steady_state()
        phosphate = 2.0 * reference.atp + reference.adp
        right = teusink.adenylate(phosphate, keq_ak=KEQ_AK)
        wrong = teusink.adenylate(phosphate, keq_ak=1.0 / KEQ_AK)
        assert float(wrong[1]) / float(right[1]) == pytest.approx(0.631, abs=0.01)
        assert float(wrong[2]) / float(right[2]) == pytest.approx(1.796, abs=0.01)
        assert (float(wrong[1]) - float(right[1])) * (float(wrong[2]) - float(right[2])) < 0.0

    @pytest.mark.parametrize("glucose", GLUCOSE_LADDER_mM)
    def test_the_two_orientations_disagree_in_SIGN_at_every_state_the_engine_occupies(
            self, engine_ladder, glucose):
        """Inverted: forward, ATP consumed. Un-inverted: backward, ATP regenerated."""
        forces = engine_ladder[glucose]["forces"]
        assert forces["inverted"] > 0.0
        assert forces["un_inverted"] < 0.0
        assert forces["inverted"] * forces["un_inverted"] < 0.0

    @pytest.mark.parametrize("glucose", GLUCOSE_LADDER_mM)
    def test_and_the_wrong_orientation_is_the_larger_of_the_two(self, engine_ladder, glucose):
        forces = engine_ladder[glucose]["forces"]
        assert abs(forces["un_inverted"] / forces["inverted"]) == pytest.approx(2.211, abs=0.02)

    def test_both_dwarf_the_driving_force_the_engine_currently_carries(self, engine_ladder):
        """As written the engine sits within 0.4% of its own implicit equilibrium, so either
        orientation converts a near-null flux into a large driven one."""
        forces = engine_ladder[50.0]["forces"]
        assert forces["as_written"] == pytest.approx(2.07e-4, rel=0.2)
        assert abs(forces["inverted"] / forces["as_written"]) > 200.0
        assert abs(forces["un_inverted"] / forces["as_written"]) > 450.0

    def test_the_engines_live_flux_agrees_in_sign_with_its_own_numerator(self, engine_ladder):
        """Ties the algebra above back to the integrated engine, so it is the same reaction."""
        run = engine_ladder[50.0]
        assert run["flux_ak"] > 0.0 and run["forces"]["as_written"] > 0.0

    def test_the_registration_states_the_trap_it_is_being_held_to(self, axis):
        assert "inverting" in axis.undriven_because
        assert "2 ADP -> ATP + AMP" in axis.undriven_because
        assert teusink.ADENYLATE_KINASE_EQUILIBRIUM_ONLY.strip()


class TestTheTrapHidesFromTheObviousCheck:
    """Why the demonstration above had to be done in the engine's own coordinates."""

    @pytest.mark.parametrize("glucose", GLUCOSE_LADDER_mM)
    def test_in_concentrations_the_two_orientations_agree_in_sign_everywhere(
            self, engine_ladder, glucose):
        """A mass-action check would pass a wiring that is inverted where it really runs."""
        concentration = engine_ladder[glucose]["concentration_forces"]
        assert concentration["inverted"] > 0.0
        assert concentration["un_inverted"] > 0.0
        assert np.sign(concentration["inverted"]) == np.sign(concentration["un_inverted"])

    def test_because_the_pool_sits_outside_the_band_the_two_can_differ_in(self, engine_ladder):
        """They can only disagree between KeqAK and 1/KeqAK; concentrations are well above."""
        quotient = engine_ladder[50.0]["concentration_forces"]["quotient"]
        assert quotient == pytest.approx(7.24, rel=0.05)
        assert quotient > 1.0 / KEQ_AK > KEQ_AK

    def test_and_the_engines_own_coordinates_put_it_squarely_inside_that_band(
            self, engine_ladder):
        quotient = engine_ladder[50.0]["forces"]["quotient"]
        assert quotient == pytest.approx(1.002, abs=0.01)
        assert KEQ_AK < quotient < 1.0 / KEQ_AK

    @pytest.mark.parametrize("glucose", GLUCOSE_LADDER_mM)
    def test_the_same_pool_quotes_two_quotients_seven_fold_apart(self, engine_ladder, glucose):
        """A THIRD refusal, beyond the two the registration names: there is no single number
        that is the engine's adenylate kinase Keq, so KeqAK has nothing to replace."""
        run = engine_ladder[glucose]
        assert run["concentration_forces"]["quotient"] / run["forces"]["quotient"] > 7.0

    def test_the_missing_bridge_between_them_is_another_engine_prior(self):
        values = engine.EngineParameters.prior().values
        assert values["energy_k"].source.startswith("PRIOR, not a published constant")
        assert values["energy_k"].units == "mmol/gDW"


class TestNoneOfThisLicensesWiringTheAxis:
    """The comparison is a price, not a permission. Nothing here retires a refusal."""

    def test_the_sources_own_gate_still_refuses(self):
        assert "REFUSED" in teusink.gate().summary()

    @pytest.mark.parametrize("name,param", [
        ("teusink.stress_to_atp_demand", teusink.STRESS_ATP_DEMAND),
        ("teusink.cytosol_l_per_gdw", teusink.CYTOSOLIC_WATER_BRIDGE),
        ("teusink.respiratory_atp_yield", teusink.RESPIRATORY_ATP_YIELD),
        ("teusink.growth_dilution", teusink.GROWTH_DILUTION)])
    def test_all_four_absent_quantities_still_raise_on_use(self, name, param):
        assert name in teusink.provenance()["refusals"]
        with pytest.raises(Exception):
            float(param)

    def test_the_axis_is_still_registered_and_undriven_after_every_run_above(
            self, axis, engine_ladder):
        assert axis.readers == () and axis.driven is False
        assert all("REGISTERED AND UNDRIVEN" in " ".join(engine_ladder[g]["validity"].unsupported)
                   for g in GLUCOSE_LADDER_mM)

    def test_the_panel_arm_this_axis_cannot_answer_is_still_named(self):
        assert "antimycin_A" in teusink.PANEL_LIMIT
        assert "does not close the atp module" in teusink.PANEL_LIMIT
