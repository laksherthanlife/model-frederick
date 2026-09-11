"""Every layer the vision names is reachable from `predict_product`, and says which it is.

`docs/DISTANCE_TO_THE_VISION.md` opened with a table of five named layers of which one was
in the prediction chain, and the reason that table could be written at all is that somebody
went and read `predict.py`'s imports. Before that, this module's own docstring claimed the
FBA audit was in the chain while `grep -n '^from' src/ystwin/predict.py` returned no
`ystwin.fba` -- a false claim that no test could see, because no test asked the question.

These tests ask it. They are not about arithmetic; every number here is checked elsewhere.
They pin three things a refactor can quietly break:

1. **Every layer is reachable.** One call reaches metabolism, the flux law, the stress
   panel, the GEM audit, E-Flux and the thermodynamic gate, and a seventh -- the latent
   stress branch -- is computed and attached. A layer that stops being reachable fails here
   rather than in a document nobody re-derives.

2. **Only the layers that may move the number do.** The audit, the gate and the latent
   branch are checks and records. The content must be bit-for-bit identical with them and
   without them, and `==` on the float is the assertion, not `approx`. This is what keeps
   an unvalidated branch from laundering itself into a headline number: if wiring the
   latent branch ever starts changing the answer, that is a decision someone must make
   deliberately, and this test is where they will be told they made it.

3. **The chemistry the gate reads matches the model it is gated against.** Both shipped
   specs declare a stoichiometry *and* a reaction id for the same step; the duplication
   only buys anything if something compares them, and `check_steps_against` is that.

Nothing here needs a GEM except the class that says so, and those skip without one.
"""

from __future__ import annotations

import math

import pytest

from ystwin import paths
from ystwin.bridge.latent_bridge import LatentState
from ystwin.bridge.thermodynamic import ThermodynamicData
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.flux import FluxCalibration
from ystwin.pathway.proteome import enzyme_content_from_mass_fraction
from ystwin.pathway.solve import NodeKinetics
from ystwin.pathway.spec import Node, PathwaySpec, load_pathway
from ystwin.pathway.thermo_gate import Feasibility, ThermodynamicallyBlocked
from ystwin.predict import (
    Environment,
    Genotype,
    LayerState,
    LegacyProductDiagnostic,
    predict_product,
)

# The five the goal names, plus the two this module splits out of them: the flux law is
# what "metabolism" actually predicts from, and E-Flux is reported separately from the GEM
# because it runs inside the audit and is inert for a different reason.
EXPECTED_LAYERS = (
    "metabolism",
    "flux law",
    "stress panel",
    # Added 2026-09-02. The chain was MEASURED to have exactly one route from the environment
    # to the product -- mu -- and yeast disagrees: Kocharin 2013 moved PHB content 3.82x
    # between glucose and ethanol at an IDENTICAL growth rate, one genotype, one vessel,
    # 64 standard deviations apart. This layer is REPORTED for beta_carotene, which has not
    # measured its own coefficient, and IN_CHAIN for phb, which has. It is not INERT, and the
    # distinction is the whole point: INERT would record the absence of a coefficient as the
    # absence of an effect.
    "environment",
    # Added 2026-08-31 with `Genotype.cassette`, and added because the architecture review
    # caught the field being accepted with no status line the same day it was written --
    # the exact failure this list exists to prevent.
    "cassette dosage",
    # Added 2026-09-02: `vmax = kcat * [E]` audited against the fitted ceiling. AUDITS and
    # never sets the number -- it is scored against nothing while the fitted capacity is
    # scored leave-one-strain-out on six states. It is here because it reaches the 63x
    # ceiling refutation from the enzyme side, which a paper had to supply before.
    "enzyme capacity",
    "GEM / FBA",
    "regulation (E-Flux)",
    "thermodynamics",
    "latent stress state",
)

# Cofactor and precursor levels the solve does not pin. Order-of-magnitude placeholders,
# and that is all they need to be: no claim in this file rests on their values, only on
# the gate reaching a verdict rather than a coverage gap.
CAROTENOID_BACKGROUND = {
    "s_0189": 1e-5,   # GGPP, the precursor
    "s_0633": 1e-4,   # diphosphate
    "s_0687": 1e-4,   # FAD
    "s_0689": 1e-4,   # FADH2
    "s_0794": 1e-7,   # proton, pH 7
}

# From scripts/thiolase_threshold.py. The two measured cytosolic pools, which the lead in
# docs/DISTANCE_TO_THE_VISION.md section 4 held to sit on opposite sides of the threshold.
# They do not: the threshold is 19 mM and both are below it. Kept because they are still the
# two real concentrations, and gating at each is still how the chain is exercised.
ACETYL_COA_GLUCOSE_M = 10e-6
ACETYL_COA_ETHANOL_M = 425e-6
PHB_BACKGROUND = {"s_0367": 1e-6, "s_0529": 100e-6}


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def thermo():
    return ThermodynamicData.load()


@pytest.fixture(scope="module")
def phb_spec():
    return load_pathway("phb")


@pytest.fixture(scope="module")
def phb_calibration():
    """Illustrative, and it does not matter: nothing here reads the flux's magnitude.

    `phb` carries no entry in `calibrations.py` deliberately -- Kocharin 2012 measures no
    expression, so there is no scalar to ship. What this fixture supplies is a well-formed
    calibration so the chain runs, not a fitted one.
    """
    return FluxCalibration(alpha=9.93e-2, entry_enzyme="PhaA", loso_rmse_log=0.2,
                           loso_skill=0.0, n_states=2, source="illustrative, not fitted")


def _predict(spec, **kwargs):
    return predict_product(spec, Genotype(1.0, "b-car4"), Environment(growth_rate_setpoint_per_h=0.18),
                           BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, **kwargs)


def _legacy(spec, **kwargs):
    result = _predict(spec, mode="legacy", **kwargs)
    assert isinstance(result, LegacyProductDiagnostic)
    assert result.supported is False
    assert result.calculation.mode == "legacy"
    return result


class TestEveryLayerIsAccountedFor:
    def test_a_bare_call_names_all_nine(self, spec):
        assert tuple(layer.name for layer in _predict(spec).layers) == EXPECTED_LAYERS

    def test_a_layer_nobody_supplied_says_what_would_run_it(self, spec):
        for layer in _predict(spec).layers:
            if layer.state == LayerState.NOT_RUN:
                assert "pass " in layer.detail, layer.name

    def test_exactly_the_layers_that_predict_are_marked_as_predicting(self, spec):
        setting = {layer.name for layer in _predict(spec).layers if layer.sets_the_number}

        assert setting == {"metabolism", "flux law"}

    def test_the_stress_panel_audits_in_a_chemostat_rather_than_being_inert(self, spec):
        """It read INERT until 2026-08-31 and that understated it.

        `INERT` means "ran, and provably could not have changed anything". This layer can
        end the call: past the lethal dose it raises `SetpointUnreachable` instead of returning a
        number. A layer with a refusal channel is `AUDITS`, which is the enum's own word for
        "ran, can refute, never moves the number".
        """
        got = _predict(spec, thermo=None)

        assert got.layer("stress panel").state == LayerState.AUDITS

    def test_and_that_refusal_channel_is_real(self, spec):
        """The reason it is not INERT, exercised rather than asserted."""
        from ystwin.predict import SetpointUnreachable

        with pytest.raises(SetpointUnreachable):
            predict_product(
                spec, Genotype(1.0, "b-car4"),
                Environment(growth_rate_setpoint_per_h=0.18, stressor="DTT", dose=10.0),
                BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS)

    def test_and_is_in_the_chain_in_batch_where_it_sets_the_growth_rate(self, spec):
        """The one condition under which stress reaches the number at all: in batch mu is
        what the cells can do, so a stressor moves it and every pool is solved at it."""
        loose = {"lycopene": NodeKinetics(
            vmax_per_growth=BETA_CAROTENE_KINETICS["lycopene"].vmax_per_growth,
            km=BETA_CAROTENE_KINETICS["lycopene"].km)}
        got = predict_product(spec, Genotype(1.0), Environment(stressor="DTT", dose=1.0),
                              BETA_CAROTENE_FLUX, loose)

        assert got.layer("stress panel").state == LayerState.IN_CHAIN

    def test_an_unknown_layer_is_a_key_error_naming_the_ones_there_are(self, spec):
        with pytest.raises(KeyError, match="thermodynamics"):
            _predict(spec).layer("regulation")

    def test_the_report_has_one_line_per_layer(self, spec):
        got = _predict(spec)

        assert len(got.layer_report().splitlines()) == len(EXPECTED_LAYERS)


class TestTheGateRunsWithoutAModel:
    """The point of declaring a stoichiometry beside the reaction id: no SBML parse."""

    def test_the_thermodynamic_layer_audits_rather_than_sitting_out(self, spec, thermo):
        got = _predict(spec, thermo=thermo, thermo_background_m=CAROTENOID_BACKGROUND)

        assert got.layer("thermodynamics").state == LayerState.AUDITS

    def test_every_declared_step_is_gated_and_none_is_left_out(self, spec, thermo):
        report = _predict(spec, thermo=thermo,
                          thermo_background_m=CAROTENOID_BACKGROUND).thermodynamics[0]

        assert report.ungated == ()

    def test_the_carotenoid_steps_come_back_cannot_say_rather_than_runs(self, spec, thermo):
        """The tables carry no formation energy for a C40 carotenoid, and reporting that
        as feasible would turn a coverage gap into a clearance. This is the honest failure
        and it is worth pinning: a change that starts returning `runs` here has invented an
        energy rather than found one."""
        report = _predict(spec, thermo=thermo,
                          thermo_background_m=CAROTENOID_BACKGROUND).thermodynamics[0]

        assert {s.feasibility for s in report.steps} == {Feasibility.CANNOT_SAY}

    def test_one_report_per_volume_because_the_conversion_is_a_convention(self, spec, thermo):
        got = _predict(spec, thermo=thermo, thermo_background_m=CAROTENOID_BACKGROUND,
                       thermo_volumes=(1.0, 2.0, 2.7))

        assert [r.cytosolic_volume_ml_per_gdcw for r in got.thermodynamics] == [1.0, 2.0, 2.7]

    def test_a_spec_with_no_declared_chemistry_refuses_an_unusable_thermo_input(self, thermo):
        bare = PathwaySpec(
            product="x", organism="S. cerevisiae", entry_enzyme="CrtE",
            precursor_metabolite="s_0189",
            nodes=(Node("a", "passthrough"), Node("x", "passthrough",
                                                  molar_mass_g_per_mol=100.0)))
        with pytest.raises(ValueError, match="declares no step chemistry"):
            predict_product(bare, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.18),
                            BETA_CAROTENE_FLUX, {}, thermo=thermo)


class TestTheThiolaseThresholdReachesThePredictionChain:
    """The mechanism section 4 of DISTANCE_TO_THE_VISION.md identifies, now inside the
    chain rather than only inside a standalone script.

    `phb`'s entry step is yeast-GEM's cytosolic thiolase r_0103, dGr'0 = +38.07 kJ/mol --
    uphill, so it has a concentration threshold, and at 19 mM acetyl-CoA neither feed is
    anywhere near it.

    **What is pinned here is the WIRING, and that is why this class survived the
    refutation.** It read +15.6 kJ/mol and `RUNS` on ethanol until 2026-08-30, when the
    energy was found to rest on a kcal/kJ error in `bridge/thermodynamic.py`. The thiolase
    lead is gone with it -- see `tests/test_thiolase_threshold.py` -- but the thing this
    class exists to check is that `predict_product` reaches the gate at all, gets a verdict
    rather than `cannot_say`, and reports the ungated steps honestly. All three still hold.
    A gate is not less wired because the answer it returns changed.
    """

    def _gate(self, phb_spec, phb_calibration, thermo, acetyl_coa_m):
        got = predict_product(
            phb_spec, Genotype(1.0, "SCKK006"), Environment(growth_rate_setpoint_per_h=0.1),
            phb_calibration, {}, thermo=thermo, mode="legacy",
            thermo_background_m={**PHB_BACKGROUND, "s_0373": acetyl_coa_m})
        assert isinstance(got, LegacyProductDiagnostic)
        assert got.supported is False
        return got.calculation.thermodynamics[0].step("acetoacetyl_coa")

    def test_the_entry_step_is_covered_by_the_tables(self, phb_spec, phb_calibration, thermo):
        step = self._gate(phb_spec, phb_calibration, thermo, ACETYL_COA_ETHANOL_M)

        assert step.feasibility != Feasibility.CANNOT_SAY
        assert step.dg0_kj_per_mol == pytest.approx(38.07, abs=0.2)

    def test_on_glucose_the_step_cannot_run(self, phb_spec, phb_calibration, thermo):
        step = self._gate(phb_spec, phb_calibration, thermo, ACETYL_COA_GLUCOSE_M)

        assert step.feasibility == Feasibility.CANNOT_RUN

    def test_on_ethanol_it_cannot_run_either(self, phb_spec, phb_calibration, thermo):
        """Asserted `RUNS` until 2026-08-30. Both feeds are below a 19 mM threshold, so the
        two no longer differ by feed -- which refutes the lead and leaves the wiring intact.
        What this still pins is that the chain reaches a real verdict on real chemistry."""
        step = self._gate(phb_spec, phb_calibration, thermo, ACETYL_COA_ETHANOL_M)

        assert step.feasibility == Feasibility.CANNOT_RUN

    def test_the_two_thioesters_below_it_are_reported_as_ungated(
            self, phb_spec, phb_calibration, thermo):
        """Heterologous steps with no GEM reaction and no declared stoichiometry. Saying so
        is the difference between a gate that cleared a pathway and one that looked at a
        third of it."""
        got = predict_product(
            phb_spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.1), phb_calibration, {},
            thermo=thermo, thermo_background_m={**PHB_BACKGROUND, "s_0373": ACETYL_COA_ETHANOL_M},
            mode="legacy")

        assert isinstance(got, LegacyProductDiagnostic)
        assert got.supported is False
        assert got.calculation.thermodynamics[0].ungated == ("r3_hydroxybutyryl_coa", "phb")

    def test_a_blocked_step_can_be_made_to_refuse_rather_than_annotate(
            self, phb_spec, phb_calibration, thermo):
        with pytest.raises(ThermodynamicallyBlocked):
            predict_product(
                phb_spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.1), phb_calibration, {},
                thermo=thermo,
                thermo_background_m={**PHB_BACKGROUND, "s_0373": ACETYL_COA_GLUCOSE_M},
                refuse_thermodynamically_blocked=True)

    def test_reports_only_in_explicit_legacy_mode(self, phb_spec, phb_calibration, thermo):
        got = predict_product(
            phb_spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.1), phb_calibration, {},
            thermo=thermo, mode="legacy",
            thermo_background_m={**PHB_BACKGROUND, "s_0373": ACETYL_COA_GLUCOSE_M})

        assert isinstance(got, LegacyProductDiagnostic)
        assert got.supported is False
        assert any("CONTRADICTS ITSELF" in note for note in got.calculation.notes)

    @pytest.mark.parametrize("acetyl_coa_m", [ACETYL_COA_GLUCOSE_M, ACETYL_COA_ETHANOL_M])
    def test_default_refuses_a_blocked_step_without_an_opt_in_flag(
            self, phb_spec, phb_calibration, thermo, acetyl_coa_m):
        with pytest.raises(ThermodynamicallyBlocked):
            predict_product(
                phb_spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.1),
                phb_calibration, {}, thermo=thermo,
                thermo_background_m={**PHB_BACKGROUND, "s_0373": acetyl_coa_m})


class TestWiringALayerInDoesNotMoveTheNumber:
    """`==`, not `approx`. The whole argument for attaching an unvalidated branch rather
    than deleting it is that it cannot reach the answer, and an assertion with a tolerance
    in it would let a small contribution through unnoticed."""

    def test_the_thermodynamic_gate_changes_nothing(self, spec, thermo):
        bare = _predict(spec)
        gated = _predict(spec, thermo=thermo, thermo_background_m=CAROTENOID_BACKGROUND)

        assert gated.content_mmol_per_gdcw == bare.content_mmol_per_gdcw

    def test_nor_does_the_latent_branch(self, spec):
        bare = _predict(spec)
        latent = _legacy(spec, latent=LatentState(0.02, "UPRE2", 4.0)).calculation

        assert latent.content_mmol_per_gdcw == bare.content_mmol_per_gdcw

    def test_nor_do_both_together(self, spec, thermo):
        bare = _predict(spec)
        both = _legacy(spec, thermo=thermo, thermo_background_m=CAROTENOID_BACKGROUND,
                       latent=LatentState(0.02, "UPRE2", 4.0)).calculation

        assert both.content_mmol_per_gdcw == bare.content_mmol_per_gdcw
        assert both.flux.flux_mmol_per_gdcw_h == bare.flux.flux_mmol_per_gdcw_h

    def test_the_latent_branch_is_carried_and_marked_unvalidated(self, spec):
        got = _legacy(spec, latent=LatentState(0.02, "UPRE2", 4.0)).calculation

        assert got.layer("latent stress state").state == LayerState.REPORTED
        assert "REFUSED by G4" in got.layer("latent stress state").detail

    def test_and_its_maintenance_demand_is_a_real_number_not_a_placeholder(self, spec):
        got = _legacy(spec, latent=LatentState(0.02, "UPRE2", 4.0)).calculation

        assert got.latent_constraints.atp_maintenance > got.latent_constraints.resting_maintenance

    def test_a_note_says_the_number_is_what_it_would_be_without_it(self, spec):
        got = _legacy(spec, latent=LatentState(0.02, "UPRE2", 4.0)).calculation

        assert any("COMPUTED AND NOT USED" in note for note in got.notes)

    def test_default_refuses_an_unvalidated_latent_prediction(self, spec):
        with pytest.raises(ValueError, match="latent"):
            _predict(spec, latent=LatentState(0.02, "UPRE2", 4.0))


class TestTheDeclaredChemistryIsCheckedAgainstTheModel:
    """Skipped without a GEM, because there is nothing to compare against."""

    @pytest.fixture(scope="class")
    def model(self):
        import warnings

        path = paths.yeast_gem()
        if path is None or not path.is_file():
            pytest.skip("yeast-GEM not present; set YSTWIN_YEAST_GEM")
        from ystwin.fba.carotenoid import add_beta_carotene_pathway
        from ystwin.fba.solver import load_model

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            loaded, _ = load_model(path)
        return add_beta_carotene_pathway(loaded)

    def test_every_carotenoid_step_agrees_with_what_carotenoid_py_builds(self, spec, model):
        assert spec.check_steps_against(model) == ("phytoene", "lycopene", "beta_carotene")

    def test_the_phb_thiolase_agrees_with_yeast_gems_r_0103(self, phb_spec, model):
        assert phb_spec.check_steps_against(model) == ("acetoacetyl_coa",)

    def test_a_drifted_coefficient_is_refused_by_name(self, model):
        drifted = PathwaySpec(
            product="beta_carotene", organism="S. cerevisiae", entry_enzyme="CrtE",
            precursor_metabolite="s_0189",
            nodes=(Node("lycopene", "passthrough", metabolite="lycopene_c", reaction="CRTI",
                        stoichiometry={"lycopene_c": 1.0, "phytoene_c": -1.0}),
                   Node("beta_carotene", "passthrough", molar_mass_g_per_mol=536.87)))

        with pytest.raises(ValueError, match="different reactions"):
            drifted.check_steps_against(model)

    def test_a_step_the_model_does_not_carry_is_skipped_not_failed(self, model):
        """A heterologous reaction that is not installed yet is a coverage gap, not a
        mismatch, and treating it as one would make the check unusable on any spec whose
        pathway has to be added to the host first -- which is every spec here, before
        `add_beta_carotene_pathway` runs."""
        uninstalled = PathwaySpec(
            product="beta_carotene", organism="S. cerevisiae", entry_enzyme="CrtE",
            precursor_metabolite="s_0189",
            nodes=(Node("lycopene", "passthrough", metabolite="lycopene_c",
                        reaction="CRTI_NOT_INSTALLED_YET",
                        stoichiometry={"lycopene_c": 1.0, "phytoene_c": -1.0}),
                   Node("beta_carotene", "passthrough", molar_mass_g_per_mol=536.87)))

        assert uninstalled.check_steps_against(model) == ()

    def test_the_whole_chain_runs_with_every_layer_at_once(self, spec, model, thermo):
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID

        got = predict_product(
            spec, Genotype(1.0, "b-car4"),
            Environment(growth_rate_setpoint_per_h=0.18, stressor="DTT", dose=1.0),
            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
            audit_model=model, audit_reaction=PRODUCT_DEMAND_ID, audit_glucose_uptake=10.0,
            thermo=thermo, thermo_model=model, thermo_background_m=CAROTENOID_BACKGROUND,
            latent=LatentState(0.02, "UPRE2", 4.0), mode="legacy",
            # The ninth layer, added 2026-09-02. kcat is ERG1's, the SLOWEST enzyme in the
            # host model's own isoprenoid branch, read out of the vendored ecYeastGEM rather
            # than asserted -- so this is the most conservative claim the audit can make and
            # it still disagrees with the fitted ceiling.
            # The content is a MASS fraction converted by the module that owns that
            # convention; passing a fraction straight in is what the 2026-09-04 unit fix
            # removed.
            enzyme_capacity={
                "enzyme": "crtYB", "kcat_per_s": 0.076,
                "enzyme_mmol_per_gdcw": enzyme_content_from_mass_fraction(
                    0.005, 74736.0)})
        assert isinstance(got, LegacyProductDiagnostic)
        assert got.supported is False
        got = got.calculation
        states = {layer.name: layer.state for layer in got.layers}

        assert LayerState.NOT_RUN not in states.values(), states
        assert states["GEM / FBA"] == LayerState.AUDITS
        assert states["enzyme capacity"] == LayerState.AUDITS
        assert states["regulation (E-Flux)"] == LayerState.INERT
        assert math.isfinite(got.content_mg_per_gdcw)

    def test_and_the_gem_audit_still_cannot_move_the_number(self, spec, model, thermo):
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID

        bare = _predict(spec)
        audited = _predict(spec, audit_model=model, audit_reaction=PRODUCT_DEMAND_ID,
                           audit_glucose_uptake=10.0)

        assert audited.content_mmol_per_gdcw == bare.content_mmol_per_gdcw

    def test_product_demand_audit_receives_terminal_production_not_total_entry_flux(self, spec, model):
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID

        bare = _predict(spec)
        terminal_rate = bare.solution.terminal.flux_in
        entry_rate = bare.flux.flux_mmol_per_gdcw_h
        assert terminal_rate < entry_rate
        assert terminal_rate == bare.rate_mmol_per_gdcw_h
        # This bound admits the solved product rate, but not the different entry
        # rate (which also supplies the accumulated lycopene intermediate).
        with model:
            model.reactions.get_by_id(PRODUCT_DEMAND_ID).upper_bound = (terminal_rate + entry_rate) / 2
            audited = _predict(spec, audit_model=model, audit_reaction=PRODUCT_DEMAND_ID,
                               audit_glucose_uptake=10.0)

        assert audited.flux_audit.predicted_flux == terminal_rate
        assert audited.flux_audit.within_envelope
        assert audited.content_mmol_per_gdcw == bare.content_mmol_per_gdcw
        assert audited.flux.flux_mmol_per_gdcw_h == entry_rate


class TestTheSpecRefusesChemistryItCannotGate:
    def test_a_step_without_the_nodes_own_species_is_refused(self):
        with pytest.raises(ValueError, match="no `metabolite`"):
            Node("lycopene", "saturating", reaction="CRTI")

    def test_a_step_written_backwards_is_refused(self):
        with pytest.raises(ValueError, match="direction that consumes it"):
            Node("lycopene", "saturating", metabolite="lycopene_c",
                 stoichiometry={"lycopene_c": -1.0, "phytoene_c": 1.0})

    def test_a_species_that_is_not_in_its_own_step_is_refused(self):
        with pytest.raises(ValueError, match="does not appear in its own step"):
            Node("lycopene", "saturating", metabolite="lycopene_c",
                 stoichiometry={"phytoene_c": -1.0, "lycopen_c": 1.0})

    def test_two_nodes_declared_to_be_the_same_species_are_refused(self):
        with pytest.raises(ValueError, match="both declared to be"):
            PathwaySpec(
                product="beta_carotene", organism="S. cerevisiae", entry_enzyme="CrtE",
                precursor_metabolite="s_0189",
                nodes=(Node("lycopene", "passthrough", metabolite="lycopene_c"),
                       Node("beta_carotene", "passthrough", metabolite="lycopene_c",
                            molar_mass_g_per_mol=536.87)))

    def test_a_node_declared_to_be_the_precursor_is_refused(self):
        """The precursor is what the pathway is pulled OFF, not a link in it. Declaring a
        node to be it makes the entry step consume and produce the same species, and Q
        would then carry that species twice with opposite signs."""
        with pytest.raises(ValueError, match="also the precursor"):
            PathwaySpec(
                product="beta_carotene", organism="S. cerevisiae", entry_enzyme="CrtE",
                precursor_metabolite="s_0189",
                nodes=(Node("ggpp", "passthrough", metabolite="s_0189"),
                       Node("beta_carotene", "passthrough", molar_mass_g_per_mol=536.87)))

    def test_an_empty_stoichiometry_is_refused_rather_than_treated_as_absent(self):
        with pytest.raises(ValueError, match="is empty"):
            Node("lycopene", "saturating", metabolite="lycopene_c", stoichiometry={})


class TestTheShippedSpecsDeclareWhatTheGateNeeds:
    def test_beta_carotene_gates_every_node(self, spec):
        assert set(spec.thermo_steps) == {n.name for n in spec.nodes}

    def test_phb_gates_the_one_step_the_tables_cover(self, phb_spec):
        assert set(phb_spec.thermo_steps) == {"acetoacetyl_coa"}

    def test_both_are_gateable(self, spec, phb_spec):
        assert spec.gateable and phb_spec.gateable

    def test_every_gated_node_also_declares_its_species(self, spec, phb_spec):
        for candidate in (spec, phb_spec):
            assert set(candidate.thermo_steps) <= set(candidate.thermo_metabolites)
