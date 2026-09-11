"""`predict_product(mode="mechanistic")`: the 109-state engine joined to the product call.

RED FILE. Every test that names `mode="mechanistic"` fails today, because today
`predict.py` line 1363 says ``mode must be 'empirical' or 'legacy'`` and there is no
`MechanisticRun` to hand it. That is the missing feature, not a typo.

The one class of test here that PASSES today is deliberate and is not dressed up as TDD:
:class:`TestTheEmpiricalDefaultDoesNotMove` is the RATCHET. It pins today's empirical
numbers as literals so that the change proves it moved no figure, and it must stay green
from the first run to the last.

Four things this file will not let the implementation do.

1. **Invent the number.** The mechanistic content must be the engine's OWN
   ``truth["content.beta_carotene"]``, compared with ``==`` against a run this file
   performs itself through `mech.engine.simulate`. Not "close to it".
2. **Return a tautology.** The seam is only worth having if it moves the number, so an
   axis the empirical route is provably blind to -- temperature at a held setpoint, where
   empirical returns 0.001994873775535095 at both 30 C and 37 C -- must move the
   mechanistic content by a real margin. 1.79x is what the engine does; the assertion
   demands 1.5x so it fails loudly if the axis is quietly dropped rather than passing on
   a float wobble.
3. **Acquire evidence it does not have.** 160 of the 173 provenance parameters behind
   this number are priors. ``validity.biological_validation`` is False and ``supported``
   is False, and the returned object is NOT a `ProductPrediction`, so nothing that quotes
   a supported claim can quote it.
4. **Default an unmeasured input.** No stressor dose becomes a medium composition, no
   read time is interpolated onto a grid the trajectory does not contain, and no
   protocol, inoculum or tolerance is supplied by the library.

Numbers quoted in this file were measured against the tree at d2be2b1 with
`/private/tmp/ystwin-quality-venv-20260908/bin/python`, not copied from a document.
"""

from __future__ import annotations

import pytest

from ystwin import predict as P
from ystwin.generator.context import CultureContext
from ystwin.mech import contracts, engine
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.spec import load_pathway
from ystwin.predict import (
    Environment,
    Genotype,
    LayerState,
    ProductPrediction,
    predict_product,
)

MEDIUM = {"glucose": 50.0, "nitrogen": 20.0, "oxygen": 0.2, "osmolyte": 250.0}
"""The declared reference medium. Every mM here is a caller declaration, not a default the
library supplies -- `MechanisticRun` has no medium of its own and must not grow one."""

TIMES_H = (0.0, 2.0, 4.0, 6.0)
SETPOINT = 0.18
"""Inside `BETA_CAROTENE_KINETICS["lycopene"]`'s fitted window [0.101, 0.2543], and the
setpoint every empirical pin below is held at."""

EMPIRICAL_CONTENT_MMOL_PER_GDCW = 0.001994873775535095
EMPIRICAL_CONTENT_MG_PER_GDCW = 1.0709878838715265
EMPIRICAL_RATE_MMOL_PER_GDCW_H = 0.00035907727959631715
"""Today's empirical answer at ``Environment(growth_rate_setpoint_per_h=0.18)``, to the
last bit. `test_predict_mechanistic.py::HELD_SETPOINT_CONTENT` carries the mg number too;
both are pinned so a mechanistic seam cannot move either without this file going red."""

EMPIRICAL_LAYERS = (
    ("metabolism", LayerState.IN_CHAIN, True),
    ("flux law", LayerState.IN_CHAIN, True),
    ("stress panel", LayerState.AUDITS, False),
    ("environment", LayerState.REPORTED, False),
    ("cassette dosage", LayerState.INERT, False),
    ("enzyme capacity", LayerState.NOT_RUN, False),
    ("GEM / FBA", LayerState.NOT_RUN, False),
    ("regulation (E-Flux)", LayerState.NOT_RUN, False),
    ("thermodynamics", LayerState.NOT_RUN, False),
    ("latent stress state", LayerState.NOT_RUN, False),
)


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


def _protocol(*, temperature_c: float = 30.0, ph: float = 5.0,
              times_h: tuple[float, ...] = TIMES_H) -> contracts.Protocol:
    return contracts.Protocol(
        times_h,
        (contracts.Control(0.0, temperature_c=temperature_c, ph=ph,
                           oxygen_transfer_per_h=30.0, oxygen_saturation_mM=0.2),),
        source="declared in-silico protocol for the mechanistic product seam",
    )


def _engine_inputs(*, temperature_c: float = 30.0, ph: float = 5.0, medium=None,
                   copies: dict[str, int] | None = None,
                   times_h: tuple[float, ...] = TIMES_H):
    """The four contracts the engine needs, all explicit. Nothing here is defaulted."""
    parameters = engine.EngineParameters.prior()
    genotype = contracts.Genotype.prior()
    if copies:
        genotype = genotype.with_copies(**copies)
    protocol = _protocol(temperature_c=temperature_c, ph=ph, times_h=times_h)
    initial = engine.initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=0.5,
                                medium_mM=dict(medium or MEDIUM))
    return protocol, genotype, initial, parameters


def _mechanistic_run(**kwargs):
    """The new input contract. `AttributeError` here IS the red: it does not exist yet."""
    read_time_h = kwargs.pop("read_time_h", None)
    protocol, genotype, initial, parameters = _engine_inputs(**kwargs)
    return P.MechanisticRun(protocol=protocol, genotype=genotype, initial_state=initial,
                            parameters=parameters, read_time_h=read_time_h)


def _mechanistic(spec, *, environment=None, **kwargs):
    return predict_product(
        spec, Genotype(1.0, "b-car4"), environment or Environment(),
        BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
        mode="mechanistic", mechanism=_mechanistic_run(**kwargs))


def _engine_truth(**kwargs) -> float:
    """The same simulation run straight through `mech.engine`, bypassing predict.py."""
    protocol, genotype, initial, parameters = _engine_inputs(**kwargs)
    result = engine.simulate(protocol, genotype, initial, parameters)
    return float(result.truth["content.beta_carotene"][-1])


def _empirical(spec, context=None, **kwargs):
    return predict_product(spec, Genotype(1.0, "b-car4"),
                           Environment(context or CultureContext(),
                                       growth_rate_setpoint_per_h=SETPOINT),
                           BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, **kwargs)


class TestTheNumberComesFromTheEnginesOwnPool:
    """(a) The content is the engine's beta-carotene pool, in mmol/gDCW, unconverted."""

    def test_the_content_is_the_engines_own_beta_carotene_trace(self, spec):
        """``==``, not ``approx``. A conversion layer would show up in the last bits, and
        `engine.py:865` already divides internal mmol by structural biomass gDW -- the
        same quantity and unit `ProductPrediction.content_mmol_per_gdcw` returns."""
        got = _mechanistic(spec)

        assert got.content_mmol_per_gdcw == _engine_truth()

    def test_the_content_is_a_real_positive_number_and_not_a_zeroed_stub(self, spec):
        got = _mechanistic(spec)

        assert got.content_mmol_per_gdcw > 0.0

    def test_the_milligram_number_is_the_specs_molar_mass_and_not_an_assumed_one(self, spec):
        got = _mechanistic(spec)

        assert got.content_mg_per_gdcw == (
            got.content_mmol_per_gdcw * spec.node("beta_carotene").molar_mass_g_per_mol)

    def test_the_whole_trajectory_comes_back_not_just_the_scalar(self, spec):
        """The caller must be able to refute the number, which needs the run behind it."""
        got = _mechanistic(spec)

        assert got.result.truth["content.beta_carotene"][-1] == got.content_mmol_per_gdcw
        assert tuple(got.result.times_h) == TIMES_H

    def test_the_rate_is_the_engines_production_flux_not_content_times_mu(self, spec):
        """``flux.lcy / biomass_gdw``. Off steady state the identity ``content * mu`` is
        wrong here by about 2.6x, so the flux is read and the identity is not used."""
        got = _mechanistic(spec)
        truth = got.result.truth

        assert got.production_rate_mmol_per_gdcw_h == (
            float(truth["flux.lcy"][-1]) / float(truth["biomass_gdw"][-1]))
        assert got.production_rate_mmol_per_gdcw_h != (
            got.content_mmol_per_gdcw * float(truth["growth_per_h"][-1]))


class TestTheSeamActuallyMovesTheNumber:
    """(b) An axis the empirical route cannot see, moving the mechanistic content."""

    def test_the_empirical_route_is_bit_blind_to_temperature_at_a_held_setpoint(self, spec):
        """The control for the next test. Both are 0.001994873775535095 exactly, so any
        temperature response the mechanistic route shows is the seam and nothing else."""
        cold = _empirical(spec, CultureContext(temperature_c=30.0))
        hot = _empirical(spec, CultureContext(temperature_c=37.0))

        assert cold.content_mmol_per_gdcw == hot.content_mmol_per_gdcw
        assert cold.content_mmol_per_gdcw == EMPIRICAL_CONTENT_MMOL_PER_GDCW

    def test_temperature_moves_the_mechanistic_content_by_a_real_margin(self, spec):
        """Measured 5.4958e-05 at 30 C against 3.0710e-05 at 37 C: 1.79x. The threshold is
        1.5x so a dropped axis fails loudly rather than passing on a float wobble."""
        cold = _mechanistic(spec, temperature_c=30.0).content_mmol_per_gdcw
        hot = _mechanistic(spec, temperature_c=37.0).content_mmol_per_gdcw

        assert hot < cold
        assert cold / hot > 1.5

    def test_the_carbon_feed_moves_the_mechanistic_content_by_a_real_margin(self, spec):
        """A second, independent axis. 50 mM glucose gives 5.4958e-05 and 5 mM gives
        7.6011e-05 -- less biomass, more per gDW -- a 1.38x move the empirical route,
        which refuses a non-default `context.glucose_g_per_L` outright, cannot express."""
        rich = _mechanistic(spec).content_mmol_per_gdcw
        lean = _mechanistic(spec, medium={**MEDIUM, "glucose": 5.0}).content_mmol_per_gdcw

        assert lean != rich
        assert max(lean, rich) / min(lean, rich) > 1.25

    def test_the_cassette_copies_move_it_too_because_the_engine_reads_copies(self, spec):
        """crtYB 1 -> 3 moves it 5.4958e-05 -> 1.3287e-04, 2.42x. This is the question
        `pathway/capacity.py` refuses for lack of a cyclase coefficient; the engine does
        not need one, because it integrates the gene rather than scaling a capacity."""
        one = _mechanistic(spec).content_mmol_per_gdcw
        three = _mechanistic(spec, copies={"crtyb": 3}).content_mmol_per_gdcw

        assert three / one > 2.0

    def test_medium_ph_is_inert_rather_than_in_chain_when_no_acetate_is_present(self, spec):
        """The reason the environment rows are decided per RUN and not per NAME. pH 4 and
        pH 5 are bit-identical without acetate, so claiming pH IN_CHAIN would be false."""
        five = _mechanistic(spec, ph=5.0)
        four = _mechanistic(spec, ph=4.0)

        assert four.content_mmol_per_gdcw == five.content_mmol_per_gdcw
        assert five.layer("engine medium pH").state == LayerState.INERT


class TestTheEmpiricalDefaultDoesNotMove:
    """(c) THE RATCHET. Green today, green after. Not a red test, and not TDD -- this is
    the proof that no existing figure moved, and it is stated as such."""

    def test_the_default_mode_is_still_empirical(self, spec):
        assert _empirical(spec).mode == "empirical"

    def test_the_pinned_empirical_numbers_are_unchanged_to_the_last_bit(self, spec):
        got = _empirical(spec)

        assert got.content_mmol_per_gdcw == EMPIRICAL_CONTENT_MMOL_PER_GDCW
        assert got.content_mg_per_gdcw == EMPIRICAL_CONTENT_MG_PER_GDCW
        assert got.rate_mmol_per_gdcw_h == EMPIRICAL_RATE_MMOL_PER_GDCW_H

    def test_naming_the_default_explicitly_returns_the_identical_object(self, spec):
        """Whole-dataclass ``==``, which `ProductPrediction` supports."""
        assert _empirical(spec) == _empirical(spec, mode="empirical")

    def test_the_empirical_layer_report_is_unchanged_row_for_row(self, spec):
        got = _empirical(spec)

        assert tuple((row.name, row.state, row.sets_the_number) for row in got.layers) == EMPIRICAL_LAYERS
        assert len(got.notes) == 4

    def test_no_engine_row_leaks_into_an_empirical_call(self, spec):
        assert not [row for row in _empirical(spec).layers if row.name.startswith("engine ")]

    def test_the_new_argument_cannot_alter_an_empirical_call_even_by_accident(self, spec):
        """`mechanism=None` must be exactly the bare call, and a mechanism handed to the
        empirical route must be refused rather than silently switching modes."""
        assert _empirical(spec) == _empirical(spec, mechanism=None)
        with pytest.raises(ValueError, match="mechanistic"):
            _empirical(spec, mechanism=_mechanistic_run())


class TestTheResultCannotBeQuotedAsSupported:
    """(d) The evidence axis stays false however wired the dataflow is."""

    def test_the_engines_validity_is_carried_and_says_not_experimentally_validated(self, spec):
        got = _mechanistic(spec)

        assert got.validity.biological_validation is False
        assert got.validity.parameter_basis == "prior-conditional; not experimentally validated"

    def test_supported_is_false_and_cannot_be_constructed_otherwise(self, spec):
        got = _mechanistic(spec)

        assert got.supported is False
        with pytest.raises((TypeError, AttributeError)):
            type(got)(**{**vars(got), "supported": True})

    def test_it_is_not_a_product_prediction_so_no_claim_consumer_can_quote_it(self, spec):
        """`scripts/register_prediction.py` writes `data/current_claims.json` from a
        `ProductPrediction`. The type is the third lock, after `supported` and validity."""
        got = _mechanistic(spec)

        assert not isinstance(got, ProductPrediction)
        assert got.mode == "mechanistic"

    def test_the_prior_parameter_count_is_computed_from_the_run_not_hardcoded(self, spec):
        got = _mechanistic(spec)

        assert got.prior_parameters == len(got.result.prior_parameters)
        assert got.total_parameters == len(got.result.provenance)
        assert got.prior_parameters > 0.9 * got.total_parameters

    def test_the_summary_says_prior_conditional_before_it_says_any_number(self, spec):
        assert _mechanistic(spec).summary().startswith("MECHANISTIC")


class TestTheMechanisticLayersReportTheirState:
    """(d, second half) Dataflow claims, each falsifiable from this run's own traces."""

    @pytest.mark.parametrize("name", [
        "engine metabolism", "engine expression", "engine vessel",
    ])
    def test_the_rows_that_set_the_number_say_so(self, spec, name):
        row = _mechanistic(spec).layer(name)

        assert row.state == LayerState.IN_CHAIN
        assert row.sets_the_number is True

    def test_the_conservation_row_audits_and_never_moves_the_number(self, spec):
        row = _mechanistic(spec).layer("engine conservation")

        assert row.state == LayerState.AUDITS
        assert row.sets_the_number is False

    def test_the_validity_row_is_reported_and_names_the_evidence_gap(self, spec):
        row = _mechanistic(spec).layer("engine validity")

        assert row.state == LayerState.REPORTED
        assert row.sets_the_number is False

    def test_eflux_stays_not_run_and_says_what_would_run_it(self, spec):
        """The reference-condition contract is declared; the gene-to-reaction map in the
        installed GEM is not held, so the row stays NOT_RUN and names what is missing."""
        row = _mechanistic(spec).layer("regulation (E-Flux)")

        assert row.state == LayerState.NOT_RUN
        assert "pass " in row.detail

    def test_the_skipped_empirical_stack_is_visible_rather_than_absent(self, spec):
        """A mechanistic report stays comparable to an empirical one line by line."""
        got = _mechanistic(spec)

        for name in ("metabolism", "flux law", "stress panel", "GEM / FBA"):
            assert got.layer(name).state == LayerState.NOT_RUN

    def test_every_not_run_row_still_says_what_would_run_it(self, spec):
        for row in _mechanistic(spec).layers:
            if row.state == LayerState.NOT_RUN:
                assert "pass " in row.detail or "mode=" in row.detail, row.name


class TestAnUnmeasuredInputIsRefusedNotDefaulted:
    """(e) Every place the seam would have to invent a number, it says so instead."""

    def test_mechanistic_mode_without_a_mechanism_names_what_is_missing(self, spec):
        """No protocol, inoculum or medium is defaulted, so the mode cannot run bare."""
        with pytest.raises(ValueError, match="MechanisticRun"):
            predict_product(spec, Genotype(1.0, "b-car4"), Environment(),
                            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS, mode="mechanistic")

    def test_a_stressor_dose_is_refused_because_no_calibration_maps_it_to_a_medium(self, spec):
        """The sharpest refusal in the file. Nothing maps "DTT 10" onto a temperature, a
        pH or a feed composition, so the Environment's culture claim is rejected rather
        than reconciled with the Protocol's."""
        with pytest.raises(ValueError, match="dose"):
            _mechanistic(spec, environment=Environment(stressor="DTT", dose=10.0))

    def test_a_growth_setpoint_is_refused_because_the_protocol_already_sets_the_vessel(self, spec):
        with pytest.raises(ValueError):
            _mechanistic(spec, environment=Environment(growth_rate_setpoint_per_h=SETPOINT))

    def test_a_read_time_off_the_protocol_grid_is_refused_rather_than_interpolated(self, spec):
        """The trajectory exists only at ``protocol.times_h``. Interpolating would return
        a number the returned run does not contain."""
        with pytest.raises(ValueError, match="times_h"):
            _mechanistic(spec, read_time_h=3.0)

    def test_a_read_time_on_the_grid_is_honoured(self, spec):
        got = _mechanistic(spec, read_time_h=4.0)

        assert got.read_time_h == 4.0
        assert got.content_mmol_per_gdcw == float(
            got.result.truth["content.beta_carotene"][TIMES_H.index(4.0)])

    def test_the_engine_contracts_have_no_defaults_to_fall_back_on(self, spec):
        with pytest.raises(TypeError):
            P.MechanisticRun()

    def test_a_legacy_input_type_is_refused_rather_than_completed_with_priors(self, spec):
        """`mech/adapters.py` already refuses this class of conversion; the seam reuses
        that refusal instead of growing a second, weaker one."""
        protocol, genotype, initial, parameters = _engine_inputs()
        with pytest.raises(ValueError, match="legacy inputs are not converted"):
            predict_product(
                spec, Genotype(1.0, "b-car4"), Environment(), BETA_CAROTENE_FLUX,
                BETA_CAROTENE_KINETICS, mode="mechanistic",
                mechanism=P.MechanisticRun(protocol=protocol, genotype=Genotype(1.0),
                                           initial_state=initial, parameters=parameters))

    @pytest.mark.parametrize("argument", [
        "audit_model", "thermo", "latent", "enzyme_capacity", "mech",
    ])
    def test_an_argument_that_audits_a_pathway_solution_is_refused(self, spec, argument):
        """Mechanistic mode produces no `PathwaySolution`, so these have nothing to audit.
        Refused by name rather than ignored, which would be a silent no-op."""
        with pytest.raises(ValueError, match=argument):
            predict_product(
                spec, Genotype(1.0, "b-car4"), Environment(), BETA_CAROTENE_FLUX,
                BETA_CAROTENE_KINETICS, mode="mechanistic",
                mechanism=_mechanistic_run(), **{argument: object()})

    def test_an_engine_refusal_propagates_rather_than_becoming_a_number(self, spec):
        """The engine's temperature domain, reached through the product call. A numeric
        fallback here would be the seam inventing biology the engine refused to state."""
        with pytest.raises((ValueError, contracts.ScientificRefusal)):
            _mechanistic(spec, temperature_c=60.0)
