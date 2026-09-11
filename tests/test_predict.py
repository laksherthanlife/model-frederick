"""The assembled chain: environment and genotype in, product out.

This file was rewritten when the chain was assembled, and the two things it now cannot test
are the point of the rewrite.

It can no longer pass a pathway flux in, because there is no longer a way to. `Environment`
used to carry `pathway_flux`, and both entry points obtained it as
``q_lycopene + q_betacarotene`` from the very state being predicted -- the model measured
the product in order to predict it. Flux now comes from the genotype.

And it can no longer assert that a stressor lowers the product rate in a chemostat, because
that behaviour was wrong. The old chain set ``growth = D`` and then multiplied by a
viability fraction, reporting mu < D at a steady state where mu = D by construction. What a
stressor decides in a chemostat is whether a steady state exists at all.
"""

from __future__ import annotations

import pytest

from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, SetpointUnreachable, predict_product


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def product_validation_module():
    import importlib.util
    import pathlib

    path = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "predict_product.py"
    loader = importlib.util.spec_from_file_location("predict_product", path)
    module = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(module)
    return module


def run(spec, expression=1.0, **environment):
    return predict_product(spec, Genotype(expression), Environment(**environment),
                           BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS)


class TestItRunsTheWholeChain:
    def test_a_genotype_and_an_environment_become_a_product_rate(self, spec):
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert got.rate_mmol_per_gdcw_h > 0

    def test_the_flux_is_predicted_and_not_supplied(self, spec):
        """The whole rewrite in one assertion. There is no way to hand a flux in."""
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert got.flux.flux_mmol_per_gdcw_h > 0

    def test_every_intermediate_comes_out_rather_than_going_in(self, spec):
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert set(got.intermediates()) == {"phytoene", "lycopene"}

    def test_more_entry_expression_gives_more_product(self, spec):
        low = run(spec, expression=0.3, growth_rate_setpoint_per_h=0.15)
        high = run(spec, expression=1.0, growth_rate_setpoint_per_h=0.15)

        assert high.rate_mmol_per_gdcw_h > low.rate_mmol_per_gdcw_h

    def test_the_content_uses_the_molar_mass_the_spec_declares(self, spec):
        got = run(spec, growth_rate_setpoint_per_h=0.15)
        declared = spec.node("beta_carotene").molar_mass_g_per_mol

        assert got.content_mg_per_gdcw == pytest.approx(
            got.content_mmol_per_gdcw * declared)

    def test_the_rate_is_the_content_washed_out_by_growth(self, spec):
        """At steady state a non-secreted product leaves only by dilution, so these are
        the same statement and any drift between them is an arithmetic slip."""
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert got.rate_mmol_per_gdcw_h == pytest.approx(
            got.content_mmol_per_gdcw * 0.15)


class TestASetpointIsASetpoint:
    """mu = mu_set is the definition of the held state, and the previous chain broke it.

    It set ``growth = D`` and then multiplied by a viability fraction, so a stressed culture
    reported mu below the rate the vessel was holding -- a state that does not exist. A
    stressor moves residual substrate and biomass.

    THE VESSEL CHANGED UNDER THIS CLASS AND THE ARITHMETIC DID NOT. What was a chemostat
    dilution rate is now an exponential fed-batch setpoint, and every assertion about mu
    below is unchanged, which is the substantive claim: the pools law divides by mu either
    way. What DID change is the failure at the top end. A chemostat pushed past its critical
    dilution rate EMPTIES; a fed-batch has no effluent and instead ACCUMULATES the carbon the
    culture can no longer take up. Both are "there is no steady state here", and the refusal
    is now `SetpointUnreachable`.
    """

    def test_an_unstressed_culture_grows_at_the_setpoint(self, spec):
        assert run(spec, growth_rate_setpoint_per_h=0.17).growth_rate_per_h == pytest.approx(0.17)

    def test_a_survivable_dose_leaves_the_growth_rate_alone(self, spec):
        stressed = run(spec, growth_rate_setpoint_per_h=0.15, stressor="DTT", dose=1.0)

        assert stressed.growth_rate_per_h == pytest.approx(0.15)

    def test_and_therefore_leaves_the_product_alone(self, spec):
        """The consequence, stated so it cannot be quietly reintroduced. Stress has no route
        to the product under a held growth rate in this architecture, and that is a claim
        about held growth rates rather than a gap: the feed fixes the only channel there
        was. It was true of the pump and it is equally true of the setpoint."""
        plain = run(spec, growth_rate_setpoint_per_h=0.15)
        stressed = run(spec, growth_rate_setpoint_per_h=0.15, stressor="DTT", dose=1.0)

        assert stressed.rate_mmol_per_gdcw_h == pytest.approx(plain.rate_mmol_per_gdcw_h)

    def test_and_it_says_so_rather_than_letting_the_silence_speak(self, spec):
        """The test above establishes that stress cannot move the product in a chemostat.
        This one is about what that looks like to a reader.

        A user comparing a stressed and an unstressed prediction sees two identical numbers.
        Without a note, the honest arithmetic reads as the unearned claim "stress does not
        affect production". Two routes are plausible and neither is modelled: a
        maintenance-ATP burden (`bridge/latent_bridge.py`, refused by G4) and flux
        redirection (E-Flux scales upper bounds, and an upper bound cannot lower a flux that
        is below its ceiling). "No route modelled" is the claim; "no effect" is not."""
        stressed = run(spec, growth_rate_setpoint_per_h=0.15, stressor="DTT", dose=1.0)

        assert any("no route modelled" in note for note in stressed.notes)

    def test_nor_on_an_unstressed_chemostat(self, spec):
        plain = run(spec, growth_rate_setpoint_per_h=0.15)

        assert not any("no route modelled" in note for note in plain.notes)

    def test_a_dose_that_pushes_mu_max_below_the_setpoint_is_refused(self, spec):
        with pytest.raises(SetpointUnreachable, match="carbon it cannot take up accumulates"):
            run(spec, growth_rate_setpoint_per_h=0.30, stressor="DTT", dose=5.0)

    def test_the_refusal_message_names_both_rates(self, spec):
        """A refusal a reader cannot act on is an obstacle. Which rate was too low, and
        against what, is the whole content."""
        with pytest.raises(SetpointUnreachable,
                           match=r"maximum growth rate of .*below the setpoint"):
            run(spec, growth_rate_setpoint_per_h=0.30, stressor="DTT", dose=5.0)

    def test_in_batch_a_stressor_does_slow_growth(self, spec):
        """The other half, tested one layer down because the chain itself now refuses batch.

        In batch mu is what the cells can do rather than what a pump imposes, so a dose that
        is inert in a chemostat is not inert here -- but batch mu is around 0.4 /h, outside
        the range the branch kinetics were fitted over, so `predict_product` refuses before
        it can show that. The growth layer is still right and is checked directly."""
        from ystwin.generator.context import CultureContext, context_growth_rate
        from ystwin.generator.stress_panel import growth_rate as panel_growth

        unstressed = panel_growth("DTT", 0.0)
        stressed = panel_growth("DTT", 1.0)
        maximum = context_growth_rate(CultureContext())

        assert maximum * (stressed / unstressed) < maximum


class TestBatchIsRefusedAndThatIsTheHonestAnswer:
    """The most common question a user has -- "how much do I get in a shake flask?" -- is
    the one this model cannot answer.

    Batch mu is around 0.4 /h. The branch kinetics were fitted over [0.101, 0.254], a
    2.5-fold window that batch sits well outside, and the growth exponent inside that window
    already has a 95% interval of [0.53, 1.78]. Extrapolating a law that loose by 60% beyond
    its range is not a small error.

    This was silently answered until the range guard was restored -- the generic solver had
    dropped a refusal the hand-written one had. `docs/WHAT_IS_LEFT.md` lists a batch
    calibration state as the fix, and it is a measurement rather than an architecture change.
    """

    def test_a_batch_environment_is_refused(self, spec):
        with pytest.raises(ValueError, match="outside the range"):
            run(spec)

    def test_the_refusal_says_what_would_lift_it(self, spec):
        with pytest.raises(ValueError, match="Measure at this rate"):
            run(spec)

    def test_a_chemostat_inside_the_window_still_answers(self, spec):
        """So the refusal is about the rate and not about the mode."""
        assert run(spec, growth_rate_setpoint_per_h=0.18).rate_mmol_per_gdcw_h > 0

    def test_the_third_rate_you_would_measure_is_refused_too(self, spec):
        """D = 0.05 and D = 0.35 are the rates worth measuring next, because a point between
        the anchors moves the growth-exponent interval from 0.46 to 0.47 and a point outside
        halves it. The model refuses both, which is exactly why they are informative."""
        for rate in (0.05, 0.35):
            with pytest.raises(ValueError, match="outside the range"):
                run(spec, growth_rate_setpoint_per_h=rate)


class TestItRefusesRatherThanInventing:
    def test_zero_entry_expression_is_not_a_small_flux(self, spec):
        with pytest.raises(ValueError, match="does not carry the pathway"):
            run(spec, expression=0.0, growth_rate_setpoint_per_h=0.15)

    def test_an_unknown_stressor_is_refused(self, spec):
        with pytest.raises(KeyError, match="no stressor"):
            run(spec, growth_rate_setpoint_per_h=0.15, stressor="not-a-stressor", dose=1.0)

    def test_a_dose_with_no_stressor_is_refused(self, spec):
        with pytest.raises(ValueError, match="no stressor to apply it to"):
            run(spec, growth_rate_setpoint_per_h=0.15, dose=1.0)

    def test_a_non_positive_dilution_rate_is_refused(self, spec):
        with pytest.raises(ValueError, match="growth_rate_setpoint_per_h must be positive"):
            run(spec, growth_rate_setpoint_per_h=0.0)

    def test_a_calibration_for_another_gene_is_refused(self, spec):
        """A units error that would otherwise surface as a wrong number. The scalar is
        relative to one gene's normalisation and means nothing against another's."""
        import dataclasses

        wrong = dataclasses.replace(BETA_CAROTENE_FLUX, entry_enzyme="CrtI")

        with pytest.raises(ValueError, match="fitted on 'CrtI'"):
            predict_product(spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.15),
                            wrong, BETA_CAROTENE_KINETICS)


class TestHistoricallySelectedFixedCrtEComparison:
    """Historical fixed-gene comparison, not a prespecified or default nested validation.
    The original numerical bounds are retained only for this explicitly selected model."""

    @pytest.fixture(scope="class")
    def fixed_comparison(self, product_validation_module, spec):
        return product_validation_module.score_forward(
            spec, product_validation_module.measured(), mode="fixed_crte_comparison", draws=0)

    def test_all_six_fixed_comparison_states_are_retained_and_labelled(self, fixed_comparison):
        assert len(fixed_comparison) == 6
        assert set(fixed_comparison.validation_mode) == {"fixed_crte_comparison"}
        assert set(fixed_comparison.selected_model) == {"CrtE:rate:saturating"}
        assert fixed_comparison.strain.nunique() == 3

    def test_fixed_comparison_product_median_error_is_under_a_fifth(self, fixed_comparison):
        fold = fixed_comparison.product_predicted / fixed_comparison.product_measured

        assert float((fold - 1).abs().median()) < 0.20

    def test_fixed_comparison_has_no_state_out_by_more_than_half(self, fixed_comparison):
        fold = fixed_comparison.product_predicted / fixed_comparison.product_measured

        assert float(fold.max()) < 1.5 and float(fold.min()) > 1 / 1.5

    def test_fixed_comparison_intermediate_is_predicted_too(self, fixed_comparison):
        """Both pools, because product alone leaves the flux scalar unidentifiable over
        five orders of magnitude -- see PathwaySpec.calibratable."""
        fold = fixed_comparison.lycopene_predicted / fixed_comparison.lycopene_measured

        assert float((fold - 1).abs().median()) < 0.50

    def test_fixed_comparison_flux_tracks_the_measured_one(self, fixed_comparison):
        """The layer that replaced the circular input, scored on its own terms."""
        fold = fixed_comparison.flux_predicted / fixed_comparison.flux_measured

        assert float((fold - 1).abs().median()) < 0.25


class TestDefaultNestedProductValidation:
    @pytest.fixture(scope="class")
    def nested(self, product_validation_module, spec):
        return product_validation_module.score_forward(
            spec, product_validation_module.measured(), draws=0)

    def test_default_evaluates_three_independent_strains_and_all_six_conditions(self, nested):
        assert set(nested.validation_mode) == {"nested"}
        assert nested.strain.nunique() == 3
        assert nested.state.nunique() == len(nested) == 6
        assert set(nested.n_independent_strains) == {3}
        assert set(nested.n_condition_predictions) == {6}
        assert set(nested.n_channels) == {2}
        for strain, fold in nested.groupby("strain"):
            assert set(fold.training_strains) == {"|".join(sorted(set(nested.strain) - {strain}))}
            assert set(fold.n_training_strains) == {2}

    def test_default_reports_the_observed_errors_not_the_favorable_fixed_comparison(self, nested):
        product = nested.product_predicted / nested.product_measured
        lycopene = nested.lycopene_predicted / nested.lycopene_measured
        entry = nested.flux_predicted / nested.flux_measured
        observed = {
            "median_product_error": float((product - 1).abs().median()),
            "maximum_product_ratio": float(product.max()),
            "median_lycopene_error": float((lycopene - 1).abs().median()),
            "median_flux_error": float((entry - 1).abs().median()),
        }
        assert observed == pytest.approx({
            "median_product_error": 0.23505445933856173,
            "maximum_product_ratio": 2.1256714982848735,
            "median_lycopene_error": 0.6718116472596215,
            "median_flux_error": 0.41814001794225114,
        }, rel=0.0, abs=1e-6)
        assert set(nested.selected_gene) == {"CrtI", "CrtYB"}

    def test_default_retains_an_unsupported_heldout_rate_in_the_denominator(
            self, product_validation_module, spec):
        import numpy as np

        from ystwin.pathway.flux import summarize_product_validation

        states = product_validation_module.measured()
        states.loc[states.condition == "4D025", "mu_per_h"] = 0.5
        result = product_validation_module.score_forward(spec, states, draws=0)
        failed = result[result.state == "4D025"].iloc[0]

        assert len(result) == 6
        assert result.strain.nunique() == 3
        assert failed.status == "failed"
        assert "outside the range" in failed.failure_reason
        assert np.isnan(failed.product_predicted)
        assert np.isinf(failed.joint_log_mse)
        summary = summarize_product_validation(result, comparisons=False).iloc[0]
        assert summary.n_independent_strains == 3
        assert summary.n_condition_predictions == 6
        assert summary.n_failed_predictions == int(result.status.ne("ok").sum())
        assert summary.n_failed_predictions >= 1
        assert np.isinf(summary.joint_rmse_log)

    def test_default_retains_all_conditions_when_one_outcome_channel_is_absent(
            self, product_validation_module, spec):
        import numpy as np

        from ystwin.pathway.flux import summarize_product_validation

        states = product_validation_module.measured().drop(columns="q_lycopene")
        result = product_validation_module.score_forward(spec, states, draws=0)
        summary = summarize_product_validation(result, comparisons=False).iloc[0]

        assert len(result) == 6
        assert result.strain.nunique() == 3
        assert result.status.ne("ok").all()
        assert result.failure_reason.str.len().gt(0).all()
        assert summary.n_condition_predictions == summary.n_failed_predictions == 6
        assert summary.n_independent_strains == 3
        assert np.isinf(summary.joint_rmse_log)


class TestTheContextAuditsReachabilityAndSaysSo:
    """`Environment.context` is an accepted input that reaches nothing in a chemostat.

    Found by an adversarial pass: two calls differing only in `carbon_source` returned
    bit-identical numbers AND identical notes, so nothing distinguished "the model says the
    carbon source does not matter" from "the model cannot see the carbon source".

    Kocharin's PHB chemostats measure the difference directly -- at D = 0.05 the feed alone
    moves the product flux 4.24x with mu held identical, which is the largest single effect
    in the only environment series this repository has. So the discard is a real gap and the
    silence about it was the part worth fixing.

    The number is still identical, deliberately. Inventing a carbon-source term to make it
    differ is exactly what `tests/test_phb_environment.py` shows the data cannot support.
    """

    def test_two_carbon_sources_agree_only_where_both_can_hold_the_setpoint(self, spec):
        """Unchanged, and it must stay unchanged until a term is FITTED rather than added."""
        from ystwin.generator.context import CultureContext

        def go(source, setpoint=0.12):
            return predict_product(
                spec, Genotype(1.0),
                Environment(context=CultureContext(carbon_source=source), growth_rate_setpoint_per_h=setpoint),
                BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS)

        assert go("glucose").content_mg_per_gdcw == go("ethanol").content_mg_per_gdcw
        assert go("glucose", 0.15).growth_rate_per_h == 0.15
        with pytest.raises(SetpointUnreachable):
            go("ethanol", 0.15)

    def test_the_note_names_the_reachability_check(self, spec):
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert any("AUDITS REACHABILITY" in note for note in got.notes)

    def test_the_note_carries_the_actual_values_so_two_calls_differ(self, spec):
        """A constant caveat is one nobody reads. Naming the discarded values means a diff
        of two predictions shows what changed even though the number did not."""
        from ystwin.generator.context import CultureContext

        def notes(source):
            return predict_product(
                spec, Genotype(1.0),
                Environment(context=CultureContext(carbon_source=source), growth_rate_setpoint_per_h=0.12),
                BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS).notes

        assert notes("glucose") != notes("ethanol")

    def test_it_points_at_the_measurement_that_makes_it_a_gap(self, spec):
        """4.24x is not a hedge, it is a number from a real experiment, and the note carries
        it so nobody has to take the caveat on trust."""
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert any("4.24x" in note for note in got.notes)


class TestThePhysiologyDescribesADifferentCulture:
    """`chemostat_physiology` reads van Hoek 1998's GLUCOSE-LIMITED chemostat, and the six
    states this pathway was calibrated on were glucose-EXCESS.

    Measured against Elizondo's own columns: glucose uptake is 4.1-6.4x the table's value at
    the same mu, and every state makes 6.0-16.1 mmol/gDCW/h of ethanol -- including at
    mu = 0.101, where a glucose-limited culture makes none at all.

    The consequence reaches the public API. `ProductPrediction.fermentative` is derived from
    mu alone against a critical rate of 0.28, so it reports False for every calibration
    state while those cultures were visibly fermenting. The flag is kept because it is right
    for the cultures van Hoek measured and two consumers depend on it; what was missing was
    anything saying it is wrong for these.
    """

    def test_the_mismatch_is_real_and_large(self):
        import pandas as pd

        from ystwin import paths
        from ystwin.generator.culture import chemostat_physiology

        path = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"
        if not path.exists():
            pytest.skip(f"{path} not present")
        measured = pd.read_csv(path, sep="\t")

        ratios = [abs(r.q_glucose) / chemostat_physiology(r.mu_per_h).glucose_uptake
                  for r in measured.itertuples()]

        assert min(ratios) > 4.0

    def test_every_calibration_state_was_fermenting(self):
        """Which is the categorical half. The table says zero ethanol below the critical
        rate; all six states produce it, at both rates."""
        import pandas as pd

        from ystwin import paths

        path = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"
        if not path.exists():
            pytest.skip(f"{path} not present")
        measured = pd.read_csv(path, sep="\t")

        assert (measured.q_ethanol > 2.0).all()

    def test_but_the_flag_says_respiratory_at_the_low_rate(self):
        """The flag is not wrong about van Hoek's cultures. It is wrong about these, and
        both facts have to be visible at once."""
        from ystwin.generator.culture import chemostat_physiology

        assert not chemostat_physiology(0.101).fermentative

    def test_and_every_prediction_now_says_which_culture_it_describes(self, spec):
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert any("GLUCOSE-LIMITED" in note for note in got.notes)

    def test_the_note_carries_the_size_of_the_mismatch(self, spec):
        """A caveat without a number is a hedge. 4.1-6.4x is a measurement."""
        got = run(spec, growth_rate_setpoint_per_h=0.15)

        assert any("6.4x" in note for note in got.notes)


class TestTheGEMAuditIsReachableAndCannotChangeTheAnswer:
    """Item 2 on the ranked list, and the reason it was on the list at all.

    For a while `predict.py`'s docstring said the FBA layer was in the chain when
    `grep '^from' src/ystwin/predict.py` returned no `ystwin.fba`. The audit existed, was
    tested, and was reached only from a script behind a flag.

    Two honest options were available -- wire it, or stop claiming it -- and the reason it
    had not been wired is real: loading a GSMM costs an 18-second SBML parse, and the audit
    CANNOT CHANGE THE NUMBER. An upper bound cannot make a flux mandatory, and capping every
    pathway reaction at the measured magnitude leaves the FVA floor at exactly zero. Making
    every prediction pay 18 seconds to learn nothing new would be a poor trade.

    So it takes an ALREADY-LOADED model. The caller pays the parse once and reuses it; a
    caller who does not want the audit pays nothing; and the docstring is true either way.
    """

    @pytest.fixture(scope="class")
    def audited_model(self):
        import warnings

        from ystwin import paths
        from ystwin.fba.carotenoid import add_beta_carotene_pathway
        from ystwin.fba.solver import load_model

        gem = paths.yeast_gem()
        if gem is None or not gem.is_file():
            pytest.skip("Yeast9 not present; set YSTWIN_YEAST_GEM")
        warnings.simplefilter("ignore")
        model, _ = load_model(gem)
        return add_beta_carotene_pathway(model)

    def test_without_a_model_nothing_is_audited_and_nothing_is_parsed(self, spec):
        got = run(spec, growth_rate_setpoint_per_h=0.18)

        assert not any("GEM audit" in note for note in got.notes)

    def test_with_a_model_the_audit_appears(self, spec, audited_model):
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID

        got = predict_product(
            spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.18),
            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
            audit_model=audited_model, audit_reaction=PRODUCT_DEMAND_ID,
            audit_glucose_uptake=10.0)

        assert any("GEM audit" in note for note in got.notes)

    def test_and_the_predicted_number_is_bit_for_bit_unchanged(self, spec, audited_model):
        """The property that makes the audit an audit. If it moved the number it would be a
        second predictor, and the whole argument for demoting FBA is that it cannot be."""
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID

        plain = run(spec, growth_rate_setpoint_per_h=0.18)
        audited = predict_product(
            spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.18),
            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
            audit_model=audited_model, audit_reaction=PRODUCT_DEMAND_ID,
            audit_glucose_uptake=10.0)

        assert audited.content_mg_per_gdcw == plain.content_mg_per_gdcw

    def test_a_model_without_its_reaction_id_is_refused(self, spec, audited_model):
        """This function cannot know which of a GSMM's several thousand reactions is the
        product, and guessing would be worse than asking."""
        with pytest.raises(ValueError, match="audit_reaction"):
            predict_product(spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.18),
                            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
                            audit_model=audited_model)

    def test_and_so_is_one_without_an_uptake(self, spec, audited_model):
        """Required rather than defaulted because the growth cost the audit reports swings
        sevenfold with it -- 0.09% at a bound of 10, 0.64% at 1.5, same flux."""
        from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID

        with pytest.raises(ValueError, match="audit_glucose_uptake"):
            predict_product(spec, Genotype(1.0), Environment(growth_rate_setpoint_per_h=0.18),
                            BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
                            audit_model=audited_model, audit_reaction=PRODUCT_DEMAND_ID)

    def test_the_docstring_no_longer_claims_more_than_the_code_does(self):
        """The defect this whole class exists for. The diagram row must say opt-in, and the
        import must be real."""
        import pathlib

        source = (pathlib.Path(__file__).resolve().parents[1] / "src" / "ystwin"
                  / "predict.py").read_text()

        assert "AUDITS, opt-in" in source
        assert "from .fba.audit import audit_predicted_flux" in source


class TestHistoricalContentBandComparisons:
    """The band omitted the term that dominates it, and the docstring said so in prose
    while the code kept returning the narrow number.

    Error budget for predicted content, each parameter across its own published 95%
    interval: the entry-flux scalar is 0.4-6.1% of the log-variance, the cyclase capacity
    71.5-96.2%. `content_interval` propagated only the first. These tests pin that supplying
    `kinetic_calibration` fixes it, that the two parameters are propagated JOINTLY rather
    than added, and the size of what remains wrong.
    """

    @staticmethod
    def _band(mu, calibration=None):
        from ystwin.kinetic.carotenoid import ELIZONDO2025
        from ystwin.pathway import calibrations
        from ystwin.pathway.spec import load_pathway
        from ystwin.predict import Environment, Genotype, predict_product

        spec = load_pathway("beta_carotene")
        got = predict_product(
            spec, Genotype(1.0, "b-car4"),
            Environment(growth_rate_setpoint_per_h=mu),
            calibrations.BETA_CAROTENE_FLUX, calibrations.BETA_CAROTENE_KINETICS)
        kwargs = {"kinetic_calibration": ELIZONDO2025} if calibration else {}
        low, high = got.content_interval(
            spec, calibrations.BETA_CAROTENE_KINETICS, **kwargs)
        return high / low

    @pytest.mark.parametrize("mu", [0.101, 0.150, 0.200, 0.254])
    def test_the_kinetics_widen_the_band_at_every_growth_rate(self, mu):
        assert self._band(mu, calibration=True) > self._band(mu)

    def test_the_omission_was_worst_where_the_repo_predicts(self):
        """9.2x too narrow as log width at Elizondo's lower calibration state. The band was
        least honest exactly where it was used most."""
        import math

        narrow = math.log(self._band(0.101))
        honest = math.log(self._band(0.101, calibration=True))

        assert honest / narrow > 9.0

    def test_the_two_parameters_are_propagated_jointly_and_not_added(self):
        """log capacity and log km correlate at +0.852 and push content in OPPOSITE
        directions, so adding the variances overstates. If this ever equals the independent
        propagation, the covariance has stopped being used."""
        from ystwin.kinetic.carotenoid import ELIZONDO2025

        gradient = (1.0, -0.5)
        joint = ELIZONDO2025.joint_log_sd(gradient)
        independent = ELIZONDO2025.independent_log_sd(gradient)

        assert joint < independent
        assert ELIZONDO2025.log_parameter_correlation > 0.8

    def test_it_refuses_a_pathway_it_cannot_attribute_the_covariance_to(self):
        """A band widened around the wrong node is arithmetic on an unrelated parameter."""
        from ystwin.kinetic.carotenoid import ELIZONDO2025
        from ystwin.pathway import calibrations
        from ystwin.pathway.spec import load_pathway
        from ystwin.predict import Environment, Genotype, predict_product

        spec = load_pathway("beta_carotene")
        got = predict_product(
            spec, Genotype(1.0, "b-car4"),
            Environment(growth_rate_setpoint_per_h=0.15),
            calibrations.BETA_CAROTENE_FLUX, calibrations.BETA_CAROTENE_KINETICS)

        with pytest.raises(ValueError, match="saturating node"):
            got.content_interval(spec, {}, kinetic_calibration=ELIZONDO2025)

    def test_the_default_is_unchanged_so_no_committed_number_moved_silently(self):
        """Omitting the calibration returns exactly what it always did."""
        assert self._band(0.150) == pytest.approx(1.0654, abs=5e-4)

    def test_documentation_separates_entry_only_comparison_from_current_joint_validation(self):
        from ystwin.predict import ProductPrediction

        text = ProductPrediction.content_interval.__doc__
        assert "entry-only" in text
        assert "score_product_validation" in text
        assert "product_prediction_intervals" in text
        assert "interval_status" in text
        for obsolete_column in ("share_alpha", "share_crtyb_capacity", "flux_only_too_narrow_by",
                                "with_kinetics_over_truth", "independent_over_joint"):
            assert f"column={obsolete_column}" not in text
