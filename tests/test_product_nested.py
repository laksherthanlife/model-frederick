"""Nested selection, exact conversions and conditional uncertainty on three strain groups."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

import ystwin.pathway.flux as flux
from ystwin.pathway.solve import solve_pathway
from ystwin.pathway.spec import load_pathway


@pytest.fixture(scope="module")
def spec():
    return load_pathway("beta_carotene")


@pytest.fixture(scope="module")
def states():
    return flux.carotenoid_measurements()


@pytest.fixture(scope="module")
def scored(spec, states):
    return flux.score_product_validation(spec, states)


@pytest.fixture(scope="module")
def intervals(spec, states):
    return flux.score_product_validation(spec, states, draws=24, seed=11)


def script(name):
    location = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    loader = importlib.util.spec_from_file_location(name, location)
    module = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(module)
    return module


def own_fold(frame, strain, columns):
    return frame.loc[frame.strain == strain, columns].reset_index(drop=True)


class TestNestedGeneAndModelSelection:
    def test_explicit_cross_product_contains_all_fourteen_genes_and_model_families(self):
        candidates = flux.product_candidates()
        assert len(candidates) == 58
        assert {c.gene for c in candidates if c.gene} == set(flux.CAROTENOID_GENES)
        assert {c.entry_law for c in candidates if c.gene} == {"rate", "content"}
        assert {c.branch for c in candidates if c.gene} == {"partition", "saturating"}
        assert {c.name for c in candidates if not c.gene} == {"constant_rate", "constant_content"}

    def test_three_independent_units_produce_six_two_channel_predictions(self, scored):
        assert len(scored) == 6
        assert scored.state.nunique() == 6
        assert scored.strain.nunique() == 3
        assert set(scored.n_independent_strains) == {3}
        assert set(scored.n_condition_predictions) == {6}
        assert set(scored.n_channels) == {2}
        assert scored.status.eq("ok").all()
        assert (scored[["product_predicted", "lycopene_rate_predicted"]] > 0).all().all()

    def test_every_outer_fold_scores_every_candidate_on_both_inner_strains(self, scored):
        selection = pd.DataFrame(scored.attrs["selection_scores"])
        assert len(selection) == 3 * 58
        assert set(selection.groupby("outer_strain").size()) == {58}
        assert set(selection.inner_strains) == {2}
        assert set(selection.inner_predictions) == {4}
        assert selection.inner_failed_predictions.eq(0).all()

    def test_selected_model_minimizes_only_its_training_inner_score(self, scored):
        selection = pd.DataFrame(scored.attrs["selection_scores"])
        for strain, group in scored.groupby("strain"):
            ranking = selection[selection.outer_strain == strain].sort_values(
                ["inner_joint_log_mse", "n_parameters", "candidate_order"])
            assert set(group.selected_model) == {ranking.iloc[0].candidate}
            assert set(group.training_strains) == {"|".join(sorted(set(scored.strain) - {strain}))}

    def test_real_selection_is_not_the_globally_ranked_crte_model(self, scored):
        by_strain = scored.groupby("strain").selected_model.first().to_dict()
        assert by_strain == {
            "b-car2": "CrtI:content:saturating",
            "b-car3": "CrtYB:content:saturating",
            "b-car4": "CrtYB:rate:saturating",
        }
        summary = flux.summarize_product_validation(scored).set_index("model")
        assert summary.loc["nested", "joint_rmse_log"] == pytest.approx(0.570295, abs=2e-6)
        assert summary.loc["fixed_crte_comparison", "joint_rmse_log"] == pytest.approx(0.291220, abs=2e-6)
        assert summary.loc["nested", "joint_rmse_log"] > summary.loc["fixed_crte_comparison", "joint_rmse_log"]

    @pytest.mark.parametrize("held_out", ["b-car2", "b-car3", "b-car4"])
    def test_poisoned_outer_outcomes_do_not_change_own_selection_prediction_or_band(
            self, spec, states, intervals, held_out):
        poisoned = states.copy()
        held = poisoned.strain == held_out
        for channel, factor in (("q_betacarotene", 70.0), ("q_lycopene", 0.03)):
            for suffix in ("", "_lo95", "_hi95"):
                poisoned.loc[held, channel + suffix] *= factor
        poisoned.loc[held, "flux"] = 999.0
        rerun = flux.score_product_validation(spec, poisoned, draws=24, seed=11)
        columns = ["selected_model", "selected_gene", "selected_entry_law", "selected_branch",
                   "inner_rmse_log", "alpha", "partition", "capacity", "km", "flux_predicted",
                   "product_predicted", "lycopene_rate_predicted", "product_rate_low",
                   "product_rate_high", "lycopene_rate_low", "lycopene_rate_high",
                   "joint_product_rate_low", "joint_lycopene_rate_high", "interval_status"]
        pd.testing.assert_frame_equal(own_fold(intervals, held_out, columns),
                                      own_fold(rerun, held_out, columns))
        before = pd.DataFrame(intervals.attrs["selection_scores"])
        after = pd.DataFrame(rerun.attrs["selection_scores"])
        pd.testing.assert_frame_equal(before[before.outer_strain == held_out].reset_index(drop=True),
                                      after[after.outer_strain == held_out].reset_index(drop=True))

    def test_heldout_inputs_can_change_predictions_but_not_selection(self, spec, states, scored):
        changed = states.copy()
        held = changed.strain == "b-car4"
        changed.loc[held, "CrtYB"] *= 1.5
        rerun = flux.score_product_validation(spec, changed)
        assert set(rerun.loc[held, "selected_model"]) == set(scored.loc[held, "selected_model"])
        assert not np.allclose(rerun.loc[held, "product_predicted"], scored.loc[held, "product_predicted"])

    def test_the_default_script_does_not_hide_a_public_prior_refusal(self, spec, states, scored, monkeypatch):
        module = script("predict_product")

        def refused(*args, **kwargs):
            raise AssertionError("validation must not use an unsupported public physiology prior")

        monkeypatch.setattr(module, "predict_product", refused)
        rerun = module.score_forward(spec, states, draws=0)
        pd.testing.assert_frame_equal(scored, rerun)

    def test_kinetic_entry_wiring_uses_the_same_nested_validator(self, states, scored):
        rerun = script("fit_carotenoid_kinetics").score_entry_models(states)
        pd.testing.assert_frame_equal(scored, rerun)

    def test_fewer_than_three_strains_cannot_support_nested_selection(self, spec, states):
        with pytest.raises(ValueError, match="at least three strain"):
            flux.score_product_validation(spec, states[states.strain != "b-car4"])


class TestDenominatorsAndTrainingBaselines:
    def test_every_comparison_keeps_all_conditions_and_both_channels(self, scored):
        comparisons = pd.DataFrame(scored.attrs["comparison_predictions"])
        assert set(comparisons.validation_mode) == {
            "nested", "constant_rate", "constant_content", "entry_partition",
            "saturating_branch", "fixed_crte_comparison",
        }
        assert set(comparisons.groupby("validation_mode").size()) == {6}
        assert set(comparisons.groupby("validation_mode").strain.nunique()) == {3}

    @pytest.mark.parametrize("name", ["constant_rate", "constant_content"])
    def test_each_baseline_uses_only_training_channel_means(self, states, scored, name):
        comparisons = pd.DataFrame(scored.attrs["comparison_predictions"])
        baseline = comparisons[comparisons.validation_mode == name]
        for row in baseline.itertuples():
            train = states[states.strain != row.strain]
            for channel, prediction in (("q_betacarotene", row.product_predicted),
                                         ("q_lycopene", row.lycopene_rate_predicted)):
                target = train[channel] / train.mu_per_h if name == "constant_content" else train[channel]
                expected = np.exp(np.log(target).mean())
                if name == "constant_content":
                    expected *= row.mu_per_h
                assert prediction == pytest.approx(expected, rel=1e-13)

    def test_missing_candidate_is_recorded_not_filtered(self, spec, states):
        result = flux.score_product_validation(spec, states, genes=(*flux.CAROTENOID_GENES, "missing"))
        selection = pd.DataFrame(result.attrs["selection_scores"])
        missing = selection[selection.gene == "missing"]
        assert len(missing) == 3 * 2 * 2
        assert missing.inner_failed_predictions.eq(4).all()
        assert np.isinf(missing.inner_joint_log_mse).all()
        assert len(result) == 6
        assert set(selection.groupby("outer_strain").size()) == {62}

    def test_ranking_keeps_a_nonpositive_gene_instead_of_silently_skipping_it(self, states):
        changed = states.copy()
        changed.loc[0, "CrtE"] = 0.0
        ranked = script("fit_pathway_flux").rank_genes(changed)
        assert len(ranked) == 14
        row = ranked[ranked.gene == "CrtE"].iloc[0]
        assert row.status == "unsupported"
        assert row.n_condition_predictions == 6
        assert np.isinf(row.loso_rmse_log)

    def test_refused_outer_rate_stays_in_the_score_denominator(self, spec, states):
        changed = states.copy()
        changed.loc[changed.condition == "4D025", "mu_per_h"] = 0.5
        result = flux.score_product_validation(spec, changed, mode="fixed_crte_comparison")
        row = result[result.state == "4D025"].iloc[0]
        assert row.status == "failed"
        assert "outside the range" in row.failure_reason
        assert len(result) == 6
        summary = flux.summarize_product_validation(result).set_index("model")
        fixed = summary.loc["fixed_crte_comparison"]
        assert fixed.n_condition_predictions == 6
        assert fixed.n_failed_predictions >= 1
        assert np.isinf(fixed.joint_rmse_log)

    def test_missing_outcome_retains_all_unscorable_cases(self, spec, states):
        result = flux.score_product_validation(spec, states.drop(columns="q_lycopene"))
        assert len(result) == 6
        assert result.status.ne("ok").all()
        assert np.isinf(result.joint_log_mse).all()
        summary = flux.summarize_product_validation(result)
        assert summary.n_condition_predictions.eq(6).all()
        assert summary.n_failed_predictions.eq(6).all()

    def test_loss_is_strain_balanced_not_pseudoreplicated(self, scored):
        changed = scored.copy()
        changed.attrs = {}
        changed["product_log_error"] = [1, 1, 2, 2, 3, 3]
        changed["lycopene_log_error"] = changed.product_log_error
        changed["joint_log_mse"] = changed.product_log_error ** 2
        extra = changed[changed.strain == "b-car2"].copy()
        extra["state"] += "_repeat"
        expanded = pd.concat([changed, extra], ignore_index=True)
        got = flux.summarize_product_validation(expanded, comparisons=False).iloc[0]
        assert got.joint_rmse_log == pytest.approx(np.sqrt((1 + 4 + 9) / 3))
        assert got.n_condition_predictions == 8
        assert got.n_independent_strains == 3


class TestRateContentAndSolverIdentity:
    def test_both_channels_and_intervals_convert_with_the_same_measured_mu(self, intervals):
        for name, content in (("product", "product_content_predicted"),
                               ("lycopene", "lycopene_predicted")):
            rate = "product_predicted" if name == "product" else "lycopene_rate_predicted"
            np.testing.assert_allclose(intervals[content] * intervals.mu_per_h, intervals[rate], rtol=1e-14)
            for side in ("low", "high"):
                np.testing.assert_allclose(intervals[f"{name}_content_{side}"] * intervals.mu_per_h,
                                           intervals[f"{name}_rate_{side}"], rtol=1e-14)
        np.testing.assert_allclose(intervals.product_predicted / intervals.product_measured,
                                   intervals.product_content_predicted / intervals.product_content_measured,
                                   rtol=1e-14)
        np.testing.assert_allclose(intervals.flux_predicted,
                                   intervals.product_predicted + intervals.lycopene_rate_predicted, rtol=1e-14)

    def test_rate_and_content_entry_models_are_distinct_laws(self, spec, states):
        train = states[states.strain != "b-car4"]
        inputs = pd.DataFrame({"mu_per_h": [0.1, 0.25], "CrtE": [0.5, 0.5]})
        for law, ratio in (("rate", 1.0), ("content", 2.5)):
            model = flux.fit_product_candidate(spec, train, flux.ProductCandidate("CrtE", law, "partition"))
            got = flux.predict_product_candidate(spec, model, inputs).sum(axis=1)
            assert got[1] / got[0] == pytest.approx(ratio)
            target = train.q_betacarotene + train.q_lycopene
            if law == "content":
                target /= train.mu_per_h
            assert model.alpha == pytest.approx(np.exp(np.log(target / train.CrtE).mean()))

    def test_saturating_predictions_reuse_the_existing_pathway_solve(self, spec, states):
        train = states[states.strain != "b-car4"]
        test = states[states.strain == "b-car4"]
        candidate = flux.ProductCandidate("CrtE", "rate", "saturating")
        model = flux.fit_product_candidate(spec, train, candidate)
        got = flux.predict_product_candidate(spec, model, test)
        for index, row in enumerate(test.itertuples()):
            solution = solve_pathway(spec, model.alpha * row.CrtE, row.mu_per_h, model.kinetics)
            assert got[index, 0] == solution.terminal.content_mmol_per_gdcw * row.mu_per_h
            assert got[index, 1] == solution.node("lycopene").content_mmol_per_gdcw * row.mu_per_h
        assert model.kinetics["lycopene"].growth_rate_range == (train.mu_per_h.min(), train.mu_per_h.max())

    @pytest.mark.parametrize("invalid", [0.0, -1.0, np.nan, np.inf])
    def test_nonfinite_or_nonpositive_mu_is_not_a_prediction(self, spec, states, invalid):
        model = flux.fit_product_candidate(spec, states, flux.ProductCandidate(None, "constant_rate", "none"))
        with pytest.raises(ValueError, match="positive and finite"):
            flux.predict_product_candidate(spec, model, pd.DataFrame({"mu_per_h": [invalid]}))


class TestConditionalTwoStageUncertainty:
    def test_nominal_targets_are_separate_from_actual_coverage(self, intervals):
        summary = flux.summarize_product_validation(intervals).set_index("model").loc["nested"]
        assert summary.nominal_channel_coverage == 0.95
        assert summary.nominal_strain_joint_coverage == 0.95
        assert summary.product_covered_conditions == int(intervals.product_inside_band.sum())
        assert summary.lycopene_covered_conditions == int(intervals.lycopene_inside_band.sum())
        assert summary.joint_covered_conditions == int(intervals.joint_inside_band.sum())
        assert summary.joint_covered_strains == int(intervals.groupby("strain").joint_inside_band.all().sum())
        assert summary.product_empirical_coverage == summary.product_covered_conditions / 6
        assert summary.joint_empirical_strain_coverage == summary.joint_covered_strains / 3
        assert intervals.interval_assumptions.str.contains("Three strains cannot establish").all()
        assert intervals.interval_assumptions.str.contains("model-selection uncertainty").all()
        assert intervals.joint_family_size.eq(4).all()

    def test_joint_strain_rectangle_is_at_least_as_wide_as_each_marginal(self, intervals):
        for channel in ("product", "lycopene"):
            assert (intervals[f"joint_{channel}_rate_low"] <= intervals[f"{channel}_rate_low"]).all()
            assert (intervals[f"joint_{channel}_rate_high"] >= intervals[f"{channel}_rate_high"]).all()

    def test_source_interval_map_recovers_both_asymmetric_endpoints(self, states):
        for column in ("CrtE", "q_betacarotene", "q_lycopene"):
            low_sd, high_sd = flux._interval_scales(states, column)
            z = norm.ppf(0.975)
            np.testing.assert_allclose(np.exp(np.log(states[column]) - z * low_sd), states[column + "_lo95"])
            np.testing.assert_allclose(np.exp(np.log(states[column]) + z * high_sd), states[column + "_hi95"])

    def test_resampling_preserves_pairs_and_refits_entry_and_kinetics_on_the_same_draw(
            self, spec, states, monkeypatch):
        train = states[states.strain != "b-car4"].copy()
        train["source_strain"] = train.strain
        test = states[states.strain == "b-car4"][["condition", "strain", "mu_per_h", "CrtE", "CrtE_lo95", "CrtE_hi95"]]
        original = flux.fit_product_candidate
        recorded = []

        def inspect(spec, sampled, candidate, **kwargs):
            assert sampled.strain.nunique() == 2
            assert set(sampled.groupby("strain").size()) == {2}
            assert sampled.groupby("strain").source_strain.nunique().eq(1).all()
            assert sampled.groupby("strain").mu_per_h.nunique().eq(2).all()
            model = original(spec, sampled, candidate, **kwargs)
            total = sampled.q_betacarotene + sampled.q_lycopene
            assert model.alpha == pytest.approx(np.exp(np.log(total / sampled.CrtE).mean()))
            recorded.append((model.alpha, model.kinetics["lycopene"].vmax_per_growth,
                             model.kinetics["lycopene"].km))
            return model

        monkeypatch.setattr(flux, "fit_product_candidate", inspect)
        bands = flux.product_prediction_intervals(
            spec, train, test, flux.ProductCandidate("CrtE", "rate", "saturating"), draws=20, seed=9)
        assert len(recorded) == 20
        assert np.all(np.std(np.asarray(recorded), axis=0) > 0)
        assert bands.interval_status.eq("ok").all()
        assert bands.interval_successful_draws.eq(20).all()

    @pytest.mark.parametrize("component", ["resample_strains", "expression_error", "observation_error"])
    def test_each_declared_uncertainty_source_actually_reaches_the_prediction(self, spec, states, component):
        train = states[states.strain != "b-car4"]
        test = states[states.strain == "b-car4"]
        candidate = flux.ProductCandidate("CrtE", "rate", "saturating")
        settings = dict(resample_strains=False, expression_error=False, observation_error=False)
        fixed = flux.product_prediction_intervals(spec, train, test, candidate, draws=16, **settings)
        settings[component] = True
        uncertain = flux.product_prediction_intervals(spec, train, test, candidate, draws=16, **settings)
        for channel in ("product", "lycopene"):
            np.testing.assert_allclose(fixed[f"{channel}_rate_high"], fixed[f"{channel}_rate_low"], rtol=1e-14)
            assert (uncertain[f"{channel}_rate_high"] > uncertain[f"{channel}_rate_low"]).all()

    def test_missing_source_interval_is_unavailable_not_a_silent_zero_uncertainty(self, spec, states):
        changed = states.drop(columns="q_betacarotene_lo95")
        result = flux.score_product_validation(spec, changed, draws=12)
        assert len(result) == 6
        assert result.interval_status.eq("failed").all()
        assert result.interval_failed_draws.eq(12).all()
        assert result.product_rate_low.isna().all()
        summary = flux.summarize_product_validation(result).set_index("model").loc["nested"]
        assert summary.n_condition_predictions == 6
        assert summary.n_intervals_available == 0
        assert summary.product_empirical_coverage == 0
        assert summary.joint_empirical_strain_coverage == 0

    def test_failed_draws_invalidate_bands_without_dropping_draws_or_conditions(self, spec, states, monkeypatch):
        original = flux.fit_product_candidate
        calls = 0

        def fail_once(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 3:
                raise RuntimeError("deliberate failed bootstrap fit")
            return original(*args, **kwargs)

        monkeypatch.setattr(flux, "fit_product_candidate", fail_once)
        result = flux.product_prediction_intervals(
            spec, states[states.strain != "b-car4"], states[states.strain == "b-car4"],
            flux.ProductCandidate("CrtE", "rate", "saturating"), draws=8)
        assert len(result) == 2
        assert result.interval_draws.eq(8).all()
        assert result.interval_successful_draws.eq(7).all()
        assert result.interval_failed_draws.eq(1).all()
        assert result.interval_status.eq("failed").all()
        assert result.product_rate_low.isna().all()

    def test_a_missing_heldout_channel_does_not_erase_the_other_channels_coverage(
            self, spec, states, intervals):
        changed = states.copy()
        changed.loc[changed.strain == "b-car2", "q_lycopene"] = np.nan
        result = flux.score_product_validation(spec, changed, draws=24, seed=11)
        before = own_fold(intervals, "b-car2", ["product_inside_band", "product_predicted"])
        after = own_fold(result, "b-car2", ["product_inside_band", "product_predicted"])
        pd.testing.assert_frame_equal(before, after)
        own = result[result.strain == "b-car2"]
        assert own.status.eq("unscorable").all()
        assert not own.lycopene_inside_band.any()
        assert not own.joint_inside_band.any()
        assert np.isfinite(own.product_log_error).all()
        assert np.isinf(own.lycopene_log_error).all()

    @pytest.mark.parametrize("settings", [{"draws": -1}, {"draws": 1.5}, {"draws": True},
                                         {"nominal_coverage": np.nan}, {"nominal_coverage": 1.0}])
    def test_invalid_interval_configuration_is_rejected_even_if_all_models_fail(self, spec, states, settings):
        with pytest.raises(ValueError, match="draws|nominal coverage"):
            flux.score_product_validation(spec, states.drop(columns="q_lycopene"), **settings)

    def test_zero_draws_is_explicitly_not_a_coverage_estimate(self, scored):
        assert scored.interval_status.eq("not_requested").all()
        summary = flux.summarize_product_validation(scored)
        assert summary.product_empirical_coverage.isna().all()
        assert summary.joint_empirical_strain_coverage.isna().all()
