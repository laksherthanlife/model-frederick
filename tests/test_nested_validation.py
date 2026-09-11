from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_latent_hidden_entries_do_not_set_training_centres(monkeypatch):
    from ystwin.analysis import latent

    data = np.random.default_rng(3).normal(size=(30, 4))
    data[0] += 100.0
    centred = data - data.mean(axis=0)
    original = latent._factorise
    means = []

    def inspect(values, observed, n_states):
        means.append(np.array([values[observed[:, c], c].mean()
                               for c in range(values.shape[1])]))
        return original(values, observed, n_states)

    monkeypatch.setattr(latent, "_factorise", inspect)
    latent._heldout_error(centred, np.isfinite(centred), 1, folds=3, seed=7)
    assert len(means) == 3
    assert np.max(np.abs(means)) < 1e-10


def test_stressor_transfer_excludes_codosings_of_the_heldout_agent():
    from ystwin.analysis.transfer import transfer_test
    from ystwin.generator.panel_experiment import panel_dataset

    data = panel_dataset(reporters=["STRE-general", "UPRE-ER", "TRX2-oxidative"],
                         stressors=["DTT", "H2O2", "heat"], doses=(0.5, 1.0),
                         replicates=3, combinations=[("DTT", "H2O2")], seed=4)
    result = transfer_test(data, "DTT", 1, observed=(0, 1))
    assert "DTT" not in result.trained_on
    assert "DTT+H2O2" not in result.trained_on
    assert set(result.trained_on) == {"H2O2", "heat"}


def test_a_fold_over_its_own_control_is_identically_one():
    from ystwin.analysis.uncertainty import fold_change

    frame = pd.DataFrame({
        "plate": np.repeat(["a", "b", "c", "d"], 3),
        "well": list("ABC") * 4,
        "construct": ["R"] * 12,
        "dose_mM": np.zeros(12),
        "activity_late": np.arange(1.0, 13.0),
    })
    result = fold_change(frame, "R", 0.0, n_resamples=200)
    assert (result.point, result.low, result.high) == (1.0, 1.0, 1.0)
    assert result.n_wells == 12


def test_missing_p_values_do_not_become_significant_or_shrink_the_family():
    from ystwin.analysis.multiplicity import benjamini_hochberg, holm

    for correction in (holm, benjamini_hochberg):
        adjusted = correction([0.01, np.nan])
        assert adjusted[0] == pytest.approx(0.02)
        assert np.isnan(adjusted[1])


def test_power_batch_noise_preserves_the_true_induction(monkeypatch):
    from ystwin.analysis import power
    from ystwin.generator.design import stressor_only_series
    from ystwin.generator.literature import parameters_from_literature

    params = parameters_from_literature("UPRE2")
    conditions = stressor_only_series(params, doses=(0.0, 0.1, 1.0))
    ratios = []
    original = power._observe

    def observe(condition, p, *args):
        ratios.append(p.promoter_peak / p.promoter_basal)
        return original(condition, p, *args)

    monkeypatch.setattr(power, "_observe", observe)
    power._one_trial(params, conditions, 1.0, 4, np.random.default_rng(1),
                     0.01, 0.4, True, well_cv=0.2)
    assert len(ratios) == 12
    assert ratios == pytest.approx(np.ones(12))


def test_branch_mass_balance_survives_small_positive_entry_flux():
    from ystwin.kinetic.carotenoid import solve_branch_from_flux

    flux = 1e-30
    state = solve_branch_from_flux(flux, 0.15)
    assert state.lycopene_content > 0.0
    assert state.beta_carotene_rate + 0.15 * state.lycopene_content == pytest.approx(
        flux, rel=1e-12, abs=0.0)


@pytest.fixture(scope="module")
def product_inputs():
    from ystwin.pathway.flux import carotenoid_measurements
    from ystwin.pathway.spec import load_pathway

    return load_pathway("beta_carotene"), carotenoid_measurements()


@pytest.fixture(scope="module")
def product_validation(product_inputs):
    from ystwin.pathway.flux import score_product_validation

    return score_product_validation(*product_inputs, draws=24, seed=9)


def test_joint_score_cannot_hide_a_bad_product_behind_the_total_flux(product_validation):
    from ystwin.pathway.flux import _log_errors

    truth = np.array([[1.0, 9.0], [9.0, 1.0]])
    predicted = truth[:, ::-1]
    assert np.allclose(predicted.sum(axis=1), truth.sum(axis=1))
    assert np.mean(_log_errors(predicted, truth) ** 2) == pytest.approx(np.log(9) ** 2)
    frame = product_validation
    assert frame.joint_log_mse.to_numpy() == pytest.approx(
        (frame.product_log_error.to_numpy() ** 2 + frame.lycopene_log_error.to_numpy() ** 2) / 2)


def test_product_evidence_counts_three_strains_six_conditions_and_two_paired_peaks(product_validation):
    from ystwin.pathway.flux import summarize_product_validation

    scores = summarize_product_validation(product_validation)
    assert (scores.n_independent_strains == 3).all()
    assert (scores.n_condition_predictions == 6).all()
    assert (scores.n_channels == 2).all()
    assert {"constant_rate", "constant_content", "entry_partition", "saturating_branch", "nested"} <= set(scores.model)
    assert product_validation.state.nunique() == len(product_validation) == 6


def test_product_selection_refits_all_candidates_inside_outer_training_only(product_validation):
    scores = pd.DataFrame(product_validation.attrs["selection_scores"])
    assert len(scores) == 3 * 58
    assert (scores.inner_strains == 2).all()
    assert (scores.inner_predictions == 4).all()
    for strain, frame in product_validation.groupby("strain"):
        assert set(frame.training_strains) == {"|".join(sorted(set(product_validation.strain) - {strain}))}


def test_outer_product_labels_cannot_change_choices_parameters_predictions_or_intervals(
        product_inputs, product_validation):
    from ystwin.pathway.flux import score_product_validation

    spec, states = product_inputs
    changed = states.copy()
    mask = changed.strain == "b-car2"
    for column in ("q_betacarotene", "q_lycopene"):
        for suffix in ("", "_lo95", "_hi95"):
            changed.loc[mask, column + suffix] *= 2
    result = score_product_validation(spec, changed, draws=24, seed=9)
    columns = ["state", "selected_gene", "selected_model", "selected_entry_law", "selected_branch",
               "alpha", "capacity", "km", "product_predicted", "lycopene_rate_predicted",
               "product_rate_low", "product_rate_high", "lycopene_rate_low", "lycopene_rate_high"]
    a = product_validation.query("strain == 'b-car2'")[columns].reset_index(drop=True)
    b = result.query("strain == 'b-car2'")[columns].reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)


def test_matched_entry_ablation_changes_only_the_saturation_branch(product_inputs, product_validation):
    from ystwin.analysis.validation import matched_product_branches

    spec, states = product_inputs
    predictions = matched_product_branches(spec, states, product_validation)
    full = predictions.query("branch == 'saturating'").set_index("condition")
    ablated = predictions.query("branch == 'partition'").set_index("condition")
    assert len(full) == len(ablated) == 6
    assert (full.selected_gene == ablated.selected_gene).all()
    assert (full.selected_entry_law == ablated.selected_entry_law).all()
    assert full.alpha.to_numpy() == pytest.approx(ablated.alpha.to_numpy())
    assert full.flux_predicted.to_numpy() == pytest.approx(ablated.flux_predicted.to_numpy())
    assert not np.allclose(full.product_predicted, ablated.product_predicted)
    for frame in (full, ablated):
        assert (frame.product_predicted + frame.lycopene_rate_predicted).to_numpy() == pytest.approx(frame.flux_predicted.to_numpy())


def test_reported_coverage_is_actual_paired_coverage_not_interval_width(product_validation):
    from ystwin.pathway.flux import summarize_product_validation

    frame = product_validation
    actual = ((frame.product_measured >= frame.joint_product_rate_low)
              & (frame.product_measured <= frame.joint_product_rate_high)
              & (frame.lycopene_rate_measured >= frame.joint_lycopene_rate_low)
              & (frame.lycopene_rate_measured <= frame.joint_lycopene_rate_high))
    coverage = summarize_product_validation(frame).query("model == 'nested'").iloc[0]
    assert coverage.n_intervals_available == 6
    assert coverage.n_independent_strains == 3
    assert coverage.joint_empirical_coverage == pytest.approx(actual.mean())
    assert coverage.nominal_strain_joint_coverage == 0.95
    assert frame.interval_assumptions.str.contains("cannot establish").all()


def test_nested_module_selection_is_blind_to_outer_targets():
    from dataclasses import replace
    from ystwin.analysis.validation import validate_modules
    from ystwin.generator.panel_experiment import panel_dataset

    data = panel_dataset(reporters=["STRE-general", "UPRE-ER", "TRX2-oxidative"],
                         stressors=["DTT", "H2O2", "heat", "NaCl"],
                         doses=(0.1, 0.5, 1.0, 2.0), replicates=3, seed=7)
    first = validate_modules(data, "DTT", max_states=2)
    targets = data.modules.copy()
    targets[data.mask("DTT")] += 1000.0
    second = validate_modules(replace(data, modules=targets), "DTT", max_states=2)
    assert first.n_states == second.n_states
    assert first.inner_scores == second.inner_scores
    assert first.prediction == pytest.approx(second.prediction)
    assert first.baseline == pytest.approx(second.baseline)


def test_group_firewall_keeps_all_replicates_together():
    from ystwin.analysis.validation import group_folds

    groups = np.array(["a", "a", "b", "b", "c", "c"])
    folds = group_folds(groups)
    assert len(folds) == 3
    for fold in folds:
        assert set(groups[fold.train]).isdisjoint(groups[fold.test])
        assert len(fold.test) == 2


def test_fold_indices_are_immutable_snapshots():
    from ystwin.analysis.validation import ValidationFold

    train = np.array([0, 1])
    fold = ValidationFold("c", train, np.array([2]), np.array([], dtype=int))
    train[0] = 2
    assert fold.train.tolist() == [0, 1]
    with pytest.raises(ValueError):
        fold.train[0] = 2
    with pytest.raises(ValueError, match="disjoint"):
        ValidationFold("c", np.array([0, 2]), np.array([2]), np.array([], dtype=int))


def test_module_validation_reports_worse_than_baseline_instead_of_clipping_it():
    from ystwin.analysis.validation import training_mean_skill

    score = training_mean_skill([[1.0], [2.0]], [[5.0], [6.0]], [[0.0], [1.0]])
    assert score[0] == pytest.approx(-11.8)


def test_matched_branch_refusals_keep_the_full_condition_denominator(product_inputs, product_validation):
    from ystwin.analysis.validation import matched_product_branches

    spec, states = product_inputs
    changed = product_validation.copy()
    changed["selected_gene"] = None
    changed["selected_model"] = "constant_rate"
    changed["selected_entry_law"] = "constant_rate"
    rows = matched_product_branches(spec, states, changed)
    assert len(rows) == 12
    assert rows.status.eq("failed").all()
    assert np.isinf(rows.joint_log_mse).all()


def _assert_matched_condition_denominators(rows, states):
    expected = {(branch, row.condition, row.strain)
                for branch in ("partition", "saturating") for row in states.itertuples()}
    assert len(rows) == len(expected) == 12
    assert set(rows[["branch", "condition", "strain"]].itertuples(index=False, name=None)) == expected
    assert rows.groupby("branch").size().eq(6).all()
    assert rows.groupby("branch").strain.nunique().eq(3).all()
    assert rows.groupby(["branch", "strain"]).size().eq(2).all()


@pytest.mark.parametrize("failure_stage, failed_conditions", [
    pytest.param("fit", ("4D01", "4D025"), id="shared-fit"),
    pytest.param("prediction", ("4D01",), id="first-sibling"),
    pytest.param("prediction", ("4D025",), id="last-sibling"),
])
@pytest.mark.parametrize("branch", ["partition", "saturating"])
@pytest.mark.parametrize("error_type", [ValueError, RuntimeError, FloatingPointError, OverflowError])
def test_matched_branch_failure_scope_retains_other_predictions(
        product_inputs, product_validation, monkeypatch, failure_stage, failed_conditions,
        branch, error_type):
    from ystwin.analysis.validation import matched_product_branches
    from ystwin.pathway import flux

    spec, states = product_inputs
    reference = matched_product_branches(spec, states, product_validation)
    original_fit, original_predict = flux.fit_product_candidate, flux.predict_product_candidate
    fitting_calls, prediction_calls = [], []
    reason = f"deliberate {failure_stage} refusal"

    def fit(spec, train, candidate):
        fitting_calls.append((tuple(train.condition), candidate.branch))
        if (failure_stage == "fit" and candidate.branch == branch
                and "b-car4" not in set(train.strain)):
            raise error_type(reason)
        return original_fit(spec, train, candidate)

    def predict(spec, model, inputs):
        prediction_calls.append((model.candidate.branch, tuple(inputs.condition)))
        if (failure_stage == "prediction" and model.candidate.branch == branch
                and inputs.condition.isin(failed_conditions).any()):
            raise error_type(reason)
        return original_predict(spec, model, inputs)

    monkeypatch.setattr(flux, "fit_product_candidate", fit)
    monkeypatch.setattr(flux, "predict_product_candidate", predict)
    rows = matched_product_branches(spec, states, product_validation)
    _assert_matched_condition_denominators(rows, states)
    failed = rows.branch.eq(branch) & rows.condition.isin(failed_conditions)
    assert rows.loc[failed, "status"].eq("failed").all()
    assert rows.loc[failed, "failure_reason"].eq(reason).all()
    assert rows.loc[failed, ["alpha", "flux_predicted", "product_predicted",
                             "lycopene_rate_predicted"]].isna().all().all()
    assert np.isinf(rows.loc[failed, "joint_log_mse"]).all()
    assert np.isinf(rows.loc[rows.branch == branch, "joint_log_mse"].to_numpy().mean())
    pd.testing.assert_frame_equal(rows.loc[~failed], reference.loc[~failed])
    expected_fits = {(tuple(states.loc[states.strain != strain, "condition"]), law)
                     for strain in states.strain.unique() for law in ("partition", "saturating")}
    assert len(fitting_calls) == len(expected_fits) == 6
    assert set(fitting_calls) == expected_fits
    if failure_stage == "prediction":
        assert all(len(conditions) == 1 for _, conditions in prediction_calls)
    expected_predictions = rows.loc[~failed] if failure_stage == "fit" else rows
    assert sorted((law, condition) for law, conditions in prediction_calls for condition in conditions) == sorted(
        expected_predictions[["branch", "condition"]].itertuples(index=False, name=None))


@pytest.mark.parametrize("failed_condition", ["4D01", "4D025"])
def test_matched_branch_real_growth_refusal_keeps_its_valid_sibling(
        product_inputs, product_validation, failed_condition):
    from ystwin.analysis.validation import matched_product_branches

    spec, states = product_inputs
    reference = matched_product_branches(spec, states, product_validation)
    changed = states.copy()
    changed.loc[changed.condition == failed_condition, "mu_per_h"] = 0.5
    rows = matched_product_branches(spec, changed, product_validation)
    _assert_matched_condition_denominators(rows, states)
    failed = rows.branch.eq("saturating") & rows.condition.eq(failed_condition)
    assert rows.loc[failed, "status"].eq("failed").all()
    assert rows.loc[failed, "failure_reason"].str.contains("outside the range").all()
    assert np.isinf(rows.loc[failed, "joint_log_mse"]).all()
    sibling = rows.branch.eq("saturating") & rows.strain.eq("b-car4") & ~failed
    pd.testing.assert_frame_equal(rows.loc[sibling], reference.loc[sibling])
    assert rows.loc[~failed, "status"].eq("ok").all()


@pytest.mark.parametrize("stage", ["fit", "predict"])
@pytest.mark.parametrize("error_type", [KeyError, TypeError])
def test_matched_branch_unexpected_errors_still_propagate(
        product_inputs, product_validation, monkeypatch, stage, error_type):
    from ystwin.analysis.validation import matched_product_branches
    from ystwin.pathway import flux

    failure = error_type("unexpected implementation error")

    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(flux, f"{stage}_product_candidate", fail)
    with pytest.raises(error_type) as exc:
        matched_product_branches(*product_inputs, product_validation)
    assert exc.value is failure


@pytest.mark.parametrize("held_out", ["b-car2", "b-car3", "b-car4"])
def test_matched_branch_fits_and_predictions_are_blind_to_outer_outcomes(
        product_inputs, product_validation, monkeypatch, held_out):
    from ystwin.analysis.validation import matched_product_branches
    from ystwin.pathway import flux

    spec, states = product_inputs
    original_fit, original_predict = flux.fit_product_candidate, flux.predict_product_candidate
    fits = []

    def fit(spec, train, candidate):
        held = set(states.strain) - set(train.strain)
        assert len(held) == 1
        assert len(train) == 4 and train.strain.nunique() == 2
        selected = product_validation.loc[product_validation.strain.isin(held)].iloc[0]
        assert (candidate.gene, candidate.entry_law) == (selected.selected_gene, selected.selected_entry_law)
        model = original_fit(spec, train, candidate)
        fits.append((held.pop(), train.copy(), model))
        return model

    def predict(spec, model, inputs):
        held, train, _ = next(record for record in fits if record[2] is model)
        assert set(inputs.strain) == {held}
        assert set(inputs.condition).isdisjoint(train.condition)
        assert list(inputs.columns) == ["condition", "strain", "mu_per_h", model.candidate.gene]
        return original_predict(spec, model, inputs)

    monkeypatch.setattr(flux, "fit_product_candidate", fit)
    monkeypatch.setattr(flux, "predict_product_candidate", predict)
    reference = matched_product_branches(spec, states, product_validation)
    changed = states.copy()
    for channel, factor in (("q_betacarotene", 70.0), ("q_lycopene", 0.03)):
        for suffix in ("", "_lo95", "_hi95"):
            changed.loc[changed.strain == held_out, channel + suffix] *= factor
    changed.loc[changed.strain == held_out, "flux"] = 999.0
    rows = matched_product_branches(spec, changed, product_validation)
    _assert_matched_condition_denominators(rows, states)
    assert len(fits) == 12
    expected_fits = sorted((strain, branch) for strain in states.strain.unique()
                           for branch in ("partition", "saturating"))
    for run in (fits[:6], fits[6:]):
        assert sorted((held, model.candidate.branch) for held, _, model in run) == expected_fits
    own_fits = [(train, model) for held, train, model in fits if held == held_out]
    for (before_train, before_model), (after_train, after_model) in zip(own_fits[:2], own_fits[2:]):
        pd.testing.assert_frame_equal(before_train, after_train)
        assert before_model == after_model
    own = rows.strain.eq(held_out)
    columns = rows.columns.difference(["product_measured", "lycopene_rate_measured", "joint_log_mse"])
    pd.testing.assert_frame_equal(rows.loc[own, columns], reference.loc[own, columns])
    assert not np.allclose(rows.loc[own, "joint_log_mse"], reference.loc[own, "joint_log_mse"])


def test_heldout_combination_is_not_retained_inside_a_larger_codose():
    from ystwin.analysis.validation import stressor_fold
    from ystwin.generator.panel_experiment import PanelDataset

    labels = np.array(["DTT", "H2O2", "heat", "DTT+H2O2", "DTT+H2O2+heat"])
    data = PanelDataset(np.ones((5, 2)), labels, np.zeros(5), np.zeros((5, 24)), ["a", "b"])
    fold = stressor_fold(data, "DTT+H2O2")
    assert data.labels[fold.excluded].tolist() == ["DTT+H2O2+heat"]
    assert data.labels[fold.train].tolist() == ["DTT", "H2O2", "heat"]
