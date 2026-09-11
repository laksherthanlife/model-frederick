"""Separate target refits are different hypotheses, not error from unit conversion.

`scripts/score_targets.py` refits target = alpha * expression independently for each
row. Its ordering describes those competing fixed-gene models. Deriving content,
yield or titre from the SAME rate prediction and measured denominator leaves every
relative/log error unchanged. The fixed CrtE ranking is exploratory, not nested
validation of choosing the gene or the target law.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "score_targets.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("score_targets", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def table(script):
    return script.score(script._states())


class TestTheSpecificRateIsTheThingExpressionPredicts:
    def test_the_specific_rate_ranks_first(self, table):
        assert table.iloc[0].target.strip().startswith("q ")

    def test_the_separately_refitted_titre_ranks_last(self, table):
        """Fitting titre directly is not deriving it from the fitted rate."""
        assert table.iloc[-1].target.strip().startswith("titre")

    def test_every_separate_target_refit_is_worse_than_the_rate_refit(self, table):
        """Each row fits a different proportionality, rather than converting one prediction."""
        rate = table[table.target.str.strip().str.startswith("q ")].iloc[0]

        for row in table.itertuples():
            if row.target.strip().startswith("q "):
                continue
            assert row.loso_skill < rate.loso_skill, row.target

    def test_separate_refitting_changes_the_error_by_about_two_and_a_half_fold(self, table):
        best, worst = table.iloc[0], table.iloc[-1]

        ratio = worst.typical_fold_error / best.typical_fold_error
        assert 2.0 < ratio < 3.5

    def test_yield_on_glucose_beats_content(self, table):
        """The independently fitted yield law beats the independently fitted content law."""
        by_target = {r.target.strip(): r for r in table.itertuples()}
        yield_row = by_target["yield on glucose (mol/mol)"]
        content_row = by_target["content (q/mu, mmol/gDCW)"]

        assert yield_row.loso_skill > content_row.loso_skill


class TestThePermutationNullIsTheRightOne:
    """Permuting strain LABELS leaves the leave-one-out partition identical, so all six
    assignments return the observed skill and the test measures nothing. The null has to
    permute which expression PROFILE belongs to which strain's fluxes."""

    def test_every_target_sits_at_the_three_strain_floor(self, table):
        """1/6 is the best rank attainable and it is not significance."""
        assert set(table.permutation_rank) == {"1/6"}
        assert all(abs(p - 1 / 6) < 1e-9 for p in table.permutation_p)

    def test_the_null_actually_varies(self, script):
        """If the permutation did nothing, every draw would equal the observed value and the
        rank would be 6/6 rather than 1/6. This is the check that the null is a null."""
        import itertools

        frame = script._states()
        strains = sorted(frame.strain.unique())
        profile = {(r.strain, r.rate_tag): r.rel_expression for r in frame.itertuples()}
        from ystwin.pathway.flux import fit_flux_law

        skills = set()
        for order in itertools.permutations(strains):
            mapping = dict(zip(strains, order))
            shuffled = [profile[(mapping[r.strain], r.rate_tag)] for r in frame.itertuples()]
            skills.add(round(fit_flux_law(shuffled, list(frame.q_total),
                                          list(frame.strain), "CrtE", "p").loso_skill, 6))

        assert len(skills) > 1, "permuting the profiles changed nothing; the null is inert"


class TestExactConversionsPreserveTheSamePredictionError:
    def test_derived_targets_have_identical_relative_and_log_errors(self, script):
        import numpy as np

        from ystwin.pathway.flux import fit_flux_law

        frame = script._states()
        predicted = np.empty(len(frame))
        for strain in sorted(frame.strain.unique()):
            train = frame[frame.strain != strain]
            held = frame.strain == strain
            fit = fit_flux_law(list(train.rel_expression), list(train.q_total),
                               list(train.strain), "CrtE", "training-only comparison")
            predicted[held] = fit.alpha * frame.loc[held, "rel_expression"]
        rate_ratio = predicted / frame.q_total.to_numpy()
        for target in script.candidate_targets(frame).values():
            measured_conversion = target.to_numpy() / frame.q_total.to_numpy()
            derived = predicted * measured_conversion
            ratio = derived / target.to_numpy()
            np.testing.assert_allclose(ratio, rate_ratio, rtol=1e-14)
            np.testing.assert_allclose(np.log(ratio), np.log(rate_ratio), atol=1e-14)

    def test_fitting_content_directly_is_not_converting_a_rate_fit(self, script):
        import numpy as np

        from ystwin.pathway.flux import fit_flux_law

        frame = script._states()
        fit_rate = fit_flux_law(list(frame.rel_expression), list(frame.q_total),
                               list(frame.strain), "CrtE", "rate")
        fit_content = fit_flux_law(list(frame.rel_expression), list(frame.q_total / frame.mu_per_h),
                                  list(frame.strain), "CrtE", "content")
        derived = fit_rate.alpha * frame.rel_expression / frame.mu_per_h
        fitted = fit_content.alpha * frame.rel_expression
        assert not np.allclose(derived, fitted, rtol=1e-3)
        assert fit_rate.loso_rmse_log != pytest.approx(fit_content.loso_rmse_log)
