"""The cross-family script: its verdict, argmax, memo key and reporting schema.

``scripts/run_cross_family.py`` writes ``cross_family_transfer.csv`` and
``cross_family_optimism.csv`` and prints the verdict ``docs/CROSS_FAMILY.md`` quotes.
Its decisions and configuration/target metadata are testable without rerunning the scoring:

* ``verdict`` picks between three outcomes the script states *before* it sees the numbers.
  Which one it prints is the result.
* ``best_solution`` is an exhaustive argmax whose only tie-break is the incoming candidate.
  A tie resolved the other way silently changes the reported configuration.
* ``family_score`` memoises on a key that has to be a complete identity. A key missing one
  argument would return a ``--scale`` number for a default-budget question.

The scoring itself runs the generator and is exercised by ``tests/test_families.py`` and
``tests/test_optimism.py``; nothing here re-runs it. The module is imported with
``YSTWIN_OUTPUTS`` redirected, because ``OUT`` binds at import time.
"""

from __future__ import annotations

import importlib.util
import json
import os
import pathlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis.nulls import NullResult, skill_score_if_defined
from ystwin.analysis.optimism import OptimismBound

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_cross_family.py"
_TRANSFER_SCRIPT = _REPO / "scripts" / "run_transfer.py"


@pytest.fixture(scope="module")
def script(tmp_path_factory):
    """The script as a module, with its output directory redirected away from outputs/."""
    redirect = tmp_path_factory.mktemp("cross_family_outputs")
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(redirect)
    try:
        spec = importlib.util.spec_from_file_location("run_cross_family", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


def _families(**scores) -> pd.DataFrame:
    """A families table with only the columns ``verdict`` reads."""
    return pd.DataFrame([
        {"family": name, "transfer_r2": r2, "beats_null": beats}
        for name, (r2, beats) in scores.items()])


class TestTheVerdictReportsTheOutcomeThatObtained:
    """The three outcomes are stated in the script's own prose before any number is
    computed, precisely so the reading cannot be chosen afterwards. ``verdict`` is where
    that promise is kept or broken."""

    def test_no_family_beating_its_null_is_the_third_outcome(self, script):
        """Non-rejection detects no superiority; it proves neither equality nor absence."""
        families = _families(baseline=(0.03, False), edge_dropped=(0.02, False),
                             biphasic=(0.01, False))

        lines = script.verdict(families, 0.03, 0.02, 0.01, 0.03)
        text = " ".join(lines).lower()

        assert lines[0].startswith("THIRD OUTCOME")
        assert "no detected superiority" in text
        assert "does not establish equivalence or absence" in text
        assert "as much transfer information" not in text
        assert "there is no effect" not in text

    def test_the_baseline_alone_beating_its_null_is_the_failure_outcome(self, script):
        """Detected superiority on the baseline alone does not validate biological transfer."""
        families = _families(baseline=(0.20, True), edge_dropped=(0.02, False),
                             biphasic=(0.01, False))

        lines = script.verdict(families, 0.20, 0.02, 0.05, 0.20)

        assert lines[0].startswith("FAILURE")

    def test_a_variant_beating_its_null_too_is_the_success_outcome(self, script):
        families = _families(baseline=(0.20, True), edge_dropped=(0.18, True),
                             biphasic=(0.01, False))

        lines = script.verdict(families, 0.19, 0.10, 0.05, 0.20)

        assert lines[0].startswith("SUCCESS")

    def test_the_worst_variant_is_named_and_the_baseline_is_not_eligible_for_it(self,
                                                                               script):
        """The spread is reported as baseline against worst *variant*. If the baseline
        could be its own worst variant the sentence would compare a number with itself."""
        families = _families(baseline=(0.01, False), edge_dropped=(0.30, False),
                             biphasic=(0.20, False))

        text = " ".join(script.verdict(families, 0.2, 0.2, 0.01, 0.01))

        assert "(biphasic)" in text

    def test_a_large_upper_bound_is_not_evidence_of_large_loss(self, script):
        """An upper confidence bound is not a lower bound on any future loss."""
        families = _families(baseline=(0.03, False), edge_dropped=(0.02, False))

        text = " ".join(script.verdict(families, 0.03, 0.02, ucbog=0.10, baseline_r2=0.03))

        assert "upper bound alone" in text
        assert "does not establish that the true gap is large" in text
        assert "cannot be defended on this evidence" not in text
        assert "cost exceeds" not in text

    def test_a_bound_exactly_equal_to_the_score_does_not_trigger_the_large_bound_note(self,
                                                                                   script):
        families = _families(baseline=(0.03, False), edge_dropped=(0.02, False))

        text = " ".join(script.verdict(families, 0.03, 0.02, ucbog=0.03, baseline_r2=0.03))

        assert "upper bound alone" not in text

    def test_split_difference_is_not_selected_pipeline_optimism(self, script):
        families = _families(baseline=(0.03, False), edge_dropped=(0.02, False))

        text = " ".join(script.verdict(families, 0.03, 0.02, ucbog=0.01, baseline_r2=0.03))

        assert "Descriptive split difference" in text
        assert "Selected-pipeline optimism: pending" in text
        assert "fixed-candidate" in text
        assert "not each future loss" in text
        assert "Selection optimism, chosen-on minus held-out" not in text

    def test_a_zero_baseline_does_not_divide_by_it(self, script):
        """A zero baseline has no finite UCBOG-to-score ratio."""
        families = _families(baseline=(0.0, False), edge_dropped=(0.0, False))

        text = " ".join(script.verdict(families, 0.0, 0.0, ucbog=0.01, baseline_r2=0.0))

        assert "For scale only" not in text


class TestTheSearchIsExhaustiveAndItsOnlyTieBreakIsTheCandidate:
    """``best_solution`` chooses the latent dimension and the co-dosing design that get
    written into ``cross_family_optimism.csv``. SPOTA passes the candidate in so a reference
    cannot report a negative gap; here that argument does nothing except settle exact
    ties, and a change that made it do more would break the bound's meaning."""

    @staticmethod
    def _flat(script, monkeypatch, value=0.5):
        """Every configuration scores the same, so every configuration is a winner."""
        monkeypatch.setattr(script, "family_score",
                            lambda family, design, seed, n_states, doses, replicates: value)

    def test_every_dimension_and_design_is_scored(self, script, monkeypatch):
        seen = []

        def record(family, design, seed, n_states, doses, replicates):
            seen.append((n_states, design))
            return 0.0

        monkeypatch.setattr(script, "family_score", record)
        script.best_solution(["one-family"], 0, (1.0,), 1, initial=None)

        assert set(seen) == {(k, d) for k in script.K_GRID for d in script.DESIGNS}

    def test_the_highest_scoring_configuration_wins(self, script, monkeypatch):
        monkeypatch.setattr(
            script, "family_score",
            lambda family, design, seed, n_states, doses, replicates:
                1.0 if (n_states, design) == (4, "recommended pairs") else 0.1)

        assert script.best_solution([object()], 0, (1.0,), 1, initial=None) == (
            4, "recommended pairs")

    def test_an_exact_tie_goes_to_the_incoming_candidate(self, script, monkeypatch):
        """A reference that starts from the candidate and ties with it must return the
        candidate, or the gap it reports is a difference between two equal scores dressed
        as a different configuration."""
        self._flat(script, monkeypatch)
        candidate = (4, "recommended pairs")

        assert script.best_solution([object()], 0, (1.0,), 1, initial=candidate) == candidate

    def test_an_exact_tie_with_no_candidate_is_settled_deterministically(self, script,
                                                                        monkeypatch):
        """With nothing to break the tie the first configuration in the grid wins, every
        time. A dict iteration order that varied would make the reported candidate a
        property of the run."""
        self._flat(script, monkeypatch)
        first = (script.K_GRID[0], next(iter(script.DESIGNS)))

        assert script.best_solution([object()], 0, (1.0,), 1, initial=None) == first

    def test_a_candidate_that_did_not_win_is_not_returned(self, script, monkeypatch):
        """``initial`` is a tie-break and nothing else."""
        monkeypatch.setattr(
            script, "family_score",
            lambda family, design, seed, n_states, doses, replicates:
                1.0 if n_states == 3 else 0.1)

        assert script.best_solution([object()], 0, (1.0,), 1,
                                    initial=(2, "blocked"))[0] == 3


class TestTheObjectiveIsAnAverageAndNotAWorstCase:
    """Peng et al.'s risk-neutral expectation over the draw. A minimum would report the
    hardest structure and would claim a robustness guarantee domain randomization does not
    give, so the difference is a claim and not a detail."""

    def test_the_score_is_the_mean_over_the_families(self, script, monkeypatch):
        scores = iter([0.1, 0.9])
        monkeypatch.setattr(
            script, "family_score",
            lambda family, design, seed, n_states, doses, replicates: next(scores))

        got = script.objective((3, "blocked"), [object(), object()], 0, (1.0,), 1)

        assert got == pytest.approx(0.5)


class TestTheMemoKeyIsACompleteIdentity:
    """A quarter of SPOTA's sweeps are repeats, which is why the cache exists. Every
    argument that reaches the generator has to be in the key: one left out and the cache
    answers a different question with a stored number, silently and only under ``--scale``
    or a second seed."""

    @pytest.fixture
    def counted(self, script, monkeypatch):
        calls = []

        monkeypatch.setattr(script, "_SCORES", {})
        monkeypatch.setattr(script, "dataset",
                            lambda family, design, seed, doses, replicates: (
                                family, design, seed, doses, replicates))
        monkeypatch.setattr(script, "median_transfer",
                            lambda data, n_states: calls.append((data, n_states)) or 0.0)
        return calls

    @staticmethod
    def _family(name="baseline", draw=0):
        from ystwin.generator.families import build_family

        return build_family(name, seed=draw)

    def test_the_same_question_is_computed_once(self, script, counted):
        family = self._family()
        for _ in range(3):
            script.family_score(family, "blocked", 0, 3, (0.5, 1.0), 3)

        assert len(counted) == 1

    @pytest.mark.parametrize("changed", [
        {"design": "recommended pairs"},
        {"seed": 1},
        {"n_states": 4},
        {"doses": (0.25, 0.5, 1.0)},
        {"replicates": 6},
    ])
    def test_changing_any_argument_recomputes(self, script, counted, changed):
        family = self._family()
        base = dict(design="blocked", seed=0, n_states=3, doses=(0.5, 1.0), replicates=3)
        script.family_score(family, **base)
        script.family_score(family, **{**base, **changed})

        assert len(counted) == 2

    def test_two_draws_of_the_same_family_are_different_questions(self, script, counted):
        """``build_family`` is deterministic in (name, draw), so the draw is what separates
        two samples of one perturbation. Keying on the name alone would collapse the whole
        reference set onto its first member."""
        base = dict(design="blocked", seed=0, n_states=3, doses=(0.5, 1.0), replicates=3)
        script.family_score(self._family(draw=0), **base)
        script.family_score(self._family(draw=1), **base)

        assert len(counted) == 2


class TestTheComparisonKeepsItsOrigin:
    """The scripts share panel, reader and revealed-channel definitions, not their
    selection procedure. Matching these inputs must not imply matching outer scores."""

    @pytest.fixture(scope="class")
    def transfer(self, tmp_path_factory):
        redirect = tmp_path_factory.mktemp("transfer_outputs")
        previous = os.environ.get("YSTWIN_OUTPUTS")
        os.environ["YSTWIN_OUTPUTS"] = str(redirect)
        try:
            spec = importlib.util.spec_from_file_location("run_transfer", _TRANSFER_SCRIPT)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        finally:
            if previous is None:
                os.environ.pop("YSTWIN_OUTPUTS", None)
            else:
                os.environ["YSTWIN_OUTPUTS"] = previous
        return module

    def test_the_panel_is_the_one_run_transfer_scores(self, script, transfer):
        assert script.PANEL == transfer.PANEL

    def test_the_readers_are_the_ones_run_transfer_reads(self, script, transfer):
        assert script.READERS == transfer.READERS

    def test_the_revealed_channels_are_the_ones_run_transfer_reveals(self, script,
                                                                     transfer):
        """``run_transfer.py`` passes ``(0, 1, 2, 3)`` to every ``leave_one_stressor_out``
        call. Revealing a different set here would change what "the same test" means."""
        assert script.OBSERVED == (0, 1, 2, 3)


class TestTheHeldOutFamiliesAreReallyHeldOut:
    """The M-open split of CIRCULARITY.md §8. It is written as one tuple and one
    comprehension over the registry, so adding a family to the registry silently adds it to
    the held-out side -- which is correct -- but adding one to ``TRAINING`` without it
    existing would silently shrink the evaluation."""

    def test_the_two_sides_do_not_overlap(self, script):
        assert not set(script.TRAINING) & set(script.HELD_OUT)

    def test_between_them_they_cover_the_registry(self, script):
        from ystwin.generator.families import family_names

        assert set(script.TRAINING) | set(script.HELD_OUT) == set(family_names())

    def test_every_training_family_exists(self, script):
        from ystwin.generator.families import family_names

        assert set(script.TRAINING) <= set(family_names())

    def test_something_is_actually_held_out(self, script):
        assert script.HELD_OUT

    def test_the_baseline_is_on_the_training_side(self, script):
        """Every delta in the table is measured against the baseline row, and the split's
        two means are reported side by side. A baseline on the held-out side would make
        the held-out mean partly a measurement of the reference point."""
        assert "baseline" in script.TRAINING


@pytest.fixture
def report(script, monkeypatch, tmp_path, capsys):
    scores = dict(zip(script.family_names(), (0.04, -0.05, 0.03, 0.02, -0.06, 0.025, 0.01)))
    calls = SimpleNamespace(datasets=[], scores=[], searches=[], optimism=[])
    monkeypatch.setattr(script, "OUT", tmp_path)

    def dataset(family, design, seed, doses, replicates):
        data = SimpleNamespace(family=family, design=design, seed=seed,
                               doses=doses, replicates=replicates)
        calls.datasets.append(data)
        return data

    def transfer(data, n_states, observed, only):
        calls.scores.append((data.family.name, n_states, observed, only))
        return pd.DataFrame({"r2": [scores[data.family.name]], "oracle_r2": [0.2],
                             "alignment": [0.5]})

    def null(statistic, data, surrogate, n_draws, greater_is_better, label, seed):
        value = statistic(data)
        draws = np.full(n_draws, np.nan if label == "loading_noise" else value + 0.1)
        return NullResult(value, draws, greater_is_better, label)

    def objective(solution, families, seed, doses, replicates):
        return float(np.mean([scores[family.name] for family in families]))

    monkeypatch.setattr(script, "dataset", dataset)
    monkeypatch.setattr(script, "leave_one_stressor_out", transfer)
    monkeypatch.setattr(script, "compare_to_null", null)
    monkeypatch.setattr(script, "objective", objective)

    def run(chosen=(3, "recommended pairs"), candidate=(2, "recommended pairs"), scale=False):
        def best_solution(families, seed, doses, replicates, initial):
            calls.searches.append(([family.name for family in families], seed,
                                   doses, replicates, initial))
            return chosen

        bounds = []

        def estimate(problem, **kwargs):
            calls.optimism.append((problem, kwargs))
            raw = np.zeros(kwargs["n_reference"])
            raw[:3] = (-0.00123456789, 0.004, 0.008)
            gaps = np.clip(raw, 0.0, None)
            bound = OptimismBound(
                bound=0.01, mean_gap=float(gaps.mean()), gaps=gaps, raw_gaps=raw,
                candidate=candidate, alpha=kwargs["alpha"], n_bootstrap=kwargs["n_bootstrap"],
                n_candidate_domains=kwargs["n_candidate_domains"],
                n_reference_domains=kwargs["n_reference_domains"],
                target="fixture candidate population optimality gap",
                method="fixture basic bootstrap upper bound",
                coverage_status="fixture nominal approximation",
            )
            bounds.append(bound)
            return bound

        monkeypatch.setattr(script, "best_solution", best_solution)
        monkeypatch.setattr(script, "estimate_optimism", estimate)
        script.main(scale=scale)
        return SimpleNamespace(
            transfer=pd.read_csv(tmp_path / "cross_family_transfer.csv"),
            optimism=pd.read_csv(tmp_path / "cross_family_optimism.csv"),
            text=capsys.readouterr().out, calls=calls, scores=scores,
            chosen=chosen, candidate=candidate, bound=bounds[0],
        )

    return run


class TestTheReportSeparatesConfigurationsAndTargets:
    @pytest.mark.parametrize("scale", [False, True])
    def test_csv_configuration_comes_from_scored_arguments_and_actual_solutions(self, script,
                                                                              report, scale):
        result = report(scale=scale)
        transfer = result.transfer
        summary = result.optimism.iloc[0]
        required = {"model", "configuration", "target", "evaluation_unit", "selection_scope",
                    "oracle_role"}
        assert required <= set(transfer.columns)
        assert transfer.model.eq("ystwin.analysis.transfer.leave_one_stressor_out").all()
        assert summary.model == transfer.model.iloc[0]
        assert transfer.configuration.nunique() == 1
        fixed = json.loads(transfer.configuration.iloc[0])
        selected = json.loads(summary.selected_configuration)
        candidate = json.loads(summary.candidate_configuration)
        assert fixed["n_states"] == 3
        assert selected["n_states"] == result.chosen[0]
        assert candidate["n_states"] == result.candidate[0] == summary.candidate_states
        assert selected != candidate
        assert json.loads(summary.transfer_configuration) == fixed
        for config, solution in ((fixed, (3, "recommended pairs")),
                                 (selected, result.chosen), (candidate, result.candidate)):
            assert config["combinations"] == [list(pair) for pair in script.DESIGNS[solution[1]]]
            assert config["doses_ec50"] == list(result.calls.datasets[0].doses)
            assert config["replicates"] == result.calls.datasets[0].replicates
            assert config["observed"] == list(script.OBSERVED)
            assert config["reporters"] == script.READERS
            assert config["stressors"] == script.PANEL
            assert config["noise_cv"] == script.MEASURED_ACTIVITY_CV
            assert config["growth_rate_se_per_h"] == script.MEASURED_GROWTH_RATE_SE
        assert all(states == fixed["n_states"] for _, states, _, _ in result.calls.scores)
        assert all(data.design == tuple(script.RECOMMENDED_DESIGN)
                   for data in result.calls.datasets)
        assert all(observed == script.OBSERVED and only == script.PANEL
                   for _, _, observed, only in result.calls.scores)

    def test_design_contents_are_not_inferred_from_their_display_names(self, script,
                                                                     monkeypatch, report):
        monkeypatch.setattr(script, "DESIGNS", {"blocked": (("DTT", "MMS"),),
                                                "recommended pairs": ()})
        result = report(chosen=(4, "blocked"), candidate=(2, "recommended pairs"))
        summary = result.optimism.iloc[0]
        assert json.loads(summary.selected_configuration)["combinations"] == [["DTT", "MMS"]]
        assert json.loads(summary.candidate_configuration)["combinations"] == []
        assert json.loads(summary.selected_configuration)["n_states"] == 4
        assert summary.candidate_design == "recommended pairs"

    def test_split_and_bound_have_distinct_units_targets_and_selection_scopes(self, script, report):
        result = report()
        summary = result.optimism.iloc[0]
        assert "split_evaluation_unit" in summary.index
        assert "named synthetic family" in summary.split_evaluation_unit
        assert "not iid" in summary.split_evaluation_unit
        assert "independent reference family set" in summary.bound_resampling_unit
        assert "descriptive" in summary.split_target
        assert summary.bound_target == result.bound.target
        assert summary.bound_method == result.bound.method
        assert summary.bound_coverage_status == result.bound.coverage_warning
        assert "pending" in summary.selected_pipeline_optimism_status
        assert "independent evaluation" in summary.selected_pipeline_optimism_status
        assert summary.score_target == result.transfer.target.iloc[0]
        scope = json.loads(summary.split_selection_scope)
        assert scope == {"training_families": list(script.TRAINING),
                         "held_out_families": list(script.HELD_OUT),
                         "family_seed": 0, "data_seed": 0}
        scope = json.loads(summary.bound_selection_scope)
        assert scope["sampler"] == "ystwin.generator.families.sample_family"
        assert scope["family_registry"] == list(script.family_names())
        assert scope["family_weights"] == "uniform"
        assert scope["seed"] == result.calls.optimism[0][1]["seed"] == 0
        assert scope["candidate_fixed_for_references"] is True
        assert json.loads(summary.configuration_grid) == {
            "n_states": list(script.K_GRID),
            "designs": {name: [list(pair) for pair in pairs]
                        for name, pairs in script.DESIGNS.items()},
        }
        assert result.calls.searches[0][0] == list(script.TRAINING)
        assert result.calls.optimism[0][1] == {
            "n_candidate_domains": 4, "n_reference_domains": 2, "n_reference": 8,
            "alpha": 0.05, "n_bootstrap": 1000, "seed": 0,
        }

    def test_runtime_and_docstring_do_not_promote_non_rejection_or_the_separate_bound(self,
                                                                                   script, report):
        result = report()
        text = " ".join(result.text.split())
        assert "there is no effect" not in text
        assert "as much transfer information" not in text
        assert "baseline row is that script's number" not in text
        assert "in-sample optimism of the selection" not in text
        assert "fixed-candidate" in text
        assert "Selected-pipeline optimism: pending" in text
        assert "not biological validation" in text
        assert "raw gap samples" in text
        assert "-0.00123456789" in text
        assert "separate" in script.__doc__
        assert "pending" in script.__doc__

    def test_negative_scores_non_rejections_and_failed_null_rows_are_unchanged(self, script, report):
        result = report()
        expected = []
        for name, score in result.scores.items():
            null_median = np.nan if name == "loading_noise" else score + 0.1
            expected.append({
                "family": name, "perturbs": script.build_family(name, seed=0).perturbs,
                "transfer_r2": score, "oracle_r2": 0.2, "alignment": 0.5,
                "null_median": null_median,
                "p_value": np.nan if name == "loading_noise" else 1.0,
                "beats_null": False,
                "skill_over_null": skill_score_if_defined(score, null_median, perfect=1.0),
                "delta_vs_baseline": score - result.scores["baseline"],
            })
        expected = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(result.transfer[expected.columns], expected)
        assert result.transfer.loc[result.transfer.family == "loading_noise", "p_value"].isna().all()
        summary = result.optimism.iloc[0]
        assert summary.ucbog == result.bound.bound
        assert summary.mean_gap == pytest.approx(result.bound.mean_gap)
        assert summary.clipped_fraction == result.bound.clipped_fraction
        assert summary.baseline_transfer_r2 == result.scores["baseline"]
        assert summary.chosen_on_mean == pytest.approx(
            np.mean([result.scores[name] for name in script.TRAINING]))
        assert summary.held_out_mean == pytest.approx(
            np.mean([result.scores[name] for name in script.HELD_OUT]))
