"""The transfer script: the panel it scores, and the table it writes.

``scripts/run_transfer.py`` writes ``transfer_by_stressor.csv`` and ``module_recovery.csv``
and prints the per-design comparison ``docs/FINDINGS.md`` quotes. The script is almost all
``main``, so what is testable without a ten-minute run is the shape of the comparison
rather than its numbers: which stressors are scored, which are only ever co-dosed, that the
table written to disk is the table that was printed, and that the p-value floor the
docstring quotes is the one the draw count gives.

Loaded with ``YSTWIN_OUTPUTS`` redirected, because ``OUT`` binds at import time and
``main`` writes into it unconditionally.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import importlib.util
import itertools
import os
import pathlib
from queue import SimpleQueue
from threading import Barrier, get_ident
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_transfer.py"


@pytest.fixture(scope="module")
def script(tmp_path_factory):
    """The script as a module, with its output directory redirected away from outputs/."""
    redirect = tmp_path_factory.mktemp("transfer_outputs")
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(redirect)
    try:
        spec = importlib.util.spec_from_file_location("run_transfer", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


@pytest.fixture(scope="module")
def four_stressor_panel(script):
    """A cheap stand-in for the twelve-agent panel: four agents, four channels."""
    from ystwin.generator.panel_experiment import panel_dataset

    return panel_dataset(
        reporters=script.READERS[:4], stressors=script.PANEL[:4], doses=(0.5, 1.0),
        replicates=2, noise_cv=0.14, seed=0)


class TestOnlyThePanelIsScored:
    """``transfer_table`` exists to pass ``only=PANEL`` to every sweep. A sweep costs two
    latent fits per stressor and the median over its rows is the headline number, so
    widening the set silently changes both the runtime and the statistic."""

    def test_the_sweep_is_restricted_to_the_panel(self, script, monkeypatch):
        seen = {}
        monkeypatch.setattr(
            script, "leave_one_stressor_out",
            lambda data, n_states, observed, only: seen.update(
                n_states=n_states, observed=observed, only=only) or pd.DataFrame())

        script.transfer_table("data", n_states=3, observed=(0, 1, 2, 3))

        assert seen["only"] == script.PANEL
        assert seen == {"n_states": 3, "observed": (0, 1, 2, 3), "only": script.PANEL}

    def test_the_table_written_to_disk_is_the_table_that_was_printed(
            self, script, four_stressor_panel):
        """``main`` prints one sweep and then recomputes the same sweep to write the CSV.
        That is only sound because the sweep is deterministic; if it ever stopped being so,
        the published table and the printed one would differ with nothing to show for it."""
        first = script.transfer_table(four_stressor_panel, 2, (0, 1))
        second = script.transfer_table(four_stressor_panel, 2, (0, 1))

        pd.testing.assert_frame_equal(first, second)


class TestThePanelAndTheReadersAreRealNames:
    """Both are literals, and a typo in either surfaces only when the generator is called
    -- which in ``--scale`` is minutes into the run."""

    def test_every_panel_stressor_exists(self, script):
        from ystwin.generator.stress_panel import STRESSORS

        assert set(script.PANEL) <= set(STRESSORS)

    def test_every_reader_exists(self, script):
        from ystwin.generator.stress_panel import REPORTERS

        assert set(script.READERS) <= set(REPORTERS)

    def test_no_stressor_or_reader_is_listed_twice(self, script):
        """A repeat would add a duplicate channel column, and the revealed indices
        ``(0, 1, 2, 3)`` are positional."""
        assert len(set(script.PANEL)) == len(script.PANEL)
        assert len(set(script.READERS)) == len(script.READERS)

    def test_the_revealed_channels_exist_and_leave_something_to_predict(self, script):
        """``transfer_test`` refuses a sweep that reveals every channel. Four revealed of
        eight readers is what makes the number a prediction rather than a copy."""
        assert max((0, 1, 2, 3)) < len(script.READERS)
        assert len(script.READERS) > 4

    def test_the_buildable_panel_can_hold_out_a_channel_at_a_time(self, script):
        """``main`` hides each build channel in turn and reveals the rest. With one channel
        there would be nothing to reveal, and the ``min(k, len(seen))`` guard would ask for
        zero states."""
        from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD

        assert len(RECOMMENDED_BUILD.stress_reporters) >= 2


class TestTheDesignComparisonIsBetweenComparableThings:
    """The three designs are labelled with their own sizes read off the designs rather than
    typed, because both PANEL and RECOMMENDED_DESIGN have grown since they were last typed
    by hand. These are the invariants that keep the labels honest."""

    def test_all_pairs_is_every_unordered_pair_of_the_panel(self, script):
        assert script.ALL_PAIRS == list(itertools.combinations(script.PANEL, 2))

    def test_no_pair_doses_one_agent_against_itself(self, script):
        assert all(a != b for a, b in script.ALL_PAIRS)

    def test_the_pair_count_is_the_binomial_one(self, script):
        n = len(script.PANEL)

        assert len(script.ALL_PAIRS) == n * (n - 1) // 2

    def test_every_recommended_pair_is_a_real_stressor_pair(self, script):
        """``panel_dataset`` validates combination members against the whole registry, not
        against the stressors argument, so a name that is not a stressor raises but a name
        that is merely outside the panel does not."""
        from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN
        from ystwin.generator.stress_panel import STRESSORS

        members = {agent for pair in RECOMMENDED_DESIGN for agent in pair}
        assert members <= set(STRESSORS)


class TestThePValueFloorIsTheOneTheDrawCountBuys:
    """``NULL_DRAWS`` is quoted in the docstring as a floor of 1/41. A draw count changed
    without the docstring would leave a p-value reported against a floor it cannot reach."""

    def test_forty_draws_give_the_quoted_floor(self, script):
        assert script.NULL_DRAWS == 40
        assert 1 / (script.NULL_DRAWS + 1) == pytest.approx(0.024, abs=0.001)


@pytest.mark.integration
class TestTheWrittenTables:
    """What ``main`` leaves on disk. Magnitudes move whenever the generator is touched, so
    only the columns and the coverage are pinned."""

    @pytest.fixture(scope="class")
    def transfer_table(self):
        from ystwin import paths

        path = paths.outputs_dir() / "transfer_by_stressor.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not present; run scripts/run_transfer.py")
        return pd.read_csv(path)

    @pytest.fixture(scope="class")
    def recovery_table(self):
        from ystwin import paths

        path = paths.outputs_dir() / "module_recovery.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not present; run scripts/run_transfer.py")
        return pd.read_csv(path)

    def test_every_panel_stressor_has_a_row(self, script, transfer_table):
        assert set(transfer_table.held_out) == set(script.PANEL)

    def test_the_oracle_column_is_there_beside_the_transfer_one(self, transfer_table):
        """The oracle is a model allowed to see the held-out stressor. Without it a low
        transfer number cannot be told apart from a low ceiling."""
        assert {"r2", "oracle_r2", "alignment"} <= set(transfer_table.columns)

    def test_the_recovery_table_carries_its_own_shuffled_control(self, recovery_table):
        """A recovery number with no shuffled column beside it is not interpretable, and
        the control is a column of the same table so the two cannot drift apart."""
        assert "shuffled control" in recovery_table.columns

    def test_every_module_has_a_recovery_row(self, recovery_table):
        from ystwin.generator.stress_panel import MODULES

        assert set(recovery_table.module) == set(MODULES)


def test_invalid_null_budget_is_refused_before_generating_or_writing(script, tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(script, "panel_dataset", lambda **kwargs: calls.append(1))
    output = tmp_path / "not_created"
    with pytest.raises(ValueError, match="at least"):
        script.main(False, output_dir=output, draws=script._MIN_DRAWS - 1)
    assert not calls
    assert not output.exists()


@pytest.mark.parametrize("workers", [1, 3])
def test_transfer_only_cli_writes_the_same_nested_table(script, tmp_path, monkeypatch, workers):
    from ystwin.generator.panel_experiment import panel_dataset

    data = panel_dataset(reporters=script.READERS, stressors=script.PANEL[:4],
                         doses=(0.5, 1.0), replicates=2, noise_cv=0.14, seed=3)
    expected = script.transfer_table(data, None, (0, 1, 2, 3))
    monkeypatch.setattr(script, "panel_dataset", lambda **kwargs: data)
    script.main(False, output_dir=tmp_path, transfer_only=True, workers=workers)
    written = pd.read_csv(tmp_path / "transfer_by_stressor.csv")
    pd.testing.assert_frame_equal(expected, written)
    assert (tmp_path / "transfer_by_stressor.csv").read_bytes() == expected.to_csv(index=False).encode()
    assert not (tmp_path / "module_recovery.csv").exists()


@pytest.fixture(scope="module")
def synthetic_folds(script):
    from ystwin.generator.panel_experiment import PanelDataset

    rng = np.random.default_rng(203)
    labels = np.repeat(["heat", "DTT+H2O2", "DTT", "NaCl", "H2O2", "sorbitol"], 4)
    readings = rng.normal(size=(len(labels), 4))
    modules = rng.normal(size=(len(labels), len(script.MODULES)))
    data = PanelDataset(readings, labels, np.zeros(len(labels)), modules, script.READERS[:4])
    for values in (data.readings, data.labels, data.doses, data.modules):
        values.setflags(write=False)
    return data


@pytest.mark.parametrize("width", [None, 1, 2])
def test_parallel_single_fold_calls_match_every_sweep_column(script, synthetic_folds, width):
    observed = (1, 0)
    expected = script.leave_one_stressor_out(synthetic_folds, width, observed, only=script.PANEL)
    serial = script.transfer_table(synthetic_folds, width, observed)
    parallel = script.transfer_table(synthetic_folds, width, observed, workers=3)

    pd.testing.assert_frame_equal(expected, serial, check_exact=True)
    pd.testing.assert_frame_equal(expected, parallel, check_exact=True)
    assert expected.to_csv(index=False).encode() == parallel.to_csv(index=False).encode()
    assert list(parallel.columns) == [
        "held_out", "alignment", "r2", "oracle_r2", "n_states", "n_train", "n_test",
        "n_excluded", "selection", "baseline", "oracle_role",
    ]
    assert parallel.held_out.tolist() == ["DTT", "H2O2", "NaCl", "heat"]
    assert parallel.n_test.eq(4).all()
    paired = parallel.held_out.isin(["DTT", "H2O2"])
    assert parallel.loc[paired, "n_excluded"].eq(4).all()
    assert parallel.loc[paired, "n_train"].eq(16).all()
    assert parallel.loc[~paired, "n_excluded"].eq(0).all()
    assert parallel.loc[~paired, "n_train"].eq(20).all()
    assert parallel.oracle_role.eq("in-sample diagnostic, not a held-out prediction").all()
    assert parallel.baseline.eq("outer training channel mean").all()
    selection = ("nested leave-stressor-out on outer training only" if width is None
                 else "fixed width supplied before outer scoring")
    assert parallel.selection.eq(selection).all()


def test_parallel_transfer_preserves_nonfinite_scores(script, synthetic_folds):
    data = replace(synthetic_folds, readings=np.zeros_like(synthetic_folds.readings))
    expected = script.leave_one_stressor_out(data, 2, (0, 1), only=script.PANEL)
    parallel = script.transfer_table(data, 2, (0, 1), workers=3)
    pd.testing.assert_frame_equal(expected, parallel, check_exact=True)
    assert expected.to_csv(index=False).encode() == parallel.to_csv(index=False).encode()
    assert parallel[["r2", "oracle_r2"]].isna().all().all()
    assert parallel.alignment.eq(0.0).all()


@pytest.mark.parametrize("failure,width,observed,error", [
    ("missing_reading", 2, (0, 1), ValueError),
    ("missing_reading", None, (0, 1), ValueError),
    ("constant_readings", None, (0, 1), RuntimeWarning),
    ("too_many_states", 3, (0, 1), ValueError),
    ("all_channels_revealed", 2, (0, 1, 2, 3), ValueError),
])
def test_parallel_transfer_preserves_refusals(script, synthetic_folds, failure, width, observed, error):
    data = synthetic_folds
    if failure == "missing_reading":
        readings = data.readings.copy()
        readings[0, 0] = np.nan
        data = replace(data, readings=readings)
    elif failure == "constant_readings":
        data = replace(data, readings=np.zeros_like(data.readings))
    with pytest.raises(error) as expected:
        script.leave_one_stressor_out(data, width, observed, only=script.PANEL)
    for workers in (1, 3):
        with pytest.raises(error) as refused:
            script.transfer_table(data, width, observed, workers=workers)
        assert type(refused.value) is type(expected.value)
        assert str(refused.value) == str(expected.value)


@pytest.mark.parametrize("panel", [[], ["absent"], ["heat"], ["heat", "DTT", "DTT", "absent"]])
def test_parallel_transfer_preserves_requested_membership_and_order(
        script, synthetic_folds, monkeypatch, panel):
    monkeypatch.setattr(script, "PANEL", panel)
    expected = script.leave_one_stressor_out(synthetic_folds, 2, (0, 1), only=panel)
    parallel = script.transfer_table(synthetic_folds, 2, (0, 1), workers=3)
    pd.testing.assert_frame_equal(expected, parallel, check_exact=True)
    assert expected.to_csv(index=False).encode() == parallel.to_csv(index=False).encode()


def test_parallel_transfer_is_bounded_in_process_and_keeps_order(
        script, synthetic_folds, monkeypatch):
    panel = ["heat", "DTT", "H2O2", "absent"]
    monkeypatch.setattr(script, "PANEL", panel)
    leave_one = script.leave_one_stressor_out
    expected = leave_one(synthetic_folds, 2, (0, 1), only=panel)
    barrier = Barrier(3, timeout=10)
    completed = SimpleQueue()
    bounds = []

    def executor(*, max_workers):
        bounds.append(max_workers)
        return ThreadPoolExecutor(max_workers=max_workers)

    def inspect(data, n_states, observed, only):
        assert len(only) == 1
        barrier.wait()
        result = leave_one(data, n_states=n_states, observed=observed, only=only)
        completed.put((only[0], os.getpid(), get_ident()))
        return result

    monkeypatch.setattr(script, "ThreadPoolExecutor", executor)
    monkeypatch.setattr(script, "leave_one_stressor_out", inspect)
    parallel = script.transfer_table(synthetic_folds, 2, (0, 1), workers=99)
    results = [completed.get_nowait() for _ in range(3)]
    assert bounds == [3]
    assert {held for held, _, _ in results} == {"DTT", "H2O2", "heat"}
    assert {pid for _, pid, _ in results} == {os.getpid()}
    assert len({thread for _, _, thread in results}) == 3
    assert get_ident() not in {thread for _, _, thread in results}
    pd.testing.assert_frame_equal(expected, parallel, check_exact=True)


@pytest.mark.parametrize("workers", [0, -1, 1.5, True, np.bool_(True), "3"])
def test_invalid_workers_are_refused_before_generating_or_writing(script, tmp_path, monkeypatch, workers):
    calls = []
    monkeypatch.setattr(script, "panel_dataset", lambda **kwargs: calls.append(1))
    monkeypatch.setattr(script, "leave_one_stressor_out", lambda *args, **kwargs: calls.append(1))
    output = tmp_path / "not_created"
    with pytest.raises(ValueError, match="workers must be a positive integer"):
        script.main(False, output_dir=output, workers=workers)
    with pytest.raises(ValueError, match="workers must be a positive integer"):
        script.transfer_table(None, 2, (0, 1), workers=workers)
    assert not calls
    assert not output.exists()


def test_default_transfer_does_not_start_a_thread_pool(script, synthetic_folds, monkeypatch):
    def forbidden(**kwargs):
        pytest.fail("the default execution path must remain serial")

    monkeypatch.setattr(script, "ThreadPoolExecutor", forbidden)
    script.transfer_table(synthetic_folds, 2, (0, 1))


def test_main_routes_workers_through_design_null_and_build_sweeps(
        script, synthetic_folds, tmp_path, monkeypatch):
    calls, budgets = [], []

    def dataset(**kwargs):
        assert "workers" not in kwargs
        return synthetic_folds

    def table(data, n_states, observed, *, workers):
        calls.append((n_states, observed, workers))
        return pd.DataFrame({"held_out": ["DTT", "H2O2"], "r2": [0.2, 0.4],
                             "oracle_r2": [0.3, 0.5]})

    def compare(statistic, data, surrogate, *, n_draws, greater_is_better, label, seed):
        budgets.append((n_draws, seed))
        return SimpleNamespace(observed=statistic(data), null_median=0.0,
                               p_value=1.0, beats_null=False)

    monkeypatch.setattr(script, "panel_dataset", dataset)
    monkeypatch.setattr(script, "select_dimension", lambda centred, max_states: 2)
    monkeypatch.setattr(script, "select_dimension_by_transfer", lambda data, max_states, observed: 2)
    monkeypatch.setattr(script, "transfer_table", table)
    monkeypatch.setattr(script, "compare_to_null", compare)
    monkeypatch.setattr(script, "module_recovery", lambda data, n_states, shuffle=False: dict.fromkeys(script.MODULES, 0.0))
    script.main(False, output_dir=tmp_path, workers=3)
    build_size = len(script.RECOMMENDED_BUILD.stress_reporters)
    assert calls == ([(None, (0, 1, 2, 3), 3)] * 9
                     + [(None, tuple(i for i in range(build_size) if i != hidden), 3)
                        for hidden in range(build_size)])
    assert budgets == [(40, 0), (40, 0)]
    assert (tmp_path / "transfer_by_stressor.csv").exists()
    assert (tmp_path / "module_recovery.csv").exists()
