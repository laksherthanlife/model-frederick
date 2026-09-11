"""The training script: what it refuses to train on, and how it counts a recovered module.

``scripts/run_training.py`` writes the four ``training_*`` tables and the three saved
models, so a silent change to its gate or to its 0.25 recovery threshold moves published
numbers with nothing to point at. The properties pinned here are the ones that would move
them without failing anything else:

* the gate refuses a plate that could not have been measured, by raising rather than by
  training anyway on a warning;
* the reference for the posterior-predictive gate is absent rather than invented when the
  measured table is not on disk;
* ``transfer_table`` counts a module as recovered strictly above 0.25 and averages only
  over the modules the held-out stressor actually drives.

The module is loaded with ``YSTWIN_OUTPUTS`` pointed at a temporary directory, because
``OUT`` is bound at import time and a test that imported it otherwise would be one stray
``main()`` away from overwriting the tracked tables in ``outputs/``.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import os
import pathlib

import numpy as np
import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_training.py"


@pytest.fixture(scope="module")
def script(tmp_path_factory):
    """The script as a module, with its output directory redirected away from outputs/."""
    redirect = tmp_path_factory.mktemp("training_outputs")
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(redirect)
    try:
        spec = importlib.util.spec_from_file_location("run_training", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


@pytest.fixture(scope="module")
def small_panel(script):
    """One cheap plate: four channels over the whole stressor set, two doses."""
    from ystwin.generator.stress_panel import STRESSORS

    from ystwin.generator.panel_experiment import panel_dataset

    return panel_dataset(
        reporters=["STRE-general", "UPRE-ER", "TRX2-oxidative", "HSE-heat"],
        stressors=list(STRESSORS), doses=(0.5, 1.0), replicates=2,
        noise_cv=script.MEASURED_ACTIVITY_CV, seed=0)


class TestTheGateRefusesRatherThanTrainingAnyway:
    """``gate`` is the only thing standing between a generator change and a model fitted to
    readings no culture could produce. It has to stop the run, not annotate it."""

    def test_a_reading_below_the_constitutive_floor_stops_the_run(self, script,
                                                                  small_panel):
        """A negative reading means a blank subtraction gone wrong or a perturbation that
        broke the floor. Either way the plate is not evidence, and training on it would
        produce a saved model with no way to tell."""
        readings = small_panel.readings.copy()
        readings[0, 0] = -1.0
        broken = dataclasses.replace(small_panel, readings=readings)

        with pytest.raises(SystemExit) as refusal:
            script.gate(broken, None, "a plate with a negative well")

        assert "refused" in str(refusal.value)

    def test_a_plate_that_could_have_been_measured_passes(self, script, small_panel):
        """The counterpart, so the test above is about the negative well and not about the
        gate refusing everything."""
        script.gate(small_panel, None, "an ordinary plate")

    def test_a_refused_plate_leaves_no_table_and_no_saved_model(self, script,
                                                                small_panel,
                                                                monkeypatch):
        """The gate runs before anything is written, and this is what says so. Four tables
        and three ``.npz`` models come out of one run; a refusal that arrived after even
        one of them would leave a file that reads exactly like a good one."""
        readings = small_panel.readings.copy()
        readings[0, 0] = -1.0
        monkeypatch.setattr(script, "build",
                            lambda *args, **kwargs: dataclasses.replace(
                                small_panel, readings=readings))

        with pytest.raises(SystemExit):
            script.main(quick=True)

        assert not list(pathlib.Path(script.OUT).iterdir())

    def test_a_reference_the_matched_slice_does_not_resemble_stops_the_run(self, script,
                                                                          small_panel):
        """The posterior-predictive gate compares a matched one-promoter, one-agent slice
        against real readings. A reference whose spread the simulator cannot reproduce is a
        refusal, not a warning -- that gate is the only one that looks at measured data."""
        implausible = np.array([[1.0], [1.0], [1000.0]])

        with pytest.raises(SystemExit) as refusal:
            script.gate(small_panel, implausible, "a plate against an alien reference")

        assert "refused" in str(refusal.value)


class TestTheGateReferenceIsAbsentRatherThanInvented:
    """``reference_readings`` reads the real characterisation table. When it is not on disk
    the answer is "no reference", never a stand-in: a fabricated reference would make the
    posterior-predictive gate pass against the generator's own assumptions, which is the
    one thing that gate exists not to do."""

    def test_no_measured_table_yields_no_reference(self, script, tmp_path, monkeypatch):
        monkeypatch.setattr(script, "MEASURED", tmp_path / "not_written_yet.csv")

        assert script.reference_readings() is None

    def test_the_reference_is_the_two_upre_constructs_paired_by_plate_and_dose(
            self, script, tmp_path, monkeypatch):
        """One row per (plate, dose) with both constructs present. A third construct must
        not widen it and an unpaired dose must not enter it half-filled, because the gate
        compares the spread of this block against a simulated one of the same shape."""
        measured = tmp_path / "sensor_characterisation.csv"
        pd.DataFrame([
            {"plate": "20260722", "construct": "UPRE1", "dose_mM": 0.5,
             "activity_late": 900.0},
            {"plate": "20260722", "construct": "UPRE2", "dose_mM": 0.5,
             "activity_late": 800.0},
            {"plate": "20260722", "construct": "TRX2", "dose_mM": 0.5,
             "activity_late": 700.0},
            # 1.0 mM has UPRE1 only, so the pair is incomplete and the row is dropped.
            {"plate": "20260722", "construct": "UPRE1", "dose_mM": 1.0,
             "activity_late": 1200.0},
        ]).to_csv(measured, index=False)
        monkeypatch.setattr(script, "MEASURED", measured)

        reference = script.reference_readings()

        assert reference.shape == (1, 2)
        assert sorted(reference[0]) == [800.0, 900.0]


class TestTheRecoveryThresholdIsStrict:
    """``modules_recovered`` is a count over ``v > 0.25`` and it is quoted in
    ``training_transfer_*.csv``. A module sitting exactly on the threshold must not be
    counted, and a change from ``>`` to ``>=`` would move every row in those tables."""

    @staticmethod
    def _with_scores(script, monkeypatch, scores):
        monkeypatch.setattr(script, "module_transfer",
                            lambda data, name, n_states: dict(scores))
        return script.transfer_table(data=None, n_states=3, tag="stub", quiet=True)

    def test_a_module_exactly_on_the_threshold_is_not_recovered(self, script, monkeypatch):
        from ystwin.generator.stress_panel import MODULES

        scores = {m: 0.0 for m in MODULES}
        scores["UPR"] = 0.25

        frame = self._with_scores(script, monkeypatch, scores)

        assert set(frame.modules_recovered) == {0}

    def test_a_module_just_above_the_threshold_is_recovered(self, script, monkeypatch):
        from ystwin.generator.stress_panel import MODULES

        scores = {m: 0.0 for m in MODULES}
        scores["UPR"] = 0.2500001

        frame = self._with_scores(script, monkeypatch, scores)

        assert set(frame.modules_recovered) == {1}

    def test_the_mean_is_taken_over_the_modules_the_stressor_drives(self, script,
                                                                   monkeypatch):
        """``mean_over_driven`` is the headline column. Averaging over all 24 modules
        instead of the held-out agent's own targets would dilute it toward zero with the
        modules nothing in that row is asking about."""
        from ystwin.generator.stress_panel import MODULES, STRESSORS

        scores = {m: 0.0 for m in MODULES}
        for module in STRESSORS["DTT"].targets:
            scores[module] = 0.8

        frame = self._with_scores(script, monkeypatch, scores)

        row = frame[frame.held_out == "DTT"].iloc[0]
        assert row.mean_over_driven == pytest.approx(0.8)

    def test_every_stressor_names_targets_the_scorer_will_have_scores_for(self, script):
        """``mean_over_driven`` averages ``[scored[m] for m in scored if m in targets]``.
        A target name that is not a module makes that list empty, and the mean of an empty
        slice is a warning and a NaN rather than a refusal -- so the invariant that keeps
        the column meaningful is checked here instead."""
        from ystwin.generator.stress_panel import MODULES, STRESSORS

        for name, stressor in STRESSORS.items():
            unknown = [m for m in stressor.targets if m not in MODULES]
            assert not unknown, f"{name} targets {unknown}, which are not modules"

    def test_one_row_per_stressor_in_a_fixed_order(self, script, monkeypatch):
        """The table is written to CSV, so its row order is part of the artifact."""
        from ystwin.generator.stress_panel import MODULES, STRESSORS

        frame = self._with_scores(script, monkeypatch, {m: 0.5 for m in MODULES})

        assert frame.held_out.tolist() == sorted(STRESSORS)


class TestWhatAPlateGivesForFree:
    """``deployable`` is what the model is actually fitted on. Both steps it performs were
    once absent, and each is invisible in the outputs if it goes missing again."""

    def test_growth_arrives_as_one_extra_channel(self, script, small_panel):
        """Growth rate comes off the optical density, so it costs no fluorophore. A fit
        that lost it would still run, on one channel less."""
        fitted = script.deployable(small_panel)

        assert fitted.readings.shape[1] == small_panel.readings.shape[1] + 1

    def test_replicates_are_pooled_before_fitting(self, script, small_panel):
        """The latent model gives every well its own state, so unpooled replicates buy
        parameters instead of averaging noise."""
        fitted = script.deployable(small_panel)

        assert fitted.readings.shape[0] < small_panel.readings.shape[0]
        assert len(set(fitted.labels)) == len(set(small_panel.labels))


class TestTheTrainingCornerIsLeftDeliberately:
    """Everything before this script trained in exponential glucose at 30 degrees, where a
    raised general stress response always means a stressor. The context list is what stops
    that, and it is a plain module constant that a refactor could quietly shorten."""

    def test_the_contexts_are_distinct_and_include_the_old_corner(self, script):
        from ystwin.generator.context import CultureContext

        assert CultureContext() in script.CONTEXTS
        assert len(set(script.CONTEXTS)) == len(script.CONTEXTS)

    def test_more_than_one_growth_phase_is_trained_on(self, script):
        """A stationary plate reads a raised ESR with no stressor present. If every context
        were exponential the model could not tell the two apart."""
        phases = {c.growth_phase for c in script.CONTEXTS}

        assert len(phases) > 1

    def test_the_wide_panel_adds_the_ratiometric_sensor_without_repeating_a_promoter(
            self, script):
        """``WIDE`` is a list, so a sensor that were also a transcriptional reporter would
        appear twice and be read twice."""
        from ystwin.generator.stress_panel import REPORTERS, Kind

        assert len(set(script.WIDE)) == len(script.WIDE)
        assert REPORTERS["roGFP2-Grx1"].kind is Kind.RATIOMETRIC


@pytest.mark.integration
class TestTheWrittenTablesSayWhichModeProducedThem:
    """``--quick`` runs four doses and three replicates against six and six, which moves the
    numbers enough to change conclusions, and a table with no mode recorded was quoted as a
    full run once already. Every tracked training table carries the stamp or it is not
    interpretable."""

    @pytest.fixture(params=["training_transfer_wide", "training_transfer_build",
                            "training_transfer_three_sensor", "training_recovery"])
    def table(self, request):
        from ystwin import paths

        path = paths.outputs_dir() / f"{request.param}.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not present; run scripts/run_training.py")
        return pd.read_csv(path)

    def test_the_run_mode_is_stamped_on_every_row(self, table):
        assert set(table.columns) >= {"run_mode", "n_doses", "replicates"}
        assert set(table.run_mode) <= {"quick", "full"}

    def test_the_stamp_is_one_value_per_table(self, table):
        """A table concatenated from two runs would carry two modes and could not be
        quoted as either."""
        for column in ("run_mode", "n_doses", "replicates"):
            assert table[column].nunique() == 1

    def test_the_dose_count_agrees_with_the_mode_it_claims(self, table):
        """The stamp is written from the same ``doses`` tuple the run used, so a mismatch
        means the columns were filled in from somewhere else."""
        expected = {"quick": 4, "full": 6}[table.run_mode.iloc[0]]

        assert int(table.n_doses.iloc[0]) == expected
