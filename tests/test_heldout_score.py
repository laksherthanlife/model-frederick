"""The forward-prediction script, on constructed ladders where the answer is known.

``scripts/run_heldout_score.py`` produces the only predicted-then-measured number in the
project, and `docs/superseded/heldout-interpolation.md` records three retracted versions
of it. Every retraction was a property of the *procedure* rather than of the data: the
wrong functional form, one partition reported as a result, and partitions that differed in
training size as well as in which dose was held out. Those are the properties pinned here.

The tests fit real ``DoseFit`` objects to constructed ladders rather than checking
arithmetic, because each retraction turned on what the procedure did with a fit, not on a
sum. Only the last class touches the measured table, and it skips when that table is
absent.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_heldout_score.py"

LADDER = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0)
"""A ladder with a falling limb, so ``fit_dose_response`` can identify both limbs."""


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("run_heldout_score", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _activity(script, dose: float, ec50: float = 0.13) -> float:
    """One noiseless point of the repo's own biphasic form."""
    return float(script._biphasic_model(dose, ec50, 3.0, 2000.0, 300.0))


def _inventory(script, ladders) -> pd.DataFrame:
    """A characterisation-shaped frame: three wells per dose on each of two plates.

    ``stressor`` has to agree with ``qpcr.STRESSOR_FOR_CONSTRUCT`` or the split module
    refuses the frame before any of this is exercised.
    """
    rows = []
    for construct, (doses, ec50) in ladders.items():
        for plate in ("20260722", "20260803"):
            for dose in doses:
                for well in ("A1", "A2", "A3"):
                    rows.append({
                        "plate": plate, "construct": construct, "stressor": "DTT",
                        "dose_mM": dose, script.TARGET: _activity(script, dose, ec50),
                        "well": f"{plate[-2:]}{construct[-1]}{well}",
                    })
    return pd.DataFrame(rows)


def _manifest(assignments: dict, kind: str = "interpolation", seed: int = 0) -> pd.DataFrame:
    """One manifest row per (plate, construct, dose) group, as ``make_splits.py`` writes."""
    return pd.DataFrame([
        {"split_kind": kind, "seed": seed, "plate_key": plate, "construct": construct,
         "dose_mM": dose, "assignment": side, "partition_hash": "0123456789abcdef"}
        for (plate, construct, dose), side in assignments.items()
    ])


def _one_construct_split(script, held_out_dose: float, held_out_activity: float,
                         train_doses=(0.0, 0.1, 0.2, 0.5, 1.0, 2.0)):
    """A manifest and readings pair holding out a single well of a single construct."""
    assignments, readings = {}, []
    for dose in train_doses:
        assignments[("20260722", "UPRE1", dose)] = "train"
        readings.append({"plate_key": "20260722", "construct": "UPRE1", "dose_mM": dose,
                         script.TARGET: _activity(script, dose)})
    assignments[("20260722", "UPRE1", held_out_dose)] = "test"
    readings.append({"plate_key": "20260722", "construct": "UPRE1",
                     "dose_mM": held_out_dose, script.TARGET: held_out_activity})
    return _manifest(assignments), pd.DataFrame(readings)


class TestAPerfectBaselineDoesNotKillTheReport:
    """``skill_score`` raises when the baseline is already perfect, and this report scores
    five partitions in one pass. One lucky nearest-dose baseline must cost that partition's
    skill number, not every other row in the table."""

    def test_the_scorer_holds_the_guarded_skill_score(self, script):
        """The guarded variant returns None where the strict one raises ValueError."""
        assert script.skill_score_if_defined(5.0, 0.0, perfect=0.0) is None

    def test_a_perfect_baseline_leaves_the_other_numbers_standing(self, script):
        """Held-out activity equal to the nearest training well makes the baseline exact.
        The skill is then undefined, and the RMSEs still have to be reported."""
        manifest, readings = _one_construct_split(
            script, held_out_dose=0.6, held_out_activity=_activity(script, 0.5))

        result = script._score("interpolation", manifest, readings, seed=0)

        assert result["rmse_nearest_dose"] == 0.0
        assert result["skill_biphasic"] is None
        assert np.isfinite(result["rmse_biphasic"])


class TestIdentifiabilityIsARealPartition:
    """``DoseFit.identifiable`` is the calibration module's own verdict, with the reason
    attached. The scorer ignored it once, and scoring refused fits alongside vouched ones
    is what turned a refusal into an apparent failure of the model."""

    @staticmethod
    def _fits(script):
        from ystwin.generator.panel_calibration import fit_dose_response

        doses = np.asarray(LADDER, dtype=float)
        biphasic = np.asarray([_activity(script, d) for d in LADDER])
        # Saturating with no falling limb: the lethal dose is unconstrained, so the
        # module refuses the fit rather than reporting the bound the optimiser ran to.
        saturating = 300.0 + 2000.0 * doses / (0.13 + doses)
        return {"vouched": fit_dose_response(doses, biphasic),
                "refused": fit_dose_response(doses, saturating)}

    def test_a_ladder_with_no_falling_limb_is_refused(self, script):
        fits = self._fits(script)

        assert fits["refused"].identifiable is False
        assert fits["vouched"].identifiable is True

    def test_the_split_partitions_the_fits_rather_than_keeping_them_all(self, script):
        good, refused = script._split_by_identifiability(self._fits(script))

        assert (good, refused) == ({"vouched"}, {"refused"})

    def test_a_fit_that_does_not_vouch_for_itself_is_refused(self, script):
        """The lookup defaults to False, so an object with no verdict cannot slip through
        into the headline number by silence."""
        good, refused = script._split_by_identifiability({"unknown": object()})

        assert (good, refused) == (set(), {"unknown"})


class TestTheInterpolationSweepExhaustsTheInteriorRungs:
    """v3 of this result sampled three seeds that withheld two, three and four doses, so
    the partitions differed in training size as well as in which dose was held out, and the
    spread across them attributed neither. The sweep exists to exhaust the interior one rung
    at a time; these are the assertions that fail if it ever samples again."""

    @staticmethod
    def _balanced(script):
        return _inventory(script, {"UPRE1": (LADDER, 0.13), "UPRE2": (LADDER, 0.15)})

    def test_every_interior_rung_is_scored_exactly_once(self, script):
        inventory = self._balanced(script)
        interior = sorted(inventory.dose_mM.unique())[1:-1]

        frame = script._exhaustive_interpolation(inventory)

        assert frame.held_out_dose.tolist() == interior

    def test_the_ends_of_the_ladder_are_never_held_out(self, script):
        """Withholding the top rung is the extrapolation split, and it licenses a different
        claim. Withholding the bottom one leaves no basal to anchor the fit."""
        inventory = self._balanced(script)

        frame = script._exhaustive_interpolation(inventory)

        assert set(frame.held_out_dose).isdisjoint({inventory.dose_mM.min(),
                                                    inventory.dose_mM.max()})

    def test_each_partition_withholds_one_whole_rung_and_nothing_else(self, script):
        """The comparability the retraction was about: on a balanced ladder every partition
        tests exactly the wells at its own held-out dose, so training sizes match and a
        difference in skill is about the dose rather than about how much data was left."""
        inventory = self._balanced(script)

        frame = script._exhaustive_interpolation(inventory)

        at_dose = [int((inventory.dose_mM == d).sum()) for d in frame.held_out_dose]
        assert frame.n_test.tolist() == at_dose

    def test_a_rung_the_split_module_refuses_is_recorded_rather_than_dropped(self, script):
        """A rung that is one construct's top dose cannot be an interpolation fold. Dropping
        such a row silently would make a sampled sweep indistinguishable from an exhausted
        one -- so the refusal is carried, with its reason."""
        ragged = _inventory(script, {"UPRE1": (LADDER, 0.13),
                                     "UPRE2": (LADDER[:-1], 0.15)})

        frame = script._exhaustive_interpolation(ragged)

        refused = frame[frame.skill.isna()]
        assert refused.held_out_dose.tolist() == [2.0]
        assert "not available" in refused.note.iloc[0]


class TestTheBaselineSeesOnlyTrainingRows:
    """"Carry the nearest measured dose forward" is the baseline the whole claim is stated
    against. A baseline that could see the held-out well would score zero error and turn
    every skill negative; one pooled across constructs would be a different baseline than
    the one reported."""

    def test_the_nearest_dose_baseline_cannot_see_the_held_out_activity(self, script):
        """0.3 mM sits nearer 0.2 than 0.5, so the baseline is the 0.2 mM training well --
        and the absurd held-out value must appear only as error, never as the prediction."""
        manifest, readings = _one_construct_split(
            script, held_out_dose=0.3, held_out_activity=12345.0)

        result = script._score("interpolation", manifest, readings, seed=0)

        assert result["rmse_nearest_dose"] == pytest.approx(
            abs(_activity(script, 0.2) - 12345.0))

    def test_the_nearest_dose_comes_from_the_same_construct(self, script):
        """A second construct measured at exactly the held-out dose is a zero-distance
        neighbour if the baseline pools. It must be ignored: the fit is per construct, so
        the baseline it is scored against has to be too."""
        manifest, readings = _one_construct_split(
            script, held_out_dose=0.3, held_out_activity=12345.0)
        manifest = pd.concat([manifest, _manifest(
            {("20260722", "UPRE2", d): "train" for d in (0.0, 0.1, 0.3, 0.5, 1.0, 2.0)})])
        readings = pd.concat([readings, pd.DataFrame([
            {"plate_key": "20260722", "construct": "UPRE2", "dose_mM": d,
             script.TARGET: 99999.0} for d in (0.0, 0.1, 0.3, 0.5, 1.0, 2.0)])])

        result = script._score("interpolation", manifest, readings, seed=0)

        assert result["rmse_nearest_dose"] == pytest.approx(
            abs(_activity(script, 0.2) - 12345.0))


class TestAnUnscoreableSplitNamesItsCause:
    """"0 test rows" reads as a data gap. Every real cause here is structural, and a report
    that says which one it hit is the difference between a missing result and a known one."""

    def test_a_held_out_plate_absent_from_the_readings_is_named(self, script):
        """Two of the four plates have no reporter-channel blank, so characterisation never
        produced rows for them -- permanently, not as a join fault."""
        manifest = _manifest({("20260803", "UPRE1", 0.5): "test",
                              ("20260722", "UPRE1", 0.1): "train"})
        readings = pd.DataFrame([{"plate_key": "20260722", "construct": "UPRE1",
                                  "dose_mM": 0.1, script.TARGET: 900.0}])

        note = script._score("interpolation", manifest, readings, seed=0)["note"]

        assert "20260803" in note and "absent from the readings table" in note

    def test_a_split_that_withholds_whole_constructs_says_which(self, script):
        """The model class fits per construct, so a heldout_construct split leaves every
        test row unpredictable. That is a limit of the model, not a gap in the data."""
        manifest, readings = _one_construct_split(
            script, held_out_dose=0.6, held_out_activity=1000.0)
        manifest = manifest[manifest.assignment == "train"]
        readings = readings[readings.dose_mM != 0.6]
        manifest = pd.concat([manifest, _manifest(
            {("20260722", "UPRE2", d): "test" for d in (0.1, 0.2, 0.5)})])
        readings = pd.concat([readings, pd.DataFrame([
            {"plate_key": "20260722", "construct": "UPRE2", "dose_mM": d,
             script.TARGET: _activity(script, d, 0.2)} for d in (0.1, 0.2, 0.5)])])

        note = script._score("interpolation", manifest, readings, seed=0)["note"]

        assert "UPRE2" in note and "no overlap" in note

    def test_a_training_side_too_small_to_fit_says_so(self, script):
        """Below six training rows neither form can be fitted. Naming the threshold is what
        distinguishes it from the two causes above."""
        manifest, readings = _one_construct_split(
            script, held_out_dose=0.3, held_out_activity=900.0,
            train_doses=(0.0, 0.1, 0.5))

        note = script._score("interpolation", manifest, readings, seed=0)["note"]

        assert "below the 6 needed" in note


@pytest.mark.integration
class TestTheSweepOnTheMeasuredLadder:
    """The published table comes from the real characterisation output. Its magnitudes move
    whenever a plate is re-read, so only the sweep's coverage is pinned here."""

    @pytest.fixture(scope="class")
    def readings(self):
        from ystwin import paths

        path = paths.outputs_dir() / "sensor_characterisation.csv"
        if not path.exists():
            pytest.skip(f"characterisation table not present: run "
                        f"scripts/run_sensor_characterisation.py to write {path.name}")
        return pd.read_csv(path)

    def test_every_interior_rung_of_the_measured_ladder_is_accounted_for(self, script,
                                                                        readings):
        """Complete rather than sampled, on the actual ladder -- including the rung the
        split module refuses, which is carried as a note."""
        interior = sorted(readings.dose_mM.dropna().unique())[1:-1]

        frame = script._exhaustive_interpolation(readings)

        assert frame.held_out_dose.tolist() == interior


def test_baseline_and_model_choices_are_blind_to_outer_test_activity(script):
    manifest, readings = _one_construct_split(
        script, held_out_dose=0.3, held_out_activity=1000.0)
    first = script._score("interpolation", manifest, readings, seed=0)
    changed = readings.copy()
    changed.loc[changed.dose_mM == 0.3, script.TARGET] = 1e8
    second = script._score("interpolation", manifest, changed, seed=0)
    assert first["selected_baselines"] == second["selected_baselines"]
    assert first["selected_models"] == second["selected_models"]
    assert first["baseline_selection"] == "inner leave-dose-out on outer training only"


def test_nearest_dose_averages_training_wells_instead_of_selecting_one(script):
    manifest, readings = _one_construct_split(
        script, held_out_dose=0.3, held_out_activity=1000.0)
    added = readings[readings.dose_mM == 0.2].copy()
    added[script.TARGET] += 1000.0
    readings = pd.concat([readings, added], ignore_index=True)
    result = script._score("interpolation", manifest, readings, seed=0)
    expected = abs(_activity(script, 0.2) + 500.0 - 1000.0)
    assert result["rmse_nearest_dose"] == pytest.approx(expected)


def test_product_target_delegates_to_the_shared_nested_validator(script, tmp_path):
    from ystwin.pathway.flux import carotenoid_measurements, score_product_validation
    from ystwin.pathway.spec import load_pathway

    expected = score_product_validation(load_pathway("beta_carotene"), carotenoid_measurements(), draws=0)
    actual = script._score_products(tmp_path, n_draws=0)
    pd.testing.assert_frame_equal(actual, expected)
    assert (tmp_path / "product_prediction.csv").exists()
    assert (tmp_path / "product_prediction_selection.csv").exists()
    comparison = pd.read_csv(tmp_path / "product_matched_entry_comparisons.csv")
    assert set(comparison.groupby("branch").size()) == {6}
    assert set(comparison.groupby("branch").strain.nunique()) == {3}


def test_uncertainty_target_writes_actual_coverage_and_simulation_assumptions(script, tmp_path):
    table = script._score_uncertainty(tmp_path, n_trials=3, n_resamples=40, confidence=0.90)
    written = pd.read_csv(tmp_path / "fold_bootstrap_coverage.csv")
    assert len(table) == len(written) == 6
    assert set(written.nominal_coverage) == {0.90}
    assert np.allclose(written.actual_coverage, written.n_covered / written.n_issued)
    assert written.coverage_scope.str.contains("not real-data calibration").all()


def test_conflicting_manifest_assignments_are_not_silently_overwritten(script):
    manifest, readings = _one_construct_split(script, 0.3, 1000.0)
    conflict = manifest.iloc[[0]].copy()
    conflict["assignment"] = "test"
    with pytest.raises(ValueError, match="conflicting"):
        script._score("interpolation", pd.concat([manifest, conflict]), readings, seed=0)
