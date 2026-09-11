"""Runner regressions use fresh temporary tables, not historical reported calls."""

from __future__ import annotations

import importlib.util
import itertools
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

from ystwin.analysis.multiplicity import DoseResponse, benjamini_hochberg, holm
from ystwin.analysis.uncertainty import fold_change

_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "fold_multiplicity.py"


@pytest.fixture
def script():
    spec = importlib.util.spec_from_file_location("fold_multiplicity_runner", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _panel(construct="R", ratios=(1.2, 1.3, 1.4, 1.5), doses=(1.0,)):
    rows = []
    for plate, ratio in enumerate(ratios):
        for dose in (0.0, *doses):
            for well, factor in enumerate((0.85, 1.0, 1.15)):
                rows.append({
                    "plate": f"p{plate}", "construct": construct, "stressor": "DTT",
                    "dose_mM": dose, "well": f"d{dose}_w{well}",
                    "activity_late": factor * (1.0 if dose == 0 else ratio),
                })
    return pd.DataFrame(rows)


def _planned_family():
    frame = _panel(doses=(1.0, 2.0))
    unavailable = frame[frame.dose_mM == 1.0]
    return pd.concat([frame, *[
        unavailable.assign(dose_mM=float(dose), well=lambda d: "missing_" + d.well,
                           activity_late=-1.0)
        for dose in range(3, 25)
    ]], ignore_index=True)


@pytest.fixture
def run(script, monkeypatch, tmp_path):
    count = itertools.count()

    def invoke(frame, *, n_resamples=400, permutation_p=None, extra_args=()):
        directory = tmp_path / f"run-{next(count)}"
        directory.mkdir()
        frame.to_csv(directory / "sensor_characterisation.csv", index=False)
        monkeypatch.setattr(script.paths, "outputs_dir", lambda: directory)
        monkeypatch.setattr(script, "N_RESAMPLES", n_resamples)
        monkeypatch.setattr(script, "N_PERMUTATIONS", 20)
        monkeypatch.setattr(sys, "argv", [str(_SCRIPT), *map(str, extra_args)])

        def trend(readings, construct, **kwargs):
            subset = readings[readings.construct == construct]
            return DoseResponse(construct, "DTT", subset.dose_mM.nunique(),
                                subset.plate.nunique(), 0.5,
                                (permutation_p or {}).get(construct, 0.25),
                                kwargs["n_permutations"])

        monkeypatch.setattr(script, "dose_response_permutation", trend)
        assert script.main() == 0
        return pd.read_csv(directory / "fold_multiplicity.csv")

    return invoke


def _folds(table):
    return table[table.measurement == "per_dose_fold"].reset_index(drop=True)


def test_holm_and_bh_receive_the_entire_planned_family(script, run, monkeypatch):
    calls = []

    def record(name, correct):
        def wrapped(values, family_size=None):
            calls.append((name, family_size))
            return correct(values, family_size=family_size)
        return wrapped

    monkeypatch.setattr(script, "holm", record("holm", holm))
    monkeypatch.setattr(script, "benjamini_hochberg", record("bh", benjamini_hochberg))
    table = _folds(run(_planned_family(), n_resamples=1200))

    assert ("holm", 24) in calls
    assert ("bh", 24) in calls
    np.testing.assert_allclose(table.p_holm, holm(list(table.p_exact_sign), family_size=24),
                               equal_nan=True)
    np.testing.assert_allclose(table.p_bh,
                               benjamini_hochberg(list(table.p_exact_sign), family_size=24),
                               equal_nan=True)


def test_unestimable_cells_remain_visible_and_do_not_shrink_the_family(run):
    table = _folds(run(_planned_family(), n_resamples=1200))

    assert len(table) == 24
    assert set(table.family_size) == {24}
    assert table.estimable.sum() == 2
    refused = table[~table.estimable]
    assert refused.p_exact_sign.isna().all()
    assert refused.p_holm.isna().all()
    assert refused.p_bh.isna().all()
    assert refused.interval_method.str.contains("refused").all()


def test_sensor_holm_is_not_a_mislabeled_bonferroni_product(run):
    values = {"A": 0.001, "B": 0.01, "C": 0.02}
    frame = pd.concat([_panel(construct=name) for name in values], ignore_index=True)
    trends = run(frame, permutation_p=values).query(
        "measurement == 'dose_response_permutation'").sort_values("construct")

    np.testing.assert_allclose(trends.p_holm, holm(list(values.values()), family_size=3))
    np.testing.assert_allclose(trends.p_permutation, list(values.values()))
    np.testing.assert_allclose(trends.p_bh,
                               benjamini_hochberg(list(values.values()), family_size=3))
    assert set(trends.family_size) == {3}


def test_two_stage_family_intervals_are_not_refused_by_cluster_only_atoms(run):
    table = _folds(run(_planned_family(), n_resamples=1200))
    estimable = table[table.dose_mM <= 2]

    assert estimable.family_wise_resolvable.all()
    assert np.isfinite(estimable[["fw_low", "fw_high"]]).all().all()
    assert set(estimable.resampling) == {"two_stage"}
    assert estimable.coverage_status.str.contains("approximate").all()
    assert (estimable.fw_low <= estimable.low).all()
    assert (estimable.fw_high >= estimable.high).all()


def test_cluster_only_policy_is_still_reported_for_data_without_well_ids(run):
    table = _folds(run(_planned_family().drop(columns="well"), n_resamples=1200))
    estimable = table[table.dose_mM <= 2]

    assert not estimable.family_wise_resolvable.any()
    assert estimable.fw_low.isna().all() and estimable.fw_high.isna().all()
    assert set(estimable.resampling) == {"cluster_only"}
    assert estimable.family_wise_method.str.contains("cluster-only").all()


def test_monte_carlo_tail_resolution_is_checked_separately_from_plate_count(run):
    table = _folds(run(_panel(doses=(1.0, 2.0)), n_resamples=40))

    assert table.estimable.all()
    assert not table.family_wise_resolvable.any()
    assert table.fw_low.isna().all() and table.fw_high.isna().all()
    assert table.family_wise_method.str.contains("Monte Carlo").all()


def test_excluded_plate_cannot_reenter_the_interval_or_the_sign_test(run):
    frame = _panel()
    invalid = frame[frame.plate == "p0"].copy()
    invalid["plate"] = "excluded"
    invalid.loc[invalid.dose_mM == 1.0, "activity_late"] = [-1000.0, -1000.0, 100.0]
    original = _folds(run(frame)).iloc[0]
    combined = _folds(run(pd.concat([frame, invalid], ignore_index=True))).iloc[0]

    np.testing.assert_allclose([combined.fold, combined.low, combined.high],
                               [original.fold, original.low, original.high], rtol=1e-13)
    assert combined.n_plates == 5
    assert combined.n_plates_contributing == original.n_plates_contributing == 4
    assert combined.p_exact_sign == original.p_exact_sign
    assert combined.exact_sign_floor == original.exact_sign_floor


def test_repeated_summary_of_one_well_is_not_a_new_bootstrap_unit(run):
    frame = _panel()
    duplicate = frame[(frame.plate == "p0") & (frame.dose_mM == 1.0)].iloc[[0]]
    original = _folds(run(frame)).iloc[0]
    combined = _folds(run(pd.concat([frame, duplicate], ignore_index=True))).iloc[0]

    np.testing.assert_allclose([combined.fold, combined.low, combined.high],
                               [original.fold, original.low, original.high], rtol=1e-13)
    assert combined.n_wells == original.n_wells


def test_sign_test_uses_only_non_tied_usable_plates(run):
    frame = _panel(ratios=(2.0, 1.0, 1.0, 1.0))
    frame["activity_late"] = 1.0
    frame.loc[(frame.plate == "p0") & (frame.dose_mM == 1.0), "activity_late"] = 2.0
    row = _folds(run(frame)).iloc[0]

    assert row.p_exact_sign == 1.0
    assert row.n_plates_sign == 1
    assert row.plates_tied_one == 3
    assert row.exact_sign_floor == 1.0


def test_all_refused_folds_produce_rows_not_an_empty_table_crash(run):
    table = _folds(run(_panel(ratios=(-1.0, -2.0, -3.0, -4.0))))

    assert len(table) == 1
    assert not table.estimable.any()
    assert table.p_exact_sign.isna().all()
    assert not table.clears.any()
    assert not table.clears_family_wise.any()


def test_reported_intervals_come_from_the_shared_fold_change_api(run):
    frame = _panel()
    expected = fold_change(frame, "R", 1.0, n_resamples=400)
    row = _folds(run(frame)).iloc[0]

    np.testing.assert_allclose([row.fold, row.low, row.high],
                               [expected.point, expected.low, expected.high], rtol=1e-13)
    assert row.n_plates_contributing == expected.n_plates_contributing
    assert row.interval_method == expected.method
    assert row.coverage_status == expected.coverage_status


def test_empirical_bootstrap_tail_is_not_exported_as_a_confirmatory_p_value(run, capsys):
    table = run(_panel())
    output = capsys.readouterr().out

    assert "p_bootstrap" not in table.columns
    assert "p_permutation" in table.columns
    assert "not a confirmatory p-value" in output
    assert "not a resolution bound for" in output


def test_explicit_output_directory_does_not_redirect_the_readings(script, monkeypatch, tmp_path):
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.mkdir()
    destination.mkdir()
    _panel().to_csv(source / "sensor_characterisation.csv", index=False)
    monkeypatch.setattr(script.paths, "outputs_dir", lambda: source)
    monkeypatch.setattr(script, "N_RESAMPLES", 100)
    monkeypatch.setattr(script, "N_PERMUTATIONS", 20)
    monkeypatch.setattr(sys, "argv", [str(_SCRIPT), "--output-dir", str(destination)])

    assert script.main() == 0
    assert (destination / "fold_multiplicity.csv").is_file()
    assert not (source / "fold_multiplicity.csv").exists()


def test_fewer_than_three_contributing_plates_refuses_the_interval(run):
    row = _folds(run(_panel(ratios=(1.2, 1.3, -1.0, -2.0)))).iloc[0]

    assert row.n_plates == 4 and row.n_plates_contributing == 2
    assert not row.estimable and not row.family_wise_resolvable
    assert np.isnan(row.low) and np.isnan(row.high)
    assert row.n_plates_sign == 2
    assert row.p_exact_sign == pytest.approx(0.5)


def test_bootstrap_tail_budget_refusal_is_not_turned_into_a_fold_call(run):
    frame = _panel()
    frame.loc[frame.dose_mM > 0, "activity_late"] = [-10.0, -10.0, 21.0] * 4
    row = _folds(run(frame)).iloc[0]

    assert row.n_plates_contributing == 4
    assert np.isfinite(row.fold)
    assert not row.estimable and not row.clears
    assert not row.family_wise_resolvable and not row.clears_family_wise
    assert "tail budget" in row.interval_method


def test_cluster_only_resolution_is_per_cell_not_the_largest_plate_count(run):
    frame = _panel(ratios=(1.2, 1.3, 1.4, 1.5, 1.6), doses=(1.0, 2.0))
    frame = frame[~((frame.dose_mM == 1.0) & frame.plate.isin(["p3", "p4"]))]
    table = _folds(run(frame.drop(columns="well"))).set_index("dose_mM")

    assert table.loc[1.0, "n_plates_contributing"] == 3
    assert not table.loc[1.0, "estimable"]
    assert not table.loc[1.0, "family_wise_resolvable"]
    assert table.loc[2.0, "n_plates_contributing"] == 5
    assert table.loc[2.0, "estimable"] and table.loc[2.0, "family_wise_resolvable"]


def test_missing_control_and_undefined_trend_stay_in_their_families(script):
    frame = _panel()
    no_control = _panel(construct="no_control").query("dose_mM > 0")
    table = script.analyse(pd.concat([frame, no_control], ignore_index=True),
                           n_resamples=100, n_permutations=20)
    refused = table[table.construct == "no_control"]
    trends = table[table.measurement == "dose_response_permutation"]

    assert len(refused) == 2
    assert set(refused.family_size) == {2}
    assert refused.p_holm.isna().all() and refused.p_bh.isna().all()
    assert np.isnan(refused.loc[refused.measurement == "per_dose_fold", "fold"]).all()
    assert trends.loc[trends.construct == "no_control", "p_permutation"].isna().all()
    assert trends.loc[trends.construct == "R", "estimable"].all()
    assert not trends.loc[trends.construct == "no_control", "estimable"].any()
    np.testing.assert_allclose(trends.p_holm, holm(list(trends.p_permutation), family_size=2),
                               equal_nan=True)


def test_explicit_readings_and_destination_do_not_use_the_default_directory(script, monkeypatch,
                                                                           tmp_path):
    source = tmp_path / "input.csv"
    destination = tmp_path / "result"
    _panel().to_csv(source, index=False)

    def forbidden_default():
        pytest.fail("explicit input and output must not consult the default outputs directory")

    monkeypatch.setattr(script.paths, "outputs_dir", forbidden_default)
    assert script.main(["--readings", str(source), "--output-dir", str(destination),
                        "--resamples", "100", "--permutations", "20"]) == 0
    assert (destination / "fold_multiplicity.csv").is_file()
