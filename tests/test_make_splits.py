"""The split manifest, and the one way writing it can destroy what it replaces.

``scripts/make_splits.py`` builds two halves: a simulated panel that needs nothing, and a
real biosensor inventory read off the NewProtocol exports. When the exports are not
reachable the real half is simply absent -- and the run used to write the panel-only
manifest over the tracked one and exit 0, losing 714 of 1650 rows while printing "wrote".

The loss was silent at the point it happened and loud somewhere else entirely:
``run_heldout_score.py`` selects ``dataset == "real_biosensor"`` and then reads
``manifest.plate``, so it died with ``'DataFrame' object has no attribute 'plate'`` -- a
column error, three scripts away from a missing environment variable.
"""

from __future__ import annotations

import importlib.util
import os
import pathlib
import subprocess
import sys
from dataclasses import replace

import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "make_splits.py"
_REAL_TABLE = _REPO / "outputs" / "sensor_characterisation.csv"
_REAL_KEY = ["plate", "construct", "stressor", "dose_mM", "well"]
_GROUP_KEY = ["plate", "construct", "dose_mM"]


@pytest.fixture
def script(tmp_path, monkeypatch):
    """The script as a module, with ``OUT`` bound to a scratch directory.

    ``OUT`` binds at import time, so the redirect has to precede the import; the module is
    re-imported per test because each test gives it a different manifest to look at.
    """
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
    spec = importlib.util.spec_from_file_location("make_splits", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _manifest(datasets: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"dataset": datasets, "split_kind": ["interpolation"] * len(datasets),
                         "n_rows": [1] * len(datasets)})


class TestItWillNotReplaceRealRowsWithNone:
    def test_a_manifest_holding_real_rows_stops_the_run(self, script):
        _manifest(["simulated_panel", "real_biosensor"]).to_csv(
            script.OUT / "split_manifest.csv", index=False)

        with pytest.raises(SystemExit, match="YSTWIN_PLATES"):
            script._refuse_to_shrink_the_manifest(panel_only=False)

    def test_the_refusal_names_the_file_it_is_protecting(self, script):
        _manifest(["real_biosensor"]).to_csv(script.OUT / "split_manifest.csv", index=False)

        with pytest.raises(SystemExit, match="split_manifest.csv"):
            script._refuse_to_shrink_the_manifest(panel_only=False)


class TestExplicitPanelOnlyKeepsTheDocumentedSimulatedPath:
    """Historical rationale for the simulated-only path:

    The guard has to stay narrow. A gate that fired whenever the plates were missing
    would make a fresh checkout unable to build the simulated half at all, which is a
    working path this repository documents.

    The path remains supported, but the caller must now explicitly request --panel-only.
    Destination contents cannot declare the population a fresh full-scope run requested.
    """

    def test_a_fresh_checkout_with_no_manifest_is_allowed_with_panel_only(self, script):
        assert not (script.OUT / "split_manifest.csv").exists()

        script._refuse_to_shrink_the_manifest(panel_only=True)

    def test_a_manifest_that_never_had_real_rows_is_allowed_with_panel_only(self, script):
        _manifest(["simulated_panel", "simulated_panel"]).to_csv(
            script.OUT / "split_manifest.csv", index=False)

        script._refuse_to_shrink_the_manifest(panel_only=True)

    def test_panel_only_says_the_shrink_is_intended(self, script):
        _manifest(["simulated_panel", "real_biosensor"]).to_csv(
            script.OUT / "split_manifest.csv", index=False)

        script._refuse_to_shrink_the_manifest(panel_only=True)

    def test_an_empty_file_is_not_treated_as_a_loss_with_panel_only(self, script):
        """A zero-byte manifest is what a run that died mid-write leaves. There is nothing
        in it to protect, and raising on it would strand the very run that fixes it."""
        (script.OUT / "split_manifest.csv").write_text("")

        script._refuse_to_shrink_the_manifest(panel_only=True)

    def test_a_manifest_with_no_dataset_column_is_not_treated_as_a_loss_with_panel_only(self, script):
        """Same reasoning for a file this script did not write. It cannot be read as
        carrying real rows, so it cannot be read as something being lost."""
        pd.DataFrame({"something_else": [1]}).to_csv(
            script.OUT / "split_manifest.csv", index=False)

        script._refuse_to_shrink_the_manifest(panel_only=True)


class TestTheGuardIsWiredIntoTheRunAndNotJustDefined:
    """A guard nothing calls is a guard that does not exist, and this one is reached only
    on the branch where ``paths.biosensor_plates()`` returns ``None``."""

    def test_a_run_without_plates_exits_non_zero_and_writes_nothing(self, script,
                                                                    monkeypatch):
        _manifest(["simulated_panel", "real_biosensor"]).to_csv(
            script.OUT / "split_manifest.csv", index=False)
        monkeypatch.setattr(script.paths, "biosensor_plates", lambda: None)
        monkeypatch.setattr(script.sys, "argv", ["make_splits.py"])

        with pytest.raises(SystemExit) as caught:
            script.main()

        assert caught.value.code != 0
        assert list(pd.read_csv(script.OUT / "split_manifest.csv").dataset) == [
            "simulated_panel", "real_biosensor"]


def test_fresh_full_run_without_real_input_refuses(tmp_path):
    output = tmp_path / "fresh-output"
    result = subprocess.run(
        [sys.executable, "-B", str(_SCRIPT)], cwd=_REPO,
        env={**os.environ, "YSTWIN_OUTPUTS": str(output),
             "YSTWIN_PLATES": str(tmp_path / "missing-plates")},
        capture_output=True, text=True, check=False,
    )

    assert result.returncode != 0, result.stdout
    assert "--real-table" in result.stderr
    assert "--panel-only" in result.stderr
    assert not (output / "split_manifest.csv").exists()
    assert not (output / "split_manifest_summary.csv").exists()


@pytest.fixture
def public_readings():
    return pd.read_csv(_REAL_TABLE, dtype={"plate": str})


def test_full_run_uses_explicit_public_table(script, monkeypatch, public_readings):
    from ystwin.plate.layout import plate_key

    monkeypatch.setattr(script.paths, "biosensor_plates",
                        lambda: pytest.fail("an explicit table must not resolve private plates"))
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT), "--real-table", str(_REAL_TABLE)])
    script.main()

    manifest = pd.read_csv(script.OUT / "split_manifest.csv", dtype={"plate": str})
    summary = pd.read_csv(script.OUT / "split_manifest_summary.csv")
    assert set(manifest.dataset) == set(summary.dataset) == {script.REAL, script.PANEL}
    assert list(manifest.columns) == [
        "dataset", "split_kind", "seed", "group_key", "within", "axis", "stratum",
        "group_id", "construct", "stressor", "plate", "dose_mM", "axis_value",
        "assignment", "n_rows", "held_out", "partition_hash",
    ]
    expected = public_readings.assign(plate=public_readings.plate.map(plate_key))
    sizes = expected.groupby(_GROUP_KEY).size()
    real = manifest[manifest.dataset == script.REAL]
    kinds = {"interpolation", "extrapolation", "heldout_construct", "heldout_replicate"}
    assert set(real.split_kind) == kinds
    for kind in kinds:
        for seed in script.SEEDS:
            block = real[(real.split_kind == kind) & (real.seed == seed)]
            assert set(block[_GROUP_KEY].itertuples(index=False, name=None)) == set(sizes.index)
            for row in block.itertuples():
                assert row.n_rows == sizes.loc[(row.plate, row.construct, row.dose_mM)]
            partition = script.make_split(expected, kind, seed=seed, dataset=script.REAL,
                                          group_key=_GROUP_KEY)
            assert set(block.partition_hash) == {partition.hash}
            assert dict(zip(block.group_id, block.assignment)) == partition.assignment
    print(f"full CLI: assignments={len(manifest)}, summaries={len(summary)}, "
          f"populations={manifest.groupby('dataset').size().to_dict()}, "
          f"schema={list(manifest.columns)}")


def test_panel_only_never_resolves_real_inputs(script, monkeypatch):
    monkeypatch.setattr(script.paths, "biosensor_plates",
                        lambda: pytest.fail("--panel-only must not resolve private plates"))
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT), "--panel-only"])
    script.main()

    manifest = pd.read_csv(script.OUT / "split_manifest.csv")
    summary = pd.read_csv(script.OUT / "split_manifest_summary.csv")
    assert set(manifest.dataset) == set(summary.dataset) == {script.PANEL}


@pytest.mark.parametrize("contents", [
    None, "", "something_else\n1\n", "dataset\nsimulated_panel\n", "dataset\nreal_biosensor\n",
])
def test_full_scope_requires_real_input_for_every_destination_state(script, monkeypatch, contents):
    manifest = script.OUT / "split_manifest.csv"
    summary = script.OUT / "split_manifest_summary.csv"
    if contents is not None:
        manifest.write_text(contents)
    summary.write_text("existing summary\n")
    monkeypatch.setattr(script.paths, "biosensor_plates", lambda: None)
    monkeypatch.setattr(script, "panel_inventory", lambda: pytest.fail("validate real inputs first"))
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT)])

    with pytest.raises(SystemExit, match="full scope requires real_biosensor"):
        script.main()

    assert (manifest.read_text() if manifest.exists() else None) == contents
    assert summary.read_text() == "existing summary\n"


def test_missing_explicit_real_table_never_falls_back_or_overwrites(script, monkeypatch):
    for name in ("split_manifest.csv", "split_manifest_summary.csv"):
        (script.OUT / name).write_text(f"retained {name}\n")
    monkeypatch.setattr(script.paths, "biosensor_plates",
                        lambda: pytest.fail("a missing explicit input must not fall back"))
    monkeypatch.setattr(script, "panel_inventory", lambda: pytest.fail("validate real inputs first"))
    missing = script.OUT / "missing-real.csv"
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT), "--real-table", str(missing)])

    with pytest.raises(SystemExit, match="required --real-table input"):
        script.main()

    for name in ("split_manifest.csv", "split_manifest_summary.csv"):
        assert (script.OUT / name).read_text() == f"retained {name}\n"


def test_real_table_and_panel_only_are_conflicting_scopes(script, monkeypatch):
    monkeypatch.setattr(script, "real_table_inventory", lambda _: pytest.fail("resolve scope first"))
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT), "--panel-only", "--real-table", str(_REAL_TABLE)])

    with pytest.raises(SystemExit) as caught:
        script.main()

    assert caught.value.code == 2
    assert not (script.OUT / "split_manifest.csv").exists()
    assert not (script.OUT / "split_manifest_summary.csv").exists()


def test_public_inventory_preserves_each_well_and_uses_canonical_plate_ids(script, public_readings):
    from ystwin.plate.layout import plate_key

    inventory = script.real_table_inventory(_REAL_TABLE)
    expected = public_readings[_REAL_KEY].assign(plate=public_readings.plate.map(plate_key))
    pd.testing.assert_frame_equal(inventory[_REAL_KEY], expected)
    expected_plates = {plate_key(export.name) for export in
                       script.replay.exports(source_set=script.REAL_SOURCE_SET)}
    assert set(inventory.plate) == expected_plates
    assert set(inventory.loc[inventory.well == "F9", "plate"]) == expected_plates
    assert inventory.plate_has_reporter_blank.all()


def test_public_identity_validation_does_not_select_a_raw_numeric_channel(script, monkeypatch, public_readings):
    from ystwin.plate.layout import plate_key

    monkeypatch.setattr(script.replay.SynergyRun, "raw_channel",
                        lambda *args: pytest.fail("well identity must not select or correct raw values"))
    inventory = script.real_table_inventory(_REAL_TABLE)
    expected = public_readings[_REAL_KEY].assign(plate=public_readings.plate.map(plate_key))
    pd.testing.assert_frame_equal(inventory[_REAL_KEY], expected)


@pytest.mark.parametrize("derived", [False, True])
def test_reporter_well_coverage_uses_declared_nonderived_blocks(script, monkeypatch, derived):
    original = script.replay.load_run
    first_export = script.replay.exports(source_set=script.REAL_SOURCE_SET)[0].name

    def changed_coverage(export, **kwargs):
        run = original(export, **kwargs)
        if export != first_export:
            return run
        blocks = list(run.blocks)
        reporter = next(i for i, block in enumerate(blocks) if block.fluorophore == "mCitrine")
        block = blocks[reporter]
        blocks[reporter] = replace(block, data=block.data.drop(columns="F9"), derived=derived)
        return replace(run, blocks=tuple(blocks))

    monkeypatch.setattr(script.replay, "load_run", changed_coverage)
    if derived:
        inventory = script.real_table_inventory(_REAL_TABLE)
        assert "F9" in set(inventory.well)
    else:
        with pytest.raises(SystemExit, match="ambiguous reporter well coverage"):
            script.real_table_inventory(_REAL_TABLE)


def test_measurement_refusals_do_not_remove_well_identities(script, public_readings):
    public_readings["activity_late"] = float("nan")
    source = script.OUT / "unavailable-estimates.csv"
    public_readings.to_csv(source, index=False)

    pd.testing.assert_frame_equal(script.real_table_inventory(source),
                                  script.real_table_inventory(_REAL_TABLE))


@pytest.mark.parametrize("missing_column", _REAL_KEY)
def test_real_table_requires_every_identity_column(script, public_readings, missing_column):
    source = script.OUT / "missing-column.csv"
    public_readings.drop(columns=missing_column).to_csv(source, index=False)

    with pytest.raises(SystemExit, match=f"missing required columns.*{missing_column}"):
        script.real_table_inventory(source)


@pytest.mark.parametrize("column,value,reason", [
    ("plate", None, "missing real well identity"),
    ("plate", "unrecorded-plate", "unrecorded plate identities"),
    ("construct", " ", "empty construct identities"),
    ("well", None, "missing real well identity"),
    ("dose_mM", "not-a-dose", "numeric dose_mM"),
    ("dose_mM", float("nan"), "missing real well identity"),
    ("dose_mM", float("inf"), "finite, nonnegative dose_mM"),
    ("dose_mM", -1.0, "finite, nonnegative dose_mM"),
])
def test_invalid_real_identities_are_refused(script, public_readings, column, value, reason):
    public_readings[column] = public_readings[column].astype(object)
    public_readings.loc[0, column] = value
    source = script.OUT / "invalid-identity.csv"
    public_readings.to_csv(source, index=False)

    with pytest.raises(SystemExit, match=reason):
        script.real_table_inventory(source)


@pytest.mark.parametrize("contents", ["", "plate,construct,stressor,dose_mM,well\n"])
def test_empty_real_table_is_not_a_panel_only_request(script, contents):
    source = script.OUT / "empty-real.csv"
    source.write_text(contents)

    with pytest.raises(SystemExit, match="invalid --real-table|no real well identities"):
        script.real_table_inventory(source)


@pytest.mark.parametrize("population", ["plate", "construct", "F9"])
def test_missing_real_subpopulation_refuses_before_writing(script, monkeypatch, public_readings, population):
    if population == "plate":
        keep = public_readings.plate != public_readings.plate.iloc[0]
    elif population == "construct":
        keep = public_readings.construct != "NativeYap1"
    else:
        keep = public_readings.well != "F9"
    source = script.OUT / "incomplete-real.csv"
    public_readings[keep].to_csv(source, index=False)
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT), "--real-table", str(source)])
    monkeypatch.setattr(script, "panel_inventory", lambda: pytest.fail("validate real inputs first"))

    with pytest.raises(SystemExit, match="missing well identities"):
        script.main()

    assert not (script.OUT / "split_manifest.csv").exists()
    assert not (script.OUT / "split_manifest_summary.csv").exists()


@pytest.mark.parametrize("column,value", [
    ("well", "H9"), ("stressor", "DTT"), ("construct", "UPRE1"), ("dose_mM", 0.25),
    ("plate", "20260821"),
])
def test_same_size_wrong_real_population_is_not_accepted(script, public_readings, column, value):
    public_readings.loc[0, column] = value
    source = script.OUT / "substituted-identity.csv"
    public_readings.to_csv(source, index=False)

    with pytest.raises(SystemExit, match="unexpected well identities"):
        script.real_table_inventory(source)


def test_two_export_aliases_of_one_well_are_not_two_replicates(script, public_readings):
    from ystwin.plate.layout import plate_key

    duplicate = public_readings.iloc[[0]].copy()
    duplicate["plate"] = duplicate.plate.map(plate_key)
    source = script.OUT / "duplicate-well.csv"
    pd.concat([public_readings, duplicate], ignore_index=True).to_csv(source, index=False)

    with pytest.raises(SystemExit, match="repeats physical plate/well identities"):
        script.real_table_inventory(source)


def test_current_source_set_selects_exports_not_dates_or_counts(script, monkeypatch):
    manifest = script.replay.load_manifest()
    expected = list(script.replay.exports(manifest=manifest, source_set=script.REAL_SOURCE_SET))
    visited = []
    original = script.replay.load_run

    def load_run(export, **kwargs):
        visited.append(pathlib.Path(export))
        return original(export, **kwargs)

    monkeypatch.setattr(script.replay, "load_run", load_run)
    script.real_table_inventory(_REAL_TABLE)

    assert visited == expected
    assert set(visited).isdisjoint(script.replay.exports(manifest=manifest, source_set="newprotocol_superseded"))
    assert set(visited).isdisjoint(script.replay.exports(manifest=manifest, source_set="july_2026_07"))


@pytest.mark.parametrize("defect,reason", [
    ("no_source_set", "missing columns.*source_set"),
    ("no_current_population", "no source_set='newprotocol' exports"),
    ("unrecorded_plate", "no recorded plate identity"),
    ("duplicate_current_plate", "multiple current exports declare plate"),
])
def test_missing_or_ambiguous_manifest_scope_refuses(script, monkeypatch, defect, reason):
    manifest = script.replay.load_manifest()
    if defect == "no_source_set":
        manifest = manifest.drop(columns="source_set")
    elif defect == "no_current_population":
        manifest = manifest[manifest.source_set != script.REAL_SOURCE_SET]
    elif defect == "unrecorded_plate":
        manifest.loc[manifest.export == manifest.export.iloc[0], "export"] = "unrecorded.xlsx"
    else:
        duplicate = manifest[manifest.export == manifest.export.iloc[0]].copy()
        duplicate["export"] = duplicate.export.map(lambda name: f"{pathlib.Path(name).stem}_copy.xlsx")
        manifest = pd.concat([manifest, duplicate], ignore_index=True)
    monkeypatch.setattr(script.replay, "load_manifest", lambda: manifest)

    with pytest.raises(SystemExit, match=reason):
        script.real_table_inventory(_REAL_TABLE)


def test_missing_public_scope_manifest_refuses(script, monkeypatch):
    monkeypatch.setattr(script.replay, "committed_dir", lambda: script.OUT / "no-public-plates")

    with pytest.raises(SystemExit, match="required real-population manifest is unavailable"):
        script.real_table_inventory(_REAL_TABLE)


def test_missing_required_public_readings_refuse(script, monkeypatch):
    def unavailable(*args, **kwargs):
        raise FileNotFoundError("required-reporter.csv")

    monkeypatch.setattr(script.replay, "load_run", unavailable)

    with pytest.raises(SystemExit, match="required real-population readings.*required-reporter.csv"):
        script.real_table_inventory(_REAL_TABLE)


def test_partial_raw_inventory_also_refuses_full_scope(script, monkeypatch):
    inventory = script.real_table_inventory(_REAL_TABLE)
    incomplete = inventory[inventory.construct != "NativeYap1"]
    monkeypatch.setattr(script.paths, "biosensor_plates", lambda: script.OUT / "selected-raw-input")
    monkeypatch.setattr(script, "real_inventory", lambda _: incomplete)
    monkeypatch.setattr(script, "panel_inventory", lambda: pytest.fail("validate real inputs first"))
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT)])

    with pytest.raises(SystemExit, match="missing well identities"):
        script.main()

    assert not (script.OUT / "split_manifest.csv").exists()
    assert not (script.OUT / "split_manifest_summary.csv").exists()


def test_refused_real_splits_cannot_be_reported_as_full_scope_success(script, monkeypatch):
    original = script.attempt

    def only_panel(splits, refusals, frame, kind, dataset, *args, **kwargs):
        if dataset == script.REAL:
            return None
        return original(splits, refusals, frame, kind, dataset, *args, **kwargs)

    monkeypatch.setattr(script, "attempt", only_panel)
    monkeypatch.setattr(script.sys, "argv", [str(_SCRIPT), "--real-table", str(_REAL_TABLE)])

    with pytest.raises(SystemExit, match="no supported splits for requested populations.*real_biosensor"):
        script.main()

    assert not (script.OUT / "split_manifest.csv").exists()
    assert not (script.OUT / "split_manifest_summary.csv").exists()
