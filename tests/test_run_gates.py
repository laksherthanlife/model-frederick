from __future__ import annotations

from dataclasses import replace
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ystwin.plate import replay
from ystwin.plate.layout import PlateLayout, RecordedPlate
from ystwin.plate.synergy import KineticBlock, SynergyRun


@pytest.fixture
def script(monkeypatch, tmp_path):
    monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
    filename = Path(__file__).resolve().parents[1] / "scripts/run_gates.py"
    spec = importlib.util.spec_from_file_location("run_gates_test", filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def recorded_run(*, failed_culture=True):
    t = np.linspace(0, 4, 25)
    growing = 0.1 * np.exp(0.3 * t)
    second = np.full_like(t, 0.005) if failed_culture else growing * 0.9
    od = pd.DataFrame({"A1": growing + 0.1, "A2": second + 0.1,
                       "H1": np.full_like(t, 0.1), "H2": np.full_like(t, 0.04)}, index=t)
    rfu = pd.DataFrame({"A1": growing * 5000 + 100, "A2": second * 5000 + 100,
                        "H1": np.full_like(t, 100), "H2": np.full_like(t, 20)}, index=t)
    recorded = RecordedPlate(PlateLayout({"UPRE1": (1, 2)}, ("A",)), ("H1",), {"DTT": (0.0,)})
    return SynergyRun(Path("opaque-name.xlsx"), (
        KineticBlock("density", "OD600", "600", "first", False, False, od, []),
        KineticBlock("reporter", "mCitrine", "480,530", "second", False, False, rfu, []),
    ), recorded_plate=recorded)


def test_recorded_failed_cultures_stay_but_unrecorded_wells_do_not(script, monkeypatch, tmp_path):
    run = recorded_run()
    monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: run)

    result = script.run(run.source)
    table = pd.read_csv(tmp_path / "g1_opaque-name.csv")

    assert result["n_cultures"] == 2
    assert set(table.well) == {"A1", "A2"}
    assert not table.set_index("well").loc["A2", "passed"]
    assert "above_blank" in table.set_index("well").loc["A2", "failures"]
    assert result["unrecorded_wells"] == ["H2"]
    assert result["blank_wells"] == ["H1"]


def test_current_xpt_has_84_recorded_cultures_not_93(script, monkeypatch, tmp_path):
    name = "20260804_ER_Oxidative_Replicate4.xpt"
    run = replay.load_run(name)
    monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: run)

    result = script.run(Path(name))
    table = pd.read_csv(tmp_path / "g1_20260804_ER_Oxidative_Replicate4.csv")

    assert result["n_cultures"] == 84
    assert result["n_pass"] == 55
    assert set(table.well) == set(run.recorded_plate.culture_wells)
    assert int((~table.passed).sum()) == 29
    assert set(result["unrecorded_wells"]) == {f"H{i}" for i in range(4, 13)}


def test_missing_recorded_culture_is_an_explicit_refusal_not_a_smaller_cohort(script, monkeypatch, tmp_path):
    run = recorded_run()
    blocks = tuple(replace(block, data=block.data.drop(columns="A2")) for block in run.blocks)
    run = replace(run, blocks=blocks)
    monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: run)

    result = script.run(run.source)

    assert "A2" in result["note"]
    assert result["n_cultures"] == 2
    assert result["missing_culture_wells"] == ["A2"]
    assert not list(tmp_path.glob("g1_*.csv"))


def test_declared_corrected_block_is_not_subtracted_again(script, monkeypatch, tmp_path, capsys):
    run = recorded_run()
    raw = run.blocks[1]
    corrected = replace(raw, channel="opaque-processed", blank_subtracted=True,
                        data=raw.data + 250.0)
    run = replace(run, blocks=(run.blocks[0], corrected, raw))
    monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: run)

    result = script.run(run.source)

    assert result["n_cultures"] == 2
    assert "reporter=reporter" in capsys.readouterr().out


def test_ambiguous_public_raw_declarations_refuse_before_writing_gates(script, monkeypatch, tmp_path):
    name = "20260803_ER&oxidativestress_Replicate3.xlsx"
    manifest = replay.load_manifest()
    manifest.loc[manifest.export == name, "blank_subtracted"] = "False"
    run = replay.load_run(name, manifest=manifest)
    monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: run)

    result = script.run(Path(name))

    assert "ambiguous raw" in result["note"]
    assert not list(tmp_path.glob("g1_*.csv"))
    assert not list(tmp_path.glob("d2_g1passed_*.csv"))


def test_d2_refusals_are_retained_with_the_attempted_denominator(script, monkeypatch):
    run = recorded_run(failed_culture=False)
    monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: run)
    original = script.plate_dilution_report

    def exclude_one(*args, **kwargs):
        report = original(*args, **kwargs)
        return replace(report, per_well=report.per_well.iloc[:1].copy(),
                       excluded=pd.DataFrame([{"well": "A2", "reason": "unidentifiable test trace"}]))

    monkeypatch.setattr(script, "plate_dilution_report", exclude_one)
    result = script.run(run.source)

    assert result["n_pass"] == result["n_d2_attempted"] == 2
    assert result["n_d2_scored"] == 1
    assert result["d2_exclusions"] == [{"well": "A2", "reason": "unidentifiable test trace"}]


def test_manifest_source_sets_select_current_representations_only(script):
    manifest = replay.load_manifest()
    selected, excluded = script._select_exports(replay.exports(manifest=manifest), manifest)
    expected = set(manifest.loc[manifest.source_set.isin(["newprotocol", "july_2026_07"]), "export"])

    assert {path.name for path in selected} == expected
    assert len(selected) == 6
    assert "20260804_ER_Oxidative_Replicate4.xpt" in expected
    assert "20260804_ER&OxidativeStress_Replicate4.xlsx" not in expected
    assert {entry["file"]: entry["status"] for entry in excluded} == {
        "20260804_ER&OxidativeStress_Replicate4.xlsx": "superseded",
        "20260814_BY4741_0mM_OD600 Data.xlsx": "outside_scope",
    }


def test_source_selection_uses_classification_not_export_spelling(script):
    manifest = pd.DataFrame([
        {"export": "looks_current.xlsx", "source_set": "newprotocol_superseded", "source_sha256": "a" * 64},
        {"export": "looks_old.xlsx", "source_set": "newprotocol", "source_sha256": "b" * 64},
    ])
    selected, excluded = script._select_exports([Path(name) for name in manifest.export], manifest)

    assert selected == [Path("looks_old.xlsx")]
    assert excluded[0]["file"] == "looks_current.xlsx"
    assert excluded[0]["status"] == "superseded"


def test_duplicate_current_source_identity_is_not_two_replicates(script):
    manifest = pd.DataFrame([
        {"export": name, "source_set": "newprotocol", "source_sha256": "a" * 64}
        for name in ("one.xlsx", "another_name.xlsx")
    ])
    with pytest.raises(ValueError, match="same source"):
        script._select_exports([Path(name) for name in manifest.export], manifest)


def test_conflicting_manifest_classifications_are_refused(script):
    manifest = pd.DataFrame([
        {"export": "one.xlsx", "source_set": scope, "source_sha256": "a" * 64}
        for scope in ("newprotocol", "newprotocol_superseded")
    ])
    with pytest.raises(ValueError, match="conflicting"):
        script._select_exports([Path("one.xlsx")], manifest)


def test_main_records_refusals_and_supersession_and_returns_failure(script, monkeypatch, tmp_path):
    manifest = pd.DataFrame([
        {"export": name, "source_set": scope, "source_sha256": digest * 64}
        for name, scope, digest in (("ok.xlsx", "newprotocol", "a"),
                                    ("blocked.xlsx", "july_2026_07", "b"),
                                    ("old.xlsx", "newprotocol_superseded", "c"))
    ])
    monkeypatch.setattr(script.replay, "load_manifest", lambda: manifest)
    monkeypatch.setattr(script.replay, "install", lambda target: [Path(name) for name in manifest.export])
    monkeypatch.setattr(script.paths, "igem_results", lambda: pytest.fail("private resolver called"))
    monkeypatch.setattr(script.paths, "biosensor_plates", lambda: pytest.fail("private resolver called"))
    attempted = []
    reason = "unknown correction-state requiring an explicit processing declaration"

    def run(path):
        attempted.append(path.name)
        return ({"file": path.name, "note": reason} if path.name == "blocked.xlsx"
                else {"file": path.name, "n_cultures": 84, "n_pass": 55})

    monkeypatch.setattr(script, "run", run)
    assert script.main(["--committed"]) == 2
    assert attempted == ["ok.xlsx", "blocked.xlsx"]
    report = json.loads((tmp_path / "gates_manifest.json").read_text())
    entries = {entry["file"]: entry for entry in report["exports"]}
    assert entries["ok.xlsx"]["status"] == "processed"
    assert entries["blocked.xlsx"]["status"] == "refused"
    assert entries["blocked.xlsx"]["note"] == reason
    assert entries["old.xlsx"]["status"] == "superseded"
    assert report["complete"] is False
    assert report["n_sources_processed"] == 1
    assert report["n_sources_refused"] == 1


def test_all_current_public_gate_sources_resolve_without_reclassifying_or_double_subtracting(
        script, monkeypatch, tmp_path):
    monkeypatch.setattr(script.paths, "igem_results", lambda: pytest.fail("private resolver called"))
    monkeypatch.setattr(script.paths, "biosensor_plates", lambda: pytest.fail("private resolver called"))

    assert script.main(["--committed"]) == 0

    report = json.loads((tmp_path / "gates_manifest.json").read_text())
    assert report["complete"] is True
    assert report["n_sources_processed"] == 6
    assert report["n_sources_refused"] == 0
    entries = {entry["file"]: entry for entry in report["exports"]}
    august3 = entries["20260803_ER&oxidativestress_Replicate3.xlsx"]
    assert august3["n_cultures"] == 84
    assert august3["blank_wells"] == ["H1", "H2", "H3"]
    assert entries["20260804_ER&OxidativeStress_Replicate4.xlsx"]["status"] == "superseded"
    assert entries["20260804_ER_Oxidative_Replicate4.xpt"]["status"] == "processed"
