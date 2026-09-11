from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.analysis.hog_data import HogDataset, load_hog_data


S4_LOW = "wild type pfk26Δ27Δ 0.4M"
S4_HIGH = "wild type pfk26Δ27Δ 0.8M"
REQUIRED_COLUMNS = {
    "source_doi", "supplement_id", "source_sheet", "source_cell", "source_time_cell",
    "experiment_id", "measurement_block", "replicate_label", "strain_id", "genotype",
    "background", "medium", "nacl_molar", "time_native", "time_unit", "time_unit_evidence",
    "time_relative_s", "time_model_s", "observable_id", "value", "raw_value", "unit",
    "raw_unit", "sd", "measurement_provenance", "normalization", "original_fit_status",
    "split",
}


@pytest.fixture(scope="module")
def raw_directory():
    pytest.importorskip("xlrd")
    directory = paths.data_dir() / "hog2013"
    for name in ("observations.xls", "western_blots.xls", "methods.pdf", "sources.json"):
        assert (directory / name).is_file(), name
    return directory


@pytest.fixture(scope="module")
def dataset(raw_directory):
    return load_hog_data(raw_directory, preliminary_time_unit="minutes")


@pytest.fixture
def raw_copy(raw_directory, tmp_path):
    for name in ("observations.xls", "western_blots.xls", "methods.pdf", "sources.json"):
        (tmp_path / name).write_bytes((raw_directory / name).read_bytes())
    return tmp_path


def observation(dataset, sheet, cell):
    rows = dataset.observations.loc[
        (dataset.observations.source_sheet == sheet)
        & (dataset.observations.source_cell == cell)
    ]
    assert len(rows) == 1
    return rows.iloc[0]


def replace_cell(monkeypatch, sheet, column, row, value):
    original = pd.ExcelFile.parse

    def parse(self, sheet_name=0, *args, **kwargs):
        frame = original(self, sheet_name, *args, **kwargs)
        if sheet_name == sheet:
            frame = frame.copy()
            frame.iat[row - 1, ord(column) - ord("A")] = value
        return frame

    monkeypatch.setattr(pd.ExcelFile, "parse", parse)


def test_requires_explicit_preliminary_time_unit():
    with pytest.raises(TypeError, match="preliminary_time_unit"):
        load_hog_data()
    with pytest.raises(TypeError):
        load_hog_data(None, "minutes")


@pytest.mark.parametrize("unit", [None, "", "min", "seconds", "Minutes"])
def test_rejects_unapproved_time_unit_before_reading_files(tmp_path, unit):
    with pytest.raises(ValueError, match="preliminary_time_unit.*minutes"):
        load_hog_data(tmp_path, preliminary_time_unit=unit)


def test_default_directory_uses_paths_and_custom_directory_accepts_string(
    raw_directory, dataset, monkeypatch
):
    monkeypatch.setattr(paths, "data_dir", lambda: raw_directory.parent)
    default = load_hog_data(preliminary_time_unit="minutes")
    explicit = load_hog_data(str(raw_directory), preliminary_time_unit="minutes")
    pd.testing.assert_frame_equal(default.observations, dataset.observations)
    pd.testing.assert_frame_equal(explicit.observations, dataset.observations)


def test_tidy_contract_and_frozen_counts(dataset):
    assert isinstance(dataset, HogDataset)
    assert REQUIRED_COLUMNS <= set(dataset.observations.columns)
    assert len(dataset.observations) == 407
    assert str(dataset.observations.strain_id.dtype) == "Int64"
    assert dataset.observations.source_doi.eq("10.1371/journal.pcbi.1003084").all()
    assert dataset.observations.groupby("supplement_id").size().to_dict() == {
        "s002": 317, "s004": 90,
    }
    assert {key: len(frame) for key, frame in dataset.split().items()} == {
        "fit": 13, "heldout_dose": 38, "heldout_assay": 52, "unused": 304,
    }


def test_known_reported_molar_values_and_source_coordinates(dataset):
    hog = observation(dataset, "Hog1PP_ALL", "C28")
    assert hog.value == 3.6344917648231455e-7
    assert hog.raw_value == hog.value
    assert hog.sd == 5.6649698724306056e-8
    assert hog.source_time_cell == "B28"
    assert hog.sd_source_cell == "D28"
    assert hog.observable_id == "Hog1PP_measured"
    assert hog.unit == "mol/L"
    assert hog.raw_unit == "mol/l"
    assert hog.strain_id == 202
    assert hog.genotype == "wild_type"
    assert hog.background == "W303"
    assert hog.medium == "YPD"
    assert hog.nacl_molar == 0.4
    gpd = observation(dataset, "Gpd1_ALL", "C21")
    assert gpd.value == 7.574314387968357e-8
    assert gpd.observable_id == "Gpd1_measured"
    glycerol = observation(dataset, "glycerol_i", "N7")
    assert glycerol.value == 0.08781879368758393
    assert glycerol.observable_id == "glycerol_measured"
    assert glycerol.source_time_cell == "M7"
    assert np.isnan(glycerol.sd)
    assert glycerol.time_model_s == 3900
    assert glycerol.time_relative_s == 300
    external = observation(dataset, "glycerol_e", "N7")
    assert external.value == 0.0026527961776522968
    assert external.observable_id == "glycerol_e"


def test_time_bases_are_explicit_and_not_silently_equated(dataset):
    s2 = dataset.observations.query("supplement_id == 's002'")
    assert s2.time_unit.eq("seconds").all()
    assert s2.time_unit_evidence.eq("explicit_workbook_header").all()
    np.testing.assert_array_equal(s2.time_model_s, s2.time_native)
    np.testing.assert_array_equal(s2.time_relative_s, s2.time_native - 3600)
    assert observation(dataset, "Hog1PP_ALL", "C26").time_relative_s == -900
    s4 = dataset.observations.query("supplement_id == 's004'")
    assert s4.time_unit.eq("minutes").all()
    assert s4.time_unit_evidence.eq("caller_supplied_inference_unlabelled_workbook").all()
    np.testing.assert_array_equal(s4.time_relative_s, s4.time_native * 60)
    np.testing.assert_array_equal(s4.time_model_s, s4.time_native * 60 + 3600)
    point = observation(dataset, S4_HIGH, "E25")
    assert point.time_native == 2
    assert point.time_relative_s == 120
    assert point.time_model_s == 3720
    assert point.value == 0.5358560927876885
    evidence = dataset.metadata["preliminary_time_unit"]
    assert evidence["unit"] == "minutes"
    assert evidence["explicit_workbook_label"] is False
    assert evidence["basis"]


def test_percent_conversion_is_only_applied_to_the_two_identified_blocks(dataset):
    frame = dataset.observations
    percent = (
        ((frame.source_sheet == S4_LOW) & (frame.measurement_block == "Measurement3"))
        | ((frame.source_sheet == S4_HIGH) & (frame.measurement_block == "Measurement4"))
    )
    assert percent.sum() == 26
    np.testing.assert_array_equal(frame.loc[percent, "value"], frame.loc[percent, "raw_value"] / 100)
    np.testing.assert_array_equal(frame.loc[~percent, "value"], frame.loc[~percent, "raw_value"])
    assert frame.loc[percent, "raw_unit"].eq("percent_of_peak").all()
    assert frame.loc[frame.supplement_id == "s004", "unit"].eq("relative_intensity").all()
    assert observation(dataset, S4_LOW, "D46").value == 1
    assert observation(dataset, S4_HIGH, "E85").value == 1


def test_missing_measurements_are_not_zeroes_or_interpolated_points(dataset):
    frame = dataset.observations
    assert frame.value.notna().all()
    assert frame.query("genotype == 'gpd1_del' and observable_id == 'glycerol_measured'").empty
    assert not ((frame.source_sheet == S4_HIGH) & (frame.source_cell == "E36")).any()
    assert observation(dataset, S4_HIGH, "E35").value == 0
    assert frame.loc[frame.supplement_id == "s004", "sd"].isna().all()
    assert frame.loc[frame.source_sheet.isin(["glycerol_i", "glycerol_e"]), "sd"].isna().all()
    assert frame.groupby("source_sheet").size().to_dict() == {
        "Hog1PP_ALL": 65, "Gpd1_ALL": 65, "glycerol_i": 84, "glycerol_e": 103,
        S4_LOW: 52, S4_HIGH: 38,
    }


def test_does_not_import_duplicate_or_inferred_targets(dataset):
    frame = dataset.observations
    assert not frame.duplicated(["supplement_id", "source_sheet", "source_cell"]).any()
    assert set(frame.observable_id) == {
        "Hog1PP_measured", "Gpd1_measured", "glycerol_measured", "glycerol_e",
    }
    for sheet, first, last in [("Hog1PP_ALL", 26, 38), ("Gpd1_ALL", 19, 31)]:
        rows = frame.loc[frame.source_sheet == sheet, "source_cell"].str.extract(r"(\d+)$")[0].astype(int)
        assert rows.between(first, last).all()
    assert not frame.source_sheet.isin(["Numbers", "Numbers 2", "wild type gpd2Δ gpd1Δ"]).any()
    assert not (
        (frame.source_sheet == S4_HIGH) & (frame.measurement_block == "Measurement3")
    ).any()


def test_frozen_split_holds_whole_groups_and_only_fits_wt_hog1(dataset):
    partitions = dataset.split()
    fit = partitions["fit"]
    assert fit.supplement_id.eq("s002").all()
    assert fit.observable_id.eq("Hog1PP_measured").all()
    for name, dose in [("fit", 0.4), ("heldout_dose", 0.8), ("heldout_assay", 0.4)]:
        frame = partitions[name]
        assert frame.strain_id.eq(202).all()
        assert frame.genotype.eq("wild_type").all()
        assert frame.background.eq("W303").all()
        assert frame.medium.eq("YPD").all()
        assert frame.nacl_molar.eq(dose).all()
    groups = {key: set(frame.experiment_id) for key, frame in partitions.items()}
    indices = {key: set(frame.index) for key, frame in partitions.items()}
    for left in partitions:
        for right in partitions:
            if left != right:
                assert groups[left].isdisjoint(groups[right])
                assert indices[left].isdisjoint(indices[right])
    assert set.union(*indices.values()) == set(dataset.observations.index)
    assert dataset.observations.groupby("experiment_id").split.nunique().eq(1).all()
    assert partitions["heldout_dose"].experiment_id.nunique() == 3
    assert partitions["heldout_assay"].experiment_id.nunique() == 3
    low_first = partitions["heldout_assay"].query("measurement_block == 'Measurement1'")
    assert len(low_first) == 26
    assert low_first.experiment_id.nunique() == 1
    assert set(low_first.replicate_label) == {"(1)", "(2)"}


def test_split_returns_copies(dataset):
    fit = dataset.split()["fit"]
    index = fit.index[0]
    original = dataset.observations.loc[index, "value"]
    fit.loc[index, "value"] = -999
    assert dataset.observations.loc[index, "value"] == original


def test_unspecified_mutant_backgrounds_and_control_dose_are_not_invented(dataset):
    frame = dataset.observations
    assert frame.loc[frame.genotype != "wild_type", "background"].isna().all()
    controls = frame.query("source_sheet == 'glycerol_i' and nacl_molar == 0")
    assert len(controls) == 14
    assert controls.genotype.eq("wild_type").all()
    assert set(frame.genotype) == {
        "wild_type", "pfk2627_del", "gpd1_del", "hog1_del", "fps1_open", "hog1_att", "stl1_del",
    }
    assert frame.loc[frame.genotype == "hog1_att", "strain_id"].isna().all()


def test_metadata_distinguishes_calibration_provenance_and_limited_coverage(dataset):
    metadata = dataset.metadata
    assert metadata["protein_units"]["direct_absolute_concentration"] is False
    assert metadata["protein_units"]["hog1_peak_fraction_assumed"] == 0.9
    assert metadata["protein_units"]["model_unit_conversion_applied"] is False
    assert metadata["published_fit"]["s002"]["includes_wild_type_and_mutants"] is True
    assert metadata["published_fit"]["s004"]["used_for_parameter_fitting"] is False
    assert metadata["published_fit"]["s004"]["part_of_model_development"] is True
    assert metadata["replication"]["biological_independence_established"] is False
    assert metadata["split_policy"]["heldout_fitting_allowed"] is False
    assert {"gpd1mRNA_measured", "relVM"} <= set(metadata["missing_observables"])
    exclusions = metadata["not_imported"]
    assert any(item["source_sheet"] == "wild type gpd2Δ gpd1Δ" for item in exclusions)
    assert any(item["source_sheet"] == S4_LOW for item in exclusions)
    assert any(item["source_sheet"] == S4_HIGH for item in exclusions)
    assert all(item["reason"] for item in exclusions)
    assert metadata["sources"]["observations"]["sha256"] == (
        "87854ff96bf22943f1347abdf5636e261ddfcdec041ea78fe498d62d924a967e"
    )
    json.dumps(metadata, allow_nan=False)


@pytest.mark.parametrize("filename", ["observations.xls", "western_blots.xls", "methods.pdf"])
def test_rejects_changed_file_hashes(raw_copy, filename):
    target = raw_copy / filename
    target.write_bytes(target.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="SHA-256.*" + re.escape(filename)):
        load_hog_data(raw_copy, preliminary_time_unit="minutes")


@pytest.mark.parametrize("field,value", [
    ("sha256", "0" * 64), ("filename", "different.xls"), ("url", "https://example.com/data"),
])
def test_manifest_cannot_relabel_the_verified_source(raw_copy, field, value):
    source = raw_copy / "sources.json"
    manifest = json.loads(source.read_text())
    manifest["assets"]["observations"][field] = value
    source.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="manifest.*observations"):
        load_hog_data(raw_copy, preliminary_time_unit="minutes")


def test_rejects_wrong_publication_identity(raw_copy):
    source = raw_copy / "sources.json"
    manifest = json.loads(source.read_text())
    manifest["publication"]["doi"] = "different-doi"
    source.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="publication"):
        load_hog_data(raw_copy, preliminary_time_unit="minutes")


@pytest.mark.parametrize("sheet,column,row,value", [
    ("Hog1PP_ALL", "C", 21, "norm intensity"),
    ("Hog1PP_ALL", "C", 23, 444),
    ("Hog1PP_ALL", "C", 24, "stdev"),
    ("Hog1PP_ALL", "B", 25, "time (min)"),
    ("Gpd1_ALL", "C", 17, "norm int"),
    ("Gpd1_ALL", "C", 18, "avg 444"),
    ("glycerol_i", "N", 1, "Wild Type1"),
    ("glycerol_i", "M", 2, "time (min)"),
    (S4_LOW, "A", 4, "YNB"),
    (S4_LOW, "A", 7, "0.8M NaCl"),
    (S4_LOW, "A", 12, "wt 689"),
    (S4_LOW, "E", 1, "double mutant (1)"),
    (S4_LOW, "F", 1, "wild type (3)"),
    (S4_LOW, "D", 20, "pfk26Δ"),
    (S4_LOW, "D", 41, "pfk26Δ27Δ"),
    (S4_HIGH, "E", 23, "pfk26d"),
    (S4_HIGH, "E", 42, 689),
    (S4_HIGH, "E", 82, "pfk26d27d"),
])
def test_rejects_layout_or_identity_changes_even_with_verified_input_bytes(
    raw_directory, monkeypatch, sheet, column, row, value
):
    replace_cell(monkeypatch, sheet, column, row, value)
    with pytest.raises(ValueError, match=re.escape(f"{sheet}!{column}{row}")):
        load_hog_data(raw_directory, preliminary_time_unit="minutes")


@pytest.mark.parametrize("sheet,column,row,value,expected", [
    (S4_LOW, "E", 2, -0.125, -0.125),
    (S4_HIGH, "E", 83, -2.0, -0.02),
])
def test_negative_values_survive_without_clipping_or_adaptive_scaling(
    raw_directory, monkeypatch, sheet, column, row, value, expected
):
    replace_cell(monkeypatch, sheet, column, row, value)
    result = load_hog_data(raw_directory, preliminary_time_unit="minutes")
    point = observation(result, sheet, f"{column}{row}")
    assert point.raw_value == value
    assert point.value == expected


@pytest.mark.parametrize("value", ["not a measurement", np.inf])
def test_invalid_numeric_observations_are_rejected(raw_directory, monkeypatch, value):
    replace_cell(monkeypatch, S4_HIGH, "E", 25, value)
    with pytest.raises(ValueError, match=re.escape(f"{S4_HIGH}!E25")):
        load_hog_data(raw_directory, preliminary_time_unit="minutes")
