from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from ystwin.analysis.parameter_evidence import (
    PROGRAMMES,
    EvidenceGap,
    audit_local_sources,
    evidence_report,
    freeze_evidence,
    load_parameter_evidence,
    load_reference_model,
    model_parameter_inventory,
    native_asset_inventory,
    source_path,
    validate_inventory,
    verify_evidence_freeze,
)
from ystwin.mech.params import ParameterUncertainty, RefusedValue, Tag

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def inventory():
    return load_parameter_evidence()


def test_every_programme_has_checked_model_leads_and_explicit_blockers(inventory):
    payload = inventory.to_dict()
    assert set(payload["programmes"]) == set(PROGRAMMES)
    for name in PROGRAMMES:
        record = inventory.programme(name)
        assert record["scope"] and record["blockers"] and record["model_ids"]
        assert all(name in payload["models"][key]["programmes"] for key in record["model_ids"])
    assert payload["models"]["jalihal2021"]["parameters"] == 128
    assert payload["models"]["jalihal2021"]["model_id"] == ""
    assert payload["sources"]["jalihal2021"]["identity"]["commit"] == "3e93dba11d4e7950ff89519e40f32ba5766e7f16"
    assert payload["sources"]["christiano2014"]["verification"]["level"] == "indexed_identity"


@pytest.mark.parametrize("name,remote_path,byte_count,sha256,git_blob", [
    ("sbml-translator.py", "src/sbml-translator.py", 5197,
     "a54860f612ad9d66f3f01f0a7d98a1ce19ce5fbd35c1f336de9622301fbfacae", "9608f93ddf47dbd2361399c6099c111bf95099a2"),
    ("simulators.py", "nutrient_signaling/simulators.py", 11688,
     "f4bc27b2d79c7c79dada509082d4e6e7fc584edce1ab32e5d7203842a95f03b7", "714e1f7200280acb6ebb30b2c00e4883d4c61f8e"),
    ("timecourse.py", "nutrient_signaling/timecourse.py", 8085,
     "2255037810367c994f277c43c3869eb72c71796570217183518db3b1d56f28f9", "633f20db43d1114ff519a952ea6ecf292d9654e1"),
])
def test_jalihal_python_evidence_is_nonimportable_and_keeps_upstream_identity(
    inventory, name, remote_path, byte_count, sha256, git_blob,
):
    source = inventory.to_dict()["sources"]["jalihal2021"]
    assert source["license"]["status"] == "unverified"
    assert "spdx" not in source["license"]
    assets = {asset["path"]: asset for asset in source["local_artifacts"]}
    original = f"data/native_reference_models/jalihal2021/{name}"
    relative = f"{original}.source"
    assert assets[relative] == {
        "path": relative, "sha256": sha256, "bytes": byte_count,
        "acquisition": {"kind": "github", "path": remote_path, "git_blob": git_blob},
    }
    assert original not in assets and not (ROOT / original).exists()
    assert (ROOT / relative).is_file()
    assert importlib.util.spec_from_file_location("jalihal_upstream_source", ROOT / relative) is None


def test_hog_and_native_nutrient_sources_are_available_but_missing_programmes_fail_closed(inventory):
    assert len(inventory.require_reference("osmotic_hog")) == 6
    for programme in ("carbon_pka_snf1", "nitrogen_tor"):
        assert inventory.require_reference(programme) == ("jalihal2021",)
    for programme in set(PROGRAMMES) - {"osmotic_hog", "carbon_pka_snf1", "nitrogen_tor"}:
        with pytest.raises(EvidenceGap, match=programme):
            inventory.require_reference(programme)
    assert len(load_reference_model("jalihal2021", inventory).species_ids) == 25
    with pytest.raises(EvidenceGap, match="MATLAB supplement"):
        load_reference_model("zheng2016", inventory)
    with pytest.raises(ValueError, match="unknown"):
        inventory.require_reference("invented_programme")


def test_model_identities_include_corrections_and_all_bytes_are_checked(inventory):
    expected = {
        "hog2013_wt": "PetelenzKuehn_osmoadaptation_WT",
        "hog2013_hog1_del": "PetelenzKuehn_osmoadaptation_hog1D",
        "hog2013_hog1_att": "PetelenzKuehn_osmoadaptation_HOG1att",
        "hog2013_fps1_open": "PetelenzKuehn_osmoadaptation_fps1D1",
        "hog2013_gpd1_del": "PetelenzKuehn_osmoadaptation_gpd1D",
        "hog2013_pfk2627_del": "PetelenzKuehn_osmoadaptation_pfk2627D",
    }
    for key, identity in expected.items():
        model = load_reference_model(key, inventory)
        assert model.model_id == identity
        assert len(model.species_ids) == 29
        assert model.metadata["learned_parameters"] == {}
    for check in audit_local_sources(inventory):
        if check["source_id"] != "gasch2000_local" or check["actual_sha256"] is not None:
            assert check["verified"], check
    assert inventory.to_dict()["sources"]["hog2013"]["identity"]["correction_doi"] == "10.1371/journal.pcbi.1003663"


def test_all_global_sbml_parameters_and_assignments_are_inventoried_without_promoting_them(inventory):
    result = model_parameter_inventory("hog2013_wt", inventory)
    assert len(result["parameters"]) + len(result["assignments"]) == 126
    assert len(result["parameters"]) == 98
    assert all(p["status"] == "prior" and p["units"] == "unspecified" for p in result["parameters"])
    assert all(p["uncertainty"]["kind"] == "not_reported" for p in result["parameters"])
    assert all(p["unresolved_directions"] for p in result["parameters"])
    assert "gpd1mRNA_measured" in result["assignments"]
    assert len(result["reaction_formulas"]) == 58
    assert "Inline equation constants" in result["gap"]


def test_parameter_statuses_scales_and_refusals_preserve_existing_param_semantics(inventory):
    parameters = {p.parameter.name: p for p in inventory.parameters}
    assert {p.status for p in parameters.values()} == {"measured", "fitted", "prior", "refused"}
    assert {p.scale for p in parameters.values()} == {"biological", "assay", "physical"}
    assert len(Tag.ALL) == 7
    assert parameters["heat.hsp70_hsf1_k2"].status == "fitted"
    assert parameters["heat.hsp70_hsf1_k2"].parameter.tag == Tag.ASSERTED
    assert parameters["heat.hsp70_hsf1_k2"].parameter.is_free
    assert parameters["oxidative.membrane_entry_alpha"].status == "prior"
    assert parameters["native.glucose_uptake_D0p10"].status == "measured"
    for key in ("upr.dtt_to_er_client", "ph.permeability_lactic_exact", "assay.od_to_dry_mass", "assay.reporter_gain"):
        assert parameters[key].parameter.value is None
        with pytest.raises(RefusedValue):
            float(parameters[key].parameter)
    assert parameters["assay.reporter_gain"].scale == "assay"
    assert parameters["native.tdh3_abundance"].identifiability == "composite_only"


def test_primary_sd_is_not_silently_relabelled_as_a_confidence_interval(inventory):
    parameter = next(p for p in inventory.parameters if p.parameter.name == "ph.permeability_acetic")
    assert parameter.uncertainty.kind == "standard_deviation"
    assert parameter.uncertainty.value == 2e-6
    assert parameter.uncertainty.units == "cm/s"
    assert parameter.uncertainty.confidence is None
    assert parameter.uncertainty.interval is None
    assert parameter.parameter.ci95 is None
    assert parameter.to_dict()["uncertainty"]["kind"] == "standard_deviation"


@pytest.mark.parametrize("kwargs", [
    {"kind": "standard_deviation", "value": 1, "confidence": 0.95},
    {"kind": "standard_deviation", "value": -1},
    {"kind": "standard_deviation", "value": True},
    {"kind": "standard_error", "value": float("nan")},
    {"kind": "not_reported", "value": 0},
    {"kind": "not_reported", "interval": (0, 1)},
    {"kind": "confidence_interval", "interval": (0, 1)},
    {"kind": "confidence_interval", "interval": (2, 1), "confidence": 0.95},
    {"kind": "sensitivity_range", "interval": (0, 1), "confidence": 0.95},
])
def test_uncertainty_rejects_semantic_shortcuts(kwargs):
    with pytest.raises(ValueError):
        ParameterUncertainty(units="s^-1", explanation="test source", **kwargs)


def test_uncertainty_roundtrip_preserves_kind_units_and_level():
    uncertainty = ParameterUncertainty("confidence_interval", "h^-1", "reported interval", interval=[1, 2], confidence=0.90)
    assert ParameterUncertainty(**uncertainty.to_dict()) == uncertainty
    assert uncertainty.confidence == 0.90


@pytest.mark.parametrize("parameter_id", ["heat.hsp70_hsf1_k2", "oxidative.membrane_entry_alpha"])
def test_a_primary_table_does_not_make_its_fitted_or_chosen_constants_measured(inventory, parameter_id):
    payload = inventory.to_dict()
    row = next(row for row in payload["parameters"] if row["id"] == parameter_id)
    row["status"] = "measured"
    with pytest.raises(ValueError, match="checked source measurement"):
        validate_inventory(payload)


@pytest.mark.parametrize("field,value", [("value", 1.2), ("units", "RFU"), ("locator", "somewhere else")])
def test_measured_value_must_match_the_exact_checked_source_tuple(inventory, field, value):
    payload = inventory.to_dict()
    row = next(row for row in payload["parameters"] if row["id"] == "native.glucose_uptake_D0p10")
    row[field] = value
    if field == "units":
        row["uncertainty"]["units"] = value
    with pytest.raises(ValueError, match="checked source measurement"):
        validate_inventory(payload)


def test_firewall_and_unresolved_directions_are_not_optional(inventory):
    payload = inventory.to_dict()
    payload["firewall"]["product_outcomes_allowed"] = True
    with pytest.raises(ValueError, match="firewall"):
        validate_inventory(payload)
    payload = inventory.to_dict()
    payload["parameters"][0]["unresolved_directions"] = []
    with pytest.raises(ValueError, match="unresolved directions"):
        validate_inventory(payload)
    payload = inventory.to_dict()
    payload["native_data"][0]["artifact_paths"].append("data/holdout_transfer/labels.json")
    with pytest.raises(ValueError, match="explicitly pinned"):
        validate_inventory(payload)


def test_inventory_is_detached_and_json_roundtrippable(inventory, tmp_path):
    payload = inventory.to_dict()
    payload["parameters"][0]["value"] = 123
    assert inventory.to_dict()["parameters"][0]["value"] != 123
    path = tmp_path / "evidence.json"
    path.write_text(json.dumps(inventory.to_dict()), encoding="utf-8")
    assert load_parameter_evidence(path).sha256 == inventory.sha256
    path.write_text('{"schema_version":1,"schema_version":1}', encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        load_parameter_evidence(path)


def test_missing_and_changed_source_bytes_are_refused(inventory, tmp_path):
    source = inventory.to_dict()["models"]["hog2013_wt"]
    path = tmp_path / source["path"]
    path.parent.mkdir(parents=True)
    path.write_bytes((ROOT / source["path"]).read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256"):
        load_reference_model("hog2013_wt", inventory, root=tmp_path)
    with pytest.raises(FileNotFoundError):
        load_reference_model("hog2013_hog1_del", inventory, root=tmp_path)
    payload = inventory.to_dict()
    payload["models"]["hog2013_wt"]["model_id"] = "mismatched_caption"
    with pytest.raises(ValueError, match="identity mismatch"):
        load_reference_model("hog2013_wt", validate_inventory(payload))


def test_freeze_reads_only_native_reference_sources_and_detects_changed_inventory(inventory):
    frozen = freeze_evidence(inventory, source_ids=["hog2013"], product_outcomes_seen=False)
    assert len(frozen.to_dict()["assets"]) == 6
    assert all(asset["path"].endswith(".xml") for asset in frozen.to_dict()["assets"])
    verify_evidence_freeze(frozen, inventory)
    for selected in (["hog2013_assays"], ["granados2018_metadata"], ["product_outcomes"]):
        with pytest.raises(ValueError, match="allowlisted"):
            freeze_evidence(inventory, source_ids=selected, product_outcomes_seen=False)
    for flag in (True, 0, None):
        with pytest.raises(ValueError, match="no-product-outcomes"):
            freeze_evidence(inventory, source_ids=["hog2013"], product_outcomes_seen=flag)
    changed = inventory.to_dict()
    changed["scope"] += " changed"
    with pytest.raises(ValueError, match="freeze"):
        verify_evidence_freeze(frozen, validate_inventory(changed))


def test_source_paths_cannot_escape_or_read_unregistered_reserved_assets(tmp_path):
    for relative in ("../outside", "/absolute", "data/../outside", "data//x", "data\\x"):
        with pytest.raises(ValueError):
            source_path(tmp_path, relative)
    link = tmp_path / "link"
    link.symlink_to(tmp_path.parent, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        source_path(tmp_path, "link/outside")


def test_native_data_inventory_covers_support_assays_but_never_sealed_response_paths(inventory):
    kinds = {record["kind"] for record in inventory.to_dict()["native_data"]}
    assert {"exchange", "expression", "protein_abundance", "turnover", "cell_mass", "reporter"} <= kinds
    assets = native_asset_inventory(inventory)
    assert assets and all("/sealed/" not in asset["path"] for asset in assets)
    assert all(asset["source_id"] != "product_outcomes" for asset in assets)
    report = evidence_report(inventory)
    assert report["full_biological_coverage"] is False
    assert report["independent_real_validation_established"] is False
