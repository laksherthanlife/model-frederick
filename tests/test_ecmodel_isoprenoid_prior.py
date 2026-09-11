from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

from ystwin import paths
from ystwin.pathway import enzyme_capacity, proteome

SCRIPT = paths.REPO_ROOT / "scripts" / "ecmodel_isoprenoid_prior.py"
SBML = "{http://www.sbml.org/sbml/level3/version1/core}"
FBC = "{http://www.sbml.org/sbml/level3/version1/fbc/version2}"
INVENTORY_FIELDS = (
    "reaction_id", "declared_gpr_genes", "gpr_genes", "branch", "protein_species_id",
    "protein_coefficient_h", "kcat_per_s", "is_reverse_arm", "included_in_forward_span",
    "association_basis", "interpretation", "source_model_sha256", "source_sbml_sha256",
    "boundary_basis",
)
SWEEP_FIELDS = (
    "enzyme", "kcat_per_s", "mass_fraction_of_protein", "enzyme_molar_mass_g_per_mol",
    "total_protein_g_per_gdcw", "growth_rate_per_h", "enzyme_mmol_per_gdcw",
    "vmax_mmol_per_gdcw_h", "capacity_mmol_per_gdcw", "reference_content_source",
    "reference_content_mg_per_gdcw", "reference_product_molar_mass_g_per_mol",
    "reference_content_mmol_per_gdcw", "fold_over_reference_content",
    "reference_capacity_source", "reference_capacity_mmol_per_gdcw",
    "fold_over_reference_capacity", "reference_basis", "interpretation", "kcat_basis",
    "source_model_sha256", "source_sbml_sha256", "boundary_basis",
)


@pytest.fixture(scope="module")
def prior():
    spec = importlib.util.spec_from_file_location("ecmodel_isoprenoid_prior_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop(spec.name, None)


def _model_xml():
    root = ET.Element(SBML + "sbml", {"level": "3", "version": "1"})
    model = ET.SubElement(root, SBML + "model", {"id": "synthetic_gecko2"})
    definitions = ET.SubElement(model, SBML + "listOfUnitDefinitions")
    definition = ET.SubElement(definitions, SBML + "unitDefinition", {
        "id": "mmol_per_gDW_per_hr",
    })
    units = ET.SubElement(definition, SBML + "listOfUnits")
    for kind, exponent, scale, multiplier in (
            ("mole", "1", "-3", "1"), ("gram", "-1", "0", "1"),
            ("second", "-1", "0", "3600")):
        ET.SubElement(units, SBML + "unit", {
            "kind": kind, "exponent": exponent, "scale": scale, "multiplier": multiplier,
        })
    species = ET.SubElement(model, SBML + "listOfSpecies")
    for name in ("prot_erg1", "prot_hmg1", "prot_aux", "prot_pool", "metab_prot_decoy"):
        ET.SubElement(species, SBML + "species", {"id": name})
    gene_products = ET.SubElement(model, FBC + "listOfGeneProducts")
    for gene in ("ERG1", "HMG1", "AUX", "OTHER"):
        ET.SubElement(gene_products, FBC + "geneProduct", {
            FBC + "id": "G_" + gene, FBC + "label": gene,
        })
    reactions = ET.SubElement(model, SBML + "listOfReactions")
    cases = (
        ("r_forward", ("ERG1",), (("prot_erg1", 1 / 7200),)),
        ("r_complex", ("HMG1", "AUX"), (("prot_hmg1", 1 / 3600), ("prot_aux", 1 / 14400))),
        ("r_slow", ("ERG1",), (("prot_erg1", 1 / 360),)),
        ("r_forward_REV", ("ERG1",), (("prot_erg1", 1 / 3.6e9),)),
        ("r_no_turnover", ("ERG1",), (("prot_pool", 1), ("metab_prot_decoy", 1))),
        ("r_other_gene", ("OTHER",), (("prot_aux", 1),)),
        ("arm_pooling", ("ERG1",), (("prot_erg1", 1e-20),)),
        ("draw_prot_erg1", ("ERG1",), (("prot_erg1", 1e-20),)),
        ("r_no_gpr", (), (("prot_erg1", 1e-20),)),
    )
    for reaction_id, genes, coefficients in cases:
        reaction = ET.SubElement(reactions, SBML + "reaction", {
            "id": reaction_id, "reversible": "false",
        })
        reactants = ET.SubElement(reaction, SBML + "listOfReactants")
        for name, coefficient in coefficients:
            ET.SubElement(reactants, SBML + "speciesReference", {
                "species": name, "stoichiometry": str(coefficient), "constant": "true",
            })
        if genes:
            association = ET.SubElement(reaction, FBC + "geneProductAssociation")
            if len(genes) > 1:
                association = ET.SubElement(association, FBC + "and")
            for gene in genes:
                ET.SubElement(association, FBC + "geneProductRef", {FBC + "geneProduct": "G_" + gene})
    return ET.tostring(root)


@pytest.fixture
def model_path(tmp_path):
    path = tmp_path / "model.xml.gz"
    path.write_bytes(gzip.compress(_model_xml(), mtime=0))
    return path


def _csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return tuple(reader.fieldnames), list(reader)


def _cli(tmp_path, *args, model=None):
    env = {
        **os.environ,
        "YSTWIN_EC_YEAST_GEM": str(model or tmp_path / "absent.xml"),
        "YSTWIN_OUTPUTS": str(tmp_path / "output"),
    }
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)], cwd=tmp_path, env=env,
        text=True, capture_output=True, timeout=30,
    )


def test_import_has_no_data_access_or_solver_side_effects(tmp_path):
    code = (
        "import runpy, sys; "
        f"runpy.run_path({str(SCRIPT)!r}, run_name='import_only'); "
        "assert not {'cobra', 'ystwin.predict', 'ystwin.kinetic.carotenoid', "
        "'ystwin.pathway.published_cassettes'} & sys.modules.keys()"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=tmp_path,
        env={**os.environ, "YSTWIN_EC_YEAST_GEM": str(tmp_path / "absent.xml"),
             "YSTWIN_OUTPUTS": str(tmp_path / "output")},
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert not (tmp_path / "output").exists()


def test_help_is_usable_without_model_or_output_directory(tmp_path):
    result = _cli(tmp_path, "--help")
    assert result.returncode == 0, result.stderr
    for flag in ("--model", "--output-dir", "--kcats-per-s", "--mass-fractions",
                 "--growth-rate-per-h", "--enzyme-molar-mass-g-per-mol",
                 "--reference-content-mg-per-gdcw", "--reference-capacity-mmol-per-gdcw"):
        assert flag in result.stdout
    assert "prior" in result.stdout.lower()
    assert not (tmp_path / "output").exists()


def test_inventory_uses_all_coefficients_and_reaction_level_gpr_membership(prior, model_path):
    inventory = prior.parse_kcats(model_path)
    assert inventory.reactions_with_gpr == 6
    assert inventory.reactions_with_kcat == 5
    assert len(inventory.rows) == 5
    assert inventory.model_id == "synthetic_gecko2"
    assert inventory.source_model_sha256 == hashlib.sha256(model_path.read_bytes()).hexdigest()
    assert inventory.source_sbml_sha256 == hashlib.sha256(_model_xml()).hexdigest()
    complex_rows = [row for row in inventory.rows if row["reaction_id"] == "r_complex"]
    assert {row["protein_species_id"] for row in complex_rows} == {"prot_hmg1", "prot_aux"}
    assert {row["kcat_per_s"] for row in complex_rows} == {1.0, 4.0}
    assert all(row["declared_gpr_genes"] == "HMG1" for row in complex_rows)
    assert all(row["gpr_genes"] == "AUX;HMG1" for row in complex_rows)
    assert all(row["branch"] == "MVA" for row in complex_rows)
    assert next(row for row in inventory.rows if row["reaction_id"] == "r_forward")["branch"] == "sterol"
    for row in inventory.rows:
        assert tuple(row) == INVENTORY_FIELDS
        assert row["association_basis"] == "reaction_GPR_membership_not_gene_to_protein_assignment"
        assert row["interpretation"] == "model_encoded_prior_not_measurement"
        assert row["protein_coefficient_h"] * row["kcat_per_s"] * 3600 == pytest.approx(1.0)
        assert row["source_model_sha256"] == inventory.source_model_sha256
        assert row["source_sbml_sha256"] == inventory.source_sbml_sha256
        assert row["boundary_basis"] == prior.BOUNDARY_BASIS


def test_plain_and_gzip_inputs_preserve_sbml_identity(prior, model_path, tmp_path):
    plain = tmp_path / "plain.xml"
    plain.write_bytes(_model_xml())
    compressed = prior.parse_kcats(model_path)
    uncompressed = prior.parse_kcats(plain)
    assert uncompressed.source_sbml_sha256 == compressed.source_sbml_sha256
    assert uncompressed.source_model_sha256 != compressed.source_model_sha256
    assert [row["kcat_per_s"] for row in uncompressed.rows] == [row["kcat_per_s"] for row in compressed.rows]


def test_reverse_is_reported_but_never_sets_forward_span(prior, model_path):
    inventory = prior.parse_kcats(model_path)
    low, high = prior.forward_span(inventory.rows)
    assert low["kcat_per_s"] == pytest.approx(0.1)
    assert high["kcat_per_s"] == 4.0
    reverse = [row for row in inventory.rows if row["is_reverse_arm"]]
    assert len(reverse) == 1
    assert reverse[0]["kcat_per_s"] == pytest.approx(1e6)
    assert not reverse[0]["included_in_forward_span"]
    with pytest.raises(ValueError, match="forward"):
        prior.forward_span(reverse)
    with pytest.raises(ValueError, match="forward"):
        prior.forward_span([])


@pytest.mark.parametrize("bad", ["0", "-1", "nan", "inf", "1e-320", "1e308", "invalid"])
def test_invalid_protein_coefficients_fail_instead_of_disappearing(prior, tmp_path, bad):
    root = ET.fromstring(_model_xml())
    ref = root.find(f".//{SBML}reaction[@id='r_other_gene']/{SBML}listOfReactants/{SBML}speciesReference")
    ref.set("stoichiometry", bad)
    path = tmp_path / "invalid.xml"
    path.write_bytes(ET.tostring(root))
    with pytest.raises(ValueError, match="r_other_gene.*prot_aux"):
        prior.parse_kcats(path)


@pytest.mark.parametrize("defect", ["unknown_gene", "unknown_protein", "duplicate_protein", "reversible", "time_units", "missing_unit_kind", "missing_units", "missing_model", "no_gpr", "no_branch"])
def test_structural_or_unit_ambiguity_fails_closed(prior, tmp_path, defect):
    root = ET.fromstring(_model_xml())
    reaction = root.find(f".//{SBML}reaction[@id='r_forward']")
    if defect == "unknown_gene":
        reaction.find(f".//{FBC}geneProductRef").set(FBC + "geneProduct", "missing")
    elif defect == "unknown_protein":
        reaction.find(f".//{SBML}speciesReference").set("species", "prot_unknown")
    elif defect == "duplicate_protein":
        refs = reaction.find(SBML + "listOfReactants")
        refs.append(ET.fromstring(ET.tostring(refs[0])))
    elif defect == "reversible":
        reaction.set("reversible", "true")
    elif defect == "time_units":
        root.find(f".//{SBML}unit[@kind='second']").set("multiplier", "1")
    elif defect == "missing_unit_kind":
        root.find(f".//{SBML}unit[@kind='second']").attrib.pop("kind")
    elif defect == "missing_units":
        model = root.find(SBML + "model")
        model.remove(model.find(SBML + "listOfUnitDefinitions"))
    elif defect == "missing_model":
        root.remove(root.find(SBML + "model"))
    elif defect == "no_gpr":
        for element in root.iter(SBML + "reaction"):
            association = element.find(FBC + "geneProductAssociation")
            if association is not None:
                element.remove(association)
    elif defect == "no_branch":
        for gene in root.iter(FBC + "geneProduct"):
            gene.set(FBC + "label", "OTHER")
    path = tmp_path / "invalid.xml"
    path.write_bytes(ET.tostring(root))
    with pytest.raises(ValueError):
        prior.parse_kcats(path)


def test_mass_content_and_hour_conversion_use_current_interfaces(prior, monkeypatch):
    calls = {"mass": [], "capacity": []}
    original_mass = proteome.enzyme_content_from_mass_fraction
    original_capacity = enzyme_capacity.derived_capacity

    def mass(*args, **kwargs):
        calls["mass"].append((args, kwargs))
        return original_mass(*args, **kwargs)

    def capacity(*args, **kwargs):
        calls["capacity"].append((args, kwargs))
        return original_capacity(*args, **kwargs)

    monkeypatch.setattr(proteome, "enzyme_content_from_mass_fraction", mass)
    monkeypatch.setattr(enzyme_capacity, "derived_capacity", capacity)
    row, = prior.sweep_capacities(
        [2.0], [0.01], enzyme="synthetic_enzyme", enzyme_molar_mass_g_per_mol=100_000,
        total_protein_g_per_gdcw=0.5, growth_rate_per_h=0.2,
    )
    assert calls["mass"] and calls["capacity"]
    assert row["enzyme_mmol_per_gdcw"] == pytest.approx(5e-5)
    assert row["vmax_mmol_per_gdcw_h"] == pytest.approx(0.36)
    assert row["capacity_mmol_per_gdcw"] == pytest.approx(1.8)
    assert row["interpretation"] == "conditional_sensitivity_not_product_prediction"
    assert row["reference_content_mmol_per_gdcw"] is None
    assert row["fold_over_reference_capacity"] is None
    molar = proteome.enzyme_content_mmol_per_gdcw(10_000, 0.5)
    assert row["enzyme_mmol_per_gdcw"] / molar == pytest.approx(0.5)


@pytest.mark.parametrize("change, factor", [
    ({"kcats_per_s": [4.0]}, 2.0), ({"mass_fractions": [0.02]}, 2.0),
    ({"enzyme_molar_mass_g_per_mol": 200_000}, 0.5),
    ({"growth_rate_per_h": 0.4}, 0.5), ({"total_protein_g_per_gdcw": 0.25}, 0.5),
])
def test_capacity_physical_scaling_invariants(prior, change, factor):
    kwargs = dict(kcats_per_s=[2.0], mass_fractions=[0.01],
                  enzyme_molar_mass_g_per_mol=100_000, total_protein_g_per_gdcw=0.5,
                  growth_rate_per_h=0.2)
    base, = prior.sweep_capacities(**kwargs)
    changed, = prior.sweep_capacities(**(kwargs | change))
    assert changed["capacity_mmol_per_gdcw"] == pytest.approx(base["capacity_mmol_per_gdcw"] * factor)


@pytest.mark.parametrize("field", ["kcats_per_s", "mass_fractions", "growth_rate_per_h", "enzyme_molar_mass_g_per_mol", "total_protein_g_per_gdcw"])
@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_sweep_rejects_nonphysical_or_nonfinite_inputs(prior, field, bad):
    kwargs = dict(kcats_per_s=[2.0], mass_fractions=[0.01])
    kwargs[field] = [bad] if field in ("kcats_per_s", "mass_fractions") else bad
    with pytest.raises(ValueError, match=field):
        prior.sweep_capacities(**kwargs)


@pytest.mark.parametrize("kwargs", [
    {"kcats_per_s": []}, {"mass_fractions": []}, {"mass_fractions": [1.01]},
    {"total_protein_g_per_gdcw": 1.01}, {"enzyme": " "},
    {"kcats_per_s": [1e308]}, {"growth_rate_per_h": 1e-320},
    {"enzyme_molar_mass_g_per_mol": 1e308, "mass_fractions": [1e-300]},
])
def test_sweep_rejects_empty_axes_and_nonfinite_derived_results(prior, kwargs):
    with pytest.raises(ValueError):
        prior.sweep_capacities(**(dict(kcats_per_s=[2.0], mass_fractions=[0.01]) | kwargs))


def test_reference_measurements_only_change_comparison_columns(prior):
    content = prior.ContentReference("synthetic assay", 10.0, 500.0)
    capacity = prior.CapacityReference("synthetic fitted comparator", 0.004)
    kwargs = dict(kcats_per_s=[2.0], mass_fractions=[0.01],
                  enzyme_molar_mass_g_per_mol=100_000, total_protein_g_per_gdcw=0.5,
                  growth_rate_per_h=0.2)
    no_reference, = prior.sweep_capacities(**kwargs)
    reference, = prior.sweep_capacities(**kwargs, content_reference=content, capacity_reference=capacity)
    assert reference["reference_content_mmol_per_gdcw"] == pytest.approx(0.02)
    assert reference["fold_over_reference_content"] == pytest.approx(90.0)
    assert reference["fold_over_reference_capacity"] == pytest.approx(450.0)
    assert reference["reference_basis"] == "user_declared_comparison_only_not_calibration"
    for field in ("enzyme_mmol_per_gdcw", "vmax_mmol_per_gdcw_h", "capacity_mmol_per_gdcw"):
        assert reference[field] == no_reference[field]


@pytest.mark.parametrize("args", [
    ("--kcats-per-s", "nan"), ("--mass-fractions", "2"), ("--growth-rate-per-h", "0"),
    ("--reference-content-mg-per-gdcw", "10"),
    ("--reference-capacity-mmol-per-gdcw", "0.01"),
    ("--reference-capacity-mmol-per-gdcw", "inf", "--reference-capacity-source", "synthetic"),
])
def test_bad_cli_inputs_fail_without_creating_output(tmp_path, model_path, args):
    result = _cli(tmp_path, *args, model=model_path)
    assert result.returncode == 2, result.stdout + result.stderr
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("kind", ["missing", "bad_xml", "bad_gzip", "bad_deflate", "truncated_gzip"])
def test_bad_models_fail_without_fallback_or_partial_output(tmp_path, model_path, kind):
    invalid = tmp_path / "invalid.xml.gz"
    if kind == "bad_xml":
        invalid.write_bytes(b"not SBML")
    elif kind == "bad_gzip":
        invalid.write_bytes(b"\x1f\x8bnot gzip")
    elif kind == "bad_deflate":
        invalid.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03\x07" + b"\x00" * 8)
    elif kind == "truncated_gzip":
        invalid.write_bytes(model_path.read_bytes()[:40])
    result = _cli(tmp_path, "--model", invalid, model=model_path)
    assert result.returncode == 2
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "output").exists()


def test_missing_environment_override_does_not_read_vendored_fallback(tmp_path):
    result = _cli(tmp_path)
    assert result.returncode == 2
    assert "YSTWIN_EC_YEAST_GEM" in result.stderr
    assert not (tmp_path / "output").exists()


def test_cli_exports_schema_provenance_and_all_declared_inputs(tmp_path, model_path):
    output = tmp_path / "chosen-output"
    result = _cli(
        tmp_path, "--model", model_path, "--output-dir", output,
        "--kcats-per-s", "2", "4", "--mass-fractions", "0.01", "0.02",
        "--enzyme", "synthetic_enzyme", "--enzyme-molar-mass-g-per-mol", "100000",
        "--total-protein-g-per-gdcw", "0.5", "--growth-rate-per-h", "0.2",
        "--reference-content-mg-per-gdcw", "10", "--reference-product-molar-mass-g-per-mol", "500",
        "--reference-content-source", "synthetic assay", "--reference-capacity-mmol-per-gdcw", "0.004",
        "--reference-capacity-source", "synthetic fitted comparator",
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "output").exists()
    assert {p.name for p in output.iterdir()} == {
        "ecmodel_isoprenoid_kcats.csv", "enzyme_capacity_sweep.csv", "ecmodel_isoprenoid_prior.json",
    }
    inventory_fields, inventory = _csv(output / "ecmodel_isoprenoid_kcats.csv")
    sweep_fields, sweep = _csv(output / "enzyme_capacity_sweep.csv")
    assert inventory_fields == INVENTORY_FIELDS
    assert sweep_fields == SWEEP_FIELDS
    assert len(inventory) == 5
    assert len(sweep) == 4
    report = json.loads((output / "ecmodel_isoprenoid_prior.json").read_text())
    assert report["schema_version"] == 1
    assert report["calibration_performed"] is False
    assert report["solver_used"] is False
    assert report["source_model"]["sha256"] == hashlib.sha256(model_path.read_bytes()).hexdigest()
    assert report["source_model"]["sbml_sha256"] == hashlib.sha256(_model_xml()).hexdigest()
    assert report["inventory"]["reactions_with_gpr"] == 6
    assert report["inventory"]["reactions_with_kcat"] == 5
    assert report["inventory"]["coverage_fraction"] == pytest.approx(5 / 6)
    assert set(report["inventory"]["declared_genes"]) == set(enzyme_capacity.ISOPRENOID_BRANCH_GENES)
    assert "ERG6" in report["inventory"]["genes_without_forward_turnovers"]
    assert report["sweep"]["kcat_basis"] == "user_declared_sensitivity_axis_not_measurement"
    assert report["sweep"]["kcats_per_s"] == [2.0, 4.0]
    assert report["sweep"]["mass_fractions_of_protein"] == [0.01, 0.02]
    assert report["sweep"]["content_reference"]["source"] == "synthetic assay"
    assert report["sweep"]["capacity_reference"]["source"] == "synthetic fitted comparator"
    for source in report["source_code"]:
        assert source["sha256"] == hashlib.sha256(pathlib.Path(source["path"]).read_bytes()).hexdigest()
    assert len(report["source_code"]) == 3
    for name, metadata in report["outputs"].items():
        assert metadata["sha256"] == hashlib.sha256((output / name).read_bytes()).hexdigest()
    for row in sweep:
        assert row["source_model_sha256"] == report["source_model"]["sha256"]
        assert row["boundary_basis"] == report["boundary_basis"]
        expected_enzyme = float(row["mass_fraction_of_protein"]) * 0.5 / 100_000 * 1000
        assert float(row["enzyme_mmol_per_gdcw"]) == pytest.approx(expected_enzyme)
        expected_capacity = expected_enzyme * float(row["kcat_per_s"]) * 3600 / 0.2
        assert float(row["capacity_mmol_per_gdcw"]) == pytest.approx(expected_capacity)
    assert "not measurement" in result.stdout.lower()


def test_default_sweep_uses_observed_forward_coefficients_and_no_product_outcomes(tmp_path, model_path):
    result = _cli(tmp_path, model=model_path)
    assert result.returncode == 0, result.stderr
    output = tmp_path / "output"
    report = json.loads((output / "ecmodel_isoprenoid_prior.json").read_text())
    assert report["sweep"]["kcats_per_s"] == pytest.approx([0.1, 1.0, 2.0, 4.0])
    assert report["sweep"]["kcat_basis"] == "forward_model_encoded_turnovers_as_sensitivity_priors"
    assert report["sweep"]["content_reference"] is None
    assert report["sweep"]["capacity_reference"] is None
    _, sweep = _csv(output / "enzyme_capacity_sweep.csv")
    assert len(sweep) == 4 * len(report["sweep"]["mass_fractions_of_protein"])
    assert all(row["fold_over_reference_content"] == row["fold_over_reference_capacity"] == "" for row in sweep)


def test_output_io_failure_is_reported_without_traceback(tmp_path, model_path):
    output = tmp_path / "occupied"
    output.write_text("keep me")
    result = _cli(tmp_path, "--output-dir", output, model=model_path)
    assert result.returncode == 2
    assert "Traceback" not in result.stderr
    assert output.read_text() == "keep me"


@pytest.mark.integration
def test_vendored_gecko2_bounded_smoke(tmp_path):
    model = paths.data_dir() / "gem" / "ecYeastGEM_batch.xml.gz"
    if not model.is_file():
        pytest.skip("vendored GECKO2 model is absent")
    result = _cli(tmp_path, "--model", model, "--mass-fractions", "0.001", "0.005")
    assert result.returncode == 0, result.stderr
    output = tmp_path / "output"
    report = json.loads((output / "ecmodel_isoprenoid_prior.json").read_text())
    assert report["source_model"]["sbml_sha256"] == "3a9d97a68e2408ea243114794c027b11d70f253ecf6a758f4bcd5f2de50170b8"
    assert report["inventory"]["reactions_with_gpr"] == 4258
    assert report["inventory"]["reactions_with_kcat"] == 3856
    assert report["inventory"]["genes_without_forward_turnovers"] == []
    low = report["inventory"]["forward_span"]["minimum"]
    high = report["inventory"]["forward_span"]["maximum"]
    assert low["declared_gpr_genes"] == "ERG6"
    assert high["declared_gpr_genes"] == "IDI1"
    assert low["kcat_per_s"] == pytest.approx(0.011, rel=1e-4)
    assert high["kcat_per_s"] == pytest.approx(29900, rel=1e-4)
    _, rows = _csv(output / "ecmodel_isoprenoid_kcats.csv")
    assert len(rows) == 39
    assert sum(row["is_reverse_arm"] == "True" for row in rows) == 3
    complex_rows = [row for row in rows if row["reaction_id"] == "r_0317No1"]
    assert len(complex_rows) == 2
    assert all(row["gpr_genes"] == "ERG11;NCP1" for row in complex_rows)
    _, sweep = _csv(output / "enzyme_capacity_sweep.csv")
    assert len(sweep) == 2 * len(report["sweep"]["kcats_per_s"])
