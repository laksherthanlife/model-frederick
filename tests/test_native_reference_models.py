from __future__ import annotations

import ast
import base64
import hashlib
import importlib.util
import io
import json
import math
import operator
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace
import warnings
import xml.etree.ElementTree as ET
import zipfile

import libsbml
import numpy as np
import pytest

from ystwin.analysis.parameter_evidence import (
    JALIHAL_EXPORT_ADAPTER,
    JALIHAL_NORMALIZED_INPUTS,
    JALIHAL_SBML_SHA256,
    EvidenceGap,
    audit_local_sources,
    freeze_evidence,
    load_parameter_evidence,
    load_reference_model,
    load_source_asset,
    model_parameter_inventory,
    simulate_native_reference,
    validate_inventory,
)
from ystwin.mech.kinetic_sbml import KineticModel, UnsupportedSBMLError

ROOT = Path(__file__).resolve().parents[1]
ASSETS = "data/native_reference_models"

#: jalihal2021 is gitignored on purpose: unlike the CC-BY article bodies, its assets come from a
#: pinned GitHub commit whose redistribution rights are unverified, so a clean clone lacks them.
_RESTRICTED = ROOT / ASSETS / "jalihal2021"
_RESTRICTED_REASON = (f"{ASSETS}/jalihal2021 is not redistributed; fetch it with "
                      "'python scripts/fetch_native_reference_models.py' to run these checks")
requires_restricted_assets = pytest.mark.skipif(not _RESTRICTED.is_dir(), reason=_RESTRICTED_REASON)


@pytest.fixture(scope="module")
def inventory():
    return load_parameter_evidence(ROOT / "data/parameter_evidence.json")


@pytest.fixture(scope="module")
def fetcher():
    spec = importlib.util.spec_from_file_location("native_source_fetcher", ROOT / "scripts/fetch_native_reference_models.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _source(inventory, source, name):
    return load_source_asset(source, f"{ASSETS}/{source}/{name}", inventory, root=ROOT)


def _table(payload, converter=str):
    return {name: converter(value.strip()) for name, value in
            (line.split("\t", 1) for line in payload.decode("utf-8").splitlines())}


@pytest.fixture(scope="module")
def nutrient_model(inventory):
    if not _RESTRICTED.is_dir():
        pytest.skip(_RESTRICTED_REASON)
    return load_reference_model("jalihal2021", inventory, root=ROOT)


@requires_restricted_assets
def test_every_new_asset_is_pinned_and_audits_do_not_double_count_the_sbml(inventory, fetcher):
    entries = fetcher.asset_entries(inventory)
    assert len(entries) == 37
    assert {e["source_id"] for e in entries} == {
        "jalihal2021", "zheng2016", "goulev2017", "pincus2010", "ke2013", "krakowiak2018",
    }
    for entry in entries:
        payload = (ROOT / entry["path"]).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == entry["sha256"]
        fetcher.verify_payload(payload, entry)
        assert not entry["path"].lower().endswith((".md", "/readme"))
    checks = audit_local_sources(inventory, ROOT)
    declared = inventory.to_dict()
    expected = {(source_id, asset["path"]) for source_id, source in declared["sources"].items()
                for asset in source["local_artifacts"]}
    expected.update((model["source_id"], model["path"]) for model in declared["models"].values()
                    if model["availability"] == "local_verified")
    assert ("native_reporter_plates", "data/plates/processing_provenance.json") in expected
    assert {(check["source_id"], check["path"]) for check in checks} == expected
    assert len(checks) == len(expected)
    assert all(a["verified"] for a in checks)
    github = [a for a in checks if a.get("acquisition", {}).get("kind") == "github"]
    assert len(github) == 7
    assert all(a["actual_git_blob"] == a["acquisition"]["git_blob"] for a in github)
    frozen = freeze_evidence(inventory, source_ids=["jalihal2021"], product_outcomes_seen=False, root=ROOT)
    assert len(frozen.to_dict()["assets"]) == 8


def test_source_loading_is_explicitly_allowlisted_and_rejects_changed_bytes(inventory, tmp_path):
    with pytest.raises(ValueError, match="registered"):
        load_source_asset("jalihal2021", "data/holdout_transfer/labels.json", inventory, root=ROOT)
    record = inventory.to_dict()["models"]["jalihal2021"]
    path = tmp_path / record["path"]
    path.parent.mkdir(parents=True)
    path.write_bytes((ROOT / record["path"]).read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256"):
        load_reference_model("jalihal2021", inventory, root=tmp_path)


def test_native_export_issues_are_audited_not_silently_fixed_globally(inventory, nutrient_model):
    path = ROOT / inventory.to_dict()["models"]["jalihal2021"]["path"]
    original = path.read_bytes()
    document = libsbml.readSBMLFromString(original.decode())
    source = document.getModel()
    assert source.getId() == ""
    assert (source.getNumSpecies(), source.getNumParameters(), source.getNumReactions()) == (25, 128, 25)
    assert all(not r.getProduct(0).isSetStoichiometry() for r in source.getListOfReactions())
    with pytest.raises(UnsupportedSBMLError, match="stoichiometry"):
        KineticModel.from_sbml(path)
    # Once missing stoichiometry is supplied, the generic loader still refuses min.
    for reaction in source.getListOfReactions():
        reaction.getProduct(0).setStoichiometry(1.0)
    with pytest.raises(UnsupportedSBMLError, match="math node"):
        KineticModel(document, path, JALIHAL_SBML_SHA256)
    audit = nutrient_model.metadata["native_reference"]
    assert audit["adapter"] == JALIHAL_EXPORT_ADAPTER
    assert audit["binary_min_nodes_lowered"] == 5
    assert len(audit["unit_product_coefficients_supplied"]) == 25
    assert audit["source_bytes_modified"] is False
    assert audit["execution_sbml_sha256"] != JALIHAL_SBML_SHA256
    assert "min(" in audit["source_reaction_formulas"]["Gcn4_r"]
    assert "piecewise(" in nutrient_model.metadata["reaction_formulas"]["Gcn4_r"]
    assert path.read_bytes() == original


def test_all_source_parameters_and_initial_conditions_match_the_correct_upstream_directory(inventory, nutrient_model):
    pars = _table(_source(inventory, "jalihal2021", "parameters.txt"), float)
    ics = _table(_source(inventory, "jalihal2021", "initialconditions.txt"), float)
    assert pars == nutrient_model.parameter_values
    np.testing.assert_array_equal(nutrient_model.initial_state(), [ics[s] for s in nutrient_model.species_ids])
    result = model_parameter_inventory("jalihal2021", inventory, root=ROOT)
    assert len(result["parameters"]) == 128
    assert all(p["status"] == "prior" and p["units"] == "unspecified" for p in result["parameters"])
    assert all(p["uncertainty"]["kind"] == "not_reported" for p in result["parameters"])
    assert {p["id"] for p in result["parameters"] if p["scale"] == "normalized_input"} == set(JALIHAL_NORMALIZED_INPUTS)
    assert all(p["unresolved_directions"] for p in result["parameters"])
    assert "Sch9 * PKA" in result["reaction_formulas"]["Dot6_r"]
    assert "min(" in result["reaction_formulas"]["Protein_r"]


def _source_expression(node, values):
    # Evaluate only the arithmetic/calls in the separately pinned upstream ODE text.
    # Do not import or exec the source simulator (which can compile code/write files).
    operations = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
    functions = {"min": min, "tRNA": min, "pRib": min,
                 "shs": lambda sigma, total: 1 / (1 + math.exp(-sigma * total))}
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return node.value
    if isinstance(node, ast.Name):
        return values[node.id]
    if isinstance(node, ast.BinOp) and type(node.op) in operations:
        return operations[type(node.op)](_source_expression(node.left, values), _source_expression(node.right, values))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_source_expression(node.operand, values)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in functions and not node.keywords:
        return functions[node.func.id](*[_source_expression(arg, values) for arg in node.args])
    raise AssertionError(f"unexpected source expression: {ast.dump(node)}")


def test_source_rhs_and_raw_roadrunner_rates_agree_on_all_species_and_min_branches(inventory, nutrient_model):
    roadrunner = pytest.importorskip("roadrunner")
    reference = roadrunner.RoadRunner(_source(inventory, "jalihal2021", "nutrient-signaling-sbml.xml").decode())
    simulator = ast.parse(_source(inventory, "jalihal2021", "simulators.py.source"))
    definitions = next(ast.literal_eval(n.value) for n in ast.walk(simulator)
                       if isinstance(n, ast.keyword) and n.arg == "fnspecs")
    assert definitions["shs"] == (["sig", "summation"], "1/(1+e^(-sig*summation))")
    assert definitions["tRNA"][1] == "min(tRNA_tot, AmAc)"
    assert definitions["pRib"][1] == "min(rib_comp,init_factor)"
    expressions = {name: ast.parse(expr, mode="eval").body for name, expr in
                   _table(_source(inventory, "jalihal2021", "variables.txt")).items()}
    assert tuple(expressions) == nutrient_model.species_ids
    rng = np.random.default_rng(31271)
    states = [nutrient_model.initial_state(), *rng.uniform(0, 1, (12, 25))]
    for glutamine in (0, 0.5, nutrient_model.parameter_values["tRNA_total"], 3.0):
        states.append(nutrient_model.initial_state({"Glutamine": glutamine, "Rib": 0.5, "eIF": 0.5, "Gcn2": 0.5}))
    reaction_order = [list(reference.model.getReactionIds()).index(r) for r in nutrient_model.reaction_ids]
    species_order = [list(reference.model.getFloatingSpeciesIds()).index(s) for s in nutrient_model.species_ids]
    for state in states:
        values = {**nutrient_model.parameter_values, **dict(zip(nutrient_model.species_ids, state, strict=True))}
        source_rhs = [_source_expression(expr, values) for expr in expressions.values()]
        for name, value in zip(nutrient_model.species_ids, state, strict=True):
            reference[name] = float(value)
        np.testing.assert_allclose(nutrient_model.rhs(0, state), source_rhs, rtol=3e-13, atol=3e-13)
        np.testing.assert_allclose(nutrient_model.reaction_rates(0, state), reference.getReactionRates()[reaction_order],
                                   rtol=3e-13, atol=3e-13)
        np.testing.assert_allclose(nutrient_model.rhs(0, state), reference.getRatesOfChange()[species_order],
                                   rtol=3e-13, atol=3e-13)


INPUT_CASES = [
    {"Carbon": 0.5, "ATP": 0.5, "Glutamine_ext": 0.0, "NH4": 0.0, "Proline": 0.0},
    {"Carbon": 1.0, "ATP": 1.0, "Glutamine_ext": 1.0, "NH4": 0.0, "Proline": 0.0},
    {"Carbon": 0.0, "ATP": 0.0, "Glutamine_ext": 1.0, "NH4": 0.0, "Proline": 0.0},
    {"Carbon": 1.0, "ATP": 1.0, "Glutamine_ext": 0.0, "NH4": 1.0, "Proline": 0.0},
    {"Carbon": 1.0, "ATP": 1.0, "Glutamine_ext": 0.0, "NH4": 0.0, "Proline": 1.0},
]


@pytest.mark.parametrize("inputs", INPUT_CASES, ids=["source_default", "glutamine", "carbon_absent", "ammonium", "proline"])
def test_native_trajectories_match_independent_raw_roadrunner_without_time_or_dose_conversion(inventory, nutrient_model, inputs):
    roadrunner = pytest.importorskip("roadrunner")
    times = np.linspace(0, 60, 121)
    trajectory = simulate_native_reference("jalihal2021", times, normalized_inputs=inputs, inventory=inventory, root=ROOT)
    reference = roadrunner.RoadRunner(_source(inventory, "jalihal2021", "nutrient-signaling-sbml.xml").decode())
    # Rich-input late-time oscillations amplify solver error; these settings were
    # refined against both Radau and tighter BDF, rather than widening agreement bounds.
    reference.integrator.relative_tolerance = 3e-14
    reference.integrator.absolute_tolerance = 3e-16
    reference.integrator.maximum_num_steps = 100000
    for name, value in inputs.items():
        reference[name] = value
    expected = np.asarray(reference.simulate(0, 60, 121, selections=["time", *nutrient_model.species_ids]))
    np.testing.assert_array_equal(trajectory.times_native, expected[:, 0])
    np.testing.assert_allclose(trajectory.states, expected[:, 1:], rtol=3e-6, atol=1e-7)
    assert np.isfinite(trajectory.states).all() and np.isfinite(trajectory.reaction_rates).all()
    assert {"PKA", "Snf1", "TORC1", "Sch9"} <= trajectory.variables.keys()
    assert not hasattr(trajectory, "times_s")
    assert trajectory.metadata["normalized_inputs"] == inputs
    assert trajectory.metadata["physical_seconds_per_native_unit"] is None
    assert trajectory.metadata["physical_input_calibration"] is None
    assert trajectory.metadata["independent_real_validation_established"] is False
    for name, value in inputs.items():
        np.testing.assert_array_equal(trajectory.variables[name], np.full(times.shape, value))


def test_a_carbon_shift_reuses_the_native_state_without_an_implicit_atp_map(inventory):
    pre = simulate_native_reference("jalihal2021", [0, 60], normalized_inputs=INPUT_CASES[1], inventory=inventory, root=ROOT)
    inputs = {**INPUT_CASES[1], "Carbon": 0.0}  # ATP deliberately remains 1, not silently coupled to Carbon.
    post = simulate_native_reference("jalihal2021", [0, 0.1, 1, 10], normalized_inputs=inputs,
                                     initial_state=pre.states[-1], inventory=inventory, root=ROOT)
    np.testing.assert_array_equal(post.states[0], pre.states[-1])
    assert np.all(post.variables["ATP"] == 1.0)
    assert np.all(post.variables["Carbon"] == 0.0)
    assert not np.array_equal(post.states[-1], post.states[0])


@pytest.mark.parametrize("patch", [{"Carbon": 20.0}, {"NH4": -1.0}, {"ATP": True}, {"Proline": float("nan")},
                                    {"Glutamine_ext": "1"}, {"glucose_g_L": 2.0}])
def test_physical_or_implicit_nutrient_inputs_are_refused(inventory, patch):
    with pytest.raises(ValueError, match="normalized"):
        simulate_native_reference("jalihal2021", [0, 1], normalized_inputs={**INPUT_CASES[0], **patch}, inventory=inventory)
    with pytest.raises(ValueError, match="all five"):
        simulate_native_reference("jalihal2021", [0, 1], normalized_inputs={"Carbon": 0.5}, inventory=inventory)


def test_source_seconds_minutes_discrepancy_is_preserved_from_authoritative_code(inventory, nutrient_model):
    # The archived exporter uses legacy non-raw regex literals. Parse its AST under
    # the historical warning policy; do not execute it or change its source bytes.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        translator = ast.parse(_source(inventory, "jalihal2021", "sbml-translator.py.source"))
    declarations = [n for n in ast.walk(translator) if isinstance(n, ast.Call)
                    and isinstance(n.func, ast.Attribute) and n.func.attr == "setTimeUnits"]
    assert len(declarations) == 1 and ast.literal_eval(declarations[0].args[0]) == "second"
    code = ast.parse(_source(inventory, "jalihal2021", "timecourse.py.source"))
    simulate = next(n for n in ast.walk(code) if isinstance(n, ast.FunctionDef) and n.name == "simulate")
    expected = ast.dump(ast.parse("tmax = tmax / 60.").body[0], include_attributes=False)
    assert any(ast.dump(n, include_attributes=False) == expected for n in ast.walk(simulate))
    audit = nutrient_model.metadata["native_reference"]["unit_audit"]
    assert audit["time"]["status"] == "conflicting_declarations"
    assert audit["states"]["absolute_units_verified"] is False
    assert set(nutrient_model.metadata["units"]["parameters"].values()) == {"unspecified"}
    assert nutrient_model.metadata["units"]["physical_calibration"] == {"time": None, "states": None, "nutrient_inputs": None}
    changed = inventory.to_dict()
    changed["models"]["jalihal2021"]["unit_audit"]["time"]["physical_seconds_per_native_unit"] = 60
    with pytest.raises(ValueError, match="unresolved"):
        validate_inventory(changed)
    changed = inventory.to_dict()
    changed["models"]["jalihal2021"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="immutable"):
        validate_inventory(changed)


def test_primary_article_resolves_zheng_exact_archive_and_cc0_license(inventory):
    article = json.loads(_source(inventory, "zheng2016", "elife-18638-v3.json"))
    assert article["version"] == 3 and article["doi"] == "10.7554/eLife.18638"
    code = next(a for a in article["additionalFiles"] if a["doi"] == "10.7554/eLife.18638.025")
    assert code["uri"] == "https://cdn.elifesciences.org/articles/18638/elife-18638-code1-v3.zip"
    assert article["copyright"]["license"] == "CC0-1.0"
    primary = ET.fromstring(_source(inventory, "zheng2016", "PMC5127643.xml"))
    media = primary.find(".//supplementary-material[@id='SD7-data']/media")
    assert media.get("{http://www.w3.org/1999/xlink}href") == "elife-18638-code1.zip"
    archive = _source(inventory, "zheng2016", "elife-18638-code1-v3.zip")
    with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
        names = [name for name in bundle.namelist() if "/" not in name and name.endswith(".m")]
        assert len(names) == 7
        for name in names:
            assert bundle.read(name) == _source(inventory, "zheng2016", name)
        # The binary resource forks are retained in the archive, never decoded or extracted.
        assert "__MACOSX/._Fig1_IPMS.m" in bundle.namelist()
    assert not (ROOT / ASSETS / "zheng2016/__MACOSX").exists()
    assert inventory.to_dict()["models"]["zheng2016"]["availability"] == "local_source_only"
    with pytest.raises(EvidenceGap, match="MATLAB"):
        load_reference_model("zheng2016", inventory)


def test_actual_goulev_equations_are_in_primary_text_not_the_strain_supplement(inventory):
    article = ET.fromstring(_source(inventory, "goulev2017", "PMC5438251.xml"))
    for identifier in ("equ6", "equ7", "equ8"):
        assert article.find(f".//disp-formula[@id='{identifier}']") is not None
    for identifier in ("s4-13-2", "s4-13-6", "tblu3", "tblu4"):
        assert article.find(f".//*[@id='{identifier}']") is not None
    api = json.loads(_source(inventory, "goulev2017", "elife-23971-v4.json"))
    assert api["additionalFiles"][0]["doi"] == "10.7554/eLife.23971.026"
    assert api["additionalFiles"][0]["uri"] == "https://cdn.elifesciences.org/articles/23971/elife-23971-supp1-v4.docx"
    with zipfile.ZipFile(io.BytesIO(_source(inventory, "goulev2017", "elife-23971-supp1-v4.docx"))) as bundle:
        text = " ".join(ET.fromstring(bundle.read("word/document.xml")).itertext())
    assert "Genotype" in text and "Euroscarf" in text
    assert inventory.to_dict()["models"]["goulev2017"]["availability"] == "equations_only"


def test_pincus_actual_text_s1_preserves_count_units_and_fitted_parameter_limitations(inventory):
    payload = _source(inventory, "pincus2010", "pbio.1000415.s016.pdf")
    assert payload.startswith(b"%PDF")
    pdftotext = shutil.which("pdftotext")
    if pdftotext is None:
        # The byte audit is still required; no equation execution is claimed.
        assert inventory.to_dict()["models"]["pincus2010"]["availability"] == "equations_only"
        return
    result = subprocess.run([pdftotext, "-layout", "-", "-"], input=payload, capture_output=True, check=True, timeout=30)
    text = result.stdout.decode()
    assert "number of molecules rather than concentrations" in text
    assert "There are seven such parameters" in text
    assert "molecular number of DTT within the ER" in text
    assert "ODE15s" in text and "mol−2 s−1" in text
    record = inventory.to_dict()["models"]["pincus2010"]
    assert record["availability"] == "equations_only"
    assert "ambiguities" in record["gap"]


def test_krakowiak_republishes_the_zheng_model_file_and_a_different_transcription_gain(inventory):
    # The same bytes under two DOIs are one source, not two agreeing sources.
    assert _source(inventory, "krakowiak2018", "elife-31668-code1-v3.m") == _source(
        inventory, "zheng2016", "titration_YFP_FB.m")
    assert "titration_YFP_noFB" in _source(inventory, "krakowiak2018", "elife-31668-code2-v3.m").decode()
    driver = _source(inventory, "krakowiak2018", "elife-31668-code4-v3.m").decode("utf-8", "replace")
    assert "beta = 1.7783/5" in driver and "UPo_39 = 10.5141" in driver
    article = json.loads(_source(inventory, "krakowiak2018", "elife-31668-v3.json"))
    assert article["version"] == 3 and article["copyright"]["license"] == "CC-BY-4.0"
    primary = _source(inventory, "krakowiak2018", "PMC5809143.xml").decode("utf-8", "replace")
    for printed in ("Previous Paper model values", "1.778", "0.3557", "10.51"):
        assert printed in primary
    assert inventory.to_dict()["models"]["krakowiak2018"]["availability"] == "local_source_only"
    with pytest.raises(EvidenceGap, match="byte-identical"):
        load_reference_model("krakowiak2018", inventory)


def test_ke2013_prints_a_complete_ode_system_with_units_and_no_buffering_term(inventory):
    record = inventory.to_dict()["models"]["ke2013"]
    assert _source(inventory, "ke2013", "pcbi.1002879.s009.pdf").startswith(b"%PDF")
    pdftotext = shutil.which("pdftotext")
    if pdftotext is None:
        # The byte audit still holds; no equation transcription is claimed either way.
        assert record["availability"] == "equations_only"
        return

    def text(name):
        payload = _source(inventory, "ke2013", name)
        result = subprocess.run([pdftotext, "-layout", "-", "-"], input=payload,
                                capture_output=True, check=True, timeout=30)
        return result.stdout.decode()

    system, transporters, physiological = (
        text(f"pcbi.1002879.s{part}.pdf") for part in ("009", "010", "012"))
    assert "ODE system for the integrative model" in system and "J Pma1" in system
    assert "kPma1" in transporters and "1e-18 mol/(s*mM2)" in transporters
    assert "adjusted during the model integration step" in transporters
    assert "Intracellular negative" in physiological
    # Free protons, no buffering capacity: the absent term must stay refused, not estimated.
    for supplement in ("009", "012", "015"):
        assert "buffer" not in text(f"pcbi.1002879.s{supplement}.pdf").lower()
    assert record["availability"] == "equations_only"
    with pytest.raises(EvidenceGap, match="no executable artifact"):
        load_reference_model("ke2013", inventory)


def test_zheng_native_runtime_or_explicit_source_only_status(inventory, tmp_path):
    octave = shutil.which("octave")
    if octave is None:
        # Absence is a source-only result, not a silently successful MATLAB replay.
        assert inventory.to_dict()["models"]["zheng2016"]["availability"] == "local_source_only"
        with pytest.raises(EvidenceGap, match="MATLAB"):
            load_reference_model("zheng2016", inventory)
        return
    for name in ("titration_YFP_FB.m", "titration_YFP_FB_decoy.m", "titration_YFP_FB_Phosph.m"):
        _source(inventory, "zheng2016", name)
    folder = str(ROOT / ASSETS / "zheng2016").replace("'", "''")
    code = f"""
addpath('{folder}');
k1=166.8; k2=2.783; k3=k1; k4=0.0464; k5=4.642e-7; beta=1.7783; Kd=0.0022;
y=[0.8;0.0005;3;0.0015;0.25;4];
r=titration_YFP_FB(0,y,k1,k2,k3,k4,k5,beta,Kd);
rd=titration_YFP_FB_decoy(0,[y;0.002;0.003],k1,k2,k3,k4,k5,beta,Kd);
rp=titration_YFP_FB_Phosph(30,y,k1,k2,k3,k4,k5,[beta,2*beta,beta],Kd,[0,30,60]);
times=linspace(0,60,61); cold=[1;0;0;1/500;0;0]; hot=[1;0;10;1/500;0;0];
[t,c]=ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd),times,cold);
[t,h]=ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd),times,hot);
[t,d]=ode23s(@(t,y)titration_YFP_FB_decoy(t,y,k1,k2,k3,k4,k5,beta,Kd),times,[hot;0.002;0.003]);
assert(all(isfinite([c(:);h(:);d(:)])));
assert(max(abs(h(:,2)+h(:,4)-1/500))<1e-10);
assert(max(abs(d(:,7)+d(:,8)-0.005))<1e-10);
result=struct('runtime',version,'rate',r,'decoy_rate',rd,'phosph_rate',rp,'cold_final',c(end,:),'hot_final',h(end,:),'min_hot_complex',min(h(:,4)),'samples',length(t));
printf('NATIVE_RESULT:%s\\n',jsonencode(result));
"""
    result = subprocess.run([octave, "--no-gui", "--quiet", "--no-history", "--no-init-file", "--no-site-file", "--eval", code],
                            cwd=tmp_path, capture_output=True, text=True, check=True, timeout=120)
    data = json.loads(next(line.removeprefix("NATIVE_RESULT:") for line in result.stdout.splitlines()
                           if line.startswith("NATIVE_RESULT:")))
    rate, decoy, phosph = (np.asarray(data[key]).ravel() for key in ("rate", "decoy_rate", "phosph_rate"))
    assert data["runtime"] and data["samples"] == 61
    assert rate[1] + rate[3] == pytest.approx(0, abs=1e-12)
    assert decoy[6] + decoy[7] == pytest.approx(0, abs=1e-12)
    # Preserve the actual source's double k5 contribution; do not silently "correct" it.
    assert rate[2] + rate[4] == pytest.approx(-2 * 4.642e-7 * 0.25, abs=1e-11)
    assert rate[0] + rate[3] + rate[4] == pytest.approx(rate[5], abs=1e-11)
    np.testing.assert_allclose(phosph[[1, 2, 3, 4]], rate[[1, 2, 3, 4]], rtol=1e-12)
    assert phosph[5] == pytest.approx(2 * rate[5])
    assert data["min_hot_complex"] < data["cold_final"][3]
    assert data["hot_final"][5] > data["cold_final"][5]


def test_fetch_verify_only_and_existing_assets_never_use_network_or_restamp(inventory, fetcher, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("existing pinned bytes must not trigger a network request")
    monkeypatch.setattr(fetcher, "download_https", forbidden)
    monkeypatch.setattr(fetcher, "download_github", forbidden)
    entries = fetcher.asset_entries(inventory)
    before = {e["path"]: (ROOT / e["path"]).stat().st_mtime_ns for e in entries}
    for fetch in (False, True):
        report = fetcher.acquire(inventory, root=ROOT, fetch=fetch)
        assert report["passed"]
        assert all(row["status"] == "verified_existing" for row in report["assets"])
        assert report["independent_real_validation_established"] is False
    assert before == {e["path"]: (ROOT / e["path"]).stat().st_mtime_ns for e in entries}


@pytest.mark.parametrize("selected", [[], ["garcia_salcedo2014"], ["not_a_registered_source"]])
def test_empty_or_unregistered_source_selection_does_not_claim_verification(inventory, fetcher, selected):
    with pytest.raises(ValueError):
        fetcher.asset_entries(inventory, selected)


def test_acquisition_rejects_bad_responses_before_creating_any_file(inventory, fetcher, monkeypatch, tmp_path):
    monkeypatch.setattr(fetcher, "download_https", lambda url, size: b"x" * size)
    report = fetcher.acquire(inventory, root=tmp_path, source_ids=["pincus2010"], fetch=True)
    assert not report["passed"] and len(report["assets"]) == 1
    assert "SHA-256 mismatch" in report["assets"][0]["error"]
    assert not (tmp_path / ASSETS).exists()


def test_acquisition_creates_only_pinned_paths_and_never_replaces_conflicting_bytes(inventory, fetcher, monkeypatch, tmp_path):
    entries = fetcher.asset_entries(inventory, ["pincus2010"])
    payloads = {e["acquisition"]["url"]: (ROOT / e["path"]).read_bytes() for e in entries}
    monkeypatch.setattr(fetcher, "download_https", lambda url, size: payloads[url])
    report = fetcher.acquire(inventory, root=tmp_path, source_ids=["pincus2010"], fetch=True)
    assert report["passed"]
    path = tmp_path / entries[0]["path"]
    path.write_bytes(b"do not overwrite this conflicting local file")
    report = fetcher.acquire(inventory, root=tmp_path, source_ids=["pincus2010"], fetch=True)
    assert not report["passed"]
    assert path.read_bytes() == b"do not overwrite this conflicting local file"


@pytest.mark.parametrize("relative", ["../escape", "data/hog2013/model_wt.xml", "data/native_reference_models/../old",
                                      "data/native_reference_models", "/absolute", "data//native_reference_models/x"])
def test_source_fetch_destinations_cannot_escape_the_owned_new_asset_directory(fetcher, tmp_path, relative):
    with pytest.raises(ValueError):
        fetcher.destination(tmp_path, relative)


def test_source_fetch_rejects_symlinks_and_https_downgrades(fetcher, tmp_path):
    (tmp_path / "data").symlink_to(ROOT / "data", target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        fetcher.destination(tmp_path, f"{ASSETS}/test.xml")
    with pytest.raises(ValueError, match="HTTPS"):
        fetcher._url("http://example.org/source.xml")
    with pytest.raises(ValueError, match="credential-free"):
        fetcher._url("https://user:password@example.org/source.xml")


def test_jalihal_fetch_retains_remote_paths_and_stores_nonimportable_source_data(inventory, fetcher, monkeypatch, tmp_path):
    entries = fetcher.asset_entries(inventory, ["jalihal2021"])
    payloads = {entry["path"]: (ROOT / entry["path"]).read_bytes() for entry in entries}
    github = {
        f"repos/{entry['repository']}/contents/{entry['acquisition']['path']}?ref={entry['commit']}": entry
        for entry in entries if entry["acquisition"]["kind"] == "github"
    }
    https = {entry["acquisition"]["url"]: entry for entry in entries if entry["acquisition"]["kind"] == "https"}
    requested = []

    def run(command, **kwargs):
        assert command[:4] == ["gh", "api", "--method", "GET"]
        entry = github[command[4]]
        requested.append(command[4])
        response = {
            "type": "file", "path": entry["acquisition"]["path"], "size": entry["bytes"],
            "sha": entry["acquisition"]["git_blob"], "encoding": "base64",
            "content": base64.b64encode(payloads[entry["path"]]).decode(),
        }
        return SimpleNamespace(stdout=json.dumps(response).encode())

    def download_https(url, size):
        entry = https[url]
        assert size == entry["bytes"]
        return payloads[entry["path"]]

    monkeypatch.setattr(fetcher.subprocess, "run", run)
    monkeypatch.setattr(fetcher, "download_https", download_https)
    report = fetcher.acquire(inventory, root=tmp_path, source_ids=["jalihal2021"], fetch=True)
    assert report["passed"] and len(report["assets"]) == 8
    assert all(row["status"] == "acquired_verified" for row in report["assets"])
    assert requested == list(github)
    for entry in entries:
        assert (tmp_path / entry["path"]).read_bytes() == payloads[entry["path"]]
    for name in ("sbml-translator.py", "simulators.py", "timecourse.py"):
        original = tmp_path / ASSETS / "jalihal2021" / name
        local = original.with_suffix(".py.source")
        assert local.is_file() and not original.exists()
        assert importlib.util.spec_from_file_location("jalihal_upstream_source", local) is None


def test_github_fetch_uses_gh_and_verifies_commit_path_blob_before_accepting_payload(inventory, fetcher, monkeypatch):
    entry = next(e for e in fetcher.asset_entries(inventory) if e["path"].endswith("nutrient-signaling-sbml.xml"))
    payload = (ROOT / entry["path"]).read_bytes()
    response = {"type": "file", "path": entry["acquisition"]["path"], "size": len(payload),
                "sha": entry["acquisition"]["git_blob"], "encoding": "base64", "content": base64.b64encode(payload).decode()}
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(stdout=json.dumps(response).encode())
    monkeypatch.setattr(fetcher.subprocess, "run", run)
    assert fetcher.download_github(entry) == payload
    assert commands[0][:4] == ["gh", "api", "--method", "GET"]
    assert commands[0][-1].endswith("?ref=3e93dba11d4e7950ff89519e40f32ba5766e7f16")
    response["sha"] = "0" * 40
    with pytest.raises(fetcher.AcquisitionError, match="identity"):
        fetcher.download_github(entry)
    bad_blob = {**entry, "acquisition": {**entry["acquisition"], "git_blob": "0" * 40}}
    with pytest.raises(ValueError, match="Git blob"):
        fetcher.verify_payload(payload, bad_blob)
