"""Offline source identity/content checks, not a COSMIC or yeast biology validation.

These tests read the retained HTTP bodies, never the network or private inputs. The
MATLAB kernel is inert text: it is neither imported nor executed. An unreviewed
CSV/workbook/archive/model needs its own source identity and schema before admission;
renaming an HTML response or accepting HTTP 200 alone is not source acquisition.
"""

from __future__ import annotations

import base64
import hashlib
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from datetime import datetime
from pathlib import Path, PurePosixPath
from urllib.parse import parse_qs, urlsplit

import pytest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/cosmic_reference"
DOI = "10.1016/j.ymben.2024.02.012"
PREPRINT_DOI = "10.1101/2023.09.13.557646"
PMID = "38387677"
TITLE = "COSMIC-dFBA: A novel multi-scale hybrid framework for bioprocess modeling"
COMMIT = "3d60d3a2ae67943d770d7f8fbbb8bcc16f391a86"
KERNEL_SHA256 = "b82441d85b96c044e82a94f3bc27ab9e1c339ba9e5db0183133a8e7d21dd9b10"
ELSEVIER = "http://www.elsevier.com/xml/svapi/article/dtd"
PRISM = "http://prismstandard.org/namespaces/basic/2.0/"
XML_ROOTS = {
    "pubmed": "PubmedArticleSet",
    "pubmed_pmc": "eLinkResult",
    "pubmed_linkout": "eLinkResult",
    "elsevier_metadata": f"{{{ELSEVIER}}}full-text-retrieval-response",
}
JSON_IDS = {
    "europepmc_journal",
    "pmc_idconv",
    "crossref",
    "europepmc_preprint",
    "biorxiv_details",
    "openalex_journal",
    "openalex_preprint",
    "datacite_related",
    "github_tree",
    "github_branches",
    "github_releases",
    "github_tags",
    "github_issues",
}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _invalid_constant(value):
    raise ValueError(f"non-JSON constant: {value}")


def _json(body):
    return json.loads(body, object_pairs_hook=_unique_object, parse_constant=_invalid_constant)


@pytest.fixture(scope="module")
def manifest():
    return _json((ROOT / "data/cosmic_sources.json").read_bytes())


@pytest.fixture(scope="module")
def responses(manifest):
    return {row["id"]: row for row in manifest["acquisition_audit"]["retained_responses"]}


def _reference_bytes(relative_path):
    relative = PurePosixPath(relative_path)
    assert not relative.is_absolute()
    assert relative.parts[:2] == ("data", "cosmic_reference")
    assert len(relative.parts) == 3 and ".." not in relative.parts
    path = ROOT / relative
    assert not path.is_symlink() and path.resolve().parent == REFERENCE.resolve()
    return path.read_bytes()


def _body(record):
    stored = _reference_bytes(record["path"])
    if "storage_encoding" not in record:
        return stored
    assert record["storage_encoding"] == "base64"
    assert len(stored) == record["stored_bytes"], record["id"]
    assert hashlib.sha256(stored).hexdigest() == record["stored_sha256"], record["id"]
    original = base64.b64decode(b"".join(stored.splitlines()), validate=True)
    assert base64.encodebytes(original) == stored
    copy = record["readable_copy"]
    readable = _reference_bytes(copy["path"])
    assert copy["transformation"] == "CRLF_to_LF_only"
    assert readable == original.replace(b"\r\n", b"\n")
    assert len(readable) == copy["bytes"]
    assert hashlib.sha256(readable).hexdigest() == copy["sha256"]
    return original


def _parse_recorded_asset(record, body):
    """Fail closed on this audited corpus, not a generic table/archive downloader."""
    if record["http_status"] != 200:
        raise ValueError("HTTP error is not a scientific asset")
    media_type = record["content_type"].partition(";")[0].strip().lower()
    format_name = record["format"]
    media_types = {
        "json": {"application/json"},
        "xml": {"application/xml", "text/xml"},
        "matlab_text": {"application/vnd.github.raw+json"},
        "text": {"application/vnd.github.raw+json"},
    }
    if format_name not in media_types:
        raise ValueError("unreviewed source format requires a source-specific schema")
    if media_type not in media_types[format_name]:
        raise ValueError("response media type does not match the reviewed source format")
    if format_name == "json":
        assert record["id"] in JSON_IDS
        parsed = _json(body)
        assert isinstance(parsed, (dict, list))
        return parsed
    if format_name == "xml":
        parsed = ET.fromstring(body)
        if parsed.tag != XML_ROOTS[record["id"]]:
            raise ValueError("wrong XML root: an HTML/error page is not the declared source")
        return parsed
    text = body.decode("utf-8")
    if format_name == "matlab_text":
        assert record["id"] == "github_kernel"
        assert text.startswith("function [dFBA_results] = COSMIC_dFBA(")
    else:
        assert record["id"] == "github_license"
        assert text.startswith("Creative Commons Legal Code\n\nCC0 1.0 Universal\n")
    return text


def _parsed(responses, name):
    record = responses[name]
    return _parse_recorded_asset(record, _body(record))


def _check_primary_identity(pubmed, epmc, crossref):
    assert len(pubmed.findall("PubmedArticle")) == 1
    citation = pubmed.find("PubmedArticle/MedlineCitation")
    assert citation.findtext("PMID") == PMID
    assert citation.findtext("Article/ArticleTitle") == TITLE + "."
    ids = {
        node.attrib["IdType"]: node.text
        for node in pubmed.findall("PubmedArticle/PubmedData/ArticleIdList/ArticleId")
    }
    assert ids == {"pubmed": PMID, "doi": DOI, "pii": "S1096-7176(24)00028-4"}
    assert epmc["request"]["queryString"] == f"EXT_ID:{PMID} AND SRC:MED"
    assert epmc["hitCount"] == 1 and len(epmc["resultList"]["result"]) == 1
    row = epmc["resultList"]["result"][0]
    assert (row["source"], row["id"], row["pmid"], row["doi"]) == ("MED", PMID, PMID, DOI)
    assert row["title"] == TITLE + "."
    work = crossref["message"]
    assert work["DOI"] == DOI and work["type"] == "journal-article"
    assert work["title"] == [TITLE]
    assert work["volume"] == "82" and work["page"] == "183-192"


def _pmc_identity_links(root):
    return [
        node.text
        for linkset in root.findall("LinkSet/LinkSetDb")
        if linkset.findtext("LinkName") == "pubmed_pmc" and linkset.findtext("DbTo") == "pmc"
        for node in linkset.findall("Link/Id")
    ]


def test_retained_corpus_is_complete_with_no_importable_source(manifest, responses):
    rows = manifest["acquisition_audit"]["retained_responses"]
    assert len(rows) == len(responses) == 19
    assert set(responses) == JSON_IDS | set(XML_ROOTS) | {"github_kernel", "github_license"}
    expected_paths = {ROOT / row["path"] for row in rows}
    assert len(expected_paths) == len(rows)
    copies = {ROOT / row["readable_copy"]["path"] for row in rows if "readable_copy" in row}
    assert len(copies) == 2 and expected_paths.isdisjoint(copies)
    expected_paths |= copies
    assert {p for p in REFERENCE.rglob("*") if p.is_file()} == expected_paths
    assert all(path.suffix in {".txt", ".json", ".xml", ".base64"} for path in expected_paths)


@pytest.mark.parametrize(
    "record_id", sorted(JSON_IDS | set(XML_ROOTS) | {"github_kernel", "github_license"})
)
def test_each_retained_response_has_exact_bytes_format_and_provenance(
    manifest, responses, record_id
):
    row = responses[record_id]
    body = _body(row)
    assert len(body) == row["bytes"]
    assert hashlib.sha256(body).hexdigest() == row["sha256"]
    assert row["method"] == "GET" and row["source_version"] and row["discovery"]
    assert row["license_ref"] in manifest["acquisition_audit"]["licenses"]
    assert datetime.fromisoformat(row["retrieved_at_utc"]).utcoffset().total_seconds() == 0
    url = urlsplit(row["url"])
    assert url.scheme == "https" and url.username is None and url.password is None
    assert row["client"] == ("gh api" if url.hostname == "api.github.com" else "curl")
    _parse_recorded_asset(row, body)


def test_exact_primary_identity_not_a_cited_work_or_similar_title(manifest, responses):
    assert (manifest["paper"]["doi"], manifest["paper"]["pmid"]) == (DOI, PMID)
    assert manifest["paper"]["title"] == TITLE
    pubmed = _parsed(responses, "pubmed")
    _check_primary_identity(
        pubmed, _parsed(responses, "europepmc_journal"), _parsed(responses, "crossref")
    )
    authors = pubmed.findall("PubmedArticle/MedlineCitation/Article/AuthorList/Author")
    assert len(authors) == 12
    assert authors[0].findtext("LastName") == "Gopalakrishnan"
    assert authors[-1].findtext("LastName") == "Lewis"


@pytest.mark.parametrize(
    "field,value", [("id", "wrong-id"), ("source", "PPR"), ("doi", "wrong-doi")]
)
def test_wrong_structured_identity_is_rejected_even_with_the_same_title(responses, field, value):
    epmc = deepcopy(_parsed(responses, "europepmc_journal"))
    epmc["resultList"]["result"][0][field] = value
    with pytest.raises(AssertionError):
        _check_primary_identity(_parsed(responses, "pubmed"), epmc, _parsed(responses, "crossref"))


def test_doi_in_reference_list_does_not_identify_the_primary_work(responses):
    crossref = deepcopy(_parsed(responses, "crossref"))
    crossref["message"]["DOI"] = "not-the-requested-primary-work"
    crossref["message"]["reference"].append({"DOI": DOI})
    with pytest.raises(AssertionError):
        _check_primary_identity(
            _parsed(responses, "pubmed"), _parsed(responses, "europepmc_journal"), crossref
        )


def test_actual_pmc_metadata_has_no_identity_link_or_hosted_supplement(manifest, responses):
    audit = manifest["paper"]["open_access_audit"]
    assert manifest["paper"]["pmcid"] is None
    for name in ("journal", "preprint"):
        row = _parsed(responses, "europepmc_" + name)["resultList"]["result"][0]
        assert "pmcid" not in row
        for flag, value in audit[f"europepmc_{name}_flags"].items():
            assert row[flag] == value == "N"
    converted = _parsed(responses, "pmc_idconv")
    assert converted["status"] == "ok"  # API success is not record resolution.
    assert converted["records"] == [
        {
            "pmid": int(PMID),
            "requested-id": PMID,
            "status": "error",
            "errmsg": "Identifier not found in PMC",
        }
    ]
    root = _parsed(responses, "pubmed_pmc")
    assert root.findtext("LinkSet/IdList/Id") == PMID
    assert _pmc_identity_links(root) == []
    assert parse_qs(urlsplit(responses["pubmed_pmc"]["url"]).query)["linkname"] == ["pubmed_pmc"]
    assert audit["pmc_identity_link_count"] == 0
    linkout = _parsed(responses, "pubmed_linkout").findall("LinkSet/IdUrlList/IdUrlSet/ObjUrl")
    assert len(linkout) == 1
    assert linkout[0].findtext("Attribute") == "subscription/membership/fee required"


def test_citing_pmc_articles_are_not_the_papers_pmc_identity(responses):
    root = deepcopy(_parsed(responses, "pubmed_pmc"))
    linkset = ET.SubElement(root.find("LinkSet"), "LinkSetDb")
    ET.SubElement(linkset, "DbTo").text = "pmc"
    ET.SubElement(linkset, "LinkName").text = "pubmed_pmc_refs"
    ET.SubElement(ET.SubElement(linkset, "Link"), "Id").text = "synthetic-citing-id"
    assert _pmc_identity_links(root) == []


def test_preprint_version_and_requested_links_are_declared_not_guessed(manifest, responses):
    work = _parsed(responses, "crossref")["message"]
    preprint = _parsed(responses, "biorxiv_details")["collection"]
    assert len(preprint) == 1
    row = preprint[0]
    assert (row["doi"], row["published"], row["version"], row["date"]) == (
        PREPRINT_DOI,
        DOI,
        "1",
        "2023-09-17",
    )
    assert work["relation"]["has-preprint"] == [
        {"id-type": "doi", "id": PREPRINT_DOI, "asserted-by": "object"}
    ]
    epmc_preprint = _parsed(responses, "europepmc_preprint")["resultList"]["result"][0]
    assert (epmc_preprint["source"], epmc_preprint["id"], epmc_preprint["doi"]) == (
        "PPR",
        "PPR726164",
        PREPRINT_DOI,
    )
    relation = epmc_preprint["commentCorrectionList"]["commentCorrection"][0]
    assert (relation["type"], relation["id"], relation["source"]) == ("Preprint of", PMID, "MED")
    attempts = {row["id"]: row for row in manifest["acquisition_audit"]["noncontent_responses"]}
    assert attempts["biorxiv_jats"]["url"] == row["jatsxml"]
    assert attempts["publisher_linkinghub"]["url"] == work["resource"]["primary"]["URL"]
    declared = {link["URL"] for link in work["link"]}
    assert responses["elsevier_metadata"]["url"] in declared
    assert attempts["publisher_plain"]["url"] in declared
    for name in ("europepmc_journal", "europepmc_preprint"):
        links = _parsed(responses, name)["resultList"]["result"][0]["fullTextUrlList"][
            "fullTextUrl"
        ]
        assert attempts["biorxiv_pdf"]["url"] in {
            link["url"] for link in links if link["documentStyle"] == "pdf"
        }


def test_publisher_http_200_xml_is_only_coredata_not_full_text(responses):
    root = _parsed(responses, "elsevier_metadata")
    assert [child.tag for child in root] == [f"{{{ELSEVIER}}}coredata"]
    assert root.findtext(f"{{{ELSEVIER}}}coredata/{{{PRISM}}}doi") == DOI
    assert root.findtext(f"{{{ELSEVIER}}}coredata/{{{ELSEVIER}}}openaccess") == "0"
    assert root.findtext(f"{{{ELSEVIER}}}coredata/{{{ELSEVIER}}}openaccessArticle") == "false"
    assert responses["elsevier_metadata"]["role"] == "publisher_coredata_only"


def test_failed_routes_are_logged_with_actual_status_and_not_saved_as_scientific_files(manifest):
    attempts = {row["id"]: row for row in manifest["acquisition_audit"]["noncontent_responses"]}
    assert {name: row["http_status"] for name, row in attempts.items()} == {
        "publisher_linkinghub": 200,
        "publisher_sciencedirect": 403,
        "publisher_plain": 400,
        "publisher_full_view": 401,
        "biorxiv_jats": 403,
        "biorxiv_pdf": 403,
    }
    for name, row in attempts.items():
        assert row["path"] is None and row["method"] == "GET"
        assert row["bytes"] > 0 and len(bytes.fromhex(row["sha256"])) == 32
        if "body_verbatim" in row:
            body = row["body_verbatim"].encode("utf-8")
            assert len(body) == row["bytes"]
            assert hashlib.sha256(body).hexdigest() == row["sha256"]
            error = ET.fromstring(body)
            assert error.tag == "service-error"
            expected = "AUTHENTICATION_ERROR" if name == "publisher_full_view" else "INVALID_INPUT"
            assert error.findtext("status/statusCode") == expected
        else:
            assert row["content_type"].partition(";")[0] == "text/html"


@pytest.mark.parametrize("format_name", ["csv", "xlsx", "zip", "tar.gz", "sbml"])
def test_unreviewed_formats_cannot_smuggle_html_as_source_data(format_name):
    fake = {"http_status": 200, "format": format_name, "content_type": "application/octet-stream"}
    with pytest.raises(ValueError, match="unreviewed source format"):
        _parse_recorded_asset(fake, b"<html><body>download blocked</body></html>")


@pytest.mark.parametrize(
    "status,content_type", [(403, "text/xml"), (200, "text/html"), (200, "text/xml")]
)
def test_status_mime_and_xml_root_all_matter_for_html_rejection(responses, status, content_type):
    row = {**responses["pubmed"], "http_status": status, "content_type": content_type}
    with pytest.raises(ValueError):
        _parse_recorded_asset(row, b"<html><body>not PubMed XML</body></html>")


def test_json_rejects_html_duplicate_identity_keys_and_nonfinite_values(responses):
    row = responses["europepmc_journal"]
    for body in (b"<html>not JSON</html>", b'{"doi":"wrong","doi":"right"}', b'{"value":NaN}'):
        with pytest.raises(ValueError):
            _parse_recorded_asset(row, body)


def test_public_tree_blob_identities_license_and_source_version(manifest, responses):
    implementation = manifest["public_implementation"]
    assert implementation["commit"] == COMMIT
    tree = _parsed(responses, "github_tree")
    assert tree["sha"] == COMMIT and tree["truncated"] is False
    entries = {entry["path"]: entry for entry in tree["tree"]}
    assert set(entries) == {"COSMIC_dFBA.m", "LICENSE"}
    for filename, record_id in (("COSMIC_dFBA.m", "github_kernel"), ("LICENSE", "github_license")):
        body = _body(responses[record_id])
        entry = entries[filename]
        assert entry["type"] == "blob" and entry["size"] == len(body)
        assert (
            hashlib.sha1(b"blob " + str(len(body)).encode() + b"\0" + body).hexdigest()
            == entry["sha"]
        )
        assert parse_qs(urlsplit(responses[record_id]["url"]).query)["ref"] == [COMMIT]
    kernel = _body(responses["github_kernel"])
    assert len(kernel) == implementation["file_bytes"] == 5864
    assert hashlib.sha256(kernel).hexdigest() == implementation["sha256"] == KERNEL_SHA256
    assert kernel.count(b"\r\n") == kernel.count(b"\n") == 205
    assert responses["github_kernel"]["path"] == implementation["local_path"]
    assert implementation["license"]["spdx_id"] == "CC0-1.0"
    assert responses["github_license"]["sha256"] == implementation["license"]["sha256"]
    assert _parsed(responses, "biorxiv_details")["collection"][0]["license"] == "cc_no"
    licenses = manifest["acquisition_audit"]["licenses"]
    assert licenses["preprint_metadata"]["spdx_id"] is None
    assert licenses["journal_metadata"]["spdx_id"] is None


def test_repository_and_doi_discovery_snapshots_do_not_contain_training_assets(responses):
    branches = _parsed(responses, "github_branches")
    assert {branch["name"] for branch in branches} == {"main", "Documentation"}
    assert {branch["commit"]["sha"] for branch in branches} == {COMMIT}
    for record_id in ("github_releases", "github_tags", "github_issues"):
        assert _parsed(responses, record_id) == []
    datacite = _parsed(responses, "datacite_related")
    assert datacite["data"] == [] and datacite["meta"]["total"] == 0
    journal = _parsed(responses, "openalex_journal")
    preprint = _parsed(responses, "openalex_preprint")
    assert journal["doi"] == "https://doi.org/" + DOI
    assert journal["open_access"]["any_repository_has_fulltext"] is False
    assert preprint["doi"] == "https://doi.org/" + PREPRINT_DOI
    assert len(preprint["locations"]) == 1
    pdf = preprint["locations"][0]["pdf_url"]
    assert pdf == preprint["open_access"]["oa_url"]
    assert urlsplit(pdf).hostname == "www.biorxiv.org"


def test_method_audit_is_anchored_to_actual_kernel_lines_not_biological_truth(manifest, responses):
    audit = manifest["source_method_audit"]
    assert audit["verbatim_source_path"] == responses["github_kernel"]["path"]
    assert audit["source_path"] == responses["github_kernel"]["readable_copy"]["path"]
    lines = _body(responses["github_kernel"]).decode("utf-8").splitlines()
    for section in audit.values():
        if isinstance(section, dict):
            assert section["evidence"]
            for evidence in section["evidence"]:
                start, end = evidence["lines"]
                assert 1 <= start <= end <= len(lines)
                assert evidence["literal"] in "\n".join(lines[start - 1 : end])
    assert [
        line.split("=", 1)[-1].strip().split("(")[0]
        for line in lines
        if line.startswith("function ")
    ] == ["COSMIC_dFBA", "simulate_model", "phase_progress", "evaluate_fluxes", "eval_reactor"]
    assert audit["mixture_weighting"]["state_count"] == 2
    assert audit["mixture_weighting"]["cell_to_biomass_conversion_verified"] is False
    assert audit["classifier_features"]["feature_names"] is None
    assert audit["classifier_features"]["fitted_coefficients_available"] is False
    assert audit["uptake_kinetics"]["parameter_fields"] == [
        "kinetic.Vm_growth",
        "kinetic.Vm_prod",
        "kinetic.Km",
    ]
    training = audit["train_fold_parameter_estimation"]
    assert training["number_of_fitted_parameters"] is None
    assert training["original_train_fold_estimation_verified"] is False
    assert training["training_routine_available"] is False
    assert training["od_only_identifiability_established"] is False
    assert len(training["unverified_parameter_families"]) == 12


def test_local_acquisition_records_original_identities_and_duplicate_aliases(manifest):
    local = manifest["local_acquisition"]
    assert local["source"] == "explicit_user_supplied_files"
    assert local["full_text_acquired"] and local["supplementary_methods_acquired"]
    assert local["normalized_process_data_acquired"]
    assert not local["complete_execution_package_available"]
    assert not local["originals_copied_to_repository"]
    assert not local["automatic_redistribution_authorized"]
    assets = local["assets"]
    assert len(assets) == 5
    assert sum(1 + len(asset["aliases"]) for asset in assets.values()) == 7
    expected = {
        "paper": (5670274, "4808da8e085626be136cf86954b62892a39b4798a25a0c7adfe413c92887198b"),
        "supplementary_methods": (1019024, "4ba691b66644f85261ae24748861b172cce7114a20a795a914304bf5fa8720f9"),
        "supplementary_tables": (4156209, "a68dc8bf86b87e6b9d4f732cd3413948bdb567bf3f85fdabe937a202bb7a7daa"),
    }
    for name, (size, checksum) in expected.items():
        assert (assets[name]["bytes"], assets[name]["sha256"]) == (size, checksum)
    for asset in assets.values():
        assert not PurePosixPath(asset["filename"]).is_absolute()
        assert len(bytes.fromhex(asset["sha256"])) == 32
        if asset["aliases"]:
            assert asset["aliases_verified_byte_identical"] is True
    for name in ("github_commit_archive", "github_branch_archive"):
        assert assets[name]["kernel_sha256"] == KERNEL_SHA256
        assert assets[name]["training_routine_in_archive"] is False


def test_local_method_audit_keeps_lda_and_inferred_labels_distinct_from_pca(manifest):
    methods = manifest["paper_and_supplement_audit"]
    assert methods["classifier"]["projection"] == "linear_discriminant_analysis"
    assert methods["exploratory_pca"]["role"] == "computed_flux_visualization"
    assert methods["phase_labels"]["independently_measured"] is False
    assert methods["phase_labels"]["method"] == "two_state_variance_weighted_nonlinear_regression"
    assert methods["training_partition"]["documented_heldout_culture_split"] is False
    assert methods["reproduction"]["exact_source_model_available"] is False
    assert methods["reproduction"]["fitted_classifier_coefficients_available"] is False
    assert methods["reproduction"]["fitted_uptake_kinetics_available"] is False
    assert methods["reproduction"]["cell_mass_conversion_verified"] is False
    assert methods["classifier"]["missing_workbook_features"] == ["oxygen_partial_pressure", "temperature"]


def test_partial_recovery_does_not_promote_inferred_cho_data_to_yeast_validation(manifest):
    data = manifest["supplementary_and_training_data"]
    assert data["original_training_inputs_available"] is False
    assert {record["sheet"] for record in data["acquired_training_tables"]} == {"ST2", "ST3", "ST4"}
    assert data["acquired_model_files"] == []
    assert data["normalized_process_measurements_available"] is True
    assert data["inferred_phase_labels_available"] is True
    assert data["declared_supplement_links_recovered"] == []
    assert data["missing"]
    host = manifest["paper"]["host_compatibility"]
    assert host["source_host_verified"] is True
    assert host["source_organism"] == "Chinese hamster"
    assert host["source_cell_type"] == "ovary (CHO)" and host["source_cell_line"] is None
    assert host["matched_yeast_host"] is False
    assert host["independent_measured_yeast_state_labels"] is False
    boundary = manifest["local_implementation_boundary"]
    assert boundary["biological_validation"] is False and boundary["demo_data_kind"] == "synthetic"
    for flag in (
        "scientific_outputs_regenerated",
        "scientific_outputs_restamped",
        "full_text_acquired",
        "training_data_acquired",
    ):
        assert manifest["acquisition_audit"][flag] is False
