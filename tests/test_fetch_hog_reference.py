from __future__ import annotations

import hashlib
import json

import libsbml
import pytest

from scripts import fetch_hog_reference as fetcher


def sbml(model_id):
    document = libsbml.SBMLDocument(3, 2)
    document.createModel().setId(model_id)
    return libsbml.writeSBMLToString(document).encode()


def test_source_fetch_records_identity_and_reuses_only_verified_bytes(tmp_path):
    asset = fetcher.ASSETS["model_wt"]
    payload = sbml(asset.model_id)
    calls = []

    def download(url):
        calls.append(url)
        return payload, url

    result = fetcher.fetch_assets(tmp_path, ["model_wt"], download)
    record = result["assets"]["model_wt"]
    assert record["model_id"] == asset.model_id
    assert record["sha256"] == hashlib.sha256(payload).hexdigest()
    assert result["license"]["id"] == "CC-BY-4.0"
    fetcher.fetch_assets(tmp_path, ["model_wt"], download)
    assert calls == [asset.url]
    assert json.loads((tmp_path / "sources.json").read_text())["source_assumptions"]["published_fit"]


def test_wrong_model_identity_is_not_accepted_as_a_genotype_label(tmp_path):
    with pytest.raises(ValueError, match="identity mismatch"):
        fetcher.fetch_assets(tmp_path, ["model_wt"], lambda url: (sbml("wrong_genotype"), url))
    assert not (tmp_path / "model_wt.xml").exists()


def test_existing_artifacts_cannot_be_silently_replaced(tmp_path):
    asset = fetcher.ASSETS["model_wt"]
    payload = sbml(asset.model_id)
    fetcher.fetch_assets(tmp_path, ["model_wt"], lambda url: (payload, url))
    path = tmp_path / "model_wt.xml"
    path.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum changed"):
        fetcher.fetch_assets(tmp_path, ["model_wt"], lambda url: (payload, url))
    assert path.read_bytes() == b"changed"


def test_unknown_existing_file_and_html_workbook_are_refused(tmp_path):
    (tmp_path / "observations.xls").write_bytes(b"unmanifested")
    with pytest.raises(ValueError, match="unmanifested"):
        fetcher.fetch_assets(tmp_path, ["observations"])
    with pytest.raises(ValueError, match="Excel workbook"):
        fetcher._validate(fetcher.ASSETS["observations"], b"<html>not data</html>")
    with pytest.raises(ValueError, match="unknown source"):
        fetcher.fetch_assets(tmp_path, ["not_registered"])


def test_mutant_sources_follow_the_published_correction():
    assert fetcher.ASSETS["model_hog1_del"].supplement == "10.1371/journal.pcbi.1003663.s001"
    assert fetcher.ASSETS["model_hog1_att"].supplement == "10.1371/journal.pcbi.1003663.s002"
    assert fetcher.ASSETS["model_fps1_open"].supplement == "10.1371/journal.pcbi.1003663.s003"
