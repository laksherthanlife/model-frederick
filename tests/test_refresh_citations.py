from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]


def test_missing_only_resolves_new_ids_and_preserves_existing_titles(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("citation_refresh", ROOT / "scripts/refresh_citations.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "TABLES", tmp_path)
    monkeypatch.setattr(module, "REPO", tmp_path)
    table = tmp_path / "doi_titles.csv"
    table.write_bytes(b'doi,registry,title,journal,year,cited_in\r\nold,crossref,"First\nsecond",J,2000,original\r\n')
    monkeypatch.setattr(module, "_cited", lambda: {"doi": {"old": {"retained"}, "new": {"source:2"}}})
    requested = []

    def resolve(ids):
        requested.extend(ids)
        return [{"doi": identity, "registry": "crossref", "title": "New title", "journal": "", "year": "None"}
                for identity in ids]

    monkeypatch.setitem(module.RESOLVERS, "doi", resolve)
    monkeypatch.setattr(sys, "argv", ["refresh_citations.py", "--kind", "doi", "--missing-only"])
    assert module.main() == 0
    assert requested == ["new"]
    with table.open(newline="") as handle:
        records = {row["doi"]: row for row in csv.DictReader(handle)}
    assert records["old"]["title"] == "First\nsecond"
    assert records["old"]["cited_in"] == "original"
    assert records["new"]["cited_in"] == "source:2"
    assert records["new"]["year"] == ""
    assert b"\r\n" not in table.read_bytes()
