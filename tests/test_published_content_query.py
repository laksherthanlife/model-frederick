from __future__ import annotations

import csv
import math
from dataclasses import replace

import pytest

import ystwin.pathway.published_cassettes as module

_TABLE = module._TABLE


@pytest.fixture(autouse=True)
def clear_reader_cache():
    reader = module.published_cassettes
    reader.cache_clear()
    yield
    reader.cache_clear()


@pytest.fixture
def source_rows():
    if not _TABLE.is_file():
        pytest.skip("the survey is not vendored")
    with _TABLE.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture
def use_table(tmp_path, monkeypatch):
    def use(rows):
        path = tmp_path / "published_batch_titres.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        monkeypatch.setattr(module, "_TABLE", path)
        module.published_cassettes.cache_clear()

    return use


def _row(**changes):
    return {
        "pmid": "synthetic",
        "reference": "Independent synthetic source",
        "journal": "Synthetic journal",
        "strain": "beta-carotene and lycopene producer",
        "carbon_source": "not stated",
        "mode": "shake flask",
        "hours": "72",
        "content_mg_per_gdcw": "8",
        "titre_mg_per_l": "8",
        "biomass_g_per_l": "1",
        "bioreactor": "no",
        "measurand": "beta_carotene",
        "dfba_ready": "no",
        "titre_as_reported": "8 mg/L total carotenoids; lycopene in another condition",
        "notes": "A caption naming total carotenoids and beta-carotene is not an assay class",
        "copies_crtE": "",
        "copies_crtYB": "",
        "copies_crtI": "",
        "promoter": "",
        "integration": "",
        "copy_evidence": "",
    } | changes


def _legacy_record():
    return module.PublishedCassette(
        "legacy", "Legacy source", "strain", "glucose", frozenset({"glucose"}),
        True, {"crtE": 0, "crtYB": 2, "crtI": "2u"}, "TDH3p", "genomic",
        "Original copy evidence", 8.0,
    )


class TestCurrentCorpusContentQuery:
    def test_reader_retains_measurand_units_and_source_context(self, source_rows):
        cassettes = module.published_cassettes()

        assert len(cassettes) == len(source_rows) == 29
        assert {row["measurand"] for row in source_rows} == {
            "beta_carotene", "both_reported", "total_carotenoid", "unverified", "other",
        }
        for cassette, row in zip(cassettes, source_rows, strict=True):
            for field in (
                "pmid", "reference", "journal", "strain", "measurand", "mode", "hours",
                "bioreactor", "titre_as_reported", "notes", "copy_evidence", "integration",
            ):
                assert getattr(cassette, field) == row[field].strip(), (cassette.pmid, field)
            raw_content = row["content_mg_per_gdcw"]
            assert cassette.content_mg_per_gdcw == (float(raw_content) if raw_content else None)

    def test_archived_specific_guard_uses_current_reader_records(self, source_rows):
        cassettes = module.published_cassettes()
        above = module.rows_above_content(1.2483)

        assert isinstance(above, tuple)
        # A snapshot of the survey, so it moves when a row is worked up.
        # Fathi 2021 and Bu 2022 joined on 2026-09-09; see their notes for the assay evidence.
        assert [row.reference for row in above] == [
            "Verwaal 2007", "Lange 2011", "Xie 2014", "Fathi 2021", "Bu 2022", "Arhar 2024",
            "Bubphasawan 2025",
        ]
        assert all(any(row is original for original in cassettes) for row in above)
        hplc = [row for row in above if "HPLC" in row.notes]
        assert sorted(row.reference for row in hplc) == [
            "Arhar 2024", "Bu 2022", "Bubphasawan 2025", "Fathi 2021", "Lange 2011", "Xie 2014",
        ]
        assert max(row.content_mg_per_gdcw for row in hplc) / 1.2483 > 60.0
        assert not all(row.bioreactor == "yes" for row in above)
        assert all(
            any(description in row.mode.lower()
                for description in ("flask", "fed-batch", "not stated", "unverified"))
            for row in above
        )

    def test_opt_out_retains_total_and_unverified_values_without_reclassifying(self, source_rows):
        specific = module.rows_above_content(1.2483)
        everything = module.rows_above_content(1.2483, specific_only=False)

        # 17 since Sun 2020 gained a content (11.4 mg/gDCW) whose measurand is still a
        # 453 nm sum -- the case this opt-out exists to keep visible without promoting.
        assert len(everything) == 17 > len(specific)
        assert {row.measurand for row in specific} == {"beta_carotene", "both_reported"}
        assert {row.measurand for row in everything} == {
            "beta_carotene", "both_reported", "total_carotenoid", "unverified",
        }
        assert next(row for row in everything if row.pmid == "27423881").content_mg_per_gdcw == 25.52
        assert next(row for row in everything if row.pmid == "33332529").content_mg_per_gdcw == 46.5

    def test_both_reported_uses_specific_content_not_the_total_or_volumetric_titre(self, source_rows):
        row = next(row for row in module.rows_above_content(7) if row.pmid == "23860829")

        assert row.measurand == "both_reported"
        assert row.content_mg_per_gdcw == 7.41
        assert "11 mg/g DCW total carotenoids (72.57 mg/L)" in row.titre_as_reported
        assert row not in module.rows_above_content(7.41)
        assert row not in module.rows_above_content(10, specific_only=False)

    def test_lycopene_caption_is_not_a_total_or_a_verified_beta_carotene_measurement(self, source_rows):
        row = next(row for row in module.published_cassettes() if row.pmid == "32236768")

        assert row.measurand == "other"
        assert row.content_mg_per_gdcw is None
        assert "33.1 mg/g CDW and 1.48 g/L" in row.titre_as_reported
        assert "lycopene, not beta-carotene" in row.titre_as_reported
        assert row not in module.rows_above_content(0, specific_only=False)

    @pytest.mark.parametrize("specific_only", [True, False])
    def test_missing_dry_mass_contents_are_not_derived_from_other_columns(self, source_rows, specific_only):
        missing = {row["pmid"] for row in source_rows if not row["content_mg_per_gdcw"]}
        above = module.rows_above_content(0, specific_only=specific_only)

        assert {"33062054", "38213763", "32236768"} <= missing
        assert missing.isdisjoint(row.pmid for row in above)


class TestStructuredMeasurandControlsSelection:
    def test_shared_literals_and_misleading_captions_do_not_replace_the_measurand(self, use_table):
        measurands = (
            "beta_carotene", "both_reported", "total_carotenoid", "unverified", "other",
            "beta_carotene",
        )
        use_table([
            _row(pmid=str(index), reference=f"Independent source {index}", measurand=measurand)
            for index, measurand in enumerate(measurands)
        ])
        original = module.published_cassettes()

        assert [row.pmid for row in module.rows_above_content(7)] == ["0", "1", "5"]
        assert module.rows_above_content(7, specific_only=False) == original
        assert all(row.content_mg_per_gdcw == 8.0 for row in original)
        assert all(row.titre_as_reported == original[0].titre_as_reported for row in original)

    @pytest.mark.parametrize("measurand", ["", "lycopene", "beta_carotene_equivalents", "unknown"])
    def test_unclassified_or_unsupported_measurands_are_not_promoted(self, use_table, measurand):
        use_table([_row(measurand=measurand)])
        row, = module.published_cassettes()

        assert row.measurand == (measurand or None)
        assert module.rows_above_content(0) == ()
        assert module.rows_above_content(0, specific_only=False) == (row,)

    def test_legacy_constructor_keeps_cassette_semantics_but_does_not_claim_specificity(self, monkeypatch):
        row = _legacy_record()
        monkeypatch.setattr(module, "published_cassettes", lambda: (row,))

        assert row.measurand is None
        assert row.states_a_numeric_crtyb_dosage
        assert row.copies == {"crtE": 0, "crtYB": 2, "crtI": "2u"}
        assert row.copy_evidence == "Original copy evidence"
        assert module.rows_above_content(0) == ()
        assert module.rows_above_content(0, specific_only=False) == (row,)

    @pytest.mark.parametrize("measurand, selected", [
        ("beta_carotene", True), ("both_reported", True), ("total_carotenoid", False),
        ("unverified", False), ("other", False),
    ])
    def test_constructed_records_use_the_same_explicit_source_classes(self, monkeypatch, measurand, selected):
        row = replace(_legacy_record(), measurand=measurand)
        monkeypatch.setattr(module, "published_cassettes", lambda: (row,))

        assert module.rows_above_content(7) == ((row,) if selected else ())
        assert module.rows_above_content(7, specific_only=False) == (row,)


class TestContentDomainAndThresholds:
    @pytest.mark.parametrize("threshold", [-1.0, -1e-12, math.nan, math.inf, -math.inf])
    @pytest.mark.parametrize("specific_only", [True, False])
    def test_invalid_thresholds_fail_before_reading(self, monkeypatch, threshold, specific_only):
        def unexpected_load():
            pytest.fail("invalid thresholds must fail before the survey is read")

        monkeypatch.setattr(module, "published_cassettes", unexpected_load)
        with pytest.raises(ValueError, match="threshold_mg_per_gdcw.*finite.*nonnegative"):
            module.rows_above_content(threshold, specific_only=specific_only)

    def test_strict_comparison_keeps_zero_distinct_from_missing(self, use_table):
        use_table([
            _row(pmid="missing", content_mg_per_gdcw="", titre_mg_per_l="300", biomass_g_per_l="2"),
            _row(pmid="zero", content_mg_per_gdcw="0"),
            _row(pmid="equal", content_mg_per_gdcw="8"),
            _row(pmid="above", content_mg_per_gdcw="8.000001"),
        ])
        missing, zero, equal, above = module.published_cassettes()

        assert missing.content_mg_per_gdcw is None
        assert zero.content_mg_per_gdcw == 0.0
        assert module.rows_above_content(0) == (equal, above)
        assert module.rows_above_content(8) == (above,)
        assert module.rows_above_content(9) == ()

    @pytest.mark.parametrize("content", [
        "nan", "inf", "-inf", "1e9999", "-0.1", "UNVERIFIED", "unknown",
        "8 mg/gDCW", "8000 ug/gDCW", "8 mmol/gDCW", "8 mg/g wet weight", "8 mg/L",
    ])
    def test_invalid_or_unit_bearing_content_cells_fail_at_load(self, use_table, content):
        use_table([_row(content_mg_per_gdcw=content)])

        with pytest.raises(module.CassetteDomainError, match="row 2.*content_mg_per_gdcw.*mg/gDCW"):
            module.published_cassettes()

    @pytest.mark.parametrize("wrong_column", [
        "content_ug_per_gdcw", "content_mg_per_l", "content_mmol_per_gdcw", "content_mg_per_gwetweight",
    ])
    def test_different_quantity_headers_are_not_reinterpreted(self, use_table, wrong_column):
        row = _row()
        row[wrong_column] = row.pop("content_mg_per_gdcw")
        use_table([row])

        with pytest.raises(module.CassetteDomainError, match="content_mg_per_gdcw"):
            module.rows_above_content(0)

    def test_missing_measurand_column_is_not_inferred_from_a_caption(self, use_table):
        row = _row()
        del row["measurand"]
        use_table([row])

        with pytest.raises(module.CassetteDomainError, match="measurand"):
            module.rows_above_content(0)

    @pytest.mark.parametrize("content", [None, math.nan, math.inf, -math.inf, -1.0])
    @pytest.mark.parametrize("specific_only", [True, False])
    def test_nonmeasurements_on_constructed_records_do_not_pass_the_query(self, monkeypatch, content, specific_only):
        row = replace(_legacy_record(), measurand="beta_carotene", content_mg_per_gdcw=content)
        monkeypatch.setattr(module, "published_cassettes", lambda: (row,))

        assert module.rows_above_content(0, specific_only=specific_only) == ()

    def test_port_is_exported_and_specific_only_is_keyword_only(self):
        assert "rows_above_content" in module.__all__
        with pytest.raises(TypeError):
            module.rows_above_content(0, False)
