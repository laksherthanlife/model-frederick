"""The fetcher for the external-validation datasets, exercised without a network.

``scripts/fetch_external_validation_data.py`` decides what ``run_external_validation.py``
will find on disk. Two of its properties are worth more than the download itself:

* a failed or malformed fetch must leave nothing behind. A truncated cache file is
  indistinguishable from a complete one on the next run, and the analysis would then be
  computed on a partial compendium with no error anywhere;
* the factor list has to match the module-to-factor table the analysis uses. A factor
  requested here but absent there is a wasted request; one named there and missing here
  makes ``load_regulons`` raise partway through a run.

``urlopen`` is replaced in every test that could reach the network, and the two cache
directories are redirected to a temporary one so the vendored ``data/external`` is never
touched.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
import urllib.error

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "fetch_external_validation_data.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("fetch_external_validation_data",
                                                  _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def sandboxed(script, tmp_path, monkeypatch):
    """Both cache directories under tmp, and no network for anyone who forgets to stub."""
    monkeypatch.setattr(script, "GASCH", tmp_path / "gasch2000")
    monkeypatch.setattr(script, "SGD", tmp_path / "sgd")
    monkeypatch.setattr(script.time, "sleep", lambda seconds: None)

    def refuse(*args, **kwargs):
        raise AssertionError("a test reached the network")

    monkeypatch.setattr(script.urllib.request, "urlopen", refuse)
    return tmp_path


class _Response:
    def __init__(self, payload: bytes):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self) -> bytes:
        return self._payload


def _serves(script, monkeypatch, payload: bytes, seen=None):
    def urlopen(url, timeout=None):
        if seen is not None:
            seen.append(url)
        return _Response(payload)

    monkeypatch.setattr(script.urllib.request, "urlopen", urlopen)


def _fails(script, monkeypatch, error=None):
    def urlopen(url, timeout=None):
        raise error or urllib.error.URLError("no route to host")

    monkeypatch.setattr(script.urllib.request, "urlopen", urlopen)


class TestAFailedFetchLeavesNothingBehind:
    """The next run trusts whatever is on disk. Anything written by a fetch that did not
    complete would be trusted too, and the compendium would be silently short."""

    def test_a_refused_connection_writes_no_file(self, script, monkeypatch, tmp_path):
        _fails(script, monkeypatch)
        destination = tmp_path / "GSE18-GPL51_series_matrix.txt.gz"

        assert script._get("https://example.invalid/x", destination) is False
        assert not destination.exists()

    def test_a_timeout_is_reported_as_a_failure_rather_than_raised(self, script,
                                                                  monkeypatch, tmp_path):
        """One platform failing must not abort the other eight; the return value carries
        the failure up to the exit status instead."""
        _fails(script, monkeypatch, TimeoutError("timed out"))

        assert script._get("https://example.invalid/x", tmp_path / "a.gz") is False

    def test_a_successful_fetch_writes_the_payload_and_makes_its_directory(
            self, script, monkeypatch, tmp_path):
        _serves(script, monkeypatch, b"payload")
        destination = tmp_path / "nested" / "a.gz"

        assert script._get("https://example.invalid/x", destination) is True
        assert destination.read_bytes() == b"payload"

    def test_a_file_already_on_disk_is_not_fetched_again(self, script, tmp_path):
        """The autouse fixture refuses every request, so this passes only if no request is
        made. Re-downloading a 60 MB matrix on every run is the thing this avoids."""
        destination = tmp_path / "a.gz"
        destination.write_bytes(b"cached")

        assert script._get("https://example.invalid/x", destination) is True
        assert destination.read_bytes() == b"cached"

    def test_an_empty_file_on_disk_is_fetched_again(self, script, monkeypatch, tmp_path):
        """A zero-byte file is what an interrupted write leaves. It is not a cache hit."""
        _serves(script, monkeypatch, b"payload")
        destination = tmp_path / "a.gz"
        destination.write_bytes(b"")

        assert script._get("https://example.invalid/x", destination) is True
        assert destination.read_bytes() == b"payload"


class TestARecordFileThatIsNotJsonIsRemoved:
    """SGD is a public API with no documented rate limit, and an error page comes back with
    a 200. Keeping one would cache an HTML document under ``MSN2.json`` and the next run
    would report it as present."""

    @pytest.fixture(autouse=True)
    def one_factor(self, script, monkeypatch):
        monkeypatch.setattr(script, "FACTORS", ("MSN2",))

    def test_a_non_json_payload_is_deleted_and_reported_as_a_failure(self, script,
                                                                     monkeypatch):
        _serves(script, monkeypatch, b"<html>service unavailable</html>")

        assert script.fetch_sgd() is False
        assert not (script.SGD / "MSN2.json").exists()

    def test_a_json_payload_is_kept(self, script, monkeypatch):
        _serves(script, monkeypatch, json.dumps([{"regulation_of": "transcription"}
                                                 ]).encode())

        assert script.fetch_sgd() is True
        assert json.loads((script.SGD / "MSN2.json").read_text())


class TestTheGaschRequestsCoverTheSeries:
    """GSE18 is printed on nine array batches and the stressors do not line up with them,
    so dropping one loses timepoints out of the middle of a time course."""

    def test_every_platform_gets_a_series_matrix(self, script, monkeypatch):
        seen = []
        _serves(script, monkeypatch, b"x", seen)

        script.fetch_gasch()

        matrices = [url for url in seen if "series_matrix" in url]
        assert len(matrices) == len(script.PLATFORMS)

    def test_only_the_platforms_with_a_published_table_get_an_annotation(self, script,
                                                                        monkeypatch):
        """GPL63 carries one array and GEO publishes no annotation for it, so its probes
        cannot be mapped to genes. Requesting it would 404 and fail the whole run."""
        seen = []
        _serves(script, monkeypatch, b"x", seen)

        script.fetch_gasch()

        annotations = [url for url in seen if url.endswith(".annot.gz")]
        assert len(annotations) == len(script.PLATFORMS) - len(script.NO_ANNOTATION)
        assert not any(f"GPL{gpl}.annot" in url
                       for url in annotations for gpl in script.NO_ANNOTATION)

    def test_a_single_failed_platform_fails_the_fetch(self, script, monkeypatch):
        _fails(script, monkeypatch)

        assert script.fetch_gasch() is False

    def test_the_unannotated_platforms_are_platforms(self, script):
        assert set(script.NO_ANNOTATION) <= set(script.PLATFORMS)


class TestTheFactorListIsTheAnalysisFactorList:
    """The comment says "every transcription factor named by a MODULES entry". That claim
    is checkable, and it is the difference between a complete regulon map and a run that
    dies on a missing file."""

    def test_the_requested_factors_are_exactly_the_ones_the_modules_name(self, script):
        from ystwin.bridge.regulation import MODULE_FACTORS

        named = {factor for factors in MODULE_FACTORS.values() for factor in factors}
        assert set(script.FACTORS) == named

    def test_no_factor_is_requested_twice(self, script):
        assert len(set(script.FACTORS)) == len(script.FACTORS)

    def test_the_copper_factor_is_requested_under_the_name_sgd_uses(self, script):
        """ACE1 is requested as CUP2, which is the name SGD and the GEO platform tables
        both use. Under the other name the file would 404 and copper would have no
        regulon."""
        assert "CUP2" in script.FACTORS
        assert "ACE1" not in script.FACTORS


class TestTheExitStatusSaysWhetherTheDataIsThere:
    """A pipeline runs this before the analysis. A zero exit with half the compendium
    downloaded would send the analysis into a partial dataset."""

    def test_a_complete_run_exits_zero(self, script, monkeypatch):
        monkeypatch.setattr(script, "FACTORS", ("MSN2",))
        _serves(script, monkeypatch, json.dumps([]).encode())

        assert script.main(["--skip-gasch"]) == 0

    def test_a_failed_fetch_exits_non_zero(self, script, monkeypatch):
        monkeypatch.setattr(script, "FACTORS", ("MSN2",))
        _fails(script, monkeypatch)

        assert script.main(["--skip-gasch"]) == 1

    def test_skipping_both_halves_fetches_nothing(self, script):
        """The autouse fixture refuses every request, so reaching one would fail here."""
        assert script.main(["--skip-gasch", "--skip-sgd"]) == 0

    def test_a_gasch_failure_alone_is_enough_to_fail(self, script, monkeypatch):
        _fails(script, monkeypatch)

        assert script.main(["--skip-sgd"]) == 1


class TestTheUrlTemplatesAddressTheRightThing:
    """Three format strings, each with one placeholder. A mismatched one produces a URL
    that 404s for every platform, which reads as the whole dataset being unavailable."""

    def test_the_series_matrix_url_names_its_platform(self, script):
        url = script.GEO_MATRIX.format(gpl="51")

        assert url.endswith("GSE18-GPL51_series_matrix.txt.gz")
        assert "{gpl}" not in url

    def test_the_annotation_url_names_its_platform_twice(self, script):
        """The path holds the platform directory and the file name, and both come from the
        same placeholder."""
        url = script.GEO_ANNOT.format(gpl="51")

        assert url.endswith("GPL51/annot/GPL51.annot.gz")

    def test_the_regulation_url_names_its_gene(self, script):
        url = script.SGD_REGULATION.format(gene="MSN2")

        assert url.endswith("/locus/MSN2/regulation_details")
