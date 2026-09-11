"""Replaying the committed plate text, and the seam that lets any script do it.

`data/plates` holds the four NewProtocol exports as text because the workbooks themselves
carry a named private individual in their document properties and cannot be committed.
That made the *measurements* reproducible by anyone. It did not make the *scripts*
reproducible: the replay lived inside one script as a private shim, so `plate_readings.py`
ran from the committed text and every other plate-reading script exited on a missing
`YSTWIN_PLATES` with the numbers sitting on disk beside it.

`ystwin.plate.replay` is that shim moved into the library. These tests are about the seam
-- what `install` rebinds, what it leaves alone, and the distinctions the pipeline above it
depends on being preserved.
"""

from __future__ import annotations

import pathlib
import types

import pandas as pd
import pytest

from ystwin.plate import replay


@pytest.fixture(scope="module")
def committed():
    """The tracked text, or a skip. It is committed, so the skip should never fire."""
    src = replay.committed_dir()
    if not (src / replay.MANIFEST).exists():
        pytest.skip(f"{src / replay.MANIFEST} is tracked and absent; this is a bad checkout")
    return src


class TestTheSeamRebindsWhatAModuleActuallyHas:
    """`install` is the whole mechanism. It has to rebind both readers where both exist,
    rebind one where one exists, and never invent a name on a module that has neither --
    `run_calibration_nis.py` reads kinetics and never reads a dose sheet."""

    def test_it_rebinds_both_readers_when_both_are_present(self, committed):
        module = types.SimpleNamespace(read_synergy_kinetic="original",
                                       read_dose_response="original")

        replay.install(module)

        assert module.read_synergy_kinetic != "original"
        assert module.read_dose_response != "original"

    def test_it_rebinds_only_the_reader_a_module_has(self, committed):
        module = types.SimpleNamespace(read_synergy_kinetic="original")

        replay.install(module)

        assert module.read_synergy_kinetic != "original"
        assert not hasattr(module, "read_dose_response")

    def test_it_returns_the_export_names_a_caller_iterates_instead_of_a_directory(
            self, committed):
        module = types.SimpleNamespace(read_synergy_kinetic=None)

        names = replay.install(module)

        assert names and all(isinstance(n, pathlib.Path) for n in names)
        # `.xpt` joined `.xlsx` when 20260804 was re-exported from the raw instrument
        # file to recover its blanks. The seam keys on the name, so the source format is
        # the reader's business and not this module's -- but it is pinned rather than
        # left open, so a third format has to be a decision somebody makes.
        assert {n.suffix for n in names} <= {".xlsx", ".xpt"}

    def test_the_names_it_returns_are_names_and_not_files(self, committed):
        """The paths are never opened. A test that asserted they existed would be asserting
        the workbooks are present, which is the thing this module exists because they are
        not."""
        assert not any(name.exists() for name in replay.exports())

    def test_a_rebound_reader_takes_the_path_the_pipeline_would_have_opened(self,
                                                                            committed):
        """The pipeline passes whole paths, so the replacement has to accept one and key on
        its name alone -- including a path with a directory in front of it that does not
        exist here."""
        module = types.SimpleNamespace(read_synergy_kinetic=None)
        names = replay.install(module)

        run = module.read_synergy_kinetic(pathlib.Path("/nowhere/at/all") / names[0].name)

        assert run.blocks


class TestTheRebuiltRunCarriesWhatTheManifestSaysItDoes:
    def test_every_export_rebuilds_with_the_blocks_the_manifest_lists(self, committed):
        manifest = replay.load_manifest(committed)

        for name in manifest.export.drop_duplicates():
            expected = list(manifest[manifest.export == name].channel)

            assert [b.channel for b in replay.load_run(name, committed).blocks] == expected

    def test_the_block_order_is_the_manifest_order(self, committed):
        """The manifest preserves exported channel identities and alignment traversal;
        ordering does not choose among raw candidates."""
        manifest = replay.load_manifest(committed)
        name = manifest.export.iloc[0]

        rebuilt = replay.load_run(name, committed)

        rows = manifest[manifest.export == name]
        assert [(b.channel, b.sheet) for b in rebuilt.blocks] == list(
            zip(rows.channel, rows.sheet))

    def test_the_time_index_is_named_what_the_reader_names_it(self, committed):
        run = replay.load_run(replay.exports(committed)[0].name, committed)

        assert all(b.data.index.name == "time_h" for b in run.blocks)


class TestTheDistinctionsAboveThisLayerDependOn:
    def test_an_export_with_no_committed_blocks_raises_rather_than_returning_empty(
            self, committed):
        with pytest.raises(ValueError, match="no committed blocks"):
            replay.load_run("never_plated.xlsx", committed)

    def test_a_plate_with_no_dose_sheet_raises_rather_than_returning_empty(self,
                                                                          committed):
        """``collect()`` tells "this plate has no dose ladder of its own" from "this plate
        is unusable" by catching ``ValueError``. Returning an empty frame would merge the
        two and quietly drop a usable plate."""
        with pytest.raises(ValueError, match="no per-construct dose-response sheet"):
            replay.load_doses("never_plated.xlsx", committed)

    def test_at_least_one_committed_export_genuinely_has_no_dose_sheet(self, committed):
        """The refusal above is only worth having if the case is real. It is: not every
        export carries a per-construct sheet, which is why the pipeline catches rather
        than assumes."""
        manifest = replay.load_manifest(committed)
        without = [name for name in manifest.export.drop_duplicates()
                   if not {f for f in manifest[manifest.export == name].dose_response_file
                           if f}]

        assert without

    def test_a_missing_manifest_names_both_ways_of_supplying_the_plates(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="YSTWIN_PLATES"):
            replay.load_manifest(tmp_path)


class TestTheCommittedTextIsReadAsTextAndNotGuessedAt:
    def test_the_optics_label_stays_a_string(self, committed):
        """``optics`` holds ``600``, a wavelength label rather than a number. Inferred as
        an int it would not compare equal to the block it came from."""
        manifest = replay.load_manifest(committed)

        assert manifest.optics.map(type).eq(str).all()

    def test_a_blank_field_is_an_empty_string_and_not_a_nan(self, committed):
        """``dose_response_file`` is empty for exports without one, and the membership test
        that reads it (``if f``) treats an empty string as absent and a NaN as present."""
        manifest = replay.load_manifest(committed)

        assert not manifest.dose_response_file.isna().any()

    def test_the_time_string_round_trips_through_the_expression_synergy_uses(self,
                                                                             committed):
        """``h + m/60 + s/3600`` and ``round(h*3600 + m*60 + s)/3600`` are the same real
        number and different floats. The committed table reproduces bit for bit or not at
        all, so the expression is the claim."""
        h, m, s = 0.0, 1.0, 1.0

        assert replay._from_hms("0:01:01") == h + m / 60.0 + s / 3600.0
        assert replay._from_hms("0:01:01") != round(h * 3600 + m * 60 + s) / 3600.0


class TestItReplaysTheNumbersAndNotAnApproximationOfThem:
    def test_a_rebuilt_block_equals_the_text_it_was_built_from(self, committed):
        """Exact equality, not ``allclose``. The text is stored at a precision that
        round-trips, so any difference at all is a defect rather than a rounding."""
        manifest = replay.load_manifest(committed)
        row = manifest.iloc[0]
        stored = pd.read_csv(committed / row.readings_file)

        block = replay.load_run(row.export, committed).blocks[0]

        wells = [c for c in stored.columns if c not in {"elapsed_hms", "temperature_c"}]
        assert list(block.data.columns) == wells
        assert (block.data.to_numpy() == stored[wells].to_numpy(dtype=float)).all()


class TestASwapMustNotBreakTheSourceItWasNotFor:
    """The subtlest thing this module does, and the one that broke first.

    `run_gates.py` reads two plate sets, and when this was written only one of them was
    committed as text. A swap that redirected *everything* to the replay turned a fix for
    that one into a break in the other -- every July export came back "no committed blocks
    for ...", and the script exited 0 having gated half of a comparison whose whole point
    is both halves. Both sets are committed now, which removes the occasion but not the
    rule: a workbook the manifest does not cover must still reach the real reader.

    Whether the manifest covers a name is the authoritative test, and it is the only thing
    that actually separates the two cases.
    """

    def test_a_covered_export_goes_to_the_replay(self, committed):
        original_calls = []
        module = types.SimpleNamespace(
            read_synergy_kinetic=lambda path: original_calls.append(path))
        names = replay.install(module)

        run = module.read_synergy_kinetic(names[0])

        assert original_calls == []
        assert run.blocks

    def test_an_uncovered_export_still_reaches_the_original_reader(self, committed):
        original_calls = []
        module = types.SimpleNamespace(
            read_synergy_kinetic=lambda path: original_calls.append(path) or "from workbook")
        replay.install(module)

        got = module.read_synergy_kinetic(pathlib.Path("a_plate_never_exported_here.xlsx"))

        assert got == "from workbook"
        assert len(original_calls) == 1

    def test_an_uncovered_export_with_no_original_reader_is_refused_clearly(self,
                                                                            committed):
        """When the name was never bound there is nothing to fall back to, and the message
        has to say which export it could not find rather than raising an AttributeError
        about a None."""
        module = types.SimpleNamespace(read_synergy_kinetic=None)
        replay.install(module)

        with pytest.raises(ValueError, match="no committed blocks for 'unknown.xlsx'"):
            module.read_synergy_kinetic(pathlib.Path("unknown.xlsx"))

    def test_installing_into_globals_works_the_same_way(self, committed):
        """Scripts pass ``globals()`` rather than ``sys.modules[__name__]``, because a
        script loaded through ``spec_from_file_location`` -- which is how every script test
        here loads one -- is never registered in ``sys.modules``."""
        namespace = {"read_synergy_kinetic": lambda path: "from workbook"}

        names = replay.install(namespace)

        assert namespace["read_synergy_kinetic"](names[0]).blocks
        assert namespace["read_synergy_kinetic"]("other.xlsx") == "from workbook"

    def test_it_does_not_create_a_name_the_namespace_did_not_have(self, committed):
        namespace = {"read_synergy_kinetic": lambda path: None}

        replay.install(namespace)

        assert "read_dose_response" not in namespace


@pytest.mark.parametrize("field", ["blank_subtracted", "derived"])
@pytest.mark.parametrize("invalid", ["", "false", "0", "unknown"])
def test_replay_does_not_turn_invalid_boolean_metadata_into_false(field, invalid, committed):
    manifest = replay.load_manifest(committed)
    name = manifest.export.iloc[0]
    manifest.loc[manifest.export == name, field] = invalid

    with pytest.raises(ValueError, match="boolean"):
        replay.load_run(name, committed, manifest)


def test_public_august3_processing_declarations_match_the_verified_source_mapping(committed):
    import numpy as np

    run = replay.load_run("20260803_ER&oxidativestress_Replicate3.xlsx", committed)
    assert {b.channel: b.blank_subtracted for b in run.blocks} == {
        "OD600[1]": False, "OD600[2]": True, "mCitrine[1]": True, "mCitrine[2]": False,
    }
    assert run.raw_channel("OD600").channel == "OD600[1]"
    assert run.raw_channel("mCitrine").channel == "mCitrine[2]"
    processed = run.channel("mCitrine[1]")
    original = run.channel("mCitrine[2]")
    offsets = original.data.to_numpy() - processed.data.to_numpy()
    blank = original.data[["H1", "H2", "H3"]].mean(axis=1).to_numpy()

    assert np.max(np.ptp(offsets, axis=1)) == 0.0
    assert np.max(np.abs(offsets - blank[:, None])) <= 1 / 3 + 1e-10


def test_all_current_public_plate_sources_have_unique_declared_raw_channels(committed):
    manifest = replay.load_manifest(committed)
    current = manifest[manifest.source_set != "newprotocol_superseded"]
    for name, rows in current.groupby("export", sort=False):
        run = replay.load_run(name, committed, manifest)
        for fluorophore in rows.fluorophore.unique():
            assert run.raw_channel(fluorophore).blank_subtracted is False


def test_two_declared_raw_blocks_are_ambiguous_even_if_their_signs_differ(committed):
    name = "20260803_ER&oxidativestress_Replicate3.xlsx"
    manifest = replay.load_manifest(committed)
    manifest.loc[manifest.export == name, "blank_subtracted"] = "False"
    run = replay.load_run(name, committed, manifest)

    with pytest.raises(ValueError, match="ambiguous raw"):
        run.raw_channel("mCitrine")
    assert all(b.blank_subtracted is False for b in run.blocks)


def test_processing_provenance_binds_declarations_to_unchanged_public_csvs(committed):
    import hashlib
    import json

    evidence = json.loads((committed / "processing_provenance.json").read_text())
    manifest = replay.load_manifest(committed)
    changed = set()
    for record in evidence["records"]:
        rows = manifest[manifest.export == record["export"]].set_index("channel")
        assert set(rows.source_sha256) == {record["source_sha256"]}
        assert set(rows.source_bytes) == {str(record["source_bytes"])}
        for block in record["blocks"]:
            row = rows.loc[block["channel"]]
            assert row.readings_file == block["readings_file"]
            assert row.blank_subtracted == str(block["blank_subtracted"])
            content = (committed / block["readings_file"]).read_bytes()
            assert hashlib.sha256(content).hexdigest() == block["readings_sha256"]
            if "reference_readings_file" in block:
                reference = (committed / block["reference_readings_file"]).read_bytes()
                assert hashlib.sha256(reference).hexdigest() == block["reference_readings_sha256"]
            if block["blank_subtracted"] != block["previous_blank_subtracted"]:
                changed.add((record["export"], block["channel"]))
    assert changed == {
        ("20260803_ER&oxidativestress_Replicate3.xlsx", "OD600[2]"),
        ("20260803_ER&oxidativestress_Replicate3.xlsx", "mCitrine[1]"),
        ("20260804_ER&OxidativeStress_Replicate4.xlsx", "OD600"),
        ("20260804_ER&OxidativeStress_Replicate4.xlsx", "mCitrine"),
    }


def test_known_source_declarations_refuse_duplicate_block_locators(monkeypatch):
    row = {"source_sha256": "a" * 64, "sheet": "opaque", "fluorophore": "OD600",
           "optics": "600", "blank_subtracted": "False"}
    monkeypatch.setattr(replay, "load_manifest", lambda: pd.DataFrame([row, row]))

    with pytest.raises(ValueError, match="ambiguous correction-state"):
        replay._processing_states_for_source("a" * 64)


def test_a_missing_manifest_leaves_source_processing_unknown(monkeypatch):
    def absent():
        raise FileNotFoundError("no manifest")

    monkeypatch.setattr(replay, "load_manifest", absent)
    assert replay._processing_states_for_source("a" * 64) == {}


def test_replay_reader_options_cannot_override_known_processing(committed):
    namespace = {"read_synergy_kinetic": None}
    replay.install(namespace, committed)
    name = "20260803_ER&oxidativestress_Replicate3.xlsx"
    read = namespace["read_synergy_kinetic"]

    run = read(name, processing_states={"OD600[1]": False, "OD600[2]": True})
    assert run.raw_channel("OD600").channel == "OD600[1]"
    with pytest.raises(ValueError, match="conflicts with source metadata"):
        read(name, processing_states={"OD600[2]": False})
