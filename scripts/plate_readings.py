"""The biosensor plates as committed text, and the dose response derived from it.

Nothing in this repository could produce a real-data number from a fresh clone. Every
tracked table in ``outputs/`` was a claim about what the code does rather than something a
stranger could re-derive, because the four NewProtocol Synergy exports the whole sensor
result rests on live outside the checkout and resolve through ``YSTWIN_PLATES``.

The obvious fix -- commit the four ``.xlsx`` -- is the wrong one, and not for the reason
you would guess. They are only 989 kB together, which is small enough to track. But every
one of them carries a named private individual in ``docProps/core.xml``
(``cp:lastModifiedBy``), and committing the workbook publishes that name into git history
permanently, where removing it costs a history rewrite. Personal data is not ours to
redistribute, and laundering the metadata out would leave a binary that no longer matches
the lab's own copy while still diffing to nothing.

So what is committed here is the **measurement matrix and nothing else**: every value the
reader wrote, as text, at a precision that round-trips bit for bit. That is smaller than
the workbooks (roughly a third), it diffs, and it carries no author, no revision history
and no embedded metadata -- only optical densities, relative fluorescence, elapsed times
and temperatures. The ``sha256`` of each source workbook is recorded in the manifest so a
holder of the original can prove these came from it without the original being here.

Fidelity is checked rather than asserted. ``--export`` reconstructs every block from the
text it just wrote and refuses to leave the export in place unless each value, each
elapsed time and each well ordering is *identical* -- ``==`` on the float, not close to it.
Elapsed times are stored as ``H:MM:SS`` and rebuilt with the same ``h + m/60 + s/3600``
expression ``plate/synergy.py`` uses, because ``round(t*3600)/3600`` is a different
floating-point expression and lands on a different float.

Two sources, one pipeline. Both modes call the real ``collect()`` from
``run_sensor_characterisation.py`` -- this script does not reimplement the analysis, it
only supplies the readers ``collect()`` reads through. Which source was used is printed as
the first and last thing you see, because a fallback that reads like a fresh derivation is
worse than no fallback: it turns "I re-derived this" into a claim nobody checked.

Usage:
    python scripts/plate_readings.py                 # reproduce; raw plates if present
    python scripts/plate_readings.py --source committed   # force the committed text
    python scripts/plate_readings.py --export        # rebuild data/plates from the plates
    python scripts/plate_readings.py --write         # also write outputs/ tables
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.plate import gen5, replay
from ystwin.plate.synergy import KineticBlock

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PLATES_DIR = paths.data_dir() / "plates"
MANIFEST = replay.MANIFEST  # one name for the file, in the module that reads it

MANIFEST_COLUMNS = [
    "export", "source_set", "channel", "fluorophore", "optics", "sheet", "derived",
    "blank_subtracted", "n_times", "n_wells", "readings_file",
    "dose_response_file", "source_sha256", "source_bytes",
]

#: Which raw plate set an export came from, keyed by the resolver that found it. Recorded
#: because two scripts need to tell the sets apart and the alternative is parsing the date
#: out of the filename -- a heuristic that is wrong the first time a plate is renamed, and
#: that would have to be repeated in every consumer. `run_d2.py` reads the July set only;
#: `run_gates.py` reads both and reports the comparison between them.
SOURCE_SETS = {
    "biosensor_plates": "newprotocol",
    "igem_results": "july_2026_07",
    "gen5_xpt_dir": "newprotocol",
}

INSTRUMENT_FILES = {
    "20260804_ER_Oxidative_Replicate4.xpt": "newprotocol",
}
"""Plates committed from the instrument's own ``.xpt`` rather than from an ``.xlsx``.

Named one by one, and short on purpose. ``paths.gen5_xpt_dir()`` holds every experiment
this reader has ever run -- AFL circuits, ATP sensors, June bring-up plates, twenty-two
files. Globbing it would commit twenty tables nobody asked for and, worse, hand them all to
``collect()``, which reads whatever ``data/plates`` covers. So a file earns a line here by
being a plate the analysis needs and cannot get any other way.

``20260804`` is the one that qualifies. Every surviving ``.xlsx`` of it went through Gen5's
blank-subtraction transform -- its media blanks H1-H3 land within 0.002 of zero -- so the
background is gone from the export and the plate sat outside the panel at n=3. The ``.xpt``
still has it: 0.087/0.090/0.088 OD and 347/349/354 RFU at the first read. Subtracting the
per-timepoint mean of H1:H3 from the archive reproduces the committed export exactly in
fluorescence (2175 of 2175 readings, once Gen5's rounding to whole RFU is applied) and to
5e-4 in absorbance, which is all a three-decimal export can say about a four-decimal
reading. ``tests/test_gen5_xpt.py`` holds that comparison.

Both files stay committed. The subtracted export is the evidence of what was lost and is
what those tests compare against; the archive text is what the analysis reads. Which of the
two is the biological replicate is decided in ``run_sensor_characterisation.collect()``, by
the plate they share, and not here.
"""


def _slug(text: str) -> str:
    """The filename convention already used by ``run_gates.py`` and ``run_d2.py``."""
    return text.replace(" ", "_").replace("&", "and")


def _channel_slug(channel: str) -> str:
    """``mCitrine[2]`` -> ``mCitrine_2``, so a repeat read stays distinguishable."""
    return channel.replace("[", "_").replace("]", "")


def _hms(hours: float) -> str:
    """Elapsed hours -> ``H:MM:SS``.

    Every timestamp in these exports is a whole number of seconds, which is checked
    here rather than assumed: a fractional second would make the stored string lossy
    and the round-trip check downstream would fail without saying why.
    """
    total = hours * 3600.0
    seconds = int(round(total))
    if abs(total - seconds) > 1e-6:
        raise ValueError(
            f"elapsed time {hours!r} h is not a whole number of seconds "
            f"({total!r} s); H:MM:SS cannot store it losslessly"
        )
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def _from_hms(text: str) -> float:
    """``H:MM:SS`` -> hours, using ``plate/synergy.py::_to_hours``'s own expression.

    The expression matters more than the value. ``h + m / 60 + s / 3600`` and
    ``round(h * 3600 + m * 60 + s) / 3600`` are the same real number and different
    floats, and the committed table is reproduced byte for byte or not at all.
    """
    h, m, s = (float(part) for part in text.split(":"))
    return h + m / 60.0 + s / 3600.0


def _by_suffix(workbook_reader, archive_reader):
    """One reader over both plate formats, dispatching on the file's own extension.

    Not a heuristic on the name: the extension *is* the format here, and the two formats
    are not two spellings of one thing. A ``.xlsx`` is what Gen5 exported -- whatever the
    operator had on screen, at three decimals, through whatever transform was selected. A
    ``.xpt`` is the experiment file Gen5 itself wrote. Handing both to one reader is what
    lets ``export()`` and ``verify()`` stay a single code path across the two.
    """
    def read(path):
        path = pathlib.Path(path)
        if path.suffix.lower() == ".xpt":
            return archive_reader(path)
        return workbook_reader(path)
    return read


def _no_dose_sheet(path) -> None:
    """A ``.xpt`` has no per-construct sheets, and says so the way ``collect()`` listens.

    Those sheets are Gen5's per-construct plots, assembled when a workbook is exported;
    the archive holds plate reads. ``collect()`` distinguishes "this plate has no dose
    ladder of its own" from "this plate is unusable" by catching :class:`ValueError`, so
    the absence is raised rather than returned empty -- and for a plate in the logbook
    registry that costs nothing but the cross-check, since the ladder and the layout are
    recorded. On ``20260804`` the two agree: the ladder in its own subtracted export is
    DTT 0/0.1/0.2/0.5/1/2/5 and H2O2 0/0.1/0.2/0.5/1/2/4, which is the recorded one, and
    recovering the layout from the archive returns the recorded columns.
    """
    raise ValueError(
        f"{pathlib.Path(path).name} is a Gen5 instrument file: it holds plate reads and "
        "no per-construct derived sheet")


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


# --------------------------------------------------------------------------- export


def _write_block(block: KineticBlock, dest: pathlib.Path, stem: str) -> str:
    """One block -> one wide CSV whose header *is* the well ordering.

    Per block rather than per export, because ``20260804`` reads its two channels in
    two different well orders, and that ordering reaches the committed output: it sets
    the order of ``shared`` in ``prepare()`` and therefore the row order of
    ``sensor_characterisation.csv``. A union header with blanks would lose it.
    """
    name = f"{stem}__{_channel_slug(block.channel)}.csv"
    missing = int(block.data.isna().to_numpy().sum())
    if missing:
        # `plate/gen5.py` turns the instrument's over-range flag into NaN, which is the
        # honest expression of "no value reported" everywhere EXCEPT here: committed text
        # is the record that a reading was taken, and a blank cell in it is indistinguishable
        # from a well that was never plated. Nothing in the corpus reaches this -- the
        # committed channels are the gain-75 ones and none of them overflows -- so this is
        # a door held shut rather than one that has been pushed.
        raise ValueError(
            f"{block.channel} on {stem}: {missing} of {block.data.size} readings are "
            "missing, which for a .xpt means the instrument flagged them over-range. "
            "Committed text records measurements; export a channel that has them")
    frame = block.data.copy()
    frame.insert(0, "elapsed_hms", [_hms(float(t)) for t in block.data.index])
    if block.temperature_c:
        if len(block.temperature_c) != len(frame):
            raise ValueError(
                f"{block.channel} on {stem}: {len(block.temperature_c)} temperatures "
                f"for {len(frame)} timepoints"
            )
        frame.insert(1, "temperature_c", list(block.temperature_c))
    frame.to_csv(dest / name, index=False)
    return name


def _write_dose_response(frame: pd.DataFrame, dest: pathlib.Path, stem: str) -> str:
    """The per-construct sheets, which carry the dose ladder and feed layout recovery.

    Kept because dropping them changes the answer rather than merely the provenance:
    without them ``collect()`` takes the dose ladder from the logbook registry instead
    of the plate, and the two are not guaranteed to agree.
    """
    name = f"{stem}__dose_response.csv"
    out = pd.DataFrame({
        "construct": frame.construct,
        "dose_mM": frame.dose_mM,
        "elapsed_min": frame.time_h * 60.0,
        "signal": frame.signal,
    })
    out.to_csv(dest / name, index=False)
    return name


def export(exports: list[pathlib.Path], dest: pathlib.Path,
           read_run, read_doses, source_set: dict | None = None) -> pd.DataFrame:
    """Write every export's measurement matrix as text, and its block metadata.

    Stale files are removed by name from the previous manifest rather than by globbing
    ``*.csv``, so a renamed export cannot leave an orphan behind and a mistyped
    ``--dir`` cannot delete somebody's unrelated tables.
    """
    dest.mkdir(parents=True, exist_ok=True)
    source_set = source_set or {}

    # READ EVERYTHING FIRST, THEN CLEAR. The previous version cleared the old files at the
    # top and then read the workbooks, so anything that failed mid-loop -- one .xlsx in the
    # directory that is not a kinetic export was enough -- left the committed text deleted
    # and the manifest gone. The tracked measurements are the only wet-lab data in this
    # repository, and an export that crashes must leave them exactly as it found them.
    #
    # Nothing is written until every workbook has been read, and the clear happens between
    # the reads and the writes rather than before both.
    read: list[tuple[pathlib.Path, str, str, int, object, object, str]] = []
    for path in exports:
        try:
            run = read_run(path)
        except ValueError as exc:
            # Not every .xlsx beside the exports is one. A plate map, a summary, a one-off
            # carries no kinetic block and is not a measurement to commit.
            print(f"  skipped {path.name}: {exc}")
            continue
        undeclared = [block.channel for block in run.blocks
                      if not isinstance(block.blank_subtracted, (bool, np.bool_))]
        if undeclared:
            raise ValueError(
                f"{path.name}: correction-state is undeclared for {undeclared}; supply "
                "explicit processing_states before exporting. Existing text is unchanged")
        try:
            doses = read_doses(path)
        except ValueError:
            doses = None
        read.append((path, _slug(path.stem), _sha256(path), path.stat().st_size, run, doses,
                     source_set.get(path, "")))
    if not read:
        raise SystemExit(
            f"no export among {len(exports)} carried a kinetic block; nothing was written "
            f"and {dest} is unchanged")

    # Clear only what THIS run is rewriting, and keep every other export's rows.
    #
    # The previous version cleared every file the manifest named, which made an export
    # all-or-nothing: the two wet-lab plate sets live in two directories, usually only one
    # is configured, and exporting from one deleted the other's committed text. Scoping the
    # clear to the exports in hand makes a partial export safe by construction rather than
    # forbidden -- and orphans, which the wholesale clear existed to prevent, are still
    # impossible, because a renamed export's old files are named by the rows being replaced.
    rewriting = {path.name for path, *_ in read}
    # A RENAME is not a different export, and the manifest already records the only thing
    # that can tell the two apart: the workbook's sha256. An entry whose name is absent
    # from this run but whose digest is present has been renamed, and its old files are
    # orphans -- a measurement nothing points at, which the next reader cannot distinguish
    # from a current one. An entry whose digest is also absent is a plate set this machine
    # simply cannot see, and its files stay exactly where they are.
    digests = {digest for _, _, digest, *_ in read}
    kept: list[dict] = []
    if (dest / MANIFEST).exists():
        previous = pd.read_csv(dest / MANIFEST, dtype=str, keep_default_na=False)
        superseded = previous.export.isin(rewriting) | previous.source_sha256.isin(digests)
        replaced = previous[superseded]
        stale = set(replaced.readings_file) | {f for f in replaced.dose_response_file if f}
        for name in sorted(stale | {MANIFEST}):
            (dest / name).unlink(missing_ok=True)
        renamed = set(replaced.export) - rewriting
        if renamed:
            print(f"  renamed since the last export, old tables removed: "
                  f"{', '.join(sorted(renamed))}")
        kept = previous[~superseded].to_dict("records")
        if kept:
            print(f"  keeping {len({r['export'] for r in kept})} export(s) already "
                  f"committed here and not in this run")
            # A kept row carries whatever it already had. Rows written before `source_set`
            # existed have none, and that is reported rather than guessed: filling it in
            # from "everything already here must be the other set" would be right today
            # and wrong the first time a third plate set is committed. Re-exporting with
            # that set's variable set is what fills them, and this names them so somebody
            # can.
            unlabelled = sorted({r["export"] for r in kept if not r.get("source_set")})
            if unlabelled:
                print(f"  {len(unlabelled)} kept export(s) predate the source_set column "
                      f"and carry none: {', '.join(unlabelled)}. Re-export with their own "
                      f"data variable set to fill it.")

    rows: list[dict] = list(kept)
    for path, stem, digest, size, run, doses, which in read:
        dose_file = "" if doses is None else _write_dose_response(doses, dest, stem)
        for block in run.blocks:
            rows.append({
                "export": path.name,
                "source_set": which,
                "channel": block.channel,
                "fluorophore": block.fluorophore,
                "optics": block.optics,
                "sheet": block.sheet,
                "derived": block.derived,
                "blank_subtracted": block.blank_subtracted,
                "n_times": len(block.data),
                "n_wells": len(block.data.columns),
                "readings_file": _write_block(block, dest, stem),
                "dose_response_file": dose_file,
                "source_sha256": digest,
                "source_bytes": size,
            })
    manifest = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    manifest.to_csv(dest / MANIFEST, index=False)
    return manifest


# --------------------------------------------------------------------------- load


# The three loaders that used to live here are now ``ystwin.plate.replay``. They were
# moved out unchanged: a shim that lets exactly one script run from the committed text is
# not reproducibility, it is one reproducible script, and every other plate-reading script
# was still exiting on a missing YSTWIN_PLATES with the numbers sitting in data/plates.


# --------------------------------------------------------------------------- verify


def _identical(left: pd.DataFrame, right: pd.DataFrame) -> str | None:
    """Exact equality, or the first difference found. Not ``allclose``."""
    if list(left.columns) != list(right.columns):
        return f"well ordering differs: {list(left.columns)[:6]} vs {list(right.columns)[:6]}"
    if left.index.name != right.index.name:
        return f"index name differs: {left.index.name!r} vs {right.index.name!r}"
    lt, rt = np.asarray(left.index, dtype=float), np.asarray(right.index, dtype=float)
    if lt.shape != rt.shape or not np.array_equal(lt, rt):
        return "elapsed times differ"
    lv, rv = left.to_numpy(dtype=float), right.to_numpy(dtype=float)
    if lv.shape != rv.shape:
        return f"shape differs: {lv.shape} vs {rv.shape}"
    bad = np.argwhere(lv != rv)
    if bad.size:
        i, j = bad[0]
        return (f"value differs at {left.index[i]!r}/{left.columns[j]}: "
                f"{lv[i, j]!r} vs {rv[i, j]!r}")
    return None


def verify(exports: list[pathlib.Path], src: pathlib.Path, read_run) -> bool:
    """Every committed block must equal the workbook it came from, exactly.

    Checks the exports the manifest covers and no others. A workbook sitting beside them
    that carries no kinetic block -- a plate map, a summary -- was never exported and has
    nothing to compare against; treating it as a failure would make the check depend on
    what else happens to be in the directory.
    """
    manifest = replay.load_manifest(src)
    covered = set(manifest.export)
    exports = [p for p in exports if p.name in covered]
    ok = True
    for path in exports:
        original = read_run(path)
        rebuilt = replay.load_run(path.name, src, manifest)
        if [b.channel for b in original.blocks] != [b.channel for b in rebuilt.blocks]:
            print(f"  FAIL {path.name}: block order or labels differ")
            ok = False
            continue
        for a, b in zip(original.blocks, rebuilt.blocks):
            problem = _identical(a.data, b.data)
            meta = [
                f"fluorophore {a.fluorophore!r} vs {b.fluorophore!r}"
                if a.fluorophore != b.fluorophore else "",
                f"optics {a.optics!r} vs {b.optics!r}" if a.optics != b.optics else "",
                f"sheet {a.sheet!r} vs {b.sheet!r}" if a.sheet != b.sheet else "",
                "derived flag" if a.derived != b.derived else "",
                "blank_subtracted flag"
                if a.blank_subtracted != b.blank_subtracted else "",
                f"temperatures {len(a.temperature_c)} vs {len(b.temperature_c)}"
                if list(a.temperature_c) != list(b.temperature_c) else "",
            ]
            problems = [m for m in [problem, *meta] if m]
            if problems:
                ok = False
                for message in problems:
                    print(f"  FAIL {path.name} {a.channel}: {message}")
            else:
                print(f"  ok   {path.name[:44]:46s} {a.channel:12s} "
                      f"{a.data.shape[0]}x{a.data.shape[1]} exact")
    return ok


# --------------------------------------------------------------------------- reproduce


def _load_pipeline():
    """Import ``run_sensor_characterisation.py`` as a module.

    The analysis is not reimplemented here. This script's whole job is to hand that
    module a different reader, so the numbers below come from the same ``collect()``
    that produced the committed table.
    """
    script = pathlib.Path(__file__).resolve().parent / "run_sensor_characterisation.py"
    spec = importlib.util.spec_from_file_location("run_sensor_characterisation", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _install_csv_readers(module, src: pathlib.Path) -> list[pathlib.Path]:
    """Point the pipeline's two readers at the committed text.

    ``collect()`` reaches the workbooks through exactly two module-level names, and
    swapping them is the smallest seam available without editing ``src/``: the plate
    handling, the blank logic, the layout recovery and the inversion all stay put. The
    paths handed back do not exist on disk and are never opened -- ``collect()`` uses
    them only for ``.name``, ``.stem`` and the logbook lookups keyed on the filename.
    """
    return replay.install(module, src)


def _banner(source: str, detail: str) -> None:
    print("=" * 92)
    print(f"SOURCE: {source}")
    print(f"        {detail}")
    print("=" * 92)


def compare_to_committed(derived: pd.DataFrame) -> None:
    """Say whether this run reproduces the tracked table, and how far off if not."""
    committed = REPO_ROOT / "outputs" / "sensor_characterisation.csv"
    if not committed.exists():
        print(f"\n  {committed} is not present; nothing to compare against")
        return
    reference = pd.read_csv(committed)
    print(f"\n  against {committed.relative_to(REPO_ROOT)} "
          f"({len(reference)} rows committed, {len(derived)} derived)")
    if len(reference) != len(derived) or list(reference.columns) != list(derived.columns):
        print("  DIFFERS in shape or columns -- not the same table")
        return
    keys = ["plate", "construct", "stressor", "dose_mM", "well"]
    merged = reference.merge(derived, on=keys, suffixes=("_ref", "_new"), how="outer")
    if len(merged) != len(reference):
        print("  DIFFERS: the key set does not line up")
        return
    worst, worst_col = 0.0, ""
    for column in [c for c in reference.columns if c not in keys]:
        delta = np.abs(merged[f"{column}_ref"] - merged[f"{column}_new"])
        finite = delta[np.isfinite(delta)]
        if len(finite) and float(finite.max()) > worst:
            worst, worst_col = float(finite.max()), column
    if worst == 0.0:
        print("  IDENTICAL: every value equal to the last bit")
    else:
        print(f"  largest absolute difference {worst:.3g} in {worst_col}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--export", action="store_true",
                        help="rebuild data/plates from the raw exports (needs YSTWIN_PLATES)")
    parser.add_argument("--verify", action="store_true",
                        help="check the committed text against the raw exports and stop")
    parser.add_argument("--source", choices=("auto", "raw", "committed"), default="auto",
                        help="which plate source to reproduce from (default: auto)")
    parser.add_argument("--write", action="store_true",
                        help="write the outputs/ tables as well as printing them")
    parser.add_argument("--window", action="store_true",
                        help="also sweep the late window, which re-reads every plate six "
                             "times and reproduces late_window_sensitivity.csv")
    parser.add_argument("--dir", type=pathlib.Path, default=PLATES_DIR,
                        help=f"committed plate text (default: {PLATES_DIR})")
    args = parser.parse_args()

    raw_dir = paths.biosensor_plates()
    raw_exports = sorted(raw_dir.glob("*.xlsx")) if raw_dir else []

    # THE INSTRUMENT FILES JOIN THE NEWPROTOCOL SET, because they are the same plates.
    # Named individually rather than globbed -- see INSTRUMENT_FILES for why -- and a name
    # that is not on this machine is reported rather than raised: `export()` is partial-safe
    # and a missing source leaves that plate's committed text exactly as it found it.
    xpt_dir = paths.gen5_xpt_dir()
    instrument = [xpt_dir / name for name in sorted(INSTRUMENT_FILES)] if xpt_dir else []
    absent = [path.name for path in instrument if not path.exists()]
    instrument = [path for path in instrument if path.exists()]
    if absent:
        print(f"  instrument files named but not under {xpt_dir}: {', '.join(absent)}")
    raw_exports += instrument

    # EXPORT COVERS EVERY RAW PLATE SET; REPRODUCTION COVERS ONLY THE NEWPROTOCOL ONE, and
    # the two lists are kept apart deliberately.
    #
    # Both sets have the same problem and deserve the same answer: the July workbooks carry
    # a named private individual in `cp:lastModifiedBy` exactly as the NewProtocol ones do,
    # so neither can be committed and both can have their numbers committed. Leaving the
    # July pair out was what made `paths.igem_results()` resolve to a sibling directory by
    # convention, which `audit_reproducibility.py` failed on for months and was right to.
    #
    # What must NOT happen is the two sets merging into one dose ladder. `collect()` below
    # is the sensor-characterisation pipeline and it is a statement about the NewProtocol
    # replicates; handing it six exports instead of four would silently change every fold
    # in `outputs/`. So `raw_exports` stays as it was and `exportable` is a separate list
    # that only `--export` and `--verify` read. `export()` is already partial-safe -- it
    # rewrites the rows for the exports in hand and keeps every other -- so committing one
    # set does not disturb the other.
    exportable = list(raw_exports)
    source_set = {path: SOURCE_SETS["biosensor_plates"] for path in raw_exports}
    source_set.update({path: INSTRUMENT_FILES[path.name] for path in instrument})
    july_dir = paths.igem_results()
    if july_dir is not None:
        july = sorted(july_dir.glob("*.xlsx"))
        exportable += july
        source_set.update({path: SOURCE_SETS["igem_results"] for path in july})

    if args.export or args.verify:
        if not exportable:
            raise SystemExit(
                "no raw exports to build or check the committed text against. Set "
                "YSTWIN_PLATES to the NewProtocol replicates, or YSTWIN_IGEM_RESULTS to "
                "the two 2026-07 ER exports, or both."
            )
        module = _load_pipeline()
        read_run = _by_suffix(module.read_synergy_kinetic, gen5.read_xpt_run)
        read_doses = _by_suffix(module.read_dose_response, _no_dose_sheet)
        if args.export:
            manifest = export(exportable, args.dir, read_run, read_doses, source_set)
            written = sorted(args.dir.glob("*.csv"))
            total = sum(p.stat().st_size for p in written)
            source_total = int(manifest.drop_duplicates("export").source_bytes.astype(int).sum())
            shown = (args.dir.relative_to(REPO_ROOT)
                     if args.dir.resolve().is_relative_to(REPO_ROOT) else args.dir)
            print(f"wrote {len(written)} files, {total / 1024:.0f} kB, into {shown}")
            print(f"  from {manifest.export.nunique()} exports totalling "
                  f"{source_total / 1024:.0f} kB of .xlsx "
                  f"({100 * total / source_total:.0f}% of the binary)")
        print("\nverifying the committed text against the workbooks, value by value:")
        if not verify(exportable, args.dir, read_run):
            raise SystemExit("\nthe committed text does not match the workbooks; "
                             "do not commit this export")
        print("\nevery block round-trips exactly.")
        return

    use_raw = args.source == "raw" or (args.source == "auto" and bool(raw_exports))
    if args.source == "raw" and not raw_exports:
        raise SystemExit("--source raw was asked for and YSTWIN_PLATES resolves to nothing")

    module = _load_pipeline()
    if use_raw:
        # The same dispatch the export path uses, so "reproduce from raw" covers the same
        # plates as "reproduce from the committed text" rather than a subset of them.
        module.read_synergy_kinetic = _by_suffix(module.read_synergy_kinetic,
                                                 gen5.read_xpt_run)
        module.read_dose_response = _by_suffix(module.read_dose_response, _no_dose_sheet)
        _banner("the raw plate files (a fresh derivation)",
                f"{len(raw_exports) - len(instrument)} .xlsx under {raw_dir}"
                + (f" and {len(instrument)} .xpt under {xpt_dir}" if instrument else ""))
        exports = raw_exports
    else:
        exports = _install_csv_readers(module, args.dir)
        _banner("the committed text in data/plates (NOT the workbooks)",
                f"{len(exports)} exports rebuilt from CSV; the workbooks are absent, so "
                f"this is a replay of committed measurements")

    frame, skipped = module.collect(exports)
    if frame.empty:
        raise SystemExit("no export yielded a usable dose ladder")

    module.report_dose_response(frame)
    sweep = module.report_autofluorescence(frame)
    # The window sweep re-runs collect() at six late-window fractions, so it reads the
    # plates six more times. It works through the committed text like everything else,
    # which is the point: the choice of window is the largest unforced assumption in the
    # dose response, and a clone that cannot sweep it has to take that choice on trust.
    window = module.report_late_window_sensitivity(exports) if args.window else pd.DataFrame()
    compare_to_committed(frame)

    if args.write:
        out = paths.outputs_dir()
        frame.to_csv(out / "sensor_characterisation.csv", index=False)
        sweep.to_csv(out / "autofluorescence_sensitivity.csv", index=False)
        print(f"\n  wrote {out / 'sensor_characterisation.csv'} ({len(frame)} well rows)")
        print(f"  wrote {out / 'autofluorescence_sensitivity.csv'} ({len(sweep)} rows)")
        if not window.empty:
            window.to_csv(out / "late_window_sensitivity.csv", index=False)
            print(f"  wrote {out / 'late_window_sensitivity.csv'} ({len(window)} rows)")
    else:
        print("\n  nothing written; pass --write to update outputs/")

    for reason, names in skipped.items():
        for name in names:
            print(f"  dropped: {name[:56]:<58} ({reason})")
    print()
    _banner("the raw plate files (a fresh derivation)" if use_raw
            else "the committed text in data/plates (NOT the workbooks)",
            "repeated here so it cannot be lost above the table")


if __name__ == "__main__":
    main()
