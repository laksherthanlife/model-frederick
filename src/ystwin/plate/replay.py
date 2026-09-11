"""Replay the committed plate text in place of the Synergy workbooks.

The four NewProtocol exports are `.xlsx` files that are **not** in this repository and
cannot be: their document properties carry a named private individual. What is committed
is `data/plates/` -- the same numbers as text, at a precision that round-trips bit for bit,
with a manifest saying which file was which block of which export.

That makes the measurements reproducible by anyone who clones this. What it did not do,
until this module existed, was make the *scripts* reproducible: the replay lived inside
``scripts/plate_readings.py`` as a private shim, so that one script could run from the
committed text and every other script that reads a plate still demanded the workbooks and
exited when they were absent. A stranger could reproduce one table.

So the seam moved here. :func:`install` swaps the two module-level readers that every
plate-reading script reaches the workbooks through, and the rest of the pipeline -- plate
handling, blanks, layout recovery, the inversion -- is untouched and unaware. The paths
handed back do not exist on disk and are never opened; ``collect()`` uses them only for
``.name``, ``.stem`` and the logbook lookups keyed on the filename.

The replay is a replay and says so. Nothing here re-derives anything from the workbooks;
it reconstitutes what a verified export of them produced. ``plate_readings.py --export``
is what checks the two against each other, value by value.
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping

import pandas as pd

from .. import paths
from .layout import RecordedPlate
from .synergy import KineticBlock, SynergyRun, _declare_processing_states

__all__ = [
    "MANIFEST",
    "committed_dir",
    "exports",
    "install",
    "load_doses",
    "load_manifest",
    "load_run",
]

MANIFEST = "manifest.csv"


def committed_dir() -> pathlib.Path:
    """Where the committed plate text lives. Tracked, so it is always present."""
    return paths.data_dir() / "plates"


def _from_hms(text: str) -> float:
    """``H:MM:SS`` -> hours, using ``plate/synergy.py::_to_hours``'s own expression.

    The expression matters more than the value. ``h + m / 60 + s / 3600`` and
    ``round(h * 3600 + m * 60 + s) / 3600`` are the same real number and different
    floats, and the committed table is reproduced byte for byte or not at all.
    """
    h, m, s = (float(part) for part in text.split(":"))
    return h + m / 60.0 + s / 3600.0


def load_manifest(src: pathlib.Path | None = None) -> pd.DataFrame:
    """The committed manifest, or a refusal naming both ways out.

    Read entirely as strings. ``optics`` holds ``600``, which is a wavelength label rather
    than a number, and an inferred int would not compare equal to the block it came from.
    """
    src = committed_dir() if src is None else src
    path = src / MANIFEST
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. Either set YSTWIN_PLATES to the raw Synergy exports, "
            f"or rebuild the committed text with scripts/plate_readings.py --export on a "
            f"machine that has them."
        )
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def exports(src: pathlib.Path | None = None,
            manifest: pd.DataFrame | None = None,
            source_set: str | None = None) -> list[pathlib.Path]:
    """The export names the committed text covers, in manifest order.

    These are names, not files. Nothing opens them.

    Args:
        src: Where the committed text lives.
        manifest: A manifest already loaded, to avoid re-reading it.
        source_set: Restrict to one raw plate set -- ``"newprotocol"`` or
            ``"july_2026_07"``. ``None`` returns every export.

            The filter exists because two scripts read different halves. `run_gates.py`
            wants both, and the comparison between them is its whole point. `run_d2.py`
            wanted only the July pair, and until both sets were committed it got that for
            free by resolving one directory. Once the replay covered everything, "whatever
            the resolver found" silently became "all seven", which would have written five
            new tracked tables nobody asked for.

            Keyed on a manifest column rather than on the date in the filename: the set an
            export belongs to is a fact about where it came from, `plate_readings.py`
            records it at export time when that is still known, and a filename heuristic
            would be wrong the first time a plate is renamed.
    """
    manifest = load_manifest(src) if manifest is None else manifest
    if source_set is not None:
        if "source_set" not in manifest.columns:
            raise ValueError(
                "this manifest predates the source_set column, so it cannot say which raw "
                "plate set an export came from. Re-run `plate_readings.py --export`")
        manifest = manifest[manifest.source_set == source_set]
    return [pathlib.Path(name) for name in manifest.export.drop_duplicates()]


def _manifest_boolean(value, field: str) -> bool:
    if type(value) is bool:
        return value
    if value not in ("True", "False"):
        raise ValueError(f"manifest {field} must be a boolean True or False, got {value!r}")
    return value == "True"


def _processing_states_for_source(source_sha256: str) -> dict[tuple[str, str, str], bool]:
    """Declared states by exact source hash and block locator, never by naming patterns."""
    try:
        manifest = load_manifest()
    except FileNotFoundError:
        return {}
    rows = manifest[manifest.source_sha256 == source_sha256]
    locator = ["sheet", "fluorophore", "optics"]
    if rows.duplicated(locator).any():
        raise ValueError(f"ambiguous correction-state block locators for source SHA256 {source_sha256}")
    return {(row.sheet, row.fluorophore, row.optics):
            _manifest_boolean(row.blank_subtracted, "blank_subtracted")
            for row in rows.itertuples()}


def load_run(export_name: str, src: pathlib.Path | None = None,
             manifest: pd.DataFrame | None = None, *,
             recorded_plate: RecordedPlate | None = None, prefer_raw: bool = True,
             processing_states: Mapping[str, bool] | None = None) -> SynergyRun:
    """Rebuild one export's :class:`SynergyRun` from the committed text.

    Block order preserves exported channel identities and the traversal in ``aligned``;
    it never resolves ambiguous raw channels. Optional processing declarations must agree
    with the manifest, not override it.
    """
    if not prefer_raw:
        raise ValueError("prefer_raw=False is unavailable in committed replay; derived blocks were not exported")
    src = committed_dir() if src is None else src
    manifest = load_manifest(src) if manifest is None else manifest
    blocks = []
    for _, row in manifest[manifest.export == export_name].iterrows():
        frame = pd.read_csv(src / row.readings_file)
        times = [_from_hms(value) for value in frame.elapsed_hms]
        temperatures: list[float] = []
        if "temperature_c" in frame.columns:
            temperatures = [float(v) for v in frame.temperature_c]
        wells = [c for c in frame.columns if c not in {"elapsed_hms", "temperature_c"}]
        data = frame[wells].astype(float)
        data.index = pd.Index(times, name="time_h")  # the name synergy.py gives it
        blocks.append(KineticBlock(
            channel=row.channel,
            fluorophore=row.fluorophore,
            optics=row.optics,
            sheet=row.sheet,
            derived=_manifest_boolean(row.derived, "derived"),
            blank_subtracted=_manifest_boolean(row.blank_subtracted, "blank_subtracted"),
            data=data,
            temperature_c=temperatures,
        ))
    if not blocks:
        raise ValueError(f"no committed blocks for {export_name!r}")
    blocks = _declare_processing_states(blocks, processing_states)
    return SynergyRun(source=pathlib.Path(export_name), blocks=tuple(blocks),
                      recorded_plate=recorded_plate)


def load_doses(export_name: str, src: pathlib.Path | None = None,
               manifest: pd.DataFrame | None = None) -> pd.DataFrame:
    """Rebuild one export's per-construct sheets, or refuse exactly as the reader does.

    ``collect()`` distinguishes "this plate has no dose ladder of its own" from "this plate
    is unusable" by catching :class:`ValueError` from ``read_dose_response``, so the
    absence has to be raised and not returned empty.
    """
    src = committed_dir() if src is None else src
    manifest = load_manifest(src) if manifest is None else manifest
    files = {f for f in manifest[manifest.export == export_name].dose_response_file if f}
    if not files:
        raise ValueError(f"no per-construct dose-response sheet committed for {export_name}")
    if len(files) != 1:
        raise ValueError(f"conflicting dose-response files for {export_name!r}: {sorted(files)}")
    frame = pd.read_csv(src / next(iter(files)))
    return pd.DataFrame({
        "construct": frame.construct,
        "dose_mM": frame.dose_mM.astype(float),
        "time_h": frame.elapsed_min.astype(float) / 60.0,
        "signal": frame.signal.astype(float),
    })


def install(target, src: pathlib.Path | None = None) -> list[pathlib.Path]:
    """Point one module's plate readers at the committed text; return the export names.

    Two names are swapped, because two are all any script reaches the workbooks through.
    A module that binds only one of them keeps the other, which is why this rebinds what is
    there rather than asserting both are: ``run_calibration_nis.py`` reads kinetics and
    never reads a dose sheet.

    Args:
        target: Either a module object, or the caller's ``globals()``. **Pass
            ``globals()`` from a script.** ``sys.modules[__name__]`` looks equivalent and
            is not: a script loaded through ``importlib.util.spec_from_file_location`` --
            which is how every script test in this repository loads one -- is never
            registered in ``sys.modules``, so that lookup raises ``KeyError`` under test
            and works when run from the command line. A rebinding that depends on how the
            module was loaded is a rebinding that is tested in a different world from the
            one it runs in.
        src: Where the committed text lives. Defaults to ``data/plates``.

    Returns:
        The export names covered, for a caller to iterate in place of a directory listing.
    """
    src = committed_dir() if src is None else src
    manifest = load_manifest(src)
    covered = set(manifest.export)

    def _dispatch(committed_reader, original):
        """Replay only the exports the manifest covers; hand the rest to the real reader.

        ``run_gates.py`` reads two plate sets -- the July Synergy exports and the
        NewProtocol replicates -- and only the second is committed as text. A swap that
        redirected *everything* made the first fail with "no committed blocks", turning a
        fix for one source into a break in the other. Whether the manifest covers a name is
        the authoritative test, and it is the one thing that actually distinguishes them.
        """
        def read(path, *args, **kwargs):
            name = pathlib.Path(path).name
            if name in covered:
                return committed_reader(name, *args, **kwargs)
            if original is None:
                raise ValueError(f"no committed blocks for {name!r}")
            return original(path, *args, **kwargs)
        return read

    readers = {
        "read_synergy_kinetic": lambda name, prefer_raw=True, **kwargs: load_run(
            name, src, manifest, prefer_raw=prefer_raw, **kwargs),
        "read_dose_response": lambda name: load_doses(name, src, manifest),
    }
    for name, committed_reader in readers.items():
        if isinstance(target, dict):
            if name in target:
                target[name] = _dispatch(committed_reader, target[name])
        elif hasattr(target, name):
            setattr(target, name, _dispatch(committed_reader, getattr(target, name)))
    return exports(src, manifest)
