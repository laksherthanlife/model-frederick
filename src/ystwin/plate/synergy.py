"""Reader for BioTek/Agilent Synergy H1 kinetic plate exports.

The exporter writes one block per optical channel, each with a ``Time`` column, a
temperature column whose header also carries the channel name and optics
(``T<deg> mCitrine:480,530``), then one column per well. Blocks may sit behind a
variable preamble of instrument metadata, so the header row is located by content
rather than by a fixed offset.

Some exports split the plate down the page instead of across it: one block per group
of plate columns, stacked vertically, and only the first carries the ``Time`` column.
Two of the four biosensor replicates are written that way, and reading only the first
block silently returns 21 wells of 87 -- which then looks like a plate with no blank.
``_continuations`` reattaches those groups to the time axis of the block above them,
and refuses rather than guesses when the row counts disagree.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import io
import pathlib
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from .layout import RecordedPlate, recorded_for_export, recorded_well_roles

WELL_RE = re.compile(r"^([A-H])(\d{1,2})$")
_TIME_COL = "time"

# Sheet-name spellings mapped to one channel each.
_SHEET_CHANNEL_ALIASES = (
    ("od600", "OD600"),
    ("od", "OD600"),
    ("absorbance", "OD600"),
    ("mcitrine", "mCitrine"),
    ("citrine", "mCitrine"),
    ("gfp", "GFP"),
    ("rfp", "RFP"),
    ("mcherry", "mCherry"),
    ("yfp", "YFP"),
)


def _channel_from_sheet_name(sheet: str) -> str | None:
    """``'OD600 (Raw)'`` -> ``'OD600'``; ``'NEW Raw mCitrine'`` -> ``'mCitrine'``.

    Matched longest-alias-first so that ``od600`` wins over the bare ``od``.
    """
    text = re.sub(r"[^a-z0-9]+", " ", str(sheet).lower())
    tokens = set(text.split())
    for alias, channel in sorted(_SHEET_CHANNEL_ALIASES, key=lambda kv: -len(kv[0])):
        if alias in tokens:
            return channel
    return None


def _is_well(label: object) -> bool:
    return isinstance(label, str) and WELL_RE.match(label.strip()) is not None


def _parse_channel_header(header: object) -> tuple[str, str] | None:
    """``'T° mCitrine:480,530[2]'`` -> ``('mCitrine', '480,530')``.

    The trailing ``[n]`` marks a repeat read of the same fluorophore -- usually a
    second gain -- so it is stripped from the optics here and restored as a
    distinguishing suffix on the channel label by the caller.
    """
    if not isinstance(header, str):
        return None
    raw = header.strip()
    if raw.lower() == _TIME_COL:
        return None
    # Strip the T prefix only when a degree sign follows it.
    text = re.sub(r"^T\s*[^A-Za-z0-9:\s]\s*", "", raw)
    name, _, optics = text.partition(":")
    name = name.strip()
    if not name or name.lower() == _TIME_COL:
        return None
    optics = re.sub(r"\[\d+\]\s*$", "", optics.strip()).strip()
    return name, optics


def _to_hours(value: object) -> float:
    """Elapsed-time cell -> hours. Handles ``HH:MM:SS``, ``time``, ``timedelta``, day fraction."""
    if isinstance(value, _dt.timedelta):
        return value.total_seconds() / 3600.0
    if isinstance(value, _dt.time):
        return value.hour + value.minute / 60.0 + value.second / 3600.0
    if isinstance(value, _dt.datetime):
        return value.hour + value.minute / 60.0 + value.second / 3600.0
    if isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool):
        return float(value) * 24.0  # Excel serial day fraction
    if isinstance(value, str):
        parts = value.strip().split(":")
        if len(parts) == 3:
            h, m, s = (float(p) for p in parts)
            return h + m / 60.0 + s / 3600.0
        if len(parts) == 2:
            m, s = (float(p) for p in parts)
            return m / 60.0 + s / 3600.0
    raise ValueError(f"unparseable elapsed-time cell: {value!r}")


@dataclass(frozen=True)
class KineticBlock:
    """One optical channel's kinetic trace across wells."""

    channel: str  # unique label: 'mCitrine' when read once, 'mCitrine[1]'/'[2]' when repeated
    fluorophore: str  # base name as exported, without any repeat-read suffix
    optics: str  # excitation,emission or absorbance wavelength, for calibration provenance
    sheet: str
    derived: bool  # True when summary columns sit among the wells (a "for plotting" sheet)
    blank_subtracted: bool | None  # declared processing state; None means unknown
    data: pd.DataFrame  # index = elapsed hours, columns = well labels
    temperature_c: list[float]

    @property
    def times_h(self) -> list[float]:
        return [float(t) for t in self.data.index]

    @property
    def wells(self) -> list[str]:
        return list(self.data.columns)


@dataclass(frozen=True)
class SynergyRun:
    """All kinetic channels recovered from one exported workbook."""

    source: pathlib.Path
    blocks: tuple[KineticBlock, ...]
    recorded_plate: RecordedPlate | None = None

    def __post_init__(self) -> None:
        source = pathlib.Path(self.source)
        object.__setattr__(self, "source", source)
        recorded = (recorded_for_export(source.name) if self.recorded_plate is None
                    else self.recorded_plate)
        object.__setattr__(self, "recorded_plate", recorded)
        if recorded is not None:
            blocks = []
            for block in self.blocks:
                data = block.data.copy(deep=False)
                data.attrs["well_roles"] = recorded.well_roles
                blocks.append(replace(block, data=data))
            object.__setattr__(self, "blocks", tuple(blocks))

    @property
    def channel_names(self) -> list[str]:
        return [b.channel for b in self.blocks]

    def channel(self, name: str) -> KineticBlock:
        """Look up one channel by its unique label.

        A bare fluorophore name is rejected when the run holds several reads of it:
        those differ in gain, and picking one silently would let an optics change
        pass as biology.
        """
        exact = [b for b in self.blocks if b.channel == name]
        if len(exact) == 1:
            return exact[0]
        repeats = [b.channel for b in self.blocks if b.fluorophore == name]
        if len(repeats) > 1:
            raise KeyError(
                f"{name!r} was read {len(repeats)} times in this run "
                f"({', '.join(repeats)}); ask for one explicitly -- they differ in gain"
            )
        raise KeyError(f"no channel {name!r}; have {self.channel_names}")

    def raw_channel(self, fluorophore: str) -> KineticBlock:
        """The uniquely declared unblanked block for ``fluorophore``.

        Raises:
            KeyError: if every block for the fluorophore is blank-subtracted, or if
                the fluorophore is absent.
            ValueError: if processing is undeclared or more than one block is raw.
        """
        candidates = [b for b in self.blocks if b.fluorophore == fluorophore]
        if not candidates:
            raise KeyError(
                f"no channel for {fluorophore!r}; have "
                f"{sorted({b.fluorophore for b in self.blocks})}"
            )
        undeclared = [b.channel for b in candidates
                      if not isinstance(b.blank_subtracted, (bool, np.bool_))]
        if undeclared:
            raise ValueError(
                f"unknown or ambiguous correction-state for {undeclared}; supply an explicit "
                "boolean processing_states declaration before selecting a raw channel"
            )
        raw = [b for b in candidates if not b.blank_subtracted]
        if not raw:
            raise KeyError(
                f"every {fluorophore!r} block in this workbook is already "
                "blank-subtracted; use it directly with od_blank=0 rather than "
                "subtracting a blank a second time"
            )
        if len(raw) != 1:
            raise ValueError(
                f"ambiguous raw {fluorophore!r} channels: {[b.channel for b in raw]}; "
                "well count, block order and sheet name cannot establish correction "
                "state or which instrument read is intended"
            )
        return raw[0]

    def aligned(self, reference: str) -> pd.DataFrame:
        """Interpolate every channel onto the reference channel's time grid.

        Channels are read sequentially, so their timestamps are offset by tens of
        seconds. Downstream state estimation needs one clock; linear interpolation
        over an interval far shorter than the sampling period is the honest choice.
        """
        grid = np.asarray(self.channel(reference).data.index, dtype=float)
        frames = {}
        for block in self.blocks:
            src_t = np.asarray(block.data.index, dtype=float)
            for well in block.wells:
                y = np.asarray(block.data[well], dtype=float)
                frames[(block.channel, well)] = np.interp(grid, src_t, y)
        out = pd.DataFrame(frames, index=pd.Index(grid, name="time_h"))
        out.columns = pd.MultiIndex.from_tuples(out.columns, names=["channel", "well"])
        roles = recorded_well_roles(*(block.data for block in self.blocks))
        if roles is not None:
            out.attrs["well_roles"] = roles
        return out.sort_index()


def _well_header_rows(raw: pd.DataFrame) -> list[tuple[int, list[int], bool]]:
    """Every row that heads a well block, with its well columns and whether it has Time."""
    n_rows, n_cols = raw.shape
    found = []
    for r in range(n_rows):
        cols = [c for c in range(n_cols) if _is_well(raw.iat[r, c])]
        if not cols:
            continue
        first = raw.iat[r, 0]
        has_time = isinstance(first, str) and first.strip().lower() == _TIME_COL
        found.append((r, cols, has_time))
    return found


def _continuations(raw: pd.DataFrame, start: int, n_times: int,
                   taken: set[str]) -> list[tuple[list[str], list[list[float]]]]:
    """Well groups stacked below a time-headed block, sharing its time axis.

    A continuation has well headers but no ``Time`` cell, and its rows carry no time of
    their own, so they are read positionally against the block above. Two conditions have
    to hold or the group is skipped rather than merged: its wells must be new -- an
    overlapping group is the raw/blank-subtracted pairing, a different thing entirely --
    and it must carry exactly as many rows as the block it would join. Merging a shorter
    group would pair readings with the wrong timepoints, which no downstream check would
    catch.
    """
    out = []
    for row, cols, has_time in _well_header_rows(raw):
        if row <= start or has_time:
            continue
        names = [str(raw.iat[row, c]).strip() for c in cols]
        if set(names) & taken:
            break
        values = []
        for rr in range(row + 1, min(row + 1 + n_times, raw.shape[0])):
            vals = [pd.to_numeric(raw.iat[rr, c], errors="coerce") for c in cols]
            if all(pd.isna(v) for v in vals):
                break
            values.append(vals)
        if len(values) != n_times:
            continue
        out.append((names, values))
        taken |= set(names)
    return out


def _extract_blocks(raw: pd.DataFrame, sheet: str) -> list[KineticBlock]:
    blocks: list[KineticBlock] = []
    n_rows, n_cols = raw.shape
    for r in range(n_rows):
        first = raw.iat[r, 0]
        if not (isinstance(first, str) and first.strip().lower() == _TIME_COL):
            continue
        header = [raw.iat[r, c] for c in range(n_cols)]
        well_cols = [c for c in range(n_cols) if _is_well(header[c])]
        if not well_cols:
            continue
        temp_cols = [c for c in range(1, min(well_cols)) if _parse_channel_header(header[c])]
        if temp_cols:
            temp_col = temp_cols[0]
            fluorophore, optics = _parse_channel_header(header[temp_col])
        else:
            # No channel column: take the name from the sheet, or refuse.
            temp_col = None
            fluorophore = _channel_from_sheet_name(sheet)
            if fluorophore is None:
                continue
            optics = ""

        times, temps, rows = [], [], []
        for rr in range(r + 1, n_rows):
            cell = raw.iat[rr, 0]
            if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                break
            if isinstance(cell, str) and cell.strip().lower() == _TIME_COL:
                break
            try:
                t = _to_hours(cell)
            except ValueError:
                break
            values = [pd.to_numeric(raw.iat[rr, c], errors="coerce") for c in well_cols]
            # Padding rows carry a placeholder time and no wells; key off the absent wells.
            if all(pd.isna(v) for v in values):
                break
            times.append(t)
            if temp_col is not None:
                temps.append(pd.to_numeric(raw.iat[rr, temp_col], errors="coerce"))
            rows.append(values)
        if not times:
            continue
        frame = pd.DataFrame(
            rows,
            index=pd.Index(times, name="time_h"),
            columns=[str(header[c]).strip() for c in well_cols],
            dtype=float,
        )
        # Reattach any column groups stacked below this block on the same time axis.
        for names, values in _continuations(raw, r, len(times), set(frame.columns)):
            frame = pd.concat(
                [frame, pd.DataFrame(values, index=frame.index, columns=names, dtype=float)],
                axis=1)
        # Summary columns among the wells mark a derived, for-plotting sheet.
        span = range(min(well_cols), max(well_cols) + 1)
        derived = any(c not in well_cols for c in span)

        blocks.append(
            KineticBlock(
                channel=fluorophore,  # provisional; disambiguated once all blocks are read
                fluorophore=fluorophore,
                optics=optics,
                sheet=sheet,
                derived=derived,
                blank_subtracted=None,
                data=frame,
                temperature_c=[float(t) for t in temps],
            )
        )
    return blocks


def _declare_processing_states(blocks: list[KineticBlock],
                               states: Mapping[str, bool] | None) -> list[KineticBlock]:
    if states is None:
        return blocks
    if not isinstance(states, Mapping):
        raise ValueError("processing_states must map channel labels to boolean correction-state declarations")
    unknown = set(states) - {block.channel for block in blocks}
    if unknown:
        raise ValueError(f"processing_states names unknown channels: {sorted(unknown, key=str)}")
    out = []
    for block in blocks:
        if block.channel not in states:
            out.append(block)
            continue
        state = states[block.channel]
        if not isinstance(state, (bool, np.bool_)):
            raise ValueError(f"correction-state declaration must be boolean for {block.channel!r}")
        if block.blank_subtracted is not None and block.blank_subtracted != state:
            raise ValueError(f"correction-state declaration conflicts with source metadata for {block.channel!r}")
        out.append(replace(block, blank_subtracted=bool(state)))
    return out


def read_synergy_kinetic(path: str | pathlib.Path, prefer_raw: bool = True, *,
                         recorded_plate: RecordedPlate | None = None,
                         processing_states: Mapping[str, bool] | None = None) -> SynergyRun:
    """Read every kinetic channel block in a Synergy H1 workbook.

    Args:
        path: Workbook to read.
        prefer_raw: Drop derived ("for plotting") blocks for any channel that also
            has a non-derived block. This does not declare blank-subtraction state.
            A channel whose only block is derived is kept either way.
        processing_states: Explicit ``channel label -> blank_subtracted`` booleans.
            Otherwise the committed manifest's declarations are used only for an exact
            source SHA256 and block locator match. Unmatched blocks remain unknown and
            cannot be selected by ``raw_channel``. Neither reading signs nor sheet names
            establish processing; explicit declarations cannot override known metadata.
    """
    from .replay import _processing_states_for_source

    path = pathlib.Path(path)
    contents = path.read_bytes()
    sheets = pd.read_excel(io.BytesIO(contents), sheet_name=None, header=None)
    blocks: list[KineticBlock] = []
    for sheet, raw in sheets.items():
        if raw.empty:
            continue
        blocks.extend(_extract_blocks(raw, sheet))
    if not blocks:
        raise ValueError(f"no kinetic channel blocks found in {path}")
    if prefer_raw:
        has_nonderived = {b.fluorophore for b in blocks if not b.derived}
        blocks = [b for b in blocks if not b.derived or b.fluorophore not in has_nonderived]
    declared = _processing_states_for_source(hashlib.sha256(contents).hexdigest())
    locators = [(b.sheet, b.fluorophore, b.optics) for b in blocks]
    blocks = [replace(b, blank_subtracted=declared.get(locator)
                      if locators.count(locator) == 1 else None)
              for b, locator in zip(blocks, locators)]
    blocks = _declare_processing_states(label_repeat_reads(blocks), processing_states)
    return SynergyRun(source=path, blocks=tuple(blocks), recorded_plate=recorded_plate)


def label_repeat_reads(blocks: list[KineticBlock]) -> list[KineticBlock]:
    """Suffix repeated reads of one fluorophore as ``name[1]``, ``name[2]``, ...

    Public because the ``.xlsx`` is no longer the only door plate data comes through:
    ``plate/gen5.py`` builds the same :class:`KineticBlock` tuple from the instrument's own
    ``.xpt``, and the two must label a repeat read identically or the committed text of one
    source would not compare against the other. One convention, in one place, rather than a
    second copy of it that drifts.
    """
    counts: dict[str, int] = {}
    for block in blocks:
        counts[block.fluorophore] = counts.get(block.fluorophore, 0) + 1
    seen: dict[str, int] = {}
    out = []
    for block in blocks:
        if counts[block.fluorophore] == 1:
            out.append(block)
            continue
        seen[block.fluorophore] = seen.get(block.fluorophore, 0) + 1
        out.append(replace(block, channel=f"{block.fluorophore}[{seen[block.fluorophore]}]"))
    return out
