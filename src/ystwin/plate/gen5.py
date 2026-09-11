"""Reader for BioTek/Agilent Gen5 ``.xpt`` experiment files -- the reader's own format.

Why this exists
---------------
``synergy.py`` reads what Gen5 *exports*: an ``.xlsx`` whose numbers are whatever the
operator had on screen when they pressed Export. That is lossy in two ways that have
already cost this project data.

1. Absorbance is exported at three decimals. The instrument records four.
2. If the operator had a blank-subtraction transform selected, the export carries the
   transformed numbers and the raw ones are gone. Plate ``20260804`` reached
   ``data/plates`` that way: in the committed text its media blanks H1-H3 span -0.004 to
   +0.004 OD where the rest of the plate spans 0.09 to 0.61, so the background could not
   be recovered from the export and the biosensor panel is ``n=3`` rather than ``n=4``
   (``docs/DATA_INVENTORY.md``, which calls re-exporting this one file the cheapest
   remaining gain in the project).

The ``.xpt`` is the experiment file Gen5 itself writes, and it holds the untransformed
readings. Nothing here re-exports anything; this module only reads.

The container
-------------
An ``.xpt`` is a Microsoft Compound File Binary document (MS-CFB, the OLE2 container
also used by ``.doc`` and ``.xls``). :class:`_CompoundFile` implements the subset needed
to pull whole streams out of one, so that reading a plate costs no dependency beyond the
standard library and numpy. ``tests/test_gen5_xpt.py`` checks it stream-for-stream
against ``olefile`` on every ``.xpt`` present, when ``olefile`` is installed.

Two streams matter:

``Contents``
    Instrument and protocol metadata. Only the protocol file name is read from it.
``SUBSETS/1/DATA``
    Every reading, as an MFC ``CArchive`` serialisation of Gen5's assay document.

The payload
-----------
The ``CArchive`` format is not documented by the vendor, so every offset below was
measured from the files in ``docs/DATA_INVENTORY.md`` and then *checked against numbers
this repository already trusts*: ``data/plates/`` holds each plate as committed text,
exported from the matching ``.xlsx`` and verified value-for-value by
``plate_readings.py --export``. Four of those plates also have their ``.xpt`` on this
machine, and reading it reproduces the committed text -- exactly for fluorescence (20,115
readings, zero disagreements), and to the 0.0005 that a three-decimal export can say
about a four-decimal reading for absorbance. Those comparisons are the tests in
``tests/test_gen5_xpt.py``, and they are the only reason to believe any of the constants
in this module.

What the archive holds, in order:

* one *plate data set* per exported channel, each beginning with a length-prefixed name
  such as ``OD600:600`` or ``mCitrine:480,530`` (with a ``[2]``/``[3]`` suffix when the
  same fluorophore was read more than once in a step, at a second gain);
* for each of those, a **well block**: a nine-field header, then
  ``n_times * n_rows * n_cols`` fixed-size records in reading order -- time-major,
  and within a timepoint row-major across the plate (A1..A12, B1..B12, ... H12) --
  followed by an **elapsed-time table** of one entry per timepoint;
* then one *temperature data set* per physical read step, named ``T<degree> `` plus the
  channel name, holding one temperature per timepoint.

Every well of the plate is present in the ``.xpt`` whether or not the operator included
it in the export: the committed 20260728 text has 87 columns, the file has 96.

What this module refuses to do
------------------------------
Guess. When a structure is present but does not validate -- a record stride that does not
hold, an elapsed table whose length disagrees with the block header -- the reader raises
:class:`Gen5FormatError` naming what failed. When a structure is legitimately absent --
an endpoint read has no elapsed axis, two 2026-06 files logged no temperature -- the
corresponding attribute is ``None`` rather than a fabricated array.
"""

from __future__ import annotations

import datetime as _dt
import pathlib
import re
import struct
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:  # `synergy` is imported inside the functions below to avoid a cycle;
    from .synergy import SynergyRun  # this makes their return annotation resolvable.

__all__ = [
    "OVERFLOW",
    "Gen5Channel",
    "Gen5FormatError",
    "Gen5Read",
    "as_synergy_run",
    "read_xpt",
    "read_xpt_run",
    "well_names",
]

OVERFLOW = -99999.0
"""What Gen5 writes in a fluorescence well whose reading passed the top of scale.

A property of the file format, which is why it lives here. `calib/gain_linearity.py` is
where the evidence for it is written down: 11,149 of the 303,552 fluorescence readings in
the 22-file corpus are exactly this, no reading anywhere lies between 99,998 and it, and
the flagged readings' surviving partners predict at least 98.5% of the 100,000 RFU ceiling.

It is a flag, not a measurement, and averaging one into a reporter trace produces a number
that looks like a reading. :meth:`Gen5Channel.frame` refuses to hand back a channel
containing it unless the caller says what to do with it.
"""


class Gen5FormatError(ValueError):
    """The file is not a Gen5 ``.xpt``, or it is one shaped in a way this reader cannot verify."""


# ---------------------------------------------------------------------------
# the container: Microsoft Compound File Binary
# ---------------------------------------------------------------------------

_CFB_SIGNATURE = bytes.fromhex("d0cf11e0a1b11ae1")
"""MS-CFB v1.0 section 2.2: the eight bytes every compound file starts with."""

_FREESECT = 0xFFFFFFFF
"""MS-CFB section 2.2: a sector-chain entry meaning "unallocated"."""

_ENDOFCHAIN = 0xFFFFFFFE
"""MS-CFB section 2.2: the last sector of a chain."""

_DIR_ENTRY_SIZE = 128
"""MS-CFB section 2.6.1. Fixed for both major versions."""

_MAX_CHAIN = 1 << 22
"""Sector-chain walk limit: ~2 GB of 512-byte sectors.

A corrupt FAT can point a chain at itself. Without a limit that is an infinite loop
inside a test suite; with one it is a :class:`Gen5FormatError` naming the stream.
"""


class _CompoundFile:
    """Whole streams out of an MS-CFB document, by path.

    Only the read side, and only what an ``.xpt`` needs: sector chains through the FAT,
    short streams through the mini FAT, and a directory walk that flattens the red-black
    sibling tree into ``'SUBSETS/1/DATA'``-style paths. No writing, no locking, no
    property sets.

    Args:
        data: The whole file. These are 0.1-2.3 MB; streaming them would buy nothing and
            would make the sector arithmetic harder to check.

    Raises:
        Gen5FormatError: if the signature, version, or sector geometry is not one this
            reader can address.
    """

    def __init__(self, data: bytes) -> None:
        if len(data) < 512 or data[:8] != _CFB_SIGNATURE:
            raise Gen5FormatError("not a compound file: signature missing")
        self._data = data
        (major,) = struct.unpack_from("<H", data, 26)
        (byte_order,) = struct.unpack_from("<H", data, 28)
        (sector_shift,) = struct.unpack_from("<H", data, 30)
        (mini_shift,) = struct.unpack_from("<H", data, 32)
        if byte_order != 0xFFFE:
            raise Gen5FormatError(f"unsupported byte order 0x{byte_order:04x}")
        if (major, sector_shift) not in ((3, 9), (4, 12)):
            raise Gen5FormatError(
                f"unsupported compound-file geometry: version {major}, sector shift "
                f"{sector_shift}")
        self._sector_size = 1 << sector_shift
        self._mini_size = 1 << mini_shift
        (self._n_fat,) = struct.unpack_from("<I", data, 44)
        (first_dir,) = struct.unpack_from("<I", data, 48)
        (self._mini_cutoff,) = struct.unpack_from("<I", data, 56)
        (first_minifat,) = struct.unpack_from("<I", data, 60)
        (n_minifat,) = struct.unpack_from("<I", data, 64)
        (first_difat,) = struct.unpack_from("<I", data, 68)
        (n_difat,) = struct.unpack_from("<I", data, 72)

        self._fat = self._read_fat(first_difat, n_difat)
        self._minifat = self._read_uint32_chain(first_minifat, n_minifat)
        self._entries = self._read_directory(first_dir)
        root = self._entries[0]
        self._mini_stream = self._read_chain(root["start"], root["size"], mini=False)
        self.paths = self._walk(root["child"], "")

    # -- sectors ------------------------------------------------------------

    def _sector(self, index: int) -> bytes:
        """MS-CFB section 2.1: sector *n* begins one sector-size in, past the header."""
        start = (index + 1) * self._sector_size
        chunk = self._data[start:start + self._sector_size]
        if len(chunk) != self._sector_size:
            raise Gen5FormatError(f"sector {index} runs past the end of the file")
        return chunk

    def _read_fat(self, first_difat: int, n_difat: int) -> list[int]:
        """The file allocation table, following the DIFAT past its first 109 entries.

        The extension loop is written from MS-CFB section 2.5 and is **not** exercised by
        anything here: the largest ``.xpt`` on this machine is 2.26 MB, which needs 35 FAT
        sectors against the 109 the header carries inline. Said rather than left implied,
        because a branch no file reaches is a branch no test covers.
        """
        difat = list(struct.unpack_from("<109I", self._data, 76))
        sector, guard = first_difat, 0
        while sector not in (_ENDOFCHAIN, _FREESECT) and guard < n_difat:
            block = self._sector(sector)
            per = self._sector_size // 4 - 1
            difat.extend(struct.unpack_from(f"<{per}I", block, 0))
            (sector,) = struct.unpack_from("<I", block, per * 4)
            guard += 1
        fat: list[int] = []
        for index in difat[:self._n_fat]:
            if index in (_FREESECT, _ENDOFCHAIN):
                continue
            block = self._sector(index)
            fat.extend(struct.unpack_from(f"<{self._sector_size // 4}I", block, 0))
        return fat

    def _read_uint32_chain(self, first: int, count: int) -> list[int]:
        """A FAT-chained run of ``count`` sectors read as little-endian uint32 (the mini FAT)."""
        out: list[int] = []
        sector, seen = first, 0
        while sector not in (_ENDOFCHAIN, _FREESECT) and seen < count:
            block = self._sector(sector)
            out.extend(struct.unpack_from(f"<{self._sector_size // 4}I", block, 0))
            sector = self._next(sector)
            seen += 1
        return out

    def _next(self, sector: int) -> int:
        if sector >= len(self._fat):
            raise Gen5FormatError(f"sector {sector} is outside the allocation table")
        return self._fat[sector]

    def _read_chain(self, start: int, size: int, *, mini: bool) -> bytes:
        """Follow one chain and return exactly ``size`` bytes of it."""
        out = bytearray()
        sector, steps = start, 0
        while sector not in (_ENDOFCHAIN, _FREESECT) and len(out) < size:
            if mini:
                offset = sector * self._mini_size
                out.extend(self._mini_stream[offset:offset + self._mini_size])
                if sector >= len(self._minifat):
                    raise Gen5FormatError(f"mini sector {sector} is outside the mini FAT")
                sector = self._minifat[sector]
            else:
                out.extend(self._sector(sector))
                sector = self._next(sector)
            steps += 1
            if steps > _MAX_CHAIN:
                raise Gen5FormatError("sector chain does not terminate")
        if len(out) < size:
            raise Gen5FormatError(f"stream ends early: wanted {size} bytes, chain gave {len(out)}")
        return bytes(out[:size])

    # -- directory ----------------------------------------------------------

    def _read_directory(self, first: int) -> list[dict]:
        raw = bytearray()
        sector, steps = first, 0
        while sector not in (_ENDOFCHAIN, _FREESECT):
            raw.extend(self._sector(sector))
            sector = self._next(sector)
            steps += 1
            if steps > _MAX_CHAIN:
                raise Gen5FormatError("directory chain does not terminate")
        entries = []
        for offset in range(0, len(raw), _DIR_ENTRY_SIZE):
            block = bytes(raw[offset:offset + _DIR_ENTRY_SIZE])
            if len(block) < _DIR_ENTRY_SIZE:
                break
            (name_len,) = struct.unpack_from("<H", block, 64)
            name = block[:max(name_len - 2, 0)].decode("utf-16-le", errors="replace")
            left, right, child = struct.unpack_from("<3i", block, 68)
            (start,) = struct.unpack_from("<I", block, 116)
            (size,) = struct.unpack_from("<Q", block, 120)
            entries.append({
                "name": name, "type": block[66], "left": left, "right": right,
                "child": child, "start": start, "size": size,
            })
        if not entries:
            raise Gen5FormatError("compound file has no directory")
        return entries

    def _walk(self, index: int, prefix: str) -> dict[str, dict]:
        """Flatten the sibling tree under ``index`` into ``path -> entry``."""
        out: dict[str, dict] = {}
        stack = [index]
        seen: set[int] = set()
        while stack:
            node = stack.pop()
            if node < 0 or node in seen or node >= len(self._entries):
                continue
            seen.add(node)
            entry = self._entries[node]
            path = f"{prefix}{entry['name']}"
            stack.extend((entry["left"], entry["right"]))
            if entry["type"] == 1:  # storage
                out.update(self._walk(entry["child"], path + "/"))
            elif entry["type"] == 2:  # stream
                out[path] = entry
        return out

    # -- public -------------------------------------------------------------

    def read(self, path: str) -> bytes:
        """One stream by path, e.g. ``'SUBSETS/1/DATA'``.

        Raises:
            Gen5FormatError: if the document has no such stream.
        """
        entry = self.paths.get(path)
        if entry is None:
            raise Gen5FormatError(
                f"no stream {path!r}; the document holds {sorted(self.paths)}")
        mini = entry["size"] < self._mini_cutoff
        return self._read_chain(entry["start"], entry["size"], mini=mini)


# ---------------------------------------------------------------------------
# the payload: Gen5's serialised assay document
# ---------------------------------------------------------------------------

_CHANNEL_ANCHOR = b"\x05\x00\x00\x00\x01\x00"
"""The six bytes that immediately precede a plate data set's name.

Measured, not documented: in all 22 ``.xpt`` files on hand the byte-length-prefixed
channel name (``\\x09OD600:600``) follows exactly this run, and the count of anchors that
resolve to a printable non-empty name equals the count of well blocks in every file.
One file (``20260713``) contains a fourth anchor whose name is empty; requiring a
non-empty printable name is what excludes it.
"""

_WELL_BLOCK_ANCHOR = b"\x18\x00\x00\x00\x07\x00"
"""Record size (24) and record kind (7) -- the header of a well-readings block.

Followed by seven uint32: ``n_times``, 1, ``n_rows``, ``n_cols``, 1, 1, 1. The three
trailing ones are constant across every file here, so they are checked rather than
interpreted.
"""

_RECORD_SIZE = 24
"""One well reading: float64 value, eight bytes of padding, an eight-byte token.

The padding is *not* a second value -- in some files it is uninitialised memory. The
token is constant within a file and differs between files, so it identifies nothing, but
it does pin the stride: :func:`_read_well_block` requires every record in a block to
carry the same token, which is what makes a wrong stride impossible to miss.
"""

_WELL_BLOCK_HEADER_SIZE = 34
"""``_WELL_BLOCK_ANCHOR`` (6) plus its seven uint32 fields (28)."""

_ELAPSED_ENTRY_SIZE = 16
"""One elapsed-time entry: uint32 zero, uint32 milliseconds, uint64 zero.

The first and last fields are zero in every kinetic file here and are checked, so a
future file that uses them fails loudly instead of being read as if they were padding.
"""

_TEMPERATURE_KIND = 2
"""The record kind of a temperature series, against 7 for well readings."""

_TEMPERATURE_RECORD_SIZE = 32
"""One temperature record: float64 Celsius, then 24 bytes this reader does not interpret.

Wider than a well record and, unlike it, with no constant token to check the stride
against. What checks it instead is the committed text: all 566 temperatures across the
four proven plates match ``data/plates`` exactly, at one decimal.
"""

_TEMPERATURE_HEADER_GAP = 53
"""Bytes between the end of a temperature data set's name and its record header.

Constant at 53 across the three name lengths the corpus contains -- ``OD600:600``,
``RFP:530,580``/``RFP:580,610`` and ``mCitrine:480,530``/``mCitrine:510,530``, so
temperature names of 12, 14 and 19 characters -- which is why this is an offset and not a
search. The header it lands on is then verified byte-for-byte before anything is read
through it, so a file where the gap differs reports no temperature rather than nonsense.
"""

_TIMESTAMP_GAP = 36
"""Bytes between the end of a channel's name and its ``__time64_t`` read-start stamp.

Three float64 (24), an eight-byte identifier, and two uint16 tags.
"""

_EPOCH_WINDOW = (_dt.datetime(1990, 1, 1, tzinfo=_dt.UTC), _dt.datetime(2100, 1, 1, tzinfo=_dt.UTC))
"""The range a decoded read-start stamp must fall in to be reported at all.

Not a plausibility filter on the science -- a filter on the decoding. A stamp outside it
means the field is not the ``__time64_t`` this reader takes it for in that file, and the
right answer is then ``None`` rather than a date.
"""

_DEGREE = "\N{DEGREE SIGN}"

_REPEAT_SUFFIX = re.compile(r"\[(\d+)\]$")
"""``mCitrine:480,530[2]`` -- the same fluorophore read again in one step, at a second gain.

Gen5 keeps one temperature series per *read step*, not per exported channel, so the
suffix is stripped when looking the temperature up. The repeat read happens at the same
instant as the first, so it is the same temperature, and this is also how the ``.xlsx``
export presents it.
"""


def well_names(rows: int, cols: int) -> tuple[str, ...]:
    """Plate well labels in the order the archive stores them: row-major, A1 first.

    Args:
        rows: Plate rows, 8 for a 96-well plate.
        cols: Plate columns, 12 for a 96-well plate.

    Returns:
        ``('A1', 'A2', ..., 'A12', 'B1', ..., 'H12')``.

    Raises:
        Gen5FormatError: for a geometry with no letter for its last row.
    """
    if not 1 <= rows <= 26 or cols < 1:
        raise Gen5FormatError(f"cannot label a {rows} x {cols} plate")
    return tuple(f"{chr(ord('A') + r)}{c + 1}" for r in range(rows) for c in range(cols))


def _cstring(buffer: bytes, offset: int) -> tuple[str, int]:
    """An MFC byte-length-prefixed string at ``offset`` -> ``(text, offset_after)``.

    Gen5 writes these in Latin-1: the temperature data sets are named with a degree sign
    (``T\N{DEGREE SIGN} OD600:600``), which is 0xB0 in Latin-1 and invalid UTF-8.
    """
    length = buffer[offset]
    end = offset + 1 + length
    return buffer[offset + 1:end].decode("latin-1"), end


def _find_all(buffer: bytes, needle: bytes) -> list[int]:
    out, index = [], buffer.find(needle)
    while index >= 0:
        out.append(index)
        index = buffer.find(needle, index + 1)
    return out


def _channel_names(archive: bytes) -> list[tuple[int, str]]:
    """``(offset_of_name, name)`` for every plate data set, in file order."""
    found = []
    for anchor in _find_all(archive, _CHANNEL_ANCHOR):
        start = anchor + len(_CHANNEL_ANCHOR)
        if start >= len(archive):
            continue
        length = archive[start]
        text = archive[start + 1:start + 1 + length]
        if length and len(text) == length and all(32 <= b < 127 for b in text):
            found.append((start, text.decode("ascii")))
    return found


def _read_start(archive: bytes, name_end: int) -> _dt.datetime | None:
    """The channel's read-start time, or ``None`` if the field does not decode as one."""
    offset = name_end + _TIMESTAMP_GAP
    if offset + 8 > len(archive):
        return None
    (stamp,) = struct.unpack_from("<q", archive, offset)
    if stamp <= 0:
        return None
    try:
        moment = _dt.datetime.fromtimestamp(stamp, _dt.UTC)
    except (OSError, OverflowError, ValueError):
        return None
    return moment if _EPOCH_WINDOW[0] <= moment <= _EPOCH_WINDOW[1] else None


def _read_well_block(archive: bytes, anchor: int) -> tuple[np.ndarray, int, int, int]:
    """One well block -> ``(values[n_times, n_wells], rows, cols, offset_after)``.

    Raises:
        Gen5FormatError: if the block runs past the stream, or if the per-record token is
            not constant across it -- which is what a wrong stride looks like.
    """
    fields = struct.unpack_from("<7I", archive, anchor + len(_WELL_BLOCK_ANCHOR))
    n_times, one_a, rows, cols, one_b, one_c, one_d = fields
    if (one_a, one_b, one_c, one_d) != (1, 1, 1, 1):
        raise Gen5FormatError(
            f"well block at {anchor} has unexpected header fields {fields}; this reader "
            "has only ever seen ones in those four positions and will not guess")
    n_wells = rows * cols
    start = anchor + _WELL_BLOCK_HEADER_SIZE
    span = n_times * n_wells * _RECORD_SIZE
    if n_times <= 0 or n_wells <= 0 or start + span > len(archive):
        raise Gen5FormatError(
            f"well block at {anchor} claims {n_times} x {rows} x {cols} records, which "
            f"does not fit in {len(archive)} bytes of archive")
    raw = np.frombuffer(archive, dtype=np.uint8, count=span, offset=start)
    grid = raw.reshape(n_times * n_wells, _RECORD_SIZE)
    tokens = np.unique(grid[:, 16:], axis=0)
    if len(tokens) != 1:
        raise Gen5FormatError(
            f"well block at {anchor}: the per-record token is not constant across "
            f"{n_times * n_wells} records ({len(tokens)} distinct), so the {_RECORD_SIZE}-byte "
            "stride does not hold for this file")
    values = grid[:, :8].copy().view("<f8").reshape(n_times, n_wells)
    return values, rows, cols, start + span


def _read_elapsed_table(archive: bytes, offset: int, n_times: int) -> np.ndarray | None:
    """Milliseconds since the run started, one slot per planned timepoint.

    ``None`` when the block has no elapsed table at all, which is how an endpoint read --
    one timepoint, no kinetic axis -- is written: length zero.

    The returned array is as long as the block's *planned* step count and may end in
    zeros; :func:`_completed_reads` is what decides how many of those slots were filled.

    Raises:
        Gen5FormatError: on any length other than zero or ``n_times``, or if the fields
            this reader treats as padding are not zero.
    """
    zero, one, count = struct.unpack_from("<HHI", archive, offset)
    if (zero, one) != (0, 1):
        raise Gen5FormatError(
            f"elapsed table at {offset} starts with {(zero, one)} rather than (0, 1)")
    if count == 0:
        return None
    if count != n_times:
        raise Gen5FormatError(
            f"elapsed table at {offset} holds {count} entries for a block of {n_times} "
            "timepoints")
    body = offset + 8
    span = count * _ELAPSED_ENTRY_SIZE
    if body + span > len(archive):
        raise Gen5FormatError(f"elapsed table at {offset} runs past the end of the archive")
    entries = np.frombuffer(archive, dtype="<u4", count=count * 4, offset=body).reshape(count, 4)
    if entries[:, 0].any() or entries[:, 2].any() or entries[:, 3].any():
        raise Gen5FormatError(
            f"elapsed table at {offset} uses fields this reader treats as padding; "
            "its milliseconds cannot be trusted")
    return entries[:, 1].astype(np.int64)


def _completed_reads(
    channel: str,
    values: np.ndarray,
    elapsed: np.ndarray | None,
    temperature: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    """Drop the slots a run that stopped early never filled.

    A well block is sized by the *protocol*, not by what happened: a 24-hour, 145-step
    protocol stopped after 123 reads still writes 145 rows, with the last 22 left as
    exact zeros in the readings, in the elapsed table, and in the temperature series.
    Reporting those rows would hand a caller 22 timepoints at OD 0.000. Six of the 22
    files here ended that way, and in every one the two counts agree exactly -- the number
    of non-zero elapsed entries equals the number of rows with any non-zero reading -- which
    is what makes the truncation a reading of the file rather than a repair of it. For
    ``20260708_ER_stress_1st.xpt`` it gives 123 reads, and ``data/plates`` holds that
    plate's export with 123 rows.

    Two files (``Experiment1``, ``RFP-read_20-min-interval_24h``, both 2026-06) instead
    write a *full-length* table of zeros over readings that are all present. A 49-point
    run did not happen at t=0 forty-nine times, so those get ``None`` for elapsed time and
    keep every row: the file carries no time axis, rather than the run having stopped
    before its first read.

    Returns:
        ``(values, elapsed, temperature)`` cut to the reads that happened.

    Raises:
        Gen5FormatError: if the elapsed times are not strictly increasing over the filled
            slots, or if the readings disagree with them about how many reads there were.
    """
    n_times = values.shape[0]
    if elapsed is None:
        return values, None, temperature
    stamped = np.flatnonzero(elapsed)
    if stamped.size == 0:
        return values, None, temperature
    n_reads = int(stamped[-1]) + 1
    if not np.all(np.diff(elapsed[:n_reads]) > 0):
        raise Gen5FormatError(
            f"channel {channel!r}: elapsed times are not strictly increasing over the "
            f"first {n_reads} slots; two reads cannot share a millisecond, so this field "
            "is not elapsed time in this file")
    if n_reads == n_times:
        return values, elapsed, temperature
    filled = np.flatnonzero(np.abs(values).sum(axis=1) > 0)
    n_filled = int(filled[-1]) + 1 if filled.size else 0
    if n_filled != n_reads:
        raise Gen5FormatError(
            f"channel {channel!r}: {n_reads} timepoints carry an elapsed time but "
            f"{n_filled} rows carry a reading. This reader will not choose between them")
    cut_temperature = None if temperature is None else temperature[:n_reads]
    return values[:n_reads], elapsed[:n_reads], cut_temperature


def _read_temperatures(archive: bytes, channel: str, n_times: int) -> np.ndarray | None:
    """The temperature series belonging to ``channel``, or ``None`` if the file has none.

    Gen5 writes a ``T<degree> <channel>`` data set per read step, and -- in the longer
    runs -- also an all-zero placeholder of the same name whose header ends in 0x00010001
    instead of 1. Requiring the header to be exactly ``(2, n_times, 1, 1)`` is what tells
    the measured series from the placeholder; two 2026-06 files contain only the
    placeholder and correctly yield ``None``.

    Raises:
        Gen5FormatError: if several candidate series validate and disagree.
    """
    base = _REPEAT_SUFFIX.sub("", channel)
    name = f"T{_DEGREE} {base}".encode("latin-1")
    key = bytes([len(name)]) + name
    wanted = struct.pack("<HIII", _TEMPERATURE_KIND, n_times, 1, 1)
    span = n_times * _TEMPERATURE_RECORD_SIZE
    series = []
    for hit in _find_all(archive, key):
        header = hit + len(key) + _TEMPERATURE_HEADER_GAP
        body = header + len(wanted)
        if body + span > len(archive):
            continue
        if archive[header:body] != wanted:
            continue
        grid = np.frombuffer(archive, dtype=np.uint8, count=span, offset=body)
        grid = grid.reshape(n_times, _TEMPERATURE_RECORD_SIZE)
        series.append(grid[:, :8].copy().view("<f8").reshape(n_times))
    if not series:
        return None
    for other in series[1:]:
        if not np.array_equal(series[0], other):
            raise Gen5FormatError(
                f"channel {channel!r} has several temperature series that disagree")
    return series[0]


def _protocol_name(contents: bytes) -> str | None:
    """The ``.prt`` protocol's file name, without the directory it sat in.

    The archive stores the full path on the instrument PC. Only the base name is
    surfaced: the directory names an account on someone's machine and says nothing about
    the measurement, and this repository keeps such paths out of tracked artefacts.
    """
    for end in (index + 4 for index in _find_all(contents, b".prt")):
        for length in range(5, 256):
            start = end - length
            if start < 1 or contents[start - 1] != length:
                continue
            text = contents[start:end]
            if all(32 <= b < 127 for b in text):
                return text.replace(b"\\", b"/").rsplit(b"/", 1)[-1].decode("ascii")
    return None


# ---------------------------------------------------------------------------
# what a caller gets
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Gen5Channel:
    """One optical channel's readings, as the instrument recorded them.

    Nothing here is blank-subtracted, pathlength-corrected, or rounded: these are the
    numbers Gen5 stored, which for absorbance carry one more decimal than its ``.xlsx``
    export does.
    """

    name: str
    """The archive's own label, e.g. ``'mCitrine:480,530'`` or ``'mCitrine:480,530[2]'``."""

    fluorophore: str
    """``'OD600'``, ``'mCitrine'``, ``'RFP'`` -- the part before the colon."""

    optics: str
    """``'600'`` for absorbance, ``'480,530'`` excitation,emission for fluorescence."""

    repeat: int
    """1, or *n* for the ``[n]`` suffix Gen5 gives a second read of one fluorophore."""

    wells: tuple[str, ...]
    """Column labels for :attr:`values`, row-major from A1.

    Every well of the plate, not only the ones an export happened to include.
    """

    values: np.ndarray
    """``(n_times, n_wells)`` float64, in the units of :attr:`fluorophore` (OD or RFU)."""

    elapsed_ms: np.ndarray | None
    """Milliseconds from the run's start to each read, or ``None`` for an endpoint read."""

    temperature_c: np.ndarray | None
    """Incubator temperature at each read, or ``None`` when the file logged none."""

    started_at: _dt.datetime | None
    """When this channel's first read began, UTC, or ``None`` if the stamp did not decode."""

    rows: int
    cols: int

    @property
    def n_times(self) -> int:
        return int(self.values.shape[0])

    def elapsed_hours(self) -> np.ndarray:
        """Elapsed time in hours.

        Raises:
            Gen5FormatError: for an endpoint read, which has no elapsed axis. Callers that
                should skip must test :attr:`elapsed_ms` for ``None`` instead.
        """
        if self.elapsed_ms is None:
            raise Gen5FormatError(f"channel {self.name!r} is an endpoint read: no elapsed axis")
        return self.elapsed_ms / 3_600_000.0

    def elapsed_hms(self) -> tuple[str, ...]:
        """Elapsed time as ``'0:08:18'`` -- the spelling ``data/plates`` uses.

        Raises:
            Gen5FormatError: for an endpoint read, as :meth:`elapsed_hours`.
        """
        if self.elapsed_ms is None:
            raise Gen5FormatError(f"channel {self.name!r} is an endpoint read: no elapsed axis")
        out = []
        for milliseconds in self.elapsed_ms.tolist():
            seconds = milliseconds // 1000
            out.append(f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}")
        return tuple(out)

    def overflowed(self) -> np.ndarray:
        """Which readings are :data:`OVERFLOW` rather than measurements, as a boolean mask."""
        return np.asarray(self.values, dtype=float) <= OVERFLOW + 1.0

    def frame(self, on_overflow: str = "raise") -> pd.DataFrame:
        """Readings as a frame: one row per timepoint, one column per well.

        Indexed by elapsed hours when there is an elapsed axis, and by read number when
        there is not.

        Args:
            on_overflow: What to do when the channel contains :data:`OVERFLOW`.
                ``'raise'`` (the default) refuses, ``'nan'`` substitutes ``NaN``, and
                ``'keep'`` returns the sentinel as the float the file holds.

        Raises:
            Gen5FormatError: when the channel contains overflow readings and
                ``on_overflow='raise'``.
            ValueError: for an unknown ``on_overflow``.

        The default refuses because the failure it prevents is silent. A sentinel that
        reaches an average is not a wrong reading, it is a number with no relationship to
        the plate, and nothing downstream can tell it from a measurement. No channel any
        caller in this repository reads contains one -- the guard is for the next caller.
        """
        if on_overflow not in ("raise", "nan", "keep"):
            raise ValueError(
                f"on_overflow must be 'raise', 'nan' or 'keep', got {on_overflow!r}")
        values = np.asarray(self.values, dtype=float)
        flagged = self.overflowed()
        if flagged.any():
            if on_overflow == "raise":
                raise Gen5FormatError(
                    f"channel {self.name!r} has {int(flagged.sum())} of {flagged.size} "
                    f"readings at the overflow flag {OVERFLOW:.0f}, which are not "
                    "measurements. Pass on_overflow='nan' to mask them or 'keep' to take "
                    "them as they are; do not average them")
            if on_overflow == "nan":
                values = np.where(flagged, np.nan, values)
        if self.elapsed_ms is None:
            index = pd.Index(range(self.n_times), name="read")
        else:
            index = pd.Index(self.elapsed_hours(), name="elapsed_h")
        return pd.DataFrame(values, index=index, columns=list(self.wells))


@dataclass(frozen=True)
class Gen5Read:
    """Everything one ``.xpt`` holds that this reader can verify."""

    source: pathlib.Path
    protocol: str | None
    """The ``.prt`` file name the run was made with, e.g. ``'4h-10min_mCitrine_OD600.prt'``."""

    channels: tuple[Gen5Channel, ...]

    def channel(self, name: str) -> Gen5Channel:
        """One channel by its archive name.

        Raises:
            KeyError: naming the channels the file does have.
        """
        for channel in self.channels:
            if channel.name == name:
                return channel
        raise KeyError(f"{name!r} not in {[c.name for c in self.channels]}")


def read_xpt(path: str | pathlib.Path) -> Gen5Read:
    """Read a Gen5 ``.xpt`` experiment file.

    Args:
        path: The ``.xpt``. Read whole into memory; these are 0.1-2.3 MB.

    Returns:
        Every channel the archive holds, in file order, with all 96 wells present whether
        or not the operator's ``.xlsx`` export included them.

    Raises:
        Gen5FormatError: if the container, the channel list, or any block fails to
            validate. This reader refuses rather than returning numbers it cannot place.
    """
    path = pathlib.Path(path)
    document = _CompoundFile(path.read_bytes())
    archive = document.read("SUBSETS/1/DATA")
    try:
        contents = document.read("Contents")
    except Gen5FormatError:
        contents = b""

    names = _channel_names(archive)
    anchors = _find_all(archive, _WELL_BLOCK_ANCHOR)
    if not anchors:
        raise Gen5FormatError(f"{path.name}: no well-readings block found")
    if len(names) < len(anchors):
        raise Gen5FormatError(
            f"{path.name}: {len(anchors)} well blocks but only {len(names)} channel names")

    channels = []
    for anchor in anchors:
        preceding = [(offset, name) for offset, name in names if offset < anchor]
        if not preceding:
            raise Gen5FormatError(
                f"{path.name}: well block at {anchor} has no channel name before it")
        name_offset, name = preceding[-1]
        values, rows, cols, after = _read_well_block(archive, anchor)
        elapsed = _read_elapsed_table(archive, after, values.shape[0])
        temperature = _read_temperatures(archive, name, values.shape[0])
        values, elapsed, temperature = _completed_reads(name, values, elapsed, temperature)
        fluorophore, _, optics = name.partition(":")
        match = _REPEAT_SUFFIX.search(optics)
        channels.append(Gen5Channel(
            name=name,
            fluorophore=fluorophore,
            optics=_REPEAT_SUFFIX.sub("", optics),
            repeat=int(match.group(1)) if match else 1,
            wells=well_names(rows, cols),
            values=values,
            elapsed_ms=elapsed,
            temperature_c=temperature,
            started_at=_read_start(archive, name_offset + 1 + len(name)),
            rows=rows,
            cols=cols,
        ))
    return Gen5Read(source=path, protocol=_protocol_name(contents), channels=tuple(channels))


# ---------------------------------------------------------------------------
# the instrument file in the shape the rest of the pipeline speaks
# ---------------------------------------------------------------------------


def as_synergy_run(read: Gen5Read) -> "SynergyRun":
    """A :class:`Gen5Read` as the :class:`~ystwin.plate.synergy.SynergyRun` every script wants.

    Everything downstream of a plate -- ``prepare()``, ``aligned()``, the layout recovery,
    ``plate_readings.py --export`` -- takes a ``SynergyRun``. Until this existed, the only
    way to build one was to read an ``.xlsx``, so a plate whose every export is
    blank-subtracted could be *read* by this module and still not be *used* by anything.
    ``20260804`` was exactly that plate.

    Nothing is converted, rescaled or repaired here. This is a change of container, with
    one exception that is a change of *expression* and not of value.

    **The time axis is decoded from ``H:MM:SS``, not divided out of the milliseconds.**
    ``elapsed_ms / 3_600_000`` and ``h + m/60 + s/3600`` are the same real number and
    different floats: on ``20260804`` they disagree by one ulp at 16 of the 50 timestamps.
    ``data/plates`` stores elapsed time as ``H:MM:SS`` and ``plate_readings.py --export``
    refuses to leave an export in place unless the rebuilt index is ``==`` the reader's, so
    a run built on the millisecond expression could never be committed. The whole-second
    check above is what makes the spelling lossless rather than merely conventional.

    Two fields have no counterpart in a ``.xpt`` and are filled from what the archive does
    hold rather than left blank or invented:

    ``sheet``
        A ``.xpt`` has no sheets. The locator that plays the same role -- "which part of
        the file this block came from" -- is the archive's own data-set name, so that is
        what goes in: ``'OD600:600'``, ``'mCitrine:480,530[2]'``. It is not a sheet name
        and does not pretend to be one; it is the name the instrument gave the read.

    ``blank_subtracted``
        ``False``: this reader decodes the archive's untransformed well-readings data sets,
        not Gen5's export transforms. The paired-source checks in ``tests/test_gen5_xpt.py``
        verify that format contract. Finite negative readings remain readings and do not
        change their processing state; only :data:`OVERFLOW` is masked as missing.

    **The overflow flag becomes NaN, and that is a change of expression.** It is the one
    place this function departs from "nothing is converted", and it is not a repair: ``NaN``
    is what "the instrument reported no value here" means to every consumer downstream,
    where -99999 means it to none of them. ``aligned()`` averaged the sentinel into a
    per-well mean of -62,622 RFU before this. The count is not hidden -- :meth:`
    Gen5Channel.overflowed` still reports it on the channel this was built from, and
    :meth:`Gen5Channel.frame` still refuses by default.

    Args:
        read: One parsed instrument file.

    Returns:
        A run whose blocks are in archive order, carrying every well the file holds.

    Raises:
        Gen5FormatError: if any channel is an endpoint read. A ``SynergyRun`` is indexed by
            elapsed hours throughout, and a read with no time axis has nothing to index by;
            substituting the read number would put a plate on a fabricated clock.
    """
    from .synergy import KineticBlock, SynergyRun, _to_hours, label_repeat_reads

    blocks = []
    for channel in read.channels:
        if channel.elapsed_ms is None:
            raise Gen5FormatError(
                f"{read.source.name}: channel {channel.name!r} is an endpoint read with no "
                "elapsed axis, and a SynergyRun is indexed by elapsed hours")
        stray = [int(ms) for ms in channel.elapsed_ms.tolist() if int(ms) % 1000]
        if stray:
            raise Gen5FormatError(
                f"{read.source.name}: channel {channel.name!r} has {len(stray)} elapsed "
                f"stamps that are not a whole number of seconds (first {stray[0]} ms). "
                "H:MM:SS cannot carry them, so this run cannot be committed as text")
        flagged = channel.overflowed()
        values = np.asarray(channel.values, dtype=float)
        readings = np.where(flagged, np.nan, values)
        frame = pd.DataFrame(
            readings,
            index=pd.Index([_to_hours(t) for t in channel.elapsed_hms()], name="time_h"),
            columns=list(channel.wells),
            dtype=float,
        )
        blocks.append(KineticBlock(
            channel=channel.fluorophore,  # provisional, as in synergy.py
            fluorophore=channel.fluorophore,
            optics=channel.optics,
            sheet=channel.name,
            derived=False,  # the archive holds plate reads; the plot sheets are made on export
            blank_subtracted=False,
            data=frame,
            temperature_c=([] if channel.temperature_c is None
                           else [float(t) for t in channel.temperature_c]),
        ))
    return SynergyRun(source=read.source, blocks=tuple(label_repeat_reads(blocks)))


def read_xpt_run(path: str | pathlib.Path) -> "SynergyRun":
    """:func:`read_xpt` then :func:`as_synergy_run`, with the signature every reader has.

    ``scripts/plate_readings.py`` swaps readers by name; this is the ``.xpt`` one, and it
    matches ``plate/synergy.py::read_synergy_kinetic`` argument for argument so the export
    and the verification can call the same thing.
    """
    return as_synergy_run(read_xpt(path))
