"""Summarise one or more Gen5 ``.xpt`` experiment files -- the reader's own format.

The ``.xlsx`` a Gen5 operator exports carries whatever transform was on screen at the
time, rounded to what was on screen. The ``.xpt`` carries the readings. This prints what
is in one, so that "is the raw background still recoverable from this plate?" can be
answered by looking rather than by remembering.

``--blanks`` is the reason the script exists. Give it the wells that held medium only and
it prints their readings next to the plate's occupied wells, so a subtracted export shows
up immediately: real blanks read ~0.09 OD and a few hundred RFU, whereas an export that
has already had its background removed reads within 0.002 of zero.

Usage:
    python scripts/read_gen5_xpt.py FILE.xpt [FILE.xpt ...]
    python scripts/read_gen5_xpt.py FILE.xpt --blanks H1,H2,H3
    python scripts/read_gen5_xpt.py FILE.xpt --csv DIR     # one CSV per channel

``--csv`` writes the same column layout as ``data/plates`` -- ``elapsed_hms``,
``temperature_c``, then one column per well -- but with every well of the plate and the
absorbance decimal the ``.xlsx`` export drops. It writes only where it is told to; it
never touches ``data/plates`` on its own.
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin.plate.gen5 import Gen5Channel, Gen5FormatError, Gen5Read, read_xpt


def _rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def _describe(channel: Gen5Channel) -> str:
    if channel.elapsed_ms is None:
        span = "endpoint (no elapsed axis)"
    else:
        hours = channel.elapsed_hours()
        span = f"{channel.elapsed_hms()[0]} to {channel.elapsed_hms()[-1]} ({hours[-1]:.2f} h)"
    if channel.temperature_c is None:
        temperature = "no temperature logged"
    else:
        temperature = (f"{channel.temperature_c.min():.1f}-{channel.temperature_c.max():.1f} "
                       f"{chr(0xB0)}C")
    return (f"  {channel.name:22s} {channel.n_times:4d} reads x {len(channel.wells)} wells | "
            f"{span} | {temperature}\n"
            f"  {'':22s} values {channel.values.min():.4g} to {channel.values.max():.4g}"
            f" | first read {channel.started_at or 'unstamped'}")


def _report_blanks(read: Gen5Read, blanks: list[str]) -> None:
    """Print the named wells against the rest of the plate, per channel."""
    _rule(f"blank wells {', '.join(blanks)}")
    for channel in read.channels:
        index = {well: i for i, well in enumerate(channel.wells)}
        missing = [w for w in blanks if w not in index]
        if missing:
            print(f"  {channel.name}: no such wells on this plate: {', '.join(missing)}")
            continue
        columns = [index[w] for w in blanks]
        blank = channel.values[:, columns]
        occupied = np.delete(channel.values, columns, axis=1)
        print(f"  {channel.name}")
        print(f"    blank wells      first read {np.array2string(blank[0], precision=4)}"
              f"  last read {np.array2string(blank[-1], precision=4)}")
        print(f"    blank mean       {blank.mean():.4f}   range over the run "
              f"{blank.min():.4f} to {blank.max():.4f}")
        print(f"    rest of plate    mean {occupied.mean():.4f}, "
              f"min {occupied.min():.4f}, max {occupied.max():.4f}")
        ratio = abs(occupied.mean()) / abs(blank.mean()) if blank.mean() else float("inf")
        floor = 1e-3 * max(abs(occupied.mean()), 1e-9)
        verdict = ("background is present" if abs(blank.mean()) > floor else
                   "blank reads ~0: this channel looks already background-subtracted")
        print(f"    plate / blank    {ratio:.1f}x -- {verdict}")


def _write_csv(read: Gen5Read, out_dir: pathlib.Path) -> list[pathlib.Path]:
    """One CSV per channel, in the ``data/plates`` column layout."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for channel in read.channels:
        safe = channel.name.replace(":", "_").replace(",", "-")
        path = out_dir / f"{read.source.stem}__{safe}.csv"
        header = ["elapsed_hms"]
        if channel.temperature_c is not None:
            header.append("temperature_c")
        header.extend(channel.wells)
        times = channel.elapsed_hms() if channel.elapsed_ms is not None else \
            tuple(str(i) for i in range(channel.n_times))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            for row in range(channel.n_times):
                line = [times[row]]
                if channel.temperature_c is not None:
                    line.append(f"{channel.temperature_c[row]:.1f}")
                line.extend(f"{v:g}" for v in channel.values[row])
                writer.writerow(line)
        written.append(path)
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="+", type=pathlib.Path, help="Gen5 .xpt experiment files")
    parser.add_argument("--blanks", default="",
                        help="comma-separated wells that held medium only, e.g. H1,H2,H3")
    parser.add_argument("--csv", type=pathlib.Path, default=None,
                        help="directory to write one CSV per channel into")
    args = parser.parse_args(argv)

    blanks = [w.strip().upper() for w in args.blanks.split(",") if w.strip()]
    failures = 0
    for path in args.files:
        if not path.exists():
            print(f"{path}: not found", file=sys.stderr)
            failures += 1
            continue
        try:
            read = read_xpt(path)
        except Gen5FormatError as exc:
            print(f"{path.name}: {exc}", file=sys.stderr)
            failures += 1
            continue
        _rule(path.name)
        print(f"  protocol {read.protocol or 'not recorded in this file'}")
        for channel in read.channels:
            print(_describe(channel))
        if blanks:
            _report_blanks(read, blanks)
        if args.csv is not None:
            for written in _write_csv(read, args.csv):
                print(f"  wrote {written}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
