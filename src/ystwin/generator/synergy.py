"""Read a Synergy H1 kinetic export: plate layout plus one table per read.

The reader writes each read as its own block -- a name, a Time header carrying the well
IDs, then one row per timepoint -- and prefixes blank-subtracted copies with "Blank". The
layout block maps wells to sample names, which is the only place the plate design is
recorded.
"""
import re
import pandas as pd

__all__ = ["read_export"]

WELL = re.compile(r"^[A-H](?:[1-9]|1[0-2])$")
CLOCK = re.compile(r"^\d+:\d\d:\d\d$")


def _clock(v):
    if isinstance(v, str) and CLOCK.match(v.strip()):
        h, m, s = (int(x) for x in v.split(":"))
        return h + m / 60 + s / 3600
    if hasattr(v, "hour"):
        return v.hour + v.minute / 60 + v.second / 3600
    return None


def read_export(path, sheet=0):
    with pd.ExcelFile(path) as book:
        df = book.parse(sheet, header=None)
    layout, tables, name = {}, {}, None
    for i in range(len(df)):
        a = str(df.iat[i, 0]).strip() if pd.notna(df.iat[i, 0]) else ""
        b = str(df.iat[i, 1]).strip() if pd.notna(df.iat[i, 1]) else ""
        if a == "Layout":
            for r in range(i + 1, i + 12):
                row = str(df.iat[r, 1]).strip() if pd.notna(df.iat[r, 1]) else ""
                if row not in list("ABCDEFGH"):
                    continue
                for j in range(12):
                    v = df.iat[r, 2 + j]
                    if pd.notna(v) and str(v).strip():
                        layout[f"{row}{j + 1}"] = str(v).strip()
        elif a and not b and ":" in a:
            name = a
        elif b == "Time" and name:
            cols = {int(c): str(df.iat[i, c]).strip() for c in range(2, df.shape[1])
                    if pd.notna(df.iat[i, c]) and WELL.match(str(df.iat[i, c]).strip())}
            times, vals = [], []
            for r in range(i + 1, len(df)):
                t = _clock(df.iat[r, 1])
                if t is None:
                    break
                times.append(t)
                vals.append([pd.to_numeric(df.iat[r, c], errors="coerce") for c in cols])
            tables[name] = pd.DataFrame(vals, columns=list(cols.values()),
                                        index=pd.Index(times, name="hour")).astype(float)
            name = None
    return layout, tables
