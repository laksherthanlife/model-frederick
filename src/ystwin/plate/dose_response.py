"""Reader for the per-construct dose-response sheets.

One sheet per construct, normalised fluorescence against time, one column per dose. This is the
reporter response as the team reads it, which is what G4 must be tested against.
"""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd

from ..qpcr import TARGET_FOR_CONSTRUCT, _normalise_construct, _parse_dose

__all__ = ["endpoint_response", "read_dose_response"]

_TIME_HEADER = re.compile(r"^\s*time\s*(\(.*\))?\s*$", re.IGNORECASE)


def read_dose_response(path: str | pathlib.Path) -> pd.DataFrame:
    """Read every per-construct dose-response sheet in a workbook.

    Args:
        path: Workbook to read. Sheets whose name is not a known construct are
            skipped, so the OD and raw fluorescence sheets pass through untouched.

    Returns:
        Columns ``construct``, ``dose_mM``, ``time_h``, ``signal``.
    """
    path = pathlib.Path(path)
    records: list[dict] = []
    for sheet, raw in pd.read_excel(path, sheet_name=None, header=None).items():
        construct = _normalise_construct(sheet)
        if construct is None or raw.empty:
            continue
        records.extend(_extract_dose_block(raw, construct))
    if not records:
        raise ValueError(
            f"no per-construct dose-response sheet found in {path.name}; expected a "
            f"sheet named for one of {sorted(TARGET_FOR_CONSTRUCT)}"
        )
    return pd.DataFrame(records)


def _extract_dose_block(raw: pd.DataFrame, construct: str) -> list[dict]:
    n_rows, n_cols = raw.shape
    for r in range(n_rows):
        time_col = next(
            (c for c in range(n_cols)
             if isinstance(raw.iat[r, c], str) and _TIME_HEADER.match(raw.iat[r, c])),
            None,
        )
        if time_col is None:
            continue
        dose_cols = {}
        for c in range(time_col + 1, n_cols):
            dose = _parse_dose(raw.iat[r, c])
            if dose is not None:
                dose_cols[c] = dose
        if not dose_cols:
            continue

        records = []
        for rr in range(r + 1, n_rows):
            minutes = pd.to_numeric(raw.iat[rr, time_col], errors="coerce")
            values = {c: pd.to_numeric(raw.iat[rr, c], errors="coerce") for c in dose_cols}
            if all(pd.isna(v) for v in values.values()):
                break
            if pd.isna(minutes):
                continue
            for c, dose in dose_cols.items():
                if pd.isna(values[c]):
                    continue
                records.append({
                    "construct": construct,
                    "dose_mM": float(dose),
                    "time_h": float(minutes) / 60.0,
                    "signal": float(values[c]),
                })
        if records:
            return records
    return []


def endpoint_response(
    tidy: pd.DataFrame,
    fraction: float = 0.25,
    relative_to_dose: float | None = None,
) -> pd.DataFrame:
    """Summarise each construct and dose by its late-timecourse signal.

    Args:
        tidy: Frame from :func:`read_dose_response`.
        fraction: Share of the timecourse, taken from the end, to average over.
            Averaging beats a single last point, which is the noisiest reading.
        relative_to_dose: If given, express every dose relative to this one, so the
            output is a fold change comparable with a qPCR anchor.
    """
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    rows = []
    for (construct, dose), group in tidy.groupby(["construct", "dose_mM"]):
        ordered = group.sort_values("time_h")
        n = max(1, int(round(len(ordered) * fraction)))
        rows.append({
            "construct": construct,
            "dose_mM": dose,
            "response": float(ordered.signal.tail(n).mean()),
        })
    out = pd.DataFrame(rows)
    if relative_to_dose is not None:
        baseline = out[np.isclose(out.dose_mM, relative_to_dose)][["construct", "response"]]
        if baseline.empty:
            raise ValueError(f"dose {relative_to_dose} is not present to normalise against")
        out = out.merge(baseline.rename(columns={"response": "_base"}), on="construct")
        out["response"] = out.response / out._base
        out = out.drop(columns="_base")
    return out.sort_values(["construct", "dose_mM"]).reset_index(drop=True)
