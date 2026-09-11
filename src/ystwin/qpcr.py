"""Relative expression from Cq -- the independent anchor for G4.

Hac1 for the ER reporters, TRX2 for the oxidative ones, both against UBC. A ratio of two
transcripts in the same sample carries no dilution term, which is what makes it independent of
the growth confound. Fold change follows Livak, generalised to efficiency below one.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = [
    "STRESSOR_FOR_CONSTRUCT",
    "TARGET_FOR_CONSTRUCT",
    "assert_single_stressor",
    "REFERENCE_TARGET",
    "QPCR_PLATE",
    "QpcrPlate",
    "delta_delta_cq",
    "gdna_share",
    "read_annotated_cq",
    "read_positional_cq",
    "rt_minus_margin",
]

REFERENCE_TARGET = "UBC"
"""Label the qPCR exports use for the normalising transcript.

The string is the export's, not a choice made here, so it stays as written -- but it is
worth flagging that **"UBC" is not a gene name in this organism.** *S. cerevisiae* has
UBC1 through UBC13 plus UBI4, and which one was actually amplified is fixed by the
ordered oligos rather than by anything recoverable from the data. The identity matters
in two ways: UBC6 is annotated to ER-associated protein catabolism, which is the very
process the ER reporters are challenged to induce and therefore a poor choice of
invariant reference under DTT; and UBC4 is one of the few yeast genes carrying an
intron, which changes what a no-RT control means for it.

Two consequences worth acting on rather than noting. Read the identity off the primer
order and record it here. And a single reference transcript is below current practice
regardless -- the standard is the geometric mean of two or three validated genes,
because a lone normaliser that itself responds to the treatment divides the signal out.
See docs/research/UPR_ANCHOR.md for a recommended set.
"""

# Which endogenous transcript reports on the same stress each construct senses.
TARGET_FOR_CONSTRUCT = {
    "UPRE1": "Hac1",
    "UPRE2": "Hac1",
    "NativeYap1": "TRX2",
    "AlteredYap1": "TRX2",
}

# The two sensor pairs are challenged with different agents, so a dose on one axis
# says nothing about the other.
STRESSOR_FOR_CONSTRUCT = {
    "UPRE1": "DTT",
    "UPRE2": "DTT",
    "NativeYap1": "H2O2",
    "AlteredYap1": "H2O2",
}


def assert_single_stressor(constructs) -> str:
    """Return the stressor shared by ``constructs``, or refuse to pool them.

    Raises:
        ValueError: if the constructs were challenged with different agents, or if
            any of them has no declared stressor.
    """
    unknown = [c for c in constructs if c not in STRESSOR_FOR_CONSTRUCT]
    if unknown:
        raise ValueError(f"no stressor declared for {sorted(unknown)}")
    agents = {STRESSOR_FOR_CONSTRUCT[c] for c in constructs}
    if len(agents) > 1:
        detail = ", ".join(f"{c} ({STRESSOR_FOR_CONSTRUCT[c]})" for c in constructs)
        raise ValueError(
            f"these constructs were challenged with different stressors -- {detail}; "
            "their dose axes are not comparable and must not be pooled"
        )
    return agents.pop()

_MIN_RT_MINUS_MARGIN = 3.0
"""Smallest RT- separation, in cycles, treated as a quantitative measurement.

Three cycles is 12.5% genomic DNA. The widely quoted alternative is five cycles, which
is 3.1% -- but that number is a software default, not a guideline requirement, and it
began life as a *no-template-control* rule. MIQE 2009, its 2010 precis, and MIQE 2.0
(2025) all require only that the +RT/-RT comparison be reported and that the
experimenter state their own tolerance.
"""

_MAX_CORRECTABLE_GDNA_SHARE = 0.60
"""Above this contamination fraction, no published method can correct the reading.

ValidPrime (Laurell et al., Nucleic Acids Research 2012;40(7):e51) is the only
validated genomic-DNA correction for RT-qPCR, and it states three times that
correction holds only "as long as the DNA contribution to the total signal is <60%",
above which its own software refuses. Sixty percent is a margin of about 0.74 cycles.

The distinction this draws is worth having in the data rather than in a footnote: a
reading that merely fails the QC threshold might still be recoverable, while one past
this ceiling is not a measurement of transcript at all and no analysis will make it
one. See docs/research/G4_STATISTICS.md for the error propagation behind the number.
"""


def gdna_share(margin_cycles: np.ndarray, efficiency: float = 1.0) -> np.ndarray:
    """Fraction of the +RT signal that is genomic DNA, from the RT- margin.

    The margin is reported in cycles because that is what the instrument prints, but
    cycles are a logarithmic proxy for the quantity that decides whether a reading
    means anything: how much of it is contamination. Since amplification is
    exponential, a margin converts exactly --

        share = (1 + efficiency) ** (-margin)

    -- and the conversion is worth doing explicitly, because the cycle scale badly
    misleads intuition near zero. Three cycles sounds much worse than five and is
    12.5% against 3.1%; 0.6 cycles sounds only a little worse than 1.0 and is 66%
    against 50%, which straddles the limit of what can be corrected at all.

    Args:
        margin_cycles: ``Cq(-RT) - Cq(+RT)``, per reading. NaN propagates.
        efficiency: Amplification efficiency as a fraction, so 1.0 means perfect
            doubling per cycle. Measure it from a dilution series rather than
            assuming it: at efficiencies in the plausible 0.93-1.05 range a single
            reported margin admits a wide range of true shares.

    Returns:
        Contamination fraction in [0, 1]; NaN where the margin is unknown.
    """
    base = 1.0 + float(efficiency)
    if base <= 1.0:
        raise ValueError("efficiency must be positive; 1.0 is perfect doubling")
    return np.asarray(base, dtype=float) ** (-np.asarray(margin_cycles, dtype=float))


def delta_delta_cq(
    tidy: pd.DataFrame,
    reference_target: str = REFERENCE_TARGET,
    control_dose: float = 0.0,
    efficiency: float = 1.0,
) -> pd.DataFrame:
    """Relative expression per construct and dose, against a reference gene and dose.

    Args:
        tidy: Long frame with ``construct``, ``target``, ``dose_mM``, ``cq``.
            Technical replicates are averaged.
        reference_target: Housekeeping transcript used to normalise input amount.
        control_dose: Dose treated as the untreated baseline.
        efficiency: Amplification efficiency in (0, 1]. One means each cycle
            doubles; lower values mean a cycle of Cq buys less than a doubling, so
            the same ddCq corresponds to a smaller true fold change.

    Returns:
        One row per construct and dose with ``delta_cq``, ``delta_delta_cq`` and
        ``fold_change``.
    """
    if not 0 < efficiency <= 1:
        raise ValueError("efficiency must be in (0, 1]")
    required = {"construct", "target", "dose_mM", "cq"}
    missing = required - set(tidy.columns)
    if missing:
        raise ValueError(f"tidy frame is missing columns {sorted(missing)}")

    mean_cq = (
        tidy.groupby(["construct", "target", "dose_mM"], as_index=False)["cq"].mean()
    )
    reference = mean_cq[mean_cq.target == reference_target]
    targets = mean_cq[mean_cq.target != reference_target]
    if reference.empty:
        raise ValueError(
            f"no {reference_target!r} measurement in this run; relative expression "
            "cannot be computed without the reference gene"
        )

    merged = targets.merge(
        reference[["construct", "dose_mM", "cq"]].rename(columns={"cq": "cq_reference"}),
        on=["construct", "dose_mM"], how="left",
    )
    orphan = merged[merged.cq_reference.isna()]
    if not orphan.empty:
        pairs = sorted({(r.construct, r.dose_mM) for r in orphan.itertuples()})
        raise ValueError(
            f"no {reference_target!r} measurement for {pairs}; each sample must "
            "carry its own reference reading"
        )
    merged["delta_cq"] = merged.cq - merged.cq_reference

    baseline = merged[np.isclose(merged.dose_mM, control_dose)]
    if baseline.empty:
        raise ValueError(
            f"control dose {control_dose} is absent from this run; there is nothing "
            "to express the treated samples relative to"
        )
    merged = merged.merge(
        baseline[["construct", "target", "delta_cq"]].rename(
            columns={"delta_cq": "delta_cq_control"}
        ),
        on=["construct", "target"], how="left",
    )
    merged["delta_delta_cq"] = merged.delta_cq - merged.delta_cq_control
    merged["fold_change"] = (1.0 + efficiency) ** (-merged.delta_delta_cq)
    return merged.drop(columns=["delta_cq_control"]).sort_values(
        ["construct", "dose_mM"]
    ).reset_index(drop=True)


def rt_minus_margin(
    tidy: pd.DataFrame,
    min_cycles: float = _MIN_RT_MINUS_MARGIN,
    efficiency: float = 1.0,
) -> pd.DataFrame:
    """Distance between the no-RT control and the sample, per reading.

    Reverse transcriptase is omitted from the control, so anything that still
    amplifies is genomic DNA rather than transcript. A margin of fewer than a few
    cycles means part of the measured "expression" is contamination.

    Three columns come back rather than one, because "failed QC" collapses two
    situations that call for different actions. A reading below ``min_cycles`` is
    contaminated and should not be quoted as expression, but may still be
    recoverable. A reading whose contamination exceeds
    :data:`_MAX_CORRECTABLE_GDNA_SHARE` is past the limit of any published
    correction: it is not a measurement of transcript, and reporting it as one that
    merely failed a threshold overstates what is there. ``correctable`` marks the
    difference so a caller can say which it has.

    Args:
        tidy: Long frame carrying ``cq`` and ``rt_minus_cq``.
        min_cycles: Smallest acceptable separation. A missing control fails, since
            an unrun control is not a passed one.
        efficiency: Amplification efficiency for the cycles-to-share conversion.

    Returns:
        The input frame with ``margin_cycles``, ``gdna_share``, ``passed`` and
        ``correctable``. ``correctable`` is False where the control is missing --
        an unknown contamination level is not a correctable one.
    """
    out = tidy.copy()
    out["margin_cycles"] = out.rt_minus_cq - out.cq
    out["gdna_share"] = gdna_share(out.margin_cycles.to_numpy(), efficiency)
    out["passed"] = out.margin_cycles.notna() & (out.margin_cycles >= min_cycles)
    out["correctable"] = out.gdna_share.notna() & (
        out.gdna_share <= _MAX_CORRECTABLE_GDNA_SHARE
    )
    return out


# --- reader for the hand-annotated Cq sheets -------------------------------

import pathlib  # noqa: E402
import re  # noqa: E402

_CONSTRUCT_ALIASES = {
    "upre1": "UPRE1",
    "upre2": "UPRE2",
    "nativeyap1": "NativeYap1",
    "alteredyap1": "AlteredYap1",
}
_DOSE_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+)\s*mm\s*$", re.IGNORECASE)


def _normalise_construct(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    return _CONSTRUCT_ALIASES.get(re.sub(r"[^a-z0-9]+", "", value.lower()))


def _parse_dose(value: object) -> float | None:
    # NaN is a float, so an empty cell would otherwise read as a dose.
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return None if np.isnan(value) else float(value)
    if isinstance(value, str):
        match = _DOSE_RE.match(value)
        if match:
            return float(match.group(1))
    return None


def _resolve_target(label: object, construct: str) -> str | None:
    """``'Hac1 / TRX2'`` names two genes; which one applies depends on the block."""
    if not isinstance(label, str):
        return None
    options = [p.strip() for p in label.split("/") if p.strip()]
    if len(options) == 1:
        return options[0]
    expected = TARGET_FOR_CONSTRUCT.get(construct)
    for option in options:
        if expected and option.lower() == expected.lower():
            return option
    return expected


def read_annotated_cq(path: str | pathlib.Path) -> pd.DataFrame:
    """Read a hand-annotated Cq sheet into the long form ``delta_delta_cq`` expects.

    The construct header row is located by content, because the preamble depth
    differs between replicates.

    Args:
        path: Workbook to read. The first sheet holding a construct header is used.

    Returns:
        Columns ``construct``, ``target``, ``dose_mM``, ``tech_rep``, ``cq``,
        ``rt_minus_cq`` -- one row per technical replicate.
    """
    path = pathlib.Path(path)
    for raw in pd.read_excel(path, sheet_name=None, header=None).values():
        if raw.empty:
            continue
        records = _extract_annotated(raw)
        if records:
            return pd.DataFrame(records)
    raise ValueError(
        f"no construct header row found in {path.name}; this reader expects the "
        "annotated analysis layout, not the raw instrument export"
    )


def _extract_annotated(raw: pd.DataFrame) -> list[dict]:
    n_rows, n_cols = raw.shape
    header_row, blocks = None, {}
    for r in range(n_rows):
        found = {}
        for c in range(n_cols):
            name = _normalise_construct(raw.iat[r, c])
            if name:
                found[name] = c
        if len(found) >= 2:
            header_row, blocks = r, found
            break
    if header_row is None:
        return []

    records: list[dict] = []
    current_target: object = None
    for r in range(header_row + 1, n_rows):
        dose = None
        for c in range(min(blocks.values())):
            dose = _parse_dose(raw.iat[r, c])
            if dose is not None:
                dose_col = c
                break
        if dose is None:
            continue
        # Target is written once per group; a gene name is never one character.
        for c in range(dose_col + 1, min(blocks.values())):
            cell = raw.iat[r, c]
            if not isinstance(cell, str):
                continue
            text = cell.strip()
            if len(text) < 2 or _parse_dose(text) is not None:
                continue
            if _normalise_construct(text):
                continue
            current_target = text
            break
        if current_target is None:
            continue

        for construct, col in blocks.items():
            target = _resolve_target(current_target, construct)
            if target is None:
                continue
            rt_minus = pd.to_numeric(raw.iat[r, col + 3], errors="coerce") \
                if col + 3 < n_cols else np.nan
            for offset in (0, 1):
                if col + offset >= n_cols:
                    continue
                cq = pd.to_numeric(raw.iat[r, col + offset], errors="coerce")
                if pd.isna(cq):
                    continue
                records.append({
                    "construct": construct,
                    "target": target,
                    "dose_mM": float(dose),
                    "tech_rep": offset + 1,
                    "cq": float(cq),
                    "rt_minus_cq": float(rt_minus) if pd.notna(rt_minus) else np.nan,
                })
    return records


# --- recovering an unannotated export from the recorded plate layout ----------


@dataclass(frozen=True)
class QpcrPlate:
    """The qPCR plate layout as the logbook records it, identical across all runs.

    Args:
        construct_columns: Construct -> its three plate columns.
        roles: What each column within a block is, in order.
        reference_rows: Row -> dose, for the reference-gene half of the plate.
        target_rows: Row -> dose, for the stress-gene half.
        empty_rows: Rows the layout leaves unused, checked against the data.
    """

    construct_columns: dict[str, tuple[int, int, int]]
    roles: tuple[str, ...]
    reference_rows: dict[str, float]
    target_rows: dict[str, float]
    empty_rows: tuple[str, ...] = ("G", "H")


QPCR_PLATE = QpcrPlate(
    construct_columns={
        "UPRE1": (1, 2, 3),
        "UPRE2": (4, 5, 6),
        "NativeYap1": (7, 8, 9),
        "AlteredYap1": (10, 11, 12),
    },
    roles=("tech1", "tech2", "rt_minus"),
    reference_rows={"A": 0.0, "B": 0.2, "C": 0.5},
    target_rows={"D": 0.0, "E": 0.2, "F": 0.5},
)


def read_positional_cq(path, plate: QpcrPlate = QPCR_PLATE) -> pd.DataFrame:
    """Annotate a bare Cq export from the recorded plate layout.

    The 2026-07-24 run was exported with Target and Sample blank, so it has no
    analysis sheet. The layout is recorded and identical across runs, so position
    supplies what the export omitted.

    Applying a layout by position is only safe if the data agree with it, so three
    things are checked before anything is returned: the rows the layout leaves empty
    really are empty, the reference gene really does amplify earlier than its target,
    and the export is not one that already carries its own annotation.

    Args:
        path: Bio-Rad style export with ``Well`` and ``Cq`` columns.
        plate: Layout to apply.

    Raises:
        ValueError: if the export is already annotated, or the data contradict the
            layout.
    """
    raw = pd.read_excel(path, sheet_name=0, header=None)
    header = [str(x).strip() for x in raw.iloc[0]]
    body = raw.iloc[1:].copy()
    body.columns = header

    def column(name):
        for c in body.columns:
            if c.lower() == name:
                return c
        raise ValueError(f"export has no {name!r} column; this is not a Cq export")

    well_col, cq_col = column("well"), column("cq")
    if "target" in [h.lower() for h in header]:
        target_col = column("target")
        if body[target_col].notna().any():
            raise ValueError(
                "this export is already annotated with targets; read it with "
                "read_annotated_cq rather than applying a layout by position"
            )

    frame = body[[well_col, cq_col]].rename(columns={well_col: "well", cq_col: "cq"})
    frame = frame.dropna(subset=["well"])
    frame["row"] = frame.well.astype(str).str[0]
    frame["col"] = pd.to_numeric(frame.well.astype(str).str[1:], errors="coerce")
    frame["cq"] = pd.to_numeric(frame.cq, errors="coerce")
    grid = frame.pivot_table(index="row", columns="col", values="cq")

    occupied = [r for r in plate.empty_rows if r in grid.index and grid.loc[r].notna().any()]
    if occupied:
        raise ValueError(
            f"the layout says rows G/H are unused but {occupied} carry readings; "
            "this plate was not run with the recorded layout"
        )

    records = []
    for construct, cols in plate.construct_columns.items():
        stress_gene = TARGET_FOR_CONSTRUCT[construct]
        for gene, rows in ((REFERENCE_TARGET, plate.reference_rows),
                           (stress_gene, plate.target_rows)):
            for row, dose in rows.items():
                if row not in grid.index:
                    continue
                values = {role: grid.loc[row, col] if col in grid.columns else np.nan
                          for role, col in zip(plate.roles, cols)}
                rt_minus = values.get("rt_minus", np.nan)
                for rep, role in enumerate(("tech1", "tech2"), start=1):
                    cq = values.get(role, np.nan)
                    if pd.isna(cq):
                        continue
                    records.append({
                        "construct": construct, "target": gene, "dose_mM": float(dose),
                        "tech_rep": rep, "cq": float(cq),
                        "rt_minus_cq": float(rt_minus) if pd.notna(rt_minus) else np.nan,
                    })

    tidy = pd.DataFrame(records)
    reference = tidy[tidy.target == REFERENCE_TARGET].cq.mean()
    target = tidy[tidy.target != REFERENCE_TARGET].cq.mean()
    if not reference < target:
        raise ValueError(
            f"the reference gene ({reference:.1f}) does not amplify earlier than the "
            f"stress genes ({target:.1f}); the rows are not laid out as recorded"
        )
    return tidy
