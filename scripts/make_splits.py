"""Write the committed split manifest: which group is held out, under which claim.

Nothing in this repository has a named partition, so no two held-out numbers were scored
on the same one and none can be checked. This writes ``outputs/split_manifest.csv`` --
one row per group, carrying its assignment, the seed, the group key and a hash of the
whole partition -- plus ``outputs/split_manifest_summary.csv`` with the counts per kind.
A result cites a row; a reader re-runs this and compares hashes.

The inventories are built from the data itself, not from a hand-written list of
conditions. The real one is read off the plate exports, so the fluorescence coverage that
``docs/DATA_INVENTORY.md`` documents -- replicates 2 and 4 read columns 1-3 only -- shows
up in the counts rather than having to be remembered. The simulated one comes from
``generator/panel_experiment.py``.

Splits that the data cannot support are reported as refusals with the count behind them,
and are absent from the manifest. That absence is the result: see ``docs/SPLITS.md``.

The coverage description above records the original exports. The current public
``data/plates/manifest.csv`` selects the required NewProtocol population. Full scope is
always real biosensor plus simulated panel; a missing input is a refusal, not permission
to change that scope. ``--panel-only`` is the explicit simulated-only request.

Usage: python scripts/make_splits.py
Public replay: python scripts/make_splits.py --real-table outputs/sensor_characterisation.csv
"""
from __future__ import annotations

import argparse
import math
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.splits import (
    SPLIT_KINDS,
    SplitNotPossible,
    feasibility,
    make_split,
    split_manifest,
    summarise_manifest,
)
from ystwin.generator.panel_experiment import panel_dataset
from ystwin.plate import replay
from ystwin.plate.layout import RECORDED_PLATES, blanks_for_export, plate_key
from ystwin.plate.synergy import read_synergy_kinetic
from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

OUT = paths.outputs_dir()
SEED = 0

SEEDS = (0, 1, 2)
"""Seeds to emit partitions for.

Only some kinds are seed-sensitive. `extrapolation` is deterministic -- the top dose is
the top dose -- so every seed gives the same partition and the same hash. `interpolation`
picks which interior rung to withhold, so it genuinely differs, and a score from one seed
is a score against one held-out dose.
"""

REAL = "real_biosensor"
PANEL = "simulated_panel"
REAL_SOURCE_SET = "newprotocol"
REAL_KEY = ("plate", "construct", "stressor", "dose_mM", "well")

# Which combination the simulated panel co-applies. DTT and H2O2 are the two agents the
# wet lab actually used, so the one interaction worth simulating is theirs.
COMBINATION = ("DTT", "H2O2")


# ---------------------------------------------------------------------------
# inventories
# ---------------------------------------------------------------------------


def _recorded_rows(date: str, read_wells: set[str], has_blank: bool) -> list[dict]:
    """Place only observed reporter wells using the recorded construct/dose layout."""
    recorded = RECORDED_PLATES[date]
    rows = []
    for construct in sorted(recorded.layout.construct_columns):
        stressor = STRESSOR_FOR_CONSTRUCT[construct]
        for index, dose in enumerate(recorded.doses_mM[stressor]):
            for well in recorded.layout.wells(construct, index):
                if well in read_wells:
                    rows.append({
                        "plate": date,
                        "construct": construct,
                        "stressor": stressor,
                        "dose_mM": float(dose),
                        "well": well,
                        "plate_has_reporter_blank": has_blank,
                    })
    return rows


def real_inventory(directory: pathlib.Path) -> pd.DataFrame:
    """One row per well that was actually read in the reporter channel.

    Coverage is read from the export rather than assumed from the protocol. Two of the
    four NewProtocol plates were configured to fluoresce columns 1-3 only, which under the
    recorded layout is UPRE1, so those plates contribute a full UPRE1 ladder and nothing
    for the other three constructs. A split built on the protocol's nominal design would
    claim four biological replicates everywhere and be wrong for three constructs.

    Args:
        directory: Where the NewProtocol exports live.

    Returns:
        Columns ``plate``, ``construct``, ``stressor``, ``dose_mM``, ``well`` and
        ``plate_has_reporter_blank`` -- the last saying whether this plate's recorded
        blank wells were themselves read in the reporter channel, which is what decides
        whether the plate can be background-corrected at all.
    """
    rows = []
    for path in sorted(directory.glob("*.xlsx")):
        date = plate_key(path.name)
        recorded = RECORDED_PLATES.get(date)
        if recorded is None:
            continue
        try:
            run = read_synergy_kinetic(path)
            reporter = run.raw_channel("mCitrine")
        except (KeyError, ValueError):
            # An OD-only export -- the reporter-free BY4741 plates are like this. Nothing
            # to split: the channel every result rests on was not read.
            continue
        read_wells = set(reporter.wells)
        blanks = blanks_for_export(path.name) or ()
        has_blank = bool(read_wells.intersection(blanks))
        rows.extend(_recorded_rows(date, read_wells, has_blank))
    if not rows:
        raise SystemExit(
            f"no export in {directory} was read in the mCitrine channel; there is "
            f"nothing to partition and {OUT / 'split_manifest.csv'} was left untouched"
        )
    return pd.DataFrame(rows)


def _committed_real_inventory() -> pd.DataFrame:
    """The required real population, identified by source_set and observed well IDs.

    Superseded exports and other plate populations are not interchangeable with the
    current NewProtocol replicates. The manifest selects the exports; their reporter
    columns and the logbook place the wells, rather than a nominal well or plate count.
    This is an identity check, not a raw-fluorescence selection or correction. Multiple
    non-derived reporter blocks must agree on well coverage; their numerical correction
    state remains the upstream sensor producer's responsibility.
    """
    try:
        manifest = replay.load_manifest()
    except FileNotFoundError as exc:
        raise SystemExit(f"required real-population manifest is unavailable: {exc}") from exc
    missing_columns = {"source_set", "export"} - set(manifest.columns)
    if missing_columns:
        raise SystemExit(f"real-population manifest is missing columns {sorted(missing_columns)}")
    exports = replay.exports(manifest=manifest, source_set=REAL_SOURCE_SET)
    if not exports:
        raise SystemExit(f"real-population manifest has no source_set={REAL_SOURCE_SET!r} exports")
    rows = []
    seen = set()
    for export in exports:
        date = plate_key(export.name)
        if date is None:
            raise SystemExit(f"no recorded plate identity for required export {export.name!r}")
        if date in seen:
            raise SystemExit(f"multiple current exports declare plate {date}; resolve the manifest identity")
        seen.add(date)
        try:
            run = replay.load_run(export.name, manifest=manifest)
        except FileNotFoundError as exc:
            raise SystemExit(f"required real-population readings are unavailable: {exc}") from exc
        reporters = [block for block in run.blocks
                     if block.fluorophore == "mCitrine" and not block.derived]
        if not reporters:
            raise SystemExit(f"required export {export.name!r} has no non-derived mCitrine well inventory")
        wells = set(reporters[0].wells)
        if any(set(block.wells) != wells for block in reporters[1:]):
            raise SystemExit(f"required export {export.name!r} has ambiguous reporter well coverage; "
                             "declare the intended channel population before making splits")
        has_blank = bool(wells.intersection(RECORDED_PLATES[date].blank_wells))
        observed = _recorded_rows(date, wells, has_blank)
        if not observed:
            raise SystemExit(f"required export {export.name!r} has no recorded reporter culture wells")
        rows.extend(observed)
    return pd.DataFrame(rows)


def _require_real_population(real: pd.DataFrame, expected: pd.DataFrame, source) -> None:
    """Refuse identity loss or substitution, including equal-sized wrong inventories."""
    duplicated = real.duplicated(["plate", "well"], keep=False)
    if duplicated.any():
        identities = real.loc[duplicated, ["plate", "well"]].to_dict("records")
        raise SystemExit(f"{source} repeats physical plate/well identities: {identities[:5]}")
    actual_ids = set(real[list(REAL_KEY)].itertuples(index=False, name=None))
    expected_ids = set(expected[list(REAL_KEY)].itertuples(index=False, name=None))
    missing = sorted(expected_ids - actual_ids)
    unexpected = sorted(actual_ids - expected_ids)
    if missing or unexpected:
        raise SystemExit(
            f"{source} does not cover the required source_set={REAL_SOURCE_SET!r} "
            f"real population: missing well identities {missing[:5]}; "
            f"unexpected well identities {unexpected[:5]}. "
            "Supply the complete real input; no split outputs were written."
        )


def real_table_inventory(source: pathlib.Path) -> pd.DataFrame:
    """Read an explicit sensor_characterisation table without changing its population.

    Measurements are not used to choose a split: every input well, including wells with
    unavailable numerical estimates, must remain. Canonical plate keys come from the
    shared logbook lookup, preserving the existing group IDs and partition hashes.
    Reporter-blank coverage comes from the public exports, not an assumption about CSV
    row counts or about which plates used to carry a blank.
    """
    if not source.is_file():
        raise SystemExit(f"required --real-table input is not a file: {source}; no split outputs were written")
    try:
        table = pd.read_csv(source, dtype={key: str for key in REAL_KEY if key != "dose_mM"})
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
        raise SystemExit(f"invalid --real-table input {source}: {exc}") from exc
    missing_columns = set(REAL_KEY) - set(table.columns)
    if missing_columns:
        raise SystemExit(f"--real-table {source} is missing required columns {sorted(missing_columns)}")
    real = table[list(REAL_KEY)].copy()
    if real.empty:
        raise SystemExit(f"--real-table {source} has no real well identities")
    if real.isna().any().any():
        raise SystemExit(f"--real-table {source} has missing real well identity values")
    for column in (key for key in REAL_KEY if key != "dose_mM"):
        if real[column].str.strip().eq("").any():
            raise SystemExit(f"--real-table {source} has empty {column} identities")
    try:
        real["dose_mM"] = pd.to_numeric(real.dose_mM, errors="raise")
    except ValueError as exc:
        raise SystemExit(f"--real-table {source} needs numeric dose_mM identities") from exc
    if not real.dose_mM.map(math.isfinite).all() or real.dose_mM.lt(0).any():
        raise SystemExit(f"--real-table {source} needs finite, nonnegative dose_mM identities")
    canonical = real.plate.map(plate_key)
    if canonical.isna().any():
        raise SystemExit(f"--real-table {source} has unrecorded plate identities: "
                         f"{sorted(set(real.loc[canonical.isna(), 'plate']))}")
    real["plate"] = canonical
    expected = _committed_real_inventory()
    _require_real_population(real, expected, source)
    return real.merge(expected, on=list(REAL_KEY), how="left", validate="one_to_one")


def panel_inventory() -> pd.DataFrame:
    """One row per simulated well, with the same column vocabulary as the real frame.

    The panel has no plate and no construct: replicates are wells of one simulated plate
    and the readings are channels, not promoters. That absence is why the replicate and
    construct splits are unavailable here and the stressor and combination splits are
    unavailable on the real data -- the two datasets support disjoint claims, which is
    worth seeing in one table.

    Doses are in each stressor's own units rather than mM throughout. The column keeps
    the shared name because every split is drawn inside one ladder, so no two ladders are
    ever compared.
    """
    data = panel_dataset(replicates=3, seed=SEED, combinations=[COMBINATION])
    return pd.DataFrame({
        "stressor": [str(label) for label in data.labels],
        "dose_mM": [float(dose) for dose in data.doses],
        "well": list(data.wells),
    })


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------


def attempt(splits, refusals, frame, kind, dataset, scope="", **kwargs):
    """Build one split, or record why the data cannot support it.

    Args:
        splits: Collected splits, appended to on success.
        refusals: Collected refusals, appended to on failure.
        frame: The inventory to partition.
        kind: Split kind to attempt.
        dataset: Name recorded in the manifest.
        scope: What ``frame`` covers, for the refusal line -- a kind can be impossible
            over the whole dataset and available over part of it, and a report that does
            not say which is being talked about is worse than no report.
        **kwargs: Passed through to :func:`ystwin.analysis.splits.make_split`.
    """
    try:
        for seed in SEEDS:
            splits.append(make_split(frame, kind, seed=seed, dataset=dataset, **kwargs))
        split = splits[-1]
    except SplitNotPossible as exc:
        refusals.append({"dataset": dataset, "split_kind": kind, "scope": scope,
                         "reason": str(exc)})
        print(f"  {kind:22s} REFUSED{f' over {scope}' if scope else ''}")
        print(f"    {exc}")
        return None
    print(f"  {split.summary()}")
    return split


def report_feasibility(frame, kind, dataset, **kwargs):
    """Print the per-stratum counts behind a refusal, which is the citable part."""
    try:
        table = feasibility(frame, kind, seed=SEED, **kwargs)
    except SplitNotPossible as exc:
        print(f"  {kind:22s} not applicable to {dataset}: {exc}".rstrip())
        return None
    for _, row in table.iterrows():
        verdict = "available" if row.feasible else "BLOCKED"
        stratum = row.stratum or "(all)"
        print(f"    {stratum:24s} n={row.n_units:<3d} {verdict}"
              + (f" -- {row.reason}" if row.reason else ""))
    return table


# ---------------------------------------------------------------------------


def _refuse_to_shrink_the_manifest(panel_only: bool) -> None:
    """Require full-scope inputs even when the destination has no previous manifest.

    Historical rationale for the original destination-only guard, now superseded:
    Stop before overwriting a manifest that has real rows with one that has none.

    Without the plate exports the real half of the inventory is simply absent, and the
    run then writes a panel-only manifest over the tracked one and exits 0. That is the
    worst kind of failure a tracked artefact can have: 714 of 1650 rows disappear, the
    console says "wrote", and the next script to read the file dies somewhere unrelated
    with ``'DataFrame' object has no attribute 'plate'``.

    The rule is narrow on purpose. Refuse only when writing would *destroy* something --
    an existing manifest carrying real rows. A fresh checkout with no manifest, or one
    that never had real rows, can build the simulated half freely.

    The current contract is the requested scope, not the destination's contents. A
    fresh staging directory must not turn a full run into a successful panel-only run.
    The simulated half remains available, but only with an explicit --panel-only.

    Args:
        panel_only: Set by ``--panel-only`` to say the shrink is intended.
    """
    if panel_only:
        return
    raise SystemExit(
        f"full scope requires {REAL} inputs, but the plate exports are not reachable; "
        f"{OUT / 'split_manifest.csv'} and its summary were left untouched. "
        "Pass --real-table outputs/sensor_characterisation.csv for the public inventory, "
        "or set YSTWIN_PLATES to the complete NewProtocol exports. "
        "Use --panel-only only when a simulated-only manifest is explicitly intended."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--panel-only", action="store_true",
                        help="explicitly request only the simulated population; never read real inputs")
    source.add_argument("--real-table", type=pathlib.Path,
                        help="sensor_characterisation CSV covering the complete current NewProtocol "
                             "well population; required for a public full-scope run without raw exports")
    args = parser.parse_args()

    real = None
    real_source = args.real_table
    if args.real_table is not None:
        real = real_table_inventory(args.real_table)
    elif not args.panel_only:
        real_source = paths.biosensor_plates()
        if real_source is None:
            _refuse_to_shrink_the_manifest(args.panel_only)
        real = real_inventory(real_source)
        _require_real_population(real, _committed_real_inventory(), real_source)

    splits = []
    refusals: list[dict] = []

    print("=" * 88)
    print("Simulated panel  (generator/panel_experiment.py)")
    print("=" * 88)
    panel = panel_inventory()
    print(f"  {len(panel)} wells, {panel.stressor.nunique()} treatments, "
          f"{panel.groupby('stressor').dose_mM.nunique().max()} rungs per ladder")
    attempt(splits, refusals, panel, "interpolation", PANEL,
            group_key=("stressor", "dose_mM"), within=("stressor",))
    attempt(splits, refusals, panel, "extrapolation", PANEL,
            group_key=("stressor", "dose_mM"), within=("stressor",))
    attempt(splits, refusals, panel, "heldout_stressor", PANEL,
            group_key=("stressor",), validation_units=1)
    attempt(splits, refusals, panel, "heldout_combination", PANEL,
            group_key=("stressor",))
    print("  heldout_replicate / heldout_construct: not applicable -- a simulated panel")
    print("    has no plate and no promoter, so neither axis exists to hold out.")

    if real is not None:
        print("\n" + "=" * 88)
        print(f"Real biosensor plates  ({real_source})")
        print("=" * 88)
        coverage = real.groupby("construct").plate.nunique().sort_index()
        print("  wells read in the reporter channel, by construct and plate:")
        print("   " + real.pivot_table(index="construct", columns="plate",
                                       values="well", aggfunc="count",
                                       fill_value=0).to_string().replace("\n", "\n   "))
        print(f"  biological replicates with fluorescence: {dict(coverage)}")

        attempt(splits, refusals, real, "interpolation", REAL,
                group_key=("plate", "construct", "dose_mM"))
        attempt(splits, refusals, real, "extrapolation", REAL,
                group_key=("plate", "construct", "dose_mM"))
        attempt(splits, refusals, real, "heldout_construct", REAL,
                group_key=("plate", "construct", "dose_mM"))
        attempt(splits, refusals, real, "heldout_replicate", REAL,
                scope="all four constructs",
                group_key=("plate", "construct", "dose_mM"))
        attempt(splits, refusals, real, "heldout_stressor", REAL,
                scope="the whole plate set",
                group_key=("plate", "construct", "stressor", "dose_mM"))

        print("\n  Replicate-level feasibility, per construct")
        print("  all plates read in the reporter channel:")
        replicate_key = ("plate", "construct", "dose_mM")
        table = report_feasibility(real, "heldout_replicate", REAL,
                                   group_key=replicate_key)
        corrigible = real[real.plate_has_reporter_blank]
        print("  plates whose blanks were also read, which is the set every current "
              "result uses:")
        report_feasibility(corrigible, "heldout_replicate", REAL,
                           group_key=replicate_key)

        # Historical rationale; feasibility and blank coverage now come from the input.
        # One construct can support the split, on plates whose background was never
        # measured. That partition is worth having in the manifest -- it is the only
        # forward-looking one the real data offers -- but the caveat travels with it.
        usable = [] if table is None else [
            row.stratum.split("=", 1)[1] for _, row in table.iterrows() if row.feasible
        ]
        if usable:
            print(f"\n  Building the replicate split for {usable} alone.")
            without_blanks = sorted(set(real.loc[
                real.construct.isin(usable) & ~real.plate_has_reporter_blank, "plate"]))
            if without_blanks:
                print(f"  Plates {without_blanks} carry no blank in the reporter channel; "
                      "a score on this partition requires a stated background bound.")
            else:
                print("  Every required plate has a recorded blank read in the reporter channel.")
            attempt(splits, refusals, real[real.construct.isin(usable)],
                    "heldout_replicate", REAL, group_key=replicate_key)

    manifest = split_manifest(splits)
    requested = {PANEL} if args.panel_only else {PANEL, REAL}
    missing_populations = requested - set(manifest.dataset)
    if missing_populations:
        raise SystemExit(f"no supported splits for requested populations {sorted(missing_populations)}; "
                         "no split outputs were written")
    summary = summarise_manifest(manifest)
    manifest.to_csv(OUT / "split_manifest.csv", index=False)
    summary.to_csv(OUT / "split_manifest_summary.csv", index=False)

    print("\n" + "=" * 88)
    print(f"Manifest: {len(manifest)} groups across {len(splits)} splits, seed {SEED}")
    print("=" * 88)
    shown = summary[["dataset", "split_kind", "assignment", "n_groups", "n_rows"]]
    print(shown.to_string(index=False))
    print("\n  partition hashes")
    for split in splits:
        print(f"    {split.dataset:16s} {split.kind:22s} {split.hash}")
    if refusals:
        print("\n  refused")
        for refusal in refusals:
            scope = f"  (over {refusal['scope']})" if refusal["scope"] else ""
            print(f"    {refusal['dataset']:16s} {refusal['split_kind']:22s}{scope}")
    unused = sorted(set(SPLIT_KINDS) - {s.kind for s in splits})
    if unused:
        print(f"  no partition at all in this manifest for: {unused}")
    print("  see docs/SPLITS.md for what each refusal costs and what would lift it.")
    print(f"\n  wrote {OUT / 'split_manifest.csv'}")
    print(f"  wrote {OUT / 'split_manifest_summary.csv'}")


if __name__ == "__main__":
    main()
