"""Held-out splits that respect the grouping this data actually has.

Every held-out number in this repository is currently produced by a split written inline
in whichever script needed one. No two of them partition the data the same way, none is
recorded anywhere, and none can be re-created by a reader. A held-out R2 quoted against an
unnamed partition is not a claim about generalisation; it is a claim about an arrangement
of rows that no longer exists. This module makes the partition an artefact: named, seeded,
counted, hashed, and written to ``outputs/split_manifest.csv`` so a result can cite the
exact one it was scored on.

The kinds here are derived from what these two datasets can support, not copied from
another project's list. Two constraints do most of the deriving.

**The unit of replication is the plate.** Wells on one plate share an inoculum, a medium
batch, a reader and a position in an incubator, so a split that puts wells from one plate
on both sides has trained on the batch effect it is about to be scored against.
``analysis/uncertainty.py`` refuses that mistake on the analysis side and
``analysis/power.py`` on the design side; this is the matching refusal for evaluation. It
is enforced structurally rather than by convention: the caller names the group key, and a
split whose groups straddle the axis being held out is refused rather than permitted.

**The two sensor pairs were challenged with different agents.** UPRE1 and UPRE2 saw DTT,
NativeYap1 and AlteredYap1 saw H2O2, and ``qpcr.assert_single_stressor`` already refuses to
pool their dose axes. A construct-level split inherits that refusal, and the frame's own
construct-to-stressor labelling is checked against ``qpcr.STRESSOR_FOR_CONSTRUCT`` rather
than trusted, so a mislabelled inventory is caught before it becomes a partition.

The kinds, and what a score on each licenses:

``interpolation``
    Hold out one *interior* rung of a dose ladder. Licenses "the response is smooth
    between the doses we ran". The lowest rung is never the held-out one -- predicting
    below the training range is extrapolation downward wearing the wrong name, and on the
    real ladder the lowest rung is the 0 mM control that is the denominator of every fold
    change, so withholding it leaves a training set that cannot express the quantity being
    predicted. Four rungs are required: holding out the middle of three leaves two, and
    two points cannot distinguish a curve from a line, so the split would not test the
    thing it is named for.

``extrapolation``
    Hold out the *top* dose. Licenses "the response continues above the doses we ran",
    which is a different and much harder claim, and the one the dose-response work in this
    repository has the most reason to doubt: the top of every ladder is where growth
    arrest, dilution artefacts and the unmeasured autofluorescence term all concentrate
    (see ``docs/DATA_INVENTORY.md``). Conflating it with interpolation is the specific
    error this module exists to prevent, so it is a separate kind with a separate hash and
    it is deterministic -- the top dose is the top dose, and no seed moves it.

``heldout_replicate``
    Hold out one whole biological replicate, which here means one plate. Licenses "this
    would have worked on a plate we had not run", the only genuinely forward-looking claim
    available on real data. The training side must retain at least
    ``uncertainty.MIN_PLATES_FOR_INTERVAL`` plates, because a claim about generalising
    across batches rests on having measured batch variation, and this package's own answer
    to how few clusters admit a statement is three. That makes four plates the minimum for
    the split, and on the current dataset it is available for no construct -- see
    ``docs/SPLITS.md``.

``heldout_construct``
    Train on UPRE1, predict UPRE2: same stressor, different promoter. Licenses "the model
    is about the pathway, not about the reporter it was fitted to". Valid only *within* a
    stressor. One training construct is enough here, unlike one training plate, because
    this is a named one-to-one transfer claim whose scope is stated by naming both
    constructs -- while a single training plate leaves the very quantity the replicate
    split is about, batch spread, unmeasured.

``heldout_stressor``
    Hold out a whole stressor. The existing leave-one-stressor-out of
    ``analysis/transfer.py``, given a manifest. Refused when constructs are perfectly
    nested inside stressors, which is exactly the real biosensor design: holding out H2O2
    there also holds out NativeYap1 and AlteredYap1 entirely, so nothing in the result
    separates "unseen stressor" from "unseen reporter". Two stressors are not enough
    either -- a latent state is a subspace, and one training stressor is one ray. When
    the design co-applies agents, holding out a stressor takes every treatment containing
    it: leaving ``DTT+H2O2`` in training while withholding ``H2O2`` withholds a label and
    not an agent, and the leak does not show up in the counts.

``heldout_combination``
    Hold out a co-applied stressor pair while both its components remain in training.
    Licenses "the interaction is predicted from the parts", which no single-agent split
    tests. Refused unless every component is present as a single in the training set,
    since otherwise it is quietly an extrapolation. Simulated panel only today: no real
    plate here co-applied two agents.

Deliberately absent: a plain random split of wells into train and test. It leaks the plate
and it leaks the dose group, it scores near-perfectly for that reason, and it licenses
nothing at all. ``validation`` is not a kind either but an *assignment* -- further units
moved out of training by the same grouping discipline, so that model selection never
touches the held-out set.

Column vocabulary is fixed rather than parameterised: ``construct``, ``stressor``,
``plate``, ``dose_mM``. A manifest is only citable if every producer of one means the same
thing by its columns, and a rename parameter is how two tables end up disagreeing about
what a group is.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..qpcr import STRESSOR_FOR_CONSTRUCT, assert_single_stressor
from .uncertainty import MIN_PLATES_FOR_INTERVAL

__all__ = [
    "CONSTRUCT_COLUMN",
    "DOSE_COLUMN",
    "PLATE_COLUMN",
    "SPLIT_KINDS",
    "STRESSOR_COLUMN",
    "Split",
    "SplitKind",
    "SplitNotPossible",
    "TEST",
    "TRAIN",
    "VALIDATION",
    "feasibility",
    "make_split",
    "partition_hash",
    "split_manifest",
    "summarise_manifest",
]

CONSTRUCT_COLUMN = "construct"
STRESSOR_COLUMN = "stressor"
PLATE_COLUMN = "plate"
DOSE_COLUMN = "dose_mM"

TRAIN = "train"
VALIDATION = "validation"
TEST = "test"

_COMBINATION_SEPARATOR = "+"
"""How ``generator/panel_experiment.py`` writes a co-applied treatment's label.

Read from the generator rather than guessed: it builds the label as
``"+".join(members)``, so the separator is a fact about the producer of the data and not
a pattern hoped for in a string.
"""


class SplitNotPossible(ValueError):
    """This dataset cannot support the requested split, and here is the count that says so.

    A subclass of :class:`ValueError` so that a caller catching either sees it. It is
    raised rather than worked around: a split that quietly shrinks to whatever the data
    allows produces a number under a name that no longer describes it.
    """


# ---------------------------------------------------------------------------
# the kinds
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SplitKind:
    """One way of holding data out, and the conditions under which it means anything.

    Args:
        name: How the kind is written in the manifest and cited in a result.
        axis: Column whose values are withheld. One value of it goes to test.
        within: Columns defining the stratum the split happens inside. Each stratum is
            split independently, because a dose ladder belongs to a construct and a
            plate count belongs to a construct -- a split pooled across them would hold
            out a dose from one construct and nothing from another.
        rule: How the held-out value is chosen. ``"top"`` takes the largest,
            ``"interior"`` draws from strictly between the extremes, ``"random"`` draws
            from all of them, ``"single"`` draws from the single-agent labels and takes
            every co-applied treatment containing the chosen agent with it, and
            ``"combination"`` draws from the co-applied labels only.
        min_units: Distinct axis values a stratum needs before the split is offered.
        min_train_units: Distinct axis values that must remain in training.
        claim: What a score on this split licenses, in one sentence.
    """

    name: str
    axis: str
    within: tuple[str, ...]
    rule: str
    min_units: int
    min_train_units: int
    claim: str

    @property
    def deterministic(self) -> bool:
        """Whether the seed moves the partition at all.

        False for ``extrapolation``: the top dose is the top dose. Two runs of it with
        different seeds must hash identically, and a kind that quietly randomised it
        would be an interpolation split under an extrapolation label.
        """
        return self.rule == "top"


SPLIT_KINDS: dict[str, SplitKind] = {
    "interpolation": SplitKind(
        name="interpolation",
        axis=DOSE_COLUMN,
        within=(CONSTRUCT_COLUMN,),
        rule="interior",
        min_units=4,
        min_train_units=3,
        claim="the dose response is smooth between the rungs that were run",
    ),
    "extrapolation": SplitKind(
        name="extrapolation",
        axis=DOSE_COLUMN,
        within=(CONSTRUCT_COLUMN,),
        rule="top",
        min_units=3,
        min_train_units=2,
        claim="the dose response continues above the top rung that was run",
    ),
    "heldout_replicate": SplitKind(
        name="heldout_replicate",
        axis=PLATE_COLUMN,
        within=(CONSTRUCT_COLUMN,),
        rule="random",
        min_units=MIN_PLATES_FOR_INTERVAL + 1,
        min_train_units=MIN_PLATES_FOR_INTERVAL,
        claim="the result would have held on a plate that was never run",
    ),
    "heldout_construct": SplitKind(
        name="heldout_construct",
        axis=CONSTRUCT_COLUMN,
        within=(STRESSOR_COLUMN,),
        rule="random",
        min_units=2,
        min_train_units=1,
        claim="the model is about the pathway rather than the promoter it was fitted to",
    ),
    "heldout_stressor": SplitKind(
        name="heldout_stressor",
        axis=STRESSOR_COLUMN,
        within=(),
        rule="single",
        min_units=3,
        min_train_units=2,
        claim="a latent stress state reaches a stressor it never saw",
    ),
    "heldout_combination": SplitKind(
        name="heldout_combination",
        axis=STRESSOR_COLUMN,
        within=(),
        rule="combination",
        min_units=1,
        min_train_units=2,
        claim="the interaction between two agents is predicted from the single agents",
    ),
}


# ---------------------------------------------------------------------------
# identity: group ids, stratum ids, and the hash over the partition
# ---------------------------------------------------------------------------


def _format(value: object) -> str:
    """One value as it appears in a group id.

    Floats go through ``%.12g`` so that a dose read from two different exports lands on
    the same id, and so that ``0.5`` never becomes ``0.5000000000000001`` in a hash.
    """
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.12g}"
    if isinstance(value, (np.integer,)):
        return str(int(value))
    return str(value)


def _identifier(columns: Sequence[str], values: Sequence[object]) -> str:
    """``construct=UPRE1|dose_mM=0.5`` -- readable, sortable, and unambiguous."""
    return "|".join(f"{c}={_format(v)}" for c, v in zip(columns, values))


def partition_hash(assignment: Mapping[str, str]) -> str:
    """A stable digest of a group-to-assignment mapping.

    Two runs that produce this digest partitioned the data identically, whatever seed,
    machine or pandas version produced them. The mapping is sorted by group id before
    hashing, so nothing here depends on a dict's insertion order or on the order the
    rows arrived in -- which is the whole point, since a hash that moved with row order
    would prove nothing.

    Args:
        assignment: Group id to ``"train"``, ``"validation"`` or ``"test"``.

    Returns:
        A 64-character SHA-256 hex digest. Quoting its first 12 characters is enough to
        distinguish partitions in prose; the manifest carries all of it.
    """
    payload = "\n".join(f"{group}\t{assignment[group]}" for group in sorted(assignment))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _stratum_rng(seed: int, kind: str, stratum: str) -> np.random.Generator:
    """A generator whose stream depends on the stratum but not on its neighbours.

    Drawing every stratum's fold from one sequential generator would make UPRE2's
    held-out dose change because a fifth construct was added to the frame. Deriving the
    stream from ``(seed, kind, stratum)`` keeps each stratum's fold stable under changes
    elsewhere, which is what lets a manifest be regenerated after the dataset grows.
    """
    digest = hashlib.sha256(f"{seed}\t{kind}\t{stratum}".encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


# ---------------------------------------------------------------------------
# the split
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Split:
    """A partition of groups, with everything needed to cite and re-create it.

    Args:
        kind: Which :class:`SplitKind` produced it.
        seed: The seed it was drawn with. Recorded even for a deterministic kind, so
            that a manifest row is self-describing.
        group_key: Columns whose distinct combinations are the indivisible units.
        within: Columns that defined the strata.
        dataset: Name of the dataset this partitions, for the manifest.
        groups: One row per group -- its id, stratum, axis value, assignment and row
            count -- alongside the group key columns themselves.
    """

    kind: str
    seed: int
    group_key: tuple[str, ...]
    within: tuple[str, ...]
    dataset: str
    groups: pd.DataFrame = field(repr=False)

    @property
    def spec(self) -> SplitKind:
        return SPLIT_KINDS[self.kind]

    @property
    def assignment(self) -> dict[str, str]:
        """Group id to assignment, in sorted order."""
        pairs = zip(self.groups["group_id"], self.groups["assignment"])
        return {group: value for group, value in sorted(pairs)}

    @property
    def hash(self) -> str:
        """:func:`partition_hash` of this split's mapping."""
        return partition_hash(self.assignment)

    def held_out(self) -> dict[str, tuple]:
        """Stratum to the axis values sent to test."""
        return self._values_by_stratum(TEST)

    def validation_values(self) -> dict[str, tuple]:
        """Stratum to the axis values sent to validation. Empty when none were carved."""
        return self._values_by_stratum(VALIDATION)

    def _values_by_stratum(self, assignment: str) -> dict[str, tuple]:
        chosen = self.groups[self.groups["assignment"] == assignment]
        out: dict[str, tuple] = {}
        for stratum, block in chosen.groupby("stratum", sort=True):
            out[str(stratum)] = tuple(sorted(set(block["axis_value"])))
        return out

    def counts(self) -> dict[str, dict[str, int]]:
        """Groups and rows on each side, which is what a reader checks first."""
        out: dict[str, dict[str, int]] = {}
        for assignment in (TRAIN, VALIDATION, TEST):
            block = self.groups[self.groups["assignment"] == assignment]
            out[assignment] = {
                "groups": int(len(block)),
                "rows": int(block["n_rows"].sum()),
            }
        return out

    def assign_rows(self, frame: pd.DataFrame) -> pd.Series:
        """Each row's assignment, aligned to ``frame``'s index.

        The manifest is the artefact, but a manifest nobody can apply is a document
        rather than a split. This is how a scoring script uses one.

        Args:
            frame: The frame the split was built from, or any frame with the same group
                key columns.

        Raises:
            SplitNotPossible: if a row falls in a group this split does not cover, which
                means the frame is not the one the manifest describes.
        """
        ids = _group_ids(frame, self.group_key)
        known = self.assignment
        missing = sorted(set(ids) - set(known))
        if missing:
            raise SplitNotPossible(
                f"{len(missing)} group(s) in this frame are not in the {self.kind} split, "
                f"starting with {missing[:3]}; the manifest describes a different frame"
            )
        return pd.Series([known[i] for i in ids], index=frame.index, name="assignment")

    def manifest(self) -> pd.DataFrame:
        """One row per group, carrying its own provenance.

        Every row names the seed, the kind, the group key and the partition hash, so a
        single row lifted out of the table still says which partition it came from. That
        redundancy is deliberate: a result cites a row, not a file.
        """
        out = self.groups.copy()
        out.insert(0, "dataset", self.dataset)
        out.insert(1, "split_kind", self.kind)
        out.insert(2, "seed", self.seed)
        out.insert(3, "group_key", "|".join(self.group_key))
        out.insert(4, "within", "|".join(self.within))
        out["axis"] = self.spec.axis
        # Per stratum, not per split: a row should say what was withheld from *its* own
        # ladder. Repeating all twenty-six of the panel's held-out doses on every row
        # would treble the file and tell a reader nothing about the row they are on.
        held = self.held_out()
        out["held_out"] = [
            ", ".join(_format(v) for v in held.get(str(stratum), ()))
            for stratum in out["stratum"]
        ]
        out["partition_hash"] = self.hash
        out["axis_value"] = [_format(v) for v in out["axis_value"]]
        return out.reset_index(drop=True)

    def summary(self) -> str:
        counts = self.counts()
        held = "; ".join(
            f"{stratum} -> {', '.join(_format(v) for v in values)}" if stratum else
            ", ".join(_format(v) for v in values)
            for stratum, values in sorted(self.held_out().items())
        )
        parts = [
            f"{self.kind} (seed {self.seed}, by {'|'.join(self.group_key)}): "
            f"{counts[TRAIN]['groups']} train / "
            f"{counts[VALIDATION]['groups']} validation / "
            f"{counts[TEST]['groups']} test groups, "
            f"{counts[TRAIN]['rows']}/{counts[VALIDATION]['rows']}/{counts[TEST]['rows']} rows",
            f"    held out: {held}",
            f"    hash: {self.hash[:12]}",
        ]
        return "\n".join(parts)


# ---------------------------------------------------------------------------
# grouping, and the refusals that protect it
# ---------------------------------------------------------------------------


def _group_ids(frame: pd.DataFrame, group_key: Sequence[str]) -> list[str]:
    """Group id for every row of ``frame``, in row order."""
    values = frame[list(group_key)].to_numpy(dtype=object)
    return [_identifier(group_key, row) for row in values]


def _require_columns(frame: pd.DataFrame, needed: Sequence[str], kind: str) -> None:
    missing = [c for c in dict.fromkeys(needed) if c not in frame.columns]
    if missing:
        raise SplitNotPossible(
            f"the {kind} split needs column(s) {missing}, which this frame does not have "
            f"(it has {sorted(frame.columns)}). That is usually not a missing column but "
            f"a dataset without the axis: a simulated panel has no plate, and a "
            f"single-construct plate has no construct to hold out."
        )
    empty = [c for c in dict.fromkeys(needed) if frame[c].isna().any()]
    if empty:
        raise SplitNotPossible(
            f"column(s) {empty} contain missing values; a row whose plate, construct, "
            f"stressor or dose is unrecorded cannot be assigned to a side, and guessing "
            f"one is how a leak gets in"
        )


def _check_construct_stressor_agreement(frame: pd.DataFrame) -> None:
    """Verify the frame's own labelling against the package's declared mapping.

    ``qpcr.STRESSOR_FOR_CONSTRUCT`` is the authoritative record of which agent each
    construct was challenged with. A frame that disagrees with it has been mislabelled
    somewhere upstream, and building a split on it would bake the mislabelling into a
    committed manifest.
    """
    if CONSTRUCT_COLUMN not in frame.columns or STRESSOR_COLUMN not in frame.columns:
        return
    for construct in sorted(set(frame[CONSTRUCT_COLUMN])):
        declared = STRESSOR_FOR_CONSTRUCT.get(construct)
        if declared is None:
            continue
        seen = set(frame.loc[frame[CONSTRUCT_COLUMN] == construct, STRESSOR_COLUMN])
        if seen != {declared}:
            raise SplitNotPossible(
                f"this frame labels {construct!r} with stressor(s) {sorted(seen)}, but "
                f"qpcr.STRESSOR_FOR_CONSTRUCT records {declared!r}; fix the frame rather "
                f"than the split"
            )


def _group_table(
    frame: pd.DataFrame,
    group_key: Sequence[str],
    axis: str,
    within: Sequence[str],
    kind: str,
) -> pd.DataFrame:
    """One row per group, refusing any group that straddles the axis or a stratum.

    This is where requirement "never split within a plate" is enforced. The group key is
    the caller's declaration of what an indivisible unit is; if a unit contains two
    plates, then assigning that unit puts one plate's wells on both sides no matter which
    side it goes to, so the split is refused rather than silently permitted.
    """
    key = list(group_key)
    extra = [c for c in dict.fromkeys([axis, *within]) if c not in key]
    counts = frame.groupby(key, sort=True, dropna=False).size().rename("n_rows").reset_index()
    variants = frame[key + extra].drop_duplicates()
    if len(variants) != len(counts):
        sizes = variants.groupby(key, sort=True, dropna=False).size()
        straddling = sizes[sizes > 1]
        offenders = []
        for name in list(straddling.index)[:3]:
            values = name if isinstance(name, tuple) else (name,)
            mask = np.ones(len(variants), dtype=bool)
            for column, value in zip(key, values):
                mask &= (variants[column] == value).to_numpy()
            spread = {c: sorted(set(variants.loc[mask, c])) for c in extra}
            offenders.append(f"{_identifier(key, values)} spans {spread}")
        raise SplitNotPossible(
            f"the {kind} split holds out whole values of {axis!r}, but the group key "
            f"{key} produces group(s) containing more than one: "
            + "; ".join(offenders)
            + f". Assigning such a group puts one {axis} on both sides of the split. "
            f"Add {extra} to the group key, or split on an axis the key already nests."
        )
    table = counts.merge(variants, on=key, how="left", validate="one_to_one")
    table["group_id"] = [
        _identifier(key, row) for row in table[key].to_numpy(dtype=object)
    ]
    table["stratum"] = [
        _identifier(within, row) for row in table[list(within)].to_numpy(dtype=object)
    ] if within else ""
    table["axis_value"] = table[axis]
    ordered = [*key, "group_id", "stratum", "axis_value", "n_rows"]
    return table[ordered + [c for c in table.columns if c not in ordered]]


# ---------------------------------------------------------------------------
# choosing what to hold out
# ---------------------------------------------------------------------------


def _components(label: object) -> tuple[str, ...]:
    """Single agents inside a co-applied treatment label."""
    return tuple(str(label).split(_COMBINATION_SEPARATOR))


def _is_single(value: object) -> bool:
    return _COMBINATION_SEPARATOR not in str(value)


def _also_leaks(value: object, values: Sequence) -> tuple:
    """Other labels that would put ``value`` back in training if they stayed there.

    Holding out H2O2 while ``DTT+H2O2`` remains in training does not hold out H2O2: the
    model has seen it, at a dose, alongside another agent. A leave-one-stressor-out score
    computed that way is not a transfer score, and the leak is invisible in the counts
    because the held-out label really was withheld. Every treatment containing the agent
    goes with it.
    """
    return tuple(v for v in values if v != value and str(value) in _components(v))


def _candidates(spec: SplitKind, values: list, where: str) -> tuple[list, list, str | None]:
    """Values eligible to be held out, values that count towards the training minimum.

    Returns:
        ``(candidates, units, reason)``. ``reason`` is set when the rule has nothing to
        choose from at all, which is a different failure from having too few of them.
    """
    if spec.rule == "top":
        return [values[-1]], values, None
    if spec.rule == "interior":
        return values[1:-1], values, None
    if spec.rule == "single":
        singles = [v for v in values if _is_single(v)]
        if not singles:
            return [], singles, (
                f"every treatment{where} is a co-applied pair; there is no single agent "
                f"to hold out"
            )
        return singles, singles, None
    if spec.rule == "combination":
        combinations = [v for v in values if not _is_single(v)]
        if not combinations:
            return [], values, (
                f"no co-applied treatment{where}: every label is a single agent, so "
                f"there is no interaction to hold out"
            )
        return combinations, values, None
    return list(values), values, None


def _choose(
    spec: SplitKind,
    stratum: str,
    values: list,
    seed: int,
    validation_units: int,
    requested: object | None,
) -> tuple[tuple, tuple, str | None]:
    """Held-out values, validation values, and a refusal reason if there are none.

    Returns:
        ``(test, validation, reason)``. ``reason`` is ``None`` when the split stands;
        otherwise it says which count fell short, and ``test`` and ``validation`` are
        empty.
    """
    rng = _stratum_rng(seed, spec.name, stratum)
    where = f" for {stratum}" if stratum else ""

    candidates, units, reason = _candidates(spec, values, where)
    if reason is not None:
        return (), (), reason
    if len(units) < spec.min_units:
        return (), (), (
            f"{len(units)} distinct {spec.axis} value(s){where}; the {spec.name} split "
            f"needs at least {spec.min_units}"
        )
    if not candidates:
        return (), (), (
            f"{len(units)} distinct {spec.axis} value(s){where} but none of them is one "
            f"the {spec.name} split may hold out"
        )

    if requested is not None:
        if requested not in values:
            return (), (), (
                f"{spec.axis} {_format(requested)} is not present{where}; have "
                f"{[_format(v) for v in values]}"
            )
        if requested not in candidates:
            return (), (), (
                f"{spec.axis} {_format(requested)} cannot be the held-out value of the "
                f"{spec.name} split{where}; that split holds out "
                + ("the top rung, which is " if spec.rule == "top" else "one of ")
                + f"{[_format(v) for v in candidates]}"
            )
        primary = requested
    else:
        primary = candidates[int(rng.integers(len(candidates)))]
    test = (primary, *_also_leaks(primary, values)) if spec.rule == "single" else (primary,)

    remaining = [v for v in values if v not in test]
    if validation_units:
        if spec.rule == "top":
            # The nearest thing to a second extrapolation fold: the next rungs down,
            # in order, so tuning happens on the same kind of question as scoring.
            source = [v for v in reversed(remaining)]
        else:
            eligible = [v for v in candidates if v not in test]
            source = [eligible[i] for i in rng.permutation(len(eligible))]
        if len(source) < validation_units:
            return (), (), (
                f"{validation_units} validation unit(s) requested{where} but only "
                f"{len(source)} eligible {spec.axis} value(s) remain after the held-out "
                f"one"
            )
        picked = source[:validation_units]
        validation = tuple(picked)
        if spec.rule == "single":
            validation = tuple(dict.fromkeys(
                [v for p in picked for v in (p, *_also_leaks(p, remaining))]
            ))
    else:
        validation = ()

    trained_on = [v for v in units if v not in test and v not in validation]
    if len(trained_on) < spec.min_train_units:
        return (), (), (
            f"{len(units)} {spec.axis} value(s){where} leaves {len(trained_on)} in "
            f"training after holding out {len(test)} and reserving {len(validation)}; the "
            f"{spec.name} split needs at least {spec.min_train_units}"
        )

    if spec.rule == "combination":
        singles = {str(v) for v in trained_on if _is_single(v)}
        for held in test:
            absent = [c for c in _components(held) if c not in singles]
            if absent:
                return (), (), (
                    f"holding out {held!r} would also remove component(s) {absent} from "
                    f"the training set, which makes this an extrapolation to unseen "
                    f"agents rather than a held-out combination"
                )

    return test, validation, None


def _confounding_reason(frame: pd.DataFrame, spec: SplitKind) -> str | None:
    """Why this axis cannot be held out on this dataset, if it cannot.

    Only one case exists so far and it is the important one. On the real biosensor plates
    every construct sits under exactly one stressor, so removing a stressor removes its
    constructs with it and nothing in the resulting score separates "unseen agent" from
    "unseen reporter". Detected from the frame's own crossing rather than from the names
    in it, so it fires for any dataset built the same way.
    """
    if spec.axis != STRESSOR_COLUMN or spec.rule == "combination":
        return None
    if CONSTRUCT_COLUMN not in frame.columns:
        return None
    crossing = frame[[CONSTRUCT_COLUMN, STRESSOR_COLUMN]].drop_duplicates()
    per_construct = crossing.groupby(CONSTRUCT_COLUMN, sort=True).size()
    if (per_construct > 1).any():
        return None
    return (
        f"every construct in this frame appears under exactly one stressor "
        f"({dict(zip(crossing[CONSTRUCT_COLUMN], crossing[STRESSOR_COLUMN]))}), so "
        f"holding out a stressor also holds out its constructs entirely and the score "
        f"cannot distinguish an unseen agent from an unseen reporter; use "
        f"heldout_construct within a stressor instead"
    )


def _plan(
    frame: pd.DataFrame,
    kind: str,
    group_key: Sequence[str],
    seed: int,
    within: Sequence[str] | None,
    validation_units: int,
    hold_out: object | Mapping[str, object] | None,
) -> tuple[SplitKind, tuple[str, ...], pd.DataFrame, dict, list[tuple[str, str]]]:
    """Everything both :func:`make_split` and :func:`feasibility` need.

    Structural problems -- a missing column, a group key that straddles the axis, a frame
    that disagrees with the declared construct-to-stressor mapping -- raise here, because
    they are mistakes in the call rather than facts about the data. Shortfalls in the data
    come back as reasons, one per stratum, so both callers report the same thing.
    """
    if kind not in SPLIT_KINDS:
        raise KeyError(f"no split kind {kind!r}; have {sorted(SPLIT_KINDS)}")
    spec = SPLIT_KINDS[kind]
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise SplitNotPossible(f"the {kind} split needs a non-empty frame")
    group_key = tuple(dict.fromkeys(group_key))
    if not group_key:
        raise SplitNotPossible(
            "the group key must be given explicitly; it is the declaration of what an "
            "indivisible unit is, and a split with an implicit one cannot be audited"
        )
    strata = tuple(spec.within if within is None else tuple(within))
    if validation_units < 0:
        raise ValueError("validation_units must not be negative")

    _require_columns(frame, [spec.axis, *strata, *group_key], kind)
    _check_construct_stressor_agreement(frame)

    if spec.axis == CONSTRUCT_COLUMN:
        _refuse_pooled_constructs(frame, strata)

    table = _group_table(frame, group_key, spec.axis, strata, kind)

    confounded = _confounding_reason(frame, spec)
    if confounded is not None:
        stratum_ids = sorted(set(table["stratum"]))
        return spec, strata, table, {}, [(s, confounded) for s in stratum_ids]

    plan: dict[str, tuple[tuple, tuple]] = {}
    failures: list[tuple[str, str]] = []
    for stratum in sorted(set(table["stratum"])):
        block = table[table["stratum"] == stratum]
        values = sorted(set(block["axis_value"]))
        requested = hold_out
        if isinstance(hold_out, Mapping):
            requested = hold_out.get(stratum)
        test, validation, reason = _choose(
            spec, stratum, values, seed, validation_units, requested
        )
        if reason is not None:
            failures.append((stratum, reason))
        else:
            plan[stratum] = (test, validation)
    return spec, strata, table, plan, failures


def _refuse_pooled_constructs(frame: pd.DataFrame, strata: Sequence[str]) -> None:
    """A construct-level split must sit inside one stressor, and this proves it does.

    Delegates to :func:`qpcr.assert_single_stressor` whenever the constructs are ones the
    package knows, so the refusal and its wording come from the one place that records
    which agent challenged which construct. Constructs it does not know -- a simulated
    panel's reporters, say -- fall back to the frame's own stressor column.
    """
    if STRESSOR_COLUMN not in frame.columns:
        return
    if not strata:
        blocks = [("", frame)]
    else:
        ids = np.array([
            _identifier(strata, row)
            for row in frame[list(strata)].to_numpy(dtype=object)
        ])
        blocks = [(name, frame[ids == name]) for name in sorted(set(ids))]
    for name, block in blocks:
        constructs = sorted(set(block[CONSTRUCT_COLUMN]))
        if all(c in STRESSOR_FOR_CONSTRUCT for c in constructs):
            try:
                assert_single_stressor(constructs)
            except ValueError as exc:
                raise SplitNotPossible(
                    f"the heldout_construct split was asked to pool {constructs} in one "
                    f"stratum: {exc}. Pass within=('{STRESSOR_COLUMN}',)."
                ) from exc
        agents = sorted(set(block[STRESSOR_COLUMN]))
        if len(agents) > 1:
            raise SplitNotPossible(
                f"the heldout_construct split was asked to pool constructs challenged with "
                f"{agents} in one stratum ({name}); their dose axes are not comparable. "
                f"Pass within=('{STRESSOR_COLUMN}',)."
            )


# ---------------------------------------------------------------------------
# the two entry points
# ---------------------------------------------------------------------------


def make_split(
    frame: pd.DataFrame,
    kind: str,
    *,
    group_key: Sequence[str],
    seed: int = 0,
    within: Sequence[str] | None = None,
    validation_units: int = 0,
    hold_out: object | Mapping[str, object] | None = None,
    dataset: str = "",
) -> Split:
    """Partition ``frame``'s groups for one held-out claim, or refuse to.

    Args:
        frame: One row per observation. Must carry the kind's axis column, its stratum
            columns and every column of the group key, with no missing values in any of
            them. Column names are the fixed vocabulary at the top of this module.
        kind: A key of :data:`SPLIT_KINDS`.
        group_key: Columns whose distinct combinations are indivisible. This is the
            anti-leak declaration and it is required rather than defaulted: for a
            replicate-level claim it must nest the plate, and for a dose-level claim it
            must nest construct and dose.
        seed: Draws the fold. Recorded in the manifest. Ignored by ``extrapolation``,
            whose held-out value is fixed by the data.
        within: Stratum columns, defaulting to the kind's own. Each stratum is split
            independently, so each construct's ladder gives up its own top dose rather
            than the ladder with the largest numbers giving up all of them.
        validation_units: Further axis values moved out of training, drawn by the same
            grouping discipline so that model selection never touches the held-out set.
            Zero leaves the split as train and test only.
        hold_out: Pin the fold instead of drawing it -- one value for every stratum, or a
            mapping from stratum id to value. Refused if the value is not one the kind is
            allowed to hold out, which is how an "extrapolation" that actually withheld a
            middle dose gets caught.
        dataset: Name recorded in the manifest, so rows from the real plates and rows
            from the simulated panel stay distinguishable in one table.

    Returns:
        A :class:`Split` whose groups each appear on exactly one side.

    Raises:
        KeyError: if ``kind`` is not a known split kind.
        SplitNotPossible: if a required column is missing or holds missing values; if the
            group key produces groups that straddle the axis; if the frame's labelling
            disagrees with ``qpcr.STRESSOR_FOR_CONSTRUCT``; or if any stratum has too few
            units for the split to mean what it is named -- with the count that says so,
            for every stratum that fell short.
    """
    spec, strata, table, plan, failures = _plan(
        frame, kind, group_key, seed, within, validation_units, hold_out
    )
    if failures:
        detail = "; ".join(reason for _, reason in failures)
        raise SplitNotPossible(
            f"the {kind} split is not available on this dataset: {detail}. "
            f"{_claim_note(spec)}"
        )

    assignments = []
    for stratum, axis_value in zip(table["stratum"], table["axis_value"]):
        test, validation = plan[stratum]
        if axis_value in test:
            assignments.append(TEST)
        elif axis_value in validation:
            assignments.append(VALIDATION)
        else:
            assignments.append(TRAIN)
    groups = table.copy()
    groups["assignment"] = assignments

    split = Split(
        kind=kind,
        seed=int(seed),
        group_key=tuple(dict.fromkeys(group_key)),
        within=tuple(strata),
        dataset=dataset,
        groups=groups,
    )
    _assert_disjoint(split)
    return split


def _claim_note(spec: SplitKind) -> str:
    return f"That split would have licensed: {spec.claim}."


def _assert_disjoint(split: Split) -> None:
    """No group on two sides, and no axis value on two sides within a stratum.

    Cheap, and it runs on every split rather than only in the tests. The first condition
    is guaranteed by construction -- one assignment per group row -- and is checked anyway
    because "guaranteed by construction" is what every leaked split was before it leaked.
    """
    counts = split.groups.groupby("group_id", sort=False)["assignment"].nunique()
    straddling = sorted(counts[counts > 1].index)
    if straddling:
        raise SplitNotPossible(
            f"group(s) {straddling[:3]} were assigned to more than one side"
        )
    per_value = split.groups.groupby(["stratum", "axis_value"], sort=False)["assignment"].nunique()
    leaked = per_value[per_value > 1]
    if len(leaked):
        raise SplitNotPossible(
            f"{split.spec.axis} value(s) {sorted(set(leaked.index))[:3]} appear on more "
            f"than one side of the split, which is the leak this module exists to prevent"
        )


def feasibility(
    frame: pd.DataFrame,
    kind: str,
    *,
    group_key: Sequence[str],
    seed: int = 0,
    within: Sequence[str] | None = None,
    validation_units: int = 0,
) -> pd.DataFrame:
    """Which strata can support this split today, and why the others cannot.

    The same checks :func:`make_split` runs, reported rather than raised. This is what a
    documentation table is built from -- "n=2 blocks the replicate split for three of four
    constructs" is a claim that should come out of the data on every run, not be typed
    into a document once.

    Args:
        frame: As :func:`make_split`.
        kind: As :func:`make_split`.
        group_key: As :func:`make_split`.
        seed: As :func:`make_split`.
        within: As :func:`make_split`.
        validation_units: As :func:`make_split`.

    Returns:
        One row per stratum with ``stratum``, ``n_units``, ``n_groups``, ``n_rows``,
        ``feasible`` and ``reason``. ``reason`` is empty when the stratum is feasible.

    Raises:
        KeyError: if ``kind`` is not a known split kind.
        SplitNotPossible: for structural problems -- a missing column, a leaking group
            key -- which are mistakes in the call rather than limits of the data.
    """
    spec, _, table, plan, failures = _plan(
        frame, kind, group_key, seed, within, validation_units, None
    )
    reasons = dict(failures)
    rows = []
    for stratum in sorted(set(table["stratum"])):
        block = table[table["stratum"] == stratum]
        # The count the rule actually works with, not the count of distinct labels: a
        # leave-one-stressor-out over 25 agents plus one co-applied pair has 25 units,
        # and reporting 26 would overstate what is there to hold out.
        _, units, _ = _candidates(spec, sorted(set(block["axis_value"])), "")
        rows.append({
            "split_kind": kind,
            "stratum": stratum,
            "axis": spec.axis,
            "n_units": int(len(units)),
            "n_groups": int(len(block)),
            "n_rows": int(block["n_rows"].sum()),
            "feasible": stratum in plan,
            "reason": reasons.get(stratum, ""),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# the manifest
# ---------------------------------------------------------------------------


def split_manifest(splits: Sequence[Split]) -> pd.DataFrame:
    """Every split's groups in one table, which is the artefact that gets committed.

    Args:
        splits: Splits to record. They need not share a group key; columns one split
            does not use are left empty for its rows.

    Returns:
        One row per group per split, sorted so that a diff between two runs is legible.
    """
    if not splits:
        raise ValueError("no splits to record; an empty manifest reads as 'no held-out data'")
    frames = [split.manifest() for split in splits]
    manifest = pd.concat(frames, ignore_index=True, sort=False)
    manifest = manifest.sort_values(
        ["dataset", "split_kind", "seed", "stratum", "assignment", "group_id"],
        kind="stable",
    ).reset_index(drop=True)
    # Fixed column order, so that two manifests diff line by line rather than by whose
    # first split happened to name its group key first.
    preferred = [
        "dataset", "split_kind", "seed", "group_key", "within", "axis", "stratum",
        "group_id", CONSTRUCT_COLUMN, STRESSOR_COLUMN, PLATE_COLUMN, DOSE_COLUMN,
        "axis_value", "assignment", "n_rows", "held_out", "partition_hash",
    ]
    order = [c for c in preferred if c in manifest.columns]
    return manifest[order + [c for c in manifest.columns if c not in order]]


def summarise_manifest(manifest: pd.DataFrame) -> pd.DataFrame:
    """Counts per split kind and side, with the hash that names the partition.

    Args:
        manifest: The output of :func:`split_manifest`.

    Returns:
        One row per dataset, kind, seed and assignment.
    """
    keys = ["dataset", "split_kind", "seed", "assignment"]
    grouped = manifest.groupby(keys, sort=True)
    summary = grouped.agg(
        n_groups=("group_id", "size"),
        n_rows=("n_rows", "sum"),
        n_strata=("stratum", "nunique"),
    ).reset_index()
    carried = manifest.drop_duplicates(subset=["dataset", "split_kind", "seed"])[
        ["dataset", "split_kind", "seed", "group_key", "within", "axis",
         "partition_hash"]
    ]
    return summary.merge(carried, on=["dataset", "split_kind", "seed"], how="left")
