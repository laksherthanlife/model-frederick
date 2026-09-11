"""What a split has to guarantee before a number scored on it means anything.

The load-bearing tests here are the refusals. A split module that only ever returns a
partition is easy to write and useless: the mistakes that matter -- a plate's wells on
both sides, a "held-out top dose" that was really a middle one, a leave-one-stressor-out
that left the held-out agent in training inside a combination -- all produce a
perfectly-shaped partition and a wrong number. Each of those is pinned here, alongside the
one refusal the real dataset actually triggers: at two biological replicates there is no
replicate-level split to be had, for three of the four constructs.

Everything is synthetic and in-memory. The suite has to pass on a machine that has never
seen a plate reader.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ystwin.analysis.splits import (
    SPLIT_KINDS,
    TEST,
    TRAIN,
    VALIDATION,
    Split,
    SplitNotPossible,
    feasibility,
    make_split,
    partition_hash,
    split_manifest,
    summarise_manifest,
)
from ystwin.analysis.uncertainty import MIN_PLATES_FOR_INTERVAL

DTT_DOSES = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
H2O2_DOSES = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0)
REAL_CONSTRUCTS = {
    "UPRE1": ("DTT", DTT_DOSES),
    "UPRE2": ("DTT", DTT_DOSES),
    "NativeYap1": ("H2O2", H2O2_DOSES),
    "AlteredYap1": ("H2O2", H2O2_DOSES),
}

DOSE_KEY = ("plate", "construct", "dose_mM")


def biosensor_frame(plates_per_construct=None, replicates=3):
    """A well-level frame shaped like the NewProtocol plates.

    Args:
        plates_per_construct: Construct to the plates it was read on. Defaults to the
            coverage the real exports have: four for UPRE1, two for the rest.
        replicates: Technical replicate wells per construct, dose and plate.
    """
    if plates_per_construct is None:
        plates_per_construct = {
            "UPRE1": ("P1", "P2", "P3", "P4"),
            "UPRE2": ("P1", "P3"),
            "NativeYap1": ("P1", "P3"),
            "AlteredYap1": ("P1", "P3"),
        }
    rows = []
    for construct, plates in plates_per_construct.items():
        stressor, doses = REAL_CONSTRUCTS[construct]
        for plate in plates:
            for dose in doses:
                for well in range(replicates):
                    rows.append({
                        "plate": plate,
                        "construct": construct,
                        "stressor": stressor,
                        "dose_mM": float(dose),
                        "well": f"{plate}-{construct}-{dose}-{well}",
                    })
    return pd.DataFrame(rows)


def panel_frame(stressors=("DTT", "H2O2", "NaCl", "heat"), doses=(1.0, 2.0, 4.0, 8.0, 16.0),
                combinations=("DTT+H2O2",), replicates=3):
    """A treatment-level frame shaped like the simulated panel: no plate, no construct."""
    rows = []
    for label in (*stressors, *combinations):
        for dose in doses:
            for well in range(replicates):
                rows.append({
                    "stressor": label,
                    "dose_mM": float(dose),
                    "well": f"{label}-{dose}-{well}",
                })
    return pd.DataFrame(rows)


def sides(split: Split) -> dict[str, set[str]]:
    """Group ids on each side, which is what every leak test looks at."""
    return {
        side: set(split.groups.loc[split.groups.assignment == side, "group_id"])
        for side in (TRAIN, VALIDATION, TEST)
    }


# ---------------------------------------------------------------------------
# no group on both sides, and nothing lost on the way
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind,frame,key,extra", [
    ("interpolation", biosensor_frame(), DOSE_KEY, {}),
    ("extrapolation", biosensor_frame(), DOSE_KEY, {}),
    ("heldout_construct", biosensor_frame(), DOSE_KEY, {}),
    ("interpolation", panel_frame(), ("stressor", "dose_mM"), {"within": ("stressor",)}),
    ("extrapolation", panel_frame(), ("stressor", "dose_mM"), {"within": ("stressor",)}),
    ("heldout_stressor", panel_frame(), ("stressor",), {}),
    ("heldout_combination", panel_frame(), ("stressor",), {}),
])
def test_no_group_appears_on_both_sides(kind, frame, key, extra):
    split = make_split(frame, kind, group_key=key, seed=1, **extra)
    by_side = sides(split)
    assert by_side[TRAIN] and by_side[TEST]
    assert not by_side[TRAIN] & by_side[TEST]
    assert not by_side[TRAIN] & by_side[VALIDATION]
    assert not by_side[TEST] & by_side[VALIDATION]


@pytest.mark.parametrize("kind,frame,key,extra", [
    ("interpolation", biosensor_frame(), DOSE_KEY, {}),
    ("extrapolation", biosensor_frame(), DOSE_KEY, {}),
    ("heldout_construct", biosensor_frame(), DOSE_KEY, {}),
    ("heldout_stressor", panel_frame(), ("stressor",), {}),
])
def test_counts_sum_to_the_input(kind, frame, key, extra):
    """A row that fell out of the partition is a row nobody notices is missing."""
    split = make_split(frame, kind, group_key=key, seed=1, **extra)
    counts = split.counts()
    assert sum(side["rows"] for side in counts.values()) == len(frame)
    assert sum(side["groups"] for side in counts.values()) == len(split.groups)


def test_every_row_is_assigned_exactly_once():
    frame = biosensor_frame()
    split = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=1)
    assigned = split.assign_rows(frame)
    assert len(assigned) == len(frame)
    assert set(assigned) <= {TRAIN, TEST}
    # A dose either is the held-out one for its construct or it is not; no well of a
    # held-out condition survives into training.
    for construct, held in split.held_out().items():
        name = construct.split("=", 1)[1]
        mask = (frame.construct == name) & frame.dose_mM.isin(held)
        assert set(assigned[mask]) == {TEST}
        assert TEST not in set(assigned[(frame.construct == name) & ~mask])


# ---------------------------------------------------------------------------
# the same seed gives the same partition, and the hash proves it
# ---------------------------------------------------------------------------


def test_same_seed_gives_the_same_partition_and_hash():
    frame = biosensor_frame()
    first = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=7)
    second = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=7)
    assert first.assignment == second.assignment
    assert first.hash == second.hash


def test_a_different_seed_gives_a_different_partition():
    frame = biosensor_frame()
    seeds = {make_split(frame, "interpolation", group_key=DOSE_KEY, seed=s).hash
             for s in range(12)}
    assert len(seeds) > 1


def test_the_hash_does_not_depend_on_row_order():
    """A digest that moved with row order would prove nothing about the partition."""
    frame = biosensor_frame()
    shuffled = frame.sample(frac=1.0, random_state=3).reset_index(drop=True)
    a = make_split(frame, "heldout_construct", group_key=DOSE_KEY, seed=2)
    b = make_split(shuffled, "heldout_construct", group_key=DOSE_KEY, seed=2)
    assert a.hash == b.hash


def test_the_hash_moves_when_one_group_changes_side():
    frame = biosensor_frame()
    split = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=1)
    tampered = dict(split.assignment)
    victim = next(g for g, side in tampered.items() if side == TRAIN)
    tampered[victim] = TEST
    assert partition_hash(tampered) != split.hash


def test_hashing_is_insensitive_to_mapping_insertion_order():
    forward = {"a": TRAIN, "b": TEST, "c": TRAIN}
    backward = {"c": TRAIN, "b": TEST, "a": TRAIN}
    assert partition_hash(forward) == partition_hash(backward)


def test_one_stratums_fold_does_not_move_when_another_construct_is_added():
    """Adding a construct must not silently re-draw the folds already published."""
    two = biosensor_frame({"UPRE1": ("P1", "P3"), "UPRE2": ("P1", "P3")})
    three = biosensor_frame(
        {"UPRE1": ("P1", "P3"), "UPRE2": ("P1", "P3"), "NativeYap1": ("P1", "P3")}
    )
    small = make_split(two, "interpolation", group_key=DOSE_KEY, seed=5)
    large = make_split(three, "interpolation", group_key=DOSE_KEY, seed=5)
    assert small.held_out()["construct=UPRE1"] == large.held_out()["construct=UPRE1"]


# ---------------------------------------------------------------------------
# extrapolation is the top dose, and is not interpolation
# ---------------------------------------------------------------------------


def test_extrapolation_holds_out_the_top_dose_of_each_ladder():
    """Per construct, because DTT tops out at 5 mM and H2O2 at 4."""
    frame = biosensor_frame()
    split = make_split(frame, "extrapolation", group_key=DOSE_KEY, seed=0)
    assert split.held_out() == {
        "construct=AlteredYap1": (4.0,),
        "construct=NativeYap1": (4.0,),
        "construct=UPRE1": (5.0,),
        "construct=UPRE2": (5.0,),
    }


def test_extrapolation_is_not_moved_by_the_seed():
    frame = biosensor_frame()
    hashes = {make_split(frame, "extrapolation", group_key=DOSE_KEY, seed=s).hash
              for s in range(8)}
    assert len(hashes) == 1
    assert SPLIT_KINDS["extrapolation"].deterministic


def test_extrapolation_refuses_to_hold_out_a_middle_dose():
    """The error this module exists to prevent, asked for explicitly."""
    frame = biosensor_frame()
    with pytest.raises(SplitNotPossible, match="cannot be the held-out value"):
        make_split(frame, "extrapolation", group_key=DOSE_KEY, hold_out=0.5)


def test_interpolation_and_extrapolation_are_scored_on_different_partitions():
    frame = biosensor_frame()
    interior = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=0)
    top = make_split(frame, "extrapolation", group_key=DOSE_KEY, seed=0)
    assert interior.hash != top.hash


def test_interpolation_never_holds_out_the_lowest_or_the_top_rung():
    """The lowest rung is the 0 mM control, and below the training range is not between it."""
    frame = biosensor_frame()
    for seed in range(20):
        split = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=seed)
        for stratum, held in split.held_out().items():
            name = stratum.split("=", 1)[1]
            ladder = sorted(REAL_CONSTRUCTS[name][1])
            for dose in held:
                assert ladder[0] < dose < ladder[-1]


def test_interpolation_refuses_a_three_rung_ladder():
    """Holding out the middle of three leaves a line, which tests nothing about curvature."""
    short = biosensor_frame({"UPRE1": ("P1", "P3")})
    short = short[short.dose_mM.isin([0.0, 1.0, 5.0])]
    with pytest.raises(SplitNotPossible, match="needs at least 4"):
        make_split(short, "interpolation", group_key=DOSE_KEY)


# ---------------------------------------------------------------------------
# grouping is respected, and a leaking key is refused rather than permitted
# ---------------------------------------------------------------------------


def test_a_replicate_split_refuses_a_group_key_that_straddles_plates():
    """Grouping by construct and dose puts every plate's wells in every group."""
    frame = biosensor_frame({"UPRE1": ("P1", "P2", "P3", "P4")})
    with pytest.raises(SplitNotPossible, match="more than one"):
        make_split(frame, "heldout_replicate", group_key=("construct", "dose_mM"))


def test_a_dose_split_refuses_a_group_key_that_straddles_doses():
    frame = biosensor_frame()
    with pytest.raises(SplitNotPossible, match="more than one"):
        make_split(frame, "interpolation", group_key=("plate", "construct"))


def test_no_plate_has_wells_on_both_sides_of_a_replicate_split():
    frame = biosensor_frame({"UPRE1": ("P1", "P2", "P3", "P4")})
    split = make_split(frame, "heldout_replicate", group_key=DOSE_KEY, seed=4)
    assigned = split.assign_rows(frame)
    per_plate = frame.assign(side=assigned).groupby("plate").side.nunique()
    assert set(per_plate) == {1}


def test_no_construct_and_dose_pair_straddles_a_dose_split():
    frame = biosensor_frame()
    split = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=4)
    assigned = split.assign_rows(frame)
    condition = frame.assign(side=assigned).groupby(["construct", "dose_mM"]).side.nunique()
    assert set(condition) == {1}


def test_a_manifest_cannot_be_applied_to_a_frame_it_does_not_describe():
    frame = biosensor_frame()
    split = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=0)
    other = biosensor_frame({"UPRE1": ("P9",)})
    with pytest.raises(SplitNotPossible, match="different frame"):
        split.assign_rows(other)


# ---------------------------------------------------------------------------
# impossible splits refuse loudly, with the count that says so
# ---------------------------------------------------------------------------


def test_two_biological_replicates_is_not_a_replicate_split():
    """The refusal the real dataset triggers, for three of its four constructs."""
    frame = biosensor_frame()
    with pytest.raises(SplitNotPossible) as raised:
        make_split(frame, "heldout_replicate", group_key=DOSE_KEY)
    message = str(raised.value)
    for construct in ("UPRE2", "NativeYap1", "AlteredYap1"):
        assert construct in message
    assert "2 distinct plate value(s)" in message
    assert "UPRE1" not in message


def test_the_replicate_minimum_is_the_packages_own_interval_floor():
    """Not a number invented here: below it, uncertainty.py refuses to state anything."""
    spec = SPLIT_KINDS["heldout_replicate"]
    assert spec.min_train_units == MIN_PLATES_FOR_INTERVAL
    assert spec.min_units == MIN_PLATES_FOR_INTERVAL + 1


def test_feasibility_reports_the_same_refusal_without_raising():
    frame = biosensor_frame()
    table = feasibility(frame, "heldout_replicate", group_key=DOSE_KEY)
    verdicts = dict(zip(table.stratum, table.feasible))
    assert verdicts["construct=UPRE1"]
    assert not verdicts["construct=UPRE2"]
    assert not verdicts["construct=NativeYap1"]
    assert not verdicts["construct=AlteredYap1"]
    assert dict(zip(table.stratum, table.n_units))["construct=UPRE1"] == 4


def test_the_feasible_construct_alone_does_give_a_replicate_split():
    frame = biosensor_frame()
    split = make_split(frame[frame.construct == "UPRE1"], "heldout_replicate",
                       group_key=DOSE_KEY, seed=0)
    assert list(split.held_out()) == ["construct=UPRE1"]
    assert split.counts()[TEST]["groups"] == len(DTT_DOSES)


def test_a_missing_axis_column_names_the_axis_rather_than_the_column():
    """A simulated panel has no plate. That is a fact about the dataset, not a typo."""
    with pytest.raises(SplitNotPossible, match="plate"):
        make_split(panel_frame(), "heldout_replicate", group_key=("stressor", "dose_mM"))


def test_a_missing_value_is_refused_rather_than_grouped():
    frame = biosensor_frame({"UPRE1": ("P1", "P2", "P3", "P4")})
    frame.loc[0, "plate"] = None
    with pytest.raises(SplitNotPossible, match="missing values"):
        make_split(frame, "heldout_replicate", group_key=DOSE_KEY)


def test_an_unknown_kind_lists_the_known_ones():
    with pytest.raises(KeyError, match="heldout_replicate"):
        make_split(biosensor_frame(), "heldout_everything", group_key=DOSE_KEY)


# ---------------------------------------------------------------------------
# the two stressors are not one axis
# ---------------------------------------------------------------------------


def test_a_construct_split_refuses_to_pool_the_two_stressors():
    """Delegates to qpcr.assert_single_stressor, whose wording comes through intact."""
    frame = biosensor_frame()
    with pytest.raises(SplitNotPossible, match="dose axes are not comparable"):
        make_split(frame, "heldout_construct", group_key=DOSE_KEY, within=())


def test_a_construct_split_within_a_stressor_holds_out_one_per_stressor():
    frame = biosensor_frame()
    split = make_split(frame, "heldout_construct", group_key=DOSE_KEY, seed=0)
    held = split.held_out()
    assert set(held) == {"stressor=DTT", "stressor=H2O2"}
    assert held["stressor=DTT"][0] in ("UPRE1", "UPRE2")
    assert held["stressor=H2O2"][0] in ("NativeYap1", "AlteredYap1")


def test_a_frame_that_mislabels_a_constructs_stressor_is_refused():
    frame = biosensor_frame()
    frame.loc[frame.construct == "UPRE1", "stressor"] = "H2O2"
    with pytest.raises(SplitNotPossible, match="STRESSOR_FOR_CONSTRUCT"):
        make_split(frame, "heldout_construct", group_key=DOSE_KEY)


def test_a_stressor_split_is_refused_when_constructs_are_nested_in_stressors():
    """Holding out H2O2 there also holds out both Yap1 reporters, so it proves nothing."""
    frame = biosensor_frame()
    with pytest.raises(SplitNotPossible, match="cannot distinguish an unseen agent"):
        make_split(frame, "heldout_stressor",
                   group_key=("plate", "construct", "stressor", "dose_mM"))


def test_two_stressors_are_not_enough_for_a_stressor_split():
    frame = panel_frame(stressors=("DTT", "H2O2"), combinations=())
    with pytest.raises(SplitNotPossible, match="needs at least 3"):
        make_split(frame, "heldout_stressor", group_key=("stressor",))


# ---------------------------------------------------------------------------
# combinations
# ---------------------------------------------------------------------------


def test_holding_out_a_stressor_takes_its_combinations_with_it():
    """Leaving DTT+H2O2 in training while withholding H2O2 withholds a label, not an agent."""
    frame = panel_frame()
    split = make_split(frame, "heldout_stressor", group_key=("stressor",),
                       seed=0, hold_out="H2O2")
    held = split.held_out()[""]
    assert "H2O2" in held
    assert "DTT+H2O2" in held
    trained = {g.split("=", 1)[1] for g in sides(split)[TRAIN]}
    assert not any("H2O2" in label for label in trained)


def test_a_combination_split_keeps_both_components_in_training():
    frame = panel_frame()
    split = make_split(frame, "heldout_combination", group_key=("stressor",), seed=0)
    assert split.held_out()[""] == ("DTT+H2O2",)
    trained = {g.split("=", 1)[1] for g in sides(split)[TRAIN]}
    assert {"DTT", "H2O2"} <= trained


def test_a_combination_split_refuses_when_a_component_was_never_run_alone():
    frame = panel_frame(stressors=("DTT", "NaCl", "heat"), combinations=("DTT+H2O2",))
    with pytest.raises(SplitNotPossible, match="extrapolation to unseen"):
        make_split(frame, "heldout_combination", group_key=("stressor",))


def test_a_combination_split_refuses_a_panel_with_no_combinations():
    frame = panel_frame(combinations=())
    with pytest.raises(SplitNotPossible, match="no interaction to hold out"):
        make_split(frame, "heldout_combination", group_key=("stressor",))


# ---------------------------------------------------------------------------
# validation is carved by group too
# ---------------------------------------------------------------------------


def test_validation_is_carved_out_of_training_and_never_out_of_test():
    frame = panel_frame(stressors=("DTT", "H2O2", "NaCl", "heat", "sorbitol"),
                        combinations=())
    split = make_split(frame, "heldout_stressor", group_key=("stressor",),
                       seed=0, validation_units=1)
    counts = split.counts()
    assert counts[VALIDATION]["groups"] == 1
    assert counts[TEST]["groups"] == 1
    assert sum(side["rows"] for side in counts.values()) == len(frame)


def test_validation_for_an_extrapolation_split_is_the_next_rung_down():
    """Tuning on the same question as scoring: the second-highest dose, not a random one."""
    frame = biosensor_frame({"UPRE1": ("P1",)})
    split = make_split(frame, "extrapolation", group_key=DOSE_KEY,
                       seed=0, validation_units=1)
    assert split.held_out()["construct=UPRE1"] == (5.0,)
    assert split.validation_values()["construct=UPRE1"] == (2.0,)


def test_too_much_validation_is_refused_rather_than_shrunk():
    frame = biosensor_frame({"UPRE1": ("P1", "P2", "P3", "P4")})
    with pytest.raises(SplitNotPossible, match="needs at least"):
        make_split(frame, "heldout_replicate", group_key=DOSE_KEY, validation_units=2)


# ---------------------------------------------------------------------------
# the manifest
# ---------------------------------------------------------------------------


def test_the_manifest_carries_its_own_provenance_on_every_row():
    frame = biosensor_frame()
    split = make_split(frame, "interpolation", group_key=DOSE_KEY, seed=3,
                       dataset="real_biosensor")
    manifest = split.manifest()
    assert len(manifest) == len(split.groups)
    for column in ("dataset", "split_kind", "seed", "group_key", "within", "axis",
                   "held_out", "partition_hash", "assignment", "n_rows"):
        assert column in manifest.columns
        assert manifest[column].nunique() >= 1
    assert set(manifest.partition_hash) == {split.hash}
    assert set(manifest.seed) == {3}
    assert set(manifest.group_key) == {"plate|construct|dose_mM"}


def test_the_manifest_records_each_rows_own_stratum():
    """A UPRE1 row says 5 mM was withheld; an AlteredYap1 row says 4."""
    frame = biosensor_frame()
    manifest = make_split(frame, "extrapolation", group_key=DOSE_KEY, seed=0).manifest()
    by_construct = manifest.groupby("construct").held_out.unique()
    assert list(by_construct["UPRE1"]) == ["5"]
    assert list(by_construct["AlteredYap1"]) == ["4"]


def test_the_manifest_summary_counts_match_the_splits():
    frame = biosensor_frame()
    splits = [
        make_split(frame, kind, group_key=DOSE_KEY, seed=0, dataset="real_biosensor")
        for kind in ("interpolation", "extrapolation", "heldout_construct")
    ]
    manifest = split_manifest(splits)
    summary = summarise_manifest(manifest)
    assert len(manifest) == sum(len(s.groups) for s in splits)
    for split in splits:
        block = summary[summary.split_kind == split.kind]
        assert int(block.n_rows.sum()) == len(frame)
        assert set(block.partition_hash) == {split.hash}


def test_an_empty_manifest_is_refused():
    """An empty table reads as 'nothing was held out', which is the opposite of the truth."""
    with pytest.raises(ValueError, match="no splits to record"):
        split_manifest([])


def test_every_declared_kind_states_the_claim_it_licenses():
    for name, spec in SPLIT_KINDS.items():
        assert spec.name == name
        assert spec.claim and spec.claim[0].islower()
        assert spec.min_train_units >= 1
        assert spec.min_units > spec.min_train_units or spec.rule == "combination"
