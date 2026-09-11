"""A second read of the same fluorophore is a different measurement, not a duplicate.

Real exports label repeat reads ``mCitrine:480,530[2]`` -- typically a second gain.
Collapsing them silently would let a gain change masquerade as biology.
"""

import numpy as np
import pytest

from ystwin.plate.synergy import read_synergy_kinetic


def test_repeat_reads_of_one_fluorophore_are_kept_as_distinct_channels(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file(second_reporter_read=True))

    assert set(run.channel_names) == {"OD600", "mCitrine[1]", "mCitrine[2]"}


def test_the_two_reads_retain_their_own_values(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file(second_reporter_read=True))

    assert run.channel("mCitrine[1]").data["A1"].to_numpy() == pytest.approx([1000, 1200, 2000])
    assert run.channel("mCitrine[2]").data["A1"].to_numpy() == pytest.approx([80, 96, 160])


def test_asking_for_the_bare_fluorophore_name_refuses_to_guess_between_reads(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file(second_reporter_read=True))

    with pytest.raises(KeyError, match=r"mCitrine\[1\].*mCitrine\[2\]"):
        run.channel("mCitrine")


def test_optics_are_preserved_for_calibration_provenance(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file(second_reporter_read=True))

    assert run.channel("mCitrine[1]").optics == "480,530"
    assert run.channel("OD600").optics == "600"


def test_alignment_keeps_both_reads_separate(make_synergy_file):
    run = read_synergy_kinetic(make_synergy_file(second_reporter_read=True))

    aligned = run.aligned(reference="OD600")
    assert {"OD600", "mCitrine[1]", "mCitrine[2]"} == set(aligned.columns.get_level_values("channel"))
    assert not np.allclose(aligned[("mCitrine[1]", "A1")], aligned[("mCitrine[2]", "A1")])
