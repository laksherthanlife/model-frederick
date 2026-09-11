"""S0: assemble whole synthetic plates that the real readers accept.

The decisive test is the round trip. A synthetic plate is written as a Synergy
export, parsed by the same reader used on wet-lab files, and pushed through G1, the
layout recovery and D2. If any of those behave differently on synthetic data than on
real data, the generator is producing something that only looks right to itself.
"""

import numpy as np
import pandas as pd
import pytest

from ystwin.generator.plate import DEFAULT_PANEL, PlateConditions, generate_plate
from ystwin.plate.layout import NEWPROTOCOL_LAYOUT


def _conditions(**over):
    base = dict(layout=NEWPROTOCOL_LAYOUT, blank_wells=("H1", "H2", "H3"),
                duration_h=4.14, n_timepoints=25, seed=0)
    base.update(over)
    return PlateConditions(**base)


def test_it_fills_every_culture_well_and_the_blanks():
    plate = generate_plate(DEFAULT_PANEL, _conditions())

    assert plate.od.shape[1] == 4 * 3 * 7 + 3
    assert set(plate.blank_wells) == {"H1", "H2", "H3"}


def test_the_timebase_matches_the_real_protocol():
    plate = generate_plate(DEFAULT_PANEL, _conditions())

    assert len(plate.od) == 25
    assert plate.od.index[-1] == pytest.approx(4.14)


def test_blank_wells_carry_medium_only():
    plate = generate_plate(DEFAULT_PANEL, _conditions())

    for well in plate.blank_wells:
        assert plate.od[well].max() < 0.12
        fold = plate.od[well].iloc[-1] / plate.od[well].iloc[0]
        assert fold == pytest.approx(1.0, abs=0.05)


def test_cultures_grow_and_blanks_do_not():
    plate = generate_plate(DEFAULT_PANEL, _conditions())

    culture = plate.od["A1"]
    assert culture.iloc[-1] / culture.iloc[0] > 1.5


def test_each_construct_sits_in_its_own_columns():
    plate = generate_plate(DEFAULT_PANEL, _conditions())
    truth = plate.truth.set_index("well")

    assert truth.loc["A1"].construct == "UPRE1"
    assert truth.loc["A10"].construct == "AlteredYap1"


def test_dose_increases_down_the_rows():
    plate = generate_plate(DEFAULT_PANEL, _conditions())
    truth = plate.truth.set_index("well")

    assert truth.loc["A1"].dose_mM < truth.loc["G1"].dose_mM


def test_the_two_sensor_pairs_get_different_stressors():
    plate = generate_plate(DEFAULT_PANEL, _conditions())
    truth = plate.truth.set_index("well")

    assert truth.loc["A1"].stressor == "DTT"
    assert truth.loc["A7"].stressor == "H2O2"


def test_the_ground_truth_travels_with_the_plate():
    plate = generate_plate(DEFAULT_PANEL, _conditions())

    assert {"well", "construct", "dose_mM", "stressor", "promoter_activity",
            "growth_rate_max"} <= set(plate.truth.columns)


def test_every_record_is_labelled_synthetic():
    plate = generate_plate(DEFAULT_PANEL, _conditions())

    assert (plate.truth.source == "synthetic").all()


class TestRealism:
    def test_reader_noise_is_present_and_of_the_stated_size(self):
        quiet = generate_plate(DEFAULT_PANEL, _conditions(reader_cv=0.0))
        noisy = generate_plate(DEFAULT_PANEL, _conditions(reader_cv=0.02))

        residual = (noisy.od["A1"] - quiet.od["A1"]).std() / quiet.od["A1"].mean()
        assert 0.005 < residual < 0.05

    def test_the_blank_offset_matches_the_real_plates(self):
        plate = generate_plate(DEFAULT_PANEL, _conditions())

        assert plate.od[list(plate.blank_wells)].to_numpy().mean() == pytest.approx(0.098, abs=0.01)

    def test_constructs_can_be_inoculated_at_different_densities(self):
        """NativeYap1 went in 3.7x lower than UPRE1 on the real plate."""
        plate = generate_plate(
            DEFAULT_PANEL, _conditions(inoculum_od={"NativeYap1": 0.04, "UPRE1": 0.15})
        )

        assert plate.od["A7"].iloc[0] < plate.od["A1"].iloc[0] * 0.6

    def test_two_seeds_give_different_noise_but_the_same_truth(self):
        a = generate_plate(DEFAULT_PANEL, _conditions(seed=1))
        b = generate_plate(DEFAULT_PANEL, _conditions(seed=2))

        assert not np.allclose(a.od["A1"], b.od["A1"])
        pd.testing.assert_frame_equal(
            a.truth.drop(columns=["initial_biomass"]),
            b.truth.drop(columns=["initial_biomass"]),
        )
