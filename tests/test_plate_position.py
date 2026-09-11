"""Where a well sits on the plate, which the generator did not know.

Edge wells evaporate faster than interior ones. The medium concentrates, so the stressor
concentrates with it, the path length shortens and the culture behaves as though it were
dosed higher than the pipette says. It is the most reproducible artifact in plate work and
the reason plate maps put controls in the interior.

It matters here for a specific reason. Measurement noise is independent per well, so
replicating averages it away; an edge effect is a property of the position, so replicates
sitting on the same edge inherit the same bias and averaging does nothing. A generator with
noise but no position effect therefore says replication buys more than it does.

Two things follow and both are checkable: a well's assignment must be deterministic, so a
design can be criticised for where it puts its controls, and a batch whose treatments are
laid out down a plate must be able to confound treatment with position.
"""

import numpy as np
import pytest

from ystwin.generator.panel_experiment import panel_dataset
from ystwin.generator.plate_layout import edge_wells, well_grid, well_positions


class TestWellGeometry:
    def test_a_96_well_plate_has_96_wells(self):
        assert len(well_grid(rows=8, columns=12)) == 96

    def test_wells_are_named_by_row_letter_and_column_number(self):
        grid = well_grid(rows=8, columns=12)

        assert grid[0] == "A1"
        assert grid[-1] == "H12"

    def test_the_edge_is_the_outer_ring(self):
        edge = edge_wells(rows=8, columns=12)

        assert "A1" in edge and "H12" in edge and "A6" in edge and "D1" in edge
        assert "D6" not in edge

    def test_a_96_well_plate_has_36_edge_wells(self):
        assert len(edge_wells(rows=8, columns=12)) == 2 * 8 + 2 * 12 - 4

    def test_positions_are_deterministic(self):
        assert well_positions(5, rows=8, columns=12) == well_positions(5, rows=8, columns=12)

    def test_it_refuses_more_samples_than_the_plate_holds(self):
        with pytest.raises(ValueError, match="only 96"):
            well_positions(97, rows=8, columns=12)


class TestTheEffectOnReadings:
    def test_a_dataset_records_which_well_each_sample_sat_in(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, edge_effect=0.15)

        assert len(data.wells) == len(data.labels)

    def test_no_effect_leaves_the_readings_alone(self):
        plain = panel_dataset(stressors=["DTT"], noise_cv=0.0, edge_effect=0.0)
        also = panel_dataset(stressors=["DTT"], noise_cv=0.0, edge_effect=0.0)

        assert plain.readings == pytest.approx(also.readings)

    def test_an_edge_well_reads_differently_from_an_interior_one(self):
        data = panel_dataset(stressors=["DTT"], noise_cv=0.0, edge_effect=0.20)
        on_edge = np.array([w in edge_wells() for w in data.wells])

        assert on_edge.any() and (~on_edge).any()
        assert not np.allclose(data.readings[on_edge].mean(axis=0),
                               data.readings[~on_edge].mean(axis=0))

    def test_replication_does_not_average_a_position_effect_away(self):
        """The reason it has to be modelled: it is not noise, and more wells do not help."""
        def edge_bias(replicates):
            data = panel_dataset(stressors=["DTT", "H2O2"], noise_cv=0.0, edge_effect=0.20,
                                 replicates=replicates)
            on_edge = np.array([w in edge_wells() for w in data.wells])
            return abs(data.readings[on_edge].mean() - data.readings[~on_edge].mean())

        assert edge_bias(4) > 0.2 * edge_bias(2)

    def test_a_small_run_can_sit_entirely_in_the_outer_ring(self):
        """Filling row by row puts a short experiment along the top row, which is all edge.
        Worth being able to see, because it is what happens when nobody plans the map."""
        data = panel_dataset(stressors=["DTT"], doses=(0.1, 0.2), replicates=2,
                             noise_cv=0.0, edge_effect=0.20)
        on_edge = [w in edge_wells() for w in data.wells]

        assert all(on_edge)

    def test_the_true_module_activities_are_untouched(self):
        """A position effect is a measurement artifact, not different biology."""
        plain = panel_dataset(stressors=["DTT"], noise_cv=0.0, edge_effect=0.0)
        edged = panel_dataset(stressors=["DTT"], noise_cv=0.0, edge_effect=0.20)

        assert edged.modules == pytest.approx(plain.modules)
