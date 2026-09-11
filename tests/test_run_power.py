"""The power script: the designs it compares, and the shape of the table it writes.

``scripts/run_power.py`` writes ``power_analysis.csv``, which is where the replicate counts
in the protocol come from. The finding in it is an *ordering* -- dosing one agent at a time
confounds stress with growth, and crossing the dose ladder with nutrient levels breaks the
confound -- so the ordering is what is pinned here rather than any particular power figure.

The simulation itself is exercised by ``tests/test_power.py`` against ``analysis/power``.
Running ``main`` costs minutes at the published ``SIMS``, so the written table is checked
where it lies rather than regenerated, and skipped when it is absent.
"""

from __future__ import annotations

import importlib.util
import os
import pathlib

import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_power.py"


@pytest.fixture(scope="module")
def script(tmp_path_factory):
    """The script as a module, with its output directory redirected away from outputs/."""
    redirect = tmp_path_factory.mktemp("power_outputs")
    previous = os.environ.get("YSTWIN_OUTPUTS")
    os.environ["YSTWIN_OUTPUTS"] = str(redirect)
    try:
        spec = importlib.util.spec_from_file_location("run_power", _SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            os.environ.pop("YSTWIN_OUTPUTS", None)
        else:
            os.environ["YSTWIN_OUTPUTS"] = previous
    return module


@pytest.fixture(scope="module")
def designs(script):
    """The three designs the script compares, built exactly as ``main`` builds them."""
    from ystwin.generator.design import decoupling_grid, stressor_only_series
    from ystwin.generator.literature import parameters_from_literature

    params = parameters_from_literature("UPRE2")
    return {
        "dose only (as run)": stressor_only_series(params, doses=script.DOSES),
        "dose x 2 nutrients": decoupling_grid(params, doses=script.DOSES,
                                              nutrient_factors=(1.0, 0.5)),
        "dose x 3 nutrients": decoupling_grid(params, doses=script.DOSES,
                                              nutrient_factors=(1.0, 0.6, 0.35)),
    }


class TestTheDesignsDifferInTheWayTheClaimNeeds:
    """"1.0 = stress and growth are the same axis". The whole table exists to show that the
    design as run sits at that ceiling and that crossing nutrients moves it off. If the
    three designs ever stopped differing, every power figure below would still compute and
    the comparison would mean nothing."""

    def test_dosing_one_agent_at_a_time_confounds_stress_with_growth(self, script,
                                                                     designs):
        from ystwin.generator.design import collinearity

        assert collinearity(designs["dose only (as run)"]) > 0.95

    def test_crossing_nutrient_levels_breaks_the_confound(self, script, designs):
        from ystwin.generator.design import collinearity

        assert collinearity(designs["dose x 2 nutrients"]) < 0.8

    def test_more_nutrient_levels_do_not_make_it_worse(self, script, designs):
        """The third level buys a little more separation and six more conditions. Both
        halves of that trade are in the table, so both are pinned."""
        from ystwin.generator.design import collinearity

        two = collinearity(designs["dose x 2 nutrients"])
        three = collinearity(designs["dose x 3 nutrients"])

        assert three <= two
        assert len(designs["dose x 3 nutrients"]) > len(designs["dose x 2 nutrients"])

    def test_the_dose_ladder_has_an_unstressed_rung(self, script):
        """A slope is fitted through the ladder, and without a zero rung its intercept is
        an extrapolation rather than a measured control."""
        assert min(script.DOSES) == 0.0

    def test_the_effect_sizes_are_folds_above_one(self, script):
        """A "fold" below one would be a reduction, and the power computed for it would
        answer a different question from the one the column heading claims."""
        assert all(fold > 1.0 for fold in script.FOLDS)
        assert list(script.FOLDS) == sorted(script.FOLDS)


@pytest.mark.integration
class TestTheWrittenTable:
    """``power_analysis.csv`` carries three different questions in one frame, each with its
    own columns. That is what makes it readable only through the ``question`` column, and
    what makes a silently renamed question unreadable."""

    @pytest.fixture(scope="class")
    def table(self):
        from ystwin import paths

        path = paths.outputs_dir() / "power_analysis.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not present; run scripts/run_power.py")
        return pd.read_csv(path)

    def test_all_three_questions_are_present(self, table):
        assert set(table.question) == {"detect", "attribute", "chemostat"}

    def test_each_question_is_asked_once_per_design_and_effect(self, table):
        """A duplicated row would be averaged by anything reading this table, and the two
        copies need not have come from the same run."""
        assert not table.duplicated(subset=["design", "question", "fold"]).any()

    def test_every_design_is_asked_the_easy_question_and_the_real_one(self, table):
        """Replicates to detect a slope, and power to attribute it with growth in the
        model. Reporting one without the other is the misreading the script exists to
        prevent."""
        detect = set(table.loc[table.question == "detect", "design"])
        attribute = set(table.loc[table.question == "attribute", "design"])

        assert detect == attribute
        assert len(detect) == 3

    def test_every_interval_brackets_its_own_estimate(self, table):
        """A power figure is a proportion over a finite number of simulations and is
        reported with its interval. An interval that does not contain the estimate is
        arithmetically impossible, so it would mean the columns had been paired wrongly."""
        for value, low, high in (("power_n6", "power_n6_low", "power_n6_high"),
                                 ("power_plate_n4", "power_plate_n4_low",
                                  "power_plate_n4_high"),
                                 ("power_chemostat_n4", "power_chemostat_n4_low",
                                  "power_chemostat_n4_high")):
            rows = table.dropna(subset=[value, low, high])
            assert len(rows), f"no rows carry {value}"
            assert (rows[low] <= rows[value]).all(), value
            assert (rows[value] <= rows[high]).all(), value

    def test_the_table_records_the_simulation_count_it_was_run_at(self, script, table):
        """The interval width is a function of it, so a table produced at a different
        ``SIMS`` from the one in the script is a table whose intervals cannot be
        reproduced."""
        recorded = table.n_simulations.dropna().unique()

        assert list(recorded) == [script.SIMS]

    def test_the_table_records_the_well_level_noise_it_assumed(self, script, table):
        """Power against an assumed CV is not a measurement of power. The assumption
        travels in the table so it cannot be quoted without it."""
        from ystwin.analysis.power import DEFAULT_WELL_CV

        recorded = table.well_cv.dropna().unique()

        assert list(recorded) == [DEFAULT_WELL_CV]

    def test_a_design_that_needs_more_replicates_than_the_ceiling_says_so(self, table):
        """``replicates_needed`` returns ``None`` above ``max_replicates``, which lands in
        the table as an empty cell rather than as the ceiling. Reading the ceiling as an
        answer would say sixteen plates suffice where nothing tried did."""
        detect = table[table.question == "detect"]

        assert detect.replicates.isna().any() or (detect.replicates <= 16).all()
