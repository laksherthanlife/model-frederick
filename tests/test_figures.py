"""What can actually break in a figure, which is never the shape of a curve.

A test that asserted a bar's height would fail on every deliberate restyle and pass on
every wrong number, so nothing here looks at pixels. What is tested is the contract the
figures make with the rest of the project:

* a missing or empty input is refused with a message naming the script that writes it,
  rather than rendering an empty axis that reads as "nothing passed";
* the filenames are a pure function of the figure, so a rebuild overwrites instead of
  accumulating;
* the palette refuses to cycle, because a cycled colour is two categories wearing one
  colour, which is worse than no figure;
* the refusals in the data survive into the drawing -- a fold change with no estimable
  interval and a dose whose recovered activity went negative must both be reachable as
  data, not silently smoothed over.

Everything runs on small synthetic frames under ``tmp_path``. Nothing here reads
``outputs/``; the two tests that do are marked ``integration`` and skip when it is absent.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest

from ystwin.viz import figures

OUTPUTS = pathlib.Path(__file__).resolve().parents[1] / "outputs"


# ---------------------------------------------------------------------------
# synthetic inputs, shaped like the real exports and nothing like their size
# ---------------------------------------------------------------------------


@pytest.fixture
def g1_frame() -> pd.DataFrame:
    """One plate: two passes, two single-criterion failures, one double."""
    return pd.DataFrame({
        "well": ["A1", "A2", "A3", "A4", "A5"],
        "passed": [True, True, False, False, False],
        "failures": ["", "", "linear_range", "sustained_decline",
                     "linear_range,sustained_decline"],
    })


@pytest.fixture
def sensor_frame() -> pd.DataFrame:
    """Two plates, one construct, three doses -- the last one uninvertible.

    The 2 mM row carries a negative recovered activity on both plates, which is what a
    dying culture produces and what the figure has to refuse rather than plot.
    """
    rows = []
    for plate, scale in (("plateA", 1.0), ("plateB", 1.1)):
        for dose, naive, activity in ((0.0, 100.0, 50.0), (0.5, 140.0, 60.0),
                                      (2.0, -20.0, -8.0)):
            rows.append({
                "plate": plate, "construct": "TESTER", "stressor": "DTT",
                "dose_mM": dose, "naive_late": naive * scale,
                "activity_late": activity * scale,
            })
    return pd.DataFrame(rows)


@pytest.fixture
def rt_minus_frame() -> pd.DataFrame:
    """Two constructs against their own anchors, plus the shared UBC reference."""
    rows = []
    for construct, anchor, margin in (("UPRE1", "Hac1", 0.13), ("NativeYap1", "TRX2", 7.8)):
        for target, value in ((anchor, margin), ("UBC", 3.4)):
            for replicate in range(3):
                rows.append({
                    "construct": construct, "target": target, "dose_mM": 0.0,
                    "tech_rep": replicate + 1, "margin_cycles": value + 0.1 * replicate,
                    "passed": value + 0.1 * replicate >= 3.0,
                })
    return pd.DataFrame(rows)


@pytest.fixture
def power_frame() -> pd.DataFrame:
    return pd.DataFrame({
        "arm": ["one reporter"] * 2 + ["ratio"] * 2,
        "arm_label": ["one reporter"] * 2 + ["ratio"] * 2,
        "induction_fold": [1.2, 1.5, 1.2, 1.5],
        "n_replicates": [6, 6, 6, 6],
        "power": [0.0, 0.05, 1.0, 1.0],
        "power_low": [0.0, 0.02, 0.97, 0.97],
        "power_high": [0.03, 0.11, 1.0, 1.0],
        "n_simulations": [120] * 4,
        "collinearity": [0.99, 0.99, 0.99, 0.99],
        "well_cv": [0.1] * 4,
    })


# ---------------------------------------------------------------------------
# the palette refuses rather than cycles
# ---------------------------------------------------------------------------


class TestPalette:
    def test_gives_distinct_colours_for_every_series_it_accepts(self):
        for n in range(1, len(figures.PALETTE) + 1):
            colours = figures.series_colours(n)
            assert len(colours) == n
            assert len(set(colours)) == n, "two series would share a colour"

    def test_is_stable_across_calls(self):
        assert figures.series_colours(3) == figures.series_colours(5)[:3]

    def test_refuses_more_series_than_it_can_separate(self):
        with pytest.raises(ValueError, match="colour-blind-safe palette"):
            figures.series_colours(len(figures.PALETTE) + 1)

    def test_refuses_a_figure_with_no_series(self):
        with pytest.raises(ValueError, match="at least one series"):
            figures.series_colours(0)

    def test_never_puts_pure_red_against_pure_green(self):
        """The one distinction a colour-blind reader cannot make."""
        assert "#FF0000" not in figures.PALETTE
        assert "#00FF00" not in figures.PALETTE


# ---------------------------------------------------------------------------
# a missing input is an error, not an empty figure
# ---------------------------------------------------------------------------


class TestReadTable:
    def test_names_the_script_that_writes_a_missing_table(self, tmp_path):
        with pytest.raises(FileNotFoundError) as exc:
            figures.read_table(tmp_path / "absent.csv", "the widget table",
                               "python3 scripts/run_widgets.py")
        assert "the widget table" in str(exc.value)
        assert "scripts/run_widgets.py" in str(exc.value)

    def test_refuses_a_zero_byte_file(self, tmp_path):
        path = tmp_path / "empty.csv"
        path.write_text("")
        with pytest.raises(ValueError, match="is empty"):
            figures.read_table(path, "the widget table", "python3 scripts/run_widgets.py")

    def test_refuses_a_header_with_no_rows(self, tmp_path):
        path = tmp_path / "headers.csv"
        path.write_text("well,passed\n")
        with pytest.raises(ValueError, match="no rows"):
            figures.read_table(path, "the widget table", "python3 scripts/run_widgets.py")

    def test_reads_a_real_table(self, tmp_path, g1_frame):
        path = tmp_path / "g1_x.csv"
        g1_frame.to_csv(path, index=False)
        assert len(figures.read_table(path, "G1", "run_gates.py")) == len(g1_frame)


# ---------------------------------------------------------------------------
# 01 - G1
# ---------------------------------------------------------------------------


class TestG1:
    def test_signatures_partition_the_wells(self, g1_frame):
        summary = figures.g1_summary({"plate one": g1_frame})
        assert summary.n_wells_with_outcome.sum() == len(g1_frame), (
            "outcomes must partition the plate; counting reasons would double-count "
            "wells that failed two criteria"
        )
        assert set(summary.plate) == {"plate one"}
        assert summary.n_pass.unique().tolist() == [2]

    def test_refuses_with_no_tables_at_all(self):
        with pytest.raises(ValueError, match="no G1 tables"):
            figures.g1_summary({})

    def test_refuses_a_frame_that_is_not_a_gate_report(self):
        with pytest.raises(ValueError, match="missing"):
            figures.g1_summary({"x": pd.DataFrame({"well": ["A1"]})})

    def test_refuses_a_plate_with_no_wells(self):
        empty = pd.DataFrame({"well": [], "passed": [], "failures": []})
        with pytest.raises(ValueError, match="no wells"):
            figures.g1_summary({"x": empty})

    def test_writes_both_formats_under_a_deterministic_stem(self, tmp_path, g1_frame):
        summary = figures.g1_summary({"plate one": g1_frame})
        svg, png = figures.g1_pass_rates(summary, tmp_path)
        assert svg.name == "01_g1_pass_rates.svg"
        assert png.name == "01_g1_pass_rates.png"
        assert svg.stat().st_size > 0 and png.stat().st_size > 0

    def test_rebuilding_overwrites_rather_than_accumulating(self, tmp_path, g1_frame):
        summary = figures.g1_summary({"plate one": g1_frame})
        first = figures.g1_pass_rates(summary, tmp_path)
        second = figures.g1_pass_rates(summary, tmp_path)
        assert first == second
        assert sorted(p.name for p in tmp_path.iterdir()) == [
            "01_g1_pass_rates.png", "01_g1_pass_rates.svg"
        ]

    def test_names_the_plates_it_has_not_read(self, tmp_path, g1_frame):
        summary = figures.g1_summary({"plate one": g1_frame})
        svg, _ = figures.g1_pass_rates(summary, tmp_path,
                                       plates_without_a_table=["2026-07-22"])
        assert "2026-07-22" in svg.read_text(encoding="utf-8")

    def test_refuses_an_empty_summary(self, tmp_path):
        with pytest.raises(ValueError, match="no rows"):
            figures.g1_pass_rates(pd.DataFrame(columns=["plate", "outcome"]), tmp_path)

    def test_refuses_more_outcomes_than_the_palette_can_separate(self, tmp_path):
        many = pd.DataFrame([
            {"plate": "p", "outcome": f"criterion_{i}", "n_wells_with_outcome": 1,
             "n_wells": 9, "n_pass": 0, "pass_rate": 0.0, "n_criteria": 1}
            for i in range(len(figures.PALETTE) + 1)
        ])
        with pytest.raises(ValueError, match="palette"):
            figures.g1_pass_rates(many, tmp_path)


# ---------------------------------------------------------------------------
# 02 - dose response
# ---------------------------------------------------------------------------


class TestDoseResponse:
    def test_reports_both_folds_and_the_refusal_at_two_plates(self, sensor_frame):
        table = figures.dose_response_table(sensor_frame)
        row = table[np.isclose(table.dose_mM, 0.5)].iloc[0]
        assert row.n_plates == 2
        assert not row.interval_estimable, "two plates cannot support a bootstrap"
        assert np.isnan(row.interval_low) and np.isnan(row.interval_high)
        assert row.naive_fold == pytest.approx(1.4)
        assert row.corrected_fold == pytest.approx(1.2)

    def test_marks_a_dose_the_inversion_could_not_invert(self, sensor_frame):
        table = figures.dose_response_table(sensor_frame)
        dying = table[np.isclose(table.dose_mM, 2.0)].iloc[0]
        assert dying.n_plates_usable_corrected == 0, (
            "a fold from a negative activity is not a small fold change"
        )

    def test_refuses_a_frame_missing_the_correction(self, sensor_frame):
        with pytest.raises(ValueError, match="activity_late"):
            figures.dose_response_table(sensor_frame.drop(columns=["activity_late"]))

    def test_refuses_an_empty_frame(self):
        empty = pd.DataFrame(columns=["plate", "construct", "stressor", "dose_mM",
                                      "naive_late", "activity_late"])
        with pytest.raises(ValueError, match="no rows"):
            figures.dose_response_table(empty)

    def test_refuses_a_construct_with_only_a_control(self, sensor_frame):
        control_only = sensor_frame[np.isclose(sensor_frame.dose_mM, 0.0)]
        with pytest.raises(ValueError, match="only its control dose"):
            figures.dose_response_table(control_only)

    def test_writes_both_formats_under_a_deterministic_stem(self, tmp_path, sensor_frame):
        table = figures.dose_response_table(sensor_frame)
        svg, png = figures.dose_response_correction(table, tmp_path)
        assert svg.name == "02_dose_response_correction.svg"
        assert png.name == "02_dose_response_correction.png"
        assert svg.stat().st_size > 0 and png.stat().st_size > 0

    def test_says_on_its_face_that_no_interval_exists(self, tmp_path, sensor_frame):
        svg, _ = figures.dose_response_correction(
            figures.dose_response_table(sensor_frame), tmp_path)
        text = svg.read_text(encoding="utf-8")
        assert "NO INTERVAL IS ESTIMABLE" in text or "no interval" in text.lower()

    def test_refuses_an_empty_fold_table(self, tmp_path):
        columns = ["construct", "stressor", "dose_mM", "naive_fold", "corrected_fold",
                   "n_plates", "n_plates_usable_corrected", "n_plates_usable_naive"]
        with pytest.raises(ValueError, match="no rows"):
            figures.dose_response_correction(pd.DataFrame(columns=columns), tmp_path)


# ---------------------------------------------------------------------------
# 03 - the decoupling grid
# ---------------------------------------------------------------------------


class TestDecouplingGrid:
    def test_crossing_nutrients_lowers_the_collinearity(self):
        designs = figures.decoupling_designs()
        scores = designs.groupby("design").collinearity.first()
        assert scores["dose only (as run)"] > 0.9
        assert scores["dose x 3 nutrients"] < scores["dose x 2 nutrients"]
        assert scores["dose x 2 nutrients"] < scores["dose only (as run)"]

    def test_is_deterministic(self):
        pd.testing.assert_frame_equal(figures.decoupling_designs(),
                                      figures.decoupling_designs())

    def test_writes_both_formats_under_a_deterministic_stem(self, tmp_path):
        svg, png = figures.decoupling_grid_scatter(figures.decoupling_designs(), tmp_path)
        assert svg.name == "03_decoupling_grid.svg"
        assert png.name == "03_decoupling_grid.png"
        assert svg.stat().st_size > 0 and png.stat().st_size > 0

    def test_refuses_an_empty_design(self, tmp_path):
        columns = ["design", "dose_mM", "nutrient_factor", "growth_rate", "collinearity"]
        with pytest.raises(ValueError, match="no rows"):
            figures.decoupling_grid_scatter(pd.DataFrame(columns=columns), tmp_path)

    def test_refuses_a_design_missing_its_growth_axis(self, tmp_path):
        designs = figures.decoupling_designs().drop(columns=["growth_rate"])
        with pytest.raises(ValueError, match="growth_rate"):
            figures.decoupling_grid_scatter(designs, tmp_path)


# ---------------------------------------------------------------------------
# 04 - attribution power
# ---------------------------------------------------------------------------


class TestAttributionPower:
    def test_writes_both_formats_under_a_deterministic_stem(self, tmp_path, power_frame):
        svg, png = figures.attribution_power(power_frame, tmp_path)
        assert svg.name == "04_attribution_power.svg"
        assert png.name == "04_attribution_power.png"
        assert svg.stat().st_size > 0 and png.stat().st_size > 0

    def test_refuses_a_table_with_no_interval(self, tmp_path, power_frame):
        with pytest.raises(ValueError, match="power_low"):
            figures.attribution_power(power_frame.drop(columns=["power_low"]), tmp_path)

    def test_refuses_an_empty_table(self, tmp_path):
        columns = ["arm", "induction_fold", "power", "power_low", "power_high",
                   "n_simulations"]
        with pytest.raises(ValueError, match="no rows"):
            figures.attribution_power(pd.DataFrame(columns=columns), tmp_path)

    def test_refuses_a_simulation_that_would_run_no_campaigns(self):
        with pytest.raises(ValueError, match="n_simulations"):
            figures.attribution_power_table(n_simulations=0)

    @pytest.mark.slow
    def test_carries_a_wilson_interval_that_brackets_its_point(self):
        table = figures.attribution_power_table(n_simulations=12, folds=(1.5,))
        assert len(table) == 3, "one row per assay arm"
        assert (table.power_low <= table.power).all()
        assert (table.power <= table.power_high).all()
        assert (table.power_high > table.power_low).all(), (
            "a bare point estimate is a number whose last digit is noise"
        )


# ---------------------------------------------------------------------------
# 05 - G4 contamination
# ---------------------------------------------------------------------------


class TestG4Contamination:
    def test_converts_cycles_to_a_share_the_ceiling_can_be_read_against(self, rt_minus_frame):
        table = figures.g4_contamination_table(rt_minus_frame)
        er = table[(table.construct == "UPRE1") & (table.role == "anchor")]
        assert er.gdna_share.median() > figures.CORRECTABLE_GDNA_CEILING
        assert not er.correctable.any()
        oxidative = table[(table.construct == "NativeYap1") & (table.role == "anchor")]
        assert oxidative.gdna_share.median() < 0.01
        assert oxidative.correctable.all()

    def test_names_the_reference_gene_as_a_reference(self, rt_minus_frame):
        table = figures.g4_contamination_table(rt_minus_frame)
        assert set(table[table.target == "UBC"].role) == {"reference"}

    def test_refuses_a_table_without_the_margin(self, rt_minus_frame):
        with pytest.raises(ValueError, match="margin_cycles"):
            figures.g4_contamination_table(rt_minus_frame.drop(columns=["margin_cycles"]))

    def test_refuses_an_empty_table(self):
        empty = pd.DataFrame(columns=["construct", "target", "margin_cycles"])
        with pytest.raises(ValueError, match="no readings"):
            figures.g4_contamination_table(empty)

    def test_refuses_a_verdict_table_it_cannot_join(self, rt_minus_frame):
        with pytest.raises(ValueError, match="verdict"):
            figures.g4_contamination_table(
                rt_minus_frame, pd.DataFrame({"construct": ["UPRE1"], "call": ["x"]}))

    def test_writes_both_formats_under_a_deterministic_stem(self, tmp_path, rt_minus_frame):
        table = figures.g4_contamination_table(rt_minus_frame)
        svg, png = figures.g4_contamination(table, tmp_path)
        assert svg.name == "05_g4_contamination.svg"
        assert png.name == "05_g4_contamination.png"
        assert svg.stat().st_size > 0 and png.stat().st_size > 0

    def test_refuses_when_no_target_is_an_anchor_or_a_reference(self, tmp_path,
                                                                rt_minus_frame):
        table = figures.g4_contamination_table(rt_minus_frame)
        table["role"] = "other"
        with pytest.raises(ValueError, match="no anchor or reference"):
            figures.g4_contamination(table, tmp_path)

    def test_refuses_an_empty_frame(self, tmp_path):
        columns = ["construct", "target", "gdna_share", "role"]
        with pytest.raises(ValueError, match="no readings"):
            figures.g4_contamination(pd.DataFrame(columns=columns), tmp_path)


# ---------------------------------------------------------------------------
# against the committed tables, when they are present
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestAgainstCommittedOutputs:
    def test_the_intervals_are_estimable_now_that_three_plates_parse(self):
        """This asserted the opposite until 20260728 was recovered, and fired as designed
        when it stopped being true. At n=2 a cluster bootstrap has three distinct
        resamples and `fold_change` refuses; at n=3 it does not."""
        path = OUTPUTS / "sensor_characterisation.csv"
        if not path.exists():
            pytest.skip(f"{path} is not present")
        table = figures.dose_response_table(pd.read_csv(path))

        assert table.interval_estimable.any()

    def test_the_er_constructs_induce_at_intermediate_dtt(self):
        """The first non-INCONCLUSIVE real-data verdict in the project. Both UPRE
        constructs clear 1.0 at 0.2, 0.5 and 1.0 mM after dilution correction, and the
        late-window sweep keeps the sign on all six.

        The dose set is asserted and the construct set is asserted per construct, because
        the two moved apart when `20260804`'s instrument file took the panel to n=4: the
        UPRE pair is still exactly 0.2/0.5/1.0 mM, and `AlteredYap1` gained a single call at
        0.5 mM H2O2 -- the first per-dose call either oxidative sensor has produced. A test
        that says "only the ER constructs resolve" would now be asserting the old panel.
        """
        path = OUTPUTS / "sensor_characterisation.csv"
        if not path.exists():
            pytest.skip(f"{path} is not present")
        table = figures.dose_response_table(pd.read_csv(path))
        induced = table[table.interval_estimable & (table.interval_low > 1.0)]

        er = induced[induced.stressor == "DTT"]
        assert set(er.construct) == {"UPRE1", "UPRE2"}
        assert sorted(er.dose_mM.unique()) == [0.2, 0.5, 1.0]
        for construct in ("UPRE1", "UPRE2"):
            assert sorted(er[er.construct == construct].dose_mM) == [0.2, 0.5, 1.0]

        oxidative = induced[induced.stressor == "H2O2"]
        assert list(oxidative.construct) == ["AlteredYap1"]
        assert list(oxidative.dose_mM) == [0.5]

    def test_the_er_anchor_is_still_past_the_correction_ceiling(self):
        path = OUTPUTS / "g4_rt_minus_qc.csv"
        if not path.exists():
            pytest.skip(f"{path} is not present")
        table = figures.g4_contamination_table(pd.read_csv(path))
        anchors = (table[table.role == "anchor"].groupby("construct").gdna_share.median())
        assert anchors["UPRE1"] > figures.CORRECTABLE_GDNA_CEILING
        assert anchors["UPRE2"] > figures.CORRECTABLE_GDNA_CEILING
        assert anchors["NativeYap1"] < 0.01
        assert anchors["AlteredYap1"] < 0.01
