"""What decides whether a figure can be drawn at all, which is never the drawing.

``scripts/make_figures.py`` is the only place that knows which table feeds which figure, so
the failures it owns are the ones a rendered figure hides: an input table that is not there
(a figure with no series in it renders perfectly and says nothing), a plate that was run and
has no gate report, and more failure signatures than the colour-blind-safe palette can
separate.

The last one is the interesting one. ``figures.g1_pass_rates`` refuses to cycle colours,
because a cycled colour is two failure modes wearing one colour. So this script groups the
rarest signatures until what is left fits -- and a grouping that lost or double-counted a
well would leave a stacked bar that no longer reaches the plate's well count, which looks
exactly like a plate with missing wells. Both properties are pinned here: it fits, and it
is still a partition.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest
from matplotlib.text import Text

from ystwin.plate.layout import RECORDED_PLATES
from ystwin.viz import figures

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "make_figures.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("make_figures", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _g1_frame(signatures: dict[str, int], n_pass: int = 4) -> pd.DataFrame:
    """One plate's gate export: ``n_pass`` passing wells and the named failures."""
    rows = [{"well": f"P{i}", "passed": True, "failures": ""} for i in range(n_pass)]
    for signature, count in signatures.items():
        rows += [{"well": f"{signature}{i}", "passed": False, "failures": signature}
                 for i in range(count)]
    return pd.DataFrame(rows)


def _crowded(n_signatures: int = 9) -> pd.DataFrame:
    """A summary with more failure signatures than the palette can colour.

    Counts fall with the index so the pooling has an unambiguous order to work in: a tie
    at the boundary would make which signature keeps its own colour a property of the sort,
    not of the plate.
    """
    signatures = {f"criterion_{i}": n_signatures - i for i in range(n_signatures)}
    return figures.g1_summary({"plate one": _g1_frame(signatures)})


class TestAPlateLabelCarriesTheDateAndStaysOffItsNeighbour:
    def test_the_date_leads_and_the_lab_filename_follows_on_its_own_line(self, script):
        assert script._plate_label("20260701_ER_preliminary_(RAW)") == \
            "2026-07-01\nER preliminary"

    def test_a_tag_that_is_only_a_date_is_only_a_date(self, script):
        assert script._plate_label("20260722") == "2026-07-22"

    def test_a_long_filename_is_wrapped_rather_than_run_into_the_next_label(self, script):
        """Two plate names overlapping on an axis is worse than either being on three
        lines, so the wrap is part of the label and not a styling choice."""
        label = script._plate_label("20260722_ERandOxidativeStress_NewProtocol_ANALYSED")

        assert label.startswith("2026-07-22\n")
        assert max(len(line) for line in label.splitlines()[1:]) <= 24

    def test_a_stem_with_no_date_in_it_is_passed_through_unchanged(self, script):
        assert script._plate_label("pilot_run") == "pilot_run"


class TestAPlateThatWasRunAndHasNoGateReportIsNamed:
    def test_the_recorded_plates_with_no_table_come_back_as_dates(self, script, tmp_path):
        """The logbook-confirmed list decides what is missing, not a list written into the
        script -- so a plate added to the registry appears in the footnote by itself."""
        (tmp_path / "g1_20260722_something.csv").write_text("well\n")

        missing = script._plates_without_a_table(tmp_path)

        assert "2026-07-22" not in missing
        assert len(missing) == len(RECORDED_PLATES) - 1
        assert all(len(date) == 10 and date.count("-") == 2 for date in missing)

    def test_and_nothing_is_named_once_every_plate_has_one(self, script, tmp_path):
        for date in RECORDED_PLATES:
            (tmp_path / f"g1_{date}_export.csv").write_text("well\n")

        assert script._plates_without_a_table(tmp_path) == []

    def test_no_gate_table_at_all_refuses_and_names_the_script_that_writes_them(
            self, script, tmp_path):
        """Figure 01 *is* the gate report. Drawing it from nothing would render an empty
        axis, which reads as a plate on which nothing passed."""
        with pytest.raises(SystemExit, match="run_gates.py"):
            script._g1_tables(tmp_path)


class TestPoolingRareSignaturesKeepsTheBarHonest:
    def test_a_summary_that_already_fits_the_palette_is_returned_untouched(self, script):
        summary = figures.g1_summary({"plate one": _g1_frame({"linear_range": 3,
                                                              "sustained_decline": 2})})

        assert script._pooled_to_palette(summary) is summary

    def test_a_crowded_summary_comes_back_inside_the_palette(self, script):
        crowded = _crowded()
        assert crowded.outcome.nunique() > len(figures.PALETTE)

        pooled = script._pooled_to_palette(crowded)

        assert pooled.outcome.nunique() <= len(figures.PALETTE)

    def test_and_the_figure_that_refused_the_crowded_one_draws_the_pooled_one(self, script,
                                                                             tmp_path):
        """The end of the chain, and the reason this function exists. Asserting only the
        outcome count would pass a grouping that left one colour too many for the palette
        after ``passed`` is counted."""
        crowded = _crowded()
        with pytest.raises(ValueError, match="palette"):
            figures.g1_pass_rates(crowded, tmp_path)

        svg, png = figures.g1_pass_rates(script._pooled_to_palette(crowded), tmp_path)

        assert svg.exists() and png.exists()

    def test_every_well_still_appears_in_exactly_one_outcome(self, script):
        """A stacked bar is only honest if the outcomes partition the wells. Pooling by
        summing counts can double-count or drop, and either leaves a bar that no longer
        reaches the plate's well count."""
        pooled = script._pooled_to_palette(_crowded())

        for _, plate in pooled.groupby("plate"):
            assert int(plate.n_wells_with_outcome.sum()) == int(plate.n_wells.iloc[0])

    def test_the_commonest_signatures_keep_their_own_colour(self, script):
        pooled = script._pooled_to_palette(_crowded())

        assert "criterion_0" in set(pooled.outcome)
        assert "criterion_8" not in set(pooled.outcome)

    def test_the_pooled_row_says_how_many_signatures_went_into_it(self, script):
        """Otherwise the figure implies there was one rare failure mode where there were
        several, which is a claim about the plates rather than about the palette."""
        crowded = _crowded(n_signatures=9)

        pooled = script._pooled_to_palette(crowded)

        label = next(o for o in set(pooled.outcome) if o.startswith("other ("))
        rolled_up = crowded.outcome.nunique() - pooled.outcome.nunique() + 1
        assert f"other ({rolled_up} rarer signatures)" == label

    def test_a_plate_carrying_none_of_the_rare_signatures_gets_no_pooled_row(self, script):
        """The pooled outcome is per plate, so a plate that only ever failed the common way
        must not gain a bar segment of height zero."""
        crowded = {f"criterion_{i}": 9 - i for i in range(9)}
        summary = figures.g1_summary({
            "crowded": _g1_frame(crowded),
            "clean": _g1_frame({"criterion_0": 2}),
        })

        pooled = script._pooled_to_palette(summary)

        clean = pooled[pooled.plate == "clean"]
        assert set(clean.outcome) == {"passed", "criterion_0"}

    def test_pooling_twice_changes_nothing_further(self, script):
        """The result has to be a fixed point: a second pass over an already-pooled summary
        would otherwise fold the kept signatures into the pooled label as well."""
        once = script._pooled_to_palette(_crowded())

        twice = script._pooled_to_palette(once)

        assert twice is once


@pytest.fixture(scope="module")
def story_tables():
    names = (
        "sensor_characterisation.csv", "heldout_interpolation_by_dose.csv",
        "external_redundancy.csv", "external_nulls.csv", "external_nulls_independent.csv",
        "cross_family_transfer.csv", "late_window_sensitivity.csv",
        "autofluorescence_sensitivity.csv", "g4_rt_minus_qc.csv", "g4_verdicts.csv",
    )
    tables = {name: pd.read_csv(_REPO / "outputs" / name) for name in names}
    tables["folds"] = figures.dose_response_table(tables["sensor_characterisation.csv"])
    return tables


@pytest.fixture
def render_story(script, story_tables, monkeypatch, tmp_path):
    """Exercise the actual recipe; stub only unrelated reports and statistical work."""
    captured = {}

    def pair(out_dir, stem):
        out_dir.mkdir(parents=True, exist_ok=True)
        paths = (out_dir / f"{stem}.svg", out_dir / f"{stem}.png")
        for path in paths:
            path.write_bytes(b"test save seam")
        return paths

    def save(fig, out_dir, stem):
        if stem in captured:
            figures.plt.close(captured[stem])
        captured[stem] = fig
        return pair(out_dir, stem)

    def unrelated_report(_table, out_dir, **_kwargs):
        return pair(out_dir, "unrelated_report")

    monkeypatch.setattr(figures, "save_figure", save)
    monkeypatch.setattr(script, "_g1_tables", lambda _: {})
    monkeypatch.setattr(script, "_plates_without_a_table", lambda _: [])
    monkeypatch.setattr(script, "_pooled_to_palette", lambda frame: frame)
    monkeypatch.setattr(figures, "g1_summary", lambda _: pd.DataFrame({"dummy": [0]}))
    monkeypatch.setattr(figures, "decoupling_designs", lambda: pd.DataFrame({"dummy": [0]}))
    monkeypatch.setattr(figures, "attribution_power_table",
                        lambda **_: pd.DataFrame({"dummy": [0]}))
    monkeypatch.setattr(figures, "g4_contamination_table",
                        lambda *_: pd.DataFrame({"dummy": [0]}))
    for name in ("g1_pass_rates", "decoupling_grid_scatter", "attribution_power",
                 "g4_contamination"):
        monkeypatch.setattr(figures, name, unrelated_report)
    monkeypatch.setattr(figures, "identifiable_ec50", lambda _: pd.DataFrame({
        "construct": ["UPRE1"], "ec50": [0.25], "identifiable": [True],
    }))
    monkeypatch.setattr(script.sys, "argv", [
        str(_SCRIPT), "--figures-dir", str(tmp_path / "report"),
        "--story-dir", str(tmp_path / "story"),
    ])

    def render(**overrides):
        tables = {name: frame.copy(deep=True) for name, frame in story_tables.items()}
        tables.update(overrides)
        before = {name: frame.copy(deep=True) for name, frame in tables.items()}
        monkeypatch.setattr(figures, "read_table", lambda path, *_: tables[path.name])
        monkeypatch.setattr(figures, "dose_response_table", lambda _: tables["folds"])
        script.main()
        for name, frame in tables.items():
            pd.testing.assert_frame_equal(frame, before[name])
        pd.testing.assert_frame_equal(
            pd.read_csv(tmp_path / "report" / "dose_response_folds.csv"), tables["folds"])
        return captured

    yield render
    for fig in captured.values():
        figures.plt.close(fig)


def _figure_text(fig):
    return " ".join(" ".join(artist.get_text() for artist in fig.findobj(Text)).split())


def _marker_count(fig, label):
    return sum(len(artist.get_offsets()) for ax in fig.axes for artist in ax.collections
               if artist.get_label() == label)


class TestAuthoritativeStoryPresentation:
    def test_all_24_conditions_and_five_ci_refusals_survive_the_recipe(
            self, render_story, story_tables):
        folds = story_tables["folds"]
        assert len(folds) == 24
        assert int(folds.interval_estimable.sum()) == 19
        assert int(folds.corrected_fold.isna().sum()) == 3
        drawn = render_story()
        for stem in ("fig01_dose_response_correction", "02_dose_response_correction"):
            fig = drawn[stem]
            text = _figure_text(fig)
            for phrase in ("24 dose conditions", "19 estimable", "7 exclude 1.0",
                           "12 span 1.0", "5 CI REFUSED", "2 finite folds", "3 missing folds"):
                assert phrase in text
            expected = {"CI excludes 1.0": 7, "CI spans 1.0": 12,
                        "CI REFUSED (finite fold)": 2, "CI REFUSED (fold missing)": 3}
            assert {label: _marker_count(fig, label) for label in expected} == expected
            assert all(len(ax.get_xticks()) == 6 for ax in fig.axes)
            assert "pointwise, not family-wise" in text
            for obsolete in ("3 plates x 3", "culture is dying", "two of the resolved calls",
                             "Every hollow marker is INCONCLUSIVE", "inductions vanish"):
                assert obsolete not in text

    def test_coverage_uses_plate_and_well_identity_per_condition_not_a_maximum(
            self, render_story, story_tables):
        sensor = story_tables["sensor_characterisation.csv"].copy()
        at_dose = (sensor.construct == "UPRE1") & (sensor.dose_mM == 0.2)
        plates = sorted(sensor.loc[at_dose, "plate"].unique())
        sensor = sensor[~(at_dose & (sensor.plate == plates[-1]))].copy()
        one_well = sensor.index[(sensor.construct == "UPRE1") & (sensor.dose_mM == 0.2)][0]
        sensor = sensor.drop(index=one_well)
        sensor["n_points"] = 9999
        folds = story_tables["folds"].copy()
        mask = (folds.construct == "UPRE1") & (folds.dose_mM == 0.2)
        folds.loc[mask, ["n_plates", "n_plates_usable_naive",
                         "n_plates_usable_corrected"]] = 3
        fig = render_story(**{"sensor_characterisation.csv": sensor, "folds": folds})[
            "fig01_dose_response_correction"]
        panel = next(ax for ax in fig.axes if ax.get_title().startswith("UPRE1"))
        ticks = [tick.get_text() for tick in panel.get_xticklabels()]
        assert "0.2\n3/8\n3/3" in ticks
        assert "0.1\n4/12\n4/4" in ticks
        text = _figure_text(fig)
        assert "3-4 recorded plates" in text and "8-12 recorded wells" in text
        assert "2-3 wells per plate" in text
        assert "P/W: recorded plates/wells" in text
        assert "N/C: usable naive/corrected plate folds" in text
        assert "9999" not in text

    @pytest.mark.parametrize("column,value", [
        ("interval_low", np.nan), ("interval_high", np.nan), ("interval_estimable", False),
    ])
    def test_a_finite_point_or_stale_verdict_never_supplies_missing_or_refused_bounds(
            self, render_story, story_tables, column, value):
        folds = story_tables["folds"].copy()
        mask = (folds.construct == "UPRE1") & (folds.dose_mM == 0.2)
        folds.loc[mask, column] = value
        fig = render_story(folds=folds)["fig01_dose_response_correction"]
        text = _figure_text(fig)
        assert "18 estimable" in text and "6 CI REFUSED" in text
        assert _marker_count(fig, "CI REFUSED (finite fold)") == 3
        assert _marker_count(fig, "CI excludes 1.0") == 6

    def test_fig02_draws_the_fixed_form_column_and_labels_its_selected_baseline(
            self, render_story, story_tables):
        fig = render_story()["fig02_heldout_skill_vs_saturation"]
        text = _figure_text(fig)
        assert "Fixed-form biphasic diagnostic" in text
        assert "1 - rmse_biphasic / rmse_selected_baseline" in text
        assert "not skill_selected_model" in text
        assert "inner leave-dose-out on outer training only" in text
        assert "nearest_dose = AlteredYap1, NativeYap1" in text
        assert "train_mean = UPRE1, UPRE2" in text
        assert "all fitted constructs; scoring uses only the vouched cohort" in text
        line = next(line for line in fig.axes[0].lines
                    if line.get_label() == "fixed biphasic skill")
        source = story_tables["heldout_interpolation_by_dose.csv"]
        np.testing.assert_allclose(line.get_ydata(), source.skill, equal_nan=True)
        assert not np.allclose(line.get_ydata()[:5], source.skill_selected_model.iloc[:5])
        for obsolete in ("skill against the nearest-dose baseline", "beats persistence",
                         "wins near", "loses in saturation", "outside its EC50 range should"):
            assert obsolete not in text

    def test_fig02_keeps_all_six_requests_and_separates_refused_inventory_from_test_counts(
            self, render_story):
        fig = render_story()["fig02_heldout_skill_vs_saturation"]
        text = _figure_text(fig)
        for phrase in ("6 requested doses", "5 scored", "1 REFUSED", "264 recorded wells",
                       "144/240 test wells scored", "96 excluded", "4 mM: 24 recorded wells",
                       "split counts not reported", "outer-training identifiable biphasic fits"):
            assert phrase in text
        assert _marker_count(fig, "REFUSED (not a score)") == 1
        np.testing.assert_allclose(fig.axes[-1].get_xticks(), np.array([.1, .2, .5, 1, 2, 4]) / .25)
        assert "24/48 scored; 24 excluded" in text
        assert "36/48 scored; 12 excluded" in text

    def test_fig02_analytic_curve_is_conditional_illustration_not_evidence(self, render_story):
        fig = render_story()["fig02_heldout_skill_vs_saturation"]
        text = _figure_text(fig)
        assert "illustration, conditional on full-data K" in text
        assert "not out-of-sample evidence or causal proof" in text
        assert "does not establish a crossover or an extrapolation rule" in text
        assert "has no free parameter" not in text

    def test_fig03_preserves_nonrejections_and_zero_peer_recovery_without_equivalence(
            self, render_story):
        fig = render_story()["fig03_transfer_honesty"]
        text = _figure_text(fig)
        for p_value in ("0.169", "0.398", "0.940", "1.000", "0.429"):
            assert f"p = {p_value}" in text
        assert text.count("does not reject") >= 5
        assert "40% recovery with 0 peers (5 pairs, panel target lists)" in text
        assert "descriptive association" in text
        assert "Non-rejection is not equivalence or proof of no generalisation" in text
        assert "Known-loadings control and fitted-basis transfer ask different questions" in text
        assert "peer-count association (Spearman)" in text
        assert "Transfer uses own-target scores on arrays and held-out channel R2 in simulation" in text
        for obsolete in ("not otherwise", "what is known works", "covers twice",
                         "Transfer is interpolation across redundant coverage"):
            assert obsolete not in text

    def test_fig04_is_only_the_selected_seven_of_nineteen_across_the_declared_sweeps(
            self, render_story, story_tables):
        fig = render_story()["fig04_interval_robustness"]
        text = _figure_text(fig)
        for phrase in ("7/7 selected pointwise calls", "selected 7 of 19 estimable",
                       "77/77", "not family-wise acceptance", "only the declared sweep settings"):
            assert phrase in text
        assert "none of them is an artefact" not in text
        assert "two of the resolved calls" not in text
        calls = story_tables["folds"].query("interval_estimable and interval_low > 1")
        for ax, name, knob in zip(fig.axes, ("late_window_sensitivity.csv",
                                           "autofluorescence_sensitivity.csv"),
                                 ("late_fraction", "autofluorescence_fraction")):
            source = story_tables[name]
            for call in calls.itertuples():
                line = next(line for line in ax.lines
                            if line.get_label() == f"{call.construct} {call.dose_mM:g} mM")
                expected = source[(source.construct == call.construct)
                                  & (source.dose_mM == call.dose_mM)].sort_values(knob)
                np.testing.assert_allclose(line.get_xdata(), expected[knob])
                np.testing.assert_allclose(line.get_ydata(), expected.low)

    def test_fig04_a_missing_sweep_bound_is_not_silently_counted_as_retained(
            self, render_story, story_tables):
        sweep = story_tables["late_window_sensitivity.csv"].copy()
        index = sweep.index[(sweep.construct == "UPRE1") & (sweep.dose_mM == 0.2)][0]
        sweep.loc[index, "low"] = np.nan
        fig = render_story(**{"late_window_sensitivity.csv": sweep})[
            "fig04_interval_robustness"]
        text = _figure_text(fig)
        assert "6/7 selected pointwise calls" in text
        assert "76/77" in text
        assert _marker_count(fig, "CI REFUSED / missing sweep row") == 1

    def test_story_captions_fit_the_canvas_and_do_not_overlap_the_panels(self, render_story):
        for stem, fig in render_story().items():
            if not stem.startswith("fig0"):
                continue
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            caption = fig.texts[-1].get_window_extent(renderer)
            assert caption.x0 >= fig.bbox.x0 and caption.x1 <= fig.bbox.x1, stem
            assert caption.y0 >= fig.bbox.y0, stem
            assert caption.y1 < min(ax.get_tightbbox(renderer).y0 for ax in fig.axes), stem
