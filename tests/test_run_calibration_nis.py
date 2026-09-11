"""Whether the honesty test is itself honest: the band, the whiteness statistic, the priors.

``scripts/run_calibration_nis.py`` is the only check in the project that asks whether the
posterior's *width* is right rather than whether its centre is, and it needs no ground truth
to do it. That makes it easy to get quietly wrong in three ways, and all three are pinned
here.

The acceptance band has to tighten with the run length. The average of K chi-squared(1) draws
is chi-squared(K)/K, so quoting the single-sample band [0.001, 5.02] for a 25-step average
would let almost anything pass and the verdict would mean nothing.

The priors have to be in the filter's own units. The reporter state is a concentration *per
gDCW* and ``observe_rfu`` multiplies it by biomass, so centring the prior on raw RFU
under-predicts every reading by the biomass factor -- about thirtyfold on these plates -- and
the filter is then scored on a units error rather than on its honesty. The script's own
comment names this as the per-OD versus per-gDCW boundary that has no converter in the
package.

And the priors have to come from the well's opening readings rather than from its outcome,
because a prior fitted to what the well did is not a prior.

Manifest and selection tests replay committed plate text. No test writes result artifacts:
the CLI uses --no-write, and the separate output-boundary test captures writes in memory.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
import sys

import numpy as np
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "run_calibration_nis.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("run_calibration_nis", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TIMES = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
OD = np.array([0.10, 0.12, 0.15, 0.19, 0.24, 0.30])
RFU = np.array([900.0, 1000.0, 1100.0, 1400.0, 1800.0, 2300.0])


class TestTheAcceptanceBandIsForTheAverageAndNotForOneDraw:
    def test_a_consistent_filter_sits_inside_the_band_at_every_run_length(self, script):
        """A time-averaged NIS of one is what a correctly specified filter produces, so a
        band that excluded it would fail every honest filter."""
        for k in (1, 3, 10, 25, 200):
            low, high = script._acceptance_band(k)
            assert low < 1.0 < high

    def test_the_band_tightens_as_the_run_lengthens(self, script):
        narrow = script._acceptance_band(25)
        wide = script._acceptance_band(1)

        assert wide[0] < narrow[0] and narrow[1] < wide[1]

    def test_the_single_sample_band_is_the_one_that_would_pass_anything(self, script):
        """The failure the docstring names, as a number: the one-draw band spans a factor of
        five thousand and the 25-step band a factor of three, so quoting the first for a
        25-step average buys a pass for a filter that is wrong by three orders of
        magnitude."""
        low, high = script._acceptance_band(1)
        assert (low, high) == pytest.approx((0.000982, 5.0239), rel=1e-3)
        assert high / low > 1000.0

        low, high = script._acceptance_band(25)
        assert high / low < 4.0

    def test_a_stricter_level_widens_the_band_rather_than_narrowing_it(self, script):
        """alpha is the probability of rejecting a consistent filter, so a smaller alpha has
        to accept more. A band that moved the other way would reject harder the more
        cautious it was told to be."""
        cautious = script._acceptance_band(25, alpha=0.01)
        default = script._acceptance_band(25, alpha=0.05)

        assert cautious[0] < default[0] and default[1] < cautious[1]


class TestTheWhitenessStatistic:
    def test_a_run_too_short_to_have_a_lag_one_pair_is_not_scored(self, script):
        assert np.isnan(script._lag1_autocorrelation([1.0, -1.0]))

    def test_a_constant_run_has_no_autocorrelation_to_report(self, script):
        """Zero variance is not zero autocorrelation. Returning 0.0 here would count a
        degenerate run as white and let it into the "percent of wells white" figure."""
        assert np.isnan(script._lag1_autocorrelation([2.0, 2.0, 2.0, 2.0]))

    def test_the_denominator_is_the_whole_sum_of_squares(self, script):
        """Pinned as an exact number because the common variant -- dividing by the sum over
        the overlapping window only -- gives a different value on the same innovations, and
        the +/-2/sqrt(K) bound the verdict uses is calibrated for this one."""
        assert script._lag1_autocorrelation([1.0, -1.0, 1.0, -1.0]) == pytest.approx(-0.75)

    def test_innovations_that_alternate_are_reported_as_strongly_negative(self, script):
        assert script._lag1_autocorrelation([1.0, -1.0] * 8) < -0.8

    def test_innovations_that_drift_are_reported_as_strongly_positive(self, script):
        """The shape that means the dynamics are missing a term rather than the noise being
        mis-scaled, which is the whole reason this statistic is not optional."""
        assert script._lag1_autocorrelation(np.arange(16.0)) > 0.7

    def test_the_statistic_stays_inside_its_own_bounds(self, script):
        rng = np.random.default_rng(0)
        for _ in range(20):
            value = script._lag1_autocorrelation(rng.normal(size=40))
            assert -1.0 <= value <= 1.0

    def test_a_non_finite_step_is_not_bridged_when_forming_lag_one_pairs(self, script):
        """Dropping a row must not relabel a lag-two pair as lag one."""
        assert script._lag1_autocorrelation([1.0, -1.0, np.nan, 1.0, -1.0]) == \
            pytest.approx(-0.5)


class TestThePriorsAreInTheFiltersOwnUnits:
    def test_the_reporter_prior_is_a_concentration_per_dry_weight(self, script):
        priors = script._priors_for(OD, RFU, TIMES)

        opening_od = float(np.mean(OD[:3]))
        opening_rfu = float(np.mean(RFU[:3]))
        assert priors.reporter[0] == pytest.approx(
            opening_rfu / (opening_od * script.GDCW_PER_OD))

    def test_and_therefore_is_not_the_raw_reading(self, script):
        """The units error this crossing exists to prevent: centred on raw RFU the filter
        under-predicts every reading by the biomass factor, about thirtyfold here, and the
        NIS then measures the mistake instead of the posterior."""
        priors = script._priors_for(OD, RFU, TIMES)

        opening_rfu = float(np.mean(RFU[:3]))
        assert priors.reporter[0] / opening_rfu > 10.0

    def test_the_biomass_prior_is_the_optical_density_converted_the_same_way(self, script):
        priors = script._priors_for(OD, RFU, TIMES)

        assert priors.biomass[0] == pytest.approx(
            float(np.mean(OD[:3])) * script.GDCW_PER_OD)

    def test_every_prior_is_given_a_spread_rather_than_a_point(self, script):
        """A prior with no width hands the filter the answer and makes the NIS a statement
        about arithmetic."""
        priors = script._priors_for(OD, RFU, TIMES)

        for centre, spread in (priors.biomass, priors.reporter, priors.promoter_activity,
                               priors.growth_rate):
            assert spread > 0.0
            assert np.isfinite(centre)

    def test_the_priors_come_from_the_opening_readings_and_not_from_the_outcome(self, script):
        """Tuning the priors against what the well went on to do would make the honesty test
        circular. Only the first three timepoints may reach them."""
        from_opening = script._priors_for(OD, RFU, TIMES)

        exploded = RFU.copy()
        exploded[3:] *= 50.0
        assert script._priors_for(OD, exploded, TIMES).reporter == from_opening.reporter

    def test_later_od_readings_cannot_change_the_opening_growth_prior(self, script):
        from_opening = script._priors_for(OD, RFU, TIMES)
        altered = OD.copy()
        altered[3:] *= np.array([2.0, 5.0, 20.0])
        assert script._priors_for(altered, RFU, TIMES).growth_rate == from_opening.growth_rate

    def test_the_growth_prior_is_held_inside_a_plausible_range(self, script):
        """A well whose opening optical density falls gives a negative estimate, and a prior
        centred there would put the filter into a regime no culture is in."""
        falling = np.array([0.30, 0.28, 0.26, 0.24, 0.22, 0.20])

        priors = script._priors_for(falling, RFU, TIMES)

        assert 0.05 <= priors.growth_rate[0] <= 0.6


class TestAWellWithNoSignalIsExcludedAndNotClamped:
    def test_a_reporter_at_or_below_background_is_refused_by_name(self, script):
        """Both Yap1 constructs fall below background above 1 mM, where the culture is
        dying. Clamping to a small positive number would invent a signal and then test the
        filter's honesty about it."""
        dead = np.array([-5.0, -3.0, -4.0, 10.0, 20.0, 30.0])

        with pytest.raises(script.BelowBackground, match="opening reporter"):
            script._priors_for(OD, dead, TIMES)

    def test_a_well_with_no_biomass_is_refused_the_same_way(self, script):
        with pytest.raises(script.BelowBackground, match="OD"):
            script._priors_for(np.zeros_like(OD), RFU, TIMES)

    def test_the_refusal_is_a_value_error_so_a_caller_can_count_it(self, script):
        """The script counts these wells and reports the count. A bare exception class the
        loop did not catch would end the run instead."""
        assert issubclass(script.BelowBackground, ValueError)


class TestNoWorkbooksMeansTheCommittedTextAndNotARefusal:
    """This report is a real-data result, and it used to be unavailable to anyone without
    the Synergy workbooks -- which are not committed and cannot be, because their document
    properties name a private individual. `data/plates` holds the same numbers as text, and
    the run now replays them rather than exiting 2.

    That is a change in what the script *can* do, not in what it will claim. The source is
    printed either way, and a replay is labelled as one."""

    def test_absent_workbooks_replay_the_committed_text_instead_of_exiting(
            self, script, monkeypatch, capsys):
        monkeypatch.setattr(script.paths, "biosensor_plates", lambda: None)
        monkeypatch.setattr(sys, "argv", ["run_calibration_nis.py", "--max-wells", "2",
                                          "--particles", "20", "--source", "auto", "--no-write"])

        assert script.main() == 0

        assert "data/plates" in capsys.readouterr().out

    def test_the_source_is_stated_and_not_left_to_be_inferred(self, script, monkeypatch, capsys):
        """The reader used is stated, and an artifact-free run creates no directories."""
        def forbidden(*args, **kwargs):
            raise AssertionError("--no-write reached the output boundary")

        monkeypatch.setattr(script, "_write_outputs", forbidden)
        assert script.main(["--max-wells", "2", "--particles", "20", "--no-write"]) == 0
        output = capsys.readouterr().out
        assert "NOT the workbooks" in output
        assert "no tables or manifest written" in output


class TestTheCalibrationSpecificationIsExecutable:
    def test_manifest_records_exact_inputs_parameters_selection_and_initialized_priors(self, script):
        args = script._parse_args(["--max-wells", "2", "--particles", "20", "--no-write",
                                   "--activity-walk", "0.21", "--growth-deceleration", "0.1"])
        innovations, summary, channels, ledger, manifest = script._collect(args)
        assert len(innovations) > 0
        assert len(summary) == 4
        assert len(channels) == 2
        assert (ledger.status == "included").sum() == 2
        assert "smoke_test_cap" in set(ledger.reason)
        assert manifest["resolved_source"] == "committed"
        assert manifest["configuration"]["activity_walk"] == 0.21
        assert manifest["configuration"]["growth_deceleration"] == 0.1
        assert manifest["configuration"]["positive_noise_mode"] == "mean_preserving"
        assert manifest["schema_version"] == 4
        assert manifest["configuration"]["observation_noise_mode"] == "predicted_scale"
        assert manifest["configuration"]["observation_scale_floor"] == 1e-6
        assert manifest["configuration"]["reading_policy"] == "finite_gaussian"
        noise = manifest["measurement_noise"]
        assert noise["mode"] == "predicted_scale"
        assert noise["family"] == "conditionally_independent_additive_gaussian"
        assert noise["relative_sigma"] == {"od": 0.02, "rfu": 0.03}
        assert noise["observation_scale_floor"] == 1e-6
        assert noise["floor_role"] == "numerical reference-scale floor in each channel's units; not an instrument noise measurement"
        assert noise["prior_predictive"] is True
        assert noise["predictive_measurement_variance"] == "sum_i prior_weight_i * sigma_i**2"
        assert manifest["reading_policy"]["name"] == "finite_gaussian"
        assert manifest["reading_policy"]["finite_values"] == "retain signed OD and RFU, including zero and negative readings"
        assert manifest["reading_policy"]["nan"] == "omit only the missing channel at its original step"
        assert manifest["reading_policy"]["infinity"] == "refuse well with field and original step indices"
        assert len(manifest["initialized_wells"]) == 2
        assert all(row["priors"]["activity_walk"] == 0.21 for row in manifest["initialized_wells"])
        assert all(len(row["sha256"]) == 64 for row in manifest["inputs"] + manifest["code"])
        assert manifest["outputs"] == []
        json.dumps(manifest, allow_nan=False)

    def test_exact_plate_and_well_selection_and_exclusions_are_recorded(self, script):
        export = script.replay.exports(source_set="newprotocol")[0].name
        args = script._parse_args(["--plate", export, "--well", "A1", "--well", "A2",
                                   "--exclude-well", f"{export}::A2", "--particles", "20", "--no-write"])
        innovations, _, _, ledger, manifest = script._collect(args)
        assert set(innovations.plate) == {export}
        assert set(innovations.well) == {"A1"}
        assert manifest["counts"]["included_wells"] == 1
        excluded = ledger[(ledger.plate == export) & (ledger.well == "A2")]
        assert excluded.reason.iloc[0] == "well_selection"

    def test_raw_reader_dispatch_uses_the_recorded_file_format_without_rebinding(self, script, monkeypatch):
        raw = pathlib.Path("/not-created/raw-plates")
        original_glob = pathlib.Path.glob
        monkeypatch.setattr(pathlib.Path, "glob", lambda self, pattern: iter([raw / "present.xpt"])
                            if self == raw else original_glob(self, pattern))
        monkeypatch.setattr(script, "read_xpt_run", lambda path: "gen5")
        monkeypatch.setattr(script, "read_synergy_kinetic", lambda path: "xlsx")
        args = script._parse_args(["--source", "workbooks", "--plates-dir", str(raw), "--no-write"])
        exports, reader, _, committed, _ = script._resolve_sources(args)
        assert not committed
        assert reader(next(p for p in exports if p.suffix == ".xpt")) == "gen5"
        assert reader(next(p for p in exports if p.suffix == ".xlsx")) == "xlsx"
        assert script.read_synergy_kinetic(None) == "xlsx"

    def test_unknown_selectors_are_refused_rather_than_silently_ignored(self, script):
        args = script._parse_args(["--plate", "unrecorded.xlsx", "--no-write"])
        with pytest.raises(ValueError, match="unknown export"):
            script._collect(args)
        args = script._parse_args(["--well", "NOT_A_WELL", "--no-write"])
        with pytest.raises(ValueError, match="matches no"):
            script._collect(args)

    def test_config_can_replay_a_manifest_without_writing_a_configuration_file(self, script, monkeypatch):
        payload = {"configuration": {"particles": 31, "activity_walk": 0.4,
                                      "positive_noise_mode": "median_preserving_legacy",
                                      "observation_noise_mode": "observed_scale_legacy",
                                      "observation_scale_floor": 1e-6,
                                      "reading_policy": "positive_od_legacy",
                                      "no_write": True, "output_dir": "/not-created/nis"}}
        original = pathlib.Path.read_text

        def read(path, *args, **kwargs):
            if path.name == "memory-config.json":
                return json.dumps(payload)
            return original(path, *args, **kwargs)

        monkeypatch.setattr(pathlib.Path, "read_text", read)
        args = script._parse_args(["--config", "memory-config.json", "--particles", "41"])
        assert args.particles == 41
        assert args.activity_walk == 0.4
        assert args.positive_noise_mode == "median_preserving_legacy"
        assert args.observation_noise_mode == "observed_scale_legacy"
        assert args.observation_scale_floor == 1e-6
        assert args.reading_policy == "positive_od_legacy"
        assert args.no_write
        assert args.output_dir == pathlib.Path("/not-created/nis")

    def test_legacy_observation_noise_is_explicit_in_manifest_and_cli(self, script, capsys):
        options = ["--max-wells", "1", "--particles", "20", "--no-write",
                   "--observation-noise-mode", "observed_scale_legacy",
                   "--positive-noise-mode", "median_preserving_legacy",
                   "--reading-policy", "positive_od_legacy"]
        _, _, _, _, manifest = script._collect(script._parse_args(options))
        assert manifest["measurement_noise"]["mode"] == "observed_scale_legacy"
        assert manifest["measurement_noise"]["family"] == "outcome_dependent_legacy_pseudo_likelihood"
        assert manifest["measurement_noise"]["prior_predictive"] is False
        assert manifest["reading_policy"]["name"] == "positive_od_legacy"
        assert manifest["initialized_wells"][0]["priors"]["observation_noise_mode"] == "observed_scale_legacy"
        assert script.main(options) == 0
        output = capsys.readouterr().out
        assert "NON-PREQUENTIAL" in output
        assert "Observation noise: observed_scale_legacy" in output
        assert "Reading policy: positive_od_legacy" in output

    def test_output_boundary_honours_output_dir_for_every_table_and_manifest(self, script, monkeypatch):
        import pandas as pd

        args = script._parse_args(["--output-dir", "/not-created/nis"])
        created, written = [], {}
        monkeypatch.setattr(pathlib.Path, "mkdir", lambda self, **kwargs: created.append(self))
        monkeypatch.setattr(pd.DataFrame, "to_csv", lambda self, path, **kwargs: written.update({path: "csv"}))
        monkeypatch.setattr(pathlib.Path, "write_text", lambda self, text, **kwargs: written.update({self: text}))
        table = pd.DataFrame({"value": [1]})
        script._write_outputs(args, table, table, table, table, {"schema_version": 2})
        assert created == [args.output_dir]
        assert {path.name for path in written} == set(script._OUTPUT_NAMES)
        assert all(path.parent == args.output_dir for path in written)
        assert json.loads(written[args.output_dir / "nis_manifest.json"]) == {"schema_version": 2}

    @pytest.mark.parametrize("options", [["--max-wells", "0"], ["--particles", "0"],
                                         ["--coverage", "nan"], ["--min-ess", "-1"]])
    def test_invalid_configurations_are_refused(self, script, options):
        with pytest.raises(SystemExit):
            script._parse_args(options)


class TestTheRunnerUsesTheSameFilterItDiagnoses:
    def test_explicit_optics_are_inverted_and_passed_through(self, script):
        optics = script.ReporterOptics(7.0, 20.0, 3.0, 0.8)
        priors = script._priors_for(OD, RFU, TIMES, optics=optics, carotenoid=1.5)
        expected = ((RFU[:3].mean() - optics.background)
                    / (OD[:3].mean() * script.GDCW_PER_OD * np.exp(-0.8 * 1.5))
                    - optics.autofluorescence) / optics.gain
        assert priors.reporter[0] == pytest.approx(expected)
        frame = script._run_well(TIMES, OD, RFU, 20, 0, optics=optics, carotenoid=1.5)
        assert frame.attrs["priors"]["reporter"][0] == pytest.approx(expected)

    def test_loss_and_noise_overrides_reach_the_filter(self, script):
        frame = script._run_well(TIMES, OD, RFU, 20, 0,
                                 prior_overrides={"k_deg": 0.7, "rfu_rel_sigma": 0.2})
        reporter = frame[frame.channel == "rfu"]
        np.testing.assert_allclose(reporter.measurement_variance,
                                   0.2**2 * (reporter.state_variance + reporter.predicted**2))
        assert frame.attrs["priors"]["k_deg"] == 0.7
        assert frame.attrs["priors"]["rfu_rel_sigma"] == 0.2
        assert frame.attrs["priors"]["observation_noise_mode"] == "predicted_scale"
        default = script._priors_for(OD, RFU, TIMES)
        changed = frame.attrs["priors"]
        assert changed["promoter_activity"][0] > default.promoter_activity[0]

    def test_current_readings_cannot_change_the_runner_prior_moments(self, script):
        altered_od, altered_rfu = OD.copy(), RFU.copy()
        altered_od[4], altered_rfu[4] = -0.1, 10000.0
        original = script._run_well(TIMES, OD, RFU, 20, 0)
        altered = script._run_well(TIMES, altered_od, altered_rfu, 20, 0)
        columns = ["predicted", "state_variance", "measurement_variance"]
        np.testing.assert_array_equal(original.loc[original.step_index == 4, columns],
                                      altered.loc[altered.step_index == 4, columns])
        assert original.attrs["priors"] == altered.attrs["priors"]

    def test_initialization_is_labelled_and_not_in_the_scoring_population(self, script):
        frame = script._run_well(TIMES, OD, RFU, 20, 0)
        assert not frame[frame.step_index < 3].scored.any()
        assert frame[frame.step_index >= 3].scored.all()
        frame["plate"], frame["well"] = "plate", "A1"
        summary = script._summarise(frame)
        assert set(summary.k_steps) == {3}
        assert set(summary.initialization_steps) == {3}
        assert set(summary.nominal_coverage) == {0.95}
        required = {"prior_ess", "posterior_ess", "post_resample_ess", "resampled",
                    "prior_unique_particles", "unique_particles", "unique_ancestors"}
        assert required <= set(frame)

    def test_missing_rfu_does_not_discard_a_valid_od_reading(self, script):
        gapped = RFU.copy()
        gapped[4] = np.nan
        frame = script._run_well(TIMES, OD, gapped, 20, 0)
        assert set(frame[frame.step_index == 4].channel) == {"od"}
        assert frame.attrs["missing_rfu_steps"] == 1

    def test_each_well_gets_its_own_reference_band_not_the_median_step_count(self, script):
        import pandas as pd

        short = script._run_well(TIMES, OD, RFU, 20, 0)
        long = script._run_well(np.arange(10) / 2, np.exp(np.arange(10) / 10),
                                100 * np.exp(np.arange(10) / 10), 20, 0)
        short["plate"], short["well"] = "p1", "A1"
        long["plate"], long["well"] = "p2", "A1"
        summary = script._summarise(pd.concat([short, long], ignore_index=True))
        for _, row in summary.iterrows():
            assert (row.band_low, row.band_high) == pytest.approx(script._acceptance_band(row.k_steps))

    def test_aggregation_changes_the_estimand_not_the_median_column(self, script):
        import pandas as pd

        frame = script._run_well(TIMES, OD, RFU, 20, 0)
        frame["plate"], frame["well"] = "p1", "A1"
        summary = script._summarise(frame)
        a = summary.iloc[[0]].copy()
        b, c = a.copy(), a.copy()
        a["mean_nis"], b["mean_nis"], c["mean_nis"] = 1.0, 3.0, 20.0
        b["well"], c["plate"] = "A2", "p2"
        population = pd.concat([a, b, c], ignore_index=True)
        med = script._channel_rows(population, "well_median")[0]
        mean = script._channel_rows(population, "well_mean")[0]
        plate = script._channel_rows(population, "plate_median")[0]
        assert med["median_mean_nis"] == mean["median_mean_nis"] == plate["median_mean_nis"] == 3.0
        assert med["aggregate_mean_nis"] == 3.0
        assert mean["aggregate_mean_nis"] == 8.0
        assert plate["aggregate_mean_nis"] == 11.0
        assert med["verdict"] == "INCONCLUSIVE"


class TestSignedGaussianReadingPolicy:
    @staticmethod
    def _collect_trace(script, monkeypatch, od, rfu, extra=()):
        from types import SimpleNamespace
        import pandas as pd

        export = script.replay.exports(source_set="newprotocol")[0]
        manifest = script.replay.load_manifest()
        frames = {"OD600": pd.DataFrame({"A1": od, "B1": np.zeros(len(TIMES))}, index=TIMES),
                  "mCitrine": pd.DataFrame({"A1": rfu, "B1": np.zeros(len(TIMES))}, index=TIMES)}
        run = SimpleNamespace(raw_channel=lambda name: SimpleNamespace(channel=name),
                              aligned=lambda channel: frames)
        monkeypatch.setattr(script, "_resolve_sources", lambda args:
                            ([export], lambda path: run, manifest, True, "synthetic policy test"))
        args = script._parse_args(["--particles", "20", "--no-write", "--blank-well", "B1", *extra])
        return script._collect(args)

    @pytest.mark.parametrize("channel", ["od", "rfu"])
    @pytest.mark.parametrize("value", [-0.1, 0.0])
    def test_finite_nonpositive_readings_are_conditioned_on_not_relabelled_missing(self, script, channel, value):
        od, rfu = OD.copy(), RFU.copy()
        (od if channel == "od" else rfu)[4] = value
        frame = script._run_well(TIMES, od, rfu, 20, 0)
        reading = frame[(frame.step_index == 4) & (frame.channel == channel)]
        assert len(reading) == 1
        assert reading.observed.iloc[0] == value
        assert reading.scored.iloc[0]
        assert frame.attrs[f"missing_{channel}_steps"] == 0
        assert frame.attrs[f"nonpositive_{channel}_steps"] == 1
        assert frame.attrs["reading_policy"] == "finite_gaussian"

    @pytest.mark.parametrize("channel", ["od", "rfu"])
    def test_nan_omits_only_that_channel_and_records_the_original_step(self, script, channel):
        od, rfu = OD.copy(), RFU.copy()
        (od if channel == "od" else rfu)[4] = np.nan
        frame = script._run_well(TIMES, od, rfu, 20, 0)
        assert set(frame[frame.step_index == 4].channel) == ({"od", "rfu"} - {channel})
        assert frame.attrs[f"missing_{channel}_steps"] == 1
        assert frame.attrs[f"missing_{channel}_step_indices"] == [4]

    @pytest.mark.parametrize("channel", ["od", "rfu"])
    @pytest.mark.parametrize("value", [np.inf, -np.inf])
    @pytest.mark.parametrize("index", [0, 4])
    def test_infinity_is_an_invalid_reading_not_a_missing_channel(self, script, channel, value, index):
        od, rfu = OD.copy(), RFU.copy()
        (od if channel == "od" else rfu)[index] = value
        with pytest.raises(ValueError, match=f"{channel}.*infinite"):
            script._run_well(TIMES, od, rfu, 20, 0)

    @pytest.mark.parametrize("value", [-0.1, 0.0, np.nan])
    def test_full_collection_keeps_a_well_with_valid_opening_and_later_nonpositive_or_missing_od(
            self, script, monkeypatch, value):
        od = OD.copy()
        od[4] = value
        frame, _, _, ledger, manifest = self._collect_trace(script, monkeypatch, od, RFU)
        assert manifest["counts"]["included_wells"] == 1
        assert ledger[ledger.well == "A1"].status.iloc[0] == "included"
        assert "rfu" in set(frame[frame.step_index == 4].channel)
        reading = frame[(frame.step_index == 4) & (frame.channel == "od")]
        assert len(reading) == (0 if np.isnan(value) else 1)
        if not reading.empty:
            assert reading.observed.iloc[0] == value

    @pytest.mark.parametrize("channel", ["od", "rfu"])
    @pytest.mark.parametrize("aggregation", ["well_median", "well_mean", "plate_median"])
    def test_a_channel_with_no_scored_readings_stays_missing_and_inconclusive(
            self, script, monkeypatch, channel, aggregation):
        od, rfu = OD.copy(), RFU.copy()
        (od if channel == "od" else rfu)[3:] = np.nan
        _, summary, channels, _, manifest = self._collect_trace(
            script, monkeypatch, od, rfu, extra=["--aggregation", aggregation])
        assert manifest["counts"]["included_wells"] == 1
        per_well = summary[summary.channel == channel].iloc[0]
        assert per_well.k_steps == 0
        assert "no_scored_predictions" in per_well.exclusion_reason
        result = channels[channels.channel == channel].iloc[0]
        assert np.isnan(result.median_mean_nis)
        assert np.isnan(result.aggregate_mean_nis)
        assert np.isnan(result.median_lag1_autocorr)
        assert result.verdict == "INCONCLUSIVE"

    @pytest.mark.parametrize("channel", ["od", "rfu"])
    def test_infinite_reading_exclusion_names_its_field_and_is_not_an_initialization_failure(
            self, script, monkeypatch, channel):
        od, rfu = OD.copy(), RFU.copy()
        (od if channel == "od" else rfu)[4] = np.inf
        _, _, _, ledger, manifest = self._collect_trace(script, monkeypatch, od, rfu)
        refusal = ledger[ledger.well == "A1"].iloc[0]
        assert manifest["counts"]["included_wells"] == 0
        assert refusal.status == "excluded"
        assert refusal.stage == "reading"
        assert refusal.field == channel
        assert "infinite" in refusal.reason and "[4]" in refusal.reason

    @pytest.mark.parametrize("channel", ["od", "rfu"])
    def test_opening_data_requirements_have_separate_field_specific_exclusions(
            self, script, monkeypatch, channel):
        od, rfu = OD.copy(), RFU.copy()
        (od if channel == "od" else rfu)[:3] = -1.0
        _, _, _, ledger, manifest = self._collect_trace(script, monkeypatch, od, rfu)
        refusal = ledger[ledger.well == "A1"].iloc[0]
        assert manifest["counts"]["included_wells"] == 0
        assert refusal.status == "excluded"
        assert refusal.stage == "initialization"
        assert refusal.field == channel
        assert "opening" in refusal.reason

    def test_supplied_physical_priors_do_not_require_positive_opening_measurements(self, script):
        priors = script._priors_for(OD, RFU, TIMES)
        od, rfu = OD.copy(), RFU.copy()
        od[0], rfu[0], od[1] = -0.1, -2.0, np.nan
        frame = script._run_well(TIMES, od, rfu, 20, 0, priors=priors)
        assert set(frame[frame.step_index == 0].channel) == {"od", "rfu"}
        assert set(frame[frame.step_index == 1].channel) == {"rfu"}
        assert frame.attrs["priors"]["biomass"] == priors.biomass

    def test_the_legacy_positive_od_population_rule_is_an_explicit_refusal_not_a_hidden_drop(self, script):
        od = OD.copy()
        od[4] = -0.1
        with pytest.raises(ValueError, match="positive_od_legacy"):
            script._run_well(TIMES, od, RFU, 20, 0, reading_policy="positive_od_legacy")

    def test_unknown_reading_policy_is_refused(self, script):
        with pytest.raises(ValueError, match="reading_policy"):
            script._run_well(TIMES, OD, RFU, 20, 0, reading_policy="silent_drop")


class TestRecordedWellPopulation:
    EXPORT = "20260804_ER_Oxidative_Replicate4.xpt"

    @staticmethod
    def _collect_plate(script, monkeypatch, od_wells, rfu_wells, *, extra=(), export=None):
        import pandas as pd
        from ystwin.plate.synergy import KineticBlock, SynergyRun

        path = pathlib.Path(export or TestRecordedWellPopulation.EXPORT)
        recorded = script.recorded_for_export(path.name)
        run = SynergyRun(source=path, recorded_plate=recorded, blocks=tuple(
            KineticBlock(channel=name, fluorophore=name, optics="", sheet=name,
                         derived=False, blank_subtracted=False, temperature_c=[],
                         data=pd.DataFrame(wells, index=TIMES))
            for name, wells in (("OD600", od_wells), ("mCitrine", rfu_wells))))
        manifest = script.replay.load_manifest()
        monkeypatch.setattr(script, "_resolve_sources", lambda args:
                            ([path], lambda source: run, manifest, True, "synthetic recorded-population test"))
        return script._collect(script._parse_args(["--particles", "20", "--no-write", *extra]))

    def test_bright_unrecorded_well_is_not_a_culture_and_bad_recorded_culture_is_retained(
            self, script, monkeypatch):
        recorded = script.recorded_for_export(self.EXPORT)
        assert len(recorded.culture_wells) == 84
        assert "A2" in recorded.culture_wells and "H4" not in recorded.well_roles
        bad_od = OD.copy()
        bad_od[4] = np.inf
        od_wells = {"A1": OD + 0.01, "A2": bad_od + 0.01, "H4": OD * 10 + 0.01,
                    **dict.fromkeys(recorded.blank_wells, np.full(len(TIMES), 0.01))}
        rfu_wells = {"A1": RFU + 10, "A2": RFU + 10, "H4": RFU * 10 + 10,
                     **dict.fromkeys(recorded.blank_wells, np.full(len(TIMES), 10.0))}
        assert np.min(od_wells["H4"]) > 0.15 and np.min(rfu_wells["H4"]) > 0

        innovations, summary, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, od_wells, rfu_wells)

        assert set(innovations.well) == {"A1"}
        assert set(summary.well) == {"A1"}
        assert manifest["counts"]["included_wells"] == 1
        entries = ledger.set_index("well")
        outside = entries.loc["H4"]
        assert outside.status == "excluded"
        assert outside.reason == "outside_recorded_layout"
        assert outside.recorded_role == "unrecorded"
        refused = entries.loc["A2"]
        assert refused.status == "excluded" and refused.recorded_role == "culture"
        assert refused.stage == "reading" and refused.field == "od"
        assert "infinite" in refused.reason and "[4]" in refused.reason
        assert set(ledger.loc[ledger.recorded_role == "culture", "well"]) == set(recorded.culture_wells)
        missing = entries.loc["A3"]
        assert missing.status == "excluded"
        assert missing.reason == "missing_od_or_reporter_channel"
        assert not ledger.duplicated(["plate", "well"]).any()
        assert manifest["counts"]["recorded_culture_wells"] == 84
        assert manifest["counts"]["recorded_blank_wells"] == 3
        assert manifest["counts"]["outside_recorded_layout_wells"] == 1
        assert manifest["counts"]["unknown_identity_wells"] == 0
        assert manifest["population"]["plates"] == [dict(
            plate=self.EXPORT, recorded_culture_wells=list(recorded.culture_wells),
            recorded_blank_wells=list(recorded.blank_wells))]
        json.dumps(manifest, allow_nan=False)

    @staticmethod
    def _blanks(script):
        recorded = script.recorded_for_export(TestRecordedWellPopulation.EXPORT)
        return (dict.fromkeys(recorded.blank_wells, np.full(len(TIMES), 0.01)),
                dict.fromkeys(recorded.blank_wells, np.full(len(TIMES), 10.0)))

    @pytest.mark.parametrize("selector", ["H4", f"{EXPORT}::H4"])
    @pytest.mark.parametrize("blank_options", [[], ["--blank-well", "H1"],
                                               ["--blank-policy", "legacy_low_od"]])
    def test_explicit_selectors_cannot_promote_bright_outside_layout_wells(
            self, script, monkeypatch, selector, blank_options):
        od_wells, rfu_wells = self._blanks(script)
        od_wells.update(A1=OD + 0.01, H4=OD * 10 + 0.01)
        rfu_wells.update(A1=RFU + 10, H4=RFU * 10 + 10)
        innovations, summary, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, od_wells, rfu_wells, extra=["--well", selector, *blank_options])
        assert innovations.empty and summary.empty
        outside = ledger.set_index("well").loc["H4"]
        assert outside.reason == "outside_recorded_layout"
        assert outside.status == "excluded" and outside.stage == "population"
        assert outside.recorded_role == "unrecorded"
        assert manifest["counts"]["included_wells"] == 0
        assert manifest["counts"]["recorded_culture_wells"] == 84
        assert manifest["counts"]["outside_recorded_layout_wells"] == 1

    def test_missing_channels_and_absent_recorded_selectors_remain_accounted_for(
            self, script, monkeypatch):
        od_wells, rfu_wells = self._blanks(script)
        od_wells.update(A1=OD + 0.01, A2=OD + 0.01, H4=OD * 10, H6=OD * 10)
        rfu_wells.update(A1=RFU + 10, A3=RFU + 10, H5=RFU * 10, H6=RFU * 10)
        innovations, _, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, od_wells, rfu_wells, extra=["--well", f"{self.EXPORT}::A4"])
        assert innovations.empty
        entries = ledger.set_index("well")
        for well, role in (("A2", "culture"), ("A3", "culture"), ("A4", "culture"),
                           ("H4", "unrecorded"), ("H5", "unrecorded")):
            assert entries.loc[well, "reason"] == "missing_od_or_reporter_channel"
            assert entries.loc[well, "stage"] == "channels"
            assert entries.loc[well, "recorded_role"] == role
        assert entries.loc["A1", "reason"] == "well_selection"
        assert entries.loc["H6", "reason"] == "outside_recorded_layout"
        assert not ledger.duplicated(["plate", "well"]).any()
        assert manifest["counts"]["recorded_culture_wells"] == 84
        assert manifest["counts"]["outside_recorded_layout_wells"] == 3
        assert manifest["counts"]["excluded_entries"] == 87

    @pytest.mark.parametrize("failure", ["missing", "nonfinite"])
    def test_blank_failures_preserve_well_identities_and_missing_channel_rows(
            self, script, monkeypatch, failure):
        od_wells, rfu_wells = self._blanks(script)
        od_wells.update(A1=OD + 0.01, A2=OD + 0.01, H4=OD * 10, H5=OD * 10)
        rfu_wells.update(A1=RFU + 10, H4=RFU * 10)
        if failure == "missing":
            od_wells.pop("H3")
            rfu_wells.pop("H3")
            reason = "recorded_or_selected_blanks_unavailable"
        else:
            rfu_wells["H3"] = np.full(len(TIMES), np.nan)
            reason = "nonfinite_blank"
        innovations, _, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, od_wells, rfu_wells)
        assert innovations.empty
        plate_failure = ledger[ledger.well.isna()].iloc[0]
        assert plate_failure.reason == reason and plate_failure.stage == "blank"
        assert plate_failure.recorded_role is None
        entries = ledger[ledger.well.notna()].set_index("well")
        assert entries.loc["A1", "reason"] == reason
        assert entries.loc["A1", "stage"] == "blank"
        for well in ("A2", "A3", "H5"):
            assert entries.loc[well, "reason"] == "missing_od_or_reporter_channel"
        assert entries.loc["H4", "reason"] == "outside_recorded_layout"
        assert not ledger.duplicated(["plate", "well"]).any()
        assert manifest["counts"]["recorded_culture_wells"] == 84
        assert manifest["counts"]["recorded_blank_wells"] == 3
        assert manifest["counts"]["outside_recorded_layout_wells"] == 2
        assert manifest["counts"]["excluded_entries"] == 90
        json.dumps(manifest, allow_nan=False)

    def test_plate_input_failure_retains_each_known_well_without_counting_the_plate_as_a_culture(
            self, script, monkeypatch):
        path = pathlib.Path(self.EXPORT)
        recorded = script.recorded_for_export(path.name)
        manifest = script.replay.load_manifest()

        def failed_reader(path):
            raise ValueError("ambiguous declared raw reporter")

        monkeypatch.setattr(script, "_resolve_sources", lambda args:
                            ([path], failed_reader, manifest, True, "synthetic plate failure"))
        innovations, _, _, ledger, manifest = script._collect(script._parse_args(
            ["--well", "A1", "--particles", "20", "--no-write"]))
        assert innovations.empty
        assert ledger[ledger.well.isna()].reason.tolist() == ["ambiguous declared raw reporter"]
        well_entries = ledger[ledger.well.notna()]
        assert set(well_entries.well) == set(recorded.well_roles)
        assert set(well_entries.stage) == {"input"}
        assert well_entries.reason.str.startswith("plate_input_unavailable:").all()
        assert manifest["counts"]["recorded_culture_wells"] == 84
        assert manifest["counts"]["recorded_blank_wells"] == 3
        assert manifest["counts"]["excluded_entries"] == 88
        assert not ledger.duplicated(["plate", "well"]).any()

    @pytest.mark.parametrize("selected_blank", ["A2", "H4"])
    def test_explicit_blank_references_change_correction_not_recorded_identity(
            self, script, monkeypatch, selected_blank):
        od_wells, rfu_wells = self._blanks(script)
        od_wells.update(A1=OD + 0.01, A2=np.full(len(TIMES), 0.01), H4=np.full(len(TIMES), 0.01))
        rfu_wells.update(A1=RFU + 10, A2=np.full(len(TIMES), 10.0), H4=np.full(len(TIMES), 10.0))
        innovations, _, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, od_wells, rfu_wells, extra=["--blank-well", selected_blank])
        assert set(innovations.well) == {"A1"}
        entries = ledger.set_index("well")
        assert entries.loc[selected_blank, "status"] == "blank"
        assert entries.loc[selected_blank, "reason"] == "blank_reference"
        assert entries.loc["A2", "recorded_role"] == "culture"
        assert entries.loc["H4", "recorded_role"] == "unrecorded"
        for well in script.recorded_for_export(self.EXPORT).blank_wells:
            assert entries.loc[well, "reason"] == "recorded_blank_not_selected"
            assert entries.loc[well, "recorded_role"] == "blank"
        assert manifest["plate_details"][0]["blank_wells"] == [selected_blank]
        assert manifest["counts"]["recorded_culture_wells"] == 84
        assert manifest["counts"]["recorded_blank_wells"] == 3
        assert manifest["counts"]["outside_recorded_layout_wells"] == 1

    def test_legacy_blank_policy_remains_explicit_without_reclassifying_biological_identity(
            self, script, monkeypatch):
        od_wells, rfu_wells = self._blanks(script)
        od_wells.update(A1=OD + 0.01, A2=np.full(len(TIMES), 0.05), H4=OD * 10)
        rfu_wells.update(A1=RFU + 10, A2=np.full(len(TIMES), 10.0), H4=RFU * 10)
        innovations, _, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, od_wells, rfu_wells, extra=["--blank-policy", "legacy_low_od"])
        assert set(innovations.well) == {"A1"}
        entries = ledger.set_index("well")
        assert entries.loc["A2", "status"] == "blank"
        assert entries.loc["A2", "recorded_role"] == "culture"
        assert entries.loc["H4", "reason"] == "outside_recorded_layout"
        assert set(manifest["plate_details"][0]["blank_wells"]) == {"A2", "H1", "H2", "H3"}
        assert manifest["counts"]["recorded_culture_wells"] == 84

    @pytest.mark.parametrize("options", [[], ["--blank-well", "B1"],
                                        ["--blank-policy", "legacy_low_od"]])
    def test_unknown_layout_requires_the_existing_explicit_exploratory_blank_policy(
            self, script, monkeypatch, options):
        export = "unregistered-plate.xlsx"
        assert script.recorded_for_export(export) is None
        innovations, _, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, {"A1": OD + 0.01, "B1": np.full(len(TIMES), 0.01)},
            {"A1": RFU + 10, "B1": np.full(len(TIMES), 10.0)}, extra=options, export=export)
        assert manifest["counts"]["included_wells"] == bool(options)
        if options:
            assert set(innovations.well) == {"A1"}
        else:
            assert innovations.empty
            assert "recorded_or_selected_blanks_unavailable" in set(ledger.reason)
        assert set(ledger.loc[ledger.well.notna(), "recorded_role"]) == {"unknown"}
        assert manifest["counts"]["recorded_culture_wells"] == 0
        assert manifest["counts"]["outside_recorded_layout_wells"] == 0
        assert manifest["counts"]["unknown_identity_wells"] == 2
        assert manifest["population"]["plates"] == [dict(
            plate=export, recorded_culture_wells=None, recorded_blank_wells=None)]

    def test_typed_record_not_position_or_filename_defines_a_culture(self, script, monkeypatch):
        from ystwin.plate.layout import PlateLayout, RecordedPlate

        export = "opaque-export.xlsx"
        recorded = RecordedPlate(PlateLayout({"reporter": (4,)}, ("H",)), ("A1",), {})
        monkeypatch.setattr(script, "recorded_for_export", lambda name: recorded if name == export else None)
        innovations, _, _, ledger, manifest = self._collect_plate(
            script, monkeypatch, {"H4": OD + 0.01, "A1": np.full(len(TIMES), 0.01), "B2": OD * 10},
            {"H4": RFU + 10, "A1": np.full(len(TIMES), 10.0), "B2": RFU * 10}, export=export)
        assert set(innovations.well) == {"H4"}
        entries = ledger.set_index("well")
        assert entries.loc["H4", "recorded_role"] == "culture"
        assert entries.loc["A1", "recorded_role"] == "blank"
        assert entries.loc["B2", "reason"] == "outside_recorded_layout"
        assert manifest["counts"]["recorded_culture_wells"] == 1

    @pytest.mark.parametrize("selected_plate", [None, EXPORT])
    def test_public_population_is_scoped_by_plate_but_not_by_well_selection_or_smoke_cap(
            self, script, selected_plate):
        options = ["--no-write", "--particles", "20", "--max-wells", "1", "--well", "A1"]
        if selected_plate is not None:
            options += ["--plate", selected_plate]
        _, _, _, ledger, manifest = script._collect(script._parse_args(options))
        expected_exports = script.replay.exports(source_set="newprotocol")
        scoped = [p.name for p in expected_exports if selected_plate is None or p.name == selected_plate]
        assert [p["plate"] for p in manifest["population"]["plates"]] == scoped
        for export in scoped:
            recorded = script.recorded_for_export(export)
            cultures = ledger[(ledger.plate == export) & (ledger.recorded_role == "culture")]
            assert set(cultures.well) == set(recorded.culture_wells)
        outside = ledger[ledger.recorded_role == "unrecorded"]
        assert set(outside.plate) == {self.EXPORT}
        assert set(outside.well) == {f"H{i}" for i in range(4, 13)}
        assert set(outside.reason) == {"outside_recorded_layout"}
        assert set(outside.stage) == {"population"}
        excluded_plates = ledger[ledger.well.isna()]
        assert len(excluded_plates) == len(expected_exports) - len(scoped)
        assert excluded_plates.recorded_role.isna().all()
        assert not ledger.duplicated(["plate", "well"]).any()
        assert manifest["counts"]["recorded_culture_wells"] == 84 * len(scoped)
        assert manifest["counts"]["recorded_blank_wells"] == 3 * len(scoped)
        assert manifest["counts"]["outside_recorded_layout_wells"] == 9
        assert manifest["counts"]["included_wells"] == 1
        assert manifest["counts"]["well_channels"] == 2
        assert manifest["counts"]["excluded_entries"] == int((ledger.status == "excluded").sum())


class TestTheErrorDecompositionArms:
    """Whether an arm changes the one thing it claims to change, and nothing else.

    The decomposition's whole force is that each arm moves one declared quantity. If a
    particle arm also perturbed the noise, or a measurement arm also moved the priors,
    the attribution would be an artifact of the harness rather than a property of the
    filter -- and it would read exactly the same in the table.
    """

    def test_a_particle_arm_changes_the_ensemble_size_and_nothing_else(self, script):
        args = script._parse_args(["--no-write", "--particles", "20"])
        arm = script._arm_args(args, particles=80, sigma_multiplier=1.0)
        moved = {key for key, value in vars(arm).items() if vars(args)[key] != value}
        assert moved == {"particles"}
        assert arm.particles == 80

    def test_a_measurement_arm_scales_both_sigmas_and_nothing_else(self, script):
        args = script._parse_args(["--no-write", "--particles", "20"])
        arm = script._arm_args(args, particles=args.particles, sigma_multiplier=4.0)
        moved = {key for key, value in vars(arm).items() if vars(args)[key] != value}
        assert moved == {"od_rel_sigma", "rfu_rel_sigma"}
        assert arm.od_rel_sigma == pytest.approx(args.od_rel_sigma * 4.0)
        assert arm.rfu_rel_sigma == pytest.approx(args.rfu_rel_sigma * 4.0)

    def test_a_measurement_arm_leaves_the_initialized_priors_untouched(self, script):
        """The claim the arm rests on: `_priors_for` reads only k_deg out of the process
        overrides, so an inflated relative sigma cannot move the state it starts from.
        If it could, the measurement arm would be a process arm wearing its label."""
        declared = script._priors_for(OD, RFU, TIMES, prior_overrides={"od_rel_sigma": 0.02,
                                                                       "rfu_rel_sigma": 0.03})
        inflated = script._priors_for(OD, RFU, TIMES, prior_overrides={"od_rel_sigma": 0.16,
                                                                       "rfu_rel_sigma": 0.24})
        assert declared.biomass == inflated.biomass
        assert declared.reporter == inflated.reporter
        assert declared.promoter_activity == inflated.promoter_activity
        assert declared.growth_rate == inflated.growth_rate
        assert (inflated.od_rel_sigma, inflated.rfu_rel_sigma) == (0.16, 0.24)

    def test_an_arm_that_repeats_the_baseline_is_dropped_rather_than_run_twice(self, script):
        """500 in the particle sweep *is* the 500-particle baseline. Running it again
        under a second name would report Monte Carlo scatter as an effect of the knob."""
        args = script._parse_args(["--no-write", "--particles", "500",
                                   "--decompose-particles", "500", "2000",
                                   "--decompose-sigma", "1", "4"])
        plan = script._arm_plan(args)
        assert [label for label, *_ in plan] == ["baseline", "particles=2000", "sigma_x4"]
        assert [family for _, family, *_ in plan] == ["baseline", "particles", "measurement"]

    def test_the_baseline_is_first_and_carries_the_declared_configuration(self, script):
        args = script._parse_args(["--no-write", "--particles", "300"])
        label, family, particles, multiplier = script._arm_plan(args)[0]
        assert (label, family, particles, multiplier) == ("baseline", "baseline", 300, 1.0)

    @pytest.mark.parametrize("options,message", [
        (["--decompose-particles", "0"], "decompose_particles"),
        (["--decompose-particles", "-4"], "decompose_particles"),
        (["--decompose-sigma", "0"], "decompose_sigma"),
        (["--decompose-sigma", "-2"], "decompose_sigma"),
        (["--decompose-sigma", "nan"], "decompose_sigma"),
    ])
    def test_an_unusable_sweep_is_refused_before_any_arm_runs(self, script, options, message):
        with pytest.raises(SystemExit):
            script._parse_args(["--no-write", *options])

    def test_the_measurement_share_is_read_off_the_run_not_assumed(self, script):
        """The share bounds what an inflated sigma could buy, so assuming it were one
        would make every measurement-only account look cheaper than it is."""
        import pandas as pd

        innovations = pd.DataFrame({
            "channel": ["od", "od", "rfu", "rfu"],
            "state_variance": [3.0, 3.0, 0.0, 0.0],
            "measurement_variance": [1.0, 1.0, 2.0, 2.0],
            "scored": [True, True, True, False],
        })
        assert script._measurement_variance_share(innovations) == {"od": 0.25, "rfu": 1.0}

    def test_unscored_initialization_rows_are_not_in_the_share(self, script):
        import pandas as pd

        innovations = pd.DataFrame({
            "channel": ["od", "od"],
            "state_variance": [0.0, 9.0],
            "measurement_variance": [1.0, 1.0],
            "scored": [True, False],
        })
        assert script._measurement_variance_share(innovations) == {"od": 1.0}

    def test_an_empty_run_has_no_share_rather_than_a_zero_one(self, script):
        import pandas as pd

        assert script._measurement_variance_share(pd.DataFrame()) == {}


class TestTheDecompositionOverRealWells:
    """One end-to-end sweep on committed plate text, small enough to run in the suite."""

    @pytest.fixture(scope="class")
    def decomposition(self, script):
        args = script._parse_args(["--no-write", "--particles", "40", "--max-wells", "3",
                                   "--decompose", "--decompose-particles", "160",
                                   "--decompose-sigma", "4"])
        innovations, _summary, channels, _ledger, _manifest = script._collect(args)
        return script._decompose(args, innovations, channels)

    def test_every_channel_gets_every_arm_exactly_once(self, decomposition):
        counts = decomposition.groupby("channel").label.agg(["count", "nunique"])
        assert set(decomposition.channel) == {"od", "rfu"}
        assert (counts["count"] == 3).all() and (counts["nunique"] == 3).all()
        assert set(decomposition.family) == {"baseline", "particles", "measurement"}

    def test_the_baseline_row_repeats_the_run_it_was_taken_from(self, script, decomposition):
        """The decomposition must not quietly re-run the baseline under different
        settings; its baseline is the run the ordinary tables report."""
        args = script._parse_args(["--no-write", "--particles", "40", "--max-wells", "3"])
        _innovations, _summary, channels, _ledger, _manifest = script._collect(args)
        baseline = decomposition[decomposition.family == "baseline"].set_index("channel")
        for row in channels.to_dict("records"):
            assert baseline.loc[row["channel"], "mean_nis"] == pytest.approx(row["aggregate_mean_nis"])

    def test_inflating_the_declared_sigma_lowers_nis_in_both_channels(self, decomposition):
        """The direction is not in doubt -- NIS is squared innovation over declared
        variance -- and an arm that moved the other way would mean the filter's posterior
        had moved with the likelihood by more than the variance it added."""
        arm = decomposition[decomposition.family == "measurement"].set_index("channel")
        base = decomposition[decomposition.family == "baseline"].set_index("channel")
        for channel in ("od", "rfu"):
            assert arm.loc[channel, "mean_nis"] < base.loc[channel, "mean_nis"]

    def test_more_particles_raises_the_effective_sample_size_it_is_meant_to_probe(
            self, decomposition):
        """A particle arm that did not move ESS would not be probing the particle
        approximation at all, and its flat NIS would mean nothing."""
        arm = decomposition[decomposition.family == "particles"].set_index("channel")
        base = decomposition[decomposition.family == "baseline"].set_index("channel")
        for channel in ("od", "rfu"):
            assert arm.loc[channel, "min_ess"] > base.loc[channel, "min_ess"]

    def test_the_measurement_arm_reports_the_ensemble_it_also_moved(self, decomposition):
        """The one-knob claim is about the declared quantity, not about the effect. On
        this data a x4 sigma multiplies the surviving ensemble many times over, and the
        run has to publish that ratio rather than let the arm read as measurement-only."""
        arm = decomposition[decomposition.family == "measurement"].set_index("channel")
        for channel in ("od", "rfu"):
            assert arm.loc[channel, "min_ess_ratio"] > 1.0
            assert arm.loc[channel, "min_ess_ratio"] == pytest.approx(
                arm.loc[channel, "min_ess"] / decomposition.set_index(
                    ["family", "channel"]).loc[("baseline", channel), "min_ess"])

    def test_the_manifest_names_process_error_as_a_residual_and_not_an_arm(self, script):
        args = script._parse_args(["--no-write", "--particles", "40", "--decompose"])
        manifest = script._decomposition_manifest(args, script._arm_plan(args))
        assert "not an arm" in manifest["families"]["process"]
        assert {entry["family"] for entry in manifest["arms"]} == {
            "baseline", "particles", "measurement"}
        for entry in manifest["arms"]:
            assert entry["od_rel_sigma"] == pytest.approx(args.od_rel_sigma * entry["sigma_multiplier"])
