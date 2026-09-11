"""The auto feedback loop plates, and what a fourth plate did to the finding.

`outputs/afl_circuit.csv` is tracked, so the assertions about the RESULT run on any clone.
The assertions that need the instrument files skip without them, the same contract every
real-data fixture in this repository honours.

Two things are pinned here, and the second one exists because it used to be wrong.

1.  **The finding.** Rank correlation with the galactose dose is positive for the corrected
    activity in every strain on every plate, and zero or negative for naive `RFU/OD`. The
    two readouts disagree about the sign of the dose response.
2.  **The claim that broke.** The finding used to be stated as endpoint ratios -- "naive
    `RFU/OD` falls from 0% to 1.8%". On 2026-07-05 it rises, and on any plate it moves with
    how long the plate was read. The counterexample is asserted so it cannot quietly come
    back.
"""

from __future__ import annotations

import itertools
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from ystwin import paths

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

DEBUGGING_STRAINS = ("AFL", "Control", "No-Gal", "No-LacI")
DOSES = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8]
MATCHED = "first-109"

#: The two plates whose run is exactly the matched length, so truncating them changes nothing.
ALREADY_MATCHED = ("20260803_AFL-debugging_2nd", "20260804_AFL-debugging_3rd")


@pytest.fixture(scope="module")
def circuit():
    """The committed result table. Tracked, so this never skips."""
    path = paths.outputs_dir() / "afl_circuit.csv"
    if not path.exists():
        pytest.skip(f"{path} not present; run scripts/score_afl_circuit.py")
    return pd.read_csv(path)


def blocks(circuit, truncation):
    """(plate, strain) -> that block's ladder, sorted by dose."""
    frame = circuit[circuit.truncation == truncation]
    return {key: block.sort_values("galactose_percent")
            for key, block in frame.groupby(["plate", "strain"])}


class TestTheseAreTwoExperimentsNotFourReplicates:
    """2026-07-05 is the first replicate of "AFL vs control", a two-strain experiment run a
    month before the four-strain debugging plates. Pooling them would be a mistake, so the
    table keeps the series apart and this says so."""

    def test_two_series(self, circuit):
        assert set(circuit.experiment) == {"AFL-vs-control", "AFL-debugging"}

    def test_the_july_plate_carries_only_the_two_strains_that_existed_then(self, circuit):
        july = circuit[circuit.experiment == "AFL-vs-control"]

        assert set(july.plate) == {"20260705_AFL&control_1st"}
        assert set(july.strain) == {"AFL", "Control"}

    def test_the_august_plates_carry_four(self, circuit):
        august = circuit[circuit.experiment == "AFL-debugging"]

        assert august.plate.nunique() == 3
        assert set(august.strain) == set(DEBUGGING_STRAINS)

    def test_the_july_layout_is_marked_as_inferred_and_the_august_one_as_recorded(self, circuit):
        """The August notebook maps head their column blocks with strain names. The
        2026-07-05 map does not, so its assignment rests on the readings and the table has
        to say so rather than presenting it as a record."""
        source = circuit.groupby("experiment").strain_source.unique()

        assert list(source["AFL-vs-control"]) == ["inferred-from-signal"]
        assert list(source["AFL-debugging"]) == ["notebook"]

    def test_seven_galactose_doses_evenly_spaced_on_both_maps(self, circuit):
        assert sorted(circuit.galactose_percent.unique()) == DOSES

    def test_every_strain_and_dose_appears_on_every_plate_of_its_series(self, circuit):
        counts = circuit[circuit.truncation == MATCHED].groupby(
            ["experiment", "strain", "galactose_percent"]).plate.nunique()
        expected = {"AFL-vs-control": 1, "AFL-debugging": 3}

        for (experiment, _, _), n in counts.items():
            assert n == expected[experiment]

    def test_the_july_blocks_are_six_wells_and_the_august_blocks_three(self, circuit):
        by_series = circuit.groupby("experiment").n_wells.max()

        assert by_series["AFL-vs-control"] == 6
        assert by_series["AFL-debugging"] == 3


class TestTheMatchedTruncationIsATruncation:
    """`RFU/OD` at the end of a run is how much stable reporter has piled up since the
    culture stopped dividing, so plates stopped at 18.1, 20.1 and 24.1 hours cannot be
    compared as they stand. All four files record the same 10-minute grid, so the first 109
    reads of each are literally the same timepoints."""

    def test_the_two_shortest_runs_are_unchanged_by_it(self, circuit):
        """If truncation were resampling or reindexing, these would move."""
        for plate in ALREADY_MATCHED:
            plate_rows = circuit[circuit.plate == plate]
            as_run = plate_rows[plate_rows.truncation == "as-run"].sort_values(
                ["strain", "galactose_percent"])
            matched = plate_rows[plate_rows.truncation == MATCHED].sort_values(
                ["strain", "galactose_percent"])

            assert np.allclose(as_run.naive_rfu_per_od, matched.naive_rfu_per_od)
            assert np.allclose(as_run.corrected_activity, matched.corrected_activity)

    def test_it_leaves_every_plate_at_the_same_length_and_window(self, circuit):
        matched = circuit[circuit.truncation == MATCHED]

        assert matched.n_timepoints.nunique() == 1
        assert matched.run_hours.nunique() == 1
        assert matched.activity_window_h.nunique() == 1
        assert matched.late_window_start_h.nunique() == 1

    def test_and_the_plates_really_were_run_to_different_lengths(self, circuit):
        as_run = circuit[circuit.truncation == "as-run"]

        assert sorted(as_run.run_hours.unique()) == [18.124, 20.124, 24.124]


class TestTheNaiveAndCorrectedReadoutsDisagreeOnTheSignOfTheResponse:
    """The finding, in the form that survived a fourth plate: as a rank correlation with
    the dose, not as a ratio of two endpoints."""

    @pytest.mark.parametrize("truncation", ["as-run", MATCHED])
    def test_corrected_activity_ranks_with_the_dose_in_every_block(self, circuit, truncation):
        for (plate, strain), block in blocks(circuit, truncation).items():
            rho = spearmanr(block.galactose_percent, block.corrected_activity).statistic

            assert rho > 0.8, (plate, strain, rho)

    @pytest.mark.parametrize("truncation", ["as-run", MATCHED])
    def test_naive_rfu_per_od_never_does(self, circuit, truncation):
        """Zero or negative in 13 of the 14 blocks; the one exception, No-Gal on 08-02, is
        +0.14, which is one adjacent swap away from flat."""
        for (plate, strain), block in blocks(circuit, truncation).items():
            rho = spearmanr(block.galactose_percent, block.naive_rfu_per_od).statistic

            assert rho < 0.2, (plate, strain, rho)

    def test_the_july_plate_agrees_with_the_august_ones_on_this(self, circuit):
        """The point of scoring the fourth plate. It is a different experiment on a
        different week with a different comparator, and it gives the same answer."""
        july = blocks(circuit, MATCHED)[("20260705_AFL&control_1st", "Control")]

        assert spearmanr(july.galactose_percent, july.corrected_activity).statistic > 0.8
        assert spearmanr(july.galactose_percent, july.naive_rfu_per_od).statistic < 0.2

    def test_the_corrected_span_agrees_between_the_two_experiments(self, circuit, ):
        """Read at the same length, the July plate's two strains recover the same corrected
        response as the August plates' same-named strains, to within about 15%."""
        matched = circuit[circuit.truncation == MATCHED]

        def span(experiment, strain):
            block = matched[(matched.experiment == experiment) & (matched.strain == strain)]
            ladder = block.groupby("galactose_percent").corrected_activity.mean()
            return float(ladder.iloc[-1] - ladder.iloc[0])

        for strain in ("AFL", "Control"):
            july, august = span("AFL-vs-control", strain), span("AFL-debugging", strain)

            assert abs(july - august) / august < 0.15, (strain, july, august)


class TestTheEndpointRatioClaimIsFalseAndStaysBuried:
    """This file used to assert that naive `RFU/OD` FALLS from 0% to 1.8% galactose. It does
    on three plates and does not on the fourth, and on any one plate the ratio moves with
    how long the plate was read. The counterexamples are asserted so the claim cannot come
    back by accident."""

    def test_the_july_control_ratio_is_carried_by_one_dying_well(self, circuit):
        """This asserted that the July Control RISES where the August three fall, and it
        passed only because well G10 is in the mean.

        G10 is the culture `docs/research/AUTO_FEEDBACK_LOOP.md` describes as dying. Over
        the matched late window its blank-corrected OD is 1.107 against its five same-dose
        neighbours' 1.40-1.62, and its RFU/OD is 598.6 against their 241-281. Removing it
        takes the July Control endpoint ratio from 1.087 to 0.897 -- falling, on the same
        side of 1.0 as all three August Controls. Found by the 2026-08-30 audit.

        So what is pinned now is the true and weaker fact: the ratio is above 1.0 only
        marginally, and the three August Controls are all well below it. The claim that the
        two experiments DISAGREE on the naive readout is withdrawn.
        """
        ratios = {}
        for (plate, strain), block in blocks(circuit, MATCHED).items():
            if strain == "Control":
                ratios[plate] = (float(block.naive_rfu_per_od.iloc[-1])
                                 / float(block.naive_rfu_per_od.iloc[0]))

        july = ratios["20260705_AFL&control_1st"]
        august = [v for k, v in ratios.items() if k != "20260705_AFL&control_1st"]

        assert all(v < 1.0 for v in august)
        # Marginal, and one well from being below 1.0 itself.
        assert 1.0 < july < 1.15
        assert july - max(august) < 0.3

    def test_reading_longer_raises_the_naive_ratio_in_every_block_that_has_a_longer_read(
            self, circuit):
        """Two plates ran past the matched length: 20260802 to 20.1 h and 20260705 to 24.1 h.
        In all six of their strain blocks the extra hours push the naive endpoint ratio UP,
        because a stable reporter keeps accumulating in cultures the inducer slowed."""
        longer = [p for p in circuit.plate.unique() if p not in ALREADY_MATCHED]
        as_run, matched = blocks(circuit, "as-run"), blocks(circuit, MATCHED)
        compared = 0

        for plate, strain in itertools.product(longer, DEBUGGING_STRAINS):
            if (plate, strain) not in as_run:
                continue
            compared += 1
            ratio = [float(b[(plate, strain)].naive_rfu_per_od.iloc[-1])
                     / float(b[(plate, strain)].naive_rfu_per_od.iloc[0])
                     for b in (matched, as_run)]

            assert ratio[1] > ratio[0], (plate, strain, ratio)
        assert compared == 6

    def test_but_it_never_flips_the_sign_of_the_corrected_response(self, circuit):
        """Same six blocks. The extra hours shrink the corrected span -- activity decays as
        the cultures go stationary -- and never turn it negative."""
        longer = [p for p in circuit.plate.unique() if p not in ALREADY_MATCHED]
        as_run, matched = blocks(circuit, "as-run"), blocks(circuit, MATCHED)

        for plate, strain in itertools.product(longer, DEBUGGING_STRAINS):
            if (plate, strain) not in as_run:
                continue
            spans = [float(b[(plate, strain)].corrected_activity.iloc[-1])
                     - float(b[(plate, strain)].corrected_activity.iloc[0])
                     for b in (matched, as_run)]

            assert spans[0] > 0 and spans[1] > 0, (plate, strain, spans)

    def test_the_doc_says_the_claim_broke(self):
        doc = (pathlib.Path(__file__).resolve().parents[1]
               / "docs" / "research" / "AUTO_FEEDBACK_LOOP.md").read_text()

        assert "endpoint ratio" in doc.lower()
        assert "no-gal" in doc.lower()


class TestTheClosedLoopMovesLessThanTheOpenedOne:
    def test_afl_spans_less_than_no_laci(self, circuit):
        matched = circuit[(circuit.truncation == MATCHED)
                          & (circuit.experiment == "AFL-debugging")]

        def span(strain):
            ladder = matched[matched.strain == strain].groupby(
                "galactose_percent").corrected_activity.mean()
            return float(ladder.iloc[-1] - ladder.iloc[0])

        assert span("AFL") < span("No-LacI")
        assert span("AFL") < span("No-Gal")

    def test_but_no_claim_is_made_that_the_circuit_works(self):
        """Three of the four instrument files are named 'debugging'. Whether these numbers
        are the intended behaviour is the team's call, and the doc says so rather than the
        table implying it."""
        doc = (pathlib.Path(__file__).resolve().parents[1]
               / "docs" / "research" / "AUTO_FEEDBACK_LOOP.md").read_text()

        assert "Not that the circuit works" in doc

    def test_the_strain_without_a_gal_uas_responds_too(self, circuit):
        """A caveat the doc has to carry: No-Gal's promoter has no galactose UAS, and its
        corrected activity still rises further than any other strain's. So the rise is a
        response to the carbon source, not necessarily to GAL induction, and neither
        readout can separate the two."""
        matched = circuit[(circuit.truncation == MATCHED)
                          & (circuit.experiment == "AFL-debugging")]
        ladder = matched[matched.strain == "No-Gal"].groupby(
            "galactose_percent").corrected_activity.mean()

        assert ladder.iloc[-1] - ladder.iloc[0] > 0


class TestTheNegativeActivitiesAreReportedNotHidden:
    def test_some_zero_inducer_activities_are_negative(self, circuit):
        """A culture that has stopped growing with a slowly falling signal gives a negative
        derivative. Clipping them would hide the one place the estimator is visibly working
        at the edge of its assumptions."""
        zero = circuit[circuit.galactose_percent == 0.0]

        assert (zero.corrected_activity < 0).any()

    def test_and_the_highest_dose_is_positive_for_every_strain(self, circuit):
        top = circuit[(circuit.galactose_percent == 1.8)
                      & (circuit.truncation == MATCHED)].groupby(
            ["experiment", "strain"]).corrected_activity.mean()

        assert (top > 0).all()


class TestItReproducesFromTheInstrumentFiles:
    """Needs the .xpt files. Skips without them, like every real-data test here."""

    @pytest.fixture(scope="class")
    def directory(self):
        found = paths.gen5_xpt_dir()
        if found is None:
            pytest.skip("Gen5 .xpt files not present; set YSTWIN_GEN5_XPT")
        return found

    @pytest.fixture(scope="class")
    def script(self):
        import score_afl_circuit

        return score_afl_circuit

    def test_all_four_plates_share_one_time_grid_over_the_matched_window(self, directory, script):
        """What makes the matched comparison a truncation rather than an interpolation."""
        from ystwin.plate import gen5

        grids = []
        for series in script.SERIES:
            for name in series.plates:
                path = directory / name
                if not path.is_file():
                    pytest.skip(f"{name} not present")
                elapsed = gen5.read_xpt(path).channel("OD600:600").elapsed_ms
                grids.append(elapsed[:script.MATCHED_TIMEPOINTS])

        assert all(np.array_equal(grids[0], g) for g in grids)

    def test_the_sampling_interval_is_ten_minutes_whatever_the_protocol_is_called(
            self, directory, script):
        """`20260705_AFL&control_1st.xpt` names its protocol `24h-30min_...`, and its
        elapsed table says 10 minutes, the same as the August runs. The file wins."""
        from ystwin.plate import gen5

        path = directory / script.SERIES[0].plates[0]
        if not path.is_file():
            pytest.skip(f"{path.name} not present")
        run = gen5.read_xpt(path)
        minutes = np.diff(run.channel("OD600:600").elapsed_ms) / 60_000.0

        assert "30min" in run.protocol
        assert np.allclose(minutes, 10.0)

    @pytest.mark.parametrize("index", [0, 1])
    def test_the_blank_row_is_medium_not_air_on_both_maps(self, directory, script, index):
        """What confirms each notebook map against the readings: H1-H3 carry medium and
        read far above the empty wells, which is what a three-well blank looks like."""
        from ystwin.plate import gen5

        path = directory / script.SERIES[index].plates[0]
        if not path.is_file():
            pytest.skip(f"{path.name} not present")
        reporter = gen5.read_xpt(path).channel("mCitrine:480,530").frame()
        blank = reporter[list(script.BLANK_WELLS)].to_numpy()[-1].mean()
        empty = reporter[[f"H{c}" for c in range(5, 13)]].to_numpy()[-1].mean()

        assert blank > 2 * empty

    def test_the_galactose_ladder_shows_up_as_a_growth_delay_down_the_rows(
            self, directory, script):
        """The 2026-07-05 map has no strain header, so the row assignment is the only part
        of it that can be checked against the data -- and it can be: time to reach
        blank-corrected OD 1.0 is ~4 hours longer in row G than in row A, on both maps."""
        from ystwin.plate import gen5

        for series in script.SERIES:
            path = directory / series.plates[0]
            if not path.is_file():
                pytest.skip(f"{path.name} not present")
            run = gen5.read_xpt(path)
            optical = run.channel("OD600:600").frame()
            hours = run.channel("OD600:600").elapsed_hours()
            blank = optical[list(script.BLANK_WELLS)].to_numpy().mean(axis=1)

            def hours_to_od_one(row):
                out = []
                for column in range(1, 13):
                    biomass = optical[f"{row}{column}"].to_numpy() - blank
                    out.append(hours[np.argmax(biomass > 1.0)])
                return float(np.mean(out))

            assert hours_to_od_one("G") - hours_to_od_one("A") > 2.0, series.name

    def test_the_july_plate_splits_six_and_six_and_nowhere_else(self, directory, script):
        """The evidence behind `strain_source='inferred-from-signal'`. Columns 7-12 read
        above columns 1-6 at every dose, and no split inside either half does that."""
        from ystwin.plate import gen5

        path = directory / script.SERIES[0].plates[0]
        if not path.is_file():
            pytest.skip(f"{path.name} not present")
        run = gen5.read_xpt(path)
        reporter = run.channel("mCitrine:480,530").frame()
        optical = run.channel("OD600:600").frame()
        blank_rfu = reporter[list(script.BLANK_WELLS)].to_numpy()[-1].mean()
        blank_od = optical[list(script.BLANK_WELLS)].to_numpy()[-1].mean()

        def ratio(row, column):
            return ((reporter[f"{row}{column}"].iloc[-1] - blank_rfu)
                    / (optical[f"{row}{column}"].iloc[-1] - blank_od))

        for row in script.GALACTOSE_PERCENT:
            left = max(ratio(row, c) for c in range(1, 7))
            right = min(ratio(row, c) for c in range(7, 13) if not (row, c) == ("G", 10))

            assert right > left, row

    def test_rescoring_reproduces_the_committed_table(self, directory, script, circuit):
        for series in script.SERIES:
            for name in series.plates:
                path = directory / name
                if not path.is_file():
                    pytest.skip(f"{name} not present")
                fresh = script.score_plate(path, series)
                committed = circuit[(circuit.plate == path.stem)
                                    & (circuit.truncation == "as-run")]
                merged = fresh.merge(committed, on=["strain", "galactose_percent"],
                                     suffixes=("_new", "_old"))

                assert len(merged) == len(fresh)
                assert np.allclose(merged.corrected_activity_new, merged.corrected_activity_old)
                assert np.allclose(merged.naive_rfu_per_od_new, merged.naive_rfu_per_od_old)

    def test_a_run_shorter_than_the_matched_window_is_refused_not_padded(
            self, directory, script):
        path = directory / script.SERIES[0].plates[0]
        if not path.is_file():
            pytest.skip(f"{path.name} not present")

        with pytest.raises(ValueError, match="fewer than"):
            script.score_plate(path, script.SERIES[0], timepoints=1000)


#: The two `AFL-vs-control` instrument files that are NOT scored, exactly as named on disk.
UNPLACED = ("20260712_AFL&control_2nd_AND_oxidative_stress_1st.xpt",
            "20260713_AFL&Control3_ER Stress1.xpt")


def _late_ratio_grid(path, script, timepoints=None):
    """(row, column) -> late-window mean blank-corrected RFU/OD, rows A-G."""
    from ystwin.plate import gen5

    run = gen5.read_xpt(path)
    optical = run.channel("OD600:600").frame()
    reporter = run.channel("mCitrine:480,530").frame()
    if timepoints is not None:
        optical, reporter = optical.iloc[:timepoints], reporter.iloc[:timepoints]
    blank_od = optical[list(script.BLANK_WELLS)].to_numpy().mean(axis=1)
    blank_rfu = reporter[list(script.BLANK_WELLS)].to_numpy().mean(axis=1)
    late = slice(int(len(optical) * script.LATE_FRACTION), None)
    grid = {}
    for row in script.GALACTOSE_PERCENT:
        for column in range(1, 13):
            well = f"{row}{column}"
            biomass = optical[well].to_numpy() - blank_od
            signal = reporter[well].to_numpy() - blank_rfu
            grid[(row, column)] = float(np.mean((signal / biomass)[late]))
    return grid


def _separating_bipartitions(grid, columns, rows):
    """Every split of `columns` whose two sides' ranges do not overlap at ANY row.

    The recovery that placed 2026-07-05: a layout is established only when exactly one
    partition survives, because a plausible one is not evidence.
    """
    columns = list(columns)
    found = []
    for mask in range(1, 2 ** (len(columns) - 1)):
        left = [c for i, c in enumerate(columns) if (mask >> i) & 1]
        right = [c for c in columns if c not in left]
        if not left or not right:
            continue
        if all(min(grid[(r, c)] for c in left) > max(grid[(r, c)] for c in right)
               or min(grid[(r, c)] for c in right) > max(grid[(r, c)] for c in left)
               for r in rows):
            found.append((tuple(left), tuple(right)))
    return found


class TestTheTwoMidJulyPlatesStayUnplaced:
    """`20260712` and `20260713` are the other two `AFL-vs-control` files, and both were put
    through the recovery that placed `20260705`. It fails on both. Pinned so that adding
    either to `SERIES` has to argue with a test rather than with nobody.

    See `docs/research/AUTO_FEEDBACK_LOOP.md`, "2026-07-12 and 2026-07-13 — read, and
    deliberately not placed".
    """

    @pytest.fixture(scope="class")
    def directory(self):
        found = paths.gen5_xpt_dir()
        if found is None:
            pytest.skip("Gen5 .xpt files not present; set YSTWIN_GEN5_XPT")
        return found

    @pytest.fixture(scope="class")
    def script(self):
        import score_afl_circuit

        return score_afl_circuit

    def grid(self, directory, script, name):
        path = directory / name
        if not path.is_file():
            pytest.skip(f"{name} not present")
        return _late_ratio_grid(path, script, script.MATCHED_TIMEPOINTS)

    def test_the_recovery_still_returns_the_july_split_it_was_built_on(self, directory, script):
        """The control. On 2026-07-05 exactly one of the 2047 partitions separates, and it is
        the 6/6 one `SERIES` uses; without this the two negatives below mean nothing."""
        grid = self.grid(directory, script, script.SERIES[0].plates[0])

        assert _separating_bipartitions(grid, range(1, 13), script.GALACTOSE_PERCENT) == [
            ((1, 2, 3, 4, 5, 6), (7, 8, 9, 10, 11, 12))]

    def test_and_recovers_the_notebook_blocks_on_a_plate_that_has_a_map(self, directory, script):
        grid = self.grid(directory, script, script.SERIES[1].plates[0])

        assert ((1, 2, 3), (4, 5, 6)) in _separating_bipartitions(
            grid, range(1, 7), script.GALACTOSE_PERCENT)
        assert _separating_bipartitions(grid, range(7, 13), script.GALACTOSE_PERCENT) == [
            ((7, 8, 9), (10, 11, 12))]

    def test_20260712_carries_no_column_structure_at_all(self, directory, script):
        """Zero of 2047 over the whole plate and zero of 31 inside columns 1-6. There is
        nothing on this plate to adopt, whatever the file name says it holds."""
        grid = self.grid(directory, script, UNPLACED[0])

        assert _separating_bipartitions(grid, range(1, 13), script.GALACTOSE_PERCENT) == []
        assert _separating_bipartitions(grid, range(1, 7), script.GALACTOSE_PERCENT) == []

    def test_20260713_has_a_unique_split_that_is_still_not_adopted(self, directory, script):
        """It returns `{1,2,3}|{4,5,6}` and nothing else inside columns 1-6. The next test is
        why that is not enough."""
        grid = self.grid(directory, script, UNPLACED[1])

        assert _separating_bipartitions(grid, range(1, 7), script.GALACTOSE_PERCENT) == [
            ((1, 2, 3), (4, 5, 6))]

    def test_both_candidate_blocks_are_far_too_bright_to_be_afl_or_control(
            self, directory, script):
        """The criterion that fixed the 2026-07-05 assignment was that its two halves' ranges
        ARE the August AFL and Control ranges. Here the two populations do not overlap at all:
        every AFL and Control well on the four scored plates reads 70-599 late RFU/OD, and
        every well in columns 1-6 of the two unplaced plates reads 1,999-9,038 -- the dimmest
        candidate well is 3.3x the brightest AFL or Control well in the series. So the unique
        split on 20260713 cannot be those two strains, and neither plate is scored."""
        recorded = []
        for series in script.SERIES:
            for name in series.plates:
                grid = self.grid(directory, script, name)
                recorded += [grid[(r, c)] for r in script.GALACTOSE_PERCENT
                             for c in series.strain_columns["AFL"] + series.strain_columns["Control"]]
        brightest = max(recorded)

        for name in UNPLACED:
            grid = self.grid(directory, script, name)
            dimmest = min(grid[(r, c)] for r in script.GALACTOSE_PERCENT for c in range(1, 7))

            assert dimmest > 3 * brightest, (name, dimmest, brightest)

    def test_the_doc_says_why_they_are_not_scored(self):
        doc = (pathlib.Path(__file__).resolve().parents[1]
               / "docs" / "research" / "AUTO_FEEDBACK_LOOP.md").read_text()

        assert "deliberately not placed" in doc
        assert "20260713_AFL&Control3_ER Stress1.xpt" in doc


def _partial_dose_rho(frame):
    """Spearman(galactose, corrected) raw, and partialled on `mu_late`."""
    import numpy as np
    from scipy import stats

    xy = stats.spearmanr(frame.galactose_percent, frame.corrected_activity).statistic
    xz = stats.spearmanr(frame.galactose_percent, frame.mu_late).statistic
    yz = stats.spearmanr(frame.corrected_activity, frame.mu_late).statistic
    denominator = np.sqrt((1 - xz**2) * (1 - yz**2))
    return xy, ((xy - xz * yz) / denominator if denominator > 0 else float("nan"))

class TestTheGrowthRateConfoundIsVisibleRatherThanAbsent:
    """A galactose ladder is also a growth ladder, so the headline has a rival explanation.

    The 2026-08-30 audit found Spearman(corrected_activity, mu_late) at a median of 0.964
    across blocks while the only committed test asserted rho(galactose, corrected) > 0.8 --
    a threshold a pure growth-stage artefact passes identically. `mu_late` was not even a
    column, so a reader could not check. These tests pin the confound's size, not its
    absence.
    """

    def test_the_growth_rate_is_a_column(self, circuit):
        assert "mu_late" in circuit.columns
        assert circuit.mu_late.notna().all()

    def test_dose_and_growth_rate_are_strongly_associated(self, circuit):
        """The confound itself. If this were weak there would be nothing to condition on."""
        from scipy import stats

        rho = stats.spearmanr(circuit.galactose_percent, circuit.mu_late).statistic

        assert rho > 0.8

    def test_conditioning_on_growth_rate_halves_the_pooled_association(self, circuit):
        """Neither "it survives" nor "it vanishes". Both would be easier to write down."""
        raw, partial = _partial_dose_rho(circuit)

        assert raw > 0.8
        assert 0.35 < partial < 0.65
        assert partial < 0.7 * raw

    def test_but_within_a_block_it_survives(self, circuit):
        """Where the comparison is one strain on one plate, growth rate does not account
        for the ranking. Seven points per block, so this saturates -- which is why the
        doc quotes it as weaker than it looks rather than as the headline."""
        import numpy as np

        values = []
        for _, block in circuit.groupby(
                ["experiment", "plate", "truncation", "strain"]):
            if block.galactose_percent.nunique() < 4:
                continue
            _, partial = _partial_dose_rho(block)
            if np.isfinite(partial):
                values.append(partial)

        assert len(values) >= 12
        assert sum(v > 0 for v in values) >= 0.8 * len(values)

    def test_a_third_of_the_table_is_negative_and_not_only_at_zero_dose(self, circuit):
        """The doc said "at 0% galactose". It is four doses and 23 of 28 blocks."""
        negative = circuit[circuit.corrected_activity < 0]

        assert len(negative) > len(circuit) / 4
        assert negative.galactose_percent.nunique() >= 4
        assert (negative.galactose_percent > 0).any()
