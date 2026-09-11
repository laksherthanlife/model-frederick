"""The gain-100 channel nobody had read, and what it does and does not license.

`outputs/gain_linearity.csv` is tracked, so the assertions about the RESULT run on any
clone. The assertions that need the instrument files skip without them, the same contract
every real-data fixture in this repository honours.

Three things are being pinned, and the third is the one most likely to be undone by a
well-meaning future edit:

1. **The fitter can find a breakpoint.** A linearity check that never fails is not
   evidence, so the synthetic cases put a known compression in and require it back at the
   right place. Without that, "no saturation observed" on the real plates would be
   unfalsifiable.
2. **The real plates show no breakpoint** -- the reported RFU scale is a fixed multiple of
   itself from the fitting window to 99,998, and it fails by writing an overflow flag
   rather than by compressing.
3. **This does not touch `od_linear_max`.** It is a photomultiplier; `od_linear_max` is a
   photodiode. A test holds the placeholder at 1.0 so that a later reader who finds a
   measured linear range in this repository cannot quietly wire it into the wrong gate.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

from ystwin import paths
from ystwin.calib.gain_linearity import (RFU_CEILING, RFU_GAIN_EXPONENT,
                                         RFU_LINEAR_MAX_MEASURED, RFU_OVERFLOW,
                                         MCITRINE_GAIN_RATIO_75_TO_100,
                                         fit_gain_pair, is_overflow)
from ystwin.gates.g1_optical import OpticalQualityGate

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

RATIO = 8.15
"""The multiple the synthetic cases are built with, near the measured one."""


def _ladder(n: int = 4000, top: float = 12_000.0) -> np.ndarray:
    """A low-gain channel spanning four decades, deterministic and without ties."""
    return np.geomspace(1.0, top, n)


def _companion(low: np.ndarray, ratio: float = RATIO, compress_above: float | None = None,
               noise: float = 0.0) -> np.ndarray:
    """The high-gain partner, optionally compressing above a stated high-channel level."""
    high = ratio * low
    if compress_above is not None:
        over = high > compress_above
        high = np.where(over, compress_above + 0.25 * (high - compress_above), high)
    if noise:
        high = high * (1.0 + noise * np.cos(np.arange(low.size)))
    return high


class TestTheFitterCanActuallyFindABreakpoint:
    """A check that cannot fail reads as assurance and is worse than none."""

    def test_a_clean_pair_recovers_the_multiple(self):
        low = _ladder()

        fit = fit_gain_pair(low, _companion(low))

        assert fit.ratio == pytest.approx(RATIO, rel=1e-6)
        assert fit.offset == pytest.approx(0.0, abs=1e-6)

    def test_a_clean_pair_reports_no_saturation(self):
        low = _ladder()

        fit = fit_gain_pair(low, _companion(low))

        assert not fit.saturation_observed
        assert fit.worst_departure < 1e-6

    def test_a_compressing_channel_is_caught(self):
        low = _ladder(top=40_000.0)

        fit = fit_gain_pair(low, _companion(low, compress_above=50_000.0))

        assert fit.saturation_observed

    def test_the_breakpoint_lands_where_the_compression_was_put(self):
        """Not to the RFU -- to the bin. The fit reports the top of the last bin that held,
        so a compression starting at 50,000 must be called somewhere below the next
        boundary above it and not two decades away."""
        low = _ladder(top=40_000.0)

        fit = fit_gain_pair(low, _companion(low, compress_above=50_000.0))

        assert 40_000.0 <= fit.linear_max <= 70_000.0

    def test_a_gentler_tolerance_lets_more_compression_pass(self):
        low = _ladder(top=40_000.0)
        high = _companion(low, compress_above=50_000.0)

        strict = fit_gain_pair(low, high, tolerance=0.02).linear_max
        loose = fit_gain_pair(low, high, tolerance=0.30).linear_max

        assert loose > strict

    def test_the_linear_range_never_exceeds_the_largest_reading(self):
        """The doctrine `od_linearity.py` states: never extrapolated past the densest well
        actually run."""
        low = _ladder()
        high = _companion(low)

        fit = fit_gain_pair(low, high)

        assert fit.linear_max <= high.max()

    def test_a_sparse_top_decade_does_not_extend_the_linear_range(self):
        """Nine readings in the top bin cannot show a departure, so they must not be
        allowed to certify one either."""
        low = np.concatenate([np.geomspace(250.0, 2_500.0, 4000), np.full(9, 12_000.0)])
        high = _companion(low)

        fit = fit_gain_pair(low, high)

        assert fit.linear_max < 0.9 * high.max()


class TestOverflowIsAFlagAndNotAReading:
    def test_is_overflow_finds_the_sentinel(self):
        values = np.array([0.0, 42.0, 99_998.0, RFU_OVERFLOW])

        assert is_overflow(values).tolist() == [False, False, False, True]

    def test_a_real_reading_at_the_ceiling_is_not_flagged(self):
        assert not is_overflow(np.array([RFU_CEILING - 2.0])).any()

    def test_overflowed_readings_are_counted_and_excluded(self):
        low = _ladder()
        high = _companion(low)
        high[-100:] = RFU_OVERFLOW

        fit = fit_gain_pair(low, high)

        assert fit.n_overflow == 100
        assert fit.n_pairs == low.size - 100

    def test_an_overflow_does_not_drag_the_multiple_down(self):
        """The failure this guards: -99999 averaged in as if it were RFU. A ratio fitted
        over a channel with sentinels in it would come out negative or absurd."""
        low = _ladder()
        clean = fit_gain_pair(low, _companion(low)).ratio
        spiked = _companion(low)
        spiked[-500:] = RFU_OVERFLOW

        assert fit_gain_pair(low, spiked).ratio == pytest.approx(clean, rel=1e-6)


class TestItRefusesRatherThanGuesses:
    def test_mismatched_channels_raise(self):
        with pytest.raises(ValueError, match="same readings"):
            fit_gain_pair(np.ones(10), np.ones(11))

    def test_a_plate_that_never_reached_the_anchor_raises(self):
        low = np.geomspace(1.0, 100.0, 500)

        with pytest.raises(ValueError, match="anchor window"):
            fit_gain_pair(low, _companion(low))

    def test_a_descending_anchor_raises(self):
        low = _ladder()

        with pytest.raises(ValueError, match="ascending"):
            fit_gain_pair(low, _companion(low), anchor=(20_000.0, 2_000.0))

    def test_scatter_grows_with_noise(self):
        low = _ladder()
        quiet = fit_gain_pair(low, _companion(low, noise=0.005)).scatter
        loud = fit_gain_pair(low, _companion(low, noise=0.05)).scatter

        assert loud > 5 * quiet


@pytest.fixture(scope="module")
def committed():
    """The tracked result table. Skips only where the table has not been produced."""
    path = paths.outputs_dir() / "gain_linearity.csv"
    if not path.exists():
        pytest.skip(f"{path} not present; run scripts/score_gain_linearity.py")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def scored(committed):
    return committed[committed.status == "scored"]


class TestTheCommittedTableSaysWhatTheModuleClaims:
    def test_the_declared_gain_pairs_all_give_the_same_multiple(self, scored):
        ratios = scored[scored.gain_low == 75].plate_ratio.unique()

        assert len(ratios) >= 9
        assert ratios.min() > 8.0 and ratios.max() < 8.3

    def test_the_constant_is_the_mean_of_those_plates(self, scored):
        ratios = scored[scored.gain_low == 75].plate_ratio.unique()

        assert ratios.mean() == pytest.approx(MCITRINE_GAIN_RATIO_75_TO_100, abs=0.01)

    def test_the_gain_step_below_gives_the_same_exponent(self, scored):
        """The 50 -> 75 pair is a different step on a different part of the scale, and it
        has to land on the same power law or the ratio is a property of one plate."""
        low_step = scored[scored.gain_low == 50].plate_ratio.unique()

        assert len(low_step) == 1
        assert np.log(low_step[0]) / np.log(75 / 50) == pytest.approx(
            RFU_GAIN_EXPONENT, abs=0.05)

    def test_no_plate_shows_the_detector_compressing(self, scored):
        assert not scored.plate_saturation_observed.any()

    def test_the_highest_populated_bins_still_sit_on_the_line(self, scored):
        """The claim in one assertion: at 90,000 RFU and above the median of
        measured/predicted is still 1, so there is no roll-off before the ceiling."""
        top = scored[(scored.bin_low_rfu >= 90_000)
                     & (scored.bin_low_rfu < RFU_CEILING)
                     & (scored.n_pairs >= 30)]

        assert len(top) >= 3
        assert (top.departure_median - 1.0).abs().max() < 0.01

    def test_nothing_was_ever_reported_between_the_top_reading_and_the_ceiling(self, scored):
        readings = scored.high_max_rfu[scored.high_max_rfu > 0]

        assert readings.max() == RFU_LINEAR_MAX_MEASURED
        assert not ((readings > RFU_LINEAR_MAX_MEASURED) & (readings < RFU_CEILING)).any()

    def test_the_overflow_rows_carry_the_flag_and_not_a_value(self, committed):
        overflow = committed[(committed.bin_low_rfu == RFU_CEILING)
                             & (committed.n_pairs > 0)]

        assert len(overflow) >= 5
        assert (overflow.high_max_rfu == RFU_OVERFLOW).all()

    def test_the_afl_plates_are_marked_and_none_of_them_overflowed(self, scored):
        afl = scored[scored.is_afl_plate]

        assert afl.plate.nunique() == 4
        assert afl.plate_n_overflow.max() == 0

    def test_a_pair_the_data_could_not_answer_is_named_rather_than_dropped(self, committed):
        """`24h_30min-interval_mCitrine` is 100% overflow in both channels. A table that
        silently omitted it would look like a corpus with no such plate in it."""
        refused = committed[committed.status != "scored"]

        assert len(refused) >= 1
        assert refused.status.str.contains("anchor window").all()


class TestThisIsTheWrongDetectorForOdLinearMax:
    """The refusal, held in place by a test.

    Fluorescence is a photomultiplier counting emitted photons; absorbance is a photodiode
    measuring transmission against a reference, and its nonlinearity is stray light and
    pathlength. A linear range measured on one is not a linear range on the other, and the
    number below must not drift towards a fluorescence figure because one now exists.
    """

    def test_od_linear_max_is_still_the_placeholder(self):
        assert OpticalQualityGate().od_linear_max == 1.0

    def test_no_fluorescence_constant_leaked_into_the_optical_gate(self):
        gate = OpticalQualityGate()

        assert gate.od_linear_max not in (RFU_CEILING, RFU_LINEAR_MAX_MEASURED,
                                          MCITRINE_GAIN_RATIO_75_TO_100)

    def test_the_module_says_so_in_writing(self):
        from ystwin.calib import gain_linearity

        assert "od_linear_max" in gain_linearity.__doc__
        assert "different detector" in gain_linearity.__doc__


@pytest.fixture(scope="module")
def xpt_dir():
    found = paths.gen5_xpt_dir()
    if found is None:
        pytest.skip("Gen5 .xpt files not present; set YSTWIN_GEN5_XPT")
    return found


@pytest.mark.integration
class TestAgainstTheInstrumentFiles:
    def test_the_two_gains_are_52_seconds_apart_every_time(self, xpt_dir):
        """Not simultaneous, and the caveat in the module docstring depends on the figure.
        Gen5 sweeps the plate at one gain and then again at the next."""
        import score_gain_linearity as script
        from ystwin.plate import gen5

        for name in script.AFL_PLATES[1:]:
            path = xpt_dir / name
            if not path.is_file():
                pytest.skip(f"{name} not present")
            read = gen5.read_xpt(path)
            low, high = script.gain_pairs(read)[0]
            offsets = (high.elapsed_ms - low.elapsed_ms) / 1000.0

            assert set(offsets.tolist()) == {52.0}

    def test_the_reporter_channel_the_afl_analysis_reads_never_overflowed(self, xpt_dir):
        """The number that matters: a saturated reading is a silently wrong one, and there
        are none in the channel `scripts/score_afl_circuit.py` scores."""
        import score_gain_linearity as script
        from ystwin.plate import gen5

        for name in script.AFL_PLATES:
            path = xpt_dir / name
            if not path.is_file():
                pytest.skip(f"{name} not present")
            values = gen5.read_xpt(path).channel("mCitrine:480,530").values

            assert not is_overflow(values).any()
            assert values.max() < 0.1 * RFU_CEILING

    def test_no_committed_plate_holds_a_saturated_reading(self):
        """`data/plates` is the biosensor pipeline's input, and it is committed text rather
        than an instrument file, so this one runs on any clone."""
        directory = paths.REPO_ROOT / "data" / "plates"
        tables = sorted(directory.glob("*mCitrine*.csv"))
        if not tables:
            pytest.skip("no committed mCitrine plates")

        highest = -np.inf
        for path in tables:
            block = pd.read_csv(path).drop(columns=["elapsed_hms", "temperature_c"],
                                           errors="ignore").to_numpy(dtype=float)
            assert not is_overflow(block).any()
            highest = max(highest, float(np.nanmax(block)))

        assert highest < 0.5 * RFU_CEILING

    def test_rescoring_reproduces_the_committed_table(self, xpt_dir, committed):
        import score_gain_linearity as script

        name = "20260803_AFL-debugging_2nd.xpt"
        path = xpt_dir / name
        if not path.is_file():
            pytest.skip(f"{name} not present")
        fresh = script.score_plate(path)
        old = committed[committed.plate == path.stem]

        assert len(fresh) == len(old)
        assert np.allclose(fresh.plate_ratio, old.plate_ratio)
        assert np.allclose(fresh.departure_median.to_numpy(dtype=float),
                           old.departure_median.to_numpy(dtype=float), equal_nan=True)

    def test_the_overflow_flag_is_not_raised_early(self, xpt_dir):
        """What separates "the top of the scale" from "the detector gave up": ask the
        surviving partner what each flagged reading would have been. If the instrument were
        quitting below its stated range these would come in well under the ceiling."""
        import score_gain_linearity as script
        from ystwin.plate import gen5

        weakest, seen = np.inf, 0
        for name, gains in script.DECLARED_GAINS.items():
            path = xpt_dir / name
            if not path.is_file() or gains[-2:] != (75, 100):
                continue
            read = gen5.read_xpt(path)
            low, high = script._pairs(read, *[c.name for c in script.gain_pairs(read)[-1]])
            flagged = is_overflow(high) & ~is_overflow(low)
            if not flagged.any():
                continue
            seen += 1
            fit = fit_gain_pair(low, high, anchor=script.ANCHOR)
            weakest = min(weakest, float((fit.ratio * low[flagged] + fit.offset).min()))
        if not seen:
            pytest.skip("no overflowing plate present")

        assert seen >= 5
        assert weakest > 0.98 * RFU_CEILING

    def test_the_absorbance_channel_shows_no_top_of_scale(self, xpt_dir):
        """The negative that keeps `od_linear_max` where it is: OD600 does not clip in this
        corpus, so there is nothing here to measure even a *ceiling* against, let alone a
        linear range."""
        import score_gain_linearity as script

        census = script.absorbance_top_of_scale(xpt_dir)

        assert census["plates"] >= 15
        assert census["distinct_plate_maxima"] == census["plates"]
        assert census["most_readings_at_a_plate_max"] <= 5
