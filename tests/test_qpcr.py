"""Relative expression from Cq, and the QC that has to pass before it means anything.

The anchor for G4 is endogenous stress-gene expression: Hac1 for the ER reporters,
TRX2 for the oxidative ones, each against UBC as reference. Two things must hold
before a fold change is interpretable:

  * the reference gene must actually be measured for that sample, and
  * the no-reverse-transcriptase control must sit far above the +RT signal, or the
    "expression" being measured is partly genomic DNA.
"""

import numpy as np
import pandas as pd
import pytest

from ystwin.qpcr import (
    TARGET_FOR_CONSTRUCT,
    delta_delta_cq,
    rt_minus_margin,
)


def _tidy(rows):
    return pd.DataFrame(
        rows, columns=["construct", "target", "dose_mM", "tech_rep", "cq", "rt_minus_cq"]
    )


class TestDeltaDeltaCq:
    def test_one_cycle_of_delta_delta_is_a_two_fold_change(self):
        tidy = _tidy([
            ("UPRE1", "UBC", 0.0, 1, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.0, 1, 26.0, 32.0),
            ("UPRE1", "UBC", 0.5, 1, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.5, 1, 25.0, 32.0),
        ])

        out = delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)
        treated = out.set_index("dose_mM").loc[0.5]

        assert treated.delta_delta_cq == pytest.approx(-1.0)
        assert treated.fold_change == pytest.approx(2.0)

    def test_the_control_dose_is_exactly_one_fold_by_construction(self):
        tidy = _tidy([
            ("UPRE1", "UBC", 0.0, 1, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.0, 1, 26.0, 32.0),
            ("UPRE1", "UBC", 0.5, 1, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.5, 1, 24.0, 32.0),
        ])

        out = delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)

        assert out.set_index("dose_mM").loc[0.0].fold_change == pytest.approx(1.0)

    def test_a_shift_in_the_reference_gene_is_divided_out(self):
        """Doubling input material moves both genes and must not look like induction."""
        base = [
            ("UPRE1", "UBC", 0.0, 1, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.0, 1, 26.0, 32.0),
            ("UPRE1", "UBC", 0.5, 1, 19.0, 30.0),
            ("UPRE1", "Hac1", 0.5, 1, 25.0, 32.0),
        ]

        out = delta_delta_cq(_tidy(base), reference_target="UBC", control_dose=0.0)

        assert out.set_index("dose_mM").loc[0.5].fold_change == pytest.approx(1.0)

    def test_technical_replicates_are_averaged_before_the_ratio(self):
        tidy = _tidy([
            ("UPRE1", "UBC", 0.0, 1, 20.0, 30.0), ("UPRE1", "UBC", 0.0, 2, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.0, 1, 25.0, 32.0), ("UPRE1", "Hac1", 0.0, 2, 27.0, 32.0),
            ("UPRE1", "UBC", 0.5, 1, 20.0, 30.0), ("UPRE1", "UBC", 0.5, 2, 20.0, 30.0),
            ("UPRE1", "Hac1", 0.5, 1, 24.0, 32.0), ("UPRE1", "Hac1", 0.5, 2, 26.0, 32.0),
        ])

        out = delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)

        assert out.set_index("dose_mM").loc[0.5].fold_change == pytest.approx(2.0)

    def test_amplification_efficiency_below_one_is_accounted_for(self):
        tidy = _tidy([
            ("UPRE1", "UBC", 0.0, 1, 20.0, 30.0), ("UPRE1", "Hac1", 0.0, 1, 26.0, 32.0),
            ("UPRE1", "UBC", 0.5, 1, 20.0, 30.0), ("UPRE1", "Hac1", 0.5, 1, 25.0, 32.0),
        ])

        perfect = delta_delta_cq(tidy, "UBC", 0.0, efficiency=1.0)
        poor = delta_delta_cq(tidy, "UBC", 0.0, efficiency=0.8)

        assert poor.set_index("dose_mM").loc[0.5].fold_change < \
               perfect.set_index("dose_mM").loc[0.5].fold_change

    def test_each_construct_is_normalised_against_its_own_reference(self):
        tidy = _tidy([
            ("UPRE1", "UBC", 0.0, 1, 20.0, 30.0), ("UPRE1", "Hac1", 0.0, 1, 26.0, 32.0),
            ("UPRE1", "UBC", 0.5, 1, 20.0, 30.0), ("UPRE1", "Hac1", 0.5, 1, 25.0, 32.0),
            ("NativeYap1", "UBC", 0.0, 1, 22.0, 30.0), ("NativeYap1", "TRX2", 0.0, 1, 15.0, 24.0),
            ("NativeYap1", "UBC", 0.5, 1, 22.0, 30.0), ("NativeYap1", "TRX2", 0.5, 1, 13.0, 24.0),
        ])

        out = delta_delta_cq(tidy, "UBC", 0.0).set_index(["construct", "dose_mM"])

        assert out.loc[("UPRE1", 0.5)].fold_change == pytest.approx(2.0)
        assert out.loc[("NativeYap1", 0.5)].fold_change == pytest.approx(4.0)

    def test_a_sample_with_no_reference_measurement_is_refused_by_name(self):
        tidy = _tidy([
            ("UPRE1", "Hac1", 0.0, 1, 26.0, 32.0),
            ("UPRE1", "Hac1", 0.5, 1, 25.0, 32.0),
        ])

        with pytest.raises(ValueError, match="UBC"):
            delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)

    def test_a_missing_control_dose_is_refused(self):
        tidy = _tidy([
            ("UPRE1", "UBC", 0.5, 1, 20.0, 30.0), ("UPRE1", "Hac1", 0.5, 1, 25.0, 32.0),
        ])

        with pytest.raises(ValueError, match="control dose"):
            delta_delta_cq(tidy, reference_target="UBC", control_dose=0.0)


class TestNoReverseTranscriptaseControl:
    def test_a_clean_sample_shows_a_wide_margin(self):
        tidy = _tidy([("UPRE1", "Hac1", 0.0, 1, 26.0, 34.0)])

        out = rt_minus_margin(tidy)

        assert out.iloc[0].margin_cycles == pytest.approx(8.0)
        assert bool(out.iloc[0].passed)

    def test_a_narrow_margin_is_flagged_as_possible_genomic_contamination(self):
        tidy = _tidy([("UPRE1", "Hac1", 0.0, 1, 26.0, 27.5)])

        out = rt_minus_margin(tidy)

        assert not bool(out.iloc[0].passed)

    def test_the_threshold_is_adjustable(self):
        tidy = _tidy([("UPRE1", "Hac1", 0.0, 1, 26.0, 30.0)])

        assert bool(rt_minus_margin(tidy, min_cycles=3.0).iloc[0].passed)
        assert not bool(rt_minus_margin(tidy, min_cycles=5.0).iloc[0].passed)

    def test_a_missing_rt_minus_reading_is_not_silently_passed(self):
        tidy = _tidy([("UPRE1", "Hac1", 0.0, 1, 26.0, np.nan)])

        out = rt_minus_margin(tidy)

        assert not bool(out.iloc[0].passed)


def test_each_construct_is_mapped_to_the_stress_gene_it_reports_on():
    assert TARGET_FOR_CONSTRUCT["UPRE1"] == "Hac1"
    assert TARGET_FOR_CONSTRUCT["UPRE2"] == "Hac1"
    assert TARGET_FOR_CONSTRUCT["NativeYap1"] == "TRX2"
    assert TARGET_FOR_CONSTRUCT["AlteredYap1"] == "TRX2"


class TestDoseParsing:
    def test_an_empty_cell_is_not_a_dose(self):
        """NaN is a float; without an explicit check it matches a blank column."""
        from ystwin.qpcr import _parse_dose

        assert _parse_dose(np.nan) is None
        assert _parse_dose(None) is None

    def test_a_labelled_dose_is_read_as_a_number(self):
        from ystwin.qpcr import _parse_dose

        assert _parse_dose("0.2 mM") == pytest.approx(0.2)
        assert _parse_dose("0 mM") == pytest.approx(0.0)
        assert _parse_dose("5 MM") == pytest.approx(5.0)

    def test_a_bare_number_is_accepted(self):
        from ystwin.qpcr import _parse_dose

        assert _parse_dose(0.5) == pytest.approx(0.5)

    def test_unrelated_text_is_not_a_dose(self):
        from ystwin.qpcr import _parse_dose

        assert _parse_dose("UBC") is None
        assert _parse_dose("Technical Replicate 1") is None


class TestStressorMapping:
    """Two constructs sense ER stress and two sense oxidative stress, with different
    stressors. Doses are therefore not comparable across the pair boundary, and a
    plot or table that labels them all 'DTT' is wrong."""

    def test_the_er_reporters_are_challenged_with_dtt(self):
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

        assert STRESSOR_FOR_CONSTRUCT["UPRE1"] == "DTT"
        assert STRESSOR_FOR_CONSTRUCT["UPRE2"] == "DTT"

    def test_the_oxidative_reporters_are_challenged_with_hydrogen_peroxide(self):
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT

        assert STRESSOR_FOR_CONSTRUCT["NativeYap1"] == "H2O2"
        assert STRESSOR_FOR_CONSTRUCT["AlteredYap1"] == "H2O2"

    def test_every_construct_with_an_anchor_also_has_a_declared_stressor(self):
        from ystwin.qpcr import STRESSOR_FOR_CONSTRUCT, TARGET_FOR_CONSTRUCT

        assert set(STRESSOR_FOR_CONSTRUCT) == set(TARGET_FOR_CONSTRUCT)

    def test_doses_are_not_pooled_across_different_stressors(self):
        from ystwin.qpcr import assert_single_stressor

        with pytest.raises(ValueError, match="different stressors"):
            assert_single_stressor(["UPRE1", "NativeYap1"])

    def test_constructs_sharing_a_stressor_may_be_pooled(self):
        from ystwin.qpcr import assert_single_stressor

        assert assert_single_stressor(["UPRE1", "UPRE2"]) == "DTT"


class TestGdnaShare:
    """The RT- margin in the unit that decides whether a reading is a measurement.

    Cycles are what the instrument prints; contamination fraction is what determines
    whether any analysis can recover transcript from the reading. The conversion is
    exact, and the point of having it is that the cycle scale misleads badly near zero.
    """

    def test_converts_the_observed_margins_to_the_published_shares(self):
        from ystwin.qpcr import gdna_share

        # The three margins this project actually measured, and the convention.
        assert gdna_share(8.6) == pytest.approx(0.0026, abs=5e-4)   # TRX2, clean
        assert gdna_share(5.0) == pytest.approx(0.0312, abs=1e-3)   # the convention
        assert gdna_share(0.60) == pytest.approx(0.660, abs=1e-3)   # UPRE2
        assert gdna_share(0.13) == pytest.approx(0.914, abs=1e-3)   # UPRE1

    def test_the_correctability_ceiling_sits_near_three_quarters_of_a_cycle(self):
        """ValidPrime's 60% limit corresponds to about 0.74 cycles."""
        from ystwin.qpcr import gdna_share

        assert gdna_share(0.74) == pytest.approx(0.60, abs=0.01)

    def test_a_zero_margin_is_all_contamination(self):
        from ystwin.qpcr import gdna_share

        assert gdna_share(0.0) == pytest.approx(1.0)

    def test_share_falls_monotonically_with_margin(self):
        from ystwin.qpcr import gdna_share

        shares = gdna_share(np.array([0.0, 0.5, 1.0, 2.0, 5.0, 10.0]))
        assert np.all(np.diff(shares) < 0)
        assert np.all((shares >= 0) & (shares <= 1))

    def test_lower_efficiency_means_a_margin_buys_less(self):
        """At efficiency below one a cycle is less than a doubling, so the same
        margin leaves more contamination behind."""
        from ystwin.qpcr import gdna_share

        assert gdna_share(3.0, efficiency=0.9) > gdna_share(3.0, efficiency=1.0)

    def test_missing_margin_propagates_rather_than_defaulting(self):
        from ystwin.qpcr import gdna_share

        assert np.isnan(gdna_share(np.array([np.nan]))[0])

    def test_non_positive_efficiency_is_refused(self):
        from ystwin.qpcr import gdna_share

        with pytest.raises(ValueError, match="efficiency must be positive"):
            gdna_share(3.0, efficiency=0.0)


class TestRtMinusCorrectability:
    """Failing QC and being beyond correction are different findings."""

    def _frame(self, margin):
        return pd.DataFrame({"cq": [25.0], "rt_minus_cq": [25.0 + margin]})

    def test_a_clean_reading_passes_and_is_correctable(self):
        from ystwin.qpcr import rt_minus_margin

        row = rt_minus_margin(self._frame(8.6)).iloc[0]
        assert row.passed and row.correctable

    def test_a_marginal_reading_fails_qc_but_stays_correctable(self):
        """Below the QC threshold, above the ValidPrime ceiling -- recoverable."""
        from ystwin.qpcr import rt_minus_margin

        row = rt_minus_margin(self._frame(1.5)).iloc[0]
        assert not row.passed
        assert row.correctable
        assert row.gdna_share == pytest.approx(0.354, abs=1e-3)

    def test_the_observed_er_margins_are_past_correction(self):
        from ystwin.qpcr import rt_minus_margin

        for margin in (0.13, 0.60):
            row = rt_minus_margin(self._frame(margin)).iloc[0]
            assert not row.passed
            assert not row.correctable, f"{margin} cycles should be uncorrectable"

    def test_a_missing_control_is_neither_passed_nor_correctable(self):
        """An unknown contamination level is not a correctable one."""
        from ystwin.qpcr import rt_minus_margin

        frame = pd.DataFrame({"cq": [25.0], "rt_minus_cq": [np.nan]})
        row = rt_minus_margin(frame).iloc[0]
        assert not row.passed
        assert not row.correctable
