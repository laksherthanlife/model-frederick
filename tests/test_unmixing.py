"""Linear unmixing of co-expressed reporters, and whether it beats the growth confound.

Co-expression changes the problem. Reporters in one cell share mu exactly, so a ratio
between two of them cancels the dilution term:

    R_i / R_j = [k_i/(mu+kdeg_i)] / [k_j/(mu+kdeg_j)]  ->  k_i/k_j   when kdeg matches

That is an exact cancellation, not a correction -- no growth estimate is needed at all.
The question is what breaks it: mismatched maturation and stability between fluorophores,
and error introduced by the unmixing itself.
"""

import numpy as np
import pytest

from ystwin.generator.unmixing import (
    PANEL,
    FluorophoreSpec,
    mixing_matrix,
    ratio_series,
    unmix,
)


class TestThePanel:
    def test_it_holds_the_five_planned_fluorophores(self):
        assert set(PANEL) == {"CFP", "eGFP", "YFP", "mOrange2", "mApple"}

    def test_each_carries_excitation_emission_and_kinetics(self):
        for name, spec in PANEL.items():
            assert spec.excitation_nm > 0 and spec.emission_nm > spec.excitation_nm, name
            assert spec.maturation_h > 0, name

    def test_they_are_ordered_across_the_spectrum(self):
        peaks = [PANEL[n].emission_nm for n in ("CFP", "eGFP", "YFP", "mOrange2", "mApple")]
        assert peaks == sorted(peaks)


class TestMixingAndUnmixing:
    def test_the_mixing_matrix_has_one_row_per_band_and_column_per_fluorophore(self):
        names = ["CFP", "eGFP", "YFP"]
        bands = [(460, 480), (500, 520), (525, 545)]

        matrix = mixing_matrix(names, bands)

        assert matrix.shape == (3, 3)

    def test_adjacent_fluorophores_genuinely_overlap(self):
        """If they did not, unmixing would be unnecessary."""
        matrix = mixing_matrix(["eGFP", "YFP"], [(500, 520), (525, 545)])

        assert matrix[0, 1] > 0.05

    def test_unmixing_recovers_known_abundances(self):
        names = ["CFP", "eGFP", "YFP", "mOrange2", "mApple"]
        bands = [(460, 480), (500, 520), (525, 545), (560, 580), (585, 610)]
        matrix = mixing_matrix(names, bands)
        truth = np.array([1.0, 2.0, 0.5, 1.5, 0.8])

        recovered = unmix(matrix @ truth, matrix)

        assert recovered == pytest.approx(truth, rel=1e-6)

    def test_detector_noise_sets_an_error_floor(self):
        names = ["CFP", "eGFP", "YFP"]
        bands = [(460, 480), (500, 520), (525, 545)]
        matrix = mixing_matrix(names, bands)
        truth = np.array([1.0, 2.0, 0.5])
        rng = np.random.default_rng(0)

        measured = matrix @ truth
        noisy = measured * (1 + rng.normal(0, 0.02, measured.shape))
        recovered = unmix(noisy, matrix)

        assert np.max(np.abs(recovered - truth) / truth) < 0.5

    def test_more_bands_than_fluorophores_is_accepted(self):
        names = ["CFP", "eGFP"]
        bands = [(460, 480), (500, 520), (525, 545), (470, 500)]
        matrix = mixing_matrix(names, bands)

        recovered = unmix(matrix @ np.array([1.0, 2.0]), matrix)

        assert recovered == pytest.approx([1.0, 2.0], rel=1e-6)

    def test_fewer_bands_than_fluorophores_is_refused(self):
        with pytest.raises(ValueError, match="at least as many bands"):
            mixing_matrix(["CFP", "eGFP", "YFP"], [(460, 480), (500, 520)])


class TestTheRatioCancelsGrowth:
    def test_matched_reporters_give_a_growth_invariant_ratio(self):
        """The key claim: no growth estimate needed at all."""
        fast = ratio_series(stress_activity=2.0e-3, reference_activity=1.0e-3,
                            growth_rate=0.40, stress_spec=PANEL["YFP"],
                            reference_spec=PANEL["YFP"])
        slow = ratio_series(stress_activity=2.0e-3, reference_activity=1.0e-3,
                            growth_rate=0.10, stress_spec=PANEL["YFP"],
                            reference_spec=PANEL["YFP"])

        assert fast == pytest.approx(slow, rel=1e-6)
        assert fast == pytest.approx(2.0, rel=1e-6)

    def test_mismatched_maturation_breaks_the_cancellation(self):
        slow_maturing = FluorophoreSpec("slow", 500, 520, maturation_h=2.0)
        fast_maturing = FluorophoreSpec("fast", 500, 520, maturation_h=0.1)

        fast = ratio_series(2.0e-3, 1.0e-3, growth_rate=0.40,
                            stress_spec=slow_maturing, reference_spec=fast_maturing)
        slow = ratio_series(2.0e-3, 1.0e-3, growth_rate=0.10,
                            stress_spec=slow_maturing, reference_spec=fast_maturing)

        assert abs(fast - slow) / slow > 0.05

    def test_the_error_from_mismatch_grows_with_the_difference(self):
        reference = FluorophoreSpec("ref", 500, 520, maturation_h=0.2)
        small = FluorophoreSpec("a", 500, 520, maturation_h=0.4)
        large = FluorophoreSpec("b", 500, 520, maturation_h=4.0)

        def spread(spec):
            fast = ratio_series(2e-3, 1e-3, 0.40, spec, reference)
            slow = ratio_series(2e-3, 1e-3, 0.10, spec, reference)
            return abs(fast - slow) / slow

        assert spread(large) > spread(small)

    def test_a_single_reporter_still_needs_a_growth_estimate(self):
        """Contrast: without a co-expressed reference the dilution term survives."""
        from ystwin.analysis.power import _reporter_at

        fast = _reporter_at(2.0e-3, 0.40, 0.0, 100.0)
        slow = _reporter_at(2.0e-3, 0.10, 0.0, 100.0)

        assert slow == pytest.approx(fast * 4, rel=1e-6)
