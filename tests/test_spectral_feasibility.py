"""A build has to be spectrally buildable, not merely informative.

Charging promoter fusions for the growth correction they need makes ratiometric sensors
strictly better on information, and an unconstrained search then fills all five channels
with them. That build cannot be made. roGFP2, HyPer7, QUEEN, pHluorin and Peredox are all
GFP- or cpYFP-derived excitation-ratio sensors emitting around 510 nm, so they occupy one
another's channel exactly; two of them in one strain are not two measurements.

Promoter fusions carry no such constraint, because the promoter is separable from the
protein -- any fusion can be built onto whichever fluorophore is still free. So the real
design question is a mix, and the thing that decides it is the spectrum, not the score.
"""

import pytest

from ystwin.analysis.design import NoiseModel
from ystwin.analysis.sensor_selection import select_sensors, spectral_conflict
from ystwin.generator.stress_panel import REPORTERS, Kind, ratiometric_reporters

PLATE = NoiseModel(relative_cv=0.05)


class TestSlots:
    def test_every_ratiometric_sensor_claims_a_spectral_slot(self):
        assert all(REPORTERS[r].spectral_slot for r in ratiometric_reporters())

    def test_a_promoter_fusion_claims_none_because_it_can_use_any_fluorophore(self):
        assert REPORTERS["STRE-general"].spectral_slot is None

    def test_the_green_ratio_sensors_all_claim_the_same_slot(self):
        slots = {REPORTERS[r].spectral_slot for r in ("roGFP2-Grx1", "HyPer7", "QUEEN-2m")}

        assert len(slots) == 1


class TestConflictDetection:
    def test_two_sensors_in_one_slot_conflict(self):
        assert spectral_conflict(["roGFP2-Grx1", "HyPer7"])

    def test_one_sensor_beside_promoter_fusions_does_not(self):
        assert not spectral_conflict(["roGFP2-Grx1", "STRE-general", "UPRE-ER"])

    def test_promoter_fusions_never_conflict_with_each_other(self):
        assert not spectral_conflict(["STRE-general", "UPRE-ER", "HSE-heat", "FeRE-iron"])

    def test_an_empty_build_does_not_conflict(self):
        assert not spectral_conflict([])


class TestSelectionRespectsTheSpectrum:
    def test_it_never_returns_a_build_that_cannot_be_made(self):
        chosen = select_sensors(n_channels=5, noise=PLATE, min_effect=0.5, growth_cv=0.30)

        assert not spectral_conflict(chosen)

    def test_it_still_takes_the_one_ratiometric_sensor_it_can_have(self):
        """The growth penalty is real, so the free green slot should not go unused."""
        chosen = select_sensors(n_channels=5, noise=PLATE, min_effect=0.5, growth_cv=0.30)

        assert sum(REPORTERS[r].kind is Kind.RATIOMETRIC for r in chosen) == 1

    def test_the_rest_of_the_build_is_promoter_fusions(self):
        chosen = select_sensors(n_channels=5, noise=PLATE, min_effect=0.5, growth_cv=0.30)

        assert sum(REPORTERS[r].kind is Kind.TRANSCRIPTIONAL for r in chosen) == 4

    def test_without_a_growth_penalty_it_is_still_buildable(self):
        chosen = select_sensors(n_channels=5, noise=PLATE, min_effect=0.5)

        assert not spectral_conflict(chosen)

    def test_a_required_pair_that_cannot_be_built_is_refused(self):
        with pytest.raises(ValueError, match="same spectral slot"):
            select_sensors(n_channels=5, noise=PLATE, min_effect=0.5,
                           required=["roGFP2-Grx1", "HyPer7"])
