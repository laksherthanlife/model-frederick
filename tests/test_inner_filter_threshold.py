"""When the reporter's inner-filter correction is worth measuring, and when it is not.

`ReporterOptics.inner_filter_coeff` refuses until it has been calibrated by spiking
purified carotenoid into a non-producing strain -- one plate. This asks whether that plate
is worth spending, and at these titres the answer is no: the attenuation is smaller than
the well-to-well noise it would be measured against.

That is a result, not a shortcut. `inner_filter_coeff_from_extinction` computes the
coupling from the pigment's own absorbance, so the threshold can be derived instead of
guessed, and the refusal stays in place for the strains that pass it.
"""

import numpy as np
import pytest

from ystwin.analysis.power import DEFAULT_WELL_CV
from ystwin.observation import (
    BETA_CAROTENE_MW_G_PER_MOL,
    BETA_CAROTENE_PETROLEUM_ETHER_450 as PEAK,
    PigmentExtinction,
    inner_filter_coeff_from_extinction,
)

# A well as this lab reads it: 200 uL in a 96-well plate is about 0.5 cm, and the
# packaging factor is the one PigmentOptics documents.
PATH_CM, BIOMASS_G_L, FLATTENING = 0.5, 1.0, 0.6

# Beta-carotene's absorbance at mCitrine's 516 nm excitation, as a fraction of its 450 nm
# peak. ASSERTED as a bracket, not measured: the spectrum falls steeply across that gap
# and no coefficient at 516 nm is vendored here. Everything below is reported across the
# whole bracket precisely so no single value inside it is load-bearing.
SPECTRAL_RATIO_BRACKET = (0.05, 0.25)

# Elizondo & Saa 2025 (PMID 40891387), the three calibration strains.
ELIZONDO_CONTENT_MG_PER_GDCW = (0.39, 0.97)


def _coeff(ratio: float) -> float:
    extinction = PigmentExtinction(
        epsilon_per_M_per_cm=PEAK.epsilon_per_M_per_cm * ratio,
        wavelength_nm=516.0, solvent=PEAK.solvent,
        source=f"asserted bracket: {ratio:.0%} of the 450 nm peak",
    )
    return inner_filter_coeff_from_extinction(
        extinction, 516.0, PATH_CM, BIOMASS_G_L, FLATTENING)


def _loss(content_mg_per_gdcw: float, ratio: float) -> float:
    content = content_mg_per_gdcw / BETA_CAROTENE_MW_G_PER_MOL
    return 1.0 - float(np.exp(-_coeff(ratio) * content))


class TestTheEffectIsBelowTheNoiseAtTheseTitres:
    @pytest.mark.parametrize("ratio", SPECTRAL_RATIO_BRACKET)
    def test_the_worst_case_loss_is_under_the_well_cv(self, ratio):
        """Across the whole spectral bracket and the brightest calibration strain, the
        reporter loses less signal to its own pigment than one well differs from its
        neighbour. Correcting it would be correcting inside the noise."""
        worst = _loss(max(ELIZONDO_CONTENT_MG_PER_GDCW), ratio)

        assert worst < DEFAULT_WELL_CV

    def test_and_that_holds_at_the_top_of_the_bracket(self):
        assert _loss(0.97, 0.25) < DEFAULT_WELL_CV


class TestItBecomesWorthMeasuringHigherUp:
    @pytest.mark.parametrize("ratio", SPECTRAL_RATIO_BRACKET)
    def test_the_threshold_is_between_one_and_six_mg_per_gdcw(self, ratio):
        """Where the loss reaches the well CV. The bracket puts it at 1.2 mg/gDCW if the
        pigment still absorbs a quarter of its peak at 516 nm, and 6.0 if a twentieth."""
        threshold = (-np.log(1 - DEFAULT_WELL_CV) / _coeff(ratio)
                     * BETA_CAROTENE_MW_G_PER_MOL)

        assert 1.0 < threshold < 6.5

    def test_published_high_producers_are_far_past_it(self):
        """Arhar 2024, PMID 39215465, reaches 79 mg/gDCW. There the correction is not
        optional -- and it is also where the refusal starts protecting something."""
        assert _loss(79.0, 0.05) > 3 * DEFAULT_WELL_CV


class TestTheDerivationRefusesTheWrongInputs:
    def test_the_450_peak_is_refused_at_the_excitation_wavelength(self):
        """The whole point: reusing the peak would overstate the coupling several-fold."""
        with pytest.raises(ValueError, match="falls steeply"):
            inner_filter_coeff_from_extinction(PEAK, 516.0, PATH_CM, BIOMASS_G_L, FLATTENING)

    def test_a_stronger_absorber_couples_harder(self):
        assert _coeff(0.25) > _coeff(0.05)
