"""The product channel: beta-carotene read as absorbance on the plate reader we own.

The observation layer read biomass and reporter and not product, which left the last link
of the chain unmeasurable on the instrument this lab already has. Beta-carotene is a
pigment, so it is readable at 450 nm in the same well on the same run.

What these tests defend is the confound that makes that harder than it sounds. A whole-cell
reading at 450 nm is pigment absorbance **plus cell scattering**, and at the contents this
organism reaches the two terms are the same order. Crediting the whole reading to pigment
is the same class of error as the naive-versus-dilution-corrected fold already recorded in
this repository: a number that looks like a titre and is mostly turbidity. So the
whole-cell path refuses without a measured scattering ratio, and the tests pin the refusal
as hard as they pin the arithmetic.

The sign test at the bottom is the one that would cost the most if it inverted. Product
attenuates the reporter's excitation, so a producing strain reads *lower* RFU at the same
promoter activity. Get that backwards and apparent stress is inflated in exactly the
producing strains the twin exists to describe, with the error growing as the strain gets
better.
"""

import numpy as np
import pytest

from ystwin.calib.od import ODCalibration
from ystwin.observation import BETA_CAROTENE_PETROLEUM_ETHER_450 as _EPS
from ystwin.observation import (
    BETA_CAROTENE_MW_G_PER_MOL,
    BETA_CAROTENE_PETROLEUM_ETHER_450,
    PigmentExtinction,
    PigmentOptics,
    ReporterOptics,
    correct_inner_filter,
    inner_filter_coeff_from_extinction,
    observe_absorbance,
    observe_absorbance_extract,
    observe_rfu,
    product_from_absorbance,
)
from ystwin.readings import RawRFU

#: The solvent the extinction coefficient was measured in. Read off the constant rather
#: than restated, so a test cannot assert against a solvent the coefficient is not for.
_SOLVENT = _EPS.solvent

# A fully specified channel: every number that could have been assumed is stated.
CALIBRATED = PigmentOptics(
    wavelength_nm=450.0,
    extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
    path_length_cm=0.5,
    blank=0.040,
    scatter_per_od600=1.35,
    flattening=0.40,
)
# The same channel before anyone measured the scattering ratio.
UNCALIBRATED = PigmentOptics(
    wavelength_nm=450.0,
    extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
    path_length_cm=0.5,
    blank=0.040,
)
GDCW_PER_OD = 0.45


class TestTheExtinctionCoefficient:
    def test_it_names_the_solvent_it_was_measured_in(self):
        """Beta-carotene's absorptivity is a property of the pigment *in a solvent*.

        Its maximum moves by 16 nm between petroleum ether and chloroform, so a
        coefficient carried without its solvent cannot be checked and cannot be reused.
        Three carotenoid-engineering papers quote 138,900 M-1 cm-1 at 450 nm citing one
        handbook and name no solvent; that is the failure this field exists to prevent
        repeating here.
        """
        assert BETA_CAROTENE_PETROLEUM_ETHER_450.solvent.strip()
        assert BETA_CAROTENE_PETROLEUM_ETHER_450.wavelength_nm == 450.0

    def test_it_is_the_tabulated_specific_absorbance_times_the_molecular_weight(self):
        """A(1%,1cm) = 2592 in petroleum ether at 450 nm, times MW / 10.

        Pinned as arithmetic on the published quantity rather than as a bare number, so a
        future edit that changes the value has to change the derivation with it.
        """
        expected = 2592 * BETA_CAROTENE_MW_G_PER_MOL / 10.0

        assert BETA_CAROTENE_PETROLEUM_ETHER_450.epsilon_per_M_per_cm == pytest.approx(
            expected, rel=1e-3
        )

    def test_a_coefficient_with_no_solvent_is_refused_at_construction(self):
        with pytest.raises(ValueError, match="solvent"):
            PigmentExtinction(epsilon_per_M_per_cm=1.0e5, wavelength_nm=450.0,
                              solvent="  ", source="somewhere")

    def test_reading_at_a_wavelength_the_coefficient_was_not_measured_at_is_refused(self):
        """Lycopene reads at 470 nm and beta-carotene at 450, off one shared filter wheel.

        Absorbance falls steeply either side of the maximum, so pairing a 450 nm
        coefficient with a 470 nm filter is an error of unknown size that produces a
        perfectly plausible number.
        """
        with pytest.raises(ValueError, match="470.*450|450.*470"):
            PigmentOptics(wavelength_nm=470.0,
                          extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
                          path_length_cm=0.5)


class TestScatteringIsTheCrux:
    def test_a_well_with_no_product_still_reads_well_above_the_blank(self):
        """Zero product is not zero absorbance: it is blank plus cells.

        This is the whole reason the channel needs a non-producing control. A reader
        subtracting only the media blank attributes every bit of the cell term to titre.
        """
        reading = observe_absorbance(product=0.0, biomass=0.9, optics=CALIBRATED,
                                     gdcw_per_od=GDCW_PER_OD)

        assert reading == pytest.approx(0.040 + 1.35 * 0.9 / GDCW_PER_OD, rel=1e-9)
        assert reading > CALIBRATED.blank

    def test_an_empty_well_reads_only_the_media_blank(self):
        assert observe_absorbance(0.0, 0.0, CALIBRATED, gdcw_per_od=GDCW_PER_OD) == (
            pytest.approx(0.040)
        )

    def test_the_scattering_term_is_the_same_order_as_the_pigment_term(self):
        """Not a rounding correction. At the content FBA calls a ceiling, they are
        comparable, which is why this channel refuses rather than approximates.

        32.84 mg/gDCW is `EXTERNAL_CAROTENOID_BOUND.md`'s arithmetic on the D1 frontier,
        used here only as a scale for the comparison.
        """
        content = 32.84 / BETA_CAROTENE_MW_G_PER_MOL  # mmol/gDCW
        biomass = 0.9

        total = observe_absorbance(content, biomass, CALIBRATED, gdcw_per_od=GDCW_PER_OD)
        scatter_only = observe_absorbance(0.0, biomass, CALIBRATED, gdcw_per_od=GDCW_PER_OD)
        pigment = total - scatter_only
        scatter = scatter_only - CALIBRATED.blank

        assert 0.1 < pigment / scatter < 10.0

    def test_without_a_measured_scattering_ratio_the_reading_is_refused(self):
        """REFUSE rather than default. Taking the ratio as zero would credit cell
        turbidity to titre, and the error grows with biomass rather than washing out."""
        with pytest.raises(ValueError, match="scatter_per_od600"):
            observe_absorbance(1e-3, 0.9, UNCALIBRATED, gdcw_per_od=GDCW_PER_OD)

    def test_the_refusal_names_the_extracted_sample_route_out_of_it(self):
        """A refusal that does not say what to do instead gets worked around."""
        with pytest.raises(ValueError, match="observe_absorbance_extract"):
            observe_absorbance(1e-3, 0.9, UNCALIBRATED, gdcw_per_od=GDCW_PER_OD)

    def test_without_a_dry_weight_factor_the_reading_is_refused(self):
        """The scattering ratio is per OD600 and the state is in g/L; the bridge between
        them is strain- and medium-specific and is measured, never assumed."""
        blind = ODCalibration(saturation_k=0.2, top_true_od=2.0, blank=0.05)

        with pytest.raises(ValueError, match="gdcw_per_od"):
            observe_absorbance(1e-3, 0.9, CALIBRATED, calibration=blind)

    def test_the_dry_weight_factor_is_taken_from_the_od_calibration_when_measured(self):
        measured = ODCalibration(saturation_k=0.2, top_true_od=2.0, blank=0.05,
                                 gdcw_per_od=GDCW_PER_OD)

        assert observe_absorbance(1e-3, 0.9, CALIBRATED, calibration=measured) == (
            pytest.approx(observe_absorbance(1e-3, 0.9, CALIBRATED,
                                             gdcw_per_od=GDCW_PER_OD))
        )


class TestPackaging:
    def test_without_a_stated_flattening_factor_the_reading_is_refused(self):
        """Intracellular pigment absorbs less than the same quantity in solution.

        Duysens 1956 named it; applying a solution extinction coefficient to intact cells
        over-reads by a factor nobody in this project has measured. The refusal is what
        stops that factor being silently taken as 1.
        """
        no_flattening = PigmentOptics(
            wavelength_nm=450.0, extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
            path_length_cm=0.5, blank=0.040, scatter_per_od600=1.35,
        )

        with pytest.raises(ValueError, match="flattening"):
            observe_absorbance(1e-3, 0.9, no_flattening, gdcw_per_od=GDCW_PER_OD)

    def test_stating_it_as_one_is_permitted_and_is_the_upper_bound(self):
        """``flattening=1.0`` is an explicit claim, not a default, and it reads highest."""
        unpackaged = PigmentOptics(
            wavelength_nm=450.0, extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
            path_length_cm=0.5, blank=0.040, scatter_per_od600=1.35, flattening=1.0,
        )

        assert observe_absorbance(1e-3, 0.9, unpackaged, gdcw_per_od=GDCW_PER_OD) > (
            observe_absorbance(1e-3, 0.9, CALIBRATED, gdcw_per_od=GDCW_PER_OD)
        )

    def test_a_flattening_factor_above_one_is_refused(self):
        """Packaging can only reduce apparent absorbance. A fitted value above 1 is a
        scattering ratio that is too small, not a real optical effect."""
        with pytest.raises(ValueError, match="flattening"):
            PigmentOptics(wavelength_nm=450.0,
                          extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
                          path_length_cm=0.5, scatter_per_od600=1.0, flattening=1.4)


class TestTheSignalItself:
    def test_absorbance_rises_with_product_at_fixed_biomass(self):
        contents = np.linspace(0.0, 0.06, 25)
        readings = [observe_absorbance(c, 0.9, CALIBRATED, gdcw_per_od=GDCW_PER_OD)
                    for c in contents]

        assert np.all(np.diff(readings) > 0)

    def test_absorbance_rises_with_biomass_at_fixed_product(self):
        """Both terms scale with cells: scattering directly, pigment because content is
        per gram. So the channel alone cannot separate a denser culture from a richer one,
        which is why it is read beside OD600 rather than instead of it."""
        densities = np.linspace(0.1, 1.5, 25)
        readings = [observe_absorbance(1e-3, x, CALIBRATED, gdcw_per_od=GDCW_PER_OD)
                    for x in densities]

        assert np.all(np.diff(readings) > 0)

    def test_the_pigment_term_is_beer_lambert_on_the_well_concentration(self):
        content, biomass = 2.0e-3, 0.8
        concentration_m = content * biomass / 1000.0
        expected = (0.40 * BETA_CAROTENE_PETROLEUM_ETHER_450.epsilon_per_M_per_cm
                    * concentration_m * 0.5)

        total = observe_absorbance(content, biomass, CALIBRATED, gdcw_per_od=GDCW_PER_OD)
        scatter_only = observe_absorbance(0.0, biomass, CALIBRATED, gdcw_per_od=GDCW_PER_OD)

        assert total - scatter_only == pytest.approx(expected, rel=1e-9)

    def test_the_detector_saturates_when_a_ceiling_is_declared(self):
        """The product channel saturates long before OD600 does: at the contents in the
        literature the pigment term alone runs past 2 AU in a 0.5 cm path. A clipped
        reading is not a small error, it is a censored one."""
        capped = PigmentOptics(
            wavelength_nm=450.0, extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
            path_length_cm=0.5, blank=0.040, scatter_per_od600=1.35, flattening=0.40,
            detector_max=3.5,
        )

        assert observe_absorbance(0.06, 1.2, capped, gdcw_per_od=GDCW_PER_OD) == (
            pytest.approx(3.5)
        )

    def test_a_reading_below_the_ceiling_is_left_alone(self):
        capped = PigmentOptics(
            wavelength_nm=450.0, extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
            path_length_cm=0.5, blank=0.040, scatter_per_od600=1.35, flattening=0.40,
            detector_max=3.5,
        )

        assert observe_absorbance(1e-3, 0.5, capped, gdcw_per_od=GDCW_PER_OD) == (
            pytest.approx(observe_absorbance(1e-3, 0.5, CALIBRATED,
                                             gdcw_per_od=GDCW_PER_OD))
        )

    def test_a_negative_product_is_refused(self):
        with pytest.raises(ValueError, match="product"):
            observe_absorbance(-1e-4, 0.9, CALIBRATED, gdcw_per_od=GDCW_PER_OD)


class TestTheExtractedSample:
    """The path that needs no uncalibrated number. It still refuses on a solvent
    mismatch, because the extinction coefficient is only valid in its own solvent."""

    def test_an_extract_needs_neither_a_scattering_ratio_nor_a_flattening_factor(self):
        reading = observe_absorbance_extract(
            product=2.0e-3, biomass_harvested_gdcw=0.01, extract_volume_l=1.0e-3,
            optics=UNCALIBRATED, extraction_solvent=_SOLVENT,
        )

        assert reading > UNCALIBRATED.blank

    def test_it_is_beer_lambert_with_the_solution_coefficient(self):
        product, harvested, volume = 2.0e-3, 0.01, 1.0e-3
        concentration_m = product * harvested / volume / 1000.0
        expected = (0.040 + BETA_CAROTENE_PETROLEUM_ETHER_450.epsilon_per_M_per_cm
                    * concentration_m * 0.5)

        assert observe_absorbance_extract(product, harvested, volume, UNCALIBRATED, _SOLVENT) == (
            pytest.approx(expected, rel=1e-9)
        )

    def test_concentrating_the_extract_into_less_solvent_raises_the_reading(self):
        loose = observe_absorbance_extract(2e-3, 0.01, 2.0e-3, UNCALIBRATED, _SOLVENT)
        tight = observe_absorbance_extract(2e-3, 0.01, 1.0e-3, UNCALIBRATED, _SOLVENT)

        assert tight - UNCALIBRATED.blank == pytest.approx(
            2.0 * (loose - UNCALIBRATED.blank), rel=1e-9
        )

    def test_a_zero_volume_extract_is_refused_rather_than_dividing_by_zero(self):
        with pytest.raises(ValueError, match="extract_volume_l"):
            observe_absorbance_extract(2e-3, 0.01, 0.0, UNCALIBRATED, _SOLVENT)


class TestInvertingTheReading:
    """The direction the instrument is actually used in."""

    def test_it_recovers_the_content_that_produced_the_reading(self):
        content, biomass = 3.0e-3, 0.9
        reading = observe_absorbance(content, biomass, CALIBRATED,
                                     gdcw_per_od=GDCW_PER_OD)

        assert product_from_absorbance(reading, biomass, CALIBRATED,
                                       gdcw_per_od=GDCW_PER_OD) == pytest.approx(
            content, rel=1e-9
        )

    def test_ignoring_the_scattering_floor_would_have_inflated_the_titre(self):
        """The size of the error the scattering term prevents, made explicit.

        Inverting with the scattering ratio set to zero is what a naive reader does. It is
        not a small bias; at these numbers it is several-fold, and it grows with biomass.
        """
        content, biomass = 3.0e-3, 0.9
        reading = observe_absorbance(content, biomass, CALIBRATED,
                                     gdcw_per_od=GDCW_PER_OD)
        naive_optics = PigmentOptics(
            wavelength_nm=450.0, extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
            path_length_cm=0.5, blank=0.040, scatter_per_od600=0.0, flattening=0.40,
        )

        naive = product_from_absorbance(reading, biomass, naive_optics,
                                        gdcw_per_od=GDCW_PER_OD)

        assert naive > 2.0 * content

    def test_a_saturated_reading_is_refused_rather_than_inverted(self):
        """A clipped reading carries no concentration information; inverting it would
        manufacture one, and it would be the same value for every well past the ceiling."""
        capped = PigmentOptics(
            wavelength_nm=450.0, extinction=BETA_CAROTENE_PETROLEUM_ETHER_450,
            path_length_cm=0.5, blank=0.040, scatter_per_od600=1.35, flattening=0.40,
            detector_max=3.5,
        )

        with pytest.raises(ValueError, match="ceiling"):
            product_from_absorbance(3.5, 0.9, capped, gdcw_per_od=GDCW_PER_OD)

    def test_a_well_below_the_scattering_floor_returns_a_negative_rather_than_a_zero(self):
        """Clipping to zero would hide a wrong scattering ratio behind a plausible titre.

        A negative content is not a measurement; it is a statement that the non-producing
        control does not describe these cells.
        """
        floor = observe_absorbance(0.0, 0.9, CALIBRATED, gdcw_per_od=GDCW_PER_OD)

        assert product_from_absorbance(floor - 0.05, 0.9, CALIBRATED,
                                       gdcw_per_od=GDCW_PER_OD) < 0

    def test_an_empty_well_carries_no_content_to_infer(self):
        with pytest.raises(ValueError, match="biomass"):
            product_from_absorbance(0.5, 0.0, CALIBRATED, gdcw_per_od=GDCW_PER_OD)


# --------------------------------------------------------------------------- #
# the coupling back into the reporter channel
# --------------------------------------------------------------------------- #

PIGMENTED_REPORTER = ReporterOptics(gain=2.0e5, background=250.0, autofluorescence=400.0,
                                    inner_filter_coeff=180.0)
UNCALIBRATED_REPORTER = ReporterOptics(gain=2.0e5, background=250.0,
                                       autofluorescence=400.0, inner_filter_coeff=None)


class TestTheInnerFilterSign:
    def test_product_reduces_apparent_fluorescence(self):
        """The sign that costs the most if it inverts.

        Beta-carotene absorbs across 400-500 nm and mCitrine is excited near 516 nm, so
        pigment eats excitation light and the reporter reads *low*. Invert this and the
        inferred promoter activity is too high in exactly the producing strains the twin
        exists to describe, with the error growing as the strain improves -- an artefact
        that would look like product-induced stress.
        """
        clean = observe_rfu(1e-4, 0.9, carotenoid=0.0, optics=PIGMENTED_REPORTER)
        producing = observe_rfu(1e-4, 0.9, carotenoid=4e-3, optics=PIGMENTED_REPORTER)

        assert producing < clean

    def test_the_attenuation_deepens_monotonically_with_titre(self):
        titre = np.linspace(0.0, 8e-3, 30)
        rfu = [observe_rfu(1e-4, 0.9, c, PIGMENTED_REPORTER) for c in titre]

        assert np.all(np.diff(rfu) < 0)

    def test_the_correction_raises_the_reading_back_and_never_lowers_it(self):
        """The correction undoes an attenuation, so it can only move the signal up.

        A correction that lowered a pigmented well's reading would have the sign the wrong
        way round twice and still look self-consistent.
        """
        titre = np.linspace(0.0, 8e-3, 30)
        raw = np.array([observe_rfu(1e-4, 0.9, c, PIGMENTED_REPORTER) for c in titre])

        corrected = correct_inner_filter(RawRFU(raw), titre, PIGMENTED_REPORTER)

        assert np.all(corrected >= raw - 1e-9)

    def test_correcting_with_an_uncalibrated_coefficient_is_refused(self):
        """Unchanged behaviour, pinned again from the product side: an uncalibrated filter
        must not be silently treated as an absent one."""
        with pytest.raises(ValueError, match="inner_filter_coeff"):
            correct_inner_filter(RawRFU(np.array([1200.0])), np.array([4e-3]),
                                 UNCALIBRATED_REPORTER)


class TestTheCoefficientImpliedByTheProductChannel:
    """The two channels see one pigment, so one physical quantity should set both."""

    def test_it_is_the_beer_lambert_coefficient_in_natural_log_units(self):
        at_excitation = PigmentExtinction(
            epsilon_per_M_per_cm=2.0e4, wavelength_nm=516.0, solvent="petroleum ether",
            source="a hypothetical measurement at the mCitrine excitation wavelength",
        )
        expected = np.log(10.0) * 0.40 * 2.0e4 * 0.9 * 0.5 / 1000.0

        assert inner_filter_coeff_from_extinction(
            at_excitation, excitation_nm=516.0, path_length_cm=0.5, biomass=0.9,
            flattening=0.40,
        ) == pytest.approx(expected, rel=1e-12)

    def test_the_derived_coefficient_attenuates_in_the_same_direction(self):
        at_excitation = PigmentExtinction(
            epsilon_per_M_per_cm=2.0e4, wavelength_nm=516.0, solvent="petroleum ether",
            source="a hypothetical measurement at the mCitrine excitation wavelength",
        )
        coeff = inner_filter_coeff_from_extinction(
            at_excitation, 516.0, path_length_cm=0.5, biomass=0.9, flattening=0.40,
        )
        derived = ReporterOptics(gain=2.0e5, background=250.0, autofluorescence=400.0,
                                 inner_filter_coeff=coeff)

        assert observe_rfu(1e-4, 0.9, 4e-3, derived) < observe_rfu(1e-4, 0.9, 0.0, derived)

    def test_using_the_450_nm_peak_for_a_516_nm_excitation_is_refused(self):
        """Beta-carotene at 516 nm is on the far shoulder of its band, not at its peak.

        Reusing the 450 nm coefficient would overstate the coupling several-fold and
        therefore over-correct the reporter -- inflating apparent activity in producing
        strains, which is the same failure as getting the sign wrong, only quieter.
        """
        with pytest.raises(ValueError, match="516|450"):
            inner_filter_coeff_from_extinction(
                BETA_CAROTENE_PETROLEUM_ETHER_450, excitation_nm=516.0,
                path_length_cm=0.5, biomass=0.9, flattening=0.40,
            )



class TestTheSolventRefusalActuallyFires:
    """The docstring promised a ValueError on solvent mismatch and the signature had no
    solvent argument, so it could never raise. A carotenoid extinction coefficient is a
    property of the pigment IN A SOLVENT -- beta-carotene's maximum moves 16 nm between
    petroleum ether and chloroform -- so the mismatch is a wrong answer, silently."""

    def test_a_mismatched_solvent_is_refused(self):
        with pytest.raises(ValueError, match="measured in"):
            observe_absorbance_extract(2e-3, 0.01, 1e-3, UNCALIBRATED, "acetone")

    def test_the_matching_solvent_is_accepted(self):
        assert observe_absorbance_extract(2e-3, 0.01, 1e-3, UNCALIBRATED, _SOLVENT) > 0

    def test_the_comparison_ignores_case_and_padding_only(self):
        """A solvent name is free text, so trivial formatting must not refuse -- but a
        different solvent must, and nothing else may be treated as equivalent."""
        assert observe_absorbance_extract(2e-3, 0.01, 1e-3, UNCALIBRATED,
                                          f"  {_SOLVENT.upper()}  ") > 0
        with pytest.raises(ValueError):
            observe_absorbance_extract(2e-3, 0.01, 1e-3, UNCALIBRATED, _SOLVENT + " ether")
