"""Measurement models mapping twin state to instrument readings.

RFU = background + [autofluorescence + gain*R] * biomass * exp(-eps * carotenoid).

The attenuation term is a path from the quantity being predicted back into the channel
predicting it, so it is real signal that shuffled-channel controls cannot detect. It is inert
until a carotenoid-producing strain exists, and refuses to correct without a measured
coefficient.

Beta-carotene is a pigment, so the third quantity in the chain -- the product -- is readable
on the same plate reader by absorbance near 450 nm. One well, one run, biomass and reporter
and titre. What stops that being free is that a whole-cell reading at 450 nm is
*mostly cells*: the measured attenuation is pigment absorbance plus light scattering, and
scattering at 450 nm is the same order as the pigment term at any titre this organism
reaches. Myers, Curtis & Curtis 2013 (PMID 24499615) put the split plainly -- OD at a
"robust" wavelength is scattering and does not vary with culture conditions, OD at a
"sensitive" wavelength is scattering *plus* the organism's pigments -- and their remedy is
the one taken here: correlate against off-peak attenuation measured on cells that carry no
pigment. So the whole-cell path takes a measured scattering ratio and refuses without one;
the extracted-sample path needs no such number because the cells are gone.

The second confound is that intracellular pigment is not dissolved pigment. Packaging a
chromophore into particles flattens its apparent absorption spectrum relative to the same
quantity in solution (Duysens 1956, doi 10.1016/0006-3002(56)90380-8), so a solution
extinction coefficient over-reads a whole-cell sample. That is carried as an explicit
``flattening`` factor which must be stated, never inferred.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .calib.od import ODCalibration
from .readings import RawRFU, require

__all__ = [
    "BETA_CAROTENE_MW_G_PER_MOL",
    "BETA_CAROTENE_PETROLEUM_ETHER_450",
    "PigmentExtinction",
    "PigmentOptics",
    "ReporterOptics",
    "correct_inner_filter",
    "inner_filter_coeff_from_extinction",
    "observe_absorbance",
    "observe_absorbance_extract",
    "observe_od",
    "observe_rfu",
    "product_from_absorbance",
]


@dataclass(frozen=True)
class ReporterOptics:
    """Instrument- and strain-specific reporter measurement parameters.

    Args:
        gain: RFU per unit mature reporter concentration per gDCW biomass.
        background: Media-only reading, RFU.
        autofluorescence: RFU per gDCW from cells carrying no reporter. Measure it
            on an isogenic reporter-free strain, not by subtracting a guess.
        inner_filter_coeff: Attenuation coefficient for beta-carotene, per unit
            intracellular concentration. ``None`` means it has not been calibrated;
            correction is then refused rather than assumed to be zero. Measure it by
            spiking purified carotenoid into wells of a non-producing strain.
        detector_max: Reading at which the photomultiplier saturates, if known.
    """

    gain: float
    background: float
    autofluorescence: float
    inner_filter_coeff: float | None
    detector_max: float | None = None

    def __post_init__(self) -> None:
        for name in ("gain", "background", "autofluorescence"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative")
        if self.inner_filter_coeff is not None and self.inner_filter_coeff < 0:
            raise ValueError("inner_filter_coeff must be non-negative")
        if self.detector_max is not None and self.detector_max <= self.background:
            raise ValueError("detector_max must exceed the background reading")

    def attenuation(self, carotenoid: np.ndarray) -> np.ndarray:
        """Fraction of excitation surviving the cell's own pigment."""
        if not self.inner_filter_coeff:
            return np.ones_like(np.asarray(carotenoid, dtype=float))
        return np.exp(-self.inner_filter_coeff * np.asarray(carotenoid, dtype=float))


def observe_od(
    viable: float,
    dead: float,
    calibration: ODCalibration,
    gdcw_per_od: float | None = None,
) -> float:
    """Predicted raw OD reading from the biomass state.

    Dead and injured cells still scatter light, so optical density reads total
    biomass. That is precisely why OD alone cannot report viability.

    Args:
        viable: Viable biomass, g/L.
        dead: Non-viable biomass, g/L.
        calibration: Reader response, supplying saturation and blank.
        gdcw_per_od: Dry weight per linearised OD unit. Falls back to the
            calibration's own measured factor.
    """
    factor = gdcw_per_od if gdcw_per_od is not None else calibration.gdcw_per_od
    if factor is None:
        raise ValueError(
            "no gdcw_per_od available; supply it here or measure it into the calibration"
        )
    true_od = (viable + dead) / factor
    saturated = true_od / (1.0 + calibration.saturation_k * true_od)
    return float(saturated + calibration.blank)


def observe_rfu(
    reporter_mature: float,
    biomass: float,
    carotenoid: float,
    optics: ReporterOptics,
) -> float:
    """Predicted raw reporter reading from the twin state.

    Args:
        reporter_mature: Mature reporter concentration, per gDCW.
        biomass: Biomass in the optical path, g/L.
        carotenoid: Intracellular beta-carotene concentration, same basis as the
            calibrated ``inner_filter_coeff``.
        optics: Instrument and strain parameters.
    """
    per_cell = optics.autofluorescence + optics.gain * reporter_mature
    signal = per_cell * biomass * float(optics.attenuation(np.asarray(carotenoid)))
    reading = optics.background + signal
    if optics.detector_max is not None:
        reading = min(reading, optics.detector_max)
    return float(reading)


def correct_inner_filter(
    rfu: RawRFU,
    carotenoid: np.ndarray,
    optics: ReporterOptics,
) -> RawRFU:
    """Undo product-driven attenuation, returning the reading an unpigmented cell would give.

    Applied before any dilution correction or latent-state inference, so that the
    reporter channel is not partly a measurement of titre.

    Args:
        rfu: Raw readings. Raw in, raw out: the background is taken off to scale the
            attenuation and then added back, so the result is still an instrument reading
            and is typed as one.
        carotenoid: Matching intracellular beta-carotene, measured or predicted.
        optics: Must carry a calibrated ``inner_filter_coeff``.

    Raises:
        ValueError: if the coefficient has not been measured. Treating an
            uncalibrated filter as absent would silently keep the confound.
    """
    if optics.inner_filter_coeff is None:
        raise ValueError(
            "inner_filter_coeff has not been calibrated for these optics; "
            "measure it with a purified-carotenoid spike-in before correcting"
        )
    raw = require(rfu, RawRFU, name="rfu")
    above_background = raw.minus_background(optics.background).array
    return RawRFU(optics.background + above_background / optics.attenuation(carotenoid))


# --------------------------------------------------------------------------- #
# the product channel: beta-carotene by absorbance
# --------------------------------------------------------------------------- #

BETA_CAROTENE_MW_G_PER_MOL = 536.87
"""C40H56. The same figure `scripts/parked/run_d1.py` converts flux to content with."""


@dataclass(frozen=True)
class PigmentExtinction:
    """A molar extinction coefficient, with the two things that make it meaningful.

    A carotenoid extinction coefficient is not a property of the carotenoid. It is a
    property of the carotenoid *in a solvent*, at a *wavelength*, and beta-carotene's
    absorption maximum moves by 16 nm between petroleum ether and chloroform. A number
    carried without those two fields is not a measurement of anything, so they are
    required fields rather than documentation.

    Args:
        epsilon_per_M_per_cm: Molar extinction coefficient, M^-1 cm^-1.
        wavelength_nm: The wavelength it was measured at.
        solvent: The solvent it was measured in.
        source: Where it came from, resolvable.
    """

    epsilon_per_M_per_cm: float
    wavelength_nm: float
    solvent: str
    source: str

    def __post_init__(self) -> None:
        if self.epsilon_per_M_per_cm <= 0:
            raise ValueError("epsilon_per_M_per_cm must be positive")
        if self.wavelength_nm <= 0:
            raise ValueError("wavelength_nm must be positive")
        if not self.solvent.strip():
            raise ValueError(
                "an extinction coefficient with no solvent named is not usable; "
                "beta-carotene's absorptivity and its maximum both move with solvent"
            )
        if not self.source.strip():
            raise ValueError("an extinction coefficient with no source is not checkable")


BETA_CAROTENE_PETROLEUM_ETHER_450 = PigmentExtinction(
    # A(1%, 1 cm) = 2592 in petroleum ether at 450 nm, times the molecular weight over 10:
    # 2592 * 536.87 / 10 = 1.3916e5 M^-1 cm^-1.
    epsilon_per_M_per_cm=1.3916e5,
    wavelength_nm=450.0,
    solvent="petroleum ether",
    source=(
        "A(1%,1cm) = 2592 for beta-carotene in petroleum ether at 450 nm, tabulated in "
        "Britton, Liaaen-Jensen & Pfander, Carotenoids: Handbook (2004), "
        "doi 10.1007/978-3-0348-7836-4, and quoted with that solvent and wavelength in "
        "doi 10.1038/s41598-026-45956-6. Converted to molar units with MW 536.87 g/mol. "
        "Corroborated to 0.2% by carotenoid-engineering groups quoting the same handbook "
        "as 138,900 M^-1 cm^-1 at 450 nm, though naming no solvent (PMC4510654)."
    ),
)
"""Solution extinction coefficient for beta-carotene. NOT a whole-cell coefficient.

Applying this to intact cells over-reads, because pigment packaged into particles absorbs
less than the same quantity in solution -- see ``PigmentOptics.flattening``.
"""


@dataclass(frozen=True)
class PigmentOptics:
    """Instrument, strain and geometry parameters for one absorbance channel.

    Args:
        wavelength_nm: The filter being read. Must match the extinction coefficient's
            wavelength: beta-carotene's absorbance falls off steeply either side of its
            maximum, so reading at 470 nm with a 450 nm coefficient is a silent error of
            unknown size.
        extinction: Solution extinction coefficient at that wavelength.
        path_length_cm: Optical path through the well. In a 96-well plate this is set by
            the fill volume, not by the plate, so it is per-experiment.
        blank: Media-only reading at this wavelength.
        scatter_per_od600: Absorbance contributed at this wavelength by cells alone, per
            unit linearised OD600. ``None`` means it has not been measured and the
            whole-cell reading is then refused rather than credited entirely to pigment.
            Measure it on an isogenic *non-producing* strain across a dilution series,
            reading both wavelengths, and take the slope of A(lambda) on linearised OD600.
        flattening: Apparent absorbance of intracellular pigment as a fraction of the same
            quantity in solution, in (0, 1]. ``None`` refuses; ``1.0`` is the explicit
            claim that there is no packaging effect, which for a membrane-localised
            carotenoid is a strong claim and should be typed deliberately.
        detector_max: Reading at which the detector saturates, if known.
    """

    wavelength_nm: float
    extinction: PigmentExtinction
    path_length_cm: float
    blank: float = 0.0
    scatter_per_od600: float | None = None
    flattening: float | None = None
    detector_max: float | None = None

    def __post_init__(self) -> None:
        if self.path_length_cm <= 0:
            raise ValueError("path_length_cm must be positive")
        if self.blank < 0:
            raise ValueError("blank must be non-negative")
        if abs(self.wavelength_nm - self.extinction.wavelength_nm) > 1e-9:
            raise ValueError(
                f"reading at {self.wavelength_nm:.1f} nm with an extinction coefficient "
                f"measured at {self.extinction.wavelength_nm:.1f} nm; supply a coefficient "
                "for the filter actually in the reader"
            )
        if self.scatter_per_od600 is not None and self.scatter_per_od600 < 0:
            raise ValueError("scatter_per_od600 must be non-negative")
        if self.flattening is not None and not 0 < self.flattening <= 1:
            raise ValueError(
                "flattening must be in (0, 1]: packaging a chromophore into cells can only "
                "reduce its apparent absorbance, never raise it "
                "(Duysens 1956, doi 10.1016/0006-3002(56)90380-8)"
            )
        if self.detector_max is not None and self.detector_max <= self.blank:
            raise ValueError("detector_max must exceed the blank reading")

    def _absorptivity_per_mmol_per_gdcw(self) -> float:
        """Whole-cell absorbance per (mmol/gDCW of product) per (gDCW/L of biomass).

        Raises:
            ValueError: if the flattening factor has not been stated.
        """
        if self.flattening is None:
            raise ValueError(
                "flattening has not been stated for these optics; intracellular pigment "
                "absorbs less than the same quantity in solution, so applying a solution "
                "extinction coefficient to intact cells over-reads by an unknown factor. "
                "Measure it as the ratio of whole-cell to extract absorbance on one "
                "culture, or state flattening=1.0 deliberately"
            )
        # q [mmol/gDCW] * X [gDCW/L] / 1000 = c [mol/L]; A = flattening * eps * c * l.
        return (self.flattening * self.extinction.epsilon_per_M_per_cm
                * self.path_length_cm / 1000.0)

    def _scatter(self, biomass: float, gdcw_per_od: float) -> float:
        if self.scatter_per_od600 is None:
            raise ValueError(
                "scatter_per_od600 has not been measured for this strain and reader; a "
                f"whole-cell reading at {self.wavelength_nm:.0f} nm is pigment absorbance "
                "plus cell scattering and the two are the same order of magnitude, so "
                "crediting the whole reading to pigment inflates titre without limit. "
                "Measure it on an isogenic non-producer, or read an extracted sample with "
                "observe_absorbance_extract"
            )
        return self.scatter_per_od600 * biomass / gdcw_per_od


def _resolve_gdcw_per_od(calibration: ODCalibration | None, gdcw_per_od: float | None) -> float:
    factor = gdcw_per_od
    if factor is None and calibration is not None:
        factor = calibration.gdcw_per_od
    if factor is None:
        raise ValueError(
            "no gdcw_per_od available; the scattering term is defined per unit OD600 and "
            "the state is in g/L, so the conversion has to be measured, not assumed"
        )
    if factor <= 0:
        raise ValueError("gdcw_per_od must be positive")
    return float(factor)


def observe_absorbance(
    product: float,
    biomass: float,
    optics: PigmentOptics,
    calibration: ODCalibration | None = None,
    gdcw_per_od: float | None = None,
) -> float:
    """Predicted whole-cell absorbance reading from the twin state.

    ``A = blank + scatter_per_od600 * OD600 + flattening * eps * (product * biomass / 1000)
    * path_length``, clipped at the detector ceiling.

    The scattering term is not a correction, it is a comparable share of the reading. It is
    also the term that makes this channel worth having: it is measurable once, on a strain
    that carries no pathway, and then it is the same for every producing well on the plate.

    Args:
        product: Intracellular beta-carotene, mmol/gDCW -- the basis
            ``kinetic/carotenoid.py`` carries pools in.
        biomass: Total biomass in the optical path, g/L. Dead cells scatter, so this is
            viable plus dead, exactly as for :func:`observe_od`.
        optics: Instrument and strain parameters, carrying the measured scattering ratio.
        calibration: OD calibration, used only for its dry-weight factor.
        gdcw_per_od: Dry weight per linearised OD unit, overriding the calibration's.

    Raises:
        ValueError: if the scattering ratio, the flattening factor or the dry-weight
            factor has not been measured. Each of the three would otherwise be silently
            taken as a convenient value, and two of them are large.
    """
    if product < 0:
        raise ValueError("product must be non-negative")
    if biomass < 0:
        raise ValueError("biomass must be non-negative")
    factor = _resolve_gdcw_per_od(calibration, gdcw_per_od)
    # Scattering is checked first because it is the larger and the more easily forgotten
    # of the two missing measurements, so it is the one a refusal should name.
    scatter = optics._scatter(biomass, factor)
    pigment = optics._absorptivity_per_mmol_per_gdcw() * product * biomass
    reading = optics.blank + scatter + pigment
    if optics.detector_max is not None:
        reading = min(reading, optics.detector_max)
    return float(reading)


def observe_absorbance_extract(
    product: float,
    biomass_harvested_gdcw: float,
    extract_volume_l: float,
    optics: PigmentOptics,
    extraction_solvent: str,
) -> float:
    """Predicted absorbance of a solvent extract, where the cells are gone.

    This is the path that needs no uncalibrated quantity: no cells, so no scattering; no
    packaging, so the solution extinction coefficient is the right one. It costs a harvest,
    an extraction and a well, and it is the reference the whole-cell channel's
    ``scatter_per_od600`` and ``flattening`` are calibrated against.

    Args:
        product: Intracellular beta-carotene of the harvested cells, mmol/gDCW.
        biomass_harvested_gdcw: Dry weight harvested into the extract, gDCW.
        extract_volume_l: Volume the pigment was taken up into, L.
        optics: Instrument parameters. ``flattening`` and ``scatter_per_od600`` are
            deliberately not consulted -- neither applies to a cell-free solution.
        extraction_solvent: The solvent the pigment was actually taken up into. Required,
            and required to match the coefficient's own solvent.

    Raises:
        ValueError: if ``extraction_solvent`` is not the solvent the extinction
            coefficient was measured in. That is the whole reason ``solvent`` is a field.

            This argument was missing and the clause above was unenforceable: the function
            accepted a petroleum-ether coefficient for an acetone extract and returned a
            number. Beta-carotene's absorption maximum moves 16 nm between petroleum ether
            and chloroform, so the answer was wrong by whatever that mismatch costs and
            nothing said so. A documented refusal that cannot fire is worse than no
            docstring, because a reader takes the guarantee and stops checking.
    """
    if extraction_solvent.strip().casefold() != optics.extinction.solvent.strip().casefold():
        raise ValueError(
            f"extracted into {extraction_solvent!r} but the extinction coefficient was "
            f"measured in {optics.extinction.solvent!r}. A carotenoid extinction "
            "coefficient is a property of the pigment in a solvent, not of the pigment, "
            "so the two must match. Supply a coefficient measured in the solvent you used."
        )
    if product < 0 or biomass_harvested_gdcw < 0:
        raise ValueError("product and harvested biomass must be non-negative")
    if extract_volume_l <= 0:
        raise ValueError("extract_volume_l must be positive")
    concentration_m = product * biomass_harvested_gdcw / extract_volume_l / 1000.0
    reading = (optics.blank
               + optics.extinction.epsilon_per_M_per_cm * concentration_m
               * optics.path_length_cm)
    if optics.detector_max is not None:
        reading = min(reading, optics.detector_max)
    return float(reading)


def product_from_absorbance(
    reading: float,
    biomass: float,
    optics: PigmentOptics,
    calibration: ODCalibration | None = None,
    gdcw_per_od: float | None = None,
) -> float:
    """Invert :func:`observe_absorbance`: reading and biomass -> intracellular content.

    This is the direction the instrument is used in, and it is where the scattering term
    earns its keep: subtracting a measured scattering floor is the difference between a
    titre and a turbidity.

    A negative result is returned rather than clipped. It means the scattering floor
    over-subtracts, which is information -- either the non-producer the ratio came from is
    not isogenic, or this well has no product -- and clipping it to zero would hide both.

    Raises:
        ValueError: if the reading is at or above the declared detector ceiling, where the
            reading carries no concentration information and inverting it invents one.
    """
    if biomass <= 0:
        raise ValueError("cannot infer intracellular content from a well with no biomass")
    if optics.detector_max is not None and reading >= optics.detector_max:
        raise ValueError(
            f"reading {reading:.3f} is at or above the detector ceiling "
            f"{optics.detector_max:.3f}; dilute the well rather than inverting a clipped "
            "reading"
        )
    factor = _resolve_gdcw_per_od(calibration, gdcw_per_od)
    pigment = reading - optics.blank - optics._scatter(biomass, factor)
    return float(pigment / (optics._absorptivity_per_mmol_per_gdcw() * biomass))


def inner_filter_coeff_from_extinction(
    extinction: PigmentExtinction,
    excitation_nm: float,
    path_length_cm: float,
    biomass: float,
    flattening: float,
) -> float:
    """The reporter's inner-filter coefficient implied by the product's absorbance.

    ``ReporterOptics.inner_filter_coeff`` is an empirical number from a spike-in. This
    computes what it should be from first principles, which makes the coupling between the
    two channels one physical quantity rather than two unrelated fitted constants: the
    pigment that the 450 nm channel measures is the same pigment that eats the reporter's
    excitation light.

    Transmitted fraction is ``10^-A``, so with ``A = flattening * eps * q * X / 1000 * l``
    the natural-log coefficient is ``ln(10) * flattening * eps * X * l / 1000`` per
    mmol/gDCW of product.

    **This is a lower bound on the true attenuation.** It accounts for the excitation path
    only. mCitrine emits near 529 nm, still inside beta-carotene's absorption band, so the
    emitted photons are attenuated on the way out as well and the measured spike-in
    coefficient should exceed this one. A measured coefficient *below* this bound is a
    result, not a rounding difference.

    Args:
        extinction: Coefficient at the reporter's *excitation* wavelength, not at 450 nm.
        excitation_nm: The reporter's excitation wavelength, ~516 nm for mCitrine.
        path_length_cm: Optical path, cm.
        biomass: Biomass of the wells the coefficient will be applied to, g/L. The
            coefficient is per intracellular concentration but the absorbance is per well
            concentration, so it carries a biomass and is only valid near that density.
        flattening: Packaging factor, as for :class:`PigmentOptics`.

    Raises:
        ValueError: if the coefficient was not measured at the excitation wavelength.
            Beta-carotene's absorbance at 516 nm is a fraction of its 450 nm peak, so
            reusing the peak value here would overstate the coupling several-fold.
    """
    if abs(extinction.wavelength_nm - excitation_nm) > 1e-9:
        raise ValueError(
            f"extinction coefficient is for {extinction.wavelength_nm:.1f} nm but the "
            f"reporter is excited at {excitation_nm:.1f} nm; beta-carotene's absorbance "
            "falls steeply across that gap and the peak value would overstate the "
            "attenuation"
        )
    if path_length_cm <= 0 or biomass < 0:
        raise ValueError("path_length_cm must be positive and biomass non-negative")
    if not 0 < flattening <= 1:
        raise ValueError("flattening must be in (0, 1]")
    return float(math.log(10.0) * flattening * extinction.epsilon_per_M_per_cm
                 * biomass * path_length_cm / 1000.0)
