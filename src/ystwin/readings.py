"""What a reading is, and whether its blank has been taken off yet.

Every number this package reads off a plate is one of five things, and until 2026-08-31
all five were the same thing: a bare ``numpy`` array. Nothing distinguished an optical
density with the medium's absorbance still in it from one with the medium subtracted, or
either of those from a fluorescence, or a fluorescence from a fluorescence already divided
by density. The convention was carried in **parameter names**, and
``docs/ARCHITECTURE.md`` §6.3 documented it in three tables running to forty lines --
which is the length a convention reaches just before it stops being followed.

The convention was also almost right, and the exact shape of "almost" is the argument for
replacing it. Counted on the tree this replaced: **22 public functions took one of these
quantities. Six took a blank alongside it and did the subtraction themselves; the other 16
wanted it already done.** That rule is easy to state and was followed.

Five of those six gave the blank a default of ``0.0``. So a caller could omit it and
silently mean "already corrected" while the signature said "raw", and the rule had an
escape hatch in 5 of the 6 places it applied. A convention with an escape hatch and no
enforcement is one a reader trusts and a compiler does not.

**What goes wrong, concretely.** Hand ``growth.specific_growth_rate`` a raw OD and it
returns a number. The blank is roughly 0.09 on this project's plates against cultures
spanning 0.09 to 0.61, so ``d ln(OD)/dt`` computed on the uncorrected trace understates
the growth rate early and converges to the right answer late. The result is not an error,
a warning or a NaN. It is a slightly wrong growth rate, and since
``k_synth = dR/dt + (mu + k_deg) R`` divides by exactly that quantity, it becomes a
slightly wrong promoter activity, which is the number this repository exists to report.

So the state is a **type** now, and the arithmetic that changes the state lives in one
place rather than at every call site:

    RawOD ---- minus_blank(blank) ----> CorrectedOD ----+
                                                        +--> SpecificFluorescence
    RawRFU --- minus_background(bg) --> CorrectedRFU ---+

**Why five and not two.** ``CorrectedRFU`` and ``SpecificFluorescence`` are different
quantities and §6.3's third table is about confusing them: ``reporter.promoter_activity``
documents its ``reporter`` argument as a per-cell concentration, and every real caller
passes RFU divided by OD, so the function returns reader units rather than molar ones.
A type that could not tell those apart would leave the sharpest of the three documented
hazards exactly where it was.

**What this does not claim.** ``CorrectedOD(raw_array)`` is a lie a caller can still tell,
and no newtype in a language without linear types can stop them. What it stops is the
*accident* -- the wrong variable of the right shape, which is the failure that actually
happens, because both are float arrays of the same length off the same plate. The
declaration is now at the call site where the author knows the answer, instead of in a
table in a document they have not opened.

**Shape is deliberately not part of the type.** A well is an ndarray and a plate is a
``pandas.DataFrame``, and both can be raw or corrected; ``gates/g1_optical.assess_plate``
and ``diagnostics/dilution_confound.detect_blank_wells`` take frames while everything else
takes traces. The tag answers "has the blank come off", which is orthogonal to whether the
container has one column or ninety-six, and a type that fused the two would need eight
names to say what five say.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = [
    "CorrectedOD",
    "CorrectedRFU",
    "RawOD",
    "RawRFU",
    "Reading",
    "ReporterPerGDCW",
    "SpecificFluorescence",
    "require",
]


@dataclass(frozen=True)
class Reading:
    """One measured channel, in a known correction state.

    Args:
        values: The measurements. An ndarray for a single well's trace, a
            ``pandas.DataFrame`` for a whole plate. A list or tuple is accepted and
            converted; a ``Reading`` is not, because wrapping one twice is a slip rather
            than an intention and the second wrapper would silently claim a state the
            first had already contradicted.
    """

    values: np.ndarray | pd.DataFrame

    def __post_init__(self) -> None:
        if isinstance(self.values, Reading):
            raise TypeError(
                f"{type(self).__name__} was given a {type(self.values).__name__}, which is "
                "already a reading. Pass the measurements, or use the conversion method if "
                "the state needs to change")
        if not isinstance(self.values, pd.DataFrame):
            object.__setattr__(self, "values", np.asarray(self.values, dtype=float))

    @property
    def array(self) -> np.ndarray:
        """The measurements as an ndarray, whatever container they arrived in."""
        return (self.values.to_numpy(dtype=float) if isinstance(self.values, pd.DataFrame)
                else self.values)

    def __array__(self, dtype=None, copy=None) -> np.ndarray:
        """The escape hatch for arithmetic inside a consumer that has already checked the
        type -- and it honours numpy's copy protocol, which the first version did not.

        ``np.array(reading, copy=True)`` returned the reading's own buffer, so writing to
        the "copy" wrote through to the reading. A type introduced to stop one silent
        aliasing mistake should not ship another; caught by testing the hatch rather than
        the types.
        """
        values = self.array
        if copy is False:
            if dtype is not None and np.dtype(dtype) != values.dtype:
                raise ValueError(
                    f"cannot view {type(self).__name__} as {np.dtype(dtype)} without "
                    "copying; pass copy=None to allow one")
            return values
        out = values if dtype is None else values.astype(dtype, copy=False)
        return out.copy() if copy else out

    def __len__(self) -> int:
        return len(self.values)


@dataclass(frozen=True)
class RawOD(Reading):
    """Optical density as the reader reported it, with the medium's absorbance included.

    This is what comes off an instrument file and out of ``plate/synergy.py``. Nothing
    that takes a logarithm or a ratio may use it: an OD of 0.09 on this project's plates
    is an empty well, not a hundredth of a culture.
    """

    def minus_blank(self, blank: float) -> "CorrectedOD":
        """Subtract the media-only absorbance.

        The one place the subtraction happens, so that "which blank" is a question with a
        single answer per call rather than a line of arithmetic repeated at each site.

        Raises:
            ValueError: if ``blank`` is not finite. A NaN blank propagates into every
                downstream density and turns a plate into NaNs several functions later,
                where the cause is no longer visible.
        """
        if not np.isfinite(blank):
            raise ValueError(f"blank must be finite, got {blank}")
        return CorrectedOD(self.values - blank)


@dataclass(frozen=True)
class CorrectedOD(Reading):
    """Optical density with the medium subtracted: density attributable to cells.

    What ``growth.py`` differentiates, what ``calib/od.ODCalibration.linearise`` takes,
    and the denominator of :class:`SpecificFluorescence`. Strictly positive wherever a
    logarithm is taken, which the consumers check for themselves.
    """


@dataclass(frozen=True)
class RawRFU(Reading):
    """Fluorescence as the reader reported it: reporter, medium and autofluorescence.

    The medium and the cells' own glow are both in here. ``docs/FINDINGS.md`` puts the
    reporter-free per-cell signal near 249 RFU/OD against a median reporter well of about
    4,000, so this is not a small term.
    """

    def minus_background(self, background: float) -> "CorrectedRFU":
        """Subtract the media-plus-autofluorescence baseline.

        Raises:
            ValueError: if ``background`` is not finite.
        """
        if not np.isfinite(background):
            raise ValueError(f"background must be finite, got {background}")
        return CorrectedRFU(self.values - background)


@dataclass(frozen=True)
class CorrectedRFU(Reading):
    """Fluorescence with the background subtracted: total signal from the reporter.

    **Total, not per cell.** A well with twice the biomass gives twice this number at the
    same promoter activity, which is why dividing by density is a separate step with a
    separate type rather than something a consumer does inline.
    """

    def per(self, optical_density: CorrectedOD) -> "SpecificFluorescence":
        """Divide by blank-corrected density to get signal per unit density.

        Raises:
            TypeError: if handed a :class:`RawOD`. Dividing a corrected fluorescence by an
                uncorrected density is the composite of the two mistakes this module
                exists to separate, and it is the one that looks most like ordinary code.
        """
        require(optical_density, CorrectedOD, name="optical_density")
        signal, density = self.array, optical_density.array
        if signal.shape != density.shape or signal.size == 0:
            raise ValueError("fluorescence and optical density must have the same nonempty shape")
        frames = (isinstance(self.values, pd.DataFrame),
                  isinstance(optical_density.values, pd.DataFrame))
        if any(frames):
            if not all(frames):
                raise ValueError("fluorescence and optical density must both carry frame labels")
            if (not self.values.index.equals(optical_density.values.index)
                    or not self.values.columns.equals(optical_density.values.columns)
                    or not self.values.index.is_unique or not self.values.columns.is_unique):
                raise ValueError("fluorescence and optical density must have matching unique labels")
        if not np.all(np.isfinite(signal)):
            raise ValueError("fluorescence must be finite before division by optical density")
        if not np.all(np.isfinite(density) & (density > 0)):
            raise ValueError("optical density must be strictly positive and finite")
        return SpecificFluorescence(self.values / optical_density.values)


@dataclass(frozen=True)
class ReporterPerGDCW(Reading):
    """Background-corrected fluorescence per gram dry cell weight.

    What :class:`~ystwin.estimator.ParticleFilter` holds as its ``reporter`` state, because
    :class:`~ystwin.observation.ReporterOptics` states its gain as "RFU per unit mature
    reporter per gDCW". It differs from :class:`SpecificFluorescence` by exactly the factor
    ``gdcw_per_od``, and ``docs/ARCHITECTURE.md`` §6.3 recorded that the two "are never used
    together and there is no converter". There is one now:
    :meth:`SpecificFluorescence.per_gdcw`, which requires the factor rather than defaulting
    it, for the reason that method's docstring gives.
    """


@dataclass(frozen=True)
class SpecificFluorescence(Reading):
    """Background-corrected fluorescence per blank-corrected optical density.

    **Reader units, not molar**, and that is the point of giving it a name.
    ``reporter.promoter_activity`` documents its input as "per-cell mature reporter
    concentration" and every caller in this repository passes this instead, so what comes
    back is RFU/OD/h. ``estimator.ParticleFilter``'s reporter state is genuinely per gDCW
    -- it feeds ``observation.observe_rfu``, whose gain is RFU per unit reporter per gDCW
    -- and the two are never interconverted anywhere. The missing factor is
    ``gdcw_per_od``; ``generator/calibrate.py`` reports ``gain * gdcw_per_od`` as a single
    product and lists the pair as unidentifiable without a measured dry weight, which is
    the honest treatment.
    """

    def per_gdcw(self, gdcw_per_od: float) -> "ReporterPerGDCW":
        """Convert to the per-gDCW basis the particle filter works in.

        ``docs/ARCHITECTURE.md`` §6.3 ends its per-OD/per-gDCW section with "the two are
        never used together and **there is no converter** … Wiring `estimator.py` onto the
        real-data path requires reconciling this first". This is that reconciliation, and it
        is a method on the only type that can legitimately be its input rather than a free
        function anything could call.

        **The factor is required and never defaulted, because it is not identifiable from a
        plate.** ``generator/calibrate.py`` fits ``gain * gdcw_per_od`` as a single product
        and lists ``("gain", "gdcw_per_od")`` as unidentifiable unless a measured dry weight
        is supplied. A default here would invent the half of a product nobody separated, and
        every activity downstream would carry it silently. 0.42 g/L per linearised OD unit is
        this project's working value; it belongs at a call site that can cite it.

        Args:
            gdcw_per_od: Grams dry cell weight per linearised OD unit, measured. Strictly
                positive.

        Raises:
            ValueError: if the factor is not positive and finite. Zero divides the trace to
                infinity and a negative one flips its sign, and both would read as a
                reporter doing something rather than as an arithmetic slip.
        """
        if not np.isfinite(gdcw_per_od) or gdcw_per_od <= 0:
            raise ValueError(
                f"gdcw_per_od must be a positive finite number, got {gdcw_per_od}. It is a "
                "measured dry-weight factor, not a default: generator/calibrate.py reports "
                "gain x gdcw_per_od as one product and calls the pair unidentifiable "
                "without a measurement")
        return ReporterPerGDCW(self.values / float(gdcw_per_od))


def require(value, *kinds: type, name: str = "value"):
    """Return ``value`` if it is one of ``kinds``, else raise a :class:`TypeError` that says
    what to do about it.

    The message names the conversion rather than the type, because "expected CorrectedOD,
    got RawOD" tells a caller what is wrong and not what to write.

    Args:
        value: The reading to check.
        kinds: Acceptable reading types, in the order a message should list them.
        name: The parameter being checked, for the message.
    """
    if isinstance(value, kinds):
        return value
    wanted = " or ".join(k.__name__ for k in kinds)
    if isinstance(value, Reading):
        hint = _CONVERSION.get((type(value).__name__, kinds[0].__name__))
        detail = f" -- {hint}" if hint else ""
        raise TypeError(
            f"{name} must be {wanted}, got {type(value).__name__}{detail}")
    raise TypeError(
        f"{name} must be {wanted}, got a bare {type(value).__name__}. Wrap it at the point "
        f"where you know the answer: {wanted}(values). Correction state used to be carried "
        "in the parameter name, which is why this is now a type -- see ystwin.readings")


#: How to get from the reading somebody passed to the one that was wanted. Only the pairs
#: a conversion exists for; the rest fall back to naming the types alone, because inventing
#: a route between two quantities that do not convert would be worse than saying nothing.
_CONVERSION = {
    ("RawOD", "CorrectedOD"): "call .minus_blank(od_blank) on it",
    ("RawRFU", "CorrectedRFU"): "call .minus_background(rfu_background) on it",
    ("CorrectedRFU", "SpecificFluorescence"): "call .per(corrected_od) on it",
    ("CorrectedOD", "RawOD"): "this function wants the reading before the blank came off",
    ("CorrectedRFU", "RawRFU"): "this function wants the reading before the background came off",
}
