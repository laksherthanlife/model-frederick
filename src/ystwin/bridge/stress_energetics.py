"""Reporter signal -> stress state -> maintenance energy -> growth, composed in one place.

The pieces of this chain were built at different times and never joined:

    reporter.py                 fluorescence -> promoter activity, by inverting
                                dR/dt = k_synth - (mu + k_deg)*R
    latent_bridge.py            activity / full scale -> a stress fraction, then
                                NGAM = resting + fraction * scale, as a PRE-SOLVE constraint
    maintenance_calibration.py  what that scale actually measures (added 2026-09-02)
    fba/dynamic_rates.py        NGAM -> the GEM -> mu -> biomass -> titre

This module is the join, and it exists to make the gates on that chain visible at the point
of use rather than one file apart.

GATE ONE, on the reporter half. `latent_bridge` refuses unless a caller acknowledges that
G4 returned INCONCLUSIVE for all four local constructs. Read what that refusal actually
says: "Hac1 anchor fails its no-RT control in 3/3 replicates; TRX2 anchor underpowered at
n=3, needing ~9 and ~80 replicates". **Those are objections to the qPCR REFERENCE GENES and
to replicate count -- to the assay, not to the mechanism.** The distinction matters because
it says what would lift the refusal: a working anchor and more replicates, not a new theory.
It is recorded here because "G4 refused it" has been repeated in this project as though the
mechanism were dead, and it is not what G4 said.

GATE TWO, on the energetics half. The maintenance scale
`latent_bridge._MAX_STRESS_MAINTENANCE = 6.5` mmol ATP/gDW/h is an envelope by its own
admission. The measurement it stood in for puts the increment at 0.26-0.32 -- twenty times
smaller. What reaches growth, though, is the TOTAL and not the increment: 7.2 asserted
against 1.31-1.64 measured, a factor of
<!-- audit:value table=outputs/maintenance_scale.csv column=penalty_ratio_asserted_over_measured row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->4.39 wherever the response stays inside one metabolic regime.

GATE THREE, ON THE REGIME, and it is a retraction of what this docstring used to say. It
read: growth on Yeast9 is exactly linear in the NGAM bound, slope
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.004642 /h per mmol ATP/gDW/h, unchanged from glucose -1.0 to -10.0 with oxygen
unlimited -- **so that factor carries straight through and does so independently of the
carbon supply.** Every clause of that is true and the conclusion drawn from it was not. The
table it rested on held oxygen_lower_bound = -1000.0 in every row: the invariance is across
glucose while the cell stays fully respiratory, which is the one axis that leaves the slope
alone. At `fba/physiology.REFERENCE_AEROBIC_BATCH`, the operating point this repository
validates its GEM against, overflow is ON -- the reference carries ethanol_secretion 13.9
mmol/gDCW/h and the model secretes
<!-- audit:value table=outputs/maintenance_scale.csv column=ethanol_at_zero_ngam_mmol_per_gdcw_h row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->15.923 -- and the slope there is
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.010778,
<!-- audit:value table=outputs/maintenance_scale.csv column=slope_ratio_to_oxygen_unlimited row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->2.322x larger. So
:attr:`StressEnergetics.growth_penalty_per_h` takes a
:class:`~ystwin.bridge.maintenance_calibration.MaintenanceRegime` and REFUSES without one.
The NGAM interval is regime-free and is still returned; what a regime buys is the only step
that converts it into growth.

The slope is MODEL-DERIVED / IN_SILICO, not MEASURED, and this module said MEASURED. It is
the gradient of a linear program under bounds a script chose.

In fractions of growth at full-scale reporter, with oxygen unlimited:
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->37.52% of growth asserted against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->8.55% measured at
glucose -1.0, and <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->3.75% against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.85% at -10.0. At the fermenting
reference point, on more glucose than either, it is
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->21.95% against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->5.00%.
Produced by `scripts/maintenance_scale.py`.

AND THE SIGN IS UNRESOLVED, which is why this module composes but does not decide. The
measurement compares a strain that CAN mount the general stress response against one that
cannot, and the responder spends LESS. `latent_bridge` adds maintenance as its reporter
RISES. A reporter reading high means the response is ON, which by the measurement means
maintenance is lower, not higher. Those may be opposite. Resolving it needs a reporter and a
maintenance measurement on the SAME cultures, which is exactly what PMID 40181231 did with
pHSP12-mRuby2 on a plate reader -- the same instrument this repository already reads.

So: the chain is joined, the scale is measured, the slope is now a function of the regime,
the sign is an open experiment, and nothing here sets a number until the gates lift.
"""

from __future__ import annotations

from dataclasses import dataclass

from .maintenance_calibration import (
    MAINTENANCE_REGIMES,
    MEASURED_MAINTENANCE,
    MaintenanceRegime,
    Refused,
    stress_maintenance_increment,
)

__all__ = [
    "REGIME_IS_REQUIRED",
    "SIGN_IS_UNRESOLVED",
    "StressEnergetics",
    "maintenance_from_reporter",
]

#: Stated once, imported wherever the chain is used, so no caller has to rediscover it.
SIGN_IS_UNRESOLVED = (
    "the measured maintenance difference compares GSR-PRESENT against GSR-ABSENT, not "
    "stressed against unstressed, and the responder spends LESS -- so whether a rising "
    "reporter should RAISE or LOWER maintenance is not settled by it. What settles it is a "
    "reporter and a maintenance measurement on the same cultures"
)

#: How much worse the fermenting slope is than the respiratory one, COMPUTED from the two
#: regimes rather than typed, so the refusal below cannot quote a ratio that has gone stale.
_FERMENTING_SLOPE_RATIO = (
    MAINTENANCE_REGIMES["aerobic_batch_reference"].ngam_growth_slope_per_h_per_mmol_atp
    / MAINTENANCE_REGIMES["respiratory_glucose_excess"].ngam_growth_slope_per_h_per_mmol_atp)

#: The refusal carried by a growth penalty asked for without a regime. Its own constant for
#: the same reason as the one above: a caller that hits it should read the reason, not a
#: TypeError from three frames down.
REGIME_IS_REQUIRED = (
    "no maintenance regime supplied, so there is no NGAM -> growth slope. It is not a "
    f"constant: it is {_FERMENTING_SLOPE_RATIO:.2f}x larger at the fermenting operating "
    "point fba/physiology.py validates against than under unlimited oxygen, and it steps "
    "rather than sliding between them. Pass "
    "regime=maintenance_calibration.regime_for(glucose, oxygen)"
)


@dataclass(frozen=True)
class StressEnergetics:
    """One reporter reading carried through to a maintenance constraint, with its gates."""

    activity_rfu_per_od_h: float
    stress_fraction: float
    ngam_low: float
    ngam_high: float
    scale_source: str
    reporter_validated: bool
    sign_resolved: bool
    regime: MaintenanceRegime | None = None

    @property
    def usable_for_a_number(self) -> bool:
        """The two gates on the NGAM VALUE. A caller that ignores this is asserting past two
        open questions. The regime is not one of them -- it gates the conversion to growth,
        and refuses there, because an NGAM interval is a perfectly good answer without one."""
        return self.reporter_validated and self.sign_resolved

    @property
    def growth_penalty_per_h(self) -> tuple[float, float] | Refused:
        """Growth lost to this maintenance interval, in 1/h, or a :class:`Refused`.

        The refusal is the whole point of the regime argument: the slope this multiplies by
        is not a property of the cell, it is a property of the cell AND its oxygen supply.
        """
        if self.regime is None:
            return Refused(REGIME_IS_REQUIRED)
        return (self.regime.growth_penalty_per_h(self.ngam_low),
                self.regime.growth_penalty_per_h(self.ngam_high))

    def summary(self) -> str:
        band = (f"NGAM {self.ngam_low:.3f}-{self.ngam_high:.3f} mmol ATP/gDW/h"
                if self.ngam_low != self.ngam_high
                else f"NGAM {self.ngam_low:.3f} mmol ATP/gDW/h")
        gates = []
        if not self.reporter_validated:
            gates.append("reporter assay not validated (G4: anchor genes, replicate count)")
        if not self.sign_resolved:
            gates.append("sign unresolved")
        tail = f"; REFUSED FOR USE -- {'; '.join(gates)}" if gates else ""
        line = (f"activity {self.activity_rfu_per_od_h:.1f} RFU/OD/h -> stress "
                f"{self.stress_fraction:.3f} -> {band} [{self.scale_source}]")

        # The growth cost is a DERIVED step and appears only when a regime priced it. The
        # two gates above are on the NGAM value, which is a complete answer without one.
        if self.regime is not None:
            low, high = self.growth_penalty_per_h
            cost = (f"-{low:.4f} to -{high:.4f} /h" if low != high else f"-{low:.4f} /h")
            line += f" -> growth {cost} in {self.regime.name}"
        return line + tail


def maintenance_from_reporter(activity_rfu_per_od_h: float,
                              activity_full_scale: float,
                              *,
                              use_measured_scale: bool = True,
                              atp_per_glucose: float | None = None,
                              regime: MaintenanceRegime | None = None,
                              reporter_validated: bool = False,
                              sign_resolved: bool = False) -> StressEnergetics:
    """Carry a reporter reading through to an NGAM interval, and optionally to growth.

    Args:
        activity_rfu_per_od_h: Promoter activity from `reporter.py`.
        activity_full_scale: The reading treated as full induction. Required, not defaulted,
            because it is a property of the CONSTRUCT and the instrument -- `latent_bridge`
            uses 1305.0 for its four constructs and that number does not transfer to a
            different reporter on a different plate reader.
        use_measured_scale: Use the measured maintenance increment rather than
            `latent_bridge`'s envelope. Default True; passing False reproduces the old
            behaviour and should be labelled where it is used.
        atp_per_glucose: Pin the P/O and get a point instead of an interval. Asserting one
            is the caller's to make, so there is no default.
        regime: The operating point, from
            :func:`~ystwin.bridge.maintenance_calibration.regime_for`. Without it the NGAM
            interval is still returned and :attr:`StressEnergetics.growth_penalty_per_h`
            REFUSES -- there is no default regime, because a default regime is the scalar
            slope that made this argument necessary.
        reporter_validated, sign_resolved: The two gates on the NGAM value. Both default
            False, which is the true state, and both must be True before
            :attr:`StressEnergetics.usable_for_a_number`.

    Raises:
        ValueError: on a non-positive full scale, or a negative activity.
    """
    if activity_full_scale <= 0:
        raise ValueError(
            f"activity_full_scale must be positive, got {activity_full_scale}. It is the "
            f"reading treated as full induction, and it belongs to the construct and the "
            f"instrument -- there is no universal value to fall back on")
    if activity_rfu_per_od_h < 0:
        raise ValueError(f"activity must be non-negative, got {activity_rfu_per_od_h}")

    fraction = min(1.0, activity_rfu_per_od_h / activity_full_scale)

    # The resting base and the stress increment are the SAME measurement converted by the
    # SAME ratio, so a caller who pins the P/O pins both. Converting only the increment left
    # the base spanning 1.056-1.320 while the increment was a point, which is not a narrower
    # answer -- it is two different assumptions inside one number.
    resting = MEASURED_MAINTENANCE["gsr_intact"]
    if atp_per_glucose is None:
        resting_low, resting_high = resting.atp_interval()
    else:
        resting_low = resting_high = resting.as_atp(atp_per_glucose)

    if use_measured_scale:
        increment_low, increment_high = stress_maintenance_increment(atp_per_glucose)
        source = f"measured, {MEASURED_MAINTENANCE['gsr_intact'].source}"
        low = resting_low + fraction * increment_low
        high = resting_high + fraction * increment_high
    else:
        from .latent_bridge import _MAX_STRESS_MAINTENANCE, _RESTING_MAINTENANCE
        source = "asserted envelope from latent_bridge, ~20x the measured effect"
        low = high = _RESTING_MAINTENANCE + fraction * _MAX_STRESS_MAINTENANCE

    return StressEnergetics(
        activity_rfu_per_od_h=float(activity_rfu_per_od_h),
        stress_fraction=fraction, ngam_low=low, ngam_high=high,
        scale_source=source,
        reporter_validated=reporter_validated, sign_resolved=sign_resolved,
        regime=regime)
