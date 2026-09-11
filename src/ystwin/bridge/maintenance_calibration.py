"""Maintenance energy, measured, so the latent-stress branch stops resting on an envelope.

`bridge/latent_bridge.py` turns a stress reporter into an ATP maintenance constraint on the
GEM. Its scale constant, ``_MAX_STRESS_MAINTENANCE = 6.5`` mmol ATP/gDW/h, is honest about
what it is -- its own comment calls it "Tier 0 still: an envelope, not their measured value,
which is published as a bar chart". It was obtained by asking how much NGAM Yeast9 could
absorb at a plausible uptake, not by measuring what a stressed cell spends.

This module carries the measurement that envelope was standing in for.

MEASURED (Microb Biotechnol 2025, PMID 40181231, PMC11968331), aerobic glucose-limited
*S. cerevisiae* CEN.PK113-7D against a Δmsn2Δmsn4 derivative that cannot mount the general
stress response:

    parental, GSR intact      m_S = 0.066 +/- 0.032 mmol glucose/gDW/h
    Δmsn2Δmsn4, GSR absent    m_S = 0.082 +/- 0.045 mmol glucose/gDW/h

**THE ENVELOPE IS ABOUT TWENTY TIMES TOO LARGE.** At any plausible aerobic P/O the measured
total maintenance is 1.1-1.6 mmol ATP/gDW/h, so the branch's full-stress value of 7.2 sits
about fivefold above anything measured, and the stress INCREMENT it asserts (6.5) is roughly
twenty times the measured difference between having the stress response and lacking it.

TWO THINGS THIS MEASUREMENT IS NOT, and both matter before anyone uses it as a stress scale.

1. **It compares REGULATOR PRESENT against REGULATOR ABSENT, not stressed against
   unstressed.** The paper's finding is that the general stress response LOWERS the cost of
   maintaining homeostasis -- deleting Msn2/Msn4 RAISES m_S. That is a statement about what
   the GSR buys, not about what stress costs. A reporter reading high means the GSR is ON,
   which by this measurement means maintenance is LOWER than it would otherwise be. The
   latent branch adds maintenance as its reporter rises. **The two are not obviously the
   same sign, and this module does not pretend to resolve that** -- it supplies the
   magnitude and flags the sign as open.
2. **The strains carry no product.** This calibrates the host's energetics, not a pathway.

WHAT THE MAGNITUDE COSTS DEPENDS ON THE REGIME, and quoting it as a scalar is the error this
module shipped with. Measured on Yeast9 (`data/gem/yeast-GEM.xml.gz`), growth is linear in the
NGAM bound, and the slope was written here as one number --
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.004642 /h per mmol
ATP/gDW/h, "the same at glucose -1.0, -1.5 and -10.0". That is true, and it is the wrong
branch. Every row of `outputs/maintenance_scale.csv` carried oxygen_lower_bound = -1000.0, so
the famous invariance was measured across the CARBON supply while the cell stayed fully
respiratory -- ethanol
<!-- audit:value table=outputs/maintenance_scale.csv column=ethanol_at_zero_ngam_mmol_per_gdcw_h row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.000 mmol/gDCW/h in all of them. Oxygen is the axis that
moves the slope, and it was the one axis never swept.

AT THIS REPOSITORY'S OWN VALIDATED OPERATING POINT THE CELL FERMENTS.
`fba/physiology.REFERENCE_AEROBIC_BATCH` carries ethanol_secretion 13.9 mmol/gDCW/h, and under
its glucose and oxygen bounds Yeast9 secretes
<!-- audit:value table=outputs/maintenance_scale.csv column=ethanol_at_zero_ngam_mmol_per_gdcw_h row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->15.923. There a unit of NGAM costs
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.010778 /h per mmol ATP/gDW/h, which is
<!-- audit:value table=outputs/maintenance_scale.csv column=slope_ratio_to_oxygen_unlimited row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->2.322x the respiratory number this
module quoted. A fermenting cell buys ATP at a worse exchange rate, so maintenance costs it
more growth. :data:`MAINTENANCE_REGIMES` therefore carries a SLOPE PER OPERATING POINT and
:func:`regime_for` refuses any point that has not been solved, rather than interpolating one.

AND THE SLOPE IS RETAGGED. It was presented as MEASURED. It is MODEL-DERIVED / IN_SILICO: the
gradient of a linear program on a vendored SBML file under bounds a script chose. No culture
was read to obtain it. The measurement in this module is PMID 40181231's maintenance
coefficient and nothing else.

WHAT THE ABLATION SAYS, against a named floor. Replacing the regime slope with the single
respiratory one moves predicted growth at the reference operating point by
<!-- audit:value table=outputs/maintenance_scale.csv column=ablation_growth_error_at_asserted_per_h row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.044179 /h at the asserted NGAM total the latent branch actually spans, which is
<!-- audit:value table=outputs/maintenance_scale.csv column=ablation_over_growth_floor_at_asserted row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->3.776x
`generator/panel_experiment.MEASURED_GROWTH_RATE_SE`, the
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_assay_floor_per_h row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.0117 /h standard error of a
growth rate fitted off a real OD trace. **At the MEASURED maintenance total it is
<!-- audit:value table=outputs/maintenance_scale.csv column=ablation_growth_error_at_measured_per_h row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.010063 /h, which is BELOW that floor.** Said plainly: the regime matters at the
envelope's scale and is not resolvable at the measurement's, and both halves of that belong
here rather than only the half that argues for the change.

WHAT IS LINEAR, AND WHERE. Within one regime the response is exactly linear -- secant spread
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_spread_relative row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.0 with oxygen free. But the NGAM bound is itself part of the regime: raising
maintenance on a respiring cell eventually forces overflow, and the slope changes there. At
glucose -1.0 against oxygen -3.7 the model secretes no ethanol up to NGAM 7.2 and slopes at
the respiratory value; by NGAM 10.0 it ferments, and a single slope is then wrong by up to
<!-- audit:value table=outputs/maintenance_scale.csv column=slope_linearity_growth_error_per_h row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.048399 /h -- above the growth floor. Each regime below therefore carries its
own ``slope_linearity_growth_error_per_h``, and the five declared here are all rows where it
sits under the floor.

IN FRACTIONS OF GROWTH at full-scale reporter, with the asserted envelope against the measured
total. Oxygen free:
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->37.52% against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->8.55% at glucose -1.0,
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->25.01% against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->5.70% at -1.5,
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->3.75% against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.85% at -10.0.
At the fermenting reference point, on far more glucose,
<!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_asserted_percent row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->21.95% against <!-- audit:value table=outputs/maintenance_scale.csv column=growth_penalty_measured_percent row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->5.00% -- six times the penalty that
the glucose -10.0 row predicts, at the same carbon supply, purely because the oxygen is
capped.

The ratio between those two columns was reported as condition-independent at
<!-- audit:value table=outputs/maintenance_scale.csv column=penalty_ratio_asserted_over_measured row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->4.39. It is, wherever the response is linear -- it is a property of the two
NGAM values, not of the metabolism -- and it reaches
<!-- audit:value table=outputs/maintenance_scale.csv column=penalty_ratio_asserted_over_measured row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->8.491 in the row where the regime changes underneath it.

For the record, Yeast9 stops growing at NGAM
<!-- audit:value table=outputs/maintenance_scale.csv column=no_growth_ngam_mmol_atp_per_gdw_h row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->19.19 at glucose -1.0 and
<!-- audit:value table=outputs/maintenance_scale.csv column=no_growth_ngam_mmol_atp_per_gdw_h row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->191.92 at -10.0 with oxygen free, and at
<!-- audit:value table=outputs/maintenance_scale.csv column=no_growth_ngam_mmol_atp_per_gdw_h row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-2.0" -->25.73 at that same glucose under an oxygen cap of 2.0 -- a
property of BOTH supplies, and an earlier version of this docstring quoted the first of them
as though it were a constant.

All of the above is produced by `scripts/maintenance_scale.py` into
`outputs/maintenance_scale.csv`, which exists because none of it was checkable when it was
first written as prose -- and the sweep is now over oxygen as well as glucose, because a
table that varies one bound and prints the other cannot show which one it depends on.

"""

from __future__ import annotations

from dataclasses import dataclass

from ..fba.physiology import REFERENCE_AEROBIC_BATCH

__all__ = [
    "ETHANOL_DETECTION_MMOL_PER_GDCW_H",
    "MAINTENANCE_REGIMES",
    "MEASURED_MAINTENANCE",
    "MaintenanceMeasurement",
    "MaintenanceRegime",
    "REFERENCE_REGIME",
    "Refused",
    "regime_for",
    "regime_label",
    "stress_maintenance_increment",
]

#: ATP yielded per glucose in aerobic yeast. A RANGE and not a constant, because the P/P and
#: P/O ratios that set it are themselves contested -- textbook full oxidation gives 32, and
#: measured yeast values sit lower. Carried as an interval so a caller converting the glucose
#: measurement below into ATP has to carry the interval too, rather than picking a number and
#: presenting the result as exact.
ATP_PER_GLUCOSE_AEROBIC = (16.0, 20.0)

#: Below this an ethanol flux out of the GEM is solver noise rather than overflow. ASSERTED,
#: and numerical rather than biological: no published threshold is being invoked, and every
#: fermenting row of the sweep sits several orders of magnitude above it.
ETHANOL_DETECTION_MMOL_PER_GDCW_H = 1e-6

#: What every slope in :data:`MAINTENANCE_REGIMES` is, and it is not a measurement. Carried as
#: a string on each regime so the tag travels with the number instead of living in a docstring
#: the caller never opens.
IN_SILICO = ("MODEL-DERIVED / IN_SILICO: gradient of a linear program on Yeast9 "
             "(data/gem/yeast-GEM.xml.gz) under declared exchange bounds, from "
             "scripts/maintenance_scale.py -> outputs/maintenance_scale.csv. No culture was "
             "read to obtain it")


class Refused:
    """A number this repository declines to supply, which raises when anyone tries to use it.

    The idiom rather than ``None`` or a default. ``None`` propagates into arithmetic as a
    TypeError three frames away from the decision, and a default propagates silently, which
    is exactly how a single respiratory slope came to stand for every operating point.
    """

    __slots__ = ("reason",)

    def __init__(self, reason: str) -> None:
        self.reason = str(reason)

    def __float__(self) -> float:
        raise TypeError(f"REFUSED: {self.reason}")

    def __repr__(self) -> str:
        return f"REFUSED({self.reason})"

    def __bool__(self) -> bool:
        return False


def regime_label(ethanol_mmol_per_gdcw_h: float) -> str:
    """``"fermentative"`` or ``"respiratory"``, from the one flux that separates them.

    Overflow is the mechanism behind the slope difference, so the label reads the ethanol
    exchange rather than the oxygen bound -- a cell can be oxygen-capped and still respire
    everything it takes up, and at glucose -1.0 against oxygen -3.7 it does.
    """
    return ("fermentative" if ethanol_mmol_per_gdcw_h > ETHANOL_DETECTION_MMOL_PER_GDCW_H
            else "respiratory")


@dataclass(frozen=True)
class MaintenanceMeasurement:
    """One measured maintenance coefficient, in the units it was published in."""

    strain: str
    genotype: str
    m_s_mmol_glucose_per_gdw_h: float
    standard_deviation: float
    source: str

    def as_atp(self, atp_per_glucose: float) -> float:
        """Convert to mmol ATP/gDW/h, the unit the GEM's NGAM reaction is in."""
        if atp_per_glucose <= 0:
            raise ValueError(f"atp_per_glucose must be positive, got {atp_per_glucose}")
        return self.m_s_mmol_glucose_per_gdw_h * atp_per_glucose

    def atp_interval(self) -> tuple[float, float]:
        """The measurement over the P/O interval, which is the honest form."""
        low, high = ATP_PER_GLUCOSE_AEROBIC
        return self.as_atp(low), self.as_atp(high)


@dataclass(frozen=True)
class MaintenanceRegime:
    """One operating point, and what a unit of NGAM costs growth THERE.

    Every field but :attr:`name` and :attr:`provenance` is a cell of
    ``outputs/maintenance_scale.csv`` at this regime's two bounds, and
    `tests/test_s_atp_regime.py` fails if any of them drifts from the table.
    """

    name: str
    glucose_lower_bound_mmol_per_gdcw_h: float
    oxygen_lower_bound_mmol_per_gdcw_h: float
    growth_at_zero_ngam_per_h: float
    ethanol_at_zero_ngam_mmol_per_gdcw_h: float
    ngam_growth_slope_per_h_per_mmol_atp: float
    no_growth_ngam_mmol_atp_per_gdw_h: float
    slope_linearity_growth_error_per_h: float
    provenance: str = IN_SILICO

    @property
    def regime(self) -> str:
        return regime_label(self.ethanol_at_zero_ngam_mmol_per_gdcw_h)

    @property
    def ferments(self) -> bool:
        """Whether overflow is on here, which is what the slope difference is made of."""
        return self.regime == "fermentative"

    def growth_penalty_per_h(self, ngam_mmol_atp_per_gdw_h: float) -> float:
        """Growth lost to a maintenance bound, in 1/h.

        Raises:
            ValueError: on a negative bound, or on one at or past
                :attr:`no_growth_ngam_mmol_atp_per_gdw_h`, where the LP is infeasible and a
                linear extrapolation would return a penalty larger than the growth there is.
        """
        if ngam_mmol_atp_per_gdw_h < 0:
            raise ValueError(f"NGAM must be non-negative, got {ngam_mmol_atp_per_gdw_h}")
        if ngam_mmol_atp_per_gdw_h >= self.no_growth_ngam_mmol_atp_per_gdw_h:
            raise ValueError(
                f"NGAM {ngam_mmol_atp_per_gdw_h} mmol ATP/gDW/h is at or past the no-growth "
                f"threshold {self.no_growth_ngam_mmol_atp_per_gdw_h} for regime "
                f"{self.name!r}; the model is infeasible there, not merely slow")
        return self.ngam_growth_slope_per_h_per_mmol_atp * ngam_mmol_atp_per_gdw_h

    def growth_per_h(self, ngam_mmol_atp_per_gdw_h: float) -> float:
        """Growth rate under a maintenance bound, in 1/h."""
        return self.growth_at_zero_ngam_per_h - self.growth_penalty_per_h(
            ngam_mmol_atp_per_gdw_h)

    def summary(self) -> str:
        return (f"{self.name}: glucose {self.glucose_lower_bound_mmol_per_gdcw_h}, oxygen "
                f"{self.oxygen_lower_bound_mmol_per_gdcw_h} mmol/gDCW/h -> {self.regime}, "
                f"{self.ngam_growth_slope_per_h_per_mmol_atp:.6f} /h per mmol ATP/gDW/h "
                f"[{self.provenance.split(':')[0]}]")


#: The two measured states. Keyed by whether the general stress response is available.
MEASURED_MAINTENANCE = {
    "gsr_intact": MaintenanceMeasurement(
        strain="CEN.PK113-7D",
        genotype="parental",
        m_s_mmol_glucose_per_gdw_h=0.066,
        standard_deviation=0.032,
        source="PMID 40181231 (PMC11968331), Microb Biotechnol 2025"),
    "gsr_absent": MaintenanceMeasurement(
        strain="ScNtd002",
        genotype="Δmsn2 Δmsn4",
        m_s_mmol_glucose_per_gdw_h=0.082,
        standard_deviation=0.045,
        source="PMID 40181231 (PMC11968331), Microb Biotechnol 2025"),
}


#: Five solved operating points out of the sixteen `scripts/maintenance_scale.py` sweeps,
#: chosen because each one carries an argument. The two oxygen-free rows are the branch this
#: module used to quote as though it were general; `aerobic_batch_reference` is the phenotype
#: `fba/physiology.py` validates the GEM against and the one this repository actually solves
#: at; the last two walk the same glucose deeper into oxygen limitation.
#:
#: THE FIVE ARE ALL ROWS WHERE A SINGLE SLOPE IS LICENSED -- every
#: ``slope_linearity_growth_error_per_h`` here is under the growth assay's own floor. The
#: rows where it is not are in the table and deliberately not here: a regime object whose
#: slope is wrong by more than the assay can resolve is not a reduction, it is a mistake with
#: a dataclass around it.
MAINTENANCE_REGIMES = {
    regime.name: regime for regime in (
        MaintenanceRegime(
            name="respiratory_glucose_limited",
            glucose_lower_bound_mmol_per_gdcw_h=-1.0,
            oxygen_lower_bound_mmol_per_gdcw_h=-1000.0,
            growth_at_zero_ngam_per_h=0.089094,
            ethanol_at_zero_ngam_mmol_per_gdcw_h=0.0,
            ngam_growth_slope_per_h_per_mmol_atp=0.004642,
            no_growth_ngam_mmol_atp_per_gdw_h=19.19,
            slope_linearity_growth_error_per_h=0.0),
        MaintenanceRegime(
            name="respiratory_glucose_excess",
            glucose_lower_bound_mmol_per_gdcw_h=-10.0,
            oxygen_lower_bound_mmol_per_gdcw_h=-1000.0,
            growth_at_zero_ngam_per_h=0.890935,
            ethanol_at_zero_ngam_mmol_per_gdcw_h=0.0,
            ngam_growth_slope_per_h_per_mmol_atp=0.004642,
            no_growth_ngam_mmol_atp_per_gdw_h=191.92,
            slope_linearity_growth_error_per_h=0.0),
        MaintenanceRegime(
            name="aerobic_batch_reference",
            glucose_lower_bound_mmol_per_gdcw_h=-REFERENCE_AEROBIC_BATCH.glucose_uptake,
            oxygen_lower_bound_mmol_per_gdcw_h=-REFERENCE_AEROBIC_BATCH.oxygen_uptake,
            growth_at_zero_ngam_per_h=0.353529,
            ethanol_at_zero_ngam_mmol_per_gdcw_h=15.923,
            ngam_growth_slope_per_h_per_mmol_atp=0.010778,
            no_growth_ngam_mmol_atp_per_gdw_h=32.80,
            slope_linearity_growth_error_per_h=0.0),
        MaintenanceRegime(
            name="oxygen_capped_glucose_excess",
            glucose_lower_bound_mmol_per_gdcw_h=-10.0,
            oxygen_lower_bound_mmol_per_gdcw_h=-2.0,
            growth_at_zero_ngam_per_h=0.276688,
            ethanol_at_zero_ngam_mmol_per_gdcw_h=15.386,
            ngam_growth_slope_per_h_per_mmol_atp=0.010223,
            no_growth_ngam_mmol_atp_per_gdw_h=25.73,
            slope_linearity_growth_error_per_h=0.004912),
        MaintenanceRegime(
            name="oxygen_limited_glucose_excess",
            glucose_lower_bound_mmol_per_gdcw_h=-10.0,
            oxygen_lower_bound_mmol_per_gdcw_h=-1.0,
            growth_at_zero_ngam_per_h=0.236160,
            ethanol_at_zero_ngam_mmol_per_gdcw_h=16.298,
            ngam_growth_slope_per_h_per_mmol_atp=0.009759,
            no_growth_ngam_mmol_atp_per_gdw_h=22.86,
            slope_linearity_growth_error_per_h=0.000224),
    )
}

#: The default for anything in this repository that has to pick one, because it is the point
#: `fba/physiology.py` validates against and the only one with a measured phenotype behind
#: its bounds. It is NOT wired in as a function default anywhere: taking the regime is the
#: whole correction, and a default regime is a scalar with extra steps.
REFERENCE_REGIME = MAINTENANCE_REGIMES["aerobic_batch_reference"]


def regime_for(glucose_lower_bound: float, oxygen_lower_bound: float,
               tolerance: float = 1e-9) -> MaintenanceRegime:
    """The solved regime at these exchange bounds, or a refusal naming the ones that exist.

    EXACT lookup, never interpolation. The slope is not a smooth function of the bounds --
    it steps when overflow switches on -- so a value read off a line between two solved
    points would be wrong in the region that matters most, and wrong invisibly.

    Args:
        glucose_lower_bound: mmol/gDCW/h, negative for uptake.
        oxygen_lower_bound: mmol/gDCW/h, negative for uptake. ``-1000.0`` for unlimited.
        tolerance: Absolute match tolerance on both bounds.

    Raises:
        KeyError: naming every declared point and the script that adds another.
    """
    for regime in MAINTENANCE_REGIMES.values():
        if (abs(regime.glucose_lower_bound_mmol_per_gdcw_h - glucose_lower_bound) <= tolerance
                and abs(regime.oxygen_lower_bound_mmol_per_gdcw_h - oxygen_lower_bound)
                <= tolerance):
            return regime
    declared = "; ".join(
        f"{r.glucose_lower_bound_mmol_per_gdcw_h}/{r.oxygen_lower_bound_mmol_per_gdcw_h}"
        for r in MAINTENANCE_REGIMES.values())
    raise KeyError(
        f"REFUSED: no solved maintenance regime at glucose {glucose_lower_bound}, oxygen "
        f"{oxygen_lower_bound} mmol/gDCW/h. The NGAM slope steps when overflow switches on, "
        f"so it is not interpolated. Declared points (glucose/oxygen): {declared}. To add "
        f"one, put its bounds in scripts/maintenance_scale.py and rerun it")


def stress_maintenance_increment(atp_per_glucose: float | None = None) -> tuple[float, float]:
    """The measured maintenance difference between GSR present and absent, mmol ATP/gDW/h.

    Returns an INTERVAL, over the P/O range if no ratio is given and over the measurement's
    own standard deviations otherwise. An interval rather than a point because the two
    published values overlap within one standard deviation of each other -- 0.066 +/- 0.032
    against 0.082 +/- 0.045 -- so a point estimate of the difference would assert a precision
    the measurement does not have.

    **This is the number `latent_bridge._MAX_STRESS_MAINTENANCE` should be compared against,
    not silently replaced by**, for the sign reason in the module docstring. What it settles
    is the SCALE: whatever the sign, the effect is a few tenths of a mmol ATP/gDW/h, not six
    and a half. What it does NOT settle is what that costs in growth -- that is the regime's
    business, and :data:`MAINTENANCE_REGIMES` is where it lives.
    """
    difference = (MEASURED_MAINTENANCE["gsr_absent"].m_s_mmol_glucose_per_gdw_h
                  - MEASURED_MAINTENANCE["gsr_intact"].m_s_mmol_glucose_per_gdw_h)
    if atp_per_glucose is None:
        low, high = ATP_PER_GLUCOSE_AEROBIC
        return difference * low, difference * high
    return difference * atp_per_glucose, difference * atp_per_glucose
