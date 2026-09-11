"""S0: the mechanistic model synthetic cultures come from.

    mu(d)      = mu_max / (1 + (d/IC50)^m)
    k_synth(d) = basal + (peak-basal) * d^n/(EC50^n + d^n)
    dR/dt      = k_synth - (mu + k_deg) R

The third equation is why this generator is worth having: growth appears in the reporter's own
balance, so a culture whose promoter never moves still shows an apparent dose response. Data
without that term would validate the pipeline against a strawman.

The growth rate is also the carbon bookkeeping
----------------------------------------------
Those three equations treat the growth rate as a number that only dilutes the reporter, and a
yeast cell does not. How much biomass a gram of glucose buys, whether the culture makes ethanol
at all, and how much oxygen it draws are all functions of the growth rate, and all of them turn
over at one point: the critical rate. :func:`chemostat_physiology` reads that dependence out of
a measured curve instead of asserting it, which matters because every stress dose moves the
growth rate and therefore moves all of it.

The carrying capacity is therefore a function of the dose
---------------------------------------------------------
Until 2026-09-05 :func:`simulate_culture` integrated against ``params.carrying_capacity``, one
asserted number per strain that no dose could move, while :func:`substrate_limited_capacity`
sat beside it deriving the same quantity from the measured yield and was called by nothing but
its own test. The logistic now closes on the substrate budget by default: the dose lowers the
growth rate, the growth rate sets the measured biomass yield, and the yield and the glucose
charge set the capacity. That removes four asserted per-strain capacities from the default
panel and adds one declared medium quantity, :data:`STANDARD_GLUCOSE_G_PER_L`.

No state is added for the sugar, and none is needed. With ``dS/dt = -mu*X/Y`` the substrate
is eliminated exactly by conservation -- ``S = S0 - (X - X0)/Y`` -- and exhaustion lands at
``X = X0 + S0*Y``, which *is* the capacity. A glucose state would be an unmeasured ODE (no
assay in this corpus reads sugar in a well) standing in for an algebraic identity.

What it does **not** do is reproduce the deceleration a high dose actually shows on these
plates, and the sign is why. van Hoek's yield *rises* as the rate falls, so a dose that slows
a culture raises its capacity: for UPRE1 the derived K is 4.55 g/L at 0 mM and 9.77 g/L at
5 mM, where the asserted one was 3.0 at both. A logistic fitted freely to the real wells
wants the opposite -- `REVISED_BUILD_LIST.md` reports K_dosed/K_control = 0.245 at 2-5 mM --
so the fitted capacity drop is not what the measured yield curve predicts, and no glucose
charge makes it so.
`scripts/carrying_capacity_check.py` measures the residual and sweeps the charge.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from .. import paths
from ..mech.adapters import simulate_protocol
from ..reporter import ReporterKinetics, simulate_reporter

__all__ = [
    "CRITICAL_GROWTH_RATE_PER_H",
    "STANDARD_GLUCOSE_G_PER_L",
    "ChemostatPhysiology",
    "CultureParameters",
    "chemostat_physiology",
    "measured_growth_rate_range",
    "resolved_capacity",
    "simulate_culture",
    "simulate_empirical_comparison",
    "simulate_protocol",
    "substrate_limited_capacity",
]

CRITICAL_GROWTH_RATE_PER_H = 0.28
"""Specific growth rate at which aerobic glucose-limited *S. cerevisiae* starts to ferment.

**Measured.** van Hoek, van Dijken & Pronk 1998, PMID 9797269, "Effect of specific growth rate
on fermentative capacity of baker's yeast", Table 1. Ethanol is below the detection limit at
every dilution rate up to and including 0.25 /h, first appears at D = 0.28 (0.11 mmol/gDW/h)
and reaches 13.9 by D = 0.40. Strain DS28911, aerobic glucose-limited chemostat, 30 degrees,
dissolved oxygen held above 60% of air saturation.

Written down as a named constant because it is a *threshold in the growth rate*, and the growth
rate is what every stressor in this generator moves. A dose that takes an untreated culture from
0.40 /h to 0.25 /h has not merely slowed it: it has switched it from respiro-fermentative to
fully respiratory, roughly doubled its biomass yield and *raised* its oxygen uptake. Nothing in
the reporter equations knows that, which is why it is exposed here."""

STANDARD_GLUCOSE_G_PER_L = 20.0
"""Glucose charged into one well of the standard medium, g/L.

**ASSERTED, and never assayed on these plates.** 2% w/v is what
:class:`generator.context.CultureContext` already defaults to and what
``docs/ARCHITECTURE.md`` flags as "glucose is a literal, never assayed". Nothing in the
committed corpus records the charge of the four biosensor plates, so this is a declared
axis rather than a measurement, and `scripts/carrying_capacity_check.py` sweeps it rather
than trusting it.

It is a *substrate budget*, not a rate: it never enters ``mu``. It sets how much biomass the
well can make, through :func:`substrate_limited_capacity`, which is the same division of
labour ``context.py`` documents.

The yield it is multiplied by carries a regime caveat that the derived capacity inherits.
van Hoek's curve is a glucose-**limited** chemostat; a well charged with 20 g/L is
glucose-**excess** for the whole of a 4.14 h read, and ``docs/PROTOCOLS.md`` records
measured excess yields 5.6-6.0x below the carbon-limited value at the same rate. So the
derived capacity is an upper bound on what an excess well can reach, and the direction of
that bias is known even though its size here is not."""


def _hill(dose: float, half: float, exponent: float) -> float:
    """Hill saturation ``d^n / (K^n + d^n)``, computed without overflowing.

    Written directly, the intermediate ``d**n`` exceeds a float long before the
    result would: at d=5, n=1000 it overflows even though the answer is exactly 1.
    In log space the same expression is the logistic of ``n*(ln d - ln K)``, which
    saturates cleanly at both ends. Augmentation draws exponents wide enough to
    reach this, so it is not a hypothetical.
    """
    d = max(float(dose), 0.0)
    if d <= 0.0:
        return 0.0
    return float(expit(exponent * (np.log(d) - np.log(half))))


@dataclass(frozen=True)
class ChemostatPhysiology:
    """What a culture is doing metabolically at one growth rate.

    Every field is read from van Hoek 1998 Table 1 (PMID 9797269) and interpolated between
    measured dilution rates. Fluxes are mmol per gram dry weight per hour, exactly as that
    table's footnote a states them; ``biomass_yield`` is grams dry biomass per gram glucose.

    Args:
        growth_rate: The rate this was evaluated at, 1/h. In the source it is the dilution
            rate of a steady-state chemostat, which equals the specific growth rate.
        glucose_uptake: q_glucose.
        oxygen_uptake: q_O2. **Not monotonic in growth rate** -- it peaks at 7.4 at the
            critical rate and falls to 3.7 at 0.40 as the cells switch to fermentation.
        co2_production: q_CO2.
        ethanol: q_ethanol. Zero below the critical rate.
        glycerol: q_glycerol. Zero below 0.35 /h.
        acetate: q_acetate.
        pyruvate: q_pyruvate.
        biomass_yield: Y_sx, g/g. 0.45-0.49 while respiratory, collapsing to 0.20 at 0.40 /h.
        interpolated: False when the growth rate landed exactly on a measured row.
    """

    growth_rate: float
    glucose_uptake: float
    oxygen_uptake: float
    co2_production: float
    ethanol: float
    glycerol: float
    acetate: float
    pyruvate: float
    biomass_yield: float
    interpolated: bool

    @property
    def fermentative(self) -> bool:
        """Whether the culture is past the Crabtree switch.

        Recorded per row by the environment sweep because a downstream fit conditioned on
        the growth rate alone cannot see it: the same rate means two different metabolisms
        on either side of the threshold.

        Keyed on the measured critical rate, not on ``ethanol > 0``. Interpolating between
        van Hoek's non-detect at D = 0.25 and his 0.11 at D = 0.28 puts a positive ethanol
        on every rate above 0.2501, so the boolean turned True a full 0.03 /h before the
        threshold it is supposed to mark. In the default sweep that shipped 288 of 8364
        rows carrying ``fermentative=True`` beside ``above_critical_growth_rate=False``,
        both columns in the same CSV, disagreeing.

        A non-detect is not a zero to interpolate through. It says the quantity was below
        the detection limit, and a line drawn from it invents a measurement the instrument
        did not make.
        """
        return self.growth_rate > CRITICAL_GROWTH_RATE_PER_H

    @property
    def respiratory_quotient(self) -> float:
        """q_CO2 / q_O2. Near 1 while respiratory, and rises steeply once ethanol appears.

        Derived rather than tabulated, so it carries no extra provenance of its own -- but it
        is the quantity a wet lab actually watches to see the switch happen, so it is worth
        having next to the fluxes it comes from.
        """
        if self.oxygen_uptake <= 0.0:
            raise ValueError(
                f"no oxygen uptake at growth rate {self.growth_rate:g}/h, so the respiratory "
                "quotient is undefined; this dataset is aerobic throughout and should never "
                "produce it"
            )
        return self.co2_production / self.oxygen_uptake


_CHEMOSTAT_COLUMNS = (
    ("GlucoseUptake", "glucose_uptake"),
    ("O2uptake", "oxygen_uptake"),
    ("CO2production", "co2_production"),
    ("Ethanol", "ethanol"),
    ("Glycerol", "glycerol"),
    ("Acetate", "acetate"),
    ("Pyruvate", "pyruvate"),
    ("BiomassYield_g_per_g", "biomass_yield"),
)


@functools.lru_cache(maxsize=1)
def _chemostat_table() -> dict[str, np.ndarray]:
    """The vendored van Hoek table, as columns.

    Read rather than transcribed into this file on purpose. The numbers are data, the file is
    where data lives, and `docs/research/CHEMOSTAT_REFERENCE.md` records the transcription
    check that was run against the paper's own carbon balance. Copying them into source would
    put a second, unchecked copy one edit away from disagreeing with the first.
    """
    path = paths.data_dir() / "physiology" / "chemostatData_VanHoek1998.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"the van Hoek 1998 chemostat table is not at {path}; growth-rate-dependent "
            "physiology is read from it and cannot be substituted with a plausible curve"
        )
    rows = [line.split("\t") for line in path.read_text().strip().splitlines()]
    header, body = rows[0], rows[1:]
    index = {name: header.index(name) for name in
             ("Drate", *(source for source, _ in _CHEMOSTAT_COLUMNS))}
    table = {name: np.array([float(row[column]) for row in body], dtype=float)
             for name, column in index.items()}
    order = np.argsort(table["Drate"])
    return {name: values[order] for name, values in table.items()}


def measured_growth_rate_range() -> tuple[float, float]:
    """Lowest and highest growth rate the chemostat curve was measured at, 1/h."""
    rates = _chemostat_table()["Drate"]
    return float(rates[0]), float(rates[-1])


def chemostat_physiology(growth_rate: float) -> ChemostatPhysiology:
    """Measured aerobic physiology at one specific growth rate.

    Linear interpolation between the measured dilution rates, and **no extrapolation**: a rate
    outside the range van Hoek ran raises rather than returning the nearest edge. The temptation
    to clamp is strong because a stressed culture in this generator goes to zero growth, and
    clamping is exactly the wrong answer -- it would report a 0.45 g/g biomass yield for a
    culture nobody has measured, and 0.45 g/g is a *fact about a chemostat at 0.025 /h*, not
    about a poisoned well.

    Args:
        growth_rate: Specific growth rate, 1/h.

    Raises:
        ValueError: if the rate is outside the measured range. The message names the range,
            because the caller's options are to stay inside it or to measure more of it.
    """
    rate = float(growth_rate)
    rates = _chemostat_table()["Drate"]
    low, high = float(rates[0]), float(rates[-1])
    if not low <= rate <= high:
        raise ValueError(
            f"growth rate {rate:g}/h is outside the measured chemostat curve [{low:g}, "
            f"{high:g}] /h (van Hoek 1998, PMID 9797269, Table 1); this generator will not "
            "extrapolate a biomass yield or a fermentation flux beyond what was measured"
        )
    table = _chemostat_table()
    values = {field: float(np.interp(rate, rates, table[source]))
              for source, field in _CHEMOSTAT_COLUMNS}
    return ChemostatPhysiology(
        growth_rate=rate,
        interpolated=not bool(np.any(np.isclose(rates, rate, rtol=0.0, atol=1e-12))),
        **values,
    )


def substrate_limited_capacity(
    initial_biomass: float, glucose_g_per_L: float, growth_rate: float
) -> float:
    """Biomass a well can reach on the glucose it was given, g/L.

    ``initial + glucose * Y_sx(mu)``, with the yield taken from the measured curve rather than
    asserted. This is the term that carries the Crabtree switch into an observable: at 0.40 /h
    the yield is 0.20 g/g and at 0.25 /h it is 0.48, so a stressor that slows a culture across
    the threshold makes it end up **denser**, not sparser, on the same glucose. A generator
    with a fixed ``carrying_capacity`` says the opposite, and the sign of that effect is
    visible in the OD channel of every well on a plate that runs to substrate exhaustion.

    :func:`simulate_culture` integrates against this by default since 2026-09-05. Note the
    scope the paragraph above sets and the 4.14 h biosensor protocol does not meet: from OD
    0.16 a well reaches under a tenth of the derived capacity, so the term shifts that read's
    late growth rate by 0.0144 /h, only 1.23x the 0.0117 /h growth assay SE. On a 24 h read
    it is 5.64x. `scripts/carrying_capacity_check.py` is the measurement.

    Args:
        initial_biomass: Inoculum, g/L.
        glucose_g_per_L: Glucose supplied. A standard 2% w/v medium is 20.0.
        growth_rate: Realised specific growth rate, 1/h, which sets the yield.

    Raises:
        ValueError: through :func:`chemostat_physiology` if the rate is unmeasured, and
            directly for a non-positive glucose charge.
    """
    if glucose_g_per_L <= 0:
        raise ValueError(f"glucose_g_per_L must be positive, got {glucose_g_per_L}")
    if initial_biomass <= 0:
        raise ValueError(f"initial_biomass must be positive, got {initial_biomass}")
    yield_g_per_g = chemostat_physiology(growth_rate).biomass_yield
    return float(initial_biomass) + float(glucose_g_per_L) * yield_g_per_g


def resolved_capacity(
    params: "CultureParameters",
    growth_rate: float,
    initial_biomass: float,
    glucose_g_per_L: float | None = STANDARD_GLUCOSE_G_PER_L,
) -> tuple[float, str]:
    """The capacity one well is integrated against, and where that number came from.

    Two routes, and which one ran is returned beside the number because they disagree by
    more than a factor of two across the Crabtree switch and a reader cannot tell from the
    value alone. ``scripts/run_environment_sweep.py`` already recorded exactly this pair
    per row; this is the same bookkeeping moved to where the choice is actually made.

    The substrate route is the default. It falls back to the asserted capacity **only**
    where the realised rate is outside the range van Hoek measured, because that is where
    :func:`chemostat_physiology` refuses to supply a yield -- and a refusal is not a licence
    to invent one. The fallback is reported, never silent. It is not rare: the default
    panel sits inside the measured range but the real plates' own fitted rates do not, with
    a median ``mu_max`` of 0.567 /h against a measured ceiling of 0.40 /h.

    Args:
        params: Strain kinetics; supplies the asserted capacity used by the opt-out.
        growth_rate: Realised specific growth rate for this well, 1/h, after the dose and
            any nutrient factor. This is what makes the capacity dose-dependent.
        initial_biomass: Inoculum, g/L.
        glucose_g_per_L: Substrate charge. ``None`` opts out and takes
            ``params.carrying_capacity`` instead -- the right answer for a well that is
            capacity-limited rather than substrate-limited, and for any caller that has
            already resolved a capacity of its own.

    Returns:
        ``(capacity_g_per_L, source)``.
    """
    if glucose_g_per_L is None:
        return float(params.carrying_capacity), "asserted per strain, substrate route declined"
    # Raise on a bad charge or inoculum; only a refused *rate* is allowed to fall back.
    if glucose_g_per_L <= 0 or initial_biomass <= 0:
        raise ValueError(
            f"a substrate-limited well needs a positive charge and inoculum, got "
            f"glucose_g_per_L={glucose_g_per_L} and initial_biomass={initial_biomass}; "
            "pass glucose_g_per_L=None to take the asserted capacity instead"
        )
    try:
        capacity = substrate_limited_capacity(
            initial_biomass, glucose_g_per_L, growth_rate)
    except ValueError:
        return (float(params.carrying_capacity),
                f"asserted per strain, growth rate {growth_rate:g}/h outside the measured "
                f"yield curve {measured_growth_rate_range()}")
    return capacity, "measured yield, van Hoek 1998 PMID 9797269"


@dataclass(frozen=True)
class CultureParameters:
    """Per-strain kinetics. Rates in 1/h, doses in mM, biomass in g/L.

    Args:
        mu_max: Specific growth rate with no stressor.
        growth_ic50: Dose halving the growth rate.
        growth_hill: Steepness of the growth-inhibition curve.
        promoter_basal: Promoter activity with no stressor.
        promoter_peak: Promoter activity at saturating dose.
        promoter_ec50: Dose giving half-maximal induction.
        promoter_hill: Steepness of the induction curve.
        carrying_capacity: Biomass the culture saturates at. Asserted per strain, and **no
            longer what :func:`simulate_culture` integrates against by default**: the
            capacity is derived from the glucose charge and the measured biomass yield at
            the realised rate, so a dose moves it. This field is the declared opt-out --
            pass ``glucose_g_per_L=None`` for a well that is capacity-limited rather than
            substrate-limited -- and the fallback wherever the realised rate lands outside
            the yield curve van Hoek measured. The two disagree by more than a factor of
            two across the Crabtree switch.
        k_deg: Reporter loss beyond dilution. Zero for a stable FP, which is what
            D2 found for mCitrine on the real plates.
        k_mat: Chromophore maturation rate, 1/h. ``None`` means instantaneous, which is what
            every simulated plate in this repository has assumed and **is not a measurement**.
            Maturation is a real delay between transcription and signal: Nagai 2002,
            PMID 11753368, identifies chromophore oxidation as the rate-limiting step of YFP
            maturation and reports that a single F46L substitution greatly accelerates it, so
            the rate is a property of the exact variant and cannot be read off "YFP". It has
            never been measured for the mCitrine constructs on these plates. Left at ``None``
            so no unmeasured number enters by default; supply it to see what the assumption
            costs. It matters most in the short reads this project runs -- an unmatured
            reporter is invisible, and a 4.14 h read gives it little time to arrive.
    """

    mu_max: float
    growth_ic50: float
    growth_hill: float
    promoter_basal: float
    promoter_peak: float
    promoter_ec50: float
    promoter_hill: float
    carrying_capacity: float
    k_deg: float = 0.0
    k_mat: float | None = None

    def __post_init__(self) -> None:
        positive = {
            "mu_max": self.mu_max, "growth_ic50": self.growth_ic50,
            "growth_hill": self.growth_hill, "promoter_ec50": self.promoter_ec50,
            "promoter_hill": self.promoter_hill,
            "carrying_capacity": self.carrying_capacity,
        }
        bad = sorted(k for k, v in positive.items() if v <= 0)
        if bad:
            raise ValueError(f"parameters must be positive: {', '.join(bad)}")
        if self.promoter_basal < 0 or self.promoter_peak < 0 or self.k_deg < 0:
            raise ValueError("promoter activities and k_deg must be non-negative")
        if self.k_mat is not None and self.k_mat <= 0:
            raise ValueError(
                f"k_mat must be positive, or None for instantaneous maturation, got "
                f"{self.k_mat}"
            )

    def growth_rate_at(self, dose: float) -> float:
        """Unrestricted growth rate at a given stressor dose."""
        return self.mu_max * (1.0 - _hill(dose, self.growth_ic50, self.growth_hill))

    def promoter_activity_at(self, dose: float) -> float:
        """Promoter output at a given stressor dose."""
        saturation = _hill(dose, self.promoter_ec50, self.promoter_hill)
        return self.promoter_basal + (self.promoter_peak - self.promoter_basal) * saturation


def simulate_empirical_comparison(
    times_h: np.ndarray,
    dose: float,
    params: CultureParameters,
    initial_biomass: float,
    initial_reporter: float | None = None,
    nutrient_factor: float = 1.0,
    glucose_g_per_L: float | None = STANDARD_GLUCOSE_G_PER_L,
) -> dict[str, np.ndarray]:
    """Simulate one well.

    Args:
        times_h: Ascending time grid, hours.
        dose: Stressor concentration, mM.
        params: Strain kinetics.
        initial_biomass: Starting biomass, g/L.
        initial_reporter: Starting per-cell reporter. Defaults to the quasi-steady
            state for the *untreated* culture, which is what a culture carries when
            the stressor is added at t=0.
        nutrient_factor: Scales growth without touching promoter activity, standing in
            for nutrient limitation. This is the knob that decouples the two axes.
        glucose_g_per_L: Substrate charged into the well. The capacity is derived from it
            and the measured yield at the realised rate, so the dose moves the capacity.
            ``None`` is the explicit opt-out and takes ``params.carrying_capacity``
            instead; use it for a well that is capacity-limited rather than
            substrate-limited, or when the caller has already resolved a capacity.

    Returns:
        ``biomass``, ``growth_rate``, ``promoter_activity`` and ``reporter``
        (per-cell concentration), each on the supplied grid, plus the scalar
        ``carrying_capacity`` this well was integrated against and the
        ``capacity_source`` that produced it.
    """
    if not 0.0 < nutrient_factor <= 1.0:
        raise ValueError(f"nutrient_factor must be in (0, 1], got {nutrient_factor}")
    t = np.asarray(times_h, dtype=float)
    mu_unrestricted = params.growth_rate_at(dose) * nutrient_factor
    k_synth = params.promoter_activity_at(dose)
    capacity, capacity_source = resolved_capacity(
        params, mu_unrestricted, initial_biomass, glucose_g_per_L)

    def rhs(_, y):
        return [mu_unrestricted * y[0] * (1.0 - y[0] / capacity)]

    solution = solve_ivp(
        rhs, (t[0], t[-1]), [initial_biomass], t_eval=t,
        method="LSODA", rtol=1e-9, atol=1e-12,
    )
    if not solution.success:
        raise RuntimeError(f"culture integration failed: {solution.message}")
    biomass = solution.y[0]

    # Realised rate dilutes the reporter; a saturating culture dilutes less.
    growth = mu_unrestricted * (1.0 - biomass / capacity)

    if initial_reporter is None:
        untreated = params.growth_rate_at(0.0) * nutrient_factor
        initial_reporter = params.promoter_basal / max(untreated + params.k_deg, 1e-9)

    reporter = simulate_reporter(
        t, np.full_like(t, k_synth), growth,
        ReporterKinetics(k_deg=params.k_deg, k_mat=params.k_mat), r0=initial_reporter,
    )
    return {
        "biomass": biomass,
        "growth_rate": growth,
        "promoter_activity": np.full_like(t, k_synth),
        "reporter": reporter,
        "carrying_capacity": capacity,
        "capacity_source": capacity_source,
        "simulation_route": "empirical_culture_comparison",
    }


simulate_culture = simulate_empirical_comparison
