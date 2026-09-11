"""The physiological baseline a stress is added to, which every simulated culture shared.

Training varied the stressor and nothing else, so every well was exponential-phase glucose
at 30 degrees. Two of the twenty-four modules are set by the baseline before any agent is
added -- the carbon regulon is derepressed the moment glucose runs out, and the general
stress response is high in stationary phase whatever else is happening -- so a model
trained in that one corner reads a late-plate ESR as stress.

Context shifts the baseline and the growth rate, not the response. What a stressor does is
still its own; what changes is what it is added to, and how fast the culture is diluting
its reporter while being measured.

The growth rate is not the end of it
------------------------------------
Two of the terms below used to be plausible curves and are now measurements, because the
difference was large. Temperature was a symmetric line through 30 degrees, which gave a
culture at 45 degrees 40% of maximum growth when the measured maximum growth temperature is
45.4; it is now Rosso's cardinal model at Salvado's cardinals. Anoxia was the same line
extrapolated, which put an unsupplemented anaerobic glucose culture at 0.28 /h when the answer
is zero, because it cannot make a sterol without oxygen.

pH is carried here and its effect on growth is refused -- see :func:`ph_growth_factor`. The pH
effect that *is* measured is on the reporter, not on the culture, and it is
:func:`ystwin.photophysics.citrine_ph_response`, which this module re-exports rather than
owns: a fluorophore's titration curve is not a property of the culture, and `observation.py`
needs to reach it without importing the generator.

`docs/research/GENERATOR_REALISM.md` is the audit these changes came out of and lists what is
still missing.
"""

from __future__ import annotations

from dataclasses import dataclass

# Re-exported, not owned. Moved to a depth-0 leaf so `observation.py` can reach the
# measurement channel without importing the generator; the names stay for existing callers.
from ..photophysics import CITRINE_PKA, citrine_ph_response
from .culture import CRITICAL_GROWTH_RATE_PER_H, chemostat_physiology
from .stress_panel import MODULES

__all__ = [
    "ANAEROBIC_MU_MAX_PER_H",
    "CARDINAL_TEMPERATURES_C",
    "CITRINE_PKA",
    "CRITICAL_GROWTH_RATE_PER_H",
    "CultureContext",
    "baseline_activity",
    "citrine_ph_response",
    "context_growth_rate",
    "context_physiology",
    "ph_growth_factor",
    "temperature_factor",
]

_CARBON = {
    "glucose": {"mu": 0.40, "carbon": 0.0, "retrograde": 0.0, "respiring": 0.3},
    "galactose": {"mu": 0.28, "carbon": 0.05, "retrograde": 0.15, "respiring": 0.6},
    "ethanol": {"mu": 0.14, "carbon": 0.80, "retrograde": 0.40, "respiring": 1.0},
}
"""Growth rates and derepression per carbon source. Glucose represses the carbon regulon
through Mig1 and is fermented; ethanol must be respired, which needs the mitochondrion and
raises the Rtg1/Rtg3 arm with it.

Galactose derepresses this regulon barely at all, which "not glucose" makes tempting to
assume otherwise. The axis is fermentable against non-fermentable, not one sugar against
another: ADH2 is repressed by galactose as well as by glucose, a cat8 null grows normally on
galactose while failing on ethanol, lactate and acetate, and ICL1 moves more than 200-fold
only on transfer to a non-fermentative source (Schuller & Entian 1996, PMID 8710508;
Walther & Schuller 2001, doi 10.1099/00221287-147-8-2037).

This table is a static label and a batch culture is not. A well started on glucose ferments
it, then switches to the ethanol it made, and arrives at the strongest derepressing state on
the plate -- which is why a late reading of a glucose well can exceed a galactose one.

The three `mu` values are **asserted**, and `literature.py` carries the glucose one with
`verified=False` for a reason worth repeating here: van Hoek 1998 (PMID 9797269) quotes a batch
mu_max on glucose of 0.42 /h and cites that to another paper rather than measuring it, so 0.40
is a round number near a figure nobody in the chain measured. It is left alone because changing
it would move every existing result, and it is named because 0.40 sits above the 0.28 /h
critical rate and therefore decides which side of the Crabtree switch the reference culture is
on. `respiring` is asserted too and does no work at full aeration."""

_PHASE = {
    "exponential": {"mu": 1.0, "ESR": 0.0},
    "diauxic": {"mu": 0.35, "ESR": 0.45},
    "stationary": {"mu": 0.05, "ESR": 0.75},
}
"""Gasch 2000: the environmental stress response is the signature of leaving exponential
growth, and is high in stationary phase with no stressor present at all."""

_REFERENCE_TEMPERATURE = 30.0
_HEAT_PER_DEGREE = 0.06
_HYPOXIA_THRESHOLD = 0.21

CARDINAL_TEMPERATURES_C = {"min": 2.8, "opt": 32.3, "max": 45.4}
"""Where *S. cerevisiae* growth starts, peaks and stops, in degrees Celsius.

**Measured.** Salvado et al. 2011, PMID 21317255, "Temperature adaptation markedly determines
evolution within the genus Saccharomyces", Table 2. Ten *S. cerevisiae* strains from ten
independent origins, grown in YNB plus 20 g/L glucose in **96-well microplates read at OD600**
at eleven temperatures from 4 to 46 degrees -- the same medium, format and channel this
generator simulates, which is why this paper rather than a fermenter study.

`opt` and `max` are the species means the paper states in its own abstract, 32.3 and 45.4.
`min` is the mean of the same ten Table 2 rows (0.74, 4.11, 5.04, 4.38, 4.02, 0.76, 3.55, 0.43,
3.72, 1.69), which comes to 2.84; the abstract does not state a species mean for T_min, and the
paper reports no significant difference between strains in that column. Averaging the other two
columns the same way reproduces 32.27 and 45.40, which is what licenses reading `min` off the
table by the same arithmetic.

**This replaces `1 - 0.04 * |T - 30|`.** That was symmetric and linear, and yeast is neither.
The old form said 20 degrees and 40 degrees cost the same 40%; the measured curve says 20
degrees is a mild slowdown and 40 degrees is most of the way to a dead culture. The two
disagree most exactly where a heat-shock experiment is run."""

ANAEROBIC_MU_MAX_PER_H = 0.31
"""Fastest growth of *S. cerevisiae* with no oxygen, in a medium that supports it.

**Measured.** Verduyn et al. 1990, PMID 1975265, "Physiology of Saccharomyces cerevisiae in
anaerobic glucose-limited chemostat cultures": CBS 8066 in mineral medium *supplemented with
ergosterol and Tween 80*, mu_max 0.31 /h, and a maximal biomass yield of 0.10 g/g at D = 0.10.

Two facts hide in that sentence and both are load-bearing here.

**The supplements are not a detail.** Sterol and unsaturated fatty acid synthesis both need
molecular oxygen, so a defined medium without them does not support anaerobic growth at all --
Andreasen & Stier 1953, PMID 13034889, "Anaerobic nutrition of Saccharomyces cerevisiae. I.
Ergosterol requirement for growth in a defined medium". A sweep that dials oxygen to zero and
reports a growth rate has silently assumed a supplemented medium; `CultureContext` therefore
makes that assumption an explicit field rather than a default.

**The yield collapses long before the rate does.** 0.31 /h anaerobic against 0.42 /h aerobic is
a 26% loss of rate; 0.10 g/g against van Hoek's 0.48 g/g at the same D = 0.10 is a **4.8-fold**
loss of yield. Any model that represents oxygen limitation as a growth-rate multiplier alone
gets the small effect right and the large one wrong."""


def temperature_factor(temperature_c: float) -> float:
    """Growth-rate multiplier from temperature, 1.0 at the 30-degree reference.

    The cardinal temperature model with inflection of Rosso, Lobry & Flandrois 1993,
    PMID 8412234, "An unexpected correlation between cardinal temperatures of microbial growth
    highlighted by a new model", evaluated at :data:`CARDINAL_TEMPERATURES_C` and then divided
    by its own value at 30 degrees.

    The normalisation is the point of care. The carbon-source growth rates in this module were
    read from 30-degree cultures, so the temperature term has to be a *shape* that leaves 30
    degrees alone; taking the model's own mu_opt would double-count the rate and silently
    rescale every existing result. Dividing by CTMI(30) makes this exactly a shape, and the
    shape is what the citation supports.

    Returns 0.0 outside (T_min, T_max), which is what those parameters mean: Salvado defines
    T_min as the temperature below which growth is no longer observed and T_max as the one
    above which no growth occurs.
    """
    t = float(temperature_c)
    return _ctmi(t) / _ctmi(_REFERENCE_TEMPERATURE)


def _ctmi(temperature_c: float) -> float:
    """Rosso's CTMI, normalised to 1 at T_opt. Zero outside the cardinal range."""
    t = float(temperature_c)
    t_min = CARDINAL_TEMPERATURES_C["min"]
    t_opt = CARDINAL_TEMPERATURES_C["opt"]
    t_max = CARDINAL_TEMPERATURES_C["max"]
    if not t_min < t < t_max:
        return 0.0
    numerator = (t - t_max) * (t - t_min) ** 2
    denominator = (t_opt - t_min) * (
        (t_opt - t_min) * (t - t_opt) - (t_opt - t_max) * (t_opt + t_min - 2.0 * t)
    )
    return float(numerator / denominator)


_FERMENTABLE = {"glucose": True, "galactose": True, "ethanol": False}
"""Whether the culture can make ATP from this substrate with no oxygen at all.

Fermentation of a sugar to ethanol yields ATP by substrate-level phosphorylation and needs no
terminal electron acceptor; ethanol as a *substrate* has to be oxidised, so a culture on it
stops dead without oxygen. This is the axis the carbon table's `respiring` column was already
reaching for, written as the boolean it actually is."""


@dataclass(frozen=True)
class CultureContext:
    """The conditions a culture is in, before any stressor is applied.

    Args:
        carbon_source: Which of the tabulated substrates the medium carries.
        growth_phase: Where the culture is on its own growth curve.
        temperature_c: Incubator temperature. Enters through :func:`temperature_factor`.
        oxygen: Gas-phase oxygen fraction. 0.21 is air.
        anaerobic_supplements: Whether the medium carries ergosterol and an unsaturated fatty
            acid. **Not cosmetic.** Sterol and unsaturated fatty acid synthesis both consume
            molecular oxygen, so a defined medium without them supports no anaerobic growth at
            all (Andreasen & Stier 1953, PMID 13034889), and the anaerobic chemostats that
            measured a rate were supplemented (Verduyn 1990, PMID 1975265). Defaults to False
            because a plate of ordinary synthetic medium does not carry them, so a sweep that
            drops oxygen to zero gets a dead culture unless it says otherwise.
        glucose_g_per_L: Sugar charged into the well. Recorded rather than used by the growth
            rate: it sets how much biomass the well can make, through
            ``culture.substrate_limited_capacity``, not how fast. 20.0 is a standard 2% w/v.
        ph_medium: pH of the medium. **Recorded, not modelled.** See
            :func:`ph_growth_factor` for why this module refuses to turn it into a rate.
        ph_cytosolic: Intracellular pH. This is the one the reporter sees, through
            :func:`citrine_ph_response`, and it is not the medium pH -- yeast holds a cytosol
            near neutral against a medium at 4.
        strain: Label carried through to the output.
    """

    carbon_source: str = "glucose"
    growth_phase: str = "exponential"
    temperature_c: float = 30.0
    oxygen: float = 0.21
    anaerobic_supplements: bool = False
    glucose_g_per_L: float = 20.0
    ph_medium: float = 5.5
    ph_cytosolic: float = 7.0
    strain: str = "BY4741"

    def __post_init__(self) -> None:
        if self.carbon_source not in _CARBON:
            raise ValueError(
                f"no physiology for carbon source {self.carbon_source!r}; "
                f"have {sorted(_CARBON)}")
        if self.growth_phase not in _PHASE:
            raise ValueError(
                f"no physiology for growth phase {self.growth_phase!r}; have {sorted(_PHASE)}")
        if self.oxygen < 0.0:
            raise ValueError(f"oxygen fraction must be non-negative, got {self.oxygen}")
        if self.oxygen > 1.0:
            raise ValueError(
                f"oxygen is a gas-phase mole fraction, so it cannot exceed 1.0, got "
                f"{self.oxygen}. Air is 0.21 and pure oxygen is 1.0")
        if self.oxygen > _HYPOXIA_THRESHOLD:
            raise ValueError(
                f"oxygen {self.oxygen:g} is enriched above air ({_HYPOXIA_THRESHOLD:g}) and "
                "no physiology here covers enrichment. The growth rate would be right -- "
                "respiration is already saturated at air, so the aeration term returns the "
                "air value for anything at or above it -- but returning it would assert that "
                "enrichment is inert, and it is not: raised pO2 raises ROS and the oxidative "
                "burden that the stress panel exists to represent. Sweep at or below air, or "
                "add a measured enrichment response"
            )
        for name in ("ph_medium", "ph_cytosolic"):
            value = getattr(self, name)
            if not 0.0 < value < 14.0:
                raise ValueError(f"{name} must be between 0 and 14, got {value}")
        if self.glucose_g_per_L < 0.0:
            raise ValueError(
                f"glucose_g_per_L must be non-negative, got {self.glucose_g_per_L}")

    @property
    def anoxic(self) -> bool:
        """Whether there is no usable oxygen at all."""
        return self.oxygen <= 0.0


def _aerated_carbon_rate(context: CultureContext) -> float:
    """Carbon-source growth rate after aeration, before phase and temperature, 1/h.

    Two regimes, and only one of them is asserted.

    **Reduced aeration, 0 < O2 < air.** Unchanged from before: a linear shortfall in the gas
    fraction, weighted by how much of the substrate has to be respired. Nothing measures this
    and it is very probably too gradual -- yeast's affinity for oxygen is high enough that
    respiration runs almost to anoxia, so the true curve should stay flat and then fall off a
    cliff. It is left alone because no verified half-saturation constant for *S. cerevisiae*
    oxygen uptake was found in this pass, and swapping one invented shape for another buys
    nothing. **Asserted.**

    **Anoxia, O2 == 0.** Now a measurement, and it is a different number from what the
    asserted curve extrapolates to. Zero on a non-fermentable substrate, which needs no
    citation beyond the definition of respiration. On a fermentable one, zero without
    ergosterol and unsaturated fatty acid -- Andreasen & Stier 1953, PMID 13034889 -- and
    otherwise capped at :data:`ANAEROBIC_MU_MAX_PER_H`, Verduyn 1990's measured 0.31 /h. That
    is applied as a **ceiling** rather than as a value, because the measurement is glucose and
    this table also holds galactose, for which no anaerobic rate was found.

    The two regimes do not meet, and the gap is left visible on purpose. The asserted curve
    carries an unsupplemented anaerobic glucose culture to 0.28 /h; the answer is zero, because
    it cannot make a sterol. A culture that cannot grow and a culture growing at two thirds of
    maximum is not a small modelling difference, and smoothing the join would hide precisely
    the thing the endpoint was added to say. A measured oxygen half-saturation constant would
    close it properly.
    """
    aerobic = _CARBON[context.carbon_source]["mu"]
    if context.anoxic:
        if not _FERMENTABLE[context.carbon_source] or not context.anaerobic_supplements:
            return 0.0
        return min(ANAEROBIC_MU_MAX_PER_H, aerobic)
    shortfall = max(0.0, 1.0 - context.oxygen / _HYPOXIA_THRESHOLD)
    return aerobic * (1.0 - _CARBON[context.carbon_source]["respiring"] * shortfall)


def ph_growth_factor(ph: float, cardinal_ph: dict[str, float] | None = None) -> float:
    """Growth-rate multiplier from pH, and a refusal when nothing measures it.

    The model is the cardinal pH model of Rosso, Lobry, Bajard & Flandrois 1995, PMID 16534932,
    "Convenient Model To Describe the Combined Effects of Temperature and pH on Microbial
    Growth", which also supplies the structure this module relies on elsewhere: that paper's
    combined model is built "using the hypothesis that the temperature and pH effects on the
    mu_max are independent", which is the licence for multiplying the terms together rather
    than fitting an interaction.

    **No cardinal pH values for *S. cerevisiae* were found in this pass**, and one relevant
    negative result argues against guessing them: Arroyo-Lopez et al. 2009, PMID 19246112, fit
    response surfaces to *S. cerevisiae* T73, *S. kudriavzevii* and their hybrid and found pH
    significant for the hybrid and for *S. kudriavzevii* but **not** for *S. cerevisiae*, where
    temperature dominated. So the honest position is that the effect over an ordinary
    experimental range is small and unpinned, not that it is a particular curve.

    There is also a modelling reason to distrust a bare pH knob. Below about pH 5 what damages
    yeast is the *undissociated* weak acid that low pH creates -- Verduyn 1990 measured exactly
    this, with acetate and propionate added to anaerobic chemostats lowering the biomass yield
    and raising ethanol production through proton import. Medium pH and weak-acid concentration
    are not separable knobs, so a one-argument pH term would be modelling the wrong variable.

    Args:
        ph: Medium pH.
        cardinal_ph: ``{"min": ..., "opt": ..., "max": ...}``, measured for the organism and
            medium in hand. Required.

    Raises:
        ValueError: when ``cardinal_ph`` is absent. The alternative is to return 1.0, which
            would quietly assert that pH does not matter in every sweep that varies it.
    """
    if cardinal_ph is None:
        raise ValueError(
            f"no growth-rate response to pH is available for this organism (asked for pH "
            f"{ph:g}); supply measured cardinal pH values as {{'min','opt','max'}} to use the "
            "Rosso 1995 cardinal pH model (PMID 16534932). This module records pH and refuses "
            "to invent its effect on growth"
        )
    low, opt, high = cardinal_ph["min"], cardinal_ph["opt"], cardinal_ph["max"]
    if not low < opt < high:
        raise ValueError(f"cardinal pH values must satisfy min < opt < max, got {cardinal_ph}")
    value = float(ph)
    if not low < value < high:
        return 0.0
    numerator = (value - low) * (value - high)
    return float(numerator / (numerator - (value - opt) ** 2))


def context_growth_rate(context: CultureContext) -> float:
    """Unstressed growth rate in these conditions, 1/h.

    Three multiplicative terms on the aerated carbon-source rate: growth phase, temperature and
    nothing else. Multiplying rather than fitting an interaction is Rosso 1995's independence
    hypothesis (PMID 16534932), stated there for temperature and pH.

    Oxygen costs a respiring culture more than a fermenting one, which is why the carbon source
    and the aeration cannot be treated as independent knobs and are resolved together in
    :func:`_aerated_carbon_rate`.

    **pH is deliberately absent.** It is carried on the context and recorded by the sweep, and
    :func:`ph_growth_factor` refuses to turn it into a rate without measured cardinal values.
    """
    rate = _aerated_carbon_rate(context) * _PHASE[context.growth_phase]["mu"]
    rate *= temperature_factor(context.temperature_c)
    return max(rate, 0.0)


def context_physiology(context: CultureContext):
    """Measured carbon and oxygen fluxes for an unstressed culture in these conditions.

    A thin bridge from a context to ``culture.chemostat_physiology``, and it inherits that
    function's refusal: a context whose growth rate lands outside van Hoek's measured range
    raises rather than returning an extrapolated yield. Most contexts in a wide sweep will,
    which is the honest answer -- the measured curve is aerobic, glucose-limited and 30
    degrees, and a stationary-phase ethanol culture at 37 is none of those things.

    Returns:
        A ``culture.ChemostatPhysiology``.
    """
    return chemostat_physiology(context_growth_rate(context))


def baseline_activity(context: CultureContext) -> dict[str, float]:
    """Module activity these conditions produce with no stressor applied.

    Everything a real plate carries before dosing: a derepressed carbon regulon off
    glucose, the general stress response of a culture past exponential, the hypoxic arm
    under poor aeration, and whatever heat response a warm incubator provokes.
    """
    carbon = _CARBON[context.carbon_source]
    baseline = {name: 0.0 for name in MODULES}
    baseline["carbon"] = carbon["carbon"]
    baseline["retrograde"] = carbon["retrograde"]
    baseline["ESR"] = _PHASE[context.growth_phase]["ESR"]
    shortfall = max(0.0, 1.0 - context.oxygen / _HYPOXIA_THRESHOLD)
    baseline["hypoxia"] = shortfall
    warmth = max(0.0, context.temperature_c - _REFERENCE_TEMPERATURE)
    baseline["heat"] = min(1.0, _HEAT_PER_DEGREE * warmth)
    return {name: min(1.0, max(0.0, value)) for name, value in baseline.items()}
