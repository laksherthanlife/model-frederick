"""The stress channel that is big enough to matter: carbon diverted into stress effectors.

WHY THIS EXISTS, and why it is not the channel this project kept testing. Six routes from a
stress state into the GEM were measured and all failed -- NGAM maintenance, GECKO protein
pool, precursor ceilings, NADPH competition, Crabtree partitioning, cassette burden. The
largest of them moved titre 0.17-1.73%, against an entry-flux law whose own leave-one-strain
-out error is 22.2%. Every one was below the noise floor of the thing it had to move.

They were all the same shape: a SCALAR TAX. NGAM adds ATP demand; the protein pool shrinks a
budget. This module is a different shape, and it is the one the architecture actually implies.

WHAT THE GEM DOES AND DOES NOT CONTAIN, checked rather than assumed. Yeast9 carries the
METABOLIC EFFECTORS of the stress response in full -- 5 trehalose reactions, 265 glycerol,
26 glutathione, 11 thioredoxin, 4 catalase/peroxide. It contains NONE of the signalling that
decides when to use them: MSN2, MSN4, HOG1, HAC1, YAP1, SKN7 and PBS2 are all absent, because
a genome-scale METABOLIC model has no regulators.

**So the GEM has the machinery and no reason to use it.** FBA maximises growth, and trehalose,
glycerol and glutathione all cost carbon and NADPH, so a growth-maximising solve will never
make them. That missing decision is exactly what a learned latent stress state supplies: not
a tax on the objective, but a REQUIREMENT that flux flows into named effector reactions.

THE CHANNEL IS LARGE PER UNIT FLUX AND THE FLUXES ARE SMALL, and an earlier version of this
docstring got that badly wrong. Measured on Yeast9 at glucose -10.0 mmol/gDCW/h with oxygen
free and growth HELD at 0.18 /h -- pure diversion, not a growth effect -- against a baseline
GGPP supply of 0.7257 mmol/gDCW/h::

    effector        0.5 mmol/gDCW/h    2.0 mmol/gDCW/h
    trehalose            -12.8%             -51.7%
    glutathione          -14.8%             -60.6%
    glycerol              -3.5%             -14.2%

**THOSE FLUXES ARE NOT PHYSIOLOGICAL AND THE FIRST VERSION OF THIS MODULE QUOTED THEM AS IF
THEY WERE.** A literature sweep of all 24 stress modules put real numbers on them, and the
gap is one to three orders of magnitude::

    module      effector       PHYSIOLOGICAL flux   measured cost
    oxidative   glutathione    0.0047               -0.11%
    peroxide    catalase       0.48                  0.00%   <- metabolically FREE
    peroxide    Tsa1 cycle     0.48                 -0.39%
    heat        trehalose      0.0526 steady        -1.32%
    heat        trehalose      0.409 acute shock   -10.41%
    ESR         glycogen       0.222                -2.83%
    osmotic     glycerol       1.49                -10.54%
    redox       GSSG cycle     9.6                  -9.69%
    dna_damage  dNTP           0.0035               -0.08%
    metals      all            <0.002              ~-0.00%

Glutathione's own pool is 8.8 umol/gDW, so holding it against dilution at mu = 0.18 costs
0.00158 mmol/gDCW/h -- the 0.5 rung above is **316x** that, and clearing the 22.2% floor
would need 0.78, which is 90-490x the physiological synthesis flux.

**SO AT REAL FLUXES THIS CHANNEL MOSTLY DOES NOT CLEAR THE FLOOR EITHER.** The largest honest
diversions are glycerol under osmotic stress (-10.5%), an acute heat-shock trehalose transient
(-10.4%) and glutathione redox cycling (-9.7%) -- an order of magnitude above the six scalar
taxes that failed before, and still half the 22.2% they would have to beat. Reported here
rather than buried, because "ten to forty times NGAM" was true and "large enough to matter"
was not.

THEN THE REQUIRED ARGUMENT WAS SUPPLIED FROM A MEASUREMENT, AND THE CHANNEL COLLAPSED
FURTHER. This module's contract says the diversion flux must come from a measurement and be
cited. A search of ~53 papers across five product families found exactly one aerobic dataset
that measures a stress effector in absolute units on one strain across graded conditions:
Hakkaart 2020 (PMID 31654410, PMC7028085), CEN.PK113-7D in an aerobic glucose-limited
chemostat at D = 0.025 /h, 2x2 factorial pH 5.0/3.0 x CO2 0.04%/50%, glycogen and trehalose
in mg/gDCW with standard deviations. Converting content to the flux that holds it against
dilution (q = mu * content, the same arithmetic this module already uses) and pricing it
through :func:`combined_cost`::

    pH 5 reference     glycogen 0.00545, trehalose 0.00142 mmol/gDCW/h    -0.0925%
    pH 3               glycogen 0.00716, trehalose 0.00092               -0.1007%
    pH 3 + 50% CO2     glycogen 0.00469, trehalose 0.00071               -0.0683%

**The stress-induced INCREMENT is 0.008 percentage points** -- the whole channel, from the
largest verified low-pH effector response in S. cerevisiae. Against the 22.2% floor that is
2,700x too small, not 1.5x.

AND ONE NUMBER ABOVE IS RETRACTED. The glycerol row -- 1.49 mmol/gDCW/h giving -10.5% -- is
not a stress increment. It is close to the UNSTRESSED ANAEROBIC baseline (0.79-1.11 across
four independent cultures), and Jouhten 2008 (PMID 18613954) measures q_glycerol =
0.00 +/- 0.00 at every oxygenated level from 0.5% to 20.9% O2. **Aerobically there is no
glycerol flux to divert at all**, so the largest surviving diversion in the table above does
not exist under the conditions this repository predicts in.

AND THE SIGN IS THE OTHER WAY ROUND. Every published dataset that dosed a stressor on one
strain and measured an absolute biomass-normalised content found it RISING, not falling:
astaxanthin +83% with biomass falling 8.87 -> 8.00 g/L (PMID 16896607); beta-carotene
35.77 -> 117.62 ug/gdw by HPLC with biomass flat under 0.75 M NaCl (PMID 39399738); NaCl
+56%, C:N x1.51-x3.98, 20 C vs 30 C x59.

**NOT ONE OF THOSE THREE IS S. CEREVISIAE**, and that limit is load-bearing rather than
pedantic. They are the carotenogenic yeasts -- X. dendrorhous (PMID 16896607),
R. odoratus (PMID 39399738) and R. mucilaginosa (PMID 31487889) -- which make these
carotenoids NATIVELY, from native promoters under native regulation. This repository
predicts a HETEROLOGOUS cassette in S. cerevisiae, where the pathway's genes are foreign
and sit under a chosen promoter. So these three establish the sign and the mechanism for
native carotenogenesis; they do NOT establish either for a cassette, and the direction a
cassette moves is set by whatever its promoter does. The whole argument of this project is
that the HOST is constant -- which is exactly why evidence from three other hosts cannot
be quietly spent on this one.

Where anyone measured the mechanism it is
TRANSCRIPTIONAL and not metabolic: in R. mucilaginosa under H2O2, HMG1 and ERG12 transcripts
did not move -- the authors write that peroxide "has no effect on the carbon flux through the
mevalonate pathway" -- while CAR0, CAR1 and CAR2 all fell (PMID 31487889). Precursor supply
was unchanged; the pathway's own genes changed.

**So this module models the wrong mechanism.** The arithmetic in it is right -- the shared
budget, the walls, the additivity are all measured and reproduce -- but stress reaches a
product by changing how much PATHWAY ENZYME there is, which a metabolic model cannot
represent, and not by stealing the pathway's carbon. That also explains the size mismatch
without appealing to anything else: diversion is small because the carbon at stake is small,
while up- or down-regulating the pathway moves flux directly and reaches 60-80%.

WHAT DOES CLEAR THE FLOOR, and it is not a diversion at all: **oxygen**. Restricting the oxygen
exchange moves precursor supply -67.2%, because supply is almost exactly proportional to
oxygen uptake (supply = 0.0616*qO2 + 0.0101, R^2 = 0.999989) for an ATP- and NADPH-hungry
mevalonate pathway. That is a REROUTING driven by an input the GEM already accepts directly,
so it needs no learned latent state -- which is why `predict.py` can already answer for
oxygen and why the hypoxia module is the one place environment reaches product for free.

WHAT THIS MODULE DOES NOT DO, and the limit is the same one that governs everything else
here. **It does not know how much effector flux a given environment causes.** That number is
the learned latent state, and this repository has not measured it -- `bridge/latent_bridge.py`
holds a reporter whose four constructs G4 returned INCONCLUSIVE for. So the diversion flux is
a REQUIRED ARGUMENT with no default: supply it from a measurement and cite the measurement,
or do not call this. What the module contributes is the transfer function from that flux to
the product, which is a property of the HOST and therefore the part that generalises.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "BUDGET_ZERO_GROWTH_RATE",
    "INFEASIBILITY_WALL_AT_MU_018",
    "STRESS_EFFECTORS",
    "combined_cost",
    "DiversionResult",
    "EffectorUnmeasured",
    "diversion_cost",
]

#: The metabolic outputs a stress response commits carbon to, by Yeast9 metabolite id.
#: Cytosolic, because that is where the protective pool accumulates. Deliberately a small
#: named set rather than a search: a name match would pick up 265 glycerol reactions, most of
#: which are lipid metabolism and not the osmotic response.
STRESS_EFFECTORS = {
    "trehalose": "s_1520",     # heat and desiccation protectant, Msn2/4 output
    "glycerol": "s_0765",      # osmotic response, the HOG pathway's metabolic output
    "glutathione": "s_0750",   # oxidative response, Yap1 output
}


#: The flux at which each effector ALONE drives the model infeasible at mu = 0.18, on Yeast9
#: at glucose -10.0 with oxygen free. Measured, not assumed.
#:
#: THESE ARE NOT METABOLIC CEILINGS. Maximum feasible growth is exactly LINEAR in total
#: diversion -- slope constant to five decimals across the whole range -- so a "wall" is just
#: where the held growth rate runs out of budget: ``D_max(mu) = D0 * (1 - mu/0.88769)``.
#: Trehalose affords 3.714 at mu = 0.18, 4.658 at mu = 0, and 2.94 at the stressed mu = 0.327
#: the panel itself predicts.
INFEASIBILITY_WALL_AT_MU_018 = {
    "er_membrane_pc": 0.647, "squalene": 0.963, "glutathione": 3.082, "dntp": 3.618,
    "trehalose": 3.714, "chitin": 5.152, "cysteine": 5.654, "glycogen": 7.309,
    "beta_1_3_glucan": 7.309, "myo_inositol": 7.550, "glycerol": 12.852,
    "nadph_cycle": 91.762, "atp": 152.450,
}

#: Growth rate at which the diversion budget vanishes, from the linear fit above.
BUDGET_ZERO_GROWTH_RATE = 0.88769


def combined_cost(effector_fluxes: dict, growth_rate_per_h: float = 0.18) -> float:
    """Total cost of several simultaneous stress effectors, as a fraction of precursor supply.

    ``cost = -sum(f_i / W_i)``, with the walls scaled to the growth rate.

    THE COSTS ADD, and that was measured rather than assumed -- it is the one result that
    makes a multi-module stress state computable at all. For each pair of effectors the joint
    infeasibility point along the ray ``(t*W_a, t*W_b)`` was found: a SHARED budget predicts
    ``t* = 0.5``, INDEPENDENT budgets predict 1.0. **All 45 pairs came back 0.5000 to 0.5137,
    mean 0.50255**, and ten effectors at once gave ``sum(f_i/W_i) = 1.011`` at the wall
    against 1.0 shared and 10.0 independent. It is one budget, to 1.1%.

    Sharing a biosynthetic intermediate does NOT make two effectors cheaper together, which
    was the obvious hypothesis and is false: trehalose and glycogen both come off UDP-glucose
    and still give exactly t* = 0.5000.

    The residual is a mild SUPER-additivity -- 44 of 45 pairs cost slightly more jointly than
    the sum of their solos, by 0.1-2.9% at realistic loads and at most 9% near the wall. The
    mechanism is that each solo cost is slightly CONVEX in flux (cheap at the first mmol,
    linear near the wall), so a second effector raises the marginal price of the first. ATP is
    the most convex, at 0.59 of the linear law at 2% of its wall, because the aerobic optimum
    absorbs the first ~15 mmol ATP/gDCW/h nearly free.

    Args:
        effector_fluxes: ``{name: mmol/gDCW/h}``, names from
            :data:`INFEASIBILITY_WALL_AT_MU_018`.

    Returns:
        A negative fraction. -0.10 means a 10% loss of precursor supply.

    Raises:
        KeyError: on an effector with no measured wall.
        ValueError: on a non-positive growth rate, or a bundle past the budget.
    """
    if growth_rate_per_h <= 0:
        raise ValueError(f"growth rate must be positive, got {growth_rate_per_h}")
    scale = 1.0 - growth_rate_per_h / BUDGET_ZERO_GROWTH_RATE
    if scale <= 0:
        raise ValueError(
            f"growth rate {growth_rate_per_h} leaves no diversion budget at all; the budget "
            f"vanishes at {BUDGET_ZERO_GROWTH_RATE} /h")
    total = 0.0
    for name, flux in effector_fluxes.items():
        if name not in INFEASIBILITY_WALL_AT_MU_018:
            raise KeyError(
                f"no measured wall for {name!r}; have "
                f"{sorted(INFEASIBILITY_WALL_AT_MU_018)}")
        if flux < 0:
            raise ValueError(f"flux for {name!r} must be non-negative, got {flux}")
        wall_at_018 = INFEASIBILITY_WALL_AT_MU_018[name]
        wall = wall_at_018 * scale / (1.0 - 0.18 / BUDGET_ZERO_GROWTH_RATE)
        total += flux / wall
    if total >= 1.0:
        raise ValueError(
            f"the bundle sums to {total:.3f} of the diversion budget at mu = "
            f"{growth_rate_per_h:g} /h, so there is no feasible solution. That boundary is "
            f"the HELD GROWTH RATE running out of carbon, not a metabolic ceiling -- lower "
            f"the growth rate and the same bundle becomes affordable")
    return -total


class EffectorUnmeasured(ValueError):
    """No diversion flux was supplied, and there is no defensible default.

    Its own type because the correct response is specific: measure the effector pool or its
    accumulation rate in the condition of interest, or state plainly that the environment ->
    stress-state half is unmeasured. Inventing a flux here would make the product number a
    property of this module rather than of the culture, which is the failure every refusal in
    this package exists to prevent.
    """


@dataclass(frozen=True)
class DiversionResult:
    """What a given effector flux costs the pathway, at a held growth rate."""

    effector: str
    flux_mmol_per_gdcw_h: float
    growth_rate_per_h: float
    precursor_supply_mmol_per_gdcw_h: float
    baseline_supply_mmol_per_gdcw_h: float
    conditions: str

    @property
    def relative_change(self) -> float:
        """Negative is diversion away from the product."""
        return ((self.precursor_supply_mmol_per_gdcw_h
                 - self.baseline_supply_mmol_per_gdcw_h)
                / self.baseline_supply_mmol_per_gdcw_h)

    @property
    def clears_the_calibration_floor(self) -> bool:
        """Against the entry-flux law's own 22.2% leave-one-strain-out error. Below this a
        change cannot be checked against anything, whatever its sign."""
        return abs(self.relative_change) > 0.222

    def summary(self) -> str:
        return (f"{self.effector} at {self.flux_mmol_per_gdcw_h:g} mmol/gDCW/h, mu held at "
                f"{self.growth_rate_per_h:g} /h: precursor supply "
                f"{self.precursor_supply_mmol_per_gdcw_h:.4f} against "
                f"{self.baseline_supply_mmol_per_gdcw_h:.4f} "
                f"({100 * self.relative_change:+.1f}%) [{self.conditions}]")


def diversion_cost(model, spec, effector: str, flux_mmol_per_gdcw_h: float,
                   growth_rate_per_h: float, *,
                   glucose_lower_bound: float = -10.0,
                   oxygen_lower_bound: float = -1000.0,
                   glucose_exchange: str = "r_1714",
                   oxygen_exchange: str = "r_1992",
                   biomass_reaction: str = "r_2111") -> DiversionResult:
    """What forcing ``flux`` of ``effector`` costs this pathway's precursor supply.

    THE GROWTH RATE IS HELD, which is the whole point. Every previous stress route in this
    repository reached the product through mu and was therefore indistinguishable from the
    growth-rate channel the stress panel already measures directly. Holding mu makes the
    diversion the only thing that can move the answer.

    Args:
        model: A cobra model. Not mutated -- all changes are inside ``with model:``.
        spec: A :class:`~ystwin.pathway.spec.PathwaySpec`; only ``precursor_metabolite`` and
            ``precursor_stoichiometry`` are read, so this works for a pathway with no
            calibration at all.
        flux_mmol_per_gdcw_h: How much effector the cell is making. **Required, and there is
            no default**: it is the learned latent state, and this repository has not measured
            it. Supply it from a measurement and cite that measurement.

    Raises:
        EffectorUnmeasured: on a non-positive flux.
        KeyError: on an effector outside :data:`STRESS_EFFECTORS`.
        ValueError: on a non-positive growth rate.
    """
    import cobra

    if effector not in STRESS_EFFECTORS:
        raise KeyError(
            f"no effector {effector!r}; have {sorted(STRESS_EFFECTORS)}. These are the "
            f"metabolic outputs of the stress response that Yeast9 actually carries -- the "
            f"signalling that drives them (MSN2/4, HOG1, YAP1) is not in a metabolic model")
    if flux_mmol_per_gdcw_h <= 0:
        raise EffectorUnmeasured(
            f"diversion flux must be positive, got {flux_mmol_per_gdcw_h}. There is no "
            f"default because this number IS the learned latent stress state, and it is "
            f"unmeasured here: bridge/latent_bridge.py's four reporter constructs all "
            f"returned INCONCLUSIVE at G4. Measure the effector pool or its accumulation "
            f"rate in the condition of interest, or do not call this")
    if growth_rate_per_h <= 0:
        raise ValueError(f"growth rate must be positive, got {growth_rate_per_h}")

    def supply(with_diversion: bool) -> float:
        with model:
            if with_diversion:
                demand = cobra.Reaction("YSTWIN_STRESS_DIVERSION")
                demand.lower_bound = flux_mmol_per_gdcw_h
                demand.upper_bound = 1000.0
                demand.add_metabolites(
                    {model.metabolites.get_by_id(STRESS_EFFECTORS[effector]): -1.0})
                model.add_reactions([demand])
            drain = cobra.Reaction("YSTWIN_PRECURSOR_DRAIN")
            drain.lower_bound, drain.upper_bound = 0.0, 1000.0
            drain.add_metabolites(
                {model.metabolites.get_by_id(spec.precursor_metabolite):
                 -float(spec.precursor_stoichiometry)})
            model.add_reactions([drain])
            model.reactions.get_by_id(glucose_exchange).lower_bound = glucose_lower_bound
            model.reactions.get_by_id(oxygen_exchange).lower_bound = oxygen_lower_bound
            model.reactions.get_by_id(biomass_reaction).bounds = (
                growth_rate_per_h, growth_rate_per_h)
            model.objective = drain.id
            solution = model.optimize()
            if solution.status != "optimal" or solution.objective_value is None:
                return float("nan")
            return float(solution.objective_value)

    return DiversionResult(
        effector=effector,
        flux_mmol_per_gdcw_h=float(flux_mmol_per_gdcw_h),
        growth_rate_per_h=float(growth_rate_per_h),
        precursor_supply_mmol_per_gdcw_h=supply(True),
        baseline_supply_mmol_per_gdcw_h=supply(False),
        conditions=(f"glucose {glucose_lower_bound:g}, oxygen {oxygen_lower_bound:g} "
                    f"mmol/gDCW/h, mu held"))
