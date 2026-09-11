"""The GEM as a source of the environment signal -- asked for a YIELD, never for a ceiling.

WHY THIS EXISTS, and why an earlier attempt concluded the opposite. Two workflows measured
whether stress or environment could reach titre through the GEM, and both asked it for a
CEILING: "how much product could the host supply?" That question fails, and the failure is
real -- the chain's predicted flux sits 12.8x to 195x below the ceiling so it never binds, and
scaling the ceiling by each state's own measured uptake gives a captured fraction spanning
216x, because at high growth rate the measured biomass yields collapse and the inferred uptake
explodes with them. The conclusion drawn was that the GEM cannot carry the signal.

That conclusion was wrong, and the reason is the question. **Content is not set by a ceiling,
it is set by a yield**, and yield is the one thing flux balance analysis is genuinely good at.
Asked for the maximum product obtainable per C-mol fed, at MATCHED carbon and MATCHED growth
rate, Yeast9 ranks the carbon sources correctly on the one axis where a measurement exists.
Under the SHIPPED formulation -- bare precursor sink, 11 mmol C/gDCW/h, mu = 0.05 /h, which is
the rate Kocharin measured the contrast at::

    max PHB per C-mol:   glucose 0.736    ethanol 1.000    -> ethanol <!-- audit:value table=outputs/gem_environment_predictions.csv column=ethanol_over_glucose row="product=phb; growth_rate_per_h=0.05" -->1.3584x better

which is the direction Kocharin 2013 measured (ethanol 3.82x glucose at mu = 0.05 /h). The
earlier "glucose is best in 5 of 6 regimes" came from maximising a precursor TAP with biomass
pinned, where glucose wins simply because it carries more ATP per C-mol.

*This read "glucose 0.716 ethanol 1.000 -> 1.397x" and named neither its formulation nor its
growth rate until 2026-09-04. It reproduces only under the FULL-PATHWAY formulation -- the one
this module does not ship -- and at mu = 0, where "MATCHED growth rate" is vacuous. A
paragraph whose own thesis is "a result that does not name its formulation is not a
measurement" may not violate it in its first result, so the sentence now carries both, and
the number is the shipped one.*

WHAT IT BUYS, AND WHAT IT DOES NOT. **It gets the direction and not the magnitude, and an
earlier version of this docstring overstated it.** Measured leave-one-CARBON-SOURCE-out on the
eleven Kocharin states -- predicting a feed the law has never seen::

    mu only       (the chain today)                       <!-- audit:value table=outputs/gem_environment_scores.csv column=leave_one_carbon_source_out_fold_error row="law=mu only" -->1.7237x
    log(GEM max)  (stoichiometric, 0 fitted)              <!-- audit:value table=outputs/gem_environment_scores.csv column=leave_one_carbon_source_out_fold_error row="law=log GEM max" -->1.4683x
    z*log(ratio)  (the pure carbon contrast, this module) <!-- audit:value table=outputs/gem_environment_scores.csv column=leave_one_carbon_source_out_fold_error row="law=GEM carbon contrast" -->2.1942x
    z             (empirical, 1 fitted)                   <!-- audit:value table=outputs/gem_environment_scores.csv column=leave_one_carbon_source_out_fold_error row="law=z empirical" -->1.3439x

The middle row looks like a win over mu-only and is not one:
``corr(log(GEM max), log mu) = <!-- audit:value table=outputs/gem_environment_scores.csv column=corr_with_log_mu row="law=log GEM max" -->-0.8190``. The absolute
maximum is 82% collinear with growth rate, so it beats mu by BEING mu, with a little carbon
attached. Strip that out and score the pure carbon contrast -- which is what an environment
descriptor has to be, zero on the reference feed -- and it correlates
**<!-- audit:value table=outputs/gem_environment_scores.csv column=corr_with_log_content row="law=GEM carbon contrast" -->0.3656** with content against the
empirical z's **<!-- audit:value table=outputs/gem_environment_scores.csv column=corr_with_log_content row="law=z empirical" -->0.5801**, and scores WORSE than
mu alone.

So the honest summary is narrower than "the GEM carries the environment signal":

  * **The SIGN is predicted for two of four products AT mu = 0.05 /h, needs no data, and is
    STABLE THERE.** Swept across the plausible carbon supply (4-25 mmol C/gDCW/h, spanning
    Kocharin's measured 3.99-4.52 up to the value this module originally shipped), ethanol
    beats glucose for beta-carotene over
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_low row="product=beta_carotene; growth_rate_per_h=0.05" -->1.2362 to <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_high row="product=beta_carotene; growth_rate_per_h=0.05" -->1.9505x
    and for PHB over
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_low row="product=phb; growth_rate_per_h=0.05" -->1.3138 to <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_high row="product=phb; growth_rate_per_h=0.05" -->2.3472x,
    never crossing 1.0. **Gadusol and glycogen are NOT predictions**: they span
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_low row="product=glycogen; growth_rate_per_h=0.05" -->0.8504 to <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_high row="product=glycogen; growth_rate_per_h=0.05" -->1.6323x,
    crossing 1.0 inside that range, so their sign is a property of the supply parameter and
    not of the chemistry. They also agree with each other to 0.15%, so they were never two
    predictions in any case. :func:`sign_is_stable` is what separates the two groups and
    should gate any use. That is a real, falsifiable, product-generic claim and it is what
    this module is for.

    **"AT mu = 0.05" IS A RETRACTION MADE ON 2026-09-04, and it is the SECOND instance of the
    finding commit cf34701 recorded.** That commit found a robustness claim quantified over
    nothing at all -- a supply parameter shipped at one unswept value. :func:`sign_is_stable`
    was the repair, and it reproduced the disease one level down: it swept the SUPPLY and
    left the GROWTH RATE pinned at whatever the caller passed, returned a bare
    ``(stable, low, high)`` with no record of what it had evaluated, and ``continue``d
    silently past every infeasible point. So the mu = 0.10 and mu = 0.15 rows shipped as
    predictions on
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=points_feasible row="product=phb; growth_rate_per_h=0.15" -->4 of
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=points_requested row="product=phb; growth_rate_per_h=0.15" -->6 swept points,
    every one of them near a boundary where the ratio diverges: PHB reads
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=ethanol_over_glucose row="product=phb; growth_rate_per_h=0.15" -->3.0139x
    at mu = 0.15 against
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=ethanol_over_glucose row="product=phb; growth_rate_per_h=0.05" -->1.3584x
    at mu = 0.05, and its glucose maximum has fallen to
    <!-- audit:value table=outputs/gem_environment_predictions.csv column=glucose_max_mmol_per_gdcw_h row="product=phb; growth_rate_per_h=0.15" -->0.146891
    from <!-- audit:value table=outputs/gem_environment_predictions.csv column=glucose_max_mmol_per_gdcw_h row="product=phb; growth_rate_per_h=0.05" -->1.537439.
    :class:`SignStability` now returns the sweep's own provenance and
    :data:`MIN_GROWTH_HEADROOM` gates on it, so those rows carry ``is_a_prediction=False``.
    The lesson is the one cf34701 wrote and this module then had to learn twice: sweeping the
    forgotten parameter fixes today's rows, and only RETURNING WHAT THE SWEEP COVERED makes
    the next forgotten parameter visible.
  * **The MAGNITUDE is not.** Measured is 3.82x for PHB at mu = 0.05; the GEM says 1.32x. All
    of the predictive value in the empirical law lives in the size of the effect, and the GEM
    does not supply it.

Use this to form a PRIOR on which feed to try and on which way a coefficient should point --
not as a substitute for measuring that coefficient.

THE FORMULATION IS PINNED, AND THAT IS NOT A DETAIL. The GEM's answer depends on how it is
asked. Same model, same data, different question. An earlier analysis saw this as "the
coefficient's sign reverses when the matched uptake changes" and read it as a refutation; it
is really a statement that an unpinned formulation has no defined answer. So
:class:`Formulation` records every choice, every result carries the formulation that produced
it, and comparisons between differently-formulated numbers are refused rather than silently
made.

MEASURED, at last, on 2026-09-04. At 20 mmol C/gDCW/h the bare precursor sink gives PHB an
ethanol/glucose ratio of
<!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_precursor_sink row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.05" -->1.3212 to <!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_precursor_sink row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.2" -->1.606
and the full PHA reaction carrying its NADPH cost gives
<!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_full_pathway row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.05" -->1.4001 to <!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_full_pathway row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.2" -->1.5818 --
a few percent apart, <!-- audit:value table=outputs/gem_environment_formulation.csv column=full_over_sink row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.05" -->1.0597x at the widest. The whole table
is ``outputs/gem_environment_formulation.csv``.

**THIS PARAGRAPH CLAIMED "1.09-1.25 AGAINST 1.40-2.00" UNTIL THEN, AND NEITHER HALF WAS A
NUMBER THIS REPOSITORY COULD PRODUCE.** The sink half was the value from before commit
385b6b8, which corrected the precursor drain to return its group carrier; that commit updated
two occurrences of "1.09-1.25" in this file and left two more, one of them in the
`precursor_sink_only` field comment below. The full-pathway half had never been computed at
all: :attr:`Formulation.precursor_sink_only` was declared, printed in :meth:`Formulation.label`,
copied through :func:`sign_is_stable` and used in :meth:`GemEnvironmentResult.comparable_with`,
and NEVER READ by :func:`_max_precursor` -- so a result labelled "full pathway" was
bit-identical to the sink, and `comparable_with` refused to compare two results that were the
same computation. It is now read, via the spec's declared net cofactor balance, and a spec
that has not declared one raises :class:`FullPathwayUndeclared` rather than falling back.
So the sensitivity this class exists to pin is real but SMALL, and the sentence that made it
the argument for the class was quoting a gap that did not exist.

WHAT IT STILL CANNOT DO, and none of these is hidden by the code:

  * **Magnitude.** It gets about 30% of the log effect -- 1.4-2.0x against a measured 3.82x.
  * **The reversal.** Kocharin's ethanol advantage SHRINKS with growth rate and flips sign at
    mu = 0.1718 /h; the GEM's grows. That is expected and it is the boundary between the two
    halves of the vision: the GEM says what is POSSIBLE, and a cell does not maximise a
    storage polymer. How much of the possible is taken -- the captured fraction -- is the
    part a learned regulatory state has to supply, and nothing here supplies it.
  * **Significance.** Three carbon sources put the smallest attainable permutation p at
    1/3! = 0.1667. Nothing computed here can be earned statistically on this dataset.

THE FALSIFIABLE PART. Because the descriptor is stoichiometric it can be computed for a
product with no data at all, which makes it a committed prediction rather than a fit. At
mu = 0.05 /h it says ethanol beats glucose for beta-carotene -- GGPP ratio
<!-- audit:value table=outputs/gem_environment_predictions.csv column=ethanol_over_glucose row="product=beta_carotene; growth_rate_per_h=0.05" -->1.2384 at the
default supply, <!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_low row="product=beta_carotene; growth_rate_per_h=0.05" -->1.2362-<!-- audit:value table=outputs/gem_environment_predictions.csv column=ratio_high row="product=beta_carotene; growth_rate_per_h=0.05" -->1.9505
across the swept one -- the same sign as PHB, so the planned matched-mu bench test has a real
prior and one that can be wrong. `scripts/gem_environment_predictions.py` commits the whole
table before the experiment runs.

*"GGPP ratio 1.24-1.33" here was the pre-385b6b8 value and outlived the commit that
invalidated it, in the same way and for the same reason as the formulation sensitivity above:
nothing pointed it at a cell.*
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

__all__ = [
    "DEFAULT_FORMULATION",
    "MIN_GROWTH_HEADROOM",
    "sign_is_stable",
    "PrecursorCarrierUnknown",
    "Formulation",
    "GemEnvironmentResult",
    "FullPathwayUndeclared",
    "InfeasibleState",
    "SignStability",
    "environment_descriptor",
    "max_growth_rate",
    "yield_ratio",
]

#: How far below the model's own maximum growth rate a swept point must sit before its yield
#: ratio counts toward a prediction. Dimensionless: ``1 - mu / mu_max`` on the tighter of the
#: two feeds.
#:
#: WHY A HEADROOM GATE EXISTS AT ALL. The ratio diverges as the LP approaches its feasibility
#: boundary, because both maxima collapse toward zero and the ratio becomes a ratio of two
#: small numbers. Measured on PHB at the default 11 mmol C/gDCW/h, one LP per row::
#:
#:     headroom  0.94   0.69   0.44   0.25   0.13   0.063  0.032
#:     ratio     1.313  1.358  1.519  1.707  2.145  3.014  4.730
#:
#: There is no knee -- the divergence is smooth -- so this number is a STATED CONVENTION and
#: not a discovered threshold, and it is anchored rather than chosen: 0.10 is the headroom of
#: the tightest state Kocharin 2013 ACTUALLY RAN. Their three mu = 0.05 chemostats took
#: 3.99, 4.52 and 4.23 mmol C/gDCW/h, which this model puts at headroom 0.107, 0.211 and
#: 0.161. A gate tighter than the source experiment's own operating point would refuse states
#: a real vessel has held; a looser one admits ratios the boundary is setting.
MIN_GROWTH_HEADROOM = 0.10

#: Carbon atoms per molecule, for converting a matched C-mol supply into exchange bounds.
_CARBON_ATOMS = {"glucose": 6.0, "ethanol": 2.0}

#: The GROUP CARRIER each precursor hands back when its carbon is consumed, and how much of
#: it per precursor molecule. ``None`` means the precursor carries no group.
#:
#: THIS EXISTS BECAUSE OMITTING IT IS A 15.5x ERROR AND A SILENT ONE. A pathway pulling
#: acetyl-CoA consumes the ACETYL and returns the CoA -- thiolase condenses two acetyl-CoA
#: into acetoacetyl-CoA plus free CoA. A demand written as ``{acetyl-CoA: -2}`` alone forces
#: the model to synthesise coenzyme A de novo for every drained molecule, so what it reports
#: is a COA-BIOSYNTHESIS CEILING and not a carbon yield. Measured on Yeast9 at 20 mmol
#: C/gDCW/h, glucose, mu = 0.10: 0.1749 without the carrier against 2.7130 with it.
#: UDP-D-glucose has the same shape (glycogen synthase returns UDP), while GGPP and
#: sedoheptulose-7-phosphate carry no group and are unaffected -- so the defect deflated one
#: precursor family fifteenfold and left the other untouched, destroying exactly the
#: between-product contrast this module exists to compute.
#:
#: Keyed by precursor id and deliberately EXHAUSTIVE rather than defaulted: a precursor
#: absent from this mapping raises, so a new pathway must state its carrier chemistry instead
#: of silently receiving a number that measures the wrong thing.
_PRECURSOR_CARRIER = {
    "s_0373": ("s_0529", 1.0),   # acetyl-CoA  -> coenzyme A, one per acetyl consumed
    "s_1543": ("s_1538", 1.0),   # UDP-D-glucose -> UDP, one per glucosyl transferred
    "s_0189": None,              # geranylgeranyl diphosphate: prenyl carbon, no group returned
    "s_1427": None,              # sedoheptulose 7-phosphate: sugar phosphate, no group
}


class PrecursorCarrierUnknown(ValueError):
    """This precursor's carrier chemistry has not been declared, so no yield is computed.

    Raised rather than defaulting to "no carrier", because that default is wrong by 15.5x on
    a CoA-thioester and right on a prenyl phosphate, and nothing in a metabolite id says
    which. A wrong number here is worse than a refusal: it is a silent 15x deflation of one
    precursor family against another.
    """


class InfeasibleState(ValueError):
    """The model cannot hold this growth rate on this carbon supply.

    Raised rather than returning zero, because zero is a legitimate descriptor value (a
    product the host cannot make) and infeasibility is not the same claim. Conflating them
    would let "this vessel cannot run" be reported as "this feed is bad for this product".
    """


@dataclass(frozen=True)
class Formulation:
    """Every choice that changes the GEM's answer, recorded so results are comparable.

    This class exists because the sensitivity is large enough to reverse conclusions. A
    result that does not name its formulation is not a measurement, and two results from
    different formulations are not comparable -- :meth:`GemEnvironmentResult.comparable_with`
    enforces that rather than trusting a reader to notice.
    """

    #: THE MOST LOAD-BEARING NUMBER IN THIS CLASS, and it was chosen badly. It shipped at
    #: 20.0 mmol C/gDCW/h, which is 4.7x Kocharin's OWN measured specific uptake at the growth
    #: rate the predictions were headlined at (3.99, 4.52 and 4.23 at mu = 0.05; median 10.99
    #: over all eleven states). At the measured supply, glycogen and gadusol REVERSE -- 0.86
    #: becomes 1.63 -- and the module's advertised claim that they "flip at low growth rate"
    #: turns out to exist only for roughly 8 <= S <= 25. It was never swept in the committed
    #: work. Default lowered to the measured median, and :func:`sign_is_stable` now exists so
    #: that a sign which does not survive the plausible range is not reported as a prediction.
    matched_carbon_cmol_per_gdcw_h: float = 11.0e-3
    oxygen_lower_bound: float = -1000.0
    glucose_exchange: str = "r_1714"
    ethanol_exchange: str = "r_1761"
    oxygen_exchange: str = "r_1992"
    biomass_reaction: str = "r_2111"
    #: Pull on the PRECURSOR the spec declares, with no cofactor cost. Product-generic -- it
    #: needs only `PathwaySpec.precursor_metabolite`, which every spec already carries.
    #: ``False`` builds the drain from `PathwaySpec.full_pathway_stoichiometry` instead, and
    #: a spec that has not declared one RAISES :class:`FullPathwayUndeclared`: only `phb` has,
    #: because beta-carotene's CrtI cofactor chemistry is written down nowhere here and
    #: inventing it to fill a table is what this module is built against.
    #:
    #: MEASURED ON PHB AT 20 mmol C/gDCW/h: sink
    #: <!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_precursor_sink row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.05" -->1.3212-<!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_precursor_sink row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.2" -->1.606,
    #: full pathway <!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_full_pathway row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.05" -->1.4001-<!-- audit:value table=outputs/gem_environment_formulation.csv column=ratio_full_pathway row="product=phb; supply_mmol_c_per_gdcw_h=20.0; growth_rate_per_h=0.2" -->1.5818.
    #: The direction is as expected -- ethanol is the more reduced substrate, so carrying the
    #: NADPH cost favours it further at three of the four growth rates -- but the size is a
    #: few percent, not the "1.09-1.25 -> 1.40-2.00" this comment asserted until 2026-09-04.
    #: That number described an arm nothing could compute; see the module docstring.
    precursor_sink_only: bool = True

    def label(self) -> str:
        return (f"matched {1000 * self.matched_carbon_cmol_per_gdcw_h:.1f} mmol C/gDCW/h, "
                f"O2 lb {self.oxygen_lower_bound:g}, "
                f"{'precursor sink' if self.precursor_sink_only else 'full pathway'}")


DEFAULT_FORMULATION = Formulation()


@dataclass(frozen=True)
class GemEnvironmentResult:
    """One ethanol-versus-glucose yield ratio, with everything needed to check it."""

    product: str
    precursor_metabolite: str
    growth_rate_per_h: float
    glucose_max_mmol_per_gdcw_h: float
    ethanol_max_mmol_per_gdcw_h: float
    formulation: Formulation

    @property
    def ratio(self) -> float:
        """Ethanol over glucose. Above 1 means ethanol is the better feed for this product."""
        return self.ethanol_max_mmol_per_gdcw_h / self.glucose_max_mmol_per_gdcw_h

    @property
    def favours_ethanol(self) -> bool:
        return self.ratio > 1.0

    def comparable_with(self, other: "GemEnvironmentResult") -> bool:
        """Two results are comparable only under the same formulation. See the class docstring."""
        return self.formulation == other.formulation

    def summary(self) -> str:
        return (f"{self.product} at mu = {self.growth_rate_per_h:.4g} /h: ethanol/glucose = "
                f"{self.ratio:.3f}x [{self.formulation.label()}]")


class FullPathwayUndeclared(ValueError):
    """This spec has not written down its pathway's net cofactor balance, so no full-pathway
    yield is computed for it.

    The same refusal shape as :class:`PrecursorCarrierUnknown`, and for the same reason: a
    cofactor cost that is not declared is not a cofactor cost of zero. Falling back to the
    bare precursor sink under a ``precursor_sink_only=False`` formulation would return a
    number labelled "full pathway" that IS the sink -- which is exactly the defect this
    exception was written to close, on 2026-09-04. See :attr:`Formulation.precursor_sink_only`.
    """


def _drain_stoichiometry(model, spec, formulation) -> dict:
    """The reaction the LP maximises, per mole of PRODUCT, under this formulation.

    Two formulations, two drains, and they were ONE until 2026-09-04 -- see
    :attr:`Formulation.precursor_sink_only` for what that cost.

    Raises:
        PrecursorCarrierUnknown: under the bare sink, when the precursor's carrier chemistry
            is not declared.
        FullPathwayUndeclared: under the full pathway, when the spec has not written its net
            balance down.
    """
    if not formulation.precursor_sink_only:
        declared = getattr(spec, "full_pathway_stoichiometry", None)
        if not declared:
            raise FullPathwayUndeclared(
                f"pathway {spec.product!r} has no full_pathway_stoichiometry, so its "
                f"cofactor cost is unwritten and no full-pathway yield can be computed. "
                f"Declare the NET reaction per mole of product in the spec, or ask for the "
                f"bare precursor sink instead. This refuses rather than falling back: a "
                f"result labelled 'full pathway' that silently held the sink's number is "
                f"what this exception was added to prevent")
        return {model.metabolites.get_by_id(mid): float(coefficient)
                for mid, coefficient in declared.items()}

    if spec.precursor_metabolite not in _PRECURSOR_CARRIER:
        raise PrecursorCarrierUnknown(
            f"precursor {spec.precursor_metabolite!r} for {spec.product!r} has no "
            f"declared carrier chemistry; add it to _PRECURSOR_CARRIER. Draining a "
            f"group-carrying precursor without returning its carrier measures the "
            f"carrier's biosynthesis rather than a carbon yield -- a 15.5x error on "
            f"acetyl-CoA and none at all on GGPP, so it cannot be defaulted either way")
    stoichiometry = float(spec.precursor_stoichiometry)
    drain = {model.metabolites.get_by_id(spec.precursor_metabolite): -stoichiometry}
    carrier = _PRECURSOR_CARRIER[spec.precursor_metabolite]
    if carrier is not None:
        carrier_id, per_precursor = carrier
        drain[model.metabolites.get_by_id(carrier_id)] = stoichiometry * per_precursor
    return drain


def _max_precursor(model, spec, carbon_source, growth_rate_per_h, formulation):
    """Maximum pull on the declared precursor at a fixed growth rate and matched carbon."""
    import cobra

    supply = formulation.matched_carbon_cmol_per_gdcw_h * 1000.0
    with model:
        glucose = model.reactions.get_by_id(formulation.glucose_exchange)
        ethanol = model.reactions.get_by_id(formulation.ethanol_exchange)
        glucose.lower_bound = 0.0
        ethanol.lower_bound = 0.0
        if carbon_source == "glucose":
            glucose.lower_bound = -supply / _CARBON_ATOMS["glucose"]
        elif carbon_source == "ethanol":
            ethanol.lower_bound = -supply / _CARBON_ATOMS["ethanol"]
        else:
            raise ValueError(
                f"carbon source must be 'glucose' or 'ethanol', got {carbon_source!r}. "
                f"A blend is a mixture of the two bounds and belongs in a caller that states "
                f"its own C-mol split")
        model.reactions.get_by_id(
            formulation.oxygen_exchange).lower_bound = formulation.oxygen_lower_bound
        model.reactions.get_by_id(formulation.biomass_reaction).bounds = (
            growth_rate_per_h, growth_rate_per_h)

        drain = _drain_stoichiometry(model, spec, formulation)

        sink = cobra.Reaction("YSTWIN_PRECURSOR_SINK")
        sink.lower_bound, sink.upper_bound = 0.0, 1000.0
        sink.add_metabolites(drain)
        model.add_reactions([sink])
        model.objective = sink.id
        # Infeasibility is an EXPECTED outcome here -- it is what `InfeasibleState` reports --
        # so cobra's status warning is not a defect and must not trip the suite's
        # `filterwarnings = ["error"]`. Caught around the solve alone, by message, so any
        # other warning still surfaces.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Solver status is", UserWarning)
            solution = model.optimize()
        if solution.status != "optimal" or solution.objective_value is None:
            raise InfeasibleState(
                f"{spec.product} at mu = {growth_rate_per_h:g} /h on {carbon_source} is "
                f"infeasible under {formulation.label()}; the model cannot hold that growth "
                f"rate on that carbon supply. Raise the matched carbon or lower the growth "
                f"rate -- do not read this as a zero yield, which is a different claim")
        return float(solution.objective_value)


def yield_ratio(model, spec, growth_rate_per_h: float,
                formulation: Formulation = DEFAULT_FORMULATION) -> GemEnvironmentResult:
    """Ethanol-over-glucose yield for this pathway's precursor, at a fixed growth rate.

    Args:
        model: A cobra model. NOT mutated -- every change is inside a ``with model:`` block.
        spec: A :class:`~ystwin.pathway.spec.PathwaySpec`. Only its ``precursor_metabolite``
            and ``precursor_stoichiometry`` are read, which is what makes this
            product-generic: a pathway with no calibration, no kinetics and no measured state
            still has a precursor, so it still gets a descriptor.

    Raises:
        InfeasibleState: when the growth rate cannot be held on the matched carbon supply.
        ValueError: on a non-positive growth rate.
    """
    if growth_rate_per_h <= 0:
        raise ValueError(f"growth rate must be positive, got {growth_rate_per_h}")
    return GemEnvironmentResult(
        product=spec.product,
        precursor_metabolite=spec.precursor_metabolite,
        growth_rate_per_h=growth_rate_per_h,
        glucose_max_mmol_per_gdcw_h=_max_precursor(
            model, spec, "glucose", growth_rate_per_h, formulation),
        ethanol_max_mmol_per_gdcw_h=_max_precursor(
            model, spec, "ethanol", growth_rate_per_h, formulation),
        formulation=formulation)


def max_growth_rate(model, carbon_source: str,
                    formulation: Formulation = DEFAULT_FORMULATION) -> float:
    """The largest growth rate this feed and supply can hold, /h. Zero if none can.

    The feasibility boundary of every LP in this module, in the units the caller states its
    question in. It is the SAME constraint set as :func:`_max_precursor` with the biomass
    bound removed and biomass made the objective, so ``growth_rate <= max_growth_rate`` is
    exactly the condition under which a yield can be asked for at all -- a property
    `tests/test_gem_environment.py` asserts rather than assumes.
    """
    with model:
        glucose = model.reactions.get_by_id(formulation.glucose_exchange)
        ethanol = model.reactions.get_by_id(formulation.ethanol_exchange)
        glucose.lower_bound = 0.0
        ethanol.lower_bound = 0.0
        supply = formulation.matched_carbon_cmol_per_gdcw_h * 1000.0
        if carbon_source == "glucose":
            glucose.lower_bound = -supply / _CARBON_ATOMS["glucose"]
        elif carbon_source == "ethanol":
            ethanol.lower_bound = -supply / _CARBON_ATOMS["ethanol"]
        else:
            raise ValueError(
                f"carbon source must be 'glucose' or 'ethanol', got {carbon_source!r}")
        model.reactions.get_by_id(
            formulation.oxygen_exchange).lower_bound = formulation.oxygen_lower_bound
        model.objective = formulation.biomass_reaction
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Solver status is", UserWarning)
            solution = model.optimize()
        if solution.status != "optimal" or solution.objective_value is None:
            return 0.0
        return float(solution.objective_value)


def _growth_headroom(model, growth_rate_per_h, formulation) -> float:
    """``1 - mu / mu_max`` on the tighter feed, floored at zero. See :data:`MIN_GROWTH_HEADROOM`."""
    headrooms = []
    for carbon_source in ("glucose", "ethanol"):
        ceiling = max_growth_rate(model, carbon_source, formulation)
        headrooms.append(0.0 if ceiling <= 0.0
                         else max(0.0, 1.0 - growth_rate_per_h / ceiling))
    return min(headrooms)


@dataclass(frozen=True)
class SignStability:
    """What a sign test COVERED, not only what it concluded.

    Added 2026-09-04, replacing the bare ``(stable, low, high)`` tuple that
    :func:`sign_is_stable` used to return. The tuple was the second instance of the finding
    commit cf34701 recorded for the carbon supply: a robustness claim quantified over one
    parameter while a second one that moves the answer at least as much sits unswept and
    unreported. Here the second parameter is the GROWTH RATE, and it entered through the back
    door -- the old loop ``continue``d past an infeasible supply point, so "stable across
    4-25 mmol C/gDCW/h" could be, and was, a verdict reached on four of six points, every one
    of them near the boundary where the ratio diverges.

    Every field below exists so a caller can see that without guessing which parameter was
    forgotten: how many points were asked for, how many were actually evaluated, and how
    close to its own feasibility boundary the closest one sat.
    """

    #: Every evaluated ratio on the same side of 1.0. Says NOTHING about coverage: read it
    #: with :attr:`swept_the_whole_range` and :attr:`min_growth_headroom`, or through
    #: :meth:`is_a_prediction`, which combines all three.
    stable: bool
    ratio_low: float
    ratio_high: float
    points_requested: int
    points_feasible: int
    supply_range_cmol_per_gdcw_h: tuple[float, float]
    #: The smallest ``1 - mu/mu_max`` over the REQUESTED points, not the evaluated ones, so an
    #: infeasible point drives this to exactly zero rather than vanishing from the summary.
    min_growth_headroom: float
    growth_rate_per_h: float
    formulation: Formulation

    @property
    def swept_the_whole_range(self) -> bool:
        return self.points_feasible == self.points_requested

    def is_a_prediction(self,
                        min_growth_headroom: float = MIN_GROWTH_HEADROOM) -> bool:
        """The gate. A sign is a prediction only if the sweep that found it was complete.

        Three conditions, and dropping any one of them was a way this module has already been
        wrong: the sign has to hold (that is :attr:`stable`), the sweep has to have covered
        the range it names (:attr:`swept_the_whole_range`), and no point may sit so close to
        the feasibility boundary that the boundary is what set its ratio
        (:data:`MIN_GROWTH_HEADROOM`).
        """
        return (self.stable and self.swept_the_whole_range
                and self.min_growth_headroom >= min_growth_headroom)

    def summary(self) -> str:
        low, high = self.supply_range_cmol_per_gdcw_h
        return (f"mu = {self.growth_rate_per_h:.4g} /h over "
                f"{1000 * low:.1f}-{1000 * high:.1f} mmol C/gDCW/h: "
                f"{self.points_feasible}/{self.points_requested} points feasible, "
                f"ratio {self.ratio_low:.4g}-{self.ratio_high:.4g}, "
                f"min growth headroom {self.min_growth_headroom:.3f}"
                f" -> {'prediction' if self.is_a_prediction() else 'NOT a prediction'}")


def sign_is_stable(model, spec, growth_rate_per_h: float,
                   supply_range_cmol_per_gdcw_h: tuple[float, float] = (4.0e-3, 25.0e-3),
                   points: int = 6,
                   formulation: Formulation = DEFAULT_FORMULATION) -> SignStability:
    """Does this product's ethanol-vs-glucose sign survive the plausible carbon supply?

    A sign that reverses inside the range of supplies the source experiment actually ran at
    is not a prediction, it is a property of one unswept parameter. This exists because that
    is exactly what happened: the committed glycogen and gadusol rows were glucose-favoured
    at the shipped 20.0 mmol C/gDCW/h and ethanol-favoured at Kocharin's measured ~4.0, and
    the flip was reported as the module's headline falsifiable content.

    The default range spans the measured uptakes (3.99-4.52 at mu = 0.05) up to 25, which is
    where the previously-shipped value sat.

    Returns:
        A :class:`SignStability`. It carries the sweep's own provenance -- how many points
        were requested against how many were feasible, and the smallest growth headroom
        among them -- because until 2026-09-04 this returned a bare ``(stable, low, high)``
        and silently ``continue``d past every infeasible point. Ask it
        :meth:`~SignStability.is_a_prediction` rather than reading ``stable`` alone.

    Raises:
        InfeasibleState: when NO requested point is feasible. Never returned as
            ``stable=False``: "the vessel cannot run" and "the sign reverses" are different
            claims, and :class:`InfeasibleState`'s own docstring forbids conflating them.
            The caller that used to receive ``(False, nan, nan)`` here wrote it out as
            ``is_a_prediction=False`` under the banner "sign reverses inside that range".
    """
    low, high = supply_range_cmol_per_gdcw_h
    ratios = []
    headrooms = []
    for index in range(points):
        supply = low + (high - low) * index / (points - 1)
        candidate = Formulation(
            matched_carbon_cmol_per_gdcw_h=supply,
            oxygen_lower_bound=formulation.oxygen_lower_bound,
            glucose_exchange=formulation.glucose_exchange,
            ethanol_exchange=formulation.ethanol_exchange,
            oxygen_exchange=formulation.oxygen_exchange,
            biomass_reaction=formulation.biomass_reaction,
            precursor_sink_only=formulation.precursor_sink_only)
        headrooms.append(_growth_headroom(model, growth_rate_per_h, candidate))
        try:
            ratios.append(yield_ratio(model, spec, growth_rate_per_h, candidate).ratio)
        except InfeasibleState:
            continue
    if not ratios:
        raise InfeasibleState(
            f"{spec.product} at mu = {growth_rate_per_h:g} /h is infeasible at every one of "
            f"the {points} supplies in {1000 * low:g}-{1000 * high:g} mmol C/gDCW/h under "
            f"{formulation.label()}. That is not an unstable sign -- it is no sign at all, "
            f"and reporting it as one would record 'this vessel cannot run' as 'this feed is "
            f"bad for this product'")
    stable = all(r > 1.0 for r in ratios) or all(r < 1.0 for r in ratios)
    return SignStability(
        stable=stable,
        ratio_low=min(ratios),
        ratio_high=max(ratios),
        points_requested=points,
        points_feasible=len(ratios),
        supply_range_cmol_per_gdcw_h=(low, high),
        min_growth_headroom=min(headrooms),
        growth_rate_per_h=growth_rate_per_h,
        formulation=formulation)


def environment_descriptor(model, spec, ethanol_carbon_fraction: float,
                           growth_rate_per_h: float,
                           formulation: Formulation = DEFAULT_FORMULATION) -> float:
    """The descriptor to put in a flux law: log of the yield at this feed, relative to glucose.

    Zero on pure glucose by construction, so the descriptor measures the feed's departure from
    the reference rather than its absolute quality -- which is what keeps it on the same
    footing as the empirical ``z`` it replaces.

    Interpolated linearly in the ethanol carbon fraction between the two pure feeds rather
    than solved at the blend. That is a stated approximation, not an oversight: solving the
    blend directly needs a second matched-carbon split and doubles the LP count, and on the
    one blend Kocharin measured the interpolation is what this module's own
    leave-one-carbon-source-out score -- the `GEM carbon contrast` row, at
    <!-- audit:value table=outputs/gem_environment_scores.csv column=leave_one_carbon_source_out_fold_error row="law=GEM carbon contrast" -->2.1942x --
    was computed with.

    *That sentence cited "the 1.462x leave-one-carbon-source-out score" until 2026-09-04. No
    cell in `outputs/gem_environment_scores.csv` holds 1.462 under any configuration -- the
    five values are 1.7237, 1.4683, 2.1942, 1.3439 and 1.7534 -- and it does not match the
    pre-carrier-fix values either. It was deleted rather than corrected, because a number
    whose provenance cannot be recovered cannot be repaired, and the cell the sentence MEANT
    is named instead.*
    """
    import math

    if not 0.0 <= ethanol_carbon_fraction <= 1.0:
        raise ValueError(
            f"ethanol carbon fraction is in [0, 1], got {ethanol_carbon_fraction}")
    if ethanol_carbon_fraction == 0.0:
        return 0.0
    result = yield_ratio(model, spec, growth_rate_per_h, formulation)
    return ethanol_carbon_fraction * math.log(result.ratio)
