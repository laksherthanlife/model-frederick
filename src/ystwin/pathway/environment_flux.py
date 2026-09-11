"""The environment channel the chain never had: carbon source moving titre at FIXED growth rate.

WHY THIS EXISTS. Everything else in this package reaches the product through `mu`. That was
measured directly (`scripts/product_environment_sweep.py`): holding mu at 0.18 /h and varying
carbon source, temperature, oxygen and stressor over 62 environments returns ONE content
value, and re-predicting 92 batch environments at their own mu with every other axis bypassed
disagrees by 6.4e-07 mg/gDCW. Content was a strict function of mu.

That is false about yeast, and this module carries the measurement that shows it.
Everything numeric below is produced by `scripts/environment_channel.py` into
`outputs/environment_channel.csv` and `outputs/environment_law_scores.csv`.

THE MEASUREMENT. Kocharin 2013 (`data/phb/kocharin2013_chemostat_states.tsv`): ONE genotype
(SCKK006), ONE vessel, eleven aerobic chemostat steady states, four dilution rates x three
carbon feeds matched at 0.666 Cmol/L. At mu = 0.05 /h, PHB content is
<!-- audit:value table=outputs/environment_channel.csv column=glucose_mg_per_gdw row="growth_rate_per_h=0.05;carbon_source=ethanol" -->4.33 mg/gDW on glucose and
<!-- audit:value table=outputs/environment_channel.csv column=test_mg_per_gdw row="growth_rate_per_h=0.05;carbon_source=ethanol" -->16.55 on ethanol -- a
<!-- audit:value table=outputs/environment_channel.csv column=fold_vs_glucose row="growth_rate_per_h=0.05;carbon_source=ethanol" -->3.822x change at IDENTICAL growth rate, with a
separation of <!-- audit:value table=outputs/environment_channel.csv column=separation_sd row="growth_rate_per_h=0.05;carbon_source=ethanol" -->63.96 standard deviations. Genotype
fixed, vessel fixed, carbon matched, mu equal. Nothing about that is a fitted claim, and the
chain answers it with zero.

AND THE EFFECT REVERSES WITH GROWTH RATE, which is the part a first reading misses and the
reason this module fits an interaction rather than a coefficient::

    mu = 0.05  ethanol <!-- audit:value table=outputs/environment_channel.csv column=fold_vs_glucose row="growth_rate_per_h=0.05;carbon_source=ethanol" -->3.822x glucose
    mu = 0.10  ethanol <!-- audit:value table=outputs/environment_channel.csv column=fold_vs_glucose row="growth_rate_per_h=0.1;carbon_source=ethanol" -->2.413x glucose
    mu = 0.15  ethanol <!-- audit:value table=outputs/environment_channel.csv column=fold_vs_glucose row="growth_rate_per_h=0.15;carbon_source=ethanol" -->0.892x glucose  -- reversed
    mu = 0.20  mix     <!-- audit:value table=outputs/environment_channel.csv column=fold_vs_glucose row="growth_rate_per_h=0.2;carbon_source=glucose_ethanol_1to2" -->0.816x glucose

A single scalar `c` in ``exp(c*z)`` would predict one direction at every growth rate and be
wrong above the crossover. Leave-one-state-out on the eleven states, median fold error::

    content constant  <!-- audit:value table=outputs/environment_law_scores.csv column=loo_median_fold_error row="law=content constant" -->1.6452x
    mu only           <!-- audit:value table=outputs/environment_law_scores.csv column=loo_median_fold_error row="law=mu only" -->1.5052x   (the chain today)
    z only            <!-- audit:value table=outputs/environment_law_scores.csv column=loo_median_fold_error row="law=z only" -->1.3438x
    mu + z            <!-- audit:value table=outputs/environment_law_scores.csv column=loo_median_fold_error row="law=mu + z" -->1.4407x
    mu + z + mu*z     <!-- audit:value table=outputs/environment_law_scores.csv column=loo_median_fold_error row="law=mu + z + mu*z" -->1.3201x   (this module)

The interaction is not decoration: dropping it is worse than dropping mu.

**BUT DO NOT EXTRAPOLATE IT TO A FEED IT WAS NOT FITTED ON.** Scored leave-one-CARBON-SOURCE
-out -- predicting a feed never seen, which is the bar a new medium would face -- the ranking
inverts and the four-parameter law is the worst of the three::

                          k   held-out STATE      held-out FEED, median and worst
    mu only               2   1.5052x / 4.0985x   <!-- audit:value table=outputs/environment_law_scores.csv column=loco_median_fold_error row="law=mu only" -->1.7237x / <!-- audit:value table=outputs/environment_law_scores.csv column=loco_worst_fold_error row="law=mu only" -->4.4928x
    z only                2   1.3438x / 2.5825x   <!-- audit:value table=outputs/environment_law_scores.csv column=loco_median_fold_error row="law=z only" -->1.3439x / <!-- audit:value table=outputs/environment_law_scores.csv column=loco_worst_fold_error row="law=z only" -->2.3235x
    mu + z                3   1.4407x / 3.3005x   <!-- audit:value table=outputs/environment_law_scores.csv column=loco_median_fold_error row="law=mu + z" -->1.4748x / <!-- audit:value table=outputs/environment_law_scores.csv column=loco_worst_fold_error row="law=mu + z" -->5.4435x
    mu + z + mu*z         4   1.3201x / 2.5481x   <!-- audit:value table=outputs/environment_law_scores.csv column=loco_median_fold_error row="law=mu + z + mu*z" -->1.4564x / <!-- audit:value table=outputs/environment_law_scores.csv column=loco_worst_fold_error row="law=mu + z + mu*z" -->5.6051x

The worst case more than DOUBLES, 2.3235x to 5.6051x. Both facts are true at once and neither
cancels the other: the growth-rate reversal is measured and real, and eleven states cannot
estimate it well enough to survive extrapolation to a fourth feed.

*The held-out-FEED column was eight bare literals until 2026-09-04, produced by no script in
this repository -- the computation existed only inside `tests/test_environment_flux.py`, which
asserted ORDERINGS between them and never their values. They were all correct, so this is a
hazard closed rather than a number changed: `scripts/environment_channel.py` now writes
``loco_median_fold_error`` and ``loco_worst_fold_error`` per law and the eight point at cells.*

So the shipped law is the right one for glucose and ethanol, which is what it was fitted on
and what :func:`environment_factor` is asked for, and the two-parameter ``z``-only law is the
one to prefer if a fourth GLUCOSE:ETHANOL RATIO is ever run -- at which point four levels take
the permutation floor from 1/3! = 0.1667 to 1/4! = 0.0417 and the coefficients become earnable
for the first time.

**A GLYCEROL OR ACETATE FEED WOULD NOT DO THAT, and this paragraph said it would until
2026-09-04.** The prescription was written from the descriptor's docstring rather than from
its signature. :func:`carbon_descriptor` takes ``(glucose_g_per_l, ethanol_g_per_l)`` and
nothing else: a glycerol-only feed raises ``ValueError("the feed contains no carbon")``, and
glucose plus glycerol returns z = 0.0 exactly, so the z-only law would predict a glycerol
state as though it were pure glucose and :func:`environment_factor` would return exactly 1.0.
Adding those feeds needs a NEW DESCRIPTOR, not new rows -- specifically a wider
:func:`carbon_descriptor` signature and entries in :data:`_CARBON_ATOMS` and
:data:`_MOLAR_MASS_G_PER_MOL` for every species the new feed carries, plus a re-derivation of
what "ethanol carbon fraction" generalises to when there are more than two carbon species.
A wrong number is recoverable; a wrong prescription aims the next experiment, which is worse.
`tests/test_environment_flux.py` now pins the descriptor's coverage against
:data:`MEASURED_ENVIRONMENT_RESPONSES`, so widening this prose without widening the descriptor
fails.

WHAT THIS IS NOT, and both limits are structural rather than fixable by more analysis.

1. **The COEFFICIENT is not established, though the CHANNEL is.** There are three carbon
   sources, so the exchangeable unit is the carbon source and the smallest attainable
   permutation p is 1/3! = <!-- audit:value table=outputs/environment_law_scores.csv column=permutation_floor row="law=z only" -->0.1667 -- which is exactly what
   a relabelling test returns (<!-- audit:value table=outputs/environment_law_scores.csv column=permutation_p row="law=z only" -->0.1667). No
   analysis of these eleven states can earn significance for the law. What does not need
   earning is the raw contrast: 4.33 +/- 0.19 against 16.55 +/- 0.02 at the same mu is 64
   standard deviations and requires no model at all. **The existence of the channel is
   measured; its functional form is fitted and unearned.**
   Note also that holding out a whole carbon source means extrapolating z from two points to
   a third, so leave-one-carbon-source-out is not a validation of the form either.

2. **It is fitted on ONE product.** PHB is acetyl-CoA-limited and ethanol feeds cytosolic
   acetyl-CoA directly through ACS, bypassing the PDH bypass -- a PATHWAY mechanism, not a
   host one. So the coefficients here may not transfer, and this module REFUSES to apply
   them to a product that has not measured its own. `beta_carotene` has not, and no
   carotenoid row in `data/carotenoid_batch/published_batch_titres.tsv` could fit one: the
   fitted descriptor is the feed's ETHANOL CARBON FRACTION, and not one of those rows has
   ethanol in the feed, so z = 0 on every one of them and the coefficient is unidentified by
   construction. That is computed rather than remembered --
   :mod:`ystwin.pathway.published_cassettes` parses each row's medium against a declared
   carbon vocabulary in which ``ethanol`` appears, so the absence is a measurement and not a
   vocabulary gap, and `tests/test_published_cassettes.py` asserts both halves.

   *That sentence said "every carotenoid row is glucose or galactose" until 2026-09-04, which
   is false -- the table carries sucrose, xylose, olive oil, oleic acid, grape juice and
   glucose+acetate rows. The CONCLUSION survives on a different and sharper reason, and the
   sharper one is also robust to a xylose row being added, which the enumeration was not.*
   The GEM's opinion, for what little it is worth: in 5 of 6 stated regimes all four
   precursor taps rank ethanol < mix < glucose, i.e. the OPPOSITE direction to the PHB
   measurement, and the coefficient it implies reverses sign on changing the matched carbon
   uptake from 1.0 to 2.0 cmol/gDCW/h. The GEM is demonstrably wrong about the one axis
   where a measurement exists, so it is not evidence -- but it is not evidence FOR
   transfer either.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "KOCHARIN_PHB",
    "MEASURED_ENVIRONMENT_RESPONSES",
    "EnvironmentResponse",
    "EnvironmentUnmeasured",
    "carbon_descriptor",
    "environment_factor",
]

#: Carbon atoms per molecule and gram-molecular mass, for the two feeds Kocharin used.
#: Kept explicit rather than importing a formula parser: two sugars do not justify one, and a
#: wrong molar mass here would silently rescale the descriptor.
_CARBON_ATOMS = {"glucose": 6.0, "ethanol": 2.0}
_MOLAR_MASS_G_PER_MOL = {"glucose": 180.156, "ethanol": 46.068}


class EnvironmentUnmeasured(ValueError):
    """This product has no measured environment response, so the chain will not invent one.

    Its own type because the correct caller response is specific: report the channel as
    unmodelled and say what would measure it, rather than falling back to a coefficient
    borrowed from another pathway. Borrowing is exactly the error the entry-flux law already
    refuses ("applying one gene's scalar to another gene's expression").
    """


@dataclass(frozen=True)
class EnvironmentResponse:
    """A fitted environment -> content response, with the interaction that makes it real.

    ``log content = intercept + b_mu*log(mu) + b_z*z + b_muz*log(mu)*z``

    Only the two z terms live here; the mu terms belong to the pathway solve, which already
    has them. What this object contributes is the part of the answer mu cannot reach.
    """

    product: str
    b_z: float
    b_mu_z: float
    growth_rate_range: tuple[float, float]
    n_states: int
    n_carbon_sources: int
    source: str

    @property
    def sign_flip_growth_rate_per_h(self) -> float:
        """The mu at which the carbon source stops helping and starts hurting.

        ``b_z + b_muz*log(mu) = 0``. On the PHB fit this is 0.1718 /h, which sits INSIDE the
        measured window -- so a model without the interaction is wrong on real data, not
        merely at an extrapolated edge.
        """
        if self.b_mu_z == 0.0:
            return float("inf")
        return math.exp(-self.b_z / self.b_mu_z)

    @property
    def permutation_floor(self) -> float:
        """Smallest p any relabelling test could return, given the carbon sources available.

        1/k! for k carbon sources. At k = 3 this is 0.1667, so the coefficients CANNOT be
        established from this dataset however they are analysed. Carried as a property so a
        caller quoting the law also has the reason it is not significant.
        """
        return 1.0 / math.factorial(self.n_carbon_sources)

    def log_factor(self, ethanol_carbon_fraction: float, growth_rate_per_h: float) -> float:
        if not 0.0 <= ethanol_carbon_fraction <= 1.0:
            raise ValueError(
                f"ethanol carbon fraction is in [0, 1], got {ethanol_carbon_fraction}")
        if growth_rate_per_h <= 0.0:
            raise ValueError(f"growth rate must be positive, got {growth_rate_per_h}")
        return ethanol_carbon_fraction * (self.b_z
                                          + self.b_mu_z * math.log(growth_rate_per_h))

    def summary(self) -> str:
        low, high = self.growth_rate_range
        return (f"{self.product}: d log content/dz = {self.b_z:+.3f} "
                f"{self.b_mu_z:+.3f}*log(mu), sign flip at mu = "
                f"{self.sign_flip_growth_rate_per_h:.4f} /h, fitted on {self.n_states} states "
                f"across {self.n_carbon_sources} carbon sources "
                f"(permutation floor p = {self.permutation_floor:.4f})")


#: Fitted on the eleven Kocharin states by ordinary least squares in log content, with mu and
#: the interaction present. The two coefficients are checked against
#: `outputs/environment_law_scores.csv` by `tests/test_environment_flux.py` rather than by an
#: audit marker: a marker is an HTML comment and only parses inside a docstring, so a
#: constant that must stay executable is pinned by a test instead. Reproduced independently twice during the audit that produced this
#: module; the leave-one-state-out median fold error of the full model is 1.320x.
KOCHARIN_PHB = EnvironmentResponse(
    product="phb",
    b_z=-2.2021,
    b_mu_z=-1.2502,
    growth_rate_range=(0.05, 0.20),
    n_states=11,
    n_carbon_sources=3,
    source="Kocharin 2013, data/phb/kocharin2013_chemostat_states.tsv; "
           "one genotype (SCKK006), one vessel, carbon matched at 0.666 Cmol/L",
)

#: Keyed by product. A product absent from this mapping has no measured response, and
#: :func:`environment_factor` refuses rather than borrowing one.
MEASURED_ENVIRONMENT_RESPONSES = {"phb": KOCHARIN_PHB}


def carbon_descriptor(glucose_g_per_l: float, ethanol_g_per_l: float) -> float:
    """The ethanol share of feed CARBON, which is the descriptor the response is fitted on.

    Carbon rather than mass or molarity, because the Kocharin feeds were matched on carbon
    (0.666 Cmol/L) and matching is what makes the contrast a carbon-QUALITY comparison
    instead of a carbon-QUANTITY one. A descriptor on grams would confound the two.
    """
    if glucose_g_per_l < 0 or ethanol_g_per_l < 0:
        raise ValueError(
            f"feed concentrations must be non-negative, got {glucose_g_per_l} and "
            f"{ethanol_g_per_l}")
    glucose_c = (glucose_g_per_l / _MOLAR_MASS_G_PER_MOL["glucose"]
                 * _CARBON_ATOMS["glucose"])
    ethanol_c = (ethanol_g_per_l / _MOLAR_MASS_G_PER_MOL["ethanol"]
                 * _CARBON_ATOMS["ethanol"])
    total = glucose_c + ethanol_c
    if total <= 0:
        raise ValueError("the feed contains no carbon; there is no descriptor to compute")
    return ethanol_c / total


def environment_factor(product: str,
                       ethanol_carbon_fraction: float,
                       growth_rate_per_h: float) -> float:
    """Fold change in content from the carbon source, at this growth rate. 1.0 means no effect.

    Raises:
        EnvironmentUnmeasured: when the product has no measured response. This is the common
            case and it is deliberate -- `beta_carotene` has none, and the correct behaviour
            is to report the channel as unmodelled rather than to apply PHB's coefficients to
            an isoprenoid pathway whose limitation is different.
    """
    response = MEASURED_ENVIRONMENT_RESPONSES.get(product)
    if response is None:
        raise EnvironmentUnmeasured(
            f"no measured environment response for {product!r}; have "
            f"{sorted(MEASURED_ENVIRONMENT_RESPONSES)}. The carbon-source channel is REAL -- "
            f"it moves PHB content 3.82x at a fixed growth rate of 0.05 /h, 64 standard "
            f"deviations, one genotype and one vessel -- but its coefficients were fitted on "
            f"an acetyl-CoA-limited pathway fed by ethanol through ACS, which is a pathway "
            f"mechanism and not a host one. Applying them here would be the same error as "
            f"applying one gene's flux scalar to another gene's expression. What lifts this: "
            f"content for {product!r} at two carbon sources reaching the SAME growth rate by "
            f"different routes -- see docs/PROTOCOLS.md")
    return math.exp(response.log_factor(ethanol_carbon_fraction, growth_rate_per_h))
