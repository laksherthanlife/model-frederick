"""Walk a declared pathway at steady state: flux and growth rate in, every pool out.

One balance per node, and it is the whole model:

    d[X]/dt = v_in - v_out(X) - mu*[X] = 0

Every rate law in :class:`~ystwin.pathway.spec.RateLaw` closes that in one line, so the
chain is solved by walking it in order rather than by integrating anything. No solver, no
tolerance, no initial guess.

    saturating    v_out = vmax*X/(km + X)   ->   mu*X^2 + (mu*km + vmax - v_in)*X - v_in*km = 0
    proportional  v_out = k*X               ->   X = v_in/(k + mu)
    passthrough   v_out = v_in              ->   X = 0
    terminal      v_out = 0                 ->   X = v_in/mu

The saturating quadratic has exactly one positive root, which is the one taken. That an
INTERMEDIATE is pinned by a mass balance rather than by an objective is the reason this
layer can predict a pool the GEM cannot: an objective bounds a flux from above, and a
balance fixes it.

The terminal line above is the DILUTED case, and it is the whole module only while every
product stays in the cell. ``X = v_in/mu`` says the product accumulates until washout removes
it as fast as the pathway makes it. A terminal node may also declare two further outlets,
each of which shifts that line and neither of which invents a new algebraic form:

    degraded   v_loss = k_deg*X          ->  X = v_in/(mu + k_deg)
    secreted   v_loss = vmax*X/(km + X)  ->  L*X^2 + (L*km + vmax - v_in)*X - v_in*km = 0

with ``L = mu + k_deg``. The second is the saturating quadratic above under a shifted linear
coefficient, so the walk gains an outlet without gaining a solver.

Both are refused at SOLVE time when their constants are missing, and that is the line this
module draws: a pathway whose chemistry it cannot describe is refused at load by
:class:`PathwaySpec`, and a pathway one measured number short is refused here, by name. A
SECRETED product used to be the first kind and is now the second.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .spec import Fate, PathwaySpec, RateLaw

__all__ = [
    "DRY_WEIGHT_G_PER_GDCW",
    "CeilingRefutation",
    "ContentCeiling",
    "ImplausibleContent",
    "NodeKinetics",
    "NodeState",
    "PathwaySolution",
    "content_ceiling",
    "solve_pathway",
]

# One gram of anything per gram of dry cell weight is the whole cell. This is arithmetic,
# not biology, and it needs no citation: a pool at this value is the cell made entirely of
# that pool, with no protein, no membrane and no nucleic acid.
#
# It is used as a REFUSAL rather than a warning and the ceiling is deliberately absurd. A
# tighter bound would be a claim about how much product a yeast can hold, which is a real
# and contested number -- Lopez 2019 reports 21 mg/gDCW total carotenoid, Ferreira 2018
# reaches 218 mg/gDCW of TAG after deleting the lipases -- and asserting one here would be a
# Tier 0 constant smuggled into a solver. A spec that knows its own ceiling can say so.
DRY_WEIGHT_G_PER_GDCW = 1.0


class ImplausibleContent(ValueError):
    """A solved pool exceeds the mass of the cell holding it.

    Raised rather than returned, because the alternative is what this replaced: at an entry
    expression of 1e4 the lycopene pool came back as 36,093 mg/gDCW -- 3,609% of dry weight
    -- with no note and no flag, as an ordinary float a caller would put in a table. The
    steady-state solve has exactly one positive root and no opinion about whether that root
    describes anything.

    It means the INPUTS are wrong, not the arithmetic. Usually an entry expression far
    outside what the flux law was fitted on, or a growth rate near zero where the washout
    term that holds every pool down goes to nothing.
    """


@dataclass(frozen=True)
class NodeKinetics:
    """Kinetic parameters for one node's consuming step.

    Args:
        vmax_per_growth: For a saturating step, ``vmax = vmax_per_growth * mu``. Carried
            per unit growth rather than absolute because that is the form the data
            supports -- see docs/research/KINETIC_FIT.md. Its 95% profile interval on the
            growth exponent is [0.53, 1.78], so the linearity is bounded, not pinned.

            ON THE STEP FEEDING A DILUTED TERMINAL NODE THIS CONSTANT IS ALSO A CEILING, and
            numerically the same one: ``mu`` cancels out of ``content <= vmax_per_growth``.
            A calibration whose ceiling a measurement has already passed should ship it as a
            :class:`ContentCeiling` carrying a :class:`CeilingRefutation`, so
            :func:`content_ceiling` can hand the refutation on to whoever reads the bound.
            A plain ``float`` here means no refutation is declared, which is the right
            default and not a claim that none exists.
        km: Half-saturation content, mmol/gDCW. Saturating steps only.
        rate_constant: First-order constant, /h. Proportional steps only.
        degradation_rate_per_h: First-order turnover by a native enzyme the pathway does
            not include, /h. **Required for, and only meaningful on, a node whose fate is
            ``DEGRADED``.** A storage compound has a third outlet besides the next step and
            growth, and omitting it is not a small error: Ferreira 2018 raises
            triacylglycerol from 129 to 218 mg/gDCW by deleting the lipases alone.

            First order is an approximation and a declared one. GPH1 and SGA1 are enzymes
            and saturate; treating their combined action as ``k*X`` is right while the pool
            sits well below their Km and wrong above it. It is used because a saturating
            outlet needs a vmax and a Km that nobody has measured for these enzymes in a
            producing strain, and inventing two constants to avoid approximating with one
            is the worse trade.
        secretion_vmax_mmol_per_gdcw_h: Export capacity of the terminal node, mmol/gDCW/h.
            **Required, with the Km, for a node whose fate is ``SECRETED``.** ABSOLUTE, not
            per unit growth -- unlike ``vmax_per_growth``, which carries a ``mu`` scaling
            because the carotenoid data supports one. No dataset supports one for an
            exporter, so this constant is a claim at the dilution rate it was measured at
            and nowhere else, which is why ``growth_rate_range`` is required beside it.
        secretion_km_mmol_per_gdcw: Content at which export runs at half capacity,
            mmol/gDCW. **Required together with the vmax**, and the pairing is not
            fastidiousness: one measured ``([X], v_sec)`` state fixes only the ratio
            ``vmax/(km + X)`` and leaves both constants free along a ridge. Two states at
            different intracellular contents invert exactly.

            SATURATING RATHER THAN FIRST ORDER, AND THIS IS A MEASURED CHOICE. The
            paragraph above defends approximating the degradation outlet as ``k*X``, on the
            grounds that inventing two constants to avoid approximating with one is the
            worse trade. That argument does not carry over here, because for export the
            curvature is not hypothetical -- it is the observation. Kastberg 2025 (PMID
            39971732) measured intracellular and secreted product on the SAME K. phaffii
            chemostat samples at D = 0.1 /h: the human insulin precursor rises 4.12 log2FC
            (17.4x) inside the cell between I1G and I1S while the secretome does not move
            at all, and Mambalgin-1 rises 4.15 log2FC (17.8x) inside against 1.8 log2FC
            (3.5x) secreted. A first-order outlet secretes in proportion to the pool BY
            CONSTRUCTION, so no value of ``k`` reproduces a 17x pool beside a flat
            supernatant. A pool piling up while its outlet does not is what a capacity
            looks like.

            WHAT THE LUMPED CONSTANT ACTUALLY COVERS, stated because it is a real
            limitation and not a caveat. Pfeffer 2011 (PMID 21703020) pulse-labelled a
            K. phaffii chemostat at the same D = 0.1 /h and split the fate of secretory
            product three ways: 58% degraded intracellularly, 35% secreted, 7% inherited by
            daughter cells, with degradation the FASTER half-life (45.8 min against 75.3).
            So the dominant sink for a secretory protein is proteolysis rather than export,
            and one saturating outlet fitted to net appearance in the supernatant absorbs
            that degradation into its own vmax. Separating them needs a node this walk does
            not have.
        growth_rate_range: The ``(low, high)`` growth rates these parameters were fitted
            over, or ``None`` to allow any rate. **Supply it.**

            This field exists because its absence was a regression. The hand-written
            carotenoid solver refused outside its fitted range; replacing it with this
            generic one dropped the guard, and the chain then answered at D = 0.05 -- half
            the lower bound -- with no warning at all. `vmax = capacity * mu` is a law
            fitted over a 2.5-fold window, and 0.05 is outside it in the direction where
            washout stops dominating and the linearity has never been checked.

            ``None`` is allowed because a rate law fitted over no stated range is a real
            situation, and forcing a fake range would be worse than recording the absence.
    """

    vmax_per_growth: float | None = None
    km: float | None = None
    rate_constant: float | None = None
    degradation_rate_per_h: float | None = None
    secretion_vmax_mmol_per_gdcw_h: float | None = None
    secretion_km_mmol_per_gdcw: float | None = None
    growth_rate_range: tuple[float, float] | None = None

    def check_growth_rate(self, node: str, growth_rate: float) -> None:
        """Refuse a growth rate outside what these parameters were fitted over.

        Raises rather than warning, and the choice is deliberate: a warning on a returned
        number is read by whoever is looking and by nobody else, and the number goes into a
        table either way. The rate the caller most wants here -- batch, mu ~ 0.4 -- is the
        one furthest outside the fit.
        """
        if self.growth_rate_range is None:
            return
        low, high = self.growth_rate_range
        if not low <= growth_rate <= high:
            raise ValueError(
                f"node {node!r}: growth rate {growth_rate:.4g} /h is outside the range its "
                f"kinetics were fitted over, [{low:.4g}, {high:.4g}]. The rate law is "
                "`vmax = capacity * mu` fitted over that window and its growth exponent has "
                "a 95% interval of [0.53, 1.78] even inside it -- extrapolating a law that "
                "loose is not a small error. Measure at this rate, or widen the range "
                "deliberately and record that you did")

    def require_saturating(self, node: str) -> tuple[float, float]:
        if self.vmax_per_growth is None or self.km is None:
            raise ValueError(
                f"node {node!r} declares rate_law 'saturating' and needs both "
                "vmax_per_growth and km; got "
                f"vmax_per_growth={self.vmax_per_growth}, km={self.km}")
        if self.vmax_per_growth <= 0 or self.km <= 0:
            raise ValueError(f"node {node!r}: vmax_per_growth and km must be positive")
        return self.vmax_per_growth, self.km

    def require_proportional(self, node: str) -> float:
        if self.rate_constant is None:
            raise ValueError(
                f"node {node!r} declares rate_law 'proportional' and needs rate_constant")
        if self.rate_constant <= 0:
            raise ValueError(f"node {node!r}: rate_constant must be positive")
        return self.rate_constant


@dataclass(frozen=True)
class NodeState:
    """One node solved: what it holds and what leaves it.

    ``flux_out`` is what the NEXT STEP takes. ``loss_flux_mmol_per_gdcw_h`` is what leaves
    the chain entirely at this node -- native degradation, export -- and it exists because
    without it :attr:`PathwaySolution.carbon_closes` cannot be computed at all.

    That was not a hypothetical gap. ``carbon_closes`` subtracted only ``mu*[X]`` and called
    the remainder an identity that "must be zero", so from the day the ``DEGRADED`` outlet
    arrived it reported that outlet's flux as a mass-balance violation: on a storage chain at
    flux 0.5, mu 0.1, k_deg 0.4 it returned 0.4, which is 80% of the entry flux, as an
    ordinary float with no flag. Nothing caught it because every test of the identity ran on
    ``DILUTED`` chains, where this field is a literal ``0.0`` and the sum is exact.
    """

    name: str
    content_mmol_per_gdcw: float
    flux_in: float
    flux_out: float
    molar_mass_g_per_mol: float | None = None
    loss_flux_mmol_per_gdcw_h: float = 0.0

    @property
    def content_mg_per_gdcw(self) -> float | None:
        """The unit the literature reports, where a molar mass was declared."""
        if self.molar_mass_g_per_mol is None:
            return None
        return self.content_mmol_per_gdcw * self.molar_mass_g_per_mol


def _degradation_loss(node, parameters) -> float:
    """The native turnover rate for a node, or a refusal naming the number that is missing.

    Zero for an ordinary ``DILUTED`` node, so the balance is unchanged and every existing
    result is bit-for-bit what it was. For a ``DEGRADED`` one the rate is required: without
    it the solve would return ``flux/mu`` and silently drop the largest outlet, which is
    exactly the wrong answer the load-time refusal existed to prevent.
    """
    if node.fate != Fate.DEGRADED:
        return 0.0
    rate = None if parameters is None else parameters.degradation_rate_per_h
    if rate is None:
        raise ValueError(
            f"node {node.name!r} is declared DEGRADED -- a native enzyme consumes it -- and "
            f"no degradation_rate_per_h was supplied for it. Solving without one returns "
            f"content = flux/mu, which omits that outlet entirely; Ferreira 2018 raises TAG "
            f"from 129 to 218 mg/gDCW by deleting the lipases alone, so the omitted term can "
            f"be the larger one. Supply NodeKinetics(degradation_rate_per_h=...) measured on "
            f"a strain that carries the degrading enzymes, or declare the node DILUTED if it "
            f"genuinely has no enzymatic outlet")
    # `math.isfinite` and not `rate != rate`. The predicate here used to read
    # `not (rate >= 0.0) or rate != rate`, which promised finiteness in its message and did
    # not test it: `inf >= 0.0` is True and `inf != inf` is False, so an infinite turnover
    # passed, gave `loss = inf`, and returned content 0.0 -- a number that looks like a
    # legitimate answer -- while `carbon_closes` reported the entire entry flux as lost. The
    # second clause was also dead: `not (nan >= 0.0)` already catches NaN.
    if not math.isfinite(rate) or rate < 0.0:
        raise ValueError(
            f"node {node.name!r}: degradation_rate_per_h must be finite and non-negative, "
            f"got {rate}")
    return float(rate)


def _secretion_capacity(node, parameters) -> tuple[float, float] | None:
    """The ``(vmax, km)`` export pair for a node, ``None`` where the node does not export.

    ``None`` for an ordinary node, so the balance is unchanged and every existing result is
    bit-for-bit what it was -- the caller branches on ``None`` rather than passing a zero
    capacity through the quadratic, because a zero vmax reaches ``flux/mu`` only up to
    rounding and this module's inertness claim is exactness, not agreement.

    For a ``SECRETED`` node both constants are required, and so is the growth-rate window.
    This is stricter than the ``DEGRADED`` outlet asks and deliberately so: a degradation
    rate is first order and dimensionless per hour, while an export capacity is an absolute
    flux that does not carry its own ``mu``, so applying one outside the dilution rate it
    was fitted at is a category error rather than an extrapolation.
    """
    if node.fate != Fate.SECRETED:
        if parameters is not None and (
                parameters.secretion_vmax_mmol_per_gdcw_h is not None
                or parameters.secretion_km_mmol_per_gdcw is not None):
            # A supplied constant the solver would drop. `tests/test_pathway_spec.py` makes
            # this argument about the loader -- a field that is read and discarded is a
            # claim the reader believes -- and `NodeKinetics` is frozen with no
            # `__post_init__`, so the point of use is the only place it can be caught.
            raise ValueError(
                f"node {node.name!r} is declared {node.fate} and carries export kinetics "
                f"(secretion_vmax_mmol_per_gdcw_h / secretion_km_mmol_per_gdcw), which "
                f"only a SECRETED node reads. Solving would silently ignore them and "
                f"return the diluted answer. Declare the node SECRETED, or drop the "
                f"constants")
        return None

    vmax = None if parameters is None else parameters.secretion_vmax_mmol_per_gdcw_h
    km = None if parameters is None else parameters.secretion_km_mmol_per_gdcw
    if vmax is None or km is None:
        raise ValueError(
            f"node {node.name!r} is declared SECRETED -- it leaves through a transporter -- "
            f"and no export capacity was supplied for it "
            f"(secretion_vmax_mmol_per_gdcw_h={vmax}, secretion_km_mmol_per_gdcw={km}). "
            f"A secreted pool has no mu*[X] dilution term doing the work: it leaves at a "
            f"rate the cell sets, not one growth sets, and solving without the capacity "
            f"returns content = flux/mu, which is the answer for a product that stays in. "
            f"Pfeffer 2011 measured export at 5.5x the dilution flux for a secretory "
            f"protein at D = 0.1 /h, so the omitted term is the larger one. BOTH constants "
            f"are needed and one measurement will not give them: a single steady state "
            f"fixes only the ratio vmax/(km + X). The experiment is paired intracellular "
            f"content and supernatant accumulation rate on the SAME culture at TWO OR MORE "
            f"expression levels, with mu recorded -- see docs/MEASUREMENTS_NEEDED.md")
    for label, value in (("secretion_vmax_mmol_per_gdcw_h", vmax),
                         ("secretion_km_mmol_per_gdcw", km)):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(
                f"node {node.name!r}: {label} must be finite and positive, got {value}")
    if parameters.growth_rate_range is None:
        raise ValueError(
            f"node {node.name!r} supplies an export capacity with no growth_rate_range. "
            f"secretion_vmax_mmol_per_gdcw_h is an ABSOLUTE flux, not a capacity per unit "
            f"growth, so unlike vmax_per_growth it does not carry its own mu and is a claim "
            f"at the dilution rate it was measured at and nowhere else. Rebnegger 2024 "
            f"measured specific productivity falling 34-fold, 0.64 to 0.019 mg/gDCW/h, "
            f"between mu = 0.17 and mu < 0.0025. Supply the (low, high) window the pair was "
            f"fitted over")
    return float(vmax), float(km)


@dataclass(frozen=True)
class CeilingRefutation:
    """A published measurement that already exceeds the ceiling a calibration implies.

    Carried as the MEASUREMENT and not as a boolean ``refuted`` flag, so
    :attr:`ContentCeiling.refuted_fold` recomputes the verdict from the shipped constant
    every time it is asked. A boolean would go on saying "refuted" after a refit that moved
    the constant, and would never have said by how much -- which is the difference between a
    bound that is slightly optimistic and one that is wrong by 63x.

    Args:
        measured_mg_per_gdcw: The measured content, in the unit the literature reports.
        analyte: What was measured, and **not** simply the pathway's product name. A
            total-carotenoid absorbance sum and an HPLC beta-carotene peak are different
            measurements, and only the second can refute a beta-carotene ceiling: Lopez 2019
            reports 21 mg/gDCW, 17x this package's bound, and refutes nothing because a plate
            reader sums every carotenoid in the extract. Declaring the analyte is what stops
            the larger, weaker number being reached for.
        assay: How it was measured, in enough words to tell a specific assay from a lumped
            one.
        source: The paper, so the number and its provenance travel together.
    """

    measured_mg_per_gdcw: float
    analyte: str
    assay: str
    source: str

    def __post_init__(self) -> None:
        if not math.isfinite(self.measured_mg_per_gdcw) or self.measured_mg_per_gdcw <= 0.0:
            raise ValueError(
                f"measured_mg_per_gdcw must be finite and positive, got "
                f"{self.measured_mg_per_gdcw}")
        for label, text in (("analyte", self.analyte), ("assay", self.assay),
                            ("source", self.source)):
            if not text.strip():
                raise ValueError(
                    f"{label} must say something: a refutation with no {label} cannot be "
                    f"checked against the measurement it claims to be")


def _rebuild_content_ceiling(value, bounded_node, capacity_node, capacity_enzyme,
                             molar_mass_g_per_mol, refuted_by):
    """Reconstructor for :meth:`ContentCeiling.__reduce__`; see the note there."""
    return ContentCeiling(value, bounded_node=bounded_node, capacity_node=capacity_node,
                          capacity_enzyme=capacity_enzyme,
                          molar_mass_g_per_mol=molar_mass_g_per_mol, refuted_by=refuted_by)


class ContentCeiling(float):
    """A ceiling on terminal content, mmol/gDCW, carrying what it is a ceiling OF.

    A ``float`` SUBCLASS, and the choice is the whole point of the class. Everything that
    made this number safe to read -- that it bounds one fitted calibration rather than the
    organism, that the enzyme setting it is not the enzyme the flux law reads, that Arhar
    2024 measured 63x past it -- used to live in docstrings, in
    :mod:`ystwin.pathway.capacity`, and in one ``predict.py`` note that fires only within 10%
    of the bound. None of it travelled with the value, so a caller holding the number held no
    scope at all. Subclassing ``float`` attaches the scope to the number itself: every
    existing comparison, multiplication and ``is None`` test is bit-for-bit what it was, and
    the scope cannot be dropped by a caller who simply did not know to ask for it.

    A field on :class:`NodeKinetics` would not have done. ``predict.py`` validates that
    dataclass by iterating :func:`dataclasses.fields` and REFUSES any populated field it does
    not recognise as used by the node's rate law and fate, so a new field would have made the
    shipped beta-carotene calibration unloadable in both modes.

    Args:
        bounded_node: The terminal node whose content this bounds.
        capacity_node: The node whose saturating step supplies ``vmax_per_growth``. Its
            content is the cyclase's substrate, not the bounded pool.
        capacity_enzyme: The gene consuming ``capacity_node``, read off the spec. Compare it
            with ``spec.entry_enzyme``: where they differ, raising the model's genotype input
            cannot move this bound, which is the beta-carotene diagnosis in one comparison.
        molar_mass_g_per_mol: Of the bounded node, so :attr:`mg_per_gdcw` needs no lookup.
        refuted_by: The measurement that exceeds this bound, or ``None`` where none is
            declared. **A declaration about one fit, never a default** -- a pathway nobody
            has refuted must not inherit beta-carotene's refutation.
    """

    __slots__ = ("bounded_node", "capacity_node", "capacity_enzyme",
                 "molar_mass_g_per_mol", "refuted_by")

    def __new__(cls, value: float, *, bounded_node: str | None = None,
                capacity_node: str | None = None, capacity_enzyme: str | None = None,
                molar_mass_g_per_mol: float | None = None,
                refuted_by: CeilingRefutation | None = None) -> "ContentCeiling":
        self = super().__new__(cls, value)
        self.bounded_node = bounded_node
        self.capacity_node = capacity_node
        self.capacity_enzyme = capacity_enzyme
        self.molar_mass_g_per_mol = molar_mass_g_per_mol
        self.refuted_by = refuted_by
        return self

    def __reduce__(self):
        # A float subclass pickles through `float.__new__`, which skips `__new__` above and
        # leaves every slot unset. Named here so a deepcopy cannot silently strip the scope.
        return (_rebuild_content_ceiling,
                (float(self), self.bounded_node, self.capacity_node, self.capacity_enzyme,
                 self.molar_mass_g_per_mol, self.refuted_by))

    @property
    def mg_per_gdcw(self) -> float | None:
        """The bound in the unit the literature reports, where a molar mass was declared."""
        if self.molar_mass_g_per_mol is None:
            return None
        return float(self) * self.molar_mass_g_per_mol

    @property
    def legacy_only(self) -> bool:
        """``True`` when a content solved against this calibration is NOT a forecast.

        Which is exactly when a measurement has already passed the bound. A calibration whose
        hard asymptote is known to be wrong still has diagnostic uses -- it is what the fit
        asserts, and a caller must be able to ask what the model thinks in order to be told
        it is wrong -- but the number it produces is a statement about the fit and not about
        the strain, in whichever ``predict_product`` mode produced it.
        """
        return self.refuted_by is not None

    @property
    def refuted_fold(self) -> float | None:
        """How many fold the refuting measurement passes this bound, or ``None``.

        Computed rather than stored. The 63x in this package's prose is a number about two
        other numbers, and the last set of claims of that shape went stale three days after
        they were written -- see :mod:`ystwin.pathway.published_cassettes`.
        """
        ceiling = self.mg_per_gdcw
        if self.refuted_by is None or ceiling is None or ceiling <= 0.0:
            return None
        return self.refuted_by.measured_mg_per_gdcw / ceiling

    def describe(self) -> str:
        """One line saying what this is a ceiling of, for a note or a table cell."""
        where = (f"{self.bounded_node} content" if self.bounded_node else "terminal content")
        how = (f", set by {self.capacity_node}'s {self.capacity_enzyme} capacity"
               if self.capacity_node and self.capacity_enzyme else "")
        unit = (f"{self.mg_per_gdcw:.4g} mg/gDCW" if self.mg_per_gdcw is not None
                else f"{float(self):.4g} mmol/gDCW")
        if self.refuted_by is None:
            return f"ceiling on {where}{how}: {unit}"
        return (f"ceiling on {where}{how}: {unit} -- REFUTED {self.refuted_fold:.0f}x by "
                f"{self.refuted_by.measured_mg_per_gdcw:g} mg/gDCW of "
                f"{self.refuted_by.analyte} ({self.refuted_by.assay}, "
                f"{self.refuted_by.source}), so this calibration is legacy-only and a "
                f"content solved against it is a property of the fit, not of the strain")


def content_ceiling(spec: PathwaySpec,
                    kinetics: dict[str, NodeKinetics]) -> ContentCeiling | None:
    """The highest terminal content this calibration can EVER return, mmol/gDCW.

    A hard asymptote, and one nobody noticed the model was asserting until an adversarial
    pass went looking. For a terminal node fed by a saturating step:

        content = flux_out / mu    and    flux_out <= vmax = vmax_per_growth * mu

    so ``content <= vmax_per_growth`` and **mu cancels exactly**. No genotype at any growth
    rate can exceed it. For beta-carotene that is 2.3252e-3 mmol/gDCW = 1.25 mg/gDCW.

    This was the sharpest falsifiable claim in the package, and **it has been falsified.**
    It is a claim about every strain that will ever be built rather than about the six that
    were measured, and it is robust to the fit -- leave-one-out moves it between 1.11 and
    1.41 mg/gDCW -- so it is not resting on a single state. It is resting on the wrong
    enzyme.

    **Arhar 2024 (PMID 39215465) measured 79 mg/gDCW** of beta-carotene by HPLC on
    gravimetric dry weight with lycopene below detection, which is a specific assay and not
    a total-carotenoid absorbance sum. The bound is low by about **63x**. See
    :mod:`ystwin.pathway.capacity`, which carries the diagnosis: the ceiling is set by
    ``vmax_per_growth`` at the crtYB cyclase, the flux law is fitted on crtE, and raising
    the one genotype input this package accepts buys 0.06% across a hundredfold range. What
    is missing is one coefficient relating cyclase dosage to capacity, and it is missing
    because the five published cassettes that state a crtYB dosage numerically hold only two
    distinct values between them, confounded with host, promoter, medium and assay. This
    sentence read "all 28 published cassettes state their dosage qualitatively" until
    2026-09-04, which `data/carotenoid_batch/published_batch_titres.tsv` refutes; see
    :mod:`ystwin.pathway.published_cassettes`, which derives the reason from the table.

    This paragraph read "It is also 17x below Lopez 2019's 21 mg/gDCW. That is a red flag
    rather than a refutation: different strain, flask rather than chemostat, and TOTAL
    carotenoid rather than beta-carotene. Somebody should close that gap." Somebody did, on
    2026-08-30, and the answer arrived in `capacity.py` and `docs/HARD_TESTS.md` without
    reaching the function that computes the bound -- so the two files disagreed about
    whether the package's central claim was alive for a day, and this one said it was.
    Corrected by the 2026-08-31 audit.

    The function still returns the bound, and should: it is what the fitted calibration
    asserts, `predict_product` annotates any prediction within 10% of it, and a caller
    needs to be able to ask what the model thinks in order to be told it is wrong. What it
    is not is a property of the organism.

    **It no longer returns that bound as a bare number.** It returns a
    :class:`ContentCeiling`, which is a ``float`` -- every existing caller's arithmetic is
    unchanged -- carrying the node it bounds, the node and enzyme whose capacity sets it, and
    the measurement that refutes it where a calibration declares one. Until this the whole
    scope lived in prose beside the value, and prose does not travel with a float: the one
    place ``predict.py`` says any of it is a note that fires only within 10% of the bound, so
    an ordinary prediction reported the refuted capacity with nothing attached.

    The enzyme is read off the spec rather than looked up per product, which is what makes
    the diagnosis checkable instead of asserted: compare ``ceiling.capacity_enzyme`` with
    ``spec.entry_enzyme`` and beta-carotene's structural failure -- the ceiling is set at
    CrtYB, the only genotype input is CrtE -- is one comparison rather than a paragraph.

    Returns ``None`` where no ceiling exists -- a chain whose terminal node is fed by a
    passthrough or proportional step has none, because nothing saturates.
    """
    if len(spec.nodes) < 2:
        return None
    if spec.nodes[-1].fate != Fate.DILUTED:
        # The bound above rests on `mu` cancelling, and it only cancels when dilution is the
        # terminal node's ONLY outlet. With a third outlet the balance is
        # `X = v_in/(mu + k_loss)`, so `X <= vmax_per_growth * mu/(mu + k_loss)` -- strictly
        # increasing in `mu`, reaching `vmax_per_growth` only in the limit. The old return
        # would still have been a valid upper bound and would no longer have been the
        # asymptote the docstring calls it, which is the kind of half-true a caller acts on.
        # This signature has no growth rate and structurally cannot express the tight bound,
        # so it declines instead. Shipped DILUTED specs are untouched.
        return None
    feeder = spec.nodes[-2]
    if feeder.rate_law != RateLaw.SATURATING:
        return None
    parameters = kinetics.get(feeder.name)
    if parameters is None or parameters.vmax_per_growth is None:
        return None
    terminal = spec.nodes[-1]
    # `getattr` and not an isinstance branch: a calibration declares its refutation by
    # shipping a marked capacity, and an unmarked plain float must stay unmarked here.
    return ContentCeiling(
        parameters.vmax_per_growth,
        bounded_node=terminal.name,
        capacity_node=feeder.name,
        capacity_enzyme=feeder.enzyme or None,
        molar_mass_g_per_mol=terminal.molar_mass_g_per_mol,
        refuted_by=getattr(parameters.vmax_per_growth, "refuted_by", None))


@dataclass(frozen=True)
class PathwaySolution:
    """Every node of one pathway at steady state."""

    product: str
    growth_rate_per_h: float
    entry_flux: float
    nodes: tuple[NodeState, ...]

    @property
    def terminal(self) -> NodeState:
        return self.nodes[-1]

    def node(self, name: str) -> NodeState:
        for candidate in self.nodes:
            if candidate.name == name:
                return candidate
        raise KeyError(f"no node {name!r}; have {[n.name for n in self.nodes]}")

    @property
    def carbon_closes(self) -> float:
        """Entry flux minus everything leaving the chain, which must be zero.

        At steady state every mole entering leaves as somebody's dilution term OR through a
        declared outlet of its own, so this is an identity rather than a fitted check -- and
        identities are exactly what a solver should be asked to reproduce, because an
        arithmetic slip breaks them.

        The second half of that sentence is new, and it was missing for a day. This summed
        only ``mu*[X]``, which is every outlet a chain of ``DILUTED`` nodes has and not every
        outlet the solver can now model. From the moment ``DEGRADED`` began solving, a
        storage chain reported its own degradation flux here as a violation -- 0.4 against an
        entry flux of 0.5 -- while the docstring called the number an identity. The export
        outlet would have done the same thing in the same position, twice as large.

        Where every node is ``DILUTED`` this is bit-for-bit what it always returned:
        ``loss_flux_mmol_per_gdcw_h`` is a literal ``0.0`` on each one.
        """
        washed_out = sum(self.growth_rate_per_h * n.content_mmol_per_gdcw
                         for n in self.nodes)
        left_the_chain = sum(n.loss_flux_mmol_per_gdcw_h for n in self.nodes)
        return self.entry_flux - washed_out - left_the_chain


def solve_pathway(
    spec: PathwaySpec,
    entry_flux: float,
    growth_rate: float,
    kinetics: dict[str, NodeKinetics] | None = None,
) -> PathwaySolution:
    """Solve every node of ``spec`` at steady state.

    Args:
        spec: The declared pathway.
        entry_flux: Flux into the first node, mmol/gDCW/h.
        growth_rate: Specific growth rate, /h. Strictly positive: a non-growing culture has
            no washout, so every pool diverges and there is no steady state to report.
        kinetics: Parameters per node name, for nodes whose rate law needs them.

    Raises:
        ValueError: on a non-positive growth rate, a negative flux, or a node whose rate
            law needs parameters that were not supplied.
    """
    if entry_flux < 0:
        raise ValueError(f"entry_flux must be non-negative, got {entry_flux}")
    if growth_rate <= 0:
        raise ValueError(
            f"growth_rate must be positive, got {growth_rate}. Every pool here is held down "
            "by washout, so at zero growth there is no steady state to solve for")
    kinetics = kinetics or {}

    # EVERY node's declared range, checked before anything is solved.
    #
    # This used to live inside the SATURATING and PROPORTIONAL arms of the walk below, which
    # made the guard an accident of a node's rate law rather than a property of the solver:
    # the terminal and passthrough arms never reached it, so `phb` and `glycogen` -- both
    # passthrough end to end -- had no growth-rate guard at all. Checking here also means a
    # refusal happens before any arithmetic, so a spec cannot half-solve and then raise.
    for node in spec.nodes:
        parameters = kinetics.get(node.name)
        if parameters is not None:
            parameters.check_growth_rate(node.name, growth_rate)

    states: list[NodeState] = []
    flux_in = float(entry_flux)
    for index, node in enumerate(spec.nodes):
        terminal = index == len(spec.nodes) - 1
        if terminal:
            # Growth dilution, plus a native enzymatic outlet where the spec declares one.
            #
            # `DEGRADED` used to be refused by `PathwaySpec` at load, which made every
            # storage compound unrepresentable rather than under-parameterised. The refusal
            # was right about the arithmetic -- `content = flux/mu` omits the term doing
            # most of the work -- and wrong about where to put it: the missing thing is a
            # rate constant, and a missing constant is what `kinetics` is for. So the
            # refusal moved here, where it can tell the difference between "this pathway is
            # out of scope" and "this pathway needs one more measured number".
            parameters = kinetics.get(node.name)
            loss = growth_rate + _degradation_loss(node, parameters)
            export = _secretion_capacity(node, parameters)
            if export is None:
                content = flux_in / loss
                loss_flux = (loss - growth_rate) * content
            else:
                # The export outlet, and it is the arm above under a shifted coefficient.
                #
                #     v_in - loss*X - vmax*X/(km + X) = 0
                #     =>  loss*X^2 + (loss*km + vmax - v_in)*X - v_in*km = 0
                #
                # which is TERM FOR TERM the SATURATING quadratic below with `mu` replaced
                # by the whole linear loss, so a third outlet adds no new algebraic form to
                # this module -- it shifts the coefficient the existing form already
                # carries. `a > 0` and `c <= 0`, so the roots have opposite signs and there
                # is exactly one non-negative root; the discriminant is `b*b + 4*a*v_in*km`
                # and cannot go negative, which is the same argument the SATURATING arm
                # rests on without stating it.
                vmax_sec, km_sec = export
                a = loss
                b = loss * km_sec + vmax_sec - flux_in
                c = -flux_in * km_sec
                discriminant = math.sqrt(b * b - 4.0 * a * c)
                # The STABLE branch, which the SATURATING arm below does not use.
                #
                # `(-b + sqrt(D))/(2a)` subtracts nearly equal numbers whenever b > 0, and
                # b > 0 is not an edge case for this outlet -- it is `vmax + loss*km > v_in`,
                # the UNSATURATED regime, where a strain sits below its export capacity.
                # That is the common case and the one a caller most wants. At km = 1e-16 the
                # naive form is wrong by 4% and at 1e-18 it returns exactly zero.
                #
                # The SATURATING arm below is deliberately left alone: it has the same
                # weakness, but at the shipped calibration (km = 5.974e-4) the two branches
                # differ by at most 8.2e-16 relative, and rewriting it would move committed
                # numbers in their last bits to fix nothing that is biting.
                if b >= 0.0:
                    content = (-2.0 * c) / (b + discriminant)
                else:
                    content = (-b + discriminant) / (2.0 * a)
                loss_flux = ((loss - growth_rate) * content
                             + vmax_sec * content / (km_sec + content))
            flux_out = 0.0
        elif node.rate_law == RateLaw.PASSTHROUGH:
            content = 0.0
            flux_out = flux_in
            loss_flux = 0.0
        elif node.rate_law == RateLaw.PROPORTIONAL:
            k = kinetics.get(node.name, NodeKinetics()).require_proportional(node.name)
            content = flux_in / (k + growth_rate)
            flux_out = k * content
            loss_flux = 0.0
        elif node.rate_law == RateLaw.SATURATING:
            vmax_per_growth, km = kinetics.get(
                node.name, NodeKinetics()).require_saturating(node.name)
            vmax = vmax_per_growth * growth_rate
            a = growth_rate
            b = growth_rate * km + vmax - flux_in
            c = -flux_in * km
            content = (-b + math.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)
            flux_out = vmax * content / (km + content)
            loss_flux = 0.0
        else:  # unreachable: RateLaw validates on construction
            raise ValueError(f"node {node.name!r}: unknown rate law {node.rate_law!r}")

        if node.molar_mass_g_per_mol is not None:
            grams = content * node.molar_mass_g_per_mol / 1000.0
            if grams > DRY_WEIGHT_G_PER_GDCW:
                raise ImplausibleContent(
                    f"node {node.name!r} solves to {content:.4g} mmol/gDCW, which at "
                    f"{node.molar_mass_g_per_mol:g} g/mol is {grams:.4g} g per gram of dry "
                    f"cell weight -- {grams * 100:.0f}% of the cell's own mass. The "
                    f"arithmetic is right and the inputs are not: entry flux "
                    f"{entry_flux:.4g} mmol/gDCW/h at mu {growth_rate:.4g} /h. Check the "
                    "entry expression against the range the flux law was fitted on")

        states.append(NodeState(node.name, content, flux_in, flux_out,
                                node.molar_mass_g_per_mol, loss_flux))
        flux_in = flux_out

    return PathwaySolution(spec.product, float(growth_rate), float(entry_flux),
                           tuple(states))
