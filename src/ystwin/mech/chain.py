"""Environment in, mechanistic states in the middle, the existing product and reporter out.

This is the join. Phase 0 built five pieces that each work and none of which was connected
to anything: `mech/ph.py` computes a cytosolic pH nobody read, `mech/burden.py` a growth tax
nobody applied, `mech/population.py` a bearing fraction nobody multiplied by, `mech/state.py`
and `mech/integrate.py` a reduction audit nobody ran on an assembled system. This module runs
them in order and hands the result to code that already existed -- `pathway/solve.py` for the
product and `generator/culture.py` for the reporter -- rather than writing a second copy of
either.

THE CHAIN, top to bottom, with the layer name each step reports under::

    context (carbon, temperature, oxygen, phase)  ->  mu_env      environment   IN_CHAIN
    medium pH + total acid  ->  [AH]_o                             weak acid     IN_CHAIN
    [AH]_o  ->  (A_i, pH_c) integrated over the window             pH states     IN_CHAIN
    [A-]_i  ->  Pma1 ATP bill  ->  mu_max                          proton bill   IN_CHAIN
    declared copies  ->  mu_max                                    burden        IN_CHAIN
    mu = min(setpoint, mu_max)                                     vessel        IN_CHAIN
    mu  ->  (F, generations)                                       population    IN_CHAIN
    mu  ->  solve_pathway  ->  content; x F  ->  titre per cell    metabolism    IN_CHAIN
    mu, dose  ->  simulate_culture  ->  reporter                   reporter      IN_CHAIN
    pH_c(t)  ->  Citrine quench  ->  observed activity             photophysics  IN_CHAIN
    [AH]_o  ->  Yeast9 pfba delta                                  GEM depth     AUDITS

**The one structural change to the incumbent, and it is a correction rather than an
addition.** `predict.py::_growth_rate` checks whether the strain can REACH a fed-batch
setpoint only when a stressor is present. With no stressor the setpoint is returned
unchecked, so an ethanol fed-batch at mu_set = 0.18 /h is accepted although
`context_growth_rate` puts that culture's maximum at 0.14 /h. That is why seven
environments returned one content to eight decimal places: not because the environment has
no route, but because the route it has was tested on one branch and not the other. This
module applies ``mu = min(mu_set, mu_max)`` on both branches, which is the same arithmetic
`fba/fedbatch.py` already enforces at design time. The proposed one-line diff to
`predict.py` is in `docs/MECHANISTIC_LAYER.md`; this file does not make it, because
`predict.py` is the repository's entry point.

**THE ASSEMBLED FREE-SCALAR GATE FAILS, AND THAT IS THE HEADLINE NEGATIVE.** Each block
passes criterion (e) on its own -- pH 2 free against 3 targets, population 2 against 2,
burden 1 against 1. Assembled, the free scalars add and the targets do not: `growth_rate`
is one measurement counted by three blocks and `reporter_activity` is one counted by two,
so the union is **7 free scalars against 3 independent targets** -- the blocks' five plus
the two this assembly asserts, against three distinct assays. :func:`assembled_gate`
computes it and :func:`require_assembled_gate` raises. The consequence is enforced in the
API rather than written in a comment: :func:`run_chain` takes a :class:`SweepPoint` with no
default, so no caller can obtain a single number without saying which corner of the swept
space produced it, and :func:`chain_band` is the shape a result should be quoted in.

**What this module does NOT do.** It does not fit anything. It adds no constant of its own
that is not either read from another module or tagged in :data:`CHAIN_PARAMS`. It does not
put the allocation arm of `mech/burden.py` or the ethanol diagnostic of
`mech/population.py` in the chain -- both fail their own gates and both report. And it does
not claim criterion (a) on the product channel: content is measured by HPLC and this
repository has no measured floor for that assay, so the content movements below are stated
as effect sizes and scored only on the two observables whose floors are measured, growth
rate at 0.0117 /h and reporter activity at a CV of 0.146.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from ..bridge.maintenance_calibration import REFERENCE_REGIME, MaintenanceRegime
from ..generator.context import context_growth_rate
from ..generator.culture import CultureParameters, simulate_culture
from ..generator.plate import DEFAULT_PANEL, PlateConditions
from ..generator.stress_panel import (
    MODULES,
    STRESSORS,
    growth_rate as panel_growth_rate,
    module_response,
    reporter_loadings,
    transcriptional_reporters,
)
from ..predict import Environment, Genotype, LayerState, LayerStatus
from . import burden as mech_burden
from . import ph as mech_ph
from . import population as mech_population
from .ablation import GROWTH_RATE_FLOOR, REPORTER_ACTIVITY_FLOOR
from .integrate import ReducedSystem, Trajectory, reduce_for_window
from .params import FreeScalarGate, Param, ParamRegistry
from .state import PLATE_READ_4H, MechState, Window

__all__ = [
    "CHAIN_PARAMS",
    "ChainResult",
    "GeometryResult",
    "LAYER_ORDER",
    "ProductArm",
    "ReporterArm",
    "SweepPoint",
    "assembled_gate",
    "assembled_states",
    "chain_band",
    "dose_time_surface",
    "geometry_audit",
    "ode_state_census",
    "provenance",
    "ray_share",
    "reachable_reactions",
    "reduction_audit",
    "require_assembled_gate",
    "run_chain",
]


LAYER_ORDER = (
    "environment", "stress panel", "weak acid", "pH states", "proton bill", "burden", "vessel",
    "population", "metabolism", "reporter", "photophysics", "GEM depth",
    "allocation", "ethanol",
)
"""Every layer this module reports, in the order the chain runs them.

Fixed so a caller can diff two runs row by row, and so a missing row is a bug rather than a
layer that quietly did not run. `predict.py` builds its own list in call order for the same
reason; this one is declared because the mechanism has optional arms.
"""


# --------------------------------------------------------------------------------------
# The chain's own constants. There are two, and both are structural rather than physical.
# --------------------------------------------------------------------------------------

CHAIN_PARAMS = ParamRegistry("mech/chain.py::assembly")
"""What the ASSEMBLY adds. Everything physical is read from the block that owns it."""

INITIAL_BEARING_FRACTION = CHAIN_PARAMS.add(Param.asserted(
    "F0_bearing", 0.95, "fraction of cells bearing the plasmid",
    "ASSERTED at 0.95, the same value mech/population.py::simulate_population declares as "
    "its own default and for the same reason: a transformant preculture has already spent "
    "generations before the run starts, and Hohnholz 2017's industrial framing (the paper "
    "mech/population.py cites for p) counts 50-60 of them from the primary seed lot. "
    "Nothing here measures it; it is carried so the two modules cannot drift apart"))

REFERENCE_COPIES = CHAIN_PARAMS.add(Param.asserted(
    "n_copies_reference", 1.0, "cassette copies",
    "ASSERTED at 1. mech/burden.py::burdened_growth_rate predicts only the INCREMENT over "
    "the genotype the environment's growth rate was measured on, and the panel strains' "
    "declared genotype is one copy. It is a declaration and not a measurement of them"))

CHAIN_PARAMS.add_target(GROWTH_RATE_FLOOR)
CHAIN_PARAMS.add_target(REPORTER_ACTIVITY_FLOOR)


IN_CHAIN_REGISTRIES = (
    mech_ph.PH_PARAMS,
    mech_population.POPULATION_PARAMS,
    mech_burden.GROWTH_PARAMS,
    CHAIN_PARAMS,
)
"""The registries whose parameters actually reach a returned number.

`mech/burden.py::ALLOCATION_PARAMS`, `mech/population.py::ETHANOL_PARAMS` and both
``NOT_BUILT`` registries are deliberately absent: they are REPORTED, they feed nothing, and
counting their free scalars against the chain's targets would make the gate answer a
question about layers that cannot move the number.
"""


def assembled_gate(registries: Sequence[ParamRegistry] = IN_CHAIN_REGISTRIES
                   ) -> FreeScalarGate:
    """Criterion (e) across the whole assembled layer, with targets DE-DUPLICATED by name.

    The de-duplication is the entire content of this function. Each block registers the
    floors it is scored against, and `mech/ablation.py`'s ``GROWTH_RATE_FLOOR`` is one
    object registered by three of them. Summing the per-block gates counts one growth assay
    three times and reports 5 free against 6 targets, which passes; counting the assay once
    reports **5 against 3**, which does not.

    Returns:
        A :class:`~ystwin.mech.params.FreeScalarGate` naming every free scalar and every
        distinct target. It is expected to FAIL, and :func:`require_assembled_gate` is what
        turns that into a refusal at a call site.
    """
    free: list[str] = []
    targets: dict[str, None] = {}
    fitted: dict[str, None] = {}
    for registry in registries:
        for param in registry.free_scalars():
            free.append(f"{registry.piece}::{param.name}")
        for target in registry.independent_targets():
            targets[target.name] = None
        for target in registry.targets:
            if target.fitted:
                fitted[target.name] = None
    return FreeScalarGate(
        piece="mech/chain.py::assembled",
        free=tuple(free),
        targets=tuple(targets),
        fitted_targets=tuple(fitted),
    )


def require_assembled_gate(registries: Sequence[ParamRegistry] = IN_CHAIN_REGISTRIES
                           ) -> FreeScalarGate:
    """Raise unless the assembled layer has at least as many targets as free scalars.

    Nothing in this module calls it on the normal path, because the chain does not build a
    single fitted point: :func:`run_chain` demands a declared :class:`SweepPoint` and
    :func:`chain_band` returns the corners. This exists so that a caller who wants to FIT
    the assembled layer is refused at the moment they ask, with the count in the message.

    Raises:
        FreeScalarGateFailed: with the free scalars and the distinct targets named.
    """
    from .params import FreeScalarGateFailed

    result = assembled_gate(registries)
    if not result.passes:
        raise FreeScalarGateFailed(
            f"{result.summary()}. Free: {', '.join(result.free)}. Distinct targets: "
            f"{', '.join(result.targets)}. Each block passes its own gate; the union does "
            "not, because the growth assay and the plate activity CV are each ONE "
            "measurement that several blocks count as their own. Run the chain over "
            "mech.chain.chain_band and quote a band, or measure one of the free scalars")
    return result


# --------------------------------------------------------------------------------------
# The declared sweep point. No defaults: the gate above is why.
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class SweepPoint:
    """One corner of the declared sweep space, stated by the caller.

    Every field is a SWEPT parameter of a block below, and the class exists so that a number
    coming out of :func:`run_chain` cannot be quoted without the corner that produced it.
    There is no default constructor: :meth:`corners` is how a caller gets the set.

    Args:
        beta_mM_per_ph: Cytosolic buffering capacity, from
            `mech/ph.py::CYTOSOLIC_BUFFER_CAPACITY`.
        cytosolic_volume_ml_per_gdcw: Cell-water volume, from
            `mech/ph.py::CYTOSOLIC_VOLUME_ML_PER_GDCW`.
        p_per_division: Segregational loss per division, from
            `mech/population.py::LOSS_PER_DIVISION`.
        burden: Growth penalty of the bearing population, from
            `mech/population.py::BURDEN`.
    """

    beta_mM_per_ph: float
    cytosolic_volume_ml_per_gdcw: float
    p_per_division: float
    burden: float

    def __post_init__(self) -> None:
        for value, param in (
                (self.beta_mM_per_ph, mech_ph.CYTOSOLIC_BUFFER_CAPACITY),
                (self.cytosolic_volume_ml_per_gdcw, mech_ph.CYTOSOLIC_VOLUME_ML_PER_GDCW),
                (self.p_per_division, mech_population.LOSS_PER_DIVISION),
                (self.burden, mech_population.BURDEN)):
            low, high = param.bounds
            if not low <= float(value) <= high:
                raise ValueError(
                    f"{param.name} is SWEPT over [{low:g}, {high:g}] and {value:g} is "
                    f"outside it. A point off the declared axis is an ASSERTED value with a "
                    "sweep's provenance attached")

    @classmethod
    def corners(cls) -> tuple["SweepPoint", ...]:
        """The 16 corners of the four declared axes, in a fixed order.

        Corners rather than a grid: every quantity the chain computes is monotone in each
        axis separately -- the ATP bill rises with volume and falls with buffering, the
        bearing fraction falls with loss and with burden -- so the extreme values of any
        output are attained at a corner and an interior grid adds cost without adding range.
        """
        axes = [mech_ph.CYTOSOLIC_BUFFER_CAPACITY.bounds,
                mech_ph.CYTOSOLIC_VOLUME_ML_PER_GDCW.bounds,
                mech_population.LOSS_PER_DIVISION.bounds,
                mech_population.BURDEN.bounds]
        out = []
        for index in range(16):
            picks = [axis[(index >> shift) & 1] for shift, axis in enumerate(axes)]
            out.append(cls(*picks))
        return tuple(out)

    @classmethod
    def midpoint(cls) -> "SweepPoint":
        """The centre of the declared axes. For a smoke test, never for a quoted number."""
        return cls(*[0.5 * (low + high) for low, high in (
            mech_ph.CYTOSOLIC_BUFFER_CAPACITY.bounds,
            mech_ph.CYTOSOLIC_VOLUME_ML_PER_GDCW.bounds,
            mech_population.LOSS_PER_DIVISION.bounds,
            mech_population.BURDEN.bounds)])

    def label(self) -> str:
        return (f"beta={self.beta_mM_per_ph:g} V={self.cytosolic_volume_ml_per_gdcw:g} "
                f"p={self.p_per_division:.4g} b={self.burden:.4g}")


# --------------------------------------------------------------------------------------
# The assembled state vector and its tau/T audit
# --------------------------------------------------------------------------------------

def assembled_states(window: Window, *, ah_out_mM: float, sweep: SweepPoint,
                     acid_name: str = "acetic") -> MechState:
    """Every mechanistic state the chain integrates, in chain order.

    Five: ``A_i`` and ``pH_c`` from `mech/ph.py`, ``product_fraction`` from
    `mech/burden.py`, ``plasmid_bearing`` and ``generations`` from `mech/population.py`.
    The biomass and reporter states of `generator/culture.py` are NOT here -- they belong to
    the incumbent generator, they are integrated by it, and :func:`ode_state_census` counts
    them separately so the two claims stay apart.

    Args:
        window: The declared window; supplies the measured growth band every time constant
            is computed against.
        ah_out_mM: Undissociated acid outside. The pH time constants depend on the dose.
        sweep: The declared corner; ``p_per_division`` and ``burden`` set the bearing
            fraction's time constant.
        acid_name: Which acid, for ``k_entry``.
    """
    p = max(float(sweep.p_per_division), 1e-9)
    return MechState(
        mech_ph.ph_states(window, ah_out_mM=max(float(ah_out_mM), 1e-9),
                          acid_name=acid_name).variables
        + mech_burden.BURDEN_STATES.variables
        + mech_population.population_states(window, p_bounds=(p, p),
                                            burden=sweep.burden).variables)


def reduction_audit(window: Window, *, ah_out_mM: float, sweep: SweepPoint,
                    acid_name: str = "acetic") -> ReducedSystem:
    """Criterion (d) for the ASSEMBLED state, in one window. The audit, not a claim.

    The result differs between the two windows and that difference is the reason `Window` is
    an argument here rather than a runtime detail: nothing that survives on the 4.14 h plate
    read survives on the 5 d fed-batch, and the state that dominates the fed-batch is frozen
    on the plate. :meth:`ReducedSystem.audit_table` prints both.
    """
    return reduce_for_window(
        assembled_states(window, ah_out_mM=ah_out_mM, sweep=sweep, acid_name=acid_name),
        window)


def ode_state_census(window: Window, *, ah_out_mM: float, sweep: SweepPoint,
                     with_maturation: bool = True) -> dict[str, int]:
    """How many differential equations the chain actually integrates, by owner.

    Written as a census rather than a constant because the answer depends on the window and
    on whether the reporter carries a maturation state, and because the number this replaces
    -- "one or two states in the old stress path" -- is the claim the whole build is scored
    against.

    Returns:
        ``{"mechanism", "incumbent generator", "total", "kept_after_reduction"}``.
    """
    states = assembled_states(window, ah_out_mM=ah_out_mM, sweep=sweep)
    reduced = reduce_for_window(states, window)
    incumbent = 1 + (2 if with_maturation else 1)
    return {
        "mechanism": len(states),
        "incumbent generator": incumbent,
        "total": len(states) + incumbent,
        "kept_after_reduction": len(reduced.integrated),
    }


# --------------------------------------------------------------------------------------
# The result
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class ProductArm:
    """What the metabolism layer returned, at the mechanistic growth rate.

    Args:
        solution: The `pathway/solve.py` solution. The existing solver, unmodified.
        content_mg_per_gdcw: Content of the terminal node, per gram of BEARING cells.
        population_content_mg_per_gdcw: The same number times the time-averaged bearing
            fraction, which is what an assay on the whole culture would read.
        mean_bearing_fraction: The average of ``F(t)`` over the window.
        flux_mmol_per_gdcw_h: Entry flux the solve ran at.
    """

    solution: object
    content_mg_per_gdcw: float
    population_content_mg_per_gdcw: float
    mean_bearing_fraction: float
    flux_mmol_per_gdcw_h: float


@dataclass(frozen=True)
class ReporterArm:
    """What the reporter and photophysics layers returned.

    Args:
        times_h: The read grid.
        activity: Per-cell mature reporter from `generator/culture.py::simulate_culture`,
            before the quench. The incumbent's own output.
        observed: ``activity * citrine_quench(pH_c(t))`` -- what the instrument sees.
        quench: The multiplier, on the same grid.
        nutrient_factor: The fraction of its unstressed growth the mechanism left the
            culture, which is the knob the mechanism enters `simulate_culture` through.
        construct: Which panel construct was simulated.
    """

    times_h: np.ndarray
    activity: np.ndarray
    observed: np.ndarray
    quench: np.ndarray
    nutrient_factor: float
    construct: str

    @property
    def endpoint_fold(self) -> float:
        """Observed signal at the end of the read over its own start. Dimensionless."""
        return float(self.observed[-1] / self.observed[0])


@dataclass(frozen=True)
class ChainResult:
    """One run of the chain at one declared sweep corner, with every layer accounted for.

    ``layers`` is the field to read first, and it carries the same
    :class:`~ystwin.predict.LayerStatus` objects `predict.py` returns, so a caller can print
    both with one function. The rest is the arithmetic those layers produced.

    Args:
        environment: The condition, as given.
        sweep: The corner. Present because no number here is meaningful without it.
        window: The declared window.
        mu_env_per_h: What the context alone says the culture can do.
        mu_max_per_h: After the acid bill and the declared cassette burden.
        growth_rate_per_h: ``min(setpoint, mu_max)`` -- the rate everything downstream ran
            at.
        feed_sets_mu: False when the setpoint is above what the strain can reach, which is
            the state `predict.py` raises on and this module reports instead.
        ah_out_mM: Undissociated acid outside, or 0.0 with no acid.
        ph_states: The ``(A_i, pH_c)`` trajectory, or None.
        proton_bill: The Pma1 ATP bill, or None.
        population: The ``(F, generations)`` trajectory.
        product: The metabolism arm, or None when no pathway was supplied.
        reporter: The reporter arm, or None when no construct was supplied.
        layers: One row per entry in :data:`LAYER_ORDER`.
        notes: Anything a reader would otherwise have to derive.
    """

    environment: Environment
    sweep: SweepPoint
    window: Window
    mu_env_per_h: float
    mu_max_per_h: float
    growth_rate_per_h: float
    feed_sets_mu: bool
    ah_out_mM: float
    ph_states: Trajectory | None
    proton_bill: object | None
    population: Trajectory
    product: ProductArm | None
    reporter: ReporterArm | None
    layers: tuple[LayerStatus, ...]
    notes: tuple[str, ...] = ()

    def layer(self, name: str) -> LayerStatus:
        for candidate in self.layers:
            if candidate.name == name:
                return candidate
        raise KeyError(f"no layer {name!r}; have {[layer.name for layer in self.layers]}")

    def layer_report(self) -> str:
        """Every layer, one per line. The same shape as `ProductPrediction.layer_report`."""
        return "\n".join(layer.summary() for layer in self.layers)

    def in_chain(self) -> tuple[str, ...]:
        """The layers the returned numbers actually depend on."""
        return tuple(layer.name for layer in self.layers
                     if layer.state == LayerState.IN_CHAIN)

    @property
    def bill_is_physical(self) -> bool:
        """Whether the pump-holds arm implies an anion pool a cell could actually contain.

        **Read this before quoting a growth effect from an acid.** The bill holds pH_c at
        rest, so ``[A-]_i = [AH]_o * 10^(pH_c - pKa)`` is 210-fold the outside undissociated
        concentration, and it passes `mech/ph.py`'s own 500 mM osmotic flag above about 3.7
        mM total acetate at medium pH 4.5. At the largest dose where it does NOT, the growth
        effect is 0.00395 /h = 0.34x the growth assay's floor. Every larger number this
        layer reports on that arm is arithmetic on an unreachable state, which is why
        `mech/ph.py` flags rather than clips and why this property exists rather than a
        silent cap.
        """
        return True if self.proton_bill is None else bool(self.proton_bill.physical)

    @property
    def ph_c_final(self) -> float:
        """Cytosolic pH at the end of the window, or the resting value with no acid."""
        if self.ph_states is None:
            return float(mech_ph.RESTING_CYTOSOLIC_PH)
        return float(self.ph_states.of("pH_c")[-1])

    def summary(self) -> str:
        product = ("--" if self.product is None
                   else f"{self.product.population_content_mg_per_gdcw:.8f} mg/gDCW")
        return (f"mu {self.mu_env_per_h:.4f} -> {self.mu_max_per_h:.4f} -> "
                f"{self.growth_rate_per_h:.4f} /h | pH_c {self.ph_c_final:.3f} | "
                f"product {product} | {self.sweep.label()}")


# --------------------------------------------------------------------------------------
# The chain
# --------------------------------------------------------------------------------------

def _acid_for(environment: Environment) -> tuple[str | None, float]:
    """The mech/ph acid name and total mM this environment carries, or (None, 0.0).

    Only `stress_panel.STRESSORS`' own weak acid is mapped, and only onto an acid
    `mech/ph.py` has a MEASURED yeast permeability for. Everything else returns None, which
    is what makes the pH layer NOT_RUN rather than silently zero.
    """
    if environment.stressor == "acetic_acid" and environment.dose > 0:
        return "acetic", float(environment.dose)
    return None, 0.0


def _speciation(acid_name: str, total_mM: float, ph_ex: float) -> tuple[float, float]:
    """``([AH]_o, pKa)``. Exact Henderson-Hasselbalch, no free scalar anywhere in it."""
    pka = float(mech_ph.acid(acid_name).pka)
    return mech_ph.undissociated_outside(total_mM, ph_ex, pka), pka


def _panel_retained(environment: Environment) -> float:
    """Fraction of its unstressed growth the incumbent panel leaves at this dose.

    Exactly `predict.py::_growth_rate`'s batch arithmetic, called rather than restated so
    the two cannot drift. 1.0 with no stressor.
    """
    if environment.stressor is None:
        return 1.0
    unstressed = panel_growth_rate(environment.stressor, 0.0)
    if unstressed <= 0:
        return 0.0
    return max(0.0, panel_growth_rate(environment.stressor, environment.dose) / unstressed)


def _mu_env(environment: Environment) -> float:
    """What the context alone says the culture can do, 1/h. The incumbent, unmodified."""
    return context_growth_rate(environment.context)


def run_chain(
    environment: Environment,
    *,
    sweep: SweepPoint,
    window: Window = PLATE_READ_4H,
    genotype: Genotype | None = None,
    n_copies: float = float(REFERENCE_COPIES),
    spec=None,
    flux_calibration=None,
    kinetics: Mapping | None = None,
    construct: str | None = None,
    culture_parameters: CultureParameters | None = None,
    n_points: int = 201,
    initial_fraction: float = float(INITIAL_BEARING_FRACTION),
    regime: MaintenanceRegime = REFERENCE_REGIME,
    gem_model=None,
) -> ChainResult:
    """Run the whole chain for one environment at one declared sweep corner.

    Args:
        environment: The condition, in `predict.py`'s own vocabulary. ``context.ph_medium``
            and a ``stressor="acetic_acid"`` dose are what drive the pH arm;
            ``growth_rate_setpoint_per_h`` is a CEILING here, not an imposition.
        sweep: The corner of the four declared axes. **No default** -- see
            :func:`assembled_gate`.
        window: PLATE_READ_4H or FEDBATCH_5D, or a declared window of the caller's own.
        genotype: Passed to the product arm. Required with ``spec``.
        n_copies: Declared cassette copies for `mech/burden.py`'s growth tax. A DECLARED
            genotype axis: nothing in this repository measures it.
        spec, flux_calibration, kinetics: The product arm. Supply all three or none.
        construct: The reporter arm. A key of `generator/plate.py::DEFAULT_PANEL`, or any
            name if ``culture_parameters`` is supplied.
        culture_parameters: Kinetics for the reporter arm; defaults to the panel's row.
        n_points: Output grid for the integrations. Solver step size is its own business.
        initial_fraction: ``F`` at inoculation.
        regime: Which metabolic regime prices the ATP bill. The default is the reference
            operating point; `bridge/maintenance_calibration.py` owns the two branches and
            the fact that the slope is MODEL-DERIVED rather than measured.
        gem_model: An already-loaded cobra model carrying the acid uncoupling. Given one,
            the GEM depth layer runs and AUDITS; without it that layer is NOT_RUN. A model
            rather than a path, for the same reason `predict.py` takes one.

    Returns:
        A :class:`ChainResult`.

    Raises:
        ValueError: if the product arm is partly supplied, or if the growth rate the
            mechanism lands on is outside what the pathway kinetics were fitted over --
            that refusal comes from `pathway/solve.py` and is left where it is.
    """
    supplied = [spec is not None, flux_calibration is not None, kinetics is not None]
    if any(supplied) and not all(supplied):
        raise ValueError(
            "the product arm needs spec, flux_calibration and kinetics together; a solve "
            "with two of the three would run against a default nobody chose")
    if spec is not None and genotype is None:
        raise ValueError("the product arm needs a genotype to read entry expression from")

    notes: list[str] = []
    layers: list[LayerStatus] = []

    # ---------------------------------------------------------------- environment
    mu_env = _mu_env(environment)
    context = environment.context
    layers.append(LayerStatus(
        "environment", LayerState.IN_CHAIN, True,
        f"{context.carbon_source}, {context.temperature_c:g} C, O2 {context.oxygen:g}, "
        f"{context.growth_phase} -> mu_env = {mu_env:.4f} /h"))

    # ---------------------------------------------------------------- stress panel
    # The incumbent's empirical dose response, kept wherever the mechanism has no entry
    # point. Where it has one it REPLACES this, because both tax the same growth rate.
    panel_retained = _panel_retained(environment)

    # ---------------------------------------------------------------- weak acid
    acid_name, total_mM = _acid_for(environment)
    ah_out = 0.0
    pka = float("nan")
    if acid_name is None:
        detail = ("pass stressor='acetic_acid' with a dose and a medium pH. mech/ph.py "
                  "carries acetic, formic and lactic and REFUSES propionic, benzoic and "
                  "sorbic by name; of those only acetic_acid is on the stress panel")
        if environment.stressor is not None:
            detail = (f"{environment.stressor} is not a weak acid with a measured yeast "
                      f"permeability, so the pH arm has no entry point for it. " + detail)
        layers.append(LayerStatus("weak acid", LayerState.NOT_RUN, False, detail))
    else:
        ah_out, pka = _speciation(acid_name, total_mM, context.ph_medium)
        layers.append(LayerStatus(
            "weak acid", LayerState.IN_CHAIN, True,
            f"{total_mM:g} mM {acid_name} at medium pH {context.ph_medium:g}, pKa {pka:g} "
            f"-> [AH]_o = {ah_out:.4g} mM. Henderson-Hasselbalch, exact, no free scalar"))

    # ---------------------------------------------------------------- pH states
    ph_trajectory: Trajectory | None = None
    bill = None
    mu_after_acid = mu_env
    if acid_name is not None and ah_out > 0.0:
        ph_trajectory = mech_ph.simulate_weak_acid(
            window, total_acid_mM=total_mM, ph_ex=context.ph_medium,
            beta_mM_per_ph=sweep.beta_mM_per_ph, acid_name=acid_name,
            growth_rate_per_h=mu_env, n_points=n_points)
        ph_final = float(ph_trajectory.of("pH_c")[-1])
        audit = mech_ph.reduction_audit(window, ah_out_mM=ah_out, acid_name=acid_name)
        layers.append(LayerStatus(
            "pH states", LayerState.IN_CHAIN, True,
            f"2 ODE states integrated by {ph_trajectory.method} in "
            f"{ph_trajectory.n_rhs_evals} rhs evals: pH_c "
            f"{float(mech_ph.RESTING_CYTOSOLIC_PH):.2f} -> {ph_final:.3f} over "
            f"{window.duration_h:g} h; tau/T verdicts "
            f"{audit.verdict('A_i').value}/{audit.verdict('pH_c').value}"))

        # ---------------------------------------------------------- proton bill
        # Two ENDS of one bracket: the bill is what holding pH_c at rest costs, the
        # trajectory is where pH_c goes with nothing held. Pma1's Vmax is NOT_BUILT.
        bill = mech_ph.proton_bill(ah_out, mu_env, sweep.cytosolic_volume_ml_per_gdcw, pka)
        mu_after_acid = mech_ph.growth_under_bill(mu_env, bill, regime=regime)
        mechanistic_retained = mu_after_acid / mu_env if mu_env > 0 else 1.0
        layers.append(LayerStatus(
            "proton bill", LayerState.IN_CHAIN, True,
            f"[A-]_i = {bill.trapped_anion_mM:.4g} mM at the resting pH -> "
            f"{bill.mmol_atp_per_gdcw_h:.4g} mmol ATP/gDCW/h at s_ATP = "
            f"{regime.ngam_growth_slope_per_h_per_mmol_atp:.6g} ({regime.name}) -> mu_max "
            f"{mu_env:.4f} -> {mu_after_acid:.4f} /h. UPPER bound on the cost: it holds "
            f"pH_c at rest, while the pH-states layer reports the LOWER bound on pH with "
            f"nothing held. A cell is between them and Pma1's Vmax is NOT_BUILT"))
        layers.append(LayerStatus(
            "stress panel", LayerState.INERT, False,
            f"the panel's empirical curve for {environment.stressor} retains "
            f"{panel_retained:.4f} of mu at this dose and the mechanism retains "
            f"{mechanistic_retained:.4f}; the mechanism REPLACES it here rather than "
            "multiplying, because both tax the same growth rate. The two are an unfitted "
            "comparison and not a calibration -- nothing here was tuned to make them agree"))
        notes.append(
            f"the two pH arms bracket rather than agree: holding pH_c at "
            f"{float(mech_ph.RESTING_CYTOSOLIC_PH):g} costs {bill.mmol_atp_per_gdcw_h:.4g} "
            f"mmol ATP/gDCW/h and no quench; holding nothing costs no ATP and takes pH_c to "
            f"{ph_final:.3f} with a quench of {mech_ph.citrine_quench(ph_final):.4f}. The "
            "chain reports the growth cost from the first and the reporter artefact from "
            "the second, so the pair is a bracket and neither is the cell")
        if not bill.physical:
            notes.append(
                f"the anion pool this dose implies, {bill.trapped_anion_mM:.4g} mM, is "
                "above the osmotic scale a cell runs at -- mech/ph.py flags it and does not "
                "clip it, and a growth rate derived from it is arithmetic rather than a "
                "prediction")
    else:
        for name in ("pH states", "proton bill"):
            layers.append(LayerStatus(
                name, LayerState.NOT_RUN, False,
                "pass an acetic_acid dose and a medium pH; with no weak acid this module's "
                "only proton load is absent and mech/ph.py refuses the basal leak (its "
                "conductance is Neurospora and its driving force needs a membrane potential)"))
        mu_after_acid = mu_env * panel_retained
        layers.append(LayerStatus(
            "stress panel", LayerState.IN_CHAIN, panel_retained < 1.0,
            f"generator/stress_panel.py's empirical curve, unchanged: "
            f"{environment.stressor or 'no stressor'} at {environment.dose:g} retains "
            f"{panel_retained:.4f} of mu -> mu_max {mu_env:.4f} -> {mu_after_acid:.4f} /h. "
            "The mechanism has no entry point for this agent, so the incumbent stands"))

    # ---------------------------------------------------------------- burden
    mu_max = mech_burden.burdened_growth_rate(mu_after_acid, n_copies,
                                              float(REFERENCE_COPIES))
    layers.append(LayerStatus(
        "burden", LayerState.IN_CHAIN, n_copies != float(REFERENCE_COPIES),
        f"{n_copies:g} declared copies at Kafri's measured "
        f"{float(mech_burden.KAFRI_GROWTH_LOSS_PER_COPY):.3g} per copy -> mu_max "
        f"{mu_after_acid:.4f} -> {mu_max:.4f} /h. DECLARED genotype axis: no environmental "
        "input, so this layer does not claim criterion 1"))
    layers.append(LayerStatus(
        "allocation", LayerState.REPORTED, False,
        f"{mech_burden.allocation_gate().summary()} -- the proteome-fraction arm is carried "
        "visible and feeds nothing"))

    # ---------------------------------------------------------------- vessel
    setpoint = environment.growth_rate_setpoint_per_h
    if setpoint is None:
        growth = mu_max
        feed_sets_mu = False
        detail = f"batch: mu is what the cells can do, {growth:.4f} /h"
    elif mu_max <= float(setpoint):
        growth = mu_max
        feed_sets_mu = False
        detail = (f"the setpoint {float(setpoint):g} /h is ABOVE mu_max = {mu_max:.4f} /h, "
                  f"so the feed no longer sets mu and the culture runs at its own maximum. "
                  f"predict.py raises SetpointUnreachable here for a stressor and does not "
                  f"check it at all for a context; this layer reports the realised rate")
        notes.append(
            f"the feed does not set mu in this condition: mu_max = {mu_max:.4f} /h against "
            f"a setpoint of {float(setpoint):g} /h. An exponential fed-batch has no effluent, "
            "so the carbon the culture cannot take up accumulates -- see fba/fedbatch.py")
    else:
        growth = float(setpoint)
        feed_sets_mu = True
        detail = (f"fed-batch: the feed holds mu = {growth:g} /h, with mu_max = "
                  f"{mu_max:.4f} /h leaving {mu_max - growth:.4f} /h of margin")
    layers.append(LayerStatus("vessel", LayerState.IN_CHAIN, True, detail))

    # ---------------------------------------------------------------- population
    model = mech_population.TwoPopulation(p_per_division=sweep.p_per_division,
                                          burden=sweep.burden)
    bearing = mech_population.simulate_population(
        model, window, mu_free_per_h=growth, initial_fraction=initial_fraction,
        n_points=n_points)
    fraction = bearing.of("plasmid_bearing")
    mean_fraction = float(bearing.time_average("plasmid_bearing"))
    layers.append(LayerStatus(
        "population", LayerState.IN_CHAIN, True,
        f"2 ODE states: F {fraction[0]:.4f} -> {fraction[-1]:.4f} over "
        f"{bearing.of('generations')[-1]:.3g} generations, time-average "
        f"{mean_fraction:.4f}"))

    # ---------------------------------------------------------------- metabolism
    product: ProductArm | None = None
    if spec is not None:
        from ..pathway.flux import predict_flux_from_expression
        from ..pathway.solve import solve_pathway

        flux = predict_flux_from_expression(genotype.entry_expression, flux_calibration)
        solution = solve_pathway(spec, flux.flux_mmol_per_gdcw_h, growth, dict(kinetics))
        content = float(solution.terminal.content_mg_per_gdcw)
        product = ProductArm(
            solution=solution,
            content_mg_per_gdcw=content,
            population_content_mg_per_gdcw=content * mean_fraction,
            mean_bearing_fraction=mean_fraction,
            flux_mmol_per_gdcw_h=float(flux.flux_mmol_per_gdcw_h),
        )
        layers.append(LayerStatus(
            "metabolism", LayerState.IN_CHAIN, True,
            f"pathway/solve.py at the mechanistic mu = {growth:.4f} /h -> {content:.8f} "
            f"mg/gDCW per bearing cell, x F = "
            f"{product.population_content_mg_per_gdcw:.8f} mg/gDCW over the culture"))
    else:
        layers.append(LayerStatus(
            "metabolism", LayerState.NOT_RUN, False,
            "pass spec, flux_calibration, kinetics and a genotype"))

    # ---------------------------------------------------------------- reporter
    reporter: ReporterArm | None = None
    if construct is not None:
        params = culture_parameters or DEFAULT_PANEL.get(construct)
        if params is None:
            raise KeyError(
                f"no culture parameters for construct {construct!r}; it is not in "
                f"generator/plate.py::DEFAULT_PANEL ({sorted(DEFAULT_PANEL)}), so pass "
                "culture_parameters explicitly")
        dose = float(environment.dose)
        # The mechanism enters simulate_culture through the knob its own docstring calls
        # "the knob that decouples the two axes". Nothing in the generator is rewritten.
        factor = _nutrient_factor(mu_max, mu_env)
        times = np.linspace(0.0, window.duration_h, n_points)
        out = simulate_culture(times, dose, params, initial_biomass=_INOCULUM_G_PER_L,
                               nutrient_factor=factor)
        activity = np.asarray(out["reporter"], dtype=float)
        if ph_trajectory is None:
            quench = np.ones_like(activity)
        else:
            ph_series = np.interp(times, ph_trajectory.t, ph_trajectory.of("pH_c"))
            quench = np.array([mech_ph.citrine_quench(p) for p in ph_series])
        reporter = ReporterArm(times_h=times, activity=activity,
                               observed=activity * quench, quench=quench,
                               nutrient_factor=factor, construct=construct)
        layers.append(LayerStatus(
            "reporter", LayerState.IN_CHAIN, True,
            f"generator/culture.py::simulate_culture on {construct} at dose {dose:g}, "
            f"nutrient_factor {factor:.4f} carrying the mechanistic growth tax"))
        if ph_trajectory is None:
            # INERT rather than IN_CHAIN: with no acid the multiplier is exactly 1, so this
            # layer ran and provably could not have moved the signal.
            layers.append(LayerStatus(
                "photophysics", LayerState.INERT, False,
                "no weak acid, so pH_c stays at its resting value and the quench is exactly "
                "1.0 -- this layer ran and could not have changed the signal"))
        else:
            layers.append(LayerStatus(
                "photophysics", LayerState.IN_CHAIN, True,
                f"Citrine quench {quench[0]:.4f} -> {quench[-1]:.4f} over the read; endpoint "
                f"fold {reporter.endpoint_fold:.4f} against "
                f"{float(activity[-1] / activity[0]):.4f} unquenched"))
    else:
        for name in ("reporter", "photophysics"):
            layers.append(LayerStatus(
                name, LayerState.NOT_RUN, False,
                "pass a construct from generator/plate.py::DEFAULT_PANEL"))

    # ---------------------------------------------------------------- GEM depth
    if gem_model is not None and acid_name is not None:
        moved = reachable_reactions(gem_model, acid_name)
        layers.append(LayerStatus(
            "GEM depth", LayerState.AUDITS, False,
            f"{len(moved)} reactions move under the acid load, across "
            f"{len({m.compartment for rid in moved for m in gem_model.reactions.get_by_id(rid).metabolites})} "
            "compartments. An audit: a flux redistribution cannot make this chain's "
            "prediction move, only disagree with it"))
    else:
        layers.append(LayerStatus(
            "GEM depth", LayerState.NOT_RUN, False,
            "pass gem_model, an already-loaded cobra model with the acid uncoupling "
            "installed by fba/stress_ph.py::add_weak_acid_uncoupling"))

    # ---------------------------------------------------------------- ethanol
    layers.append(LayerStatus(
        "ethanol", LayerState.REPORTED, False,
        "mech/population.py::EthanolDiagnostic reports accumulation and REFUSES a growth "
        "penalty; five unmeasured constants, all named in its NOT_BUILT registry"))

    ordered = tuple(sorted(layers, key=lambda row: LAYER_ORDER.index(row.name)))
    return ChainResult(
        environment=environment, sweep=sweep, window=window,
        mu_env_per_h=mu_env, mu_max_per_h=mu_max, growth_rate_per_h=growth,
        feed_sets_mu=feed_sets_mu, ah_out_mM=ah_out, ph_states=ph_trajectory,
        proton_bill=bill, population=bearing, product=product, reporter=reporter,
        layers=ordered, notes=tuple(notes))


_PLATE = PlateConditions()
_INOCULUM_G_PER_L = _PLATE.inoculum_for("UPRE1") * _PLATE.gdcw_per_od
"""Starting biomass, g/L, taken from `generator/plate.py`'s own committed conditions.

Read off the plate defaults rather than typed, so the reporter arm inoculates at the density
the real reads did and cannot drift from them.
"""


def _nutrient_factor(mu_max: float, mu_env: float) -> float:
    """The mechanism's growth tax, expressed as `simulate_culture`'s own knob.

    ``simulate_culture`` already applies the incumbent's dose response through
    ``params.growth_rate_at(dose)``; the factor returned here carries ONLY the part the
    mechanism adds, so the two are not multiplied twice. Clipped into (0, 1] because the
    knob is defined there and a mechanism that predicted a speed-up would be reporting the
    absence of a tax rather than a benefit.
    """
    if mu_env <= 0:
        return 1.0
    return float(min(1.0, max(mu_max / mu_env, 1e-9)))


def chain_band(environment: Environment, **kwargs) -> tuple[ChainResult, ...]:
    """:func:`run_chain` at every corner of :meth:`SweepPoint.corners`.

    The shape a number from this layer should be quoted in, and the reason
    :func:`run_chain` has no default sweep. Sixteen runs; the pH integration is the
    expensive one and it costs milliseconds.
    """
    if "sweep" in kwargs:
        raise TypeError("chain_band runs the corners itself; do not pass sweep")
    return tuple(run_chain(environment, sweep=corner, **kwargs)
                 for corner in SweepPoint.corners())


# --------------------------------------------------------------------------------------
# The geometry test
# --------------------------------------------------------------------------------------

def ray_share(surface) -> float:
    """Fraction of a dose x time surface's squared Frobenius norm on its leading ray.

    ``sigma_1^2 / sum(sigma_i^2)`` of the ``(dose, time)`` matrix. It is 1 exactly when
    every dose row is a scalar multiple of one common time profile -- when the surface is a
    straight line through the origin in reading-space rather than a two-dimensional object.
    That is the geometry a fitted latent state has to beat: on a rank-1 surface a random
    rotation of the latent basis reconstructs the data as well as the fitted one, which this
    repository has measured twice.

    Args:
        surface: ``(n_doses, n_times)``. Rows must not be centred -- centring would measure
            the rank of the deviations and the claim is about the surface itself.
    """
    matrix = np.asarray(surface, dtype=float)
    if matrix.ndim != 2 or min(matrix.shape) < 2:
        raise ValueError(
            f"a ray share needs a 2-D surface with at least two doses and two times, got "
            f"shape {matrix.shape}")
    singular = np.linalg.svd(matrix, compute_uv=False)
    total = float(np.sum(singular ** 2))
    if total <= 0.0:
        raise ValueError(
            "the surface is identically zero, so it has no leading ray. A reporter that "
            "never moves is a finding about the loading, not a geometry")
    return float(singular[0] ** 2 / total)


def dose_time_surface(stressor: str, promoter: str, doses: Sequence[float], *,
                      window: Window = PLATE_READ_4H,
                      sweep: SweepPoint | None = None,
                      mechanism: bool = True,
                      quench: bool = True,
                      growth_coupling: bool = True,
                      ph_ex: float = 4.5,
                      n_times: int = 25,
                      regime: MaintenanceRegime = REFERENCE_REGIME) -> np.ndarray:
    """The ``(dose, time)`` surface one transcriptional reporter records, old or new.

    ONE function for both generators, because a geometry comparison between two functions is
    a comparison of two implementations. ``mechanism=False`` is the incumbent exactly: a
    constant promoter activity from `stress_panel.module_response` into the reporter ODE at
    the dose's own growth rate. ``mechanism=True`` adds the two things the mechanism
    supplies -- a growth rate that MOVES over the read as the anion pool fills, and the
    Citrine quench at the cytosolic pH the same states produce.

    The mechanistic branch runs the NO-PUMP end of `mech/ph.py`'s bracket throughout, which
    is the internally consistent choice for a surface: nothing is exported, so the anion
    pool stays physical and the ATP tax it implies is a LOWER bound. :func:`run_chain`
    reports both ends because it reports a growth rate; this reports a shape.

    Args:
        stressor: A key of `stress_panel.STRESSORS`.
        promoter: A transcriptional reporter name, from
            `stress_panel.transcriptional_reporters`. Named ``promoter`` rather than
            ``reporter`` because `tests/test_blank_correction_contract.py` reserves the
            latter for a measured channel, and this is the name of a construct.
        doses: The ladder. At least two rungs.
        window: Declared window; its duration is the read.
        sweep: Required when ``mechanism`` is True.
        mechanism: False for the incumbent surface, True for the assembled one.
        quench: Apply the Citrine pH quench. Setting it False with ``growth_coupling``
            True is how :class:`GeometryResult`'s decomposition separates the instrument
            artefact from the cell's own dynamics, and the answer is that the quench does
            essentially all of it.
        growth_coupling: Let the filling anion pool move the growth rate over the read.
        ph_ex: Medium pH the acid arm runs at. Only reached when the stressor is the
            panel's weak acid.
        n_times: Points on the read grid.
        regime: Which branch prices the ATP bill.

    Returns:
        ``(len(doses), n_times)``.
    """
    if stressor not in STRESSORS:
        raise KeyError(f"no stressor {stressor!r}; have {sorted(STRESSORS)}")
    names = transcriptional_reporters()
    if promoter not in names:
        raise KeyError(f"no transcriptional reporter {promoter!r}; have {sorted(names)}")
    if mechanism and sweep is None:
        raise ValueError(
            "the mechanistic surface needs a declared SweepPoint: the buffering capacity "
            "and the cell-water volume are both SWEPT and neither has a defensible default")

    from ..reporter import ReporterKinetics, simulate_reporter

    column = names.index(promoter)
    loadings = reporter_loadings(names)
    order = list(MODULES)
    times = np.linspace(0.0, window.duration_h, n_times)
    kinetics = ReporterKinetics()
    acid_name = "acetic" if stressor == "acetic_acid" else None
    pka = float(mech_ph.acid(acid_name).pka) if acid_name else float("nan")
    control = panel_growth_rate(stressor, 0.0)
    # The undosed steady state, which is what a culture carries when the dose lands at t=0.
    r0 = _BASAL_ACTIVITY / (control + kinetics.k_deg)

    rows = []
    for dose in doses:
        mu = panel_growth_rate(stressor, float(dose))
        activity = np.array([module_response(stressor, float(dose))[m] for m in order])
        k_synth = float(max((activity @ loadings.T)[column], 0.0))
        if not mechanism or acid_name is None or float(dose) <= 0.0:
            trace = simulate_reporter(times, np.full_like(times, k_synth),
                                      np.full_like(times, mu), kinetics, r0=r0)
            rows.append(trace)
            continue
        states = mech_ph.simulate_weak_acid(
            window, total_acid_mM=float(dose), ph_ex=ph_ex,
            beta_mM_per_ph=sweep.beta_mM_per_ph, acid_name=acid_name,
            growth_rate_per_h=mu, n_points=max(n_times, 201))
        ph_series = np.interp(times, states.t, states.of("pH_c"))
        pool = np.interp(times, states.t, states.of("A_i"))
        anion = pool * np.array([1.0 - mech_ph.undissociated_fraction(p, pka)
                                 for p in ph_series])
        slope = regime.ngam_growth_slope_per_h_per_mmol_atp
        mu_t = (mu / (1.0 + slope * anion * sweep.cytosolic_volume_ml_per_gdcw / 1000.0)
                if growth_coupling else np.full_like(times, mu))
        trace = simulate_reporter(times, np.full_like(times, k_synth), mu_t, kinetics, r0=r0)
        if quench:
            trace = trace * np.array([mech_ph.citrine_quench(p) for p in ph_series])
        rows.append(trace)
    return np.asarray(rows, dtype=float)


_BASAL_ACTIVITY = 0.9
"""`generator/stress_panel.py::_DEFAULT_BASAL`, the constitutive floor a promoter sits on.

Restated rather than imported because it is private there. `tests/test_mech_chain.py` pins
the two together so they cannot drift; if the generator's basal moves, that test fails.
"""


@dataclass(frozen=True)
class GeometryResult:
    """One (stressor, reporter) surface's ray share, before and after the mechanism.

    Args:
        stressor, reporter: Which surface.
        incumbent: Ray share of the old generator's surface.
        mechanistic: Ray share of the assembled one.
        off_ray_gain: ``(1 - mechanistic) / (1 - incumbent)`` -- how much more of the
            surface's energy lies off the leading ray. Stated as a ratio because the
            quantity that matters is the energy a latent state could use, and that is the
            complement rather than the share.
    """

    stressor: str
    reporter: str
    incumbent: float
    mechanistic: float

    @property
    def off_ray_gain(self) -> float:
        base = 1.0 - self.incumbent
        return float("inf") if base <= 0 else (1.0 - self.mechanistic) / base

    def line(self) -> str:
        return (f"  {self.stressor:<18s} {self.reporter:<17s} {self.incumbent:.6f} -> "
                f"{self.mechanistic:.6f}   off-ray x{self.off_ray_gain:.1f}")


def geometry_audit(doses: Sequence[float], *, sweep: SweepPoint,
                   stressor: str = "acetic_acid",
                   window: Window = PLATE_READ_4H,
                   ph_ex: float = 4.5) -> tuple[GeometryResult, ...]:
    """The geometry test on every transcriptional reporter, for one stressor.

    Both surfaces come from :func:`dose_time_surface`, so the comparison is between two
    settings of one function and not between two implementations.
    """
    out = []
    for promoter in transcriptional_reporters():
        old = dose_time_surface(stressor, promoter, doses, window=window,
                                mechanism=False, ph_ex=ph_ex)
        if np.linalg.norm(old) <= 0.0:
            continue
        new = dose_time_surface(stressor, promoter, doses, window=window, sweep=sweep,
                                mechanism=True, ph_ex=ph_ex)
        out.append(GeometryResult(stressor, promoter, ray_share(old), ray_share(new)))
    return tuple(out)


# --------------------------------------------------------------------------------------
# Criterion 4: how far into the network one environmental input reaches
# --------------------------------------------------------------------------------------

def reachable_reactions(model, acid_name: str = "acetic", influx: float = 3.0
                        ) -> tuple[str, ...]:
    """Reaction ids whose parsimonious flux moves under an acid load. Does not modify.

    A thin call into `fba/stress_ph.py::proton_load_depth`, here so that the chain can
    report its own depth without a caller having to know which of the eight ``stress_*``
    modules owns the probe. The model must already carry the uncoupling
    (``add_weak_acid_uncoupling``) and should have its proton budget closed, or the load is
    disposed of at a third of an ATP through pyrophosphate excretion and under-priced.

    Args:
        model: A loaded cobra model with the acid installed.
        acid_name: A key of `fba/stress_ph.py::WEAK_ACIDS`.
        influx: mmol/gDCW/h of undissociated acid forced in.
    """
    from ..fba.stress_ph import proton_load_depth

    return tuple(rid for rid, _name, _delta in proton_load_depth(model, acid_name, influx))


def reaction_depths(model, acid_name: str = "acetic", influx: float = 3.0
                    ) -> dict[int, int]:
    """How many moved reactions sit at each graph distance from the acid's entry step.

    Depth is breadth-first on the metabolite-sharing graph, restricted to the reactions that
    actually moved, seeded at the reactions `fba/stress_ph.py` installed for this acid. It
    answers criterion 4 -- "how many reactions does pH reach, and at what depths" -- with a
    histogram rather than a count, because a count alone cannot tell a transporter cluster
    from a network-wide redistribution.

    **The exact count is not reproducible.** pFBA has alternate optima on this model, so the
    moved set is not a function of the constraints alone: over 8 repeats in one process the
    count runs 484-487, the union is 494 and the intersection 477, i.e. 3.4% of the reactions
    flip between runs. The depth CEILING (6), the compartment count (8) and the shape of the
    histogram are stable. Quote a band, never an integer.

    Returns:
        ``{depth: n_reactions}``, ascending.
    """
    import collections

    from ..fba.stress_ph import WEAK_ACIDS, acid_reaction_ids

    moved = set(reachable_reactions(model, acid_name, influx))
    seed = set(acid_reaction_ids(WEAK_ACIDS[acid_name])) & moved
    if not seed:
        raise ValueError(
            f"none of the reactions installed for {acid_name!r} moved under the load, so "
            "there is no entry point to measure depth from. Check that the uncoupling is "
            "installed on this model and that the proton budget is closed")
    neighbours: dict[str, set[str]] = collections.defaultdict(set)
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            neighbours[metabolite.id].add(reaction.id)
    depth = dict.fromkeys(seed, 0)
    frontier, seen, level = set(seed), set(seed), 0
    while frontier:
        level += 1
        nxt: set[str] = set()
        for rid in frontier:
            for metabolite in model.reactions.get_by_id(rid).metabolites:
                for other in neighbours[metabolite.id]:
                    if other in moved and other not in seen:
                        nxt.add(other)
                        seen.add(other)
                        depth[other] = level
        frontier = nxt
    return dict(sorted(collections.Counter(depth.values()).items()))


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------

def provenance(heading_level: int = 3) -> str:
    """The assembled layer's parameters and its gate, as markdown. Regenerates the doc."""
    from .params import provenance_table, tag_census_table

    hashes = "#" * heading_level
    return "\n\n".join([
        f"{hashes} Assembled free-scalar gate",
        assembled_gate().summary(),
        f"{hashes} Census across the in-chain registries",
        tag_census_table(*IN_CHAIN_REGISTRIES),
        provenance_table(*IN_CHAIN_REGISTRIES, heading_level=heading_level),
    ])
