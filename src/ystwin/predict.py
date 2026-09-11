"""Environment and genotype in, product out — the chain actually joined.

Answers one question: *given a culture environment and a strain's expression cassette, how
much product?* Not "given a measurement of the product".

    genotype (entry-enzyme expression)  ->  pathway flux        pathway.flux     PREDICTS
    environment                         ->  growth rate         generator        BOUNDS
    stressor and dose                   ->  module activities   stress_panel     CLASSIFIES
    flux + growth rate                  ->  every pool          pathway.solve    PREDICTS
    flux + GEM                          ->  feasibility, cost   fba.audit        AUDITS, opt-in
    environment                         ->  mu_max, pH_c, F     mech.chain       BOUNDS, opt-in

**"opt-in" on the last two lines is load-bearing.** Four rows run on every call. The fifth
runs only when a caller passes ``audit_model``, because loading a GSMM costs an 18-second
SBML parse and making every prediction pay it to learn nothing new would be a poor trade --
the audit cannot change the number, only say whether the network could carry it. The sixth
runs only when a caller passes ``mech``, and that one CAN change the number: it applies
``mu = min(mu_set, mu_max)`` on both branches, where this module checked reachability on the
stressor branch only. It is off by default because every leave-one-strain-out figure this
package quotes was measured without it.

For a while this docstring said the audit was in the chain when it was not, which was worse
than either honest option. The argument below exists so that the sentence is true: a caller
who wants the audit passes a model they already loaded, and one who does not pays nothing.

**The circularity is gone.** `pathway_flux` used to be a caller-supplied constant, and both
entry points obtained it as ``q_lycopene + q_beta_carotene`` from the state being predicted.
It is now predicted from relative expression of the pathway's entry enzyme -- one RT-qPCR
number, no product measurement of the strain -- at 1.22x typical fold error leave-one-strain
-out. See :mod:`ystwin.pathway.flux`.

**What the GEM does and does not do.** It audits. Capping every pathway reaction at the
measured flux magnitude leaves the FVA floor at exactly zero, so no capacity constraint
makes a heterologous flux mandatory and no FBA formulation supplies one. What it can say is
whether a predicted flux fits inside the envelope, what it costs in growth, and what the
native competition on the precursor node looks like.

**Scope: intracellular, non-secreted products.** The steady state closes on
``d[X]/dt = v_in - v_out - mu*[X] = 0``, and the ``mu*[X]`` term is growth dilution. A
secreted product leaves through a transporter at a rate the cell sets; the same arithmetic
does not describe it. :class:`~ystwin.pathway.spec.PathwaySpec` refuses those at load time
rather than returning a number that means nothing.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, fields
import math
from types import MappingProxyType

from .fba.fedbatch import SetpointAboveCapacity
from .generator.context import CultureContext, context_growth_rate
from .generator.culture import CRITICAL_GROWTH_RATE_PER_H, chemostat_physiology
from .generator.stress_panel import STRESSORS, growth_rate as stressed_growth_rate
from .generator.stress_panel import module_response
from .mech import engine as _engine
from .mech.adapters import reaction_table, simulate_protocol
from .mech.contracts import (
    EXTRACELLULAR,
    INTRACELLULAR,
    PROMOTER_UNCOUPLED_SIGNALS,
    Control,
    Event,
    Genotype as ProtocolGenotype,
    ObservationModel,
    PhysicalState,
    Protocol,
    SimulationResult,
    Validity,
)
from .mech.engine import EngineParameters, initialize
from .pathway.capacity import (
    CAPACITY_SETTING_ENZYME,
    CapacityUnmeasured,
    capacity_for,
)
from .pathway.enzyme_capacity import (
    ISOPRENOID_KCAT_RANGE_PER_S,
    derived_capacity,
)
from .pathway.flux import FluxCalibration, PredictedFlux, predict_flux_from_expression
from .pathway.solve import (
    NodeKinetics,
    PathwaySolution,
    content_ceiling,
    solve_pathway,
)
from .pathway.spec import Fate, PathwaySpec, RateLaw
from .pathway.thermo_gate import (
    CYTOSOLIC_VOLUMES_ML_PER_GDCW,
    GateReport,
    gate_across_volumes,
    require_feasible,
    volume_dependent_steps,
)

__all__ = [
    "Control",
    "EngineParameters",
    "Environment",
    "Event",
    "Genotype",
    "LayerState",
    "LayerStatus",
    "LegacyProductDiagnostic",
    "MechanisticProductPrediction",
    "MechanisticRun",
    "ObservationModel",
    "PhysicalState",
    "ProductPrediction",
    "Protocol",
    "ProtocolGenotype",
    "SetpointUnreachable",
    "SimulationResult",
    "initialize",
    "predict_product",
    "predict_protocol",
]


def predict_protocol(
    protocol: Protocol,
    genotype: ProtocolGenotype,
    initial_state: PhysicalState,
    parameters: EngineParameters,
    *,
    observation: ObservationModel | None = None,
    **solver_options,
) -> SimulationResult:
    """Run the shared physical clock without invoking empirical or FBA comparisons.

    Supply the complete physical contracts explicitly; synthetic latent coordinates
    and empirical activity folds are not molecular inventories or enzyme abundances.
    The returned ``variables`` retain units, compartments and biomass bases, and
    ``validity`` retains the engine's conditional assumptions and scientific gaps.
    Observations remain separate from simulated truth; no fitted biology is implied.
    """
    return _engine.simulate(protocol, genotype, initial_state, parameters,
                            observation=observation, **solver_options)


class LayerState:
    """What one of the vision's named layers did on a given call.

    The five states exist because "wired in" is not one thing, and collapsing them was how
    this module came to claim the FBA audit was in the chain when it was not. A layer can
    run and set the number, run and only check it, run and be carried unused, run and be
    incapable of changing anything, or not run at all -- and a caller deciding whether to
    trust a prediction needs to know which.

    ``IN_CHAIN`` -- ran, and the returned number depends on it.

    ``AUDITS`` -- ran against the returned number and can agree or refute, never move it.
    Stoichiometry and thermodynamics are both this: an upper bound cannot make a flux
    mandatory, and a driving force cannot set a rate.

    ``REPORTED`` -- ran, and its output is attached and labelled but feeds nothing. This is
    what an unvalidated branch gets. It is the state that keeps a refused layer visible
    instead of deleted.

    ``INERT`` -- ran, and provably could not have changed anything. Distinct from
    ``AUDITS`` because an audit might have refuted and did not; an inert layer never could.

    ``NOT_RUN`` -- the caller supplied nothing for it. The ``detail`` says what to pass.
    """

    IN_CHAIN = "in-chain"
    AUDITS = "audits"
    REPORTED = "reported"
    INERT = "inert"
    NOT_RUN = "not-run"
    ALL = (IN_CHAIN, AUDITS, REPORTED, INERT, NOT_RUN)


@dataclass(frozen=True)
class LayerStatus:
    """One layer's participation in one prediction.

    Args:
        name: The layer, in the vocabulary the project's goal uses.
        state: One of :class:`LayerState`.
        sets_the_number: Whether the returned content depends on this layer. **Two are
            ``True`` under a held setpoint and three in batch** -- `metabolism`, `flux law`,
            and
            `stress panel` when it sets the growth rate. This line read "exactly one" until
            2026-08-31, disagreeing with the code below it and with
            `tests/test_everything_is_wired.py`, which pins ``{"metabolism", "flux law"}``
            and was right. The docstring whose job is to stop stale wiring claims carried
            one for months.

            The count also overstates the independence. `pathway/solve.py` makes the
            terminal content a function of ``flux / mu`` ALONE -- verified to twelve
            significant figures -- so the two layers that set the number set one ratio
            between them, and a genotype that doubles expression is indistinguishable from
            an environment that halves growth.

            **That count is the EMPIRICAL route's, and only its.** ``mode="mechanistic"``
            returns a :class:`MechanisticProductPrediction` whose rows are the engine's, and
            more of them are ``True`` because the engine integrates what the empirical route
            holds fixed: `engine metabolism`, `engine expression`, `engine signalling`,
            `engine vessel`, and whichever environment axes THIS run's own traces show
            reaching the pool -- which is not the same set on every run, as medium pH going
            `INERT` without acetate shows. A `ProductPrediction` never carries those rows.
        detail: What it did, or what the caller would pass to make it run.
    """

    name: str
    state: str
    sets_the_number: bool
    detail: str

    def __post_init__(self) -> None:
        if self.state not in LayerState.ALL:
            raise ValueError(f"state must be one of {LayerState.ALL}, got {self.state!r}")

    def summary(self) -> str:
        return f"{self.name:<22} {self.state:<9} {self.detail}"


def _layer(layers: tuple[LayerStatus, ...], name: str) -> LayerStatus:
    for candidate in layers:
        if candidate.name == name:
            return candidate
    raise KeyError(f"no layer {name!r}; have {[layer.name for layer in layers]}")




class SetpointUnreachable(SetpointAboveCapacity):
    """The vessel cannot hold this growth rate: the dose leaves mu_max at or below it.

    Raised rather than returning a small number. This is the fed-batch replacement for what
    was `SetpointUnreachable`, and the physical situation is the same one seen through a different
    vessel: a chemostat below its critical dilution rate EMPTIES, while an exponential
    fed-batch has no effluent and instead ACCUMULATES the carbon the culture can no longer
    take up. Either way there is no stressed steady state, so reporting a reduced product
    rate would describe a vessel that is not in the state the number assumes.

    A subclass of `fba/fedbatch.py`'s :class:`SetpointAboveCapacity` because it IS that
    condition -- a setpoint above what the strain can do -- reached at prediction time rather
    than at design time. A caller catching the design-time refusal catches this too, which is
    the correct relationship: both say the feed does not set mu.
    """


@dataclass(frozen=True)
class Genotype:
    """What the strain carries, as the model can see it.

    Args:
        entry_expression: Relative expression of the pathway's entry enzyme, on the
            normalisation its flux calibration was fitted on.
        label: Strain name, for reporting.
        cassette: Relative dosage of every heterologous gene the strain carries, keyed by
            enzyme name, on the calibration strain's scale -- so ``1.0`` is the strain the
            capacity was fitted on. Optional, and ``None`` means "not stated", which is
            different from "one copy of everything".

    **Why `cassette` exists, and why it is mostly empty.** Until 2026-08-31 this class held
    one number and the docstring said so: "the whole heterologous cassette enters here". It
    does not, and that single number is the reason the model's sharpest claim is refuted.

    For a terminal node fed by a saturating step the solver gives
    ``content = flux_out/mu`` with ``flux_out <= vmax_per_growth * mu``, so **mu cancels**
    and content is capped at ``vmax_per_growth`` = 1.2483 mg/gDCW, whatever the genotype.
    Raising `entry_expression` a hundredfold buys 0.06%. Arhar 2024 (PMID 39215465) measured
    **79 mg/gDCW** by HPLC. The model is refuted by ~63x and no value of the one input it
    accepts can move it, because the capacity is a property of a DIFFERENT enzyme -- the
    crtYB cyclase, not the crtE the flux law is fitted on.

    So the structure now exists to say which genes a strain carries. What does NOT exist is
    the coefficient relating crtYB dosage to capacity, and the reason was corrected on
    2026-09-04: this paragraph said "28 published strains and every one is qualitative",
    which the repository's own survey refutes -- five rows of
    `data/carotenoid_batch/published_batch_titres.tsv` state a crtYB dosage as a NUMBER, with
    named evidence. They hold only TWO distinct values across five different hosts, promoters,
    media and assays, so the coefficient is unfittable through CONFOUNDING rather than through
    absent numbers. `pathway/capacity.py`'s refusal message quotes the derived reason, which
    `pathway/published_cassettes.py` recomputes from the table at the point of use.

    `pathway/capacity.py` therefore ACCEPTS a cassette and REFUSES to scale by it, naming
    the measurement. That is a smaller lie than a model whose genotype input cannot express
    the question at all.
    """

    entry_expression: float
    label: str = ""
    cassette: Mapping[str, float] | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.entry_expression) or self.entry_expression <= 0:
            raise ValueError(
                f"entry_expression must be positive and finite, got {self.entry_expression}. Zero is "
                "not a small flux, it is a strain that does not carry the pathway")
        for enzyme, dosage in (self.cassette or {}).items():
            if not math.isfinite(dosage) or dosage <= 0:
                raise ValueError(
                    f"cassette dosage for {enzyme!r} must be positive and finite, got {dosage}. A gene "
                    "the strain does not carry is absent from the mapping, not zero -- the "
                    "two mean different things and only one of them is a strain")


@dataclass(frozen=True)
class Environment:
    """One culture condition, as an experiment would specify it.

    Args:
        context: Temperature, carbon source, oxygen, pH, growth phase, strain.
        stressor: Name from ``stress_panel.STRESSORS``, or ``None``.
        dose: Stressor concentration in that stressor's own units.
        growth_rate_setpoint_per_h: The growth rate the vessel holds, /h. When given, the
            growth rate IS this. When ``None`` the growth rate comes from the context as a
            batch maximum.

            **This used to be a chemostat dilution rate and is now an exponential fed-batch
            setpoint.** The arithmetic downstream is unchanged -- the pools law divides by
            mu either way -- but three things about the experiment are not. The setpoint is
            CONVERGED ON rather than imposed instantly, so a run has to settle before it is
            measured (`fba/fedbatch.py` quantifies that: bounded by about 5.6 generations
            however wrong the assumed yield is, against the chemostat's 7.2 residence times).
            There is no effluent line, so the culture cannot wash out -- the analogous
            failure is that the feed stops setting mu and carbon accumulates, which is what
            :class:`SetpointUnreachable` now reports. And mu is MEASURED afterwards from the
            biomass series rather than asserted from a pump, which is the one thing a
            chemostat could never do; see `fba/fedbatch.estimate_growth_rate`.
    """

    context: CultureContext = field(default_factory=CultureContext)
    stressor: str | None = None
    dose: float = 0.0
    growth_rate_setpoint_per_h: float | None = None

    def __post_init__(self) -> None:
        if self.stressor is not None and self.stressor not in STRESSORS:
            raise KeyError(f"no stressor {self.stressor!r}; have {sorted(STRESSORS)}")
        if not math.isfinite(self.dose) or self.dose < 0:
            raise ValueError(f"dose must be finite and non-negative, got {self.dose}")
        if self.stressor is None and self.dose:
            raise ValueError("a dose was given with no stressor to apply it to")
        if (self.growth_rate_setpoint_per_h is not None
                and (not math.isfinite(self.growth_rate_setpoint_per_h)
                     or self.growth_rate_setpoint_per_h <= 0)):
            raise ValueError(
                f"growth_rate_setpoint_per_h must be positive and finite, got "
                f"{self.growth_rate_setpoint_per_h}")
        for name in ("temperature_c", "oxygen", "glucose_g_per_L", "ph_medium", "ph_cytosolic"):
            if not math.isfinite(getattr(self.context, name)):
                raise ValueError(f"context.{name} must be finite")


@dataclass(frozen=True)
class ProductPrediction:
    """What the chain says, with the provenance of each part.

    ``layers`` is the part worth reading first. It records what every one of the five
    layers the project's goal names did on this call, so "is that layer wired in" is a
    question the returned object answers rather than one a reader settles by grepping
    imports -- which is how the claim that the FBA audit was in the chain survived for as
    long as it did while `predict.py` did not import it.
    """

    product: str
    solution: PathwaySolution
    flux: PredictedFlux
    growth_rate_per_h: float
    fermentative: bool
    stress_modules: dict[str, float]
    mode: str = field(default="empirical", kw_only=True)
    genotype: Genotype | None = field(default=None, kw_only=True)
    environment: Environment | None = field(default=None, kw_only=True)
    notes: tuple[str, ...] = ()
    layers: tuple[LayerStatus, ...] = ()
    thermodynamics: tuple[GateReport, ...] = ()
    flux_audit: object | None = None
    latent_constraints: object | None = None
    mech_chain: object | None = None
    """The `mech/chain.py` result when ``mech`` was passed, and ``None`` when it was not.

    Carried the way ``flux_audit`` and ``latent_constraints`` are -- an opaque object rather
    than unpacked fields -- because the mechanistic run holds trajectories, a bearing
    fraction and a reporter arm that this class has no shape for. ``mech_chain.summary()``
    is one line; ``mech_chain.product.population_content_mg_per_gdcw`` is what an assay on
    the whole culture would read, against ``content_mg_per_gdcw`` here, which is per
    BEARING cell.
    """

    def layer(self, name: str) -> LayerStatus:
        return _layer(self.layers, name)

    def layer_report(self) -> str:
        """Every layer, one per line. What a caller prints to see the whole chain."""
        return "\n".join(layer.summary() for layer in self.layers)

    @property
    def content_mmol_per_gdcw(self) -> float:
        return self.solution.terminal.content_mmol_per_gdcw

    def content_interval(self, spec: PathwaySpec,
                         kinetics: dict | None = None,
                         *,
                         kinetic_calibration=None) -> tuple[float, float]:
        """Historical entry-only content-band comparison, in mmol/gDCW.

        By default this maps the selected empirical entry law's one-sigma log spread
        through ``solve_pathway`` at fixed growth and kinetics. The monotone endpoint
        mapping is exact for that conditional calculation; it does not refit the entry
        law or kinetics, propagate observation error, or establish predictive coverage.
        It is a comparison surface, not the uncertainty returned by the current nested
        product validation and not a band from :func:`predict_protocol`.

        Args:
            spec: The pathway used for the empirical comparison.
            kinetics: The same node kinetics used for that comparison.
            kinetic_calibration: Optional legacy local-sensitivity extension for one
                saturating node, with ``joint_log_sd`` and a log-parameter covariance.
                It propagates the capacity/km covariance and combines its log half-width
                with the entry term in quadrature. This linearized extension still omits
                entry/kinetic cross-covariance, model selection, and observation error;
                it must not be presented as a calibrated joint predictive interval.

        Returns:
            ``(low, high)`` in mmol/gDCW, conditional on the supplied empirical model.

        Current conditional joint uncertainty is implemented by
        :func:`ystwin.pathway.flux.product_prediction_intervals` and consumed by
        :func:`ystwin.pathway.flux.score_product_validation`. It resamples training
        strains with their paired conditions and released measurement intervals, refitting
        entry and branch parameters jointly after selection inside the outer training
        fold. Both beta-carotene and lycopene are evaluated on held-out conditions;
        failed predictions and unavailable intervals remain in the denominator.

        ``scripts/content_band_coverage.py`` now reports those actual held-out outcomes
        using ``product_inside_band``, ``lycopene_inside_band``, ``joint_inside_band``
        and ``interval_status``, with separate nominal and observed coverage summaries.
        Its former simulated-width columns do not describe this validation and are no
        longer cited here. Three independent strain groups do not by themselves establish
        nominal predictive coverage, even when a conditional interval is available.
        """
        import dataclasses as _dc
        import math as _math

        # Validate BEFORE solving. A bad call should not pay for two solves first, and the
        # solver's own "node needs parameters" error would otherwise mask this one.
        node = None
        if kinetic_calibration is not None:
            # The spec declares which node saturates; that is the authoritative signal, and
            # it is the only node this calibration describes. Refuse rather than guess if
            # the pathway has none or more than one -- a band widened around the wrong node
            # is arithmetic performed on an unrelated parameter.
            saturating = [n.name for n in spec.nodes if n.rate_law == "saturating"]
            if len(saturating) != 1 or not kinetics or saturating[0] not in kinetics:
                raise ValueError(
                    f"pathway {spec.product!r} declares {len(saturating)} saturating "
                    f"node(s) {saturating}, and the supplied kinetics cover "
                    f"{sorted(kinetics or ())}. `kinetic_calibration` describes exactly one "
                    f"saturating node, so there is nothing to attribute its covariance to")
            node = saturating[0]

        low, high = self.flux.interval()
        lo = solve_pathway(spec, low, self.growth_rate_per_h,
                           kinetics).terminal.content_mmol_per_gdcw
        hi = solve_pathway(spec, high, self.growth_rate_per_h,
                           kinetics).terminal.content_mmol_per_gdcw
        if kinetic_calibration is None:
            return (lo, hi)

        base = kinetics[node]

        def _content(capacity, km):
            bumped = dict(kinetics)
            bumped[node] = _dc.replace(base, vmax_per_growth=capacity, km=km)
            return solve_pathway(spec, self.flux.flux_mmol_per_gdcw_h,
                                 self.growth_rate_per_h,
                                 bumped).terminal.content_mmol_per_gdcw

        # Log-log sensitivities by central difference on the shipped solver, so the gradient
        # belongs to the same solve the band is drawn around rather than to a linearisation
        # of it.
        step = 0.01
        c0, k0 = base.vmax_per_growth, base.km
        centre = _math.log(_content(c0, k0))
        d_cap = (_math.log(_content(c0 * _math.exp(step), k0)) - centre) / step
        d_km = (_math.log(_content(c0, k0 * _math.exp(step))) - centre) / step
        sd = kinetic_calibration.joint_log_sd((d_cap, d_km))

        # Combine IN QUADRATURE, not multiplicatively. The two sources are independent -- a
        # held-out spread over strains and a parameter covariance from one fit -- so their
        # log half-widths add as squares. Multiplying them instead (the first version of
        # this method, for about ten minutes) is the conservative arithmetic and it
        # OVERSTATES: it returned 1.444x at mu = 0.101 against a joint Monte Carlo's 1.181x.
        # A band too wide is a smaller sin than one too narrow and still not the number.
        half_flux = 0.5 * _math.log(hi / lo)
        half = _math.hypot(half_flux, 1.959963985 * sd)
        mid = _math.sqrt(lo * hi)
        return (mid * _math.exp(-half), mid * _math.exp(half))

    @property
    def content_mg_per_gdcw(self) -> float:
        """The unit the literature reports, from the spec's own molar mass."""
        return self.solution.terminal.content_mg_per_gdcw

    @property
    def rate_mmol_per_gdcw_h(self) -> float:
        """Specific production rate. At steady state this is content times growth.

        That identity holds because dilution is the product's ONLY outlet, which is true of
        every pathway this package can currently answer for and is not true in general. A
        terminal node that also degrades or exports loses ``k_loss*[X]`` as well, so the
        rate is ``content * (mu + k_loss)`` and this property would under-report it -- for a
        secreted product by the factor ``mu/(mu + v_sec/X)``, which on Pfeffer 2011's
        chemostat is 6.4%: the number would be sixteen times too small.

        It refuses rather than guessing which number the caller meant, because for a lossy
        terminal there are two defensible answers -- everything made, and the part that
        reaches the supernatant, which is what a titre assay reads -- and they differ by the
        degraded fraction. `PathwaySolution` carries both; a scalar named "the rate" cannot.
        """
        terminal = self.solution.terminal
        if terminal.loss_flux_mmol_per_gdcw_h:
            raise ValueError(
                f"product {self.product!r} has a terminal outlet besides dilution "
                f"({terminal.loss_flux_mmol_per_gdcw_h:.4g} mmol/gDCW/h leaves the chain), "
                f"so 'content times growth' is not the production rate -- it is only the "
                f"washout part of it, and here that is "
                f"{self.growth_rate_per_h * terminal.content_mmol_per_gdcw:.4g}. Ask the "
                f"solution for the number you mean: solution.entry_flux for everything "
                f"made, or solution.terminal.loss_flux_mmol_per_gdcw_h for the part that "
                f"leaves the cell, which is what a supernatant assay reads")
        return self.content_mmol_per_gdcw * self.growth_rate_per_h

    def intermediates(self) -> dict[str, float]:
        """Every pool before the product, mmol/gDCW. Predicted, never supplied."""
        return {n.name: n.content_mmol_per_gdcw for n in self.solution.nodes[:-1]}

    def summary(self) -> str:
        low, high = self.flux.interval()
        return (f"mu {self.growth_rate_per_h:.3f} /h"
                f"{' fermentative' if self.fermentative else ' respiratory'}  ->  "
                f"flux {self.flux.flux_mmol_per_gdcw_h:.4g} [{low:.4g}, {high:.4g}]  ->  "
                f"{self.product} {self.content_mg_per_gdcw:.3g} mg/gDCW "
                f"({self.rate_mmol_per_gdcw_h:.4g} mmol/gDCW/h)")


@dataclass(frozen=True)
class LegacyProductDiagnostic:
    calculation: ProductPrediction
    mode: str = field(default="legacy", init=False)
    supported: bool = field(default=False, init=False)

    def summary(self) -> str:
        return "LEGACY DIAGNOSTIC, not a supported prediction: " + self.calculation.summary()


@dataclass(frozen=True)
class MechanisticRun:
    """Everything `mech/engine.py` needs for one call, declared by the caller.

    **Nothing here has a default, and that is the point.** `mech/adapters.py` refuses to
    convert legacy stress coordinates into physical state -- it would need calibrated
    molecular inventories and gene-specific activities that nobody has measured -- and this
    class does not grow a second, weaker version of that conversion. The protocol, the
    inoculum, the medium and the engine's own genotype are declarations.

    Args:
        protocol: The vessel programme. ``protocol.times_h`` is the ONLY grid the returned
            trajectory exists on, which is why :attr:`read_time_h` must be a member of it.
        genotype: The ENGINE genotype -- gene copies, promoters, per-gene rate constants.
            Not :class:`Genotype`, whose ``entry_expression`` is a relative mRNA level on a
            flux law's normalisation: no measured conversion from that onto a copy number
            exists in this host, the same gap `pathway/capacity.py` refuses on for cyclase
            dosage. Hand the engine its own contract instead.
        initial_state: Inoculum and medium, from :func:`ystwin.mech.engine.initialize`.
        parameters: The engine's parameter set. Its ``thermodynamics`` request is how a
            thermodynamic gate reaches a mechanistic call; ``predict_product``'s own
            ``thermo`` argument audits a `PathwaySolution` this mode never produces.
        observation: Optional reporter model. Observations stay separate from truth, and
            the content this returns always comes from truth.
        read_time_h: Which member of ``protocol.times_h`` to read. Defaults to the last.
        solver_options: ``rtol``, ``atol``, ``max_step_h``, ``method``, forwarded verbatim
            and recorded in ``result.diagnostics``. A numerical declaration, not biology --
            and there is no retry at looser tolerances anywhere below, because a run that
            only integrates once the tolerance is relaxed is a different claim.
    """

    protocol: Protocol
    genotype: ProtocolGenotype
    initial_state: PhysicalState
    parameters: EngineParameters
    observation: ObservationModel | None = None
    read_time_h: float | None = None
    solver_options: Mapping[str, object] = MappingProxyType({})

    def __post_init__(self) -> None:
        # Types are `mech/adapters.py`'s refusal to make; a copy of it here would be a
        # second rule saying the same thing, and the weaker one would win by running first.
        times = tuple(self.protocol.times_h) if isinstance(self.protocol, Protocol) else ()
        if self.read_time_h is not None and times and float(self.read_time_h) not in times:
            raise ValueError(
                f"read_time_h {float(self.read_time_h):g} is not in protocol.times_h "
                f"{tuple(float(t) for t in times)}. The trajectory exists only on that "
                "grid, so add the time to the protocol rather than interpolating a number "
                "the returned run does not contain")

    def simulate(self) -> SimulationResult:
        """Run the shared engine through the one seam that forwards contracts unchanged."""
        return simulate_protocol(self.protocol, self.genotype, self.initial_state,
                                 self.parameters, observation=self.observation,
                                 **dict(self.solver_options))

    def index_of(self, result: SimulationResult) -> int:
        """Position of the read time in the returned trajectory."""
        times = [float(time) for time in result.times_h]
        if self.read_time_h is None:
            return len(times) - 1
        return times.index(float(self.read_time_h))


@dataclass(frozen=True)
class MechanisticProductPrediction:
    """What the 109-state engine says the cell contains, and what that is not.

    Deliberately NOT a :class:`ProductPrediction` and not a subclass of one. The dataflow
    claim here is strong -- the returned number IS the engine's own pool, no conversion --
    and the EVIDENCE claim is absent: 160 of the 173 parameters behind it are priors and
    ``validity.biological_validation`` is False. `scripts/register_prediction.py` writes
    `data/current_claims.json` from a `ProductPrediction`, so making this one would put a
    prior-conditional titre one function call away from a supported claim. The type is the
    lock; ``supported`` cannot be constructed True, and the validity row says it in words.

    ``content_mmol_per_gdcw`` is ``truth["content.<species>"]`` at the read time, which
    `mech/engine.py` already reports as internal mmol over structural biomass gDW. That is
    the same quantity and the same unit `ProductPrediction.content_mmol_per_gdcw` returns,
    so nothing converts it -- but the BASIS is live structural biomass, and a published
    titre is per gravimetric dry weight. The ratio between those is not measured here, so
    ``whole_broth_content_mmol_per_gdw`` is reported beside it and neither is claimed to be
    what an extraction assay would read.
    """

    product: str
    species: str
    read_time_h: float
    content_mmol_per_gdcw: float
    content_mg_per_gdcw: float
    production_rate_mmol_per_gdcw_h: float
    consumption_rate_mmol_per_gdcw_h: float
    net_accumulation_rate_mmol_per_gdcw_h: float
    growth_rate_per_h: float
    biomass_gdw: float
    dead_biomass_gdw: float
    whole_broth_content_mmol_per_gdw: float
    result: SimulationResult
    validity: Validity
    prior_parameters: int
    total_parameters: int
    layers: tuple[LayerStatus, ...] = ()
    notes: tuple[str, ...] = ()
    mode: str = field(default="mechanistic", init=False)
    supported: bool = field(default=False, init=False)

    def layer(self, name: str) -> LayerStatus:
        return _layer(self.layers, name)

    def layer_report(self) -> str:
        """Every layer, one per line, comparable row for row with an empirical call."""
        return "\n".join(layer.summary() for layer in self.layers)

    @property
    def rate_mmol_per_gdcw_h(self) -> float:
        """Specific production rate, when there is only one number that could mean.

        Refuses for a pool the engine also drains, for the reason
        `ProductPrediction.rate_mmol_per_gdcw_h` refuses: everything made and what stays in
        the pool are two defensible answers that differ by the drained fraction, and a
        scalar called "the rate" cannot carry both. In this engine lycopene is always such
        a pool, because psy and lcy both scale with the same crtYB role.
        """
        if self.consumption_rate_mmol_per_gdcw_h:
            raise ValueError(
                f"the engine drains {self.species!r} as well as filling it "
                f"({self.consumption_rate_mmol_per_gdcw_h:.4g} mmol/gDCW/h leaves the "
                f"pool against {self.production_rate_mmol_per_gdcw_h:.4g} entering it), so "
                "'the rate' is ambiguous here. Ask for production_rate_mmol_per_gdcw_h "
                "for everything made, or net_accumulation_rate_mmol_per_gdcw_h for what "
                "stays in the pool")
        return self.production_rate_mmol_per_gdcw_h

    def summary(self) -> str:
        return (f"MECHANISTIC, prior-conditional and not experimentally validated "
                f"({self.prior_parameters}/{self.total_parameters} parameters are priors): "
                f"{self.product} {self.content_mg_per_gdcw:.4g} mg/gDCW "
                f"({self.content_mmol_per_gdcw:.4g} mmol/gDCW) at {self.read_time_h:g} h, "
                f"mu {self.growth_rate_per_h:.3f} /h, production "
                f"{self.production_rate_mmol_per_gdcw_h:.4g} mmol/gDCW/h. Not a supported "
                f"claim and not comparable to a measured titre")


def _growth_rate(environment: Environment, notes: list[str]) -> float:
    """Growth rate, respecting what a held growth rate actually is.

    When the vessel holds a setpoint the feed sets the growth rate and a stressor does NOT
    change it. That is the definition of the held state, mu = mu_set, and the previous
    version of this module broke it: it set ``growth = D`` and then multiplied by a viability
    fraction, reporting mu below the rate the vessel was holding -- a state that does not
    exist. What a stressor changes is the residual substrate and the biomass concentration --
    or, if it pushes mu_max below the setpoint, there is no held state at all.

    THE VESSEL BEHIND THIS CHANGED ON 2026-09-02 AND THE ARITHMETIC DID NOT. It was a
    chemostat, where mu = D is imposed instantly by the pump; it is now an exponential
    fed-batch, where mu CONVERGES on mu_set through the feedback
    d(ln q_S)/dt = mu_set - mu(q_S). For a steady-state solve those are the same number, and
    `fba/fedbatch.py` carries the evidence: the settled residual substrate is the chemostat's
    own Monod residual at the same rate. What differs is the failure at the top end, and only
    the name of the refusal records it -- a chemostat past its critical dilution rate EMPTIES,
    a fed-batch ACCUMULATES the carbon it can no longer take up.

    So the setpoint branch checks reachability and otherwise leaves the rate alone. The batch
    branch is where a stressor slows growth, because there mu is what the cells can do rather
    than what the feed imposes.
    """
    capacity = context_growth_rate(environment.context)
    retained = 1.0
    if environment.stressor is not None:
        unstressed = stressed_growth_rate(environment.stressor, 0.0)
        stressed = stressed_growth_rate(environment.stressor, environment.dose)
        retained = max(0.0, stressed / unstressed) if unstressed > 0 else 0.0
        capacity *= retained
    if not math.isfinite(capacity) or capacity < 0:
        raise ValueError("the environment must produce a finite non-negative growth capacity")
    if environment.growth_rate_setpoint_per_h is None:
        if environment.stressor is not None:
            notes.append(
                f"batch: {environment.stressor} at {environment.dose:g} leaves "
                f"{retained:.0%} of the unstressed growth rate")
        return capacity

    setpoint = float(environment.growth_rate_setpoint_per_h)
    context = environment.context
    condition = (f"context ({context.carbon_source}, {context.temperature_c:g} C, "
                 f"O2 {context.oxygen:g}, {context.growth_phase})")
    if environment.stressor is not None:
        condition += f" with {environment.stressor} at {environment.dose:g}"
    if capacity <= setpoint:
        raise SetpointUnreachable(
            f"{condition} leaves a maximum growth rate of {capacity:.3f} /h, at or below "
            f"the setpoint {setpoint:g} /h. An exponential feed sets mu through the negative "
            "feedback d(ln q_S)/dt = mu_set - mu(q_S), which has no fixed point once the "
            "strain cannot reach the setpoint: the culture grows at its own maximum and "
            "the carbon it cannot take up accumulates. There is no held steady state to "
            "report, and a reduced product rate here would describe a vessel filling "
            "with unconsumed substrate")
    notes.append(f"fed-batch: growth rate is held at the setpoint mu_set = {setpoint:g} /h")
    notes.append(
        f"{condition} leaves mu_max = {capacity:.3f} /h, above the setpoint, so the feed "
        "still sets mu. Context and dose audit reachability rather than changing the held rate")
    return setpoint


# ---------------------------------------------------------------------------
# the layers
#
# One function per layer, each owning its own LayerStatus. The point is not the line count:
# `LayerState` describes five ways a layer can relate to the answer, and while the statuses
# were built inline nothing owned that relationship. A new input reached `predict_product`
# and its status line was written -- or forgotten -- by whoever remembered.
# `Genotype.cassette` was forgotten, and `_cassette_layer` carries the record of it. A layer
# with a function has somewhere for its status to live.
# ---------------------------------------------------------------------------


def _metabolism_layer(solution) -> LayerStatus:
    """IN_CHAIN: the closed-form solve that produces the number."""
    # The residual is REPORTED, not branched on. `LayerState` has exactly five members --
    # IN_CHAIN / AUDITS / REPORTED / INERT / NOT_RUN -- and they describe how a layer that
    # RAN relates to the answer, not whether its arithmetic held. A broken identity is an
    # arithmetic fault and belongs in an exception, which is where `solve_pathway` now puts
    # it; inventing a sixth state to carry it would make every consumer of `LayerState.ALL`
    # wrong. What this line must not do is print a large residual beside the word "closes",
    # which it did for a day: `carbon_closes` omitted every non-dilution outlet, so a
    # storage chain reported 80% of its entry flux here as though it were closure.
    return LayerStatus(
        "metabolism", LayerState.IN_CHAIN, True,
        f"{len(solution.nodes)} nodes solved in closed form; carbon closes to "
        f"{abs(solution.carbon_closes):.1e}")


def _flux_law_layer(flux, genotype) -> LayerStatus:
    """IN_CHAIN: expression to flux -- the step that removed the circularity."""
    return LayerStatus(
        "flux law", LayerState.IN_CHAIN, True,
        f"{flux.calibration.entry_enzyme} expression {genotype.entry_expression:g}"
        f" -> {flux.flux_mmol_per_gdcw_h:.4g} mmol/gDCW/h, "
        f"{flux.calibration.typical_fold_error:.2f}x typical LOSO error")


def _stress_panel_layer(environment, modules, chain=None) -> LayerStatus:
    """IN_CHAIN in batch, AUDITS under a setpoint the feed actually holds.

    AUDITS, not INERT, and the distinction is the enum's own. `INERT` is "ran, and provably
    could not have changed anything"; this layer can END the call, because a dose past the
    lethal threshold raises `SetpointUnreachable` rather than returning a number. That is a refusal
    channel. Labelled INERT until 2026-08-31, which understated it in the one direction
    that matters.

    ``chain`` is the `mech/chain.py` result when the caller opted in. It is read for one bit
    -- whether the feed held -- because with the mechanism on a setpoint is a CEILING and a
    dose that takes mu_max below it puts the panel back in the chain. Reading that bit off
    the run rather than re-deriving it here is what keeps the two from drifting; the count
    of moved modules is `abs`, for the reason `stress_modules` is.
    """
    setpoint_holds = (environment.growth_rate_setpoint_per_h is not None
                      and (chain is None or chain.feed_sets_mu))
    batch = not setpoint_holds
    if environment.stressor is None:
        detail = "no stressor"
    elif setpoint_holds:
        moved = len([v for v in modules.values() if abs(v) > 0.01])
        detail = (f"{environment.stressor} at {environment.dose:g} moves {moved} "
                  "modules; under a fed-batch setpoint the feed fixes mu and mu is its "
                  "only route to the product")
    elif chain is not None and environment.growth_rate_setpoint_per_h is not None:
        detail = (f"{environment.stressor} at {environment.dose:g} takes mu_max to "
                  f"{chain.mu_max_per_h:.4f} /h, at or below the setpoint "
                  f"{float(environment.growth_rate_setpoint_per_h):g} /h, so the feed no "
                  "longer fixes mu and the pools are solved at the realised rate")
    else:
        detail = (f"{environment.stressor} at {environment.dose:g} sets the batch "
                  "growth rate, which the pools are solved at")
    return LayerStatus("stress panel",
                       LayerState.IN_CHAIN if batch else LayerState.AUDITS,
                       batch, detail)


def _cassette_layer(genotype, spec, ceiling, legacy=False) -> LayerStatus:
    """What `pathway/capacity.py` said about this strain's cassette. Never sets the number.

    THE NEWEST ACCEPTED INPUT, AND IT REACHED NOTHING. `Genotype.cassette` was added on
    2026-08-31 so the model could express which genes a strain carries -- and it was added
    without a status line, which is precisely the failure the `layers` mechanism exists to
    prevent. Caught by the architecture review the same day; given a function of its own on
    2026-08-31 so that the next input has an obvious place to be forgotten in.

    **And then the status line itself was the stale claim.** Until 2026-09-02 both branches
    below returned a hand-written `LayerStatus`, and the string NAMED `pathway/capacity.py`
    -- "names crtYB as the enzyme that sets the ceiling and refuses to scale by it" -- while
    this module did not import it and `capacity_for` was called from nowhere but its own
    test. A layer that reports what another module would say, without asking it, is the same
    defect as a docstring that reports what the code would do: it is right until the day the
    other module changes, and nothing here would notice. Both the state and the detail are
    now what the call returned.

    Wiring it in must NOT let the cassette move the number, and that is the whole design of
    :func:`~ystwin.pathway.capacity.capacity_for`: it ACCEPTS a cassette and REFUSES to
    scale by it, because the flux law is fitted on crtE while the ceiling is set at crtYB
    and the coefficient relating cyclase dosage to capacity does not exist. So the anchor
    passed in is this call's own ceiling, the value that comes back is the anchor unchanged,
    and what changes with the cassette is which of two things the layer reports:

    ``INERT`` -- capacity.py answered, with the calibration strain's own capacity. It ran
    and provably could not have moved anything, which is the enum's definition.

    ``REPORTED`` -- capacity.py REFUSED, and the refusal is attached verbatim. `REPORTED` is
    the enum's state for output that is "attached and labelled but feeds nothing... the
    state that keeps a refused layer visible instead of deleted", which is exactly what a
    `CapacityUnmeasured` is. The refusal is caught rather than propagated for the same
    reason `_latent_branch` attaches its G4 verdict instead of raising: the prediction is
    still computed against the calibration strain's capacity, and saying so loudly beats
    refusing to answer a question the model is otherwise able to express. What the caller
    must not be able to do is read the number without reading that.

    No sixth state was invented for the refusal, and none should be: `LayerState.ALL` is
    consumed elsewhere and a new member makes every consumer of it wrong.
    """
    stated = (", ".join(f"{gene} {dosage:g}x"
                        for gene, dosage in sorted(genotype.cassette.items()))
              if genotype.cassette else "none stated")
    if ceiling is None:
        # No anchor to ask about. `content_ceiling` returns None where nothing saturates at
        # the end of the chain, and capacity.py's argument is a measured capacity -- passing
        # a stand-in for one would be the invention this layer exists not to make.
        return LayerStatus(
            "cassette dosage", LayerState.INERT, False,
            f"cassette ({stated}): pathway {spec.product!r} has no saturating step at the "
            "end of its chain, so this calibration asserts no capacity for a dosage to be "
            "asked about")
    try:
        capacity = capacity_for(spec.product, ceiling, genotype.cassette)
    except CapacityUnmeasured as refusal:
        if not legacy:
            raise
        # Verbatim and on one line: `layer_report` is one line per layer, and paraphrasing a
        # refusal is how the paraphrase and the refusal drift apart.
        return LayerStatus(
            "cassette dosage", LayerState.REPORTED, False,
            f"cassette ({stated}): pathway/capacity.py REFUSED, and the content above is "
            f"still solved against the calibration strain's capacity -- "
            + " ".join(str(refusal).split()))
    enzyme = CAPACITY_SETTING_ENZYME.get(spec.product)
    if enzyme is None:
        why = (f"no enzyme is named as setting {spec.product!r}'s ceiling, so there is no "
               "dosage for it to read")
    elif not genotype.cassette:
        why = "no cassette was stated, so this is the strain the anchor was measured on"
    elif enzyme not in genotype.cassette:
        why = f"the cassette states no dosage for {enzyme}, the enzyme that sets the ceiling"
    else:
        why = f"{enzyme} is at the calibration strain's own dosage"
    return LayerStatus(
        "cassette dosage", LayerState.INERT, False,
        f"cassette ({stated}): pathway/capacity.py answers {capacity:.4g} mmol/gDCW, the "
        f"calibration strain's own capacity, because {why}. It scales by no dosage, and "
        "would refuse rather than invent one")


def _enzyme_capacity_layer(enzyme_capacity, ceiling, growth_rate) -> LayerStatus:
    """AUDITS: `vmax = kcat * [E]` beside the fitted capacity, and never instead of it.

    The fitted `vmax_per_growth` is scored leave-one-strain-out against six measured states;
    this one is scored against nothing, so it reports and refuses. What it CAN do is
    contradict -- and on the shipped pathway it does, in the same direction and by about the
    same factor as Arhar 2024's measurement, which it reaches from the enzyme side instead
    of from a paper.

    NOT_RUN rather than INERT when no kcat is supplied, because unlike E-Flux this layer is
    not structurally incapable of saying something -- it is simply not holding the two
    measurements it needs, and `LayerState.NOT_RUN`'s own contract is "pass X and it runs".
    """
    if enzyme_capacity is None:
        return LayerStatus(
            "enzyme capacity", LayerState.NOT_RUN, False,
            "pass enzyme_capacity={'enzyme':..., 'kcat_per_s':..., "
            "'enzyme_mmol_per_gdcw':...}; the content comes from "
            "pathway/proteome.py's enzyme_content_mmol_per_gdcw (a molar ppm) or "
            "enzyme_content_from_mass_fraction (a PRM mass fraction). It audits the fitted "
            "ceiling from kcat and enzyme content and cannot move the number")
    derived = derived_capacity(**enzyme_capacity)

    low, high = ISOPRENOID_KCAT_RANGE_PER_S
    outside = "" if low <= derived.kcat_per_s <= high else (
        f"; NOTE kcat {derived.kcat_per_s:.3g}/s is outside the min and max the vendored "
        f"ecModel gives the host's whole ergosterol branch, [{low:g}, {high:g}] /s, so this "
        f"is a claim about an enzyme unlike any of them. That interval spans six orders of "
        f"magnitude, so this NOTE catches a unit slip and nothing finer -- a kcat INSIDE it "
        f"is not thereby anchored. It read [0.076, 46.0] until 2026-09-04, when it turned "
        f"out to be two named enzymes rather than a min/max")
    if ceiling is None:
        return LayerStatus(
            "enzyme capacity", LayerState.AUDITS, False,
            derived.summary() + "; no fitted ceiling to compare against, because nothing "
            "saturates at the end of this chain" + outside)
    return LayerStatus(
        "enzyme capacity", LayerState.AUDITS, False,
        derived.summary(ceiling, growth_rate) + outside)


def _environment_layer(spec, environment, growth_rate_per_h) -> LayerStatus:
    """The carbon-source channel. REPORTED on both branches, and never IN_CHAIN.

    This layer exists because the chain was measured to have exactly ONE route from the
    environment to the product -- mu -- and yeast does not agree. Kocharin 2013 moved PHB
    content 3.82x between glucose and ethanol at an IDENTICAL growth rate, one genotype, one
    vessel, carbon matched, 64 standard deviations apart. `scripts/product_environment_sweep.py`
    shows the chain answering that with zero.

    REPORTED and not INERT where no coefficient exists. `INERT` is "ran, and provably could
    not have changed anything", and that is false here: the channel is real and large, it is
    simply unmeasured FOR THAT PRODUCT. Labelling it INERT would record the absence of a
    coefficient as the absence of an effect, which is the confusion this layer is built to
    prevent.

    REPORTED **and not IN_CHAIN** where a coefficient does exist, which is the correction of
    2026-09-04. This layer previously returned `IN_CHAIN, sets_the_number=True` for `phb` and
    printed "ethanol -> 4.679x on content" -- while `predict_product` returned a BIT-IDENTICAL
    1.7363905332 on both carbon sources. The factor was computed, formatted into the detail
    string, and never applied. `IN_CHAIN` is a contract that the returned number depends on
    the layer, and it was the one layer in nine where that was false.

    It was fixed by relabelling rather than by wiring the factor in, and the reason is not
    caution -- it is that wiring it could not have produced information:

      * `phb.toml` declares every node above the product passthrough, so `solve_pathway`
        returns `content = v_in / mu`, **which is its own input**. The spec says it outright:
        "no content this pathway can produce disagrees with any measurement at any parameter
        value". Multiplying an unfalsifiable quantity by 4.679 leaves it unfalsifiable. It
        would agree with Kocharin only by restating Kocharin.
      * The coefficient cannot be scored. Three carbon sources give a permutation floor of
        1/3! = 0.1667 and the fit ATTAINS it exactly; leave-one-carbon-source-out is
        extrapolation from two points to a third, not validation. `LayerState.REPORTED` is
        documented as "what an unvalidated branch gets", which is this branch precisely.
      * Kocharin's 3.82x is worth more as a HELD-OUT TARGET than as a fitted coefficient.
        Spending it here would leave the mechanistic route -- GEM precursor yield, which
        transfers across products, times the promoter's own measured stress response -- with
        nothing independent left to be scored against.

    What would earn `IN_CHAIN` here: a factor derived from host measurements rather than from
    this product's own titres, applied to `v_in`, and scored against Kocharin as a state the
    fit never saw.
    """
    from .pathway.environment_flux import (
        EnvironmentUnmeasured,
        MEASURED_ENVIRONMENT_RESPONSES,
        environment_factor,
    )

    # `CultureContext` names a carbon source; the response is fitted on the ethanol share of
    # feed CARBON. Only the two Kocharin feeds map, and galactose deliberately does not --
    # nothing has measured it, and interpolating a third sugar onto a two-point axis would
    # invent the number this layer exists to withhold.
    share = {"glucose": 0.0, "ethanol": 1.0}.get(environment.context.carbon_source)
    if share is None:
        return LayerStatus(
            "environment", LayerState.NOT_RUN, False,
            f"carbon source {environment.context.carbon_source!r} has no measured carbon "
            f"descriptor; the response is fitted on glucose and ethanol only")

    try:
        factor = environment_factor(spec.product, share, growth_rate_per_h)
    except EnvironmentUnmeasured:
        available = ", ".join(sorted(MEASURED_ENVIRONMENT_RESPONSES))
        return LayerStatus(
            "environment", LayerState.REPORTED, False,
            f"carbon source is NOT in the chain for {spec.product!r}: no measured response, "
            f"and the fitted one ({available}) belongs to an acetyl-CoA-limited pathway fed "
            f"by ethanol through ACS. The channel is real and large -- 3.82x on PHB at a "
            f"fixed mu of 0.05 /h -- so this number is a mu-only answer on an axis where the "
            f"truth is known to move")

    return LayerStatus(
        "environment", LayerState.REPORTED, False,
        f"carbon source {environment.context.carbon_source} -> a measured {factor:.3f}x at "
        f"mu = {growth_rate_per_h:.4g} /h, REPORTED and NOT applied to the content above. "
        f"The coefficient attains its own permutation floor exactly (3 carbon sources, "
        f"1/3! = 0.1667), so it is unscored; and {spec.product!r} solves to v_in/mu, which "
        f"is its own input, so applying it could not make the number falsifiable")


def _core_layers(solution, flux, genotype, environment, modules, spec,
                 ceiling, enzyme_capacity=None, chain=None, legacy=False) -> list[LayerStatus]:
    """The six layers that run on every call, in the order the chain runs them."""
    return [
        _metabolism_layer(solution),
        _flux_law_layer(flux, genotype),
        _stress_panel_layer(environment, modules, chain),
        _environment_layer(spec, environment, solution.growth_rate_per_h),
        _cassette_layer(genotype, spec, ceiling, legacy),
        _enzyme_capacity_layer(enzyme_capacity, ceiling,
                               solution.growth_rate_per_h),
    ]


def _gem_audit(spec, solution, modules, audit_model, audit_reaction, audit_glucose_uptake,
               layers: list[LayerStatus], notes: list[str], legacy=False):
    """The GEM audit and the E-Flux layer that runs inside it. Returns the verdict or None.

    AUDITS and INERT respectively: an upper bound cannot make a flux mandatory, and E-Flux
    scales ceilings only.
    """
    verdict = None
    if audit_model is not None:
        if audit_reaction is None or audit_glucose_uptake is None:
            raise ValueError(
                "audit_model needs audit_reaction and audit_glucose_uptake alongside it. "
                "The reaction because this function cannot know which of the model's "
                "reactions is the product, and the uptake because the growth cost the audit "
                "reports swings sevenfold with it and is meaningless without its regime")
        from .fba.audit import audit_predicted_flux

        # The caller names the terminal product reaction, not the entry reaction.
        # Entry flux also supplies intermediate dilution/losses; assigning it to
        # the terminal demand audits a different number from the solved product.
        verdict = audit_predicted_flux(
            audit_model, audit_reaction, solution.terminal.flux_in,
            glucose_uptake=audit_glucose_uptake,
            precursor_metabolite=spec.precursor_metabolite,
            module_activity=modules or None)
        if not verdict.within_envelope and not legacy:
            raise ValueError("the predicted flux is outside the GEM feasible envelope; no product prediction is supported")
        notes.append(f"GEM audit: {verdict.summary()}")
        notes.append(
            "GEM audit uses terminal production in mmol/gDCW/h on the caller-declared "
            "unit-product reaction. It is a single-flux check at the audit's growth-fraction "
            "regime, not a joint certification of intermediate losses or the supplied growth setpoint.")
        notes.extend(verdict.notes)
        if not verdict.within_envelope:
            notes.append(
                "the predicted flux is OUTSIDE what the network can carry at this growth "
                "rate and uptake, so stoichiometry refutes it -- this is the one thing the "
                "audit can say that the prediction cannot")
        layers.append(LayerStatus("GEM / FBA", LayerState.AUDITS, False, verdict.summary()))
        # E-Flux runs INSIDE the audit, which is the only place a ceiling can do anything,
        # and `regulation` records whether it did. It reports `unstressed` when no module
        # activity was supplied, and otherwise a summary in which the binding count is
        # reliably zero -- every transcriptional module a stressor touches is induced, and
        # E-Flux applies upper bounds only. INERT rather than AUDITS for that reason: an
        # audit might have refuted and did not, and this could not have.
        layers.append(LayerStatus(
            "regulation (E-Flux)", LayerState.INERT, False,
            f"{verdict.regulation} -- scaling a ceiling cannot lower a flux sitting "
            f"{verdict.headroom:.0f}x below it"))
    else:
        for name in ("GEM / FBA", "regulation (E-Flux)"):
            layers.append(LayerStatus(
                name, LayerState.NOT_RUN, False,
                "pass audit_model, audit_reaction and audit_glucose_uptake; E-Flux runs "
                "inside the audit because a ceiling is the only thing it can move"))
    return verdict


def _thermodynamic_gate(spec, solution, thermo, thermo_model, audit_model,
                        thermo_background_m, thermo_volumes,
                        refuse_thermodynamically_blocked,
                        layers: list[LayerStatus], notes: list[str]):
    """The driving-force gate across cytosolic volumes. Returns the reports, possibly empty.

    AUDITS: a driving force says whether a step may run, never how fast, so the only thing
    this can add is a contradiction.
    """
    reports: tuple[GateReport, ...] = ()
    if thermo is not None:
        steps = spec.thermo_steps
        if not steps:
            layers.append(LayerStatus(
                "thermodynamics", LayerState.NOT_RUN, False,
                f"pathway {spec.product!r} declares no step chemistry -- give a node a "
                "`metabolite` and a `reaction` or `stoichiometry` in its TOML"))
        else:
            model = thermo_model if thermo_model is not None else audit_model
            # Only when a model is at hand, and it raises rather than warns. A spec whose
            # declared chemistry has drifted from the model's would gate the old reaction
            # under the new one's name, and every verdict below would be about neither.
            crosschecked = spec.check_steps_against(model) if model is not None else ()
            reports = gate_across_volumes(
                solution, volumes=thermo_volumes, steps=steps,
                metabolites=spec.thermo_metabolites, thermo=thermo,
                background_m=thermo_background_m, model=model)
            conditional = volume_dependent_steps(reports)
            first = reports[0]
            notes.append(
                f"thermodynamic gate: {first.summary()}. It cannot move the number -- a "
                "driving force says whether a step may run, never how fast -- so what it "
                "adds is the one thing the solver could not check about its own output")
            if conditional:
                notes.append(
                    f"and {list(conditional)} change verdict across "
                    f"{list(thermo_volumes)} mL/gDCW, so those rest on the volume "
                    "convention rather than on the chemistry and are not verdicts")
            blocked = sorted({s.node for r in reports for s in r.cannot_run})
            if blocked:
                notes.append(
                    f"THE SOLVE CONTRADICTS ITSELF at {blocked}: carbon was moved through a "
                    "step whose driving force points the other way at the concentrations "
                    "this same solve produced")
            if refuse_thermodynamically_blocked:
                for report in reports:
                    require_feasible(report)
            covered = len(first.runs) + len(first.cannot_run)
            layers.append(LayerStatus(
                "thermodynamics", LayerState.AUDITS, False,
                f"{covered}/{len(steps)} declared steps carry an energy; "
                f"{len(first.cannot_say)} cannot say, {len(first.ungated)} not gated; "
                f"{len(conditional)} verdicts depend on the volume convention; "
                f"{len(crosschecked)} step(s) cross-checked against the model"))
    else:
        layers.append(LayerStatus(
            "thermodynamics", LayerState.NOT_RUN, False,
            "pass thermo=ThermodynamicData.load(), and thermo_background_m for anything "
            "the solve does not pin"))
    return reports


def _latent_branch(latent, growth, physiology, layers: list[LayerStatus],
                   notes: list[str]):
    """The maintenance-ATP branch, computed and attached rather than used.

    REPORTED: it is the only route by which stress could reach the product number, and G4
    refused the anchor that would validate it.
    """
    constraints = None
    if latent is not None:
        from .bridge.latent_bridge import latent_constraints as _latent_constraints
        from .bridge.physiology_bridge import PhysiologicalConstraints

        constraints = _latent_constraints(
            latent,
            PhysiologicalConstraints(
                growth_lower_bound=growth,
                glucose_uptake=abs(physiology.glucose_uptake),
                oxygen_uptake=abs(physiology.oxygen_uptake),
                time_h=latent.time_h),
            acknowledge_unvalidated=True)
        notes.append(
            f"latent stress branch COMPUTED AND NOT USED: {constraints.summary()}. It is "
            "the only route by which stress could reach this product number, and G4 "
            "refused the anchor that would validate it -- so it is attached and labelled "
            "rather than fed in. The content above is bit-for-bit what it would be with "
            "`latent=None`")
        layers.append(LayerStatus(
            "latent stress state", LayerState.REPORTED, False,
            f"{constraints.atp_maintenance:.3f} mmol ATP/gDW/h from "
            f"{latent.construct} -- {constraints.validation}"))
    else:
        layers.append(LayerStatus(
            "latent stress state", LayerState.NOT_RUN, False,
            "pass latent=LatentState(...); it is reported, never used -- G4 refused it"))
    return constraints


# The mech rows whose only route to this content is mu_max, so the feed can cut all five.
# `pathway/solve.py` makes the terminal content a function of (flux, mu) alone.
_MECH_ROWS_UPSTREAM_OF_MU_MAX = frozenset({
    "environment", "stress panel", "weak acid", "proton bill", "burden"})

# `vessel` picks mu = min(mu_set, mu_max) and `metabolism` solves the pools at it, so both
# reach the content on either branch of that min. No other mech row can move the number.
_MECH_ROWS_THAT_SET_THE_CONTENT = _MECH_ROWS_UPSTREAM_OF_MU_MAX | {"vessel", "metabolism"}


def _mech_layers(environment, mech, spec, genotype, flux_calibration, kinetics,
                 layers: list[LayerStatus], notes: list[str]):
    """The mechanistic chain, opt-in, and the one thing it is allowed to move.

    Runs `mech/chain.py::run_chain` on this environment and copies every row it reports into
    the ``layers`` record under a ``mech `` prefix, so a caller prints one report and gets
    both stacks. The prefix is not decoration: `mech/chain.py` names three of its layers
    `environment`, `stress panel` and `metabolism`, which are three of this module's own,
    and two rows under one name is how `layer()` starts answering the wrong question.

    **What it moves, and it is exactly one thing.** The chain's structural correction to
    this module is ``mu = min(mu_set, mu_max)`` on BOTH branches -- `_growth_rate` returns a
    setpoint unchecked when no stressor is present, so an ethanol fed-batch at mu_set = 0.18
    /h was accepted although the context puts that culture's maximum at 0.14 /h. With
    ``mech`` on, the solve runs at the chain's realised rate. Nothing else changes: the flux
    law, the spec and the kinetics are the same objects, and the chain calls the same
    `solve_pathway`.

    **The rows are copied but their `sets_the_number` is RECOMPUTED, and an `IN_CHAIN` that
    fails it is demoted to `REPORTED`**, because the two classes are answering about
    different numbers. The chain's headline is titre over the whole culture, so its
    `population` row sets it; `ProductPrediction`'s content is per BEARING cell and the
    bearing fraction does not appear in it. Copying that True verbatim would put a false
    claim in the field whose whole job is to carry a true one, and leaving the row IN_CHAIN
    beside a False would break this module's own definition of the word. What the chain
    called it is kept in the detail, so nothing is lost by the demotion.

    **And whether a row reaches the content is a fact about THE RUN, not about the row's
    name.** `mech/chain.py` already decides its own flags per run -- ``panel_retained <
    1.0`` for the stress panel, ``n_copies != REFERENCE_COPIES`` for burden -- and this
    module owes the same. When the feed holds mu at the setpoint, ``mu = min(mu_set,
    mu_max)`` returns the setpoint, mu_max never reaches the pools, and every row upstream
    of it is out of THIS chain however hard it worked: at 60 mM acetic acid, medium pH 4.5,
    setpoint 0.18 /h the content is 1.0709878838715265 mg/gDCW at [AH]_o = 0 and unchanged
    in its last digit at [AH]_o = 193.6 mM. A name allowlist alone cannot see that, and it
    printed `stress panel` False beside `weak acid` True for the same mu. The gate is
    `ChainResult.feed_sets_mu` rather than a shorter list of names because the acid arm
    genuinely DOES set the content once the feed lets go -- medium pH 4.5, 300 mM, setpoint
    0.2543 /h: mu_max 0.22142 /h, content 1.0191 against 0.9766 at the held setpoint.

    **The assembled free-scalar gate is REFUSED and the refusal is attached rather than
    raised.** `mech/chain.py::assembled_gate` counts 7 free scalars against 3 independent
    targets. It is reported instead of raising because the refused scalars provably do not
    reach this number: buffering capacity and cytosolic volume enter only the acid arm's
    growth cost, loss-per-division and burden only the bearing fraction, and the content is
    invariant to all sixteen corners of the declared sweep -- see
    `tests/test_predict_mechanistic.py`. What the caller must not be able to do is quote the
    chain's OWN population titre without the refusal, which is why it is a layer.
    """
    from .mech.chain import assembled_gate, run_chain

    result = run_chain(environment, sweep=mech, genotype=genotype, spec=spec,
                       flux_calibration=flux_calibration, kinetics=kinetics)
    for row in result.layers:
        # IN_CHAIN means "the returned number depends on it" HERE, and whether it does is a
        # fact about this run: a feed that holds mu cuts off everything upstream of mu_max.
        reaches_the_content = row.name in _MECH_ROWS_THAT_SET_THE_CONTENT and not (
            result.feed_sets_mu and row.name in _MECH_ROWS_UPSTREAM_OF_MU_MAX)
        sets = row.sets_the_number and reaches_the_content
        demoted = row.state == LayerState.IN_CHAIN and not sets
        layers.append(LayerStatus(
            f"mech {row.name}",
            LayerState.REPORTED if demoted else row.state, sets,
            row.detail + (" [mech/chain.py calls this in-chain; it does not reach the "
                          "content above]" if demoted else "")))
    layers.append(LayerStatus(
        "mech free scalars", LayerState.AUDITS, False,
        f"{assembled_gate().summary()}; the refused scalars reach the bearing fraction and "
        f"the acid arm, and the content above is identical at all 16 declared corners"))
    notes.extend(result.notes)
    setpoint = environment.growth_rate_setpoint_per_h
    if setpoint is not None and not result.feed_sets_mu:
        notes.append(
            f"MECHANISM CLIPPED THE SETPOINT: mu_max = {result.mu_max_per_h:.4f} /h against "
            f"a setpoint of {float(setpoint):g} /h, so the pools are solved at "
            f"{result.growth_rate_per_h:.4f} /h and not at the setpoint. Without `mech` this "
            "call returns the setpoint unchecked with no stressor, and raises "
            "SetpointUnreachable with one -- the mechanism reports the realised rate on "
            "both branches instead, which is the arithmetic fba/fedbatch.py already "
            "enforces at design time")
    notes.append(
        f"mech ON at sweep corner {mech.label()}: mu {result.mu_env_per_h:.4f} -> "
        f"{result.mu_max_per_h:.4f} -> {result.growth_rate_per_h:.4f} /h. Content below is "
        f"per BEARING cell; over the whole culture it is "
        f"{result.product.population_content_mg_per_gdcw:.8f} mg/gDCW at a bearing fraction "
        f"of {result.product.mean_bearing_fraction:.4f}, and "
        f"{assembled_gate().summary()}")
    return result


_EFLUX_CONTRACT = (
    "not wired, and the contract for wiring it is declared rather than approximated: pass a "
    "second MechanisticRun with the SAME Genotype under a declared reference Protocol, and "
    "scale each gene by its active mmol/gDW here over its value there. Two things are still "
    "missing and neither is invented here -- the gene-to-reaction association for the "
    "cassette in the installed GEM, and the reference-condition flux in that model's own "
    "units for the normalised ratio to scale")

_EMPIRICAL_ROWS_A_MECHANISTIC_CALL_SKIPS = (
    ("metabolism", "the empirical pathway solve, content = flux/mu on a fitted capacity; "
                   "pass mode='empirical' to run it. Here the pools are integrated instead "
                   "-- see 'engine metabolism'"),
    ("flux law", "the expression-to-flux scalar on genotype.entry_expression; pass "
                 "mode='empirical'. The engine reads gene COPIES and integrates the "
                 "cassette, so no fitted entry scalar is consulted"),
    ("stress panel", "empirical module activities from a stressor and a dose; pass "
                     "mode='empirical'. A mechanistic call refuses a dose outright, because "
                     "nothing maps one onto a medium composition"),
    ("environment", "the empirical context-to-growth-rate route; pass mode='empirical'. "
                    "Here the vessel is the Protocol -- see the 'engine <axis>' rows"),
    ("cassette dosage", "pathway/capacity.py's cyclase audit, refused for want of a "
                        "coefficient; pass mode='empirical'. The engine needs none, because "
                        "it integrates the gene -- see 'engine expression'"),
    ("enzyme capacity", "kcat times measured enzyme against the fitted ceiling; pass "
                        "mode='empirical' with enzyme_capacity="),
    ("GEM / FBA", "the stoichiometric audit of a PathwaySolution; pass mode='empirical' "
                  "with audit_model=. A mechanistic call produces no PathwaySolution"),
    ("regulation (E-Flux)", _EFLUX_CONTRACT),
    ("thermodynamics", "the gate on a PathwaySolution's concentrations; pass "
                       "mode='empirical' with thermo=. The engine's own gate is reached by "
                       "passing EngineParameters(thermodynamics=ThermodynamicRequest(...)) "
                       "instead -- see 'engine thermodynamics'"),
    ("latent stress state", "the maintenance-ATP branch, refused by G4; pass mode='legacy' "
                            "with latent="),
)

_MECHANISTIC_COMPANIONS_WITHOUT_A_SOLUTION = (
    ("audit_model", "audits a PathwaySolution's flux against a GSMM"),
    ("thermo", "gates a PathwaySolution's concentrations; the engine's own gate is reached "
               "by EngineParameters(thermodynamics=ThermodynamicRequest(...))"),
    ("latent", "constrains a PathwaySolution through the maintenance-ATP branch"),
    ("enzyme_capacity", "audits a fitted capacity the engine does not use"),
    ("mech", "routes an Environment through mech/chain.py to a PathwaySolution"),
)


def _drawing_on(table, truth, index, coordinate) -> tuple[str, ...]:
    """Reactions that consume one coordinate and carried extent at the read time."""
    return tuple(name for name, stoichiometry in table.items()
                 if stoichiometry.get(coordinate, 0.0) < 0
                 and float(truth[f"flux.{name}"][index]) != 0.0)


def _engine_layers(mechanism, result, index, species, table, producers, consumers,
                   production, biomass) -> list[LayerStatus]:
    """The engine's rows, each decided from THIS run's traces rather than from its name.

    Medium pH is why. At pH 4 and pH 5 the returned content is bit-identical when the run
    carries no acetate, so a row claiming pH IN_CHAIN by name would be false on that run.
    """
    truth = result.truth

    def at(name):
        return float(truth[name][index])

    def peak(name):
        return float(max(truth[name]))

    growth = at("growth_per_h")
    drivers = mechanism.genotype.promoter_drivers
    diagnostics = result.diagnostics
    genes = ", ".join(f"{gene.name}x{gene.copies}={at(f'gene.{gene.name}.active') / biomass:.3g}"
                      for gene in mechanism.genotype.genes)
    acetate = max(peak("external.acetate"), peak("internal.acetate"))
    carbon_sources = tuple(f"external.{name}" for name, (carbon, _) in EXTRACELLULAR.items()
                           if carbon > 0)
    taken_up = tuple(name for source in carbon_sources
                     for name in _drawing_on(table, truth, index, source))
    fed = tuple(sorted({name for control in mechanism.protocol.controls
                        for name in control.feed_mM} if any(
                            control.feed_l_h > 0 for control in mechanism.protocol.controls) else ()))
    breathing = _drawing_on(table, truth, index, "external.oxygen")
    request = mechanism.parameters.thermodynamics
    layers = [
        LayerStatus(
            "engine metabolism", LayerState.IN_CHAIN, True,
            f"the 109-state ODE sets the returned number: content.{species} = "
            f"internal/{biomass:.4g} gDW, filled by {list(producers)} and drained by "
            f"{list(consumers)}, integrated by {diagnostics['method']} at rtol "
            f"{diagnostics['rtol']:g} over {diagnostics['n_rhs_evaluations']} evaluations"),
        LayerStatus(
            "engine expression", LayerState.IN_CHAIN, True,
            f"per-gene mRNA -> protein -> active enzyme, in mmol/gDW at the read time: "
            f"{genes}. Gene copies enter the rate directly, so no expression-to-copy "
            f"conversion and no fitted capacity is consulted"),
        LayerStatus(
            "engine signalling", LayerState.IN_CHAIN if drivers else LayerState.REPORTED,
            bool(drivers),
            f"signals reaching a promoter on this genotype: {list(drivers)}; promoter "
            f"activity multiplies transcription, so it reaches the enzyme and the pool. "
            f"{list(PROMOTER_UNCOUPLED_SIGNALS)} reach metabolism and deliberately no "
            f"promoter -- see contracts.PROMOTER_COUPLING_NOT_BUILT"),
        LayerStatus(
            "engine vessel", LayerState.IN_CHAIN, True,
            f"{mechanism.protocol.mode}: biomass {biomass:.4g} gDW is the DENOMINATOR of "
            f"the content, at growth {growth:.4g} /h, death {at('death_per_h'):.4g} /h, "
            f"washout {at('washout_per_h'):.4g} /h"),
        LayerStatus(
            "engine temperature",
            LayerState.IN_CHAIN if production or growth else LayerState.INERT,
            bool(production or growth),
            f"control {mechanism.protocol.control_at(float(result.times_h[index])).temperature_c:g} C "
            f"scales every reaction rate through the thermal factor, so it reaches any "
            f"nonzero flux; production here is {production:.4g} mmol/gDCW/h"),
        LayerStatus(
            "engine oxygen", LayerState.IN_CHAIN if breathing else LayerState.INERT,
            bool(breathing),
            f"external oxygen is consumed by {list(breathing)} at the read time"
            if breathing else "no reaction drew on external oxygen at the read time"),
        LayerStatus(
            "engine feed", LayerState.IN_CHAIN if taken_up or fed else LayerState.INERT,
            bool(taken_up or fed),
            f"carbon is drawn from the medium by {list(taken_up)}"
            f"{f' and fed as {fed}' if fed else ''}"
            if taken_up or fed else "no carbon-bearing medium species was consumed or fed"),
        LayerStatus(
            "engine medium pH", LayerState.IN_CHAIN if acetate else LayerState.INERT,
            bool(acetate),
            f"medium pH enters through weak-acid speciation only, and this run carries "
            f"{acetate:.4g} mmol of acetate" if acetate else
            "medium pH enters through weak-acid speciation only, and this run carries no "
            "acetate at all -- it ran and provably could not have changed the number"),
        LayerStatus(
            "engine peroxide",
            LayerState.IN_CHAIN if peak("external.peroxide") else LayerState.INERT,
            bool(peak("external.peroxide")),
            f"peak external peroxide {peak('external.peroxide'):.4g} mmol"),
        LayerStatus(
            "engine osmolyte",
            LayerState.IN_CHAIN if peak("external.osmolyte") else LayerState.INERT,
            bool(peak("external.osmolyte")),
            f"peak external osmolyte {peak('external.osmolyte'):.4g} mmol"),
    ]
    if request is None:
        layers.append(LayerStatus(
            "engine thermodynamics", LayerState.NOT_RUN, False,
            "pass EngineParameters(thermodynamics=ThermodynamicRequest(...)) to gate the "
            "engine's own reactions"))
    else:
        constrains = request.mode == "constrain_supported"
        layers.append(LayerStatus(
            "engine thermodynamics",
            LayerState.IN_CHAIN if constrains else LayerState.AUDITS, constrains,
            f"mode {diagnostics['thermodynamic_mode']!r} on "
            f"{list(diagnostics['thermodynamic_constrained_reactions'])}; "
            f"{'the rate is rewritten inside the integrator' if constrains else 'a driving force can refute this run and never move it'}"))
    layers.append(LayerStatus(
        "engine conservation", LayerState.AUDITS, False,
        f"carbon residual {at('balance.carbon_residual'):.4g} mmol, nitrogen "
        f"{at('balance.nitrogen_residual'):.4g} mmol, live biomass "
        f"{at('balance.live_biomass_residual'):.4g} gDW -- it can refute the run and can "
        f"never move the number"))
    layers.append(LayerStatus(
        "engine observation",
        LayerState.NOT_RUN if mechanism.observation is None else LayerState.REPORTED, False,
        "pass observation=ObservationModel.prior() to attach simulated readouts; the "
        "content above comes from truth either way" if mechanism.observation is None else
        f"reporter {mechanism.observation.reporter!r} attached beside truth; the content "
        f"above still comes from truth"))
    layers.append(LayerStatus(
        "engine validity", LayerState.REPORTED, False,
        f"{result.validity.parameter_basis}; biological_validation="
        f"{result.validity.biological_validation}, "
        f"{len(result.prior_parameters)}/{len(result.provenance)} parameters are priors. "
        f"This number cannot enter data/current_claims.json as a supported claim, and the "
        f"returned object is not a ProductPrediction so nothing that writes one can quote it"))
    return layers


def _predict_mechanistic(spec, mechanism) -> MechanisticProductPrediction:
    """The engine's own pool, returned as the product number, with nothing converted."""
    if not isinstance(mechanism, MechanisticRun):
        raise ValueError(
            f"mode='mechanistic' requires a MechanisticRun, got {type(mechanism).__name__}")
    species = spec.product
    if species not in INTRACELLULAR:
        raise ValueError(
            f"the shared engine carries no intracellular {species!r}, so there is no pool "
            f"to read; it integrates {sorted(INTRACELLULAR)}")
    molar_mass = spec.node(species).molar_mass_g_per_mol
    if molar_mass is None:
        raise ValueError(
            f"pathway {spec.product!r} declares no molar_mass_g_per_mol for {species!r}, "
            f"and mg/gDCW without one would be an assumed conversion")

    result = mechanism.simulate()
    index = mechanism.index_of(result)
    truth = result.truth
    biomass = float(truth["biomass_gdw"][index])
    if biomass <= 0:
        raise ValueError(
            f"the run holds no live biomass at {float(result.times_h[index]):g} h, so "
            f"content per gDW has no basis; read the trajectory directly instead")

    # Producers and consumers come from the engine's OWN stoichiometry, so the rate is
    # reaction extent: death transfer and growth dilution are not reactions in that matrix.
    table = reaction_table(mechanism.parameters, mechanism.genotype)
    produced = consumed = 0.0
    producers, consumers = [], []
    for name, stoichiometry in table.items():
        coefficient = stoichiometry.get(f"internal.{species}", 0.0)
        if coefficient > 0:
            producers.append(name)
            produced += coefficient * float(truth[f"flux.{name}"][index])
        elif coefficient < 0:
            consumers.append(name)
            consumed += -coefficient * float(truth[f"flux.{name}"][index])
    production, consumption = produced / biomass, consumed / biomass

    content = float(truth[f"content.{species}"][index])
    dead_biomass = float(truth["dead_biomass_gdw"][index])
    broth = float(truth[f"internal.{species}"][index]) + float(truth[f"dead.{species}"][index])
    layers = _engine_layers(mechanism, result, index, species, table, tuple(producers),
                            tuple(consumers), production, biomass)
    layers.extend(LayerStatus(name, LayerState.NOT_RUN, False, detail)
                  for name, detail in _EMPIRICAL_ROWS_A_MECHANISTIC_CALL_SKIPS)
    return MechanisticProductPrediction(
        product=spec.product,
        species=species,
        read_time_h=float(result.times_h[index]),
        content_mmol_per_gdcw=content,
        content_mg_per_gdcw=content * molar_mass,
        production_rate_mmol_per_gdcw_h=production,
        consumption_rate_mmol_per_gdcw_h=consumption,
        net_accumulation_rate_mmol_per_gdcw_h=production - consumption,
        growth_rate_per_h=float(truth["growth_per_h"][index]),
        biomass_gdw=biomass,
        dead_biomass_gdw=dead_biomass,
        whole_broth_content_mmol_per_gdw=broth / (biomass + dead_biomass),
        result=result,
        validity=result.validity,
        prior_parameters=len(result.prior_parameters),
        total_parameters=len(result.provenance),
        layers=tuple(layers),
        notes=(
            f"MECHANISTIC: the number is the engine's own {species} pool, "
            f"truth['content.{species}'] at {float(result.times_h[index]):g} h, in the "
            f"mmol/gDW the engine already reports -- no conversion layer. "
            f"{len(result.prior_parameters)} of {len(result.provenance)} parameters behind "
            f"it are priors, so it is prior-conditional and NOT a supported claim",
            "the basis is LIVE STRUCTURAL biomass. A published titre is per GRAVIMETRIC dry "
            "weight and the ratio between the two is not measured here, so "
            f"whole_broth_content_mmol_per_gdw = {broth / (biomass + dead_biomass):.6g} is "
            "reported beside it and neither is claimed to be what an extraction assay reads "
            "-- that would also need a carotenoid recovery factor this repository does not hold",
            "genotype.entry_expression, the flux calibration and the node kinetics were "
            "shape-checked and then NOT read: the engine takes gene copies and its own rate "
            "constants, and no measured conversion between the two exists in this host",
        ),
    )


def _validate_prediction_inputs(spec, genotype, environment, flux_calibration, kinetics,
                                audit_model, audit_reaction, audit_glucose_uptake,
                                thermo, thermo_model, thermo_background_m, thermo_volumes,
                                refuse_thermodynamically_blocked, latent, legacy):
    if audit_model is None and (audit_reaction is not None or audit_glucose_uptake is not None):
        raise ValueError("audit_reaction and audit_glucose_uptake require audit_model")
    if audit_model is not None and (audit_reaction is None or audit_glucose_uptake is None):
        raise ValueError("audit_model requires audit_reaction and audit_glucose_uptake")
    if audit_glucose_uptake is not None and (
            not math.isfinite(audit_glucose_uptake) or audit_glucose_uptake <= 0):
        raise ValueError("audit_glucose_uptake must be finite and positive")
    if thermo is None and (thermo_model is not None or thermo_background_m is not None
                           or thermo_volumes is not None
                           or refuse_thermodynamically_blocked is not None):
        raise ValueError("thermodynamic inputs require thermo")
    if refuse_thermodynamically_blocked is not None:
        if not isinstance(refuse_thermodynamically_blocked, bool):
            raise ValueError("refuse_thermodynamically_blocked must be a boolean or None")
        if not refuse_thermodynamically_blocked and not legacy:
            raise ValueError("refuse_thermodynamically_blocked=False requires mode='legacy'; empirical predictions cannot ignore a blocked step")
    if thermo is not None and not spec.thermo_steps:
        raise ValueError(f"thermo cannot audit pathway {spec.product!r}: it declares no step chemistry")
    if thermo_volumes is not None and (
            not thermo_volumes or any(not math.isfinite(v) or v <= 0 for v in thermo_volumes)):
        raise ValueError("thermo_volumes must contain finite positive volumes")
    if not math.isfinite(flux_calibration.alpha) or flux_calibration.alpha <= 0:
        raise ValueError("flux calibration alpha must be finite and positive")
    if not math.isfinite(flux_calibration.loso_rmse_log) or flux_calibration.loso_rmse_log < 0:
        raise ValueError("flux calibration loso_rmse_log must be finite and non-negative")
    unknown = set(kinetics) - {node.name for node in spec.nodes}
    if unknown:
        raise ValueError(f"kinetics supplied for absent nodes: {sorted(unknown)}")
    for name, parameters in kinetics.items():
        node = spec.node(name)
        used = {"growth_rate_range"}
        if node is not spec.nodes[-1]:
            if node.rate_law == RateLaw.SATURATING:
                used.update(("vmax_per_growth", "km"))
            elif node.rate_law == RateLaw.PROPORTIONAL:
                used.add("rate_constant")
        if node.fate == Fate.DEGRADED:
            used.add("degradation_rate_per_h")
        elif node.fate == Fate.SECRETED:
            used.update(("secretion_vmax_mmol_per_gdcw_h", "secretion_km_mmol_per_gdcw"))
        for parameter in fields(NodeKinetics):
            value = getattr(parameters, parameter.name)
            if value is None:
                continue
            if parameter.name not in used:
                raise ValueError(f"kinetics {name}.{parameter.name} is not used by this node's declared rate law and fate")
            values = value if parameter.name == "growth_rate_range" else (value,)
            if any(not math.isfinite(v) for v in values):
                raise ValueError(f"kinetics {name}.{parameter.name} must be finite")
    if not legacy:
        if latent is not None:
            raise ValueError("latent input has no validated route to product prediction; mode='legacy' returns diagnostics only")
        defaults = CultureContext()
        for name in ("glucose_g_per_L", "ph_medium", "ph_cytosolic"):
            if getattr(environment.context, name) != getattr(defaults, name):
                raise ValueError(f"context.{name} is not modelled by the empirical product path; mode='legacy' returns diagnostics only")
        if genotype.cassette:
            enzyme = CAPACITY_SETTING_ENZYME.get(spec.product)
            unused = set(genotype.cassette) - {enzyme}
            if unused:
                raise ValueError(f"cassette dosage has no calibrated route for {sorted(unused)}; entry_expression is measured expression, not gene copy number")
            if content_ceiling(spec, kinetics) is None:
                raise CapacityUnmeasured(f"cassette dosage for {spec.product!r} has no measured capacity to audit")


def predict_product(
    spec: PathwaySpec,
    genotype: Genotype,
    environment: Environment,
    flux_calibration: FluxCalibration,
    kinetics: dict[str, NodeKinetics],
    audit_model=None,
    audit_reaction: str | None = None,
    audit_glucose_uptake: float | None = None,
    thermo=None,
    thermo_model=None,
    thermo_background_m: dict[str, float] | None = None,
    thermo_volumes: tuple[float, ...] | None = None,
    refuse_thermodynamically_blocked: bool | None = None,
    latent=None,
    enzyme_capacity: dict | None = None,
    mech=None,
    *,
    mode: str = "empirical",
    mechanism: MechanisticRun | None = None,
) -> ProductPrediction | LegacyProductDiagnostic | MechanisticProductPrediction:
    """Run the chain for one strain in one environment.

    Args:
        spec: The declared pathway. Loaded from ``data/pathways``; nothing here is
            product-specific.
        genotype: The strain's entry-enzyme expression.
        environment: The condition to predict for.
        flux_calibration: Fitted expression-to-flux scalar for this pathway and host.
        kinetics: Parameters for every node whose rate law needs them.
        audit_model: An ALREADY-LOADED cobra model carrying the installed pathway. Given
            one, the GSMM is asked whether the predicted flux fits, what it costs in growth,
            and what the native competition on the precursor node looks like. The audit
            cannot change the prediction -- an upper bound cannot make a flux mandatory, and
            capping every pathway reaction at the measured magnitude leaves the FVA floor at
            exactly zero -- so it is a check on a number rather than a source of one.

            A model rather than a path, deliberately: the caller pays the 18-second SBML
            parse once and reuses it, instead of this function paying it per call.
        audit_reaction: The terminal production reaction id in that model, using one
            product molecule per reaction flux unit. Required with ``audit_model``.
            The audit receives ``solution.terminal.flux_in``, not the total entry
            flux that also supplies intermediate dilution/losses. It checks this
            single flux at its own growth-fraction regime, not the whole solved
            pathway jointly at the supplied setpoint.
        audit_glucose_uptake: Uptake bound the audit runs under. Required with
            ``audit_model``, and required rather than defaulted because the growth cost swings
            sevenfold with it -- the same pathway at the same flux costs 0.09% at a bound of
            10 and 0.64% at 1.5. A burden number without its regime is not a number.
        thermo: A :class:`~ystwin.bridge.thermodynamic.ThermodynamicData`. Given one, every
            step the spec declares chemistry for is gated at the concentrations THIS SOLVE
            produced, across ``thermo_volumes``. Like the GEM audit it cannot move the
            number: a driving force says whether a step may run, never how fast, so the
            only thing it can add is a contradiction -- carbon moved through a step that
            points the other way at the concentrations the same solve returned.

            Loading it costs about 1.5 s and reaches disk, which is why it is a parameter
            rather than something this function fetches. ``ThermodynamicData.load()``.
        thermo_model: Cobra model used to resolve any step the spec declares as a bare
            reaction id. Defaults to ``audit_model``. Both shipped specs declare
            stoichiometries as well, so the gate runs with no model at all -- and when a
            model IS given the gate refuses a pair whose stoichiometry and reaction have
            drifted apart.
        thermo_background_m: Molar concentrations for participants the solve does not pin --
            cofactors, and the precursor the pathway is pulled off. Consulted only there:
            a background value for a node that solved to a positive content is refused
            rather than silently preferred. Without it most steps come back ``cannot_say``,
            which is the correct answer to "is this step feasible" when nobody has said
            what the cell contains.
        thermo_volumes: Cytosolic volumes, mL/gDCW, to sweep. The mmol/gDCW-to-molar
            conversion needs one and nobody here has measured it, so the gate runs at
            several and :func:`~ystwin.pathway.thermo_gate.volume_dependent_steps` names
            the verdicts that did not survive the choice. A verdict that flips across this
            range is an artefact of the convention and is reported as one.
        refuse_thermodynamically_blocked: Raise rather than report when a gated step cannot
            run at the solved concentrations. Default ``False``, because the gate's coverage
            is 44.5% of yeast-GEM's reactions and a refusal resting on a table that mostly
            says ``cannot_say`` would fire unevenly. Set it where a blocked step should
            stop a pipeline rather than annotate one.
        enzyme_capacity: Optional ``{"enzyme", "kcat_per_s", "enzyme_mmol_per_gdcw"}``. The
            content is in mmol/gDCW, and the two conversions into it live in
            :mod:`ystwin.pathway.proteome`, one per convention -- ``enzyme_content_mmol_per_gdcw``
            for a molar ppm and ``enzyme_content_from_mass_fraction`` for a PRM mass
            fraction. This dict took ``proteome_fraction``/``molar_mass_g_per_mol`` until
            2026-09-04, when the fraction's convention turned out to be unstated on both
            sides of the seam. Given one, the fitted ceiling is audited against
            ``vmax = kcat * [E]`` -- a capacity computed from properties of the PROTEIN
            rather than fitted to this product's own states. Like the GEM audit and the
            thermodynamic gate it cannot move the number, and for the same class of reason:
            it is scored against nothing, while the fitted capacity is scored
            leave-one-strain-out against six measured states.

            It exists because it would have caught a live error before the literature did.
            The shipped ceiling is 1.2483 mg/gDCW and Arhar 2024 measured 79; at any
            turnover inside the host model's own isoprenoid range this audit lands on the
            measurement and misses the fit by orders of magnitude. See
            :mod:`ystwin.pathway.enzyme_capacity` for what would promote it.
        latent: A :class:`~ystwin.bridge.latent_bridge.LatentState` -- a dilution-corrected
            promoter activity from the estimator. Given one, the maintenance-ATP branch is
            computed and attached, and **it is attached rather than used**: G4 refused the
            anchor that would validate it, and the branch would be the only route by which
            stress could reach the product number. Feeding it in would launder an
            unvalidated branch into a headline number, so it comes back as
            ``latent_constraints`` under ``LayerState.REPORTED`` with the G4 verdict on it.
            The number is bit-for-bit identical with and without it, and
            ``tests/test_everything_is_wired.py`` pins that.
        mech: A :class:`~ystwin.mech.chain.SweepPoint`. Given one, the environment is routed
            through `mech/chain.py` -- weak-acid speciation, the two pH states, the Pma1
            proton bill, the cassette burden, the vessel and the two-population split -- and
            every layer it reports is copied into ``layers`` under a ``mech `` prefix, with
            the whole run on ``mech_chain``.

            **DEFAULT OFF, and the default is the contract.** Every number this package
            quotes leave-one-strain-out was measured with the mechanism absent, so turning
            it on silently would move them without anyone deciding to. ``mech=None`` is
            bit-for-bit the incumbent, layer report included.

            **What it buys, measured rather than asserted.** One thing: ``mu = min(mu_set,
            mu_max)`` on BOTH branches. `_growth_rate` returns a setpoint unchecked when no
            stressor is present, which is why seven environments returned one content to
            eight decimal places -- not because the environment has no route but because the
            route was checked on one branch only. With ``mech`` on at mu_set = 0.18 /h the
            ethanol fed-batch is solved at its own 0.140 /h and returns 1.11791951 against
            1.07098788 mg/gDCW, a 1.0438x move, while 37 C, pH 6, near-anoxia, 1 mM H2O2 and
            60 mM acetic acid at pH 4.5 all still return 1.07098788 to eight decimals,
            because their mu_max stays above the setpoint. That is the honest size of it.

            **The acid arm cannot reach the content at any admissible setpoint.** 60 mM
            acetic acid at pH 4.5 takes mu_max from 0.400 to 0.344 /h, and the kinetics are
            fitted over [0.101, 0.2543] /h -- so no setpoint `pathway/solve.py` will accept
            is low enough for the clip to bite. Its anion pool is flagged unphysical by
            `mech/ph.py` on top of that.

            A ``SweepPoint`` has no default corner, deliberately: `mech/chain.py`'s
            assembled free-scalar gate is REFUSED at 7 free scalars against 3 independent
            targets, and a number from a refused assembly must at least carry the corner it
            came from. The refusal is attached as the ``mech free scalars`` layer. Callers
            wanting the reporter arm, another window or a declared copy number call
            :func:`~ystwin.mech.chain.run_chain` directly; this argument is the product
            channel only.
        mechanism: A :class:`MechanisticRun`, and required by ``mode="mechanistic"``. The
            109-state engine is integrated and the returned content IS its own
            ``truth["content.<product>"]``, unconverted -- so temperature, medium, feed and
            gene COPIES move the number, which the empirical route cannot express at a held
            setpoint. A :class:`MechanisticProductPrediction` comes back rather than a
            :class:`ProductPrediction`, because 160 of the 173 parameters behind that number
            are priors: the DATAFLOW claim is strong and the EVIDENCE claim is absent, and
            the type keeps the second from being borrowed from the first.

            **``mode="empirical"`` remains the default and this argument changes nothing
            about it.** Passing a mechanism to any other mode is refused rather than
            silently switching, and ``mechanism=None`` is bit-for-bit the incumbent call.

    Raises:
        SetpointUnreachable: when the feed cannot hold this setpoint at this dose. **Not
            raised with ``mech`` on**: the mechanism reports the realised rate on both
            branches instead, because refusing on the stressor branch while returning an
            unchecked setpoint on the other is not one rule.
        ThermodynamicallyBlocked: with ``refuse_thermodynamically_blocked``, when a gated
            step cannot run at the concentrations this solve produced.
        ValueError: from the layers, when the condition is outside what they measured.
    """
    if mode not in {"empirical", "legacy", "mechanistic"}:
        raise ValueError("mode must be 'empirical', 'legacy' or 'mechanistic'")
    if mechanism is not None and mode != "mechanistic":
        raise ValueError(
            f"a MechanisticRun was supplied with mode={mode!r}. It is refused rather than "
            "run, because switching modes on the presence of an argument is how a figure "
            "moves without anyone deciding to; pass mode='mechanistic' to use it")
    if mode == "mechanistic":
        if mechanism is None:
            raise ValueError(
                "mode='mechanistic' needs a MechanisticRun: the engine takes a Protocol, "
                "its own Genotype, a PhysicalState inoculum and EngineParameters, and no "
                "default protocol, medium or inoculum is invented here")
        supplied = {"audit_model": audit_model, "thermo": thermo, "latent": latent,
                    "enzyme_capacity": enzyme_capacity, "mech": mech}
        for name, role in _MECHANISTIC_COMPANIONS_WITHOUT_A_SOLUTION:
            if supplied[name] is not None:
                raise ValueError(
                    f"{name} {role}, and mode='mechanistic' produces no PathwaySolution for "
                    f"it to act on. It is refused by name rather than ignored, which would "
                    f"be a silent no-op")
        bare = Environment()
        offending = [f.name for f in fields(Environment)
                     if getattr(environment, f.name) != getattr(bare, f.name)]
        if offending:
            raise ValueError(
                f"mode='mechanistic' refuses the environment fields {offending}: the "
                "Protocol already states the culture, and nothing reconciles two claims "
                "about one vessel. No calibration maps a stressor dose in its own units "
                "onto Control.temperature_c, ph or feed_mM either. Pass Environment() and "
                "declare the condition in the MechanisticRun's Protocol")
    legacy = mode == "legacy"
    _validate_prediction_inputs(
        spec, genotype, environment, flux_calibration, kinetics, audit_model, audit_reaction,
        audit_glucose_uptake, thermo, thermo_model, thermo_background_m, thermo_volumes,
        refuse_thermodynamically_blocked, latent, legacy)
    thermo_volumes = CYTOSOLIC_VOLUMES_ML_PER_GDCW if thermo_volumes is None else thermo_volumes
    if flux_calibration.entry_enzyme != spec.entry_enzyme:
        raise ValueError(
            f"pathway {spec.product!r} reads {spec.entry_enzyme!r} but the calibration was "
            f"fitted on {flux_calibration.entry_enzyme!r}. Applying one gene's scalar to "
            "another gene's expression is a units error that would not otherwise surface")

    # The empirical inputs above were shape-checked and are not read below: the engine takes
    # gene copies and its own rate constants, and no measured conversion joins the two.
    if mode == "mechanistic":
        return _predict_mechanistic(spec, mechanism)

    notes: list[str] = []
    # The mechanistic run happens FIRST when it was asked for, because the one thing it
    # changes is the growth rate every line below is solved at.
    mech_rows: list[LayerStatus] = []
    chain = None
    if mech is None:
        growth = _growth_rate(environment, notes)
    else:
        if not legacy:
            from .mech.chain import require_assembled_gate

            require_assembled_gate()
        chain = _mech_layers(environment, mech, spec, genotype, flux_calibration, kinetics,
                             mech_rows, notes)
        if not legacy and not chain.bill_is_physical:
            raise ValueError("the mechanistic proton bill is not physical; no product prediction is supported")
        if (not legacy and environment.growth_rate_setpoint_per_h is not None
                and not chain.feed_sets_mu):
            raise SetpointUnreachable("the mechanism cannot hold the supplied growth setpoint")
        growth = chain.growth_rate_per_h

    modules = ({} if environment.stressor is None
               else module_response(environment.stressor, environment.dose))
    feed_holds = chain is None or chain.feed_sets_mu
    if (environment.stressor is not None
            and environment.growth_rate_setpoint_per_h is not None and feed_holds):
        # SAY THIS OUT LOUD. In a chemostat the pump fixes mu, and mu was the only channel
        # stress had to the product, so the number returned here is bit-for-bit the
        # unstressed one. That is the correct arithmetic and it is NOT the claim
        # "stress does not affect production" -- it is "this model has no route by which it
        # could, and does not know whether one exists".
        #
        # Two routes are plausible and neither is modelled. A maintenance-ATP burden would
        # divert carbon from the product: `bridge/latent_bridge.py` implements exactly that
        # and is REFUSED by G4, so wiring it in would launder an unvalidated branch into a
        # headline number. And regulation could redirect flux, but E-Flux scales upper
        # bounds, and an upper bound cannot lower a flux that is not at its ceiling.
        #
        # Silence here would be the worst option: a user comparing a stressed and an
        # unstressed prediction would read the identical numbers as a finding.
        notes.append(
            f"{environment.stressor} at {environment.dose:g} reaches the growth rate and "
            "nothing else, and a fed-batch feed fixes that -- so this product "
            "number is identical to the unstressed one BY CONSTRUCTION. Read it as 'no "
            "route modelled', not as 'no effect'. docs/CLAIM_BOUNDARY.md")

    flux = predict_flux_from_expression(genotype.entry_expression, flux_calibration)
    notes.extend(flux.notes)

    physiology = chemostat_physiology(growth)
    # THE PHYSIOLOGY STAMPED ON THIS RESULT DESCRIBES A DIFFERENT CULTURE FROM THE ONE THE
    # KINETICS WERE FITTED ON, and always will until a second table exists.
    #
    # `chemostat_physiology` reads van Hoek 1998's GLUCOSE-LIMITED chemostat. Elizondo's six
    # calibration states are glucose-EXCESS: measured glucose uptake is 4.1-6.4x the table's
    # value at the same mu, and they make 6.0-16.1 mmol/gDCW/h of ethanol at EVERY rate
    # including mu = 0.101, where a glucose-limited culture makes none.
    #
    # So `fermentative` -- derived from mu alone against a critical rate of 0.28 -- reports
    # False for every calibration state while those cultures were visibly fermenting. The
    # flag is kept because `generator/design.py` and the environment sweep consume it and it
    # is right for the cultures van Hoek measured; it is flagged here because it is wrong
    # for the cultures this pathway was fitted to.
    notes.append(
        "physiology (glucose uptake, ethanol, the fermentative flag) is read from van "
        "Hoek 1998's GLUCOSE-LIMITED table, and the states this pathway was calibrated on "
        "were glucose-EXCESS -- 4.1-6.4x that uptake, and fermenting at every rate. Treat "
        "the physiology as context for a different culture, not as a property of this one")

    # The chain solved this already, with the same spec, flux and kinetics -- reusing it
    # rather than re-solving is what keeps the two from being able to disagree.
    solution = (chain.product.solution if chain is not None
                else solve_pathway(spec, flux.flux_mmol_per_gdcw_h, growth, kinetics))

    if physiology.growth_rate > CRITICAL_GROWTH_RATE_PER_H:
        notes.append(
            "past the Crabtree threshold, so carbon is going to ethanol; the flux "
            "calibration states were all below it and none saw this regime")
    if environment.growth_rate_setpoint_per_h is not None and feed_holds:
        # The context is an ACCEPTED INPUT that reaches nothing. Under a setpoint the feed
        # fixes mu, and mu is the only route the context has, so carbon source, temperature,
        # oxygen and pH are all read and discarded -- two calls differing only in
        # `carbon_source` return bit-identical numbers AND identical notes.
        #
        # Kocharin's PHB chemostats measure that discard directly: at D = 0.05 the feed
        # alone moves the product flux 4.24x with mu held identical. So this is not a
        # rounding-level omission, it is the largest single effect in the only environment
        # series this repository has.
        notes.append(
            f"context ({environment.context.carbon_source}, "
            f"{environment.context.temperature_c:g} C, O2 {environment.context.oxygen:g}, "
            f"pH {environment.context.ph_medium:g}) AUDITS REACHABILITY of the held "
            "growth rate. The empirical product law uses the realised mu, not a direct "
            "carbon-source response. Kocharin's one-genotype series moves flux 4.24x at "
            "identical mu, so that direct response remains unmodelled -- see "
            "docs/EXTERNAL_PRODUCT_VALIDATION.md section 3")

    ceiling = content_ceiling(spec, kinetics)
    if ceiling is not None and solution.terminal.content_mmol_per_gdcw > 0.9 * ceiling:
        # Near an asymptote the prediction stops being about the genotype. Worth saying,
        # because the number still looks like an ordinary answer and a caller comparing two
        # strains here would read a difference of 0.4% as a real ranking.
        molar = spec.nodes[-1].molar_mass_g_per_mol or 1.0
        notes.append(
            f"within 10% of the ceiling THIS CALIBRATION imposes, {ceiling * molar:.3g} "
            "mg/gDCW. mu cancels out of content <= vmax_per_growth, so no genotype at any "
            "growth rate can exceed it WITHIN THIS MODEL -- but that is a statement about "
            "the model and NOT about the organism: Arhar 2024 (PMID 39215465) measured 79 "
            "mg/gDCW of beta-carotene by HPLC on gravimetric dry weight, about 63x this "
            "bound, by retargeting CrtYB and CrtI. Capacity was fitted flat across strains "
            "differing only in the ENTRY enzyme, and the cyclase dosage that actually sets "
            "it is not a field Genotype carries. Predictions here are about the fitted "
            "capacity rather than about the strain, and the capacity is known to be wrong "
            "as an upper bound -- see docs/HARD_TESTS.md")
    if not spec.calibratable:
        notes.append(
            f"pathway {spec.product!r} has fewer than two measurable nodes, so its flux "
            "scalar cannot be identified from its own measurements")

    # ------------------------------------------------------------------ the five layers
    #
    # Everything from here down runs beside the number rather than into it, and the
    # `layers` record is what makes that legible without reading this file. Two of them
    # take a caller-supplied object because the object is expensive -- an 18-second SBML
    # parse, a 1.5-second table load -- and paying that per call to learn nothing new would
    # be a poor trade. `NOT_RUN` says which argument would have run it.
    layers: list[LayerStatus] = _core_layers(
        solution, flux, genotype, environment, modules, spec, ceiling, enzyme_capacity,
        chain, legacy)

    verdict = _gem_audit(spec, solution, modules, audit_model, audit_reaction,
                         audit_glucose_uptake, layers, notes, legacy)

    reports = _thermodynamic_gate(
        spec, solution, thermo, thermo_model, audit_model, thermo_background_m,
        thermo_volumes, refuse_thermodynamically_blocked or not legacy, layers, notes)

    constraints = _latent_branch(latent, growth, physiology, layers, notes)
    # Last, so a bare call's report is unchanged row for row and the mechanism reads as the
    # addition it is.
    layers.extend(mech_rows)

    prediction = ProductPrediction(
        product=spec.product,
        mode=mode,
        genotype=genotype,
        environment=environment,
        solution=solution,
        flux=flux,
        growth_rate_per_h=growth,
        fermentative=physiology.fermentative,
        # abs(), because `> 0.01` was one-sided and dropped every module the panel scored
        # NEGATIVE -- acetic acid's ph = -0.696 at 60 mM among them. See `_stress_panel_layer`.
        stress_modules={k: v for k, v in modules.items() if abs(v) > 0.01},
        notes=tuple(notes),
        layers=tuple(layers),
        thermodynamics=reports,
        flux_audit=verdict,
        latent_constraints=constraints,
        mech_chain=chain,
    )
    return LegacyProductDiagnostic(prediction) if legacy else prediction
