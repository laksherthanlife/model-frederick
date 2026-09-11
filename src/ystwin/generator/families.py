"""The same panel with its structure moved, rather than its constants redrawn.

The older panel workflow starts from `stress_panel.py`: one crosstalk topology, one set of
EC50s, one loading matrix. These families challenge that structure. Muratore, Gienger &
Peters (IEEE TPAMI 43(4):1172-1183, 2021) show that a globally optimized finite-sample
objective is optimistic IN EXPECTATION relative to the population optimum under an assumed
distribution. This is not a pointwise guarantee that held-out simulator scores exceed real
biology. Mak, Morton & Wood (Oper. Res. Lett. 24(1-2):47-56, 1999, Theorem 2) give the
expected sample-optimum monotonicity result under iid sampling; a hand-built registry alone
does not establish those assumptions or bound the biological reality gap.

What a family deliberately is not is a re-draw. Perturbing the same topology's constants
asks whether the fitter tolerates noise in numbers it was handed; the question that decides
whether the transfer claim is about yeast or about this file is whether it tolerates the
topology being wrong. So each family moves one piece of structure -- an edge, the shape of a
falling limb, whether a regulon stays on -- and each states what it moves and why that
alternative is live. A family nobody would defend makes robustness cheap and meaningless.

`loading_noise` is the exception and earns its place by being the control the others are read
against: if transfer moves as far under a re-draw of unmeasured constants as under a
rewiring, then the structural families have shown nothing about structure.

Two design choices are worth stating because they are what make the numbers comparable.
Propagation is a linear solve, `(I - D W) a = D d`, rather than the sequential pass in
`stress_panel._propagate`. On the literature topology the two agree exactly -- the recorded
cascade is acyclic and `_CASCADE_ORDER` is a topological order of it -- and the solve is what
lets a family close a loop, which a single pass in a fixed order cannot represent. And the
per-module dose response is `panel_calibration._model`, the same closed form that module fits
to real ladders, imported rather than retyped so a family cannot drift away from the fitted
shape. With every scale at one it reduces to `module_response`, so the baseline family
reproduces `panel_dataset` to floating point.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from .context import CultureContext, baseline_activity, context_growth_rate
from .panel_calibration import _model as _biphasic
from .panel_experiment import _DEFAULT_RUNGS, PanelDataset
from .plate_layout import well_positions
from .stress_panel import (
    MODULES,
    REPORTERS,
    STRESSORS,
    _POOLS,
    growth_rate,
    healthy_ladder,
    module_ec50,
)

__all__ = [
    "FAMILY_BUILDERS",
    "PanelFamily",
    "baseline_family",
    "build_family",
    "family_dataset",
    "family_names",
    "module_activity",
    "sample_family",
]

_DROPPED_EDGE = ("proteasome", "heat")
"""Which cascade edge `edge_dropped` removes: Hsf1 driving RPN4, weight 0.40."""

_ADDED_EDGE = ("oxidative", "ESR", 0.25)
"""Which edge `edge_added` introduces, and at what weight: Msn2/4 driving YAP1."""

_FEEDBACK_EDGE = ("heat", "proteasome", -0.25)
"""The loop `feedback` closes: restored proteostasis capacity relieving Hsf1."""

_LIMB_SCALES = {"ESR": 1.5, "proteasome": 1.3} | {pool: 0.7 for pool in _POOLS}
"""Per-module multipliers on the falling limb for the `biphasic` family.

The direction of each is sourced and the magnitude is asserted; see the family's own
`plausible_because` for which is which.
"""

_ESR_SURVIVING_FRACTION = 0.35
"""Fraction of the general stress response still present when the plate is read."""

_LOADING_BRACKET = 2.0
"""Factor within which 95% of `loading_noise` draws of a crosstalk weight fall."""

_MAX_LOOP_GAIN = 0.95
"""Largest spectral radius of the cascade a family may carry.

At one the loop is self-sustaining and the linear solve is singular: activity would be
determined by the network rather than by the dose, which is not a stress response.
"""


@dataclass(frozen=True)
class PanelFamily:
    """One generator configuration, and the case that it is a live alternative.

    `perturbs` and `plausible_because` are not documentation of the code, they are the
    scientific content: a robustness claim is only worth as much as the plausibility of the
    configurations it survived, so a family that cannot say why its perturbation might be
    the true one has made the claim easier without making it stronger. Both are required.

    Args:
        name: Registry key.
        perturbs: What this family changes relative to the literature topology.
        plausible_because: Why that alternative is live, and which parts of the magnitude
            are sourced against which are asserted.
        drivers: Cascade edges as ``module -> {driver: weight}``, replacing every
            ``Module.driven_by`` at once. Cycles are permitted, which is the point of
            carrying the whole matrix rather than a diff.
        also_reads: Crosstalk loadings as ``reporter -> {module: weight}``. The reporter's
            own element is implicit at 1.0 and is not represented here, because it is a
            normalisation of the channel's units rather than a measurement.
        lethal_scale: Multipliers on each module's falling-limb dose. Empty means the
            shared viability factor of `stress_panel.module_response`.
        adaptation: Fraction of each module's activity still present when the plate is
            read. Empty means every regulon is sustained, which is what the panel assumes.
        draw: Seed behind any stochastic part of this configuration, so a sampled family
            is reproducible from its identity alone.
    """

    name: str
    perturbs: str
    plausible_because: str
    drivers: Mapping[str, Mapping[str, float]] = field(default_factory=dict)
    also_reads: Mapping[str, Mapping[str, float]] = field(default_factory=dict)
    lethal_scale: Mapping[str, float] = field(default_factory=dict)
    adaptation: Mapping[str, float] = field(default_factory=dict)
    draw: int = 0

    def __post_init__(self) -> None:
        if not self.perturbs or not self.plausible_because:
            raise ValueError(
                f"family {self.name!r} must say what it perturbs and why that is plausible; "
                "an undefended family makes a robustness claim easier without strengthening it"
            )
        unknown = sorted(
            {m for m in self.drivers if m not in MODULES}
            | {d for edges in self.drivers.values() for d in edges if d not in MODULES}
            | {m for m in self.lethal_scale if m not in MODULES}
            | {m for m in self.adaptation if m not in MODULES}
        )
        if unknown:
            raise KeyError(f"family {self.name!r} names modules that do not exist: {unknown}")
        absent = sorted(r for r in self.also_reads if r not in REPORTERS)
        if absent:
            raise KeyError(f"family {self.name!r} names reporters that do not exist: {absent}")
        bad_scale = sorted(m for m, s in self.lethal_scale.items() if not s > 0.0)
        if bad_scale:
            raise ValueError(
                f"family {self.name!r} gives a non-positive falling-limb scale to {bad_scale}; "
                "a module with no lethal dose never turns over"
            )
        bad_adapt = sorted(m for m, a in self.adaptation.items() if not 0.0 < a <= 1.0)
        if bad_adapt:
            raise ValueError(
                f"family {self.name!r} gives {bad_adapt} an adaptation outside (0, 1]; a "
                "regulon that adapts to more than its sustained level is an induction, not "
                "an adaptation, and belongs in the weights"
            )
        gain = float(np.max(np.abs(np.linalg.eigvals(self._closed_loop()))))
        if gain >= _MAX_LOOP_GAIN:
            raise ValueError(
                f"family {self.name!r} has cascade loop gain {gain:.3f}, at or above the "
                f"{_MAX_LOOP_GAIN} ceiling: its activity would be set by the network rather "
                "than by the dose"
            )

    def cascade_matrix(self) -> np.ndarray:
        """Cascade weights as ``(modules, modules)``, entry ``[i, j]`` for driver ``j``."""
        order = list(MODULES)
        index = {name: i for i, name in enumerate(order)}
        matrix = np.zeros((len(order), len(order)))
        for module, edges in self.drivers.items():
            for driver, weight in edges.items():
                matrix[index[module], index[driver]] = float(weight)
        return matrix

    def adaptation_vector(self) -> np.ndarray:
        """Surviving fraction per module, in module order; one where nothing adapts."""
        return np.array([float(self.adaptation.get(name, 1.0)) for name in MODULES])

    def _closed_loop(self) -> np.ndarray:
        return self.adaptation_vector()[:, None] * self.cascade_matrix()

    def loadings(self, reporters: Sequence[str]) -> np.ndarray:
        """How each reporter loads on every module, from this family's promoter content.

        The counterpart of `stress_panel.reporter_loadings`, reading the family's crosstalk
        rather than the panel's. Only promoter content appears, for the same reason: a
        module driving another is already carried by the cascade, and reading it again here
        would count one piece of regulation twice.
        """
        order = list(MODULES)
        index = {name: i for i, name in enumerate(order)}
        matrix = np.zeros((len(reporters), len(order)))
        for row, name in enumerate(reporters):
            if name not in REPORTERS:
                raise KeyError(f"no reporter {name!r}; have {sorted(REPORTERS)}")
            matrix[row, index[REPORTERS[name].module]] = 1.0
            for module, weight in self.also_reads.get(name, {}).items():
                if module not in index:
                    raise KeyError(f"reporter {name!r} reads no such module {module!r}")
                matrix[row, index[module]] += float(weight)
        return matrix


def module_activity(family: PanelFamily, doses: Mapping[str, float]) -> dict[str, float]:
    """Activation of every module under one or more agents, under this family's structure.

    The family analogue of `stress_panel.combination_response`, and the same three steps:
    each agent's direct action on its own targets, the cascade carrying those into the
    regulons they drive, and adaptation removing whatever has decayed by the time the plate
    is read. Direct contributions are summed before the cascade rather than after, which is
    identical because the cascade is linear, and cheaper because it solves once.
    """
    order = list(MODULES)
    direct = np.zeros(len(order))
    index = {name: i for i, name in enumerate(order)}
    unknown = sorted(name for name in doses if name not in STRESSORS)
    if unknown:
        raise KeyError(f"no stressor(s) {unknown}; have {sorted(STRESSORS)}")
    for stressor, dose in doses.items():
        spec = STRESSORS[stressor]
        applied = max(float(dose), 0.0)
        if applied <= 0.0:
            continue
        for module, weight in spec.targets.items():
            if weight == 0.0:
                continue
            limb = spec.lethal_dose * family.lethal_scale.get(module, 1.0)
            direct[index[module]] += float(
                _biphasic(applied, module_ec50(stressor, module), limb, weight, 0.0))
    surviving = family.adaptation_vector()
    system = np.eye(len(order)) - surviving[:, None] * family.cascade_matrix()
    return dict(zip(order, np.linalg.solve(system, surviving * direct)))


def family_dataset(
    family: PanelFamily,
    reporters: Iterable[str] | None = None,
    stressors: Iterable[str] | None = None,
    doses: Sequence[float] | None = None,
    replicates: int = 3,
    noise_cv: float = 0.05,
    growth_rate_se: float = 0.0,
    seed: int = 0,
    combinations: Iterable[Sequence[str]] = (),
) -> PanelDataset:
    """Simulate a reporter panel from one family, in the shape `panel_dataset` returns.

    Deliberately the same arithmetic in the same order as
    `panel_experiment.panel_dataset`, down to the order the noise is drawn in, so the
    baseline family's readings are that function's readings and every per-family score is
    comparable with the number already in `docs/NULL_RESULTS.md`. What it does not carry are
    the axes that function already varies -- read time, plate position, culture context,
    reporter turnover -- because those are orthogonal to structure, and a family that moved
    them too would confound "the topology was wrong" with "the plate was read early".

    Growth is not perturbed by any family. The lethal doses behind it are the panel's one
    block of measured dose parameters, read off optical density rather than inferred from a
    reporter, so a family that moved them would be perturbing a measurement.

    Args:
        family: Configuration to simulate.
        reporters: Reporter names to read; defaults to the whole panel.
        stressors: Agents to apply; defaults to all of them.
        doses: Multiples of each agent's EC50. Left unset, each agent gets its own healthy
            ladder, as in `panel_dataset`.
        replicates: Independent wells per treatment and dose.
        noise_cv: Multiplicative measurement noise on a reading.
        growth_rate_se: Standard error of the estimated growth rate, in 1/h, propagated
            through the dilution correction.
        seed: Seed for the noise.
        combinations: Agent tuples to co-apply.
    """
    reporters = list(reporters) if reporters is not None else list(REPORTERS)
    stressors = list(stressors) if stressors is not None else list(STRESSORS)
    unknown = sorted(s for s in stressors if s not in STRESSORS)
    if unknown:
        raise KeyError(f"no such stressor(s): {unknown}")
    loadings = family.loadings(reporters)
    basal = np.array([REPORTERS[r].basal for r in reporters])
    order = list(MODULES)
    rng = np.random.default_rng(seed)

    treatments = [(name, (name,)) for name in stressors]
    for group in combinations:
        members = tuple(group)
        missing = sorted(m for m in members if m not in STRESSORS)
        if missing:
            raise KeyError(f"no stressor(s) {missing}; have {sorted(STRESSORS)}")
        treatments.append(("+".join(members), members))

    context = CultureContext()
    floor = baseline_activity(context)
    unstressed = context_growth_rate(context)
    rows, labels, applied, truth, growths = [], [], [], [], []
    dropped = 0
    for label, members in treatments:
        rungs = (list(doses) if doses is not None
                 else list(healthy_ladder(members[0], _DEFAULT_RUNGS)
                           / STRESSORS[members[0]].ec50))
        for multiple in rungs:
            applied_doses = {m: multiple * STRESSORS[m].ec50 for m in members}
            activity = module_activity(family, applied_doses)
            vector = np.array([min(1.0, activity[m] + floor[m]) for m in order])
            steady = basal + loadings @ vector
            slowdown = min((growth_rate(m, applied_doses[m]) / growth_rate(m, 0.0)
                            for m in members), default=1.0)
            mu = unstressed * slowdown
            for _ in range(replicates):
                noise = rng.normal(1.0, noise_cv, size=steady.shape) if noise_cv else 1.0
                reading = steady * noise
                if growth_rate_se:
                    error = rng.normal(0.0, growth_rate_se)
                    correction = 1.0 + error / max(mu, 1e-6)
                    if correction <= 0.0:
                        dropped += 1
                        continue
                    reading = reading * correction
                rows.append(reading)
                labels.append(label)
                applied.append(multiple * STRESSORS[members[0]].ec50)
                truth.append(vector)
                growths.append(float(mu))
    wells = (well_positions(len(rows)) if len(rows) <= 96
             else [f"W{i}" for i in range(len(rows))])
    return PanelDataset(
        readings=np.array(rows),
        labels=np.array(labels),
        doses=np.array(applied),
        modules=np.array(truth),
        reporters=reporters,
        wells=wells,
        contexts=[f"{context.growth_phase}/{context.carbon_source}"] * len(rows),
        read_times=[float("inf")] * len(rows),
        growth_rates=growths,
        n_unusable=dropped,
    )


# ---------------------------------------------------------------------------
# the families
# ---------------------------------------------------------------------------


def _literature_drivers() -> dict[str, dict[str, float]]:
    return {name: dict(module.driven_by) for name, module in MODULES.items()
            if module.driven_by}


def _literature_also_reads() -> dict[str, dict[str, float]]:
    return {name: dict(reporter.also_reads) for name, reporter in REPORTERS.items()
            if reporter.also_reads}


def baseline_family(seed: int = 0) -> PanelFamily:
    """The literature topology as `stress_panel.py` states it.

    In the registry because the distribution over configurations has to include the one the
    package already fits, or the comparison has no origin: every other family's score is
    read as a displacement from this one.
    """
    return PanelFamily(
        name="baseline",
        perturbs="nothing -- the crosstalk, dose responses and loadings exactly as the "
                 "panel records them",
        plausible_because="it is the configuration every number in this package was "
                          "obtained on, so it is the origin the others are measured from "
                          "rather than a competitor to them",
        drivers=_literature_drivers(),
        also_reads=_literature_also_reads(),
        draw=seed,
    )


def _edge_dropped(seed: int = 0) -> PanelFamily:
    module, driver = _DROPPED_EDGE
    drivers = _literature_drivers()
    weight = drivers[module].pop(driver)
    return PanelFamily(
        name="edge_dropped",
        perturbs=f"removes the {driver} -> {module} cascade edge, weight {weight:.2f} "
                 "(Hsf1 transcribing RPN4)",
        plausible_because=(
            "the edge is inferred from promoter content. Hahn 2006 (PMID 16556235) confirms "
            "an HSE in the RPN4 promoter by mutagenesis, which establishes the site is "
            "functional under heat -- not that Hsf1 carries a fixed share of proteasome "
            "induction under every agent the panel routes through it. Boy-Marcotte 1999 "
            "(PMID 10411744) finds the Hsf1 and Msn2/4 regulons overlap only slightly and "
            "Ciccarelli 2023 (PMID 37467033) describes the coupling as compensatory rather "
            "than additive, so a constant 0.40 under MG132, MMS and DTT alike asserts more "
            "than either measured. It is also the consequential drop rather than the safe "
            "one: the weaker candidate is iron <- oxidative, which the panel's own source "
            "note already flags as single-source and never replicated, and removing a 0.10 "
            "edge would test almost nothing"
        ),
        drivers=drivers,
        also_reads=_literature_also_reads(),
        draw=seed,
    )


def _edge_added(seed: int = 0) -> PanelFamily:
    module, driver, weight = _ADDED_EDGE
    drivers = _literature_drivers()
    drivers.setdefault(module, {})[driver] = weight
    return PanelFamily(
        name="edge_added",
        perturbs=f"adds a {driver} -> {module} cascade edge at weight {weight:.2f} "
                 "(Msn2/4 transcribing YAP1)",
        plausible_because=(
            "the panel's own annotation for the oxidative module records 'no documented "
            "Msn2/4 transcription of YAP1' -- an absence of evidence rather than a measured "
            "zero, which is exactly the kind of edge whose absence is a modelling choice. "
            "Msn2/4 sits upstream of most of the network, so this is also the most "
            "consequential single addition available: every agent in the panel drives the "
            "ESR, and the edge carries that into the proteasome and iron arms through the "
            "cascade weights already recorded there. The weight sits between the panel's "
            "weakest recorded cascade edge (0.10) and its strongest (0.45), because an "
            "edge nobody has measured should not be given more confidence than one somebody "
            "has"
        ),
        drivers=drivers,
        also_reads=_literature_also_reads(),
        draw=seed,
    )


def _biphasic_family(seed: int = 0) -> PanelFamily:
    return PanelFamily(
        name="biphasic",
        perturbs="gives each module its own falling limb instead of the single viability "
                 "factor every module currently shares: the ESR turns over 1.5-fold later, "
                 "the proteasome 1.3-fold later, the metabolite pools 0.7-fold earlier",
        plausible_because=(
            "the baseline is already biphasic, but through a factor shared by every module, "
            "which freezes the direction of the activity vector above the peak: past the "
            "turnover every regulon falls in lockstep and the high-dose wells add no new "
            "direction. A module-specific limb does not, and that is a difference a "
            "subspace method is supposed to notice. The direction of each scale is sourced "
            "and the magnitude is not. Gasch 2000 shows the ESR being induced under exactly "
            "the conditions where global transcription is being cut back, so it is the "
            "regulon least likely to fall with the average; Rpn4 is degraded by the "
            "proteasome it induces, so impairing the proteasome stabilises the factor and "
            "holds that regulon on while the cell fails; a metabolite pool is displaced by "
            "chemistry rather than defended by transcription and goes with the cell. The "
            "measured evidence that the limb is not one number is a 14% spread in the "
            "growth-halving dose between UPRE1 and UPRE2 under the same agent (1.45 against "
            "1.65 mM DTT), which establishes that it varies without saying how far, so the "
            "magnitudes here are asserted and the family is a probe rather than an estimate"
        ),
        drivers=_literature_drivers(),
        also_reads=_literature_also_reads(),
        lethal_scale=dict(_LIMB_SCALES),
        draw=seed,
    )


def _adapting(seed: int = 0) -> PanelFamily:
    return PanelFamily(
        name="adapting",
        perturbs="reads the general stress response as a decayed transient -- "
                 f"{_ESR_SURVIVING_FRACTION:.0%} of its sustained level at read time -- "
                 "rather than as a level that stays up",
        plausible_because=(
            "Gasch 2000 measured the environmental stress response as a transient: it peaks "
            "within ten to twenty minutes and decays back toward baseline inside the hour "
            "even while the insult persists, whereas the stressor-specific regulons stay on. "
            "A plate read hours after dosing therefore reads an adapted ESR, and the panel "
            "encodes a sustained one. This is the family that bears hardest on the transfer "
            "claim, which is why it is here: the ESR is the component five of the eight "
            "reporters share, so it is the direction a low-rank model has most reason to "
            "find, and this family is the case where that direction is largely gone by the "
            "time anyone measures. The surviving fraction is asserted -- a promoter fusion "
            "averages over about 1/mu = 2.5 h, so it reports a time-average of a decaying "
            "signal rather than either endpoint, and a third is a plausible average and not "
            "a measured one"
        ),
        drivers=_literature_drivers(),
        also_reads=_literature_also_reads(),
        adaptation={"ESR": _ESR_SURVIVING_FRACTION},
        draw=seed,
    )


def _feedback(seed: int = 0) -> PanelFamily:
    module, driver, weight = _FEEDBACK_EDGE
    drivers = _literature_drivers()
    drivers.setdefault(module, {})[driver] = weight
    return PanelFamily(
        name="feedback",
        perturbs=f"closes a negative loop by adding {driver} -> {module} at {weight:.2f}, "
                 "against the heat -> proteasome edge already recorded, so the heat arm "
                 "damps itself",
        plausible_because=(
            "Hsf1 is held off by the chaperones it induces -- the canonical account of heat "
            "shock attenuation, titration of Hsf1 by the Hsp70 pool it drives up "
            "(Krakowiak et al., eLife 2018; Zheng et al., eLife 2016). The panel records the "
            "forward half of that loop, Hsf1 onto RPN4 through the HSE Hahn 2006 confirmed, "
            "and has no way to represent the return: `_propagate` resolves modules in one "
            "fixed order, so a loop cannot be written down there at all. That is a "
            "representational limit rather than a claim that no loop exists, and a family is "
            "the honest way to ask what the limit costs. The return weight is asserted; what "
            "is not asserted is the sign, which is what makes the family a perturbation of "
            "the topology rather than of a number"
        ),
        drivers=drivers,
        also_reads=_literature_also_reads(),
        draw=seed,
    )


def _loading_noise(seed: int = 0) -> PanelFamily:
    rng = np.random.default_rng(seed)
    sigma = np.log(_LOADING_BRACKET) / 1.959963984540054
    perturbed = {}
    for reporter, edges in sorted(_literature_also_reads().items()):
        perturbed[reporter] = {
            module: float(weight * np.exp(sigma * rng.normal()))
            for module, weight in sorted(edges.items())
        }
    return PanelFamily(
        name="loading_noise",
        perturbs="redraws every crosstalk loading, each within a factor of "
                 f"{_LOADING_BRACKET:.0f} at 95%, leaving the topology alone",
        plausible_because=(
            "this is the control the other families are read against, and it is here "
            "precisely because it is a re-draw and not a rewiring. The crosstalk loadings "
            "are the largest block of the roughly 103 panel parameters with no literature "
            "value: the statement behind each is qualitative -- 'KAR2 carries both a UPRE "
            "and an HSE' -- which fixes the sign and not the coefficient, so a two-fold "
            "bracket is close to the least such a sentence can mean. If transfer moves as "
            "far here as it does under a dropped edge, then what the structural families "
            "measure is sensitivity to unmeasured constants and not to structure, and the "
            "cross-family comparison says nothing about topology. The reporter's own "
            "element is left at 1.0 because it is a normalisation of the channel's units "
            "rather than a measurement, and moving it would rescale the channel"
        ),
        drivers=_literature_drivers(),
        also_reads=perturbed,
        draw=seed,
    )


FAMILY_BUILDERS: Mapping[str, Callable[[int], PanelFamily]] = {
    "baseline": baseline_family,
    "edge_dropped": _edge_dropped,
    "edge_added": _edge_added,
    "biphasic": _biphasic_family,
    "adapting": _adapting,
    "feedback": _feedback,
    "loading_noise": _loading_noise,
}
"""The configurations, keyed by name so a caller iterates rather than hardcodes.

A registry rather than a list of functions because the point of the exercise is that no
number should be reported at one configuration: a caller that has to name a family to get
one is a caller that will name `baseline` and stop.

Absent, and worth recording as absent: a family that perturbs the cascade *weights*
continuously. It is the same axis as `edge_dropped` -- dropping an edge is the limit of
shrinking it -- so it would add a draw rather than a structure, which `loading_noise`
already covers for the loadings.
"""


def family_names() -> tuple[str, ...]:
    """Registry keys, in registration order with `baseline` first."""
    return tuple(FAMILY_BUILDERS)


def build_family(name: str, seed: int = 0) -> PanelFamily:
    """One named configuration.

    Args:
        name: Registry key.
        seed: Draw behind any stochastic part. Deterministic families ignore it, so a
            caller can seed uniformly without knowing which is which.
    """
    if name not in FAMILY_BUILDERS:
        raise KeyError(f"no family {name!r}; have {sorted(FAMILY_BUILDERS)}")
    return FAMILY_BUILDERS[name](seed)


def sample_family(rng: np.random.Generator) -> PanelFamily:
    """One draw from the distribution over configurations, uniform over the registry.

    This is the `p(xi)` every bound in `analysis/optimism.py` is stated with respect to,
    so what it is matters as much as the bound: uniform over seven families is an assumed
    distribution, not a measured one, and nothing here estimates how far a real plate sits
    from it. Uniform because there is no evidence with which to weight one alternative
    topology above another -- weighting them by plausibility would be asserting exactly the
    knowledge whose absence makes the families necessary.
    """
    name = str(rng.choice(np.array(family_names())))
    return build_family(name, seed=int(rng.integers(0, 2**31 - 1)))
