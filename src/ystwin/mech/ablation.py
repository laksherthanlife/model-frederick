"""Criterion (a), made executable: an ablation compares two FITTED models, never a trace
against a frozen incumbent.

**Why this file exists, and it is a methodological repair rather than a feature.** Four of
the five block rescues in the 2026-09-05 adjudication failed the same way, and
`REVISED_BUILD_LIST.md` §5 records it: they scored a mechanistic trace against an incumbent
that had been frozen, stripped of a term, or denied its own free parameter. The clearest
case is measured -- refitting the *real* incumbent (`generator/culture.py::simulate_culture`,
logistic, free K) to the TORC1 advocate's own lagged trace gives a relative RMSE of
0.0021-0.0068, **0.015x to 0.047x the 0.146 floor**, where the advocate had reported a pass
after deleting the incumbent's capacity term and comparing against ``X0*exp(mu*t)``. That
error is one-directional: every instance of it inflates the mechanistic block's apparent
contribution. With ten blocks scored against one incumbent, an architecture that does not
write the rule down will repeat it against itself.

**The rule.**

1. Both models are refit, on the same data, in the same call. :func:`ablate` fits them; it
   does not accept two fits a caller prepared.
2. Both parameter counts are reported. A block that wins by spending a parameter has to
   show the parameter.
3. The reduced model must be a restriction of the full one. More free scalars in the
   reduced model is not an ablation, and :func:`nested_out_of_sample` refuses it by name.
4. The floor is a :class:`~ystwin.mech.params.Target`, not a number. A ``Target`` carries
   the assay, and the assay is checked against the observable's own. **Scoring an FBA
   product ceiling against a plate-reader activity CV is the specific error this prevents**
   -- it is what sank Route C2, in a repository whose own memory records the product
   ceiling being refuted by 63x.

**And the rule is enforced by measurement, not by a docstring.** A frozen model is detected
rather than declared: :func:`ablate` refits each side to a deterministically perturbed copy
of the same data and requires the parameters to move. A model that returns the same numbers
whatever it is shown is frozen however it describes itself, and raises
:class:`FrozenIncumbent`. A model with no free scalars at all is legitimate -- it is
reported as parameter-free, which is a strength, not a failure.

**The two floors this platform has measured**, both in `generator/panel_experiment.py` and
both re-exported here as ``Target`` rows with their assays attached:

| Observable | Floor | Assay |
|---|---|---|
| dilution-corrected reporter activity | 0.146 relative CV | 96-well plate reader, mCitrine 480/530, 28 matched conditions on two plates |
| specific growth rate | 0.0117 /h | slope of log OD600 against time, median of 210 real wells |

There are exactly two. Every other observable must arrive with its own named assay --
:data:`REFUSED_FLOORS` records the ones that have been asked for here and have none, with
the reason, so the next reader gets a refusal naming the missing measurement instead of a
plausible number.

**The worked instance, and its result.** The pair this module ships is the nested pair the
revised build list wants run in Phase 0 (Tier 0.3): the incumbent's own reporter algebra --
an instantaneous saturating dose response driving one first-order reporter pool, which is
`generator/stress_panel.py::module_response` composed with
`generator/kinetics.py::reporter_at_time` -- against the same thing plus the L8 maturation
state, on the committed NewProtocol dose ladders. Run leave-one-dose-out on all twelve
committed (plate, construct) blocks by `tests/test_ablation.py`, the maturation state
improves interpolating out-of-sample prediction by at most **0.01970 of training scale =
0.135x the 0.146 floor**, wins on ten of the twelve blocks and loses on two, and moves the
endpoint fold induction by at most **0.02975 = 0.204x the floor**, clearing it on none of
the twelve. That is a fail on criterion (a), and it is the first ablation in this project
scored under this rule rather than against a frozen incumbent.

Two things that result is NOT. It is not experiment 0.2/A4 -- the lag-1 innovation
autocorrelation of +0.954 lives in `estimator.py` and is a different measurement. And the
H2O2 blocks are a comparison between two models that both miss: the 2 and 4 mM peroxide
wells are the ones `generator/panel_experiment.py::PanelDataset.n_unusable` documents as
going negative after the dilution correction, and neither model describes them.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from ..generator.panel_experiment import MEASURED_GROWTH_RATE_SE, OBSERVED_ACTIVITY_CV
from ..plate import replay
from .params import FreeScalarGate, FreeScalarGateFailed, Param, ParamRegistry, Target

__all__ = [
    "ABLATION_PARAMS",
    "AblationResult",
    "FittableModel",
    "Fit",
    "FloorNotNamed",
    "FrozenIncumbent",
    "GROWTH_RATE_FLOOR",
    "MaturingReporter",
    "NestedComparison",
    "NotNested",
    "Observable",
    "REFUSED_FLOORS",
    "REPORTER_ACTIVITY_FLOOR",
    "RELATIVE_FLOOR_UNITS",
    "SaturatingReporter",
    "ablate",
    "endpoint_fold_induction",
    "floor_for",
    "free_scalar_gate",
    "nested_out_of_sample",
    "reporter_block",
    "require_free_scalar_gate",
]


# --------------------------------------------------------------------------------------
# The floors, and the refusal to invent one
# --------------------------------------------------------------------------------------

REPORTER_ACTIVITY_FLOOR = Target(
    name="reporter_activity",
    observable="dilution-corrected mCitrine promoter activity",
    assay="96-well plate reader, mCitrine 480/530, OD600-corrected",
    noise_floor=OBSERVED_ACTIVITY_CV,
    units="relative CV",
    source=("generator/panel_experiment.py::OBSERVED_ACTIVITY_CV -- total plate-to-plate CV "
            "across 28 matched conditions on the 2026-07-22 and 2026-08-03 plates"),
)
"""0.146. The only floor that may be used to score a reporter signal, and only a reporter signal."""

GROWTH_RATE_FLOOR = Target(
    name="growth_rate",
    observable="specific growth rate",
    assay="slope of log OD600 against time, 96-well plate reader",
    noise_floor=MEASURED_GROWTH_RATE_SE,
    units="1/h",
    source=("generator/panel_experiment.py::MEASURED_GROWTH_RATE_SE -- standard error of the "
            "log-OD slope, median across 210 real wells with a fittable window, "
            "IQR 0.0079 to 0.0200"),
)
"""0.0117 /h, and ABSOLUTE. It is 3% of a healthy rate and 12% of one slowed to 0.10, which is
why it must not be turned into a CV at the point of use."""

RELATIVE_FLOOR_UNITS = ("relative CV", "fraction", "dimensionless")
"""Floor units that may score a relative effect. A declared vocabulary rather than a string
sniff: an absolute comparison requires the floor's units to equal the observable's exactly."""

_FLOORS: dict[str, Target] = {
    REPORTER_ACTIVITY_FLOOR.name: REPORTER_ACTIVITY_FLOOR,
    GROWTH_RATE_FLOOR.name: GROWTH_RATE_FLOOR,
}

REFUSED_FLOORS: dict[str, str] = {
    "product_ceiling": (
        "an FBA product ceiling is an LP bound, not a reading. It has no assay and therefore "
        "no noise floor here. Route C2 scored one against the plate-reader activity CV and "
        "REVISED_BUILD_LIST.md records that as the wrong floor for the wrong observable -- in "
        "a repository whose own memory records the ceiling being refuted by 63x. To score a "
        "ceiling, name the analytical method that would measure the product and its own "
        "replicate CV"),
    "titre": (
        "no titre assay has been run on this project's strains. The 22.2% figure that gets "
        "reached for is the entry-flux law's leave-one-STRAIN-out generalisation error "
        "(bridge/stress_diversion.py), which is a cross-strain prediction error and not the "
        "within-strain repeatability of any instrument. Name the HPLC or spectrophotometric "
        "method and its replicate CV"),
    "product_content": (
        "mg/gDCW of a carotenoid read on a plate reader sums lycopene with beta-carotene and "
        "reads above the stoichiometric ceiling -- docs and the project memory both record "
        "that absorbance is not specific. An HPLC method and its own CV would close it"),
    "flux": (
        "a GEM flux is a solved quantity. 13C-MFA confidence intervals would be the floor and "
        "none has been run here"),
}
"""Observables that have been scored here, or nearly were, and have no floor on this platform.

Each entry is the message :func:`floor_for` raises. Recorded rather than left absent so that
the next reader gets the missing measurement named instead of re-deriving the borrowing."""


class FloorNotNamed(ValueError):
    """No floor was named for an observable, or the one named comes from a different assay.

    ``ValueError`` and not ``NotImplementedError`` because this is a fault in the CALL, not
    an absent measurement: the caller has to name a floor and its assay, and the fix is at
    the call site. `params.RefusedValue` is the other one, for a parameter that has no
    number at all.
    """


class FrozenIncumbent(ValueError):
    """A model presented for ablation does not respond to its data.

    The failure this module exists to prevent, caught by measurement rather than by trust:
    the model was refit to a perturbed copy of the same rows and returned the same
    parameters. Either it ignores the data, or its ``refit`` is returning constants.
    """


class NotNested(ValueError):
    """The reduced model is not a restriction of the full one."""


def floor_for(observable: Observable) -> Target:
    """The registered floor for an observable, or a refusal naming what is missing.

    Args:
        observable: What is being scored. Both its ``name`` and its ``assay`` are checked --
            a floor from a different instrument is refused even when the name matches.

    Raises:
        FloorNotNamed: if nothing is registered under that name, if the name is one of
            :data:`REFUSED_FLOORS`, or if the registered floor's assay is not the
            observable's.
    """
    name = observable.name
    if name in REFUSED_FLOORS:
        raise FloorNotNamed(
            f"observable {name!r} has no noise floor on this platform. {REFUSED_FLOORS[name]}")
    try:
        target = _FLOORS[name]
    except KeyError:
        raise FloorNotNamed(
            f"no floor is registered for observable {name!r}. This repository has measured "
            f"exactly two: {', '.join(sorted(_FLOORS))}. For anything else, NAME ITS ASSAY "
            f"and its floor and pass a mech.params.Target -- do not borrow the plate CV, "
            f"which is what REVISED_BUILD_LIST.md records as the error that sank Route C2"
        ) from None
    _require_matching_assay(observable, target)
    return target


def _require_matching_assay(observable: Observable, floor: Target) -> None:
    """The check that stops a plate-reader CV scoring a flux. Assays are declared on both
    objects, so this compares data rather than guessing from a name."""
    if observable.assay != floor.assay:
        raise FloorNotNamed(
            f"observable {observable.name!r} is measured by {observable.assay!r} and the "
            f"floor offered ({floor.name!r}, {floor.noise_floor:g} {floor.units}) was "
            f"measured on {floor.assay!r}. A floor from a different instrument is not this "
            f"observable's noise. Source of the floor offered: {floor.source}")


# --------------------------------------------------------------------------------------
# What a model has to be able to do
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Fit:
    """One model, fitted to one frame.

    Args:
        model: Name of the model that produced it.
        parameters: Fitted values, keyed by the model's own ``parameter_names``.
        n_free: How many scalars the fit searched over. Reported in every ablation, because
            a block that wins by spending a parameter has to show the parameter.
        n_rows: Rows fitted. Carried so two fits can be checked to have seen the same data.
        residual_rms: Root-mean-square residual in the response column's own units.
        converged: What the optimiser said. A non-converged fit is not an error here -- it
            is reported, because silently retrying until it converges is a search over
            starting points that nobody counted.
    """

    model: str
    parameters: Mapping[str, float]
    n_free: int
    n_rows: int
    residual_rms: float
    converged: bool


class FittableModel:
    """The interface :func:`ablate` requires. Subclass or duck-type it.

    Four members, and each exists because the protocol enforces something with it:
    ``name`` for the report, ``parameter_names`` for the count and the nesting check,
    ``response_column`` so the frozen-model probe knows what to perturb, and the two methods.

    The method is ``refit`` and not ``fit`` because that is the contract: a caller does
    not hand in a fitted model, it hands in something that will be fitted here. It must
    return a NEW :class:`Fit` from the frame it is given -- one that ignores its argument is
    exactly the frozen incumbent this module refuses, and :func:`ablate` will detect it.
    """

    name: str = ""
    parameter_names: tuple[str, ...] = ()
    response_column: str = "signal"

    def refit(self, data: pd.DataFrame) -> Fit:
        raise NotImplementedError

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        raise NotImplementedError


@dataclass(frozen=True)
class Observable:
    """A number both models predict, the data they are fitted to, and how it is scored.

    Args:
        name: Keys :func:`floor_for`. Use a registered one or bring a ``Target``.
        units: Of the summarised value.
        assay: The instrument that would measure it. Checked against the floor's assay --
            see :func:`_require_matching_assay`.
        scored: ``"relative"`` or ``"absolute"``. A relative effect is
            ``|full - reduced| / |reduced|`` and needs a dimensionless floor; an absolute
            effect is ``|full - reduced|`` and needs a floor in the observable's own units.
        data: The rows BOTH models are fitted to. One frame, so they cannot differ.
        summarise: ``(model, fit) -> float``. Takes the model as well as the fit so it can
            evaluate the prediction anywhere, not only where there is a measurement.
        description: One line for the report.
    """

    name: str
    units: str
    assay: str
    scored: str
    data: pd.DataFrame
    summarise: Callable[[FittableModel, Fit], float]
    description: str = ""

    def __post_init__(self) -> None:
        if self.scored not in ("relative", "absolute"):
            raise ValueError(
                f"observable {self.name!r}: scored must be 'relative' or 'absolute', got "
                f"{self.scored!r}")
        for label, text in (("units", self.units), ("assay", self.assay)):
            if not str(text).strip():
                raise ValueError(
                    f"observable {self.name!r}: {label} is required. Criterion (a) is scored "
                    "against a NAMED assay and an unnamed one cannot be checked against a floor")
        if len(self.data) == 0:
            raise ValueError(f"observable {self.name!r}: no rows to fit")


# --------------------------------------------------------------------------------------
# The frozen-incumbent probe
# --------------------------------------------------------------------------------------

ABLATION_PARAMS = ParamRegistry("mech/ablation.py")
"""This module's own constants. Three, all numerical rather than physical, all ASSERTED.

Deliberately NOT put through ``require_gate()``: none of them enters a prediction, so none is
a degree of freedom criterion (e) is about. They decide whether a check fires, and a check
that fires on 5% and on 8% alike has not spent a parameter. The scalars criterion (e) counts
for an ablation are the fitted model's, and :func:`free_scalar_gate` counts those."""

_PROBE_AMPLITUDE = ABLATION_PARAMS.add(Param.asserted(
    "probe_amplitude", 0.05, "fraction of the response column",
    "ASSERTED. Size of the deterministic perturbation used to prove a model responds to its "
    "data. Any value a converged optimiser can see does the job; 5% is comfortably above "
    "least_squares' own xtol/ftol of 1e-8 and comfortably below the range of every response "
    "column in this repository",
    missing="nothing -- this is a numerical probe, not a measurement"))

_PROBE_MOVED = ABLATION_PARAMS.add(Param.asserted(
    "probe_moved", 1e-6, "relative",
    "ASSERTED. A fitted parameter counts as having moved when it moves by more than this "
    "fraction of itself. Two orders above least_squares' default tolerances, so a converged "
    "refit on perturbed data clears it and a frozen one cannot",
    missing="nothing -- this is a numerical tolerance"))

_FIT_MAX_NFEV = ABLATION_PARAMS.add(Param.asserted(
    "fit_max_nfev", 20000, "function evaluations",
    "ASSERTED. Evaluation budget for one least-squares fit. Reached by no fit in "
    "tests/test_ablation.py; it exists so a pathological frame fails loudly instead of hanging",
    missing="nothing -- this is a solver budget"))


def _probe_copy(data: pd.DataFrame, column: str) -> pd.DataFrame:
    """The same rows with the response tilted by a fixed, reproducible pattern.

    Deterministic on purpose. A random perturbation would make the frozen-model check give a
    different answer on a rerun, which is the one thing a gate must not do.
    """
    if column not in data.columns:
        raise KeyError(
            f"model declares response_column={column!r}, which is not in the frame "
            f"({', '.join(map(str, data.columns))})")
    tilt = 1.0 + float(_PROBE_AMPLITUDE) * np.cos(
        2.0 * np.pi * np.arange(len(data)) / max(len(data), 1))
    probed = data.copy()
    probed[column] = np.asarray(probed[column], dtype=float) * tilt
    return probed


def _responds_to_its_data(model: FittableModel, data: pd.DataFrame, base: Fit) -> bool:
    """Refit on a perturbed copy and report whether any parameter moved."""
    if base.n_free == 0:
        return True
    alternative = model.refit(_probe_copy(data, model.response_column))
    tolerance = float(_PROBE_MOVED)
    return any(
        abs(alternative.parameters[k] - base.parameters[k])
        > tolerance * max(1.0, abs(base.parameters[k]))
        for k in base.parameters)


# --------------------------------------------------------------------------------------
# Criterion (a)
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class AblationResult:
    """What one ablation measured, and everything needed to check it was measured fairly."""

    observable: str
    units: str
    description: str
    floor: Target
    scored: str
    full_name: str
    full_n_free: int
    full_value: float
    full_parameter_free: bool
    reduced_name: str
    reduced_n_free: int
    reduced_value: float
    reduced_parameter_free: bool
    n_rows: int
    effect: float
    ratio: float

    @property
    def clears_floor(self) -> bool:
        """Criterion (a). Strictly greater: an effect equal to the floor is not above it."""
        return self.ratio > 1.0

    def report(self) -> str:
        verdict = "CLEARS THE FLOOR" if self.clears_floor else "DOES NOT CLEAR THE FLOOR"
        lines = [
            f"ABLATION  {self.observable} [{self.units}]"
            + (f" -- {self.description}" if self.description else ""),
            f"  full     {self.full_name:<24s} {self.full_n_free} free scalars"
            + ("  (parameter-free)" if self.full_parameter_free else "")
            + f"   value {self.full_value:.6g}",
            f"  reduced  {self.reduced_name:<24s} {self.reduced_n_free} free scalars"
            + ("  (parameter-free)" if self.reduced_parameter_free else "")
            + f"   value {self.reduced_value:.6g}",
            f"  both refit on the same {self.n_rows} rows; neither is frozen",
            f"  effect   {self.effect:.6g} ({self.scored})",
            f"  floor    {self.floor.noise_floor:g} {self.floor.units} -- assay: {self.floor.assay}",
            f"           {self.floor.source}",
            f"  verdict  {self.ratio:.3g}x the floor -- {verdict}",
        ]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.report()


def ablate(full_model: FittableModel, reduced_model: FittableModel,
           observable: Observable, floor: Target | None = None) -> AblationResult:
    """Score removing a piece, with both models refit at stated parameter counts.

    This is the whole rule in one call. Both models are fitted here, to
    ``observable.data``, so neither can arrive frozen or stripped; both are then probed to
    prove they respond to that data at all; and the effect is scored against a floor whose
    assay must be the observable's own.

    Args:
        full_model: The model WITH the piece.
        reduced_model: The model with the piece removed. Must not have more free scalars
            than ``full_model``.
        observable: What moves, the rows both models see, and how it is summarised.
        floor: The noise the effect must beat, as a
            :class:`~ystwin.mech.params.Target` carrying its assay. ``None`` looks the
            observable up in the registry, which refuses for anything this platform has not
            measured.

    Returns:
        An :class:`AblationResult` whose ``report()`` names the effect, both parameter
        counts, the floor and the assay the floor came from.

    Raises:
        FloorNotNamed: no floor, an unregistered observable, or a floor from another assay.
        FrozenIncumbent: either model returned the same parameters on perturbed data.
        ValueError: the reduced model has more free scalars than the full one, or the floor's
            units do not fit how the observable is scored.
    """
    if floor is None:
        floor = floor_for(observable)
    elif not isinstance(floor, Target):
        raise FloorNotNamed(
            f"observable {observable.name!r}: floor must be a mech.params.Target, got "
            f"{type(floor).__name__}. A bare number names no assay, and a floor without an "
            f"assay is what lets a plate-reader CV score an FBA quantity")
    else:
        _require_matching_assay(observable, floor)

    if observable.scored == "absolute" and floor.units != observable.units:
        raise ValueError(
            f"observable {observable.name!r} is scored absolutely in {observable.units!r} and "
            f"its floor is in {floor.units!r}. An absolute comparison needs the same units on "
            f"both sides; score it relatively, or bring a floor in the right units")
    if observable.scored == "relative" and floor.units not in RELATIVE_FLOOR_UNITS:
        raise ValueError(
            f"observable {observable.name!r} is scored relatively and its floor is in "
            f"{floor.units!r}, which is not one of {RELATIVE_FLOOR_UNITS}. A dimensional floor "
            f"cannot bound a dimensionless effect -- {GROWTH_RATE_FLOOR.name!r} in particular "
            f"is 0.0117 /h ABSOLUTE and is 3% of a healthy rate against 12% of a slowed one")

    if reduced_model.parameter_names and full_model.parameter_names:
        if len(reduced_model.parameter_names) > len(full_model.parameter_names):
            raise ValueError(
                f"{reduced_model.name!r} has {len(reduced_model.parameter_names)} free scalars "
                f"and {full_model.name!r} has {len(full_model.parameter_names)}. Removing a "
                f"piece cannot add a degree of freedom; this is not an ablation")

    data = observable.data
    full_fit = full_model.refit(data)
    reduced_fit = reduced_model.refit(data)
    for model, fit in ((full_model, full_fit), (reduced_model, reduced_fit)):
        if not _responds_to_its_data(model, data, fit):
            raise FrozenIncumbent(
                f"{model.name!r} declares {fit.n_free} free scalars and returned identical "
                f"parameters when refit on a perturbed copy of the same {len(data)} rows. A "
                f"frozen model in an ablation inflates the other side's contribution, which "
                f"is the failure REVISED_BUILD_LIST.md records in four of five block rescues")

    full_value = float(observable.summarise(full_model, full_fit))
    reduced_value = float(observable.summarise(reduced_model, reduced_fit))
    difference = abs(full_value - reduced_value)
    if observable.scored == "relative":
        if reduced_value == 0.0:
            raise ValueError(
                f"observable {observable.name!r} is scored relatively and the reduced model "
                f"predicts exactly zero, so there is no scale to be relative to")
        effect = difference / abs(reduced_value)
    else:
        effect = difference

    return AblationResult(
        observable=observable.name,
        units=observable.units,
        description=observable.description,
        floor=floor,
        scored=observable.scored,
        full_name=full_model.name,
        full_n_free=full_fit.n_free,
        full_value=full_value,
        full_parameter_free=full_fit.n_free == 0,
        reduced_name=reduced_model.name,
        reduced_n_free=reduced_fit.n_free,
        reduced_value=reduced_value,
        reduced_parameter_free=reduced_fit.n_free == 0,
        n_rows=len(data),
        effect=effect,
        ratio=effect / floor.noise_floor,
    )


# --------------------------------------------------------------------------------------
# Criterion (e), as a callable check
# --------------------------------------------------------------------------------------

def free_scalar_gate(model: FittableModel, targets: Sequence[Target],
                     piece: str | None = None) -> FreeScalarGate:
    """Count a model's fitted scalars against the independent targets it is scored on.

    Every fitted scalar is free by definition -- there is no measurement pinning it, which
    is what fitting means -- so this counts ``parameter_names`` rather than asking the model
    to grade itself. Targets marked ``fitted`` do not count, for the reason
    `ARCHITECTURE_GAPS.md` 0.3 gives: scoring on whether some point in a sweep reproduces a
    number is fitting with extra steps, and such a target in the denominator lets a sweep
    pay for itself.

    Args:
        model: The model whose scalars are counted.
        targets: Independent measurements it is scored against.
        piece: Name for the gate. Defaults to the model's.

    Returns:
        A :class:`~ystwin.mech.params.FreeScalarGate`; call ``.passes``, or use
        :func:`require_free_scalar_gate` to refuse.
    """
    return FreeScalarGate(
        piece=piece or model.name,
        free=tuple(model.parameter_names),
        targets=tuple(t.name for t in targets if not t.fitted),
        fitted_targets=tuple(t.name for t in targets if t.fitted),
    )


def require_free_scalar_gate(model: FittableModel, targets: Sequence[Target],
                             piece: str | None = None) -> FreeScalarGate:
    """:func:`free_scalar_gate`, refusing rather than reporting. Criterion (e) with teeth."""
    gate = free_scalar_gate(model, targets, piece)
    if not gate.passes:
        raise FreeScalarGateFailed(
            f"{gate.summary()}. Free: {', '.join(gate.free)}. Targets: "
            f"{', '.join(gate.targets) or '(none)'}. A piece with more free scalars than "
            f"independent targets adds mechanism without adding refutability")
    return gate


# --------------------------------------------------------------------------------------
# The nested out-of-sample comparison (Tier 0.3)
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class NestedComparison:
    """Leave-one-group-out prediction error for two nested models at stated counts.

    ``full_error`` and ``reduced_error`` are keyed by the held-out group and are RMSE over
    that group's rows, divided by the mean absolute response of the TRAINING rows of the same
    fold. Relative, so it can be compared with a CV floor; and relative to the training fold
    rather than to the held-out one, so the denominator carries no information the fit did
    not have.
    """

    full_name: str
    reduced_name: str
    full_n_free: int
    reduced_n_free: int
    group_column: str
    full_error: Mapping[float | str, float]
    reduced_error: Mapping[float | str, float]
    floor: Target
    gate: FreeScalarGate
    interior_groups: tuple = field(default_factory=tuple)

    @property
    def groups(self) -> tuple:
        return tuple(self.full_error)

    def mean_error(self, groups: Sequence | None = None) -> tuple[float, float]:
        """``(full, reduced)`` mean relative error over the named folds, all of them by default."""
        keys = tuple(self.groups if groups is None else groups)
        if not keys:
            raise ValueError("no folds selected")
        return (
            float(np.mean([self.full_error[k] for k in keys])),
            float(np.mean([self.reduced_error[k] for k in keys])),
        )

    def delta(self, groups: Sequence | None = None) -> float:
        """How much the FULL model beats the reduced one out of sample. Negative means it loses."""
        full, reduced = self.mean_error(groups)
        return reduced - full

    def ratio(self, groups: Sequence | None = None) -> float:
        return self.delta(groups) / self.floor.noise_floor

    @property
    def full_model_wins(self) -> bool:
        return self.delta() > 0.0

    @property
    def clears_floor(self) -> bool:
        return self.ratio() > 1.0

    def report(self) -> str:
        full_all, reduced_all = self.mean_error()
        lines = [
            f"NESTED OUT-OF-SAMPLE  leave-one-{self.group_column}-out, {len(self.groups)} folds",
            f"  full     {self.full_name:<24s} {self.full_n_free} free scalars   "
            f"mean relative error {full_all:.4f}",
            f"  reduced  {self.reduced_name:<24s} {self.reduced_n_free} free scalars   "
            f"mean relative error {reduced_all:.4f}",
            f"  gate     {self.gate.summary()}",
        ]
        for key in self.groups:
            lines.append(
                f"    held out {key!s:>8s}   full {self.full_error[key]:.4f}   "
                f"reduced {self.reduced_error[key]:.4f}   "
                f"delta {self.reduced_error[key] - self.full_error[key]:+.4f}")
        lines.append(
            f"  delta    {self.delta():+.4f} = {self.ratio():+.3g}x the "
            f"{self.floor.noise_floor:g} {self.floor.units} floor ({self.floor.assay})")
        if self.interior_groups:
            lines.append(
                f"  interior {self.delta(self.interior_groups):+.4f} = "
                f"{self.ratio(self.interior_groups):+.3g}x, over the "
                f"{len(self.interior_groups)} folds that are interpolation rather than "
                f"extrapolation")
        lines.append(
            "  verdict  " + ("full model wins out of sample" if self.full_model_wins
                             else "full model does NOT beat the reduced one out of sample")
            + (" and clears the floor" if self.clears_floor else "; below the floor"))
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.report()


def nested_out_of_sample(full_model: FittableModel, reduced_model: FittableModel,
                         data: pd.DataFrame, floor: Target,
                         group_column: str = "dose_mM") -> NestedComparison:
    """Leave-one-group-out comparison of two nested models, both refit on every fold.

    The comparison Tier 0.3 of `REVISED_BUILD_LIST.md` asks for: current algebra against the
    mechanistic model with all sweeps free, leave-one-dose-out, parameters counted. If the
    richer model does not beat the simpler one on data already in hand, later phases will not
    rescue it.

    Nesting is checked, not assumed: ``reduced_model.parameter_names`` must be a strict
    subset of ``full_model.parameter_names``. Two models that merely both fit the data are
    not a nested pair and their error difference is not attributable to the extra piece.

    Args:
        full_model: The richer model.
        reduced_model: Its restriction.
        data: All folds. Split by ``group_column``.
        floor: What the difference is scored against, with its assay.
        group_column: The column held out one value at a time.

    Returns:
        A :class:`NestedComparison`. Its ``interior_groups`` are the folds that are
        interpolation, i.e. every numeric group except the smallest and largest, because a
        leave-one-dose-out on a ladder necessarily makes its two end folds extrapolation.

    Raises:
        NotNested: the reduced model's parameters are not a strict subset of the full one's.
        FreeScalarGateFailed: the full model has more free scalars than there are folds.
    """
    full_names, reduced_names = set(full_model.parameter_names), set(reduced_model.parameter_names)
    if not reduced_names < full_names:
        raise NotNested(
            f"{reduced_model.name!r} has parameters {sorted(reduced_names)} and "
            f"{full_model.name!r} has {sorted(full_names)}. A nested comparison needs the "
            f"reduced model to be a strict restriction of the full one; otherwise the "
            f"difference in error is not attributable to the piece being ablated")
    if group_column not in data.columns:
        raise KeyError(
            f"no column {group_column!r} to leave out; the frame has "
            f"{', '.join(map(str, data.columns))}")

    groups = sorted(pd.unique(data[group_column]))
    targets = tuple(
        Target(name=f"{floor.name}@{group_column}={g!s}", observable=floor.observable,
               assay=floor.assay, noise_floor=floor.noise_floor, units=floor.units,
               source=f"{floor.source}; held-out fold {group_column}={g!s}")
        for g in groups)
    gate = require_free_scalar_gate(
        full_model, targets, piece=f"{full_model.name} on leave-one-{group_column}-out")

    full_error: dict = {}
    reduced_error: dict = {}
    column = full_model.response_column
    for held in groups:
        train = data[data[group_column] != held]
        test = data[data[group_column] == held]
        scale = float(np.mean(np.abs(np.asarray(train[column], dtype=float))))
        if scale == 0.0:
            raise ValueError(
                f"fold {group_column}={held!s}: the training rows have zero mean absolute "
                f"{column}, so there is no scale to report a relative error against")
        for model, store in ((full_model, full_error), (reduced_model, reduced_error)):
            fit = model.refit(train)
            residual = model.predict(test, fit.parameters) - np.asarray(
                test[column], dtype=float)
            store[held] = float(np.sqrt(np.mean(residual ** 2)) / scale)

    numeric = all(isinstance(g, (int, float, np.integer, np.floating)) for g in groups)
    interior = tuple(groups[1:-1]) if numeric and len(groups) > 2 else ()
    return NestedComparison(
        full_name=full_model.name,
        reduced_name=reduced_model.name,
        full_n_free=len(full_model.parameter_names),
        reduced_n_free=len(reduced_model.parameter_names),
        group_column=group_column,
        full_error=full_error,
        reduced_error=reduced_error,
        floor=floor,
        gate=gate,
        interior_groups=interior,
    )


# --------------------------------------------------------------------------------------
# The worked instance: the incumbent's reporter algebra, and the same plus one state
# --------------------------------------------------------------------------------------

def _dose_response_columns(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    for column in ("time_h", "dose_mM", "signal"):
        if column not in data.columns:
            raise KeyError(
                f"a reporter block needs columns time_h, dose_mM and signal; got "
                f"{', '.join(map(str, data.columns))}")
    return (np.asarray(data.time_h, dtype=float),
            np.asarray(data.dose_mM, dtype=float),
            np.asarray(data.signal, dtype=float))


def _start_and_bounds(t: np.ndarray, d: np.ndarray, y: np.ndarray) -> tuple[dict, dict, dict]:
    """Initial guess and box, both derived from the frame's own scale.

    Nothing here is a claim about yeast. Bounds set from the data rather than from absolute
    magnitudes is the same argument `pathway/flux.py` makes about its extrapolation window:
    a guard whose bounds are unrelated to the data cannot fire where it matters.
    """
    scale = float(np.max(np.abs(y))) or 1.0
    first = float(np.mean(y[t == t.min()]))
    rate = 0.5
    positive = d[d > 0]
    start = {
        "R0": first,
        "k_basal": abs(first) * rate,
        "k_max": max(scale - abs(first), abs(first)) * rate,
        "K_dose": float(np.median(positive)) if positive.size else 1.0,
        "lambda_": rate,
        "tau_mat": 0.25,
    }
    low = {"R0": -10.0 * scale, "k_basal": 0.0, "k_max": 0.0, "K_dose": 1e-3,
           "lambda_": 1e-3, "tau_mat": 0.0}
    high = {"R0": 10.0 * scale, "k_basal": 100.0 * scale, "k_max": 1000.0 * scale,
            "K_dose": 1e3, "lambda_": 20.0, "tau_mat": 10.0}
    return start, low, high


def _least_squares(residual, names, start, low, high) -> tuple[dict, bool]:
    x0 = [min(max(start[n], low[n]), high[n]) for n in names]
    solution = least_squares(
        residual, x0, bounds=([low[n] for n in names], [high[n] for n in names]),
        max_nfev=int(float(_FIT_MAX_NFEV)))
    return dict(zip(names, (float(v) for v in solution.x))), bool(solution.status > 0)


@dataclass(frozen=True)
class SaturatingReporter:
    """The incumbent, in one expression, with its own free parameters intact.

    ``k(d) = k_basal + k_max * d / (K_dose + d)`` is
    `generator/stress_panel.py::module_response`'s saturating dose term, and
    ``dR/dt = k - lambda*R`` is `generator/kinetics.py::reporter_at_time`'s balance, whose
    relaxation rate is dilution plus degradation. Dose enters instantaneously: that is the
    property this project set out to replace, and it is the property being ablated against.

    Five free scalars, all fitted, none claimed as measured. Written out here rather than
    imported so it can be refit -- the incumbent this repository has been scoring against
    was the one whose parameters were pinned, which is the whole finding.
    """

    name: str = "SaturatingReporter"
    parameter_names: tuple[str, ...] = ("R0", "k_basal", "k_max", "K_dose", "lambda_")
    response_column: str = "signal"

    def _curve(self, t, d, R0, k_basal, k_max, K_dose, lambda_):
        k = k_basal + k_max * d / (K_dose + d)
        decay = np.exp(-lambda_ * t)
        return R0 * decay + (k / lambda_) * (1.0 - decay)

    def predict(self, data: pd.DataFrame, parameters: Mapping[str, float]) -> np.ndarray:
        t, d, _ = _dose_response_columns(data)
        return self._curve(t, d, *(float(parameters[n]) for n in self.parameter_names))

    def refit(self, data: pd.DataFrame) -> Fit:
        t, d, y = _dose_response_columns(data)
        start, low, high = _start_and_bounds(t, d, y)

        def residual(x):
            return self._curve(t, d, *x) - y

        parameters, converged = _least_squares(residual, self.parameter_names, start, low, high)
        rms = float(np.sqrt(np.mean(residual([parameters[n] for n in self.parameter_names]) ** 2)))
        return Fit(self.name, parameters, len(self.parameter_names), len(data), rms, converged)


@dataclass(frozen=True)
class MaturingReporter(SaturatingReporter):
    """The incumbent plus the L8 maturation state -- one extra state, one extra scalar.

    ``dR_im/dt = k(d) - (lambda + 1/tau_mat) R_im`` and ``dR_m/dt = R_im/tau_mat - lambda R_m``,
    with only the mature pool observed. Integrated in closed form, so this is a state and not
    a solver.

    Parameterised by ``tau_mat`` rather than by a maturation rate so that the nesting is
    exact: at ``tau_mat = 0`` the expression IS :class:`SaturatingReporter`, which is a point
    of the box rather than a limit at infinity, and :func:`nested_out_of_sample` can check
    the subset relation instead of being told about it.

    ``tau_mat`` is fitted here, and it is not identified by these traces: across the twelve
    committed blocks it lands between 0.108 h and 4.15 h, a 38x spread on replicates of the
    same construct. That is the measurement standing behind the revised list's experiment
    1.2 -- measure ``k_mat`` on this instrument -- and the reason this class does not present
    the fitted value as a maturation time.
    """

    name: str = "MaturingReporter"
    parameter_names: tuple[str, ...] = (
        "R0", "k_basal", "k_max", "K_dose", "lambda_", "tau_mat")

    def _curve(self, t, d, R0, k_basal, k_max, K_dose, lambda_, tau_mat):
        if tau_mat <= 0.0:
            return SaturatingReporter._curve(self, t, d, R0, k_basal, k_max, K_dose, lambda_)
        rate = 1.0 / tau_mat
        k = k_basal + k_max * d / (K_dose + d)
        decay = np.exp(-lambda_ * t)
        return R0 * decay + (k / (lambda_ + rate)) * (
            rate * (1.0 - decay) / lambda_ - decay * (1.0 - np.exp(-rate * t)))


def reporter_block(export: str, construct: str) -> pd.DataFrame:
    """One committed (plate, construct) dose ladder as ``time_h``, ``dose_mM``, ``signal``.

    Reads the committed text through :mod:`ystwin.plate.replay`, so it needs no workbook and
    no environment variable and is the same numbers on any checkout.
    """
    frame = replay.load_doses(export)
    block = frame[frame.construct == construct]
    if block.empty:
        raise KeyError(
            f"{export} carries no construct {construct!r}; it has "
            f"{', '.join(sorted(pd.unique(frame.construct)))}")
    return block.reset_index(drop=True)


def endpoint_fold_induction(data: pd.DataFrame, dose: float | None = None) -> Observable:
    """Endpoint signal at one dose over the same model's zero-dose prediction.

    The observable a reader of one of these plates actually quotes, and the one an ablation
    of the reporter chain should move: fold induction at the last read. Dimensionless, so it
    is scored relatively against :data:`REPORTER_ACTIVITY_FLOOR`.

    Args:
        data: The block both models are fitted to.
        dose: Which rung to quote. Defaults to the highest on the ladder.
    """
    _, doses, _ = _dose_response_columns(data)
    rungs = np.unique(doses)
    top = float(rungs.max()) if dose is None else float(dose)
    if top not in set(rungs.tolist()):
        raise ValueError(f"dose {top:g} is not on this ladder ({', '.join(f'{r:g}' for r in rungs)})")
    end = float(np.max(np.asarray(data.time_h, dtype=float)))
    at = pd.DataFrame({"time_h": [end, end], "dose_mM": [top, 0.0], "signal": [0.0, 0.0]})

    def summarise(model: FittableModel, fit: Fit) -> float:
        dosed, control = model.predict(at, fit.parameters)
        if abs(control) < 1e-12:
            raise ValueError(
                f"{model.name!r} predicts a zero-dose endpoint of {control:g}, so a fold "
                f"induction has no denominator. The 2 and 4 mM peroxide wells do this: "
                f"generator/panel_experiment.py records them as unusable after the dilution "
                f"correction, and a model fitted through them can land on zero basal signal")
        return float(dosed / control)

    return Observable(
        name=REPORTER_ACTIVITY_FLOOR.name,
        units="fold over the zero-dose control",
        assay=REPORTER_ACTIVITY_FLOOR.assay,
        scored="relative",
        data=data,
        summarise=summarise,
        description=f"endpoint fold induction at {top:g} mM, read at {end:g} h",
    )
