"""Simulated reporter panels across the whole stressor set.

One stressor at one dose gives a module activity vector; a dose series scales it, so each
stressor traces a ray through module space. That geometry is the point -- whether a state
learned from some stressors describes another depends on whether the rays share a
subspace, and a dataset that hides the rays could not answer the question.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pandas is imported inside `modules_frame` -- this package does not
    import pandas as pd  # need it at import time -- so the annotation resolves here.

from ..mech.state import PLATE_GROWTH_BAND_PER_H, Window
from .context import CultureContext, baseline_activity, context_growth_rate
from .kinetics import reporter_at_time
from .plate_layout import edge_multiplier, well_positions
from .stress_panel import (
    MODULES,
    REPORTERS,
    STRESSORS,
    _propagate,
    combination_response,
    growth_rate,
    healthy_ladder,
    module_response,
    reporter_loadings,
    viability,
)

__all__ = [
    "MECHANISED_CHANNELS",
    "MEASURED_ACTIVITY_CV",
    "Mechanism",
    "pool_replicates",
    "with_growth_channel",
    "dose_time_ray_shares",
    "MEASURED_GROWTH_RATE_SE",
    "OBSERVED_ACTIVITY_CV",
    "PanelDataset",
    "panel_dataset",
]

OBSERVED_ACTIVITY_CV = 0.146
"""Total plate-to-plate CV of dilution-corrected activity, across 28 matched conditions on
the 2026-07-22 and 2026-08-03 plates. What a simulated plate has to reproduce, and what the
posterior-predictive gate checks against."""

MEASURED_GROWTH_RATE_SE = 0.0117
"""Standard error of an estimated specific growth rate, in 1/h.

Derived, not fitted. A growth rate is the slope of log optical density against time, so its
error is the standard error of that slope -- residual spread over the spread of the time
points and the root of their number, every term measurable on the trace. Median across 210
real wells with a fittable window, IQR 0.0079 to 0.0200, median window 25 points.

It is an absolute rate error, which is the whole point. The same 0.0117 is 3% of a healthy
growth rate and 12% of one slowed to 0.10, so the dose structure of the confound follows
from the parameterisation instead of being dialled in. An earlier version expressed it as a
CV times the control rate and set it by grid search, which came to 0.048 -- four times too
large, and unidentifiable besides."""

_DEFAULT_RUNGS = 5

MEASURED_ACTIVITY_CV = 0.14
"""Residual measurement CV on a reading, once the growth-correction error is accounted for.

The two are not co-fitted. This one is bounded by the observed total above, the growth term
is measured independently from the traces, and what the pair produce is then checked
against that total rather than solved for it."""


MECHANISED_CHANNELS = {
    ("H2O2", "oxidative"): "TRX2-oxidative",
    ("DTT", "UPR"): "UPRE-ER",
}
"""The ``(stressor, module)`` pairs a mech block drives, and the reporter each one reaches.

Two of the 25 x 24 stressor-by-module cells, and they are the only two whose output this
project can put beside a measurement: `mech/oxidative.py` is scored on the NativeYap1 and
AlteredYap1 H2O2 ladders and `mech/upr.py` on the UPRE1 and UPRE2 DTT ladders, which is four
of the panel's 600 reporter x stressor cells.

Nothing else is routed, and each omission has a reason rather than a schedule.
`mech/ph.py` is not a promoter-activity route in this panel. Its charge-balance repair
has landed: cytosolic pH is reducible to algebra on the audited windows, and the old dynamic
reporter-ablation claim was retracted. The separate `mech/chain.py` uses its weak-acid bounds
and Citrine quench; that optical effect is not a pH-regulated promoter here. `mech/burden.py` and
`mech/population.py` produce no promoter activity at all -- they move growth and the bearing
fraction -- so there is no channel for them to drive here. Tunicamycin is excluded although
it shares the UPR module: `mech/upr.py::TUNICAMYCIN_ENTRY` is REFUSED with no declared axis
to stand in for it, while DTT has :class:`~ystwin.mech.upr.DoseEntry`. Menadione, diamide,
cobalt and copper are excluded although they share the oxidative module: the block's dose
law is catalase and peroxiredoxin consumption of extracellular H2O2, which is not how a
redox cycler or a metal reaches Yap1.
"""

def _blocks():
    """The two mech blocks that drive a reporter channel, imported on first use.

    Deferred rather than declared at the top of this file, and for a real cycle rather than
    a preference: `mech/ablation.py` and `mech/integrate.py` both import this module's two
    MEASURED floors, so `mech/oxidative.py` and `mech/upr.py` cannot be named at import time
    here without closing the loop. `mech/state.py` has no such edge and is imported normally
    above. The floors are properties of a plate assay and belong beside the assay rather than
    in the generator; moving them is a change to files this seam does not own.
    """
    from ..mech import oxidative, upr

    return oxidative, upr


_MECHANISM_GRID_POINTS = 801
"""Points the mechanistic driver and the reporter convolution are evaluated on.

ASSERTED as a numerical setting, in the same sense as `mech/upr.py::INTERNAL_GRID_POINTS`,
and pinned by measurement rather than by argument: tests/test_panel_mechanistic.py doubles
it and checks the readings move by far less than :data:`OBSERVED_ACTIVITY_CV`. Over an 8 h
schedule it is 0.01 h of spacing against a pipetted-H2O2 half-life of 0.20-0.69 h at the
plate inoculum, so the fastest thing on the grid is resolved 20-fold.
"""


@dataclass(frozen=True)
class Mechanism:
    """One declared corner of the mechanistic path, carrying every axis the blocks refuse.

    There is no default constructor value for any of the six scalars, and that is the point:
    `mech/oxidative.py::CELLS_PER_ML_PER_OD600` and `mech/upr.py::FOLDING_RATE_PER_H` are
    SWEPT, `mech/upr.py::DTT_ENTRY` is REFUSED and reachable only through a declared
    :class:`~ystwin.mech.upr.DoseEntry`, and the remaining three are the fitted scalars of
    `mech/upr.py::PlateUpr`, two of which sit on their own box. A number produced through
    this class therefore cannot be quoted without the corner that produced it, which is the
    same discipline `mech/chain.py::SweepPoint` enforces on the product chain.

    Args:
        cells_per_ml_per_od600: A point of the SWEPT band in `mech/oxidative.py`. It sets how
            fast the well consumes the pipetted H2O2 and is the widest axis here.
        folding_rate_per_h: A point of `mech/upr.py::FOLDING_RATE_PER_H`.
        basal_load: Unstressed ER client influx as a fraction of folding capacity.
        top_load: The same quantity at ``top_dose_mM``. This IS the declared dose axis
            `mech/upr.py::DTT_ENTRY` refuses to supply.
        hill_n: Ire1 cooperativity, inside Korennykh 2009's measured bracket.
        basal_share: UPRE-independent share of the promoter rate.
        od600: Starting optical density of the well, for the H2O2 consumption law.
        top_dose_mM: The dose ``top_load`` refers to.
        grid_points: Overrides :data:`_MECHANISM_GRID_POINTS` for a convergence check.
    """

    cells_per_ml_per_od600: float
    folding_rate_per_h: float
    basal_load: float
    top_load: float
    hill_n: float
    basal_share: float
    od600: float
    top_dose_mM: float = 1.0
    grid_points: int = _MECHANISM_GRID_POINTS

    def __post_init__(self) -> None:
        oxidative, upr = _blocks()
        oxidative.CELLS_PER_ML_PER_OD600.at(float(self.cells_per_ml_per_od600))
        upr.UPRE_BASAL_SHARE.at(float(self.basal_share))
        if self.od600 <= 0.0:
            raise ValueError(
                f"od600 sets the biomass that eats the dose and must be positive, got "
                f"{self.od600!r}")
        if not 0.0 < self.basal_load < self.top_load < 1.0:
            raise ValueError(
                f"the ER load must satisfy 0 < basal_load < top_load < 1 -- at 1 the folding "
                f"machinery is exactly saturated and the client pool diverges. Got "
                f"basal_load={self.basal_load!r}, top_load={self.top_load!r}")
        if int(self.grid_points) < 3:
            raise ValueError(f"a convolution grid needs at least 3 points, got {self.grid_points!r}")
        self.upr_entry()
        self.upr_chain(float(oxidative.PLATE_MEDIAN_GROWTH_PER_H))

    def upr_entry(self):
        """The DECLARED mM-to-ER-load axis. `mech/upr.py` refuses to supply it."""
        return _blocks()[1].DoseEntry(
            top_load=float(self.top_load), top_dose_mM=float(self.top_dose_mM), declared=True,
            note="declared by generator/panel_experiment.py::Mechanism; the mM-to-lumen "
                 "mapping stays REFUSED and this is a sensitivity axis, not a measurement")

    def upr_chain(self, mu_per_h: float):
        """The six-state chain at one well's MEASURED growth rate."""
        return _blocks()[1].UprChain(
            basal_load=float(self.basal_load), hill_n=float(self.hill_n),
            basal_share=float(self.basal_share),
            folding_rate_per_h=float(self.folding_rate_per_h), mu_per_h=float(mu_per_h))

    @classmethod
    def corners(cls) -> tuple["Mechanism", ...]:
        """The four corners of the two SWEPT axes, at the blocks' own committed scalars.

        Corners rather than a grid, for the reason `mech/chain.py::SweepPoint.corners` gives:
        the H2O2 half-life is monotone in the cell density and the UPR onset is monotone in
        the folding rate, so the extremes of anything computed here are attained at a corner.

        The four non-swept scalars are `mech/upr.py::PlateUpr`'s committed fit on the
        20260722 UPRE1 fold -- ``basal_load`` 0.2424, ``top_load`` 0.99, ``hill_n`` 8.0,
        ``basal_share`` 1.8093 -- and its docstring says plainly that the middle two sit on
        their upper bounds because twelve wells at four doses do not identify four scalars.
        The inoculum is `mech/oxidative.py::PLATE_INOCULUM_OD600`, MEASURED on this project's
        own plates.
        """
        oxidative, upr = _blocks()
        cells = oxidative.CELLS_PER_ML_PER_OD600.bounds
        folding = upr.FOLDING_RATE_PER_H.bounds
        return tuple(cls(cells_per_ml_per_od600=c, folding_rate_per_h=f,
                         od600=float(oxidative.PLATE_INOCULUM_OD600), **_COMMITTED_UPR_FIT)
                     for c in cells for f in folding)

    @classmethod
    def midpoint(cls) -> "Mechanism":
        """The centre of both swept axes. For a smoke test, never for a quoted number."""
        oxidative, upr = _blocks()
        cells = oxidative.CELLS_PER_ML_PER_OD600.bounds
        folding = upr.FOLDING_RATE_PER_H.bounds
        return cls(cells_per_ml_per_od600=0.5 * (cells[0] + cells[1]),
                   folding_rate_per_h=0.5 * (folding[0] + folding[1]),
                   od600=float(oxidative.PLATE_INOCULUM_OD600), **_COMMITTED_UPR_FIT)


_COMMITTED_UPR_FIT = {
    "basal_load": 0.2424,
    "top_load": 0.99,
    "hill_n": 8.0,
    "basal_share": 1.8093,
}
"""`mech/upr.py::PlateUpr`'s committed fit on the 20260722 UPRE1 fold.

Restated rather than imported because that module carries these four in its class docstring
and in tests/test_mech_upr.py rather than as data. tests/test_panel_mechanistic.py pins them
against both so the two copies cannot drift apart, and reports what that module already
says: ``top_load`` and ``hill_n`` sit on their own upper bounds, because twelve wells at
four doses do not identify four scalars."""


def _cascade_column(module: str) -> np.ndarray:
    """Where a unit of one module's DIRECT activity lands, after `stress_panel`'s cascade.

    The mechanistic override replaces a direct activity, so it has to travel the same
    cascade the Hill it replaces travelled -- oxidative reaches the proteasome at 0.45 and
    iron at 0.1, and the proteasome carries 0.2 of that into xenobiotic. `_propagate` is
    linear in the direct vector, so one call on a unit vector is the whole column and the
    override is exact rather than a re-implementation of the cascade.
    """
    unit = {name: (1.0 if name == module else 0.0) for name in MODULES}
    carried = _propagate(unit)
    return np.array([carried[name] for name in MODULES])


def _mechanistic_direct(spec: Mechanism, stressor: str, dose: float, mu: float,
                        grid: np.ndarray) -> tuple[str, np.ndarray] | None:
    """One block's direct module activity over the read, or ``None`` if it drives nothing.

    Both arms return the SAME quantity `stress_panel.module_response` computes -- the
    stressor's target weight times an activation in [0, 1] times the panel's own viability
    envelope -- and differ from it in one thing: the activation is a function of time.
    Viability is kept because it is the plate's measured lethality and neither block models
    it, so the two paths differ only along the time axis.
    """
    module = next((m for (s, m) in MECHANISED_CHANNELS if s == stressor), None)
    if module is None or dose <= 0.0:
        return None
    weight = STRESSORS[stressor].targets[module]
    alive = viability(stressor, dose)
    oxidative, upr = _blocks()
    if stressor == "H2O2":
        k0 = oxidative.dose_decay_rate_per_h(float(spec.od600),
                                             float(spec.cells_per_ml_per_od600))
        pool = oxidative.extracellular_h2o2(grid, float(dose), k0, mu)
        activation = np.asarray(oxidative.yap1_nuclear_fraction(pool), dtype=float)
        return module, weight * activation * alive
    lethal = float(upr.DTT_LETHAL_MM)
    if dose > lethal:
        raise ValueError(
            f"{dose:g} mM DTT is above mech/upr.py's MEASURED lethal dose of {lethal:g} mM, "
            f"where growth has halved and the block declares itself inadmissible. The "
            f"algebraic path has no such refusal; pass mechanism=None to run this dose")
    window = Window(
        name="PANEL_READ", duration_h=float(grid[-1]),
        growth_rate_low_per_h=PLATE_GROWTH_BAND_PER_H[0],
        growth_rate_high_per_h=PLATE_GROWTH_BAND_PER_H[1],
        source="the panel's own read schedule, over mech/state.py's MEASURED plate band")
    chain = spec.upr_chain(mu)
    trajectory = chain.simulate(spec.upr_entry(), float(dose), window,
                                n_points=int(spec.grid_points))
    occupancy = upr.UprChain.occupancy(trajectory)
    resting = float(chain.basal()["theta0"])
    span = float(upr.HAC1_MAX_OCCUPANCY) - resting
    activation = np.clip((occupancy - resting) / span, 0.0, 1.0)
    return module, weight * np.interp(grid, trajectory.t, activation) * alive


def _relaxed(delta: np.ndarray, relaxation: float, grid: np.ndarray,
             times: np.ndarray) -> np.ndarray:
    """``relaxation * int_0^t delta(s) e^{-relaxation (t - s)} ds``, at each requested time.

    A promoter fusion obeys ``dR/dt = k - (mu + k_deg) R``, so what it reports is its own
    promoter activity through this one-pole filter and not the instantaneous value. The
    incumbent's `kinetics.reporter_at_time` is the closed form of the same filter for an
    activity that steps once and never moves; this is the same filter for one that does.
    It is applied to the DIFFERENCE from the algebraic path so that every module the
    mechanism does not touch keeps the incumbent's exact closed form.
    """
    step = float(grid[1] - grid[0])
    decay = float(np.exp(-relaxation * step))
    accumulated = np.zeros_like(delta)
    half = 0.5 * relaxation * step
    for index in range(1, len(grid)):
        accumulated[index] = (accumulated[index - 1] * decay
                              + half * (delta[index - 1] * decay + delta[index]))
    return np.array([np.interp(times, grid, accumulated[:, column])
                     for column in range(delta.shape[1])]).T


def _mechanistic_readings(spec: Mechanism, members, applied_doses, activity, floor,
                          mu: float, k_deg: float, vector: np.ndarray, basal: np.ndarray,
                          loadings: np.ndarray, integrates: np.ndarray,
                          grid: np.ndarray, times: np.ndarray) -> np.ndarray | None:
    """Readings at every requested time, or ``None`` when no block drives this treatment.

    ``None`` is not a fallback, it is the guarantee: a treatment no mech block touches never
    enters the quadrature at all and comes out of the incumbent's own closed form, so the
    mechanistic path is byte-identical to the algebraic one everywhere off
    :data:`MECHANISED_CHANNELS`. Where a block does drive, only the DIFFERENCE from the
    algebraic activity goes through :func:`_relaxed`, so the untouched modules keep the
    closed form there too and the quadrature error is confined to the driven column and its
    cascade.
    """
    order = list(MODULES)
    delta = np.zeros((len(grid), len(order)))
    driven_any = False
    for member in members:
        driven = _mechanistic_direct(spec, member, applied_doses[member], mu, grid)
        if driven is None:
            continue
        module, activation = driven
        if MODULES[module].driven_by:
            raise RuntimeError(
                f"module {module!r} is now driven by {sorted(MODULES[module].driven_by)}, so "
                f"stress_panel.module_response no longer reports its DIRECT activity and the "
                f"override in MECHANISED_CHANNELS would double-count the cascade")
        incumbent = module_response(member, applied_doses[member])[module]
        delta = delta + (activation - incumbent)[:, None] * _cascade_column(module)[None, :]
        driven_any = True
    if not driven_any:
        return None
    unclipped = np.array([activity[name] + floor[name] for name in order])
    trace = np.minimum(1.0, unclipped[None, :] + delta)
    relaxation = mu + k_deg
    filtered = (vector[None, :] * (1.0 - np.exp(-relaxation * times))[:, None]
                + _relaxed(trace - vector[None, :], relaxation, grid, times))
    instant = np.array([np.interp(times, grid, trace[:, column])
                        for column in range(trace.shape[1])]).T
    return np.where(integrates[None, :] > 0,
                    basal[None, :] + filtered @ loadings.T,
                    basal[None, :] + instant @ loadings.T)


def dose_time_ray_shares(dataset: PanelDataset) -> dict[str, float]:
    """Leading-ray share of each stressor's centred dose x time readings.

    ``sigma_1^2 / sum(sigma_i^2)`` of the SVD of one stressor's readings after subtracting
    the column mean. It is 1 when every well of that stressor differs from the mean in one
    fixed direction -- when the whole dose ladder and the whole read schedule between them
    trace a single line in reading space, which is the geometry a latent state cannot learn
    from. Centred, unlike `mech/chain.py::ray_share`, which asks the same question about one
    reporter's uncentred surface; that one measures the surface, this one measures the
    deviations a fit would see.
    """
    shares = {}
    for label in dict.fromkeys(dataset.labels.tolist()):
        block = dataset.readings[dataset.labels == label]
        if block.shape[0] < 2:
            continue
        singular = np.linalg.svd(block - block.mean(axis=0), compute_uv=False)
        total = float(np.sum(singular ** 2))
        if total <= 0.0:
            continue
        shares[label] = float(singular[0] ** 2 / total)
    return shares


@dataclass(frozen=True)
class PanelDataset:
    """Reporter readings with the labels and hidden truth behind them."""

    readings: np.ndarray
    labels: np.ndarray
    doses: np.ndarray
    modules: np.ndarray
    reporters: list[str]
    wells: list[str] = field(default_factory=list)
    contexts: list[str] = field(default_factory=list)
    read_times: list[float] = field(default_factory=list)
    growth_rates: list[float] = field(default_factory=list)
    n_unusable: int = 0
    """Wells dropped because the dilution correction could not be applied to them.

    Activity is recovered by dividing growth out, so a culture that is barely growing gives
    an estimate that can land below zero -- not a weak promoter but an arithmetic artifact,
    and exactly what the real plates showed at 2 mM peroxide. An experimenter drops those
    wells; keeping them would teach a model that promoters run backwards in stationary
    phase."""

    def mask(self, stressor: str) -> np.ndarray:
        return self.labels == stressor

    def modules_frame(self) -> "pd.DataFrame":
        """True module activities as a named table, for scoring against a prediction."""
        import pandas as pd

        return pd.DataFrame(self.modules, columns=list(MODULES))


def panel_dataset(
    reporters=None,
    stressors=None,
    doses=None,
    replicates: int = 3,
    noise_cv: float = 0.05,
    seed: int = 0,
    combinations=(),
    growth_rate_se: float = 0.0,
    k_deg: float = 0.0,
    edge_effect: float = 0.0,
    rows: int = 8,
    columns: int = 12,
    contexts=None,
    read_times_h=None,
    mechanism: "Mechanism | None" = None,
) -> PanelDataset:
    """Simulate a reporter panel over stressors and doses.

    Args:
        reporters: Reporter names to read; defaults to the whole panel.
        stressors: Stressors to apply; defaults to all of them.
        doses: Multiples of each stressor's EC50 to apply. Left unset, each agent gets a
            ladder spanning up to a fraction of its own lethal dose, so no well is dosed
            past the point where a culture stops reporting.
        replicates: Independent wells per stressor and dose.
        noise_cv: Multiplicative measurement noise on each reading.
        seed: Seed for the noise.
        combinations: Stressor tuples to co-apply, which visit the directions between the
            single-agent rays a blocked panel never leaves.
        growth_rate_se: Standard error of the estimated growth rate, in 1/h. Activity is
            recovered as ``k = dR/dt + mu*R``, so an error in mu propagates as
            ``k_hat = k * (1 + eps / (mu + k_deg))``. Because it is absolute, the same error
            is a small relative one in a healthy culture and a large one in a slowed
            culture, and the dose structure follows without being imposed. Measured at
            0.0117 on real traces; zero gives idealised activity.
        k_deg: Reporter loss beyond dilution, which damps the confound. Zero for a stable
            fluorescent protein, which is what these constructs use.
        edge_effect: Relative lift on wells in the outer ring, where evaporation
            concentrates the medium. Unlike noise it does not average out over replicates,
            because every replicate in that ring carries the same bias.
        rows: Plate rows, for assigning wells and finding the edge.
        columns: Plate columns.
        contexts: Culture conditions to run the whole design in. Each shifts the baseline
            the stressor is added to and the rate the culture grows at, so a stationary well
            carries a raised general stress response with no agent in it. Defaults to one
            exponential glucose culture, which is the corner everything was trained in.
        mechanism: A declared :class:`Mechanism` corner, or ``None`` for the algebraic path.
            OFF IS THE DEFAULT AND MUST STAY THE DEFAULT. Every committed table and every
            pinned number in this repository was computed on the algebraic path, so flipping
            the default would move all of them at once. With a corner given, the two
            ``(stressor, module)`` pairs of :data:`MECHANISED_CHANNELS` are driven by
            `mech/oxidative.py` and `mech/upr.py` instead of by
            `stress_panel.module_response`'s instantaneous Hill, and every other cell of the
            panel is untouched -- byte-identical, not merely close. It requires
            ``read_times_h``: a mechanistic reading is a time course, and the algebraic
            path's timeless steady state is exactly what the blocks say does not exist.

    Raises:
        ValueError: if ``mechanism`` is given without a finite, positive read schedule.
    """
    reporters = list(reporters or REPORTERS)
    stressors = list(stressors or STRESSORS)
    unknown = [s for s in stressors if s not in STRESSORS]
    if unknown:
        raise KeyError(f"no such stressor(s): {sorted(unknown)}")
    loadings = reporter_loadings(reporters)
    basal = np.array([REPORTERS[r].basal for r in reporters])
    order = list(MODULES)
    rng = np.random.default_rng(seed)

    treatments = [(name, (name,)) for name in stressors]
    for group in combinations:
        members = tuple(group)
        unknown = [m for m in members if m not in STRESSORS]
        if unknown:
            raise KeyError(f"no stressor(s) {sorted(unknown)}; have {sorted(STRESSORS)}")
        treatments.append(("+".join(members), members))

    contexts = list(contexts) if contexts else [CultureContext()]
    schedule = list(read_times_h) if read_times_h else [None]
    integrates = np.array([REPORTERS[r].integrates_hours for r in reporters])
    mech_times = mech_grid = None
    if mechanism is not None:
        if any(read_at is None for read_at in schedule):
            raise ValueError(
                "the mechanistic path needs read_times_h: mech/oxidative.py's dose is "
                "consumed within the read and mech/upr.py's occupancy is a pulse, so "
                "neither has the timeless steady state the algebraic path reads")
        mech_times = np.array(schedule, dtype=float)
        if not np.all(np.isfinite(mech_times)) or mech_times.min() < 0.0 or mech_times.max() <= 0.0:
            raise ValueError(f"read times must be finite and non-negative, got {schedule}")
        mech_grid = np.linspace(0.0, float(mech_times.max()), int(mechanism.grid_points))
    rows_out, labels, applied, truth, context_names, when = [], [], [], [], [], []
    growth_out = []
    dropped = 0
    for context in contexts:
        name = f"{context.growth_phase}/{context.carbon_source}"
        floor = baseline_activity(context)
        unstressed = context_growth_rate(context)
        for label, members in treatments:
            rungs = (list(doses) if doses is not None
                     else list(healthy_ladder(members[0], _DEFAULT_RUNGS)
                               / STRESSORS[members[0]].ec50))
            for multiple in rungs:
                applied_doses = {m: multiple * STRESSORS[m].ec50 for m in members}
                activity = combination_response(applied_doses)
                vector = np.array([min(1.0, activity[m] + floor[m]) for m in order])
                steady = basal + loadings @ vector
                slowdown = min((growth_rate(m, applied_doses[m]) / growth_rate(m, 0.0)
                                for m in members), default=1.0)
                mu = unstressed * slowdown
                mech_rows = (None if mechanism is None else _mechanistic_readings(
                    mechanism, members, applied_doses, activity, floor, mu, k_deg, vector,
                    basal, loadings, integrates, mech_grid, mech_times))
                for at, read_at in enumerate(schedule):
                  clean = steady if read_at is None else np.array([
                      reporter_at_time(read_at, steady[i], mu, initial=basal[i], k_deg=k_deg)
                      if integrates[i] > 0 else steady[i]
                      for i in range(len(reporters))])
                  if mech_rows is not None:
                      clean = mech_rows[at]
                  for _ in range(replicates):
                    noise = rng.normal(1.0, noise_cv, size=clean.shape) if noise_cv else 1.0
                    reading = clean * noise
                    if growth_rate_se:
                        error = rng.normal(0.0, growth_rate_se)
                        correction = 1.0 + error / max(mu + k_deg, 1e-6)
                        if correction <= 0.0:
                            dropped += 1
                            continue
                        reading = reading * correction
                    rows_out.append(reading)
                    labels.append(label)
                    applied.append(multiple * STRESSORS[members[0]].ec50)
                    truth.append(vector)
                    context_names.append(name)
                    growth_out.append(float(mu))
                    when.append(float("inf") if read_at is None else float(read_at))
    wells = well_positions(len(rows_out), rows, columns) if len(rows_out) <= rows * columns \
        else [f"W{i}" for i in range(len(rows_out))]
    readings = np.array(rows_out)
    if edge_effect:
        readings = readings * edge_multiplier(wells, edge_effect, rows, columns)[:, None]
    return PanelDataset(
        readings=readings,
        labels=np.array(labels),
        doses=np.array(applied),
        modules=np.array(truth),
        reporters=reporters,
        wells=wells,
        contexts=context_names,
        read_times=when,
        growth_rates=growth_out,
        n_unusable=dropped,
    )


def with_growth_channel(dataset: PanelDataset) -> PanelDataset:
    """Add the growth rate as an extra observed channel.

    Three stress channels cannot separate culture context from stress -- a high general
    response in stationary phase looks like a stressed exponential one. Growth tells them
    apart and costs no fluorophore, because optical density is measured anyway.
    """
    readings = np.column_stack([dataset.readings, np.asarray(dataset.growth_rates, dtype=float)])
    return replace(dataset, readings=readings,
                   reporters=list(dataset.reporters) + ["growth"])


def pool_replicates(dataset: PanelDataset) -> PanelDataset:
    """Average replicate wells of a condition before fitting.

    The latent model gives every well its own state, so replicates add parameters rather
    than averaging noise and more wells buy nothing. Pooling first turns the same plate into
    roughly twice as many recoverable modules.
    """
    groups: dict = {}
    for index, key in enumerate(zip(dataset.labels, dataset.doses,
                                    dataset.contexts or [""] * len(dataset.labels))):
        groups.setdefault(key, []).append(index)
    readings, labels, doses, modules, contexts = [], [], [], [], []
    for (label, dose, context), rows in groups.items():
        readings.append(dataset.readings[rows].mean(axis=0))
        modules.append(dataset.modules[rows].mean(axis=0))
        labels.append(label)
        doses.append(dose)
        contexts.append(context)
    return replace(dataset, readings=np.array(readings), labels=np.array(labels),
                   doses=np.array(doses), modules=np.array(modules), contexts=contexts,
                   wells=[], read_times=[], growth_rates=[])
