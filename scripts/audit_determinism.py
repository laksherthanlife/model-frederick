"""Check that seeded code gives the same answer twice, and a different one for a new seed.

Almost every number this repository reports passes through a random number generator: a
cluster bootstrap for the fold changes, a Monte-Carlo null for the skill scores, a
cross-validated split for the latent dimension, a particle filter for the state estimate,
a simulated campaign for the power analysis. All of them take a ``seed``, and the tables in
``outputs/`` are only citable if that seed does what it appears to do.

Two things have to hold, and neither is checked by the test suite, which asserts that
results fall in a range rather than that they are the same result.

**Same seed, same answer.** Otherwise a table cannot be regenerated, and a diff against a
re-run reports noise as change. The usual causes are a global ``np.random`` call somewhere
in the chain, a set iterated in hash order, or a dictionary keyed by an object whose
``id()`` leaks into the ordering.

**Different seed, different answer.** This is the check people leave out, and it is the one
that catches a seed threaded halfway. A function that accepts ``seed``, passes it to a
generator, and then never draws from it is indistinguishable from a correct one -- both
reproduce -- until someone tries to compute a spread across seeds and discovers every
replicate is the same number.

The second property is not universal, and pretending it is would make the audit lie. A
function returning a *selection* -- a latent width, an integer argmin -- can legitimately
return the same choice under every seed, because the seed moves the scores and not the
winner. Those specimens declare ``seed_sensitive=False`` with a reason, and the audit then
fails them if they turn out to be sensitive after all. The expectation is asserted in both
directions so that declaring it cannot be used to silence a finding.

Entry points are discovered by inspecting signatures across ``ystwin`` rather than listed
here, so a seeded function added tomorrow shows up -- as a coverage failure until someone
gives it a synthetic input, which is the intended nag. Anything needing real wet-lab data
or a genome-scale model is exempt by name: it cannot run on a machine without them, and
this audit has to run everywhere.

Emits ``outputs/audit_determinism.csv`` and exits non-zero if any check fails.

Usage: python scripts/audit_determinism.py [--quiet] [--seed N]
"""
from __future__ import annotations

import argparse
import dataclasses
import functools
import hashlib
import importlib
import inspect
import pathlib
import pkgutil
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from ystwin.readings import CorrectedOD, SpecificFluorescence
import pandas as pd

__all__ = [
    "REPO",
    "SPECIMENS",
    "Specimen",
    "IMPORT_FAILED",
    "check_coverage",
    "check_specimen",
    "digest",
    "discover_seeded",
]

REPO = pathlib.Path(__file__).resolve().parents[1]

#: Seeds tried when asking whether the seed matters. Several, because a coarse output --
#: a hit count out of forty simulations, an integer width -- can collide by chance for one
#: alternative seed without the seed being ignored. One difference anywhere proves the
#: seed is threaded; sameness across all of them is the finding.
ALTERNATIVE_SEEDS = (1, 2, 3)

#: Stands in for the name of a callable in a module that would not import, so a machine
#: missing an optional dependency reports the module rather than pretending it holds
#: nothing seeded. Marked in the key, not by a sentinel value, so "did not import" and
#: "signature unavailable" cannot be confused for each other.
IMPORT_FAILED = "<import failed>"


def _row(check: str, expected: object, actual: object, status: object, detail: str = "") -> dict:
    """One CSV row, in the shape ``audit_claims.py`` established.

    ``status`` may be a bool or the literal ``"SKIP"`` for a specimen that could not run
    here -- an absent optional dependency, say. A skip does not fail the audit and is not
    a pass: nothing was demonstrated.
    """
    text = status if isinstance(status, str) else ("PASS" if status else "FAIL")
    return {
        "check": check,
        "expected": str(expected),
        "actual": str(actual),
        "status": text,
        "detail": detail,
    }


# ---------------------------------------------------------------------------
# comparing two results exactly
# ---------------------------------------------------------------------------


def _canonical(obj, depth: int = 0) -> str:
    """A byte-exact, order-stable rendering of whatever a seeded function returned.

    Bit patterns rather than printed digits: two bootstrap runs that agree to fifteen
    decimals and differ in the sixteenth are *not* reproducible, and a formatted
    comparison would call them equal. NaN compares equal to NaN here, which is right --
    the question is whether the same computation happened, not whether the answer is a
    number.
    """
    if depth > 12:
        return "<deep>"
    if obj is None or isinstance(obj, (bool, int)) and not isinstance(obj, np.generic):
        return repr(obj)
    if isinstance(obj, str):
        return f"str:{obj}"
    if type(obj) is float:
        return f"f:{float.hex(obj)}"
    if isinstance(obj, np.generic):
        return f"np:{obj.dtype}:{_canonical(obj.item(), depth + 1)}"
    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            return "obj_array:" + "|".join(_canonical(v, depth + 1) for v in obj.ravel())
        return f"arr:{obj.shape}:{obj.dtype}:{obj.tobytes().hex()}"
    if isinstance(obj, pd.DataFrame):
        return ("df:" + _canonical(list(obj.columns), depth + 1)
                + _canonical(obj.index.to_numpy(), depth + 1)
                + "|".join(_canonical(obj[c].to_numpy(), depth + 1) for c in obj.columns))
    if isinstance(obj, pd.Series):
        return "series:" + _canonical(obj.to_numpy(), depth + 1)
    if isinstance(obj, pathlib.PurePath):
        return f"path:{obj}"
    if isinstance(obj, dict):
        return "{" + ",".join(f"{_canonical(k, depth + 1)}:{_canonical(v, depth + 1)}"
                              for k, v in sorted(obj.items(), key=lambda kv: repr(kv[0]))) + "}"
    if isinstance(obj, (list, tuple)):
        return "[" + ",".join(_canonical(v, depth + 1) for v in obj) + "]"
    if isinstance(obj, (set, frozenset)):
        return "{" + ",".join(sorted(_canonical(v, depth + 1) for v in obj)) + "}"
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return (f"{type(obj).__name__}("
                + ",".join(f"{f.name}={_canonical(getattr(obj, f.name), depth + 1)}"
                           for f in dataclasses.fields(obj)) + ")")
    if isinstance(obj, float):
        # A float subclass -- PowerEstimate -- carries both a value and its provenance.
        state = getattr(obj, "__dict__", {})
        return f"{type(obj).__name__}({float.hex(float(obj))},{_canonical(dict(state), depth + 1)})"
    state = getattr(obj, "__dict__", None)
    if state is not None:
        return f"{type(obj).__name__}({_canonical(dict(state), depth + 1)})"
    return f"{type(obj).__name__}:{obj!r}"


def digest(obj) -> str:
    """A short hash of :func:`_canonical`, so the CSV can show what differed without the data."""
    return hashlib.blake2b(_canonical(obj).encode("utf-8"), digest_size=8).hexdigest()


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------


def discover_seeded(package) -> dict[str, inspect.Signature]:
    """Every public callable in ``package`` that takes a ``seed`` argument.

    Discovery rather than a list, so the audit keeps covering the package as it grows and
    reports the gap when it does not. Only objects a module actually defines are
    considered, so a function re-exported through three ``__init__`` files is one entry
    and not four.

    Returns:
        ``{"ystwin.analysis.latent::fit_latent": <Signature>, ...}``.
    """
    found: dict[str, inspect.Signature] = {}
    for module_info in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
        try:
            module = importlib.import_module(module_info.name)
        except Exception:
            # An optional dependency this machine lacks. Reported by check_coverage as an
            # unreachable module rather than silently dropped.
            found[f"{module_info.name}::{IMPORT_FAILED}"] = None
            continue
        for name, obj in vars(module).items():
            if name.startswith("_"):
                continue
            if getattr(obj, "__module__", None) != module_info.name:
                continue
            if not (inspect.isfunction(obj) or inspect.isclass(obj)):
                continue
            try:
                signature = inspect.signature(obj)
            except (ValueError, TypeError):
                continue
            if "seed" in signature.parameters:
                found[f"{module_info.name}::{name}"] = signature
    return found


#: Seeded entry points this audit deliberately does not exercise, and why. Every one needs
#: something the machine may not have, so covering them would make the audit's own result
#: depend on whether the wet-lab data is present -- which is the opposite of the point.
EXEMPT: dict[str, str] = {}


def _plate_geometry():
    """The committed plates' sampling design, for the estimator-accuracy specimens."""
    from ystwin.analysis.estimator_accuracy import PLATE

    return PLATE


_PLATE_GEOMETRY = _plate_geometry()


def _sensor_readings():
    """The committed per-well activity table, for the dose-response permutation."""
    import pandas as pd

    from ystwin import paths

    return pd.read_csv(paths.outputs_dir() / "sensor_characterisation.csv")


_CROSSTALK_READINGS = _sensor_readings()


# ---------------------------------------------------------------------------
# synthetic inputs
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=None)
def _activities() -> np.ndarray:
    """A small low-rank matrix with one missing reading, for the latent fits."""
    rng = np.random.default_rng(20260826)
    states = rng.normal(size=(24, 2))
    loadings = rng.normal(size=(2, 5))
    data = states @ loadings + 0.05 * rng.normal(size=(24, 5))
    data[3, 1] = np.nan
    return data


@functools.lru_cache(maxsize=None)
def _panel():
    """A three-stressor, five-reporter panel. Built once, with its own fixed seed.

    The dataset must not move when the seed under test moves, or the two would be
    confounded and every specimen would look seed-sensitive for the wrong reason.
    """
    from ystwin.generator.panel_experiment import panel_dataset
    from ystwin.generator.stress_panel import transcriptional_reporters

    return panel_dataset(reporters=transcriptional_reporters()[:5],
                         stressors=["DTT", "H2O2", "NaCl"], doses=(0.5, 1.0),
                         replicates=2, noise_cv=0.03, seed=20260826)


@functools.lru_cache(maxsize=None)
def _readings() -> pd.DataFrame:
    """Four plates of a dosed and a control condition, the shape ``fold_change`` wants."""
    rng = np.random.default_rng(20260826)
    rows = []
    for plate in range(4):
        shift = float(np.exp(rng.normal(0.0, 0.15)))
        for dose, level in ((0.0, 1.0), (1.0, 1.4)):
            for well in range(3):
                rows.append({
                    "plate": f"P{plate}", "well": f"A{well}", "construct": "R",
                    "dose_mM": dose,
                    "activity_late": shift * level * float(np.exp(rng.normal(0.0, 0.08))),
                })
    return pd.DataFrame(rows)


@functools.lru_cache(maxsize=None)
def _trace():
    """A short reporter trace with its optical density and growth rate.

    Unpacked with `*` into `activity_uncertainty`, so the readings are typed here. This
    specimen is reached through a resolved `target(...)` handle rather than by name, which
    is why the 2026-08-31 migration's call-site rewriter could not see it and why this
    audit caught it instead of the test suite.
    """
    t = np.linspace(0.0, 8.0, 40)
    mu = np.full_like(t, 0.28)
    od = 0.1 * np.exp(0.28 * t)
    reporter = 1.2e-3 / (0.28 + 0.05) * (1.0 - np.exp(-(0.28 + 0.05) * t)) + 1e-4
    return t, SpecificFluorescence(reporter), CorrectedOD(od), mu


@functools.lru_cache(maxsize=None)
def _design():
    """Literature kinetics and a short dose series, for the power simulations."""
    from ystwin.generator.design import stressor_only_series
    from ystwin.generator.literature import parameters_from_literature

    params = parameters_from_literature("UPRE2")
    return params, stressor_only_series(params, doses=(0.0, 0.5, 1.0))


def _split_frame():
    """A manifest four plates wide, consistent with qpcr.STRESSOR_FOR_CONSTRUCT.

    The splits module cross-checks construct against stressor, so a frame invented
    freely is refused. Both DTT constructs and both peroxide ones appear, at five doses
    and two replicates, which is the smallest frame the interpolation split accepts.
    """
    import pandas as pd

    rows = []
    for plate in ("p0", "p1", "p2", "p3"):
        for construct, stressor in (("UPRE1", "DTT"), ("UPRE2", "DTT"),
                                    ("NativeYap1", "H2O2"), ("AlteredYap1", "H2O2")):
            for dose in (0.0, 0.25, 0.5, 1.0, 2.0):
                for replicate in (1, 2):
                    rows.append({"plate": plate, "construct": construct,
                                 "stressor": stressor, "dose_mM": dose,
                                 "replicate": replicate, "value": float(len(rows))})
    return pd.DataFrame(rows)


_SPLIT_KEY = ("plate", "construct", "dose_mM")


def _optimism_problem():
    """Sample a location, fit its mean, score by negative squared error.

    Deliberately not this package's latent fit. The estimator is a statement about
    sample-based optimisation in general, and tying the specimen to the thing under test
    would leave the audit unable to separate the two.
    """
    import numpy as np

    from ystwin.analysis.optimism import OptimismProblem

    return OptimismProblem(
        sample=lambda rng: float(rng.normal()),
        fit=lambda domains, seed, initial: float(np.mean(domains)),
        score=lambda solution, domains, seed: float(
            -np.mean([(solution - d) ** 2 for d in domains])),
    )


@functools.lru_cache(maxsize=None)
def _teacher_inputs():
    """A six-point hourly grid and a two-channel ramp in [0, 1], the teacher's input shape."""
    times = np.linspace(0.0, 1.0, 6)
    return times, np.column_stack([np.linspace(0.0, 0.6, 6), np.linspace(0.0, 0.3, 6)])


@functools.lru_cache(maxsize=None)
def _training_episodes():
    """Two short teacher episodes, on their own fixed seed, for the student fit.

    Built away from the seed under test so the training set cannot move when that seed does,
    which would confound the two.
    """
    from ystwin.generator.in_silico import generate_episodes

    return generate_episodes(2, seed=20260826, hours=1.0, dt=0.2)


@functools.lru_cache(maxsize=None)
def _magnitude_domain():
    """A three-item poset on the unit box, with the proposal source the sampler requires."""
    from ystwin.analysis import partial_orders as po

    source = po.Source("synthetic", "audit_determinism specimen", "constructed, not measured")
    scope = po.Scope("magnitude", "determinism-audit", "specimen quantity", "arbitrary_unit",
                     po.Calibration("determinism-audit", "known"))
    order = po.compile_order(("a", "b", "c"), scope, [po.OrderClaim("a", "b", scope, source)])
    bounds = [po.Bound(item, 0.0, 1.0, scope, source) for item in order.items]
    return po.MagnitudeDomain(order, bounds, gap=0.0, gap_source=source), source


@functools.lru_cache(maxsize=None)
def _engine_ensemble():
    """Four initialised physical states, their shared priors, and one short interval.

    Four and not two because systematic resampling of two particles cannot depend on the
    draw; short because these specimens ask about the seed, not about the trajectory.
    """
    from ystwin.mech.contracts import Control, Genotype, Protocol
    from ystwin.mech.engine import EngineParameters, initialize

    parameters, genotype = EngineParameters.prior(), Genotype.prior()
    medium = {"glucose": 50.0, "nitrogen": 20.0, "oxygen": 0.2, "osmolyte": 250.0,
              "acetate": 2.0}
    states = tuple(initialize(parameters, genotype, volume_l=1.0, biomass_gdw_l=density,
                              medium_mM=medium) for density in (0.08, 0.10, 0.13, 0.17))
    protocol = Protocol((0.0, 0.002), (Control(0.0, oxygen_transfer_per_h=20.0,
                                               oxygen_saturation_mM=0.2),))
    return parameters, genotype, states, protocol


@functools.lru_cache(maxsize=None)
def _carotenoid():
    """The committed pathway spec and steady-state tables the product validation reads.

    The measured tables rather than an invented frame, for the reason the dose-response
    specimen gives: this is the call the scripts actually make.
    """
    from ystwin.pathway.flux import carotenoid_measurements
    from ystwin.pathway.spec import load_pathway

    return load_pathway("beta_carotene"), carotenoid_measurements()


def _run_discrimination_power(discrimination_power, seed):
    """A decoupled grid, because the collinear design can only ever answer zero.

    This specimen used ``stressor_only_series``, where dose and growth rate move
    together -- exactly the design the function exists to show is unusable. Power was
    0.0 at every fold from 1.02 to 1.5, so no seed could change it, and the sensitivity
    check was reporting a property of the design rather than of the function. On
    ``decoupling_grid`` at 1.3x it sits mid-range and moves with the seed.
    """
    from ystwin.generator.design import decoupling_grid
    from ystwin.generator.literature import parameters_from_literature

    params = parameters_from_literature("UPRE2")
    grid = decoupling_grid(params, doses=(0.0, 0.5, 1.0), nutrient_factors=(1.0, 0.5))
    return discrimination_power(params, grid, induction_fold=1.3, n_replicates=3,
                                n_simulations=15, seed=seed)


_SVD_REASON = (
    "fit_latent recovers its states by SVD, which is deterministic, so the seed reaches "
    "only the cross-validated held-out error -- which this function does not use. The "
    "docstring states it, and a spread computed across seeds would read as a precise "
    "estimate while being no estimate at all."
)

_FEASIBILITY_REASON = (
    "reports which strata can support a split at all, which is a property of the frame "
    "rather than of the draw. The seed is accepted so one call can report and then make "
    "the same split."
)

_DETERMINISTIC_FIT_REASON = (
    "the student is a closed-form ridge solve, and the module says so itself: the fitted "
    "model records deterministic_fit=True among its own hyperparameters. The seed is kept "
    "as provenance and never drawn from, so the specimen digests the four learned maps and "
    "their residual covariances rather than the metadata that echoes the seed back -- "
    "digesting the echo would satisfy the sensitivity check with no draw anywhere."
)

_DECLARED_DESIGN_REASON = (
    "a reference design is a declaration, not a draw. The seed is stored and reaches a "
    "generator only at _assay_noise, and every assay in both default designs declares "
    "noise_sd=0.0, so nothing downstream moves under it either. The specimen digests the "
    "declared experiments and assays: a container that echoes its seed back would otherwise "
    "pass the sensitivity check while no seeded computation had happened at all."
)

_FAMILY_REPORTERS = ("HSE-heat", "PACE-proteasome", "CSRE-carbon")

_STATED_TOPOLOGY_REASON = (
    "the baseline family IS the literature topology as stress_panel.py states it, so its "
    "loadings are fixed by the panel rather than drawn. Of the seven builders only "
    "loading_noise perturbs the loadings themselves; the rest change topology "
    "deterministically and the seed reaches only the sampling that follows. build_family "
    "is therefore exercised on loading_noise."
)

_TOP_DOSE_REASON = (
    "the extrapolation split holds out the highest dose, which is a fact about the "
    "ladder and not a choice; there is nothing for a seed to vary."
)


# ---------------------------------------------------------------------------
# specimens
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Specimen:
    """One seeded entry point, with a synthetic input small enough to run twice.

    Args:
        target: ``module::name`` of the function this exercises. Resolved by the audit and
            handed to ``call``, so a specimen always exercises the object discovery found
            rather than a re-export that may drift away from it.
        call: Takes the resolved target and a seed, returns whatever it produced.
        seed_sensitive: Whether a different seed *must* give a different answer. False only
            for a function whose output is a discrete choice; see ``reason``.
        reason: Why insensitivity is expected. Required when ``seed_sensitive`` is False,
            and checked, so the flag cannot be set silently.
    """

    target: str
    call: object
    seed_sensitive: bool = True
    reason: str = ""

    def resolve(self):
        """The function or class named by :attr:`target`."""
        module, _, attribute = self.target.partition("::")
        return getattr(importlib.import_module(module), attribute)


def _run_filter(particle_filter, seed):
    """Drive a filter over five noiseless observations and forecast past the last one."""
    from ystwin.estimator import Observation, TwinPriors
    from ystwin.observation import ReporterOptics

    optics = ReporterOptics(gain=2.0e5, background=250.0, autofluorescence=400.0,
                            inner_filter_coeff=0.0)
    priors = TwinPriors(biomass=(0.05, 0.01), reporter=(4e-3, 1e-3),
                        promoter_activity=(1.0e-3, 5e-4), growth_rate=(0.30, 0.08),
                        k_deg=0.05)
    filt = particle_filter(priors=priors, optics=optics, gdcw_per_od=0.42, od_blank=0.09,
                           n_particles=200, seed=seed)
    posteriors = []
    for step in range(5):
        t = step * 1.5
        biomass = 0.05 * np.exp(0.28 * t)
        posteriors.append(filt.update(Observation(
            time_h=t, od=biomass / 0.42 + 0.09,
            rfu=250.0 + (400.0 + 2.0e5 * 3.6e-3) * biomass)))
    return posteriors, filt.forecast(2.0)


def _run_plate(plate_conditions, seed):
    """Generate a six-timepoint plate, which is where the seed actually lands.

    Six because the reporter fit refuses fewer than five, and the audit has to exercise
    the generator rather than its input validation.
    """
    from ystwin.generator.plate import DEFAULT_PANEL, generate_plate

    plate = generate_plate(DEFAULT_PANEL,
                           plate_conditions(n_timepoints=6, duration_h=2.0, seed=seed))
    return plate.od, plate.rfu


def _run_null(compare_to_null, seed):
    """Correlation between two channels against marginal-preserving surrogates."""
    from ystwin.analysis.nulls import matched_marginals

    data = np.nan_to_num(_activities()[:, :3])

    def statistic(values):
        return float(abs(np.corrcoef(values[:, 0], values[:, 1])[0, 1]))

    return compare_to_null(statistic, data, matched_marginals, n_draws=25, seed=seed)


def _run_ratio_power(ratio_discrimination_power, seed):
    """A ratio design at a fold small enough that the power is not pinned at 1.

    The induction fold is chosen low deliberately. At 1.5x this design detects every
    simulated campaign, the estimate saturates at 1.0, and the seed cannot change a
    quantity already at its ceiling -- which would make the sensitivity check report a
    finding about the specimen rather than about the function.
    """
    from ystwin.generator.unmixing import PANEL

    params, conditions = _design()
    return ratio_discrimination_power(
        params, conditions, induction_fold=1.1, n_replicates=3,
        stress_spec=PANEL["YFP"], reference_spec=PANEL["YFP"],
        n_simulations=15, seed=seed)


def _run_design_transfer(design_transfer, seed):
    from ystwin.generator.stress_panel import transcriptional_reporters

    return design_transfer([("DTT", "H2O2")], n_states=2,
                           reporters=transcriptional_reporters()[:5],
                           doses=(0.5, 1.0), replicates=1,
                           stressors=["DTT", "H2O2", "NaCl"], observed=(0, 1, 2), seed=seed)


def _run_panel_dataset(panel_dataset, seed):
    from ystwin.generator.stress_panel import transcriptional_reporters

    return panel_dataset(reporters=transcriptional_reporters()[:5],
                         stressors=["DTT", "H2O2"], doses=(1.0,), replicates=2, seed=seed)


def _run_sensor_model(fit_sensor_model, seed):
    """Calibrate on every committed plate but the last, the split the benchmark holds out."""
    plates = tuple(sorted(_CROSSTALK_READINGS.plate.unique()))
    return fit_sensor_model(_CROSSTALK_READINGS, training_plates=plates[:-1], seed=seed)


def _run_student(fit_student, seed):
    """The learned maps only; _DETERMINISTIC_FIT_REASON says why not the whole model."""
    student = fit_student(_training_episodes(), seed=seed)
    return (student.encoder, student.startup_encoder, student.dynamics, student.controller,
            student.encoder_residual_covariance, student.startup_residual_covariance)


def _run_simulate_teacher(simulate_teacher, seed):
    """One episode on a fixed input ramp, so only the observation noise carries the seed."""
    times, inputs = _teacher_inputs()
    return simulate_teacher(times, inputs, seed=seed, noise_cv=0.05)


def _run_observe_reporters(observe_reporters, seed):
    """Read four reporters off a fixed latent path and a fixed culture."""
    times, _ = _teacher_inputs()
    latent = np.column_stack([np.linspace(0.05, 0.5, times.size)] * 3)
    return observe_reporters(times, latent, np.full(times.size, 0.3),
                             0.1 * np.exp(0.3 * times), seed=seed, noise_cv=0.05)


def _declared_content(design):
    """What a design constructor actually built, without the seed it merely stores."""
    return design.experiments, design.assays


def _run_reference_design(reference_design, seed):
    """Round-trip a HOG recovery design through its own canonical serialisation."""
    from ystwin.analysis.mechanistic_reference import default_reference_design

    base = default_reference_design(seed=0)
    design = reference_design(base.experiments, base.assays, seed)
    return _declared_content(reference_design.from_dict(design.to_dict()))


def _run_multisource_design(multisource_design, seed):
    """Round-trip the HOG-plus-native design, whose two source kinds serialise separately."""
    from ystwin.analysis.mechanistic_reference import default_multisource_reference_design

    base = default_multisource_reference_design(seed=0)
    design = multisource_design(base.experiments, base.assays, seed)
    return _declared_content(multisource_design.from_dict(design.to_dict()))


def _run_uniform_box(sample_uniform_box, seed):
    """Eight accepted points from the unit box, with a budget the rejection loop clears."""
    domain, source = _magnitude_domain()
    return sample_uniform_box(domain, n=8, seed=seed, budget=200, source=source).points


def _run_sample_batch(sample_batch, seed):
    """The persisted report of a batch, built by the only function that constructs one."""
    from ystwin.analysis.partial_orders import sample_uniform_box

    domain, source = _magnitude_domain()
    batch = sample_uniform_box(domain, n=8, seed=seed, budget=200, source=source)
    if not isinstance(batch, sample_batch):
        raise TypeError("the uniform sampler no longer returns the discovered SampleBatch")
    return batch.report()


def _run_observation_model(observation_model, seed):
    """One short interval read through an instrument that actually draws.

    The prior optics put od_noise_sd, rfu_noise_sd and missing_probability at zero, and a
    calibration with no noise leaves the seed nothing to move.
    """
    from ystwin.mech.engine import simulate
    from ystwin.mech.params import Param

    parameters, genotype, states, protocol = _engine_ensemble()
    prior = observation_model.prior(seed=seed)
    noisy = dataclasses.replace(prior, values={**prior.values, **{
        name: Param.asserted(name, value, prior.values[name].units,
                             "synthetic determinism specimen, not a fitted calibration")
        for name, value in (("od_noise_sd", 0.02), ("rfu_noise_sd", 0.5),
                            ("missing_probability", 0.1))}})
    observed = simulate(protocol, genotype, states[0], parameters,
                        observation=noisy).observations
    return observed.values, observed.missing


def _run_physical_filter(physical_filter, seed):
    """One update whose posterior weights are uneven enough to resample.

    The seed reaches only the systematic-resampling offset, so at the shipped threshold of
    0.5 this ensemble never resamples and the specimen would report the threshold instead.
    """
    from ystwin.mech.contracts import ObservationModel

    parameters, genotype, states, protocol = _engine_ensemble()
    filt = physical_filter(states, parameters=parameters, genotype=genotype,
                           observation_model=ObservationModel.prior(seed=0),
                           measurement_noise={"od": 0.04, "rfu": 0.3},
                           weights=(1.0,) * len(states), seed=seed, resample_threshold=1.0)
    result = filt.step(protocol, {"od": 0.22, "rfu": None})
    return (result.parent_indices, result.weights, result.posterior_ess,
            result.unique_ancestors,
            tuple(dict(particle.expected_readings) for particle in result.particles))


def _run_product_intervals(product_prediction_intervals, seed):
    """Predictive bands for one held-out strain, at a draw count small enough to run twice."""
    from ystwin.pathway.flux import ProductCandidate

    spec, states = _carotenoid()
    held_out = sorted(states.strain.unique())[0]
    return product_prediction_intervals(
        spec, states[states.strain != held_out], states[states.strain == held_out],
        ProductCandidate("CrtE", "rate", "saturating"), draws=6, seed=seed)


def _run_product_validation(score_product_validation, seed):
    """One candidate rather than fifty-eight, and the draws the seed actually reaches.

    At the shipped draws=0 nothing is drawn at all, so a specimen left at the default would
    report the seed as unused when it was simply never asked for.
    """
    spec, states = _carotenoid()
    return score_product_validation(spec, states, genes=("CrtE",), entry_laws=("rate",),
                                    branches=("saturating",), draws=6, seed=seed)


_ARGMIN_REASON = (
    "returns an argmin over integer widths; the seed moves the cross-validated scores, "
    "and a stable winner under a moving score is the intended behaviour"
)

SPECIMENS: dict[str, Specimen] = {
    # The estimator's own accuracy measurement. Small replicate counts on purpose: this
    # audit asks whether the same seed gives the same answer, which twelve draws settle as
    # firmly as two hundred, and `scripts/estimator_accuracy.py` is what runs it for real.
    "dose_response_permutation": Specimen(
        "ystwin.analysis.multiplicity::dose_response_permutation",
        lambda fn, seed: fn(_CROSSTALK_READINGS, "UPRE1", max_dose=1.0,
                            n_permutations=200, seed=seed),
        seed_sensitive=False, reason=(
            "the p-value is at its floor and stays there. UPRE1's plate-averaged Spearman "
            "is +0.946, and not one of 200 within-plate reshuffles reaches it, so the "
            "count of exceedances is zero for every seed and p is 1/201 by arithmetic. "
            "That is a property of the data rather than of the draw -- and the specimen is "
            "kept on the real table rather than moved to a weaker case that would wobble, "
            "because a determinism check should exercise the call the scripts actually "
            "make")),
    "recovery_error": Specimen(
        "ystwin.analysis.estimator_accuracy::recovery_error",
        lambda fn, seed: fn(_PLATE_GEOMETRY, 1.333, n_replicates=12, seed=seed)),
    "fold_recovery": Specimen(
        "ystwin.analysis.estimator_accuracy::fold_recovery",
        lambda fn, seed: fn(_PLATE_GEOMETRY, n_replicates=12, seed=seed)),
    "fold_change": Specimen(
        "ystwin.analysis.uncertainty::fold_change",
        lambda fn, seed: fn(_readings(), "R", 1.0, n_resamples=300, seed=seed)),
    "activity_uncertainty": Specimen(
        "ystwin.analysis.uncertainty::activity_uncertainty",
        lambda fn, seed: fn(*_trace(), n_draws=30, seed=seed)),
    "compare_to_null": Specimen(
        "ystwin.analysis.nulls::compare_to_null", _run_null),
    "fit_latent": Specimen(
        "ystwin.analysis.latent::fit_latent",
        lambda fn, seed: fn(_activities(), n_states=2, seed=seed)),
    "select_dimension": Specimen(
        "ystwin.analysis.latent::select_dimension",
        lambda fn, seed: fn(_activities(), max_states=3, seed=seed),
        seed_sensitive=False, reason=_ARGMIN_REASON),
    "detection_power": Specimen(
        "ystwin.analysis.power::detection_power",
        lambda fn, seed: fn(*_design(), induction_fold=1.5, n_replicates=3,
                            n_simulations=15, seed=seed)),
    "discrimination_power": Specimen(
        "ystwin.analysis.power::discrimination_power", _run_discrimination_power),
    "power_curve": Specimen(
        "ystwin.analysis.power::power_curve",
        lambda fn, seed: fn(*_design(), induction_folds=(1.5,), replicate_counts=(3,),
                            n_simulations=15, seed=seed)),
    "replicates_needed": Specimen(
        "ystwin.analysis.power::replicates_needed",
        lambda fn, seed: fn(*_design(), induction_fold=1.5, max_replicates=4,
                            n_simulations=15, seed=seed),
        seed_sensitive=False,
        reason="returns the first replicate count clearing a power threshold; an integer "
               "search over a noisy curve is allowed to land on the same answer"),
    "ratio_discrimination_power": Specimen(
        "ystwin.analysis.power::ratio_discrimination_power", _run_ratio_power),
    "module_recovery": Specimen(
        "ystwin.analysis.recovery::module_recovery",
        lambda fn, seed: fn(_panel(), n_states=2, seed=seed),
        seed_sensitive=False, reason=_SVD_REASON),
    "basic_bootstrap_upper": Specimen(
        "ystwin.analysis.optimism::basic_bootstrap_upper",
        lambda fn, seed: fn([0.4, 0.1, 0.9, 0.3, 0.6, 0.2, 0.8, 0.5],
                            n_bootstrap=200, seed=seed)),
    "estimate_optimism": Specimen(
        "ystwin.analysis.optimism::estimate_optimism",
        lambda fn, seed: fn(_optimism_problem(), n_candidate_domains=3,
                            n_reference_domains=2, n_reference=4, n_bootstrap=100,
                            seed=seed)),
    "make_split": Specimen(
        "ystwin.analysis.splits::make_split",
        lambda fn, seed: fn(_split_frame(), "interpolation", group_key=_SPLIT_KEY,
                            seed=seed, validation_units=1).groups),
    "make_split[extrapolation]": Specimen(
        "ystwin.analysis.splits::make_split",
        lambda fn, seed: fn(_split_frame(), "extrapolation", group_key=_SPLIT_KEY,
                            seed=seed, validation_units=1).groups,
        seed_sensitive=False, reason=_TOP_DOSE_REASON),
    "Split": Specimen(
        "ystwin.analysis.splits::Split",
        lambda fn, seed: _split_frame().pipe(
            lambda f: __import__("ystwin.analysis.splits", fromlist=["make_split"])
            .make_split(f, "interpolation", group_key=_SPLIT_KEY, seed=seed,
                        validation_units=1)).groups),
    "feasibility": Specimen(
        "ystwin.analysis.splits::feasibility",
        lambda fn, seed: fn(_split_frame(), "interpolation", group_key=_SPLIT_KEY,
                            seed=seed),
        seed_sensitive=False, reason=_FEASIBILITY_REASON),
    "baseline_family": Specimen(
        "ystwin.generator.families::baseline_family",
        lambda fn, seed: fn(seed=seed).loadings(_FAMILY_REPORTERS),
        seed_sensitive=False, reason=_STATED_TOPOLOGY_REASON),
    "build_family": Specimen(
        "ystwin.generator.families::build_family",
        lambda fn, seed: fn("loading_noise", seed=seed).loadings(_FAMILY_REPORTERS)),
    "family_dataset": Specimen(
        "ystwin.generator.families::family_dataset",
        lambda fn, seed: fn(
            __import__("ystwin.generator.families", fromlist=["baseline_family"])
            .baseline_family(seed=0),
            reporters=_FAMILY_REPORTERS, stressors=("DTT", "H2O2"),
            doses=(0.5, 1.0), replicates=2, seed=seed).readings),
    "attribution_power_table": Specimen(
        "ystwin.viz.figures::attribution_power_table",
        lambda fn, seed: fn(n_simulations=12, seed=seed, folds=(1.3,), n_replicates=3)),
    "module_recovery[shuffle]": Specimen(
        "ystwin.analysis.recovery::module_recovery",
        lambda fn, seed: fn(_panel(), n_states=2, seed=seed, shuffle=True)),
    "train_stress_model": Specimen(
        "ystwin.analysis.stress_model::train_stress_model",
        lambda fn, seed: fn(_panel(), n_states=2, seed=seed)),
    "module_transfer": Specimen(
        "ystwin.analysis.stress_model::module_transfer",
        lambda fn, seed: fn(_panel(), held_out="DTT", n_states=2, seed=seed)),
    "select_width_for_modules": Specimen(
        "ystwin.analysis.stress_model::select_width_for_modules",
        lambda fn, seed: fn(_panel(), max_states=2, sample=1, seed=seed),
        seed_sensitive=False, reason=_ARGMIN_REASON),
    "transfer_test": Specimen(
        "ystwin.analysis.transfer::transfer_test",
        lambda fn, seed: fn(_panel(), held_out="DTT", n_states=2, observed=(0, 1, 2),
                            seed=seed)),
    "leave_one_stressor_out": Specimen(
        "ystwin.analysis.transfer::leave_one_stressor_out",
        lambda fn, seed: fn(_panel(), n_states=2, observed=(0, 1, 2), seed=seed)),
    "select_dimension_by_transfer": Specimen(
        "ystwin.analysis.transfer::select_dimension_by_transfer",
        lambda fn, seed: fn(_panel(), max_states=2, observed=(0, 1, 2), seed=seed),
        seed_sensitive=False, reason=_ARGMIN_REASON),
    "design_transfer": Specimen(
        "ystwin.analysis.experiment_design::design_transfer", _run_design_transfer),
    "panel_dataset": Specimen(
        "ystwin.generator.panel_experiment::panel_dataset", _run_panel_dataset),
    "PlateConditions": Specimen(
        "ystwin.generator.plate::PlateConditions", _run_plate),
    "ParticleFilter": Specimen(
        "ystwin.estimator::ParticleFilter", _run_filter),
    "fit_sensor_model": Specimen(
        "ystwin.analysis.hybrid_stress::fit_sensor_model", _run_sensor_model),
    "fit_student": Specimen(
        "ystwin.analysis.in_silico::fit_student", _run_student,
        seed_sensitive=False, reason=_DETERMINISTIC_FIT_REASON),
    "generate_episodes": Specimen(
        "ystwin.generator.in_silico::generate_episodes",
        lambda fn, seed: fn(2, seed=seed, hours=1.0, dt=0.2)),
    "simulate_teacher": Specimen(
        "ystwin.generator.in_silico::simulate_teacher", _run_simulate_teacher),
    "observe_reporters": Specimen(
        "ystwin.generator.in_silico::observe_reporters", _run_observe_reporters),
    "fold_change_coverage": Specimen(
        "ystwin.analysis.uncertainty::fold_change_coverage",
        lambda fn, seed: fn(3, n_trials=4, n_resamples=60, seed=seed)),
    "validate_modules": Specimen(
        "ystwin.analysis.validation::validate_modules",
        lambda fn, seed: fn(_panel(), "DTT", max_states=2, seed=seed)),
    "sample_uniform_box": Specimen(
        "ystwin.analysis.partial_orders::sample_uniform_box", _run_uniform_box),
    "SampleBatch": Specimen(
        "ystwin.analysis.partial_orders::SampleBatch", _run_sample_batch),
    "ReferenceDesign": Specimen(
        "ystwin.analysis.mechanistic_reference::ReferenceDesign", _run_reference_design,
        seed_sensitive=False, reason=_DECLARED_DESIGN_REASON),
    "MultiSourceReferenceDesign": Specimen(
        "ystwin.analysis.mechanistic_reference::MultiSourceReferenceDesign",
        _run_multisource_design, seed_sensitive=False, reason=_DECLARED_DESIGN_REASON),
    "default_reference_design": Specimen(
        "ystwin.analysis.mechanistic_reference::default_reference_design",
        lambda fn, seed: _declared_content(fn(seed=seed)),
        seed_sensitive=False, reason=_DECLARED_DESIGN_REASON),
    "default_multisource_reference_design": Specimen(
        "ystwin.analysis.mechanistic_reference::default_multisource_reference_design",
        lambda fn, seed: _declared_content(fn(seed=seed)),
        seed_sensitive=False, reason=_DECLARED_DESIGN_REASON),
    "ObservationModel": Specimen(
        "ystwin.mech.contracts::ObservationModel", _run_observation_model),
    "PhysicalParticleFilter": Specimen(
        "ystwin.mech.inference::PhysicalParticleFilter", _run_physical_filter),
    "product_prediction_intervals": Specimen(
        "ystwin.pathway.flux::product_prediction_intervals", _run_product_intervals),
    "score_product_validation": Specimen(
        "ystwin.pathway.flux::score_product_validation", _run_product_validation),
}


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------


def check_specimen(name: str, specimen: Specimen, seed: int = 0,
                   alternatives=ALTERNATIVE_SEEDS) -> list[dict]:
    """Run one specimen twice at ``seed``, then once at each alternative.

    Returns two rows: repeatability, and seed sensitivity against the specimen's declared
    expectation. A specimen that raises produces a failure rather than a skip, because a
    seeded entry point that cannot be called with a small synthetic input is itself a
    finding.
    """
    label = specimen.target.split("::")[-1] if name == specimen.target else name
    if not specimen.seed_sensitive and not specimen.reason:
        return [_row(f"determinism: {label}", "a declared reason", "none", False,
                     "seed_sensitive=False must state why, or it is a way of hiding a bug")]

    started = time.perf_counter()
    try:
        target = specimen.resolve()
        first = digest(specimen.call(target, seed))
        second = digest(specimen.call(target, seed))
    except Exception as exc:
        return [_row(f"determinism: {label}", "runs on a synthetic input",
                     f"{type(exc).__name__}", False, str(exc)[:200])]
    elapsed = time.perf_counter() - started

    rows = [_row(
        f"same seed, same answer: {label}", f"digest({seed}) == digest({seed})",
        "identical" if first == second else f"{first} != {second}", first == second,
        f"{elapsed:.2f} s for two runs" if first == second else
        "the same seed produced two different results; this table cannot be regenerated")]

    digests = {}
    for other in alternatives:
        try:
            digests[other] = digest(specimen.call(target, other))
        except Exception as exc:
            rows.append(_row(f"seed changes the answer: {label}", "runs on a new seed",
                             f"{type(exc).__name__}", False, str(exc)[:200]))
            return rows
    differing = sorted(s for s, d in digests.items() if d != first)
    if specimen.seed_sensitive:
        rows.append(_row(
            f"seed changes the answer: {label}",
            f"a different digest for at least one of {list(alternatives)}",
            f"{len(differing)}/{len(alternatives)} differ", bool(differing),
            "" if differing else
            "every seed gives the same answer, so the seed is accepted and not used; any "
            "spread computed across seeds here is zero by construction"))
    else:
        rows.append(_row(
            f"seed changes the answer: {label}",
            f"same answer for every seed -- {specimen.reason}",
            f"{len(differing)}/{len(alternatives)} differ", not differing,
            "" if not differing else
            "declared seed-insensitive and it is not; the declaration is now wrong"))
    return rows


def check_coverage(discovered: dict, specimens: dict[str, Specimen],
                   exempt: dict[str, str]) -> list[dict]:
    """Every discovered seeded entry point must be exercised or exempt by name.

    The nag that keeps the audit honest as the package grows. A new seeded function is
    reported here the day it appears, rather than being quietly outside the audit forever
    -- which is what a hardcoded list of targets would give.
    """
    rows = []
    unreachable = [k for k in discovered if k.endswith(f"::{IMPORT_FAILED}")]
    for name in unreachable:
        rows.append(_row(f"module importable: {name.split('::')[0]}", "imports",
                         "import failed", "SKIP",
                         "an optional dependency this machine does not have"))
    covered = {s.target for s in specimens.values()}
    missing = sorted(set(discovered) - covered - set(exempt) - set(unreachable))
    rows.append(_row(
        "every seeded entry point is exercised",
        f"{len(discovered) - len(unreachable)} discovered",
        f"{len(covered & set(discovered))} covered, {len(exempt)} exempt, {len(missing)} neither",
        not missing,
        "" if not missing else
        f"no synthetic input for: {', '.join(missing)}"))
    stale = sorted(covered - set(discovered))
    rows.append(_row(
        "no specimen targets a function that no longer exists", "none", f"{len(stale)}",
        not stale, ", ".join(stale)))
    return rows


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def collect(seed: int = 0) -> list[dict]:
    """Run every check and return the rows."""
    import ystwin

    discovered = discover_seeded(ystwin)
    rows = check_coverage(discovered, SPECIMENS, EXEMPT)
    for name, specimen in SPECIMENS.items():
        rows += check_specimen(name, specimen, seed=seed)
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true", help="only print failures")
    parser.add_argument("--seed", type=int, default=0, help="the seed to repeat")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        help="write diagnostics here; default verification is read-only")
    args = parser.parse_args(argv)

    rows = collect(seed=args.seed)
    frame = pd.DataFrame(rows)
    report = None
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        report = args.output_dir / "audit_determinism.csv"
        frame.to_csv(report, index=False)

    failures = frame[frame.status == "FAIL"]
    skips = frame[frame.status == "SKIP"]
    width = 72
    print("=" * width)
    print("Determinism audit: same seed twice, and a different seed once")
    print("=" * width)
    for _, r in frame.iterrows():
        if args.quiet and r.status == "PASS":
            continue
        mark = {"PASS": "ok  ", "SKIP": "skip"}.get(r.status, "FAIL")
        print(f"  [{mark}] {r.check}")
        if r.status != "PASS":
            print(f"         expected {r.expected}, got {r.actual}")
            if r.detail:
                print(f"         {r.detail}")

    print()
    destination = f" -> {report}" if report else " (read-only; no report written)"
    print(f"  {len(frame) - len(failures) - len(skips)}/{len(frame) - len(skips)} checks "
          f"pass ({len(skips)} skipped){destination}")
    if len(failures):
        print(f"  {len(failures)} FAILING. A seeded result here is not what it appears to be.")
    return 1 if len(failures) else 0


if __name__ == "__main__":
    raise SystemExit(main())
