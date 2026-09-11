"""Predict the flux INTO a heterologous pathway from the genotype, not from the product.

This module exists to remove a circularity. `Environment.pathway_flux` used to be a
caller-supplied constant, and both entry points obtained it as
``q_lycopene + q_beta_carotene`` from the very chemostat state being predicted -- the model
measured the product in order to predict the product. `scripts/predict_product.py` said so
in its own output ("this is a residual and not a prediction"), which was honest, and no
non-circular source existed anywhere in the tree.

**Why the GEM cannot supply it, settled rather than assumed.** Capping every reaction of
the installed pathway at exactly 1e-3, then at 1e-4 -- the measured flux magnitude -- gives
FVA ranges of [0, 1e-3] and [0, 1e-4]. The ceiling becomes exact and the FLOOR NEVER LEAVES
ZERO. A capacity constraint is an upper bound, an upper bound cannot make a flux mandatory,
and GECKO's ``v <= kcat*[E]`` under ``sum(MW*e) <= sigma*f*P`` (PMID 28779005) has the zero
vector feasible too. No amount of better parameterisation changes this. The GEM's honest job
is to audit a flux, not to supply one.

**What does supply it: an empirical expression law.** Relative expression times a fitted
scalar is an exploratory candidate, not an established causal cassette-dosage law. The
per-gene leave-one-STRAIN-out comparison against training-only geometric means gives:

    CrtE    rmse(log) 0.200   skill +0.745
    CrtYB   rmse(log) 0.229   skill +0.708
    CrtI    rmse(log) 0.503   skill +0.359
    ERG9    rmse(log) 0.782   skill +0.003

The other native genes score between -0.003 and -0.268. Ranking fourteen genes on these
same three strains does not validate the winner or rule out expression confounding.
:func:`score_product_validation` therefore selects BOTH the entry gene and the model
family inside each outer training fold, scoring beta-carotene and lycopene together.
There are three independent strain groups and six condition predictions, not six
independent samples. The fixed CrtE law remains an explicitly selected-on-this-data
comparison, not a prespecified model.

**The shipped fixed rate law has no growth-rate term.** Across three strains the flux
ratio between dilution rates is 1.09, 0.60 and 0.89. A rate-proportional entry law and a
content-proportional entry law are distinct hypotheses; nested validation compares them
without confusing either with the exact conversion content = rate / measured mu. The
saturating branch uses the existing solver and a growth window learned only from the
training fold. Its kinetic parameters and the entry scale are never borrowed from the
full-data calibration for validation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import expit
from scipy.stats import norm

from .solve import NodeKinetics, solve_pathway
from .spec import Fate, PathwaySpec, RateLaw

__all__ = [
    "CAROTENOID_GENES",
    "FluxCalibration",
    "PredictedFlux",
    "ProductCandidate",
    "FittedProductModel",
    "carotenoid_measurements",
    "fit_flux_law",
    "fit_product_candidate",
    "fit_saturating_branch",
    "predict_flux_from_expression",
    "predict_product_candidate",
    "product_candidates",
    "product_prediction_intervals",
    "score_product_validation",
    "summarize_product_validation",
]


@dataclass(frozen=True)
class FluxCalibration:
    """The one fitted scalar relating entry-enzyme expression to pathway flux.

    Args:
        alpha: Flux per unit relative expression, mmol/gDCW/h. The whole law.
        entry_enzyme: Which gene's expression the scalar is relative to. Carried so a
            calibration cannot be silently applied to a different gene's numbers.
        loso_rmse_log: Leave-one-strain-out log RMSE conditional on the named gene,
            not an error estimate for selecting that gene on the same outcomes.
        loso_skill: RMSE skill against the geometric mean of the training strains.
        n_states: How many measured states it was fitted on.
        source: Citation.
        expression_range: The ``(low, high)`` relative expressions actually seen in the
            fit. The extrapolation warning is raised against THIS rather than against a
            literal window -- an earlier version hard-coded ``0.2 <= e <= 5.0`` and
            described it in the note as "roughly 0.5-3x", while the measured CrtE range is
            [0.245, 1.000], a 4.08-fold span whose maximum is 5x below the literal ceiling.
            A guard whose bounds are unrelated to the data cannot fire where it matters.
    """

    alpha: float
    entry_enzyme: str
    loso_rmse_log: float
    loso_skill: float
    n_states: int
    source: str
    expression_range: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        if self.alpha <= 0:
            raise ValueError(f"alpha must be positive, got {self.alpha}")
        if self.n_states < 2:
            raise ValueError("a scalar cannot be fitted on fewer than two states")

    @property
    def typical_fold_error(self) -> float:
        """The LOSO error as a fold, which is the unit a reader thinks in."""
        return math.exp(self.loso_rmse_log)


@dataclass(frozen=True)
class PredictedFlux:
    """Pathway flux and where it came from."""

    flux_mmol_per_gdcw_h: float
    expression: float
    calibration: FluxCalibration
    notes: tuple[str, ...] = ()

    def interval(self) -> tuple[float, float]:
        """A one-sigma-in-log band from the leave-one-strain-out error.

        Not a confidence interval on ``alpha`` or a calibrated predictive probability.
        It is a fixed-gene held-out spread and omits gene selection, input measurement
        error and kinetic uncertainty. Product validation propagates those measurement
        and kinetic components separately via conditional two-stage strain resampling.
        """
        spread = math.exp(self.calibration.loso_rmse_log)
        return self.flux_mmol_per_gdcw_h / spread, self.flux_mmol_per_gdcw_h * spread


def predict_flux_from_expression(
    expression: float,
    calibration: FluxCalibration,
) -> PredictedFlux:
    """Pathway flux from one relative-expression number.

    Args:
        expression: Relative expression of ``calibration.entry_enzyme``, on the same
            normalisation the calibration was fitted on.
        calibration: Fitted scalar.

    Raises:
        ValueError: on a non-positive expression. Zero expression is not a small flux --
            it is a strain that does not carry the cassette, and the law is multiplicative
            in log space, so it has nothing to say about that case.
    """
    if not math.isfinite(expression) or expression <= 0:
        raise ValueError(
            f"relative expression must be positive and finite, got {expression}. Zero is not a small "
            "flux, it is a strain without the pathway, and this law is fitted in log space")
    notes: list[str] = []
    flux = calibration.alpha * expression
    window = calibration.expression_range
    if window is not None and not window[0] <= expression <= window[1]:
        fold = (expression / window[1]) if expression > window[1] else (window[0] / expression)
        notes.append(
            f"relative expression {expression:g} is outside the range the law was fitted "
            f"over, [{window[0]:.3g}, {window[1]:.3g}] -- {fold:.1f}x beyond it. The law is "
            "a proportionality and nothing checks that it stays one out here")
    elif window is None:
        notes.append(
            "this calibration declares no expression range, so nothing can say whether "
            "the input is inside what it was fitted on")
    return PredictedFlux(flux, float(expression), calibration, tuple(notes))


def fit_flux_law(
    expression: list[float],
    flux: list[float],
    strain: list[str],
    entry_enzyme: str,
    source: str,
) -> FluxCalibration:
    """Fit ``flux = alpha * expression`` in log space, scored leave-one-strain-out.

    Leave-one-STRAIN-out rather than leave-one-state-out, because the question this law has
    to answer is "what does a strain I have never measured produce". Holding out a state
    while keeping its sibling leaves the strain's own flux level in the training set, which
    is most of the answer. This function conditions on the named gene and does not
    validate choosing it from a gene sweep; that selection needs an inner validation.

    Args:
        expression: Relative expression of ``entry_enzyme``, one per state.
        flux: Measured total pathway flux, mmol/gDCW/h, one per state.
        strain: Strain label per state; states sharing a label are held out together.
        entry_enzyme: The gene ``expression`` refers to.
        source: Citation for the data.

    Raises:
        ValueError: on mismatched lengths, non-positive values, or fewer than two strains
            (with one strain there is nothing to hold out and the score would be in-sample).
    """
    if not len(expression) == len(flux) == len(strain):
        raise ValueError("expression, flux and strain must be the same length")
    if len(expression) < 2:
        raise ValueError("need at least two states")
    if any(not math.isfinite(e) or e <= 0 for e in expression) or any(
            not math.isfinite(f) or f <= 0 for f in flux):
        raise ValueError("expression and flux must be positive and finite; the fit is in log space")
    strains = sorted(set(strain))
    if len(strains) < 2:
        raise ValueError(
            f"need at least two strains to score leave-one-strain-out, got {strains}. "
            "A single-strain fit can only be scored in sample, which is not a score")

    log_e = [math.log(e) for e in expression]
    log_f = [math.log(f) for f in flux]
    alpha = math.exp(sum(f - e for e, f in zip(log_e, log_f)) / len(log_e))

    residuals: list[float] = []
    baseline_residuals: list[float] = []
    for held_out in strains:
        train = [i for i, s in enumerate(strain) if s != held_out]
        test = [i for i, s in enumerate(strain) if s == held_out]
        fitted = sum(log_f[i] - log_e[i] for i in train) / len(train)
        residuals += [fitted + log_e[i] - log_f[i] for i in test]
        training_mean = sum(log_f[i] for i in train) / len(train)
        baseline_residuals += [training_mean - log_f[i] for i in test]
    rmse = math.sqrt(sum(r * r for r in residuals) / len(residuals))

    baseline = math.sqrt(sum(r * r for r in baseline_residuals) / len(baseline_residuals))
    skill = 1.0 - rmse / baseline if baseline > 0 else float("nan")

    return FluxCalibration(alpha, entry_enzyme, rmse, skill, len(expression), source,
                           expression_range=(min(expression), max(expression)))


CAROTENOID_GENES = (
    "CrtE", "CrtYB", "CrtI", "ERG10", "ERG12", "ERG13", "ERG20", "ERG8", "ERG9",
    "HMG1", "HMG2", "BTS1", "IDI1", "MVD1",
)
"""Explicit exploratory candidate set: all fourteen released transcript measurements.

The gene list is not inferred from outcome columns or filtered using held-out values.
Its order also breaks exact inner-score ties after parameter count.
"""

_PRODUCT_CHANNELS = ("q_betacarotene", "q_lycopene")


@dataclass(frozen=True)
class ProductCandidate:
    gene: str | None
    entry_law: str
    branch: str

    def __post_init__(self) -> None:
        if self.gene is None:
            if self.entry_law not in {"constant_rate", "constant_content"} or self.branch != "none":
                raise ValueError("a no-expression baseline must be constant_rate or constant_content")
        elif not self.gene or self.entry_law not in {"rate", "content"} or self.branch not in {
                "partition", "saturating"}:
            raise ValueError("entry candidates require a gene, rate/content, and partition/saturating")

    @property
    def name(self) -> str:
        return (self.entry_law if self.gene is None
                else f"{self.gene}:{self.entry_law}:{self.branch}")

    @property
    def n_parameters(self) -> int:
        return 3 if self.branch == "saturating" else 2


@dataclass(frozen=True)
class FittedProductModel:
    candidate: ProductCandidate
    alpha: float | None = None
    constant_channels: tuple[float, float] | None = None
    partition: float | None = None
    kinetics: dict[str, NodeKinetics] = field(default_factory=dict)
    expression_range: tuple[float, float] | None = None


def product_candidates(
    genes: tuple[str, ...] = CAROTENOID_GENES,
    entry_laws: tuple[str, ...] = ("rate", "content"),
    branches: tuple[str, ...] = ("partition", "saturating"),
) -> tuple[ProductCandidate, ...]:
    """Two log-loss-optimal training baselines plus the declared gene/model cross product."""
    if len(set(genes)) != len(genes) or any(not isinstance(g, str) or not g for g in genes):
        raise ValueError("candidate genes must be distinct nonempty names")
    candidates = (
        ProductCandidate(None, "constant_content", "none"),
        ProductCandidate(None, "constant_rate", "none"),
        *(ProductCandidate(g, e, b) for g, e, b in product(genes, entry_laws, branches)),
    )
    if len(set(candidates)) != len(candidates):
        raise ValueError("candidate model families must be distinct")
    return candidates


def carotenoid_measurements(directory: Path | None = None) -> pd.DataFrame:
    """Join released rates and expression intervals without averaging or losing conditions.

    Rates and expression share the authoritative (condition, strain, dilution rate) key.
    Missing expression stays missing so an unsupported candidate remains in the evaluation.
    """
    if directory is None:
        from .. import paths

        directory = paths.data_dir() / "carotenoid"
    states = pd.read_csv(directory / "elizondo2025_steady_states.tsv", sep="\t")
    mrna = pd.read_csv(directory / "elizondo2025_relative_mrna.tsv", sep="\t")
    keys = ["condition", "strain", "mu_per_h"]
    wide = mrna.pivot(index=keys, columns="gene",
                      values=["rel_expression", "rel_lo95", "rel_hi95"])
    suffix = {"rel_expression": "", "rel_lo95": "_lo95", "rel_hi95": "_hi95"}
    wide.columns = [f"{gene}{suffix[value]}" for value, gene in wide.columns]
    frame = states.merge(wide.reset_index(), on=keys, how="left", validate="one_to_one")
    frame["flux"] = frame.q_betacarotene + frame.q_lycopene
    frame["rel_expression"] = frame["CrtE"]
    return frame


def _positive_column(frame: pd.DataFrame, column: str) -> np.ndarray:
    if column not in frame:
        raise ValueError(f"missing required column {column!r}")
    values = frame[column].to_numpy(dtype=float)
    if not len(values) or np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError(f"{column} must be positive and finite in every condition")
    return values


def _strain_weights(frame: pd.DataFrame) -> np.ndarray:
    if "strain" not in frame or frame.strain.isna().any():
        raise ValueError("every condition requires a strain group")
    counts = frame.groupby("strain", sort=False).strain.transform("size").to_numpy(dtype=float)
    return 1.0 / (counts * frame.strain.nunique())


def _training_arrays(train: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu = _positive_column(train, "mu_per_h")
    observed = np.column_stack([_positive_column(train, c) for c in _PRODUCT_CHANNELS])
    return mu, observed, _strain_weights(train)


def _require_carotenoid_spec(spec: PathwaySpec) -> None:
    if (spec.product != "beta_carotene" or spec.nodes[-1].name != "beta_carotene"
            or any(n.fate != Fate.DILUTED for n in spec.nodes)
            or [n.name for n in spec.nodes if n.rate_law == RateLaw.SATURATING] != ["lycopene"]
            or any(n.rate_law != RateLaw.PASSTHROUGH
                   for n in spec.nodes if n.name not in {"lycopene", "beta_carotene"})):
        raise ValueError("this validation requires the diluted, lumped carotenoid branch")


def fit_saturating_branch(spec: PathwaySpec, train: pd.DataFrame) -> dict[str, NodeKinetics]:
    """Refit the existing solver on both training pools; no shipped parameters or test data.

    The log-parameter starts, weights, observed total flux, and growth window all come
    from this training subset. Inner folds may contain one strain at two conditions;
    that permits a conditional fit, not an independent validation of its parameters.
    """
    _require_carotenoid_spec(spec)
    mu, observed, weights = _training_arrays(train)
    total = observed.sum(axis=1)
    window = (float(mu.min()), float(mu.max()))

    def residual(theta):
        capacity, km = np.exp(theta)
        kinetics = {"lycopene": NodeKinetics(vmax_per_growth=capacity, km=km,
                                              growth_rate_range=window)}
        predicted = []
        for flux, growth in zip(total, mu):
            solution = solve_pathway(spec, float(flux), float(growth), kinetics)
            predicted.append([solution.terminal.content_mmol_per_gdcw * growth,
                              solution.node("lycopene").content_mmol_per_gdcw * growth])
        values = np.asarray(predicted)
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("the branch fit produced a nonpositive or nonfinite pool")
        return ((np.log(values) - np.log(observed)) * np.sqrt(weights[:, None])).ravel()

    start = np.log([2 * np.exp(weights @ np.log(observed[:, 0] / mu)),
                    np.exp(weights @ np.log(observed[:, 1] / mu))])
    bounds = (np.log([1e-8, 1e-8]), np.log([1.0, 1.0]))
    fits = []
    for offset in (0.0, -1.0, 1.0):
        guess = np.clip(start + offset, bounds[0] + 1e-6, bounds[1] - 1e-6)
        fitted = least_squares(residual, guess, bounds=bounds, max_nfev=1000,
                               ftol=1e-10, xtol=1e-10, gtol=1e-10)
        if fitted.success and np.all(np.isfinite(fitted.fun)):
            fits.append(fitted)
    if not fits:
        raise RuntimeError("the saturating branch did not converge on the training fold")
    fitted = min(fits, key=lambda f: float(f.fun @ f.fun))
    capacity, km = np.exp(fitted.x)
    return {"lycopene": NodeKinetics(vmax_per_growth=float(capacity), km=float(km),
                                     growth_rate_range=window)}


def _fit_partition(train: pd.DataFrame) -> float:
    _, observed, weights = _training_arrays(train)
    fractions = observed / observed.sum(axis=1, keepdims=True)

    def residual(theta):
        log_fractions = np.array([-np.logaddexp(0.0, -theta[0]),
                                  -np.logaddexp(0.0, theta[0])])
        return ((log_fractions - np.log(fractions)) * np.sqrt(weights[:, None])).ravel()

    start = float(weights @ np.log(observed[:, 0] / observed[:, 1]))
    fitted = least_squares(residual, [np.clip(start, -19, 19)], bounds=(-20, 20),
                           ftol=1e-10, xtol=1e-10, gtol=1e-10)
    if not fitted.success:
        raise RuntimeError("the training partition did not converge")
    return float(expit(fitted.x[0]))


def fit_product_candidate(
    spec: PathwaySpec, train: pd.DataFrame, candidate: ProductCandidate,
    *, _branch_cache: dict | None = None,
) -> FittedProductModel:
    """Fit all scales on training strains only, without pretending a one-strain fit has CV.

    Rate entry fits total rate = alpha * expression; content entry fits total content
    = alpha * expression and multiplies by the supplied mu exactly once at prediction.
    The partition and saturating models see the same measured training total flux.
    """
    _require_carotenoid_spec(spec)
    mu, observed, weights = _training_arrays(train)
    if candidate.gene is None:
        target = observed / mu[:, None] if candidate.entry_law == "constant_content" else observed
        constants = np.exp(weights @ np.log(target))
        return FittedProductModel(candidate, constant_channels=tuple(map(float, constants)))
    expression = _positive_column(train, candidate.gene)
    target = observed.sum(axis=1)
    if candidate.entry_law == "content":
        target = target / mu
    alpha = float(np.exp(weights @ (np.log(target) - np.log(expression))))
    cache = {} if _branch_cache is None else _branch_cache
    if candidate.branch not in cache:
        try:
            cache[candidate.branch] = (fit_saturating_branch(spec, train)
                                       if candidate.branch == "saturating" else _fit_partition(train))
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
            cache[candidate.branch] = exc
    branch = cache[candidate.branch]
    if isinstance(branch, Exception):
        raise branch
    return FittedProductModel(
        candidate, alpha=alpha,
        partition=branch if candidate.branch == "partition" else None,
        kinetics=branch if candidate.branch == "saturating" else {},
        expression_range=(float(expression.min()), float(expression.max())))


def predict_product_candidate(
    spec: PathwaySpec, model: FittedProductModel, inputs: pd.DataFrame,
) -> np.ndarray:
    """Return (beta-carotene, lycopene) rates using only expression and measured mu.

    This is a conditional model comparison through the existing mass-balance solver,
    not a call to an unsupported public physiology or cassette prior. Outcomes and
    their intervals are neither required nor read here.
    """
    mu = _positive_column(inputs, "mu_per_h")
    candidate = model.candidate
    if candidate.gene is None:
        predicted = np.tile(model.constant_channels, (len(inputs), 1))
        if candidate.entry_law == "constant_content":
            predicted *= mu[:, None]
    else:
        flux = model.alpha * _positive_column(inputs, candidate.gene)
        if candidate.entry_law == "content":
            flux *= mu
        if candidate.branch == "partition":
            predicted = flux[:, None] * np.array([model.partition, 1 - model.partition])
        else:
            rows = []
            for f, growth in zip(flux, mu):
                solution = solve_pathway(spec, float(f), float(growth), model.kinetics)
                rows.append([solution.terminal.content_mmol_per_gdcw * growth,
                             solution.node("lycopene").content_mmol_per_gdcw * growth])
            predicted = np.asarray(rows)
    if np.any(~np.isfinite(predicted)) or np.any(predicted <= 0):
        raise ValueError("every predicted channel must be positive and finite")
    return predicted


def _input_frame(frame: pd.DataFrame, candidate: ProductCandidate | None) -> pd.DataFrame:
    columns = ["condition", "strain", "mu_per_h"]
    if candidate is not None and candidate.gene is not None:
        columns += [candidate.gene, f"{candidate.gene}_lo95", f"{candidate.gene}_hi95"]
    return frame.reindex(columns=columns).copy()


def _log_errors(predicted: np.ndarray, observed: np.ndarray) -> np.ndarray:
    errors = np.full_like(observed, np.inf, dtype=float)
    valid = (np.isfinite(predicted) & (predicted > 0)
             & np.isfinite(observed) & (observed > 0))
    errors[valid] = np.log(predicted[valid]) - np.log(observed[valid])
    return errors


def _inner_scores(spec, train, candidates):
    records = {candidate: [] for candidate in candidates}
    for inner_strain in sorted(train.strain.unique()):
        fitting = train[train.strain != inner_strain]
        validation = train[train.strain == inner_strain]
        cache = {}
        for candidate in candidates:
            try:
                model = fit_product_candidate(spec, fitting, candidate, _branch_cache=cache)
                failure = ""
            except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                model, failure = None, str(exc)
            for index in range(len(validation)):
                test = validation.iloc[[index]]
                try:
                    if model is None:
                        raise ValueError(failure)
                    predicted = predict_product_candidate(spec, model, _input_frame(test, candidate))
                    reason = ""
                except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                    predicted, reason = np.full((1, 2), np.nan), str(exc)
                observed = test.reindex(columns=_PRODUCT_CHANNELS).to_numpy(dtype=float)
                errors = _log_errors(predicted, observed)
                loss = float(np.mean(errors ** 2))
                records[candidate].append({"strain": inner_strain, "loss": loss,
                                           "failure": reason or ("invalid observation" if not np.isfinite(loss) else "")})
    rows = []
    for order, candidate in enumerate(candidates):
        frame = pd.DataFrame(records[candidate])
        losses = [float(g.loss.to_numpy().mean()) for _, g in frame.groupby("strain")]
        loss = float(np.mean(losses))
        rows.append({"candidate": candidate.name, "gene": candidate.gene,
                     "entry_law": candidate.entry_law, "branch": candidate.branch,
                     "n_parameters": candidate.n_parameters, "candidate_order": order,
                     "inner_joint_log_mse": loss, "inner_rmse_log": math.sqrt(loss),
                     "inner_strains": len(losses), "inner_predictions": len(frame),
                     "inner_failed_predictions": int(frame.failure.ne("").sum()),
                     "failure_reason": "; ".join(dict.fromkeys(f for f in frame.failure if f))})
    return rows


def _choose_candidate(scores, candidates):
    by_name = {candidate.name: candidate for candidate in candidates}
    eligible = [r for r in scores if r["candidate"] in by_name
                and math.isfinite(r["inner_joint_log_mse"])]
    if not eligible:
        return None, float("inf")
    winner = min(eligible, key=lambda r: (r["inner_joint_log_mse"], r["n_parameters"],
                                         r["candidate_order"]))
    return by_name[winner["candidate"]], winner["inner_rmse_log"]


INTERVAL_ASSUMPTIONS = (
    "Conditional on the outer-selected gene/model and measured mu; resample training strains "
    "with their paired conditions, then perturb released expression and product estimates "
    "using their asymmetric 95% intervals in log space and refit entry plus branch together. "
    "Released intervals describe estimate uncertainty, not raw replicate variation. "
    "Future expression uses its input interval; future product estimate precision is sampled "
    "from training intervals at the same mu (all training rows if absent). Measurement "
    "errors are independent conditional on strain; unavailable cross-channel, shared "
    "normalisation and assay covariance are not identified. The reference expression is "
    "conditioned fixed. No additional model-discrepancy or model-selection uncertainty is "
    "claimed. Three strains cannot establish frequentist coverage."
)


def _interval_scales(frame: pd.DataFrame, column: str) -> tuple[np.ndarray, np.ndarray]:
    centre = _positive_column(frame, column)
    low = _positive_column(frame, f"{column}_lo95")
    high = _positive_column(frame, f"{column}_hi95")
    if np.any(low > centre) or np.any(high < centre):
        raise ValueError(f"{column} source intervals must bracket the released estimates")
    z = norm.ppf(0.975)
    return (np.log(centre) - np.log(low)) / z, (np.log(high) - np.log(centre)) / z


def _perturb_intervals(frame: pd.DataFrame, columns, rng) -> pd.DataFrame:
    perturbed = frame.copy()
    for column in columns:
        low_sd, high_sd = _interval_scales(frame, column)
        z = rng.normal(size=len(frame))
        perturbed[column] = np.exp(np.log(frame[column].to_numpy(dtype=float))
                                   + z * np.where(z < 0, low_sd, high_sd))
    return perturbed


def _check_interval_options(draws: int, nominal_coverage: float) -> None:
    if isinstance(draws, bool) or not isinstance(draws, (int, np.integer)) or draws < 0:
        raise ValueError("draws must be a nonnegative integer")
    if not math.isfinite(nominal_coverage) or not 0 < nominal_coverage < 1:
        raise ValueError("nominal coverage must be between zero and one")


def product_prediction_intervals(
    spec: PathwaySpec, train: pd.DataFrame, inputs: pd.DataFrame,
    candidate: ProductCandidate, *, draws: int = 500, seed: int = 0,
    nominal_coverage: float = 0.95, resample_strains: bool = True,
    expression_error: bool = True, observation_error: bool = True,
) -> pd.DataFrame:
    """Two-stage, paired-strain conditional predictive bands for both rate channels.

    Stage one samples training strains, retaining all their conditions. Stage two draws
    uncertain expression and rate estimates from the released intervals, refits alpha
    and branch parameters on the SAME draw, and propagates future expression and assay
    uncertainty. This preserves fitted entry/kinetic covariance rather than combining
    independent marginal kinetic intervals. No held-out outcome or outcome interval is
    read. Mu and the outer-selected model are conditioned fixed.

    Marginal bands target ``nominal_coverage`` separately. A second Bonferroni rectangle
    targets that nominal coverage for BOTH channels and ALL conditions of the held-out
    strain. Quantile calibration and the unreported assay covariance are assumptions,
    not results. A failed draw invalidates its condition's band rather than silently
    conditioning on successful draws. Zero draws explicitly requests no intervals.
    """
    _check_interval_options(draws, nominal_coverage)
    if not len(inputs):
        raise ValueError("at least one held-out condition is required")
    rng = np.random.default_rng(seed)
    strains = sorted(train.strain.unique())
    samples = np.full((draws, len(inputs), 2), np.nan)
    failure_reasons = set()
    train_columns = list(_PRODUCT_CHANNELS) if observation_error else []
    if expression_error and candidate.gene is not None:
        train_columns.append(candidate.gene)
    try:
        for column in train_columns:
            _interval_scales(train, column)
        if expression_error and candidate.gene is not None:
            _interval_scales(inputs, candidate.gene)
        observation_scales = ({c: _interval_scales(train, c) for c in _PRODUCT_CHANNELS}
                              if observation_error else {})
        if not strains:
            raise ValueError("intervals require training strain groups")
    except (ValueError, KeyError) as exc:
        failure_reasons.add(str(exc))
        observation_scales = {}
    if not failure_reasons:
        for draw in range(draws):
            chosen = rng.choice(strains, size=len(strains), replace=True) if resample_strains else strains
            groups = []
            for index, strain in enumerate(chosen):
                group = train[train.strain == strain].copy()
                group["strain"] = f"bootstrap_{index}"
                groups.append(group)
            try:
                sampled_train = _perturb_intervals(pd.concat(groups, ignore_index=True),
                                                    train_columns, rng)
                model = fit_product_candidate(spec, sampled_train, candidate)
            except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                failure_reasons.add(str(exc))
                continue
            for index in range(len(inputs)):
                test = inputs.iloc[[index]]
                try:
                    if expression_error and candidate.gene is not None:
                        test = _perturb_intervals(test, [candidate.gene], rng)
                    predicted = predict_product_candidate(spec, model, test)[0]
                    if observation_error:
                        same_rate = np.flatnonzero(train.mu_per_h.to_numpy() == float(test.mu_per_h.iloc[0]))
                        source = int(rng.choice(same_rate if len(same_rate) else np.arange(len(train))))
                        z = rng.normal(size=2)
                        for channel, column in enumerate(_PRODUCT_CHANNELS):
                            low_sd, high_sd = observation_scales[column]
                            predicted[channel] *= math.exp(z[channel] * (
                                low_sd[source] if z[channel] < 0 else high_sd[source]))
                    if np.any(~np.isfinite(predicted)) or np.any(predicted <= 0):
                        raise ValueError("a bootstrap channel is nonpositive or nonfinite")
                    samples[draw, index] = predicted
                except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                    failure_reasons.add(str(exc))
    rows = []
    tail = (1 - nominal_coverage) / 2
    joint_tail = tail / (2 * len(inputs))
    for index in range(len(inputs)):
        valid = np.all(np.isfinite(samples[:, index]), axis=1)
        row = {"interval_method": "two_stage_strain_bootstrap",
               "interval_status": ("not_requested" if draws == 0 else
                                   "ok" if valid.all() else "failed"),
               "interval_draws": draws, "interval_successful_draws": int(valid.sum()),
               "interval_failed_draws": int(draws - valid.sum()),
               "interval_failure_reason": "; ".join(sorted(failure_reasons)),
               "nominal_channel_coverage": nominal_coverage,
               "nominal_strain_joint_coverage": nominal_coverage,
               "joint_family_size": 2 * len(inputs), "interval_assumptions": INTERVAL_ASSUMPTIONS}
        marginal = (np.quantile(samples[:, index], [tail, 1 - tail], axis=0)
                    if draws and valid.all() else np.full((2, 2), np.nan))
        joint = (np.quantile(samples[:, index], [joint_tail, 1 - joint_tail], axis=0)
                 if draws and valid.all() else np.full((2, 2), np.nan))
        for channel, name in enumerate(("product", "lycopene")):
            row[f"{name}_rate_low"] = float(marginal[0, channel])
            row[f"{name}_rate_high"] = float(marginal[1, channel])
            row[f"joint_{name}_rate_low"] = float(joint[0, channel])
            row[f"joint_{name}_rate_high"] = float(joint[1, channel])
        rows.append(row)
    return pd.DataFrame(rows)


def _prediction_row(spec, train, test, candidate, model, mode, inner_rmse, failure):
    state = test.iloc[0]
    try:
        if model is None:
            raise ValueError(failure or "no candidate has a complete finite inner validation score")
        predicted = predict_product_candidate(spec, model, _input_frame(test, candidate))[0]
        status, reason = "ok", ""
    except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
        predicted, status, reason = np.full(2, np.nan), "failed", str(exc)
    observed = test.reindex(columns=_PRODUCT_CHANNELS).to_numpy(dtype=float)[0]
    errors = _log_errors(predicted, observed)
    if status == "ok" and not np.all(np.isfinite(errors)):
        status, reason = "unscorable", "invalid held-out observation; prediction retained"
    mu = float(state.mu_per_h)
    divisor = mu if math.isfinite(mu) and mu > 0 else float("nan")
    molar_mass = spec.nodes[-1].molar_mass_g_per_mol
    if molar_mass is None:
        raise ValueError("product molar mass must be declared for mg/gDCW reporting")
    branch = model.kinetics.get("lycopene") if model is not None else None
    return {
        "state": state.condition, "strain": state.strain, "mu_per_h": mu,
        "validation_mode": mode, "selected_model": candidate.name if candidate else "none",
        "selected_gene": candidate.gene if candidate else None,
        "selected_entry_law": candidate.entry_law if candidate else None,
        "selected_branch": candidate.branch if candidate else None,
        "inner_rmse_log": inner_rmse, "training_strains": "|".join(sorted(train.strain.unique())),
        "n_training_strains": train.strain.nunique(), "status": status, "failure_reason": reason,
        "entry_expression": state.get(candidate.gene, np.nan) if candidate and candidate.gene else np.nan,
        "alpha": model.alpha if model is not None and model.alpha is not None else np.nan,
        "partition": model.partition if model is not None and model.partition is not None else np.nan,
        "capacity": branch.vmax_per_growth if branch else np.nan,
        "km": branch.km if branch else np.nan,
        "flux_measured": float(observed.sum()), "flux_predicted": float(predicted.sum()),
        "product_measured": observed[0], "product_predicted": predicted[0],
        "lycopene_rate_measured": observed[1], "lycopene_rate_predicted": predicted[1],
        "product_content_measured": observed[0] / divisor,
        "product_content_predicted": predicted[0] / divisor,
        "lycopene_measured": observed[1] / divisor, "lycopene_predicted": predicted[1] / divisor,
        "content_mg_per_gdcw": predicted[0] / divisor * molar_mass,
        "product_log_error": errors[0], "lycopene_log_error": errors[1],
        "joint_log_mse": float(np.mean(errors ** 2)),
    }


def _attach_interval(row, interval, molar_mass):
    row.update(interval)
    mu = row["mu_per_h"]
    divisor = mu if math.isfinite(mu) and mu > 0 else float("nan")
    for name, observed in (("product", row["product_measured"]),
                           ("lycopene", row["lycopene_rate_measured"])):
        low, high = row[f"{name}_rate_low"], row[f"{name}_rate_high"]
        row[f"{name}_content_low"] = low / divisor
        row[f"{name}_content_high"] = high / divisor
        predicted = row["product_predicted" if name == "product" else "lycopene_rate_predicted"]
        row[f"{name}_inside_band"] = bool(
            math.isfinite(predicted) and predicted > 0 and math.isfinite(observed)
            and observed > 0 and low <= observed <= high)
    row["content_low_mg_per_gdcw"] = row["product_content_low"] * molar_mass
    row["content_high_mg_per_gdcw"] = row["product_content_high"] * molar_mass
    row["measured_inside_band"] = row["product_inside_band"]
    row["joint_inside_band"] = bool(
        row["status"] == "ok"
        and row["joint_product_rate_low"] <= row["product_measured"] <= row["joint_product_rate_high"]
        and row["joint_lycopene_rate_low"] <= row["lycopene_rate_measured"] <= row["joint_lycopene_rate_high"])


def score_product_validation(
    spec: PathwaySpec, states: pd.DataFrame, *, genes: tuple[str, ...] = CAROTENOID_GENES,
    entry_laws: tuple[str, ...] = ("rate", "content"),
    branches: tuple[str, ...] = ("partition", "saturating"),
    mode: str = "nested", draws: int = 0, seed: int = 0, nominal_coverage: float = 0.95,
) -> pd.DataFrame:
    """Outer strain validation with inner strain selection of BOTH gene and model family.

    The default selection set is 58 candidates: two per-channel training baselines plus
    fourteen genes x two entry laws x two branch laws. Inner loss averages squared log
    errors equally over both channels, conditions within strain, then strains. Nothing
    is selected using outer outcomes. All parameters and training baselines are refitted
    within each fold. Unsupported predictions count as infinite loss, never missing rows.

    ``fixed_crte_comparison`` is the historically data-selected CrtE/rate/saturating
    model, explicitly NOT a prespecified validation. Companion family comparisons also
    select their gene and entry law in inner folds, rather than using a global ranking.
    Selection scores and comparison predictions are carried in DataFrame.attrs as plain
    records so CLI callers can persist the complete denominators. Predictive intervals
    are conditional on each outer-selected model; see product_prediction_intervals.
    """
    _check_interval_options(draws, nominal_coverage)
    _require_carotenoid_spec(spec)
    if states.empty or states.condition.isna().any() or states.condition.duplicated().any():
        raise ValueError("validation requires distinct, nonmissing condition identifiers")
    if states.strain.isna().any() or states.strain.nunique() < 3:
        raise ValueError("nested validation requires at least three strain groups")
    if mode not in {"nested", "fixed_crte_comparison"}:
        raise ValueError("mode must be nested or fixed_crte_comparison")
    candidates = product_candidates(genes, entry_laws, branches)
    fixed = ProductCandidate("CrtE", "rate", "saturating")
    groups = {
        "nested": candidates,
        "constant_content": (ProductCandidate(None, "constant_content", "none"),),
        "constant_rate": (ProductCandidate(None, "constant_rate", "none"),),
        "entry_partition": tuple(c for c in candidates if c.branch == "partition"),
        "saturating_branch": tuple(c for c in candidates if c.branch == "saturating"),
        "fixed_crte_comparison": (fixed,),
    }
    frame = states.sort_values(["strain", "condition"]).reset_index(drop=True)
    comparisons, selections, primary = [], [], []
    for fold, held_out in enumerate(sorted(frame.strain.unique())):
        train = frame[frame.strain != held_out].copy()
        test = frame[frame.strain == held_out].copy()
        scores = _inner_scores(spec, train, candidates)
        selections.extend({"outer_strain": held_out, **r} for r in scores)
        fitted, cache = {}, {}
        for label, subset in groups.items():
            candidate, inner_rmse = (_choose_candidate(scores, subset) if label != "fixed_crte_comparison"
                                      else (fixed, float("nan")))
            failure = ""
            if candidate is not None and candidate not in fitted:
                try:
                    fitted[candidate] = fit_product_candidate(spec, train, candidate, _branch_cache=cache)
                except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
                    fitted[candidate] = exc
            model = fitted.get(candidate)
            if isinstance(model, Exception):
                failure, model = str(model), None
            rows = [_prediction_row(spec, train, test.iloc[[i]], candidate, model, label,
                                    inner_rmse, failure) for i in range(len(test))]
            if label == mode:
                interval_candidate = candidate or ProductCandidate(None, "constant_rate", "none")
                intervals = product_prediction_intervals(
                    spec, train, _input_frame(test, interval_candidate), interval_candidate,
                    draws=draws if model is not None else 0, seed=seed + fold,
                    nominal_coverage=nominal_coverage)
                for row, interval in zip(rows, intervals.to_dict("records")):
                    if model is None and draws:
                        interval.update(interval_status="failed", interval_draws=draws,
                                        interval_successful_draws=0, interval_failed_draws=draws,
                                        interval_failure_reason=failure or "no supported fitted model")
                    _attach_interval(row, interval, spec.nodes[-1].molar_mass_g_per_mol)
                primary.extend(rows)
            comparisons.extend(rows)
    design = {"n_independent_strains": frame.strain.nunique(), "n_condition_predictions": len(frame),
              "n_channels": 2, "n_selection_candidates": len(candidates),
              "selection_metric": "strain-balanced mean squared log error over both channels"}
    for row in comparisons:
        row.update(design)
    result = pd.DataFrame(primary)
    result.attrs["selection_scores"] = selections
    result.attrs["comparison_predictions"] = comparisons
    result.attrs["design"] = design
    return result


def summarize_product_validation(scored: pd.DataFrame, *, comparisons: bool = True) -> pd.DataFrame:
    """Summarize all attempted conditions, including failures, with strains as the units.

    Skill is 1 - RMSE / baseline RMSE, reported against BOTH training-only constant
    baselines, not the best baseline selected after inspecting outer outcomes.
    Channel and simultaneous empirical coverage are descriptive counts, not evidence
    that three independent strains calibrate a nominal predictive probability.
    """
    frame = (pd.DataFrame(scored.attrs["comparison_predictions"])
             if comparisons and "comparison_predictions" in scored.attrs else scored.copy())
    rows = []
    for mode, group in frame.groupby("validation_mode", sort=False):
        row = {"model": mode, "n_independent_strains": group.strain.nunique(),
               "n_condition_predictions": len(group), "n_channels": 2,
               "n_successful_predictions": int(group.status.eq("ok").sum()),
               "n_failed_predictions": int(group.status.ne("ok").sum())}
        for name in ("product", "lycopene"):
            losses = [np.mean(g[f"{name}_log_error"].to_numpy(dtype=float) ** 2)
                      for _, g in group.groupby("strain")]
            row[f"{name}_rmse_log"] = float(np.sqrt(np.mean(losses)))
            row[f"{name}_failed_predictions"] = int(
                (~np.isfinite(group[f"{name}_log_error"].to_numpy(dtype=float))).sum())
        row["joint_rmse_log"] = float(np.sqrt(np.mean([
            np.mean(g.joint_log_mse.to_numpy(dtype=float)) for _, g in group.groupby("strain")])))
        intervals_attempted = ("interval_draws" in group and group.interval_draws.fillna(0).gt(0).any())
        row["n_intervals_available"] = (int(group.interval_status.eq("ok").sum())
                                         if "interval_status" in group else 0)
        for name in ("product", "lycopene", "joint"):
            covered = (group[f"{name}_inside_band"].eq(True)
                       if f"{name}_inside_band" in group else pd.Series(False, index=group.index))
            row[f"{name}_covered_conditions"] = int(covered.sum())
            row[f"{name}_empirical_coverage"] = float(covered.sum() / len(group)) if intervals_attempted else np.nan
            if name == "joint":
                by_strain = covered.groupby(group.strain).all()
                row["joint_covered_strains"] = int(by_strain.sum())
                row["joint_empirical_strain_coverage"] = (float(by_strain.mean())
                                                          if intervals_attempted else np.nan)
        for column in ("nominal_channel_coverage", "nominal_strain_joint_coverage"):
            row[column] = (float(group[column].iloc[0])
                           if intervals_attempted and column in group else np.nan)
        rows.append(row)
    table = pd.DataFrame(rows)
    for baseline in ("constant_rate", "constant_content"):
        reference = table.loc[table.model == baseline, "joint_rmse_log"]
        denominator = float(reference.iloc[0]) if len(reference) else np.nan
        table[f"rmse_skill_vs_{baseline}"] = (
            1 - table.joint_rmse_log / denominator if math.isfinite(denominator) and denominator > 0
            else np.nan)
    return table
