"""Prior-predictive calibration on the real plates, with a reproducible specification.

NIS measures squared signed innovations in units of predictive variance. Its time-average
has a chi-squared reference band only under independent, approximately Gaussian errors.
Coverage and signed lag-1 autocorrelation must therefore accompany it. Autocorrelation
signals misspecification, but does not by itself identify optics, growth or noise as the
cause. By default both the normalized Gaussian likelihood and predictive moments use
sigma from each particle's expected reading, never the observed outcome. Mixture variance
is Var(prediction) + E(sigma**2). Relative scales and the numerical reference-scale floor
are declared assumptions, not independent noise assays or quantities tuned to NIS.
The explicit observed_scale_legacy mode retains the former outcome-dependent weighting.

Every update records prior, posterior/pre-resampling and post-resampling ESS, a resampling
flag, exact particle diversity and surviving initial ancestors. A collapsed ensemble is
INCONCLUSIVE. The channel table retains descriptive magnitudes even when no well clears
the gate; a printed median is not evidence that the associated variance is trustworthy.

The default uses committed exports selected by the manifest's source_set, recorded blank
positions, and opening-only priors. Initialization readings are not scored by default.
Recorded culture identities, not nonblank readings, define each known plate's population.
Missing cultures and outside-layout wells remain in the ledger, with recorded identity
separate from analysis status; explicit exploratory blank policies do not rewrite identity.
Finite signed OD/RFU readings are retained under additive Gaussian noise; NaN omits only
that channel, while infinity refuses the well with its field and step indices recorded.
Positive opening data needed to initialize physical states are a separate requirement,
not a reason to censor later measurements. The explicit positive_od_legacy reading policy
restores the old full-trace positive-OD population restriction by refusing such wells.
Explicit legacy prior/blank/noise options support comparisons, not a claim to reproduce an
undocumented historical experiment. In particular, docs/research/PREDICTION.md reports
1.76 / 4.44 with 2000 particles and 167 wells, but does not supply an executable run manifest
or exact prior/well/seed specification; attributing the difference to priors is unverified.
The former estimator docstring's deceleration sweep is likewise not reproducible evidence.

--decompose asks which of the three things an inflated NIS can be it actually is. Arms
change exactly one declared quantity from the baseline and hold seed, wells, priors and
readings: a particle sweep moves the particle approximation, and scaling od_rel_sigma and
rfu_rel_sigma together moves the measurement model, since initialization reads neither.
What neither family removes is process-model error by elimination among these three and
only among these three.

One knob is not one effect, and the table says so rather than assuming otherwise. Both
knobs act on the same weight degeneracy -- a looser likelihood keeps particles alive
exactly as a larger ensemble does -- so min_ess_ratio reports how far each arm moved the
ensemble, and a measurement arm away from 1 has not isolated the measurement model. An arm
is still an aggregate over well-channels that individually fail the ESS gate; the
decomposition says where the excess is, never that the variance is sound.

Three evidence tables, an exclusion ledger and a JSON configuration/input manifest are
written only when --no-write is absent, and nis_error_decomposition.csv under
--decompose. --output-dir isolates a run from tracked results;
--config accepts the configuration section of a previous manifest. No timestamp is added.
The manifest describes this computation, not the provenance of a historical experiment.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from ystwin import paths
from ystwin.analysis.estimator_accuracy import (
    NisArm,
    nis_error_decomposition,
    prediction_metrics,
)
from ystwin.estimator import Observation, ParticleFilter, TwinPriors
from ystwin.growth import max_specific_growth_rate
from ystwin.observation import ReporterOptics
from ystwin.plate import replay
from ystwin.plate.gen5 import read_xpt_run
from ystwin.plate.layout import RecordedPlate, recorded_for_export
from ystwin.plate.synergy import read_synergy_kinetic
from ystwin.readings import CorrectedOD

# Pre-registered before looking at the result, so that a well failing the ESS floor is
# excluded for a stated reason rather than because excluding it improved the verdict.
MIN_ESS = 100.0
"""Legacy ESS floor; configurable, and recorded in every run manifest."""

ALPHA = 0.05
"""Two-sided level for the chi-squared reference band."""

GDCW_PER_OD = 0.3
"""Literature dry-weight scale, not a measurement on these plates."""

_PROCESS_FIELDS = (
    "k_deg", "activity_walk", "growth_walk", "biomass_walk", "od_rel_sigma",
    "rfu_rel_sigma", "growth_deceleration", "positive_noise_mode",
    "observation_noise_mode", "observation_scale_floor",
)
_READING_POLICIES = ("finite_gaussian", "positive_od_legacy")
_OUTPUT_NAMES = (
    "nis_innovations.csv", "nis_summary.csv", "nis_channel_summary.csv",
    "nis_exclusions.csv", "nis_manifest.json",
)

DECOMPOSITION_NAME = "nis_error_decomposition.csv"
"""Written only under --decompose, so an ordinary run's output set is unchanged."""


def _acceptance_band(k: int, alpha: float = ALPHA) -> tuple[float, float]:
    """Chi-squared(k)/k reference band for an average of k independent innovations."""
    if not np.isfinite(k) or k < 1 or int(k) != k:
        raise ValueError("k must be a positive integer")
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be in (0, 1)")
    return (
        float(stats.chi2.ppf(alpha / 2.0, df=k) / k),
        float(stats.chi2.ppf(1.0 - alpha / 2.0, df=k) / k),
    )


def _median_or_nan(values: pd.Series) -> float:
    """An unestimable statistic stays missing, including an entirely unestimable group."""
    available = values.dropna()
    return float(available.median()) if len(available) else float("nan")


def _channel_rows(summary: pd.DataFrame, aggregation: str = "well_median",
                  population: str = "all") -> list[dict]:
    """Descriptive aggregates; ESS failures remain counted even in tested-only summaries.

    The compatibility column median_mean_nis always means a median across wells.
    aggregate_mean_nis uses the explicitly selected estimator: equal-well median,
    equal-well mean, or equal-plate median of within-plate medians. These are different
    estimands. Per-well bands and whiteness bounds use each well's own scored K, never
    a pooled K or the median K. Channel medians have no chi-squared reference test.
    """
    if aggregation not in ("well_median", "well_mean", "plate_median"):
        raise ValueError("unknown aggregation")
    if population not in ("all", "tested"):
        raise ValueError("unknown aggregation population")
    rows = []
    if summary.empty:
        return rows
    for channel, block in summary.groupby("channel", sort=True):
        included = block if population == "all" else block[block.tested]
        tested = int(block.tested.sum())
        if aggregation == "well_median":
            aggregate = _median_or_nan(included.mean_nis)
        elif aggregation == "well_mean":
            aggregate = float(included.mean_nis.mean())
        else:
            aggregate = _median_or_nan(included.groupby("plate").mean_nis.agg(_median_or_nan))
        autocorr_estimable = included.lag1_autocorr.notna()
        white = float(included.white.mean()) if len(included) else float("nan")
        verdict = "DESCRIPTIVE_ONLY"
        if tested < len(block) or included.empty or not autocorr_estimable.all():
            verdict = "INCONCLUSIVE"
        elif white < 0.5:
            verdict = "AUTOCORRELATED"
        rows.append({
            "channel": channel, "well_channels": int(len(block)),
            "plates": int(block.plate.nunique()), "k_steps": _median_or_nan(included.k_steps),
            "band_low": _median_or_nan(included.band_low),
            "band_high": _median_or_nan(included.band_high),
            "median_mean_nis": _median_or_nan(included.mean_nis),
            "aggregate_mean_nis": aggregate, "aggregation": aggregation,
            "aggregation_population": population, "aggregated_well_channels": int(len(included)),
            "fraction_in_band": float(included.nis_in_band.mean()),
            "median_lag1_autocorr": _median_or_nan(included.lag1_autocorr),
            "white_bound": _median_or_nan(included.white_bound), "fraction_white": white,
            "autocorrelation_testable_well_channels": int(autocorr_estimable.sum()),
            "median_coverage": _median_or_nan(included.coverage),
            "nominal_coverage": float(block.nominal_coverage.iloc[0]),
            "tested_well_channels": tested, "median_min_ess": _median_or_nan(block.min_ess),
            "median_min_unique_particles": _median_or_nan(block.min_unique_particles),
            "median_min_unique_ancestors": _median_or_nan(block.min_unique_ancestors),
            "resampled_steps": int(block.resampled_steps.sum()), "verdict": verdict,
        })
    return rows


def _lag1_autocorrelation(values: np.ndarray) -> float:
    """Signed lag-1 autocorrelation without joining across missing original steps."""
    x = np.asarray(values, dtype=float)
    return float(prediction_metrics(x, np.zeros_like(x), np.ones_like(x))["lag1_autocorr"])


class BelowBackground(ValueError):
    """Opening data cannot initialize physical states; separate from reading validity."""

    def __init__(self, message: str, *, field: str | None = None):
        super().__init__(message)
        self.field = field


class ReadingRefusal(ValueError):
    """A field violates the explicitly selected reading policy."""

    def __init__(self, message: str, *, field: str):
        super().__init__(message)
        self.field = field


def _validate_readings(od: np.ndarray, rfu: np.ndarray, reading_policy: str) -> None:
    if reading_policy not in _READING_POLICIES:
        raise ValueError(f"unknown reading_policy: {reading_policy!r}")
    for field, values in (("od", od), ("rfu", rfu)):
        infinite = np.flatnonzero(np.isinf(values)).tolist()
        if infinite:
            raise ReadingRefusal(f"{field} contains infinite readings at step indices {infinite}",
                                 field=field)
    if reading_policy == "positive_od_legacy":
        invalid = np.flatnonzero(np.isnan(od) | (od <= 0)).tolist()
        if invalid:
            raise ReadingRefusal(
                f"positive_od_legacy requires finite positive od at every step; invalid indices {invalid}",
                field="od",
            )


def _trace_arrays(times_h, od, rfu, opening_points: int):
    t, o, r = (np.asarray(x, dtype=float) for x in (times_h, od, rfu))
    if t.ndim != 1 or o.shape != t.shape or r.shape != t.shape:
        raise ValueError("time, OD and RFU must share a one-dimensional shape")
    if opening_points < 2 or t.size < opening_points:
        raise ValueError("need at least opening_points >= 2 readings for initialization")
    if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError("times must be finite and strictly increasing")
    return t, o, r


def _priors_for(od: np.ndarray, rfu: np.ndarray, times_h: np.ndarray, *,
                opening_points: int = 3, prior_policy: str = "opening",
                prior_overrides: dict | None = None, optics: ReporterOptics | None = None,
                gdcw_per_od: float = GDCW_PER_OD, carotenoid: float = 0.0) -> TwinPriors:
    """Data-informed initialization, with its support window made explicit.

    The default growth slope uses only opening readings, not the maximum growth rate
    over the future trace. Legacy_full_trace is deliberately named: it leaks future OD
    and cannot be used as evidence for prospective calibration. Custom optics are
    inverted in their own units, and supplied degradation contributes to synthesis.
    """
    t, od, rfu = _trace_arrays(times_h, od, rfu, opening_points)
    opening = slice(0, opening_points)
    if not np.all(np.isfinite(od[opening])) or not np.all(od[opening] > 0):
        raise BelowBackground("opening OD must be finite and positive for biomass and log-growth priors",
                              field="od")
    if not np.all(np.isfinite(rfu[opening])):
        raise BelowBackground("opening reporter must be finite for reporter initialization", field="rfu")
    opening_rfu = float(np.mean(rfu[opening]))
    opening_od = float(np.mean(od[opening]))
    if not opening_rfu > 0:
        raise BelowBackground(
            f"opening reporter mean {opening_rfu:.1f} is not positive for reporter initialization",
            field="rfu",
        )
    if not np.isfinite(gdcw_per_od) or gdcw_per_od <= 0:
        raise ValueError("gdcw_per_od must be finite and positive")
    # The reporter state is a concentration PER gDCW: observe_rfu multiplies it by
    # biomass. Centring the prior on raw RFU instead would under-predict the reading by
    # exactly the biomass factor -- about thirtyfold here -- and the filter would then be
    # scored on a units error rather than on its honesty. This is the per-OD versus
    # per-gDCW boundary that has no converter in the package; crossing it is the caller's
    # job and this is the crossing.
    opening_biomass = opening_od * gdcw_per_od
    opening_per_gdcw = opening_rfu / opening_biomass
    if optics is not None:
        attenuation = float(optics.attenuation(np.asarray(carotenoid)))
        if optics.gain <= 0 or not np.isfinite(attenuation) or attenuation <= 0:
            raise ValueError("optics need positive gain and attenuation to initialize reporter")
        if optics.detector_max is not None and np.any(rfu[opening] >= optics.detector_max):
            raise ValueError("opening reporter is detector-saturated; no invertible reporter prior")
        opening_per_gdcw = ((opening_rfu - optics.background) / (opening_biomass * attenuation)
                            - optics.autofluorescence) / optics.gain
        if not np.isfinite(opening_per_gdcw) or opening_per_gdcw <= 0:
            raise BelowBackground("opening reporter after supplied photophysics is not positive", field="rfu")
    if prior_policy == "opening":
        growth = float(np.polyfit(t[opening] - t[0], np.log(od[opening]), 1)[0])
    elif prior_policy == "legacy_full_trace":
        if not np.all(np.isfinite(od)) or not np.all(od > 0):
            raise BelowBackground("legacy_full_trace growth initialization needs finite positive OD throughout",
                                  field="od")
        growth = float(max_specific_growth_rate(t, CorrectedOD(od)))
    else:
        raise ValueError("unknown prior_policy")
    if not np.isfinite(growth):
        raise ValueError("growth prior is not estimable")
    overrides = dict(prior_overrides or {})
    unknown = overrides.keys() - set(_PROCESS_FIELDS)
    if unknown:
        raise ValueError(f"unknown process/noise parameters: {sorted(unknown)}")
    growth = float(np.clip(growth, 0.05, 0.6))
    synthesis_rate = 0.3 if prior_policy == "legacy_full_trace" else growth
    synthesis_rate += overrides.get("k_deg", 0.0)
    return TwinPriors(
        biomass=(opening_biomass, opening_biomass * 0.5),
        reporter=(opening_per_gdcw, opening_per_gdcw * 0.5),
        promoter_activity=(opening_per_gdcw * synthesis_rate, opening_per_gdcw * synthesis_rate),
        growth_rate=(growth, 0.15), **overrides,
    )


def _run_well(times_h, od, rfu, n_particles: int, seed: int, *,
              prior_overrides: dict | None = None, optics: ReporterOptics | None = None,
              priors: TwinPriors | None = None, opening_points: int = 3,
              prior_policy: str = "opening", include_initialization: bool = False,
              resample_threshold: float = 0.5, gdcw_per_od: float = GDCW_PER_OD,
              carotenoid: float = 0.0, coverage: float = 0.95,
              reading_policy: str = "finite_gaussian") -> pd.DataFrame:
    """One-step prior predictions, including labelled, unscored initialization rows.

    NaN channels are omitted independently and their step indices are recorded. Finite
    signed measurements are retained by default, whereas infinity is an explicit
    reading refusal. Data-informed initialization separately needs a valid opening
    window; supplied physical priors do not impose that requirement on the readings.
    The manifest receives the actual initialized priors via DataFrame.attrs. Supplying
    fixed priors and process overrides together is refused rather than ignoring either.
    """
    times_h, od, rfu = _trace_arrays(times_h, od, rfu, opening_points)
    _validate_readings(od, rfu, reading_policy)
    optics = optics if optics is not None else ReporterOptics(
        gain=1.0, background=0.0,
        # Never measured on any plate -- see docs/DATA_INVENTORY.md. Zero is the
        # assumption the reported analyses already make, carried here so the filter is
        # tested under the same assumption rather than a kinder one.
        autofluorescence=0.0, inner_filter_coeff=None,
    )
    if priors is not None and prior_overrides:
        raise ValueError("supply priors or prior_overrides, not both")
    priors = priors if priors is not None else _priors_for(
        od, rfu, times_h, opening_points=opening_points, prior_policy=prior_policy,
        prior_overrides=prior_overrides, optics=optics, gdcw_per_od=gdcw_per_od,
        carotenoid=carotenoid,
    )
    filt = ParticleFilter(
        priors, optics, gdcw_per_od=gdcw_per_od, od_blank=0.0,
        n_particles=n_particles, seed=seed, resample_threshold=resample_threshold,
    )
    z = float(stats.norm.ppf((1.0 + coverage) / 2.0))
    rows = []
    for index, (t, o, r) in enumerate(zip(times_h, od, rfu)):
        _, innovations = filt.update_with_innovations(Observation(
            time_h=float(t), od=None if np.isnan(o) else float(o),
            rfu=None if np.isnan(r) else float(r), carotenoid=carotenoid,
        ))
        for innovation in innovations:
            spread = np.sqrt(innovation.total_variance)
            rows.append({
                "time_h": innovation.time_h, "step_index": index, "channel": innovation.channel,
                "observed": innovation.observed, "predicted": innovation.predicted,
                "state_variance": innovation.state_variance,
                "measurement_variance": innovation.measurement_variance,
                "ess": innovation.posterior_ess, "prior_ess": innovation.prior_ess,
                "posterior_ess": innovation.posterior_ess,
                "post_resample_ess": innovation.post_resample_ess, "resampled": innovation.resampled,
                "prior_unique_particles": innovation.prior_unique_particles,
                "unique_particles": innovation.unique_particles,
                "unique_ancestors": innovation.unique_ancestors,
                "nis": innovation.nis, "standardised": innovation.standardised,
                "prediction_low": innovation.predicted - z * spread,
                "prediction_high": innovation.predicted + z * spread,
                "covered": bool(abs(innovation.standardised) <= z),
                "scored": bool(include_initialization or index >= opening_points),
            })
    frame = pd.DataFrame(rows)
    frame.attrs["priors"] = asdict(priors)
    frame.attrs["reading_policy"] = reading_policy
    for channel, values in (("od", od), ("rfu", rfu)):
        missing = np.flatnonzero(np.isnan(values)).tolist()
        frame.attrs[f"missing_{channel}_steps"] = len(missing)
        frame.attrs[f"missing_{channel}_step_indices"] = missing
        frame.attrs[f"nonpositive_{channel}_steps"] = int((values <= 0).sum())
    return frame


def _summarise(innovations: pd.DataFrame, *, min_ess: float = MIN_ESS,
               min_unique_particles: int = 2, alpha: float = ALPHA,
               coverage: float = 0.95) -> pd.DataFrame:
    rows = []
    if innovations.empty:
        return pd.DataFrame(rows)
    for (plate, well, channel), group in innovations.groupby(["plate", "well", "channel"]):
        scored = group[group.scored].sort_values("step_index")
        metrics = prediction_metrics(
            scored.observed, scored.predicted, scored.state_variance + scored.measurement_variance,
            nominal_coverage=coverage, step_indices=scored.step_index,
        )
        k = metrics.pop("n_predictions")
        low, high = _acceptance_band(k, alpha) if k else (float("nan"), float("nan"))
        bound = 2.0 / np.sqrt(k) if k else float("nan")
        minimum = float(group.posterior_ess.min())
        minimum_prior = float(group.prior_ess.min())
        unique = int(group.unique_particles.min())
        unique_prior = int(group.prior_unique_particles.min())
        reasons = []
        if k == 0:
            reasons.append("no_scored_predictions")
        if minimum < min_ess or minimum_prior < min_ess:
            reasons.append("ess_below_floor")
        if min(unique, unique_prior) < min_unique_particles:
            reasons.append("particle_diversity_below_floor")
        rows.append({
            "plate": plate, "well": well, "channel": channel, "k_steps": k, **metrics,
            "initialization_steps": int((~group.scored).sum()),
            "min_ess": minimum, "min_prior_ess": minimum_prior,
            "min_post_resample_ess": float(group.post_resample_ess.min()),
            "min_unique_particles": unique, "min_prior_unique_particles": unique_prior,
            "min_unique_ancestors": int(group.unique_ancestors.min()),
            "resampled_steps": int(group.resampled.sum()), "band_low": low, "band_high": high,
            "nis_in_band": bool(low <= metrics["mean_nis"] <= high), "white_bound": bound,
            "white": bool(np.isfinite(metrics["lag1_autocorr"])
                          and abs(metrics["lag1_autocorr"]) <= bound),
            "tested": not reasons, "exclusion_reason": ";".join(reasons),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- error decomposition


def _measurement_variance_share(innovations: pd.DataFrame) -> dict[str, float]:
    """Median measurement share of predictive variance per channel, over scored rows.

    Read off the run rather than assumed, because this share is exactly what bounds
    how much an inflated measurement sigma could ever buy.
    """
    if innovations.empty:
        return {}
    scored = innovations[innovations.scored]
    total = scored.state_variance + scored.measurement_variance
    share = (scored.measurement_variance / total).where(total > 0)
    return {str(channel): float(value)
            for channel, value in share.groupby(scored.channel).median().items()}


def _arm_plan(args) -> list[tuple[str, str, int, float]]:
    """``(label, family, particles, sigma multiplier)``, baseline first, one knob each.

    A requested ensemble size or multiplier equal to the baseline's is the baseline and
    is dropped rather than run twice under a second name.
    """
    plan = [("baseline", "baseline", args.particles, 1.0)]
    for count in dict.fromkeys(args.decompose_particles):
        if count != args.particles:
            plan.append((f"particles={count}", "particles", count, 1.0))
    for multiplier in dict.fromkeys(float(m) for m in args.decompose_sigma):
        if multiplier != 1.0:
            plan.append((f"sigma_x{multiplier:g}", "measurement", args.particles, multiplier))
    return plan


def _arm_args(args, *, particles: int, sigma_multiplier: float):
    """The baseline configuration with one knob moved and nothing else.

    ``_priors_for`` reads only ``k_deg`` out of the process overrides, so scaling the
    relative sigmas moves the measurement model and leaves initialization and the walks
    exactly where the baseline put them.
    """
    return argparse.Namespace(**{
        **vars(args), "particles": particles, "no_write": True, "decompose": False,
        "od_rel_sigma": args.od_rel_sigma * sigma_multiplier,
        "rfu_rel_sigma": args.rfu_rel_sigma * sigma_multiplier,
    })


def _decompose(args, innovations: pd.DataFrame, channels: pd.DataFrame) -> pd.DataFrame:
    """Attribute the baseline's NIS excess among particles, measurement and process.

    Each arm reuses the baseline seed, wells, priors and readings and changes one
    declared quantity, so what moves is attributable to the knob named. Changing the
    ensemble size does change the random stream, so a particle arm is not a coupled
    redraw of the baseline; the comparison is between medians over many well-channels,
    not between wells.
    """
    plan = _arm_plan(args)
    shares = {"baseline": _measurement_variance_share(innovations)}
    aggregates = {"baseline": {row["channel"]: row for row in channels.to_dict("records")}}
    for label, _family, particles, multiplier in plan[1:]:
        arm_innovations, _summary, arm_channels, _ledger, _manifest = _collect(
            _arm_args(args, particles=particles, sigma_multiplier=multiplier))
        if arm_channels.empty:
            raise ValueError(f"decomposition arm {label} scored no well-channels")
        shares[label] = _measurement_variance_share(arm_innovations)
        aggregates[label] = {row["channel"]: row for row in arm_channels.to_dict("records")}
    rows = []
    for channel in sorted(aggregates["baseline"]):
        arms = []
        for label, family, particles, multiplier in plan:
            row = aggregates[label].get(channel)
            if row is None:
                raise ValueError(f"decomposition arm {label} produced no {channel} aggregate")
            arms.append(NisArm(
                label=label, family=family, n_particles=particles, sigma_multiplier=multiplier,
                mean_nis=row["aggregate_mean_nis"],
                measurement_variance_share=shares[label].get(channel, float("nan")),
                min_ess=row["median_min_ess"],
                min_unique_ancestors=row["median_min_unique_ancestors"],
                coverage=row["median_coverage"],
                well_channels=row["aggregated_well_channels"],
            ))
        rows += [{"channel": channel, **entry} for entry in nis_error_decomposition(arms)]
    return pd.DataFrame(rows)


def _decomposition_manifest(args, plan) -> dict:
    return {
        "arms": [{"label": label, "family": family, "particles": particles,
                  "sigma_multiplier": multiplier,
                  "od_rel_sigma": args.od_rel_sigma * multiplier,
                  "rfu_rel_sigma": args.rfu_rel_sigma * multiplier}
                 for label, family, particles, multiplier in plan],
        "families": {
            "particles": "ensemble size only; model, seed, wells, priors and readings held",
            "measurement": "od_rel_sigma and rfu_rel_sigma scaled together; dynamics held, "
                           "including initialization, which reads no relative sigma",
            "process": "not an arm; what neither other family removes, by elimination "
                       "among these three and only among these three",
        },
        "excess_removed": "share of (baseline aggregate NIS - 1) that the arm removes; "
                          "signed and unclipped, so an arm that makes NIS worse is negative",
        "implied_sigma_multiplier": "the multiplier that would reach NIS 1 if the whole "
                                    "excess were an under-declared measurement scale, holding "
                                    "the prediction and state variance; a lower bound on any "
                                    "measurement-only account, and a prediction the measurement "
                                    "arms then confirm or falsify",
        "min_ess_ratio": "the arm's median minimum ESS over the baseline's; a measurement arm "
                         "away from 1 has also moved the particle approximation, because a "
                         "looser likelihood keeps particles alive exactly as more of them do, "
                         "so its excess_removed is not the measurement model's alone",
        "not_established": "the arms bound contributions to the aggregate NIS; they do not "
                           "identify which process term is misspecified, they do not separate "
                           "the two families wherever min_ess_ratio departs from 1, and a "
                           "collapsed ensemble stays INCONCLUSIVE whatever the decomposition says",
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=pathlib.Path, help="JSON configuration or prior run manifest")
    parser.add_argument("--source", choices=("committed", "workbooks", "auto"), default="committed")
    parser.add_argument("--source-set", default="newprotocol")
    parser.add_argument("--plates-dir", type=pathlib.Path)
    for flag in ("plate", "well", "exclude-plate", "exclude-well", "blank-well"):
        parser.add_argument(f"--{flag}", action="append", default=[],
                            help="repeat for exact identifiers; wells may be EXPORT::WELL")
    parser.add_argument("--blank-policy", choices=("recorded", "legacy_low_od"), default="recorded")
    parser.add_argument("--reading-policy", choices=_READING_POLICIES, default="finite_gaussian",
                        help="retain finite signed readings, or explicitly require positive OD throughout; "
                             "NaN is missing and infinity is refused")
    parser.add_argument("--max-wells", type=int, default=None,
                        help="prefix cap for smoke testing only, never a representative sample")
    parser.add_argument("--particles", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--opening-points", type=int, default=3)
    parser.add_argument("--prior-policy", choices=("opening", "legacy_full_trace"), default="opening")
    parser.add_argument("--include-initialization", action="store_true")
    parser.add_argument("--resample-threshold", type=float, default=0.5)
    parser.add_argument("--min-ess", type=float, default=MIN_ESS)
    parser.add_argument("--min-unique-particles", type=int, default=2)
    parser.add_argument("--alpha", type=float, default=ALPHA)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--aggregation", choices=("well_median", "well_mean", "plate_median"),
                        default="well_median")
    parser.add_argument("--aggregation-population", choices=("all", "tested"), default="all")
    parser.add_argument("--gdcw-per-od", type=float, default=GDCW_PER_OD)
    parser.add_argument("--gain", type=float, default=1.0)
    parser.add_argument("--background", type=float, default=0.0)
    parser.add_argument("--autofluorescence", type=float, default=0.0)
    parser.add_argument("--inner-filter-coeff", type=float, default=None)
    parser.add_argument("--detector-max", type=float, default=None)
    parser.add_argument("--carotenoid", type=float, default=0.0)
    for name in _PROCESS_FIELDS:
        default = TwinPriors.__dataclass_fields__[name].default
        if name == "positive_noise_mode":
            parser.add_argument("--positive-noise-mode", default=default,
                                choices=("mean_preserving", "median_preserving_legacy"))
        elif name == "observation_noise_mode":
            parser.add_argument("--observation-noise-mode", default=default,
                                choices=("predicted_scale", "observed_scale_legacy"),
                                help="particle-expected Gaussian SD, or explicit historical outcome scaling")
        elif name == "observation_scale_floor":
            parser.add_argument("--observation-scale-floor", type=float, default=default,
                                help="numerical reference-scale floor in each channel's units, not measured noise")
        else:
            parser.add_argument("--" + name.replace("_", "-"), type=float, default=default)
    parser.add_argument("--decompose", action="store_true",
                        help="also run one-knob arms and attribute the NIS excess")
    parser.add_argument("--decompose-particles", type=int, nargs="+", default=[500, 2000, 8000],
                        help="ensemble sizes for the particle-approximation arms")
    parser.add_argument("--decompose-sigma", type=float, nargs="+", default=[2.0, 4.0, 8.0],
                        help="relative-sigma multipliers for the measurement-only arms")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path(os.environ.get("YSTWIN_OUTPUTS") or paths.REPO_ROOT / "outputs"))
    parser.add_argument("--no-write", action="store_true", help="compute and report, creating no artifacts")
    return parser


def _parse_args(argv=None):
    parser = _parser()
    preliminary, _ = parser.parse_known_args(argv)
    if preliminary.config is not None:
        payload = json.loads(preliminary.config.read_text(encoding="utf-8"))
        configuration = payload.get("configuration", payload)
        allowed = {action.dest for action in parser._actions} - {"help", "config"}
        if not isinstance(configuration, dict) or configuration.keys() - allowed:
            parser.error("config must contain only supported configuration keys")
        parser.set_defaults(**configuration)
    args = parser.parse_args(argv)
    for action in parser._actions:
        value = getattr(args, action.dest, None)
        if action.choices is not None and value not in action.choices:
            parser.error(f"invalid {action.dest}: {value!r}")
    for name in ("particles", "opening_points", "min_unique_particles"):
        value = getattr(args, name)
        if not isinstance(value, int) or value < (2 if name == "opening_points" else 1):
            parser.error(f"{name} must be a positive integer")
    if args.max_wells is not None and (not isinstance(args.max_wells, int) or args.max_wells < 1):
        parser.error("max_wells must be positive")
    if not np.isfinite(args.min_ess) or args.min_ess < 0:
        parser.error("min_ess must be finite and non-negative")
    if not np.isfinite(args.coverage) or not 0 < args.coverage < 1:
        parser.error("coverage must be in (0, 1)")
    if any(not isinstance(count, int) or count < 1 for count in args.decompose_particles):
        parser.error("decompose_particles must be positive integers")
    if any(not np.isfinite(m) or m <= 0 for m in args.decompose_sigma):
        parser.error("decompose_sigma must be positive and finite")
    _acceptance_band(1, args.alpha)
    args.output_dir = pathlib.Path(args.output_dir).expanduser().resolve()
    if args.plates_dir is not None:
        args.plates_dir = pathlib.Path(args.plates_dir).expanduser().resolve()
    return args


def _fingerprint(path: pathlib.Path) -> dict:
    return {"path": paths.display_path(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size}


def _resolve_sources(args):
    manifest = replay.load_manifest()
    if args.source_set not in set(manifest.source_set):
        raise ValueError(f"source_set {args.source_set!r} is absent from the committed manifest")
    exports = replay.exports(manifest=manifest, source_set=args.source_set)
    names = {p.name for p in exports}
    unknown = (set(args.plate) | set(args.exclude_plate)) - names
    if unknown:
        raise ValueError(f"unknown export identifiers in source_set: {sorted(unknown)}")
    plates = args.plates_dir
    if plates is None and args.source != "committed":
        plates = paths.biosensor_plates()
    if plates is not None and not (any(plates.glob("*.xlsx")) or any(plates.glob("*.xpt"))):
        # Set but empty is the same situation as unset: no workbooks to read.
        plates = None
    committed = args.source == "committed" or (args.source == "auto" and plates is None)
    if committed:
        # The workbooks are absent, which is the normal case for anyone who is not the
        # person who ran the plates -- they carry a named individual in their document
        # properties and are deliberately not committed. The numbers are committed, as
        # text, under data/plates, so this replays those rather than exiting. Nothing
        # downstream changes: `install` swaps the reader and the pipeline is unaware.
        # FILTERED, as `run_d2.py` is. An unfiltered `install` returns every committed
        # export -- the July preliminaries and the BY4741 control included -- so this script
        # read 4 plates on the machine holding the workbooks and 8 on a clone, and wrote a
        # tracked table with 1.7x the rows. That is the "reproducible on one machine" defect
        # this repository exists to refuse. Found by the 2026-08-30 audit.
        def reader(path):
            return replay.load_run(path.name, manifest=manifest)
        source = (f"the committed text in {replay.committed_dir()} "
                  f"({len(exports)} {args.source_set} exports, NOT the workbooks)")
    else:
        if plates is None:
            raise ValueError("workbooks requested but no workbook directory is available")
        def reader(path):
            return read_xpt_run(path) if path.suffix.lower() == ".xpt" else read_synergy_kinetic(path)
        exports = [plates / p.name for p in exports]
        source = f"the raw Synergy workbooks/Gen5 files in {plates}, manifest set {args.source_set}"
    return exports, reader, manifest, committed, source


def _matches_well(selectors: list[str], export: str, well: str) -> bool:
    return well in selectors or f"{export}::{well}" in selectors


def _collect(args):
    exports, reader, input_manifest, committed, source = _resolve_sources(args)
    overrides = {name: getattr(args, name) for name in _PROCESS_FIELDS}
    optics = ReporterOptics(args.gain, args.background, args.autofluorescence,
                            args.inner_filter_coeff, args.detector_max)
    ParticleFilter(TwinPriors((1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.3, 0.0), **overrides),
                   optics, args.gdcw_per_od, 0.0, n_particles=args.particles,
                   resample_threshold=args.resample_threshold, seed=args.seed)
    per_well, selections, well_priors, inputs, plate_details = [], [], [], [], []
    populations: dict[str, RecordedPlate | None] = {}
    seen_wells = set()
    for path in exports:
        export = path.name
        if (args.plate and export not in args.plate) or export in args.exclude_plate:
            selections.append(dict(plate=export, well=None, status="excluded", reason="plate_selection"))
            continue
        recorded = recorded_for_export(export)
        populations[export] = recorded
        roles = {} if recorded is None else recorded.well_roles
        seen_wells.update((export, well) for well in roles)
        if committed:
            entries = input_manifest[input_manifest.export == export]
            for name in entries.readings_file.drop_duplicates():
                inputs.append(_fingerprint(replay.committed_dir() / name))
        elif path.is_file():
            inputs.append(_fingerprint(path))
        try:
            run = reader(path)
            od_block, fl_block = run.raw_channel("OD600"), run.raw_channel("mCitrine")
            aligned = run.aligned(od_block.channel)
            od_frame, fl_frame = aligned[od_block.channel], aligned[fl_block.channel]
        except (ValueError, KeyError, FileNotFoundError) as exc:
            selections.append(dict(plate=export, well=None, status="excluded", reason=str(exc), stage="input"))
            for well in roles:
                selections.append(dict(plate=export, well=well, status="excluded",
                                       reason=f"plate_input_unavailable: {exc}", stage="input"))
            continue
        shared = [w for w in od_frame.columns if w in set(fl_frame.columns)]
        shared_set = set(shared)
        all_wells = set(od_frame.columns) | set(fl_frame.columns)
        seen_wells.update((export, well) for well in all_wells)
        # Identity survives absent channels, including wells absent from BOTH blocks.
        # A missing-channel row takes precedence over layout/blank/selection failures.
        for well in sorted((all_wells | set(roles)) - shared_set):
            selections.append(dict(plate=export, well=well, status="excluded",
                                   reason="missing_od_or_reporter_channel", stage="channels"))
        if args.blank_well:
            blanks = [w for w in shared if _matches_well(args.blank_well, export, w)]
        elif args.blank_policy == "recorded":
            blanks = [] if recorded is None else list(recorded.blank_wells)
        else:
            blanks = [w for w in shared if od_frame[w].median() < 0.15]
        # Explicit exploratory blank choices affect correction, never recorded identity.
        # Only exports with no record retain the exploratory shared-nonblank population.
        blank_set = set(blanks)
        cultures = [w for w in shared if w not in blank_set
                    and (recorded is None or roles.get(w) == "culture")]
        analysis_wells = blank_set | set(cultures)
        for well in shared:
            if well not in analysis_wells:
                selections.append(dict(
                    plate=export, well=well, status="excluded", stage="population",
                    reason="recorded_blank_not_selected" if roles.get(well) == "blank"
                    else "outside_recorded_layout"))
        blank_failure = None
        if not blanks or not blank_set <= shared_set:
            blank_failure = "recorded_or_selected_blanks_unavailable"
        else:
            od_blank = float(od_frame[blanks].to_numpy().mean())
            rfu_blank = float(fl_frame[blanks].to_numpy().mean())
            if not np.isfinite(od_blank) or not np.isfinite(rfu_blank):
                blank_failure = "nonfinite_blank"
        if blank_failure is not None:
            selections.append(dict(plate=export, well=None, status="excluded",
                                   reason=blank_failure, stage="blank"))
            for well in shared:
                if well in analysis_wells:
                    selections.append(dict(plate=export, well=well, status="excluded",
                                           reason=blank_failure, stage="blank"))
            continue
        times = od_frame.index.to_numpy(dtype=float)
        plate_details.append(dict(plate=export, od_channel=od_block.channel, rfu_channel=fl_block.channel,
                                  blank_wells=blanks, od_blank=od_blank, rfu_blank=rfu_blank,
                                  n_times=len(times), times_h=times.tolist()))
        for well in blanks:
            selections.append(dict(plate=export, well=well, status="blank", reason="blank_reference"))
        # Enumerate every eligible well. A prefix of the well list is not a sample: on a
        # plate laid out with dose along rows, the first wells are the lowest doses,
        # where the model is nearest right -- and that subset reports the opposite
        # verdict from the full plate.
        for well in cultures:
            reason, field = None, None
            if (args.well and not _matches_well(args.well, export, well)) or _matches_well(
                    args.exclude_well, export, well):
                reason = "well_selection"
            elif args.max_wells is not None and len(per_well) >= args.max_wells:
                reason = "smoke_test_cap"
            stage = "selection" if reason else None
            od = od_frame[well].to_numpy() - od_blank
            rfu = fl_frame[well].to_numpy() - rfu_blank
            if reason is None:
                try:
                    frame = _run_well(
                        times, od, rfu, args.particles, args.seed, prior_overrides=overrides,
                        optics=optics, opening_points=args.opening_points, prior_policy=args.prior_policy,
                        include_initialization=args.include_initialization,
                        resample_threshold=args.resample_threshold, gdcw_per_od=args.gdcw_per_od,
                        carotenoid=args.carotenoid, coverage=args.coverage, reading_policy=args.reading_policy,
                    )
                except ReadingRefusal as exc:
                    reason, stage, field = f"reading_refusal: {exc}", "reading", exc.field
                except BelowBackground as exc:
                    reason, stage, field = f"initialization_refusal: {exc}", "initialization", exc.field
                except ValueError as exc:
                    reason, stage = f"initialization_or_filter_refusal: {exc}", "initialization_or_filter"
                else:
                    if frame.empty:
                        reason, stage = "no_measured_channels", "filter"
                    else:
                        well_priors.append(dict(plate=export, well=well, seed=args.seed, **frame.attrs))
                        frame["well"], frame["plate"] = well, export
                        frame.attrs.clear()
                        per_well.append(frame)
            selections.append(dict(plate=export, well=well, status="excluded" if reason else "included",
                                   reason=reason or "", stage=stage, field=field))
    for selector in args.well + args.exclude_well + args.blank_well:
        if not any(_matches_well([selector], export, well) for export, well in seen_wells):
            raise ValueError(f"well selector {selector!r} matches no available or recorded well")
    # One identity per well, independent of whether it supplied blanks, failed, or ran.
    # Plate-level rows have no well identity and are not biological denominator entries.
    for entry in selections:
        recorded = populations.get(entry["plate"])
        if entry["well"] is None:
            entry["recorded_role"] = None
        elif recorded is None:
            entry["recorded_role"] = "unknown"
        else:
            entry["recorded_role"] = recorded.well_roles.get(entry["well"], "unrecorded")
    innovations = pd.concat(per_well, ignore_index=True) if per_well else pd.DataFrame()
    summary = _summarise(innovations, min_ess=args.min_ess, min_unique_particles=args.min_unique_particles,
                         alpha=args.alpha, coverage=args.coverage)
    channels = pd.DataFrame(_channel_rows(summary, args.aggregation, args.aggregation_population))
    ledger = pd.DataFrame(selections, columns=["plate", "well", "status", "reason", "stage", "field",
                                               "recorded_role"])
    configuration = {key: str(value) if isinstance(value, pathlib.Path) else value
                     for key, value in vars(args).items() if key != "config"}
    code_paths = [pathlib.Path(__file__).resolve()] + [paths.PACKAGE_ROOT / name for name in (
        "estimator.py", "analysis/estimator_accuracy.py", "observation.py", "reporter.py",
        "growth.py", "plate/replay.py", "plate/synergy.py", "plate/gen5.py", "plate/layout.py", "calib/od.py",
    )]
    prequential_noise = args.observation_noise_mode == "predicted_scale"
    finite_readings = args.reading_policy == "finite_gaussian"
    manifest = {
        "schema_version": 4, "configuration": configuration,
        "resolved_source": "committed" if committed else "workbooks", "source_description": source,
        "input_manifest": _fingerprint(replay.committed_dir() / replay.MANIFEST),
        "inputs": inputs, "code": [_fingerprint(p) for p in code_paths],
        "software": {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy")},
        "python": sys.version, "optics": asdict(optics), "plate_details": plate_details,
        "initialized_wells": well_priors, "selection": selections,
        "population": {
            "scope": "source_set after plate selection; before well selection, smoke cap or data refusals",
            "recorded": "recorded_for_export: only RecordedPlate.culture_wells are culture candidates; "
                        "selected blank references are not scored",
            "outside_layout": "never culture candidates, including when explicitly selected by --well",
            "unknown_layout": "shared nonblank wells are exploratory only; explicit blank selectors or "
                              "legacy_low_od are required without recorded blanks",
            "ledger": "one row per (export, well), retaining absent recorded wells; well=null is a plate-level "
                      "exclusion, not a biological denominator entry; recorded_role is independent of status",
            "plates": [dict(plate=export,
                            recorded_culture_wells=None if record is None else list(record.culture_wells),
                            recorded_blank_wells=None if record is None else list(record.blank_wells))
                       for export, record in populations.items()],
        },
        "measurement_noise": {
            "mode": args.observation_noise_mode,
            "family": ("conditionally_independent_additive_gaussian" if prequential_noise
                       else "outcome_dependent_legacy_pseudo_likelihood"),
            "sigma": ("rel_sigma * max(abs(particle_expected_measurement), observation_scale_floor)"
                      if prequential_noise else "rel_sigma * max(abs(observed), observation_scale_floor)"),
            "relative_sigma": {"od": args.od_rel_sigma, "rfu": args.rfu_rel_sigma},
            "observation_scale_floor": args.observation_scale_floor,
            "floor_role": "numerical reference-scale floor in each channel's units; not an instrument noise measurement",
            "log_likelihood": ("-0.5 * ((observed - predicted_i) / sigma_i)**2 - log(sigma_i) - 0.5 * log(2*pi)"
                               if prequential_noise else
                               "-0.5 * ((observed - predicted_i) / sigma)**2; normalization omitted for legacy bit preservation"),
            "predictive_measurement_variance": ("sum_i prior_weight_i * sigma_i**2" if prequential_noise
                                                else "sigma(observed)**2; outcome-dependent, not prior-predictive"),
            "prior_predictive": prequential_noise,
        },
        "reading_policy": {
            "name": args.reading_policy,
            "finite_values": ("retain signed OD and RFU, including zero and negative readings"
                              if finite_readings else "require strictly positive OD throughout; retain signed RFU"),
            "nan": ("omit only the missing channel at its original step" if finite_readings
                    else "refuse well for missing OD; omit only missing RFU at its original step"),
            "infinity": "refuse well with field and original step indices",
            "initialization": "positive physical-prior initialization is a separate requirement, not a reading filter",
        },
        "blank_correction": "constant mean across recorded/selected blank wells and all times",
        "alignment": "SynergyRun.aligned: reporter interpolated to OD timestamps",
        "initialization": "finite positive opening OD for biomass/log-growth and finite opening RFU with positive mean concentration; "
                          "legacy_full_trace additionally needs positive full-trace OD; initialization rows labelled and normally unscored",
        "ess_semantics": "ess is posterior/pre-resampling; post_resample_ess is not particle diversity",
        "ancestor_semantics": "number of surviving initial particle labels; not restored by process noise",
        "intervals": "one-step Gaussian moment approximation to the particle mixture, not exact mixture quantiles",
        "aggregation_caveat": "per-well K bands; descriptive channel aggregates have no chi-squared null",
        "prospective_prior": args.prior_policy == "opening" and not args.include_initialization,
        "smoke_test": args.max_wells is not None,
        "counts": {"included_wells": len(per_well), "well_channels": len(summary),
                   "excluded_entries": int((ledger.status == "excluded").sum()),
                   "recorded_culture_wells": int((ledger.recorded_role == "culture").sum()),
                   "recorded_blank_wells": int((ledger.recorded_role == "blank").sum()),
                   "outside_recorded_layout_wells": int((ledger.recorded_role == "unrecorded").sum()),
                   "unknown_identity_wells": int((ledger.recorded_role == "unknown").sum())},
        "outputs": [] if args.no_write else [str(args.output_dir / name) for name in _OUTPUT_NAMES],
    }
    return innovations, summary, channels, ledger, manifest


def _write_outputs(args, innovations, summary, channels, ledger, manifest, decomposition=None):
    """The only write boundary; --no-write never reaches it."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in zip(_OUTPUT_NAMES, (innovations, summary, channels, ledger)):
        table.to_csv(args.output_dir / name, index=False)
    if decomposition is not None:
        decomposition.to_csv(args.output_dir / DECOMPOSITION_NAME, index=False)
    (args.output_dir / "nis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    args = _parse_args(argv)
    try:
        innovations, summary, channels, ledger, manifest = _collect(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f"calibration refused: {exc}", file=sys.stderr)
        return 2
    print(f"SOURCE: {manifest['source_description']}")
    if summary.empty:
        print("no eligible wells found in the resolved plates", file=sys.stderr)
        return 2
    print("Filter calibration: one-step predictions, not a test of a fitted posterior centre")
    print(f"Observation noise: {args.observation_noise_mode}; {manifest['measurement_noise']['sigma']}")
    print(f"Reading policy: {args.reading_policy}; {manifest['reading_policy']['finite_values']}; "
          f"{manifest['reading_policy']['nan']}; infinity refuses the well")
    if not manifest["measurement_noise"]["prior_predictive"]:
        print("NON-PREQUENTIAL: observed_scale_legacy noise uses the outcome being scored")
    if args.max_wells:
        print(f"--max-wells={args.max_wells}: smoke test, NOT a representative result")
    if not manifest["prospective_prior"]:
        print("NON-PROSPECTIVE: future-derived priors or scored initialization were explicitly requested")
    print(f"{manifest['counts']['included_wells']} distinct (plate, well) pairs; "
          f"{summary.plate.nunique()} plates; {args.particles} particles; ESS floor {args.min_ess:g}")
    print(f"{manifest['counts']['excluded_entries']} excluded entries, recorded with reasons")
    for row in channels.to_dict("records"):
        untested = row["well_channels"] - row["tested_well_channels"]
        print(f"{row['channel']}: median mean NIS {row['median_mean_nis']:.3g}; "
              f"{args.aggregation}={row['aggregate_mean_nis']:.3g} ({args.aggregation_population}); "
              f"lag-1 {row['median_lag1_autocorr']:+.3f}; coverage {row['median_coverage']:.1%} "
              f"(nominal {args.coverage:.1%}); {row['verdict']}")
        if untested:
            # The floor is on the diversity left AFTER the reading is weighted in, which
            # is the number that can be small; before 2026-09-08 it could not be.
            print(f"  {untested}/{row['well_channels']} well-channels INCONCLUSIVE: ESS/diversity gate")
    print("Per-well bands use each well's scored K; channel aggregates are descriptive, not a chi-squared test.")
    print("Historical 1.76 / 4.44 and the old deceleration docstring sweep lack executable run provenance.")
    decomposition = None
    written = list(_OUTPUT_NAMES)
    if args.decompose:
        plan = _arm_plan(args)
        print(f"Error decomposition: {len(plan)} arms, one knob each, baseline "
              f"{args.particles} particles at the declared sigma")
        try:
            decomposition = _decompose(args, innovations, channels)
        except ValueError as exc:
            print(f"decomposition refused: {exc}", file=sys.stderr)
            return 2
        manifest["error_decomposition"] = _decomposition_manifest(args, plan)
        for channel, block in decomposition.groupby("channel", sort=True):
            base = block[block.family == "baseline"].iloc[0]
            needed = (f"a measurement-only account needs sigma x{base.implied_sigma_multiplier:.3g}"
                      if base.baseline_has_excess else "no excess over the target to apportion")
            print(f"{channel}: baseline mean NIS {base.mean_nis:.4g}; measurement is "
                  f"{base.measurement_variance_share:.1%} of predictive variance, so "
                  f"{needed}")
            for row in block[block.family != "baseline"].to_dict("records"):
                print(f"  {row['label']:<14} {row['family']:<12} mean NIS {row['mean_nis']:>9.4g}  "
                      f"excess removed {row['excess_removed']:+8.1%}  "
                      f"min ESS {row['min_ess']:.3g} (x{row['min_ess_ratio']:.3g})  "
                      f"ancestors {row['min_unique_ancestors']:.3g}")
        print("  Unremoved excess is process-model error only by elimination among these three.")
        entangled = decomposition[(decomposition.family == "measurement")
                                  & (decomposition.min_ess_ratio.sub(1.0).abs() > 0.05)]
        if not entangled.empty:
            print(f"  {len(entangled)} measurement arm-channels also moved the ensemble "
                  f"(min ESS ratio up to x{entangled.min_ess_ratio.max():.3g}); a looser "
                  "likelihood keeps particles alive, so those two families are not separated here.")
        if not args.no_write:
            written.append(DECOMPOSITION_NAME)
            manifest["outputs"] = [str(args.output_dir / name) for name in written]
    if args.no_write:
        print("--no-write: no tables or manifest written")
    else:
        _write_outputs(args, innovations, summary, channels, ledger, manifest, decomposition)
        print(f"Outputs: {args.output_dir} ({', '.join(written)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
