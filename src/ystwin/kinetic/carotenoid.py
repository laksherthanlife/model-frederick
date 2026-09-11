"""Kinetic GGPP -> beta-carotene branch (2 GGPP -crtYB-> phytoene -crtI-> lycopene -crtYB-> product).

No longer parked. Elizondo & Saa 2025 (PMID 40891387) published six chemostat steady states
on three beta-carotene CEN.PK2-1c strains, with lycopene and beta-carotene rates reported
*separately*, and that is enough to fit part of this branch. Which part, and only which
part, is the whole content of ``docs/research/KINETIC_FIT.md``:

  * **The cyclase step is identifiable.** Carotenoids are not secreted, so at chemostat
    steady state each measured rate is the growth dilution of an intracellular pool and the
    lycopene *content* is therefore a measurement, not a fitted quantity. That gives the
    cyclase both its substrate concentration and its rate, six times over.
  * **The synthase and desaturase steps are not.** Phytoene and GGPP were not measured, so
    their pools are free, and any (vmax, km) pair can carry the required flux by moving an
    unobserved pool. :func:`calibrated_kinetics` refuses to invent them rather than
    returning a plausible number.

What the identifiable half says is not what this module was built assuming. A cyclase whose
vmax does not depend on growth rate is *refuted* by the six states -- in every strain the
higher lycopene content sits with the lower beta-carotene rate, which Michaelis-Menten
forbids -- and every growth-rate-independent form scores worse than predicting the training
mean. The form that survives is a cyclase capacity proportional to growth rate, which is
the same statement as a growth-rate-independent relation between the two *contents*. See
:data:`ELIZONDO2025` and the module's tests.

crtYB is bifunctional, so its two domains share one enzyme pool; intracellular pools are
diluted by growth. Phytoene, lycopene and beta-carotene elute from one HPLC run, so their
ratios name the limiting enzyme at almost no marginal cost.
"""

from __future__ import annotations

import math

import pathlib
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from ystwin import paths

__all__ = [
    "solve_branch_from_flux",
    "BranchSolution",
    "BranchSteadyState",
    "CarotenoidKinetics",
    "CyclaseCalibration",
    "ELIZONDO2025",
    "UNIDENTIFIABLE_STEPS",
    "beta_carotene_content",
    "branch_rates",
    "calibrated_cyclase",
    "calibrated_kinetics",
    "calibration_states",
    "limiting_step",
    "predict_beta_carotene_rate",
    "simulate_branch",
    "steady_state_pools",
]

_SPECIES = ("ggpp", "phytoene", "lycopene", "beta_carotene")

UNIDENTIFIABLE_STEPS = ("psy", "crti")
"""Steps no measurement in the calibration set constrains.

Elizondo measured lycopene and beta-carotene. Phytoene and GGPP were not measured, so the
substrate concentration of the synthase and of the desaturase is unobserved in all six
states. Their rates are pinned -- the desaturase carries ``q_lycopene + q_betacarotene``
and the synthase carries that plus phytoene dilution -- but a rate without its substrate
concentration determines no (vmax, km) pair: any pair delivers that rate at some pool size.
Named here so that a caller asking for them gets a list rather than a surprise.
"""


@dataclass(frozen=True)
class CyclaseCalibration:
    """A fitted lycopene-cyclase law, with what it was fitted to and how it scored.

    **Every number here is FITTED, not looked up.** No turnover number is on record for
    phytoene synthase in any organism, nor for the *Xanthophyllomyces dendrorhous* crtE,
    crtYB or crtI these strains carry, so there was nothing to look up. Elizondo's own team
    hit the same wall and sampled a thermodynamically-constrained ensemble (ABC-GRASP)
    rather than parameterising from a database. ``docs/research/KINETIC_FIT.md`` §5 records
    that. These constants are Tier 1 under ``docs/CLAIM_BOUNDARY.md``: a deterministic fit
    to real measurements, scored out of sample, but not a registered prediction.

    The law is a saturating relation between the two intracellular *contents*, with no
    explicit growth-rate term:

        [beta-carotene] = capacity * [lycopene] / (km + [lycopene])

    and the rate the twin needs follows by multiplying by the growth rate, because at
    chemostat steady state a non-secreted product leaves only in the biomass.

    Args:
        capacity_mmol_per_gdcw: The beta-carotene content approached as lycopene saturates
            the cyclase. Equivalently ``vmax_lcy / growth_rate``.
        km_mmol_per_gdcw: Lycopene content at half that value.
        capacity_ci95, km_ci95: 95% intervals from the fit's own residual scatter, which is
            the wider of the two uncertainty estimates and so the one carried here.
        log_covariance: The fit's 2x2 covariance of ``(log capacity, log km)``, in that
            order. The two intervals above are its diagonal; this is the whole matrix.
            **The two parameters MUST NOT be propagated independently** -- see below.
        growth_rate_range: The dilution rates fitted. Outside it the law is an
            extrapolation and :func:`calibrated_cyclase` refuses.
        n_states, n_parameters: Six observations, two free parameters.
        loo_rmse_log, baseline_rmse_log: Leave-one-out error and the training-mean baseline.
        source: The paper, so the number and its provenance travel together.

    CAPACITY AND KM LIE ON THE FIT'S SATURATION RIDGE AND MUST NOT BE PROPAGATED
    INDEPENDENTLY. Their log-correlation is +0.852: a larger capacity is compensated by a
    larger km along a direction the six states barely constrain, which is the same ridge
    ``scripts/fit_carotenoid_kinetics.py``'s multi-start exists to keep the optimiser off.

    CORRECTION OF 2026-09-04, and it is why :attr:`log_covariance` exists. The fit computed
    the whole 2x2 matrix, `uncertainty()` returned its off-diagonal as ``log_correlation``,
    and ``main`` PRINTED it -- and this dataclass stored only the diagonal, as two intervals
    side by side with nothing saying they were related. Two intervals side by side are an
    invitation to add their variances, and adding them is wrong here: the two log parameters
    enter predicted beta-carotene content with OPPOSITE signs (d log content / d log capacity
    is +0.93 to +0.60 over mu = 0.101 to 0.254, d/d log km is -0.13 to -0.25), so their errors
    partly cancel. Measured through the shipped :func:`solve_branch_from_flux` at the median
    calibration flux, independent propagation overstates the parameter log-sd by 1.48x, 1.97x
    and 2.60x at those three growth rates -- a factor of 2.2 to 6.7 in variance.
    :meth:`joint_log_sd` exists so the correct propagation is the one a consumer reaches for
    first, rather than the one that takes an extra thought.
    """

    capacity_mmol_per_gdcw: float
    km_mmol_per_gdcw: float
    capacity_ci95: tuple[float, float]
    km_ci95: tuple[float, float]
    log_covariance: tuple[tuple[float, float], tuple[float, float]]
    growth_rate_range: tuple[float, float]
    n_states: int
    n_parameters: int
    loo_rmse_log: float
    baseline_rmse_log: float
    source: str

    @property
    def skill_vs_training_mean(self) -> float:
        """Fraction of the baseline's squared log error the fit removes out of sample.

        Negative means the fit is worse than predicting the mean of the training states,
        which is the outcome for every growth-rate-independent form tried.
        """
        return 1.0 - (self.loo_rmse_log / self.baseline_rmse_log) ** 2

    @property
    def log_parameter_correlation(self) -> float:
        """Correlation of ``log capacity`` and ``log km``, off :attr:`log_covariance`.

        Derived rather than stored, so it cannot disagree with the matrix it comes from.
        """
        (vaa, vab), (_, vbb) = self.log_covariance
        return vab / math.sqrt(vaa * vbb)

    def log_parameter_sd(self) -> tuple[float, float]:
        """``(sd log capacity, sd log km)`` -- the diagonal, as standard deviations."""
        (vaa, _), (_, vbb) = self.log_covariance
        return math.sqrt(vaa), math.sqrt(vbb)

    def joint_log_sd(self, gradient: tuple[float, float]) -> float:
        """Propagate both parameters TOGETHER through a log-log sensitivity.

        Args:
            gradient: ``(d log y / d log capacity, d log y / d log km)`` for whatever ``y``
                the caller is propagating into. Finite differences on the shipped solver are
                the usual way to get it; nothing here assumes a particular ``y``.

        Returns:
            ``sqrt(g' C g)``, the log standard deviation of ``y`` from these two parameters.

        This method exists because the arithmetic it replaces -- adding the two variances --
        is what a dataclass carrying two intervals and no covariance invites, and it is wrong
        by a factor of 2.2 to 6.7 in variance on this fit. See the class docstring.
        """
        (vaa, vab), (vba, vbb) = self.log_covariance
        ga, gb = gradient
        variance = ga * ga * vaa + gb * gb * vbb + ga * gb * (vab + vba)
        if variance < 0.0:  # only reachable if a caller supplies a non-PSD matrix
            raise ValueError(
                f"the log covariance is not positive semi-definite along {gradient!r}; "
                f"it came out of a fit and should not be edited by hand")
        return math.sqrt(variance)

    def independent_log_sd(self, gradient: tuple[float, float]) -> float:
        """The WRONG propagation, named so a test can measure how wrong it is.

        Kept rather than left implicit: the size of the error a consumer would make is the
        argument for :meth:`joint_log_sd`, and an argument nothing can evaluate goes stale.
        """
        (vaa, _), (_, vbb) = self.log_covariance
        ga, gb = gradient
        return math.sqrt(ga * ga * vaa + gb * gb * vbb)


ELIZONDO2025 = CyclaseCalibration(
    capacity_mmol_per_gdcw=2.3252e-3,
    km_mmol_per_gdcw=5.974e-4,
    capacity_ci95=(1.7563e-3, 3.0784e-3),
    km_ci95=(2.961e-4, 1.2051e-3),
    # The SAME `uncertainty()` call that produced the two intervals above. Carried whole
    # rather than as its diagonal -- see the class docstring's 2026-09-04 correction.
    log_covariance=((0.01021528, 0.02176597), (0.02176597, 0.06387742)),
    growth_rate_range=(0.100987987, 0.254320182),
    n_states=6,
    n_parameters=2,
    loo_rmse_log=0.1858,
    baseline_rmse_log=0.4372,
    source="PMID 40891387 (doi 10.1021/acssynbio.5c00256); fitted by scripts/fit_carotenoid_kinetics.py",
)
"""The one fitted law this dataset supports. Regenerate with ``scripts/fit_carotenoid_kinetics.py``."""


@dataclass(frozen=True)
class CarotenoidKinetics:
    """Michaelis-Menten parameters for the three enzymatic steps.

    Rates are mmol / gDCW / h; half-saturation constants are mmol / gDCW.

    Args:
        vmax_psy, km_psy: Phytoene synthase, the crtYB synthase domain.
        vmax_crti, km_crti: Phytoene desaturase.
        vmax_lcy, km_lcy: Lycopene cyclase, the crtYB cyclase domain. These two are the
            only ones the calibration set constrains, and only at a stated growth rate --
            see :func:`calibrated_cyclase`.
        crtyb_pool: Combined capacity of the two crtYB domains. When set, the two
            activities are scaled down proportionally if their unconstrained sum
            would exceed it. ``None`` treats them as independent enzymes, which is
            the wrong model for a bifunctional protein but useful as a contrast.
    """

    vmax_psy: float
    km_psy: float
    vmax_crti: float
    km_crti: float
    vmax_lcy: float
    km_lcy: float
    crtyb_pool: float | None = None

    def __post_init__(self) -> None:
        values = {
            "vmax_psy": self.vmax_psy, "km_psy": self.km_psy,
            "vmax_crti": self.vmax_crti, "km_crti": self.km_crti,
            "vmax_lcy": self.vmax_lcy, "km_lcy": self.km_lcy,
        }
        bad = [k for k, v in values.items() if v <= 0]
        if bad:
            raise ValueError(f"kinetic parameters must be positive: {', '.join(sorted(bad))}")
        if self.crtyb_pool is not None and self.crtyb_pool <= 0:
            raise ValueError("crtyb_pool must be positive, or None for independent enzymes")


@dataclass(frozen=True)
class BranchSteadyState:
    """Pools and fluxes once the branch has equilibrated at a fixed supply."""

    ggpp: float
    phytoene: float
    lycopene: float
    beta_carotene: float
    psy_flux: float
    crti_flux: float
    lcy_flux: float
    growth_rate: float

    @property
    def beta_carotene_flux(self) -> float:
        """Net rate of beta-carotene formation, mmol / gDCW / h."""
        return self.lcy_flux

    @property
    def lycopene_accumulation(self) -> float:
        """Rate at which lycopene leaves in the biomass, mmol / gDCW / h.

        This -- not :attr:`crti_flux` -- is what an HPLC of the pellet measures, because
        carotenoids are not secreted and a chemostat's only outlet for them is washout with
        the cells. It equals ``crti_flux - lcy_flux`` at steady state, so a strain that
        cyclises everything reports zero lycopene while still carrying full desaturase flux.
        Reading the measured lycopene rate as a desaturase flux instead double-counts the
        carbon that went on to beta-carotene.
        """
        return self.growth_rate * self.lycopene

    @property
    def beta_carotene_accumulation(self) -> float:
        """The measured beta-carotene rate, mmol / gDCW / h.

        Equal to :attr:`beta_carotene_flux` at steady state, and kept separate because that
        equality is a property of the steady state rather than of the pathway: away from
        steady state the flux and the washout differ, and the twin's dynamic mode reports
        both.
        """
        return self.growth_rate * self.beta_carotene


def _mm(substrate: float, vmax: float, km: float) -> float:
    s = max(substrate, 0.0)
    return vmax * s / (km + s)


def branch_rates(pools: dict[str, float], kinetics: CarotenoidKinetics) -> dict[str, float]:
    """Instantaneous enzymatic rates at the given intracellular concentrations.

    Args:
        pools: Concentrations, mmol/gDCW, keyed ``ggpp``, ``phytoene``, ``lycopene``.
        kinetics: Enzyme parameters.

    Returns:
        Rates keyed ``psy``, ``crti``, ``lcy``, mmol/gDCW/h. ``psy`` is the rate of
        phytoene formation, so it consumes twice that much GGPP.
    """
    psy = _mm(pools.get("ggpp", 0.0), kinetics.vmax_psy, kinetics.km_psy)
    crti = _mm(pools.get("phytoene", 0.0), kinetics.vmax_crti, kinetics.km_crti)
    lcy = _mm(pools.get("lycopene", 0.0), kinetics.vmax_lcy, kinetics.km_lcy)

    if kinetics.crtyb_pool is not None:
        demand = psy + lcy
        if demand > kinetics.crtyb_pool:
            # One protein cannot run both domains beyond its total capacity: they share
            # one pool, so both are scaled down together.
            scale = kinetics.crtyb_pool / demand
            psy *= scale
            lcy *= scale
    return {"psy": psy, "crti": crti, "lcy": lcy}


def _derivatives(pools_vec, supply, kinetics, growth_rate):
    pools = dict(zip(_SPECIES, pools_vec))
    r = branch_rates(pools, kinetics)
    return np.array([
        supply - 2.0 * r["psy"] - growth_rate * pools["ggpp"],
        r["psy"] - r["crti"] - growth_rate * pools["phytoene"],
        r["crti"] - r["lcy"] - growth_rate * pools["lycopene"],
        r["lcy"] - growth_rate * pools["beta_carotene"],
    ])


def simulate_branch(
    times_h: np.ndarray,
    ggpp_supply,
    kinetics: CarotenoidKinetics,
    growth_rate,
    initial_pools: dict[str, float] | None = None,
) -> dict[str, np.ndarray]:
    """Integrate the branch under a given precursor supply.

    Args:
        times_h: Ascending time grid, hours.
        ggpp_supply: GGPP supply rate, mmol/gDCW/h. Scalar or per-timepoint array;
            in the coupled twin this comes from the GSMM solve.
        kinetics: Enzyme parameters.
        growth_rate: Specific growth rate, scalar or per-timepoint array.
        initial_pools: Starting concentrations; defaults to zero.

    Returns:
        One array per species plus ``beta_carotene_flux``, on the supplied grid.
    """
    t = np.asarray(times_h, dtype=float)
    supply = np.broadcast_to(np.asarray(ggpp_supply, dtype=float), t.shape)
    mu = np.broadcast_to(np.asarray(growth_rate, dtype=float), t.shape)
    start = initial_pools or {}
    y0 = np.array([float(start.get(s, 0.0)) for s in _SPECIES])

    def rhs(x, y):
        return _derivatives(
            y, float(np.interp(x, t, supply)), kinetics, float(np.interp(x, t, mu))
        )

    sol = solve_ivp(rhs, (t[0], t[-1]), y0, t_eval=t, method="LSODA", rtol=1e-9, atol=1e-14)
    if not sol.success:
        raise RuntimeError(f"carotenoid branch integration failed: {sol.message}")

    out = {s: sol.y[i] for i, s in enumerate(_SPECIES)}
    out["beta_carotene_flux"] = np.array([
        branch_rates({k: out[k][i] for k in _SPECIES}, kinetics)["lcy"] for i in range(t.size)
    ])
    return out


def steady_state_pools(
    ggpp_supply: float,
    kinetics: CarotenoidKinetics,
    growth_rate: float,
) -> BranchSteadyState:
    """Solve the branch for its steady state at a fixed supply.

    Args:
        ggpp_supply: GGPP supply rate, mmol/gDCW/h.
        kinetics: Enzyme parameters.
        growth_rate: Specific growth rate, 1/h.
    """
    if ggpp_supply <= 0:
        return BranchSteadyState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, growth_rate)

    # Integrate to steady state rather than solving the coupled algebra: the shared
    # crtYB pool couples the two branches, so there is no closed form.
    t = np.linspace(0.0, 5000.0, 20000)
    traj = simulate_branch(t, ggpp_supply, kinetics, growth_rate)
    final = {s: float(traj[s][-1]) for s in _SPECIES}
    rates = branch_rates(final, kinetics)
    return BranchSteadyState(
        ggpp=final["ggpp"],
        phytoene=final["phytoene"],
        lycopene=final["lycopene"],
        beta_carotene=final["beta_carotene"],
        psy_flux=rates["psy"],
        crti_flux=rates["crti"],
        lcy_flux=rates["lcy"],
        growth_rate=growth_rate,
    )


def limiting_step(state: BranchSteadyState) -> str:
    """Name the enzyme whose substrate has backed up the most.

    An intermediate accumulates immediately upstream of the slow step, so the
    largest pool names the bottleneck. This is the reading an HPLC trace supports
    directly, with no additional assay.
    """
    upstream = {
        "crtYB_PSY": state.ggpp,
        "crtI": state.phytoene,
        "crtYB_LCY": state.lycopene,
    }
    return max(upstream, key=upstream.get)


# --------------------------------------------------------------------------- #
# the calibrated cyclase
# --------------------------------------------------------------------------- #


def beta_carotene_content(
    lycopene_content: float,
    calibration: CyclaseCalibration = ELIZONDO2025,
) -> float:
    """Steady-state beta-carotene content from lycopene content, both mmol/gDCW.

    The fitted law, stated in the space it was fitted in. No growth rate appears, which is
    the finding: the same two constants reproduce all six states across three strains and a
    2.5-fold range of dilution rate.

    Args:
        lycopene_content: Intracellular lycopene, mmol/gDCW. In a chemostat this is the
            measured ``q_lycopene / growth_rate``.
        calibration: Which fit to use.

    Raises:
        ValueError: on a negative content.
    """
    if lycopene_content < 0:
        raise ValueError(f"lycopene content cannot be negative, got {lycopene_content}")
    c, km = calibration.capacity_mmol_per_gdcw, calibration.km_mmol_per_gdcw
    return c * lycopene_content / (km + lycopene_content)


@dataclass(frozen=True)
class BranchSolution:
    """Both carotenoid pools, solved from the pathway flux rather than assumed.

    Args:
        lycopene_content: Steady-state lycopene, mmol/gDCW. A PREDICTION.
        beta_carotene_content: Steady-state beta-carotene, mmol/gDCW.
        beta_carotene_rate: Specific production rate, mmol/gDCW/h.
        pathway_flux: The desaturase flux that produced them, mmol/gDCW/h.
        growth_rate: The rate it was solved at, /h.
    """

    lycopene_content: float
    beta_carotene_content: float
    beta_carotene_rate: float
    pathway_flux: float
    growth_rate: float


def solve_branch_from_flux(
    pathway_flux: float,
    growth_rate: float,
    calibration: CyclaseCalibration = ELIZONDO2025,
) -> BranchSolution:
    """Solve both pools from the total carotenoid flux. Lycopene is an output here.

    An earlier version of the chain took lycopene content as an INPUT, on the reasoning
    that FBA bounds the product from above and never from below so nothing pins it. That
    was giving up too early. Lycopene is an intermediate, and an intermediate at steady
    state is pinned by a mass balance, not by an objective:

        d[lyc]/dt = v_desaturase - v_cyclase - mu*[lyc] = 0

    With the cyclase Michaelis-Menten in [lyc] that is a quadratic with one positive root,

        mu*[lyc]^2 + (mu*Km + Vmax - v)*[lyc] - v*Km = 0

    so the pool follows from the flux and the fitted cyclase alone. No PSY or crtI
    parameters are needed, which matters because those two are the steps the calibration
    set cannot identify -- see :data:`UNIDENTIFIABLE_STEPS`. Solved this way lycopene
    reproduces all six measured states to a median 5.4%.

    And the flux itself is a STRAIN constant rather than an environment one. Across the
    2.52-fold growth change in the calibration set it moves by 1.09x, 0.60x and 0.89x in
    the three strains -- Pearson(mu, flux) is -0.215. A heterologous pathway on
    constitutive promoters does not know what the dilution rate is. So one number
    characterises a strain, and it is the number a gene-expression measurement supplies.

    Args:
        pathway_flux: Total carotenoid flux through the desaturase, mmol/gDCW/h. At
            steady state this equals ``q_lycopene + q_beta_carotene``, both measurable.
        growth_rate: Specific growth rate, /h.
        calibration: Fitted cyclase parameters.

    Raises:
        ValueError: on a negative flux, or a growth rate outside the calibrated range.
    """
    if not math.isfinite(pathway_flux) or pathway_flux < 0:
        raise ValueError("pathway_flux must be finite and non-negative")
    if not math.isfinite(growth_rate) or growth_rate <= 0:
        raise ValueError("growth_rate must be finite and positive; a non-growing cell has no steady state")
    _check_growth_rate(growth_rate, calibration)

    vmax, km = calibrated_cyclase(growth_rate, calibration)
    a = growth_rate
    b = growth_rate * km + vmax - pathway_flux
    c = -pathway_flux * km
    root = math.hypot(b, 2.0 * math.sqrt(-a * c))
    lycopene = -2.0 * c / (b + root) if b >= 0.0 else (root - b) / (2.0 * a)
    cyclase_rate = vmax * lycopene / (km + lycopene)
    return BranchSolution(
        lycopene_content=float(lycopene),
        beta_carotene_content=float(cyclase_rate / growth_rate),
        beta_carotene_rate=float(cyclase_rate),
        pathway_flux=float(pathway_flux),
        growth_rate=float(growth_rate),
    )


def predict_beta_carotene_rate(
    lycopene_content: float,
    growth_rate: float,
    calibration: CyclaseCalibration = ELIZONDO2025,
) -> float:
    """Specific beta-carotene production rate, mmol/gDCW/h.

    The quantity a pellet HPLC plus a dilution rate measures, and the one the twin has to
    forecast. It is :func:`beta_carotene_content` times the growth rate, because a product
    that is never secreted leaves a chemostat only inside the cells.

    Args:
        lycopene_content: Intracellular lycopene, mmol/gDCW.
        growth_rate: Specific growth rate, 1/h. Must lie in the calibrated range.
        calibration: Which fit to use.

    Raises:
        ValueError: outside the calibrated growth-rate range. See :func:`calibrated_cyclase`.
    """
    _check_growth_rate(growth_rate, calibration)
    return growth_rate * beta_carotene_content(lycopene_content, calibration)


def _check_growth_rate(growth_rate: float, calibration: CyclaseCalibration) -> None:
    """Refuse a growth rate the fit cannot speak to.

    Two dilution rates were measured. Between them the growth-rate dependence is an
    interpolation over an interval whose endpoints are both data; outside them it is an
    extrapolation of a proportionality that exactly two points cannot distinguish from any
    other monotone function through the same pair. Extrapolating it would put a number with
    no support behind it into a product forecast, which is the failure this refusal exists
    to prevent -- and the fix is a third dilution rate, not a wider bound.
    """
    lo, hi = calibration.growth_rate_range
    if not lo <= growth_rate <= hi:
        raise ValueError(
            f"growth rate {growth_rate} /h is outside the calibrated range "
            f"[{lo}, {hi}] /h. Exactly two dilution rates were measured "
            f"({calibration.source}), so the growth-rate dependence of the cyclase capacity "
            "is an interpolation between those two and has no support outside them. Measure "
            "a third dilution rate, or take the extrapolation deliberately by widening "
            "growth_rate_range on a copy of the calibration and saying so at the call site."
        )


def calibrated_cyclase(
    growth_rate: float,
    calibration: CyclaseCalibration = ELIZONDO2025,
) -> tuple[float, float]:
    """``(vmax_lcy, km_lcy)`` for the cyclase at this growth rate, mmol/gDCW/h and mmol/gDCW.

    The kinetic reading of :data:`ELIZONDO2025`. A growth-rate-independent content relation
    is the same statement as ``vmax_lcy = capacity * growth_rate``, and that reading is the
    part of the result that should be treated with suspicion: a Michaelis-Menten cyclase
    with a fixed enzyme complement has a vmax that does *not* scale with growth rate, and
    Elizondo's own RT-qPCR has CrtYB mRNA falling by a third from the low dilution rate to
    the high one, which is the wrong direction. The fit is good and its mechanism is not
    established. ``docs/research/KINETIC_FIT.md`` §6 lays out the two candidates.

    Args:
        growth_rate: Specific growth rate, 1/h. Must lie in the calibrated range.
        calibration: Which fit to use.

    Raises:
        ValueError: outside the calibrated growth-rate range, or on a non-positive rate.
    """
    if growth_rate <= 0:
        raise ValueError(f"growth rate must be positive, got {growth_rate}")
    _check_growth_rate(growth_rate, calibration)
    return growth_rate * calibration.capacity_mmol_per_gdcw, calibration.km_mmol_per_gdcw


def calibrated_kinetics(growth_rate: float) -> CarotenoidKinetics:
    """Refuse to build a full parameter set, and name what is missing.

    Present because the obvious next call after :func:`calibrated_cyclase` is "give me the
    whole thing", and the whole thing does not exist. Elizondo measured lycopene and
    beta-carotene; phytoene and GGPP were not measured, so the synthase and the desaturase
    have a known flux and an unknown substrate concentration, and that pins no (vmax, km)
    pair -- infinitely many reproduce the data exactly. Returning plausible numbers here
    would put four unfitted constants inside a forecast that looks fitted.

    Raises:
        NotImplementedError: always, naming the steps and the measurement that would fix it.
    """
    raise NotImplementedError(
        "the calibration set does not identify the full branch. "
        f"Steps {', '.join(UNIDENTIFIABLE_STEPS)} have no fitted parameters because "
        "PMID 40891387 measured lycopene and beta-carotene but not phytoene or GGPP, so "
        "their substrate concentrations are unobserved in all six steady states and any "
        "(vmax, km) pair reproduces the data. Use calibrated_cyclase(growth_rate) for the "
        "step that is identified, and supply the upstream pair explicitly -- marking it "
        "asserted -- if you need the full CarotenoidKinetics. One phytoene HPLC peak, "
        "which the same injection already resolves at 285 nm, would close this."
    )


def calibration_states(path: pathlib.Path | None = None):
    """The six steady states the cyclase was fitted to, as a dataframe.

    Adds the four quantities the fit works in and the raw file does not carry:
    ``lycopene_content`` and ``betacarotene_content`` (mmol/gDCW, each ``q / mu``),
    ``desaturase_flux`` (their two rates summed, which is what crtI carries) and
    ``cyclase_fraction`` (the share of desaturase flux that reaches the product).

    Args:
        path: Override for the vendored TSV.

    Returns:
        One row per (strain, dilution rate).

    Raises:
        FileNotFoundError: naming the file, rather than returning an empty frame.
    """
    import pandas as pd

    tsv = path or paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"
    if not tsv.exists():
        raise FileNotFoundError(
            f"carotenoid calibration states not vendored at {tsv}. "
            "See data/carotenoid/SOURCE.md for what the file is and where it came from."
        )
    frame = pd.read_csv(tsv, sep="\t")
    frame["lycopene_content"] = frame.q_lycopene / frame.mu_per_h
    frame["betacarotene_content"] = frame.q_betacarotene / frame.mu_per_h
    frame["desaturase_flux"] = frame.q_lycopene + frame.q_betacarotene
    frame["cyclase_fraction"] = frame.q_betacarotene / frame.desaturase_flux
    return frame
