"""Branch B of g_psi: latent stress state -> ATP maintenance constraint. Research only.

G4 returned INCONCLUSIVE for all four constructs, and §9 keeps a latent-to-constraint mapping
only once it predicts an excluded anchor. So the branch label is always "latent", a refusal
string travels on every prediction, and calling it requires acknowledging that status.

The effect is a named pre-solve constraint, not a post-solve multiplier: a multiplier cannot be
infeasible, so it can never be wrong, so it never tests anything.
"""

from __future__ import annotations

from dataclasses import dataclass

import cobra
import pandas as pd

from .physiology_bridge import PhysiologicalConstraints, predict_flux

__all__ = ["LatentState", "LatentConstraints", "compare_branches", "latent_constraints"]

NGAM_REACTION = "r_4046"
_RESTING_MAINTENANCE = 0.7  # mmol ATP/gDW/h, the Yeast9 NGAM value

# The stress-driven maintenance demand at full scale, mmol ATP/gDW/h ON TOP of resting.
#
# This replaces `_MAINTENANCE_PER_ACTIVITY = 300.0`, which was not merely unfitted but
# dimensionally incoherent: it multiplied a promoter activity in RFU/OD/h and called the
# product mmol ATP/gDW/h. On the original plate snapshot the excess activity ran to ~1305
# RFU/OD/h, so it demanded an NGAM of 316,000-391,000 against a resting value of 0.7.
# Yeast9 goes infeasible at 19.19, and growth is already zero there, so the old constant
# sat about 2e4 past the point where the cell stops growing at all.
#
# The scale here is bounded by reproducing Lahtvee 2016's design in the model rather than
# read off their figure: glucose-limited chemostat at D = 0.1 /h, where a yield of ~0.5
# g/g puts q_glucose near 1.1-1.5 mmol/gDW/h, and the NGAM headroom above resting at that
# uptake is 6.5. Lahtvee is the right anchor because it holds the dilution rate fixed --
# "specific growth rate-dependent changes are eliminated" -- and concludes that raised
# maintenance ATP is what underpins the general stress response, which is the mechanism
# this branch already assumed. Tier 0 still: an envelope, not their measured value, which
# is published as a bar chart.
_MAX_STRESS_MAINTENANCE = 6.5

# Excess promoter activity treated as full scale, RFU/OD/h. The measured dynamic range on
# the four constructs, so the ratio below is dimensionless and lands in roughly [0, 1].
#
# DERIVED, NOT CHOSEN, and the derivation is pinned by a test rather than by this comment:
# the largest `activity_late` excess over that construct's own zero-dose median, taken
# across all four constructs in `outputs/sensor_characterisation.csv`. The verified
# current sensor-characterisation output gives 1457.2761833373581 RFU/OD/h.
#
# It was 1305.0 until 2026-09-04, then 1465.0 before the corrected-source replay. These
# are snapshots of a current-derived calibration, not a frozen historical model. `src/`
# must not import from `outputs/`, so the constant stays a literal;
# `tests/test_latent_maintenance.py::TestTheFullScaleIsStillWhatThePlatesSay` independently
# re-derives it from the artifact and fails on drift. Update it only from a verified
# `scripts/run_sensor_characterisation.py` output, using that same excess formula.
_ACTIVITY_FULL_SCALE = 1457.2761833373581

# The G4 outcome, recorded here so the refusal is stated wherever the branch is used.
_G4_STATUS = (
    "REFUSED by G4: all four constructs returned INCONCLUSIVE "
    "(Hac1 anchor fails its no-RT control in 3/3 replicates; TRX2 anchor "
    "underpowered at n=3, needing ~9 and ~80 replicates)"
)


@dataclass(frozen=True)
class LatentState:
    """An inferred stress state. Not a measurement.

    Args:
        promoter_activity: Dilution-corrected promoter activity from the estimator.
        construct: Which sensor produced it.
        time_h: Time of the estimate.
    """

    promoter_activity: float
    construct: str
    time_h: float


@dataclass(frozen=True)
class LatentConstraints:
    """Branch A's bounds plus a latent-driven maintenance term, labelled as such."""

    growth_lower_bound: float
    glucose_uptake: float
    oxygen_uptake: float | None
    atp_maintenance: float
    resting_maintenance: float
    time_h: float
    validation: str
    construct: str
    named_effect: str = "ATP maintenance"
    branch: str = "latent"

    def summary(self) -> str:
        return (
            f"[{self.branch}] t={self.time_h:.2f} h  {self.named_effect} "
            f"{self.atp_maintenance:.3f} mmol/gDW/h "
            f"(resting {self.resting_maintenance:.3f})  --  {self.validation}"
        )


def latent_constraints(
    latent: LatentState,
    physical: PhysiologicalConstraints,
    acknowledge_unvalidated: bool = False,
    basal_activity: float = 1.0e-3,
    max_stress_maintenance: float = _MAX_STRESS_MAINTENANCE,
    activity_full_scale: float = _ACTIVITY_FULL_SCALE,
    resting_maintenance: float = _RESTING_MAINTENANCE,
) -> LatentConstraints:
    """Add a latent-driven maintenance demand to Branch A's bounds.

    Args:
        latent: Inferred stress state.
        physical: Branch A constraints to build on.
        acknowledge_unvalidated: Must be ``True``. The gate exists so that using
            this branch is a deliberate act recorded at the call site.
        basal_activity: Promoter activity treated as unstressed.
        max_stress_maintenance: Additional maintenance at full scale, mmol ATP/gDW/h.
        activity_full_scale: Excess activity treated as full scale, RFU/OD/h. Divides the
            excess so the coefficient above is in ATP units rather than ATP-per-RFU.
        resting_maintenance: Non-growth maintenance with no stress.

    Raises:
        ValueError: unless ``acknowledge_unvalidated`` is set.
    """
    if not acknowledge_unvalidated:
        raise ValueError(
            "this branch is not validated -- " + _G4_STATUS + "; pass "
            "acknowledge_unvalidated=True to use it for exploration, and label any "
            "number it produces accordingly"
        )
    if activity_full_scale <= 0:
        raise ValueError("activity_full_scale must be positive; it is a divisor")
    excess = max(latent.promoter_activity - basal_activity, 0.0)
    # Normalised, then scaled into ATP units. Clipped at full scale because extrapolating
    # past the measured activity range is where the old constant did its damage.
    fraction = min(excess / activity_full_scale, 1.0)
    return LatentConstraints(
        growth_lower_bound=physical.growth_lower_bound,
        glucose_uptake=physical.glucose_uptake,
        oxygen_uptake=physical.oxygen_uptake,
        atp_maintenance=float(resting_maintenance + max_stress_maintenance * fraction),
        resting_maintenance=float(resting_maintenance),
        time_h=latent.time_h,
        validation=_G4_STATUS,
        construct=latent.construct,
    )


def compare_branches(
    model: cobra.Model,
    physical: PhysiologicalConstraints,
    latent: LatentConstraints,
    reactions: list[str],
    ngam_reaction: str = NGAM_REACTION,
) -> pd.DataFrame:
    """Solve both branches and return them side by side, attributed.

    One row per branch, with ``load_bearing`` marking which may be relied on. The
    whole point of running them together is to see what the latent layer would add
    without letting its numbers pass as validated.
    """
    rows = []

    physical_prediction = predict_flux(model, physical, reactions)
    rows.append({
        "branch": "physiology", "load_bearing": True,
        "validation": "measured inputs only; no latent module used",
        "atp_maintenance": _RESTING_MAINTENANCE,
        **{f"flux_{r}": v for r, v in physical_prediction.fluxes.items()},
    })

    with model as m:
        maintenance = m.reactions.get_by_id(ngam_reaction)
        maintenance.bounds = (latent.atp_maintenance, latent.atp_maintenance)
        as_physical = PhysiologicalConstraints(
            growth_lower_bound=latent.growth_lower_bound,
            glucose_uptake=latent.glucose_uptake,
            oxygen_uptake=latent.oxygen_uptake,
            time_h=latent.time_h,
            branch="latent",
        )
        latent_prediction = predict_flux(m, as_physical, reactions)
    rows.append({
        "branch": "latent", "load_bearing": False,
        "validation": latent.validation,
        "atp_maintenance": latent.atp_maintenance,
        **{f"flux_{r}": v for r, v in latent_prediction.fluxes.items()},
    })
    return pd.DataFrame(rows)
