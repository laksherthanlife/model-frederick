"""Choosing which reporters to build, given a hard cap on spectral channels.

Five fluorophores is five channels, and one is spent on the constitutive reference that
identifies growth, so four stress reporters must cover a two-dozen-module landscape whose
modules are not orthogonal. Selection is D-optimal on Fisher information, which makes the
noise of the actual reader part of the question rather than an afterthought.

A build usually has a purpose beyond spanning volume. That is expressed as a priority
weight on the modules that matter, so the objective is asked for what the build is for
instead of being overridden once it answers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

import numpy as np

from ..generator.panel_experiment import MEASURED_GROWTH_RATE_SE, OBSERVED_ACTIVITY_CV
from ..generator.stress_panel import REPORTERS, Kind, reporter_loadings
from .design import NoiseModel, design_value, identifiable_by_precision

__all__ = [
    "RECOMMENDED_BUILD",
    "THREE_SENSOR_BUILD",
    "REFERENCE_GROWTH_RATE",
    "SensorBuild",
    "growth_penalty_cv",
    "interpretable",
    "uncharacterised_ph",
    "select_sensors",
    "spectral_conflict",
]

REFERENCE_GROWTH_RATE = 0.22
"""Growth rate a build is chosen for, in 1/h.

A stressed culture, not an unstressed one: UPRE1 runs at 0.22 under 1 mM DTT, and the point
of the panel is to read cultures that are being stressed. The penalty below is relative, so
it needs a rate to be relative to, and picking the unstressed one would understate it.
"""


def growth_penalty_cv(rate_se: float | None = None, growth_rate: float | None = None) -> float:
    """Extra relative error a promoter fusion carries from dividing out its growth rate.

    Computed rather than written down. The measured quantity is an absolute rate error and
    the penalty is relative, so the two are related by the growth rate of the culture --
    and freezing the quotient as a literal leaves a constant that agrees with its inputs
    only until one of them is remeasured. A ratiometric sensor pays none of it.
    """
    rate_se = MEASURED_GROWTH_RATE_SE if rate_se is None else rate_se
    growth_rate = REFERENCE_GROWTH_RATE if growth_rate is None else growth_rate
    if not np.isfinite(growth_rate) or growth_rate <= 0:
        raise ValueError(f"growth rate must be positive and finite, got {growth_rate}")
    if not np.isfinite(rate_se) or rate_se < 0:
        raise ValueError("growth-rate standard error must be finite and non-negative")
    return float(rate_se) / float(growth_rate)


def spectral_conflict(names) -> bool:
    """Whether a build asks two sensors to share one exclusive channel.

    Only biochemical sensors claim a slot. They are single fixed proteins, and the green
    excitation-ratio family -- roGFP2, HyPer, QUEEN, pHluorin, Peredox -- all emit around
    510 nm, so two of them in one strain give one measurement, not two.
    """
    claimed = [REPORTERS[n].spectral_slot for n in names if REPORTERS[n].spectral_slot]
    return len(claimed) != len(set(claimed))


def interpretable(names) -> bool:
    """Whether a build can attribute its own readings.

    A pH-sensitive sensor needs a pH reading beside it. One pH unit moves a glutathione
    probe by 60 mV where a genuine peroxide response is 40 to 50, and in the same direction,
    so without pH the two cannot be told apart.
    """
    names = list(names)
    if not any(REPORTERS[n].ph_sensitive for n in names):
        return True
    return any(REPORTERS[n].module == "ph" for n in names)


def uncharacterised_ph(names) -> list[str]:
    """Sensors in a build that were never shown flat over the pH the cytosol reaches.

    Separate from :func:`interpretable`, and deliberately so. A pH-coupled probe is known
    to move with pH and needs a pH channel beside it. This is the other case: nobody has
    published the measurement either way, so the build carries an untested region rather
    than a known confound, and saying so is not the same as calling the sensor bad.
    """
    return [n for n in names if REPORTERS[n].ph_uncharacterised]


def select_sensors(
    n_channels: int,
    n_reference: int = 0,
    library=None,
    required=(),
    noise: NoiseModel | None = None,
    priority: dict[str, float] | None = None,
    min_effect: float = 0.5,
    growth_cv: float = 0.0,
) -> list[str]:
    """Reporter set measuring the most modules within the channel budget.

    Args:
        n_channels: Spectral channels available.
        n_reference: Channels reserved for constitutive references, which carry no
            stress information but identify growth.
        library: Reporter names to choose from; defaults to the whole panel.
        required: Reporters that must appear, whatever the score says.
        noise: Measurement noise of the reader the build is for.
        priority: Weights raising the value of particular modules.
        min_effect: Change in module activity the build has to resolve.
        growth_cv: Relative noise a promoter fusion inherits from having growth divided
            out of it. Ratiometric sensors never pay it, which is what separates the two
            kinds once loading magnitude has been normalised away.
    """
    names = list(REPORTERS if library is None else library)
    if any(not isinstance(n, (int, np.integer)) or isinstance(n, bool) or n < 0
           for n in (n_channels, n_reference)):
        raise ValueError("channel and reference counts must be non-negative integers")
    if len(set(names)) != len(names):
        raise ValueError("a reporter cannot occupy multiple library entries")
    unknown = set(names) - set(REPORTERS)
    if unknown:
        raise KeyError(f"reporters not in the declared library: {sorted(unknown)}")
    if not np.isfinite([min_effect, growth_cv]).all() or min_effect <= 0 or growth_cv < 0:
        raise ValueError("min_effect must be positive and growth_cv non-negative, both finite")
    budget = n_channels - n_reference
    if budget < 1:
        raise ValueError(f"{n_channels} channels with {n_reference} reserved leaves none for sensors")
    if budget > len(names):
        raise ValueError(f"asked for {budget} sensors but the library has only {len(names)}")
    required = list(required)
    if len(set(required)) != len(required):
        raise ValueError("required reporters must be unique")
    if len(required) > budget:
        raise ValueError("required reporters exceed the available channel budget")
    missing = [r for r in required if r not in names]
    if missing:
        raise KeyError(f"required reporters not in the library: {missing}")

    if spectral_conflict(required):
        raise ValueError(f"required reporters {required} want the same spectral slot")

    noise = noise or NoiseModel()
    free = [n for n in names if n not in required]
    best, best_value = None, (-np.inf, -np.inf)
    for extra in combinations(free, budget - len(required)):
        candidate = required + list(extra)
        if spectral_conflict(candidate) or not interpretable(candidate):
            continue
        value = design_value(reporter_loadings(candidate), noise, min_effect, priority,
                             extra_cv=_growth_penalty(candidate, growth_cv))
        if value > best_value:
            best, best_value = candidate, value
    if best is None:
        raise ValueError("no spectrally compatible, interpretable reporter set fits the channel budget")
    return list(best)


def _growth_penalty(names, growth_cv: float) -> np.ndarray:
    """Growth-correction noise, charged only to the channels that need the correction."""
    return np.array([growth_cv if REPORTERS[n].kind is Kind.TRANSCRIPTIONAL else 0.0
                     for n in names])


@dataclass(frozen=True)
class SensorBuild:
    """A concrete strain design, with the assumptions it was chosen under."""

    stress_reporters: list[str]
    reference_fluorophore: str
    noise: NoiseModel
    min_effect: float
    priority: dict[str, float] = field(default_factory=dict)
    growth_cv: float = 0.0
    selection_basis: str = "prespecified reporter loadings, noise and priorities; not held-out performance evidence"

    @property
    def n_channels(self) -> int:
        return len(self.stress_reporters) + 1

    @property
    def identifiable(self) -> list[str]:
        return identifiable_by_precision(
            reporter_loadings(self.stress_reporters), self.noise, self.min_effect,
            extra_cv=_growth_penalty(self.stress_reporters, self.growth_cv))

    @property
    def buildable(self) -> bool:
        return not spectral_conflict(self.stress_reporters)


_PLATE_NOISE = NoiseModel(relative_cv=OBSERVED_ACTIVITY_CV)
"""The noise of the plates this build is actually for, not a plausible figure for a plate.

It was 0.05, against 0.146 measured on two real plates over 28 matched conditions. The
build survives unchanged, because standard errors scale linearly in the constant and the
selector compares sets under one noise level -- but the smallest effect the build can
resolve roughly triples, so the recommendation is the same and the claim attached to it is
not."""
_AXIS_PRIORITY = {"ESR": 20.0, "UPR": 10.0, "oxidative": 10.0, "atp": 10.0}
"""What the build is for: the three stress axes it is being made to report, and the general
response that separates a stressed culture from a merely late one.

Weighting ESR alone chose a peroxide channel and no ATP channel at all, which reported the
adenylate pool at 0.01. Naming the axes swaps HyPer7 for QUEEN and PACE for TRX2, lifts ATP
to 0.67 and held-out transfer from five modules to seven, and costs peroxide 0.85 to 0.50.
The choice is stable: halving the ESR weight or more than doubling the ATP weight returns
the same four reporters."""
_MIN_EFFECT = 0.5
_GROWTH_CV = growth_penalty_cv()


RECOMMENDED_BUILD = SensorBuild(
    stress_reporters=select_sensors(
        n_channels=5, n_reference=1, noise=_PLATE_NOISE, priority=_AXIS_PRIORITY,
        min_effect=_MIN_EFFECT, growth_cv=_GROWTH_CV),
    reference_fluorophore="mOrange2",
    noise=_PLATE_NOISE,
    min_effect=_MIN_EFFECT,
    priority=_AXIS_PRIORITY,
    growth_cv=_GROWTH_CV,
)


THREE_SENSOR_BUILD = ["UPRE-ER", "TRX2-oxidative", "QUEEN-2m"]
"""The smallest build that still reports the landscape: ER, oxidative and ATP.

Three channels, one per stress axis, all three demonstrated in S. cerevisiae. The oxidative
channel has to be the transcriptional one, because HyPer7 and QUEEN are both green
excitation-ratio and cannot share a plate.

It reads three modules and recovers nine, because crosstalk and the regulon cascade carry
information about modules no channel touches. That only holds if the analysis pools
replicates before fitting and takes growth rate from optical density as a fourth input --
without both, the same plate reports two modules instead of nine.

Its limit is the axes themselves. Held out entirely, a stressor routed through ER,
oxidative or ATP transfers at six to ten modules; one that is purely osmotic, metal or TOR
transfers at one. Three sensors buy generality across the stresses they are wired to, not
across all of them.
"""
