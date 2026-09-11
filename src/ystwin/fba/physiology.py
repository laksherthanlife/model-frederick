"""Pin the GSMM to a measured phenotype before it is asked for bounds.

With a free oxygen exchange Yeast9 respires glucose completely: no ethanol, and roughly twice the
growth a Crabtree-positive yeast reaches. Every downstream bound is only meaningful inside a
regime checked against a real phenotype.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from numbers import Real

import cobra

__all__ = [
    "REFERENCE_AEROBIC_BATCH",
    "REFERENCE_GLUCOSE_LIMITED_CHEMOSTAT",
    "PhysiologyReference",
    "PhysiologyReport",
    "aerobic_batch_constraints",
    "cap_uptake",
    "chemostat_uptake_constraints",
    "validate_physiology",
]

GLUCOSE_EXCHANGE = "r_1714"
OXYGEN_EXCHANGE = "r_1992"
ETHANOL_EXCHANGE = "r_1761"
CO2_EXCHANGE = "r_1672"
BIOMASS_REACTION = "r_2111"


def _exchange_boundaries(model: cobra.Model, exchange_id: str) -> tuple[cobra.Reaction, ...]:
    reaction = model.reactions.get_by_id(exchange_id)
    if len(reaction.metabolites) != 1:
        raise ValueError(f"{exchange_id}: uptake must be a one-metabolite boundary reaction")
    metabolite, coefficient = next(iter(reaction.metabolites.items()))
    if not math.isfinite(coefficient) or coefficient == 0:
        raise ValueError(f"{exchange_id}: uptake stoichiometry must be finite and nonzero")
    reverse_id = f"{exchange_id}_REV"
    if reverse_id in model.reactions:
        reverse = model.reactions.get_by_id(reverse_id)
        if len(reverse.metabolites) == 1 and metabolite in reverse.metabolites:
            reverse_coefficient = reverse.metabolites[metabolite]
            if not math.isfinite(reverse_coefficient) or reverse_coefficient == 0:
                raise ValueError(f"{reverse_id}: uptake stoichiometry must be finite and nonzero")
            if (coefficient > 0) != (reverse_coefficient > 0):
                return reaction, reverse
    return (reaction,)


def _uptake_boundary(model: cobra.Model, exchange_id: str) -> cobra.Reaction:
    boundaries = _exchange_boundaries(model, exchange_id)
    reaction = boundaries[0]
    if len(boundaries) == 2:
        reverse = boundaries[1]
        coefficient = next(iter(reaction.metabolites.values()))
        reverse_coefficient = next(iter(reverse.metabolites.values()))
        supplies = reaction.upper_bound > 0 if coefficient > 0 else reaction.lower_bound < 0
        reverse_supplies = (reverse.upper_bound > 0 if reverse_coefficient > 0
                            else reverse.lower_bound < 0)
        if supplies and not reverse_supplies:
            return reaction
        if ((reverse_coefficient > 0 and reverse.lower_bound >= 0)
                or (reverse_coefficient < 0 and reverse.upper_bound <= 0)):
            return reverse
    return reaction


def cap_uptake(model: cobra.Model, exchange_id: str, magnitude: float) -> str:
    """Cap uptake through ``exchange_id``, whichever direction this model carries it in.

    An enzyme-constrained model is irreversibly split: ``r_1714`` becomes an export fixed
    at ``(0, 0)`` and a separate ``r_1714_REV`` carries the supply. Setting
    ``lower_bound = -magnitude`` on the export therefore caps nothing and opens a route
    that *creates* extracellular substrate -- the opposite of the intent. Measured on
    ecYeastGEM_batch: growth was identical before and after, and glucose flowed uncapped
    at 17.93 through ``r_1714_REV`` against a requested cap of 21.3.

    Reads the split from the model rather than from its filename, and returns the
    reaction it actually constrained so a caller can assert on it.
    """
    if isinstance(magnitude, bool) or not isinstance(magnitude, Real):
        raise ValueError("uptake magnitude must be a finite real number")
    try:
        amount = abs(float(magnitude))
    except OverflowError:
        raise ValueError("uptake magnitude must be finite") from None
    if not math.isfinite(amount):
        raise ValueError("uptake magnitude must be finite")
    reaction = _uptake_boundary(model, exchange_id)
    coefficient = next(iter(reaction.metabolites.values()))
    reaction_flux = amount / abs(coefficient)
    if not math.isfinite(reaction_flux):
        raise ValueError("uptake magnitude must convert to a finite reaction flux")
    changes = []
    for boundary in _exchange_boundaries(model, exchange_id):
        coefficient = next(iter(boundary.metabolites.values()))
        lower, upper = boundary.bounds
        if boundary is reaction:
            if coefficient < 0:
                lower = -reaction_flux
            else:
                upper = reaction_flux
        elif coefficient < 0:
            lower = max(0.0, lower)
        else:
            upper = min(0.0, upper)
        if lower > upper:
            raise ValueError(f"{boundary.id}: uptake magnitude conflicts with a required flux")
        changes.append((boundary, (lower, upper)))
    for boundary, bounds in changes:
        boundary.bounds = bounds
    return reaction.id


@dataclass(frozen=True)
class PhysiologyReference:
    """A measured whole-cell phenotype used to check the model."""

    name: str
    growth_rate: float
    glucose_uptake: float
    oxygen_uptake: float
    ethanol_secretion: float
    co2_secretion: float
    source: str


# van Hoek's own highest dilution rate, read from Table 1 of the cited paper.
#
# This replaces a set of numbers attributed to the same paper that are not in it. The
# strings 21.3, 27.4 and 20.4 occur ZERO times in PMID 9797269, which ran no batch culture
# at all -- every steady state in it is a glucose-limited chemostat -- and whose strain is
# DS28911, not the CEN.PK113-7D the tests claimed (that string is also absent). Checked
# directly against the full text. The other van Hoek 1998 paper, PMID 9603825, does use
# CEN.PK113-7D but reports no extracellular fluxes and does not contain those numbers
# either; they are the only two van Hoek 1998 records in PubMed.
#
# The old values were close to twice these at the same growth rate, which reads like a
# genuine glucose-excess batch phenotype from a source nobody has found rather than
# transcription drift. So they are not preserved: a number we cannot source is not a
# measurement, and this repository has no other way to say so.
#
# CONSEQUENCE, stated because it moves a published figure. With the old constants plain
# Yeast9 grew at 0.6925 /h against a measured 0.40 and the docs reported "overpredicts by
# 73%". With these it grows at 0.3460, a 13.5% UNDERprediction. The direction of the
# error reverses. The ec model is still the better choice on this phenotype, but not for
# the reason that was written down.
#
# LABEL, corrected 2026-09-09. The old identifier says "aerobic batch"; this row is a steady
# state of an aerobic glucose-LIMITED CHEMOSTAT, and PMID 9797269 ran no batch culture.
#
# yeast-GEM cannot reach this row's growth under this row's own uptakes -- 0.40 measured
# against a model maximum of 0.34598, efficiency 1.156. See docs/FINDINGS.md 2026-09-09.
REFERENCE_GLUCOSE_LIMITED_CHEMOSTAT = PhysiologyReference(
    name="glucose-limited chemostat, D = 0.40 /h",
    growth_rate=0.40,
    glucose_uptake=11.1,
    oxygen_uptake=3.7,
    ethanol_secretion=13.9,
    co2_secretion=18.9,
    source="van Hoek, van Dijken & Pronk 1998, PMID 9797269, Table 1 at D = 0.40 /h; "
           "S. cerevisiae DS28911, aerobic glucose-limited chemostat",
)

# Compatibility alias only. Twenty-odd modules, scripts and tests import the old misnomer
# and they are not this lane's files, so the wrong word survives as a name and nowhere else.
REFERENCE_AEROBIC_BATCH = REFERENCE_GLUCOSE_LIMITED_CHEMOSTAT


@dataclass(frozen=True)
class PhysiologyReport:
    """Model fluxes beside the measured phenotype they are supposed to reproduce."""

    reference: PhysiologyReference
    observed: dict[str, float]
    expected: dict[str, float]
    tolerance: float
    relative_error: dict[str, float] = field(default_factory=dict)

    @property
    def max_relative_error(self) -> float:
        return max(self.relative_error.values()) if self.relative_error else float("inf")

    @property
    def passed(self) -> bool:
        return self.max_relative_error <= self.tolerance

    def summary(self) -> str:
        lines = [f"{self.reference.name} ({'PASS' if self.passed else 'FAIL'}), "
                 f"tolerance {self.tolerance:.0%}"]
        for key in self.expected:
            lines.append(
                f"  {key:<14} model {self.observed[key]:>8.3f}  "
                f"measured {self.expected[key]:>8.3f}  "
                f"rel.err {self.relative_error[key]:>6.1%}"
            )
        lines.append(f"  source: {self.reference.source}")
        return "\n".join(lines)


def chemostat_uptake_constraints(
    model: cobra.Model,
    reference: PhysiologyReference = REFERENCE_GLUCOSE_LIMITED_CHEMOSTAT,
    glucose_uptake: float | None = None,
    oxygen_uptake: float | None = None,
) -> cobra.Model:
    """Apply the measured uptake regime to ``model`` in place.

    Args:
        model: Model to constrain. Use inside a ``with model:`` block to scope it.
        reference: Phenotype supplying the default uptake rates.
        glucose_uptake: Override, mmol/gDW/h (magnitude).
        oxygen_uptake: Override, mmol/gDW/h (magnitude). This is the constraint
            that produces overflow metabolism; leaving it free does not.
    """
    glc = reference.glucose_uptake if glucose_uptake is None else glucose_uptake
    o2 = reference.oxygen_uptake if oxygen_uptake is None else oxygen_uptake
    cap_uptake(model, GLUCOSE_EXCHANGE, glc)
    cap_uptake(model, OXYGEN_EXCHANGE, o2)
    return model


# Compatibility alias, kept for the importers outside this lane; the regime it applies is a
# chemostat steady state, not a batch culture.
aerobic_batch_constraints = chemostat_uptake_constraints


def validate_physiology(
    model: cobra.Model,
    reference: PhysiologyReference = REFERENCE_GLUCOSE_LIMITED_CHEMOSTAT,
    apply_constraints: bool = True,
    tolerance: float = 0.35,
    keys: tuple[str, ...] | None = None,
) -> PhysiologyReport:
    """Compare model fluxes with a measured phenotype. Does not modify ``model``.

    Args:
        model: Model to check.
        reference: Measured phenotype to compare against.
        apply_constraints: Apply the reference uptake regime first. Leave off for
            an enzyme-constrained model, which sets its own uptake from its
            protein pool.
        tolerance: Largest relative error still counted as a pass.
        keys: Subset of reference fluxes to score. Defaults to all of them. Use a
            subset when a model is known to reproduce one regime and not another,
            so that a pass is a claim about the part actually being relied on.
    """
    # ec models split reversible exchanges, so uptake is a positive _REV flux.
    with model as m:
        if apply_constraints:
            chemostat_uptake_constraints(m, reference)
        m.objective = BIOMASS_REACTION
        sol = m.optimize(raise_error=True)

        def boundary_rate(rid: str) -> float:
            return float(sum(
                next(iter(boundary.metabolites.values())) * sol.fluxes[boundary.id]
                for boundary in _exchange_boundaries(m, rid)
            ))

        observed = {
            "growth_rate": float(sol.objective_value),
            "glucose": boundary_rate(GLUCOSE_EXCHANGE),
            "oxygen": boundary_rate(OXYGEN_EXCHANGE),
            "ethanol": -boundary_rate(ETHANOL_EXCHANGE),
            "co2": -boundary_rate(CO2_EXCHANGE),
        }
    all_expected = {
        "growth_rate": reference.growth_rate,
        "glucose": reference.glucose_uptake,
        "oxygen": reference.oxygen_uptake,
        "ethanol": reference.ethanol_secretion,
        "co2": reference.co2_secretion,
    }
    if keys is not None:
        unknown = set(keys) - set(all_expected)
        if unknown:
            raise KeyError(
                f"no reference value for {sorted(unknown)}; "
                f"have {sorted(all_expected)}"
            )
        expected = {k: all_expected[k] for k in keys}
    else:
        expected = all_expected
    rel = {
        k: abs(observed[k] - expected[k]) / abs(expected[k]) if expected[k] else float("inf")
        for k in expected
    }
    scored = {k: observed[k] for k in expected}
    return PhysiologyReport(
        reference=reference, observed=scored, expected=expected,
        tolerance=tolerance, relative_error=rel,
    )
