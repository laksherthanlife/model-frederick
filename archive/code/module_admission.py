"""Admission check for latent modules that want to touch an FBA constraint.

A module earns a constraint only if its independent anchor is a flux the model reproduces.
Otherwise validation compares a prediction against a number the model never got right. Encoded as
a guard because a deferral in prose is forgotten the next time somebody adds a module.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..fba.physiology import PhysiologyReport

__all__ = ["ModuleDecision", "ModuleRefused", "ModuleSpec", "admit_module"]


class ModuleRefused(Exception):
    """Raised when a refused module is used anyway."""


@dataclass(frozen=True)
class ModuleSpec:
    """A candidate latent module and what it would do to the solve.

    Args:
        name: Biological name, e.g. ``oxidative_burden`` -- not "sensor 3 state".
        anchor_fluxes: Reference fluxes the module's independent anchor is read in.
        constraint: The pre-solve constraint it would set, in words.
    """

    name: str
    anchor_fluxes: tuple[str, ...]
    constraint: str


@dataclass(frozen=True)
class ModuleDecision:
    """Whether a module may be wired to a constraint, and why."""

    module: ModuleSpec
    admitted: bool
    reason: str

    def require(self) -> ModuleSpec:
        """Return the module, or raise if it was refused."""
        if not self.admitted:
            raise ModuleRefused(f"{self.module.name}: {self.reason}")
        return self.module


def admit_module(
    module: ModuleSpec,
    physiology_reports: list[PhysiologyReport],
    max_relative_error: float = 0.35,
) -> ModuleDecision:
    """Decide whether a module's anchor fluxes are ones the model can reproduce.

    Args:
        module: Candidate module.
        physiology_reports: Validation reports for the model that would carry the
            constraint, covering the fluxes in question.
        max_relative_error: Largest relative error in an anchor flux that still
            permits admission. Tighten it for a confirmatory claim.
    """
    if not module.anchor_fluxes:
        return ModuleDecision(
            module, False,
            "no anchor flux declared; a latent state with nothing to validate "
            "against cannot be mapped to a constraint",
        )

    errors: dict[str, float] = {}
    for report in physiology_reports:
        errors.update(report.relative_error)

    missing = [f for f in module.anchor_fluxes if f not in errors]
    if missing:
        return ModuleDecision(
            module, False,
            f"the model reports no value for {', '.join(sorted(missing))}, so the "
            f"anchor cannot be checked against it",
        )

    failed = {f: errors[f] for f in module.anchor_fluxes if errors[f] > max_relative_error}
    if failed:
        detail = ", ".join(f"{f} off by {errors[f]:.0%}" for f in sorted(failed))
        return ModuleDecision(
            module, False,
            f"{detail} (limit {max_relative_error:.0%}); the model does not "
            f"reproduce this module's anchor, so validating against it proves nothing",
        )

    worst = max(errors[f] for f in module.anchor_fluxes)
    return ModuleDecision(
        module, True,
        f"anchor fluxes {', '.join(module.anchor_fluxes)} reproduced to "
        f"{worst:.0%} or better; constraint '{module.constraint}' may be wired",
    )
