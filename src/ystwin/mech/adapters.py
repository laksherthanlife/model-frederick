from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import TYPE_CHECKING

from .contracts import (
    Genotype, ObservationModel, PhysicalState, Protocol, SimulationResult,
)

if TYPE_CHECKING:
    from .engine import EngineParameters


__all__ = ["reaction_table", "simulate_protocol"]


def simulate_protocol(
    protocol: Protocol,
    genotype: Genotype,
    initial_state: PhysicalState,
    parameters: EngineParameters,
    *,
    observation: ObservationModel | None = None,
    **solver_options,
) -> SimulationResult:
    """Forward complete physical inputs unchanged to the single shared engine.

    No legacy or synthetic latent-to-state conversion is defined: it would require
    calibrated molecular inventories, gene-specific expression/activities and a
    compatible observation model, not merely similarly named stress coordinates.
    The engine's unit-bearing variables, validity and scientific refusals pass through.
    """
    from . import engine

    for name, value, expected in (
        ("protocol", protocol, Protocol),
        ("genotype", genotype, Genotype),
        ("initial_state", initial_state, PhysicalState),
        ("parameters", parameters, engine.EngineParameters),
    ):
        if not isinstance(value, expected):
            raise ValueError(f"shared {name} requires {expected.__name__}; legacy inputs are not converted or completed with priors")
    if observation is not None and not isinstance(observation, ObservationModel):
        raise ValueError("shared observation requires ObservationModel or None")
    return engine.simulate(protocol, genotype, initial_state, parameters,
                           observation=observation, **solver_options)


def reaction_table(parameters: EngineParameters,
                   genotype: Genotype) -> Mapping[str, Mapping[str, float]]:
    """The engine's own stoichiometry for this genotype: reaction -> coordinate -> coefficient.

    Which reactions make and which consume a pool is a STRUCTURAL fact the engine already
    states, so a caller reading it here cannot drift from what the integrator did. The
    alternative -- a product-to-reaction name map written down somewhere else -- would go
    stale the first time a cassette or a reaction moved, and nothing would notice.

    Death transfer, growth dilution, feed and gas transfer are applied outside this matrix
    and are absent from it, so a rate computed from these coefficients is reaction extent
    and nothing else. The keys index ``SimulationResult.truth`` as ``flux.<reaction>``, in
    mmol/h, so an extent at a read time is one lookup away.

    It lives here rather than in a caller because this module is the shared engine's seam:
    the kernel is private, and one place reaching into it is better than several.
    """
    from . import engine

    if not isinstance(genotype, Genotype) or not isinstance(parameters, engine.EngineParameters):
        raise ValueError("shared reaction table requires Genotype and EngineParameters; legacy inputs are not converted or completed with priors")
    kernel = engine._Kernel(parameters, genotype)
    return MappingProxyType({name: MappingProxyType({k: float(v) for k, v in reaction.items()})
                             for name, reaction in kernel.stoichiometry.items()})
