"""Thermodynamic flux analysis, because fixing concentrations made the model infeasible.

Computing dG' from one assumed concentration for every metabolite and freezing directions
shuts 221 reactions off entirely and yeast-GEM stops solving. That is a made-up metabolome
being refused, not thermodynamics refusing the model: no single value is right for 1573
metabolites at once.

Henry 2007 (PMID 17172310) makes ``ln(concentration)`` a variable inside a physiological
window and asks whether some assignment admits a feasible flux distribution. Directionality
couples to it through one binary per reaction, which is what makes the problem
mixed-integer::

    dG'_i  =  dG'0_i + RT * sum_j S_ji * x_j        x_j = ln c_j
    v_i    <=  M * z_i                              forward flux only when z_i = 1
    v_i    >= -M * (1 - z_i)
    v_fwd  <=  M * zf                                forward flux needs zf = 1
    v_rev  <=  M * zr                                reverse flux needs zr = 1
    zf + zr <= 1                                     one direction at a time
    dG'_i  <=  M * (1 - zf) - eps                    running forward forces dG' negative
    dG'_i  >= -M * (1 - zr) + eps                    running reverse forces dG' positive

One binary per direction rather than one per reaction, so a reaction carrying no flux has
zf = zr = 0 and its dG' is unconstrained. A single binary forces every reaction to pick a
signed dG' whether or not it is running, which over-constrains the concentrations badly
enough to need a window spanning thirteen orders of magnitude before the model grows at all.

A measurement enters by tightening the window on the cofactors a sensor reads. Whether that
changes any flux is then answerable, which is the only reason to build any of it.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

from ..fba.solver import growth_or_none

from .thermodynamic import (
    GAS_CONSTANT_KJ,
    STANDARD_TEMPERATURE_K,
    ThermodynamicData,
    reaction_energy,
)

__all__ = [
    "NARROWING_TOLERANCE",
    "CONCENTRATION_CEILING_M",
    "CONCENTRATION_FLOOR_M",
    "ThermoModel",
]

CONCENTRATION_FLOOR_M = 1e-5
CONCENTRATION_CEILING_M = 2e-2
"""The window Henry 2007 uses for intracellular species. Wide on purpose: it asserts only
that a metabolite is present at some plausible level, which is the whole difference between
this and the fixed-concentration version that failed.

**This floor was briefly widened to 1e-8 and the change was wrong.** The reasoning was that
measured cytosolic acetyl-CoA in glucose-limited CEN.PK113-7D is about 10 uM -- right on the
floor -- so TMFA "could not represent the glucose state". It can: `_apply` replaces the
window on any metabolite a caller measures, says so in its own docstring, and cites cytosolic
NADH at 3-30 uM as exactly this case. The floor is a prior for metabolites nobody measured,
not a claim about the ones they did.

Widening it also had a cost that the reasoning had not weighed: a lower floor lets more
reactions run in either direction, so TMFA narrows less, and two calibrated narrowing tests
began to fail. Left at Henry's value.
"""

_BIG_M = 1e6

NARROWING_TOLERANCE = 1e-3
"""Relative width reduction below which a narrowing is not worth reporting.

Relative rather than absolute because flux widths here span zero to thousands, and a fixed
1e-9 means one thing on an exchange and another on a small internal flux. It changes none of
the results below -- the narrowings that occur are 15 to 20 percent, well clear of it."""
_EPSILON = 1e-3
_FLUX_M = 1000.0


@dataclass
class ThermoModel:
    """A metabolic model carrying thermodynamic directionality constraints."""

    model: object = field(repr=False)
    thermo: ThermodynamicData = field(repr=False)
    constrained_reactions: list[str] = field(repr=False)
    concentration_variables: dict[str, object] = field(repr=False)
    scope: list[str] = field(repr=False)
    _unconstrained_growth: float = 0.0

    @property
    def n_concentration_variables(self) -> int:
        return len(self.concentration_variables)

    @property
    def n_direction_variables(self) -> int:
        return len(self.constrained_reactions)

    @property
    def n_constrained_reactions(self) -> int:
        return len(self.constrained_reactions)

    @property
    def concentration_bounds(self) -> tuple[float, float]:
        return (float(np.log(CONCENTRATION_FLOOR_M)), float(np.log(CONCENTRATION_CEILING_M)))

    @property
    def coverage(self) -> float:
        return len(self.constrained_reactions) / max(len(self.scope), 1)

    @classmethod
    def build(cls, thermo: ThermodynamicData, subsystems=None, reactions=None, model=None,
              temperature_k: float = STANDARD_TEMPERATURE_K) -> ThermoModel:
        """Attach concentration variables and directionality coupling to a model.

        Args:
            thermo: Transformed formation energies.
            subsystems: Subsystems to constrain; None constrains everything covered. Scoping
                matters because the mixed-integer problem grows with the reaction count.
            reactions: Explicit reaction ids to constrain, taking precedence over
                subsystems. Useful for the reactions a measured cofactor could actually
                decide, which do not sit tidily inside one subsystem.
            model: A cobra model; defaults to the project's yeast GEM.
            temperature_k: Temperature the energies apply at.
        """
        from .thermodynamic import _default_model

        model = (model or _default_model()).copy()
        if reactions is not None:
            wanted = set(reactions)
            scope = [r for r in model.reactions if r.id in wanted]
        else:
            scope = [r for r in model.reactions
                     if subsystems is None or (r.subsystem or "") in set(subsystems)]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # NOT `slim_optimize() or 0.0`: an infeasible LP returns nan, nan is
            # truthy, and `nan or 0.0` is nan -- a guard that reads as one and is not.
            unconstrained = growth_or_none(model) or 0.0

        rt = GAS_CONSTANT_KJ * temperature_k
        concentrations: dict[str, object] = {}
        constrained: list[str] = []
        low, high = float(np.log(CONCENTRATION_FLOOR_M)), float(np.log(CONCENTRATION_CEILING_M))

        for reaction in scope:
            stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
            dg0 = (thermo.reaction_energy(reaction) if hasattr(thermo, "reaction_energy")
                   else reaction_energy(reaction, thermo))
            if dg0 is None:
                continue
            for metabolite in stoichiometry:
                if metabolite not in concentrations:
                    variable = model.problem.Variable(f"lnC_{metabolite}", lb=low, ub=high)
                    concentrations[metabolite] = variable
                    model.add_cons_vars([variable])

            # One binary per direction, so a reaction carrying no flux leaves dG' free.
            runs_forward = model.problem.Variable(f"zf_{reaction.id}", type="binary")
            runs_reverse = model.problem.Variable(f"zr_{reaction.id}", type="binary")
            energy = model.problem.Variable(f"dG_{reaction.id}", lb=-_BIG_M, ub=_BIG_M)
            model.add_cons_vars([runs_forward, runs_reverse, energy])

            definition = model.problem.Constraint(
                energy - rt * sum(c * concentrations[m] for m, c in stoichiometry.items()),
                lb=dg0, ub=dg0, name=f"dGdef_{reaction.id}")
            model.add_cons_vars([
                definition,
                model.problem.Constraint(runs_forward + runs_reverse, ub=1.0,
                                         name=f"onedir_{reaction.id}"),
                model.problem.Constraint(reaction.forward_variable - _FLUX_M * runs_forward,
                                         ub=0.0, name=f"fwd_{reaction.id}"),
                model.problem.Constraint(reaction.reverse_variable - _FLUX_M * runs_reverse,
                                         ub=0.0, name=f"rev_{reaction.id}"),
                model.problem.Constraint(energy + _BIG_M * runs_forward,
                                         ub=_BIG_M - _EPSILON, name=f"sgnf_{reaction.id}"),
                model.problem.Constraint(energy - _BIG_M * runs_reverse,
                                         lb=_EPSILON - _BIG_M, name=f"sgnr_{reaction.id}"),
            ])
            constrained.append(reaction.id)

        return cls(model, thermo, constrained, concentrations,
                   [r.id for r in scope], unconstrained)

    def _apply(self, measured):
        """Replace the assumed window on measured metabolites with the measured one.

        A measurement may sit outside the generic window, and when it does the measurement
        wins. The window exists because concentrations are unknown; where one is known the
        assumption has no standing against it. This is not hypothetical -- free cytosolic
        NADH in yeast runs 3 to 30 uM (Canelas 2008, PMID 18383140) against a floor of 10,
        so a large part of the measured range lies under it.

        What is still refused is an ordering or a sign that cannot be true, because those
        are mistakes rather than surprises.
        """
        low, high = self.concentration_bounds
        for metabolite, (lo, hi) in (measured or {}).items():
            if metabolite not in self.concentration_variables:
                raise KeyError(f"{metabolite!r} is not a constrained metabolite in this scope")
            if not 0 < lo <= hi:
                raise ValueError(
                    f"{metabolite!r} measured at {lo}-{hi} M is not a positive interval")
            variable = self.concentration_variables[metabolite]
            # Widen before narrowing; setting lb first fails if it passes the current ub.
            variable.lb = -_BIG_M
            variable.ub = float(np.log(hi))
            variable.lb = float(np.log(lo))
        return low, high

    def _restore(self):
        low, high = self.concentration_bounds
        for variable in self.concentration_variables.values():
            variable.lb, variable.ub = low, high

    def is_feasible(self) -> bool:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return self.model.slim_optimize() is not None and not np.isnan(
                self.model.slim_optimize())

    def growth_rate(self) -> float:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            value = self.model.slim_optimize()
        return float(value) if value is not None and not np.isnan(value) else 0.0

    def unconstrained_growth_rate(self) -> float:
        return self._unconstrained_growth

    def flux_range(self, reaction_id: str, measured=None) -> tuple[float, float]:
        """Feasible range of one flux, optionally with measured concentrations applied.

        The objective is put back afterwards. Leaving it pointing at the reaction under
        examination makes every later growth query optimise the wrong thing, which is a
        silent wrong answer rather than an error.
        """
        self._apply(measured)
        original = self.model.objective
        try:
            reaction = self.model.reactions.get_by_id(reaction_id)
            flux = reaction.forward_variable - reaction.reverse_variable
            bounds = []
            for direction in ("min", "max"):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self.model.objective = self.model.problem.Objective(flux, direction=direction)
                    value = self.model.slim_optimize()
                bounds.append(float(value) if value is not None and not np.isnan(value) else 0.0)
            return (bounds[0], bounds[1])
        finally:
            self.model.objective = original
            if measured:
                self._restore()

    def measurement_value(self, measured: dict, reactions) -> dict[str, dict]:
        """What tightening the concentration window does to each flux range.

        The question the coupling exists to answer. A measurement that leaves every range
        where it found it has told the model nothing, however well grounded its energies.

        Read it against a control, because a tight box narrows things on its own. Pinning
        two randomly chosen metabolites as tightly as the real measurement narrows one flux
        in eight on average and up to three; the measured ATP and ADP narrow none. The
        number alone does not distinguish a measurement from an arbitrary constraint.
        """
        report = {}
        for reaction_id in reactions:
            free = self.flux_range(reaction_id)
            pinned = self.flux_range(reaction_id, measured=measured)
            width_free = free[1] - free[0]
            width_pinned = pinned[1] - pinned[0]
            kept = float(width_pinned / width_free) if width_free > 0 else 1.0
            report[reaction_id] = {
                "free": free, "measured": pinned,
                "narrowed": bool(width_free > 0 and kept < 1.0 - NARROWING_TOLERANCE),
                "fraction_kept": kept,
            }
        return report
