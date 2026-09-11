"""Reaction energies from component contribution, replacing the group-contribution table.

The SEED table shipped with pytfa is group contribution, which fails where substrate and
product are structurally similar because the errors do not cancel. It put ribose-5-phosphate
isomerase at +75.4 kJ/mol with a claimed error of 1.7, made an essential reaction impossible
and blocked growth through the whole pentose phosphate pathway.

eQuilibrator is component contribution, from the group that curates BioNumbers. It handles
pH, ionic strength and magnesium natively rather than through hand-added association
constants, reports an uncertainty per reaction, and on ATP hydrolysis returns the textbook
-30.5 kJ/mol exactly.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from functools import lru_cache

__all__ = ["EquilibratorData"]

CYTOSOLIC_PH = 7.5
IONIC_STRENGTH = "0.25M"
TEMPERATURE = "303.15K"
P_MG = 3.0
"""pMg 3 is 1 mM free Mg2+, the concentration the textbook ATP hydrolysis value assumes."""


def _formula_of(compound):
    """The matched compound's formula, or ``None`` if it does not carry one."""
    try:
        return compound.formula
    except Exception:
        return None


def _same_heavy_atoms(model_formula: str | None, compound) -> bool:
    """Whether the matched compound is the same molecule the model means.

    Heavy atoms only: protonation state legitimately differs with pH, so UDP as ``C9H11...``
    against ``C9H12...`` is the same molecule and must not be refused. Everything else must
    agree exactly. A missing formula on either side is accepted, because "cannot check" and
    "checked and wrong" are different and only the second is a reason to drop a resolution.
    """
    import re

    their = _formula_of(compound)
    if not model_formula or not their:
        return True

    def heavy(formula: str) -> dict:
        # Accumulated, not assigned -- see `fba/secretion.py::_less_thiol_hydrogens`. Here
        # the input is third-party eQuilibrator data, so a non-Hill spelling is not under
        # this repository's control and a comprehension would make the guard pass a
        # mismatch it exists to catch.
        counts: dict[str, int] = {}
        for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
            if element:
                counts[element] = counts.get(element, 0) + (int(number) if number else 1)
        counts.pop("H", None)
        return counts

    return heavy(model_formula) == heavy(their)


@dataclass
class EquilibratorData:
    """Component-contribution energies for the metabolites of a yeast model."""

    contribution: object = field(repr=False)
    compounds: dict = field(repr=False)
    ph: float
    p_mg: float
    coverage: float
    ionic_strength: str = IONIC_STRENGTH
    temperature: str = TEMPERATURE
    mismatched: dict = field(default_factory=dict, repr=False)
    """Metabolite id -> ``(model formula, matched formula)`` for the resolutions REFUSED.

    A KEGG identifier resolving is not the same as the right molecule being found. These are
    the ones where it was not, and they are excluded from :attr:`compounds` entirely, so a
    reaction containing one is uncovered rather than mis-scored. See :func:`_same_heavy_atoms`.
    """

    @classmethod
    def load(cls, model=None, ph: float = CYTOSOLIC_PH, p_mg: float = P_MG,
             ionic_strength: str = IONIC_STRENGTH,
             temperature: str = TEMPERATURE) -> EquilibratorData:
        """Resolve a model's metabolites against eQuilibrator by KEGG identifier.

        Args:
            model: A cobra model; defaults to the project's yeast GEM.
            ph: Compartment pH the energies are reported at.
            p_mg: Negative log of free Mg2+; 3.0 is 1 mM.
            ionic_strength: Ionic strength with units.
            temperature: Temperature with units.
        """

        from .thermodynamic import _default_model

        model = model or _default_model()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cc = _contribution()
        compounds: dict = {}
        mismatched: dict = {}
        for metabolite in model.metabolites:
            identifier = metabolite.annotation.get("kegg.compound")
            identifier = identifier[0] if isinstance(identifier, list) else identifier
            if not identifier:
                continue
            try:
                compound = cc.get_compound(f"kegg:{identifier}")
            except Exception:
                continue
            if compound is None:
                continue
            # THE MATCH IS CHECKED AGAINST THE FORMULA, not trusted because the KEGG id
            # resolved. yeast-GEM writes a polymer as one representative repeat unit and
            # eQuilibrator holds a specific oligomer, so `s_0773` glycogen (C6H10O5 in the
            # model) resolved to a TETRASACCHARIDE, C24H42O21, and glycogen synthase came
            # back at -849.2 kJ/mol for a glycosyl transfer textbooks put near -13. Nothing
            # errored: a mass-unbalanced reaction was scored and a number returned.
            # 41 of 1,720 resolved metabolites mismatch this way, all polymers and
            # isoprenoids -- `s_0641` dodecaprenyl diphosphate, C60 in the model, matched a
            # C10. Hydrogen is ignored because protonation state legitimately differs at pH.
            if not _same_heavy_atoms(metabolite.formula, compound):
                mismatched[metabolite.id] = (metabolite.formula,
                                             _formula_of(compound) or "?")
                continue
            compounds[metabolite.id] = compound
        return cls(cc, compounds, ph, p_mg,
                   len(compounds) / max(len(model.metabolites), 1),
                   ionic_strength, temperature, mismatched)

    def _apply_conditions(self) -> None:
        """Set the conditions before each query; the underlying object is shared."""
        from equilibrator_api import Q_

        self.contribution.p_h = Q_(self.ph)
        self.contribution.ionic_strength = Q_(self.ionic_strength)
        self.contribution.temperature = Q_(self.temperature)
        self.contribution.p_mg = Q_(self.p_mg)

    def reaction_dg0(self, stoichiometry: dict[str, float], with_error: bool = False):
        """Standard reaction energy in kJ/mol, or None if any participant is uncovered.

        Returning None rather than a partial sum is deliberate: an energy computed from some
        of a reaction's metabolites belongs to a different reaction.
        """
        from equilibrator_api import Reaction

        # COEFFICIENTS ARE SUMMED PER COMPOUND, NOT ASSIGNED. eQuilibrator identifies a
        # molecule, not a molecule-in-a-compartment: ATP[c] and ATP[m] resolve to ONE
        # `Compound` object, and 427 of the objects behind this model's metabolites are
        # shared that way. A dict comprehension keyed by the compound therefore kept the
        # LAST coefficient and threw the other away, turning every transport into a
        # creation from nothing -- `r_1111`, a strict ADP/ATP antiport whose standard
        # chemistry energy is zero, came back at **-3,634.7 kJ/mol** and `gate_step`
        # reported it as `runs`.
        #
        # 663 of the 1,916 fully-covered reactions collapse this way, median error 375
        # kJ/mol. Nothing caught it because every test and every statistic quoted about this
        # backend -- including the 1,098-reaction comparison that justified promoting it to
        # primary -- filters to SINGLE-COMPARTMENT reactions, where no compound is shared.
        # Found by the 2026-08-31 code review.
        #
        # Summed, a pure antiport collapses to `{}` and gives dG0 = 0, which is right: the
        # compartment work is the separate `electrical_work` term and belongs there.
        terms: dict = {}
        try:
            for metabolite, coefficient in stoichiometry.items():
                compound = self.compounds[metabolite]
                terms[compound] = terms.get(compound, 0.0) + coefficient
        except KeyError:
            return (None, None) if with_error else None
        self._apply_conditions()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            estimate = self.contribution.standard_dg_prime(Reaction(terms))
        value = float(estimate.value.m_as("kJ/mol"))
        if not with_error:
            return value
        return value, float(estimate.error.m_as("kJ/mol"))


    def dgf(self, metabolite_id: str):
        """Refused. Component contribution has no per-compound formation energy to give.

        **This returned ``0.0`` for any covered metabolite until 2026-08-30**, described as
        "present for interface compatibility". It is not compatible: every caller that sums
        ``coefficient * dgf(m)`` -- which is what `bridge/thermodynamic.py::reaction_dg0`
        and therefore the whole thermodynamic gate does -- would get **exactly zero** for
        every reaction and read it as a real energy near equilibrium. A stub that answers a
        question it cannot answer, with a number rather than an error, is the failure mode
        this repository's audit doc names first.

        Use :meth:`reaction_dg0`. Component contribution estimates the energy of a
        *reaction* from the groups that change across it; the per-compound term is not
        identified on its own, which is precisely why it handles the reactions group
        contribution gets wrong.

        Raises:
            NotImplementedError: always, when the metabolite is covered. ``None`` where it
                is not, so a coverage probe still works.
        """
        if metabolite_id not in self.compounds:
            return None
        raise NotImplementedError(
            "component contribution has no per-compound formation energy: it estimates a "
            "REACTION from the groups that change across it. Call reaction_dg0 with the "
            f"whole stoichiometry. (Asked for {metabolite_id!r}.)")

    def covers(self, metabolite_id: str) -> bool:
        """Whether this metabolite resolves, without asking for an energy it cannot give."""
        return metabolite_id in self.compounds

    def refutation(self, metabolite_id: str):
        """The row of ``refuted_energies.tsv`` refusing this metabolite's energy, or None.

        Shared with `bridge/thermodynamic.py` rather than duplicated: the refutation is a
        statement about the COMPOUND -- that a C40 polyene is outside the group
        decomposition either method was fitted on -- and it applies to whichever estimator
        is asked. An estimator that quietly dropped the refusal when it became the primary
        backend would be the one direction a refusal must never fail in.
        """
        from .thermodynamic import _refuted_energies, declared_seed_ids

        seed = declared_seed_ids().get(metabolite_id)
        if seed is None:
            return None
        return _refuted_energies(None).get(seed)

    def reaction_energy(self, reaction, split: bool = False):
        """Reaction energy including the electrical work of moving charge across a membrane."""
        from .thermodynamic import electrical_work

        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
        chemistry = self.reaction_dg0(stoichiometry)
        if chemistry is None:
            return None

        compartments = {m.compartment for m in reaction.metabolites}
        moves: dict = {}
        for metabolite, coefficient in reaction.metabolites.items():
            charge = metabolite.charge or 0
            if charge:
                key = (metabolite.compartment, charge)
                moves[key] = moves.get(key, 0.0) + coefficient
        net = sum(c * q for (comp, q), c in moves.items() if comp != "c")
        if len(compartments) > 1 and abs(net) > 1.0:
            return None

        electrical = electrical_work(moves)
        return (chemistry, electrical) if split else chemistry + electrical

@lru_cache(maxsize=1)
def _contribution():
    from equilibrator_api import ComponentContribution

    return ComponentContribution()


class PreferredEnergies:
    """Component contribution where it reaches, group contribution where it does not.

    **Why a composite rather than a choice.** The two estimators are not equally good and
    they do not cover the same reactions. Over the 2,663 single-compartment yeast-GEM
    reactions, component contribution covers 1,289 and the ModelSEED group-contribution
    table 1,164; together they reach 1,355. And on the one reaction with an external check,
    component contribution is the one that is right: the cytosolic thiolase `r_0103` comes
    out at +24.96 +- 0.86 kJ/mol against a literature value near +26, where group
    contribution gives +38.07 +- 4.55.

    So this asks component contribution first and falls back. **Never within a reaction** --
    one estimator answers a reaction or none does. Mixing a formation energy from one method
    with a reaction energy from another produces a number belonging to neither, which is the
    same class of error as the kcal/kJ mix `bridge/thermodynamic.py::KJ_PER_KCAL` records.

    A refutation applies to the COMPOUND and is asked of both, so a species neither method
    should be trusted on stays refused whichever one reached it.
    """

    def __init__(self, component: EquilibratorData, group):
        self.component = component
        self.group = group

    @classmethod
    def load(cls, model=None, **kwargs) -> "PreferredEnergies":
        """Both backends over one model."""
        from .thermodynamic import ThermodynamicData

        return cls(EquilibratorData.load(model=model),
                   ThermodynamicData.load(model=model, **kwargs))

    def reaction_dg0(self, stoichiometry, with_error: bool = False):
        """Component contribution's answer, else group contribution's, else ``None``."""
        # Aliased: this method, the component backend's method and the group-contribution
        # FUNCTION are all called reaction_dg0, and three of one name in one scope is how a
        # call goes to the wrong one. tests/test_call_signatures.py caught exactly that.
        from .thermodynamic import reaction_dg0 as group_dg0
        from .thermodynamic import reaction_dg0_uncertainty as group_sigma

        value = self.component.reaction_dg0(stoichiometry, with_error=with_error)
        got = value[0] if with_error else value
        if got is not None:
            return value
        fallback = group_dg0(stoichiometry, self.group)
        if fallback is None:
            return (None, None) if with_error else None
        if not with_error:
            return fallback
        return fallback, group_sigma(stoichiometry, self.group)

    def covers(self, metabolite_id: str) -> bool:
        return (self.component.covers(metabolite_id)
                or self.group.dgf(metabolite_id) is not None)

    def refutation(self, metabolite_id: str):
        """Refused by either method's record is refused."""
        return (self.component.refutation(metabolite_id)
                or self.group.refutation(metabolite_id))

    def both_dg0(self, stoichiometry):
        """``(component, group)`` standard energies, or ``None`` unless BOTH answer.

        The raw pair, with no judgement attached. Whether a verdict survives the two
        disagreeing is a question about dG at real concentrations, not about dG0, so it is
        asked where dG is known -- `pathway/thermo_gate.py::_verdict_survives_refusal` --
        and not here.

        **This replaced a `disagreement()` that answered it here, on dG0, with the rule
        "the gap is smaller than the smaller estimate".** That rule was wrong and the
        counterexample is two lines of arithmetic: component -60 and group -35 pass it, and
        at a concentration ratio worth +40 kJ/mol they give -20 and +5 -- opposite verdicts
        from a check that called them robust. A threshold on dG0 cannot answer a question
        about dG.
        """
        from .thermodynamic import reaction_dg0 as group_dg0

        component = self.component.reaction_dg0(stoichiometry)
        group = group_dg0(stoichiometry, self.group)
        if component is None or group is None:
            return None
        return float(component), float(group)

    def which(self, stoichiometry) -> str:
        """Which estimator answers this reaction: ``'component'``, ``'group'`` or ``'none'``.

        Reported rather than inferred, because "the gate returned a number" does not say
        which of two methods produced it, and the two disagree by 13 kJ/mol on the one
        reaction anybody has checked.
        """
        from .thermodynamic import reaction_dg0 as group_dg0

        if self.component.reaction_dg0(stoichiometry) is not None:
            return "component"
        return "group" if group_dg0(stoichiometry, self.group) is not None else "none"
