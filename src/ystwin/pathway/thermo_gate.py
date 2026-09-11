"""Whether the pathway the solver just solved can run at the concentrations it just solved.

`pathway/solve.py` PRODUCES concentrations -- one content per node, mmol/gDCW -- and has no
opinion about whether the chemistry can proceed there, because a steady-state balance is
arithmetic on rates and no rate law can tell that a step is uphill.
`bridge/thermodynamic.py` CONSUMES concentrations: it computes dGr'0 for any reaction from
the vendored ModelSEED tables, and needs a Q before it can say anything about direction.
One module produces exactly what the other wants and nothing has ever joined them. This is
the join.

**The defect it prevents.** A step that cannot run at the concentrations the solver itself
produced is an internal contradiction, and without this it comes back as an ordinary float
in a table. yeast-GEM's cytosolic thiolase `r_0103` -- the chemical analogue of PhaA, the
entry step of the PHB pathway this package ships as `data/pathways/phb.toml` -- has
dGr'0 = +38.07 kJ/mol. Uphill, so it has a CONCENTRATION THRESHOLD, near 19 mM acetyl-CoA at
[CoA] = 100 uM and [acetoacetyl-CoA] = 1 uM, below which there is no net forward flux at all.
Measured cytosolic acetyl-CoA is about 10 uM on glucose and 425 uM on ethanol, and BOTH are
below it.

That last sentence used to read "on opposite sides of it", and the difference is this
module's most important correction. The energy was +15.6 kJ/mol and the threshold 221 uM
until 2026-08-30, when `bridge/thermodynamic.py` was found to be reading a kcal/mol table
under `pytfa`'s kJ/mol default. The straddle -- the whole thiolase lead, the only mechanism
of six that ever explained Kocharin's 4.24x -- existed only in the erroneous number.
eQuilibrator, which never touched that bug and whose +24.96 kJ/mol is the literature value
for a thiolase condensation, also puts ethanol below. `scripts/thiolase_threshold.py` carries
the provenance and the caveats.

What the gate itself lost is nothing: it still returns `cannot_run` on the glucose state and
now on the ethanol state too, and the wiring, the guards and the three states are unaffected.
A gate whose job is to refuse an infeasible step does not depend on the step being feasible
somewhere. What was lost is a biological claim, and it was lost to a check the repository
could have run against its own eQuilibrator backend at any time.

**Three states, and the third one is the point.** A step is reported as

    runs        dG < 0 at the given concentrations
    cannot_run  dG >= 0 -- the solver has contradicted itself
    cannot_say  no verdict could be reached, and :class:`Unresolved` says which of five
                reasons it was

**The reasons are different pieces of work, and two of them were being read as one.**
`cannot_say` used to be one state carrying a free-text note. `no_energy` is a coverage or a
NAMING gap and is closed by supplying an identifier: all three C40 carotenoids of
`data/pathways/beta_carotene.toml` reported it, and the tables carry all three -- what was
missing was that `fba/carotenoid.py` invents the metabolite ids, and an id invented here is
in no published alias table. `data/thermo/heterologous_metabolites.tsv` declares them and
the gap is closed.
`refuted_uncertainty` is the state underneath it, and it is a different claim: the energy
now resolves, and it still must not be used, because two estimators of it disagree by 91 to
275 kJ/mol while quoting 5 to 15. The report carries dGr'0 for such a step so a reader can
see what was withheld. `no_concentration` is the solver's business, `charge_translocation`
is a transport written without its counter-ion, and `uncertainty_exceeds_effect` is the
opt-in comparison of the driving force against
:func:`~ystwin.bridge.thermodynamic.reaction_dg0_uncertainty`. See
`docs/research/CAROTENOID_ENERGIES.md`.

:func:`~ystwin.bridge.thermodynamic.reaction_dg0` returns ``None`` rather than a partial sum
when a metabolite is uncovered, deliberately -- an energy computed from some of a reaction's
participants is a different reaction's energy. The tables cover 1840 of yeast-GEM's
4131 reactions -- 44.5%, NOT the 54.7% metabolite coverage it is easy to quote -- so
`cannot_say` is the common case rather than the exceptional one, and a gate that folded it
into `runs` would be worse than no gate: it would report a clearance it never earned.

**The unit gap is a convention, and it is required rather than defaulted.** The solver works
in mmol/gDCW and thermodynamics needs molar, and the conversion divides by a cytosolic
volume per gram dry weight that nobody here has measured.
`scripts/thiolase_threshold.py` reports across 1.0, 2.0 and 2.7 mL/gDW for exactly this
reason. So `cytosolic_volume_ml_per_gdcw` has no default anywhere in this module, and
:func:`gate_across_volumes` with :func:`volume_dependent_steps` exists to name the verdicts
that are artefacts of the choice. A verdict that flips across that range is not a verdict.

**What this does not do.** It does not change a flux, refuse a solve, or enter any published
number. It answers a question the solver could not ask. Wiring it into `solve.py` is a
separate decision with a separate cost, and the thiolase lead it was built to express rests
on two concentrations from two laboratories -- so what belongs in the chain today is the
question, not an answer that inherits that splice.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from ..bridge.thermodynamic import (
    GAS_CONSTANT_KJ,
    STANDARD_TEMPERATURE_K,
    ThermodynamicData,
    dg_prime,
    reaction_dg0,
    reaction_dg0_uncertainty,
    reaction_energy,
)
from .solve import PathwaySolution

__all__ = [
    "CYTOSOLIC_VOLUMES_ML_PER_GDCW",
    "MAXIMUM_PLAUSIBLE_MOLARITY",
    "Feasibility",
    "GateReport",
    "StepEnergy",
    "ThermodynamicallyBlocked",
    "Unresolved",
    "content_from_molar",
    "gate_across_volumes",
    "gate_pathway",
    "gate_step",
    "molar_from_content",
    "require_feasible",
    "substrate_threshold_m",
    "volume_dependent_steps",
]

MAXIMUM_PLAUSIBLE_MOLARITY = 55.5
"""Water is 55.5 M, so nothing dissolved in it is more concentrated. Arithmetic, not biology.

The same bound `solve.py` puts on a content for the same reason: a conversion with a free
unit in it will happily return 1e6 M from a volume given in the wrong unit, and every
verdict downstream is a comparison against zero that such a number wins.
"""

CYTOSOLIC_VOLUMES_ML_PER_GDCW = (1.0, 2.0, 2.7)
"""The plausible range of cytosolic volume per gram dry weight, mL/gDW.

A CONVENTION, not a measurement, and the same three values `scripts/thiolase_threshold.py`
reports its threshold across. Nothing in this repository measures it, so no single value is
offered as a default and every function that needs one takes it as a required argument. What
a result may rest on is agreement across this range, which is what
:func:`volume_dependent_steps` reports on.
"""


class Feasibility:
    """What a gated step can be said to do at the concentrations given.

    ``RUNS`` -- ``dG < 0``, so net forward flux is thermodynamically allowed. It says
    nothing about the rate, which is a different question with a different formalism.

    ``CANNOT_RUN`` -- ``dG >= 0``. Equilibrium counts here: at ``dG = 0`` there is no net
    flux in either direction, and the solver's balance assumed one.

    ``CANNOT_SAY`` -- no energy could be computed, or some participant has no concentration.
    Distinct from ``RUNS`` on purpose. The tables cover 44.5% of yeast-GEM's REACTIONS
    (1840/4131; the flattering 54.7% is metabolite coverage, a different denominator)
    metabolites, and reporting an uncovered reaction as feasible would turn a coverage gap
    into a clearance.
    """

    RUNS = "runs"
    CANNOT_RUN = "cannot_run"
    CANNOT_SAY = "cannot_say"
    ALL = (RUNS, CANNOT_RUN, CANNOT_SAY)


class Unresolved:
    """Why a step is ``cannot_say``. Five facts that used to be one state and a sentence.

    ``NO_ENERGY`` -- :func:`~ystwin.bridge.thermodynamic.reaction_dg0` returned None, so
    some participant is not in the tables *under the id it was asked for*. Two different
    fixes hide in that: a real coverage gap, which needs a measurement, and a naming gap,
    which needs an identifier. The three C40 carotenoids were the second and had been read
    as the first for as long as the spec has existed.

    ``CHARGE_TRANSLOCATION`` -- the reaction moves more net charge across a membrane than a
    carrier plausibly does in one turn, so its counter-ion is missing from the equation as
    written and the electrical term cannot be computed from it.

    ``NO_CONCENTRATION`` -- there is an energy but no Q: a participant has no concentration,
    or sits at exactly zero, which is what a ``passthrough`` node solves to and means.

    ``REFUTED_UNCERTAINTY`` -- there IS an energy, it resolves, and
    ``data/thermo/refuted_energies.tsv`` refuses it with a citation. **This is the state
    that must not collapse into ``NO_ENERGY``.** The two look identical from outside -- no
    verdict -- and they call for opposite work: an absent energy is closed by supplying an
    identifier or a number, and a refuted one is not closed by supplying anything except a
    MEASUREMENT, because the numbers on offer are extrapolations that disagree with each
    other by twenty times their quoted error. A step in this state reports its dGr'0 anyway,
    so a reader can see the size of what was withheld and check the claim.

    ``UNCERTAINTY_EXCEEDS_EFFECT`` -- a verdict was reached and then withdrawn, because
    ``|dG|`` is inside ``uncertainty_sigma`` times
    :func:`~ystwin.bridge.thermodynamic.reaction_dg0_uncertainty`. Opt-in: see
    :func:`gate_step`.
    """

    NO_ENERGY = "no_energy"
    CHARGE_TRANSLOCATION = "charge_translocation"
    NO_CONCENTRATION = "no_concentration"
    REFUTED_UNCERTAINTY = "refuted_uncertainty"
    UNCERTAINTY_EXCEEDS_EFFECT = "uncertainty_exceeds_effect"
    ALL = (NO_ENERGY, CHARGE_TRANSLOCATION, NO_CONCENTRATION, REFUTED_UNCERTAINTY,
           UNCERTAINTY_EXCEEDS_EFFECT)


class ThermodynamicallyBlocked(ValueError):
    """A solved step cannot run at the concentrations the solve produced.

    Raised by :func:`require_feasible` rather than returned, for the reason
    :class:`~ystwin.pathway.solve.ImplausibleContent` is raised: the alternative is a flag
    on a float that a caller puts in a table anyway. It means the pathway as solved is
    internally contradictory -- the balance moved carbon through a step whose driving force
    points the other way -- and not that the arithmetic is wrong.
    """


@dataclass(frozen=True)
class StepEnergy:
    """One step of a pathway, gated.

    Args:
        node: The node this step PRODUCES. A linear chain has one producing step per node,
            including the first, which is why steps are keyed this way rather than by the
            node they consume -- the entry step has no consumed node inside the chain.
        reaction_id: The GEM reaction, or ``""`` where the caller gave a bare stoichiometry.
            A bare stoichiometry carries no compartments, so no electrical work term is
            added to it; give a reaction for anything that crosses a membrane.
        feasibility: One of :class:`Feasibility`.
        dg0_kj_per_mol: Standard reaction energy, or ``None`` where none could be computed.
        dg_kj_per_mol: ``dG0 + RT ln Q`` at the concentrations used, or ``None``.
        concentrations_m: The molar concentrations Q was built from, per participant. Kept
            because the verdict is meaningless without them -- a threshold is a statement
            about all of the participants at once.
        note: Empty unless the verdict is ``CANNOT_SAY``, where it names what is missing
            and what would supply it.
        dg0_uncertainty_kj_per_mol: Standard error on ``dg0_kj_per_mol``, propagated over
            net structural cues by
            :func:`~ystwin.bridge.thermodynamic.reaction_dg0_uncertainty`, or ``None``
            where no decomposition covers the reaction. Reported next to every verdict
            rather than only next to the withdrawn ones, because a ``runs`` inside its own
            error bar is the case a reader most needs to see and the one a gate cannot see
            for them. Note that this IS kJ/mol while ``dg0_kj_per_mol``'s scale carries the
            caveat in :meth:`~ystwin.bridge.thermodynamic.ThermodynamicData.dgf`.
        reason: One of :class:`Unresolved` when ``feasibility`` is ``CANNOT_SAY``, and ``""``
            otherwise. The machine-readable half of ``note``: a caller that wants to treat
            a naming gap differently from a refuted energy compares this, and a caller that
            wants to explain it to a person prints ``note``.
    """

    node: str
    reaction_id: str
    feasibility: str
    dg0_kj_per_mol: float | None
    dg_kj_per_mol: float | None
    concentrations_m: dict[str, float] = field(default_factory=dict, repr=False)
    note: str = ""
    dg0_uncertainty_kj_per_mol: float | None = None
    reason: str = ""

    # Deliberately no `runs` boolean. Three states do not collapse into one predicate, and
    # `not step.runs` would read `cannot_say` as a failure exactly where this module exists
    # to keep them apart. Callers compare against :class:`Feasibility` and say which.


@dataclass(frozen=True)
class GateReport:
    """Every step of one solved pathway, at one cytosolic volume.

    ``ungated`` is part of the report rather than an omission from it: a gate that says
    "nothing is blocked" while covering one step of three has not cleared the pathway, and
    the only way for a caller to tell is for the report to say what it never looked at.
    """

    product: str
    cytosolic_volume_ml_per_gdcw: float
    temperature_k: float
    steps: tuple[StepEnergy, ...]
    ungated: tuple[str, ...] = ()

    def step(self, node: str) -> StepEnergy:
        for candidate in self.steps:
            if candidate.node == node:
                return candidate
        raise KeyError(f"no gated step producing {node!r}; "
                       f"have {[s.node for s in self.steps]}")

    @property
    def verdicts(self) -> dict[str, str]:
        """Node name -> feasibility, which is what comparing two volumes needs."""
        return {s.node: s.feasibility for s in self.steps}

    @property
    def runs(self) -> tuple[StepEnergy, ...]:
        return tuple(s for s in self.steps if s.feasibility == Feasibility.RUNS)

    @property
    def cannot_run(self) -> tuple[StepEnergy, ...]:
        return tuple(s for s in self.steps if s.feasibility == Feasibility.CANNOT_RUN)

    @property
    def cannot_say(self) -> tuple[StepEnergy, ...]:
        return tuple(s for s in self.steps if s.feasibility == Feasibility.CANNOT_SAY)

    def because(self, reason: str) -> tuple[StepEnergy, ...]:
        """The ``cannot_say`` steps with one :class:`Unresolved` reason.

        Raises:
            ValueError: on a reason that is not one of :attr:`Unresolved.ALL`, because a
                typo would return an empty tuple, and an empty tuple here reads as "no step
                failed that way" -- a positive claim from a query that could not have found
                one.
        """
        if reason not in Unresolved.ALL:
            raise ValueError(f"{reason!r} is not a reason a step can be unresolved for; "
                             f"have {list(Unresolved.ALL)}")
        return tuple(s for s in self.steps if s.reason == reason)

    @property
    def reasons(self) -> dict[str, str]:
        """Node name -> :class:`Unresolved` reason, for the steps that have one."""
        return {s.node: s.reason for s in self.steps if s.reason}

    def summary(self) -> str:
        """One line, with the volume in it because the verdicts are conditional on it."""
        return (f"{self.product}: {len(self.runs)} run, {len(self.cannot_run)} cannot run, "
                f"{len(self.cannot_say)} cannot say, {len(self.ungated)} not gated, at "
                f"{self.cytosolic_volume_ml_per_gdcw:g} mL/gDCW")


def molar_from_content(content_mmol_per_gdcw: float,
                       cytosolic_volume_ml_per_gdcw: float) -> float:
    """A solved content as a concentration, molar.

    The arithmetic is exact and the assumption is not: mmol/gDCW divided by mL/gDCW is
    mmol/mL, which is mol/L, so the conversion is one division and the whole of its content
    is the volume it divides by. Nothing in this repository measures that volume -- see
    :data:`CYTOSOLIC_VOLUMES_ML_PER_GDCW` -- which is why it is an argument and not a
    constant, here and everywhere below.

    Raises:
        ValueError: on a non-positive volume, which is not a small cell but no cell, and on
            a negative content, which the solver cannot produce.
    """
    if cytosolic_volume_ml_per_gdcw <= 0:
        raise ValueError(
            "cytosolic_volume_ml_per_gdcw must be positive, got "
            f"{cytosolic_volume_ml_per_gdcw}. A pool per gram of dry weight becomes a "
            f"concentration only by being given a volume to sit in; "
            f"{CYTOSOLIC_VOLUMES_ML_PER_GDCW} is the plausible range")
    if content_mmol_per_gdcw < 0:
        raise ValueError(f"content must be non-negative, got {content_mmol_per_gdcw}")
    molar = float(content_mmol_per_gdcw / cytosolic_volume_ml_per_gdcw)
    if molar > MAXIMUM_PLAUSIBLE_MOLARITY:
        raise ValueError(
            f"{content_mmol_per_gdcw:g} mmol/gDCW in {cytosolic_volume_ml_per_gdcw:g} "
            f"mL/gDCW is {molar:g} M, above the {MAXIMUM_PLAUSIBLE_MOLARITY:g} M of pure "
            "water. This is arithmetic and not biology, so it needs no citation, and the "
            "way it happens is a unit slip: passing the volume in litres rather than mL is "
            "a factor of 1000 and lands a blocked step comfortably on the other side of "
            "zero, because ln Q moves with the log of the concentration and never refuses. "
            f"{CYTOSOLIC_VOLUMES_ML_PER_GDCW} mL/gDCW is the range this module assumes")
    return molar


def content_from_molar(concentration_m: float,
                       cytosolic_volume_ml_per_gdcw: float) -> float:
    """A concentration back in the solver's units, mmol/gDCW.

    The inverse of :func:`molar_from_content`, and it exists so a threshold can be quoted in
    the units the caller's own numbers are in. A threshold in molar is not actionable by a
    layer that never sees a molarity.
    """
    if cytosolic_volume_ml_per_gdcw <= 0:
        raise ValueError(
            "cytosolic_volume_ml_per_gdcw must be positive, got "
            f"{cytosolic_volume_ml_per_gdcw}")
    if concentration_m < 0:
        raise ValueError(f"concentration must be non-negative, got {concentration_m}")
    return float(concentration_m * cytosolic_volume_ml_per_gdcw)


def gate_step(
    stoichiometry: Mapping[str, float],
    concentrations_m: Mapping[str, float],
    thermo: ThermodynamicData,
    *,
    node: str = "",
    produces: str = "",
    reaction=None,
    temperature_k: float = STANDARD_TEMPERATURE_K,
    uncertainty_sigma: float | None = None,
) -> StepEnergy:
    """Gate one reaction at one set of concentrations.

    The whole of the physics is two lines -- ``dG = dG0 + RT ln Q`` and a comparison against
    zero -- and everything else here is bookkeeping about what could not be computed.

    Args:
        stoichiometry: Coefficients keyed by GEM metabolite id, negative for substrates.
        concentrations_m: Molar, for every participant with a non-zero coefficient. A
            participant that is missing, or that is at exactly zero, gives ``CANNOT_SAY``
            rather than an exception, because a gate is asked about pathways it may not
            fully know.
        thermo: Formation energies. ``ThermodynamicData.load()`` builds them.
        node: The node this step produces, for reporting.
        produces: The GEM metabolite id of that node. Supplying it checks the reaction is
            written in the direction that MAKES the node; see :func:`_check_direction` for
            why an unchecked direction silently reverses the verdict's sign.
        reaction: The cobra reaction, when there is one. Must agree with ``stoichiometry``
            when both are given -- a mismatched pair takes dG0 from one reaction and Q from
            another, which is checked rather than trusted. Supplying it adds the electrical
            work term for anything crossing a membrane, and inherits
            :func:`~ystwin.bridge.thermodynamic.reaction_energy`'s refusal of a reaction
            whose net charge translocation is larger than a carrier plausibly does.
        temperature_k: Absolute temperature. Defaults to 30 C, the cultivation temperature
            everything in this package is calibrated at.
        uncertainty_sigma: Withdraw a verdict whose driving force is inside this many
            standard errors of the estimated energy, reporting
            ``Unresolved.UNCERTAINTY_EXCEEDS_EFFECT``. **No default, and off unless asked
            for.** Withdrawing a verdict is a claim about that verdict, and it belongs to
            whoever owns the result -- made deliberately, not smuggled in as a load-time
            default. The error is reported on every step either way, which is the part that
            does not need permission.

            **The example that used to justify this is withdrawn, and the default is now
            argued on the principle alone.** This paragraph read, until 2026-09-04, that
            `phb`'s entry step "runs at -3.3 kJ/mol on ethanol against a cue-propagated
            error of 4.5, so a one-sigma rule applied by default would silently erase it".
            That was never true as written: commit 8809344 introduced the sentence, and the
            same commit's message records the correction that killed it -- a kcal/kJ unit
            error meant thiolase r_0103 moved +15.6 -> +38.07 kJ/mol and its threshold
            221 uM -> 19 mM, "so ethanol's measured 425 uM no longer clears it". Measured
            through the shipped path today the step is +38.07 kJ/mol against sigma 4.55,
            which is 8.4 sigma UPHILL and returns ``cannot_run``. A one-sigma rule would
            erase nothing, because there is no verdict there to erase. Kept off by default
            because the principle above stands on its own; not because of that number.
            It is separate
            from ``refuted_uncertainty``, which is not a threshold and is never optional:
            there the quoted error is the thing that was refuted, so no multiple of it
            would mean anything.
    """
    reaction_id = getattr(reaction, "id", "") if reaction is not None else ""
    stoichiometry = {str(m): float(c) for m, c in stoichiometry.items()}
    _check_pair_agrees(stoichiometry, reaction)
    _check_direction(stoichiometry, produces, node, reaction_id)
    dg0, note, reason = _standard_energy(stoichiometry, thermo, reaction)
    if dg0 is None:
        return StepEnergy(node, reaction_id, Feasibility.CANNOT_SAY, None, None, {}, note,
                          None, reason)
    sigma = _backend_sigma(stoichiometry, thermo)

    # Before Q, because this is a fact about the ENERGY and not about the solve. A step
    # whose energy is refused would otherwise be reported as missing a concentration, and
    # supplying the concentration would then produce a verdict off the refused number.
    # THE DISAGREEMENT CHECK RUNS ON EVERY REACTION, not only on the three compounds
    # somebody hand-listed. `data/thermo/refuted_energies.tsv` was written after looking at
    # the carotenoids and nothing else, and the standard it applied was not the one the rest
    # of the model was held to: CRTYB_PSY was refused at a 22.2 kJ/mol disagreement -- BELOW
    # the model-wide median of 25.2 -- while `r_4710` was reported `runs` with the two
    # estimators at -963.0 and +1.4, opposite signs and 964 kJ/mol apart, because no
    # participant of it was on the list. A curated list cannot be a standard; a measurement
    # can. Over the 1,098 single-compartment reactions both estimators cover the gap runs
    # 7.3 / 25.2 / 46.8 / 81.8 kJ/mol at the 25th / 50th / 75th / 90th percentile, so
    # refusing a step for 22 would refuse most of yeast-GEM.
    pair = _both_estimates(stoichiometry, thermo)
    refused = _refuted_note(stoichiometry, thermo)
    if refused and pair is None:
        # Refused and only one estimator to ask, so there is nothing to check the verdict
        # against. This is the conservative branch and it is the default: a single method
        # agreeing with itself is not agreement.
        return StepEnergy(node, reaction_id, Feasibility.CANNOT_SAY, float(dg0), None, {},
                          refused, sigma, Unresolved.REFUTED_UNCERTAINTY)

    participants = {m: c for m, c in stoichiometry.items() if c != 0}
    used = {m: float(concentrations_m[m]) for m in participants if m in concentrations_m}
    missing = sorted(m for m in participants if m not in used)
    # Exactly zero, not merely non-positive: a NEGATIVE concentration is a caller's bug and
    # `dg_prime` raises on it below, which is the right outcome. Zero is what a passthrough
    # node legitimately solves to, and turning that into an exception would make the gate
    # unusable on half the specs this package ships.
    empty = sorted(m for m, value in used.items() if value == 0)
    if missing or empty:
        return StepEnergy(node, reaction_id, Feasibility.CANNOT_SAY, float(dg0), None, used,
                          _missing_note(missing, empty), sigma, Unresolved.NO_CONCENTRATION)

    dg = dg_prime(float(dg0), participants, used, temperature_k)
    if (uncertainty_sigma is not None and sigma is not None
            and abs(dg) <= uncertainty_sigma * sigma):
        return StepEnergy(
            node, reaction_id, Feasibility.CANNOT_SAY, float(dg0), float(dg), used,
            f"dG = {dg:+.4g} kJ/mol is inside {uncertainty_sigma:g} x {sigma:.4g} kJ/mol, "
            "the standard error of this reaction's estimated energy propagated over its net "
            "structural cues. The sign of a driving force smaller than the error on it is a "
            "property of the estimate, not of the chemistry. A measured energy, or a "
            "substrate pool far enough from the threshold that the verdict survives the "
            "error, is what would settle it", sigma, Unresolved.UNCERTAINTY_EXCEEDS_EFFECT)
    feasibility = Feasibility.RUNS if dg < 0 else Feasibility.CANNOT_RUN
    if pair is not None:
        # THE CHECK IS ON dG, NOT dG0, and that is the whole of it. Both estimators are put
        # at the SAME concentrations and asked for a verdict; if they give the same one, the
        # refused magnitude did not decide it and the verdict stands. If they differ, it did,
        # and nothing may be said. No threshold, because a threshold on dG0 cannot answer a
        # question about dG -- component -60 and group -35 pass "the gap is smaller than the
        # smaller estimate" and then disagree at a Q worth +40 kJ/mol.
        component, group = pair
        other = dg + (group - component) if abs(dg0 - component) < abs(dg0 - group) \
            else dg + (component - group)
        prefix = f"{refused}\n\n" if refused else ""
        if (dg < 0) != (other < 0):
            return StepEnergy(
                node, reaction_id, Feasibility.CANNOT_SAY, float(dg0), float(dg), used,
                f"{prefix}The two estimators do not agree on the verdict: at these "
                f"concentrations they give dG = {dg:+.4g} and {other:+.4g} kJ/mol, which are "
                "opposite directions. Group contribution and component contribution differ "
                "across this whole model -- median 25.2 kJ/mol over the 1,098 reactions both "
                "cover -- and this is a step where that difference decides the answer",
                sigma, Unresolved.REFUTED_UNCERTAINTY)
        if refused:
            return StepEnergy(
                node, reaction_id, feasibility, float(dg0), float(dg), used,
                f"{prefix}The MAGNITUDE is refused and the VERDICT is not: at these "
                f"concentrations both estimators give the same sign, dG = {dg:+.4g} against "
                f"{other:+.4g} kJ/mol. Read the direction, not the number", sigma)
    return StepEnergy(node, reaction_id, feasibility, float(dg0), float(dg), used, "", sigma)


def gate_pathway(
    solution: PathwaySolution,
    *,
    steps: Mapping[str, object],
    metabolites: Mapping[str, str],
    thermo: ThermodynamicData,
    cytosolic_volume_ml_per_gdcw: float,
    background_m: Mapping[str, float] | None = None,
    model=None,
    temperature_k: float = STANDARD_TEMPERATURE_K,
    uncertainty_sigma: float | None = None,
) -> GateReport:
    """Gate every declared step of a solved pathway at that solve's own concentrations.

    This is the function a solver can call. It does not modify the solution and does not
    raise on an infeasible step -- :func:`require_feasible` is the refusal, kept separate so
    that reporting and refusing are two decisions rather than one.

    Args:
        solution: What :func:`~ystwin.pathway.solve.solve_pathway` returned.
        steps: Node name -> the reaction that PRODUCES it. Keyed by the produced node
            because the entry step's substrate is the precursor, which is not a node in the
            chain, so keying by the consumed node would leave the entry step -- the one the
            thiolase lead is about -- unreachable. A value is a GEM reaction id (needs
            ``model``), a cobra reaction, or a bare ``{metabolite id: coefficient}`` mapping.
            A partial mapping is normal and the nodes left out are reported as ``ungated``.
        metabolites: Node name -> GEM metabolite id, for the nodes whose solved content
            should enter Q. Required rather than defaulted: which species a node IS is a
            claim, and a wrong one silently gates a different reaction.
        thermo: Formation energies.
        cytosolic_volume_ml_per_gdcw: mL of cytosol per gram dry weight, used to turn every
            solved content into a molarity. **No default** -- see
            :data:`CYTOSOLIC_VOLUMES_ML_PER_GDCW`, and prefer :func:`gate_across_volumes`.
        background_m: Molar concentrations for everything the solution does not pin --
            cofactors, the precursor the pathway is pulled off, and any node the solver
            gave no pool. Consulted only there: a background value for a node that solved
            to a positive content is refused rather than silently preferred, because the
            question this module answers is whether the SOLVER's concentrations are
            consistent, and overriding one quietly would answer a different question.
        model: A cobra model, needed only to resolve reaction ids in ``steps``.
        temperature_k: Absolute temperature.
        uncertainty_sigma: Passed to :func:`gate_step`; off unless asked for, and that
            function's docstring says why.

    Raises:
        KeyError: for a step, or a metabolite, keyed by a node the solution does not have.
        ValueError: for a reaction id with no model to look it up in, a background
            concentration that clashes with a solved one, or two nodes declared to be the
            same species.
    """
    background = dict(background_m or {})
    _check_nodes_exist(solution, metabolites)
    solved: dict[str, float] = {}
    solved_to_zero: set[str] = set()
    claimed_by: dict[str, str] = {}
    for state in solution.nodes:
        metabolite = metabolites.get(state.name)
        if metabolite is None:
            continue
        if metabolite in claimed_by:
            raise ValueError(
                f"nodes {claimed_by[metabolite]!r} and {state.name!r} are both declared to "
                f"be {metabolite!r}, so one of their contents would silently replace the "
                "other and which one depends on the order the spec lists them in. A node "
                "IS a species; give them different metabolite ids, or drop the one whose "
                "content is not the pool being gated")
        claimed_by[metabolite] = state.name
        if state.content_mmol_per_gdcw > 0:
            solved[metabolite] = molar_from_content(state.content_mmol_per_gdcw,
                                                    cytosolic_volume_ml_per_gdcw)
        else:
            solved_to_zero.add(metabolite)

    clash = sorted(set(solved) & set(background))
    if clash:
        raise ValueError(
            f"background_m gives concentrations for {clash}, which the solve already pins "
            "from its own contents. This gate asks whether the solver's concentrations are "
            "self-consistent, so preferring either silently would answer a different "
            "question. Drop them from background_m, or drop the node from `metabolites` and "
            "say plainly that its content is not being tested")
    # Three layers, and the order is the argument. A node the solver pinned at zero starts
    # as a literal zero so the report can say WHY there is no concentration rather than
    # merely that there is none; `background_m` may replace that with a measurement; a
    # positively solved content wins outright, and the clash check above is what stops
    # anything from quietly replacing one.
    concentrations = {**dict.fromkeys(solved_to_zero, 0.0), **background, **solved}

    gated: list[StepEnergy] = []
    for node, declared in steps.items():
        solution.node(node)  # raises, naming the nodes, if the key is not one of them
        stoichiometry, reaction = _resolve(node, declared, model)
        gated.append(gate_step(stoichiometry, concentrations, thermo, node=node,
                               produces=metabolites.get(node, ""), reaction=reaction,
                               temperature_k=temperature_k,
                               uncertainty_sigma=uncertainty_sigma))

    ungated = tuple(s.name for s in solution.nodes if s.name not in steps)
    return GateReport(solution.product, float(cytosolic_volume_ml_per_gdcw),
                      float(temperature_k), tuple(gated), ungated)


def gate_across_volumes(
    solution: PathwaySolution,
    *,
    volumes: tuple[float, ...] = CYTOSOLIC_VOLUMES_ML_PER_GDCW,
    **kwargs,
) -> tuple[GateReport, ...]:
    """One report per cytosolic volume, because one report per pathway would hide the choice.

    Takes the same keyword arguments as :func:`gate_pathway` except
    ``cytosolic_volume_ml_per_gdcw``, which is what it sweeps. Pair it with
    :func:`volume_dependent_steps`, which names the conclusions that did not survive.
    """
    if not volumes:
        raise ValueError("give at least one cytosolic volume to report across")
    return tuple(gate_pathway(solution, cytosolic_volume_ml_per_gdcw=volume, **kwargs)
                 for volume in volumes)


def volume_dependent_steps(reports: tuple[GateReport, ...]) -> tuple[str, ...]:
    """The steps whose verdict is not the same in every report.

    These are the ones that rest on the volume convention rather than on the chemistry, and
    naming them is the whole reason the sweep exists. `scripts/thiolase_threshold.py` makes
    the same move on the same conversion and for the same reason: what survived there was
    the ORDERING of two feeds, not the threshold's absolute value.
    """
    if len(reports) < 2:
        raise ValueError(
            f"volume dependence needs at least two volumes to compare, got {len(reports)}. "
            "A single-report sweep returns no dependent steps for the same reason a "
            "one-point line has no slope, and the empty tuple then reads as 'no verdict "
            "rests on the volume convention' -- a positive claim, produced by a check "
            f"incapable of detecting the thing it denies. {CYTOSOLIC_VOLUMES_ML_PER_GDCW} "
            "is what `gate_across_volumes` sweeps by default")
    nodes = sorted({node for report in reports for node in report.verdicts})
    return tuple(node for node in nodes
                 if len({report.verdicts.get(node) for report in reports}) > 1)


def require_feasible(report: GateReport) -> None:
    """Raise if any gated step cannot run at the solved concentrations.

    Does NOT raise on ``cannot_say``, and that is deliberate rather than lenient: refusing a
    pathway because the tables do not cover it would be a refusal with no energy behind it,
    which is a fabricated result in the shape of a caution. Those steps are in
    ``report.cannot_say``, and a caller that wants them treated as failures can say so in
    one line.

    Raises:
        ThermodynamicallyBlocked: naming every blocked step, its energies and the
            concentrations they were computed at.
    """
    if not report.steps:
        raise ThermodynamicallyBlocked(
            f"nothing was gated for {report.product!r}, so there is no feasibility here to "
            f"require: {len(report.ungated)} node(s) {list(report.ungated)} went untested. "
            "Passing on an empty report is how a gate reports success for a pathway it "
            "never looked at. Declare the steps in `steps=`, or do not call this")
    blocked = report.cannot_run
    if not blocked:
        return
    lines = []
    for step in blocked:
        levels = ", ".join(f"[{m}] = {value * 1e6:.4g} uM"
                           for m, value in sorted(step.concentrations_m.items()))
        lines.append(
            f"  {step.node!r} via {step.reaction_id or 'the given stoichiometry'}: "
            f"dGr'0 = {step.dg0_kj_per_mol:+.4g} kJ/mol, dG = {step.dg_kj_per_mol:+.4g} "
            f"kJ/mol at {levels}")
    raise ThermodynamicallyBlocked(
        f"{len(blocked)} step(s) of {report.product!r} cannot run at the concentrations "
        f"this solve produced:\n" + "\n".join(lines) + "\n"
        "The balance moved carbon through them anyway, so the solution is internally "
        "contradictory rather than merely unlikely. Three things lift it and they are "
        "different claims: a substrate pool above its threshold (`substrate_threshold_m` "
        "returns it in molar, `content_from_molar` puts it back in mmol/gDCW), measured "
        "rather than assumed levels for the cofactors in Q, or a cytosolic volume other "
        f"than the {report.cytosolic_volume_ml_per_gdcw:g} mL/gDCW assumed here -- check "
        "`volume_dependent_steps` before believing the last one")


def substrate_threshold_m(
    stoichiometry: Mapping[str, float],
    metabolite: str,
    concentrations_m: Mapping[str, float],
    dg0_kj_per_mol: float,
    *,
    temperature_k: float = STANDARD_TEMPERATURE_K,
) -> float:
    """The concentration of one participant at which dG crosses zero, molar.

    Solving ``dG0 + RT ln Q = 0`` for one participant while the rest are held::

        ln x_i = ( -dG0/RT - sum_(j != i) c_j ln x_j ) / c_i

    For a SUBSTRATE (negative coefficient) the reaction runs forward ABOVE the value
    returned; for a PRODUCT it runs forward below it. The coefficient divides, which is why
    the thiolase's threshold is a square root -- it consumes its substrate two at a time, so
    halving the pool quarters the driving term and the threshold is correspondingly sharp.

    Args:
        stoichiometry: Coefficients keyed by metabolite id.
        metabolite: The participant to solve for. Its own concentration is ignored.
        concentrations_m: Molar, for every OTHER participant with a non-zero coefficient.
        dg0_kj_per_mol: The standard energy, from :attr:`StepEnergy.dg0_kj_per_mol` or
            :func:`~ystwin.bridge.thermodynamic.reaction_dg0`.
        temperature_k: Absolute temperature.

    Raises:
        ValueError: if ``metabolite`` does not participate, or if its coefficient is zero,
            in which case its concentration does not appear in Q and no threshold exists.
        KeyError: if another participant has no concentration, so Q is undefined.
    """
    coefficient = stoichiometry.get(metabolite)
    if coefficient is None:
        raise ValueError(f"{metabolite!r} does not take part in this reaction; "
                         f"have {sorted(stoichiometry)}")
    if coefficient == 0:
        raise ValueError(
            f"{metabolite!r} has a coefficient of zero, so it does not appear in Q and no "
            "concentration of it changes the direction of anything")
    others = 0.0
    for participant, other in stoichiometry.items():
        if participant == metabolite or other == 0:
            continue
        if participant not in concentrations_m:
            raise KeyError(
                f"no concentration given for {participant!r}, so the threshold on "
                f"{metabolite!r} is undefined -- it is a statement about all of the other "
                "participants at once, not a property of this metabolite")
        value = concentrations_m[participant]
        if value <= 0:
            raise ValueError(
                f"concentration of {participant!r} must be positive, got {value}")
        others += other * np.log(value)
    driving = -dg0_kj_per_mol / (GAS_CONSTANT_KJ * temperature_k)
    return float(np.exp((driving - others) / coefficient))


def _check_pair_agrees(stoichiometry: dict[str, float], reaction) -> None:
    """Refuse a ``(stoichiometry, reaction)`` pair that describes two different reactions.

    ``gate_step`` takes dG0 from ``reaction`` when one is given and forms Q from
    ``stoichiometry`` always. Nothing made them the same reaction, so a caller holding two
    reactions could pair either energy with either quotient and get a verdict that is a
    property of neither: pairing hexokinase's dGr'0 = +36.72 kJ/mol with the thiolase
    quotient turns glucose from `cannot_run` (+15.63 kJ/mol) into `runs` (-13.40), and the
    report names the thiolase throughout. `gate_pathway` builds both from one `_resolve`
    call and so cannot make this mistake; a direct caller of the exported `gate_step` can.
    """
    if reaction is None:
        return
    theirs = {m.id: float(c) for m, c in reaction.metabolites.items()}
    if theirs == stoichiometry:
        return
    only_given = sorted(set(stoichiometry) - set(theirs))
    only_reaction = sorted(set(theirs) - set(stoichiometry))
    differing = sorted(m for m in set(theirs) & set(stoichiometry)
                       if theirs[m] != stoichiometry[m])
    raise ValueError(
        f"the stoichiometry given and reaction {reaction.id!r} are different reactions, so "
        "dG0 would come from one and Q from the other, and the verdict would be a property "
        f"of neither. Only in the mapping: {only_given}; only in the reaction: "
        f"{only_reaction}; different coefficients: {differing}. Pass the stoichiometry the "
        "reaction actually has, or pass no reaction and gate the chemistry alone")


def _check_direction(stoichiometry: dict[str, float], produces: str, node: str,
                     reaction_id: str) -> None:
    """Refuse a step written in the direction that CONSUMES the node it declares it makes.

    The sign of dG0 is the sign of the direction the reaction is written in, and an SBML
    file's direction is an authoring convention rather than a claim about the cell: yeast-
    GEM writes the thiolase r_0103 reversibly, bounds (-1000, 1000), and as written it gives
    dGr'0 = +15.606 kJ/mol -- the whole of the thiolase lead. Written the other way it is
    -15.606 and glucose clears the threshold that lead is about. Nothing downstream can
    recover the intended direction, so it is checked at the only place that knows it: the
    step is declared to PRODUCE `node`, so `node`'s species must carry a positive
    coefficient.
    """
    if not produces:
        return
    coefficient = stoichiometry.get(produces)
    if coefficient is None:
        raise ValueError(
            f"step for {node!r} is declared to produce {produces!r}, which does not appear "
            f"in {reaction_id or 'the given stoichiometry'} at all "
            f"(participants: {sorted(stoichiometry)}). A metabolite id that matches nothing "
            "is a typo, and a typo here leaves the node's solved content out of Q entirely "
            "-- `background_m` then supplies that species and the gate reports a verdict "
            "about an assumed level while saying it tested the solve's")
    if coefficient <= 0:
        raise ValueError(
            f"step for {node!r} is declared to produce {produces!r}, but it appears with "
            f"coefficient {coefficient:g} in "
            f"{reaction_id or 'the given stoichiometry'} -- that is the direction that "
            "CONSUMES it. dG0's sign is the sign of the direction as written, so gating "
            "this would report the reverse reaction's verdict under the forward one's name. "
            "Reverse the coefficients (and the sign of dG0 follows), or correct which node "
            "this step produces")


def _check_nodes_exist(solution: PathwaySolution, metabolites: Mapping[str, str]) -> None:
    """Refuse a ``metabolites`` key that is not a node of this solution.

    The one way the guarantee above can be broken quietly. A mistyped node name maps
    nothing, so the solve's own content never enters Q, `background_m` supplies that
    metabolite instead, and the clash check never fires because no solved value claimed it
    -- the gate then reports a verdict about an assumed level while saying it tested the
    solver's. On the synthetic chain that turns `runs` at 45 mM into `cannot_run` at 1 uM
    with nothing in the report to say why. Steps are already checked one loop below; this is
    the same check for the other mapping.
    """
    unknown = sorted(set(metabolites) - {state.name for state in solution.nodes})
    if not unknown:
        return
    raise KeyError(
        f"`metabolites` maps {unknown}, which no node of {solution.product!r} is named; "
        f"have {[state.name for state in solution.nodes]}. A node absent from that mapping "
        "is one whose content is not being tested, which is a legitimate thing to say and "
        "is said by leaving it out -- so a name that matches nothing is a typo, and a typo "
        "here silently swaps a solved content for a background level. Correct the name, or "
        "drop the entry")


def _backend_covers(thermo, metabolite_id: str) -> bool:
    """Whether this backend can reach the metabolite, without asking it for an energy.

    `bridge/equilibrator.py` has no per-compound energy to give and raises if asked, so a
    coverage probe written as ``thermo.dgf(m) is None`` cannot be used on it. Its ``covers``
    is the question actually being asked here.
    """
    probe = getattr(thermo, "covers", None)
    if probe is not None:
        return bool(probe(metabolite_id))
    return thermo.dgf(metabolite_id) is not None


def _backend_dg0(stoichiometry: dict[str, float], thermo):
    """The reaction's standard energy, from whichever kind of estimator this is.

    Two shapes exist and they are not interchangeable. `bridge/thermodynamic.py` estimates
    per COMPOUND -- group contribution -- so its reaction energy is a sum of formation
    energies and the module-level :func:`reaction_dg0` does that sum. Component
    contribution estimates per REACTION and has no compound term to sum; asking it for one
    is meaningless, which is why `EquilibratorData.dgf` refuses rather than returning the
    ``0.0`` it used to.

    A backend that brings its own ``reaction_dg0`` is therefore asked directly. Nothing is
    ever mixed WITHIN a reaction: one estimator answers a reaction or none does.
    """
    own = getattr(thermo, "reaction_dg0", None)
    if callable(own):
        return own(stoichiometry)
    return reaction_dg0(stoichiometry, thermo)


def _backend_sigma(stoichiometry: dict[str, float], thermo) -> float | None:
    """The standard error on that energy, from whichever estimator produced it.

    Component contribution carries its own covariance and reports an error per reaction;
    group contribution's is propagated over net structural cues by
    :func:`reaction_dg0_uncertainty`. Taking the cue-propagated error for a component
    contribution energy would pair a number with someone else's uncertainty.
    """
    own = getattr(thermo, "reaction_dg0", None)
    if callable(own):
        try:
            _, error = own(stoichiometry, with_error=True)
        except TypeError:
            return None
        return None if error is None else float(error)
    return reaction_dg0_uncertainty(stoichiometry, thermo)


def _standard_energy(stoichiometry: dict[str, float], thermo: ThermodynamicData,
                     reaction) -> tuple[float | None, str, str]:
    """``(dG0, note, reason)``, with ``None`` and a note naming what is missing.

    Two ways to have no energy and they are different facts, so they get different notes
    and different :class:`Unresolved` reasons: the tables not covering a metabolite, and
    `reaction_energy` refusing a charge translocation larger than a carrier plausibly
    performs.
    """
    chemistry = _backend_dg0(stoichiometry, thermo)
    if chemistry is None:
        uncovered = sorted(m for m in stoichiometry if not _backend_covers(thermo, m))
        return None, (
            f"no standard energy: the thermodynamic tables do not cover {uncovered}. "
            "reaction_dg0 returns None rather than a partial sum on purpose, because an "
            "energy summed over some of a reaction's metabolites is a different reaction's. "
            "Two different things produce this and they need different work: a coverage "
            "gap, which needs a measurement, and a NAMING gap, where the tables carry the "
            "species under an id nobody declared -- the case every C40 carotenoid was in, "
            "closed by data/thermo/heterologous_metabolites.tsv. Check the second before "
            "concluding the first"), Unresolved.NO_ENERGY
    if reaction is None:
        return float(chemistry), "", ""
    total = reaction_energy(reaction, thermo)
    if total is None:
        return None, (
            "no standard energy: reaction_energy refuses this one because its net charge "
            "across a membrane is larger than a carrier plausibly moves in one turn, which "
            "means the counter-ion is missing from the equation as written. Write the "
            "transport with its counter-ion, or gate the chemistry alone by passing the "
            "stoichiometry without the reaction"), Unresolved.CHARGE_TRANSLOCATION
    return float(total), "", ""


def _both_estimates(stoichiometry: dict[str, float], thermo):
    """``(component, group)`` from a backend that carries two estimators, else ``None``.

    Only `bridge/equilibrator.py::PreferredEnergies` answers. A single-estimator backend
    returns ``None`` and the caller takes the conservative branch, which is the point: one
    method agreeing with itself is not agreement, and the default behaviour of the gate is
    unchanged for every caller that does not deliberately supply two.
    """
    probe = getattr(thermo, "both_dg0", None)
    return probe(stoichiometry) if probe is not None else None


def _refuted_note(stoichiometry: dict[str, float], thermo: ThermodynamicData) -> str:
    """Why this step's energy must not be used, or ``""`` if nothing refuses it.

    The energy is present -- it is in ``dg0_kj_per_mol`` on the returned step -- and the
    note is what stops it being read as a verdict. Naming every refused participant and
    quoting the row's own reason and source, because "the gate declined" with no citation
    is the shape of a result nobody can check or overturn.
    """
    refused = [(m, thermo.refutation(m)) for m in sorted(stoichiometry)
               if thermo.refutation(m) is not None]
    if not refused:
        return ""
    lines = [f"{m} ({record.name}, {record.modelseed_id}): {record.reason} "
             f"[{record.source}]" for m, record in refused]
    return (
        "a standard energy exists for every participant and this one is refused anyway. "
        f"{len(refused)} of them carry a formation energy whose quoted uncertainty is "
        "refuted by evidence outside the estimator that quoted it, so dGr'0 is reported "
        "above but no verdict is taken from it. This is NOT the tables failing to cover "
        "the species -- they cover it. What would lift this is a MEASURED formation "
        "energy, not another estimate. " + " ".join(lines))


def _missing_note(missing: list[str], empty: list[str]) -> str:
    """Why Q could not be formed, and what would supply it."""
    parts = []
    if missing:
        parts.append(
            f"no concentration for {missing}, so Q is undefined. These are the participants "
            "no node of the solve maps onto -- cofactors and the precursor the pathway is "
            "pulled off -- and `background_m` is where they go")
    if empty:
        parts.append(
            f"{empty} sits at exactly zero, which is not a concentration -- ln Q diverges "
            "there. A node whose spec declares `passthrough` solves to exactly this, and "
            "means it: the spec asserts that intermediate holds no appreciable pool. "
            "Declare a rate law that keeps a pool, or give the node a measured level in "
            "`background_m`")
    return "; ".join(parts)


def _resolve(node: str, declared: object, model) -> tuple[dict[str, float], object]:
    """A declared step as ``(stoichiometry, reaction or None)``."""
    if isinstance(declared, str):
        if model is None:
            raise ValueError(
                f"step {node!r} names reaction {declared!r} but no model was given to look "
                "it up in. Pass `model=` a cobra model, or declare the step as a "
                "{metabolite id: coefficient} mapping -- which is also how a step gets "
                "gated without an SBML parse, and how the cheap half of this module's tests "
                "run on a machine with no GEM at all")
        reaction = model.reactions.get_by_id(declared)
        return {m.id: c for m, c in reaction.metabolites.items()}, reaction
    if hasattr(declared, "metabolites"):
        return {m.id: c for m, c in declared.metabolites.items()}, declared
    if isinstance(declared, Mapping):
        return {str(m): float(c) for m, c in declared.items()}, None
    raise TypeError(
        f"step {node!r} is a {type(declared).__name__}; give a GEM reaction id, a cobra "
        "reaction, or a {metabolite id: coefficient} mapping")
