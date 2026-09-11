"""Coupling a measured cofactor ratio to the metabolic model, as far as the evidence allows.

The obvious coupling is refuted. A sensor reads a pool, and a pool does not set a flux:
Larsson 1997 (PMID 9393686) found glycolytic flux correlates negatively with intracellular
ATP in yeast chemostats and not at all with the ATP/ADP ratio, and Vemuri 2007
(PMID 17287356) drained cytosolic NADH enough to abolish most glycerol production without
moving the critical dilution rate at all. Cytosolic redox follows the flux; it does not
govern it.

What has precedent is a thermodynamic constraint. Kummel 2006 (PMID 16788595) put the
adenylate energy charge and the NADH/NAD ratio on iND750, and Martinez 2014 (PMID 25028891)
did the same on Yeast 5, bounding Gibbs energies to constrain reaction DIRECTIONALITY --
which way a reaction may run, not how fast. Both chose bounds deliberately loose enough not
to bind, so making them bind from measurement is the part that is new.

This module supplies the energies. Whether a measured ratio actually narrows anything is a
question for the model, and one worth answering honestly either way.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from .. import paths

__all__ = [
    "COMPARTMENTS",
    "FARADAY_KJ_PER_V",
    "FREE_MAGNESIUM_M",
    "GAS_CONSTANT_KJ",
    "KJ_PER_KCAL",
    "NULL_CUE_UNCERTAINTY_KCAL",
    "RefutedEnergy",
    "STANDARD_TEMPERATURE_K",
    "ThermodynamicData",
    "dg_prime",
    "electrical_work",
    "reaction_dg0",
    "reaction_dg0_uncertainty",
    "reaction_energy",
]

FARADAY_KJ_PER_V = 96.485

KJ_PER_KCAL = 4.184
"""The vendored thermodynamic database is in kcal/mol, and it says so.

`thermo_data.thermodb`'s own ``units`` field reads ``kcal/mol``, and its ``cpd00001``
(water) carries ``deltaGf_std = -56.687``, which is -237.2 kJ/mol -- the standard value.
So every number read out of that file's ``deltaGf_err`` and out of its structural-cue table
is kcal/mol and is multiplied by this before being reported.

`docs/research/CAROTENOID_ENERGIES.md` records an analysis that was refuted three votes to
zero for missing exactly this, which is why the check is written down rather than assumed:
``pytfa``'s :class:`~pytfa.thermo.metabolite.MetaboliteThermo` DEFAULTS to
``thermo_unit='kJ/mol'``, so a caller who does not pass the unit gets the pH and ionic
strength transform in kJ on top of a tabulated energy in kcal -- a mixed scale that is
neither unit. :func:`_load` therefore passes ``thermo_unit='kcal/mol'``, which is the
parameter pytfa provides for exactly this, and converts once with this constant.

What settles it is a value outside either estimator. Transformed to pH 7 and I = 0.25 M,
water's standard formation energy is -155.66 kJ/mol (Alberty, *Thermodynamics of
Biochemical Reactions*). Reading the table in its own unit and converting gives -155.70.
Not passing the unit gives **+24.85**, and no scale on which the formation energy of water
is positive is a scale at all.
"""

NULL_CUE_UNCERTAINTY_KCAL = 2.0
"""Uncertainty to report when a reaction's net structural cues cancel exactly, kcal/mol.

``pytfa``'s own ``ThermoModel.prepare(null_error_override=2)`` uses this value and its
docstring gives the reason: "2kcal/mol is standard in estimation frameworks like GCM". A
net cue vector that cancels to zero means the group decomposition sees no difference
between the two sides, which is a statement about the decomposition's resolution and not a
claim that the energy is known exactly. Reporting 0 would let an opt-in uncertainty
threshold pass anything, which is the failure this constant exists to prevent."""

FREE_MAGNESIUM_M = 1.0e-3
"""Free cytosolic Mg2+. ATP, ADP and phosphate all bind it, and at this concentration most
ATP is present as MgATP. Binding stabilises the more highly charged species preferentially,
which is why the textbook -30.5 kJ/mol for ATP hydrolysis is quoted at 1 mM Mg and the same
reaction without it sits nearer -37."""

_LOG_K_MAGNESIUM = {"cpd00002": 4.19, "cpd00008": 3.17, "cpd00009": 1.88}
"""Mg2+ association constants for ATP, ADP and phosphate, from Alberty's tabulations at
I = 0.25 M. Only these three matter at the precision anything here is used to."""


@dataclass(frozen=True)
class Compartment:
    """Where a metabolite sits, and the two things that changes about its energy."""

    ph: float
    potential_v: float
    source: str


COMPARTMENTS = {
    "c": Compartment(7.5, 0.0,
                     "cytosol, the electrical reference. pH measured by pHluorin alongside "
                     "roGFP2 in Ayer 2013, PMID 23762325"),
    "m": Compartment(7.8, -0.16,
                     "mitochondrial matrix, alkaline and electrically negative to the "
                     "cytosol; the inner membrane potential drives oxidative phosphorylation"),
    "e": Compartment(5.0, 0.20,
                     "medium at the pH these plates are run at; the plasma membrane holds "
                     "the cytosol negative to it"),
    "v": Compartment(6.2, 0.03,
                     "vacuole, acidified by the V-ATPase and slightly positive inside"),
    "n": Compartment(7.5, 0.0, "nucleus, continuous with the cytosol through the pores"),
    "p": Compartment(8.2, 0.0, "peroxisome, alkaline and without a maintained potential"),
    "er": Compartment(7.2, 0.0, "endoplasmic reticulum"),
    "g": Compartment(6.8, 0.0, "Golgi, progressively acidified along the stack"),
}
"""Compartment pH and electrical potential relative to the cytosol.

Both terms were missing. A formation energy transformed at one pH puts every other
compartment at the wrong protonation, and a charge crossing a membrane does electrical work
of zF*dPsi whatever the chemistry does -- without which proton pumping looks impossible.
"""

GAS_CONSTANT_KJ = 8.314e-3
STANDARD_TEMPERATURE_K = 303.15
_NO_VALUE = 10000000


def _table(name: str) -> pathlib.Path:
    """One of the vendored thermodynamic tables, resolved through :mod:`ystwin.paths`.

    Not a module-level constant: the tables are a third-party download that is not
    tracked, so where they are is a question with a per-machine answer, and asking it at
    call time is what lets the answer be "nowhere, set YSTWIN_THERMO".
    """
    return paths.require(paths.thermo_dir(), "the thermodynamic tables (data/thermo)",
                         "YSTWIN_THERMO") / name


def _authored_table(name: str) -> pathlib.Path:
    """One of the tables this repository AUTHORS *about* the vendored ones.

    Resolved against the checkout rather than through :func:`ystwin.paths.thermo_dir`, and
    the difference is the point. ``aliases.tsv`` and ``thermo_data.thermodb`` are a
    third-party download that a machine may keep anywhere, which is what ``YSTWIN_THERMO``
    is for. ``heterologous_metabolites.tsv`` and ``refuted_energies.tsv`` are this
    repository's own claims -- which species an id installed by `fba/carotenoid.py` is, and
    which tabulated energies must not be used -- and they ship with the source. Resolving
    them through the override would mean a machine that points ``YSTWIN_THERMO`` at its own
    ModelSEED copy silently loses a refusal, which is the one direction a refusal must
    never fail in.
    """
    return paths.REPO_ROOT / "data" / "thermo" / name


CYTOSOLIC_PH = 7.5
"""Measured, not assumed. Ayer 2013 (PMID 23762325) reports 7.5 in the yeast cytosol using
pHluorin alongside roGFP2. Much of the 40 mV disagreement between the two E_GSH literatures
comes from one line assuming 7.0 and the other measuring 7.5, so it is a parameter here
rather than a constant buried in a transform."""

IONIC_STRENGTH_M = 0.25


@dataclass(frozen=True)
class RefutedEnergy:
    """A tabulated formation energy that resolves and must not be used anyway.

    Read from ``data/thermo/refuted_energies.tsv``, one row per ModelSEED compound. It
    carries no corrected value and there is no field for one: a row records that a
    disagreement was measured, not what the right number is.

    Args:
        modelseed_id: The compound whose energy is refused.
        name: Its name, for the report.
        reason: Why the quoted uncertainty is not believed, in the words of the row.
        source: Where that was established.
    """

    modelseed_id: str
    name: str
    reason: str
    source: str


@dataclass(frozen=True)
class ThermodynamicData:
    """Standard formation energies for the metabolites of a yeast model.

    Keyed by the model's own metabolite ids, resolved through the ModelSEED alias tables.
    Coverage is partial and the number is reported rather than hidden: a reaction needs
    every one of its metabolites to carry an energy before it can be given a constraint.

    Energies are transformed to the given pH and ionic strength, not the raw tabulated
    values. Untransformed, ATP hydrolysis comes out at -6.8 kJ/mol against a textbook -30;
    transformed it is -39.5, which is what the textbook figure becomes once the magnesium
    it assumes is taken away. Magnesium binding is NOT modelled here, and it matters most
    for exactly the adenylates.
    """

    by_metabolite: dict[str, float] = field(repr=False)
    coverage: float
    ph: float = CYTOSOLIC_PH
    ionic_strength: float = IONIC_STRENGTH_M
    rejected_on_formula: int = 0
    free_magnesium_m: float = FREE_MAGNESIUM_M
    seed_by_metabolite: dict[str, str] = field(default_factory=dict, repr=False)
    """Which ModelSEED compound each metabolite id resolved to. The audit trail for
    everything below: a refutation and a cue decomposition are both properties of the
    COMPOUND, and without this a caller cannot tell which one an energy came from."""
    cues_by_metabolite: dict[str, dict[str, float]] = field(default_factory=dict,
                                                            repr=False)
    """Structural-cue counts per metabolite, from the same database entry as the energy.
    Empty for a metabolite whose entry carries no decomposition."""
    cue_uncertainty_kcal: dict[str, float] = field(default_factory=dict,
                                                   repr=False)
    """Per-cue standard error, kcal/mol, from the database's own ``cues`` table. Kept in the
    table's units and converted at the point of use, so there is one place to look."""
    refuted: dict[str, RefutedEnergy] = field(default_factory=dict, repr=False)
    """Metabolite id -> the row of ``refuted_energies.tsv`` that refuses its energy. The
    energy is still in ``by_metabolite``: withholding it there would make a refuted energy
    indistinguishable from an absent one, which is the distinction this exists to keep."""
    declared: tuple[str, ...] = ()
    """Metabolite ids that came from ``heterologous_metabolites.tsv`` rather than from the
    model -- species this repository installs, which are in no published alias table."""

    def dgf(self, metabolite_id: str) -> float | None:
        """Standard formation energy, or None where the tables do not reach.

        **On the scale.** kJ/mol. These are ``pytfa``'s ``deltaGf_tr`` read under
        ``thermo_unit='kcal/mol'`` -- the unit the database declares -- and converted once
        by :data:`KJ_PER_KCAL`.

        This was wrong until 2026-08-30 and the error was a unit, not an estimate. Read
        under pytfa's default the transform lands in kJ on a kcal table, and water comes
        out at +24.85; read correctly it is -155.70 against Alberty's -155.66. The
        correction was not adopted for that reason alone. Scored against eQuilibrator's
        component contribution over the 1,098 single-compartment yeast-GEM reactions all
        three cover, the same reactions read the two ways give

            mixed scale    RMSE 320.1 kJ/mol   r 0.435   median |error| 35.9
            kcal -> kJ     RMSE 116.6 kJ/mol   r 0.942   median |error| 25.2

        so an independent estimator agrees with the corrected reading and not the mixed
        one. `tests/test_thermodynamic_physics.py` recomputes this rather than trusting the prose.

        Every committed thermodynamic number changed with it, which is why the change is
        recorded here and in `docs/research/CAROTENOID_ENERGIES.md` rather than folded in
        quietly. :func:`reaction_dg0_uncertainty` did not move -- it was always computed
        from the cue table in kcal/mol and converted once -- and a :class:`RefutedEnergy`
        carries no number to move.
        """
        return self.by_metabolite.get(metabolite_id)

    def cues(self, metabolite_id: str) -> dict[str, float] | None:
        """Cue counts for one metabolite, or None where there is no decomposition."""
        return self.cues_by_metabolite.get(metabolite_id) or None

    def refutation(self, metabolite_id: str) -> RefutedEnergy | None:
        """The row refusing this metabolite's energy, or None if none refuses it."""
        return self.refuted.get(metabolite_id)

    @classmethod
    def load(cls, model=None, aliases=None, thermodb=None, ph: float = CYTOSOLIC_PH,
             ionic_strength: float = IONIC_STRENGTH_M,
             free_magnesium_m: float = FREE_MAGNESIUM_M,
             heterologous=None, refuted=None) -> ThermodynamicData:
        """Build the table by mapping model metabolites onto SEED compounds.

        Args:
            model: A cobra model; defaults to the project's yeast GEM.
            aliases: ModelSEED alias table; defaults to the copy under `data/thermo`.
            thermodb: pytfa thermo database; defaults to the same directory.
            ph: Compartment pH the energies are transformed to.
            ionic_strength: Ionic strength in molar.
            free_magnesium_m: Free Mg2+, which ATP, ADP and phosphate bind. Zero switches
                the correction off, which is how the size of it gets checked.
            heterologous: Declarations for metabolites this repository installs into a host
                model, which are in no published alias table; defaults to
                ``data/thermo/heterologous_metabolites.tsv``. These are added whether or not
                ``model`` carries them, because the gate is asked about a pathway before any
                model has it installed.
            refuted: Energies that resolve and must not be used; defaults to
                ``data/thermo/refuted_energies.tsv``.
        """
        return _load(model, aliases, thermodb, ph, ionic_strength, free_magnesium_m,
                     heterologous, refuted)


def _load(model=None, aliases=None, thermodb=None, ph=CYTOSOLIC_PH,
          ionic_strength=IONIC_STRENGTH_M,
          free_magnesium_m=FREE_MAGNESIUM_M,
          heterologous=None, refuted=None) -> ThermodynamicData:
    import warnings

    import pandas as pd
    from pytfa.io import load_thermoDB
    from pytfa.thermo.metabolite import MetaboliteThermo

    if model is None:
        model = _default_model()
    alias_table = pd.read_csv(aliases or _table("aliases.tsv"), sep="\t")
    # The bundled database unpickles through a numpy call deprecated in 2.4.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        loaded = load_thermoDB(str(thermodb or _table("thermo_data.thermodb")))
    database = loaded["metabolites"]
    cue_errors = {str(name): float(entry["error"]) for name, entry in loaded["cues"].items()}

    formulas = {k: v.get("formula") for k, v in database.items()}
    usable = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for key, entry in database.items():
            if entry.get("deltaGf_std") in (None, _NO_VALUE):
                continue
            # pytfa 0.9.1 wipes __dict__ in __init__, so read deltaGf_tr not calcDGis().
            # thermo_unit is the database's, not the caller's -- see KJ_PER_KCAL. Without
            # it the pH and ionic-strength transform lands in kJ on a kcal table and the
            # result is neither unit.
            value = MetaboliteThermo(entry, pH=ph, ionicStr=ionic_strength,
                                     temperature=STANDARD_TEMPERATURE_K,
                                     thermo_unit="kcal/mol").deltaGf_tr
            if value is None or abs(value) >= _NO_VALUE:
                continue
            value *= KJ_PER_KCAL
            log_k = _LOG_K_MAGNESIUM.get(key)
            if log_k is not None and free_magnesium_m > 0:
                bound = 10.0**log_k * free_magnesium_m
                value -= GAS_CONSTANT_KJ * STANDARD_TEMPERATURE_K * np.log1p(bound)
            usable[key] = float(value)
    routes = {"kegg.compound": "KEGG", "bigg.metabolite": "BiGG"}
    lookup: dict[str, dict[str, list]] = {}
    for key, source in routes.items():
        table: dict[str, list] = {}
        rows = alias_table[alias_table.Source == source]
        for external, seed in zip(rows["External ID"], rows["ModelSEED ID"]):
            table.setdefault(external, []).append(seed)
        lookup[key] = table

    energies, rejected, seeds = {}, 0, {}
    for met in model.metabolites:
        for key, table in lookup.items():
            value = met.annotation.get(key)
            value = value[0] if isinstance(value, list) else value
            candidates = [c for c in table.get(value, []) if c in usable] if value else []
            matched = next((c for c in candidates
                            if _same_skeleton(met.formula, formulas.get(c))), None)
            if candidates and matched is None:
                rejected += 1
            if matched:
                energies[met.id] = float(usable[matched])
                seeds[met.id] = matched
                break
    # Coverage is over the MODEL's metabolites and is computed before the declared species
    # are merged in, because they are not in the model's denominator and adding them would
    # report a reach the tables do not have.
    coverage = len(energies) / max(len(model.metabolites), 1)

    declared = _declared_metabolites(heterologous, usable, formulas)
    for metabolite_id, seed in declared.items():
        energies.setdefault(metabolite_id, float(usable[seed]))
        seeds.setdefault(metabolite_id, seed)

    refusals = _refuted_energies(refuted)
    cues = {metabolite_id: {str(c): float(n)
                            for c, n in (database[seed].get("struct_cues") or {}).items()}
            for metabolite_id, seed in seeds.items()}
    return ThermodynamicData(
        energies, coverage, ph, ionic_strength, rejected, free_magnesium_m,
        seed_by_metabolite=seeds, cues_by_metabolite=cues,
        cue_uncertainty_kcal=cue_errors,
        refuted={metabolite_id: refusals[seed] for metabolite_id, seed in seeds.items()
                 if seed in refusals},
        declared=tuple(sorted(declared)))


def declared_seed_ids(heterologous=None) -> dict[str, str]:
    """``metabolite_id -> ModelSEED id`` from ``heterologous_metabolites.tsv``, unvalidated.

    The mapping alone, without the database cross-checks :func:`_declared_metabolites`
    performs, because a caller may need to know WHICH compound an installed id is without
    having loaded a thermodynamic database at all -- `bridge/equilibrator.py` asking whether
    a refutation applies to it, for instance. The validation still happens on every real
    load; this is the lookup, not a second way in.
    """
    import pandas as pd

    path = heterologous if heterologous is not None else _authored_table(
        "heterologous_metabolites.tsv")
    table = pd.read_csv(path, sep="\t", comment="#")
    return {str(row.metabolite_id): str(row.modelseed_id)
            for row in table.itertuples(index=False)}


def _declared_metabolites(heterologous, usable, formulas) -> dict[str, str]:
    """``metabolite id -> ModelSEED id`` for the species this repository installs.

    The naming gap `docs/research/CAROTENOID_ENERGIES.md` names, closed as data. The alias
    tables map ModelSEED ids onto ids from PUBLISHED models, so an id `fba/carotenoid.py`
    invents when it installs a pathway resolves to nothing and never will. Declaring it is
    the fix, and the declaration is checked rather than trusted: the row's formula must
    agree with the database entry's, which is the same guard `_same_skeleton` already
    applies to an alias hit.

    Raises:
        ValueError: naming a row whose ModelSEED id is absent from the database or carries
            no usable energy, and a row whose declared formula disagrees with it. Either is
            a typo in a hand-written table, and a typo here silently attaches one molecule's
            energy to another molecule's name.
    """
    import pandas as pd

    path = heterologous if heterologous is not None else _authored_table(
        "heterologous_metabolites.tsv")
    table = pd.read_csv(path, sep="\t", comment="#")
    out: dict[str, str] = {}
    for row in table.itertuples(index=False):
        seed = str(row.modelseed_id)
        if seed not in usable:
            raise ValueError(
                f"{path} declares "
                f"{row.metabolite_id!r} to be {seed!r}, which the thermodynamic database "
                "does not carry with a usable formation energy. A declaration that "
                "resolves to nothing is worse than no declaration: it reads as a closed "
                "gap. Check the id against data/thermo/compounds.tsv")
        if not _same_skeleton(str(row.formula), formulas.get(seed)):
            raise ValueError(
                f"declared formula {row.formula!r} for {row.metabolite_id!r} disagrees "
                f"with {seed!r}'s {formulas.get(seed)!r}. The row names one molecule and "
                "the id names another; supplying the energy anyway would attach the "
                "second's number to the first's name")
        out[str(row.metabolite_id)] = seed
    return out


def _refuted_energies(refuted) -> dict[str, RefutedEnergy]:
    """``ModelSEED id -> RefutedEnergy`` for the energies that resolve and must not be used."""
    import pandas as pd

    path = refuted if refuted is not None else _authored_table("refuted_energies.tsv")
    table = pd.read_csv(path, sep="\t", comment="#")
    return {str(row.modelseed_id): RefutedEnergy(str(row.modelseed_id), str(row.name),
                                                 str(row.reason), str(row.source))
            for row in table.itertuples(index=False)}


@lru_cache(maxsize=1)
def _default_model():
    """Yeast9, wherever this machine keeps it.

    This was an absolute path into one person's home directory, which meant every
    thermodynamic number in the README was reproducible on exactly one laptop. Resolved
    now, and a refusal that names ``YSTWIN_YEAST_GEM`` if it is nowhere -- callers that
    should skip instead of failing check :func:`ystwin.paths.yeast_gem` for ``None``.

    Raises:
        FileNotFoundError: if the model cannot be located.
    """
    import warnings


    model = paths.require(paths.yeast_gem(), "the Yeast9 GSMM (yeast-GEM.xml)",
                          "YSTWIN_YEAST_GEM")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Pinned loader, not cobra.io directly: a different simplex breaks ties
        # among alternate optima differently. See fba/solver.py.
        from ..fba.solver import load_model

        return load_model(model)[0]


def _elements(formula: str | None) -> dict[str, int] | None:
    """Element counts, hydrogen dropped so protonation state does not matter."""
    import re

    if not formula:
        return None
    counts: dict[str, int] = {}
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not element or element == "H":
            continue
        counts[element] = counts.get(element, 0) + (int(number) if number else 1)
    return counts


def _same_skeleton(model_formula, db_formula) -> bool:
    """Whether two formulas describe the same molecule up to protonation.

    An identifier can map to several database compounds and the alias tables give no way to
    choose. Without this check ribose-5-phosphate resolved to a chlorinated herbicide, and
    its formation energy went into every reaction the pentose phosphate pathway needs.
    """
    left, right = _elements(model_formula), _elements(db_formula)
    return left is not None and right is not None and left == right


def electrical_work(charge_moves: dict) -> float:
    """Work done moving charge between compartments, in kJ/mol.

    Keyed by ``(compartment, charge)`` with the reaction's coefficient as the value. A
    charge crossing a membrane does zF*dPsi of work whatever the chemistry does, and
    omitting it makes proton pumping look impossible.
    """
    total = 0.0
    for (compartment, charge), coefficient in charge_moves.items():
        potential = COMPARTMENTS[compartment].potential_v if compartment in COMPARTMENTS else 0.0
        total += coefficient * charge * FARADAY_KJ_PER_V * potential
    return float(total)


def reaction_dg0(stoichiometry: dict[str, float], thermo: ThermodynamicData) -> float | None:
    """Standard reaction energy, or None if any participant is uncovered.

    Returning None rather than a partial sum is the point: an energy computed from some of
    a reaction's metabolites is not a smaller energy, it is a different reaction's.
    """
    total = 0.0
    for metabolite, coefficient in stoichiometry.items():
        energy = thermo.dgf(metabolite)
        if energy is None:
            return None
        total += coefficient * energy
    return total


def reaction_dg0_uncertainty(stoichiometry: dict[str, float],
                             thermo: ThermodynamicData) -> float | None:
    """Standard error on a reaction's standard energy, kJ/mol, or None where it is unknown.

    **Propagated over the NET structural cues, not over the formation energies, and the
    difference is not a detail.** Group contribution estimates a formation energy by adding
    up the groups a molecule is made of, so two molecules that share a scaffold share most
    of their error, and quadrature over their individual errors counts that shared part
    twice. yeast-GEM's cytosolic thiolase is the case: acetyl-CoA, acetoacetyl-CoA and CoA
    carry tabulated formation errors of 14.8, 15.0 and 14.4 kJ/mol, which in quadrature
    give 36.2 kJ/mol for a reaction whose energy is +38.07 -- produced entirely by
    triple-counting one CoA moiety, and larger than the effect was thought to be when this
    was written, the energy then reading +15.6. Cancelling the cues first
    leaves a net vector of five groups and 4.5 kJ/mol, and 4.5 is the number that means
    something. This is what ``pytfa`` itself does (``calcDGR_cues``, whose result becomes
    ``deltaGRerr`` on every reaction of a thermodynamic model) and what Henry 2007's group
    contribution method prescribes.

    The cue table is in kcal/mol -- see :data:`KJ_PER_KCAL` -- and is converted here, which
    is why this number is in kJ/mol even though :meth:`ThermodynamicData.dgf`'s are not.

    Returns:
        The standard error in kJ/mol; :data:`NULL_CUE_UNCERTAINTY_KCAL` converted where the
        net cues cancel exactly; ``None`` where any participant has no cue decomposition,
        for the same reason :func:`reaction_dg0` returns None on an uncovered metabolite --
        an uncertainty summed over some of a reaction is a different reaction's.
    """
    net: dict[str, float] = {}
    for metabolite, coefficient in stoichiometry.items():
        cues = thermo.cues(metabolite)
        if cues is None:
            return None
        for cue, count in cues.items():
            net[cue] = net.get(cue, 0.0) + coefficient * count
    variance = 0.0
    for cue, count in net.items():
        if count == 0:
            continue
        error = thermo.cue_uncertainty_kcal.get(cue)
        if error is None:
            return None
        variance += (count * error) ** 2
    if variance == 0.0:
        return float(NULL_CUE_UNCERTAINTY_KCAL * KJ_PER_KCAL)
    return float(np.sqrt(variance) * KJ_PER_KCAL)


def dg_prime(
    dg0: float, stoichiometry: dict[str, float], concentrations: dict[str, float],
    temperature_k: float = STANDARD_TEMPERATURE_K,
) -> float:
    """Reaction energy under given concentrations: ``dG0 + RT ln Q``.

    This is where a measured ratio enters. A cofactor pair appears in Q with opposite
    signs, so their ratio moves the energy directly -- and where that crosses zero it
    decides which way the reaction may run.

    Args:
        dg0: Standard reaction energy, kJ/mol.
        stoichiometry: Coefficients keyed by metabolite id, negative for substrates.
        concentrations: Molar concentrations for every participant with a coefficient.
        temperature_k: Absolute temperature.
    """
    log_q = 0.0
    for metabolite, coefficient in stoichiometry.items():
        if coefficient == 0:
            continue
        if metabolite not in concentrations:
            raise KeyError(f"no concentration given for {metabolite!r}, so Q is undefined")
        value = concentrations[metabolite]
        if value <= 0:
            raise ValueError(f"concentration of {metabolite!r} must be positive, got {value}")
        log_q += coefficient * np.log(value)
    return float(dg0 + GAS_CONSTANT_KJ * temperature_k * log_q)


def reaction_energy(reaction, thermo: ThermodynamicData, split: bool = False):
    """Standard energy of a cobra reaction, chemistry plus electrical work.

    The electrical term is not optional for anything crossing a membrane. A proton entering
    the yeast mitochondrial matrix releases about 15 kJ/mol of it, and a chemistry-only
    energy makes proton pumping look impossible.

    Args:
        reaction: A cobra reaction.
        thermo: Formation energies.
        split: Return the two parts separately rather than their sum.
    """
    stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
    chemistry = reaction_dg0(stoichiometry, thermo)
    if chemistry is None:
        return None
    compartments = {m.compartment for m in reaction.metabolites}
    moves: dict = {}
    for metabolite, coefficient in reaction.metabolites.items():
        charge = metabolite.charge or 0
        if charge:
            key = (metabolite.compartment, charge)
            moves[key] = moves.get(key, 0.0) + coefficient

    if len(compartments) > 1 and abs(_net_charge_moved(moves)) > _MAX_NET_CHARGE:
        return None

    electrical = electrical_work(moves)
    return (chemistry, electrical) if split else chemistry + electrical


def _net_charge_moved(moves: dict) -> float:
    """Charge leaving the cytosol, summed over compartments."""
    return sum(coefficient * charge
               for (compartment, charge), coefficient in moves.items()
               if compartment != "c")


_MAX_NET_CHARGE = 1.0
"""Net charge a carrier plausibly translocates in one turn.

Real electrogenic carriers move about one -- the ADP/ATP exchanger swaps ATP(4-) out for
ADP(3-) in, net one. Where a model writes a bare uniport of a highly charged species the
counter-ion is missing from the equation, and charging the full zF*dPsi gives +77 kJ/mol for
PRPP transport and +124 for propionyl-CoA, which no carrier does. Those energies cannot be
estimated from the stoichiometry given, so they are refused rather than fabricated."""
