"""A heterologous pathway declared as data, so a new product is a file and not a diff.

Everything about beta-carotene that used to live in Python lives here as a `PathwaySpec`
loaded from `data/pathways/<product>.toml`: the chain of intracellular species, which step
each enzyme catalyses, how the flux splits at a branch, what the terminal product weighs,
and which native metabolite the pathway is pulled off.

**Scope, stated once and enforced.** The solver walks intracellular pools in one chain.
A DILUTED node loses material only to the next step and growth, so its steady-state balance
is:

    d[X]/dt = v_in - v_out - mu*[X] = 0

A SECRETED or DEGRADED intermediate breaks that walk and is refused at load time. Either
fate is accepted on a TERMINAL node, where the solver requires measured loss kinetics
rather than silently treating the pool as diluted only. Ferreira 2018 raises
triacylglycerol from 129 to 218 mg/gDCW by deleting the lipases alone, which is how large
the omitted term can be.

For those terminal fates the gross entry flux is mu times the sum of the pools PLUS the
terminal loss. Pool measurements alone therefore do not identify gross synthesis.

**What a spec must supply and what it must not.** It supplies chemistry and stoichiometry,
which are properties of the pathway. It does not supply fitted constants: those live in a
calibration table beside a citation, because a fitted number without the data behind it is
the thing this repository exists to avoid.

**Whether a pathway can be calibrated is a property of its chain, and it is announced at
load.** :class:`Calibratability` works it out from the rate laws and says so on the way
past, because the alternative -- a boolean property nobody reads -- is what let `phb` sit
in `data/pathways/` looking fittable. It is a record and not a refusal: a spec that cannot
be calibrated still exercises the solver and the environment layer, which is exactly what
`phb` is for, and refusing it would delete a working capability to prevent a mistake
nobody has made yet.
"""

from __future__ import annotations

import logging
import pathlib
import tomllib
from dataclasses import dataclass

__all__ = [
    "Calibratability",
    "Fate",
    "Node",
    "PathwaySpec",
    "RateLaw",
    "available_pathways",
    "load_pathway",
    "pathways_dir",
]

# Deliberately a log record and not ``warnings.warn``. `pyproject.toml` sets
# ``filterwarnings = ["error"]``, so a warning raised here would make `load_pathway("phb")`
# fail in every test that touches it -- turning a note into the refusal this must not be.
# A log record reaches stderr through `logging.lastResort` on any unconfigured run, which
# is where a person loading a spec is looking, and stays inert under a test that is not.
_LOG = logging.getLogger(__name__)


class RateLaw:
    """How the flux out of a node depends on that node's content.

    ``SATURATING`` is Michaelis-Menten, ``v = vmax * X / (km + X)``, and is the only law
    with a closed-form steady state for a single node, so it is the default.

    ``PROPORTIONAL`` is ``v = k * X``, first order. It is NOT a special case of the above
    that the solver can reach by taking km large: it has no ceiling, so a node under it
    cannot accumulate without bound and cannot saturate. Declaring it is a claim that the
    enzyme is far from saturation, which is checkable.

    ``PASSTHROUGH`` is ``v_out = v_in`` -- the node holds no appreciable pool. It exists so
    a spec can say "this intermediate does not accumulate" explicitly instead of leaving the
    reader to infer it from a missing node, which is the assumption that made
    ``q_lycopene + q_betacarotene`` get called "the pathway flux" without anyone checking
    whether phytoene accumulates.
    """

    SATURATING = "saturating"
    PROPORTIONAL = "proportional"
    PASSTHROUGH = "passthrough"
    ALL = (SATURATING, PROPORTIONAL, PASSTHROUGH)


class Fate:
    """What happens to a node's pool besides being consumed by the next step.

    ``DILUTED`` -- removed by growth dilution only, with no additional terminal loss.

    ``SECRETED`` -- leaves the cell through a transporter. On a TERMINAL node this is
    accepted and refused at SOLVE time for want of an export capacity, because the balance
    closes on ``mu*[X]`` and a secreted pool leaves at a rate the cell sets rather than one
    growth sets, so the missing thing is a measured constant rather than the chemistry. On an
    INTERMEDIATE it is refused HERE, and no measurement lifts that: the walk hands one flux
    from each node to the next, and a node leaking sideways breaks the walk itself.

    ``DEGRADED`` -- consumed by a native enzyme the pathway does not include. Same split, same
    reason: refused at solve time on a terminal node, at load time on an intermediate.

    This docstring said both were "**Refused**" flatly until 2026-09-02, three days after
    `solve_pathway` began handling a terminal node of either fate -- a stale claim of exactly
    the kind `scripts/audit_claims.py` exists to catch, in the one place that gate does not
    look, because it pins numbers and this was a word.

    The native storage compounds are why ``DEGRADED`` exists. Glycogen has GPH1 and SGA1,
    trehalose has NTH1/NTH2/ATH1, triacylglycerol has TGL3/4/5 -- and Ferreira 2018 gets its
    single largest increase, 129 to 218 mg/gDCW, by DELETING the lipases, which is direct
    experimental proof that the terminal node has a large enzymatic outlet. A spec that
    declared one of those as ``DILUTED`` would return ``content = flux/mu`` and silently
    omit the term doing most of the work.

    Naming the fate is what turns that from a wrong number into a refusal. `data/pathways/
    glycogen.toml` documents its own GPH1 and SGA1 reactions and is refused on the strength
    of its own notes.
    """

    DILUTED = "diluted"
    SECRETED = "secreted"
    DEGRADED = "degraded"
    ALL = (DILUTED, SECRETED, DEGRADED)


@dataclass(frozen=True)
class Node:
    """One intracellular species in the chain.

    Args:
        name: Species name, used in results and in error messages.
        rate_law: How flux leaves it. One of :class:`RateLaw`.
        fate: What else removes it. One of :class:`Fate`.
        enzyme: Gene name of the enzyme consuming it, or ``""`` for the terminal node.
        molar_mass_g_per_mol: For reporting in mg/gDCW. Required on the terminal node,
            optional elsewhere.
        measurable: Whether this species is routinely quantified. What it costs to be
            false depends on ``rate_law`` and not on how many other nodes are measurable:
            an unmeasured node that HOLDS A POOL is a free parameter in the entry flux,
            and an unmeasured ``passthrough`` node is a declared zero nobody can check.
            :class:`Calibratability` is what tells the two apart.
        metabolite: GEM metabolite id this node IS, or ``""`` where none is declared.
            Chemistry, so it belongs here: which species a node is, is a property of the
            pathway. Required before this node's step can be gated, because
            :mod:`ystwin.pathway.thermo_gate` uses it twice -- to route the solved content
            into Q, and to check the step is written in the direction that MAKES this node
            rather than the one that consumes it.
        reaction: GEM reaction id that PRODUCES this node, or ``""``. Keyed by the produced
            node rather than the consumed one because the entry step's substrate is the
            precursor, which is not a node in the chain.
        stoichiometry: ``{metabolite id: coefficient}`` for that same producing step,
            negative for substrates. Declaring it as well as ``reaction`` is not
            redundancy: it is what lets the gate run **without an SBML parse**, and where
            both are given the gate refuses a pair that describes two different reactions.
            A heterologous step has no GEM reaction until one is installed, so for those
            this is the only route.
    """

    name: str
    rate_law: str = RateLaw.SATURATING
    fate: str = Fate.DILUTED
    enzyme: str = ""
    molar_mass_g_per_mol: float | None = None
    measurable: bool = True
    metabolite: str = ""
    reaction: str = ""
    stoichiometry: dict[str, float] | None = None

    def __post_init__(self) -> None:
        if self.rate_law not in RateLaw.ALL:
            raise ValueError(
                f"node {self.name!r}: rate_law must be one of {RateLaw.ALL}, "
                f"got {self.rate_law!r}")
        if self.fate not in Fate.ALL:
            raise ValueError(
                f"node {self.name!r}: fate must be one of {Fate.ALL}, got {self.fate!r}")
        if self.molar_mass_g_per_mol is not None and self.molar_mass_g_per_mol <= 0:
            raise ValueError(
                f"node {self.name!r}: molar mass must be positive, got "
                f"{self.molar_mass_g_per_mol}")
        if (self.reaction or self.stoichiometry) and not self.metabolite:
            raise ValueError(
                f"node {self.name!r} declares a producing step but no `metabolite`. The "
                "gate needs the node's own species id to check the step is written in the "
                "direction that makes it: yeast-GEM writes the thiolase r_0103 reversibly, "
                "and as written it is +38.07 kJ/mol while the reverse is -38.07. Without the "
                "id nothing downstream can tell those apart, so the step would be reported "
                "under the forward name with the reverse verdict")
        if self.stoichiometry is not None:
            if not self.stoichiometry:
                raise ValueError(
                    f"node {self.name!r}: `stoichiometry` is empty. Omit it rather than "
                    "declaring a reaction with no participants")
            coefficient = self.stoichiometry.get(self.metabolite)
            if coefficient is None:
                raise ValueError(
                    f"node {self.name!r} is declared to be {self.metabolite!r}, which does "
                    f"not appear in its own step's stoichiometry "
                    f"({sorted(self.stoichiometry)}). A metabolite id that matches nothing "
                    "is a typo, and this one leaves the node's solved content out of Q "
                    "entirely -- the gate would then report a verdict about whatever "
                    "background level was supplied while saying it tested the solve's")
            if coefficient <= 0:
                raise ValueError(
                    f"node {self.name!r} appears with coefficient {coefficient:g} in the "
                    "step declared to PRODUCE it, which is the direction that consumes it. "
                    "dG0's sign is the sign of the direction as written, so gating this "
                    "would report the reverse reaction's verdict under the forward one's "
                    "name. Reverse the coefficients")


def _hidden_pool_mechanism(node: Node, is_terminal: bool) -> str:
    """What one unmeasured pool buys, written in the parameters THAT node actually has.

    One sentence for all three cases was wrong twice over. A saturating node hides a
    ``(vmax, km)`` ridge; a proportional node hides a single rate constant and no ridge at
    all; and a DILUTED product hides no kinetic parameter, because ``solve_pathway`` gives
    it ``v_in/mu`` regardless of its rate law. A terminal loss instead needs its own
    measured kinetics. Naming a next-step enzyme or declaring the product passthrough
    cannot eliminate either its pool or its declared loss outlet. Naming the mechanism
    is the whole point of the record.
    """
    if is_terminal and node.fate in (Fate.DEGRADED, Fate.SECRETED):
        return (
            f"{node.name!r} is the product and its pool is unmeasured. Its "
            f"{node.fate.upper()} fate also requires an independently measured terminal "
            "loss; declaring the product 'passthrough' removes neither the pool nor that "
            "outlet")
    if is_terminal:
        return (
            f"{node.name!r} is the product, and it holds no kinetic parameter to blame: the "
            "solver gives a terminal node content = v_in/mu whatever rate law it declares, "
            "so this is not a ridge but a chain with nothing measured at its end, where "
            "every content returned is proportional to an entry flux no measurement "
            "touches. Declaring the product 'passthrough' does not lift it, because the "
            "solver ignores that declaration on the last node")
    if node.rate_law == RateLaw.PROPORTIONAL:
        return (
            f"{node.name!r} is proportional, so its pool is v/k and the entry flux carries a "
            "factor of (1 + mu/k). The free constant is one rate constant and not a "
            "(vmax, km) pair, and nothing downstream bounds it: as k falls the implied entry "
            "flux rises without limit while every measured content stays put")
    return (
        f"{node.name!r} is saturating, so for any (vmax, km) whose vmax exceeds the flux it "
        "must carry, a pool of km*v/(vmax - v) delivers that flux exactly. Parameter sets "
        "orders of magnitude apart then reproduce every measured content to machine "
        "precision while implying entry fluxes hundreds of times apart -- see "
        "docs/research/KINETIC_FIT.md section 3, which reports exactly this ridge for the "
        "desaturase step")


@dataclass(frozen=True)
class Calibratability:
    """Whether a pathway's flux scalar can be fitted from its own measurements, and why not.

    This replaces a count. The property it replaces returned ``True`` when at least two
    nodes were measurable, and a count is the wrong test in BOTH directions:

    * **It says yes where the scalar is not identified.** For a DILUTED chain the entry
      flux is ``mu`` times every pool, so an unmeasured pool is a free parameter. A chain
      ``A -> B -> C`` with A and C measured and B a saturating node nobody measures has two
      measurable nodes and no identified entry flux: for any ``(vmax, km)`` with vmax above
      the flux B must carry, the pool ``km*v/(vmax - v)`` delivers it exactly, so entry
      fluxes more than 300x apart reproduce the same measured A and C to nine decimals.
      ``tests/test_pathway_spec.py`` runs that through
      :func:`~ystwin.pathway.solve.solve_pathway` rather than asserting it. It is the ridge
      ``docs/research/KINETIC_FIT.md`` section 3 already reports for the desaturase step,
      where "six such parameter sets spanning three orders of magnitude in vmax" all
      "reproduce all six states to machine precision, differing only in a pool nobody
      measured".
    * **It says no where the scalar IS identified.** `phb` is passthrough end to end, so
      its only pool is the product and the entry flux is ``mu*[PHB]`` with no parameter in
      it -- feeding the eleven Kocharin 2013 chemostat states back through the solver with
      ``kinetics={}`` returns their tabulated contents to 4.7e-6 relative at worst, and
      that residual is the six-figure rounding of the table's own columns rather than
      anything the solver did. The count called it non-calibratable for a stated reason --
      "the entry flux and the branch capacity trade off along a ridge" -- that this chain
      structurally cannot have, because it has no branch capacity. Right answer, wrong
      mechanism, and the mechanism is what a reader acts on.

    ``entry flux = mu * sum(pools)`` applies only when every node's fate is ``DILUTED``.
    A terminal ``DEGRADED`` or ``SECRETED`` fate is accepted at load but requires an
    additional loss term. The spec supplies no measured loss flux or loss kinetics, so
    those fates make ``flux_is_recoverable`` false even if every pool is measured.
    ``required_loss_parameters`` names what a matched kinetic calibration must supply;
    an independently measured loss flux can also close the bookkeeping.

    And a second thing, which belongs to the DATASET rather than to the spec: that the
    state the pools were measured in is a steady one, with mu known. That identity is the
    sum of the ``d[X]/dt = 0`` balances, so in a batch culture -- mu moving through the
    growth phases, a storage pool still filling -- ``mu * sum(pools)`` is not the entry
    flux at all. This is not hypothetical for the spec being flagged here: `phb`'s own
    GENOTYPE dataset, Kocharin 2012, is shake flasks and batch bioreactors, and
    `data/pathways/phb.toml` says so in its own notes. The eleven-state chemostat series
    beside it, Kocharin 2013, is the one this identity applies to -- and it does apply,
    exactly: feeding those states back through the solver returns their tabulated contents
    to 4.7e-6 relative at worst, which is the six-figure rounding of the table's own
    columns. A verdict of ``flux_is_recoverable`` says the chain hides no pool or terminal
    loss. It cannot say the culture was at steady state, and nothing in this module can.

    The structural test is swept rather than argued. ``tests/test_pathway_spec.py`` builds
    every chain of two to four nodes -- every rate law at every position, every node
    measured or not, 1,548 specs -- solves each one and checks the verdict against what
    the solver actually put in each node. The nodes it excludes from ``measured_pools``
    are exactly the nodes the solver leaves empty, in all 1,548. Where it says recoverable,
    mu times the measured pools returns the entry flux to 5.8e-13 relative; where it says
    not, the hidden pools take at least 5.9% of that flux with it, so nothing in the sweep
    is a borderline call that a tolerance could flip.

    So the test is structural. A node holds a pool iff it is the terminal node, regardless
    of its rate law, or its rate law is not ``passthrough``. A terminal loss adds another
    measurement requirement beyond those pools. The recovery and fitted-pool conditions
    fail differently:

    Args:
        flux_is_recoverable: Every pool is measured and the terminal fate is ``DILUTED``,
            so the entry flux is ``mu`` times measurements with no fitted parameter.
            False requires a missing pool or an independent terminal-loss measurement.
        has_a_fitted_pool: Some pool above the product is measured, so branch kinetics
            can be constrained once entry flux is known. False alone is not an
            identifiability failure: on a DILUTED all-passthrough chain an expression
            series can identify ``alpha``, although ``solve_pathway`` then simply returns
            ``content = v_in/mu``. PHB has that structure and needs an expression axis.
            This flag does not count terminal loss kinetics, which require separate
            measurements before gross entry flux can be recovered.
        measured_pools: Pool-bearing nodes with a measurement.
        unmeasured_pools: Pool-bearing nodes without one. The free parameters.
        unfalsifiable_zeros: Nodes declared ``passthrough`` AND unmeasurable. Recovery rests
            on each of these holding no pool, and nothing can ever check that. Carried
            because it is the caveat on a calibration that otherwise looks clean:
            `beta_carotene` has one, phytoene, and `beta_carotene.toml` records that
            Chen 2016 measures it at 3.99% of total carotenoid -- so the declaration
            carrying the identification is one the file itself calls known-wrong.
        reason: Which recovery or fitted-pool requirement failed and what would lift it.
        required_loss_parameters: Terminal loss kinetics absent from the spec, selected by
            the terminal node's fate. Empty only when dilution is its sole outlet.
    """

    flux_is_recoverable: bool
    has_a_fitted_pool: bool
    measured_pools: tuple[str, ...]
    unmeasured_pools: tuple[str, ...]
    unfalsifiable_zeros: tuple[str, ...]
    reason: str
    required_loss_parameters: tuple[str, ...] = ()

    @property
    def calibratable(self) -> bool:
        """Both conditions. The scalar is identified and something can be fitted to it."""
        return self.flux_is_recoverable and self.has_a_fitted_pool


@dataclass(frozen=True)
class PathwaySpec:
    """One heterologous pathway, as declared.

    Args:
        product: The terminal species name. Must be the last node.
        organism: Host the spec is written for.
        nodes: The chain, in order, first committed intermediate first.
        entry_enzyme: Gene whose expression sets the flux INTO the chain. This is the
            first committed heterologous step, and naming it is a claim about where flux
            control sits -- see ``docs/research/KINETIC_FIT.md`` section 7.
        precursor_metabolite: GEM metabolite id the pathway is pulled off, for the audit
            layer's precursor-budget check.
        precursor_stoichiometry: Moles of precursor per mole of the first node.
        full_pathway_stoichiometry: The NET reaction of the whole heterologous pathway, per
            mole of PRODUCT, as ``{GEM metabolite id: coefficient}`` with consumed species
            negative. Optional, and absent means "nobody has written this pathway's cofactor
            balance down" rather than "it has none" -- `pathway/gem_environment.py` raises
            :class:`~ystwin.pathway.gem_environment.FullPathwayUndeclared` rather than
            falling back to the bare precursor pull, for the same reason its
            ``PrecursorCarrierUnknown`` refuses a missing carrier: a cofactor cost omitted is
            not a cofactor cost of zero.
        source: Citation for the chemistry.
        notes: Anything a reader needs that the fields cannot carry.
    """

    product: str
    organism: str
    nodes: tuple[Node, ...]
    entry_enzyme: str
    precursor_metabolite: str
    precursor_stoichiometry: float = 1.0
    full_pathway_stoichiometry: dict[str, float] | None = None
    source: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.nodes:
            raise ValueError(f"pathway {self.product!r} declares no nodes")
        names = [n.name for n in self.nodes]
        if len(set(names)) != len(names):
            raise ValueError(f"pathway {self.product!r} repeats a node name: {names}")
        if names[-1] != self.product:
            raise ValueError(
                f"pathway {self.product!r}: the last node is {names[-1]!r}. The product "
                "must be the terminal node, because the solver walks the chain in order "
                "and what comes out of the end is what is reported.")
        # INTERMEDIATES only -- `self.nodes[:-1]`. A terminal node that exports is one
        # measured capacity short of a closable balance, so it is refused by `solve_pathway`
        # by name; an INTERMEDIATE that also leaves the chain is a different thing, and it
        # is refused here because no measurement fixes it. The solver walks a chain by
        # carrying one flux forward, and an intermediate losing carbon sideways breaks the
        # walk itself rather than a coefficient in it.
        #
        # DEGRADED is checked here too, and that is a repair rather than a tidy-up. When the
        # DEGRADED refusal moved from load time to solve time it landed inside the terminal
        # arm, which silently NARROWED the guard from every node to the last one: a chain
        # with a degraded INTERMEDIATE and no rate constant loaded, solved, and returned the
        # diluted answer with its largest outlet dropped and no refusal at either place.
        leaking = [n.name for n in self.nodes[:-1]
                   if n.fate in (Fate.SECRETED, Fate.DEGRADED)]
        if leaking:
            raise ValueError(
                f"pathway {self.product!r} declares the intermediate(s) {leaking} as "
                "secreted or degraded, and this walk carries one flux forward from each "
                "node to the next. An intermediate that also leaves the chain has no "
                "mu*[X] dilution term doing that work -- it goes at a rate the cell sets, "
                "not one growth sets -- so the flux reaching the next node is not the flux "
                "this one passed on, and the balance would return a number that means "
                "nothing. Unlike a terminal node, no measurement lifts this: it is a "
                "property of the chain's shape. See docs/CLAIM_BOUNDARY.md.")
        # DEGRADED is NOT refused here, and it was until 2026-09-01. The refusal read: the
        # solver gives a terminal node exactly two outlets, so it would return
        # `content = flux/mu` and omit the term doing most of the work. That was true of
        # the solver and not of the pathway -- glycogen is not outside what this arithmetic
        # can describe, it is one rate constant short of it. `solve_pathway` now carries a
        # third outlet and refuses at SOLVE time when the constant is absent, which draws
        # the line where the missing thing actually is: a spec declares chemistry, and a
        # rate constant is a fitted number that belongs in `kinetics` beside a citation.
        #
        # The consequence for a caller is smaller than it looks. `glycogen.toml` still
        # cannot produce a number, because nobody has measured GPH1/SGA1 turnover in a
        # producing strain -- it now says so with the name of the missing measurement
        # instead of with "outside what this solver describes".
        terminal = self.nodes[-1]
        if terminal.molar_mass_g_per_mol is None:
            raise ValueError(
                f"pathway {self.product!r}: the terminal node needs "
                "molar_mass_g_per_mol, because content is reported in mg/gDCW and the "
                "alternative is a magic number in the reporting layer")
        if terminal.enzyme:
            raise ValueError(
                f"pathway {self.product!r}: the terminal node {terminal.name!r} names an "
                f"enzyme {terminal.enzyme!r}, but nothing consumes the product. Drop the "
                "enzyme, or add the node it feeds.")
        if self.precursor_stoichiometry <= 0:
            raise ValueError("precursor_stoichiometry must be positive")
        if self.full_pathway_stoichiometry is not None:
            if not self.full_pathway_stoichiometry:
                raise ValueError(
                    f"pathway {self.product!r} declares an EMPTY full_pathway_stoichiometry. "
                    "Omit the key to say the balance is unwritten; an empty table says the "
                    "pathway consumes nothing, which is a different and false claim")
            if self.precursor_metabolite not in self.full_pathway_stoichiometry:
                raise ValueError(
                    f"pathway {self.product!r}: full_pathway_stoichiometry does not mention "
                    f"the precursor {self.precursor_metabolite!r}. The net reaction is "
                    "written per mole of PRODUCT and must consume the precursor the rest of "
                    "this spec is built on")
            supplied = self.full_pathway_stoichiometry[self.precursor_metabolite]
            if supplied != -self.precursor_stoichiometry:
                raise ValueError(
                    f"pathway {self.product!r}: full_pathway_stoichiometry consumes "
                    f"{-supplied:g} of {self.precursor_metabolite!r} but "
                    f"precursor_stoichiometry says {self.precursor_stoichiometry:g}. Two "
                    "declarations of the same quantity that disagree is worse than one")
        claimed: dict[str, str] = {}
        for node in self.nodes:
            if not node.metabolite:
                continue
            if node.metabolite in claimed:
                raise ValueError(
                    f"pathway {self.product!r}: nodes {claimed[node.metabolite]!r} and "
                    f"{node.name!r} are both declared to be {node.metabolite!r}. One "
                    "content would silently replace the other in the thermodynamic gate, "
                    "and which one depends on the order this file lists them in. A node "
                    "IS a species")
            claimed[node.metabolite] = node.name
        if self.precursor_metabolite in claimed:
            raise ValueError(
                f"pathway {self.product!r}: node {claimed[self.precursor_metabolite]!r} is "
                f"declared to be {self.precursor_metabolite!r}, which is also the "
                "precursor. The precursor is what the pathway is pulled OFF and is not a "
                "node in the chain, so the same id in both places means the entry step "
                "consumes and produces the same species")

    @property
    def intermediates(self) -> tuple[Node, ...]:
        """Every node before the product."""
        return self.nodes[:-1]

    @property
    def pool_bearing_nodes(self) -> tuple[Node, ...]:
        """Nodes whose steady state holds a pool, in chain order.

        The terminal node always does, whatever rate law it declares; its fate determines
        whether dilution alone or dilution plus a measured loss balances incoming flux.
        An intermediate holds a pool unless it declares ``passthrough``, whose steady
        state is exactly zero. A ``proportional`` node holds ``v_in/(k + mu)`` and a
        ``saturating`` one the positive root of its quadratic; both contain a constant
        that is not known until somebody fits it.

        Entry flux sums dilution over this set PLUS any terminal loss. Calibratability
        therefore checks both these pool measurements and the declared terminal fate.
        """
        upstream = tuple(n for n in self.intermediates if n.rate_law != RateLaw.PASSTHROUGH)
        return upstream + (self.nodes[-1],)

    @property
    def calibratability(self) -> Calibratability:
        """Whether this pathway's flux scalar can be fitted from its own measurements.

        Read :class:`Calibratability` for why this is not a count of measurable nodes.
        ``load_pathway`` logs ``reason`` when the answer is no, so a spec that cannot be
        fitted says so on the way past rather than waiting to be asked.
        """
        pools = self.pool_bearing_nodes
        measured = tuple(n.name for n in pools if n.measurable)
        unmeasured = tuple(n.name for n in pools if not n.measurable)
        zeros = tuple(n.name for n in self.intermediates
                      if n.rate_law == RateLaw.PASSTHROUGH and not n.measurable)
        terminal_node = self.nodes[-1]
        required_loss_parameters: tuple[str, ...] = ()
        entry_balance = "mu times the sum of every pool"
        if terminal_node.fate == Fate.DEGRADED:
            required_loss_parameters = ("degradation_rate_per_h",)
            entry_balance = (
                f"mu*sum(pools) + degradation_rate_per_h*[{terminal_node.name}]")
        elif terminal_node.fate == Fate.SECRETED:
            required_loss_parameters = (
                "secretion_vmax_mmol_per_gdcw_h", "secretion_km_mmol_per_gdcw",
                "growth_rate_range")
            entry_balance = (
                f"mu*sum(pools) + secretion_vmax_mmol_per_gdcw_h*[{terminal_node.name}]"
                f"/(secretion_km_mmol_per_gdcw + [{terminal_node.name}])")
        recoverable = not unmeasured and not required_loss_parameters
        fitted = any(n.measurable for n in pools[:-1])

        if unmeasured:
            terminal = terminal_node.name
            mechanisms = ". ".join(_hidden_pool_mechanism(n, n.name == terminal)
                                   for n in pools if not n.measurable)
            bounded = (f"what the measurements bound is mu times {list(measured)} and "
                       "nothing more" if measured else
                       "nothing here is measured to bound it at all")
            lifts = (
                f"measure {terminal!r}, and no declaration substitutes for it"
                if terminal in unmeasured else
                f"measure {unmeasured[0]!r}, or declare it 'passthrough' if it truly holds "
                "no pool, which is a claim about magnitude and is itself checkable")
            reason = (
                f"pathway {self.product!r} CANNOT BE CALIBRATED: {list(unmeasured)} hold a "
                f"pool that nothing measures. The entry flux is {entry_balance}, so an "
                f"unmeasured pool is a free parameter in it: {bounded}. "
                f"{mechanisms}. What lifts it: {lifts}")
        elif required_loss_parameters:
            reason = (
                f"pathway {self.product!r} CANNOT BE CALIBRATED from pool measurements "
                f"alone: {list(measured)} are measured, but their dilution flux does not "
                "include the terminal loss")
        elif not fitted:
            above = (f"every node above {self.product!r} is declared passthrough"
                     if self.intermediates else
                     f"this chain declares no node above {self.product!r} at all")
            reason = (
                f"pathway {self.product!r} CANNOT BE CALIBRATED, and NOT because its entry "
                f"flux is unidentifiable -- {above}, so the entry flux is exactly "
                f"mu*[{self.product}] and there is no branch capacity for it to trade off "
                "against. The failure is the opposite one: with no pool but the product, "
                "solve_pathway returns content = v_in/mu, which is its own input, so no "
                "content this pathway can produce disagrees with any measurement at any "
                "parameter value. A scalar nothing here can contradict is restated rather "
                "than fitted. What lifts it is a MEASURED POOL ABOVE THE PRODUCT -- a node "
                "that actually holds one, so a second passthrough node buys nothing -- and "
                "what that buys is a parameter a measurement can contradict. What it does "
                "NOT buy is a fitted flux scalar: that needs a set of states whose relative "
                f"{self.entry_enzyme} expression actually differs, which no spec can carry "
                "because it is a property of a dataset and not of the chemistry. Until "
                "then use this spec to exercise the solver and the environment layer, "
                "which is what it is good for")
        else:
            reason = (
                f"pathway {self.product!r} can be calibrated: {list(measured)} each hold a "
                "pool and each is measured, so the entry flux is mu times their sum with "
                "no fitted parameter in it")
            if zeros:
                reason += (
                    f". That rests on {list(zeros)} holding no pool at all, which is "
                    "declared and never measured -- if any of them accumulates, the entry "
                    "flux is understated by mu times that pool and nothing here would say "
                    "so")
        if required_loss_parameters:
            reason += (
                f". Terminal node {terminal_node.name!r} is declared "
                f"{terminal_node.fate.upper()}; entry flux = {entry_balance}. What lifts "
                "the loss ambiguity is an independently measured loss flux, or a matched "
                f"calibration of {', '.join(required_loss_parameters)}. Those measurements "
                "are not supplied by this chemistry spec; pool measurements and mu alone "
                "do not identify gross synthesis. Solving with supplied loss kinetics is "
                "possible, but does not make those kinetics identifiable from the pools")
        return Calibratability(recoverable, fitted, measured, unmeasured, zeros, reason,
                              required_loss_parameters=required_loss_parameters)

    @property
    def calibratable(self) -> bool:
        """Shorthand for ``self.calibratability.calibratable``.

        Kept as a bare boolean because `predict.py` reads it to decide whether to attach a
        note. Anything deciding what to DO should read :attr:`calibratability` instead:
        a missing pool, an unmeasured terminal loss and a missing expression axis require
        different measurements, and a boolean cannot tell them apart.
        """
        return self.calibratability.calibratable

    def require_calibratable(self) -> None:
        """Refuse to proceed unless this pathway's own measurements identify its scalar.

        The explicit gate for a fitting path, kept separate from loading because loading
        must not refuse. Nothing calls it yet, and that is worth saying rather than
        hiding: `scripts/fit_pathway_flux.py` never loads a spec at all -- it reads
        `data/carotenoid/*.tsv` straight into :func:`~ystwin.pathway.flux.fit_flux_law` --
        so today the load-time record is the ONLY thing standing between a
        non-calibratable spec and a fit. This is what that script should call when it
        grows a ``--pathway`` argument.

        Raises:
            ValueError: carrying :attr:`Calibratability.reason`, which names which of the
                two conditions failed and what would lift it.
        """
        verdict = self.calibratability
        if not verdict.calibratable:
            raise ValueError(verdict.reason)

    @property
    def thermo_metabolites(self) -> dict[str, str]:
        """Node name -> GEM metabolite id, for every node that declares one.

        What :func:`~ystwin.pathway.thermo_gate.gate_pathway` wants for ``metabolites``.
        Partial by design: a node with no declared species is simply not gated, and the
        gate reports it under ``ungated`` rather than treating silence as clearance.
        """
        return {n.name: n.metabolite for n in self.nodes if n.metabolite}

    @property
    def thermo_steps(self) -> dict[str, dict[str, float] | str]:
        """Node name -> the step that PRODUCES it, ready to hand to the gate.

        A declared ``stoichiometry`` wins over a bare ``reaction`` id, and the preference
        is the whole reason both fields exist: a stoichiometry needs no model, so the gate
        runs on a machine with no GEM and without paying an 18-second SBML parse.

        That preference is also why it does **not** cross-check the pair. `gate_step`
        refuses a ``(stoichiometry, reaction)`` pair that describes two different
        reactions, but it only ever sees one of them here -- ``_resolve`` on a mapping
        returns no reaction, and ``_resolve`` on a reaction id returns that reaction's own
        stoichiometry, so the pair agrees trivially either way. The comparison that
        actually bites is :meth:`check_steps_against`, and a caller holding a model should
        run it.
        """
        steps: dict[str, dict[str, float] | str] = {}
        for node in self.nodes:
            if node.stoichiometry:
                steps[node.name] = dict(node.stoichiometry)
            elif node.reaction:
                steps[node.name] = node.reaction
        return steps

    def check_steps_against(self, model) -> tuple[str, ...]:
        """Refuse a declared stoichiometry that has drifted from the model's reaction.

        Both shipped specs declare a stoichiometry *and* a reaction id for the same step.
        The duplication is deliberate -- the mapping is what lets the gate run with no SBML
        parse -- but a duplicate is only worth having if something compares the copies, and
        nothing did: whichever one the gate resolves, it agrees with itself. `beta_carotene`
        is the case that matters, because its reactions are built in Python by
        `fba/carotenoid.py` rather than read from a file, so a coefficient could change
        there and the TOML would keep gating the old chemistry under the new name.

        Returns the node names actually compared, so a caller can tell "all agree" from
        "nothing was checked" -- a heterologous step the model does not carry yet is
        skipped rather than treated as a mismatch.

        Raises:
            ValueError: naming the node, the reaction and the coefficients that differ.
        """
        checked: list[str] = []
        for node in self.nodes:
            if not node.stoichiometry or not node.reaction:
                continue
            if not model.reactions.has_id(node.reaction):
                continue
            theirs = {m.id: float(c)
                      for m, c in model.reactions.get_by_id(node.reaction).metabolites.items()}
            mine = {str(m): float(c) for m, c in node.stoichiometry.items()}
            if theirs != mine:
                only_spec = sorted(set(mine) - set(theirs))
                only_model = sorted(set(theirs) - set(mine))
                differing = sorted(m for m in set(mine) & set(theirs) if mine[m] != theirs[m])
                raise ValueError(
                    f"pathway {self.product!r}, node {node.name!r}: the stoichiometry "
                    f"declared in the spec and reaction {node.reaction!r} in the model are "
                    "different reactions. The gate would take dG0 from one and Q from the "
                    f"other. Only in the spec: {only_spec}; only in the model: "
                    f"{only_model}; different coefficients: {differing}")
            checked.append(node.name)
        return tuple(checked)

    @property
    def gateable(self) -> bool:
        """Whether any step of this chain declares the chemistry the gate needs."""
        return bool(self.thermo_steps)

    def node(self, name: str) -> Node:
        for candidate in self.nodes:
            if candidate.name == name:
                return candidate
        raise KeyError(f"pathway {self.product!r} has no node {name!r}; "
                       f"have {[n.name for n in self.nodes]}")


def pathways_dir() -> pathlib.Path:
    """Where the specs live. Tracked, because a pathway is part of the model."""
    from .. import paths

    return paths.data_dir() / "pathways"


def available_pathways(directory: pathlib.Path | None = None) -> list[str]:
    """Product names with a spec on disk, sorted."""
    directory = pathways_dir() if directory is None else directory
    if not directory.is_dir():
        return []
    return sorted(p.stem for p in directory.glob("*.toml"))


def load_pathway(name: str, directory: pathlib.Path | None = None) -> PathwaySpec:
    """Read one spec, and say on the way past if it cannot be calibrated.

    **The record is not a refusal, deliberately.** A non-calibratable spec is still a
    working capability -- `phb` is passthrough end to end and exercises the solver, the
    growth-rate guard and the environment layer on a second product and a second paper --
    so refusing it would remove something that works to prevent a mistake nobody has made.
    What it must not do is load silently, because ``calibratable`` was discoverable only by
    asking, and nobody asks a property they do not know exists.

    It goes to :mod:`logging` rather than :mod:`warnings` for a reason that is about this
    repository and not about taste: `pyproject.toml` sets ``filterwarnings = ["error"]``,
    so ``warnings.warn`` here would raise inside every test that loads `phb` -- which is
    the refusal the paragraph above rules out. An unconfigured process routes the record to
    stderr through ``logging.lastResort``, which is where somebody running a script is
    looking.

    Args:
        name: File stem under ``data/pathways``, e.g. ``"beta_carotene"``.
        directory: Override the location.

    Raises:
        FileNotFoundError: naming what is available, because a typo in a product name
            should not read like an unsupported product.
        ValueError: from :class:`PathwaySpec` or :class:`Node` on a malformed spec.
    """
    directory = pathways_dir() if directory is None else directory
    path = directory / f"{name}.toml"
    if not path.exists():
        raise FileNotFoundError(
            f"no pathway spec {path}. Available: {available_pathways(directory) or 'none'}")
    with path.open("rb") as handle:
        raw = tomllib.load(handle)

    if "pathway" not in raw:
        raise ValueError(f"{path} has no [pathway] table")
    header = dict(raw["pathway"])
    declared = raw.get("node", [])
    if not declared:
        raise ValueError(f"{path} declares no [[node]] entries")
    nodes = tuple(Node(**entry) for entry in declared)
    spec = PathwaySpec(nodes=nodes, **header)
    verdict = spec.calibratability
    if not verdict.calibratable:
        _LOG.warning("%s: %s", path.name, verdict.reason)
    return spec
