"""The missing layer: latent module activity -> the reactions each regulon controls -> bounds.

`bridge/latent_bridge.py` collapses all 24 modules into one ATP maintenance scalar. Every
module's specificity is discarded before metabolism sees anything, so nothing tells the
model *which* reactions stress turns down, and D1 (`scripts/parked/run_d1.py`) then reports
the product flux range as `[0, ceiling]` with relative width 1.000 at every growth level.
This module is the alternative: per-reaction regulation, mapped through an authoritative
gene-reaction relation rather than through an invented scalar.

Two mechanisms, and the difference between them is the whole finding.

**E-Flux** (:func:`eflux_layer`, Colijn 2009, PMID 19714220) scales the *upper* bound of
every reaction a regulon controls. It is the simplest defensible form and it is
implemented first because it needs no measurement beyond the latent state itself. It also
**cannot lift a lower bound off zero**, ever, by construction -- see
:data:`EFLUX_LIFTS_LOWER_BOUND` and the note under :class:`EFluxLayer`. Since D1's
relative width is `(max - min) / max` and its `min` is zero because nothing in the network
requires product formation, E-Flux cannot move D1's headline number however specific it
gets. It moves the ceiling. Those are different claims and the report must not merge them.

**And on this panel it does not move the ceiling either.** Driven by every one of the 25
stressors at 3x its own EC50, the number of bounds this layer tightens is **zero** on both
shipped GEMs, because every module that has a regulon is non-negative and
`scale = max(0, 1 + a) >= 1` always. The layer is therefore `LayerState.INERT`, permanently
and by construction rather than by dose -- :data:`EFLUX_LAYER_STATE` carries the label, the
measurement and the two-part reason, and :data:`EFLUX_BOUNDS_TIGHTENED_AT_3X_EC50` carries
the count. It stays here as the null arm beside :func:`task_priority_order`; it is not a live
regulation route.

**Measured task-efficiency adaptation** (:func:`task_priority_order`, inspired by
COSMIC-dFBA, PMID 38387677) does lift a floor. Its lower-bound task-mining precedent is
verified against the supplied paper §4.3 and supplementary algorithm. The original CHO
model and full fitted execution configuration are still missing; this is not their numerical
reproduction. See ``data/cosmic_sources.json``. For each candidate
metabolic task -- biomass, product, each secreted byproduct -- it takes the ratio of the
*measured* flux to that task's individually maximised FBA flux, calls that the task
efficiency, and enforces the measured flux of the most efficient task as a lower bound
before re-ranking the rest. The floor therefore comes from measurement, not from the
model. That is also its cost, and it is the reason this implementation **refuses rather
than defaults** when a task carries no measured flux: substituting a simulated or assumed
flux would choose the answer, and the narrow interval that came back would be an artefact
of the substitution. No strain here carries the beta-carotene pathway, so no measured
secretion rate exists and :func:`task_priority_order` will raise on this repository's own
product today. Naming the missing measurement is the deliverable.

Precedent, so the novelty is not overstated
-------------------------------------------
Coupling a regulatory state to a GEM's bounds is not new. **PROM** (Chandrasekaran & Price
2010, PMID 20876091) does exactly this pattern -- probabilistic regulatory state to flux
bounds -- in *E. coli* and *M. tuberculosis*, and it is sixteen years old. What is new here
is the organism and the fact that the state is inferred from live reporters rather than from
a transcriptome; the mechanism is borrowed. And the field's own systematic comparison
(Machado & Herrgard 2014, PMID 24762745) found that expression-integration methods
frequently fail to beat plain pFBA, so any claim this layer improves a prediction owes a
pFBA arm and a no-regulation arm beside it. Both are what :file:`scripts/run_regulation.py`
prints as its unregulated column.

Claim tier
----------
Everything here is **Tier 1 at best** -- a deterministic computation on models and a
binding map we adopted (`docs/CLAIM_BOUNDARY.md`). The gene-reaction relation is the GEM's
own; the module-to-gene relation is SGD's curated regulation records restricted to one
binding-motif source; the *magnitudes* relating regulon activity to enzyme capacity are
**asserted** and marked as such at :data:`_ACTIVITY_TO_CAPACITY_NOTE`. Nothing in this
repository measures them.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
import pathlib
from dataclasses import asdict, dataclass
from numbers import Real

import cobra

from .. import paths
from ..fba.solver import growth_or_none

__all__ = [
    "COLIJN_PMID",
    "COSMIC_PMID",
    "EFLUX_BOUNDS_TIGHTENED_AT_3X_EC50",
    "EFLUX_LAYER_STATE",
    "EFLUX_LIFTS_LOWER_BOUND",
    "MACISAAC_PMID",
    "MODULE_FACTORS",
    "AppliedRegulation",
    "AppliedUptake",
    "constrain_uptake",
    "exchange_coefficients",
    "specific_rate",
    "EFluxLayer",
    "InfeasibleRegulation",
    "MetabolicTask",
    "MissingMeasurement",
    "ReactionScale",
    "Regulon",
    "TaskPriority",
    "apply_eflux",
    "apply_task_bounds",
    "eflux_layer",
    "gene_scales",
    "load_regulons",
    "override_capacity",
    "reference_capacity",
    "task_priority_order",
]

# MacIsaac et al. 2006 -- "An improved map of conserved regulatory sites for Saccharomyces
# cerevisiae", BMC Bioinformatics. Conserved binding motifs over the Harbison and Lee ChIP
# data. Restricting SGD's regulation records to this one reference is not a convenience:
# the unrestricted set runs to 500+ promiscuous ChIP hits per factor, and anything carrying
# expression evidence would make the regulon partly a clustering of expression data, which
# is what this layer is supposed to be independent of. Same restriction, same reason, as
# `scripts/run_external_validation.py`.
MACISAAC_PMID = 16522208

# Colijn et al. 2009, PMID 19714220 -- "Interpreting expression data with metabolic flux
# models: predicting Mycobacterium tuberculosis mycolic acid production", PLoS Comput Biol.
# E-Flux: expression bounds reaction capacity, `|v_i| <= c * e_i`.
COLIJN_PMID = 19714220

# Gopalakrishnan et al. 2024, PMID 38387677 -- "COSMIC-dFBA: A novel multi-scale hybrid
# framework for bioprocess modeling", Metab Eng 82:183-192. The lexicographic
# task-efficiency ordering of its §4.3 is what lifts a product lower bound off zero.
COSMIC_PMID = 38387677

EFLUX_LIFTS_LOWER_BOUND = False
"""E-Flux constrains capacity, so it can never require flux. Stated as a value, not prose.

A named constant rather than a sentence in a docstring because the question "can the
regulation layer rescue D1's relative width" has exactly this answer, and a caller
reporting on the layer should be able to assert on it rather than paraphrase it.

Precisely: E-Flux is ``|v| <= c * e``, a bound on magnitude. On a reversible reaction it
does move the *numerical* lower bound, pulling a negative one up toward zero -- so the
implementation below writes lower bounds and this constant still reads ``False``. The
distinction is the one that matters for D1: the reachable minimum of ``|v|`` stays zero,
because a constraint on how much flux is permitted never says any flux is required. Only
a lower bound strictly above zero narrows ``(max - min) / max``, and no capacity
constraint can produce one. :func:`task_priority_order` can, and needs a measurement to.
"""

EFLUX_LAYER_STATE = "inert"
"""E-Flux is ``LayerState.INERT`` here, permanently, and this is the measurement behind it.

Not a judgement and not a deprecation: an inert layer *ran and provably could not have
changed anything*, which is exactly `predict.LayerState.INERT`. The string is duplicated
rather than imported because `predict` is the orchestrator and a bridge module importing it
inverts the dependency; :data:`EFLUX_LAYER_STATE` is what a report should print.

**Measured** on 2026-09-05, by driving every one of the 25 stressors in
`generator/stress_panel.py::STRESSORS` to 3x its own EC50, taking `module_response` at that
dose and counting `EFluxLayer.binding`: **0 bounds tightened**, on yeast-GEM 9.0.2 and on
ecYeastGEM_batch alike. Zero out of scales that genuinely computed, not zero because nothing
ran -- each stressor bounds 54-218 reactions on yeast-GEM and 16-80 on ecYeastGEM, and every
module at once bounds 334 and 124. `tests/test_eflux_capacities.py` re-runs it.

The reason is structural, not a property of the doses chosen, and it has two halves -- both
measured, because the first half is stated wrong in the architecture document. **Every module
that carries a transcription factor is non-negative**: the minimum of `module_response` over
:data:`MODULE_FACTORS`, over all 25 stressors, over a 31-point dose sweep from 0 to 3x EC50,
is exactly 0.0. So ``scale = max(0, 1 + a) >= 1`` on every regulon that exists, and this layer
only ever writes an upper bound at ``scale * capacity`` -- raising it. **The four modules that
do go negative -- atp (-0.381), nadh (-0.150), ph (-0.448), redox (-0.109) -- are metabolite
pools with no transcription factor**, so they have no regulon by construction and this layer
cannot see them at all; `gene_scales` raises on them unless they are declared unmappable, and
declaring them excludes them. It is *not* true that no activity in the panel is negative.

Forcing a regulon's activity negative by hand does not rescue the layer either: gene scales
add across regulons, so a two-regulon 50% repression clips an essential reaction to zero and
the model dies. It is a knockout operator, not a dial.

Kept, not deleted. It is the honest null arm beside :func:`task_priority_order`, and
Machado & Herrgard 2014, PMID 24762745 is the field's own finding that expression-integration
methods frequently fail to beat plain pFBA -- so a layer measured to change nothing is a
result, not a bug. What it must not be is offered as a live regulation route.
"""

EFLUX_BOUNDS_TIGHTENED_AT_3X_EC50 = 0
"""The count behind :data:`EFLUX_LAYER_STATE`, as a number a test can assert on."""

_ACTIVITY_TO_CAPACITY_NOTE = (
    "ASSERTED: the map from module activity a to relative enzyme capacity, "
    "scale = max(0, 1 + a), is a choice and not a measurement. Two parts of it are "
    "defensible from the panel itself -- the unstressed state is a = 0 by construction in "
    "generator/stress_panel.py::module_response, so scale(0) = 1 exactly and an inactive "
    "state must leave every bound alone; and activities of several modules on one promoter "
    "add, which is the same additive convention combination_response already uses. What is "
    "asserted is the LINEARITY and the SLOPE of 1 per unit activity. Nothing in this "
    "repository relates a regulon's activity to its enzymes' capacity, so no fitted slope "
    "exists to use instead. The clip at zero is not asserted: a negative capacity is not a "
    "quantity."
)

# Transcription factors named by each MODULES entry, as gene symbols. Identical to the
# table in `scripts/run_external_validation.py`, from the same reading of
# `generator/stress_panel.py`'s MODULES. Duplicated rather than imported because that is a
# script and this is library code, and a library importing a script inverts the
# dependency; the script should import this table instead, which is a change outside this
# module's remit and is reported rather than made.
#
# ACE1 is CUP2, the name SGD uses. The five pool modules -- redox, peroxide, atp, ph, nadh
# -- are absent BY CONSTRUCTION and not by omission: they have no transcription factor and
# therefore no regulon, because a metabolite pool is not transcribed. Asking this layer to
# regulate one of them raises.
MODULE_FACTORS: dict[str, tuple[str, ...]] = {
    "ESR": ("MSN2", "MSN4"), "UPR": ("HAC1",), "oxidative": ("YAP1", "SKN7"),
    "heat": ("HSF1",), "osmotic": ("SKO1", "HOT1"), "proteasome": ("RPN4",),
    "iron": ("AFT1", "AFT2"), "cell_wall": ("RLM1",), "dna_damage": ("RFX1",),
    "calcium": ("CRZ1",), "carbon": ("ADR1", "CAT8"), "nitrogen": ("GLN3", "GAT1"),
    "hypoxia": ("ROX1", "UPC2", "HAP1"), "copper": ("CUP2", "MAC1"), "zinc": ("ZAP1",),
    "sulfur": ("MET4", "MET31", "MET32"), "xenobiotic": ("PDR1", "PDR3"),
    "retrograde": ("RTG1", "RTG3"), "alkaline_ph": ("RIM101",),
}

_PROT_POOL_PREFIX = "prot_pool"


class MissingMeasurement(ValueError):
    """A task-efficiency ordering was asked for without the measurement it needs.

    Separate from ``ValueError`` at the call site so a caller can tell "you have not
    measured this yet" apart from "you passed nonsense". The first is the normal state of
    this repository and the message names the quantity and its units.
    """


class InfeasibleRegulation(ValueError):
    """Applying a regulation layer left the model unable to grow.

    Raised rather than returned because the failure mode this guards is a silent one: a
    layer tight enough to make the model infeasible produces `None` from
    ``slim_optimize``, and a caller that does not check reads that as a flux of zero.
    """

    def __init__(self, message, *, status="infeasible"):
        self.status = status
        super().__init__(message)


def specific_rate(value: float, unit: str = "mmol/gDW/h") -> float:
    """Convert an explicitly named specific-rate unit, never a sign convention."""
    factors = {"mmol/gDW/h": 1.0, "mmol/gDCW/h": 1.0, "umol/gDW/h": 1e-3,
               "umol/gDCW/h": 1e-3, "mmol/gDW/min": 60.0, "1/h": 1.0}
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value < 0:
        raise ValueError("specific rate must be a finite nonnegative number")
    if unit not in factors:
        raise ValueError(f"unsupported specific-rate unit {unit!r}; declare a mass/time conversion")
    result = float(value) * factors[unit]
    if not math.isfinite(result):
        raise ValueError("converted specific rate must be finite nonnegative")
    return result


def exchange_coefficients(
    model: cobra.Model, exchange_id: str, *, direction: str = "uptake",
) -> dict[str, float]:
    """Net physical boundary rate from S, including all split routes of the same metabolite.

    Positive uptake supplies the boundary metabolite; positive secretion removes it.
    Reaction names and suffixes carry no directional information. A caller must identify
    a boundary metabolite, not an intracellular transport reaction.
    """
    if direction not in {"uptake", "secretion"}:
        raise ValueError("exchange direction must be uptake or secretion")
    reaction = model.reactions.get_by_id(exchange_id)
    if len(reaction.metabolites) != 1:
        raise ValueError(f"{exchange_id}: expected a one-metabolite boundary reaction")
    metabolite = next(iter(reaction.metabolites))
    sign = 1.0 if direction == "uptake" else -1.0
    coefficients = {}
    for boundary in sorted(metabolite.reactions, key=lambda row: row.id):
        if len(boundary.metabolites) != 1:
            continue
        coefficient = float(boundary.metabolites[metabolite])
        if not math.isfinite(coefficient) or coefficient == 0:
            raise ValueError(f"{boundary.id}: invalid boundary stoichiometry")
        coefficients[boundary.id] = sign * coefficient
    return coefficients


def _expression(model, coefficients):
    return sum(coefficient * model.reactions.get_by_id(rid).flux_expression
               for rid, coefficient in coefficients.items())


def _constraint_name(prefix, coefficients):
    encoded = json.dumps(sorted(coefficients.items()), separators=(",", ":")).encode()
    return prefix + hashlib.sha256(encoded).hexdigest()[:20]


@dataclass(frozen=True)
class AppliedUptake:
    exchange: str
    magnitude: float
    unit: str
    mode: str
    source: str
    coefficients: tuple[tuple[str, float], ...]
    bounds: tuple[tuple[str, float, float], ...]

    def report(self):
        result = asdict(self)
        result["bounds"] = [(rid, None if math.isinf(lower) else lower,
                              None if math.isinf(upper) else upper) for rid, lower, upper in self.bounds]
        result["bound_null_means"] = "unbounded endpoint, not a zero bound or an unsolved flux"
        return result


def constrain_uptake(
    model: cobra.Model, exchange_id: str, magnitude: float, *,
    unit: str = "mmol/gDW/h", mode: str = "cap", source: str,
) -> AppliedUptake:
    """Set an explicit medium supply in place, preserving each model's split convention.

    Cap means net uptake <= magnitude, not a forced uptake. Fixed means a measured
    net rate, and must be explicitly requested. Existing open supply directions replace
    their medium bound; closed reverse/export directions are never opened by a guess.
    One net constraint prevents parallel imports from each receiving the entire budget.
    """
    if not isinstance(source, str) or not source.strip():
        raise MissingMeasurement("uptake constraints need an explicit measurement or assumption source")
    if mode not in {"cap", "fixed"}:
        raise ValueError("uptake mode must be cap or fixed")
    if unit == "1/h":
        raise ValueError("uptake needs an amount/mass/time unit, not a growth-rate unit")
    amount = specific_rate(magnitude, unit)
    coefficients = exchange_coefficients(model, exchange_id)
    changed = []
    for rid, coefficient in coefficients.items():
        reaction = model.reactions.get_by_id(rid)
        lower, upper = reaction.bounds
        if coefficient > 0 and upper > 0:
            upper = amount / coefficient
        elif coefficient < 0 and lower < 0:
            lower = amount / coefficient
        if lower > upper:
            raise InfeasibleRegulation(f"{rid}: measured uptake conflicts with a required flux")
        if reaction.bounds != (lower, upper):
            reaction.bounds = (lower, upper)
        changed.append((rid, float(lower), float(upper)))
    name = _constraint_name("regulation_uptake_", coefficients)
    if name in model.constraints:
        model.remove_cons_vars([model.constraints[name]])
    model.add_cons_vars([model.problem.Constraint(
        _expression(model, coefficients), lb=amount if mode == "fixed" else None,
        ub=amount, name=name,
    )])
    return AppliedUptake(exchange_id, amount, "mmol/gDW/h", mode, source,
                         tuple(coefficients.items()), tuple(changed))


def _set_floor(model, coefficients, floor):
    if len(coefficients) == 1:
        rid, coefficient = next(iter(coefficients.items()))
        reaction = model.reactions.get_by_id(rid)
        lower, upper = reaction.bounds
        if coefficient > 0:
            lower = max(lower, floor / coefficient)
        else:
            upper = min(upper, floor / coefficient)
        if lower > upper:
            raise InfeasibleRegulation(f"task floor {floor:g} contradicts {rid}'s existing capacity")
        changed = reaction.bounds != (lower, upper)
        reaction.bounds = lower, upper
        return int(changed)
    name = _constraint_name("regulation_floor_", coefficients)
    if name in model.constraints:
        previous = model.constraints[name]
        if previous.lb is not None and previous.lb >= floor:
            return 0
        model.remove_cons_vars([previous])
    model.add_cons_vars([model.problem.Constraint(_expression(model, coefficients), lb=floor, name=name)])
    return 1


def _solve_status(model, value):
    status = str(model.solver.status)
    if status == "optimal" and (value is None or not math.isfinite(value)):
        return "numerical_failure"
    return status


# --------------------------------------------------------------------------- #
# module -> gene, from SGD's curated regulation records
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Regulon:
    """The genes one module's transcription factors bind, and where that came from.

    Args:
        module: Module name, as in ``generator/stress_panel.py``'s ``MODULES``.
        factors: Gene symbols of the transcription factors, as SGD names them.
        genes: Systematic ORF names of the target genes.
        source: Provenance, carried so a regulon can be reported beside its origin.
    """

    module: str
    factors: tuple[str, ...]
    genes: tuple[str, ...]
    source: str

    def __len__(self) -> int:
        return len(self.genes)


def load_regulons(
    directory: pathlib.Path | None = None,
    pmid: int = MACISAAC_PMID,
    factors: dict[str, tuple[str, ...]] | None = None,
) -> dict[str, Regulon]:
    """Read module regulons out of SGD's curated regulation records.

    Only records in which the factor is the **regulator** are kept. SGD's files are
    per-gene rather than per-role, so ``FACTOR.json`` holds every record the factor
    appears in, including ones where it is somebody else's target: ``CAT8.json`` carries
    two MacIsaac records and CAT8 is the *target* in both. Taking ``locus2`` without
    checking ``locus1`` therefore puts the factor's own gene into its own regulon and, for
    CAT8, produces a one-gene regulon consisting entirely of CAT8. Reading the direction
    off the record is the fix.

    Args:
        directory: Directory of ``<FACTOR>.json`` files. Defaults to the vendored
            ``data/external/sgd/regulation``.
        pmid: Reference to restrict to. Defaults to the MacIsaac conserved-motif map.
        factors: Module-to-factor table. Defaults to :data:`MODULE_FACTORS`.

    Returns:
        One :class:`Regulon` per module in ``factors``. A module whose factors have no
        records at this reference gets an empty regulon rather than being dropped, so the
        gap is visible to the caller instead of absent from the result.

    Raises:
        FileNotFoundError: naming the factor whose records are missing. The regulon of a
            factor whose file is absent is not empty, it is unknown, and returning it as
            empty would silently under-regulate that module.
    """
    table = MODULE_FACTORS if factors is None else factors
    root = pathlib.Path(directory) if directory is not None else (
        paths.data_dir() / "external" / "sgd" / "regulation"
    )
    out: dict[str, Regulon] = {}
    for module, symbols in table.items():
        targets: set[str] = set()
        for symbol in symbols:
            path = root / f"{symbol}.json"
            if not path.exists():
                raise FileNotFoundError(
                    f"no SGD regulation records for {symbol!r} (module {module!r}) at "
                    f"{path}. An absent file is an unknown regulon, not an empty one; "
                    "fetch it rather than letting this module read as unregulated."
                )
            for record in json.loads(path.read_text()):
                if record.get("regulation_of") != "transcription":
                    continue
                if record.get("reference", {}).get("pubmed_id") != pmid:
                    continue
                if record.get("locus1", {}).get("display_name") != symbol:
                    continue  # the factor is the target here, not the regulator
                targets.add(record["locus2"]["format_name"])
        out[module] = Regulon(
            module=module, factors=tuple(symbols), genes=tuple(sorted(targets)),
            source=f"SGD curated regulation, restricted to PMID {pmid}",
        )
    return out


# --------------------------------------------------------------------------- #
# gene -> reaction, from the GEM's own gene-reaction rules
# --------------------------------------------------------------------------- #


def gene_scales(
    activity: dict[str, float],
    regulons: dict[str, Regulon],
    allow_unmapped: frozenset[str] | set[str] = frozenset(),
) -> dict[str, float]:
    """Relative capacity scale per gene, from signed module activity.

    Activity is signed and the sign is honoured, not absolute-valued: a module a stressor
    represses carries a negative weight in ``STRESSORS[...].targets`` and arrives here
    negative, and it must lower the bound. ``module_response`` has already multiplied the
    weight through, so the sign is in the activity.

    Activities of several modules on one gene add, following
    ``combination_response``'s convention, and the sum becomes a relative scale
    ``max(0, 1 + a)``. See :data:`_ACTIVITY_TO_CAPACITY_NOTE` for what in that is asserted.

    Args:
        activity: Module name to signed activity. Zero-activity entries are ignored, so
            passing a full ``module_response`` dict is fine.
        regulons: Output of :func:`load_regulons`.
        allow_unmapped: Modules the caller has deliberately accepted as unregulatable --
            the metabolite pools, and any module whose regulon is genuinely empty at this
            evidence level. Naming them is what stops a module silently doing nothing.

    Raises:
        ValueError: naming any module with non-zero activity that carries no regulon and
            was not listed in ``allow_unmapped``. A module that cannot reach a reaction
            and is not declared as such is the exact failure this layer exists to replace.
    """
    refused: list[str] = []
    scales: dict[str, float] = {}
    for module, value in activity.items():
        if value == 0.0:
            continue
        if module in allow_unmapped:
            continue
        regulon = regulons.get(module)
        if regulon is None:
            refused.append(
                f"{module!r}: not a transcriptional module -- it has no factor, so no "
                "regulon exists to scale (the five metabolite pools are like this)"
            )
            continue
        if not regulon.genes:
            refused.append(
                f"{module!r}: factors {list(regulon.factors)} have no target genes at "
                f"this evidence level ({regulon.source})"
            )
            continue
        for gene in regulon.genes:
            scales[gene] = scales.get(gene, 0.0) + float(value)
    if refused:
        raise ValueError(
            "these modules carry activity but cannot reach any reaction:\n  "
            + "\n  ".join(refused)
            + "\nList them in allow_unmapped to proceed with them excluded, and say so "
            "in whatever reports the result. Silently dropping them would make the layer "
            "look more specific than it is."
        )
    return {gene: max(0.0, 1.0 + total) for gene, total in scales.items()}


def _rule_scale(node: ast.expr | None, scales: dict[str, float]) -> float:
    """Fold gene scales up one gene-reaction rule. Unknown genes score 1.0.

    ``and`` takes the minimum and ``or`` takes the **maximum**. Colijn 2009 sums over
    ``or`` because it works in absolute expression levels, where two isozymes' capacities
    add. These are *relative* scales against the unstressed state, so summing would give
    two untouched isozymes a scale of 2.0 -- it would invent capacity out of the rule's
    shape. The maximum is the least-committal fold: a complex is limited by its scarcest
    subunit, and a reaction keeps its reference capacity as long as one isozyme does.

    A gene absent from ``scales`` -- either not in any regulon, or in a regulon this
    stressor does not touch -- scores 1.0 rather than 0.0. The panel says nothing about
    those genes, and treating silence as repression would let a 19-module panel constrain
    the whole proteome.
    """
    if node is None:
        return 1.0
    if isinstance(node, ast.Name):
        return scales.get(node.id, 1.0)
    if isinstance(node, ast.BoolOp):
        values = [_rule_scale(v, scales) for v in node.values]
        return min(values) if isinstance(node.op, ast.And) else max(values)
    if isinstance(node, ast.Expression):  # pragma: no cover - cobra hands us a Module
        return _rule_scale(node.body, scales)
    raise ValueError(f"unsupported node in a gene-reaction rule: {ast.dump(node)}")


def _protein_pool_cap(model: cobra.Model) -> tuple[float, str] | None:
    """The GECKO protein pool's total budget and the metabolite carrying it, if any.

    Returned rather than assumed so a plain Yeast9 model, which has no pool, takes the
    other branch in :func:`reference_capacity` instead of being special-cased by filename.
    """
    for metabolite in model.metabolites:
        if not metabolite.id.startswith(_PROT_POOL_PREFIX):
            continue
        supplies = []
        for reaction in metabolite.reactions:
            if len(reaction.metabolites) != 1:
                continue
            coefficient = reaction.metabolites[metabolite]
            capacity = max(coefficient * reaction.lower_bound, coefficient * reaction.upper_bound, 0.0)
            supplies.append(capacity)
        if supplies and all(math.isfinite(capacity) for capacity in supplies):
            return math.fsum(supplies), metabolite.id
    return None


def reference_capacity(
    reaction: cobra.Reaction,
    pool: tuple[float, str] | None,
) -> tuple[float, float] | None:
    """The reaction's own capacity ceilings, forward and reverse, unregulated.

    E-Flux needs an absolute scale, and the choice of it is the load-bearing decision in
    the whole layer -- pick it badly and the numbers are about the constant rather than
    about the regulation. So it is read off the model in both cases rather than asserted:

    * **An enzyme draw in a GECKO model.** ``draw_prot_<uniprot>`` converts
      ``prot_pool`` into one enzyme at a stoichiometric cost of that enzyme's molecular
      weight, and the pool exchange is capped at the measured total protein budget. The
      most of that enzyme the model can ever hold is therefore ``pool_cap / MW``, in
      closed form, from the model's own stoichiometry. At scale 1.0 the bound this gives
      is exactly non-binding, which is what makes an unstressed state a no-op.
    * **Any reaction with a finite upper bound.** Its own bound.

    Everything else -- and in ecYeastGEM that is every metabolic reaction, all 5746 of
    them carry ``(0, inf)`` -- has no finite reference and returns ``None``. Those are
    reported as skipped, not silently scaled: scaling infinity yields infinity, so a layer
    that pretended to constrain them would report a reach it does not have.

    **How loose this is, measured.** ``pool_cap / MW`` is the model's own hard ceiling, and
    it is a ceiling on one enzyme holding the *entire* protein budget. Since roughly a
    thousand enzymes share that budget, it over-states any single enzyme's realistic
    capacity by more than an order of magnitude: at maximum growth on ecYeastGEM the
    busiest enzyme in the H2O2 regulons draws 5.8% of its own implied ceiling and the
    median draws none of it. A relative bound against this reference is therefore inert
    until the scale falls below a few percent. That is not a reason to invent a tighter
    constant -- it is a reason for a caller to pass a *computed* one through
    ``capacities``, and for any report to say which reference it used.

    The looseness is not an artefact of this arithmetic either. Elsemman 2022, PMID
    35145105 shows the proteome constraints that actually drive yeast's metabolic strategy
    are **compartment-specific**, and Dinh & Maranas 2023, PMID 37080482 finds the Crabtree
    effect limited by *mitochondrial* proteome capacity rather than by the overall budget.
    A single global pool is, on that evidence, the wrong constraint to scale -- so a bound
    derived from it being non-binding is what the literature would predict.

    Returns:
        ``(forward, reverse)`` capacities as non-negative magnitudes, or ``None`` when the
        model implies no finite forward capacity. ``reverse`` is zero for an irreversible
        reaction, which is every enzyme-constrained one -- GECKO splits reversible
        reactions, so nothing there carries a negative bound to scale.
    """
    if pool is not None:
        cap, pool_id = pool
        for metabolite, coefficient in reaction.metabolites.items():
            if metabolite.id == pool_id and coefficient < 0:
                return cap / abs(coefficient), 0.0
    if reaction.upper_bound < float("inf"):
        reverse = -reaction.lower_bound if reaction.lower_bound > -float("inf") else 0.0
        return float(reaction.upper_bound), float(max(reverse, 0.0))
    return None


def override_capacity(
    reaction: cobra.Reaction,
    value: float | tuple[float, float],
) -> tuple[float, float]:
    """Split one entry of :func:`eflux_layer`'s ``capacities`` into forward and reverse.

    **This exists because the override used to return ``(value, 0.0)``, and that is a live
    correctness bug.** Zero reverse capacity is not "the caller did not say"; it is the
    assertion that the reaction cannot run backwards at all. On yeast-GEM 9.0.2 a caller
    handing back the layer's *own* forward capacities -- the identity, which must be a no-op
    -- locked **92 reversible reactions** forward-only and took max growth from 0.08584 /h to
    **0.00000 /h** at a module activity of 1e-12, i.e. with no regulation applied whatsoever.
    Measured 2026-09-05; `tests/test_eflux_capacities.py` pins both halves. ecYeastGEM is
    untouched by it because GECKO splits reversible reactions, so nothing there carries a
    negative bound to lose.

    A bare float overrides the FORWARD capacity only. The reverse then comes from the same
    place it comes from with no override at all -- the reaction's own lower bound -- so
    overriding a capacity can no longer change a reaction's *direction*, only its magnitude.
    A caller with both directions measured (FVA under the regime being studied gives both)
    passes a ``(forward, reverse)`` pair and states them.

    Returns:
        ``(forward, reverse)`` as non-negative magnitudes, in the same convention
        :func:`reference_capacity` returns.

    Raises:
        ValueError: on a negative capacity, which is not a quantity; or when the reverse
            capacity would be infinite, which is no reference at all -- pass the pair.
    """
    if isinstance(value, tuple):
        forward, reverse = float(value[0]), float(value[1])
    else:
        forward = float(value)
        low = float(reaction.lower_bound)
        reverse = -low if low < 0.0 else 0.0
    if forward < 0.0 or reverse < 0.0:
        raise ValueError(
            f"capacity for {reaction.id!r} must be a non-negative magnitude, "
            f"got forward={forward}, reverse={reverse}")
    if reverse == float("inf"):
        raise ValueError(
            f"{reaction.id!r} has an infinite reverse bound, so the model implies no "
            "reverse reference capacity to scale; pass a (forward, reverse) pair to say "
            "what it is. Substituting zero would silently make the reaction irreversible.")
    return forward, reverse


# --------------------------------------------------------------------------- #
# E-Flux
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ReactionScale:
    """One reaction's E-Flux bounds, and the arithmetic that produced them.

    ``lower_bound`` is negative or zero. A reversible reaction's reverse capacity scales
    by the same factor, because E-Flux bounds ``|v|``; an irreversible one has
    ``reverse_capacity == 0.0`` and its lower bound stays at zero.
    """

    reaction: str
    scale: float
    reference_capacity: float
    reverse_capacity: float
    genes: tuple[str, ...]

    @property
    def upper_bound(self) -> float:
        return self.scale * self.reference_capacity

    @property
    def lower_bound(self) -> float:
        return -self.scale * self.reverse_capacity


@dataclass(frozen=True)
class EFluxLayer:
    """Bounds E-Flux would impose, computed but not applied. Upper bounds only.

    Deliberately inert. The layer is a value so it can be inspected, counted and diffed
    before anything solves with it, and so that applying it twice is an assignment rather
    than a compounding multiplication -- :func:`apply_eflux` never reads a bound it may
    itself have written, which is what makes it idempotent.

    ``lifts_lower_bound`` is always ``False``. It is carried on the object because a
    caller asking "does regulation narrow the product range" needs it, and reading it off
    a field is safer than remembering it.

    ``layer_state`` is always ``"inert"``, for the reason and the measurement in
    :data:`EFLUX_LAYER_STATE`, and it is carried and printed for the same reason: a report
    should read the label off the layer rather than paraphrase it. ``binding`` is the
    quantity it is measured on, and on this panel it is empty at every dose.

    ``modules_unreached`` names modules that carry activity and a non-empty regulon and
    still reach no reaction in *this* model, because every gene in the regulon is absent
    from it. :func:`gene_scales` cannot refuse those -- it has not seen the model -- and
    raising here would make a full ``module_response`` dict unusable, since most stress
    regulons are not metabolic and that is normal rather than wrong. So they are named
    instead, and named in :meth:`summary`. UPR is the one that matters: Hac1's
    conserved-motif regulon is seven genes and the GSMM carries none of them, so the module
    two of the four wet-lab constructs report on cannot constrain metabolism at all.
    """

    scales: tuple[ReactionScale, ...]
    genes_reached: tuple[str, ...]
    genes_absent: tuple[str, ...]
    modules_excluded: tuple[str, ...]
    modules_unreached: tuple[str, ...]
    unbounded_skipped: int
    reference: str = "model-implied"
    lifts_lower_bound: bool = EFLUX_LIFTS_LOWER_BOUND
    layer_state: str = EFLUX_LAYER_STATE
    citation: str = f"E-Flux, Colijn 2009, PMID {COLIJN_PMID}"

    @property
    def binding(self) -> tuple[ReactionScale, ...]:
        """Reactions whose bound this layer actually tightens."""
        return tuple(s for s in self.scales if s.scale < 1.0)

    def summary(self) -> str:
        unreached = (f"; reached no reaction: {list(self.modules_unreached)}"
                     if self.modules_unreached else "")
        return (
            f"[E-Flux, {self.layer_state.upper()}] {len(self.scales)} reactions bounded "
            f"({len(self.binding)} tightened), {len(self.genes_reached)} regulon genes in "
            f"the model, {len(self.genes_absent)} absent, {self.unbounded_skipped} "
            f"reactions skipped for having no finite reference capacity; reference "
            f"{self.reference}; no lower bound raised above zero{unreached}"
        )


def eflux_layer(
    model: cobra.Model,
    activity: dict[str, float],
    regulons: dict[str, Regulon] | None = None,
    allow_unmapped: frozenset[str] | set[str] = frozenset(),
    capacities: dict[str, float | tuple[float, float]] | None = None,
) -> EFluxLayer:
    """Compute the E-Flux bounds a latent module state implies. Does not modify ``model``.

    The gene-to-reaction step is the GEM's own ``gene_reaction_rule``, read through
    ``model.genes``, and nothing here curates it: which reactions a gene controls is a
    property of the model and hand-listing it would put a second, unreviewed copy of the
    model's annotation into this file.

    Args:
        model: Model to compute against. Read only.
        activity: Signed module activity, e.g. from
            ``generator/stress_panel.py::module_response``.
        regulons: Output of :func:`load_regulons`. Loaded from the vendored SGD records
            when omitted.
        allow_unmapped: Modules accepted as unregulatable. See :func:`gene_scales`.
        capacities: Per-reaction reference capacity, overriding
            :func:`reference_capacity`. Exists because the choice of reference is the
            load-bearing decision in this layer and the model's implied ceiling is a very
            loose one -- so a caller with a computed alternative (an FVA maximum under the
            regime being studied, say) supplies it here, where the choice is visible in
            the call rather than buried as a default. A reaction listed here is bounded
            even if the model implies no finite capacity for it. A bare float overrides the
            forward capacity and leaves the reverse where the model puts it; a
            ``(forward, reverse)`` pair states both. See :func:`override_capacity` for what
            this used to do instead and what it cost.

    Raises:
        ValueError: from :func:`gene_scales`, naming any module with activity and no
            regulon; or from :func:`override_capacity` on a bad ``capacities`` entry.
    """
    override = {} if capacities is None else capacities
    table = load_regulons() if regulons is None else regulons
    scales = gene_scales(activity, table, allow_unmapped=allow_unmapped)

    present = {g.id for g in model.genes}
    reached = sorted(g for g in scales if g in present)
    absent = sorted(g for g in scales if g not in present)

    pool = _protein_pool_cap(model)
    candidates: set[str] = set()
    for gene in reached:
        candidates.update(r.id for r in model.genes.get_by_id(gene).reactions)

    unreached = tuple(sorted(
        module for module, value in activity.items()
        if value != 0.0 and module not in allow_unmapped
        and (regulon := table.get(module)) is not None and regulon.genes
        and not any(gene in present and model.genes.get_by_id(gene).reactions
                    for gene in regulon.genes)
    ))

    rows: list[ReactionScale] = []
    skipped = 0
    for rid in sorted(candidates):
        reaction = model.reactions.get_by_id(rid)
        if rid in override:
            forward, reverse = override_capacity(reaction, override[rid])
        else:
            capacity = reference_capacity(reaction, pool)
            if capacity is None:
                skipped += 1
                continue
            forward, reverse = capacity
        scale = _rule_scale(getattr(reaction.gpr, "body", None), scales)
        rows.append(ReactionScale(
            reaction=rid, scale=float(scale), reference_capacity=float(forward),
            reverse_capacity=float(reverse),
            genes=tuple(sorted(g.id for g in reaction.genes if g.id in scales)),
        ))
    return EFluxLayer(
        scales=tuple(rows), genes_reached=tuple(reached), genes_absent=tuple(absent),
        modules_excluded=tuple(sorted(m for m in allow_unmapped if activity.get(m, 0.0))),
        modules_unreached=unreached, unbounded_skipped=skipped,
        reference="caller-supplied capacities" if override else "model-implied",
    )


@dataclass(frozen=True)
class AppliedRegulation:
    """What applying a layer actually changed, and whether the model still grows."""

    bounds_changed: int
    bounds_tightened: int
    max_growth: float | None
    mechanism: str
    lifted_lower_bounds: int = 0
    status: str = "optimal"

    @property
    def feasible(self) -> bool:
        return self.status == "optimal"

    def summary(self) -> str:
        growth = "unavailable" if self.max_growth is None else f"{self.max_growth:.5f} /h"
        return (
            f"[{self.mechanism}] {self.bounds_changed} bounds set "
            f"({self.bounds_tightened} tightened, {self.lifted_lower_bounds} lower bounds "
            f"lifted); status {self.status}; max growth {growth}"
        )


def apply_eflux(
    model: cobra.Model,
    layer: EFluxLayer,
    biomass_reaction: str = "r_2111",
    acknowledge_infeasible: bool = False,
) -> AppliedRegulation:
    """Impose ``layer``'s capacity bounds on ``model`` **in place**, and check it still grows.

    Use inside a ``with model:`` block; this mutates. Bounds move as ``|v| <= scale * c``,
    so a reversible reaction's negative bound is scaled too -- and no bound ever becomes
    positive, which is why :data:`EFLUX_LIFTS_LOWER_BOUND` is ``False``.

    Applying the same layer twice is a no-op, because every bound is assigned from the
    layer's own recorded reference capacity rather than multiplied into whatever bound the
    reaction currently has.

    Args:
        model: Model to constrain in place.
        layer: Output of :func:`eflux_layer`.
        biomass_reaction: Growth reaction, used only for the feasibility check.
        acknowledge_infeasible: Return the report instead of raising when the model can no
            longer grow. Set it only when infeasibility is the thing being measured.

    Raises:
        InfeasibleRegulation: if the model can no longer grow, naming the tightest bounds.
    """
    changed = tightened = 0
    for row in layer.scales:
        reaction = model.reactions.get_by_id(row.reaction)
        new = (max(row.lower_bound, reaction.lower_bound) if reaction.lower_bound > 0 else row.lower_bound,
               min(row.upper_bound, reaction.upper_bound) if reaction.upper_bound < 0 else row.upper_bound)
        if new[0] > new[1]:
            name = _constraint_name("regulation_capacity_", {reaction.id: 1.0})
            existing = model.constraints.get(name)
            if existing is None or (existing.lb, existing.ub) != (row.lower_bound, row.upper_bound):
                if existing is not None:
                    model.remove_cons_vars([existing])
                model.add_cons_vars([model.problem.Constraint(reaction.flux_expression, lb=row.lower_bound,
                                                              ub=row.upper_bound, name=name)])
                changed += 1
                tightened += 1
            continue
        if new[1] < reaction.upper_bound or new[0] > reaction.lower_bound:
            tightened += 1
        if new != reaction.bounds:
            changed += 1
        reaction.bounds = new

    model.objective = biomass_reaction
    model.objective_direction = "max"
    # An infeasible LP returns nan and `nan <= 0.0` is False, so reading this with
    # `slim_optimize` made the refusal below unreachable -- the exact case it exists for.
    growth = growth_or_none(model)
    status = _solve_status(model, growth)
    if status != "optimal" and not acknowledge_infeasible:
        worst = sorted(layer.binding, key=lambda s: s.scale)[:5]
        raise InfeasibleRegulation(
            "this regulation layer leaves the model unable to grow "
            f"(max growth {growth!r}). Tightest bounds: "
            + ", ".join(f"{s.reaction} x{s.scale:.3f}" for s in worst)
            + ". Either the activities are past what the network tolerates or a regulon "
            "reaches an essential reaction; pass acknowledge_infeasible=True to measure "
            "that rather than be stopped by it.", status=status,
        )
    return AppliedRegulation(
        bounds_changed=changed, bounds_tightened=tightened,
        max_growth=float(growth) if status == "optimal" else None,
        mechanism=f"E-Flux (PMID {COLIJN_PMID})", status=status,
    )


# --------------------------------------------------------------------------- #
# lexicographic task efficiency -- the part that can lift a floor
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MetabolicTask:
    """One sourced flux measurement for a COSMIC-inspired task-consistency calculation.

    Args:
        name: Human name -- "biomass", "beta-carotene", "ethanol".
        reaction: Reaction id carrying the task.
        measured_flux: The **measured** specific rate, mmol/gDW/h (or 1/h for biomass).
            ``None`` means unmeasured, and every function that would use it refuses.
        source: Where the measurement came from. Required when there is one: a measured
            flux with no provenance is indistinguishable from an assumed one, and the
            whole method rests on the difference.
    """

    name: str
    reaction: str
    measured_flux: float | None = None
    source: str = ""
    unit: str = "mmol/gDW/h"
    direction: str = "forward"
    condition_id: str = ""
    host: str = ""

    def coefficients(self, model: cobra.Model) -> dict[str, float]:
        if self.direction == "forward":
            model.reactions.get_by_id(self.reaction)
            return {self.reaction: 1.0}
        return exchange_coefficients(model, self.reaction, direction=self.direction)

    def require_measurement(self) -> float:
        """The measured flux, or a refusal naming what to measure.

        Raises:
            MissingMeasurement: naming the task, the reaction and the units.
        """
        if self.measured_flux is None:
            raise MissingMeasurement(
                f"task {self.name!r} ({self.reaction}) has no measured flux. The "
                "task-efficiency ordering is measured_flux / fba_maximum, so without the "
                "numerator there is no efficiency and no priority. Supply the specific "
                "rate in mmol/gDW/h from a measured uptake/secretion balance. It is NOT "
                "substitutable by a simulated flux: the FBA maximum is the denominator, "
                "so using a model flux as the numerator makes every efficiency 1.0 and "
                "the resulting lower bound an artefact of the substitution."
            )
        if not isinstance(self.source, str) or not self.source.strip():
            raise MissingMeasurement(
                f"task {self.name!r} carries a flux of {self.measured_flux} with no "
                "source. Record where it was measured; an unsourced number here is "
                "indistinguishable from an assumption, and this method's only claim to "
                "narrowing anything is that the numerator is a measurement."
            )
        return specific_rate(self.measured_flux, self.unit)


@dataclass(frozen=True)
class TaskPriority:
    """One task's efficiency and its rank in the lexicographic order."""

    rank: int
    name: str
    reaction: str
    measured_flux: float
    fba_maximum: float
    efficiency: float | None
    coefficients: tuple[tuple[str, float], ...] = ()
    unit: str = "mmol/gDW/h"
    source: str = ""
    condition_id: str = ""
    host: str = ""
    numerical_ties: tuple[str, ...] = ()

    def summary(self) -> str:
        efficiency = "undefined (0/0)" if self.efficiency is None else f"{self.efficiency:.3f}"
        ties = f"; numerical tie broken by name: {self.numerical_ties}" if self.numerical_ties else ""
        return (
            f"{self.rank}. {self.name} ({self.reaction}): measured "
            f"{self.measured_flux:.4g} / max {self.fba_maximum:.4g} = efficiency "
            f"{efficiency}; lower bound only{ties}"
        )


def task_priority_order(
    model: cobra.Model,
    tasks: list[MetabolicTask] | tuple[MetabolicTask, ...], *,
    condition_id: str | None = None, host: str | None = None,
) -> list[TaskPriority]:
    """Rank measured-over-maximum efficiencies, imposing only measured lower bounds.

    This is a COSMIC-inspired consistency adaptation of the lower-bound procedure now
    verified in paper section 4.3 and its supplementary algorithm, not a numerical replay
    of the original CHO study. The public kernel consumes precomputed priorities and
    fractions of successive maxima; see data/cosmic_sources.json.
    At each round the largest defined ratio is selected and its measurement becomes a
    floor, leaving the upper bound untouched. A zero/zero ratio is explicitly undefined.
    Numerical ties are reported, with name order used only as a deterministic convention.

    ``model`` is modified inside a scoped context and left as it was found. Both glucose
    and oxygen must already have their condition-matched measured bounds from the caller;
    a different nutrient regime makes the denominators incomparable with measurements.

    Args:
        model: Model with the task reactions present, uptakes already constrained.
        tasks: Candidate tasks. Every one must carry a measured flux and a source.

    Raises:
        MissingMeasurement: naming the first task without a measurement or a source. This
            is the expected outcome in this repository today -- no strain carries the
            beta-carotene pathway, so no product secretion rate has been measured.
        ValueError: if a task's measured flux exceeds what the model can do at all, which
            means the model and the measurement disagree and no efficiency below 1 can be
            read as a priority.
    """
    if not tasks:
        raise ValueError("no tasks given; a priority order over nothing is not a result")
    remaining = list(tasks)
    seen_names, seen_expressions = set(), set()
    for task in remaining:
        task.require_measurement()  # refuse before solving anything
        coefficients = tuple(task.coefficients(model).items())
        if task.name in seen_names or coefficients in seen_expressions:
            raise ValueError("duplicate task name or flux expression")
        seen_names.add(task.name)
        seen_expressions.add(coefficients)
        if condition_id is not None and task.condition_id != condition_id:
            raise ValueError(f"task {task.name!r} has an unmatched measurement condition")
        if host is not None and task.host != host:
            raise ValueError(f"task {task.name!r} has an incompatible measurement host")
    if len({task.condition_id for task in remaining}) > 1:
        raise ValueError("task measurements mix condition identifiers")
    if len({task.host for task in remaining}) > 1:
        raise ValueError("task measurements mix host identities")

    order: list[TaskPriority] = []
    with model as scoped:
        while remaining:
            scored: list[TaskPriority] = []
            for task in remaining:
                coefficients = task.coefficients(scoped)
                scoped.objective = scoped.problem.Objective(_expression(scoped, coefficients), direction="max")
                maximum = growth_or_none(scoped)  # nan <= 0.0 is False; see above
                measured = task.require_measurement()
                status = _solve_status(scoped, maximum)
                if status != "optimal":
                    raise InfeasibleRegulation(
                        f"task {task.name!r} ({task.reaction}): solver status {status}; "
                        "an unsolved maximum is not zero growth or a task efficiency", status=status,
                    )
                if measured > maximum + max(1e-9, scoped.tolerance):
                    raise InfeasibleRegulation(
                        f"task {task.name!r} ({task.reaction}) was measured at "
                        f"{measured:.4g} but the model's maximum is {maximum:.4g}. "
                        "Efficiency above 1 is not a priority, it is a contradiction "
                        "between the model and the measurement."
                    )
                scored.append(TaskPriority(
                    rank=0, name=task.name, reaction=task.reaction,
                    measured_flux=measured, fba_maximum=float(maximum),
                    efficiency=float(measured / maximum) if maximum > 0 else None,
                    coefficients=tuple(coefficients.items()),
                    unit="1/h" if task.unit == "1/h" else "mmol/gDW/h",
                    source=task.source, condition_id=task.condition_id, host=task.host,
                ))
            scored.sort(key=lambda row: row.name)
            best = max(scored, key=lambda row: -math.inf if row.efficiency is None else row.efficiency)
            ties = tuple(row.name for row in scored if row.name != best.name and (
                row.efficiency == best.efficiency or (
                    row.efficiency is not None and best.efficiency is not None
                    and abs(row.efficiency - best.efficiency) <= 1e-9)))
            order.append(TaskPriority(**{**vars(best), "rank": len(order) + 1, "numerical_ties": ties}))
            _set_floor(scoped, dict(best.coefficients), best.measured_flux)
            remaining = [task for task in remaining if task.name != best.name]
        scoped.objective = scoped.problem.Objective(0, direction="max")
        value = growth_or_none(scoped)
        status = _solve_status(scoped, value)
        if status != "optimal":
            raise InfeasibleRegulation(f"the complete measured task set is {status}; no task was dropped", status=status)
    return order


def apply_task_bounds(
    model: cobra.Model,
    order: list[TaskPriority] | tuple[TaskPriority, ...],
    up_to_rank: int | None = None, *, biomass_reaction: str | None = None,
    acknowledge_infeasible: bool = False,
) -> AppliedRegulation:
    """Impose each ranked task's measured lower bound, **in place**, in priority order.

    Use inside a ``with model:`` block. Unlike :func:`apply_eflux` this does move lower
    bounds, without setting an equality or relaxing a pre-existing stricter floor.
    A task can exceed its measurement if another task or mass balance requires it.
    The feasibility check explicitly maximises biomass and preserves the prior objective.

    Args:
        model: Model to constrain in place.
        order: Output of :func:`task_priority_order`.
        up_to_rank: Bound only the tasks at or above this rank, leaving the rest free. Use
            it to ask what the higher-priority tasks alone imply about a product whose own
            rate has not been measured.
    """
    if biomass_reaction is None:
        objective = cobra.util.solver.linear_reaction_coefficients(model)
        if len(objective) != 1 or next(iter(objective.values())) != 1.0:
            raise ValueError("name the biomass_reaction explicitly; the current objective is not one growth flux")
        biomass_reaction = next(iter(objective)).id
    model.reactions.get_by_id(biomass_reaction)
    if up_to_rank is not None and (isinstance(up_to_rank, bool) or not isinstance(up_to_rank, int) or up_to_rank < 1):
        raise ValueError("up_to_rank must be a positive integer")
    lifted = changed = 0
    for priority in sorted(order, key=lambda p: p.rank):
        if up_to_rank is not None and priority.rank > up_to_rank:
            continue
        amount = specific_rate(priority.measured_flux, priority.unit)
        coefficients = dict(priority.coefficients) or {priority.reaction: 1.0}
        count = _set_floor(model, coefficients, amount)
        changed += count
        lifted += count
    with model as scoped:
        scoped.objective = biomass_reaction
        scoped.objective_direction = "max"
        growth = growth_or_none(model)  # else an infeasible solve is reported as max_growth=nan
        status = _solve_status(scoped, growth)
    if status != "optimal" and not acknowledge_infeasible:
        raise InfeasibleRegulation(f"measured task lower bounds: solver status {status}, not zero growth", status=status)
    return AppliedRegulation(
        bounds_changed=changed, bounds_tightened=changed, lifted_lower_bounds=lifted,
        max_growth=float(growth) if status == "optimal" else None, status=status,
        mechanism=f"measured task-efficiency lower bounds (COSMIC-inspired, PMID {COSMIC_PMID})",
    )
