"""A number, its units, where it came from, and whether it is allowed to be a number at all.

Seven tags, and the tag is part of the specification rather than commentary. Six of them
carry a value; one refuses to. Two of the six refuse to collapse to a single float.

    MEASURED   a number in a named *S. cerevisiae* paper, with units
    BOUNDED    a point chosen inside a bracket whose floor and ceiling are both measured
    BORROWED   measured, but not in this organism; the organism travels with the number
    DERIVED    arithmetic on MEASURED rows, and the arithmetic is shown
    ASSERTED   a number a source states without evidence, or a modelling convention
    SWEPT      a declared sensitivity axis; the model reports a band, not a point
    REFUSED    no defensible value exists, and the code raises rather than defaulting

WHY THE TAG IS MANDATORY AND NOT A DEFAULT. Seven classifier sweeps of this package,
adjudicated by re-running the code, graded its 906 substantive literals: **591 (65%) are
ASSERTED -- chosen because the arithmetic needed a value -- against 125 (14%) MEASURED and
65 (7%) DERIVED, of which 40 are unit conversions.** Every flattering grade that was tested
failed more often than it held: a defined SI constant graded MEASURED, a typed literal
graded DERIVED because a docstring explained how it *could* be computed, a chosen
experimental dose graded BORROWED because a paper had used it. A tag with a default would
be the same mistake with a type annotation on it, so :class:`Param` has no default tag and
no positional-optional escape: a value arriving here without a grade is a construction
error, in the same way a pathway node arriving without a rate law is.

WHAT float() DOES, WHICH IS THE WHOLE POINT.

* MEASURED, BOUNDED, BORROWED, DERIVED, ASSERTED -- returns the number. ASSERTED is
  allowed *because it is flagged*: the tag travels with the value into
  :func:`provenance_table` and into the free-scalar count, which is what "never silently"
  means operationally. Refusing it would be honest and would also make 65% of this
  repository unconstructible.
* SWEPT -- raises :class:`SweptValue`. A band is not a point, and the failure this
  prevents is specific: `ARCHITECTURE_GAPS.md` 0.3 found that every emergence test in the
  target architecture is scored on whether *some point in a sweep* reproduces the target,
  "which is fitting with extra steps". Taking a midpoint silently is the first step of
  that. :meth:`Param.grid` and :meth:`Param.at` are the sanctioned ways through, and both
  put the choice on the record.
* REFUSED -- raises :class:`RefusedValue` naming the reason, the measurement that would
  close it, and the citation. This is the existing house idiom, not a new one:
  `pathway/solve.py`'s ``require_saturating`` refuses a node one constant short by name,
  and `kinetic/carotenoid.py`'s ``calibrated_kinetics`` refuses the unidentified half of a
  branch and says which HPLC peak would fix it.

THE FREE-SCALAR GATE (criterion (e)). :class:`ParamRegistry` counts a piece's free scalars
against its independent targets and refuses the piece where free > targets. A scalar is
FREE when no single measurement pins it to a point -- BOUNDED, ASSERTED, SWEPT and REFUSED
-- and PINNED when one does: MEASURED and DERIVED pin a point in this organism, BORROWED
pins a point in another one and carries the flag saying so. The calibration for that rule
is `ARCHITECTURE_GAPS.md` 0.3, which counts the pH slice at **eight free scalars against
two usable targets**: beta, Jmax_pump, pKa_pma1, n, g_H and V_CYT_L_PER_GDCW all SWEPT,
Vmax_ex and K_ex both REFUSED. :data:`FREE_TAGS` reproduces that eight exactly, and
``tests/test_mech_params.py`` encodes the slice so the rule cannot drift away from the
document that motivated it.

BOUNDED IS FREE, AND THAT IS THE DEBATABLE ONE. Its value is constrained by something
measurable, so it passes criterion (b); it can still move inside its bracket without
contradicting any measurement, so it is a degree of freedom and criterion (e) counts it.
The audit's own worked example is `stress_panel.py`'s H2O2 target EC50 = 0.05 mM, a
bracket between a measured 20 uM floor (Kritsiligkou 2021, PMID 34118234) and a 0.15 mM
regulon EC50, with no HyPer7 half-max on record for yeast. Nothing measured moves if that
0.05 becomes 0.08.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ystwin import citations

__all__ = [
    "FREE_TAGS",
    "FreeScalarGate",
    "FreeScalarGateFailed",
    "NoSingleValue",
    "PINNED_TAGS",
    "Param",
    "ParamRegistry",
    "ParameterUncertainty",
    "RefusedValue",
    "SweptValue",
    "Tag",
    "Target",
    "provenance_table",
    "tag_census_table",
]


class Tag:
    """The seven grades. A plain string constant class, as :class:`ystwin.pathway.spec.Fate`
    is, so a tag survives a round trip through a TSV or a markdown table unchanged.

    ``BORROWED`` is what `ARCHITECTURE_TARGET.md` writes as MEASURED-dagger: measured, in a
    named organism that is not *S. cerevisiae*, imported only with the flag attached. The
    flag is attached here by construction -- :class:`Param` refuses a BORROWED row with no
    organism, and refuses one whose organism is the host, because that is a MEASURED row
    that has been graded down by mistake.
    """

    MEASURED = "MEASURED"
    BOUNDED = "BOUNDED"
    BORROWED = "BORROWED"
    DERIVED = "DERIVED"
    ASSERTED = "ASSERTED"
    SWEPT = "SWEPT"
    REFUSED = "REFUSED"
    ALL = (MEASURED, BOUNDED, BORROWED, DERIVED, ASSERTED, SWEPT, REFUSED)


PINNED_TAGS = (Tag.MEASURED, Tag.DERIVED, Tag.BORROWED)
"""Grades where one measurement fixes the number to a point. See the module docstring."""

FREE_TAGS = (Tag.BOUNDED, Tag.ASSERTED, Tag.SWEPT, Tag.REFUSED)
"""Grades that cost a degree of freedom. Reproduces `ARCHITECTURE_GAPS.md` 0.3's eight."""

_HOST_ORGANISM_NAMES = (
    "s. cerevisiae", "s.cerevisiae", "saccharomyces cerevisiae", "scerevisiae", "yeast",
)


class NoSingleValue(NotImplementedError):
    """``float()`` was asked for a number this parameter does not have.

    ``NotImplementedError`` rather than ``ValueError`` because that is what
    ``kinetic/carotenoid.py``'s ``calibrated_kinetics`` raises for the same act -- being
    asked for a parameter set the measurements do not identify -- and a caller that already
    catches one should catch the other.
    """


class RefusedValue(NoSingleValue):
    """A REFUSED parameter has no defensible value and will not pretend to have one."""


class SweptValue(NoSingleValue):
    """A SWEPT parameter is a band. Use :meth:`Param.grid` or :meth:`Param.at`."""


class FreeScalarGateFailed(ValueError):
    """A piece declares more free scalars than it has independent targets."""


def _finite(value, label: str, name: str) -> float:
    if value is None:
        raise ValueError(f"parameter {name!r}: {label} is required and got None")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"parameter {name!r}: {label} must be finite, got {number}")
    return number


@dataclass(frozen=True)
class Param:
    """One constant, graded.

    The positional signature is the one `ARCHITECTURE_TARGET.md` Phase 0 specifies --
    ``Param(value, units, tag, source)`` -- and everything after it is keyword-shaped
    because it is grade-specific. The seven classmethod constructors below are the
    preferred spelling: they cannot mistype a tag, and they require exactly the fields
    their grade needs.

    Args:
        value: The number, in :attr:`units`. **Must be None for SWEPT and REFUSED**, which
            is not a formality -- a swept parameter carrying a value has an attribute a
            caller can read past the refusal, and the refusal is then decorative.
        units: Required, always, including the literal ``"dimensionless"``. Three
            published time bases collide in this architecture (min^-1, s^-1, h^-1) and the
            resulting 60x has already caused errors in the source literature.
        tag: One of :class:`Tag`. No default. See the module docstring.
        source: Where the number, or its absence, comes from. Required for every grade: a
            REFUSED parameter still has to say what was searched and came up empty.
        name: Used in messages, tables and the registry. :meth:`ParamRegistry.add`
            requires it; construction does not, so the positional form above stays legal.
        bounds: ``(low, high)``. Required for BOUNDED (the measured bracket the value sits
            inside) and for SWEPT (the axis). Forbidden elsewhere.
        ci95: A 95% interval around a MEASURED or DERIVED value, containing it. Optional,
            and worth carrying wherever the source reports one.
        organism: Required for BORROWED, forbidden elsewhere.
        reason: Why no defensible value exists. Required for REFUSED.
        missing: The measurement that would close it. Required for REFUSED and for BOUNDED
            -- the audit's definition of BOUNDED is a grade that "names the floor, the
            ceiling and the absent measurement", and two out of three is a bracket with no
            way out of it.
    """

    value: float | None
    units: str
    tag: str
    source: str
    name: str = ""
    bounds: tuple[float, float] | None = None
    ci95: tuple[float, float] | None = None
    organism: str = ""
    reason: str = ""
    missing: str = ""

    def __post_init__(self) -> None:
        label = self.name or "<unnamed>"
        if self.tag not in Tag.ALL:
            raise ValueError(
                f"parameter {label!r} has tag {self.tag!r}, which is not one of "
                f"{', '.join(Tag.ALL)}. The tag is mandatory and has no default: seven "
                "classifier sweeps of this package graded 65% of its substantive numbers "
                "ASSERTED against 14% MEASURED, and an ungraded number is how that "
                "happened. Pick the grade that is true, including ASSERTED")
        if not str(self.units).strip():
            raise ValueError(
                f"parameter {label!r}: units are required. Use 'dimensionless' if it has "
                "none -- three published time bases collide in this architecture and the "
                "60x has already caused errors upstream")
        if not str(self.source).strip():
            raise ValueError(
                f"parameter {label!r}: source is required for every grade, including "
                "REFUSED, which still has to say what was searched and came up empty")

        if self.tag in (Tag.SWEPT, Tag.REFUSED):
            if self.value is not None:
                raise ValueError(
                    f"parameter {label!r} is {self.tag} and carries value {self.value!r}. "
                    "A value here is readable past the refusal, which makes the refusal "
                    "decorative; put a SWEPT parameter's range in bounds")
        else:
            object.__setattr__(self, "value", _finite(self.value, "value", label))

        if self.tag in (Tag.BOUNDED, Tag.SWEPT):
            if self.bounds is None:
                raise ValueError(
                    f"parameter {label!r} is {self.tag} and needs bounds=(low, high): "
                    + ("the measured bracket its value sits inside"
                       if self.tag == Tag.BOUNDED else "the axis it is swept over"))
            low = _finite(self.bounds[0], "bounds[0]", label)
            high = _finite(self.bounds[1], "bounds[1]", label)
            if not low < high:
                raise ValueError(
                    f"parameter {label!r}: bounds must be increasing, got ({low}, {high})")
            object.__setattr__(self, "bounds", (low, high))
            if self.tag == Tag.BOUNDED and not low <= self.value <= high:
                raise ValueError(
                    f"parameter {label!r} is BOUNDED at {self.value} {self.units} but its "
                    f"bracket is [{low}, {high}]. A value outside its own bracket is an "
                    "ASSERTED value with a citation attached to something else")
        elif self.bounds is not None:
            raise ValueError(
                f"parameter {label!r} is {self.tag} and carries bounds. Only BOUNDED and "
                "SWEPT have a range; anything else is a point or a refusal")

        if self.ci95 is not None:
            low = _finite(self.ci95[0], "ci95[0]", label)
            high = _finite(self.ci95[1], "ci95[1]", label)
            if self.value is None or not low <= self.value <= high:
                raise ValueError(
                    f"parameter {label!r}: ci95 [{low}, {high}] does not contain its value "
                    f"{self.value!r}")
            object.__setattr__(self, "ci95", (low, high))

        if self.tag == Tag.BORROWED:
            organism = str(self.organism).strip()
            if not organism:
                raise ValueError(
                    f"parameter {label!r} is BORROWED and names no organism. The flag is "
                    "the whole content of the grade: import only with it attached")
            if organism.lower().rstrip(".") in _HOST_ORGANISM_NAMES:
                raise ValueError(
                    f"parameter {label!r} is BORROWED but its organism is {organism!r}, "
                    "which is the host. That is a MEASURED row graded down by mistake")
        elif str(self.organism).strip():
            raise ValueError(
                f"parameter {label!r} is {self.tag} and names an organism. Only BORROWED "
                "carries one, because only BORROWED is a claim about a different cell")

        if self.tag == Tag.REFUSED and not str(self.reason).strip():
            raise ValueError(
                f"parameter {label!r} is REFUSED with no reason. The refusal has to say "
                "what is missing at the point where a caller hits it")
        if self.tag != Tag.REFUSED and str(self.reason).strip():
            raise ValueError(
                f"parameter {label!r} is {self.tag} and carries a refusal reason")
        if self.tag in (Tag.REFUSED, Tag.BOUNDED) and not str(self.missing).strip():
            raise ValueError(
                f"parameter {label!r} is {self.tag} and does not name the measurement that "
                "would close it. BOUNDED names the floor, the ceiling AND the absent "
                "measurement; two out of three is a bracket with no way out")

    @classmethod
    def measured(cls, name: str, value: float, units: str, source: str, *,
                 ci95: tuple[float, float] | None = None) -> Param:
        """A point measured in *S. cerevisiae*, cited at the point of use."""
        return cls(value, units, Tag.MEASURED, source, name=name, ci95=ci95)

    @classmethod
    def bounded(cls, name: str, value: float, units: str, source: str, *,
                bounds: tuple[float, float], missing: str) -> Param:
        """A point inside a measured bracket, with the measurement that would pin it."""
        return cls(value, units, Tag.BOUNDED, source, name=name, bounds=bounds,
                   missing=missing)

    @classmethod
    def borrowed(cls, name: str, value: float, units: str, source: str, *,
                 organism: str, ci95: tuple[float, float] | None = None) -> Param:
        """Measured in a named organism that is not the host. `ARCHITECTURE_TARGET`'s
        MEASURED-dagger."""
        return cls(value, units, Tag.BORROWED, source, name=name, organism=organism,
                   ci95=ci95)

    @classmethod
    def derived(cls, name: str, value: float, units: str, source: str, *,
                ci95: tuple[float, float] | None = None) -> Param:
        """Arithmetic on measured rows. ``source`` shows the arithmetic, not just the paper."""
        return cls(value, units, Tag.DERIVED, source, name=name, ci95=ci95)

    @classmethod
    def asserted(cls, name: str, value: float, units: str, source: str, *,
                 missing: str = "") -> Param:
        """A number with no evidence behind it, said out loud. ``float()`` still works."""
        return cls(value, units, Tag.ASSERTED, source, name=name, missing=missing)

    @classmethod
    def swept(cls, name: str, units: str, source: str, *, bounds: tuple[float, float],
              missing: str = "") -> Param:
        """A declared sensitivity axis. ``float()`` raises; :meth:`grid` and :meth:`at` do not."""
        return cls(None, units, Tag.SWEPT, source, name=name, bounds=bounds, missing=missing)

    @classmethod
    def refused(cls, name: str, units: str, source: str, *, reason: str,
                missing: str) -> Param:
        """The sentinel. ``float()`` raises with the reason, the missing measurement and
        the citation."""
        return cls(None, units, Tag.REFUSED, source, name=name, reason=reason,
                   missing=missing)

    def __float__(self) -> float:
        label = self.name or "<unnamed>"
        if self.tag == Tag.REFUSED:
            raise RefusedValue(
                f"parameter {label!r} ({self.units}) is REFUSED and has no float. "
                f"Reason: {self.reason} What would close it: {self.missing} "
                f"Source: {self.source}")
        if self.tag == Tag.SWEPT:
            low, high = self.bounds
            raise SweptValue(
                f"parameter {label!r} ({self.units}) is SWEPT over [{low:g}, {high:g}] and "
                f"has no single float -- the model reports a band, not a point, and a "
                f"midpoint taken here is a fitted parameter nobody counted. Use "
                f"{label or 'param'}.sweep_points(n) to sweep it or .at(x) to pin one point on the "
                f"record. Source: {self.source}"
                + (f" What would collapse the band: {self.missing}" if self.missing else ""))
        return float(self.value)

    @property
    def is_free(self) -> bool:
        """Whether this parameter costs a degree of freedom under criterion (e)."""
        return self.tag in FREE_TAGS

    @property
    def identifiers(self) -> tuple[str, ...]:
        """Literature identifiers in :attr:`source`, found by :mod:`ystwin.citations`.

        The repository's own scanner rather than a regex written here, so a source this
        reports as uncited is uncited by the same definition ``tests/test_citations.py``
        uses. Empty is not an error -- a number measured on this project's own plates cites
        a plate, not a PMID -- but :meth:`ParamRegistry.uncited` makes the set visible.
        """
        text = str(self.source)
        found = citations.scan_text(text, [1] * len(text), text, "<param source>")
        return tuple(f"{c.kind}:{c.ident}" for c in found)

    def sweep_points(self, n: int) -> tuple[float, ...]:
        """``n`` evenly spaced points across a SWEPT band, endpoints included."""
        if self.tag != Tag.SWEPT:
            raise ValueError(
                f"parameter {self.name or '<unnamed>'!r} is {self.tag}, not SWEPT; a grid "
                "over a point is a sweep that reports one answer and calls it a band")
        if n < 2:
            raise ValueError(f"a swept axis needs at least both endpoints, got n={n}")
        low, high = self.bounds
        return tuple(low + (high - low) * i / (n - 1) for i in range(n))

    def at(self, value: float) -> float:
        """One point of a SWEPT band, checked to be inside it. The explicit way through."""
        if self.tag != Tag.SWEPT:
            raise ValueError(
                f"parameter {self.name or '<unnamed>'!r} is {self.tag}, not SWEPT")
        low, high = self.bounds
        number = _finite(value, "at()", self.name or "<unnamed>")
        if not low <= number <= high:
            raise ValueError(
                f"parameter {self.name or '<unnamed>'!r} is swept over [{low:g}, {high:g}] "
                f"and {number:g} is outside it")
        return number

    def display_value(self) -> str:
        """What goes in the value column of a provenance table."""
        if self.tag == Tag.REFUSED:
            return "--"
        if self.tag == Tag.SWEPT:
            return f"{self.bounds[0]:g} - {self.bounds[1]:g}"
        if self.ci95 is not None:
            return f"{self.value:g} [{self.ci95[0]:g}, {self.ci95[1]:g}]"
        if self.tag == Tag.BOUNDED:
            return f"{self.value:g} in [{self.bounds[0]:g}, {self.bounds[1]:g}]"
        return f"{self.value:g}"

    def provenance(self) -> str:
        """What goes in the source column: the citation, plus what the grade owes."""
        parts = [str(self.source)]
        if self.organism:
            parts.append(f"organism: {self.organism}")
        if self.reason:
            parts.append(f"refused: {self.reason}")
        if self.missing:
            parts.append(f"would close it: {self.missing}")
        return " -- ".join(parts)


@dataclass(frozen=True)
class Target:
    """One independent measurement a piece is scored against. The denominator of the gate.

    Args:
        name: Unique within a registry. Two targets with one name are one target counted
            twice, which is the arithmetic the gate exists to prevent.
        observable: What is measured, in words.
        assay: The instrument or protocol. Criterion (a) requires the assay to be NAMED for
            any observable, rather than the plate-reader CV being borrowed for, say, an FBA
            product ceiling -- `REVISED_BUILD_LIST.md` records that borrowing as the error
            that sank Route C2.
        noise_floor: The assay's own noise, in :attr:`units`. Measured, and the two this
            repository has are in `generator/panel_experiment.py`:
            ``OBSERVED_ACTIVITY_CV = 0.146`` and ``MEASURED_GROWTH_RATE_SE = 0.0117`` /h.
        units: Units of the floor.
        source: Where the target value and its floor come from.
        fitted: True when the model was allowed to search for the point that matches this
            target. Such a target does not count: `ARCHITECTURE_GAPS.md` 0.3's finding is
            that scoring on whether *some point in a sweep* reproduces a number "is fitting
            with extra steps", and a fitted target in the denominator lets a sweep pay for
            itself.
    """

    name: str
    observable: str
    assay: str
    noise_floor: float
    units: str
    source: str
    fitted: bool = False

    def __post_init__(self) -> None:
        for label, text in (("name", self.name), ("observable", self.observable),
                            ("assay", self.assay), ("units", self.units),
                            ("source", self.source)):
            if not str(text).strip():
                raise ValueError(f"target {self.name!r}: {label} is required")
        floor = _finite(self.noise_floor, "noise_floor", self.name)
        if floor <= 0:
            raise ValueError(
                f"target {self.name!r}: noise_floor must be positive, got {floor}. A target "
                "with no floor cannot be passed or failed by criterion (a)")


@dataclass(frozen=True)
class FreeScalarGate:
    """The result of criterion (e) for one piece: free scalars against independent targets."""

    piece: str
    free: tuple[str, ...]
    targets: tuple[str, ...]
    fitted_targets: tuple[str, ...] = ()

    @property
    def passes(self) -> bool:
        return len(self.free) <= len(self.targets)

    def summary(self) -> str:
        verdict = "PASSES" if self.passes else "REFUSED"
        line = (f"{self.piece}: {len(self.free)} free scalars against "
                f"{len(self.targets)} independent targets -- {verdict}")
        if self.fitted_targets:
            line += (f" ({len(self.fitted_targets)} further target(s) not counted, fitted: "
                     f"{', '.join(self.fitted_targets)})")
        return line

    def __str__(self) -> str:
        return self.summary()


class ParamRegistry:
    """Every constant one piece of the model declares, and what it is scored against.

    A model that can report its own free-scalar count is the difference between criterion
    (e) being a rule and being a paragraph. Registration is explicit -- nothing scrapes the
    module -- because a scraper would silently stop seeing a constant that moved into a
    default argument, and that is precisely where the audit found several of them.
    """

    def __init__(self, piece: str) -> None:
        if not str(piece).strip():
            raise ValueError("a registry needs the name of the piece it belongs to")
        self.piece = piece
        self._params: dict[str, Param] = {}
        self._targets: dict[str, Target] = {}

    def add(self, param: Param) -> Param:
        """Register one parameter and hand it back, so a module constant can be assigned
        from the call."""
        if not param.name.strip():
            raise ValueError(
                f"registry {self.piece!r}: a parameter needs a name to be registered; an "
                "anonymous row is invisible in the provenance table and uncountable in the "
                "gate")
        if param.name in self._params:
            raise ValueError(
                f"registry {self.piece!r}: {param.name!r} is already registered. Two "
                "declarations of one constant is how a value and its copy drift apart")
        self._params[param.name] = param
        return param

    def add_target(self, target: Target) -> Target:
        if target.name in self._targets:
            raise ValueError(
                f"registry {self.piece!r}: target {target.name!r} is already registered; "
                "one measurement counted twice is the arithmetic the gate prevents")
        self._targets[target.name] = target
        return target

    @property
    def params(self) -> tuple[Param, ...]:
        return tuple(self._params.values())

    @property
    def targets(self) -> tuple[Target, ...]:
        return tuple(self._targets.values())

    def __getitem__(self, name: str) -> Param:
        try:
            return self._params[name]
        except KeyError:
            raise KeyError(
                f"registry {self.piece!r} has no parameter {name!r}; it has "
                f"{', '.join(sorted(self._params)) or '(none)'}") from None

    def __contains__(self, name: str) -> bool:
        return name in self._params

    def __len__(self) -> int:
        return len(self._params)

    def free_scalars(self) -> tuple[Param, ...]:
        """The parameters no single measurement pins to a point. See :data:`FREE_TAGS`."""
        return tuple(p for p in self._params.values() if p.is_free)

    def independent_targets(self) -> tuple[Target, ...]:
        return tuple(t for t in self._targets.values() if not t.fitted)

    def by_tag(self) -> dict[str, int]:
        """A census in :class:`Tag` order, zeros included, so a table has fixed columns."""
        counts = dict.fromkeys(Tag.ALL, 0)
        for param in self._params.values():
            counts[param.tag] += 1
        return counts

    def uncited(self) -> tuple[Param, ...]:
        """Registered parameters whose source names no literature identifier.

        Not an error and not a gate. A number measured on this project's own plates cites a
        plate; a number measured in a paper and cited by prose alone is the one worth
        seeing, because rule 3 of this build -- a number attributed to a paper that does not
        contain it -- has burned this project repeatedly and cannot be checked at all
        without an identifier to resolve.
        """
        return tuple(p for p in self._params.values()
                     if p.tag in (Tag.MEASURED, Tag.BOUNDED, Tag.BORROWED, Tag.DERIVED)
                     and not p.identifiers)

    def borrowed(self) -> tuple[Param, ...]:
        """Every number imported from another organism, with the organism attached."""
        return tuple(p for p in self._params.values() if p.tag == Tag.BORROWED)

    def refusals(self) -> tuple[Param, ...]:
        return tuple(p for p in self._params.values() if p.tag == Tag.REFUSED)

    def gate(self) -> FreeScalarGate:
        """Criterion (e), computed rather than claimed."""
        return FreeScalarGate(
            piece=self.piece,
            free=tuple(p.name for p in self.free_scalars()),
            targets=tuple(t.name for t in self.independent_targets()),
            fitted_targets=tuple(t.name for t in self._targets.values() if t.fitted),
        )

    def require_gate(self) -> FreeScalarGate:
        """Refuse to build a piece with more free scalars than independent targets."""
        result = self.gate()
        if not result.passes:
            raise FreeScalarGateFailed(
                f"{result.summary()}. Free: {', '.join(result.free)}. Targets: "
                f"{', '.join(result.targets) or '(none)'}. Adding mechanism without adding "
                "refutability is what this gate refuses; either measure one of the free "
                "scalars, or find a target the piece has not already been fitted to")
        return result


def _md_cell(text: str) -> str:
    return str(text).replace("|", r"\|").replace("\n", " ")


def provenance_table(*registries: ParamRegistry, heading_level: int = 3) -> str:
    """A markdown provenance table per registry, for a docs page.

    Markdown rather than a rendered object because the docs in this repository are
    markdown that other tests read: `tests/test_architecture_doc_matches_the_code.py`
    already checks a document against the code, and a table it can parse is one a stale
    number cannot survive in.
    """
    if not registries:
        raise ValueError("provenance_table needs at least one registry")
    hashes = "#" * max(1, int(heading_level))
    out: list[str] = []
    for registry in registries:
        out.append(f"{hashes} {registry.piece}")
        out.append("")
        out.append(f"{registry.gate().summary()}")
        out.append("")
        out.append("| Parameter | Value | Units | Tag | Free | Provenance |")
        out.append("|---|---|---|---|---|---|")
        for param in registry.params:
            out.append(
                f"| `{_md_cell(param.name)}` | {_md_cell(param.display_value())} "
                f"| {_md_cell(param.units)} | **{param.tag}** "
                f"| {'yes' if param.is_free else 'no'} "
                f"| {_md_cell(param.provenance())} |")
        if registry.targets:
            out.append("")
            out.append("| Target | Observable | Assay | Noise floor | Counted | Source |")
            out.append("|---|---|---|---|---|---|")
            for target in registry.targets:
                out.append(
                    f"| `{_md_cell(target.name)}` | {_md_cell(target.observable)} "
                    f"| {_md_cell(target.assay)} | {target.noise_floor:g} "
                    f"{_md_cell(target.units)} "
                    f"| {'no (fitted)' if target.fitted else 'yes'} "
                    f"| {_md_cell(target.source)} |")
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def tag_census_table(*registries: ParamRegistry) -> str:
    """One row per registry, one column per tag, and the asserted share.

    The shape of the scoreboard the hardcoding audit produced for the existing package --
    65% ASSERTED against 14% MEASURED across 906 substantive literals. Whether the
    mechanistic layer does better than that is a number, and this prints it.
    """
    if not registries:
        raise ValueError("tag_census_table needs at least one registry")
    header = "| Piece | " + " | ".join(Tag.ALL) + " | Free | ASSERTED share |"
    rows = [header, "|" + "---|" * (len(Tag.ALL) + 3)]
    for registry in registries:
        counts = registry.by_tag()
        total = len(registry)
        share = "--" if not total else f"{counts[Tag.ASSERTED] / total:.0%}"
        rows.append(
            f"| {_md_cell(registry.piece)} | "
            + " | ".join(str(counts[tag]) for tag in Tag.ALL)
            + f" | {len(registry.free_scalars())} | {share} |")
    return "\n".join(rows) + "\n"


@dataclass(frozen=True)
class ParameterUncertainty:
    kind: str
    units: str
    explanation: str
    value: float | None = None
    interval: tuple[float, float] | None = None
    confidence: float | None = None

    def __post_init__(self) -> None:
        kinds = (
            "not_reported", "standard_deviation", "standard_error",
            "confidence_interval", "sensitivity_range", "unresolved", "exact",
        )
        if self.kind not in kinds:
            raise ValueError(f"unknown uncertainty kind {self.kind!r}")
        if not isinstance(self.units, str) or not self.units.strip():
            raise ValueError("uncertainty units are required")
        if not isinstance(self.explanation, str) or not self.explanation.strip():
            raise ValueError("uncertainty explanation is required")
        if self.kind in ("standard_deviation", "standard_error"):
            if isinstance(self.value, bool):
                raise ValueError("uncertainty value must be numeric, not boolean")
            number = _finite(self.value, "uncertainty", self.kind)
            if number < 0:
                raise ValueError("uncertainty value must be nonnegative")
            object.__setattr__(self, "value", number)
        elif self.value is not None:
            raise ValueError(f"{self.kind} cannot carry a standard uncertainty value")
        if self.kind in ("confidence_interval", "sensitivity_range"):
            if self.interval is None or len(self.interval) != 2:
                raise ValueError(f"{self.kind} requires a two-endpoint interval")
            if any(isinstance(value, bool) for value in self.interval):
                raise ValueError("uncertainty interval must be numeric, not boolean")
            low, high = (_finite(value, "interval", self.kind) for value in self.interval)
            if low > high:
                raise ValueError("uncertainty interval must be ordered")
            object.__setattr__(self, "interval", (low, high))
        elif self.interval is not None:
            raise ValueError(f"{self.kind} cannot carry an interval")
        if self.kind == "confidence_interval":
            if isinstance(self.confidence, bool):
                raise ValueError("confidence must be numeric, not boolean")
            confidence = _finite(self.confidence, "confidence", self.kind)
            if not 0 < confidence < 1:
                raise ValueError("confidence must be strictly between zero and one")
            object.__setattr__(self, "confidence", confidence)
        elif self.confidence is not None:
            raise ValueError("only a confidence interval has a confidence level; an SD is not a CI")

    def to_dict(self) -> dict:
        return {
            "kind": self.kind, "units": self.units, "explanation": self.explanation,
            "value": self.value,
            "interval": None if self.interval is None else list(self.interval),
            "confidence": self.confidence,
        }
