from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from enum import Enum
from types import MappingProxyType
from typing import Any

import numpy as np

from .params import Param, ParamRegistry, Tag


class ScientificRefusal(NotImplementedError):
    pass


class Compartment(str, Enum):
    MEDIUM = "medium"
    CELL_WATER = "cell_water"
    ER = "er"
    REGULATORY = "regulatory"
    DEAD = "dead"
    REACTOR = "reactor"
    LEDGER = "ledger"
    INSTRUMENT = "instrument"


@dataclass(frozen=True)
class Variable:
    units: str
    compartment: Compartment
    basis: str
    carbon: float = 0.0
    nitrogen: float = 0.0


EXTRACELLULAR = MappingProxyType({
    "glucose": (6.0, 0.0), "ethanol": (2.0, 0.0), "nitrogen": (0.0, 1.0),
    "oxygen": (0.0, 0.0), "glycerol": (3.0, 0.0), "peroxide": (0.0, 0.0),
    "acetate": (2.0, 0.0), "osmolyte": (0.0, 0.0),
})
INTRACELLULAR = MappingProxyType({
    "glucose": (6.0, 0.0), "triose": (3.0, 0.0), "pyruvate": (3.0, 0.0),
    "acetyl": (2.0, 0.0), "pentose": (5.0, 0.0), "glycogen": (6.0, 0.0),
    "glycerol": (3.0, 0.0), "ipp": (5.0, 0.0), "ggpp": (20.0, 0.0),
    "phytoene": (40.0, 0.0), "lycopene": (40.0, 0.0), "beta_carotene": (40.0, 0.0),
    "native_isoprenoid": (5.0, 0.0), "atp": (10.0, 5.0), "adp": (10.0, 5.0),
    "amp": (10.0, 5.0), "camp": (10.0, 5.0), "nad": (21.0, 7.0),
    "nadh": (21.0, 7.0), "nadp": (21.0, 7.0), "nadph": (21.0, 7.0),
    "ros": (0.0, 0.0), "acetate": (2.0, 0.0),
})
SIGNALS = ("pka", "snf1", "torc1", "yap1", "hsf1", "hog1", "upre", "acid")
DYNAMIC_SIGNALS = SIGNALS[:6]
UPR_VARIABLES = ("er_client", "hac1_unspliced", "hac1_spliced", "hac1_protein",
                 "kar2_mrna", "chaperone")
ROLES = ("crte", "crti", "crtyb", "gpd1", "hsp70", "antioxidant", "reporter", "inert")
LEDGERS = ("carbon_in", "carbon_out", "nitrogen_in", "nitrogen_out", "co2",
           "biomass_formed", "biomass_died", "biomass_washed_out", "oxygen_transferred",
           "atp_produced", "atp_consumed", "nadph_produced", "nadph_consumed",
           "expression_carbon_dead", "expression_nitrogen_dead", "protons_pumped", "acetate_exported")


def finite(value: float, name: str, *, minimum: float | None = None,
           maximum: float | None = None, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_, str, bytes)) or np.iscomplexobj(value):
        raise ValueError(f"{name} must be a real finite number, not {value!r}")
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a real finite number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    if positive and number <= 0:
        raise ValueError(f"{name} must be positive")
    if minimum is not None and number < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    if maximum is not None and number > maximum:
        raise ValueError(f"{name} must be at most {maximum}")
    return number


def frozen_mapping(values: Mapping) -> Mapping:
    return MappingProxyType(dict(values))


def amount_mapping(values: Mapping[str, float], allowed: Sequence[str] | Mapping,
                   name: str, *, fill: bool = False) -> Mapping[str, float]:
    if not isinstance(values, Mapping):
        raise ValueError(f"{name} must be a named mapping")
    extra = set(values) - set(allowed)
    if extra:
        raise ScientificRefusal(f"unsupported {name}: {sorted(extra)}; supported: {tuple(allowed)}")
    result = dict.fromkeys(allowed, 0.0) if fill else {}
    result.update({key: finite(value, f"{name}.{key}", minimum=0.0)
                   for key, value in values.items()})
    return frozen_mapping(result)


def parameter_value(param: Param, units: str, *, minimum: float = 0.0,
                    maximum: float | None = None, positive: bool = False) -> float:
    if not isinstance(param, Param):
        raise ValueError("biological and observation parameters require a provenance-bearing Param")
    if param.units != units:
        raise ValueError(f"{param.name}: expected units {units!r}, received {param.units!r}")
    return finite(float(param), param.name, minimum=minimum, maximum=maximum, positive=positive)


def prior(name: str, value: float, units: str, mechanism: str) -> Param:
    return Param.asserted(name, value, units,
                          f"PRIOR, not a published constant: shared-engine {mechanism}",
                          missing="strain-, medium-, and time-resolved independent calibration")


@dataclass(frozen=True)
class Control:
    time_h: float
    temperature_c: float = 30.0
    ph: float = 5.0
    feed_l_h: float = 0.0
    outflow_l_h: float = 0.0
    feed_mM: Mapping[str, float] = field(default_factory=dict)
    oxygen_transfer_per_h: float = 0.0
    oxygen_saturation_mM: float = 0.0
    folding_inhibition: float = 0.0

    def __post_init__(self) -> None:
        for name in ("time_h", "feed_l_h", "outflow_l_h", "oxygen_transfer_per_h",
                     "oxygen_saturation_mM"):
            object.__setattr__(self, name, finite(getattr(self, name), name, minimum=0.0))
        object.__setattr__(self, "temperature_c", finite(self.temperature_c, "temperature_c",
                                                        minimum=-273.15))
        object.__setattr__(self, "ph", finite(self.ph, "ph", positive=True, maximum=14.0))
        object.__setattr__(self, "folding_inhibition", finite(
            self.folding_inhibition, "folding_inhibition", minimum=0.0, maximum=1.0))
        object.__setattr__(self, "feed_mM", amount_mapping(self.feed_mM, EXTRACELLULAR, "feed_mM"))


@dataclass(frozen=True)
class Event:
    name: str
    time_h: float
    add_volume_l: float = 0.0
    add_mmol: Mapping[str, float] = field(default_factory=dict)
    withdraw_l: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("events require a stable nonempty name for restart-safe receipts")
        for name in ("time_h", "add_volume_l", "withdraw_l"):
            object.__setattr__(self, name, finite(getattr(self, name), name, minimum=0.0))
        object.__setattr__(self, "add_mmol", amount_mapping(self.add_mmol, EXTRACELLULAR, "add_mmol"))


@dataclass(frozen=True)
class Protocol:
    times_h: Sequence[float]
    controls: tuple[Control, ...]
    mode: str = "batch"
    events: tuple[Event, ...] = ()
    source: str = "declared in-silico protocol, not an experimental record"

    def __post_init__(self) -> None:
        times = tuple(finite(t, "times_h", minimum=0.0) for t in self.times_h)
        if len(times) < 2 or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("times_h requires at least two strictly increasing absolute times")
        if self.mode not in ("batch", "fedbatch", "chemostat"):
            raise ScientificRefusal(f"unsupported vessel mode {self.mode!r}")
        controls = tuple(self.controls)
        if not controls or any(not isinstance(c, Control) for c in controls):
            raise ValueError("at least one typed Control is required")
        if controls[0].time_h > times[0] or controls[-1].time_h > times[-1]:
            raise ValueError("controls must cover the start and not extend beyond the run")
        if any(b.time_h <= a.time_h for a, b in zip(controls, controls[1:])):
            raise ValueError("control times must be strictly increasing; no ambiguous ties")
        for c in controls:
            if self.mode == "batch" and (c.feed_l_h != 0 or c.outflow_l_h != 0):
                raise ValueError("batch mode cannot have continuous feed or outflow")
            if self.mode == "fedbatch" and c.outflow_l_h != 0:
                raise ValueError("fedbatch outflow is sampling through an Event, not continuous")
            if self.mode == "chemostat" and c.feed_l_h != c.outflow_l_h:
                raise ValueError("chemostat feed and outflow must be equal")
        events = tuple(self.events)
        if any(not isinstance(e, Event) for e in events):
            raise ValueError("events must be Event instances")
        if len({e.name for e in events}) != len(events):
            raise ValueError("event names must be unique")
        if any(not times[0] <= e.time_h <= times[-1] for e in events):
            raise ValueError("event times must lie inside the observation interval")
        if any(b.time_h < a.time_h for a, b in zip(events, events[1:])):
            raise ValueError("events must be time-ordered; ties execute in declared order")
        if not str(self.source).strip():
            raise ValueError("protocol source is required")
        object.__setattr__(self, "times_h", times)
        object.__setattr__(self, "controls", controls)
        object.__setattr__(self, "events", events)

    def control_at(self, time_h: float) -> Control:
        if time_h < self.controls[0].time_h:
            raise ValueError("no control at the requested time")
        return next(c for c in reversed(self.controls) if c.time_h <= time_h)

    @property
    def boundaries_h(self) -> tuple[float, ...]:
        start, stop = self.times_h[0], self.times_h[-1]
        return tuple(sorted({start, stop, *(c.time_h for c in self.controls if c.time_h >= start),
                             *(e.time_h for e in self.events)}))


PROMOTER_UNCOUPLED_SIGNALS = ("pka", "snf1")
"""The two dynamic signals that reach metabolism and deliberately reach no promoter.

They are not absent from the model: :mod:`mech.signalling` integrates both, and
:data:`mech.signalling.CARBON_NITROGEN_FLUX_COUPLING` names the rate laws that read them.
What is absent is a *transcriptional* coefficient, and that absence is a refusal with a
named missing measurement in :data:`PROMOTER_COUPLING_NOT_BUILT`, not an oversight.
:attr:`Genotype.promoter_drivers` is what makes the disconnection checkable on any
genotype, and :func:`sourced_promoter_response` is what would retire the refusal.
"""

QUANTIFIED_UNCERTAINTY = ("standard_deviation", "standard_error", "confidence_interval",
                          "sensitivity_range")
"""The `ParameterUncertainty` kinds that are an uncertainty. `not_reported` and
`unresolved` are the other two, and they are what every carbon and nitrogen row carries."""


def sourced_promoter_response(record: Mapping[str, Any]) -> bool:
    """Would one ``data/parameter_evidence.json`` row support a wired promoter coefficient?

    This repository's bar for a coefficient is a source with units and an uncertainty, so
    that is the whole predicate: declared units, and an uncertainty that is a number rather
    than a note saying none was reported. It is deliberately not a promoter-specific test --
    no carbon or nitrogen row clears even this weaker bar, so nothing sharper is needed to
    justify the refusal, and the day a row does clear it the refusal has to be revisited.
    """
    if not isinstance(record, Mapping):
        raise ValueError("an evidence row is a mapping of the schema-1 parameter fields")
    uncertainty = record.get("uncertainty")
    return (record.get("units_status") == "source_declared"
            and isinstance(uncertainty, Mapping)
            and uncertainty.get("kind") in QUANTIFIED_UNCERTAINTY)


PROMOTER_COUPLING_NOT_BUILT = ParamRegistry("mech/contracts.py::promoter-coupling-not-built")
"""The transcriptional coefficients work package 7 asked for and no source supports.

On `mech/oxidative.py::NOT_BUILT`'s precedent: a refusal that nothing reads is a record,
not a degree of freedom. The record is the point here -- the carbon programme's states are
*already* coupled, to catabolite-repression fluxes, so the only thing missing is the edge
into a promoter, and inventing that edge is what this registry exists to prevent."""

_JALIHAL = ("Jalihal, Kraikivski, Murali & Tyson 2021, Mol Biol Cell 32:ar24, PMID 34495680, "
            "doi 10.1091/mbc.E20-02-0117; local SBML sha256 a5416c7785df8001044602f2c47cfa08"
            "4109396a0ee165b36e844d35b00fb1ed, audited in data/parameter_evidence.json")
_ODUIBHIR = ("O'Duibhir et al. 2014, Mol Syst Biol 10:732, doi 10.15252/msb.20145172, cited in "
             "docs/research/HYBRID_METHODS.md")
_NO_PROMOTER_NODE = (
    "the source stops one edge short of the coefficient wanted here. Of its 25 species, PKA, "
    "Snf1, Mig1, Dot6, Gis1 and Gln3 are normalized regulator activities, not promoter "
    "occupancies. Exactly one node is transcription-shaped -- Rib' = k_transcription * "
    "(1 - Dot6) - k_mRNA_degr * Rib, at 0.235950807701 and 0.0730004044196 -- and it is the "
    "wrong coefficient three times over: its target is the ribosomal-biogenesis regulon and "
    "not any promoter this project ships; PKA enters Dot6 only as the product Sch9*PKA, so "
    "no PKA-specific weight is separable from it; and both constants carry a rate unit that "
    "the inventory records as unresolved between the source's minutes and the exporter's "
    "seconds, which is a 60x ambiguity in their magnitude. Neither they nor the model's one "
    "regulator-to-target edge (w_gln1_gln3 = 0.519802074308, Gln3 to glutamine synthetase, "
    "and on the nitrogen arm) appears in data/parameter_evidence.json at all, so neither has "
    "been through this repository's source audit. The inventory separately blocks the two "
    "conversions any such import would need: glucose to the normalized Carbon/ATP input, and "
    "normalized state to absolute kinase activity. The expression-compendium route is closed "
    "too -- " + _ODUIBHIR + " shows a large share of an apparently stress-specific signature "
    "is a growth-rate artefact, so a coefficient regressed on batch arrays would encode "
    "'slow', not 'this kinase'")

PROMOTER_COUPLING_NOT_BUILT.add(Param.refused(
    "promoter_response.pka", "dimensionless promoter occupancy per unit PKA activity",
    _JALIHAL + ", row carbon.pka.w_pka_camp = 102.109826697",
    reason="the only carbon-side PKA coefficient this repository owns is a weight inside "
           "Jalihal's normalized network, and it fails both halves of the standard: its "
           "units_status is unspecified_in_source and its uncertainty is not_reported "
           "against a fit the source itself describes as an ensemble of acceptable "
           "parameter sets, so a single point is not identified. " + _NO_PROMOTER_NODE,
    missing="a cAMP or PKA dose series on a promoter this project actually reads, at fixed "
            "growth rate, with PKA activity assayed in the same cells -- occupancy per unit "
            "activity, with an interval. A chemostat or a retentostat removes the growth "
            "confound that a batch series cannot"))

PROMOTER_COUPLING_NOT_BUILT.add(Param.refused(
    "promoter_response.snf1", "dimensionless promoter occupancy per unit Snf1 activity",
    _JALIHAL + ", row carbon.snf1.w_mig_snf = 1.21479427525; Garcia-Salcedo et al. 2014, "
    "FEBS J, PMID 24529170, doi 10.1111/febs.12753, inventory_only with no parameter table "
    "checked, and the Persson 2022 review (PMC8916112) counts 24 alternative fitted small "
    "models of the same Snf1/Mig1 step",
    reason="w_mig_snf is a Snf1-to-Mig1 weight, not a Mig1-to-promoter occupancy, and it "
           "carries unspecified_in_source units with not_reported uncertainty. Twenty-four "
           "alternative published fits of the step it describes is the definition of an "
           "unidentified coefficient. Mig1 is also a terminal node of the source network: "
           "no equation in the model reads it, so the Snf1 arm dead-ends there in the "
           "source exactly as it does here, and importing it would import the dead end. "
           + _NO_PROMOTER_NODE,
    missing="a glucose ladder at fixed growth rate on a Snf1-dependent promoter in this "
            "host, with Snf1 phosphorylation and Mig1 localization measured in the same "
            "cells -- occupancy per unit activity, with an interval"))


MEDIUM_SPECIES_NOT_CARRIED = ParamRegistry("mech/contracts.py::medium-species-not-carried")
"""Medium components a transcribed stress axis doses and :data:`EXTRACELLULAR` does not carry.

Adding a species to EXTRACELLULAR is one line of Python and a whole experiment: every row
carries a carbon and a nitrogen count, joins the conservation audit, and needs a transport
law before any flux may read it. So an axis whose dose has no species is registered UNDRIVEN
against a row here, rather than given a medium this repository invented for it.
"""

MEDIUM_SPECIES_NOT_CARRIED.add(Param.refused(
    "medium.calcium", "mM",
    "generator/stress_panel.py STRESSORS['calcium_chloride'], lethal dose 200 mM, calcium "
    "weight 1.00; refused from the source side in mech/calcium_ke2013.py "
    "NOT_DRIVEN['extracellular_calcium_mM']",
    reason="EXTRACELLULAR is exactly eight species and none is calcium. The gap is not only "
           "this contract's: Ke 2013 Text S1 Eqn 3.3, the only calcium ODE transcribed in this "
           "repository, has no extracellular-calcium term at all -- its cytosolic pool takes a "
           "constant basal influx plus sodium and external-pH arms -- so a medium calcium "
           "species added here would have nothing to drive. Both halves are missing at once",
    missing="a published parameterised S. cerevisiae calcium model whose influx reads medium "
            "calcium, AND the transport law plus elemental counts a new EXTRACELLULAR row "
            "needs. Cui & Kaandorp 2006 is the only candidate for the first and does not "
            "supply the second"))

MEDIUM_SPECIES_NOT_CARRIED.add(Param.refused(
    "medium.sodium", "mM",
    "generator/stress_panel.py STRESSORS['NaCl'], calcium weight 0.25; refused from the source "
    "side in mech/calcium_ke2013.py NOT_DRIVEN['external_sodium_mM'] and, one compartment in, "
    "NOT_DRIVEN['cytosolic_sodium_mM']",
    reason="EXTRACELLULAR has no sodium, so a NaCl dose cannot be keyed into Control.feed_mM, "
           "and INTRACELLULAR has no sodium either. 'osmolyte' is deliberately ion-agnostic: "
           "reading it as sodium would assert an ionic fraction of the osmolyte pool that "
           "nothing here measures, which is precisely the kind of proxy this row exists to "
           "refuse. Ke's Eqn 3.3 sodium arms therefore run at his own unstressed baselines "
           "(cytosolic 100 mM, external 5 mM) as CONSTANTS, not as signals",
    missing="a medium sodium species and an intracellular sodium state in the shared vector, "
            "with a transport law joining them. ph_ke2013 carries the transport laws already; "
            "the two species are what is absent"))

MEDIUM_SPECIES_NOT_CARRIED.add(Param.refused(
    "medium.respiratory_inhibitor", "mM",
    "generator/stress_panel.py STRESSORS['antimycin_A'], retrograde weight +0.90 and atp "
    "weight -0.85; refused from the source side in mech/atp_teusink2000.py PANEL_LIMIT",
    reason="EXTRACELLULAR carries no xenobiotic of any kind, and the two sources that would "
           "read one cannot: Jalihal 2021's deposit contains no Rtg2, Mks1, Bmh1 or Bmh2 and "
           "reaches its rho0 condition by lowering a glutamine input by hand rather than by a "
           "mechanism, and Teusink 2000 has no oxygen species and no respiratory chain to "
           "inhibit. A medium species alone would drive nothing",
    missing="a parameterised dynamic model of mitochondrial-dysfunction signalling to Rtg1/3, "
            "plus a respiratory arm for the inhibitor to act on. The panel's single heaviest "
            "retrograde stressor stays forfeited until both exist"))

MEDIUM_SPECIES_NOT_CARRIED.add(Param.refused(
    "medium.cell_wall_damaging_agent", "ug/mL congo red or mM caffeine",
    "generator/stress_panel.py STRESSORS['congo_red'] (50 ug/mL) and STRESSORS['caffeine'] "
    "(10 mM), both routed to the cell_wall module; refused from the source side in "
    "mech/cell_wall_talemi2016.py NOT_DRIVEN['cell_wall_damage_dose'] and, on this same "
    "contract, NOT_DRIVEN['medium.cell_wall_damaging_agent']",
    reason="EXTRACELLULAR is exactly eight species and none is a wall-damaging agent, and the "
           "gap is worse than a missing row: congo red and calcofluor white are not "
           "metabolised and not transported, they bind chitin and glucan AT THE WALL, and this "
           "contract has no cell-wall compartment for them to act in. The source side is "
           "missing too. Talemi 2016, the only transcribed CWI model here, has no damage input "
           "at all -- its single stimulus is external osmolarity and its single Slt2 driver is "
           "cell volume, through v9 = k9*Vos, which activates Slt2 when the cell SWELLS. So a "
           "damaging agent wired to it would also invert its sign",
    missing="a cell-wall compartment or a wall-binding term, the elemental counts and transport "
            "law a new EXTRACELLULAR row needs, and a published parameterised model whose Slt2 "
            "activation reads the dose. Talemi's own calcofluor module is the nearest candidate "
            "and is refused separately at its source, being absent from the deposit"))


AXIS_PROMOTER_NOT_BUILT = ParamRegistry("mech/contracts.py::axis-promoter-not-built")
"""Promoter coefficients the newly transcribed axes would need, and no source supplies.

Separate from :data:`PROMOTER_COUPLING_NOT_BUILT` on purpose. That registry is the carbon
programme's refusal and is pinned one-to-one against :data:`PROMOTER_UNCOUPLED_SIGNALS`,
which is a statement about the six integrated signals. These rows are about axes that are
registered and NOT integrated, so folding them into that registry would misreport their
status as 'a signal that reaches metabolism but no promoter' when the truth is weaker: they
reach nothing yet, and the missing edge is named here so the difference stays visible.
"""

AXIS_PROMOTER_NOT_BUILT.add(Param.refused(
    "promoter_response.crz1_cdre", "dimensionless promoter occupancy per mM Crz1p",
    "Ke 2013, PLoS Comput Biol 9(1):e1002879, Text S1 Eqn 3.6 and the Table S3 ENA1 block; "
    "the panel declares this module against a CDRE (Stathopoulos & Cyert 1997, PMID 9407035)",
    reason="Crz1p is where Ke's calcium block ends as an executable object. Its only printed "
           "output arm is ENA1, and that arm is untranscribable twice over: it sits downstream "
           "of the Nrg1p ODE, whose KmNrg1,pH appears in Eqn 3.5 and in no table, and Table S3 "
           "prints KmENA1,Nrg1 as '1.143.0*1e-4', which is not a number. No CDRE gain exists "
           "anywhere in this source, so the panel's declared reporter has no coefficient",
    missing="a published parameterised CDRE (PMC1/CMK2) transcription model, or a Crz1p dose "
            "response on a CDRE reporter in the working strain, with an interval"))

AXIS_PROMOTER_NOT_BUILT.add(Param.refused(
    "promoter_response.rlm1_box", "dimensionless promoter occupancy per uM Slt2PP",
    "generator/stress_panel.py declares the cell_wall module as 'Slt2/Mpk1 -> Rlm1' against an "
    "'RLM1 box', citing Jung & Levin 1999 PMID 10594829; refused from the source side in "
    "mech/cell_wall_talemi2016.py NOT_DRIVEN['promoter_response.rlm1_box']",
    reason="Talemi 2016 MODEL1606100000 ends at phosphorylated Slt2 as a KINASE and has no "
           "output arm of any kind: the strings Rlm1, RLM1 and 'transcription' occur zero times "
           "in the pinned SBML, in the article body and in the supplementary text. Jung & Levin "
           "1999 establishes that Rlm1 carries most CWI transcription and reports no "
           "coefficient; it is not a parameterised model. The panel's declared reporter has no "
           "gain anywhere in this axis's evidence",
    missing="a published parameterised RLM1-box transcription model, or an Slt2 dose response "
            "on an RLM1-box reporter in the working strain, with an interval -- and the sign "
            "settled with it, since the only Slt2 activation this source models is hypo-osmotic "
            "and the panel's stressors are wall-damaging"))


@dataclass(frozen=True)
class AxisRegistration:
    """One transcribed stress axis, its states, its environment input, and why it is undriven.

    The point of the type is that "registered but undriven" has to cost something to declare.
    An axis with no reader must name the reason in :attr:`undriven_because`, every medium
    species it says is missing must be a registered refusal AND genuinely absent from
    :data:`EXTRACELLULAR`, and an axis the panel declares against a promoter must cite a
    refused coefficient. None of those can be satisfied by silence.
    """

    axis: str
    module: str
    source: str
    states: tuple[tuple[str, str, Compartment, str], ...]
    environment_input: str
    missing_medium: tuple[str, ...]
    panel_promoter: str
    promoter_refusals: tuple[str, ...]
    readers: tuple[str, ...]
    undriven_because: str

    def __post_init__(self) -> None:
        for name in ("axis", "module", "source", "environment_input", "panel_promoter"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"axis registration requires a nonempty {name}")
        states = tuple(self.states)
        if not states:
            raise ValueError(f"axis {self.axis!r} registers no state variable")
        if len({row[0] for row in states}) != len(states):
            raise ValueError(f"axis {self.axis!r} registers a state name twice")
        for row in states:
            if len(row) != 4 or not isinstance(row[2], Compartment):
                raise ValueError(f"axis {self.axis!r} states are (name, units, Compartment, meaning)")
            if not all(isinstance(field_, str) and field_.strip() for field_ in (row[0], row[1], row[3])):
                raise ValueError(f"axis {self.axis!r} states need a name, units and a meaning")
        for qualified in self.missing_medium:
            row = refused_row(qualified)
            species = row.name.split(".", 1)[1]
            if species in EXTRACELLULAR:
                raise ScientificRefusal(
                    f"axis {self.axis!r} claims EXTRACELLULAR lacks {species!r} and it does not")
        if self.panel_promoter != "none" and not self.promoter_refusals:
            raise ScientificRefusal(
                f"axis {self.axis!r} is declared against the {self.panel_promoter} promoter with "
                "no refused coefficient; a promoter an axis cannot reach is a refusal with a "
                "named missing measurement, not an omission")
        for qualified in self.promoter_refusals:
            refused_row(qualified)
        if not self.readers and len(self.undriven_because.split()) < 12:
            raise ScientificRefusal(
                f"axis {self.axis!r} is registered with no flux or promoter that reads it and no "
                "stated reason; an axis reported as undriven must say what is missing, or it is "
                "being reported as dead")
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "missing_medium", tuple(self.missing_medium))
        object.__setattr__(self, "promoter_refusals", tuple(self.promoter_refusals))
        object.__setattr__(self, "readers", tuple(self.readers))

    @property
    def driven(self) -> bool:
        return bool(self.readers)

    @property
    def state_names(self) -> tuple[str, ...]:
        """The registered names, namespaced the way the integrated states are."""
        return tuple(f"axis.{self.axis}.{row[0]}" for row in self.states)

    def variables(self) -> Mapping[str, Variable]:
        """The registration in the engine's own vocabulary, with a basis that cannot integrate.

        ``registered, not integrated`` is deliberately not one of the bases `_Kernel.physical`
        accepts, so one of these rows joining the state vector by accident is a hard error.
        """
        return frozen_mapping({
            f"axis.{self.axis}.{name}": Variable(units, compartment, "registered, not integrated")
            for name, units, compartment, _ in self.states})


_AXIS_REFUSAL_REGISTRIES = MappingProxyType({
    "PROMOTER_COUPLING_NOT_BUILT": PROMOTER_COUPLING_NOT_BUILT,
    "AXIS_PROMOTER_NOT_BUILT": AXIS_PROMOTER_NOT_BUILT,
    "MEDIUM_SPECIES_NOT_CARRIED": MEDIUM_SPECIES_NOT_CARRIED,
})


def refused_row(qualified: str) -> Param:
    """Resolve ``REGISTRY::name`` to a REFUSED parameter, or raise.

    An axis cites its refusals by string so the citation is checkable at import: a row that
    was renamed, or promoted out of REFUSED, fails here instead of leaving an axis quietly
    claiming a refusal that no longer says what it used to.
    """
    registry_name, _, name = str(qualified).partition("::")
    if registry_name not in _AXIS_REFUSAL_REGISTRIES or not name:
        raise ValueError(f"{qualified!r} must be REGISTRY::name, one of {sorted(_AXIS_REFUSAL_REGISTRIES)}")
    param = _AXIS_REFUSAL_REGISTRIES[registry_name][name]
    if param.tag != Tag.REFUSED:
        raise ScientificRefusal(f"{qualified} is cited as a refusal and is tagged {param.tag}")
    return param


@dataclass(frozen=True)
class Promoter:
    basal: Param
    responses: Mapping[str, Param] = field(default_factory=dict)

    def __post_init__(self) -> None:
        parameter_value(self.basal, "dimensionless", maximum=1.0)
        extra = set(self.responses) - set(SIGNALS)
        if extra:
            raise ScientificRefusal(f"unsupported transcription factor inputs: {sorted(extra)}")
        for param in self.responses.values():
            parameter_value(param, "dimensionless", minimum=-1.0, maximum=1.0)
        object.__setattr__(self, "responses", frozen_mapping(self.responses))

    def activity(self, signals: Mapping[str, float]) -> float:
        return float(np.clip(float(self.basal) + sum(float(p) * signals[name]
                                                    for name, p in self.responses.items()), 0.0, 1.0))


@dataclass(frozen=True)
class Gene:
    name: str
    role: str
    copies: int
    protein_aa: int
    transcript_nt: int
    promoter: Promoter
    transcription_mmol_gdw_h: Param
    translation_per_h: Param
    mrna_decay_per_h: Param
    protein_decay_per_h: Param
    source: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name or "." in self.name:
            raise ValueError("gene names must be nonempty and cannot contain '.'")
        if self.role not in ROLES:
            raise ScientificRefusal(f"unsupported gene role {self.role!r}")
        for name in ("copies", "protein_aa", "transcript_nt"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{name} must be an integer")
            if value < (0 if name == "copies" else 1):
                raise ValueError(f"invalid {name}")
        if self.transcript_nt < 3 * self.protein_aa:
            raise ValueError("transcript_nt cannot be shorter than the coding sequence")
        if not isinstance(self.promoter, Promoter) or not str(self.source).strip():
            raise ValueError("a typed promoter and gene provenance are required")
        parameter_value(self.transcription_mmol_gdw_h, "mmol/gDW/h")
        for param in (self.translation_per_h, self.mrna_decay_per_h, self.protein_decay_per_h):
            parameter_value(param, "1/h")

    @property
    def parameters(self) -> tuple[Param, ...]:
        return (self.promoter.basal, *self.promoter.responses.values(),
                self.transcription_mmol_gdw_h, self.translation_per_h,
                self.mrna_decay_per_h, self.protein_decay_per_h)


@dataclass(frozen=True)
class Genotype:
    genes: tuple[Gene, ...]
    signal_activity: Mapping[str, float] = field(default_factory=dict)
    source: str = "declared genotype; no copy-number inference"
    inheritance: str = "stable"

    def __post_init__(self) -> None:
        genes = tuple(self.genes)
        if any(not isinstance(g, Gene) for g in genes):
            raise ValueError("genes must be typed Gene instances")
        if len({g.name for g in genes}) != len(genes):
            raise ValueError("gene names must be unique")
        if self.inheritance != "stable":
            raise ScientificRefusal("episomal segregation/copy-number control requires a multi-population calibration; only stable cassettes are supported")
        activity = amount_mapping(self.signal_activity, SIGNALS, "signal_activity")
        if any(x > 1 for x in activity.values()):
            raise ValueError("signal_activity is a residual activity fraction in [0, 1]")
        if not str(self.source).strip():
            raise ValueError("genotype source is required")
        object.__setattr__(self, "genes", genes)
        object.__setattr__(self, "signal_activity", activity)

    @property
    def promoter_drivers(self) -> tuple[str, ...]:
        """Which signals actually reach a promoter here, in :data:`SIGNALS` order.

        The complement is the checkable half of :data:`PROMOTER_UNCOUPLED_SIGNALS`: read it
        off the genotype rather than grepping for driver names, because a cassette added
        later moves this and a grep would not notice.
        """
        driven = {name for gene in self.genes for name in gene.promoter.responses}
        return tuple(name for name in SIGNALS if name in driven)

    def with_copies(self, **copies: int) -> Genotype:
        extra = set(copies) - {g.name for g in self.genes}
        if extra:
            raise ValueError(f"unknown genes: {sorted(extra)}")
        return replace(self, genes=tuple(replace(g, copies=copies.get(g.name, g.copies))
                                         for g in self.genes))

    @classmethod
    def prior(cls, *, product: bool = True, reporter: bool = True) -> Genotype:
        definitions = [("hsp70", "hsf1"), ("antioxidant", "yap1"), ("gpd1", "hog1")]
        if product:
            definitions.extend((name, "torc1") for name in ("crte", "crti", "crtyb"))
        if reporter:
            definitions.append(("reporter", "upre"))
        genes = []
        for name, driver in definitions:
            promoter = Promoter(prior(f"{name}.basal", 0.15, "dimensionless", "basal promoter occupancy"),
                                {driver: prior(f"{name}.{driver}", 0.85, "dimensionless",
                                               "linear clipped TF-to-promoter response")})
            genes.append(Gene(
                name, name, 1, 300, 900, promoter,
                prior(f"{name}.transcription", 1e-9, "mmol/gDW/h", "transcription capacity"),
                prior(f"{name}.translation", 100.0, "1/h", "translation capacity per transcript"),
                prior(f"{name}.mrna_decay", 2.0, "1/h", "mRNA decay"),
                prior(f"{name}.protein_decay", 0.03, "1/h", "protein turnover"),
                "PRIOR synthetic cassette, not a sequence annotation or a measured strain",
            ))
        return cls(tuple(genes), source="PRIOR illustrative stable-cassette genotype")


@dataclass(frozen=True)
class ExpressionState:
    mrna_mmol: float = 0.0
    protein_mmol: float = 0.0
    active_mmol: float = 0.0

    def __post_init__(self) -> None:
        for name in ("mrna_mmol", "protein_mmol", "active_mmol"):
            object.__setattr__(self, name, finite(getattr(self, name), name, minimum=0.0))
        if self.active_mmol > self.protein_mmol:
            raise ValueError("active protein is a subset of total protein")


@dataclass(frozen=True)
class PhysicalState:
    time_h: float
    volume_l: float
    biomass_gdw: float
    extracellular_mmol: Mapping[str, float]
    intracellular_mmol: Mapping[str, float] = field(default_factory=dict)
    dead_biomass_gdw: float = 0.0
    dead_intracellular_mmol: Mapping[str, float] = field(default_factory=dict)
    expression: Mapping[str, ExpressionState] = field(default_factory=dict)
    signals: Mapping[str, float] = field(default_factory=dict)
    upr: Mapping[str, float] = field(default_factory=dict)
    cell_volume_ratio: float = 1.0
    unfolded_fraction: float = 0.0
    ledgers: Mapping[str, float] = field(default_factory=dict)
    event_receipts: tuple[Event, ...] = ()

    def __post_init__(self) -> None:
        for name in ("time_h", "biomass_gdw", "dead_biomass_gdw"):
            object.__setattr__(self, name, finite(getattr(self, name), name, minimum=0.0))
        for name in ("volume_l", "cell_volume_ratio"):
            object.__setattr__(self, name, finite(getattr(self, name), name, positive=True))
        object.__setattr__(self, "unfolded_fraction", finite(
            self.unfolded_fraction, "unfolded_fraction", minimum=0.0, maximum=1.0))
        for name, allowed in (("extracellular_mmol", EXTRACELLULAR),
                              ("intracellular_mmol", INTRACELLULAR),
                              ("dead_intracellular_mmol", INTRACELLULAR),
                              ("signals", DYNAMIC_SIGNALS), ("upr", UPR_VARIABLES)):
            object.__setattr__(self, name, amount_mapping(getattr(self, name), allowed, name, fill=True))
        if any(x > 1 for x in self.signals.values()):
            raise ValueError("signal states are fractions in [0, 1]")
        if any(not isinstance(x, ExpressionState) for x in self.expression.values()):
            raise ValueError("expression values must be ExpressionState instances")
        if set(self.ledgers) - set(LEDGERS):
            raise ValueError("unknown conservation ledger")
        ledger = dict.fromkeys(LEDGERS, 0.0)
        ledger.update({k: finite(v, k) for k, v in self.ledgers.items()})
        object.__setattr__(self, "ledgers", frozen_mapping(ledger))
        object.__setattr__(self, "expression", frozen_mapping(self.expression))
        receipts = tuple(self.event_receipts)
        if any(not isinstance(e, Event) or e.time_h > self.time_h for e in receipts):
            raise ValueError("event receipts must be events at or before the state's clock")
        if len({e.name for e in receipts}) != len(receipts):
            raise ValueError("duplicate event receipts")
        object.__setattr__(self, "event_receipts", receipts)
        if self.biomass_gdw == 0 and (any(self.intracellular_mmol.values()) or
                                     any(e.protein_mmol or e.mrna_mmol for e in self.expression.values())):
            raise ValueError("live intracellular inventories require positive live biomass")

    @classmethod
    def from_concentrations(cls, *, volume_l: float, biomass_gdw_l: float,
                            medium_mM: Mapping[str, float], time_h: float = 0.0) -> PhysicalState:
        volume = finite(volume_l, "volume_l", positive=True)
        density = finite(biomass_gdw_l, "biomass_gdw_l", minimum=0.0)
        medium = amount_mapping(medium_mM, EXTRACELLULAR, "medium_mM")
        return cls(time_h, volume, density * volume, {k: v * volume for k, v in medium.items()})


OBSERVATION_UNITS = MappingProxyType({
    "gDW_per_OD_l": "gDW/L/OD", "od_saturation": "1/OD", "od_blank": "OD",
    "rfu_gain": "RFU L/mmol", "rfu_blank": "RFU", "autofluorescence": "RFU L/gDW",
    "inner_filter": "gDW/mmol", "maturation": "1/h", "maturation_oxygen_k": "mM",
    "od_noise_sd": "OD", "rfu_noise_sd": "RFU", "missing_probability": "dimensionless",
})


@dataclass(frozen=True)
class ObservationModel:
    reporter: str
    values: Mapping[str, Param]
    seed: int = 0
    initial_mature_mmol: float = 0.0

    def __post_init__(self) -> None:
        missing = set(OBSERVATION_UNITS) - set(self.values)
        if missing:
            raise ScientificRefusal(f"missing optical/maturation calibration: {sorted(missing)}")
        if set(self.values) - set(OBSERVATION_UNITS):
            raise ValueError("unknown observation parameters")
        for name, units in OBSERVATION_UNITS.items():
            parameter_value(self.values[name], units,
                            positive=name in ("gDW_per_OD_l", "rfu_gain", "maturation_oxygen_k"),
                            maximum=1.0 if name == "missing_probability" else None)
        if isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        finite(self.initial_mature_mmol, "initial_mature_mmol", minimum=0.0)
        object.__setattr__(self, "values", frozen_mapping(self.values))

    @classmethod
    def prior(cls, reporter: str = "reporter", *, seed: int = 0) -> ObservationModel:
        numbers = (0.5, 0.1, 0.0, 1e9, 0.0, 5.0, 10.0, 2.0, 0.01, 0.0, 0.0, 0.0)
        return cls(reporter, {name: prior(name, value, units, "uncalibrated illustrative optics")
                              for (name, units), value in zip(OBSERVATION_UNITS.items(), numbers, strict=True)},
                   seed=seed)


@dataclass(frozen=True)
class AcidExportParameters:
    vmax: Param
    km_anion: Param
    atp_per_proton: Param
    atp_per_anion: Param
    source: str

    def __post_init__(self) -> None:
        parameter_value(self.vmax, "mmol/gDW/h")
        parameter_value(self.km_anion, "mM", positive=True)
        parameter_value(self.atp_per_proton, "mmol ATP/mmol H+", positive=True)
        parameter_value(self.atp_per_anion, "mmol ATP/mmol anion")
        if not str(self.source).strip():
            raise ValueError("paired acid-export/proton-pumping provenance is required")

    @classmethod
    def prior(cls) -> AcidExportParameters:
        return cls(
            prior("acid_export_vmax", 0.5, "mmol/gDW/h", "paired anion/proton export capacity; no measured Vmax in mech/ph.py::NOT_BUILT"),
            prior("acid_export_km", 10.0, "mM", "anion export affinity; no measured Km in mech/ph.py::NOT_BUILT"),
            prior("atp_per_proton", 1.0, "mmol ATP/mmol H+", "mech/ph.py::proton_bill glucose-fed 1 H+/ATP convention; coupling outside that context is not measured"),
            prior("atp_per_anion", 1.0, "mmol ATP/mmol anion", "unidentified anion-export ATP cost; not assigned to Tpo2, Tpo3 or Pdr12"),
            "PRIOR paired export of one acetate anion and one proton, locally electroneutral; "
            "not an independently regulated Pma1 or identified transporter model. Missing capacity, "
            "affinity, coupling and regulation are recorded in data/parameter_evidence.json, programmes.ph.",
        )


@dataclass(frozen=True)
class ThermodynamicRequest:
    data: Any
    source: str
    require_complete: bool = False
    mode: str = "audit"
    reactions: tuple[str, ...] = ("lcy",)
    metabolite_ids: Mapping[str, str] = field(default_factory=lambda: {
        "lycopene": "lycopene_c", "beta_carotene": "betacarotene_c",
    })
    uncertainty_sigma: float | None = None
    activity_coefficients: Mapping[str, Param] = field(default_factory=dict)

    def __post_init__(self) -> None:
        from ..bridge.thermodynamic import ThermodynamicData

        if not isinstance(self.data, ThermodynamicData) or not str(self.source).strip():
            raise ValueError("thermodynamic data and its source are required")
        if not isinstance(self.require_complete, bool):
            raise ValueError("require_complete must be a bool")
        if self.mode not in ("audit", "constrain_supported"):
            raise ValueError("thermodynamic mode must be audit or constrain_supported")
        reactions = tuple(self.reactions)
        if not reactions or len(set(reactions)) != len(reactions):
            raise ValueError("thermodynamic reactions must be nonempty and unique")
        required = {"lcy": {"lycopene", "beta_carotene"},
                    "adenylate_kinase": {"atp", "amp", "adp"}}
        unsupported = set(reactions) - set(required)
        if unsupported:
            raise ScientificRefusal(f"unsupported thermodynamic constraints {sorted(unsupported)}: complete molecular stoichiometry, compartment-matched activities and reliable transformed energies are required")
        needed = set().union(*(required[name] for name in reactions))
        missing = needed - set(self.metabolite_ids)
        if missing:
            raise ScientificRefusal(f"explicit thermodynamic metabolite identity mapping required for {sorted(missing)}")
        identifiers = tuple(self.metabolite_ids.values())
        if (set(self.metabolite_ids) - set().union(*required.values())
                or any(not isinstance(x, str) or not x.strip() for x in identifiers)
                or len(set(identifiers)) != len(identifiers)):
            raise ValueError("thermodynamic metabolite mappings must be known, nonempty and one-to-one")
        if self.uncertainty_sigma is not None:
            finite(self.uncertainty_sigma, "uncertainty_sigma", minimum=0.0)
        if self.mode == "constrain_supported" and needed - set(self.activity_coefficients):
            raise ScientificRefusal("thermodynamic constraints require explicit activity_coefficients: free/accessible activity is not identified by total intracellular content, especially for carotenoids")
        if set(self.activity_coefficients) - needed:
            raise ValueError("activity coefficients must belong to the requested molecular reactions")
        for param in self.activity_coefficients.values():
            parameter_value(param, "dimensionless", positive=True)
        object.__setattr__(self, "reactions", reactions)
        object.__setattr__(self, "metabolite_ids", frozen_mapping(self.metabolite_ids))
        object.__setattr__(self, "activity_coefficients", frozen_mapping(self.activity_coefficients))


@dataclass(frozen=True)
class PoolPartition:
    integrated: tuple[str, ...]
    fast_algebraic: Mapping[str, str]
    retained_fast_pool_candidates: tuple[str, ...]
    observation_only_dynamic: tuple[str, ...]

    def __post_init__(self) -> None:
        if set(self.integrated).intersection(self.fast_algebraic):
            raise ValueError("an independent pool cannot have both dynamic and algebraic owners")
        if not set(self.retained_fast_pool_candidates) <= set(self.integrated):
            raise ValueError("fast pool candidates must remain in the integrated state")
        object.__setattr__(self, "fast_algebraic", frozen_mapping(self.fast_algebraic))


@dataclass(frozen=True)
class Validity:
    parameter_basis: str
    assumptions: tuple[str, ...]
    unsupported: tuple[str, ...]
    thermodynamic_covered: tuple[str, ...] = ()
    thermodynamic_uncovered: tuple[str, ...] = ()
    biological_validation: bool = False


@dataclass(frozen=True)
class Observations:
    values: Mapping[str, np.ndarray]
    missing: Mapping[str, np.ndarray]
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        for name in ("values", "missing"):
            out = {}
            for key, value in getattr(self, name).items():
                arr = np.array(value, copy=True)
                arr.setflags(write=False)
                out[key] = arr
            object.__setattr__(self, name, frozen_mapping(out))
        object.__setattr__(self, "metadata", frozen_mapping(self.metadata))


@dataclass(frozen=True)
class SimulationResult:
    times_h: np.ndarray
    truth: Mapping[str, np.ndarray]
    variables: Mapping[str, Variable]
    final_state: PhysicalState
    validity: Validity
    provenance: Mapping[str, Param]
    diagnostics: Mapping[str, Any]
    observations: Observations | None = None

    def __post_init__(self) -> None:
        times = np.array(self.times_h, dtype=float, copy=True)
        times.setflags(write=False)
        out = {}
        if set(self.truth) != set(self.variables):
            raise ValueError("every truth variable needs a unit/compartment/basis contract")
        for key, value in self.truth.items():
            arr = np.array(value, dtype=float, copy=True)
            if arr.shape != times.shape or not np.isfinite(arr).all():
                raise ValueError(f"invalid truth trace {key}")
            arr.setflags(write=False)
            out[key] = arr
        object.__setattr__(self, "times_h", times)
        object.__setattr__(self, "truth", frozen_mapping(out))
        for name in ("variables", "provenance", "diagnostics"):
            object.__setattr__(self, name, frozen_mapping(getattr(self, name)))

    def trace(self, name: str) -> np.ndarray:
        return self.truth[name]

    @property
    def prior_parameters(self) -> tuple[str, ...]:
        return tuple(name for name, p in self.provenance.items()
                     if p.tag in (Tag.ASSERTED, Tag.BOUNDED, Tag.SWEPT))
