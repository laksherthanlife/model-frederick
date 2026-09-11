from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from ..calib.od import ODCalibration
from ..fba.stress_proteostasis import ATP_EQUIVALENTS_PER_PEPTIDE_BOND
from ..kinetic.carotenoid import CarotenoidKinetics, branch_rates
from ..observation import ReporterOptics, observe_od, observe_rfu
from ..bridge.thermodynamic import GAS_CONSTANT_KJ, STANDARD_TEMPERATURE_K
from ..pathway.thermo_gate import Feasibility, StepEnergy, Unresolved, gate_step
from ..photophysics import CITRINE_PKA, citrine_ph_response
from . import oxidative, ph, upr
from .contracts import (
    DYNAMIC_SIGNALS, EXTRACELLULAR, INTRACELLULAR, LEDGERS, UPR_VARIABLES,
    AcidExportParameters, Compartment, ExpressionState, Genotype, ObservationModel,
    Observations, PhysicalState, PoolPartition,
    Protocol, ScientificRefusal, SimulationResult, ThermodynamicRequest, Validity, Variable,
    finite, frozen_mapping, parameter_value, prior,
)
from .integrate import PINNED_METHOD
from .params import Param
from .signalling import TRANSCRIBED_AXES, Signalling, saturation
from .state import TIMESCALE_CATALOGUE


_PRIORS = {
    "mu_max": (0.35, "1/h"), "reference_growth": (0.3, "1/h"),
    "cell_water_l_gdw": (0.002, "L/gDW"), "cells_gdw": (5e10, "cells/gDW"),
    "biomass_carbon": (40.0, "mmol/gDW"), "biomass_nitrogen": (7.0, "mmol/gDW"),
    "growth_atp": (60.0, "mmol/gDW"), "growth_nadph": (6.0, "mmol/gDW"),
    "maintenance_atp": (1.0, "mmol/gDW/h"), "growth_basal_fraction": (0.1, "dimensionless"),
    "temperature_reference": (30.0, "degC"), "temperature_width": (12.0, "degC"),
    "heat_scale": (8.0, "degC"), "pka_basal_fraction": (0.2, "dimensionless"),
    "glucose_transport": (12.0, "mmol/gDW/h"), "ethanol_oxidation": (4.0, "mmol/gDW/h"),
    "upper_glycolysis": (12.0, "mmol/gDW/h"), "lower_glycolysis": (24.0, "mmol/gDW/h"),
    "ppp": (4.0, "mmol/gDW/h"), "pentose_recycle": (4.0, "mmol/gDW/h"),
    "pdh": (15.0, "mmol/gDW/h"), "tca": (10.0, "mmol/gDW/h"),
    "respiration": (30.0, "mmol/gDW/h"), "fermentation": (30.0, "mmol/gDW/h"),
    "adenylate_kinase": (20.0, "mmol/gDW/h"), "respiratory_atp_yield": (1.5, "dimensionless"),
    "carbon_k": (0.05, "mmol/gDW"), "glucose_k": (1.0, "mM"),
    "nitrogen_k": (0.1, "mM"), "oxygen_k": (0.005, "mM"),
    "energy_k": (0.0002, "mmol/gDW"), "redox_k": (0.0002, "mmol/gDW"),
    "ethanol_k": (10.0, "mM"), "ethanol_repression_k": (2.0, "mM"),
    "ethanol_toxicity_k": (500.0, "mM"), "death_baseline": (0.002, "1/h"),
    "death_starvation": (0.02, "1/h"), "death_stress": (0.2, "1/h"),
    "glycogen_synthesis": (2.0, "mmol/gDW/h"), "glycogen_mobilization": (1.0, "mmol/gDW/h"),
    "glycerol_kcat": (1e6, "1/h"), "glycerol_permeability": (0.5, "1/h"),
    "ipp_synthesis": (0.02, "mmol/gDW/h"), "native_isoprenoid": (0.004, "mmol/gDW/h"),
    "native_k": (0.0001, "mmol/gDW"), "crte_kcat": (1000.0, "1/h"),
    "psy_kcat": (400.0, "1/h"), "crti_kcat": (600.0, "1/h"), "lcy_kcat": (400.0, "1/h"),
    "km_crte": (0.0001, "mmol/gDW"), "km_psy": (0.00001, "mmol/gDW"),
    "km_crti": (0.00001, "mmol/gDW"), "km_lcy": (0.00001, "mmol/gDW"),
    "ros_leak_fraction": (0.01, "dimensionless"), "ros_k": (0.0001, "mmol/gDW"),
    "detox_kcat": (1e6, "1/h"), "basal_detox": (0.2, "mmol/gDW/h"),
    "catalase": (0.05, "1/h"), "camp_synthesis": (0.005, "mmol/gDW/h"),
    "camp_hydrolysis": (10.0, "1/h"), "camp_k": (0.0001, "mmol/gDW"),
    "pka_tau": (0.03, "h"), "snf1_tau": (0.1, "h"), "torc1_tau": (0.1, "h"),
    "yap1_tau": (2.0 / 60.0, "h"), "hsf1_tau": (TIMESCALE_CATALOGUE["hsf1_free"].tau_h[0], "h"),
    "hsf_chaperone_k": (1e-7, "mmol/gDW"), "hsf_client_k": (0.02, "dimensionless"),
    "native_unfolding": (0.01, "1/h"), "heat_unfolding": (0.5, "1/h"),
    "ros_unfolding": (0.1, "1/h"), "folding": (3.0, "1/h"),
    "hsp70_fold_k": (1e-7, "mmol/gDW"), "folding_atp": (30.0, "mmol/gDW"),
    "protein_refolding_atp": (10.0, "dimensionless"), "mrna_c_per_nt": (9.5, "dimensionless"),
    "mrna_n_per_nt": (3.75, "dimensionless"), "protein_c_per_aa": (5.0, "dimensionless"),
    "protein_n_per_aa": (1.3, "dimensionless"), "transcription_atp": (2.0, "dimensionless"),
    "translation_atp": (ATP_EQUIVALENTS_PER_PEPTIDE_BOND, "dimensionless"),
    "translation_nadph": (1.0, "dimensionless"),
    "initial_atp": (0.002, "mmol/gDW"), "initial_adp": (0.001, "mmol/gDW"),
    "initial_amp": (0.0002, "mmol/gDW"), "initial_nad": (0.001, "mmol/gDW"),
    "initial_nadh": (0.0003, "mmol/gDW"), "initial_nadp": (0.0003, "mmol/gDW"),
    "initial_nadph": (0.001, "mmol/gDW"), "initial_carbon": (0.01, "mmol/gDW"),
    "initial_protein": (1e-7, "mmol/gDW"), "initial_mrna": (1e-10, "mmol/gDW"),
    "upr_basal_load": (0.2, "dimensionless"), "upr_hill": (4.5, "dimensionless"),
    "upr_basal_share": (0.1, "dimensionless"), "upr_folding": (4.0, "1/h"),
    "upr_atp": (0.05, "mmol/gDW/h"), "er_expression_load": (1e5, "gDW/mmol"),
    "hog_pool": (1e-7, "mmol/gDW"), "buffer_capacity": (150.0, "mM/pH"),
    "resting_ph": (float(ph.RESTING_CYTOSOLIC_PH), "dimensionless"),
    "acid_signal_scale": (1.0, "dimensionless"),
    "yap1_k": (float(oxidative.K_EX_MM), "mM"),
}
_POSITIVE = {
    "reference_growth", "cell_water_l_gdw", "cells_gdw", "biomass_carbon", "biomass_nitrogen",
    "temperature_width", "heat_scale", "carbon_k", "glucose_k", "nitrogen_k", "oxygen_k",
    "energy_k", "redox_k", "ethanol_k", "ethanol_repression_k", "ethanol_toxicity_k",
    "native_k", "km_crte", "km_psy", "km_crti", "km_lcy", "ros_k", "camp_k",
    "pka_tau", "snf1_tau", "torc1_tau", "yap1_tau", "hsf1_tau", "hsf_chaperone_k",
    "hsf_client_k", "hsp70_fold_k", "buffer_capacity", "resting_ph", "acid_signal_scale", "yap1_k",
    "upr_basal_load", "upr_folding", "upr_hill",
}
_FRACTIONS = {"ros_leak_fraction", "growth_basal_fraction", "pka_basal_fraction", "upr_basal_load"}

REGISTERED_AXIS_VARIABLES = MappingProxyType({
    name: variable
    for axis in TRANSCRIBED_AXES.values()
    for name, variable in axis.variables().items()
})
"""The transcribed stress axes' states, named and typed exactly as the integrated ones are.

They are registered HERE and not in ``_Kernel.variables``, which is the state vector, because
not one of them is driven: :data:`mech.signalling.TRANSCRIBED_AXES` carries the reason per
axis. Keeping them in the engine's own vocabulary is what makes the gap checkable rather than
merely stated -- ``_Kernel`` refuses a name that collides with a real state, and refuses an
axis that claims a rate law this engine does not have.
"""


@dataclass(frozen=True)
class EngineParameters:
    values: Mapping[str, Param]
    allow_prior: bool = False
    thermodynamics: ThermodynamicRequest | None = None
    acid_export: AcidExportParameters | None = None
    _resolved: Mapping[str, float] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        missing = set(_PRIORS) - set(self.values)
        if missing:
            raise ScientificRefusal(f"missing biological parameters: {sorted(missing)}; choose priors explicitly or supply measurements")
        if set(self.values) - set(_PRIORS):
            raise ValueError(f"unknown engine parameters: {sorted(set(self.values) - set(_PRIORS))}")
        values = {name: parameter_value(self.values[name], units, positive=name in _POSITIVE,
                                        maximum=1.0 if name in _FRACTIONS else None)
                  for name, (_, units) in _PRIORS.items()}
        if not isinstance(self.allow_prior, bool):
            raise ValueError("allow_prior must be a bool")
        if not self.allow_prior:
            raise ScientificRefusal("this connected model contains unvalidated structural priors; explicit allow_prior=True is required even with measured rate constants")
        if values["upr_basal_load"] >= 1:
            raise ValueError("upr_basal_load must be below one")
        upr.FOLDING_RATE_PER_H.at(values["upr_folding"])
        if not upr.IRE1_HILL_N.bounds[0] <= values["upr_hill"] <= upr.IRE1_HILL_N.bounds[1]:
            raise ValueError("upr_hill is outside the source mechanism's supported bracket")
        if values["resting_ph"] >= 14:
            raise ValueError("resting_ph must be below 14")
        if self.thermodynamics is not None and not isinstance(self.thermodynamics, ThermodynamicRequest):
            raise ValueError("thermodynamics requires a ThermodynamicRequest")
        if self.acid_export is not None and not isinstance(self.acid_export, AcidExportParameters):
            raise ValueError("acid_export requires explicit AcidExportParameters")
        object.__setattr__(self, "values", frozen_mapping(self.values))
        object.__setattr__(self, "_resolved", frozen_mapping(values))

    @classmethod
    def prior(cls, *, acid_export: AcidExportParameters | None = None) -> EngineParameters:
        values = {name: prior(name, value, units, name.replace("_", " "))
                  for name, (value, units) in _PRIORS.items()}
        values["resting_ph"] = ph.RESTING_CYTOSOLIC_PH
        values["yap1_k"] = replace(oxidative.K_EX_MM, units="mM")
        values["hsf1_tau"] = Param.asserted(
            "hsf1_tau", TIMESCALE_CATALOGUE["hsf1_free"].tau_h[0], "h",
            TIMESCALE_CATALOGUE["hsf1_free"].source)
        values["translation_atp"] = Param.asserted(
            "translation_atp", ATP_EQUIVALENTS_PER_PEPTIDE_BOND, "dimensionless",
            "Stoichiometric convention imported from fba/stress_proteostasis.py::ATP_EQUIVALENTS_PER_PEPTIDE_BOND; Stouthamer 1973 PMID 4148026, four ATP equivalents per residue")
        return cls(values, allow_prior=True, acid_export=acid_export)

    def with_values(self, **values: Param) -> EngineParameters:
        return replace(self, values={**self.values, **values})


def initialize(parameters: EngineParameters, genotype: Genotype, *, volume_l: float,
               biomass_gdw_l: float, medium_mM: Mapping[str, float], time_h: float = 0.0) -> PhysicalState:
    state = PhysicalState.from_concentrations(volume_l=volume_l, biomass_gdw_l=biomass_gdw_l,
                                             medium_mM=medium_mM, time_h=time_h)
    if state.biomass_gdw == 0:
        return state
    p = parameters._resolved
    signalling = Signalling(p, genotype)
    pools = dict(state.intracellular_mmol)
    for name in ("atp", "adp", "amp", "nad", "nadh", "nadp", "nadph"):
        pools[name] = p[f"initial_{name}"] * state.biomass_gdw
    for name in ("glucose", "triose", "pyruvate", "acetyl", "pentose"):
        pools[name] = p["initial_carbon"] * state.biomass_gdw
    expression = {
        gene.name: ExpressionState(p["initial_mrna"] * state.biomass_gdw * gene.copies,
                                   p["initial_protein"] * state.biomass_gdw * gene.copies,
                                   p["initial_protein"] * state.biomass_gdw * gene.copies)
        for gene in genotype.genes
    }
    signals = dict(state.signals)
    signals["hog1"] = (signalling.hog.reference["Hog1PP"] / signalling.hog.total
                       * genotype.signal_activity.get("hog1", 1.0))
    return replace(state, intracellular_mmol=pools, expression=expression, signals=signals,
                   upr={name: signalling.upr_basal[name] for name in UPR_VARIABLES})


class _Kernel:
    def __init__(self, parameters: EngineParameters, genotype: Genotype):
        self.parameters = parameters
        self.p = parameters._resolved
        self.genotype = genotype
        self.signalling = Signalling(self.p, genotype)
        self.variables = {
            "volume_l": Variable("L", Compartment.REACTOR, "volume"),
            "biomass_gdw": Variable("gDW", Compartment.REACTOR, "structural biomass", self.p["biomass_carbon"], self.p["biomass_nitrogen"]),
            "dead_biomass_gdw": Variable("gDW", Compartment.DEAD, "structural biomass", self.p["biomass_carbon"], self.p["biomass_nitrogen"]),
            "cell_volume_ratio": Variable("dimensionless", Compartment.CELL_WATER, "relative osmotically accessible volume"),
            "unfolded_fraction": Variable("dimensionless", Compartment.REGULATORY, "subset of native proteome"),
        }
        for prefix, species, compartment in (("external", EXTRACELLULAR, Compartment.MEDIUM),
                                              ("internal", INTRACELLULAR, Compartment.CELL_WATER),
                                              ("dead", INTRACELLULAR, Compartment.DEAD)):
            for name, (carbon, nitrogen) in species.items():
                self.variables[f"{prefix}.{name}"] = Variable("mmol", compartment, "amount", carbon, nitrogen)
        for name in DYNAMIC_SIGNALS:
            self.variables[f"signal.{name}"] = Variable("dimensionless", Compartment.REGULATORY, "active fraction")
        for name in UPR_VARIABLES:
            self.variables[f"upr.{name}"] = Variable("relative", Compartment.ER, "source-normalized native allocation; not additional biomass")
        for gene in genotype.genes:
            for kind, carbon, nitrogen in (
                ("mrna", gene.transcript_nt * self.p["mrna_c_per_nt"], gene.transcript_nt * self.p["mrna_n_per_nt"]),
                ("protein", gene.protein_aa * self.p["protein_c_per_aa"], gene.protein_aa * self.p["protein_n_per_aa"]),
                ("active", 0.0, 0.0),
            ):
                self.variables[f"gene.{gene.name}.{kind}"] = Variable(
                    "mmol", Compartment.CELL_WATER, "subset of protein" if kind == "active" else "amount", carbon, nitrogen)
        for name in LEDGERS:
            units = "gDW" if name in ("biomass_formed", "biomass_died", "biomass_washed_out") else "mmol"
            self.variables[f"ledger.{name}"] = Variable(units, Compartment.LEDGER, "cumulative extent")
        self.variables["ledger.expression_carbon_dead"] = Variable("mmol", Compartment.DEAD, "carbon equivalents", 1.0)
        self.variables["ledger.expression_nitrogen_dead"] = Variable("mmol", Compartment.DEAD, "nitrogen equivalents", 0.0, 1.0)
        self.names = tuple(self.variables)
        self.index = {name: i for i, name in enumerate(self.names)}
        self.carbon = np.array([v.carbon for v in self.variables.values()])
        self.nitrogen = np.array([v.nitrogen for v in self.variables.values()])
        self.physical = np.array([i for i, v in enumerate(self.variables.values())
                                  if v.basis in ("amount", "structural biomass", "subset of protein", "carbon equivalents", "nitrogen equivalents")])
        self.death_pairs = tuple((self.index[f"internal.{name}"], self.index[f"dead.{name}"])
                                 for name in INTRACELLULAR)
        self.expression = tuple((g, *(self.index[f"gene.{g.name}.{kind}"] for kind in ("mrna", "protein", "active")))
                                for g in genotype.genes)
        self.upr_slice = np.array([self.index[f"upr.{name}"] for name in UPR_VARIABLES])
        self.ext_indices = np.array([self.index[f"external.{name}"] for name in EXTRACELLULAR])
        self.in_indices = np.array([self.index[f"internal.{name}"] for name in INTRACELLULAR])
        self.stoichiometry = self._reactions()
        self.reaction_names = tuple(self.stoichiometry)
        collisions = sorted(set(REGISTERED_AXIS_VARIABLES) & set(self.names))
        if collisions:
            raise RuntimeError(f"registered axis states collide with integrated states: {collisions}")
        for axis in TRANSCRIBED_AXES.values():
            unknown = sorted(set(axis.readers) - set(self.reaction_names))
            if unknown:
                raise ScientificRefusal(
                    f"axis {axis.axis!r} claims rate laws this engine does not have: {unknown}")
        self.matrix = np.zeros((len(self.names), len(self.reaction_names)))
        for j, reaction in enumerate(self.stoichiometry.values()):
            for name, coefficient in reaction.items():
                self.matrix[self.index[name], j] = coefficient
        carbon = self.carbon.copy()
        carbon[self.index["ledger.co2"]] = 1.0
        if not np.allclose(carbon @ self.matrix, 0.0, atol=1e-10, rtol=0):
            raise RuntimeError("reaction definitions violate carbon conservation")
        if not np.allclose(self.nitrogen @ self.matrix, 0.0, atol=1e-10, rtol=0):
            raise RuntimeError("reaction definitions violate nitrogen conservation")
        self.branch_shape = CarotenoidKinetics(1.0, self.p["km_psy"], 1.0, self.p["km_crti"],
                                               1.0, self.p["km_lcy"], crtyb_pool=1.0)
        self.acid_k = ph.entry_rate_constant_per_h("acetic")
        self.acid_pka = float(ph.PKA_ACETIC)

    def _reactions(self):
        p = self.p
        reactions = {
            "glucose_transport": {"external.glucose": -1, "internal.glucose": 1},
            "upper_glycolysis": {"internal.glucose": -1, "internal.triose": 2, "internal.atp": -2, "internal.adp": 2},
            "lower_glycolysis": {"internal.triose": -1, "internal.pyruvate": 1, "internal.nad": -1, "internal.nadh": 1, "internal.adp": -2, "internal.atp": 2},
            "ppp": {"internal.glucose": -1, "internal.pentose": 1, "ledger.co2": 1, "internal.nadp": -2, "internal.nadph": 2, "internal.atp": -1, "internal.adp": 1},
            "pentose_recycle": {"internal.pentose": -1, "internal.triose": 5 / 3},
            "pdh": {"internal.pyruvate": -1, "internal.acetyl": 1, "ledger.co2": 1, "internal.nad": -1, "internal.nadh": 1},
            "tca": {"internal.acetyl": -1, "ledger.co2": 2, "internal.nad": -4, "internal.nadh": 4, "internal.adp": -1, "internal.atp": 1},
            "fermentation": {"internal.pyruvate": -1, "external.ethanol": 1, "ledger.co2": 1, "internal.nadh": -1, "internal.nad": 1},
            "ethanol_oxidation": {"external.ethanol": -1, "internal.acetyl": 1, "internal.nad": -2, "internal.nadh": 2, "internal.atp": -2, "internal.adp": 2},
            "respiration": {"internal.nadh": -1, "internal.nad": 1, "external.oxygen": -0.5, "internal.adp": -p["respiratory_atp_yield"], "internal.atp": p["respiratory_atp_yield"]},
            "ros_leak": {"internal.nadh": -1, "internal.nad": 1, "external.oxygen": -1, "internal.ros": 1},
            "detox": {"internal.ros": -1, "internal.nadph": -1, "internal.nadp": 1},
            "catalase": {"internal.ros": -1, "external.oxygen": 0.5},
            "peroxide_entry": {"external.peroxide": -1, "internal.ros": 1},
            "glycogen_synthesis": {"internal.glucose": -1, "internal.glycogen": 1, "internal.atp": -1, "internal.adp": 1},
            "glycogen_mobilization": {"internal.glycogen": -1, "internal.glucose": 1},
            "glycerol_synthesis": {"internal.triose": -1, "internal.glycerol": 1, "internal.nadh": -1, "internal.nad": 1},
            "glycerol_export": {"internal.glycerol": -1, "external.glycerol": 1},
            "acid_transport": {"external.acetate": -1, "internal.acetate": 1},
            "ipp_synthesis": {"internal.acetyl": -3, "internal.ipp": 1, "ledger.co2": 1, "internal.atp": -3, "internal.adp": 3, "internal.nadph": -2, "internal.nadp": 2},
            "native_isoprenoid": {"internal.ipp": -1, "internal.native_isoprenoid": 1},
            "crte": {"internal.ipp": -4, "internal.ggpp": 1},
            "psy": {"internal.ggpp": -2, "internal.phytoene": 1},
            "crti": {"internal.phytoene": -1, "internal.lycopene": 1, "external.oxygen": -4, "internal.ros": 4},
            "lcy": {"internal.lycopene": -1, "internal.beta_carotene": 1},
            "camp_synthesis": {"internal.atp": -1, "internal.camp": 1},
            "camp_hydrolysis": {"internal.camp": -1, "internal.amp": 1},
            "adenylate_kinase": {"internal.atp": -1, "internal.amp": -1, "internal.adp": 2},
            "maintenance": {"internal.atp": -1, "internal.adp": 1},
            "regulatory_cost": {"internal.atp": -1, "internal.adp": 1},
        }
        cofactors = {"adp": p["initial_atp"] + p["initial_adp"] + p["initial_amp"],
                     "nad": p["initial_nad"] + p["initial_nadh"],
                     "nadp": p["initial_nadp"] + p["initial_nadph"]}
        cofactor_c = sum(INTRACELLULAR[k][0] * v for k, v in cofactors.items())
        cofactor_n = sum(INTRACELLULAR[k][1] * v for k, v in cofactors.items())
        reactions["growth"] = {
            "biomass_gdw": 1, "internal.acetyl": -(p["biomass_carbon"] + cofactor_c) / 2,
            "external.nitrogen": -(p["biomass_nitrogen"] + cofactor_n),
            "internal.atp": -p["growth_atp"], "internal.adp": p["growth_atp"] + cofactors["adp"],
            "internal.nadph": -p["growth_nadph"], "internal.nadp": p["growth_nadph"] + cofactors["nadp"],
            "internal.nad": cofactors["nad"],
        }
        export = self.parameters.acid_export
        if export is not None:
            atp_cost = ph.proton_pumping_atp(1.0, atp_per_proton=float(export.atp_per_proton)) + float(export.atp_per_anion)
            reactions["acid_export"] = {
                "internal.acetate": -1, "external.acetate": 1,
                "internal.atp": -atp_cost, "internal.adp": atp_cost,
            }
        for gene in self.genotype.genes:
            for kind, cost, redox in (("mrna", p["transcription_atp"] * gene.transcript_nt, 0.0),
                                      ("protein", p["translation_atp"] * gene.protein_aa,
                                       p["translation_nadph"] * gene.protein_aa)):
                key = f"gene.{gene.name}.{kind}"
                var = self.variables[key]
                reactions[f"{gene.name}.{kind}_synthesis"] = {
                    key: 1, "internal.acetyl": -var.carbon / 2, "external.nitrogen": -var.nitrogen,
                    "internal.atp": -cost, "internal.adp": cost,
                    "internal.nadph": -redox, "internal.nadp": redox,
                }
                reactions[f"{gene.name}.{kind}_decay"] = {
                    key: -1, "internal.acetyl": var.carbon / 2, "external.nitrogen": var.nitrogen,
                }
            reactions[f"{gene.name}.folding"] = {
                f"gene.{gene.name}.active": 1,
                "internal.atp": -p["protein_refolding_atp"], "internal.adp": p["protein_refolding_atp"],
            }
        return reactions

    def pack(self, state: PhysicalState):
        unknown = set(state.expression) - {g.name for g in self.genotype.genes}
        if unknown:
            raise ValueError(f"initial expression contains unknown genes {sorted(unknown)}")
        values = {k: getattr(state, k) for k in ("volume_l", "biomass_gdw", "dead_biomass_gdw",
                                               "cell_volume_ratio", "unfolded_fraction")}
        for prefix, mapping in (("external", state.extracellular_mmol), ("internal", state.intracellular_mmol),
                                ("dead", state.dead_intracellular_mmol), ("signal", state.signals),
                                ("upr", state.upr), ("ledger", state.ledgers)):
            values.update({f"{prefix}.{name}": value for name, value in mapping.items()})
        for gene in self.genotype.genes:
            expression = state.expression.get(gene.name, ExpressionState())
            values.update({f"gene.{gene.name}.{kind}": getattr(expression, f"{kind}_mmol")
                           for kind in ("mrna", "protein", "active")})
        return np.array([values[name] for name in self.names])

    def unpack(self, time_h, vector, receipts):
        y = dict(zip(self.names, vector, strict=True))
        return PhysicalState(
            time_h, y["volume_l"], y["biomass_gdw"],
            {n: y[f"external.{n}"] for n in EXTRACELLULAR},
            {n: y[f"internal.{n}"] for n in INTRACELLULAR}, y["dead_biomass_gdw"],
            {n: y[f"dead.{n}"] for n in INTRACELLULAR},
            {g.name: ExpressionState(*(y[f"gene.{g.name}.{kind}"] for kind in ("mrna", "protein", "active")))
             for g in self.genotype.genes},
            {n: y[f"signal.{n}"] for n in DYNAMIC_SIGNALS},
            {n: y[f"upr.{n}"] for n in UPR_VARIABLES}, y["cell_volume_ratio"],
            y["unfolded_fraction"], {n: y[f"ledger.{n}"] for n in LEDGERS}, tuple(receipts),
        )

    def inventories(self, y):
        return (math.fsum(float(c * value) for c, value in zip(self.carbon, y, strict=True)),
                math.fsum(float(n * value) for n, value in zip(self.nitrogen, y, strict=True)))

    def cytosolic_ph(self, y):
        biomass = max(y[self.index["biomass_gdw"]], 0.0)
        acid = max(y[self.index["internal.acetate"]], 0.0)
        if biomass == 0 or acid == 0:
            return self.p["resting_ph"]
        acid_reference_mM = acid / (biomass * self.p["cell_water_l_gdw"])
        return brentq(lambda x: ph.charge_imbalance_mM(
            acid_reference_mM, x, self.p["buffer_capacity"], self.acid_pka,
            resting_ph=self.p["resting_ph"]), 1e-12, self.p["resting_ph"], xtol=1e-12)

    def thermodynamic_report(self, name, y, control):
        request = self.parameters.thermodynamics
        temperature = control.temperature_c + 273.15
        if not math.isclose(temperature, STANDARD_TEMPERATURE_K, abs_tol=1e-8, rel_tol=0):
            return StepEnergy(name, "", Feasibility.CANNOT_SAY, None, None, {},
                              "temperature-matched transformed energies or reaction enthalpies are missing", None, Unresolved.NO_ENERGY)
        if name == "adenylate_kinase" and not math.isclose(
                self.cytosolic_ph(y), request.data.ph, abs_tol=1e-8, rel_tol=0):
            return StepEnergy(name, "", Feasibility.CANNOT_SAY, None, None, {},
                              "cytosolic pH does not match the transformed energy table; a pH-resolved transformation is required", None, Unresolved.NO_ENERGY)
        water_l = max(y[self.index["biomass_gdw"]], 0.0) * self.p["cell_water_l_gdw"] * max(y[self.index["cell_volume_ratio"]], 0.0)
        stoichiometry, activities = {}, {}
        for key, coefficient in self.stoichiometry[name].items():
            if coefficient == 0:
                continue
            compartment, species = key.split(".", 1)
            if compartment != "internal" or species not in request.metabolite_ids:
                raise ScientificRefusal(f"{name} lacks a complete intracellular molecular identity mapping")
            metabolite = request.metabolite_ids[species]
            stoichiometry[metabolite] = coefficient
            concentration = max(y[self.index[key]], 0.0) / water_l / 1000.0 if water_l > 0 else 0.0
            coefficient = float(request.activity_coefficients[species]) if species in request.activity_coefficients else 1.0
            activities[metabolite] = concentration * coefficient
        return gate_step(stoichiometry, activities, request.data, node=name,
                         temperature_k=temperature, uncertainty_sigma=request.uncertainty_sigma)

    def constrain_thermodynamic_rate(self, name, rate, y, control):
        report = self.thermodynamic_report(name, y, control)
        if report.feasibility == Feasibility.CANNOT_SAY:
            if report.reason == Unresolved.NO_CONCENTRATION and report.dg0_kj_per_mol is not None:
                if not math.isfinite(report.dg0_kj_per_mol):
                    raise ScientificRefusal(f"nonfinite thermodynamic energy for {name}")
                if rate == 0:
                    return 0.0
                direction = 1.0 if rate > 0 else -1.0
                reactants, products = [], []
                request = self.parameters.thermodynamics
                for key, coefficient in self.stoichiometry[name].items():
                    metabolite = request.metabolite_ids[key.split(".", 1)[1]]
                    value = report.concentrations_m[metabolite]
                    (reactants if coefficient * direction < 0 else products).append(value)
                if all(value > 0 for value in reactants) and any(value == 0 for value in products):
                    return rate
            raise ScientificRefusal(f"thermodynamic energy/activities for {name} cannot constrain its flux: {report.reason}; {report.note}")
        if rate == 0:
            return 0.0
        dg = float(report.dg_kj_per_mol) * (1.0 if rate > 0 else -1.0)
        if not math.isfinite(dg):
            raise ScientificRefusal(f"nonfinite thermodynamic energy for {name}")
        if dg >= 0:
            return 0.0
        return rate * (-math.expm1(dg / (GAS_CONSTANT_KJ * (control.temperature_c + 273.15))))

    def physiology(self, y, control):
        y = np.maximum(y, 0.0)
        p, idx = self.p, self.index
        volume, biomass = y[idx["volume_l"]], y[idx["biomass_gdw"]]
        external = dict(zip(EXTRACELLULAR, y[self.ext_indices] / volume, strict=True))
        content = dict(zip(INTRACELLULAR, y[self.in_indices] / biomass if biomass > 0
                           else np.zeros(len(INTRACELLULAR)), strict=True))
        signals = {name: float(np.clip(y[idx[f"signal.{name}"]], 0.0, 1.0)) for name in DYNAMIC_SIGNALS}
        roles = dict.fromkeys(("crte", "crti", "crtyb", "gpd1", "hsp70", "antioxidant"), 0.0)
        for gene, _, _, active in self.expression:
            if gene.role in roles and biomass > 0:
                roles[gene.role] += y[active] / biomass
        energy = saturation(content["atp"], p["energy_k"])
        redox = saturation(content["nadph"], p["redox_k"])
        unfolded = float(np.clip(y[idx["unfolded_fraction"]], 0.0, 1.0))
        thermal = math.exp(-((control.temperature_c - p["temperature_reference"]) / p["temperature_width"]) ** 2)
        heat = max(control.temperature_c - p["temperature_reference"], 0.0) / p["heat_scale"]
        ph_c = self.cytosolic_ph(y)
        ph_factor = min(1.0, ph.orij_growth_per_h(ph_c) / ph.orij_growth_per_h(p["resting_ph"]))
        growth = (p["mu_max"] * saturation(content["acetyl"], p["carbon_k"])
                  * saturation(external["nitrogen"], p["nitrogen_k"]) * energy * redox
                  * (p["growth_basal_fraction"] + (1.0 - p["growth_basal_fraction"]) * signals["torc1"])
                  * thermal * ph_factor * (1.0 - unfolded)) if biomass > 0 else 0.0
        death = (p["death_baseline"] + p["death_starvation"] * (1.0 - energy)
                 + p["death_stress"] * (saturation(content["ros"], p["ros_k"])
                                        + saturation(external["ethanol"], p["ethanol_toxicity_k"])
                                        + heat ** 2 + unfolded)) if biomass > 0 else 0.0
        return external, content, signals, roles, energy, redox, thermal, growth, death, ph_c

    def evaluate(self, time_h, state, control):
        y = np.maximum(state, 0.0)
        p, idx = self.p, self.index
        volume, biomass = y[idx["volume_l"]], y[idx["biomass_gdw"]]
        e, c, signals, roles, energy, redox, thermal, growth, death, ph_c = self.physiology(y, control)
        rates = dict.fromkeys(self.reaction_names, 0.0)
        dy = np.zeros(len(self.names))
        drivers = {**signals, "upre": 0.0, "acid": 0.0}
        if biomass > 0:
            sat_c = {name: saturation(value, p["carbon_k"]) for name, value in c.items()}
            nad = saturation(c["nad"], p["redox_k"])
            nadh = saturation(c["nadh"], p["redox_k"])
            nadp = saturation(c["nadp"], p["redox_k"])
            adp = saturation(c["adp"], p["energy_k"])
            nitrogen = saturation(e["nitrogen"], p["nitrogen_k"])
            oxygen = saturation(e["oxygen"], p["oxygen_k"])
            g_signal = saturation(e["glucose"], p["glucose_k"])
            rate_scale = biomass * thermal
            def flux(name, factor):
                return p[name] * factor * rate_scale
            respiration = flux("respiration", nadh * oxygen * adp)
            expression_content = sum(y[protein] for _, _, protein, _ in self.expression) / biomass
            signal_rates = self.signalling.evaluate(
                time_h, signals={name: state[idx[f"signal.{name}"]] for name in DYNAMIC_SIGNALS},
                upr_state=y[self.upr_slice],
                cell_volume_ratio=y[idx["cell_volume_ratio"]], unfolded=y[idx["unfolded_fraction"]],
                extracellular=e, intracellular_content=c, hsp70_content=roles["hsp70"],
                expression_content=expression_content, growth=growth, energy=energy, redox=redox,
                ph_c=ph_c, control=control,
            )
            drivers = signal_rates.drivers
            for name, value in signal_rates.derivatives.items():
                dy[idx[f"signal.{name}"]] = value
            dy[self.upr_slice] = signal_rates.upr_derivatives
            dy[idx["cell_volume_ratio"]] = signal_rates.volume_rate
            dy[idx["unfolded_fraction"]] = signal_rates.unfolded_rate - growth * y[idx["unfolded_fraction"]]
            rates.update({
                "glucose_transport": flux("glucose_transport", g_signal * (p["pka_basal_fraction"] + (1.0 - p["pka_basal_fraction"]) * signals["pka"])),
                "upper_glycolysis": flux("upper_glycolysis", sat_c["glucose"] * energy),
                "lower_glycolysis": flux("lower_glycolysis", sat_c["triose"] * nad * adp),
                "ppp": flux("ppp", sat_c["glucose"] * nadp * energy),
                "pentose_recycle": flux("pentose_recycle", sat_c["pentose"]),
                "pdh": flux("pdh", sat_c["pyruvate"] * nad),
                "tca": flux("tca", sat_c["acetyl"] * nad * adp),
                "fermentation": flux("fermentation", sat_c["pyruvate"] * nadh),
                "ethanol_oxidation": flux("ethanol_oxidation", saturation(e["ethanol"], p["ethanol_k"]) * nad * energy * signals["snf1"] * (1.0 - saturation(e["glucose"], p["ethanol_repression_k"]))),
                "respiration": respiration * (1.0 - p["ros_leak_fraction"]),
                "ros_leak": respiration * p["ros_leak_fraction"],
                "detox": (p["basal_detox"] + p["detox_kcat"] * roles["antioxidant"]) * saturation(c["ros"], p["ros_k"]) * redox * rate_scale,
                "catalase": p["catalase"] * y[idx["internal.ros"]],
                "peroxide_entry": oxidative.dose_decay_from_density_per_h(p["cells_gdw"] * biomass / (1000.0 * volume)) * y[idx["external.peroxide"]],
                "glycogen_synthesis": flux("glycogen_synthesis", sat_c["glucose"] * energy * (1.0 - signals["torc1"])),
                "glycogen_mobilization": flux("glycogen_mobilization", sat_c["glycogen"] * signals["snf1"]),
                "glycerol_synthesis": p["glycerol_kcat"] * roles["gpd1"] * sat_c["triose"] * nadh * rate_scale,
                "ipp_synthesis": flux("ipp_synthesis", sat_c["acetyl"] * energy * redox),
                "native_isoprenoid": flux("native_isoprenoid", saturation(c["ipp"], p["native_k"])),
                "crte": p["crte_kcat"] * roles["crte"] * saturation(c["ipp"], p["km_crte"]) * rate_scale,
                "camp_synthesis": flux("camp_synthesis", g_signal * energy),
                "camp_hydrolysis": p["camp_hydrolysis"] * y[idx["internal.camp"]],
                "adenylate_kinase": flux("adenylate_kinase", energy * saturation(c["amp"], p["energy_k"]) - adp * adp),
                "maintenance": p["maintenance_atp"] * energy * biomass,
                "regulatory_cost": biomass * (signal_rates.native_repair_atp + 2.0 * p["hog_pool"] * signal_rates.hog_phosphorylation),
                "growth": growth * biomass,
            })
            shape = branch_rates(c, self.branch_shape)
            rates.update(psy=shape["psy"] * p["psy_kcat"] * roles["crtyb"] * rate_scale,
                         crti=shape["crti"] * p["crti_kcat"] * roles["crti"] * oxygen * rate_scale,
                         lcy=shape["lcy"] * p["lcy_kcat"] * roles["crtyb"] * rate_scale)
            ratio = y[idx["cell_volume_ratio"]]
            water = biomass * p["cell_water_l_gdw"] * ratio
            rates["glycerol_export"] = p["glycerol_permeability"] * (1.0 - signals["hog1"]) * (
                y[idx["internal.glycerol"]] - water * e["glycerol"])
            acid_reference = c["acetate"] / p["cell_water_l_gdw"]
            acid_rhs = ph.weak_acid_rhs(
                time_h, (acid_reference, ph_c),
                ah_out_mM=ratio * ph.undissociated_outside(e["acetate"], control.ph, self.acid_pka),
                k_entry_per_h=self.acid_k, pka=self.acid_pka, beta_mM_per_ph=p["buffer_capacity"],
                growth_rate_per_h=growth, resting_ph=p["resting_ph"],
            )
            rates["acid_transport"] = (acid_rhs[0] + growth * acid_reference) * biomass * p["cell_water_l_gdw"]
            export = self.parameters.acid_export
            if export is not None:
                anion_mM = acid_reference / ratio * (1.0 - ph.undissociated_fraction(ph_c, self.acid_pka))
                rates["acid_export"] = float(export.vmax) * biomass * saturation(anion_mM, float(export.km_anion)) * energy
            expression_supply = energy * redox * nitrogen * sat_c["acetyl"] * thermal
            hsp = saturation(roles["hsp70"], p["hsp70_fold_k"])
            unfold_rate = (p["native_unfolding"] + p["heat_unfolding"] * (
                max(control.temperature_c - p["temperature_reference"], 0.0) / p["heat_scale"]) ** 2
                + p["ros_unfolding"] * saturation(c["ros"], p["ros_k"]))
            for gene, mrna, protein, active in self.expression:
                promoter = gene.promoter.activity(drivers)
                rates[f"{gene.name}.mrna_synthesis"] = float(gene.transcription_mmol_gdw_h) * gene.copies * promoter * biomass * expression_supply
                rates[f"{gene.name}.mrna_decay"] = float(gene.mrna_decay_per_h) * y[mrna]
                rates[f"{gene.name}.protein_synthesis"] = float(gene.translation_per_h) * y[mrna] * expression_supply
                rates[f"{gene.name}.protein_decay"] = float(gene.protein_decay_per_h) * y[protein]
                rates[f"{gene.name}.folding"] = p["folding"] * energy * hsp * max(y[protein] - y[active], 0.0)
                dy[active] -= (float(gene.protein_decay_per_h) + unfold_rate) * y[active]
        request = self.parameters.thermodynamics
        if request is not None and request.mode == "constrain_supported":
            for name in request.reactions:
                rates[name] = self.constrain_thermodynamic_rate(name, rates[name], y, control)
        rate_vector = np.array([rates[name] for name in self.reaction_names])
        dy += self.matrix @ rate_vector
        outflow = control.outflow_l_h / volume
        dy[idx["volume_l"]] = control.feed_l_h - control.outflow_l_h
        dy[self.physical] -= outflow * y[self.physical]
        for name in EXTRACELLULAR:
            dy[idx[f"external.{name}"]] += control.feed_l_h * control.feed_mM.get(name, 0.0)
        oxygen_transfer = control.oxygen_transfer_per_h * volume * (control.oxygen_saturation_mM - e["oxygen"])
        dy[idx["external.oxygen"]] += oxygen_transfer
        live_death = death * y[idx["biomass_gdw"]]
        dy[idx["biomass_gdw"]] -= live_death
        dy[idx["dead_biomass_gdw"]] += live_death
        for live, dead in self.death_pairs:
            transfer = death * y[live]
            dy[live] -= transfer
            dy[dead] += transfer
        for _, mrna, protein, active in self.expression:
            for i in (mrna, protein, active):
                dy[i] -= death * y[i]
            dy[idx["ledger.expression_carbon_dead"]] += death * (self.carbon[mrna] * y[mrna] + self.carbon[protein] * y[protein])
            dy[idx["ledger.expression_nitrogen_dead"]] += death * (self.nitrogen[mrna] * y[mrna] + self.nitrogen[protein] * y[protein])
        c_in = sum(EXTRACELLULAR[name][0] * value for name, value in control.feed_mM.items()) * control.feed_l_h
        n_in = sum(EXTRACELLULAR[name][1] * value for name, value in control.feed_mM.items()) * control.feed_l_h
        carbon, nitrogen = self.inventories(y)
        for name, value in (("carbon_in", c_in), ("nitrogen_in", n_in),
                            ("carbon_out", outflow * carbon), ("nitrogen_out", outflow * nitrogen),
                            ("biomass_formed", growth * biomass), ("biomass_died", live_death),
                            ("biomass_washed_out", outflow * biomass), ("oxygen_transferred", oxygen_transfer)):
            dy[idx[f"ledger.{name}"]] += value
        for species in ("atp", "nadph"):
            changes = self.matrix[idx[f"internal.{species}"]] * rate_vector
            dy[idx[f"ledger.{species}_produced"]] = np.sum(np.maximum(changes, 0.0))
            dy[idx[f"ledger.{species}_consumed"]] = np.sum(np.maximum(-changes, 0.0))
        exported = rates.get("acid_export", 0.0)
        dy[idx["ledger.protons_pumped"]] = exported
        dy[idx["ledger.acetate_exported"]] = exported
        diagnostics = {"growth_per_h": growth, "death_per_h": death, "washout_per_h": outflow,
                       "net_live_per_h": growth - death - outflow, "ph_c": ph_c,
                       "oxygen_transfer_mmol_h": oxygen_transfer, "energy_availability": energy,
                       "osmolarity_mOsm": sum(v for k, v in e.items() if k != "oxygen")}
        if self.parameters.acid_export is not None:
            export = self.parameters.acid_export
            diagnostics["ph_pump_atp_mmol_h"] = ph.proton_pumping_atp(exported, atp_per_proton=float(export.atp_per_proton))
            diagnostics["anion_export_atp_mmol_h"] = exported * float(export.atp_per_anion)
        diagnostics.update({f"driver.{name}": value for name, value in drivers.items()})
        diagnostics.update({f"flux.{name}": value for name, value in rates.items()})
        for gene, _, _, _ in self.expression:
            diagnostics[f"promoter.{gene.name}"] = gene.promoter.activity(drivers)
        return dy, diagnostics

    def apply_event(self, y, event):
        state = y.copy()
        volume = state[self.index["volume_l"]]
        if event.withdraw_l >= volume:
            raise ValueError(f"event {event.name!r} empties the reactor before its addition")
        fraction = event.withdraw_l / volume
        carbon, nitrogen = self.inventories(state)
        state[self.index["ledger.carbon_out"]] += fraction * carbon
        state[self.index["ledger.nitrogen_out"]] += fraction * nitrogen
        state[self.index["ledger.biomass_washed_out"]] += fraction * state[self.index["biomass_gdw"]]
        state[self.physical] *= 1.0 - fraction
        state[self.index["volume_l"]] = volume - event.withdraw_l + event.add_volume_l
        for name, amount in event.add_mmol.items():
            state[self.index[f"external.{name}"]] += amount
            state[self.index["ledger.carbon_in"]] += EXTRACELLULAR[name][0] * amount
            state[self.index["ledger.nitrogen_in"]] += EXTRACELLULAR[name][1] * amount
        return state, 1.0 - fraction


def _reachable_coordinates(kernel, y, control):
    physical = set(kernel.physical)
    present = {i for i, value in enumerate(y) if value > 0 or i not in physical}
    present.update(kernel.index[f"external.{name}"] for name, concentration in control.feed_mM.items()
                   if concentration > 0 and control.feed_l_h > 0)
    if control.oxygen_transfer_per_h > 0 and control.oxygen_saturation_mM > 0:
        present.add(kernel.index["external.oxygen"])
    edges = []
    catalysts = {}
    disabled = set()
    for name, role in (("crte", "crte"), ("psy", "crtyb"), ("crti", "crti"),
                       ("lcy", "crtyb"), ("glycerol_synthesis", "gpd1")):
        catalysts[name] = [{active} for gene, _, _, active in kernel.expression if gene.role == role]
    for gene, mrna, protein, _ in kernel.expression:
        if gene.copies == 0 or float(gene.transcription_mmol_gdw_h) == 0:
            disabled.add(f"{gene.name}.mrna_synthesis")
        if float(gene.translation_per_h) == 0:
            disabled.add(f"{gene.name}.protein_synthesis")
        catalysts[f"{gene.name}.protein_synthesis"] = [{mrna}]
        catalysts[f"{gene.name}.folding"] = [{protein}]
    for name, reaction in kernel.stoichiometry.items():
        if name in disabled:
            continue
        consumes = {kernel.index[key] for key, coefficient in reaction.items() if coefficient < 0}
        produces = {kernel.index[key] for key, coefficient in reaction.items() if coefficient > 0}
        for catalyst in catalysts.get(name, [set()]):
            edges.append((consumes | catalyst | {kernel.index["biomass_gdw"]}, produces))
        if name in ("glycerol_export", "acid_transport", "adenylate_kinase"):
            edges.append((produces | {kernel.index["biomass_gdw"]}, consumes))
    edges.extend(({live}, {dead}) for live, dead in kernel.death_pairs)
    edges.append(({kernel.index["biomass_gdw"]}, {kernel.index["dead_biomass_gdw"]}))
    for _, mrna, protein, _ in kernel.expression:
        for i in (mrna, protein):
            edges.append(({i}, {kernel.index["ledger.expression_carbon_dead"],
                                kernel.index["ledger.expression_nitrogen_dead"]}))
    while True:
        old = len(present)
        for consumes, produces in edges:
            if consumes <= present:
                present.update(produces)
        if len(present) == old:
            break
    return np.array(sorted(present), dtype=int)


@dataclass(frozen=True)
class _DenseState:
    solution: object
    reference: np.ndarray
    coordinates: np.ndarray

    def restore(self, values):
        values = np.asarray(values)
        restored = self.reference.copy() if values.ndim == 1 else np.repeat(
            self.reference[:, None], values.shape[1], axis=1)
        restored[self.coordinates] = values
        return restored

    def __call__(self, time_h):
        return self.restore(self.solution(time_h))


def _integrate_interval(kernel, start, end, y, control, method, rtol, atol, max_step_h):
    coordinates = _reachable_coordinates(kernel, y, control)
    restore = _DenseState(None, y.copy(), coordinates).restore
    tolerance = atol[coordinates] if isinstance(atol, np.ndarray) else atol
    solution = solve_ivp(lambda t, reduced: kernel.evaluate(t, restore(reduced), control)[0][coordinates],
                         (start, end), y[coordinates], method=method, rtol=rtol, atol=tolerance,
                         max_step=max_step_h, dense_output=True)
    if not solution.success:
        raise RuntimeError(f"shared integration failed on [{start}, {end}]: {solution.message}")
    return _DenseState(solution.sol, y.copy(), coordinates), solution.nfev


def _tolerances(kernel, initial, atol):
    if isinstance(atol, Mapping):
        if set(atol) != set(kernel.names):
            raise ValueError("atol mapping must name every integrated state exactly once")
        return np.array([finite(atol[name], name, positive=True) for name in kernel.names])
    if atol is not None:
        return finite(atol, "atol", positive=True)
    scale = max(initial.volume_l, initial.biomass_gdw, 1e-8)
    return np.array([1e-17 * scale if name.startswith("gene.") else
                     1e-10 if variable.compartment in (Compartment.REGULATORY, Compartment.ER) else
                     1e-12 * scale for name, variable in kernel.variables.items()])


def _regulatory_roundoff(kernel, values, atol, corrections):
    fraction_names = ("unfolded_fraction", *(f"signal.{name}" for name in DYNAMIC_SIGNALS))
    for name in (*fraction_names, *(f"upr.{name}" for name in UPR_VARIABLES)):
        i = kernel.index[name]
        upper = 1.0 if name in fraction_names else np.inf
        current = values[i]
        canonical = np.clip(current, 0.0, upper)
        error = float(np.max(np.abs(canonical - current)))
        tolerance = atol[i] if isinstance(atol, np.ndarray) else atol
        if error > 10.0 * tolerance:
            raise RuntimeError(f"regulatory state {name} exceeds its physical bounds by {error:g}; tighten tolerances")
        if error:
            if kernel.carbon[i] or kernel.nitrogen[i] or i in kernel.physical:
                raise RuntimeError("refusing a mass-changing regulatory normalization")
            corrections[name] = max(corrections.get(name, 0.0), error)
            values[i] = canonical


def _preflight(protocol, initial, parameters, observation):
    if initial.time_h != protocol.times_h[0]:
        raise ValueError("initial_state.time_h must equal the protocol's first absolute time")
    for control in protocol.controls:
        if not 20.0 <= control.temperature_c <= 40.0:
            raise ScientificRefusal("temperature is outside the declared prior operating window [20, 40] degC")
        if not 3.0 <= control.ph <= 8.0:
            raise ScientificRefusal("controlled medium pH is outside the declared prior window [3, 8]")
    if parameters.thermodynamics is not None and parameters.thermodynamics.require_complete:
        raise ScientificRefusal("complete thermodynamic feasibility cannot be certified: only LCY and adenylate kinase have supported molecular mappings; central-carbon lumping, transport work, phosphate/CoA, activity calibration and energy uncertainty remain unresolved")
    receipts = {e.name: e for e in initial.event_receipts}
    for event in protocol.events:
        if event.name in receipts and event != receipts[event.name]:
            raise ValueError(f"event {event.name!r} differs from its recorded receipt")
    volume = initial.volume_l
    boundaries = protocol.boundaries_h
    for i, time_h in enumerate(boundaries):
        for event in protocol.events:
            if event.time_h == time_h and event.name not in receipts:
                if event.withdraw_l >= volume:
                    raise ValueError(f"event {event.name!r} empties the reactor")
                volume += event.add_volume_l - event.withdraw_l
        if i < len(boundaries) - 1:
            control = protocol.control_at(time_h)
            volume += (control.feed_l_h - control.outflow_l_h) * (boundaries[i + 1] - time_h)
            if volume <= 0:
                raise ValueError("the protocol empties the reactor")
    if observation is not None and not isinstance(observation, ObservationModel):
        raise ValueError("observation must be an ObservationModel")


def _parameter_provenance(kernel, observation):
    out = dict(kernel.parameters.values)
    out.update(kernel.signalling.hog.provenance)
    for prefix, params in (
        ("upr_source", (upr.HAC1_MRNA_DECAY_PER_H, upr.HAC1_SPLICING_PER_H,
                        upr.HAC1_PROTEIN_DECAY_PER_H, upr.KAR2_PROTEIN_DECAY_PER_H,
                        upr.KAR2_MRNA_DECAY_PER_H, upr.HAC1_ISOFORM_DECAY_RATIO,
                        upr.HAC1_MAX_OCCUPANCY)),
        ("ph_source", (ph.PKA_ACETIC, ph.P_AH_ACETIC, ph.SURFACE_TO_VOLUME_PER_CM,
                       ph.ORIJ_INTERCEPT, ph.ORIJ_SLOPE)),
        ("peroxide_source", (oxidative.K_REF_PER_H, oxidative.X_REF_CELLS_PER_ML)),
    ):
        out.update({f"{prefix}.{param.name}": param for param in params})
    for gene in kernel.genotype.genes:
        out.update({f"gene.{gene.name}.{i}.{param.name}": param for i, param in enumerate(gene.parameters)})
    if kernel.parameters.acid_export is not None:
        export = kernel.parameters.acid_export
        out.update({f"acid_export.{name}": getattr(export, name)
                    for name in ("vmax", "km_anion", "atp_per_proton", "atp_per_anion")})
    if kernel.parameters.thermodynamics is not None:
        out.update({f"thermodynamic_activity.{name}": param
                    for name, param in kernel.parameters.thermodynamics.activity_coefficients.items()})
    if observation is not None:
        out.update({f"observation.{name}": param for name, param in observation.values.items()})
        out["observation.citrine_pka"] = Param.measured(
            "citrine_pka", CITRINE_PKA, "dimensionless",
            "photophysics.py::CITRINE_PKA; Griesbeck et al. 2001 PMID 11387331")
    return out


def _thermodynamic_audit(kernel, times, states, protocol):
    request = kernel.parameters.thermodynamics
    if request is None:
        return (), ()
    reports = tuple(kernel.thermodynamic_report(name, y, protocol.control_at(float(time)))
                    for time, y in zip(times, states.T, strict=True) for name in request.reactions)
    return reports, request.reactions


def simulate(protocol: Protocol, genotype: Genotype, initial_state: PhysicalState,
             parameters: EngineParameters, *, observation: ObservationModel | None = None,
             rtol: float = 1e-6, atol: float | Mapping[str, float] | None = None,
             max_step_h: float = 0.1, method: str = PINNED_METHOD) -> SimulationResult:
    if not isinstance(protocol, Protocol) or not isinstance(genotype, Genotype):
        raise ValueError("typed Protocol and Genotype are required")
    if not isinstance(initial_state, PhysicalState) or not isinstance(parameters, EngineParameters):
        raise ValueError("typed PhysicalState and EngineParameters are required")
    finite(rtol, "rtol", positive=True)
    finite(max_step_h, "max_step_h", positive=True)
    if method not in ("BDF", "Radau"):
        raise ValueError("the shared stiff clock supports BDF or Radau, without silent fallback")
    _preflight(protocol, initial_state, parameters, observation)
    kernel = _Kernel(parameters, genotype)
    y0 = kernel.pack(initial_state)
    y = y0.copy()
    tolerances = _tolerances(kernel, initial_state, atol)
    receipts = list(initial_state.event_receipts)
    receipt_names = {e.name for e in receipts}
    boundaries = protocol.boundaries_h
    boundary_states, event_factors, segments = {}, {}, []
    regulatory_corrections = {}
    nfev = 0
    for i, start in enumerate(boundaries):
        event_factors[start] = 1.0
        for event in protocol.events:
            if event.time_h == start and event.name not in receipt_names:
                y, factor = kernel.apply_event(y, event)
                event_factors[start] *= factor
                receipts.append(event)
                receipt_names.add(event.name)
        boundary_states[start] = y.copy()
        if i == len(boundaries) - 1:
            break
        end = boundaries[i + 1]
        control = protocol.control_at(start)
        dense, evaluations = _integrate_interval(kernel, start, end, y, control, method,
                                                 rtol, tolerances, max_step_h)
        nfev += evaluations
        segments.append((start, end, dense, control))
        y = dense(end)
        _regulatory_roundoff(kernel, y, tolerances, regulatory_corrections)
    times = np.asarray(protocol.times_h)
    states = np.empty((len(kernel.names), len(times)))
    for j, time in enumerate(times):
        if time in boundary_states:
            states[:, j] = boundary_states[time]
        else:
            segment = next(s for s in segments if s[0] < time < s[1])
            states[:, j] = segment[2](time)
    if not np.isfinite(states).all():
        raise RuntimeError("nonfinite shared trajectory")
    _regulatory_roundoff(kernel, states, tolerances, regulatory_corrections)
    minimum = float(np.min(states[kernel.physical]))
    if minimum < 0:
        raise RuntimeError(f"negative physical inventory {minimum:g}; tighten tolerances; no mass-changing projection was applied")
    truth = {name: states[i] for name, i in kernel.index.items()}
    variables = dict(kernel.variables)
    instantaneous = [kernel.evaluate(float(time), states[:, j], protocol.control_at(float(time)))[1]
                     for j, time in enumerate(times)]
    for name in instantaneous[0]:
        truth[name] = np.array([row[name] for row in instantaneous])
        if name == "flux.growth":
            units, basis = "gDW/h", "formation rate"
        elif name.startswith("flux.") or name in ("oxygen_transfer_mmol_h", "ph_pump_atp_mmol_h", "anion_export_atp_mmol_h"):
            units, basis = "mmol/h", "reaction extent rate"
        elif name in ("growth_per_h", "death_per_h", "washout_per_h", "net_live_per_h"):
            units, basis = "1/h", "specific rate"
        elif name == "osmolarity_mOsm":
            units, basis = "mOsm/L", "ideal dilute osmotic particle approximation"
        else:
            units, basis = "dimensionless", "algebraic driver"
        variables[name] = Variable(units, Compartment.REGULATORY, basis)
    for name in EXTRACELLULAR:
        truth[f"medium_mM.{name}"] = truth[f"external.{name}"] / truth["volume_l"]
        variables[f"medium_mM.{name}"] = Variable("mM", Compartment.MEDIUM, "concentration")
    for name in INTRACELLULAR:
        truth[f"content.{name}"] = np.divide(truth[f"internal.{name}"], truth["biomass_gdw"],
                                            out=np.zeros_like(times), where=truth["biomass_gdw"] > 0)
        variables[f"content.{name}"] = Variable("mmol/gDW", Compartment.CELL_WATER, "content per structural biomass")
    neutral_fraction = np.array([ph.undissociated_fraction(value, kernel.acid_pka) for value in truth["ph_c"]])
    truth["acid.intracellular_neutral_mmol"] = truth["internal.acetate"] * neutral_fraction
    truth["acid.intracellular_anion_mmol"] = truth["internal.acetate"] * (1.0 - neutral_fraction)
    truth["acid.buffered_protons_mmol"] = (truth["biomass_gdw"] * kernel.p["cell_water_l_gdw"]
                                           * kernel.p["buffer_capacity"] * (kernel.p["resting_ph"] - truth["ph_c"]))
    for name in ("acid.intracellular_neutral_mmol", "acid.intracellular_anion_mmol", "acid.buffered_protons_mmol"):
        variables[name] = Variable("mmol", Compartment.CELL_WATER, "algebraic partition, not an additional inventory")
    variables["ph_c"] = Variable("dimensionless", Compartment.CELL_WATER, "fast acid-buffer equilibrium")
    partition = PoolPartition(
        integrated=kernel.names,
        fast_algebraic={
            "ph_c": "mech/ph.py::charge_imbalance_mM, acid-buffer equilibrium; independent electrical charge dynamics are not claimed",
            "acid.intracellular_neutral_mmol": "mech/ph.py::undissociated_fraction, partition of the integrated acetate pool",
            "acid.intracellular_anion_mmol": "complement of the same neutral partition, never counted as additional carbon",
            "acid.buffered_protons_mmol": "acid-derived buffer proton equivalents, not the complete metabolic proton balance",
            "driver.upre": "UprChain Hac1 DNA occupancy; Hac1 and chaperone pools remain dynamic",
        },
        retained_fast_pool_candidates=tuple(f"internal.{name}" for name in
                                           ("atp", "adp", "amp", "camp", "nad", "nadh", "nadp", "nadph", "acetate")),
        observation_only_dynamic=("mature_reporter_mmol",) if observation is not None else (),
    )
    carbon, nitrogen = np.array([kernel.inventories(column) for column in states.T]).T
    c0, n0 = kernel.inventories(y0)
    for name, inventory, baseline in (("carbon", carbon, c0), ("nitrogen", nitrogen, n0)):
        consumed = truth["ledger.co2"] - y0[kernel.index["ledger.co2"]] if name == "carbon" else 0.0
        residual = (inventory + consumed - baseline
                    - (truth[f"ledger.{name}_in"] - y0[kernel.index[f"ledger.{name}_in"]])
                    + (truth[f"ledger.{name}_out"] - y0[kernel.index[f"ledger.{name}_out"]]))
        truth[f"balance.{name}_residual"] = residual
        variables[f"balance.{name}_residual"] = Variable("mmol", Compartment.LEDGER, "elemental balance residual")
    truth["balance.live_biomass_residual"] = (
        truth["biomass_gdw"] - initial_state.biomass_gdw
        - (truth["ledger.biomass_formed"] - initial_state.ledgers["biomass_formed"])
        + (truth["ledger.biomass_died"] - initial_state.ledgers["biomass_died"])
        + (truth["ledger.biomass_washed_out"] - initial_state.ledgers["biomass_washed_out"]))
    variables["balance.live_biomass_residual"] = Variable("gDW", Compartment.LEDGER, "birth-death-washout residual")
    reports, covered = _thermodynamic_audit(kernel, times, states, protocol)
    thermodynamics = parameters.thermodynamics
    constrained = thermodynamics.reactions if thermodynamics is not None and thermodynamics.mode == "constrain_supported" else ()
    ph_model = "paired_acid_export_prior" if parameters.acid_export is not None else "no_pump_bound"
    validity = Validity(
        parameter_basis="prior-conditional; not experimentally validated",
        assumptions=(
            "One well-mixed reactor clock in hours; event samples are right-continuous, withdrawal precedes addition, simultaneous events follow declared order.",
            "Biomass is structural/catalytic gDW excluding explicitly inventoried small molecules and recombinant expression; cell-water volume and cell count per gDW are explicit priors.",
            "Carbon and nitrogen are conserved through all reactions, death, flow, and events. ATP and pyridine nucleotides are explicit coupled pools, not external capacity multipliers.",
            "External mM is nominal broth-equivalent concentration; water displacement, crowding, and acid/base titrant inventories are not modeled. No physical amount is projected; sub-tolerance native regulatory roundoff is reported separately.",
            "Central-carbon moieties, amino-acid/nucleotide compositions, pentose recycling, native isoprenoid sink, and temperature response are reduced stoichiometric/kinetic priors, not a genome-scale metabolic or atom-tracing model.",
            "HOG reuses only source v16f/v16r/vVos on the shared state and clock; no source OD fit, source glucose balance, native protein molarity, or absolute native glycerol flux is transferred.",
            "UPR reference synthesis remains fixed while actual growth changes; native normalized UPR and unfolded-proteome pools are allocations within structural biomass, with an explicit prior ATP bill.",
            "pH uses the existing acid-buffer relation; medium pH is controlled, and cell-water changes scale acid and buffer concentrations together. With explicit AcidExportParameters, acetate anion and one acid-derived proton are co-exported and the matrix pays their declared ATP costs once; otherwise this is the no-pump bound.",
            "Yap1 uses the existing external response; conversion of imported and respiratory peroxide to a whole-cell oxidant-equivalent pool is a prior, not an identified cytosolic gradient.",
            "Fast acid protonation/buffering and UPRE binding are algebraic partitions, not duplicate inventories. ATP/ADP/AMP/cAMP, pyridine nucleotides and total acetate remain dynamic in the shared stiff solve; no kinetic pool is silently set to steady state.",
            "Thermodynamic mode is explicit: audit leaves rates unchanged; constrain_supported applies continuous direction/affinity factors to the named molecular steps using caller-declared activities, refusing missing or refuted energies. Neither mode certifies the uncovered network or biological validity.",
        ),
        unsupported=(
            "Absolute DTT/tunicamycin entry: supply a declared folding_inhibition intervention, not an inferred drug dose.",
            "Other carbon/nitrogen species, weak acids other than acetate, dynamic medium titration, independently regulated Pma1/anion export, membrane voltage/leak, osmolyte ion-specific toxicity, and thermal excursions outside [20, 40] degC.",
            "Plasmid segregation, copy-number control, heterogeneous populations, secretion, whole-proteome sequence composition, and products other than the CrtE/CrtI/CrtYB branch.",
            "Complete thermodynamic feasibility, phosphate/CoA inventories, charge/electron/oxygen-atom conservation, net metabolic/ATP protonation balance, membrane electrical work, and calibrated genome-scale fluxes. The acid export pairing closes only the transport charges, not the complete chemical charge balance.",
            "Absolute optical predictions without calibration; reporter fluorescence after cell death and detector-dependent missingness are not identified.",
            "The four transcribed stress axes in signalling.TRANSCRIBED_AXES -- " + ", ".join(sorted(TRANSCRIBED_AXES)) + " -- are REGISTERED AND UNDRIVEN: no state of theirs is integrated and no rate law or promoter reads one, so nothing in this result is evidence about calcium, cAMP/PKA input, published glycolytic ATP or Slt2/CWI. Each registration names the refusal that would have to be retired first. cell_wall_slt2 is undriven for a different kind of reason than the other three: its input and its output consumer both already exist, and integrating it would represent this engine's osmotic response twice.",
        ), thermodynamic_covered=covered,
        thermodynamic_uncovered=tuple(name for name in kernel.reaction_names if name not in covered),
    )
    observed = _observe(kernel, times, states, segments, event_factors, observation, rtol, max_step_h) if observation is not None else None
    matrix = kernel.matrix.copy()
    matrix.setflags(write=False)
    diagnostics = {
        "method": method, "n_rhs_evaluations": nfev, "rtol": rtol,
        "event_boundaries_h": boundaries, "events_applied": tuple(e.name for e in receipts),
        "protocol_source": protocol.source, "genotype_source": genotype.source,
        "pool_partition": partition, "ph_model": ph_model,
        "acid_export_source": parameters.acid_export.source if parameters.acid_export else "no active export model supplied",
        "thermodynamic_constrained_reactions": constrained,
        "thermodynamic_mode": thermodynamics.mode if thermodynamics is not None else "not_requested",
        "thermodynamic_rate_law": ("existing signed kinetic rate times max(0, 1-exp(direction*dG/RT)); direction/affinity constraint only, not identification of reverse kinetics; exact zero-product limits use no pseudocount"
                                   if constrained else "not applied; audit-only or no thermodynamics never changes a flux"),
        "scientific_blockers": frozen_mapping({
            "independent_pma1": "Requires Pma1 abundance/activity, pH-dependent regulation, net H+/ATP coupling in the actual carbon regime, membrane voltage and proton-leak measurements; paired acid export is not this model.",
            "anion_export_identity": "Requires transporter/substrate identity, Vmax, Km and ATP/ion coupling; the supplied export capacity, affinity and cost are declared priors, not identified Tpo2/Tpo3/Pdr12 kinetics.",
            "pump_thermodynamics": "Requires membrane electrical work and state-matched ATP/ADP/phosphate/proton activities and full pump/counter-transport stoichiometry; no pump thermodynamic certificate is supplied.",
            "complete_element_charge_redox_balance": "Requires molecular formulas, phosphate/CoA/water/proton inventories, speciation and compartmental electron carriers beyond the current carbon/nitrogen moieties and tracked cofactors.",
            "metabolic_fast_pool_elimination": "Requires an error-controlled conservation-preserving reduction with pool-specific kinetics and driver timescales; fast ATP and redox candidates therefore remain dynamic, not silently set to steady state.",
        }),
        "minimum_physical_inventory": minimum, "amount_projection": "none",
        "state_projection": "native regulatory roundoff only" if regulatory_corrections else "none",
        "regulatory_roundoff_corrections": frozen_mapping(regulatory_corrections),
        "structurally_zero_by_interval": tuple((start, end, tuple(name for i, name in enumerate(kernel.names)
                                                                 if i not in dense.coordinates))
                                                for start, end, dense, _ in segments),
        "zero_elimination": "stoichiometric and catalyst reachability, not a timescale approximation",
        "stoichiometric_state_names": kernel.names, "reaction_names": kernel.reaction_names,
        "native_hog": kernel.signalling.hog.metadata,
        "stoichiometric_matrix": matrix, "thermodynamic_reports": reports,
        "thermodynamic_scope": "covered names were submitted to the gate; constraints apply only to thermodynamic_constrained_reactions. Audit-only and per-sample cannot_say never certify flux feasibility. Full-network certification is refused.",
        "thermodynamic_report_coordinates": tuple((float(time), name) for time in times for name in covered),
        "thermodynamic_verdict_counts": frozen_mapping({verdict: sum(r.feasibility == verdict for r in reports)
                                                        for verdict in Feasibility.ALL}),
        "thermodynamic_source": parameters.thermodynamics.source if parameters.thermodynamics else "not supplied",
        "thermodynamically_blocked_sample_indices": tuple(sorted({i // len(covered) for i, report in enumerate(reports)
                                                                  if report.feasibility == Feasibility.CANNOT_RUN})),
    }
    return SimulationResult(times, truth, variables, kernel.unpack(float(times[-1]), states[:, -1], receipts),
                            validity, _parameter_provenance(kernel, observation), diagnostics, observed)


def _observe(kernel, times, states, segments, event_factors, observation, rtol, max_step_h):
    genes = [gene for gene in kernel.genotype.genes if gene.name == observation.reporter]
    if len(genes) != 1 or genes[0].role != "reporter":
        raise ValueError("observation.reporter must name a reporter gene in this genotype")
    gene = genes[0]
    idx = kernel.index
    protein_index = idx[f"gene.{gene.name}.protein"]
    active_index = idx[f"gene.{gene.name}.active"]
    p = {name: float(param) for name, param in observation.values.items()}
    mature = observation.initial_mature_mmol
    start_factor = event_factors[float(times[0])]
    if mature * start_factor > states[protein_index, 0]:
        raise ValueError("initial mature reporter cannot exceed the initial total reporter")
    mature *= start_factor
    boundary_mature = {float(times[0]): mature}
    trajectories = []
    for start, end, dense, control in segments:
        def rhs(time, value):
            y = dense(time)
            e, _, _, _, _, _, _, _, death, _ = kernel.physiology(y, control)
            kmat = p["maturation"] * saturation(e["oxygen"], p["maturation_oxygen_k"])
            loss = float(gene.protein_decay_per_h) + death + control.outflow_l_h / y[idx["volume_l"]]
            return [kmat * (max(y[protein_index], 0.0) - value[0]) - loss * value[0]]
        solution = solve_ivp(rhs, (start, end), [mature], method=PINNED_METHOD,
                             rtol=rtol, atol=1e-18 * max(1.0, states[protein_index, 0]),
                             max_step=max_step_h, dense_output=True)
        if not solution.success:
            raise RuntimeError(f"observation maturation failed: {solution.message}")
        trajectories.append((start, end, solution.sol))
        mature = float(solution.y[0, -1]) * event_factors[end]
        boundary_mature[end] = mature
    maturation = np.array([boundary_mature[float(t)] if float(t) in boundary_mature else
                           next(dense for start, end, dense in trajectories if start < t < end)(t)[0]
                           for t in times])
    od_calibration = ODCalibration(p["od_saturation"], 1.0, p["od_blank"], p["gDW_per_OD_l"])
    optics = ReporterOptics(p["rfu_gain"], p["rfu_blank"], p["autofluorescence"], p["inner_filter"])
    od, rfu = [], []
    for j in range(len(times)):
        y = states[:, j]
        volume, biomass = y[idx["volume_l"]], y[idx["biomass_gdw"]]
        live_density = biomass / volume
        od.append(observe_od(live_density, y[idx["dead_biomass_gdw"]] / volume, od_calibration))
        protein = y[protein_index]
        folded_fraction = y[active_index] / protein if protein > 0 else 0.0
        reporter_content = maturation[j] / biomass if biomass > 0 else 0.0
        pigment_content = y[idx["internal.beta_carotene"]] / biomass if biomass > 0 else 0.0
        quench = citrine_ph_response(kernel.cytosolic_ph(y))
        rfu.append(observe_rfu(reporter_content * folded_fraction * quench, live_density,
                               pigment_content, optics))
    rng = np.random.default_rng(observation.seed)
    values, missing = {}, {}
    for channel, expected in (("od", od), ("rfu", rfu)):
        noisy = np.asarray(expected) + rng.normal(0.0, p[f"{channel}_noise_sd"], len(times))
        mask = rng.random(len(times)) < p["missing_probability"]
        values[channel] = np.where(mask, np.nan, noisy)
        missing[channel] = mask
    maturation.setflags(write=False)
    return Observations(values, missing, {
        "reporter": gene.name, "seed": observation.seed, "units": MappingProxyType({"od": "OD", "rfu": "RFU"}),
        "mature_mmol": maturation, "final_mature_mmol": mature,
        "truth_feedback": False, "calibration_basis": tuple(sorted({param.tag for param in observation.values.values()})),
        "missingness": "independent per-channel Bernoulli MCAR; no imputation",
        "death": "live reporter only; fluorescence from dead cells is unsupported",
        "measurement_basis": "concentration at a fixed calibrated optical path, not total-well inventory; path-length changes require a separate calibration",
        "maturation": "one-way integration over the same absolute protocol intervals, not the read grid; oxygen-gated prior rate",
    })
