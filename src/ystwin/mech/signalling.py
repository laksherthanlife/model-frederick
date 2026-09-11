"""The regulatory states, and where each one is allowed to go.

Six dynamic signals are integrated here. Four reach a promoter in the shipped genotype:
``hog1``, ``hsf1``, ``yap1`` and ``torc1`` -- the last through a coefficient that is a
declared PRIOR, not a sourced one. The carbon pair ``pka`` and ``snf1`` reach no promoter
at all, and that is a refusal with two named missing measurements
(:data:`contracts.PROMOTER_COUPLING_NOT_BUILT`), not an oversight.

They are not inert, which is the part a "drives no promoter" reading misses:
:data:`CARBON_NITROGEN_FLUX_COUPLING` names the rate laws that read them, so the carbon
programme is coupled to metabolism, and only to metabolism.

Four further axes are TRANSCRIBED but not integrated. :data:`TRANSCRIBED_AXES` registers each
one's states, its environment input, and the single named refusal that keeps it undriven, and
:data:`TRANSCRIBED_AXIS_FLUX_COUPLING` applies the same rule to them that the check below
applies to the six: an axis that names no flux has to say why, or it is being reported as dead.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType

import numpy as np

from . import (
    atp_teusink2000, calcium_ke2013, carbon_williamson2009, cell_wall_talemi2016, oxidative, upr,
)
from .contracts import (
    PROMOTER_UNCOUPLED_SIGNALS, SIGNALS, AxisRegistration, Compartment, ScientificRefusal,
)
from .hog import HogModel
from .params import Param

CARBON_NITROGEN_FLUX_COUPLING = MappingProxyType({
    "pka": ("glucose_transport",),
    "snf1": ("ethanol_oxidation", "glycogen_mobilization"),
    "torc1": ("growth", "glycogen_synthesis"),
})
"""Every `mech/engine.py` rate law that reads a carbon or nitrogen state directly.

Declared where the states are made, and held true by ``tests/test_shared_engine.py``: each
flux has to move when its own signal's genotype gain moves, so this map cannot drift away
from the rate law without a test failing. ``pka`` and ``snf1`` appear here and in no
promoter, which is the whole content of :data:`contracts.PROMOTER_UNCOUPLED_SIGNALS`.
"""

if set(CARBON_NITROGEN_FLUX_COUPLING) - set(SIGNALS) or not set(
        PROMOTER_UNCOUPLED_SIGNALS) <= set(CARBON_NITROGEN_FLUX_COUPLING):
    raise ScientificRefusal("a signal declared uncoupled from transcription must still name "
                            "the fluxes that read it, or it is being reported as dead")


_CAMP_COMPARTMENTS = {"Glucout": Compartment.MEDIUM, "Glucin": Compartment.CELL_WATER,
                      "cAMP": Compartment.CELL_WATER}

TRANSCRIBED_AXES = MappingProxyType({axis.axis: axis for axis in (
    AxisRegistration(
        axis="calcium",
        module="mech/calcium_ke2013.py",
        source=("Ke R, Ingram PJ, Haynes K 2013, PLoS Comput Biol 9(1):e1002879, "
                "doi 10.1371/journal.pcbi.1002879, CC BY. Text S1 Eqns 3.3-3.4, already "
                "vendored and sha256-pinned as data/native_reference_models/ke2013 and "
                "already transcribed as mech/ph_ke2013.py::calcineurin_rhs"),
        states=calcium_ke2013.STATE_VARIABLES,
        environment_input=(
            "Control.ph. This is the one new axis whose driver the contract already carries: "
            "Eqn 3.3's alkaline arm is +kCa,pH*(pH_ext - 6.5), and Control.ph is a real "
            "control the preflight admits over [3, 8]. It is also linear and unbounded below, "
            "so below pH 6.3327 the equation has no non-negative rest point at all and "
            "CalciumAxis.check_ph raises rather than clipping -- and every protocol in this "
            "repository runs at pH 4.0-5.0, Control.ph defaulting to 5.0"),
        missing_medium=("MEDIUM_SPECIES_NOT_CARRIED::medium.calcium",
                        "MEDIUM_SPECIES_NOT_CARRIED::medium.sodium"),
        panel_promoter="CDRE",
        promoter_refusals=("AXIS_PROMOTER_NOT_BUILT::promoter_response.crz1_cdre",),
        readers=(),
        undriven_because=(
            "Nothing can read it. Crz1p is the axis's output and its only printed output arm "
            "has no transcribable coefficient, so integrating these four states would add a "
            "regulator that no flux and no promoter consumes -- the failure the check above "
            "refuses one level in. The driver is real but narrow: the axis is undefined below "
            "pH 6.3327 and unvalidated below 6.5, so it is silent across this repository's "
            "whole operating range, and the two panel stressors that would open it "
            "(calcium_chloride at weight 1.00, NaCl at 0.25) need medium species EXTRACELLULAR "
            "does not carry. Clamping pH upward to keep the axis alive would invent the one "
            "constant Ke does not print"),
    ),
    AxisRegistration(
        axis="carbon_camp_input",
        module="mech/carbon_williamson2009.py",
        source=("Williamson T, Schwartz J-M, Kell DB, Stateva L 2009, BMC Syst Biol 3:70, "
                "doi 10.1186/1752-0509-3-70, CC BY 2.0. The 'Complete cAMP Model', vendored "
                "as data/native_reference_models/williamson2009 with the SBML pinned at "
                "sha256 e608d6651dd3332c0a2f146bdd9402c77ee07f7c6358284bef911da5bcd1b782"),
        states=tuple((name, "mM", _CAMP_COMPARTMENTS.get(name, Compartment.REGULATORY),
                      "Williamson 2009 Additional file 1 listOfSpecies")
                     for name in carbon_williamson2009.DYNAMIC_SPECIES),
        environment_input=(
            "external.glucose (mM) -> Glucout and internal.glucose -> Glucin, both of which "
            "EXTRACELLULAR and INTRACELLULAR already carry. Nothing is missing on the input "
            "side; what is missing is on the output side and in the units"),
        missing_medium=(),
        panel_promoter="CSRE",
        promoter_refusals=("PROMOTER_COUPLING_NOT_BUILT::promoter_response.pka",
                           "PROMOTER_COUPLING_NOT_BUILT::promoter_response.snf1"),
        readers=(),
        undriven_because=(
            "This axis does NOT close the panel's carbon module and is registered under its own "
            "name so that it cannot be read as though it did: the panel declares carbon as "
            "Snf1 -> Adr1/Cat8 -> CSRE, and the strings Snf1, Adr1, Cat8, Mig1, CSRE and Msn2 "
            "occur zero times in either pinned Williamson artifact. What it transcribes is the "
            "cAMP/PKA INPUT side, which terminates in pka -- already integrated, already in "
            "PROMOTER_UNCOUPLED_SIGNALS, and refused at the promoter by the two rows cited "
            "above, which this axis does not retire. Adopting it as the engine's cAMP law is "
            "separately blocked: its mM-to-mmol/gDW bridge is refused twice over "
            "(williamson.cell_volume_litres and williamson.absolute_abundance_scale), and "
            "running it through cell_water_l_gdw would launder one assertion through another"),
    ),
    AxisRegistration(
        axis="atp_glycolysis",
        module="mech/atp_teusink2000.py",
        source=("Teusink et al. 2000, Eur J Biochem 267:5313-29, "
                "doi 10.1046/j.1432-1327.2000.01527.x, as BioModels BIOMD0000000064, CC0 1.0, "
                "vendored at sha256 "
                "d6a75df939d6707fdb904b44c1f530c8dc3600fd44267ff8e8f38dac560b4361"),
        states=tuple((name, "mM", Compartment.CELL_WATER,
                      "Teusink 2000 BIOMD0000000064 listOfSpecies; 'P' is the AXP phosphate "
                      "pool that ATP, ADP and AMP are read off by the deposit's own "
                      "assignment rules")
                     for name in atp_teusink2000.SPECIES_ORDER),
        environment_input=(
            "GLCo -> external.glucose (mM), which EXTRACELLULAR carries. It is the one panel "
            "stressor this axis answers (glucose_starvation, weight -0.60) and it works down "
            "to a validated floor of 2.0 mM, below which the source's own turbo-design "
            "instability stops it settling and steady_state raises rather than returning an "
            "unconverged number"),
        missing_medium=("MEDIUM_SPECIES_NOT_CARRIED::medium.respiratory_inhibitor",),
        panel_promoter="none",
        promoter_refusals=(),
        readers=(),
        undriven_because=(
            "This is a prior-replacement, not a new axis, and the replacement is blocked at the "
            "units: teusink.cytosol_l_per_gdw is REFUSED, so the deposit's mM cannot become the "
            "engine's mmol/gDW except through cell_water_l_gdw, which is an engine PRIOR and not "
            "a Teusink number. Two further mismatches would have to be settled first, and both "
            "are recorded rather than assumed: the deposit's KeqAK = 0.45 is written for "
            "2 ADP -> ATP + AMP and engine.py writes that reaction reversed, so wiring them "
            "together without inverting inverts the adenylate charge response; and the engine "
            "holds a fourth adenylate reservoir, cAMP, that the source has no counterpart for, "
            "so SUM_P cannot be mapped onto atp + adp + amp alone. The panel declares atp a "
            "metabolite pool with no promoter, so no transcriptional refusal is owed here"),
    ),
    AxisRegistration(
        axis="cell_wall_slt2",
        module="mech/cell_wall_talemi2016.py",
        source=("Talemi SR et al. 2016, Sci Rep 6:30950, doi 10.1038/srep30950, as BioModels "
                "MODEL1606100000, CC0 1.0, vendored at sha256 "
                "c4937f14c2d8e36b13a449b5ac9771d4e17426543298d7a2c9ca9c4df835e5dc"),
        states=cell_wall_talemi2016.STATE_VARIABLES,
        environment_input=(
            "external.osmolyte (mM), which EXTRACELLULAR already carries and which this engine "
            "ALREADY doses -- the default medium runs it at 250 mM. This is the first axis "
            "whose input and whose output consumer both exist, and measuring what that would "
            "cost is the whole content of undriven_because. The dose side is refused anyway at "
            "external_osmolarity_basal: Talemi's cen is an ADDITION on top of a declared "
            "ce0 = 260000 uM total medium osmolarity, and this repository declares no total for "
            "its media, so there is no number to add a dose to"),
        missing_medium=("MEDIUM_SPECIES_NOT_CARRIED::medium.cell_wall_damaging_agent",),
        panel_promoter="RLM1 box",
        promoter_refusals=("AXIS_PROMOTER_NOT_BUILT::promoter_response.rlm1_box",),
        readers=(),
        undriven_because=(
            "IT WOULD DOUBLE-COUNT THE OSMOTIC RESPONSE, and that was measured rather than "
            "argued. Five of the twelve states below are already integrated here. Three come "
            "from a DIFFERENT deposit -- Petelenz-Kurdziel 2013, run by NativeHogResponse: Vos "
            "is cell_volume_ratio and Hog1/Hog1PP are the two sides of signal.hog1. The other "
            "two are this engine's own stoichiometry: Glyin is internal.glycerol, made by "
            "glycerol_synthesis, and Glyex is external.glycerol, reached by glycerol_export -- "
            "which this engine already gates on (1 - hog1), the same Fps1 closure Talemi writes "
            "as v5/v6 over Fps1 and Fps1P, so the adaptation loop is duplicated in its feedback "
            "and not only in its states. Both blocks read the same medium "
            "osmolarity, so integrating both makes one osmotic dose drive two Hog1 pools and "
            "two cell volumes in one cell. That they are the same biology is not an assumption: "
            "both close the same loop -- glycerol accumulates, export shuts, volume recovers -- "
            "and run at the same +750 mM step from each model's own rest state they agree that "
            "Hog1 saturates, peaking at 0.890 of the pool (engine, 0.10 h) against 0.781 "
            "(Talemi, 0.05 h). That they cannot both be right is not an assumption either: they "
            "disagree on the adaptation and go on disagreeing. Talemi is adapted by 0.5 h and "
            "flat thereafter at Hog1PP 0.134 and volume 0.968; the engine is still at 0.760 and "
            "0.595 at 6 h, having recovered a fifth of its volume excursion where Talemi "
            "recovered 89% of a much smaller one in half an hour. That is 5.7x in Hog1 and 12.8x "
            "in volume excursion, permanently, on one dose. Integrating both would put two "
            "irreconcilable answers to one question in one state vector with nothing to "
            "arbitrate them, which is the degeneracy this repository has already paid for once. "
            "Wiring only the Slt2 arm from the engine's volume avoids the duplicated states and "
            "is refused separately and quantitatively at osmotically_active_cell_volume_bridge: "
            "the arm's rest point at 6 h reads 2.75% Slt2 phosphorylation off the engine's "
            "volume against 14.80% off Talemi's own, a 5.4x error that never relaxes, against "
            "the source's own published resting baseline of 24.6% -- so the axis would report a "
            "cell permanently at maximal CWI shutdown after one osmotic shock. And 0.595 is far "
            "outside the 0.90-1.10 band over which the source's sign was checked at all. "
            "Nothing is lost by leaving it registered: the axis "
            "answers neither panel stressor -- congo red and caffeine have no medium species and "
            "no source-side dose -- and supplies no promoter, so the engine's osmotic behaviour "
            "is strictly better off with one representation than with two"),
    ),
)})
"""Every axis transcribed against a published source, and its exact status in this engine.

All four are REGISTERED AND UNDRIVEN, each for a different reason, and the differences are
the content: calcium has its driver and no consumer, carbon has both ends and no unit bridge,
atp has its dose and a refused unit bridge, and cell_wall_slt2 -- the only one whose driver and
whose output consumer BOTH already exist -- is refused because integrating it would represent
this engine's osmotic response twice, which its registration measures rather than asserts.
Nothing here is a placeholder for work that was skipped: each entry cites the refusal that
would have to be retired first, and :func:`~mech.contracts.refused_row` resolves every one of
those at import, so a refusal that is renamed or quietly promoted fails here.
"""

TRANSCRIBED_AXIS_FLUX_COUPLING = MappingProxyType({
    "calcium": (),
    "carbon_camp_input": (),
    "atp_glycolysis": (),
    "cell_wall_slt2": (),
})
"""The same declaration :data:`CARBON_NITROGEN_FLUX_COUPLING` makes, for the new axes.

Empty is a claim, not an absence: it says no rate law in `mech/engine.py` reads these states,
which is checked against the registrations below and against the real reaction names in
`_Kernel.__init__`. An axis wired later has to appear in both places or neither.
"""

if set(TRANSCRIBED_AXIS_FLUX_COUPLING) != set(TRANSCRIBED_AXES) or any(
        TRANSCRIBED_AXIS_FLUX_COUPLING[name] != axis.readers
        for name, axis in TRANSCRIBED_AXES.items()):
    raise ScientificRefusal("the same rule one level out: a transcribed axis names the fluxes "
                            "that read it here and in its own registration, and the two must "
                            "agree -- an axis that names none carries its reason instead")


def saturation(amount: float, half: float) -> float:
    positive = max(float(amount), 0.0)
    return positive / (half + positive)


@lru_cache(maxsize=1)
def source_hog() -> HogModel:
    try:
        return HogModel.from_source()
    except (OSError, ValueError) as exc:
        raise ScientificRefusal("the checksum-audited HOG source is required; no substitute HOG curve is installed") from exc


class NativeHogResponse:
    def __init__(self, source: HogModel | None = None):
        self.source = source_hog() if source is None else source
        if not isinstance(self.source, HogModel):
            raise ValueError("source must be a HogModel")
        self.model = self.source.kinetic_model
        self.initial = self.model.initial_state()
        self.indices = {name: i for i, name in enumerate(self.model.species_ids)}
        self.reactions = {name: self.model.reaction_ids.index(name)
                          for name in ("v16f", "v16r", "vVos")}
        self.reference = self.model.evaluate(0.0, self.initial)
        self.total = self.reference["Hog1"] + self.reference["Hog1PP"]
        metadata = self.model.metadata
        source = ("Imported source definition, not an independent measurement: "
                  "Petelenz-Kurdziel et al. 2013 doi:10.1371/journal.pcbi.1003084; "
                  f"audited SBML sha256 {metadata['sha256']}")
        provenance = {
            f"hog_source.{name}": Param.asserted(f"hog_source.{name}", value, "source-native", source)
            for name, value in self.model.parameter_values.items()
            if name in ("kv16f_1", "kv16f_2", "kv16f_3", "kv16r_1", "vV_1", "vV_2", "vV_R", "vV_T")
        }
        provenance.update({
            f"hog_reference.{name}": Param.asserted(f"hog_reference.{name}", self.reference[name],
                                                    "source-native", source)
            for name in ("intra", "cellvol", "cin", "Hog1", "Hog1PP", "OsmoE")
        })
        self.provenance = MappingProxyType(provenance)
        self.metadata = MappingProxyType({
            "model_id": metadata["model_id"], "sha256": metadata["sha256"],
            "selected_reaction_formulas": MappingProxyType({name: metadata["reaction_formulas"][name]
                                                            for name in self.reactions}),
            "selected_assignment_formulas": MappingProxyType({name: metadata["assignment_formulas"][name]
                                                              for name in ("Vm", "CellSurface", "Turgor", "OsmoE")}),
            "mapping": "shared instantaneous osmolarity replaces the source step; volume is a ratio and Hog1 is a fraction, so no native protein molarity or native absolute glycerol flux is assumed",
        })

    def rates(self, time_h: float, *, fraction: float, volume_ratio: float,
              osmolarity_M: float, glycerol_external_M: float, glycerol_internal_M: float,
              temperature_c: float) -> tuple[float, float, float]:
        ratio = max(float(volume_ratio), 1e-8)
        fraction = float(np.clip(fraction, 0.0, 1.0))
        y = self.initial.copy()
        updates = {
            "cellvol": self.reference["cellvol"] * ratio,
            "Hog1": self.total * (1.0 - fraction) / ratio,
            "Hog1PP": self.total * fraction / ratio,
            "cin": self.reference["cin"] / ratio,
            "glycerol_i": max(glycerol_internal_M, 0.0),
            "glycerol_e": max(glycerol_external_M, 0.0),
        }
        for name, value in updates.items():
            y[self.indices[name]] = value
        time_s = 3600.0 * time_h
        parameters = {
            "t_stress": time_s - 5.0,
            "parameter_97": osmolarity_M - self.reference["OsmoE"],
            "vV_T": temperature_c + 273.15,
        }
        rates = self.model.reaction_rates(time_s, y, parameters)
        factor = 3600.0 * ratio / (self.reference["intra"] * self.total)
        activation = factor * rates[self.reactions["v16f"]]
        deactivation = factor * rates[self.reactions["v16r"]]
        volume_rate = (3600.0 * rates[self.reactions["vVos"]]
                       / (self.reference["intra"] * self.reference["cellvol"]))
        return float(activation), float(deactivation), float(volume_rate)


@dataclass(frozen=True)
class SignalRates:
    derivatives: dict[str, float]
    upr_derivatives: np.ndarray
    drivers: dict[str, float]
    volume_rate: float
    unfolded_rate: float
    native_repair_atp: float
    hog_phosphorylation: float


class Signalling:
    def __init__(self, values, genotype, *, hog: HogModel | None = None):
        self.p = values
        self.genotype = genotype
        self.hog = NativeHogResponse(hog)
        self.upr = upr.UprChain(values["upr_basal_load"], values["upr_hill"],
                                values["upr_basal_share"], values["upr_folding"],
                                values["reference_growth"])
        self.upr_basal = self.upr.basal()
        self.upr_indices = {name: i for i, name in enumerate(upr.UPR_STATES.names)}

    def evaluate(self, time_h, *, signals, upr_state, cell_volume_ratio, unfolded,
                 extracellular, intracellular_content, hsp70_content, expression_content,
                 growth, energy, redox, ph_c, control) -> SignalRates:
        p = self.p
        gain = self.genotype.signal_activity
        thermal_load = max(control.temperature_c - p["temperature_reference"], 0.0) / p["heat_scale"]
        hsp = saturation(hsp70_content, p["hsp70_fold_k"])
        stress_ros = saturation(intracellular_content["ros"], p["ros_k"])
        unfolded = float(np.clip(unfolded, 0.0, 1.0))
        unfold = (p["native_unfolding"] + p["heat_unfolding"] * thermal_load ** 2
                  + p["ros_unfolding"] * stress_ros) * (1.0 - unfolded)
        repair = p["folding"] * energy * hsp * unfolded
        free_hsp = hsp70_content / (1.0 + unfolded / p["hsf_client_k"])
        targets = {
            "pka": saturation(intracellular_content["camp"], p["camp_k"]),
            "snf1": 1.0 - saturation(extracellular["glucose"], p["glucose_k"]) * energy,
            "torc1": saturation(extracellular["nitrogen"], p["nitrogen_k"]) * energy * (1.0 - unfolded),
            "yap1": float(oxidative.yap1_nuclear_fraction(
                extracellular["peroxide"] + p["yap1_k"] * intracellular_content["ros"] / p["ros_k"],
                p["yap1_k"])),
            "hsf1": p["hsf_chaperone_k"] / (p["hsf_chaperone_k"] + free_hsp),
        }
        derivatives = {name: (gain.get(name, 1.0) * target - signals[name]) / p[f"{name}_tau"]
                       for name, target in targets.items()}
        cell_water = p["cell_water_l_gdw"] * max(cell_volume_ratio, 1e-8)
        osmolarity = sum(value for name, value in extracellular.items()
                         if name not in ("oxygen", "glycerol")) / 1000.0
        on, off, volume_rate = self.hog.rates(
            time_h, fraction=signals["hog1"], volume_ratio=cell_volume_ratio,
            osmolarity_M=osmolarity, glycerol_external_M=extracellular["glycerol"] / 1000.0,
            glycerol_internal_M=intracellular_content["glycerol"] / cell_water / 1000.0,
            temperature_c=control.temperature_c,
        )
        phosphorylation = gain.get("hog1", 1.0) * on * energy
        derivatives["hog1"] = phosphorylation - off
        native_synthesis = energy * redox * saturation(extracellular["nitrogen"], p["nitrogen_k"])
        influx = (p["upr_basal_load"] + p["er_expression_load"] * expression_content) * native_synthesis
        upr_derivatives = np.array(self.upr.rhs(
            influx, control.folding_inhibition, growth_rate_per_h=growth,
            synthesis_scale=native_synthesis, folding_scale=energy,
        )(time_h, upr_state))
        hac1 = max(upr_state[self.upr_indices["hac1_protein"]], 0.0)
        drivers = {name: float(np.clip(value, 0.0, 1.0)) for name, value in signals.items()}
        drivers["upre"] = gain.get("upre", 1.0) * hac1 / (1.0 + hac1)
        drivers["acid"] = gain.get("acid", 1.0) * float(np.clip(
            (p["resting_ph"] - ph_c) / p["acid_signal_scale"], 0.0, 1.0))
        native_cost = p["folding_atp"] * repair + p["upr_atp"] * native_synthesis * (
            max(upr_state[self.upr_indices["kar2_mrna"]], 0.0)
            + max(upr_state[self.upr_indices["chaperone"]], 0.0))
        return SignalRates(derivatives, upr_derivatives, drivers, volume_rate,
                           unfold - repair, native_cost, phosphorylation)
