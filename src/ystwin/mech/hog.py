from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np

from .. import paths
from .kinetic_sbml import KineticModel, Trajectory

__all__ = ["HogModel", "HogProtocol", "HogGlycerolBalance", "hog_glycerol_balance", "hog_metabolic_signals"]

_SOURCE_MODEL_ID = "PetelenzKuehn_osmoadaptation_WT"
_SOURCE_SHA256 = "7da99398c53185503022a4dc2ec54b7666280f4d461cf96a7ccc762d56963c19"
_SOURCE_NACL_MOLAR = 0.4
_SOURCE_MIXING_SECONDS = 5.0
_PROTOCOL_PARAMETERS = frozenset({"parameter_97", "t_stress"})


@dataclass(frozen=True)
class HogProtocol:
    nacl_molar: float
    shock_time_s: float = 3600.0

    def __post_init__(self):
        for name, value in asdict(self).items():
            if isinstance(value, bool) or not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be a finite nonnegative physical input")


@dataclass(frozen=True)
class HogModel:
    kinetic_model: KineticModel
    source_metadata: dict

    @classmethod
    def from_source(cls, directory=None):
        root = paths.data_dir() / "hog2013" if directory is None else Path(directory)
        manifest = json.loads((root / "sources.json").read_text())
        record = manifest["assets"]["model_wt"]
        if (manifest["publication"]["doi"] != "10.1371/journal.pcbi.1003084"
                or record["filename"] != "model_wt.xml" or record["sha256"] != _SOURCE_SHA256):
            raise ValueError("HOG source manifest does not identify the audited wild-type model")
        path = root / "model_wt.xml"
        if hashlib.sha256(path.read_bytes()).hexdigest() != _SOURCE_SHA256:
            raise ValueError("HOG model checksum does not match the primary source")
        model = KineticModel.from_sbml(path)
        if model.model_id != _SOURCE_MODEL_ID:
            raise ValueError("HOG model identity is not the published wild type")
        return cls(model, manifest)

    @property
    def parameter_values(self):
        return self.kinetic_model.parameter_values

    def simulate(self, times_s, *, protocol: HogProtocol, parameters=None,
                 initial_state=None, rtol=1e-7, atol=1e-10) -> Trajectory:
        if not isinstance(protocol, HogProtocol):
            raise ValueError("an explicit HogProtocol is required")
        overrides = {} if parameters is None else dict(parameters)
        if _PROTOCOL_PARAMETERS.intersection(overrides):
            raise ValueError("input protocol parameters must be supplied through HogProtocol")
        source = self.parameter_values
        overrides.update(
            parameter_97=source["parameter_97"] * protocol.nacl_molar / _SOURCE_NACL_MOLAR,
            t_stress=protocol.shock_time_s,
        )
        trajectory = self.kinetic_model.simulate(
            times_s, parameters=overrides, initial_state=initial_state,
            breakpoints=(protocol.shock_time_s, protocol.shock_time_s + _SOURCE_MIXING_SECONDS),
            rtol=rtol, atol=atol)
        metadata = {
            **trajectory.metadata, "hog_protocol": asdict(protocol),
            "source_mixing_seconds": _SOURCE_MIXING_SECONDS,
            "dose_mapping": "source convention: 0.4 M NaCl corresponds to an added model osmolarity of 0.8; linear scaling is an explicit extrapolation",
            "source_model_sha256": _SOURCE_SHA256,
            "parameter_overrides": {} if parameters is None else dict(parameters),
            "source_assumptions": self.source_metadata.get("source_assumptions", {}),
            "biological_validation": False,
            "quantity_basis": "source-model native variables; use source measurement transforms and explicit observation scales, not silently assumed molarity",
        }
        return replace(trajectory, metadata=metadata)


def hog_metabolic_signals(trajectory: Trajectory, control: Trajectory) -> dict[str, np.ndarray]:
    if not np.array_equal(trajectory.times_s, control.times_s):
        raise ValueError("stressed and control trajectories must share an observation grid")
    baseline = control.variables["Gpd1_measured"]
    current = trajectory.variables["Gpd1_measured"]
    if (not np.isfinite(baseline).all() or not np.isfinite(current).all()
            or np.any(baseline <= 0) or np.any(current < 0)):
        raise ValueError("Gpd1 requires a positive control and finite nonnegative protein amount")
    signals = {"time_s": trajectory.times_s.copy(), "gpd1_capacity_multiplier": current / baseline}
    total_hog = trajectory.variables["Hog1"] + trajectory.variables["Hog1PP"]
    if np.any(total_hog <= 0):
        raise ValueError("phosphorylated fraction requires positive modeled total Hog1")
    signals["modeled_hog1_phosphorylated_fraction"] = trajectory.variables["Hog1PP"] / total_hog
    signals["intracellular_glycerol_native"] = trajectory.variables["glycerol_measured"].copy()
    signals["extracellular_glycerol_native"] = trajectory.variables["glycerol_e"].copy()
    return signals


@dataclass(frozen=True)
class HogGlycerolBalance:
    times_s: np.ndarray
    reference_rates: dict[str, np.ndarray]
    extracellular_rates: dict[str, np.ndarray]
    metadata: dict

    def __post_init__(self):
        times = np.array(self.times_s, dtype=float, copy=True)
        times.setflags(write=False)
        object.__setattr__(self, "times_s", times)
        for name in ("reference_rates", "extracellular_rates"):
            rates = {}
            for key, value in getattr(self, name).items():
                array = np.array(value, dtype=float, copy=True)
                array.setflags(write=False)
                rates[key] = array
            object.__setattr__(self, name, rates)
        object.__setattr__(self, "metadata", deepcopy(self.metadata))

    def carbon_commitment(self) -> np.ndarray:
        synthesis = self.reference_rates["synthesis"]
        glucose = self.reference_rates["glucose_uptake"]
        components = (synthesis, self.reference_rates["synthesis_v6"], self.reference_rates["synthesis_v6b"])
        if (any(rate.shape != self.times_s.shape or not np.isfinite(rate).all() or np.any(rate < 0)
                for rate in components) or glucose.shape != self.times_s.shape
                or not np.isfinite(glucose).all() or np.any(glucose <= 0)):
            raise ValueError("carbon commitment requires nonnegative irreversible synthesis and positive glucose uptake")
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            fraction = 0.5 * (synthesis / glucose)
        if not np.isfinite(fraction).all() or np.any(fraction > 1):
            raise ValueError("synthesis exceeds current glucose carbon; carbon-credit/remobilization mapping is unsupported")
        return fraction


def hog_glycerol_balance(model: HogModel, trajectory: Trajectory) -> HogGlycerolBalance:
    if not isinstance(model, HogModel) or not isinstance(trajectory, Trajectory):
        raise ValueError("native balances require a HogModel and its source Trajectory")
    kinetic = model.kinetic_model
    source = kinetic.metadata
    recorded = trajectory.metadata
    if (source["model_id"] != _SOURCE_MODEL_ID or source["sha256"] != _SOURCE_SHA256
            or recorded.get("source_model_sha256") != _SOURCE_SHA256
            or recorded.get("sha256") != _SOURCE_SHA256
            or recorded.get("model_id") != _SOURCE_MODEL_ID
            or tuple(recorded.get("species_ids", ())) != kinetic.species_ids
            or tuple(recorded.get("reaction_ids", ())) != kinetic.reaction_ids):
        raise ValueError("native glycerol roles require the audited WT source model and trajectory")
    try:
        protocol = HogProtocol(**recorded["hog_protocol"])
        overrides = dict(recorded["parameter_overrides"])
    except (KeyError, TypeError) as exc:
        raise ValueError("source trajectory must record its effective HOG protocol and parameters") from exc
    if _PROTOCOL_PARAMETERS.intersection(overrides):
        raise ValueError("source trajectory protocol parameters must be recorded through HogProtocol")
    overrides.update(
        parameter_97=model.parameter_values["parameter_97"] * protocol.nacl_molar / _SOURCE_NACL_MOLAR,
        t_stress=protocol.shock_time_s,
    )
    times = np.asarray(trajectory.times_s, dtype=float)
    states = np.asarray(trajectory.states, dtype=float)
    if (times.ndim != 1 or not times.size or not np.isfinite(times).all()
            or np.any(times < 0) or not np.all(np.diff(times) > 0)
            or states.shape != (len(times), len(kinetic.species_ids)) or not np.isfinite(states).all()):
        raise ValueError("source trajectory requires finite ordered times and matching finite states")
    contexts = [kinetic.evaluate(t, state, overrides) for t, state in zip(times, states, strict=True)]

    def variable(name):
        return np.array([context[name] for context in contexts])

    volume, reference_volume = variable("cellvol"), variable("cellvol_init")
    intra, extra = variable("intra"), variable("extra")
    if np.any(volume <= 0) or np.any(reference_volume <= 0):
        raise ValueError("native physical cell volume and reference volume must be positive")
    derivatives = np.stack([kinetic.rhs(t, state, overrides) for t, state in zip(times, states, strict=True)])
    reactions = np.stack([kinetic.reaction_rates(t, state, overrides) for t, state in zip(times, states, strict=True)])

    def rhs(name):
        return derivatives[:, kinetic.species_ids.index(name)]

    def rate(name):
        return reactions[:, kinetic.reaction_ids.index(name)]

    factor = volume / (reference_volume * intra)
    synthesis = factor * (rate("v6") + rate("v6b"))
    passive, active = factor * rate("v13a"), factor * rate("v13b")
    outward = np.where(passive >= 0, passive, 0.) + np.where(active < 0, -active, 0.)
    inward = np.where(active >= 0, active, 0.) + np.where(passive < 0, -passive, 0.)
    chain_rule = variable("glycerol_i") * rhs("cellvol") / reference_volume
    retention = (volume * rhs("glycerol_i") + variable("glycerol_i") * rhs("cellvol")) / reference_volume
    external_transport = (rate("v13a") - rate("v13b")) / extra
    batch_correction = (rate("v13aBatch") - rate("v13bBatch")) / extra
    external_glucose = (rate("v1") + rate("v1Batch")) / extra
    reference_rates = {
        "synthesis": synthesis,
        "synthesis_v6": factor * rate("v6"),
        "synthesis_v6b": factor * rate("v6b"),
        "retention": retention,
        "outward_transport": outward,
        "inward_transport": inward,
        "passive_transport": passive,
        "active_import": active,
        "glucose_uptake": factor * rate("v1"),
        "volume_dilution": -factor * rate("vVglyci"),
        "volume_chain_rule": chain_rule,
        "balance_residual": retention - (synthesis - outward + inward),
    }
    extracellular_rates = {
        "net_transport": external_transport,
        "batch_correction": batch_correction,
        "net_accumulation": rhs("glycerol_e"),
        "balance_residual": rhs("glycerol_e") - (external_transport + batch_correction),
        "glucose_uptake": external_glucose,
        "glucose_batch_correction": rate("v1Batch") / extra,
        "glucose_balance_residual": rhs("glucose_e") + external_glucose,
    }
    metadata = {
        "source_model_id": _SOURCE_MODEL_ID,
        "source_model_sha256": _SOURCE_SHA256,
        "hog_protocol": asdict(protocol),
        "parameter_overrides": deepcopy(recorded["parameter_overrides"]),
        "reference_rate_unit": "mol/L_reference/s (source metabolite convention)",
        "extracellular_rate_unit": "mol/L_external/s (source metabolite convention)",
        "declared_sbml_volume_unit": source["units"]["definitions"]["volume"],
        "constant_sbml_compartments": source["units"]["compartments"],
        "unit_caveat": "SBML declares millilitres and defaults substance to mole, but source HPLC glycerol is mol/L; no global 1000-fold rescaling of mixed native variables is justified. Protein native values are not established molarity.",
        "physical_volume": "cellvol is a dynamic physical-volume species, not the constant SBML intra compartment; cellvol_init is the observation reference volume",
        "derivative_method": "source rhs chain rule at each supplied state, reevaluating source laws with effective protocol and parameters; no finite differences",
        "reference_equation": "d(glycerol_i*cellvol/cellvol_init)/dt = (cellvol*rhs_glycerol_i + glycerol_i*rhs_cellvol)/cellvol_init",
        "reaction_basis": "multiply each native reaction extent rate by cellvol/(cellvol_init*intra) for the same reference basis, including v1; do not divide synthesis by extracellular batch depletion",
        "balance_equation": "retention = synthesis - outward_transport + inward_transport; vVglyci dilution cancels the volume chain-rule term",
        "extracellular_equation": "rhs_glycerol_e = (v13a - v13b + v13aBatch - v13bBatch)/extra; -rhs_glucose_e = (v1 + v1Batch)/extra",
        "reaction_roles": {
            "new_glycerol_synthesis": ("v6", "v6b"),
            "signed_passive_outward_transport": "v13a",
            "active_inward_transport": "v13b",
            "glucose_uptake": "v1",
            "physical_volume": "vVos",
            "explicit_glycerol_volume_dilution": "vVglyci",
            "extracellular_batch_corrections": ("v13aBatch", "v13bBatch", "v1Batch"),
        },
        "transport_convention": "gross directed contributions of modeled reaction paths; signed passive and active rates are preserved. v13a is itself a net diffusive law, so unidirectional molecular turnover is not identified. Sign splitting conserves even tiny regulator roundoff; no source state clipping.",
        "growth_convention": "OD and cellnum are imposed experimental fits, not GEM growth; explicit batch corrections are external only; no added mu dilution",
        "native_training_scope": "glucose_only",
        "source_medium": "W303 YPD; source extracellular glycerol is allowed and its uptake is not new synthesis",
        "medium_dependence": "the native coefficient can depend on medium, time, protocol and parameters; transfer beyond the glucose-only training scope requires separate native evidence",
        "carbon_commitment_equation": "3*synthesis/(6*glucose_uptake) = (v6+v6b)/(2*v1); instantaneous total new synthesis, not retained or secreted flux",
        "carbon_origin_identified": False,
        "source_trehalose_remobilization_present": bool(np.any(rate("v3") < 0)),
        "carbon_scope": "an operational synthesis/current-glucose ratio, not atom tracing; source transient precursor and storage pools exist. No external glycerol credit or inventory release is transferred to the steady-state GEM. Ratios outside [0,1] and nonpositive glucose are refused.",
        "absolute_flux_identified": False,
        "required_absolute_calibration": "independently measured source-reference accessible volume in L_reference/gDW at the modeled state, dry-mass basis, provenance and uncertainty; no default is supplied by native physiology/stress_pools",
        "absolute_conversion_equation": "q_mmol_gDW_h = rate_mol_Lreference_s * calibrated_Lreference_gDW * 1000 * 3600",
        "fate_specific_gem_flux_identified": False,
        "biological_validation": False,
    }
    return HogGlycerolBalance(times, reference_rates, extracellular_rates, metadata)
