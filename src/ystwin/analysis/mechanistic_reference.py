from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from numbers import Real
from pathlib import Path

import libsbml
import numpy as np
import scipy

from ystwin.analysis.parameter_evidence import (
    EvidenceGap,
    EvidenceInventory,
    FrozenEvidence,
    JALIHAL_EXPORT_ADAPTER,
    JALIHAL_NORMALIZED_INPUTS,
    NativeReferenceTrajectory,
    PROGRAMMES,
    canonical_json,
    content_sha256,
    freeze_evidence,
    load_parameter_evidence,
    load_reference_model,
    repository_root,
    simulate_native_reference,
    source_path,
    verify_evidence_freeze,
)
from ystwin.mech.kinetic_sbml import KineticModel, Trajectory

STRUCTURES = ("published", "no_hog1_gpd1_edge")
HOLDOUT_AXES = ("environment", "genotype", "exposure_history", "structure")
OBSERVABLES = (
    "Hog1PP_measured", "Gpd1_measured", "glycerol_measured", "glycerol_e",
    "gpd1mRNA_measured", "relVM",
)
SOURCE_RAMP_SECONDS = 5.0
SOURCE_OSMOTIC_PARTICLES_PER_NACL = 2.0
NATIVE_TIME_UNIT = "source_native_unresolved"
NATIVE_HOLDOUT_AXES = ("input_combination", "exposure_history", "parameter", "structure")
NATIVE_OBSERVABLES = {
    "carbon_pka_snf1": ("PKA", "Snf1", "Mig1"),
    "nitrogen_tor": ("TORC1", "Sch9"),
}
NATIVE_STRUCTURAL_TARGETS = {
    "w_pka_camp": "PKA_r", "w_mig_snf": "Mig1_r", "w_sch9_torc": "Sch9_r",
}
NATIVE_SOLVER = {"method": "Radau", "rtol": 1e-10, "atol": 1e-12}
_IMPLEMENTATION_FILES = (
    "src/ystwin/analysis/mechanistic_reference.py",
    "src/ystwin/analysis/parameter_evidence.py",
    "src/ystwin/mech/kinetic_sbml.py", "src/ystwin/mech/params.py",
)


def _number(value, label: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite real number")
    number = float(value)
    if nonnegative and number < 0:
        raise ValueError(f"{label} must be nonnegative")
    return number


def _name(value, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be nonempty text")
    return value


@dataclass(frozen=True)
class SaltChange:
    time_s: float
    nacl_molar: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_s", _number(self.time_s, "salt-change time", nonnegative=True))
        object.__setattr__(self, "nacl_molar", _number(self.nacl_molar, "NaCl concentration", nonnegative=True))


@dataclass(frozen=True)
class ReferenceExperiment:
    experiment_id: str
    model_key: str
    times_s: tuple[float, ...]
    salt_history: tuple[SaltChange, ...]
    split: str
    structure: str = "published"
    holdout_axis: str | None = None

    def __post_init__(self) -> None:
        _name(self.experiment_id, "experiment id")
        _name(self.model_key, "model key")
        times = tuple(_number(value, "sample time", nonnegative=True) for value in self.times_s)
        if len(times) < 2 or times[0] != 0 or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("sample times must start at zero and strictly increase, with at least two samples")
        history = tuple(self.salt_history)
        if not history or not all(isinstance(change, SaltChange) for change in history):
            raise ValueError("salt history requires explicit SaltChange records")
        if any(b.time_s < a.time_s + SOURCE_RAMP_SECONDS for a, b in zip(history, history[1:])):
            raise ValueError("salt changes must be ordered and their five-second ramps must not overlap")
        if any(change.time_s + SOURCE_RAMP_SECONDS > times[-1] for change in history):
            raise ValueError("every salt ramp must finish within the observation horizon")
        if self.structure not in STRUCTURES or self.split not in ("train", "test"):
            raise ValueError("unknown structure or split")
        if (self.split == "train" and self.holdout_axis is not None) or (
            self.split == "test" and self.holdout_axis not in HOLDOUT_AXES
        ):
            raise ValueError("test experiments require a declared holdout axis; training experiments do not")
        object.__setattr__(self, "times_s", times)
        object.__setattr__(self, "salt_history", history)


@dataclass(frozen=True)
class AssayScale:
    observable_id: str
    gain: float
    background: float
    noise_sd: float
    unit: str = "synthetic_assay_unit"

    def __post_init__(self) -> None:
        if self.observable_id not in OBSERVABLES:
            raise ValueError("observable is not in the native reference allowlist")
        _name(self.unit, "assay unit")
        for field in ("gain", "background", "noise_sd"):
            object.__setattr__(self, field, _number(getattr(self, field), field, nonnegative=field != "background"))
        if self.gain <= 0:
            raise ValueError("assay gain must be positive")


@dataclass(frozen=True)
class ReferenceDesign:
    experiments: tuple[ReferenceExperiment, ...]
    assays: tuple[AssayScale, ...]
    seed: int

    def __post_init__(self) -> None:
        experiments, assays = tuple(self.experiments), tuple(self.assays)
        if not experiments or not all(isinstance(e, ReferenceExperiment) for e in experiments):
            raise ValueError("design requires explicit experiments")
        if not assays or not all(isinstance(a, AssayScale) for a in assays):
            raise ValueError("design requires explicit, separate assay scales")
        if len({e.experiment_id for e in experiments}) != len(experiments):
            raise ValueError("duplicate experiment id")
        if len({a.observable_id for a in assays}) != len(assays):
            raise ValueError("duplicate assay observable")
        if {e.split for e in experiments} != {"train", "test"}:
            raise ValueError("a recovery design needs train and test experiments")
        groups = {}
        for experiment in experiments:
            key = (experiment.model_key, experiment.structure, experiment.salt_history)
            if key in groups and groups[key] != experiment.split:
                raise ValueError("the same biological experiment cannot cross train/test boundaries")
            groups[key] = experiment.split
        if type(self.seed) is not int or self.seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        object.__setattr__(self, "experiments", experiments)
        object.__setattr__(self, "assays", assays)

    def to_dict(self) -> dict:
        return json.loads(canonical_json(asdict(self)))

    @classmethod
    def from_dict(cls, payload: dict) -> ReferenceDesign:
        if not isinstance(payload, dict) or set(payload) != {"experiments", "assays", "seed"}:
            raise ValueError("unexpected reference design fields; learner coefficients/outcomes are not inputs")
        try:
            experiments = []
            for record in payload["experiments"]:
                record = dict(record)
                record["salt_history"] = tuple(SaltChange(**change) for change in record["salt_history"])
                experiments.append(ReferenceExperiment(**record))
            return cls(tuple(experiments), tuple(AssayScale(**assay) for assay in payload["assays"]), payload["seed"])
        except (KeyError, TypeError) as exc:
            raise ValueError(f"invalid reference design: {exc}") from exc


@dataclass(frozen=True)
class NutrientInputs:
    Carbon: float
    ATP: float
    Glutamine_ext: float
    NH4: float
    Proline: float

    def __post_init__(self) -> None:
        for name in JALIHAL_NORMALIZED_INPUTS:
            value = _number(getattr(self, name), f"normalized source input {name}")
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be a normalized source input in [0, 1], not a physical dose")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class NutrientChange:
    time_native: float
    inputs: NutrientInputs

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_native", _number(self.time_native, "native change time", nonnegative=True))
        if not isinstance(self.inputs, NutrientInputs):
            raise ValueError("declare all five normalized source axes using NutrientInputs; no physical feature conversion")


@dataclass(frozen=True)
class SourceCounterfactual:
    parameter_id: str
    factor: float
    kind: str = "parameter_scale"

    def __post_init__(self) -> None:
        _name(self.parameter_id, "source parameter id")
        object.__setattr__(self, "factor", _number(self.factor, "source-coordinate multiplier", nonnegative=True))
        if self.kind not in ("parameter_scale", "structural_zero"):
            raise EvidenceGap("only source parameter_scale or structural_zero counterfactuals are supported, not gene interventions")
        if self.factor == 1 or (self.kind == "structural_zero" and self.factor != 0):
            raise ValueError("counterfactuals must change the source coordinate; structural_zero requires factor zero")


@dataclass(frozen=True)
class NativeReferenceExperiment:
    experiment_id: str
    model_key: str
    times_native: tuple[float, ...]
    input_history: tuple[NutrientChange, ...]
    split: str
    holdout_axis: str | None = None
    counterfactuals: tuple[SourceCounterfactual, ...] = ()
    time_unit: str = NATIVE_TIME_UNIT

    def __post_init__(self) -> None:
        _name(self.experiment_id, "experiment id")
        _name(self.model_key, "model key")
        if self.time_unit != NATIVE_TIME_UNIT:
            raise EvidenceGap("native source clock is unresolved: use times_native, not seconds, minutes or hours")
        times = tuple(_number(value, "native sample time", nonnegative=True) for value in self.times_native)
        if len(times) < 2 or times[0] != 0 or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("native sample times must start at zero and strictly increase, with at least two samples")
        history = tuple(self.input_history)
        if not history or not all(isinstance(change, NutrientChange) for change in history):
            raise ValueError("native history requires explicit NutrientChange records")
        if history[0].time_native != 0 or any(b.time_native <= a.time_native for a, b in zip(history, history[1:])):
            raise ValueError("native input changes must start at zero and strictly increase")
        if history[-1].time_native >= times[-1]:
            raise ValueError("every native input change must precede the observation horizon")
        if any(a.inputs == b.inputs for a, b in zip(history, history[1:])):
            raise ValueError("redundant native input changes cannot define a new held-out experiment")
        counterfactuals = tuple(self.counterfactuals)
        if not all(isinstance(change, SourceCounterfactual) for change in counterfactuals):
            raise ValueError("counterfactuals require explicit SourceCounterfactual records")
        if len({change.parameter_id for change in counterfactuals}) != len(counterfactuals):
            raise ValueError("duplicate source counterfactual coordinate")
        if self.split not in ("train", "test") or (self.split == "train" and self.holdout_axis is not None):
            raise ValueError("unknown split or training holdout axis")
        if self.split == "test":
            if self.holdout_axis not in NATIVE_HOLDOUT_AXES:
                raise EvidenceGap("declare a native input_combination, exposure_history, parameter or structure holdout; biological genotypes are unsupported")
            if self.holdout_axis == "input_combination" and len(history) != 1:
                raise ValueError("an input-combination holdout is a constant-input experiment; use exposure_history for changes")
            if self.holdout_axis == "exposure_history" and len(history) < 2:
                raise ValueError("an exposure-history challenge requires a change after time zero")
            if self.holdout_axis in ("parameter", "structure") and not counterfactuals:
                raise ValueError("parameter/structure holdouts require explicit source counterfactuals")
            if (self.holdout_axis == "structure") != any(c.kind == "structural_zero" for c in counterfactuals):
                raise ValueError("structural-zero challenges require the structure holdout axis")
            if counterfactuals and self.holdout_axis not in ("parameter", "structure"):
                raise ValueError("counterfactual challenges must be scored separately from input/history holdouts")
        object.__setattr__(self, "times_native", times)
        object.__setattr__(self, "input_history", history)
        object.__setattr__(self, "counterfactuals", tuple(sorted(counterfactuals, key=lambda change: change.parameter_id)))


@dataclass(frozen=True)
class NativeAssayScale:
    model_key: str
    programme: str
    observable_id: str
    gain: float = 1.0
    background: float = 0.0
    noise_sd: float = 0.0
    unit: str = "synthetic_assay_unit"

    def __post_init__(self) -> None:
        for name in ("model_key", "programme", "observable_id"):
            _name(getattr(self, name), name)
        if self.unit != "synthetic_assay_unit":
            raise EvidenceGap("native states have no verified physical assay calibration; use synthetic_assay_unit")
        for name in ("gain", "background", "noise_sd"):
            object.__setattr__(self, name, _number(getattr(self, name), name, nonnegative=name != "background"))
        if self.gain <= 0:
            raise ValueError("assay gain must be positive")


@dataclass(frozen=True)
class MultiSourceReferenceDesign:
    experiments: tuple[ReferenceExperiment | NativeReferenceExperiment, ...]
    assays: tuple[AssayScale | NativeAssayScale, ...]
    seed: int

    def __post_init__(self) -> None:
        experiments, assays = tuple(self.experiments), tuple(self.assays)
        if not experiments or not all(isinstance(e, (ReferenceExperiment, NativeReferenceExperiment)) for e in experiments):
            raise ValueError("multi-source design requires explicitly typed source experiments")
        if not assays or not all(isinstance(a, (AssayScale, NativeAssayScale)) for a in assays):
            raise ValueError("multi-source design requires explicitly typed source assay scales")
        if len({e.experiment_id for e in experiments}) != len(experiments):
            raise ValueError("duplicate experiment id across sources")
        hog = tuple(e for e in experiments if isinstance(e, ReferenceExperiment))
        native = tuple(e for e in experiments if isinstance(e, NativeReferenceExperiment))
        if not hog or not native:
            raise ValueError("a multi-source design requires both HOG and native nutrient experiments")
        ReferenceDesign(hog, tuple(a for a in assays if isinstance(a, AssayScale)), self.seed)
        native_assays = tuple(a for a in assays if isinstance(a, NativeAssayScale))
        keys = {(a.model_key, a.observable_id) for a in native_assays}
        if len(keys) != len(native_assays):
            raise ValueError("duplicate native assay observable")
        if {a.model_key for a in native_assays} != {e.model_key for e in native}:
            raise ValueError("every native source requires its own assay declarations, without unused assays")
        groups = {}
        for experiment in native:
            key = (experiment.model_key, experiment.input_history,
                   tuple((c.parameter_id, c.factor) for c in experiment.counterfactuals))
            if key in groups and groups[key] != experiment.split:
                raise ValueError("the same native experiment cannot cross train/test boundaries, even with different sample times or labels")
            groups[key] = experiment.split
        for model_key in {e.model_key for e in native}:
            if {e.split for e in native if e.model_key == model_key} != {"train", "test"}:
                raise ValueError("every native source requires separate training and test experiments")
        object.__setattr__(self, "experiments", experiments)
        object.__setattr__(self, "assays", assays)

    def to_dict(self) -> dict:
        return json.loads(canonical_json({
            "design_kind": "multi_source_native_v1", "seed": self.seed,
            "experiments": [{"kind": "hog" if isinstance(e, ReferenceExperiment) else "source_native", **asdict(e)}
                            for e in self.experiments],
            "assays": [{"kind": "hog" if isinstance(a, AssayScale) else "source_native", **asdict(a)} for a in self.assays],
        }))

    @classmethod
    def from_dict(cls, payload: dict) -> MultiSourceReferenceDesign:
        if not isinstance(payload, dict) or set(payload) != {"design_kind", "experiments", "assays", "seed"} or payload["design_kind"] != "multi_source_native_v1":
            raise ValueError("unexpected multi-source design fields; learner coefficients, gene outcomes and physical feature conversions are not inputs")
        try:
            experiments, assays = [], []
            for record in payload["experiments"]:
                record = dict(record)
                kind = record.pop("kind")
                if kind == "hog":
                    record["salt_history"] = tuple(SaltChange(**change) for change in record["salt_history"])
                    experiments.append(ReferenceExperiment(**record))
                elif kind == "source_native":
                    record["input_history"] = tuple(_nutrient_change_from_dict(change) for change in record["input_history"])
                    record["counterfactuals"] = tuple(SourceCounterfactual(**change) for change in record.get("counterfactuals", ()))
                    experiments.append(NativeReferenceExperiment(**record))
                else:
                    raise ValueError("unknown declared experiment kind")
            for record in payload["assays"]:
                record = dict(record)
                kind = record.pop("kind")
                if kind not in ("hog", "source_native"):
                    raise ValueError("unknown declared assay kind")
                assays.append((AssayScale if kind == "hog" else NativeAssayScale)(**record))
            return cls(tuple(experiments), tuple(assays), payload["seed"])
        except (KeyError, TypeError) as exc:
            raise ValueError(f"invalid explicit source design: {exc}; all five normalized inputs are required, not physical or product inputs") from exc


def _nutrient_change_from_dict(change) -> NutrientChange:
    if not isinstance(change, dict) or set(change) != {"time_native", "inputs"}:
        raise ValueError("native changes require exactly time_native and all five explicit normalized inputs")
    return NutrientChange(change["time_native"], NutrientInputs(**change["inputs"]))


def reference_design_from_dict(payload: dict) -> ReferenceDesign | MultiSourceReferenceDesign:
    if isinstance(payload, dict) and "design_kind" in payload:
        return MultiSourceReferenceDesign.from_dict(payload)
    return ReferenceDesign.from_dict(payload)


def default_reference_design(*, seed: int = 0) -> ReferenceDesign:
    times = (0, 1800, 3300, 3600, 3605, 3660, 3900, 4200, 4800, 5400, 6000, 6600, 7200)
    step = (SaltChange(3600, 0.4),)
    experiments = [
        ReferenceExperiment("train_control", "hog2013_wt", times, (SaltChange(3600, 0),), "train"),
        ReferenceExperiment("train_step", "hog2013_wt", times, step, "train"),
        ReferenceExperiment("test_dose", "hog2013_wt", times, (SaltChange(3600, 0.8),), "test", holdout_axis="environment"),
        ReferenceExperiment("test_history", "hog2013_wt", times,
                            (*step, SaltChange(4200, 0), SaltChange(4800, 0.4)), "test", holdout_axis="exposure_history"),
        ReferenceExperiment("test_structure", "hog2013_wt", times, step, "test",
                            structure="no_hog1_gpd1_edge", holdout_axis="structure"),
    ]
    for genotype in ("hog1_del", "hog1_att", "fps1_open", "gpd1_del", "pfk2627_del"):
        experiments.append(ReferenceExperiment(
            f"test_{genotype}", f"hog2013_{genotype}", times, step, "test", holdout_axis="genotype",
        ))
    return ReferenceDesign(
        tuple(experiments),
        tuple(AssayScale(name, gain=1.0, background=0.0, noise_sd=0.0) for name in OBSERVABLES[:4]),
        seed,
    )


def default_multisource_reference_design(*, seed: int = 0) -> MultiSourceReferenceDesign:
    hog = default_reference_design(seed=seed)
    times = (0, 0.1, 1, 5, 10, 20, 30, 40, 45, 60)
    baseline = NutrientInputs(0.5, 0.5, 0, 0, 0)
    experiments = [NativeReferenceExperiment(
        "nutrient_train_control", "jalihal2021", times, (NutrientChange(0, baseline),), "train",
    )]
    for axis in JALIHAL_NORMALIZED_INPUTS:
        varied = NutrientInputs(**{**asdict(baseline), axis: 1.0})
        experiments.extend((
            NativeReferenceExperiment(f"nutrient_train_{axis}", "jalihal2021", times,
                                      (NutrientChange(0, varied),), "train"),
            NativeReferenceExperiment(f"nutrient_test_history_{axis}", "jalihal2021", times,
                                      (NutrientChange(0, baseline), NutrientChange(20, varied), NutrientChange(40, baseline)),
                                      "test", holdout_axis="exposure_history"),
        ))
    for index, values in enumerate(((1, 0.25, 0.5, 0.25, 0.75), (0.25, 1, 0.75, 0.5, 0.25))):
        experiments.append(NativeReferenceExperiment(
            f"nutrient_test_combination_{index}", "jalihal2021", times,
            (NutrientChange(0, NutrientInputs(*values)),), "test", holdout_axis="input_combination",
        ))
    rich = NutrientInputs(0.5, 0.5, 1, 0, 0)
    experiments.append(NativeReferenceExperiment(
        "nutrient_test_parameter", "jalihal2021", times, (NutrientChange(0, rich),), "test",
        holdout_axis="parameter", counterfactuals=(SourceCounterfactual("gammasnf", 0.5),),
    ))
    for index, parameter_id in enumerate(NATIVE_STRUCTURAL_TARGETS):
        experiments.append(NativeReferenceExperiment(
            f"nutrient_test_structure_{index}", "jalihal2021", times, (NutrientChange(0, rich),), "test",
            holdout_axis="structure", counterfactuals=(SourceCounterfactual(parameter_id, 0, "structural_zero"),),
        ))
    assays = tuple(NativeAssayScale("jalihal2021", programme, observable)
                   for programme, observables in NATIVE_OBSERVABLES.items() for observable in observables)
    return MultiSourceReferenceDesign((*hog.experiments, *experiments), (*hog.assays, *assays), seed)


def _implementation_identity() -> dict:
    root = repository_root()
    return {
        "files": {relative: hashlib.sha256(source_path(root, relative).read_bytes()).hexdigest()
                  for relative in _IMPLEMENTATION_FILES},
        "numpy": np.__version__, "scipy": scipy.__version__, "libsbml": libsbml.getLibSBMLDottedVersion(),
    }


@dataclass(frozen=True)
class FrozenReference:
    _json: str

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self._json.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict:
        return {"reference_id": self.sha256, **json.loads(self._json)}

    @classmethod
    def from_dict(cls, payload: dict) -> FrozenReference:
        if not isinstance(payload, dict) or type(payload.get("schema_version")) is not int or payload["schema_version"] not in (1, 2):
            raise ValueError("invalid frozen reference schema")
        fields = {"schema_version", "reference_id", "design", "evidence", "implementation", "shared_assumptions", "scope"}
        if payload["schema_version"] == 2:
            fields |= {"source_contracts", "input_split", "source_counterfactuals"}
        if set(payload) != fields:
            raise ValueError("invalid frozen reference schema")
        body = {key: value for key, value in payload.items() if key != "reference_id"}
        if content_sha256(body) != payload["reference_id"]:
            raise ValueError("frozen reference content digest mismatch")
        design_type = ReferenceDesign if payload["schema_version"] == 1 else MultiSourceReferenceDesign
        design_type.from_dict(body["design"])
        return cls(canonical_json(body))


def freeze_reference(
    design: ReferenceDesign | MultiSourceReferenceDesign, *, product_outcomes_seen: bool,
    inventory: EvidenceInventory | None = None, root=None,
) -> FrozenReference:
    if not isinstance(design, (ReferenceDesign, MultiSourceReferenceDesign)):
        raise ValueError("a validated ReferenceDesign or MultiSourceReferenceDesign is required")
    if product_outcomes_seen is not False:
        raise ValueError("freeze requires an explicit no-product-outcomes declaration for this run")
    inventory = load_parameter_evidence() if inventory is None else inventory
    if isinstance(design, MultiSourceReferenceDesign):
        return _freeze_multisource_reference(design, inventory, root=root)
    data = inventory.to_dict()
    allowed = inventory.require_reference("osmotic_hog")
    if any(experiment.model_key not in allowed for experiment in design.experiments):
        raise EvidenceGap("only the pinned local HOG mechanisms accept ReferenceDesign; declare MultiSourceReferenceDesign for native nutrient sources")
    evidence = freeze_evidence(inventory, source_ids=("hog2013",), product_outcomes_seen=product_outcomes_seen, root=root)
    return FrozenReference(canonical_json({
        "schema_version": 1, "design": design.to_dict(), "evidence": evidence.to_dict(),
        "implementation": _implementation_identity(), "shared_assumptions": data["shared_assumptions"],
        "scope": "synthetic conditional reference; freeze before generation; not full biology or independent real validation",
    }))


def verify_reference_freeze(
    frozen: FrozenReference, inventory: EvidenceInventory, *, root=None,
) -> ReferenceDesign | MultiSourceReferenceDesign:
    if not isinstance(frozen, FrozenReference):
        raise ValueError("generation/scoring requires a FrozenReference, not unfrozen parameters")
    payload = FrozenReference.from_dict(frozen.to_dict()).to_dict()
    design = reference_design_from_dict(payload["design"])
    verify_evidence_freeze(FrozenEvidence(canonical_json(payload["evidence"])), inventory, root=root)
    current = freeze_reference(design, product_outcomes_seen=False, inventory=inventory, root=root)
    if current.sha256 != frozen.sha256:
        raise ValueError("reference design, implementation, assumptions or evidence changed after freeze")
    return design


def _math_nodes(node):
    if node is None:
        raise EvidenceGap("the verified source expression could not be inspected")
    yield node
    for index in range(node.getNumChildren()):
        yield from _math_nodes(node.getChild(index))


def _source_counterfactuals(experiment: NativeReferenceExperiment, model: KineticModel) -> list[dict]:
    values, metadata = model.parameter_values, model.metadata
    result = []
    for change in experiment.counterfactuals:
        name = change.parameter_id
        if name in JALIHAL_NORMALIZED_INPUTS:
            raise EvidenceGap("nutrient coordinates are imposed inputs, not hidden parameter interventions; declare the input history")
        if name not in values or metadata["parameter_definitions"][name]["assignment_target"]:
            raise EvidenceGap(f"{name}: counterfactual target is not a verified independent source parameter")
        reaction = None
        if change.kind == "structural_zero":
            reaction = NATIVE_STRUCTURAL_TARGETS.get(name)
            if reaction is None:
                raise EvidenceGap("structural_zero supports only the explicitly audited source weights: " + ", ".join(NATIVE_STRUCTURAL_TARGETS))
            formula = metadata["native_reference"]["source_reaction_formulas"][reaction]
            symbols = {node.getName() for node in _math_nodes(libsbml.parseL3Formula(formula)) if node.getType() == libsbml.AST_NAME}
            if name not in symbols:
                raise EvidenceGap("the declared source structural coordinate no longer occurs in its audited equation")
        altered = _number(values[name] * change.factor, "counterfactual source value")
        if altered == values[name]:
            raise ValueError("counterfactual multiplier does not change the actual source value")
        result.append({
            **asdict(change), "source_value": values[name], "reference_value": altered,
            "unit": metadata["units"]["parameters"][name], "factor_unit": "dimensionless",
            "source_reaction_id": reaction,
            "status": "declared_source_counterfactual_not_biological_gene_intervention",
        })
    return result


def _freeze_multisource_reference(design: MultiSourceReferenceDesign, inventory: EvidenceInventory, *, root=None) -> FrozenReference:
    data = inventory.to_dict()
    allowed_hog = inventory.require_reference("osmotic_hog")
    contracts, interventions, selected_sources, native_models = {}, {}, set(), {}
    for experiment in design.experiments:
        record = data["models"].get(experiment.model_key)
        if isinstance(experiment, ReferenceExperiment):
            if experiment.model_key not in allowed_hog:
                raise EvidenceGap("HOG salt-history experiments require a pinned HOG source model")
            contracts[experiment.model_key] = {
                "source_id": record["source_id"], "source_sha256": record["sha256"],
                "programmes": record["programmes"], "times_field": "times_s", "row_time_field": "time_s",
                "time_unit": "s", "physical_seconds_per_native_unit": 1.0,
                "input_axes": ["nacl_molar"], "input_scope": "source YPD/NaCl convention, not arbitrary media or full environment",
            }
        else:
            if (record is None or record["availability"] != "local_verified" or experiment.model_key != "jalihal2021"
                    or record.get("execution_adapter") != JALIHAL_EXPORT_ADAPTER):
                gap = "unknown source" if record is None else record["gap"]
                raise EvidenceGap(f"{experiment.model_key}: no audited native nutrient dispatch; {gap}")
            if experiment.model_key not in native_models:
                model = load_reference_model(experiment.model_key, inventory, root=root)
                formulas = model.metadata["native_reference"]["source_reaction_formulas"].values()
                if any(node.getType() == libsbml.AST_NAME_TIME for formula in formulas
                       for node in _math_nodes(libsbml.parseL3Formula(formula))):
                    raise EvidenceGap("native history dispatch requires the verified autonomous source; absolute-time equations cannot be reset")
                native_models[experiment.model_key] = model
            model = native_models[experiment.model_key]
            interventions[experiment.experiment_id] = _source_counterfactuals(experiment, model)
            contracts[experiment.model_key] = {
                "source_id": record["source_id"], "source_sha256": record["sha256"],
                "programmes": record["programmes"], "times_field": "times_native", "row_time_field": "time_native",
                "time_unit": NATIVE_TIME_UNIT, "physical_seconds_per_native_unit": None,
                "unit_audit": record["unit_audit"], "input_axes": list(JALIHAL_NORMALIZED_INPUTS),
                "physical_input_calibration": None, "state_scale": "source_normalized_not_absolute_activity",
                "solver": dict(NATIVE_SOLVER), "clock_policy": "autonomous ODE segments carry state; absolute native history never resets state",
                "numerical_scope": record["numerical_audit"]["sensitivity_caveat"],
            }
        selected_sources.add(record["source_id"])
    for assay in design.assays:
        if isinstance(assay, NativeAssayScale):
            record = data["models"][assay.model_key]
            if assay.programme not in record["programmes"] or assay.observable_id not in NATIVE_OBSERVABLES.get(assay.programme, ()):
                raise EvidenceGap("unsupported native programme/observable; nutrient coverage is not whole-environment biological coverage")
            if assay.observable_id not in native_models[assay.model_key].species_ids:
                raise EvidenceGap("declared native observable is absent from the verified source states")
    evidence = freeze_evidence(inventory, source_ids=sorted(selected_sources), product_outcomes_seen=False, root=root)
    input_split = {
        "unit": "whole_source_experiment_not_random_timepoints", "frozen_before_generation": True,
        "train": [e.experiment_id for e in design.experiments if e.split == "train"],
        "test": [e.experiment_id for e in design.experiments if e.split == "test"],
        "experiments_sha256": content_sha256(design.to_dict()["experiments"]),
    }
    supported = {programme for contract in contracts.values() for programme in contract["programmes"]}
    scope = {
        "evaluation_kind": "synthetic_conditional_multisource_recovery",
        "simulated_programmes": sorted(supported),
        "scored_programmes": sorted({"osmotic_hog", *(a.programme for a in design.assays if isinstance(a, NativeAssayScale))}),
        "unsupported_programmes": {name: data["programmes"][name]["blockers"] for name in PROGRAMMES if name not in supported},
        "unsupported_axes": ["physical_glucose_or_nitrogen_dose", "temperature", "peroxide", "pH", "oxygen", "biological_gene_interventions", "product_inputs_or_outcomes"],
        "source_coupling": "separate source-native challenges; no invented HOG/nutrient cross-programme kinetics",
        "full_environment_coverage": False, "full_biology_claim": False, "independent_real_validation": False,
        "genotype_scope": "HOG source labels remain published conditional fits; nutrient counterfactuals are source coordinates, not gene deletions or alleles",
        "physical_feature_conversions": False,
    }
    return FrozenReference(canonical_json({
        "schema_version": 2, "design": design.to_dict(), "evidence": evidence.to_dict(),
        "implementation": _implementation_identity(), "shared_assumptions": data["shared_assumptions"],
        "source_contracts": contracts, "input_split": input_split, "source_counterfactuals": interventions, "scope": scope,
    }))


def _input_formula(history: Sequence[SaltChange]) -> str:
    terms, previous = [], 0.0
    for change in history:
        delta = SOURCE_OSMOTIC_PARTICLES_PER_NACL * (change.nacl_molar - previous)
        start, end = change.time_s, change.time_s + SOURCE_RAMP_SECONDS
        ramp = f"piecewise(0, time < {start:.17g}, (time - {start:.17g}) / {SOURCE_RAMP_SECONDS:.17g}, time < {end:.17g}, 1)"
        terms.append(f"({delta:.17g} / parameter_97) * {ramp}")
        previous = change.nacl_molar
    return " + ".join(terms)


def _simulate(experiment: ReferenceExperiment, inventory: EvidenceInventory, *, root=None) -> tuple[Trajectory, dict]:
    original = load_reference_model(experiment.model_key, inventory, root=root)
    record = inventory.to_dict()["models"][experiment.model_key]
    root = repository_root() if root is None else root
    path = source_path(root, record["path"])
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != record["sha256"]:
        raise ValueError("source bytes changed during reference construction")
    document = libsbml.readSBMLFromString(payload.decode("utf-8-sig"))
    model = document.getModel()
    rule = model.getRuleByVariable("input")
    if rule is None or original.parameter_values.get("parameter_97", 0) <= 0:
        raise EvidenceGap("source lacks the audited salt-input convention")
    formula = _input_formula(experiment.salt_history)
    if rule.setMath(libsbml.parseL3Formula(formula)) != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise ValueError("cannot construct reference input history")
    interventions = []
    if experiment.structure == "no_hog1_gpd1_edge":
        model.getParameter("kv17f_1").setValue(0.0)
        interventions.append({
            "target": "kv17f_1", "source_value": original.parameter_values["kv17f_1"], "reference_value": 0.0,
            "meaning": "remove only the Hog1-dependent additive GPD1 transcription arm; retain AOG and basal synthesis",
            "status": "declared_counterfactual_not_measured_mutant",
        })
    transformed = libsbml.writeSBMLToString(document).encode("utf-8")
    transformed_sha = hashlib.sha256(transformed).hexdigest()
    independent = KineticModel(document, path, transformed_sha)
    knots = sorted({
        *(change.time_s for change in experiment.salt_history),
        *(change.time_s + SOURCE_RAMP_SECONDS for change in experiment.salt_history),
        4800.0, original.parameter_values["kv22_Hog1D_t"],
    })
    trajectory = independent.simulate(experiment.times_s, breakpoints=knots)
    provenance = {
        "model_key": experiment.model_key, "source_path": record["path"], "model_id": independent.model_id,
        "source_sha256": record["sha256"], "transformed_sbml_sha256": transformed_sha,
        "source_parameter_values": original.parameter_values,
        "reference_parameter_values": independent.parameter_values,
        "input_formula": formula, "structural_interventions": interventions,
        "source_clock_knots_s": [4800.0, original.parameter_values["kv22_Hog1D_t"]],
        "clock_policy": "absolute source clock is retained, including empirical mutant terms; histories never reset state",
        "units": independent.metadata["units"],
        "shared_constants": {
            "salt_ramp_s": SOURCE_RAMP_SECONDS,
            "osmotic_particles_per_nacl": SOURCE_OSMOTIC_PARTICLES_PER_NACL,
            "source": "HOG input/OsmoE assignment rules: 0.8 osmotic increment for 0.4 M NaCl over five seconds",
        },
        "source_initial_conditions": trajectory.states[0].tolist(),
        "imported_learner_coefficients": False, "parameter_uncertainty": "source point prior; no covariance supplied",
        "identifiability": {
            "status": "not_established",
            "unresolved_directions": ["native kinetic vector versus assay/abundance scale", "empirical OD/cell number versus cell mass"],
        },
    }
    return trajectory, provenance


def _salt_at(times: np.ndarray, history: Sequence[SaltChange]) -> np.ndarray:
    result, previous = np.zeros_like(times), 0.0
    for change in history:
        ramp = np.minimum(1.0, np.maximum(0.0, (times - change.time_s) / SOURCE_RAMP_SECONDS))
        result += (change.nacl_molar - previous) * ramp
        previous = change.nacl_molar
    return result


def _assay_noise(seed: int, experiment_id: str, observable_id: str, count: int, sd: float) -> np.ndarray:
    key = canonical_json([seed, experiment_id, observable_id]).encode("utf-8")
    stream = int.from_bytes(hashlib.sha256(key).digest()[:16], "big")
    return np.random.default_rng(stream).normal(0.0, sd, count)


@dataclass(frozen=True)
class MechanisticReference:
    frozen: FrozenReference
    _learner_json: str
    _truth_json: str

    def learner_inputs(self) -> dict:
        return json.loads(self._learner_json)

    def evaluator_truth(self) -> dict:
        return json.loads(self._truth_json)


def _generate_hog_rows(design: ReferenceDesign, inventory: EvidenceInventory, *, root=None) -> tuple:
    source_models = inventory.to_dict()["models"]
    training, queries, targets, experiments = [], [], [], []
    for experiment in design.experiments:
        trajectory, provenance = _simulate(experiment, inventory, root=root)
        salt = _salt_at(trajectory.times_s, experiment.salt_history)
        conditions = {
            "experiment_id": experiment.experiment_id,
            "genotype": source_models[experiment.model_key]["genotype"],
            "salt_history": [asdict(change) for change in experiment.salt_history],
            "medium": "source_YPD_context", "time_unit": "s",
        }
        for assay in design.assays:
            latent = trajectory.variables[assay.observable_id]
            expectation = assay.gain * latent + assay.background
            observed = expectation + _assay_noise(design.seed, experiment.experiment_id, assay.observable_id, len(latent), assay.noise_sd)
            for index, time_s in enumerate(trajectory.times_s):
                row = {
                    **conditions, "row_id": f"{experiment.experiment_id}:{assay.observable_id}:{index}",
                    "time_s": float(time_s), "nacl_molar": float(salt[index]),
                    "observable_id": assay.observable_id, "unit": assay.unit,
                    "known_od_input": float(trajectory.variables["OD"][index]),
                }
                if experiment.split == "train":
                    training.append({**row, "value": float(observed[index]), "noise_sd": assay.noise_sd})
                else:
                    queries.append(row)
                    targets.append({
                        **row, "value": float(observed[index]), "expectation": float(expectation[index]),
                        "holdout_axis": experiment.holdout_axis,
                    })
        experiments.append({
            "design": asdict(experiment), "provenance": provenance,
            "species_ids": list(trajectory.metadata["species_ids"]),
            "reaction_ids": list(trajectory.metadata["reaction_ids"]),
            "times_s": trajectory.times_s.tolist(), "states": trajectory.states.tolist(),
            "reaction_rates": trajectory.reaction_rates.tolist(),
        })
    return training, queries, targets, experiments


def generate_reference(frozen: FrozenReference, *, inventory: EvidenceInventory | None = None, root=None) -> MechanisticReference:
    inventory = load_parameter_evidence() if inventory is None else inventory
    design = verify_reference_freeze(frozen, inventory, root=root)
    if isinstance(design, MultiSourceReferenceDesign):
        return _generate_multisource_reference(frozen, design, inventory, root=root)
    training, queries, targets, experiments = _generate_hog_rows(design, inventory, root=root)
    inputs = {
        "schema_version": 1, "reference_id": frozen.sha256, "domain": "synthetic_native_reference",
        "training_observations": training, "prediction_queries": queries,
        "input_scope": "only conditions, imposed OD and training assay observations; no latent truth, source coefficients, assay gain or test outcomes",
    }
    truth = {
        "schema_version": 1, "reference_id": frozen.sha256, "evaluation_kind": "synthetic_recovery",
        "targets": targets, "experiments": experiments, "assay_scales": [asdict(a) for a in design.assays],
        "assay_scope": "declared synthetic observation transforms; not measured molecular calibrations",
        "shared_assumptions": inventory.to_dict()["shared_assumptions"],
        "full_biology_claim": False, "independent_real_validation": False,
        "learner_inputs_sha256": content_sha256(inputs),
    }
    return MechanisticReference(frozen, canonical_json(inputs), canonical_json(truth))


def _simulate_native_history(
    experiment: NativeReferenceExperiment, inventory: EvidenceInventory, *, root=None,
) -> tuple[NativeReferenceTrajectory, dict]:
    model = load_reference_model(experiment.model_key, inventory, root=root)
    metadata = model.metadata
    interventions = _source_counterfactuals(experiment, model)
    overrides = {change["parameter_id"]: change["reference_value"] for change in interventions}
    times = np.asarray(experiment.times_native)
    states = np.empty((len(times), len(model.species_ids)))
    rates = np.empty((len(times), len(model.reaction_ids)))
    variables, segments, state = {}, [], None
    for index, change in enumerate(experiment.input_history):
        last = index == len(experiment.input_history) - 1
        start = change.time_native
        end = times[-1] if last else experiment.input_history[index + 1].time_native
        selected = np.flatnonzero((times >= start) & ((times <= end) if last else (times < end)))
        local_times = sorted({0.0, float(end - start), *(float(times[i] - start) for i in selected)})
        inputs = asdict(change.inputs)
        if overrides:
            trajectory = model.simulate(local_times, parameters={**inputs, **overrides}, initial_state=state, **NATIVE_SOLVER)
        else:
            trajectory = simulate_native_reference(
                experiment.model_key, local_times, normalized_inputs=inputs, initial_state=state,
                inventory=inventory, root=root, **NATIVE_SOLVER,
            )
        if not variables:
            variables = {name: np.empty(len(times)) for name in trajectory.variables}
        positions = {time: i for i, time in enumerate(local_times)}
        take = [positions[float(times[i] - start)] for i in selected]
        states[selected] = trajectory.states[take]
        rates[selected] = trajectory.reaction_rates[take]
        for name in variables:
            variables[name][selected] = trajectory.variables[name][take]
        state = trajectory.states[-1].copy()
        segments.append({
            "start_native": float(start), "end_native": float(end), "normalized_inputs": inputs,
            "parameter_overrides": overrides, "initial_state": trajectory.states[0].tolist(),
            "terminal_state": state.tolist(),
        })
    provenance = {
        "model_key": experiment.model_key, "model_id": model.model_id, "source_id": "jalihal2021",
        "source_path": inventory.to_dict()["models"][experiment.model_key]["path"],
        "source_sha256": metadata["sha256"], "native_reference": metadata["native_reference"],
        "source_parameter_values": model.parameter_values,
        "reference_parameter_values": {name: overrides.get(name, value) for name, value in model.parameter_values.items()
                                       if name not in JALIHAL_NORMALIZED_INPUTS},
        "reference_parameter_scope": "non-input source coordinates only; all five imposed inputs are tracked separately in each history segment",
        "source_counterfactuals": interventions,
        "structural_interventions": [change for change in interventions if change["kind"] == "structural_zero"],
        "counterfactual_scope": "source-coordinate challenges, not measured gene deletions, alleles or biological mutation predictions",
        "units": {
            "time": NATIVE_TIME_UNIT, "physical_seconds_per_native_unit": None,
            "states": "source_normalized_not_absolute_concentration_or_kinase_activity",
            "reaction_rates": "source_normalized_derivative_per_native_numeric_unit",
            "parameters": metadata["units"]["parameters"], "physical_input_calibration": None,
        },
        "clock_policy": "verified autonomous source ODE; piecewise constant inputs carry terminal state, never reset state; report absolute native clock",
        "boundary_policy": "states continuous; imposed inputs and reaction rates are right-continuous at changes, including unobserved knots",
        "segments": segments, "solver": dict(NATIVE_SOLVER), "source_initial_conditions": model.initial_state().tolist(),
        "execution_interface": "verified source-native model.simulate parameter override API" if overrides else "simulate_native_reference",
        "imported_learner_coefficients": False, "parameter_uncertainty": "source fitted point/ensemble, no per-parameter covariance supplied",
        "identifiability": {
            "status": "not_established", "covariance": None,
            "unresolved_directions": ["source-reported acceptable-fit/sloppy directions", "native kinetic weights versus state and assay gain",
                                      "source minute/second clock ambiguity and uncalibrated nutrient inputs"],
        },
    }
    return NativeReferenceTrajectory(times.copy(), states, variables, rates, {
        "time_unit": NATIVE_TIME_UNIT, "species_ids": list(model.species_ids), "reaction_ids": list(model.reaction_ids),
    }), provenance


def _generate_multisource_reference(
    frozen: FrozenReference, design: MultiSourceReferenceDesign, inventory: EvidenceInventory, *, root=None,
) -> MechanisticReference:
    hog = ReferenceDesign(
        tuple(e for e in design.experiments if isinstance(e, ReferenceExperiment)),
        tuple(a for a in design.assays if isinstance(a, AssayScale)), design.seed,
    )
    training, queries, targets, experiments = _generate_hog_rows(hog, inventory, root=root)
    for rows in (training, queries, targets):
        for row in rows:
            row.update({"source_id": "hog2013", "programme": "osmotic_hog"})
    for experiment in experiments:
        experiment.update({"source_id": "hog2013", "programmes": ["osmotic_hog"]})
    frozen_payload = frozen.to_dict()
    for experiment in design.experiments:
        if not isinstance(experiment, NativeReferenceExperiment):
            continue
        trajectory, provenance = _simulate_native_history(experiment, inventory, root=root)
        if provenance["source_counterfactuals"] != frozen_payload["source_counterfactuals"][experiment.experiment_id]:
            raise ValueError("source counterfactuals changed after freeze")
        conditions = {
            "experiment_id": experiment.experiment_id, "source_id": provenance["source_id"], "model_key": experiment.model_key,
            "input_history": [asdict(change) for change in experiment.input_history],
            "time_unit": NATIVE_TIME_UNIT, "input_unit": "source_normalized_0_1",
            "source_condition_id": content_sha256([experiment.model_key, [(c.parameter_id, c.factor) for c in experiment.counterfactuals]]),
        }
        changes = np.searchsorted([c.time_native for c in experiment.input_history], trajectory.times_native, side="right") - 1
        for assay in design.assays:
            if not isinstance(assay, NativeAssayScale) or assay.model_key != experiment.model_key:
                continue
            latent = trajectory.variables[assay.observable_id]
            expectation = assay.gain * latent + assay.background
            observed = expectation + _assay_noise(design.seed, experiment.experiment_id, assay.observable_id, len(latent), assay.noise_sd)
            for index, time_native in enumerate(trajectory.times_native):
                row = {
                    **conditions, "programme": assay.programme,
                    "row_id": f"{experiment.experiment_id}:{assay.observable_id}:{index}",
                    "time_native": float(time_native), "normalized_inputs": asdict(experiment.input_history[changes[index]].inputs),
                    "observable_id": assay.observable_id, "unit": assay.unit,
                }
                if experiment.split == "train":
                    training.append({**row, "value": float(observed[index]), "noise_sd": assay.noise_sd})
                else:
                    queries.append(row)
                    targets.append({**row, "value": float(observed[index]), "expectation": float(expectation[index]),
                                    "holdout_axis": experiment.holdout_axis})
        experiments.append({
            "design": asdict(experiment), "provenance": provenance, "source_id": provenance["source_id"],
            "programmes": frozen_payload["source_contracts"][experiment.model_key]["programmes"],
            "species_ids": trajectory.metadata["species_ids"], "reaction_ids": trajectory.metadata["reaction_ids"],
            "times_native": trajectory.times_native.tolist(), "time_unit": NATIVE_TIME_UNIT,
            "states": trajectory.states.tolist(), "reaction_rates": trajectory.reaction_rates.tolist(),
        })
    inputs = {
        "schema_version": 2, "reference_id": frozen.sha256, "domain": "synthetic_native_reference",
        "training_observations": training, "prediction_queries": queries,
        "input_scope": "only source-labelled conditions, imposed inputs and frozen training assays; opaque source conditions are not gene labels; no test values, latent truth, coefficients, counterfactual targets or assay gains",
    }
    truth = {
        "schema_version": 2, "reference_id": frozen.sha256, "evaluation_kind": "synthetic_recovery",
        "targets": targets, "experiments": experiments, "assay_scales": design.to_dict()["assays"],
        "assay_scope": "declared synthetic observation transforms, not physical molecular calibrations",
        "shared_assumptions": inventory.to_dict()["shared_assumptions"], "scope": frozen_payload["scope"],
        "source_contracts": frozen_payload["source_contracts"], "input_split": frozen_payload["input_split"],
        "full_biology_claim": False, "independent_real_validation": False,
        "learner_inputs_sha256": content_sha256(inputs), "training_observations_sha256": content_sha256(training),
        "prediction_queries_sha256": content_sha256(queries),
    }
    return MechanisticReference(frozen, canonical_json(inputs), canonical_json(truth))


def _score_rows(targets: list[dict], predictions: Sequence[dict]) -> dict:
    if not isinstance(predictions, (list, tuple)):
        raise ValueError("predictions require a list of row_id/prediction/unit records")
    expected = {row["row_id"]: row for row in targets}
    supplied = {}
    for record in predictions:
        if not isinstance(record, dict) or set(record) != {"row_id", "prediction", "unit"}:
            raise ValueError("prediction records require exactly row_id, prediction and unit")
        key = record["row_id"]
        if key not in expected or key in supplied:
            raise ValueError("unknown or duplicate prediction row_id")
        if record["unit"] != expected[key]["unit"]:
            raise ValueError("prediction and target units differ")
        supplied[key] = None if record["prediction"] is None else _number(record["prediction"], "prediction")
    rows = []
    for target in targets:
        prediction = supplied.get(target["row_id"])
        rows.append({
            **target, "prediction": prediction, "scored": prediction is not None,
            "residual": None if prediction is None else prediction - target["value"],
        })
    metrics = []
    groups = sorted({(r["experiment_id"], r["observable_id"], r["unit"]) for r in rows})
    for experiment, observable, unit in groups:
        selected = [r for r in rows if (r["experiment_id"], r["observable_id"], r["unit"]) == (experiment, observable, unit)]
        errors = [r["residual"] for r in selected if r["scored"]]
        metrics.append({
            "experiment_id": experiment, "observable_id": observable, "unit": unit,
            **{key: selected[0][key] for key in ("source_id", "programme", "time_unit") if key in selected[0]},
            "holdout_axis": selected[0].get("holdout_axis"),
            "n_total": len(selected), "n_scored": len(errors), "coverage": len(errors) / len(selected),
            "rmse": math.hypot(*errors) / math.sqrt(len(errors)) if errors else None,
            "mae": math.fsum(abs(error) / len(errors) for error in errors) if errors else None,
        })
    count = sum(row["scored"] for row in rows)
    return {
        "n_total": len(rows), "n_scored": count,
        "coverage": count / len(rows) if rows else 0.0, "metrics": metrics, "rows": rows,
        "aggregation": "no pooling of different observable units; whole-experiment metrics",
    }


def _programme_scores(score: dict, evaluation_kind: str) -> dict:
    result = {}
    for programme in sorted({row.get("programme", "osmotic_hog") for row in score["rows"]}):
        rows = [row for row in score["rows"] if row.get("programme", "osmotic_hog") == programme]
        count = sum(row["scored"] for row in rows)
        result[programme] = {
            "evaluation_kind": evaluation_kind, "n_total": len(rows), "n_scored": count, "coverage": count / len(rows),
            "metrics": [metric for metric in score["metrics"] if metric.get("programme", "osmotic_hog") == programme],
            "aggregation": "no pooled RMSE across source clocks, observables or assay units",
        }
    return result


def score_synthetic_reference(reference: MechanisticReference, predictions: Sequence[dict]) -> dict:
    if not isinstance(reference, MechanisticReference):
        raise ValueError("synthetic scoring requires generated reference truth")
    truth = reference.evaluator_truth()
    score = _score_rows(truth["targets"], predictions)
    return {
        "evaluation_kind": "synthetic_recovery", "reference_id": reference.frozen.sha256,
        "independent_real_validation": False, "full_biology_claim": False,
        "identifiability": "prediction recovery does not identify the generating parameter vector",
        "programme_scores": _programme_scores(score, "synthetic_recovery"),
        "scope": reference.frozen.to_dict()["scope"], **score,
    }


def score_parameter_recovery(reference: MechanisticReference, estimates: Sequence[dict]) -> dict:
    if not isinstance(reference, MechanisticReference) or not isinstance(estimates, (list, tuple)):
        raise ValueError("parameter recovery needs generated synthetic truth and explicit estimate records")
    expected = {}
    for experiment in reference.evaluator_truth()["experiments"]:
        if experiment["design"]["split"] != "test":
            continue
        provenance = experiment["provenance"]
        for name, value in provenance["reference_parameter_values"].items():
            key = (experiment["design"]["experiment_id"], name)
            expected[key] = {
                "reference_value": value, "unit": provenance["units"]["parameters"][name],
                "identifiability": provenance["identifiability"],
                **{field: experiment[field] for field in ("source_id", "programmes") if field in experiment},
            }
    supplied = {}
    for record in estimates:
        if not isinstance(record, dict) or set(record) != {"experiment_id", "parameter_id", "estimate", "unit"}:
            raise ValueError("estimates require exactly experiment_id, parameter_id, estimate and unit")
        key = (record["experiment_id"], record["parameter_id"])
        if key not in expected or key in supplied:
            raise ValueError("unknown or duplicate source-coordinate parameter estimate")
        if record["unit"] != expected[key]["unit"]:
            raise ValueError("parameter recovery units must match source coordinates, including unspecified units")
        supplied[key] = None if record["estimate"] is None else _number(record["estimate"], "parameter estimate")
    rows = []
    for (experiment_id, parameter_id), truth in expected.items():
        estimate = supplied.get((experiment_id, parameter_id))
        rows.append({
            "experiment_id": experiment_id, "parameter_id": parameter_id, **truth,
            "estimate": estimate, "scored": estimate is not None,
            "point_error": None if estimate is None else estimate - truth["reference_value"],
            "resolved_uncertainty": None,
        })
    count = sum(row["scored"] for row in rows)
    return {
        "reference_id": reference.frozen.sha256, "evaluation_kind": "synthetic_parameter_coordinate_diagnostic",
        "n_total": len(rows), "n_scored": count, "coverage": count / len(rows) if rows else 0.0,
        "rows": rows, "unique_parameter_recovery_established": False,
        "independent_real_validation": False,
        "scope": "point differences in explicitly matched source coordinates only; unresolved directions and missing uncertainty remain unresolved; unlike-unit errors are not pooled",
    }


def run_synthetic_benchmark(reference: MechanisticReference, learner: Callable[[dict], Sequence[dict]]) -> dict:
    if not callable(learner):
        raise ValueError("learner must accept only the learner-input projection")
    predictions = learner(reference.learner_inputs())
    return score_synthetic_reference(reference, predictions)


def training_mean_baseline(inputs: Mapping) -> list[dict]:
    if inputs.get("domain") != "synthetic_native_reference":
        raise ValueError("baseline accepts only synthetic native learner inputs")

    def coordinate(row):
        return row.get("source_id"), row.get("programme"), row["observable_id"], row["unit"]

    groups = {}
    for row in inputs["training_observations"]:
        groups.setdefault(coordinate(row), []).append(row["value"])
    return [
        {"row_id": row["row_id"], "prediction": float(np.mean(groups[coordinate(row)])) if coordinate(row) in groups else None,
         "unit": row["unit"]}
        for row in inputs["prediction_queries"]
    ]


def sensitivity_identifiability(sensitivity, parameter_ids: Sequence[str], parameter_units: Sequence[str], *, relative_tolerance: float = 1e-8) -> dict:
    matrix = np.asarray(sensitivity)
    if matrix.dtype.kind not in "fiu" or matrix.ndim != 2 or not np.isfinite(matrix).all() or not matrix.shape[0]:
        raise ValueError("sensitivity must be a nonempty finite real matrix")
    names = tuple(_name(name, "parameter id") for name in parameter_ids)
    units = tuple(_name(unit, "parameter unit") for unit in parameter_units)
    if not names or len(set(names)) != len(names) or matrix.shape[1] != len(names) or len(units) != len(names):
        raise ValueError("sensitivity columns require unique names and explicit units")
    tolerance = _number(relative_tolerance, "relative rank tolerance")
    if not 0 < tolerance < 1:
        raise ValueError("rank tolerance must be between zero and one")
    _, singular, vectors = np.linalg.svd(matrix.astype(float), full_matrices=True)
    threshold = float(singular[0]) * tolerance if singular.size else 0.0
    rank = int(np.count_nonzero(singular > threshold))
    nullspace = vectors[rank:]
    return {
        "parameter_ids": list(names), "parameter_units": list(units), "rank": rank,
        "singular_values": singular.tolist(), "relative_tolerance": tolerance,
        "unresolved_dimension": len(names) - rank,
        "unresolved_directions": [dict(zip(names, vector.tolist(), strict=True)) for vector in nullspace],
        "covariance": None,
        "scope": "local numerical sensitivity only; retain nullspace, never report a pseudoinverse as resolved uncertainty",
    }


def evaluate_native_reference(
    frozen: FrozenReference, *, evaluation_kind: str,
    inventory: EvidenceInventory | None = None, root=None,
) -> dict:
    if evaluation_kind != "published_fit_replay":
        raise EvidenceGap("independent real validation is not established; only explicitly requested published_fit_replay is available")
    inventory = load_parameter_evidence() if inventory is None else inventory
    verify_reference_freeze(frozen, inventory, root=root)
    from ystwin.analysis.hog_data import load_native_hog_data

    root = repository_root() if root is None else Path(root)
    dataset = load_native_hog_data(source_path(root, "data/hog2013"))
    observations = dataset.observations
    selected = observations.loc[observations.genotype.eq("wild_type") & observations.nacl_molar.eq(0.4)].copy()
    grid = tuple(sorted({0.0, *selected.time_model_s.astype(float)}))
    experiment = ReferenceExperiment("native_wt_replay", "hog2013_wt", grid, (SaltChange(3600, 0.4),), "test", holdout_axis="environment")
    trajectory, provenance = _simulate(experiment, inventory, root=root)
    time_index = {float(time): index for index, time in enumerate(trajectory.times_s)}
    targets, predictions, refused = [], [], []
    for row in selected.to_dict("records"):
        row_id = f"{row['supplement_id']}:{row['source_sheet']}:{row['source_cell']}"
        target = {
            "row_id": row_id, "experiment_id": row["experiment_id"],
            "observable_id": row["observable_id"], "time_s": float(row["time_model_s"]),
            "time_unit": "s", "source_id": "hog2013", "programme": "osmotic_hog",
            "value": float(row["value"]), "unit": row["unit"],
            "source_sheet": row["source_sheet"], "source_cell": row["source_cell"],
            "original_fit_status": row["original_fit_status"],
            "sd": None if not math.isfinite(row["sd"]) else float(row["sd"]),
        }
        targets.append(target)
        if row["observable_id"] in ("glycerol_measured", "glycerol_e"):
            if row["unit"] != "mol/L" or row["measurement_provenance"] != "processed_hplc_measurement":
                raise ValueError("native HPLC observation convention changed")
            predictions.append({
                "row_id": row_id, "unit": row["unit"],
                "prediction": float(trajectory.variables[row["observable_id"]][time_index[target["time_s"]]]),
            })
        else:
            refused.append({"row_id": row_id, "reason": "source protein-native units do not identify a molar or absolute phosphorylation scale; no fitted gain substituted"})
    score = _score_rows(targets, predictions)
    refused_programmes = {
        programme: {"reason": "no real observation/time/assay alignment is supplied by this source-native challenge generator",
                    "blockers": inventory.programme(programme)["blockers"]}
        for contract in frozen.to_dict().get("source_contracts", {}).values()
        for programme in contract["programmes"] if programme != "osmotic_hog"
    }
    return {
        "evaluation_kind": "real_native_published_fit_replay", "reference_id": frozen.sha256,
        "independent_real_validation": False, "synthetic_recovery": False,
        "observation_mapping": "HPLC glycerol source-native mol/L convention only; no global SBML volume rescaling or dry-mass conversion",
        "no_refit": True, "refused_observations": refused, "provenance": provenance,
        "blockers": inventory.programme("osmotic_hog")["blockers"],
        "refused_programmes": refused_programmes, "programme_scores": _programme_scores(score, "real_native_published_fit_replay"),
        **score,
    }
