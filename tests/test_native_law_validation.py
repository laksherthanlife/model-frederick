from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import hashlib
import hmac
import json
import math
from pathlib import Path
from statistics import NormalDist

import pytest

from ystwin.analysis.native_law_validation import (
    Authority,
    EvidenceBundle,
    SealedArtifact,
    ValidationTrust,
    assess_native_law,
)


CANDIDATE = "granados2018_sfp1_carbon_localization"
PROTOCOL = "native-law-v2-preregistered-01"
MODELS = (
    "candidate", "training_constant", "affine_inputs", "clock_only",
    "source_template", "simple_mechanism", "dependency_ablated", "input_reassigned",
)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(value).hexdigest()


class Fixture:
    def __init__(self, root):
        self.root = root
        self.authorities = (
            Authority("fixture-registry", b"r" * 32, frozenset({"registry"}), "fixture"),
            Authority("fixture-custodian", b"d" * 32, frozenset({"dataset", "observations", "support"}), "fixture"),
            Authority("fixture-lineage", b"l" * 32, frozenset({"lineage"}), "fixture"),
            Authority("fixture-runner", b"x" * 32, frozenset({"freeze", "predictions", "controls"}), "fixture"),
        )
        repo = Path(__file__).resolve().parents[1]
        for name in ("protocol", "claim_domain"):
            path = f"data/native_law_v2/{name}.json"
            target = root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((repo / path).read_bytes())
        self.protocol_sha = digest((root / "data/native_law_v2/protocol.json").read_bytes())
        self.domain_sha = digest((root / "data/native_law_v2/claim_domain.json").read_bytes())
        self.bodies = {}
        self.refs = {}
        self.build()
        self.seal()

    def plain(self, name, value):
        path = f"artifacts/{name}.json"
        target = self.root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        raw = canonical(value)
        target.write_bytes(raw)
        return {"path": path, "sha256": digest(raw)}

    def sign(self, kind, body):
        authority = next(a for a in self.authorities if kind in a.kinds)
        envelope = {
            "schema_version": 1, "issuer": authority.issuer, "kind": kind,
            "receipt_id": f"fixture-{kind}", "payload": body,
        }
        envelope["mac"] = hmac.new(authority.secret, canonical(envelope), hashlib.sha256).hexdigest()
        ref = self.plain(f"receipt-{kind}", envelope)
        return SealedArtifact(**ref)

    def build(self):
        context = {
            "species": "Saccharomyces cerevisiae", "strain": "synthetic-test-strain",
            "protein": "Sfp1", "temperature_K": 303.15, "medium": "synthetic-test-medium",
            "platform": "synthetic-test-instrument",
        }
        targets = [f"experiment-{i:02}" for i in range(12)]
        targets.sort(key=lambda g: digest((PROTOCOL + CANDIDATE + g).encode()))
        self.final, self.dev, self.train = targets[:6], targets[6:9], targets[9:]
        self.sham = [f"sham-{i}" for i in range(2)]
        self.support = [f"support-{i}" for i in range(6)]
        groups = []
        for role, ids in (("train", self.train), ("development", self.dev), ("final", self.final),
                          ("sham", self.sham), ("support", self.support)):
            for i, gid in enumerate(ids):
                groups.append({
                    "group_id": gid, "role": role, "biological_ids": [f"culture:{gid}"],
                    "dependency_ids": [f"day:{gid}"], "context": context,
                    "acquisition_clock_s": float(i // 2 * 600),
                    "inputs": {"initial_glucose_percent": 2.0,
                               "final_glucose_percent": 2.0 if role == "sham" else 0.1 + 0.1 * (i % 2),
                               "stimulus_relative_time_s": 60.0},
                })
        input_units = {"initial_glucose_percent": "percent_w_v", "final_glucose_percent": "percent_w_v",
                       "stimulus_relative_time_s": "s"}
        for group in groups:
            group["input_locators"] = {name: {"source_id": "native-source", "locator": f"{group['group_id']}/{name}", "unit": unit}
                                       for name, unit in input_units.items()}
            group["input_schedule"] = [{"time_s": -300.0, "value": 2.0, "quantity": "glucose", "unit": "percent_w_v"},
                                       {"time_s": 0.0, "value": group["inputs"]["final_glucose_percent"],
                                        "quantity": "glucose", "unit": "percent_w_v"}]
        self.groups = groups
        records, observations, support_records, contrasts = [], [], [], []
        for group in groups:
            gid = group["group_id"]
            rid = f"measurement:{gid}"
            if group["role"] == "support":
                support_records.append({
                    "record_id": rid, "group_id": gid, "source_id": "marker-source",
                    "low_locator": f"{gid}/before", "high_locator": f"{gid}/after",
                    "physical_observation_ids": [f"marker-before:{gid}", f"marker-after:{gid}"],
                    "unit": "dimensionless_ratio",
                })
                contrasts.append({"record_id": rid, "low": 1.0, "high": 1.5,
                                  "low_sd": 0.03, "high_sd": 0.04, "covariance": 0.0})
                continue
            value = float(self.train.index(gid) + 1) if gid in self.train else 2.0
            records.append({
                "record_id": rid, "group_id": gid, "source_id": "native-source",
                "locator": f"{gid}/target", "physical_observation_id": f"physical:{gid}",
                "observable": "Sfp1_max5_over_median", "source_unit": "dimensionless_ratio",
                "status": "quantified", "time_s": 60.0,
            })
            observations.append({"record_id": rid, "raw_value": value, "raw_sd": 0.05,
                                 "numerator": value * 100, "denominator": 100.0,
                                 "baseline_raw_value": 2.0 if gid in self.sham else None,
                                 "baseline_raw_sd": 0.05 if gid in self.sham else None})
        source = self.plain("source-fixture", {"origin": "synthetic_software_fixture", "rows": observations,
                                               "groups": groups, "unit": "dimensionless_ratio", "input_units": input_units})
        marker = self.plain("marker-fixture", {"origin": "synthetic_software_fixture", "rows": contrasts, "unit": "dimensionless_ratio"})
        target_decoder = {"kind": "json_records", "table_pointer": "/rows", "key_pointer": "/record_id",
                          "fields": {key: f"/{key}" for key in observations[0]}, "unit_pointer": "/unit", "origin_pointer": "/origin"}
        marker_decoder = {"kind": "json_records", "table_pointer": "/rows", "key_pointer": "/record_id",
                          "fields": {key: f"/{key}" for key in contrasts[0]}, "unit_pointer": "/unit", "origin_pointer": "/origin"}
        input_decoder = {"kind": "json_records", "table_pointer": "/groups", "key_pointer": "/group_id",
                         "fields": {**{key: f"/inputs/{key}" for key in input_units},
                                    **{key: f"/{key}" for key in ("group_id", "biological_ids", "dependency_ids", "context", "acquisition_clock_s", "input_schedule")}},
                         "units_pointer": "/input_units", "origin_pointer": "/origin"}
        method = self.plain("methods-fixture", {"origin": "synthetic_software_fixture", "method": "ratio"})
        adapter = {"kind": "ratio", "factor": 1.0, "offset": 0.0,
                   "source_unit": "dimensionless_ratio", "canonical_unit": "dimensionless_ratio"}
        self.bodies["dataset"] = {
            "candidate_id": CANDIDATE, "source_doi": "10.1073/pnas.1716659115",
            "context": context, "groups": groups, "records": records, "support_records": support_records,
            "input_units": input_units,
            "sources": [
                {"source_id": "native-source", "artifact": source, "methods": method,
                 "kind": "native_imaging", "accession": "synthetic-fixture-native", "assay_id": "imaging",
                 "decoder": target_decoder, "input_decoder": input_decoder},
                {"source_id": "marker-source", "artifact": marker, "methods": method,
                 "kind": "nuclear_marker_imaging", "accession": "synthetic-fixture-marker", "assay_id": "independent-marker",
                 "decoder": marker_decoder},
            ],
            "measurement": {"origin": "experimental_observation", "observable": "Sfp1_max5_over_median",
                            "unit": "dimensionless_ratio", "uncertainty_kind": "standard_uncertainty",
                            "adapter": adapter, "calibration_source_id": "native-source"},
            "physical": {"kind": "localization_ratio", "lower": 0.0, "upper": None,
                         "source_id": "native-source"},
        }
        self.bodies["observations"] = {"candidate_id": CANDIDATE, "rows": observations}
        models = []
        code = self.plain("inference-code-fixture", {"origin": "synthetic_software_fixture", "revision": 1})
        source_model = self.plain("source-model-fixture", {"origin": "synthetic_software_fixture", "source_parameters": {"tau_s": 60.0}})
        source_spec = {"artifact": source_model, "source_reference": "synthetic-fixture-only-not-a-real-citation",
                       "source_locator": "equation-and-parameters", "parameter_policy": "published_fixed_except_observation_gain"}
        dependencies = [code, source_model]
        features = list(groups[0]["inputs"])
        for name in MODELS:
            inputs = ["acquisition_clock_s"] if name == "clock_only" else features
            if name == "training_constant":
                inputs = []
            if name == "dependency_ablated":
                inputs = ["initial_glucose_percent", "stimulus_relative_time_s"]
            form = "adaptive" if name == "candidate" else name
            count = 2 if name == "candidate" else 1
            model = {"model_id": name, "form_id": form, "feature_names": inputs,
                     "parameter_count": count, "parameters": [1.0] * count,
                     "functional_form": [form, *inputs], "implementation_sha256": code["sha256"]}
            if name == "source_template":
                model["source_artifact_sha256"] = source_model["sha256"]
            ref = self.plain(f"model-{name}", model)
            dependencies.append(ref)
            models.append({"model_id": name, "artifact": ref, "form_id": form,
                           "feature_names": inputs, "parameter_count": count, "functional_form": model["functional_form"],
                           "fit_groups": self.train, "optimizer_budget": 100})
        loser = self.plain("model-development-loser", {
            "model_id": "development-loser", "form_id": "linear_response", "feature_names": features,
            "parameter_count": 1, "parameters": [1.0], "implementation_sha256": code["sha256"],
            "functional_form": ["affine_response", *features],
        })
        dependencies.append(loser)
        library = [{"form_id": "adaptive", "feature_names": features, "parameter_count": 2, "functional_form": ["adaptive", *features]},
                   {"form_id": "linear_response", "feature_names": features, "parameter_count": 1, "functional_form": ["affine_response", *features]}]
        model_library = {"library": library, "source_template_specification": source_spec,
                         "comparators": [{key: model[key] for key in ("model_id", "form_id", "feature_names", "parameter_count", "functional_form", "optimizer_budget")}
                                         for model in models if model["model_id"] != "candidate"]}
        attempts = [
            {"attempt_id": "attempt-1", "form_id": "adaptive", "checkpoint": models[0]["artifact"],
             "status": "ok", "seed": 7, "evaluations": 10,
             "development_predictions": {f"measurement:{g}": 2.0 for g in self.dev}},
            {"attempt_id": "attempt-2", "form_id": "linear_response", "checkpoint": loser,
             "status": "ok", "seed": 7, "evaluations": 10,
             "development_predictions": {f"measurement:{g}": 2.3 for g in self.dev}},
        ]
        normalization = {"method": "training_group_iqr", "fit_groups": self.train, "offset": 0.0,
                         "scale": 1.0, "training_group_means": {g: float(i + 1) for i, g in enumerate(self.train)}}
        support_adapter = {"factor": 1.0, "offset": 0.0, "unit": "dimensionless_ratio",
                           "fit_groups": self.train}
        error = {"kind": "gaussian_prediction_interval", "nominal_coverage": 0.9,
                 "fit_groups": self.train + self.dev, "prediction_sd": 0.04,
                 "measurement_sd": 0.05, "support_prediction_sd": 0.03}
        configs = {}
        for name, body in (("normalization", normalization), ("observation_adapter", adapter),
                           ("support_adapter", support_adapter), ("error_model", error),
                           ("model_library", model_library), ("selection_log", {"library": library, "attempts": attempts})):
            configs[name] = self.plain(name, body)
            dependencies.append(configs[name])
        self.bodies["freeze"] = {
            "candidate_id": CANDIDATE, "run_id": "fixture-run-1", "context": context,
            "models": models, "library": library, "attempts": attempts, "selected_attempt_id": "attempt-1",
            "normalization": normalization, "observation_adapter": adapter,
            "support_adapter": support_adapter, "error_model": error, "config_artifacts": configs,
            "code_artifacts": [code], "dependencies": dependencies, "source_template_specification": source_spec,
            "input_domain": {
                "initial_glucose_percent": {"unit": "percent_w_v", "lower": 2.0, "upper": 2.0},
                "final_glucose_percent": {"unit": "percent_w_v", "lower": 0.1, "upper": 2.0},
                "stimulus_relative_time_s": {"unit": "s", "lower": 0.0, "upper": 3600.0},
            },
            "gauge_scope": "observable_only",
            "unidentified_directions": ["latent abundance/activity/observation gain scale"],
        }
        training_source = self.plain("training-release-fixture", {"groups": self.train + self.dev})
        exposed = [v for g in groups if g["role"] in {"train", "development"}
                   for v in [g["group_id"], *g["biological_ids"], *g["dependency_ids"]]]
        nodes = [{"artifact": ref, "role": "code" if ref == code else "model_dependency",
                  "parents": [] if ref == code else [training_source["sha256"]],
                  "exposed_dependency_ids": [], "source_id": "implementation", "locators": [ref["path"]]}
                 for ref in dependencies]
        nodes.append({"artifact": training_source, "role": "measurement_release", "parents": [],
                      "exposed_dependency_ids": exposed, "source_id": "native-source", "locators": ["train/development"]})
        self.bodies["lineage"] = {"candidate_id": CANDIDATE, "roots": [r["sha256"] for r in dependencies], "nodes": nodes}
        rows = []
        halfwidth = NormalDist().inv_cdf(0.95) * math.hypot(error["prediction_sd"], error["measurement_sd"])
        for record in records:
            if record["group_id"] not in self.final + self.sham:
                continue
            for model in MODELS:
                value = 2.02 if model == "candidate" else 3.0
                rows.append({"record_id": record["record_id"], "model_id": model, "value": value,
                             "lower": value - halfwidth, "upper": value + halfwidth,
                             "prediction_sd": 0.04, "measurement_sd": 0.05,
                             "last_input_time_s": 0.0})
        self.bodies["predictions"] = {
            "candidate_id": CANDIDATE, "run_id": "fixture-run-1", "rows": rows,
            "support_predictions": [{"record_id": r["record_id"], "value": 0.48, "prediction_sd": 0.03}
                                    for r in support_records],
        }
        self.bodies["support"] = {"candidate_id": CANDIDATE, "kind": "independent_nuclear_marker_observation", "rows": contrasts}
        candidate_points = {r["record_id"]: r["value"] for r in rows if r["model_id"] == "candidate"}
        final_observations = [r for r in observations if r["record_id"] in candidate_points]
        self.bodies["controls"] = {
            "candidate_id": CANDIDATE, "run_id": "fixture-run-1",
            "input_reassignment": {g: self.final[(i + 1) % len(self.final)] for i, g in enumerate(self.final)},
            "absolute_clock_shift": {
                "seconds": 600.0, "predictions": candidate_points,
                "clocks": [{"group_id": g["group_id"], "before": g["acquisition_clock_s"],
                            "after": g["acquisition_clock_s"] + 600.0,
                            "inputs_before": g["inputs"], "inputs_after": g["inputs"]}
                           for g in groups if g["group_id"] in self.final + self.sham],
            },
            "source_template_recovery": {
                "checkpoint_sha256": next(m["artifact"]["sha256"] for m in models if m["model_id"] == "source_template"),
                "label_origin": "engineering_source_template",
                "teacher": {r["record_id"]: r["value"] for r in rows if r["model_id"] == "source_template"},
                "recovered": {r["record_id"]: r["value"] for r in rows if r["model_id"] == "source_template"},
            },
            "normalization_poison": {
                "factor": 7.0, "predictions": candidate_points,
                "artifacts_before": dependencies, "artifacts_after": dependencies,
                "poisoned_values": {r["record_id"]: r["raw_value"] * 7 for r in final_observations},
            },
            "measurement_scale_gauge": {
                "factor": 2.0, "predictions": candidate_points,
                "rows": [{"record_id": r["record_id"], "raw_value": r["raw_value"], "raw_sd": r["raw_sd"],
                          "numerator": r["numerator"] * 2, "denominator": r["denominator"] * 2,
                          "adapter_factor": 1.0} for r in final_observations],
            },
            "unit_conversion": {
                "factor": 1000.0, "predictions": candidate_points,
                "rows": [{"record_id": r["record_id"], "raw_value": r["raw_value"] * 1000,
                          "raw_sd": r["raw_sd"] * 1000, "adapter_factor": 0.001}
                         for r in final_observations],
            },
            "physical_consistency": {"rows": [{"record_id": r["record_id"], "model_id": r["model_id"],
                                                "value": r["value"], "last_input_time_s": r["last_input_time_s"]}
                                               for r in rows]},
        }

    def seal(self):
        for kind in ("dataset", "lineage", "freeze", "observations", "predictions", "support", "controls"):
            self.refs[kind] = self.sign(kind, self.bodies[kind])
        def ids(role):
            return [g["group_id"] for g in self.groups if g["role"] in role]

        events = [
            {"event": "protocol_registered", "protocol_sha256": self.protocol_sha, "claim_domain_sha256": self.domain_sha},
            {"event": "dataset_sealed", "sha256": self.refs["dataset"].sha256},
            {"event": "split_registered", "partitions": {g["group_id"]: g["role"] for g in self.groups}},
            {"event": "model_library_registered", "sha256": self.bodies["freeze"]["config_artifacts"]["model_library"]["sha256"]},
            {"event": "development_opened", "group_ids": ids({"train", "development"})},
            {"event": "fit_started"},
            *[{"event": "development_attempt", "attempt_id": a["attempt_id"]} for a in self.bodies["freeze"]["attempts"]],
            {"event": "selection_frozen", "sha256": self.refs["freeze"].sha256},
            {"event": "predictions_sealed", "sha256": self.refs["predictions"].sha256},
            {"event": "final_outcomes_opened", "group_ids": ids({"final", "sham"})},
            {"event": "support_outcomes_opened", "group_ids": ids({"support"})},
            {"event": "evaluation_sealed", "artifacts": {k: self.refs[k].sha256 for k in ("observations", "support", "controls")}},
        ]
        for i, event in enumerate(events, 1):
            event["sequence"] = i
        history = []
        for event in events:
            if event["event"] in {"development_opened", "final_outcomes_opened", "support_outcomes_opened"}:
                deps = [v for g in self.groups if g["group_id"] in event["group_ids"]
                        for v in [g["group_id"], *g["biological_ids"], *g["dependency_ids"]]]
                history.append({"protocol_id": PROTOCOL, "candidate_id": CANDIDATE, "run_id": "fixture-run-1",
                                "event_sequence": event["sequence"], "purpose": event["event"],
                                "dependency_ids": deps})
        if "registry" not in self.bodies:
            self.bodies["registry"] = {}
        prior = self.bodies["registry"].get("additional_history", [])
        self.bodies["registry"].update({
            "protocol_id": PROTOCOL, "candidate_id": CANDIDATE, "run_id": "fixture-run-1",
            "protocol_sha256": self.protocol_sha, "claim_domain_sha256": self.domain_sha,
            "artifacts": {k: r.sha256 for k, r in self.refs.items() if k != "registry"},
            "events": events, "access_history": history + prior,
        })
        self.refs["registry"] = self.sign("registry", self.bodies["registry"])
        self.bundle = EvidenceBundle(**self.refs)
        self.trust = ValidationTrust(self.authorities, {f"{PROTOCOL}/{CANDIDATE}": self.refs["registry"].sha256})

    def assess(self):
        return assess_native_law(self.bundle, root=self.root, trust=self.trust)


@pytest.fixture
def evidence(tmp_path):
    return Fixture(tmp_path)


def reasons(result):
    return " ".join(reason for gate in result.gates for reason in gate.reasons)


def test_empty_evidence_reports_missing_without_authorizing(tmp_path):
    result = assess_native_law(EvidenceBundle(), root=tmp_path, trust=ValidationTrust((), {}))
    assert not result.authorized
    assert result.status == "missing_evidence"
    assert "registry" in reasons(result)


def test_complete_software_fixture_exercises_scoped_passing_path_but_is_not_biology(evidence):
    result = evidence.assess()
    assert all(g.state == "pass" for g in result.gates), reasons(result)
    assert result.status == "fixture_only_not_biology"
    assert not result.authorized
    assert result.grade == "independently_validated_scoped_empirical_native_response"
    quality = next(g.metrics for g in result.gates if g.name == "prediction_quality")
    assert quality["comparisons"]["clock_only"]["p_value"] == pytest.approx(1 / 64)
    assert quality["mean_group_scaled_mae"] == pytest.approx(0.02)
    assert "Sfp1" in result.scope


def test_recomputed_artifact_digest_does_not_replace_authority_mac(evidence):
    path = evidence.root / evidence.bundle.dataset.path
    envelope = json.loads(path.read_bytes())
    envelope["payload"]["measurement"]["origin"] = "generated_teacher"
    raw = canonical(envelope)
    path.write_bytes(raw)
    evidence.bundle = replace(evidence.bundle, dataset=SealedArtifact(evidence.bundle.dataset.path, digest(raw)))
    result = evidence.assess()
    assert not result.authorized
    assert "authentication" in reasons(result)


def test_caller_cannot_select_old_registry_head(evidence):
    evidence.trust = replace(evidence.trust, registry_heads={f"{PROTOCOL}/{CANDIDATE}": "0" * 64})
    result = evidence.assess()
    assert "current external registry head" in reasons(result)


def test_tampered_checkpoint_refused_even_with_intact_receipts(evidence):
    ref = evidence.bodies["freeze"]["models"][0]["artifact"]
    (evidence.root / ref["path"]).write_text('{"changed":true}')
    assert "artifact digest" in reasons(evidence.assess())


def test_signed_metadata_pass_flags_do_not_supply_measurements(evidence):
    evidence.bodies["observations"] = {"candidate_id": CANDIDATE, "passed": True}
    evidence.seal()
    result = evidence.assess()
    assert not result.authorized
    assert any(g.state == "missing" for g in result.gates)


def test_whole_biological_dependency_alias_cannot_cross_groups(evidence):
    groups = evidence.bodies["dataset"]["groups"]
    a = next(g for g in groups if g["role"] == "train")
    b = next(g for g in groups if g["role"] == "final")
    b["dependency_ids"].append(a["dependency_ids"][0])
    evidence.seal()
    assert "biological dependency" in reasons(evidence.assess())


def test_duplicate_physical_observation_is_not_an_independent_row(evidence):
    records = evidence.bodies["dataset"]["records"]
    records[1]["physical_observation_id"] = records[0]["physical_observation_id"]
    evidence.seal()
    assert "physical observation" in reasons(evidence.assess())


def test_upstream_fit_exposure_cannot_be_erased_by_later_split(evidence):
    evidence.bodies["lineage"]["nodes"][-1]["exposed_dependency_ids"].append(f"culture:{evidence.final[0]}")
    evidence.seal()
    assert "ancestral" in reasons(evidence.assess())


def test_failed_final_from_another_protocol_consumes_same_biological_holdout(evidence):
    evidence.bodies["registry"]["additional_history"] = [{
        "protocol_id": "older-protocol", "candidate_id": "older-candidate", "run_id": "failed-run",
        "event_sequence": 1, "purpose": "final_outcomes_opened",
        "dependency_ids": [f"culture:{evidence.final[0]}"],
    }]
    evidence.seal()
    assert "consumed" in reasons(evidence.assess())


def test_freeze_must_precede_first_final_access(evidence):
    body = evidence.bodies["registry"]
    a = next(e for e in body["events"] if e["event"] == "selection_frozen")
    b = next(e for e in body["events"] if e["event"] == "final_outcomes_opened")
    a["sequence"], b["sequence"] = b["sequence"], a["sequence"]
    ref = evidence.sign("registry", body)
    evidence.bundle = replace(evidence.bundle, registry=ref)
    evidence.trust = replace(evidence.trust, registry_heads={f"{PROTOCOL}/{CANDIDATE}": ref.sha256})
    assert "sequence" in reasons(evidence.assess()) or "before" in reasons(evidence.assess())


def test_normalization_cannot_use_reserved_outcomes(evidence):
    evidence.bodies["freeze"]["normalization"]["fit_groups"].append(evidence.final[0])
    evidence.seal()
    assert "normalization" in reasons(evidence.assess())


def test_partial_final_predictions_are_not_dropped(evidence):
    evidence.bodies["predictions"]["rows"].pop()
    evidence.seal()
    assert "complete" in reasons(evidence.assess())


def test_good_scores_without_independent_support_do_not_authorize(evidence):
    evidence.bodies["support"]["rows"].pop()
    evidence.seal()
    assert "support" in reasons(evidence.assess())


def test_localization_may_not_be_relabelled_phosphorylation(evidence):
    evidence.bodies["dataset"]["measurement"]["observable"] = "Hog1_phosphorylation"
    evidence.seal()
    assert "observable" in reasons(evidence.assess())


def test_absolute_clock_confounding_is_not_resolved_by_invariance_replay(evidence):
    for group in evidence.groups:
        group["acquisition_clock_s"] = 0.0
    evidence.seal()
    assert "clock" in reasons(evidence.assess())


def test_poisoned_reserved_labels_must_leave_actual_predictions_unchanged(evidence):
    poison = deepcopy(evidence.bodies["controls"]["normalization_poison"])
    poison["predictions"][f"measurement:{evidence.final[0]}"] += 0.1
    evidence.bodies["controls"]["normalization_poison"] = poison
    evidence.seal()
    assert "normalization_poison" in reasons(evidence.assess())


def test_unconverted_uncertainty_fails_unit_replay(evidence):
    evidence.bodies["controls"]["unit_conversion"]["rows"][0]["raw_sd"] /= 1000
    evidence.seal()
    assert "unit_conversion" in reasons(evidence.assess())


def test_broad_feasibility_envelope_is_not_a_calibrated_prediction_interval(evidence):
    for row in evidence.bodies["predictions"]["rows"]:
        if row["model_id"] == "candidate":
            row["lower"], row["upper"] = 0.0, 100.0
    evidence.seal()
    assert "interval" in reasons(evidence.assess())


def test_support_same_signal_is_not_independent_measurement(evidence):
    evidence.bodies["dataset"]["support_records"][0]["source_id"] = "native-source"
    evidence.seal()
    assert "support" in reasons(evidence.assess())


def test_role_separation_is_not_supplied_by_receipt_flags(evidence):
    authorities = tuple(replace(a, secret=b"s" * 32) for a in evidence.authorities)
    evidence.trust = replace(evidence.trust, authorities=authorities)
    assert "independent authority" in reasons(evidence.assess())


def test_no_dataset_arrays_or_model_fits_are_needed_for_missing_evidence(tmp_path):
    result = assess_native_law(EvidenceBundle(), root=tmp_path, trust=ValidationTrust((), {}))
    assert result.status == "missing_evidence"
    assert not list(tmp_path.iterdir())


def test_malformed_typed_reference_returns_refusal_not_an_uncaught_exception(evidence):
    evidence.bundle = replace(evidence.bundle, observations=True)
    result = evidence.assess()
    assert not result.authorized
    assert "sealed artifact reference" in reasons(result)


def test_failed_final_comparison_retains_every_experiment_and_comparator(evidence):
    for row in evidence.bodies["predictions"]["rows"]:
        if row["model_id"] == "clock_only":
            row["value"] -= 0.98
            row["lower"] -= 0.98
            row["upper"] -= 0.98
    evidence.seal()
    result = evidence.assess()
    gate = next(g for g in result.gates if g.name == "prediction_quality")
    assert gate.state == "fail"
    assert set(gate.metrics["comparisons"]) == set(MODELS) - {"candidate"}
    assert all(len(losses) == 6 for losses in gate.metrics["group_losses"].values())
    assert gate.metrics["comparisons"]["clock_only"]["p_value"] == 1.0


def test_an_input_unit_name_is_not_an_authoritative_unit_adapter(evidence):
    evidence.bodies["freeze"]["input_domain"]["final_glucose_percent"]["unit"] = "seconds"
    evidence.seal()
    assert "input units" in reasons(evidence.assess())


def test_no_change_label_does_not_make_a_stimulated_group_a_sham(evidence):
    group = next(g for g in evidence.groups if g["role"] == "sham")
    group["inputs"]["final_glucose_percent"] = 0.1
    evidence.seal()
    assert "sham" in reasons(evidence.assess())


def test_identical_registered_forms_with_different_names_do_not_establish_form_learning(evidence):
    library = evidence.bodies["freeze"]["library"]
    library[1]["functional_form"] = library[0]["functional_form"]
    evidence.seal()
    assert "functional" in reasons(evidence.assess())


def test_source_template_requires_a_fixed_source_parameter_policy(evidence):
    evidence.bodies["freeze"]["source_template_specification"]["parameter_policy"] = "fit_all_native_parameters"
    evidence.seal()
    assert "source template" in reasons(evidence.assess())


def test_signed_unknown_ancestor_cannot_be_ignored(evidence):
    evidence.bodies["lineage"]["nodes"][0]["parents"] = ["f" * 64]
    evidence.seal()
    assert "unknown leaf" in reasons(evidence.assess())


def test_symlink_cannot_replace_an_authoritative_receipt(evidence):
    path = evidence.root / evidence.bundle.observations.path
    target = path.with_name("alternate-receipt.json")
    path.rename(target)
    path.symlink_to(target.name)
    assert "symlink" in reasons(evidence.assess())


def test_duplicate_json_key_is_rejected_before_authentication(evidence):
    ref = evidence.bundle.registry
    path = evidence.root / ref.path
    raw = path.read_bytes().replace(b'"schema_version":1', b'"schema_version":1,"schema_version":1')
    path.write_bytes(raw)
    evidence.bundle = replace(evidence.bundle, registry=SealedArtifact(ref.path, digest(raw)))
    assert "duplicate JSON field" in reasons(evidence.assess())


def test_unsupported_measurement_gain_cannot_hide_a_normalization_leak(evidence):
    evidence.bodies["dataset"]["measurement"]["adapter"]["offset"] = 2.0
    evidence.seal()
    assert "normalization" in reasons(evidence.assess())


def test_nonfinite_prediction_is_not_silently_removed(evidence):
    ref = evidence.bundle.predictions
    path = evidence.root / ref.path
    raw = path.read_bytes().replace(b'"value":2.02', b'"value":NaN', 1)
    path.write_bytes(raw)
    evidence.bundle = replace(evidence.bundle, predictions=SealedArtifact(ref.path, digest(raw)))
    result = evidence.assess()
    assert not result.authorized
    assert "registry" in reasons(result) or "nonfinite" in reasons(result)


def test_actual_source_values_are_recomputed_not_inferred_from_matching_hashes(evidence):
    source = evidence.bodies["dataset"]["sources"][0]
    raw = json.loads((evidence.root / source["artifact"]["path"]).read_bytes())
    raw["rows"][0]["raw_value"] += 100.0
    raw["rows"][0]["numerator"] = raw["rows"][0]["raw_value"] * raw["rows"][0]["denominator"]
    source["artifact"] = evidence.plain("source-fixture", raw)
    evidence.seal()
    result = evidence.assess()
    gate = next(g for g in result.gates if g.name == "source_bound_measurements")
    assert gate.state == "fail"
    assert "source/projection" in reasons(result)


def test_unsupported_source_decoder_is_missing_evidence_not_source_verification(evidence):
    evidence.bodies["dataset"]["sources"][0]["decoder"]["kind"] = "unimplemented_binary_format"
    evidence.seal()
    result = evidence.assess()
    gate = next(g for g in result.gates if g.name == "source_bound_measurements")
    assert gate.state == "missing"
    assert not result.authorized


def test_support_contrast_is_recomputed_from_its_independent_source(evidence):
    evidence.bodies["support"]["rows"][0]["high"] += 0.01
    evidence.seal()
    result = evidence.assess()
    assert "support source/projection" in reasons(result)


def test_source_unit_cannot_be_changed_only_in_a_receipt(evidence):
    source = evidence.bodies["dataset"]["sources"][0]
    raw = json.loads((evidence.root / source["artifact"]["path"]).read_bytes())
    raw["unit"] = "mmol/L"
    source["artifact"] = evidence.plain("source-fixture", raw)
    evidence.seal()
    assert "source/projection unit" in reasons(evidence.assess())


def test_relabelling_authority_domain_does_not_turn_synthetic_source_into_biology(evidence):
    evidence.trust = replace(evidence.trust, authorities=tuple(replace(authority, domain="real") for authority in evidence.authorities))
    result = evidence.assess()
    assert not result.authorized
    assert "source origin" in reasons(result)


def test_source_clock_is_not_replaceable_by_authenticated_metadata(evidence):
    for group in evidence.groups:
        group["acquisition_clock_s"] += 100.0
    evidence.seal()
    result = evidence.assess()
    assert "experiment identity/design source/projection" in reasons(result)
