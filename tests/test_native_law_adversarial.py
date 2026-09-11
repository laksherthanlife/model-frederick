"""Independent caller-seam tests; every numeric fixture here is engineering-only.

No reserved observations, old product calibrations, or synthetic biological
validation claims are used. These tests call the real learners/verifier, not
mocked fitting or scoring implementations.
"""

from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import math
from statistics import NormalDist

import numpy as np
import pytest


def _rate_case():
    """A small identifiable engineering rate panel with explicit raw identities."""
    groups = {f"culture-{i}": "development" if i % 4 == 0 else "train" for i in range(8)}
    records = []
    for i, (group, role) in enumerate(groups.items()):
        for j in range(2):
            enzyme = 0.8 + 0.17 * ((3 * i + j) % 7)
            substrate = 0.2 + 0.8 * i + 0.3 * j
            values = {"E": enzyme, "S": substrate, "v": 3.0 * enzyme * substrate}
            record_id = f"engineering-{i}-{j}"
            records.append({
                "record_id": record_id,
                "group_id": group,
                "role": role,
                "grouping": {"culture": group, "experiment": f"experiment-{i}"},
                "values": values,
                "provenance": {
                    name: {
                        "evidence_type": "synthetic_engineering",
                        "source_id": "adversarial-engineering-not-biology",
                        "measurement_id": f"{record_id}:{name}",
                    }
                    for name in values
                },
            })
    measurements = {
        "schema_version": "native_rate_measurements.v1",
        "evidence_kind": "synthetic_engineering",
        "reaction": {
            "reaction_id": "engineering-reaction",
            "direction": 1,
            "stoichiometry": {"A": -1.0, "B": 1.0},
            "source_id": "adversarial-engineering-not-biology",
        },
        "target": "v",
        "quantities": {
            "E": {"physical_type": "enzyme_amount_per_dry_mass", "unit": "umol/gDW", "entity_id": "engineering-enzyme"},
            "S": {"physical_type": "metabolite_concentration", "unit": "mM", "entity_id": "A"},
            "v": {"physical_type": "specific_flux", "unit": "umol/(gDW*s)", "entity_id": "engineering-reaction"},
        },
        "records": records,
    }
    training = [name for name, role in groups.items() if role == "train"]
    development = [name for name, role in groups.items() if role == "development"]
    protocol = {
        "schema_version": "native_rate_learning_protocol.v1",
        "protocol_id": "adversarial-engineering-not-biology",
        "evidence_kind": "synthetic_engineering",
        "approval": {"status": "synthetic_engineering_only"},
        "permitted_groups": groups,
        "independence_factors": ["culture", "experiment"],
        "selection": {
            "mode": "development_groups", "metric": "equal_group_scaled_mae",
            "tie_tolerance": 1e-12, "refit_policy": "train_only",
        },
        "folds": [
            {
                "fold_id": f"train-fold-{i}",
                "fit_groups": [name for name in training if name not in training[i::3]],
                "validation_groups": training[i::3],
            }
            for i in range(3)
        ] + [{"fold_id": "development", "fit_groups": training, "validation_groups": development}],
        "candidates": [
            {"candidate_id": "enzyme", "family": "enzyme_only", "enzyme": "E"},
            {"candidate_id": "mass", "family": "mass_action", "enzyme": "E", "substrate": "S"},
        ],
        "settings": {
            "max_nfev": 200,
            "profile_points": 9,
            "profile_confidence": 0.95,
            "log_shape_bounds": [-8.0, 8.0],
            "coefficient_bounds": [-4.0, 4.0],
            "start_values": [-2.0, 0.0, 2.0],
            "condition_limit": 1000000.0,
            "constant_target_rtol": 1e-12,
            "profile_floor": 1e-10,
        },
    }
    return measurements, protocol


def _response_experiment(name, dose, value, partition="train", times=(0.0,)):
    from ystwin.analysis.native_response_learning import NativeExperiment, PhysicalHistory, Readout

    return NativeExperiment(
        experiment_id=name,
        biological_replicate_id=f"biological-{name}",
        source_groups=(f"physical-movie-{name}",),
        partition=partition,
        history=PhysicalHistory(
            quantity="solute_concentration", dose_unit="M", time_unit="s",
            start_time=-1.0, initial_dose=0.0, changes=((0.0, dose),),
        ),
        times=times,
        values=(tuple(value for _ in times),),
        readout=Readout(
            observable_id="engineering-ratio", context_id="engineering-assay-only", unit="dimensionless",
            quantity_kind="nuclear_localization_proxy", normalization="raw_max5_over_median",
            normalization_reference="synthetic engineering fixture; never biological evidence",
        ),
    )


def _response_recipe():
    from ystwin.analysis.native_response_learning import SelectionRecipe

    return SelectionRecipe(
        tau_seconds=(2.0,), ridge=(0.0,), max_state_terms=1,
        time_only_degrees=(1,), support_policy="flag",
    )


def test_rate_target_cannot_be_laundered_as_enzyme_via_a_different_column():
    """The same raw measurement is still the outcome after a column rename."""
    from ystwin.analysis.native_rate_learning import ContractError, fit_rate_laws

    measurements, protocol = _rate_case()
    control = fit_rate_laws(measurements, protocol)
    assert control["status"] == "fitted"
    assert control["selection"]["selected_candidate_id"] == "mass"
    leaked = deepcopy(measurements)
    for record in leaked["records"]:
        record["values"]["E"] = record["values"]["v"]
        record["provenance"]["E"] = deepcopy(record["provenance"]["v"])
    try:
        result = fit_rate_laws(leaked, protocol)
    except ContractError:
        return
    assert result["model"] is None, (
        "A phenotype/outcome measurement was reused as protein abundance under a new "
        "column name. Reject shared target/feature source identity, not just feature == target."
    )


def test_response_winner_obeys_registered_group_mae_not_outlier_driven_mse():
    """protocol.json selection requires group MAE; an outlier reverses MSE's winner."""
    from ystwin.analysis.native_response_learning import NativeResponseModel, fit_native_response

    train = [_response_experiment(f"train-{i}", dose, 1.0 + dose)
             for i, dose in enumerate((0.2, 0.4, 0.6))]
    dev = [_response_experiment(f"dev-{i}", 0.5, 1.4 if i < 4 else 2.4, "development")
           for i in range(5)]
    result = fit_native_response(train, dev, recipe=_response_recipe())
    winner = NativeResponseModel.from_dict(result["model"])
    constant = NativeResponseModel.from_dict(result["baselines"]["constant"]["model"])

    def group_mae(model):
        return float(np.mean([
            np.mean(np.abs(np.asarray(model.predict(e.history, e.times)["values"]) - e.values[0]))
            for e in dev
        ]))

    assert group_mae(constant) == pytest.approx(0.2)
    assert group_mae(winner) <= group_mae(constant) + 1e-12, (
        "The selected response loses to its frozen constant on the preregistered "
        "equal-group MAE. Do not substitute MSE/complexity penalties for protocol selection."
    )


def test_response_prediction_prefix_does_not_depend_on_query_suffix_or_clock_origin():
    """A real fitted predictor must have a causal, batch-independent caller seam."""
    from ystwin.analysis.native_response_learning import NativeResponseModel, fit_native_response

    train = [_response_experiment(f"train-{i}", dose, 1.0 + dose, times=(0.0, 1.0, 2.0))
             for i, dose in enumerate((0.2, 0.4, 0.6))]
    dev = [_response_experiment("dev", 0.5, 1.5, "development", times=(0.0, 1.0, 2.0))]
    model = NativeResponseModel.from_dict(fit_native_response(train, dev, recipe=_response_recipe())["model"])
    history = dev[0].history
    prefix = (0.0, 0.5, 1.0)
    expected = model.predict(history, prefix)["values"]
    extended = replace(history, changes=history.changes + ((100.0, 999.0),))
    assert model.predict(extended, prefix + (101.0,))["values"][:len(prefix)] == pytest.approx(expected)
    shift = 600.0
    shifted = replace(history, start_time=history.start_time + shift,
                      changes=tuple((t + shift, u) for t, u in history.changes))
    assert model.predict(shifted, tuple(t + shift for t in prefix))["values"] == pytest.approx(expected)


@pytest.fixture
def engineering_evidence(tmp_path):
    # Reuse only the sibling test's signed-file scaffolding. Every assessment,
    # numerical mutation and assertion below exercises the real public verifier.
    # All authority domains remain 'fixture'; no production authority is minted.
    from test_native_law_validation import Fixture

    case = Fixture(tmp_path)
    baseline = case.assess()
    assert all(gate.state == "pass" for gate in baseline.gates), [
        (gate.name, gate.reasons) for gate in baseline.gates if gate.state != "pass"
    ]
    assert baseline.status == "fixture_only_not_biology"
    assert not baseline.authorized
    return case


def _gate(assessment, name):
    return next(gate for gate in assessment.gates if gate.name == name)


def _refresh_config(case, name):
    """Rebind changed fixture data, so failures cannot be stale-hash artifacts."""
    freeze = case.bodies["freeze"]
    value = ({"library": freeze["library"], "attempts": freeze["attempts"]}
             if name == "selection_log" else freeze[name])
    old = freeze["config_artifacts"][name]
    new = case.plain(name, value)
    freeze["config_artifacts"][name] = new
    freeze["dependencies"][:] = [new if ref == old else ref for ref in freeze["dependencies"]]
    lineage = case.bodies["lineage"]
    lineage["roots"] = [new["sha256"] if value == old["sha256"] else value for value in lineage["roots"]]
    for node in lineage["nodes"]:
        if node["artifact"] == old:
            node["artifact"] = new
        node["parents"] = [new["sha256"] if value == old["sha256"] else value for value in node["parents"]]
    poison = case.bodies["controls"]["normalization_poison"]
    poison["artifacts_before"] = deepcopy(freeze["dependencies"])
    poison["artifacts_after"] = deepcopy(freeze["dependencies"])


def _sync_prediction_controls(case):
    """Keep arithmetic witnesses consistent while testing actual score defects."""
    rows = case.bodies["predictions"]["rows"]
    points = {row["record_id"]: row["value"] for row in rows if row["model_id"] == "candidate"}
    controls = case.bodies["controls"]
    for name in ("absolute_clock_shift", "normalization_poison", "measurement_scale_gauge", "unit_conversion"):
        controls[name]["predictions"] = dict(points)
    controls["physical_consistency"]["rows"] = [
        {key: row[key] for key in ("record_id", "model_id", "value", "last_input_time_s")}
        for row in rows
    ]
    template = {row["record_id"]: row["value"] for row in rows if row["model_id"] == "source_template"}
    controls["source_template_recovery"]["teacher"] = dict(template)
    controls["source_template_recovery"]["recovered"] = dict(template)


@pytest.mark.parametrize("source_id, receipt_kind, field, gate_name", [
    ("native-source", "observations", "raw_value", "source_bound_measurements"),
    ("marker-source", "support", "high", "independent_matched_support"),
])
def test_verifier_binds_numeric_observations_to_referenced_source_not_only_hashes(
    engineering_evidence, source_id, receipt_kind, field, gate_name,
):
    """Authenticated stale target/support projections must be source-bound."""
    case = engineering_evidence
    source = next(source for source in case.bodies["dataset"]["sources"] if source["source_id"] == source_id)
    raw = json.loads((case.root / source["artifact"]["path"]).read_bytes())
    # The fixture source is an actual readable table, not an unreadable dummy ref.
    rid = raw["rows"][0]["record_id"]
    reported = next(row for row in case.bodies[receipt_kind]["rows"] if row["record_id"] == rid)
    assert raw["rows"][0][field] == reported[field]
    raw["rows"][0][field] += 123.0
    if receipt_kind == "observations":
        raw["rows"][0]["numerator"] = raw["rows"][0]["raw_value"] * raw["rows"][0]["denominator"]
    source["artifact"] = case.plain(f"adversarial-{receipt_kind}-source", raw)
    case.seal()
    result = case.assess()
    assert not result.authorized  # fixture scope is never real biological evidence
    assert _gate(result, gate_name).state != "pass", (
        "The signed quantitative projection contradicts its referenced source table, "
        "yet the source-bound gate passes. Reopen typed source locators and compare "
        "values/units/uncertainties, or report missing source-adapter evidence."
    )


def test_verifier_recomputed_reference_hash_does_not_authenticate_fabricated_evidence(engineering_evidence):
    from ystwin.analysis.native_law_validation import SealedArtifact

    case = engineering_evidence
    path = case.root / case.bundle.observations.path
    forged = json.loads(path.read_bytes())
    forged["payload"]["rows"][0]["raw_value"] += 1.0
    forged["payload"].update(passed=True, authorized=True, evidence_status="verified")
    # Retain the old MAC while recomputing all caller-visible identity hashes.
    raw = json.dumps(forged, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    path.write_bytes(raw)
    changed = SealedArtifact(case.bundle.observations.path, hashlib.sha256(raw).hexdigest())
    registry = case.bodies["registry"]
    registry["artifacts"]["observations"] = changed.sha256
    next(event for event in registry["events"] if event["event"] == "evaluation_sealed")["artifacts"]["observations"] = changed.sha256
    head = case.sign("registry", registry)
    case.bundle = replace(case.bundle, observations=changed, registry=head)
    key = f"{registry['protocol_id']}/{registry['candidate_id']}"
    case.trust = replace(case.trust, registry_heads={key: head.sha256})
    result = case.assess()
    assert not result.authorized
    gate = _gate(result, "observations_authentication")
    assert gate.state == "fail" and any("authentication failed" in reason for reason in gate.reasons)


def test_verifier_recomputes_development_winner_despite_excellent_final_scores(engineering_evidence):
    case = engineering_evidence
    attempts = case.bodies["freeze"]["attempts"]
    assert case.bodies["freeze"]["selected_attempt_id"] == attempts[0]["attempt_id"]
    attempts[0]["development_predictions"] = {rid: 4.0 for rid in attempts[0]["development_predictions"]}
    _refresh_config(case, "selection_log")
    case.seal()
    result = case.assess()
    assert _gate(result, "selection_observation_and_scope_freeze").state == "pass"
    assert _gate(result, "prediction_quality").state == "pass"
    gate = _gate(result, "recomputed_development_selection_and_normalization")
    assert gate.state == "fail" and any("development-only selection" in reason for reason in gate.reasons)
    assert not result.authorized


def test_verifier_full_interval_coverage_cannot_replace_point_accuracy(engineering_evidence):
    case = engineering_evidence
    error = case.bodies["freeze"]["error_model"]
    error["prediction_sd"] = 0.6
    _refresh_config(case, "error_model")
    final_ids = {record["record_id"] for record in case.bodies["dataset"]["records"] if record["group_id"] in case.final}
    halfwidth = NormalDist().inv_cdf(0.95) * math.hypot(error["prediction_sd"], error["measurement_sd"])
    for row in case.bodies["predictions"]["rows"]:
        if row["model_id"] == "candidate" and row["record_id"] in final_ids:
            row["value"] = 2.8
        row["prediction_sd"] = error["prediction_sd"]
        row["lower"], row["upper"] = row["value"] - halfwidth, row["value"] + halfwidth
    candidate = [row for row in case.bodies["predictions"]["rows"] if row["model_id"] == "candidate" and row["record_id"] in final_ids]
    assert all(row["lower"] < 2.0 < row["upper"] for row in candidate)
    assert 2 * halfwidth < 2.0  # even the registered width budget is satisfied
    assert np.mean([abs(row["value"] - 2.0) for row in candidate]) == pytest.approx(0.8)
    _sync_prediction_controls(case)
    case.seal()
    result = case.assess()
    assert _gate(result, "complete_presealed_predictions").state == "pass"
    gate = _gate(result, "prediction_quality")
    assert gate.state == "fail" and any("candidate final mean error" in reason for reason in gate.reasons)
    assert not result.authorized


@pytest.mark.parametrize("comparator", ["dependency_ablated", "input_reassigned", "source_template"])
def test_verifier_no_dependency_or_shuffle_advantage_cannot_earn_credit(engineering_evidence, comparator):
    case = engineering_evidence
    rows = case.bodies["predictions"]["rows"]
    candidate = {row["record_id"]: row for row in rows if row["model_id"] == "candidate"}
    for row in rows:
        if row["model_id"] == comparator:
            for field in ("value", "lower", "upper"):
                row[field] = candidate[row["record_id"]][field]
    _sync_prediction_controls(case)
    case.seal()
    result = case.assess()
    assert _gate(result, "negative_controls_gauge_and_physics").state == "pass"
    gate = _gate(result, "prediction_quality")
    assert gate.state == "fail" and any(comparator in reason for reason in gate.reasons)
    assert not result.authorized


@pytest.mark.parametrize("mutation", ["drop", "replace_with_zero"])
def test_verifier_retains_censored_source_rows_without_imputing_accuracy(engineering_evidence, mutation):
    case = engineering_evidence
    record = deepcopy(next(row for row in case.bodies["dataset"]["records"] if row["group_id"] in case.final))
    original_id = record["record_id"]
    record.update(record_id=f"{original_id}:censored", physical_observation_id=f"censored:{original_id}",
                  locator=f"{record['group_id']}/censored-target", status="censored")
    observation = deepcopy(next(row for row in case.bodies["observations"]["rows"] if row["record_id"] == original_id))
    observation.update(record_id=record["record_id"], raw_value=None, raw_sd=None, numerator=None, denominator=None)
    case.bodies["dataset"]["records"].append(record)
    case.bodies["observations"]["rows"].append(observation)
    source = next(row for row in case.bodies["dataset"]["sources"] if row["source_id"] == "native-source")
    raw = json.loads((case.root / source["artifact"]["path"]).read_bytes())
    raw["rows"].append(deepcopy(observation))
    source["artifact"] = case.plain("source-fixture", raw)
    case.seal()
    baseline = case.assess()
    assert all(gate.state == "pass" for gate in baseline.gates), baseline.gates
    if mutation == "drop":
        case.bodies["observations"]["rows"].remove(observation)
    else:
        observation["raw_value"] = 0.0
        observation["raw_sd"] = 0.01
    case.seal()
    result = case.assess()
    gate = _gate(result, "source_bound_measurements")
    assert gate.state == "fail", "Missing/censored source rows cannot vanish or become exact zeros."
    assert not result.authorized
