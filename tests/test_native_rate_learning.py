import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from ystwin.analysis.native_rate_learning import (
    ContractError,
    content_sha256,
    fit_rate_laws,
    predict_rate_law,
)


def engineering_case(family="saturation", group_count=8, repeats=3):
    rng = np.random.default_rng(73021)
    records = []
    groups = {f"culture-{i}": "development" if i % 4 == 0 else "train" for i in range(group_count)}
    substrates = np.geomspace(0.04, 60.0, group_count)
    for i, (group, role) in enumerate(groups.items()):
        for j in range(repeats):
            enzyme = float(rng.uniform(0.6, 1.8))
            substrate = float(substrates[i] * (0.85 + 0.15 * j))
            regulator = float(rng.uniform(0.2, 12.0))
            if family == "mass_action":
                rate = 4.0 * enzyme * substrate
            elif family == "enzyme_only":
                rate = 7.0 * enzyme
            elif family == "activation":
                rate = 7.0 * enzyme * substrate / (2.5 + substrate) * regulator / (1.8 + regulator)
            elif family == "inhibition":
                rate = 7.0 * enzyme * substrate / (2.5 + substrate) * 1.8 / (1.8 + regulator)
            else:
                rate = 7.0 * enzyme * substrate / (2.5 + substrate)
            values = {"E": enzyme, "S": substrate, "R": regulator, "v": rate}
            record_id = f"synthetic-{i}-{j}"
            records.append(
                {
                    "record_id": record_id,
                    "group_id": group,
                    "role": role,
                    "grouping": {"culture": group, "experiment": f"synthetic-experiment-{i}"},
                    "values": values,
                    "provenance": {
                        name: {
                            "evidence_type": "synthetic_engineering",
                            "source_id": "engineering-fixture-not-biological-evidence",
                            "measurement_id": f"{record_id}:{name}",
                            "source_observation_ids": [],
                        }
                        for name in values
                    },
                }
            )
    measurements = {
        "schema_version": "native_rate_measurements.v1",
        "evidence_kind": "synthetic_engineering",
        "reaction": {
            "reaction_id": "synthetic-reaction",
            "direction": 1,
            "stoichiometry": {"synthetic-A": -1.0, "synthetic-B": 1.0},
            "source_id": "engineering-fixture-not-biological-evidence",
        },
        "target": "v",
        "quantities": {
            "E": {
                "physical_type": "enzyme_amount_per_dry_mass",
                "unit": "umol/gDW",
                "entity_id": "synthetic-enzyme",
            },
            "S": {
                "physical_type": "metabolite_concentration",
                "unit": "mM",
                "entity_id": "synthetic-A",
            },
            "R": {
                "physical_type": "metabolite_concentration",
                "unit": "mM",
                "entity_id": "synthetic-regulator",
            },
            "v": {
                "physical_type": "specific_flux",
                "unit": "umol/(gDW*s)",
                "entity_id": "synthetic-reaction",
            },
        },
        "records": records,
    }
    group_names = list(groups)
    folds = []
    for i in range(4):
        validation = group_names[i::4]
        folds.append(
            {
                "fold_id": f"group-fold-{i}",
                "fit_groups": [group for group in group_names if group not in validation],
                "validation_groups": validation,
            }
        )
    protocol = {
        "schema_version": "native_rate_learning_protocol.v1",
        "protocol_id": "synthetic-engineering-protocol",
        "evidence_kind": "synthetic_engineering",
        "approval": {"status": "synthetic_engineering_only"},
        "permitted_groups": groups,
        "independence_factors": ["culture", "experiment"],
        "selection": {
            "mode": "group_cross_validation",
            "metric": "equal_group_scaled_mae",
            "tie_tolerance": 1e-12,
            "refit_policy": "all_permitted_groups",
        },
        "folds": folds,
        "candidates": [
            {"candidate_id": "enzyme", "family": "enzyme_only", "enzyme": "E"},
            {"candidate_id": "saturation", "family": "saturation", "enzyme": "E", "substrate": "S"},
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


def feature_records(measurements):
    target = measurements["target"]
    records = []
    for record in measurements["records"]:
        records.append(
            {
                "record_id": record["record_id"],
                "group_id": record["group_id"],
                "values": {name: value for name, value in record["values"].items() if name != target},
                "provenance": {name: value for name, value in record["provenance"].items() if name != target},
            }
        )
    quantities = {name: spec for name, spec in measurements["quantities"].items() if name != target}
    return records, quantities


def candidate_report(result, candidate_id):
    return next(item for item in result["candidates"] if item["candidate_id"] == candidate_id)


def test_synthetic_saturation_recovery_is_engineering_only():
    measurements, protocol = engineering_case()
    result = fit_rate_laws(measurements, protocol)
    assert result["status"] == "fitted"
    assert result["claim_scope"] == "synthetic_engineering_only"
    assert result["selection"]["selected_candidate_id"] == "saturation"
    assert {item["candidate_id"] for item in result["candidates"]} == {"enzyme", "saturation"}
    parameters = result["model"]["parameters"]
    assert parameters["coefficient"]["value"] == pytest.approx(7.0, rel=2e-5)
    assert parameters["K_substrate"]["value"] == pytest.approx(0.0025, rel=2e-5)
    assert parameters["K_substrate"]["unit"] == "mol/L"
    assert parameters["K_substrate"]["profile"]["status"] == "bounded_on_grid"
    assert candidate_report(result, "saturation")["cv_score"] < 1e-12
    assert candidate_report(result, "enzyme")["cv_score"] > 0.01
    assert result["model"]["given_invariants"]["stoichiometry"] == "provided_not_learned"
    assert len(result["predictions"]["out_of_fold"]) == len(measurements["records"])


def test_json_round_trip_keeps_predictions_parameters_and_stoichiometry():
    measurements, protocol = engineering_case()
    result = fit_rate_laws(measurements, protocol)
    restored = json.loads(json.dumps(result, allow_nan=False))
    records, quantities = feature_records(measurements)
    predictions = predict_rate_law(restored["model"], records, quantities)
    assert predictions == predict_rate_law(result["model"], records, quantities)
    assert result["input_sha256"] == content_sha256(measurements)
    for record, prediction in zip(measurements["records"], predictions, strict=True):
        assert prediction["predicted_rate"] == pytest.approx(record["values"]["v"], rel=2e-5)
        assert prediction["model_sha256"] == result["model"]["model_sha256"]
        assert prediction["stoichiometric_contributions"]["synthetic-A"] == -prediction["predicted_rate"]
        assert prediction["stoichiometric_contributions"]["synthetic-B"] == prediction["predicted_rate"]
        assert prediction["uncertainty"]["scope"] == "conditional_parameter_only_not_calibrated_coverage"


@pytest.mark.parametrize("role", ["test", "reserved", "calibration"])
def test_held_out_roles_are_rejected_before_numeric_values_are_read(role):
    measurements, protocol = engineering_case()
    measurements["records"][0]["role"] = role
    measurements["records"][0]["values"]["v"] = float("nan")
    with pytest.raises(ContractError, match="train/dev"):
        fit_rate_laws(measurements, protocol)


def test_unpermitted_groups_and_mixed_roles_are_rejected():
    measurements, protocol = engineering_case()
    measurements["records"][0]["group_id"] = "unapproved-culture"
    with pytest.raises(ContractError, match="permitted"):
        fit_rate_laws(measurements, protocol)
    measurements, protocol = engineering_case()
    measurements["records"][0]["role"] = "train"
    with pytest.raises(ContractError, match="role"):
        fit_rate_laws(measurements, protocol)


@pytest.mark.parametrize("factor", ["culture", "experiment"])
def test_aliasing_independent_units_across_groups_is_rejected(factor):
    measurements, protocol = engineering_case()
    measurements["records"][3]["grouping"][factor] = measurements["records"][0]["grouping"][factor]
    with pytest.raises(ContractError, match="independence"):
        fit_rate_laws(measurements, protocol)


def test_shared_measurement_provenance_cannot_cross_groups():
    measurements, protocol = engineering_case()
    measurements["records"][3]["provenance"]["E"] = copy.deepcopy(measurements["records"][0]["provenance"]["E"])
    with pytest.raises(ContractError, match="measurement.*groups"):
        fit_rate_laws(measurements, protocol)


def test_crossvalidation_requires_disjoint_complete_whole_group_folds():
    measurements, protocol = engineering_case()
    protocol["folds"][0]["fit_groups"].append(protocol["folds"][0]["validation_groups"][0])
    with pytest.raises(ContractError, match="disjoint"):
        fit_rate_laws(measurements, protocol)
    measurements, protocol = engineering_case()
    protocol["folds"] = protocol["folds"][:-1]
    with pytest.raises(ContractError, match="exactly once"):
        fit_rate_laws(measurements, protocol)


def test_validation_measurements_never_affect_fold_normalization_or_fit():
    measurements, protocol = engineering_case()
    first = fit_rate_laws(measurements, protocol)
    changed = copy.deepcopy(measurements)
    validation = set(protocol["folds"][0]["validation_groups"])
    for record in changed["records"]:
        if record["group_id"] in validation:
            record["values"]["S"] *= 100.0
            record["values"]["E"] *= 10.0
            record["values"]["v"] *= 20.0
    second = fit_rate_laws(changed, protocol)
    a = candidate_report(first, "saturation")["folds"][0]
    b = candidate_report(second, "saturation")["folds"][0]
    assert a["normalization"] == b["normalization"]
    assert a["parameters"] == b["parameters"]
    assert a["fit_groups"] == b["fit_groups"]
    assert set(a["fit_groups"]).isdisjoint(a["validation_groups"])


def test_explicit_physical_unit_conversion_preserves_canonical_parameters():
    measurements, protocol = engineering_case()
    original = fit_rate_laws(measurements, protocol)
    changed = copy.deepcopy(measurements)
    changed["quantities"]["E"]["unit"] = "nmol/gDW"
    changed["quantities"]["S"]["unit"] = "umol/L"
    changed["quantities"]["v"]["unit"] = "mmol/(gDW*h)"
    for record in changed["records"]:
        record["values"]["E"] *= 1000.0
        record["values"]["S"] *= 1000.0
        record["values"]["v"] *= 3.6
    converted = fit_rate_laws(changed, protocol)
    for parameter in ("coefficient", "K_substrate"):
        assert converted["model"]["parameters"][parameter]["value"] == pytest.approx(
            original["model"]["parameters"][parameter]["value"], rel=1e-6
        )
    records, quantities = feature_records(changed)
    predictions = predict_rate_law(original["model"], records, quantities)
    for record, prediction in zip(measurements["records"], predictions, strict=True):
        assert prediction["predicted_rate_canonical"] == pytest.approx(record["values"]["v"] * 1e-6, rel=2e-5)


@pytest.mark.parametrize("unit", ["mM", "unknown-protein-units"])
def test_incompatible_or_undefined_units_are_not_guessed(unit):
    measurements, protocol = engineering_case()
    measurements["quantities"]["E"]["unit"] = unit
    with pytest.raises(ContractError, match="unit"):
        fit_rate_laws(measurements, protocol)


def test_relative_enzyme_scaling_retains_only_an_effective_coefficient():
    measurements, protocol = engineering_case()
    measurements["quantities"]["E"] = {
        "physical_type": "relative_enzyme",
        "unit": "relative",
        "entity_id": "synthetic-enzyme",
        "reference_id": "synthetic-protein-reference-A",
    }
    first = fit_rate_laws(measurements, protocol)
    changed = copy.deepcopy(measurements)
    changed["quantities"]["E"]["reference_id"] = "synthetic-protein-reference-B"
    for record in changed["records"]:
        record["values"]["E"] *= 11.0
    second = fit_rate_laws(changed, protocol)
    a = first["model"]["parameters"]["coefficient"]
    b = second["model"]["parameters"]["coefficient"]
    assert a["interpretation"] == b["interpretation"] == "relative_effective_coefficient"
    assert b["value"] == pytest.approx(a["value"] / 11.0, rel=1e-6)
    assert second["model"]["parameters"]["K_substrate"]["value"] == pytest.approx(0.0025, rel=2e-5)
    records, quantities = feature_records(changed)
    with pytest.raises(ContractError, match="reference"):
        predict_rate_law(first["model"], records, quantities)


def test_relative_metabolite_cannot_produce_an_absolute_half_saturation():
    measurements, protocol = engineering_case()
    measurements["quantities"]["S"] = {
        "physical_type": "relative_metabolite",
        "unit": "relative",
        "entity_id": "synthetic-A",
        "reference_id": "synthetic-metabolite-reference",
    }
    result = fit_rate_laws(measurements, protocol)
    parameter = result["model"]["parameters"]["K_substrate"]
    assert parameter["value"] == pytest.approx(2.5, rel=2e-5)
    assert parameter["interpretation"] == "relative_effective_half_saturation"
    assert parameter["unit"] == "relative"


def test_constant_substrate_refuses_separate_capacity_and_half_saturation():
    measurements, protocol = engineering_case()
    for record in measurements["records"]:
        record["values"]["S"] = 2.0
        record["values"]["v"] = 3.0 * record["values"]["E"]
    result = fit_rate_laws(measurements, protocol)
    assert result["selection"]["selected_candidate_id"] == "enzyme"
    saturation = candidate_report(result, "saturation")
    assert saturation["status"] == "failed"
    assert saturation["model"] is None
    assert saturation["parameters"]["K_substrate"]["value"] is None
    assert saturation["parameters"]["coefficient"]["value"] is None
    assert any(failure["code"] == "unidentifiable" for failure in saturation["failures"])


def test_collinear_regulators_do_not_get_arbitrary_distinct_constants():
    measurements, protocol = engineering_case(group_count=12)
    protocol["candidates"].append(
        {"candidate_id": "activation", "family": "activation", "enzyme": "E", "substrate": "S", "regulator": "R"}
    )
    for record in measurements["records"]:
        enzyme, substrate = record["values"]["E"], record["values"]["S"]
        record["values"]["R"] = substrate
        record["values"]["v"] = 12.0 * enzyme * substrate / (1.0 + substrate) * substrate / (5.0 + substrate)
    result = fit_rate_laws(measurements, protocol)
    activation = candidate_report(result, "activation")
    assert activation["status"] == "failed"
    assert activation["model"] is None
    assert activation["parameters"]["K_regulator"]["value"] is None


def test_constant_target_has_no_spurious_r_squared_or_parameter_recovery():
    measurements, protocol = engineering_case()
    for record in measurements["records"]:
        record["values"]["v"] = 2.0
    result = fit_rate_laws(measurements, protocol)
    assert result["status"] == "refused_constant_target"
    assert result["model"] is None
    assert result["selection"]["selected_candidate_id"] is None
    assert result["baseline"]["r_squared"] is None
    assert all(item["status"] == "not_run_constant_target" for item in result["candidates"])
    json.dumps(result, allow_nan=False)


def test_group_block_target_shuffle_destroys_synthetic_mass_action_recovery():
    measurements, protocol = engineering_case("mass_action", group_count=24, repeats=2)
    protocol["candidates"][1] = {"candidate_id": "mass", "family": "mass_action", "enzyme": "E", "substrate": "S"}
    true_result = fit_rate_laws(measurements, protocol)
    shuffled = copy.deepcopy(measurements)
    permutation = np.random.default_rng(1981).permutation(24)
    for destination, source in enumerate(permutation):
        for slot in range(2):
            shuffled["records"][destination * 2 + slot]["values"]["v"] = measurements["records"][int(source) * 2 + slot]["values"]["v"]
            shuffled["records"][destination * 2 + slot]["provenance"]["v"]["source_id"] = "engineering-group-block-null-shuffle"
    null_result = fit_rate_laws(shuffled, protocol)
    assert candidate_report(true_result, "mass")["cv_score"] < 1e-15
    assert candidate_report(null_result, "mass")["cv_score"] > 0.1
    assert true_result["input_sha256"] != null_result["input_sha256"]
    assert null_result["claim_scope"] == "synthetic_engineering_only"


@pytest.mark.parametrize("family", ["activation", "inhibition"])
def test_regulatory_families_learn_synthetic_constants_without_fixed_regulators(family):
    measurements, protocol = engineering_case(family, group_count=12)
    protocol["candidates"].append(
        {"candidate_id": family, "family": family, "enzyme": "E", "substrate": "S", "regulator": "R"}
    )
    result = fit_rate_laws(measurements, protocol)
    assert result["selection"]["selected_candidate_id"] == family
    assert result["model"]["parameters"]["K_regulator"]["value"] == pytest.approx(0.0018, rel=2e-4)


def test_flexible_comparator_learns_both_coefficient_signs():
    measurements, protocol = engineering_case(group_count=12)
    protocol["candidates"] = [{"candidate_id": "flex", "family": "log_linear", "enzyme": "E", "covariates": ["S", "R"]}]
    for record in measurements["records"]:
        enzyme, substrate, regulator = (record["values"][name] for name in ("E", "S", "R"))
        record["values"]["v"] = float(enzyme * np.exp(0.4 + 0.8 * np.log1p(substrate / 10.0) - 0.7 * np.log1p(regulator / 6.0)))
    result = fit_rate_laws(measurements, protocol)
    assert result["status"] == "fitted"
    assert result["model"]["parameters"]["beta:S"]["value"] > 0
    assert result["model"]["parameters"]["beta:R"]["value"] < 0
    assert result["model"]["parameters"]["beta:S"]["interpretation"] == "predictive_coefficient_not_kinetic_constant"


def test_unknown_and_outcome_using_candidates_are_reported_not_silently_dropped():
    measurements, protocol = engineering_case()
    protocol["candidates"].extend(
        [
            {"candidate_id": "unknown", "family": "unapproved_family", "enzyme": "E"},
            {"candidate_id": "target-leak", "family": "enzyme_only", "enzyme": "v"},
        ]
    )
    result = fit_rate_laws(measurements, protocol)
    for candidate_id in ("unknown", "target-leak"):
        report = candidate_report(result, candidate_id)
        assert report["status"] == "failed"
        assert report["model"] is None
        assert report["failures"]
    assert candidate_report(result, "target-leak")["failures"][0]["code"] == "target_as_feature"


def test_predictions_reject_outcomes_tampering_and_missing_features():
    measurements, protocol = engineering_case()
    result = fit_rate_laws(measurements, protocol)
    with pytest.raises(ContractError, match="feature-only"):
        predict_rate_law(result["model"], measurements["records"], measurements["quantities"])
    records, quantities = feature_records(measurements)
    changed = copy.deepcopy(result["model"])
    changed["parameters"]["coefficient"]["value"] *= 2.0
    with pytest.raises(ContractError, match="digest"):
        predict_rate_law(changed, records, quantities)
    records[0]["values"].pop("S")
    with pytest.raises(ContractError, match="missing"):
        predict_rate_law(result["model"], records, quantities)


def test_declared_direction_and_zero_enzyme_preserve_given_stoichiometric_meaning():
    measurements, protocol = engineering_case()
    measurements["reaction"]["direction"] = -1
    measurements["quantities"]["S"]["entity_id"] = "synthetic-B"
    for record in measurements["records"]:
        record["values"]["v"] *= -1.0
    result = fit_rate_laws(measurements, protocol)
    records, quantities = feature_records(measurements)
    records[0]["values"]["E"] = 0.0
    predictions = predict_rate_law(result["model"], records, quantities)
    assert predictions[0]["predicted_rate"] == 0.0
    assert "E" in predictions[0]["outside_training_range"]
    assert all(item["predicted_rate"] <= 0 for item in predictions)
    assert predictions[1]["stoichiometric_contributions"]["synthetic-A"] > 0
    measurements["reaction"]["direction"] = 1
    with pytest.raises(ContractError, match="direction"):
        fit_rate_laws(measurements, protocol)


def test_unapproved_native_training_is_blocked_before_optimization():
    measurements, protocol = engineering_case()
    measurements["evidence_kind"] = "native_measurement"
    protocol["evidence_kind"] = "native_measurement"
    with pytest.raises(ContractError, match="approval"):
        fit_rate_laws(measurements, protocol)


def test_inputs_are_not_mutated_and_protocol_settings_are_digest_bound():
    measurements, protocol = engineering_case()
    original_measurements, original_protocol = copy.deepcopy(measurements), copy.deepcopy(protocol)
    result = fit_rate_laws(measurements, protocol)
    assert measurements == original_measurements
    assert protocol == original_protocol
    assert result["protocol_sha256"] == content_sha256(protocol)
    assert result["model"]["protocol_sha256"] == result["protocol_sha256"]
    assert result["selection"]["tie_break"] == "fewer_parameters_then_model_id"


def test_optimizer_budget_cannot_be_unbounded():
    measurements, protocol = engineering_case()
    protocol["settings"]["max_nfev"] = 1000000000
    with pytest.raises(ContractError, match="budget"):
        fit_rate_laws(measurements, protocol)


def development_protocol(protocol):
    protocol = copy.deepcopy(protocol)
    training = [group for group, role in protocol["permitted_groups"].items() if role == "train"]
    development = [group for group, role in protocol["permitted_groups"].items() if role == "development"]
    protocol["selection"]["mode"] = "development_groups"
    protocol["selection"]["refit_policy"] = "train_only"
    protocol["folds"] = [
        {"fold_id": f"training-cv-{i}", "fit_groups": [group for group in training if group not in training[i::3]], "validation_groups": training[i::3]}
        for i in range(3)
    ] + [{"fold_id": "development-selection", "fit_groups": training, "validation_groups": development}]
    return protocol


def test_development_only_selection_keeps_all_fits_and_normalizers_train_only():
    measurements, protocol = engineering_case(group_count=12)
    protocol = development_protocol(protocol)
    first = fit_rate_laws(measurements, protocol)
    assert first["selection"]["selected_candidate_id"] == "saturation"
    train = {group for group, role in protocol["permitted_groups"].items() if role == "train"}
    assert set(first["model"]["normalization"]["fit_groups"]) == train
    for report in first["candidates"]:
        for fold in report["folds"]:
            assert set(fold["fit_groups"]) <= train
        development = next(fold for fold in report["folds"] if fold["fold_id"] == "development-selection")
        assert report["selection_score"] == pytest.approx(np.mean([loss["scaled_mae"] for loss in development["group_losses"]]))
    changed = copy.deepcopy(measurements)
    for record in changed["records"]:
        if record["role"] == "development":
            record["values"]["E"] *= 3.0
            record["values"]["v"] *= 7.0
    second = fit_rate_laws(changed, protocol)
    for candidate_id in ("enzyme", "saturation"):
        assert candidate_report(first, candidate_id)["model"] == candidate_report(second, candidate_id)["model"]


def test_development_groups_cannot_enter_training_cv():
    measurements, protocol = engineering_case(group_count=12)
    protocol = development_protocol(protocol)
    dev = next(group for group, role in protocol["permitted_groups"].items() if role == "development")
    protocol["folds"][0]["fit_groups"].append(dev)
    with pytest.raises(ContractError, match="development groups cannot"):
        fit_rate_laws(measurements, protocol)


def test_linear_regime_retains_mass_action_combination_not_unbounded_saturation_constants():
    measurements, protocol = engineering_case("mass_action")
    protocol["candidates"].append({"candidate_id": "mass", "family": "mass_action", "enzyme": "E", "substrate": "S"})
    result = fit_rate_laws(measurements, protocol)
    assert result["selection"]["selected_candidate_id"] == "mass"
    saturation = candidate_report(result, "saturation")
    assert saturation["parameters"]["coefficient"]["value"] is None
    assert saturation["parameters"]["K_substrate"]["value"] is None
    assert result["model"]["parameters"]["coefficient"]["value"] == pytest.approx(4000.0)


def test_identical_laws_use_lexical_tie_break_not_attempt_order():
    measurements, protocol = engineering_case("enzyme_only")
    protocol["candidates"] = [
        {"candidate_id": "z-model", "family": "enzyme_only", "enzyme": "E"},
        {"candidate_id": "a-model", "family": "enzyme_only", "enzyme": "E"},
    ]
    result = fit_rate_laws(measurements, protocol)
    assert result["selection"]["selected_candidate_id"] == "a-model"


def test_untyped_extra_projection_fields_and_unknown_conversion_fields_are_rejected():
    measurements, protocol = engineering_case()
    measurements["untyped_outcomes"] = [1.0]
    with pytest.raises(ContractError, match="projection fields"):
        fit_rate_laws(measurements, protocol)
    measurements, protocol = engineering_case()
    measurements["quantities"]["E"]["custom_conversion"] = 2.0
    with pytest.raises(ContractError, match="quantity specification fields"):
        fit_rate_laws(measurements, protocol)


def test_malformed_candidate_is_a_reported_failure():
    measurements, protocol = engineering_case()
    protocol["candidates"].append({"candidate_id": "malformed", "family": [], "enzyme": "E"})
    result = fit_rate_laws(measurements, protocol)
    assert candidate_report(result, "malformed")["status"] == "failed"
    assert candidate_report(result, "malformed")["failures"][0]["code"] == "unknown_family"


def test_selected_substrate_must_be_consumed_under_provided_stoichiometry():
    measurements, protocol = engineering_case()
    measurements["quantities"]["S"]["entity_id"] = "synthetic-B"
    result = fit_rate_laws(measurements, protocol)
    report = candidate_report(result, "saturation")
    assert report["status"] == "failed"
    assert report["failures"][0]["code"] == "substrate_stoichiometry"
    assert result["selection"]["selected_candidate_id"] == "enzyme"


def test_source_standard_uncertainties_are_converted_but_not_claimed_as_prediction_intervals():
    measurements, protocol = engineering_case()
    for record in measurements["records"]:
        record["standard_uncertainties"] = {"E": 0.1, "S": 0.2, "R": 0.3, "v": 0.4}
    result = fit_rate_laws(measurements, protocol)
    values = result["measurement_uncertainty"]["records"][0]["canonical_standard_uncertainties"]
    assert values["S"] == pytest.approx(0.0002)
    assert values["v"] == pytest.approx(0.4e-6)
    assert result["predictions"]["refit"][0]["uncertainty"]["measurement_plus_prediction_interval"] is None


@pytest.mark.parametrize("uncertainty", [{"untyped": 0.1}, {"E": None}, {"E": -0.1}])
def test_source_uncertainty_mapping_is_explicit_finite_and_nonnegative(uncertainty):
    measurements, protocol = engineering_case()
    measurements["records"][0]["standard_uncertainties"] = uncertainty
    with pytest.raises(ContractError, match="uncertainty"):
        fit_rate_laws(measurements, protocol)


def invoke_engineering_cli(arguments):
    root = Path(__file__).resolve().parents[1]
    return subprocess.run(
        [sys.executable, str(root / "scripts" / "run_native_rate_learning.py"), *arguments],
        cwd=root,
        env={**os.environ, "PYTHONPATH": str(root / "src"), "PYTHONDONTWRITEBYTECODE": "1"},
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )


def test_cli_blocks_unapproved_protocol_before_opening_a_projection(tmp_path):
    protocol_path = tmp_path / "protocol.json"
    protocol_path.write_text(json.dumps({"schema_version": 1, "protocol_id": "synthetic-blocked-protocol", "status": "not_approved"}))
    run = invoke_engineering_cli(["--protocol", str(protocol_path), "--measurements", str(tmp_path / "must-not-open.json")])
    assert run.returncode == 2
    result = json.loads(run.stdout)
    assert result["status"] == "blocked"
    assert "learning_contract" in result["reason"]
    assert "must-not-open" not in result["reason"]


def test_cli_engineering_mode_is_explicit_and_emits_complete_json(tmp_path):
    measurements, protocol = engineering_case()
    measurement_path, protocol_path = tmp_path / "synthetic-measurements.json", tmp_path / "protocol.json"
    measurement_path.write_text(json.dumps(measurements))
    protocol_path.write_text(json.dumps(protocol))
    arguments = ["--protocol", str(protocol_path), "--measurements", str(measurement_path)]
    refused = invoke_engineering_cli(arguments)
    assert refused.returncode == 2
    accepted = invoke_engineering_cli([*arguments, "--engineering-only"])
    assert accepted.returncode == 0, accepted.stderr
    result = json.loads(accepted.stdout)
    assert result["claim_scope"] == "synthetic_engineering_only"
    assert result["execution"]["measurement_file_sha256"]
    assert result["model"]["parameters"]["K_substrate"]["value"] == pytest.approx(0.0025, rel=2e-5)


def measurement_reference(origin):
    return {"source_id": origin["source_id"], "measurement_id": origin["measurement_id"]}


def synthetic_source_node(source, measurement, parents=(), derivation_id=None):
    result = {
        "source_id": source,
        "measurement_id": measurement,
        "evidence_type": "synthetic_engineering",
        "source_observation_ids": list(parents),
    }
    if derivation_id is not None:
        result["derivation_id"] = derivation_id
    return result


def derive_from(record, quantity, parents, inventory=()):
    origin = record["provenance"][quantity]
    origin["source_observation_ids"] = list(parents)
    origin["derivation_id"] = f"{record['record_id']}:{quantity}:derivation-instance"
    origin["ancestry"] = copy.deepcopy(list(inventory))


@pytest.mark.parametrize("feature", ["E", "S", "R"])
@pytest.mark.parametrize("source_slot", [0, 1])
def test_target_feature_ancestry_rejects_direct_alias_in_any_same_group_record(feature, source_slot):
    measurements, protocol = engineering_case()
    protocol["candidates"].append({"candidate_id": "regulatory", "family": "activation", "enzyme": "E", "substrate": "S", "regulator": "R"})
    target_record, feature_record = measurements["records"][source_slot], measurements["records"][0]
    feature_record["values"][feature] = target_record["values"]["v"]
    feature_record["provenance"][feature] = copy.deepcopy(target_record["provenance"]["v"])
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_is_checked_before_values_or_optimization(monkeypatch):
    from ystwin.analysis import native_rate_learning

    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    record["provenance"]["E"] = copy.deepcopy(record["provenance"]["v"])
    record["values"]["v"] = float("nan")

    def forbidden_fit(*args, **kwargs):
        pytest.fail("measurement ancestry must be checked before any fitting")

    monkeypatch.setattr(native_rate_learning, "_fit_one", forbidden_fit)
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_survives_renamed_target_enzyme_and_group_aliases():
    measurements, protocol = engineering_case()
    replacements = {"E": "protein_sensor", "v": "observed_flux"}
    measurements["target"] = "observed_flux"
    for old, new in replacements.items():
        measurements["quantities"][new] = measurements["quantities"].pop(old)
    aliases = {group: f"alias-{group}" for group in protocol["permitted_groups"]}
    protocol["permitted_groups"] = {aliases[group]: role for group, role in protocol["permitted_groups"].items()}
    for fold in protocol["folds"]:
        for field in ("fit_groups", "validation_groups"):
            fold[field] = [aliases[group] for group in fold[field]]
    for candidate in protocol["candidates"]:
        candidate["enzyme"] = "protein_sensor"
    for record in measurements["records"]:
        record["group_id"] = aliases[record["group_id"]]
        for field in ("values", "provenance"):
            for old, new in replacements.items():
                record[field][new] = record[field].pop(old)
        record["values"]["protein_sensor"] = record["values"]["observed_flux"]
        record["provenance"]["protein_sensor"] = copy.deepcopy(record["provenance"]["observed_flux"])
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


@pytest.mark.parametrize("identity_kind", ["source_location", "physical_measurement"])
def test_target_feature_ancestry_uses_physical_identity_across_source_aliases(identity_kind):
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    record["provenance"]["E"]["source_id"] = "synthetic-renamed-export"
    record["provenance"]["E"]["measurement_id"] = "synthetic-renamed-reading"
    for name in ("E", "v"):
        if identity_kind == "source_location":
            record["provenance"][name].update({"source_artifact_sha256": "1" * 64, "source_locator": "synthetic-sheet:cell-7"})
        else:
            record["provenance"][name]["physical_measurement_id"] = "synthetic-physical-observation-7"
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_resolves_transitive_physical_aliases():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    record["provenance"]["v"]["physical_measurement_id"] = "synthetic-global-observation"
    location = {"source_artifact_sha256": "2" * 64, "source_locator": "synthetic-cell-13"}
    record["provenance"]["E"].update(location)
    bridge = synthetic_source_node("synthetic-identity-registry", "alias-bridge")
    bridge.update(location, physical_measurement_id="synthetic-global-observation")
    record["provenance"]["E"]["ancestry"] = [bridge]
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_rejects_a_cycle_hidden_by_physical_aliases():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    source = synthetic_source_node("synthetic-alias-source", "aliased-input")
    source["physical_measurement_id"] = "synthetic-same-input-output"
    record["provenance"]["E"]["physical_measurement_id"] = "synthetic-same-input-output"
    derive_from(record, "E", [measurement_reference(source)], [source])
    with pytest.raises(ContractError, match="ancestry.*cycle"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_resolves_legacy_local_parent_references():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    derive_from(record, "E", [record["provenance"]["v"]["measurement_id"]])
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_resolves_indirect_cross_source_derivations():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    intermediate = synthetic_source_node(
        "synthetic-intermediate-export", "intermediate-0", [measurement_reference(record["provenance"]["v"])], "synthetic-intermediate-event"
    )
    derive_from(record, "E", [measurement_reference(intermediate)], [intermediate])
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_is_symmetric_when_the_target_is_derived_from_a_feature():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    derive_from(record, "v", [measurement_reference(record["provenance"]["E"])])
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_rejects_a_shared_raw_ancestor_of_distinct_derivations():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    shared = synthetic_source_node("synthetic-raw-source", "shared-detector-reading")
    for name in ("E", "v"):
        derive_from(record, name, [measurement_reference(shared)], [shared])
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_rejects_a_shared_derivation_instance():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    for name in ("E", "v"):
        raw = synthetic_source_node("synthetic-raw-source", f"reading-{name}")
        derive_from(record, name, [measurement_reference(raw)], [raw])
        record["provenance"][name]["derivation_id"] = "same-synthetic-output-event"
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_allows_independent_identical_numeric_measurements():
    measurements, protocol = engineering_case()
    protocol["candidates"] = [{"candidate_id": "enzyme", "family": "enzyme_only", "enzyme": "E"}]
    for record in measurements["records"]:
        record["values"]["E"] = record["values"]["v"]
    result = fit_rate_laws(measurements, protocol)
    assert result["status"] == "fitted"
    assert result["selection"]["selected_candidate_id"] == "enzyme"
    assert result["measurement_ancestry"]["resolution"] == "declared_graph_closed"
    assert result["measurement_ancestry"]["authority"] == "declared_provenance_not_external_source_authorization"


def test_target_feature_ancestry_allows_disjoint_derivations_using_the_same_method():
    measurements, protocol = engineering_case()
    protocol["candidates"] = [{"candidate_id": "enzyme", "family": "enzyme_only", "enzyme": "E"}]
    for record in measurements["records"]:
        record["values"]["E"] = record["values"]["v"]
        for name in ("E", "v"):
            raw = synthetic_source_node("synthetic-instrument", f"{record['record_id']}:{name}:raw")
            derive_from(record, name, [measurement_reference(raw)], [raw])
            record["provenance"][name]["derivation_method"] = "same-synthetic-calibration-method"
    result = fit_rate_laws(measurements, protocol)
    assert result["status"] == "fitted"
    assert result["measurement_ancestry"]["resolution"] == "declared_graph_closed"


@pytest.mark.parametrize("missing", ["dangling_reference", "undeclared_ancestor", "missing_parent_list", "unknown_ancestry_shape"])
def test_target_feature_ancestry_does_not_upgrade_unknown_derived_ancestry(missing):
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    root = synthetic_source_node("synthetic-raw-source", "not-in-the-projection")
    derive_from(record, "E", [measurement_reference(root)], [root])
    if missing == "dangling_reference":
        record["provenance"]["E"]["ancestry"] = []
    elif missing == "undeclared_ancestor":
        record["provenance"]["E"]["ancestry"][0].pop("source_observation_ids")
    elif missing == "missing_parent_list":
        record["provenance"]["E"].pop("source_observation_ids")
    else:
        record["provenance"]["E"]["ancestry"] = {"status": "unknown"}
    with pytest.raises(ContractError, match="ancestry"):
        fit_rate_laws(measurements, protocol)


@pytest.mark.parametrize("defect", ["unknown", "cycle"])
def test_target_feature_ancestry_inventory_cannot_hide_unresolved_or_cyclic_nodes(defect):
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    first = synthetic_source_node("synthetic-inventory", "first")
    if defect == "unknown":
        first.pop("source_observation_ids")
        inventory = [first]
    else:
        second = synthetic_source_node("synthetic-inventory", "second", [measurement_reference(first)], "second-event")
        first["source_observation_ids"] = [measurement_reference(second)]
        first["derivation_id"] = "first-event"
        inventory = [first, second]
    record["provenance"]["E"]["ancestry"] = inventory
    with pytest.raises(ContractError, match="ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_rejects_cycles():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    intermediate = synthetic_source_node("synthetic-intermediate", "cycle", [measurement_reference(record["provenance"]["E"])], "synthetic-cycle-event")
    derive_from(record, "E", [measurement_reference(intermediate)], [intermediate])
    with pytest.raises(ContractError, match="ancestry.*cycle"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_rejects_inconsistent_definitions_of_a_node():
    measurements, protocol = engineering_case()
    record = measurements["records"][0]
    raw = synthetic_source_node("synthetic-source", "conflicting-observation")
    conflicting = copy.deepcopy(raw)
    conflicting["source_observation_ids"] = [{"source_id": "synthetic-source", "measurement_id": "other"}]
    conflicting["derivation_id"] = "different-event"
    derive_from(record, "E", [measurement_reference(raw)], [raw, conflicting])
    with pytest.raises(ContractError, match="conflicting.*ancestry"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_preserves_whole_group_ancestor_separation():
    measurements, protocol = engineering_case()
    shared = synthetic_source_node("synthetic-source", "shared-enzyme-assay")
    for record in (measurements["records"][0], measurements["records"][3]):
        derive_from(record, "E", [measurement_reference(shared)], [shared])
    with pytest.raises(ContractError, match="measurement.*groups"):
        fit_rate_laws(measurements, protocol)


def test_target_feature_ancestry_legacy_synthetic_roots_are_not_claimed_closed():
    measurements, protocol = engineering_case()
    for record in measurements["records"]:
        for origin in record["provenance"].values():
            origin.pop("source_observation_ids")
    result = fit_rate_laws(measurements, protocol)
    assert result["status"] == "fitted"
    assert result["measurement_ancestry"]["resolution"] == "synthetic_identity_only_unresolved"
    assert result["measurement_ancestry"]["unresolved_synthetic_roots"]
    assert result["claim_scope"] == "synthetic_engineering_only"


@pytest.mark.parametrize("indirect", [False, True])
def test_target_feature_ancestry_is_enforced_again_on_prediction_inputs(indirect):
    measurements, protocol = engineering_case()
    fitted = fit_rate_laws(measurements, protocol)
    records, quantities = feature_records(measurements)
    target = copy.deepcopy(measurements["records"][0]["provenance"]["v"])
    if indirect:
        derive_from(records[0], "E", [measurement_reference(target)], [target])
    else:
        records[0]["provenance"]["E"] = target
    with pytest.raises(ContractError, match="target.*feature.*ancestry"):
        predict_rate_law(fitted["model"], records, quantities)


def test_target_feature_ancestry_cannot_be_dropped_from_a_reserialized_model():
    measurements, protocol = engineering_case()
    fitted = fit_rate_laws(measurements, protocol)
    model = copy.deepcopy(fitted["model"])
    model.pop("measurement_ancestry")
    model["model_sha256"] = content_sha256({key: value for key, value in model.items() if key != "model_sha256"})
    records, quantities = feature_records(measurements)
    with pytest.raises(ContractError, match="measurement ancestry"):
        predict_rate_law(model, records, quantities)
