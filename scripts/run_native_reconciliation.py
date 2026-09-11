from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import gzip
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import sys
import time

import cobra
from cobra.core.formula import Formula, elements_and_molecular_weights

from ystwin.analysis.training_freeze import ArtifactRef, InputSpec, TrainingContract, validate_checkpoint
from ystwin.fba.native_obligations import NativeCarbonObligation, NativeObligationRoles, impose_native_obligation
from ystwin.fba.native_reconciliation import (
    EnergyRoles,
    FluxObservable,
    biomass_flux_equivalence,
    bound_metadata,
    condition_constraints,
    ec_native_observables,
    energy_diagnostic,
    minimum_relaxation,
    native_conditions,
    solve_native,
)
from ystwin.fba.solver import configure


FROZEN_CHECKPOINT = "outputs/native_training_run_02/checkpoint.json"
FROZEN_SHA256 = "384175e19b2451f542ab08e347ee9d9021e22772fbc205e0c404bbc28416ec94"
MODEL_PATH = "data/gem/ecYeastGEM_batch.xml.gz"
ENERGY_ROLES = EnergyRoles("r_4046", "r_4041", "s_0434[c]", "s_0803[c]", "s_0394[c]", "s_0794[c]", "s_1322[c]")
OBLIGATION_ROLES = NativeObligationRoles("s_0765[c]", "s_0565[e]", ("r_1714", "r_1714_REV"),
                                         ("r_1808",), ("prot_pool[c]",))


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path, payload):
    data = json.dumps(payload, sort_keys=True, indent=2, allow_nan=False).encode() + b"\n"
    with path.open("xb") as handle:
        handle.write(gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)


def _frozen_native_context(root, *, artifact_source="original", manifest_path=None, expected_sha256=None):
    root = Path(root).resolve(strict=True)
    if artifact_source == "portable":
        from ystwin.analysis.portable_replay import load_portable_native_checkpoint

        payload, receipt = load_portable_native_checkpoint(
            root, manifest_path=manifest_path, expected_sha256=expected_sha256)
        config = payload["model"]["config"]
        if config["host_model_path"] != MODEL_PATH:
            raise ValueError("native checkpoint host differs from the declared experiment host")
        return config["native_obligation_scenarios"], {
            **receipt, "path": receipt["checkpoint_reference"]["path"],
            "sha256": receipt["checkpoint_reference"]["sha256"],
            "all_frozen_inputs_and_code_verified": True,
            "original_checkpoint_bytes_verified": False,
            "use": "Only preserved native obligation endpoints are reused; public-export integrity is checked without fitting, outcome-based selection or an original-checkpoint-byte claim",
        }
    if artifact_source != "original":
        raise ValueError("artifact_source must be 'original' or 'portable'")
    if manifest_path is not None or expected_sha256 is not None:
        raise ValueError("portable manifest arguments require artifact_source='portable'")
    path = root / FROZEN_CHECKPOINT
    if _digest(path) != FROZEN_SHA256:
        raise ValueError("frozen native checkpoint digest mismatch")
    payload = json.loads(path.read_text())
    manifest = payload["training_manifest"]
    contract = TrainingContract(root, *(
        tuple(InputSpec(e["path"], e["sha256"], e["role"], e["source"]) for e in manifest[section])
        for section in ("inputs", "code")
    ))
    verified = validate_checkpoint(ArtifactRef(path, FROZEN_SHA256), contract=contract)
    config = verified["model"]["config"]
    if config["host_model_path"] != MODEL_PATH:
        raise ValueError("native checkpoint host differs from the declared experiment host")
    return config["native_obligation_scenarios"], {
        "path": FROZEN_CHECKPOINT, "sha256": FROZEN_SHA256,
        "all_frozen_inputs_and_code_verified": True,
        "input_count": len(manifest["inputs"]), "code_count": len(manifest["code"]),
        "use": "Only frozen native obligation endpoints are reused; no fitting and no outcome or prediction artifacts are read",
    }


def _model_identity(model):
    structure = {
        "reactions": {r.id: {"bounds": [bound_metadata(v) for v in r.bounds],
                             "stoichiometry": {m.id: float(v) for m, v in r.metabolites.items()}}
                      for r in model.reactions},
        "constraints": {c.name: [bound_metadata(c.lb), bound_metadata(c.ub), str(c.expression)]
                        for c in model.constraints},
        "objective": [str(model.objective.expression), model.objective.direction],
    }
    return hashlib.sha256(json.dumps(structure, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _readouts(result, condition, observables):
    values = {key: None for key in condition["observations"]}
    if result["status"] == "optimal":
        solved = {key: math.fsum(c * result["fluxes"][rid] for rid, c in observable.coefficients.items())
                  for key, observable in observables.items()}
        values.update({key: value for key, value in solved.items() if key in values})
        glucose = solved["GlucoseUptake"]
        if glucose > 0:
            values["BiomassYield_g_per_g"] = solved["growth"] / (glucose * Formula("C6H12O6").weight / 1000.0)
            carbon = math.fsum(solved[key] * coefficient for key, coefficient in (
                ("CO2production", 1.0), ("Ethanol", 2.0), ("Glycerol", 3.0),
                ("Acetate", 2.0), ("Pyruvate", 3.0),
            ))
            biomass_carbon = 1000.0 * 0.48 * solved["growth"] / elements_and_molecular_weights["C"]
            values["CarbonRecovery_pct"] = 100.0 * (carbon + biomass_carbon) / (6.0 * glucose)
    return {
        key: {"value": values[key], "observed_value": record["value"],
              "reported_value": record["reported_value"], "unit": record["unit"],
              "observation_status": record["observation_status"],
              "residual": (values[key] - record["value"] if values[key] is not None and record["value"] is not None else None),
              "interpretation": ("derived diagnostic, not an additional independent constraint" if key in ("BiomassYield_g_per_g", "CarbonRecovery_pct")
                                 else "a feasible witness, not identified realized flux allocation")}
        for key, record in condition["observations"].items()
    }


def _native_structure(model):
    def reaction_record(reaction):
        return {"id": reaction.id, "bounds": [bound_metadata(v) for v in reaction.bounds],
                "stoichiometry": {m.id: float(c) for m, c in reaction.metabolites.items()},
                "species": {m.id: {"formula": m.formula, "annotation": m.annotation, "charge": m.charge}
                            for m in reaction.metabolites}}

    supply = []
    for reaction in model.boundary:
        if any((c > 0 and reaction.upper_bound > 0) or (c < 0 and reaction.lower_bound < 0)
               for c in reaction.metabolites.values()):
            supply.append(reaction_record(reaction))
    equivalence = biomass_flux_equivalence(model, growth_reaction_id="r_2111",
                                           assembly_reaction_id="r_4041", biomass_metabolite_id="s_0450[c]")
    if equivalence["assembly_flux_per_growth_flux"] != 1.0:
        raise ValueError("native biomass assembly and growth mass normalization differ")
    return {
        "model_id": model.id, "n_reactions": len(model.reactions), "n_metabolites": len(model.metabolites),
        "positive_lower_bounds": [reaction_record(r) for r in model.reactions if r.lower_bound > 0],
        "energy_roles": asdict(ENERGY_ROLES),
        "energy_reactions": [reaction_record(model.reactions.get_by_id(rid)) for rid in ("r_4041", "r_4046")],
        "active_boundary_supplies": supply,
        "protein_pool_unit": "g protein/gDW; native GECKO2 allocation resource, not a chemical metabolite flux",
        "biomass_assembly_equals_growth": True, "biomass_equivalence_evidence": equivalence,
        "strain_assignment": "The native SBML is a generic EC host, not a demonstrated DS28911-specific reconstruction",
        "medium_comparison": {
            "source": "van Hoek 1998 Materials and Methods, page 4227",
            "observed": "DS28911, defined mineral medium with vitamins, 30 C, pH 5.0; glucose feed 7.5 g/L",
            "gem": "Only glucose supplies chemical carbon in the native SBML; inorganic imports and protein-pool supply are listed explicitly",
            "unresolved": "No measured DS28911 biomass-composition/protein-budget mapping; no medium or dry-mass correction is inferred or applied",
        },
    }


def _energy_allowance(host, condition, observables, timeout_s):
    original_gam = -host.reactions.get_by_id("r_4041").get_coefficient("s_0434[c]")
    original_ngam = host.reactions.get_by_id("r_4046").lower_bound
    growth = condition["growth_rate_per_h"]
    original_total = original_gam * growth + original_ngam
    model, provenance = energy_diagnostic(
        host, roles=ENERGY_ROLES, gam_mmol_atp_per_gdw=0.0,
        source="Stoichiometric energy-budget diagnostic, bounded by the original total ATP demand; no fitted parameter",
    )
    model.reactions.get_by_id("r_4046").bounds = (0.0, original_total)
    provenance["ngam_bounds_mmol_atp_per_gdw_h"] = [0.0, original_total]
    objective = FluxObservable("total_maintenance", "mmol ATP/gDW/h", {"r_4046": 1.0},
                               "Complete native ATP hydrolysis vector used as a bounded energy-demand allowance")
    constraints = condition_constraints(condition, observables, mode="uptake_caps", impose_growth=True)
    result = solve_native(model, constraints=constraints, objective=objective, timeout_s=timeout_s)
    allowance = result["objective_value"] if result["status"] == "optimal" else None
    ngam_limit = allowance - original_gam * growth if allowance is not None else None
    gam_limit = (allowance - original_ngam) / growth if allowance is not None else None
    result["energy_diagnostic"] = {
        "provenance": provenance, "original_total_mmol_atp_per_gdw_h": original_total,
        "compatible_total_within_original_budget_mmol_atp_per_gdw_h": allowance,
        "minimum_original_energy_reduction_mmol_atp_per_gdw_h": original_total - allowance if allowance is not None else None,
        "compatible_ngam_at_original_gam_mmol_atp_per_gdw_h": ngam_limit if ngam_limit is not None and ngam_limit >= 0 else None,
        "compatible_gam_at_original_ngam_mmol_atp_per_gdw": gam_limit if gam_limit is not None and gam_limit >= 0 else None,
        "conditional_necessary_upper_bound": "GAM*growth + NGAM <= compatible total, with GAM>=0 and NGAM>=0. The lower compatible total was not minimized; sufficiency throughout this half-plane is not claimed.",
        "parameter_semantics": "GAM here denotes the total biomass-assembly ATP hydrolysis coefficient. Historical GECKO scaleBioMass adds a polymerization component to its fitted GAM; the total is not a separately measured physiological maintenance coefficient.",
        "scope": "Original nutrient capacity observations and all other native constraints; not simultaneous agreement with every measured exchange",
        "calibration": "Not adopted. Feasibility bounds are not identified measurements of either GAM or NGAM.",
    }
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Native-only GEM feasibility and bounded discrepancy diagnostics; no parameter fitting")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout-s", type=int, default=15, choices=range(1, 31))
    parser.add_argument("--growth-rates", type=float, nargs="+")
    parser.add_argument("--skip-obligations", action="store_true")
    parser.add_argument("--artifact-source", choices=("original", "portable"), default="original")
    parser.add_argument("--portable-manifest", type=Path)
    parser.add_argument("--portable-manifest-sha256")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    outputs = (root / "outputs").resolve()
    relative = output.relative_to(outputs) if output.is_relative_to(outputs) else Path()
    if not relative.parts or not relative.parts[0].startswith("native_reconciliation_"):
        raise ValueError("experiment output must be within a fresh outputs/native_reconciliation_* directory")
    if not output.parent.is_dir():
        raise ValueError("experiment output parent must already exist")
    if output.exists():
        raise ValueError("refusing to overwrite an existing experiment directory")
    conditions = native_conditions()
    if args.growth_rates:
        if not set(args.growth_rates) <= {c["growth_rate_per_h"] for c in conditions}:
            raise ValueError("only verified native Table 1 growth rates can be selected")
        conditions = [c for c in conditions if c["growth_rate_per_h"] in args.growth_rates]
    artifact_options = {"artifact_source": args.artifact_source, "manifest_path": args.portable_manifest,
                        "expected_sha256": args.portable_manifest_sha256}
    obligations, frozen = _frozen_native_context(root, **artifact_options)
    model_sha256 = _digest(root / MODEL_PATH)
    host = cobra.io.read_sbml_model(str(root / MODEL_PATH))
    configure(host)
    host_identity = _model_identity(host)
    observables = ec_native_observables(host)
    structure = _native_structure(host)
    output.mkdir()
    reports = output / "cases"
    reports.mkdir()
    _json(output / "native_conditions.json", conditions)
    _json(output / "native_structure.json", structure)
    _json(output / "protocol.json", {
        "version": "native_reconciliation_diagnostic_v1", "calibrated": False,
        "frozen_checkpoint": frozen, "host_path": MODEL_PATH, "host_sha256": model_sha256,
        "native_source_precision": conditions[0]["precision_provenance"],
        "environment": {"python": sys.version, "packages": {name: importlib.metadata.version(name)
                        for name in ("cobra", "numpy", "scipy", "pandas", "swiglpk", "python-libsbml")}},
        "implementation_sha256": {path: _digest(root / path) for path in (
            "src/ystwin/fba/native_reconciliation.py", "scripts/run_native_reconciliation.py",
            "tests/test_native_reconciliation.py")},
        "per_solve_timeout_s": args.timeout_s,
        "maximum_case_count": 11 * len(conditions) + 9,
        "infeasibility_policy": "Keep failed/infeasible cases with null physical results; no zero filling, original-bound relaxation, or solver-policy change",
        "rounding_policy": "Primary printed half-digits under nearest rounding, not confidence intervals; source zeros unquantified and never fixed to zero",
        "energy_policy": "Only explicit counterfactual copies may change energy terms; no production parameter is calibrated or adopted",
        "obligation_policy": "All frozen native dose/window endpoints, without selecting a favorable one; no HOG mol/L-to-gDW conversion",
        "data_split": "All ten native conditions were already development evidence, including interpolation checks; none is a fresh blind validation set",
        "independent_validation": {
            "used": False,
            "reason": "No parameter calibration is adopted; no final external holdout was opened",
            "candidate_requirements": "Independently generated native chemostat fluxes, strain-specific dry biomass/protein composition, and maintenance evidence with an audited GEM/GECKO fitting lineage",
            "overlap_audit": "van Hoek Table 1 is already used in native training. HOG2013 s002 overlaps the published/local kinetic fit. Historical GECKO commit adac1095a670848036d9e62bc85a0d048c501ec0 fits GAM using the same van Hoek rows at D=0.025,0.05,0.10,0.15,0.20; these are not independent prior-validation rows. The exact historical EC SBML build dependency SHA remains unverified.",
        },
    })
    rows = []

    def record(name, condition, result, metadata=None):
        growth = condition["growth_rate_per_h"]
        case_id = f"D_{growth:g}__{name}"
        result["case_id"] = case_id
        result["condition_id"] = condition["condition_id"]
        result["case_metadata"] = metadata or {}
        result["native_readouts"] = _readouts(result, condition, observables)
        actual = None
        if result["status"] == "optimal":
            actual = math.fsum(c * result["fluxes"][rid] for rid, c in observables["growth"].coefficients.items())
        energy = result.get("energy_diagnostic", {})
        protein = result["fluxes"]["prot_pool_exchange"] if result["fluxes"] is not None else None
        rows.append({
            "case_id": case_id, "growth_per_h": growth, "case": name, "status": result["status"],
            "achievable_or_imposed_growth_per_h": actual,
            "growth_shortfall_per_h": max(0.0, growth - actual) if actual is not None else None,
            "protein_g_per_gdw": protein,
            "protein_original_budget_g_per_gdw": host.reactions.get_by_id("prot_pool_exchange").upper_bound,
            "minimum_discrepancy": result.get("minimum_relaxation"),
            "discrepancy_unit": result.get("relaxation_unit"),
            "minimum_energy_reduction_mmol_atp_per_gdw_h": energy.get("minimum_original_energy_reduction_mmol_atp_per_gdw_h"),
            "compatible_ngam_mmol_atp_per_gdw_h": energy.get("compatible_ngam_at_original_gam_mmol_atp_per_gdw_h"),
            "compatible_gam_mmol_atp_per_gdw": energy.get("compatible_gam_at_original_ngam_mmol_atp_per_gdw"),
            "max_mass_residual": result["audit"]["summary"]["max_abs_mass_balance"] if result["audit"] else None,
            "max_constraint_violation": result["audit"]["summary"]["max_constraint_violation"] if result["audit"] else None,
            "reason": result["reason"], "report": f"cases/{case_id}.json.gz",
        })
        _json(reports / f"{case_id}.json.gz", result)
        print(json.dumps({k: rows[-1][k] for k in ("case_id", "status", "achievable_or_imposed_growth_per_h", "minimum_discrepancy")}, allow_nan=False), flush=True)

    start = time.monotonic()
    for condition in conditions:
        for mode in ("uptake_caps", "rounded_uptake_caps", "quantified_point", "printed_rounding"):
            fixed_growth = mode in ("quantified_point", "printed_rounding")
            constraints = condition_constraints(condition, observables, mode=mode, impose_growth=fixed_growth)
            record(f"baseline_{mode}", condition, solve_native(
                host, constraints=constraints, objective=observables["growth"], timeout_s=args.timeout_s,
            ))
        record("energy_allowance", condition, _energy_allowance(host, condition, observables, args.timeout_s))
        larger_pool = host.copy()
        pool = larger_pool.reactions.get_by_id("prot_pool_exchange")
        if {m.id: c for m, c in pool.metabolites.items()} != {"prot_pool[c]": 1.0}:
            raise ValueError("declared non-molecular protein resource stoichiometry mismatch")
        pool.upper_bound *= 2.0
        record("protein_budget_double_diagnostic", condition, solve_native(
            larger_pool, constraints=condition_constraints(condition, observables, mode="uptake_caps"),
            objective=observables["growth"], timeout_s=args.timeout_s,
        ), {"calibrated": False, "changed_bound": ["prot_pool_exchange", "upper", pool.upper_bound],
            "unit": "g protein/gDW", "source": "Counterfactual 2x original budget, not an inferred biological parameter"})
        exact = condition_constraints(condition, observables, mode="quantified_point", impose_growth=True)
        scales = {c.observable.key: (condition["observations"][c.observable.key]["rounding_half_width"],) * 2
                  for c in exact if c.observable.key != "growth"}
        record("minimum_joint_exchange_rounding_units", condition, minimum_relaxation(
            host, constraints=exact, relaxations=scales, max_relaxation=20.0,
            relaxation_unit="multiples of the printed half-last-digit", timeout_s=args.timeout_s,
        ))
        if condition["growth_rate_per_h"] in (0.025, 0.1, 0.4):
            caps = condition_constraints(condition, observables, mode="uptake_caps", impose_growth=True)
            for key in ("GlucoseUptake", "O2uptake"):
                record(f"minimum_{key}_cap_increase", condition, minimum_relaxation(
                    host, constraints=caps, relaxations={key: (0.0, 1.0)},
                    max_relaxation=max(1.0, condition["observations"][key]["value"]),
                    relaxation_unit="mmol/gDW/h", timeout_s=args.timeout_s,
                ))
            record("minimum_relative_uptake_cap_increase", condition, minimum_relaxation(
                host, constraints=caps,
                relaxations={key: (0.0, condition["observations"][key]["value"])
                             for key in ("GlucoseUptake", "O2uptake")},
                max_relaxation=1.0, relaxation_unit="fraction of each reported uptake cap", timeout_s=args.timeout_s,
            ), {"scope": "Capacity-envelope sensitivity only; not an identified dry-mass conversion"})
    if not args.skip_obligations:
        for scenario in obligations.values():
            for endpoint in ("minimum", "maximum"):
                name = f"frozen_hog_NaCl_{scenario['nacl_molar']:g}_{endpoint}"
                obligation = NativeCarbonObligation(
                    scenario[endpoint], scenario["source"], scenario["source_medium"],
                    metadata={"native_training_scope": "glucose_only", "frozen_checkpoint_sha256": FROZEN_SHA256,
                              "source_nacl_molar": scenario["nacl_molar"], "endpoint": endpoint,
                              "scope": scenario["scope"]},
                )
                try:
                    installation_host = host.copy()
                    installation_host.solver.configuration.timeout = args.timeout_s
                    conditioned, report = impose_native_obligation(
                        installation_host, roles=OBLIGATION_ROLES, obligation=obligation,
                        demand_reaction_id="native_reconciliation_glycerol",
                    )
                    metadata = asdict(report)
                    bounds = metadata["metadata"]["applied_boundary_input_bounds"]
                    metadata["metadata"]["applied_boundary_input_bounds"] = {
                        rid: [bound_metadata(v) for v in pair] for rid, pair in bounds.items()
                    }
                    for condition in conditions:
                        result = solve_native(
                            conditioned, constraints=condition_constraints(condition, observables, mode="uptake_caps"),
                            objective=observables["growth"], timeout_s=args.timeout_s,
                        )
                        if result["status"] == "optimal":
                            flux = result["fluxes"]
                            result["native_commitment"] = {
                                "required_mmol_per_gdw_h": report.native_moles_per_glucose_mole * math.fsum(
                                    c * flux[rid] for rid, c in report.glucose_flux_coefficients.items()),
                                "existing_export_mmol_per_gdw_h": flux["r_1808"],
                                "new_demand_mmol_per_gdw_h": flux[report.demand_reaction_id],
                                "total_terminal_mmol_per_gdw_h": math.fsum(c * flux[rid] for rid, c in report.terminal_flux_coefficients.items()),
                                "accounting": "Existing export contributes to the same total terminal requirement; it is not charged a second time",
                            }
                        record(name, condition, result, metadata)
                except (ValueError, RuntimeError) as exc:
                    for condition in conditions:
                        record(name, condition, {"status": "failed", "reason": str(exc), "objective_value": None,
                                                "fluxes": None, "observables": None, "audit": None},
                               {"stage": "native obligation installation", "fraction": scenario[endpoint]})
    if _model_identity(host) != host_identity or _digest(root / MODEL_PATH) != model_sha256:
        raise RuntimeError("original native host changed during diagnostic experiments")
    _, after = _frozen_native_context(root, **artifact_options)
    if after != frozen:
        raise RuntimeError("frozen native inputs changed during reconciliation")
    with (output / "case_summary.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    _json(output / "case_summary.json", rows)
    _json(output / "completion.json", {
        "status": "completed", "n_cases": len(rows), "elapsed_s": time.monotonic() - start,
        "status_counts": {status: sum(r["status"] == status for r in rows) for status in ("optimal", "infeasible", "failed")},
        "source_model_unchanged": True, "frozen_checkpoint": after,
        "calibration_adopted": False,
        "conclusion_scope": "Conditional native consistency diagnostics, not a unique biological law or an independent validation of a new calibrated model",
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
