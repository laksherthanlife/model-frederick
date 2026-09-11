"""Fit and solve explicit state policies, writing only to a requested empty directory.

Usage:
    python scripts/run_state_allocation.py --synthetic-demo --output-dir /tmp/allocation-demo
    python scripts/run_state_allocation.py --config experiment.json --output-dir /tmp/allocation-real

A schema_version=1 JSON config supplies model (SBML path, relative to config),
identity (ModelIdentity fields), features (Feature fields), states (StatePolicy
fields), observations (StateObservation fields, including split and group_id),
fit (projection_dimension, ridge, basis), and baselines. Every source is a Source
object with kind/reference/method. Environmental values have explicit units.
Uptake feature names refer to the feature table. Baselines must explicitly name
mechanistic_state and fixed_objective_state; fixed_transition is optional and
supplies states, feature, midpoint, width, basis, source. No label or kinetic table
is imputed. Optional strain_transfer is a reason, not a cross-strain validation.

The demo's labels and capacities are synthetic engineering fixtures. Neither its
scores nor a successful fit establish learned yeast biology or reproduce COSMIC.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import cobra

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ystwin.bridge.regulation import MetabolicTask
from ystwin.fba.allocation import (
    Environment, Feature, FixedObjectivePolicy, FixedTransition, FluxCapacity, MechanisticPolicy,
    MixturePolicy, ModelIdentity, ObjectiveTask, ProteinBudget, Source, StateObservation, StatePolicy,
    UptakePolicy, fit_state_distribution, solve_allocation,
)
from ystwin.fba.solver import configure, load_model


def _source(row):
    return Source(**row)


def _environment(row):
    return Environment(**{**row, "source": _source(row["source"])})


def _with_source(cls, row, **overrides):
    return cls(**{**row, "source": _source(row["source"]), **overrides})


def load_configuration(path):
    raw = path.read_bytes()
    document = json.loads(raw)
    if document.get("schema_version") != 1:
        raise ValueError("config schema_version must be 1")
    model_path = Path(document["model"])
    if not model_path.is_absolute():
        model_path = path.resolve().parent / model_path
    model, settings = load_model(model_path)
    identity = _with_source(ModelIdentity, document["identity"])
    features = tuple(Feature(**row) for row in document["features"])
    by_feature = {feature.name: feature for feature in features}
    states = []
    for row in document["states"]:
        uptakes = tuple(_with_source(UptakePolicy, uptake,
                                    feature=by_feature[uptake["feature"]] if uptake.get("feature") else None)
                        for uptake in row.get("uptakes", []))
        capacities = tuple(_with_source(FluxCapacity, capacity) for capacity in row.get("enzyme_capacities", []))
        tasks = tuple(_with_source(ObjectiveTask, task) for task in row.get("tasks", []))
        measured = tuple(MetabolicTask(**task) for task in row.get("measured_tasks", []))
        budget = _with_source(ProteinBudget, row["protein_budget"]) if row.get("protein_budget") else None
        mass_source = _source(row["dry_mass_source"]) if row.get("dry_mass_source") else None
        states.append(_with_source(StatePolicy, row, uptakes=uptakes, enzyme_capacities=capacities,
                                   tasks=tasks, measured_tasks=measured, protein_budget=budget,
                                   dry_mass_source=mass_source))
    observations = tuple(_with_source(StateObservation, row, environment=_environment(row["environment"]))
                         for row in document["observations"])
    baselines = dict(document["baselines"])
    if baselines.get("fixed_transition") is not None:
        transition = baselines["fixed_transition"]
        baselines["fixed_transition"] = _with_source(FixedTransition, transition,
                                                     feature=by_feature[transition["feature"]],
                                                     states=tuple(transition["states"]))
    metadata = {"config_sha256": hashlib.sha256(raw).hexdigest(), "model_path": str(model_path.resolve()),
                "solver": str(settings), "data_kind": "+".join(sorted({row.source.kind for row in observations})),
                "strain_transfer": document.get("strain_transfer")}
    return model, identity, features, tuple(states), observations, document["fit"], baselines, metadata


def synthetic_demo():
    source = Source("synthetic", "run_state_allocation analytic toy", "constructed labels/capacities; no measured biology")
    model = cobra.Model("allocation_demo")
    metabolites = {name: cobra.Metabolite(name, compartment=compartment)
                   for name, compartment in (("g_e", "e"), ("g_c", "c"), ("o_e", "e"), ("o_c", "c"), ("p_c", "c"))}
    definitions = (
        ("EX_G", {"g_e": -2.0}, (-1000.0, 1000.0)),
        ("EX_O", {"o_e": 1.0}, (0.0, 1000.0)),
        ("TG", {"g_e": -1.0, "g_c": 1.0}, (0.0, 1000.0)),
        ("TO", {"o_e": -1.0, "o_c": 1.0}, (0.0, 1000.0)),
        ("BIO", {"g_c": -1.0, "o_c": -1.0}, (0.0, 1000.0)),
        ("MAKE_P", {"g_c": -1.0, "p_c": 1.0}, (0.0, 1000.0)),
        ("PRODUCT", {"p_c": -1.0}, (0.0, 1000.0)),
    )
    for name, stoichiometry, bounds in definitions:
        reaction = cobra.Reaction(name, lower_bound=bounds[0], upper_bound=bounds[1])
        reaction.add_metabolites({metabolites[key]: value for key, value in stoichiometry.items()})
        model.add_reactions([reaction])
    model.objective = "BIO"
    settings = configure(model)
    identity = ModelIdentity(model.id, "synthetic organism", "toy strain", "analytic-v1", source)
    features = (Feature("oxygen", "mmol/L"), Feature("growth", "1/h", "growth"),
                Feature("elapsed", "h", "time"))
    base = StatePolicy(
        "off", source, biomass_reaction="BIO",
        uptakes=(UptakePolicy("EX_G", 10.0, source),
                 UptakePolicy("EX_O", 10.0, source, feature=features[0], half_saturation=0.1)),
        enzyme_capacities=(FluxCapacity("BIO", 2.0, 2.0, source, unit="1/h"),
                           FluxCapacity("MAKE_P", 0.0, 0.0, source, unit="mmol/gDW/h")),
        tasks=(ObjectiveTask("PRODUCT", 1.0, source),),
        dry_mass_g_per_cell=1e-12, dry_mass_source=source,
    )
    production = replace(base, name="on", dry_mass_g_per_cell=3e-12,
                         uptakes=(UptakePolicy("EX_G", 12.0, source),
                                  UptakePolicy("EX_O", 8.0, source, feature=features[0], half_saturation=0.1)),
                         enzyme_capacities=(base.enzyme_capacities[0],
                                            FluxCapacity("MAKE_P", 0.0, 8.0, source, unit="mmol/gDW/h")))
    observations = []
    for split, data in (("train", ((1.0, 0.05), (2.0, 0.1), (8.0, 0.9), (9.0, 0.95))),
                        ("test", ((3.0, 0.2), (7.0, 0.8)))):
        for index, (oxygen, fraction) in enumerate(data):
            identifier = f"{split}-{index}"
            environment = Environment(identifier, identity.organism, identity.strain,
                                      {"oxygen": oxygen, "growth": 2.0, "elapsed": 12.0},
                                      {"oxygen": "mmol/L", "growth": "1/h", "elapsed": "h"}, source)
            observations.append(StateObservation(identifier, f"culture-{identifier}", split, environment,
                                                 {"off": 1.0 - fraction, "on": fraction}, source))
    fit = {"projection_dimension": 1, "ridge": 0.01, "basis": "cell"}
    baselines = {"mechanistic_state": "on", "fixed_objective_state": "on",
                 "fixed_transition": FixedTransition(("off", "on"), features[2], 12.0, 1.0, source, "cell")}
    metadata = {"config_sha256": None, "model_path": None, "solver": str(settings), "data_kind": "synthetic",
                "strain_transfer": None}
    return model, identity, features, (base, production), tuple(observations), fit, baselines, metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--config", type=Path)
    choice.add_argument("--synthetic-demo", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise SystemExit("output directory must be empty; prior results will not be overwritten")
    try:
        model, identity, features, states, observations, fit, baselines, metadata = (
            synthetic_demo() if args.synthetic_demo else load_configuration(args.config)
        )
        fitted = fit_state_distribution(observations, features=features, states=[state.name for state in states], **fit)
        by_name = {state.name: state for state in states}
        policies = [MixturePolicy(states, fitted),
                    MechanisticPolicy(by_name[baselines["mechanistic_state"]]),
                    FixedObjectivePolicy(by_name[baselines["fixed_objective_state"]])]
        if baselines.get("fixed_transition") is not None:
            policies.append(MixturePolicy(states, baselines["fixed_transition"]))
        results, summary = [], []
        for row in observations:
            for policy in policies:
                allocation = solve_allocation(model, identity, row.environment, policy,
                                               strain_transfer=metadata["strain_transfer"])
                results.append({"sample_id": row.sample_id, "group_id": row.group_id, "split": row.split,
                                "allocation": allocation.report()})
                summary.append({"sample_id": row.sample_id, "group_id": row.group_id, "split": row.split,
                                "policy": allocation.policy, "status": allocation.status,
                                "growth_rate_per_h": allocation.growth_rate,
                                "cell_population_growth_rate_per_h": allocation.cell_population_growth_rate,
                                "mass_balance_residual": allocation.mass_balance_residual,
                                "biomass_fractions": json.dumps(allocation.biomass_fractions, sort_keys=True)})
        payload = {**metadata, "schema_version": 1, "biological_validation": False,
                   "claim": "explicit policy solves and a fit to sourced labels; not independent biological validation",
                   "cosmic_evidence": "data/cosmic_sources.json; original fitted tables not available",
                   "evaluation": "train/test labels are not pooled for fitting; no held-out label or flux score is claimed",
                   "results": results}
        encoded = json.dumps(payload, indent=2, allow_nan=False)
        trained = json.dumps(fitted.report(), indent=2, allow_nan=False)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        raise SystemExit(str(exc)) from exc
    output.mkdir(parents=True, exist_ok=True)
    (output / "state_allocation.json").write_text(encoded + "\n")
    (output / "state_distribution.json").write_text(trained + "\n")
    with (output / "state_allocation.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f"wrote {output}; {len(results)} solves; data_kind={metadata['data_kind']}; biological_validation=false")


if __name__ == "__main__":
    main()
