from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pathlib
import textwrap
from collections import Counter
from dataclasses import asdict, fields

import pandas as pd

from ystwin import paths
from ystwin.fba.solver import load_model
from ystwin.fba.stress_pathways import COVERAGE_PATH, load_coverage, model_gene_coverage
from ystwin.generator.in_silico import CONTROL_NAMES, INPUT_NAMES, LATENT_NAMES, OBSERVATION_NAMES, TeacherParameters
from ystwin.generator.panel_experiment import MECHANISED_CHANNELS
from ystwin.generator.stress_panel import METABOLITE_POOLS, MODULES, STRESSORS

MAP_PATH = paths.data_dir() / "gem" / "stress_response_map.json"


def validate_map(mapping, items):
    if mapping.get("schema_version") != 1:
        raise ValueError("unsupported stress-map schema")
    known_items = {item.item for item in items}
    seen = set()
    required = {"id", "title", "input", "signal", "regulators", "effectors", "catalogue_modules",
                "coverage_items", "implementation_refs", "learner_aggregates", "scope_note"}
    for branch in mapping["branches"]:
        if required - branch.keys() or branch["id"] in seen:
            raise ValueError("incomplete or duplicate stress-map branch")
        seen.add(branch["id"])
        if set(branch["coverage_items"]) - known_items:
            raise ValueError(f"unknown coverage item in {branch['id']}")
        if set(branch["catalogue_modules"]) - MODULES.keys():
            raise ValueError(f"unknown catalogue programme in {branch['id']}")
        if set(branch["learner_aggregates"]) - set(LATENT_NAMES):
            raise ValueError(f"unknown learned aggregate in {branch['id']}")
        for module in branch["implementation_refs"]:
            try:
                found = importlib.util.find_spec(f"ystwin.{module}")
            except ModuleNotFoundError:
                found = None
            if found is None:
                raise ValueError(f"unknown implementation reference {module}")
    return mapping


def teacher_assumptions(parameters=None):
    parameters = TeacherParameters() if parameters is None else parameters
    groups = {
        "synthetic_state_law": (
            "family", "tau_u_h", "tau_o_h", "tau_b_h", "cross_uo_per_h", "burden_weights",
            "saturation_half", "saturation_power"),
        "virtual_sensor_design": (
            "reporter_basal_synthesis_per_h", "reporter_loadings", "reporter_kinetics"),
        "synthetic_culture": ("assay_growth_rate_per_h", "initial_cell_density"),
        "numerical_setting": ("solver_rtol", "solver_atol", "solver_max_step_h", "solver_method"),
        "synthetic_noise_default": ("default_noise_cv", "noise_model"),
        "frozen_synthetic_control_law": (
            "control_feature_names", "growth_cost_coefficients", "enzyme_cost_coefficients",
            "ngam_coefficients", "allocation_logit_coefficients"),
    }
    roles = {name: role for role, names in groups.items() for name in names}
    declared = {field.name: field for field in fields(parameters)}
    if roles.keys() != declared.keys():
        raise ValueError("teacher parameter inventory must classify every field explicitly")
    values = asdict(parameters)
    return pd.DataFrame([
        {"parameter": name, "role": roles[name], "value": json.dumps(value),
         "constructor_input": declared[name].init,
         "source": "declared reference-teacher assumption, not fitted biological evidence",
         "student_relation": "student fits separate weights from generated labels; this value is not learned here"}
        for name, value in values.items()
    ])


def build_report(models, mapping, items):
    validate_map(mapping, items)
    catalogue = {item.item: item for item in items}
    rows, branches, model_summary = [], [], {}
    for role, model in models.items():
        observed = {row.item: row for row in model_gene_coverage(model, items)}
        model_summary[role] = {
            "model_id": model.id, "genes": len(model.genes), "reactions": len(model.reactions),
            "metabolites": len(model.metabolites),
            "catalogue_items": len(items),
            "items_with_gpr_associations": sum(bool(row.reaction_ids) for row in observed.values()),
            "items_by_reference_diagram_layer": dict(Counter(item.layer for item in items)),
        }
        for key, association in observed.items():
            item = catalogue[key]
            rows.append({
                "model_role": role, "model_id": model.id, "item": key, "label": item.label,
                "reference_diagram_layer": item.layer, "historical_verdict": item.verdict,
                "declared_orfs": "|".join(association.declared_orfs),
                "present_orfs": "|".join(association.present_orfs),
                "absent_orfs": "|".join(association.absent_orfs),
                "actual_gpr_reactions": "|".join(association.reaction_ids),
                "evidence_scope": association.evidence_scope,
                "source": item.source,
            })
        for branch in mapping["branches"]:
            queried = [observed[key] for key in branch["coverage_items"]]
            present = sorted({orf for row in queried for orf in row.present_orfs})
            absent = sorted({orf for row in queried for orf in row.absent_orfs})
            reactions = sorted({reaction for row in queried for reaction in row.reaction_ids})
            mechanisms = [f"{stressor}:{programme}:{channel}"
                          for (stressor, programme), channel in MECHANISED_CHANNELS.items()
                          if programme in branch["catalogue_modules"]]
            branches.append({
                "model_role": role, "model_id": model.id, "branch": branch["id"],
                "title": branch["title"], "queried_items": len(queried),
                "present_orfs": "|".join(present), "absent_orfs": "|".join(absent),
                "actual_gpr_reactions": "|".join(reactions),
                "catalogue_programmes": "|".join(branch["catalogue_modules"]),
                "optional_panel_mechanisms": "|".join(mechanisms),
                "separate_implementation_refs": "|".join(branch["implementation_refs"]),
                "current_synthetic_aggregate_only": "|".join(branch["learner_aggregates"]),
                "scope_note": branch["scope_note"],
            })
    summary = {
        "schema_version": 1,
        "purpose": "coverage and assumption audit; no metabolic optimization or biological validation",
        "models": model_summary,
        "older_panel": {
            "catalogue_programmes": list(MODULES), "catalogue_count": len(MODULES),
            "transcriptional_programme_count": len(MODULES.keys() - METABOLITE_POOLS),
            "metabolite_pools": sorted(METABOLITE_POOLS),
            "stressor_count": len(STRESSORS), "stressors": list(STRESSORS),
            "optional_mechanistic_reporter_routes": [
                {"stressor": key[0], "programme": key[1], "reporter": channel}
                for key, channel in MECHANISED_CHANNELS.items()
            ],
            "meaning": "catalogued phenomenological programmes are not all dynamical or identifiable states",
        },
        "current_student": {
            "imposed_inputs": list(INPUT_NAMES),
            "latent_coordinates": list(LATENT_NAMES), "control_outputs": list(CONTROL_NAMES),
            "observation_channels": list(OBSERVATION_NAMES),
            "state_semantics": "supervised synthetic aggregate coordinates, not identified molecular cascades",
            "teacher_source": TeacherParameters().provenance["parameter_source"],
            "product_labels_used_for_training": False,
            "shared_metabolic_backend_in_reference_comparison": True,
            "independent_biological_validation": False,
            "not_learned": ["GEM stoichiometry", "all stress pathways", "gene-specific regulation",
                            "real sensor-to-pathway calibration", "real stress-to-enzyme control laws"],
        },
        "interpretation": mapping["interpretation"],
        "missing_quantitative_layers": mapping["missing_quantitative_layers"],
    }
    return summary, pd.DataFrame(rows), pd.DataFrame(branches)


def render_map(mapping, branch_table, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    columns = 3
    nrows = (len(mapping["branches"]) + columns - 1) // columns
    figure = plt.figure(figsize=(25, 7 + nrows * 5))
    axis = figure.add_axes((0, 0, 1, 1))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.text(0.025, 0.985, mapping["title"], fontsize=22, weight="bold", va="top")
    axis.text(0.025, 0.957,
              "BIOLOGICAL MAP, NOT A CLAIM THAT THE CURRENT MODEL IMPLEMENTS THESE CASCADES\n"
              "Within each branch: input -> sensing/signaling -> transcriptional regulation -> metabolic or population effects.",
              fontsize=12, va="top", color="#374151")
    top, bottom = 0.915, 0.18
    height = (top - bottom) / nrows
    role = "ec_predictor" if "ec_predictor" in set(branch_table.model_role) else branch_table.model_role.iloc[0]
    evidence = branch_table[branch_table.model_role == role].set_index("branch")
    for index, branch in enumerate(mapping["branches"]):
        row, column = divmod(index, columns)
        x, y = 0.025 + column * 0.325, top - (row + 1) * height
        width, box_height = 0.307, height - 0.012
        axis.add_patch(FancyBboxPatch((x, y), width, box_height,
                                      boxstyle="round,pad=0.005,rounding_size=0.004",
                                      facecolor="#f6f8fa", edgecolor="#aeb9c4", linewidth=1))
        axis.text(x + 0.012, y + box_height - 0.011, branch["title"],
                  fontsize=13, weight="bold", va="top", color="#142c40")
        body = []
        for name, key in (("Input", "input"), ("Signal", "signal"), ("Regulation", "regulators"), ("Effects", "effectors")):
            body.extend(textwrap.wrap(f"{name}: {branch[key]}", width=60))
            body.append("")
        axis.text(x + 0.012, y + box_height - 0.026, "\n".join(body), fontsize=10.2,
                  va="top", linespacing=1.22, color="#243645")
        observed = evidence.loc[branch["id"]]
        association = "some queried genes have GPR associations" if observed.actual_gpr_reactions else "no queried GPR associations"
        proxy = ", ".join(branch["learner_aggregates"]) or "none"
        footer = f"EC evidence: {association}.\nCurrent learned aggregate proxy: {proxy}.\nNot evidence that this signaling branch is implemented."
        axis.text(x + 0.012, y + 0.011, footer, fontsize=9.5, va="bottom", color="#623d80")
    axis.text(0.025, 0.157, "SHARED FEEDBACKS AND MISSING STATE", fontsize=14, weight="bold", va="top")
    feedback = "\n".join("\n".join(textwrap.wrap(line, width=185)) for line in mapping["shared_feedbacks"])
    axis.text(0.025, 0.142, feedback, fontsize=11, va="top", linespacing=1.3)
    axis.add_patch(FancyBboxPatch((0.02, 0.018), 0.955, 0.068,
                                  boxstyle="round,pad=0.004,rounding_size=0.004",
                                  facecolor="#eef0ff", edgecolor="#656caf"))
    axis.text(0.033, 0.075, "CURRENT EXECUTABLE LEARNING LOOP", fontsize=13, weight="bold", va="top")
    axis.text(0.033, 0.057,
              f"Virtual reporter history -> synthetic {' / '.join(LATENT_NAMES)} aggregates -> {len(CONTROL_NAMES)} global controls -> supplied GEM -> product pools\n"
              "These aggregate states do not instantiate the detailed pathways above. Shared-teacher recovery is an engineering check, not biological validation.",
              fontsize=12, va="top", linespacing=1.4)
    figure.savefig(output / "yeast_stress_expanded.svg", bbox_inches="tight")
    figure.savefig(output / "yeast_stress_expanded.png", dpi=130, bbox_inches="tight")
    plt.close(figure)


def render_overview(output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    figure, axis = plt.subplots(figsize=(18, 12))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.text(0.02, 0.985, "A yeast stress model needs more than a GEM", fontsize=22, weight="bold", va="top")
    axis.text(0.02, 0.94, "Conceptual architecture: arrows identify required couplings, not fitted quantitative laws.", fontsize=12)

    def box(x, y, title, text, color, width=0.27, height=0.17):
        axis.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.008",
                                      facecolor=color, edgecolor="#405064", linewidth=1.2))
        axis.text(x + 0.014, y + height - 0.021, title, fontsize=12, weight="bold", va="top")
        axis.text(x + 0.014, y + height - 0.055, text, fontsize=10.5, va="top", linespacing=1.35)

    def arrow(start, end, label="", color="#526174", curve=0, dashed=False):
        axis.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                                       connectionstyle=f"arc3,rad={curve}", color=color,
                                       linewidth=1.5, linestyle="--" if dashed else "-"))
        if label:
            axis.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.015,
                      label, fontsize=9, ha="center", color=color,
                      bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5})

    box(0.03, 0.72, "ENVIRONMENT AND GENOTYPE",
        "Heat, pH, oxidants, osmolarity\nNutrients, oxygen, toxicants\nExpression, secretion, plasmids", "#f3f4f6")
    box(0.365, 0.72, "PHYSICAL CELL STATE",
        "Misfolding, ROS, pH, turgor\nATP/redox pools; membrane state\nDamage, repair and organelles", "#fff4d9")
    box(0.70, 0.72, "SIGNALING AND REGULATION",
        "HOG / CWI / UPR / HSF / ESR\nSnf1 / PKA / TOR / Rim / Crz1\nActivation, localization, repression", "#fff4d9")
    box(0.70, 0.425, "PROTEOME AND ENZYME STATE",
        "Transcripts are not enzymes\nTranslation, turnover, activity\nCapacities and explicit resource costs", "#fff4d9")
    box(0.365, 0.425, "METABOLIC MODEL / GEM",
        "Stoichiometry and GPR associations\nFlux feasibility, uptake, enzyme budgets\nNot the signaling network above", "#e4f1e6")
    box(0.03, 0.425, "PROCESS AND POPULATION",
        "Biomass, substrate, product inventories\nViability, producer fraction, aggregation\nOnly a subset is in the current loop", "#f3f4f6")
    arrow((0.30, 0.805), (0.365, 0.805))
    arrow((0.635, 0.805), (0.70, 0.805))
    arrow((0.835, 0.72), (0.835, 0.595), "expression / activity")
    arrow((0.70, 0.51), (0.635, 0.51))
    arrow((0.365, 0.51), (0.30, 0.51))
    arrow((0.47, 0.595), (0.47, 0.72), "energy / redox feedback", curve=0.12)
    arrow((0.16, 0.595), (0.40, 0.72), "growth / demand feedback", curve=0.18)

    axis.text(0.02, 0.35, "CURRENT SYNTHETIC QUALIFICATION LANE", fontsize=13, weight="bold", color="#4c418c")
    box(0.03, 0.115, "OBSERVATION PHYSICS",
        "State-dependent synthetic reporters\nReporter kinetics and growth dilution\nReal pH/optical confounds are not fitted", "#eef0ff")
    box(0.365, 0.115, "LEARNED SYNTHETIC STATE",
        "Causal synthesis-activity estimates\nUPR / oxidative / burden aggregates\nNot all molecular pathways", "#eef0ff")
    box(0.70, 0.115, "LEARNED GLOBAL CONTROLS",
        "Growth retention; protein budget\nATP maintenance; allocation fraction\nReference law is a synthetic assumption", "#eef0ff")
    arrow((0.16, 0.425), (0.16, 0.285), "density / reporter context")
    arrow((0.38, 0.72), (0.28, 0.285), curve=0.04)
    arrow((0.30, 0.20), (0.365, 0.20), color="#6255aa")
    arrow((0.635, 0.20), (0.70, 0.20), color="#6255aa")
    arrow((0.835, 0.285), (0.58, 0.425), "current GEM coupling", color="#6255aa", dashed=True)
    axis.text(0.02, 0.04,
              "Shared-teacher agreement verifies this lower lane under supplied assumptions.\n"
              "Independent mechanisms, interventions and biological observations are needed to validate the upper network.",
              fontsize=11, color="#374151")
    figure.savefig(output / "yeast_stress_architecture.svg", bbox_inches="tight")
    figure.savefig(output / "yeast_stress_architecture.png", dpi=150, bbox_inches="tight")
    plt.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Audit stress-map scope against actual GEM gene associations")
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--map", type=pathlib.Path, default=MAP_PATH)
    parser.add_argument("--coverage", type=pathlib.Path, default=COVERAGE_PATH)
    parser.add_argument("--no-render", action="store_true")
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise SystemExit(f"Refusing to overwrite existing audit: {output}")
    items = load_coverage(args.coverage)
    mapping = json.loads(args.map.read_text())
    model_paths = {
        "ec_predictor": paths.require(paths.ec_yeast_gem(), "EC model", "YSTWIN_EC_YEAST_GEM"),
        "yeast9_reference": paths.require(paths.yeast_gem(), "yeast-GEM", "YSTWIN_YEAST_GEM"),
    }
    models = {name: load_model(path)[0] for name, path in model_paths.items()}
    summary, associations, branches = build_report(models, mapping, items)
    summary["input_hashes"] = {
        **{name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in model_paths.items()},
        "coverage_table": hashlib.sha256(args.coverage.read_bytes()).hexdigest(),
        "conceptual_map": hashlib.sha256(args.map.read_bytes()).hexdigest(),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "scope.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    associations.to_csv(output / "gene_associations.csv", index=False)
    branches.to_csv(output / "pathway_scope.csv", index=False)
    teacher_assumptions().to_csv(output / "teacher_assumptions.csv", index=False)
    if not args.no_render:
        render_map(mapping, branches, output)
        render_overview(output)
    print(json.dumps(summary["models"], indent=2))
    print(f"Catalogue: {len(MODULES)} programmes/pools, {len(STRESSORS)} stressors")
    print(f"Optional mechanistic reporter routes: {len(MECHANISED_CHANNELS)}")
    print(f"Current synthetic state coordinates: {list(LATENT_NAMES)}")
    print(f"Stress-map audit: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
