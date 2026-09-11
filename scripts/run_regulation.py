"""D3: does per-reaction regulation narrow the product flux range, where one ATP scalar did not?

D1 collapsed the whole 24-module latent state into a single non-growth ATP maintenance
term (`bridge/latent_bridge.py`, `_MAINTENANCE_PER_ACTIVITY = 300.0`, the project's own
audit calling it "a plausible scale, not a fitted one") and then reported the beta-carotene
range as `[0, ceiling]` with relative width 1.000 at every growth level. Its conclusion --
FBA cannot forecast product -- may have been about the coupling rather than about FBA.

This script replaces the scalar with a regulation layer: module activity -> the genes each
regulon binds (SGD, MacIsaac 2006 conserved motifs) -> the reactions those genes control
(the GEM's own gene-reaction rules) -> bounds on those reactions. Then it re-runs
`product_flux_range` at D1's own growth fractions, with and without the layer, and prints
the two numbers side by side.

Five sections, and section 3 is the answer:

  1. What the mapping actually reaches. Printed first, because a narrowing computed on a
     mapping that reaches nothing is a number about nothing.
  2. Which direction the panel can even push. `STRESSORS[...].targets` carries signed
     weights and the sign is honoured, so this is a measurement of the panel, not an
     assumption about it.
  3. E-Flux (Colijn 2009, PMID 19714220) against the unregulated model, at D1's growth
     fractions: ceiling, floor, relative width.
  4. A counterfactual repression probe. Section 2 determines whether the panel supplies a
     repressive regulon; if it does not, the mechanism is still worth exercising, so the
     activities are sign-flipped and clearly labelled as a counterfactual. It answers
     "would this layer bind if the panel could repress", which is a different question
     from "does it bind today" and must not be reported as the same one.
  5. A measured task-efficiency adaptation inspired by COSMIC-dFBA (PMID 38387677).
     The exact paper training procedure is unverified; data/cosmic_sources.json records
     the accessible public kernel and missing tables. The entire measured task set is
     retained, including biomass, and unmeasured product rates are refused.

Usage: python scripts/run_regulation.py --tasks-only --output-dir /tmp/regulation-audit
Supply --strain-transfer with an explicit reason to compare the S288C model with DS28911.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import os
import pathlib
import sys

import cobra
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.bridge.regulation import (
    COLIJN_PMID,
    COSMIC_PMID,
    EFLUX_LIFTS_LOWER_BOUND,
    InfeasibleRegulation,
    MetabolicTask,
    MissingMeasurement,
    apply_eflux,
    apply_task_bounds,
    constrain_uptake,
    eflux_layer,
    load_regulons,
    task_priority_order,
)
from ystwin.fba.allocation import Environment, ModelIdentity, Source, require_matched_models, validate_host
from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID, add_beta_carotene_pathway
from ystwin.fba.fva import product_flux_range
from ystwin.fba.physiology import REFERENCE_AEROBIC_BATCH
from ystwin.fba.solver import FVA_PROCESSES, growth_or_none, load_model

from ystwin.generator.stress_panel import MODULES, module_response

OUT = None
BIOMASS = "r_2111"
ETHANOL = "r_1761"
CO2 = "r_1672"

# D1's own grid, so the before/after columns are comparable to the published table.
GROWTH_FRACTIONS = (1.0, 0.99, 0.95, 0.90, 0.75, 0.50, 0.25)

# The five metabolite pools have no transcription factor and therefore no regulon; calcium
# has a factor (Crz1) but zero target genes at the MacIsaac evidence level. Both are
# declared here rather than filtered silently, which is what `eflux_layer` insists on.
POOLS = frozenset({"redox", "peroxide", "atp", "ph", "nadh"})
UNMAPPED = POOLS | {"calcium"}

# One agent per mechanism the panel covers, dosed at its own EC50 so the activities are
# mid-scale rather than saturated.
PROBES = (("H2O2", 0.15), ("DTT", 1.0), ("menadione", 100.0), ("tunicamycin", 1.0))

# A refused task set is a result. This records WHY, so the refusal is readable as a finding
# rather than as a gate that someone should quietly widen.
COSMIC_SCOPED_OUT = """
SCOPED-OUT HOST-MODEL INCOMPATIBILITY -- a correct outcome, recorded, not hidden.
  Measured 0.400 /h against a model maximum of 0.346 /h, efficiency 1.156, with glucose
  11.1 and oxygen 3.7 both pinned at the same measured condition (van Hoek 1998 Table 1,
  D = 0.40, PMID 9797269). Nothing below is relaxed, tuned or dropped to remove it.

What refuses. yeast-GEM carries no overflow (proteome-allocation) mechanism, and D = 0.40
  is the most fermentative row in the table: 41.7% of substrate carbon leaves as ethanol on
  3.7 mmol O2/gDW/h. Handed the same two uptakes the model respires instead, and secretes
  MORE ethanol (16.031 against 13.9) and MORE CO2 (20.424 against 18.9) while growing less.
  Run upstream fitGAM.m's own protocol instead -- fix growth, leave oxygen free, minimise
  glucose -- and the model reproduces this same strain's seven respiratory rows to a mean
  ratio of 1.009 (0.935-1.057), then collapses to 0.542 (0.408-0.756) on the three
  fermentative ones. The sub-critical refusals (1.007 at D = 0.200, 1.019 at D = 0.150) are
  that residual meeting an inequality with zero tolerance, not a disagreement about
  physiology. Deliberately NOT the explanation: the row is carbon-closed at 97.9%, and
  restoring the glycerol/acetate/pyruvate columns leaves it infeasible either way.

Where the number actually lives. r_2111 is a bare drain (s_0450 -->) with no ATP term at
  all. Growth-associated maintenance is 55.3 mmol ATP/gDW in r_4041; non-growth maintenance
  is r_4046, pinned at (0.7, 0.7). Zeroing r_4046 moves the efficiency only 1.156 -> 1.131,
  so maintenance is not the lever.

Provenance of 55.3, traced upstream rather than inferred. yeast-GEM's history.md records
  "Fitted GAM to chemostat data (PR #159)" at v8.3.1, and the shipped coefficient went
  61.9779 (v8.0.0-v8.3.0) -> 55.4 (v8.3.1) -> 55.3 (v8.3.3 through v9.0.2, this model).
  code/otherChanges/fitGAM.m fits it against data/physiology/chemostatData_VanHoek1998.tsv,
  which upstream is the four SUB-CRITICAL rows of this same table (D = 0.025/0.05/0.10/0.15).
  The 40.8 that v8.3.4's modelCorrections.m attributes to van Hoek was never the shipped
  value -- fitGAM overwrote it, and that script sits in code/.deprecated/ by v9.0.2. So 55.3
  is a documented refit, and this refusal is a property of GAM 55.3 read with oxygen pinned.
  It is NOT a units, sign, row-selection or uptake-mode error (fixed and cap both return
  0.3459844408), and it is not the paper's carbon recovery. The one gap upstream leaves open
  is that no GAM refit is recorded after v8.3.1, so 55.3 crossed the later
  biomass-composition curation unchanged; scaleBioMass.m rescales components and never
  touches GAM.

What a GAM change would and would not do, recorded because it is true and REFUSED as a
  repair. No GAM value clears the table: five rows (0.025, 0.100, 0.150, 0.200, 0.250) stay
  above 1 all the way down to GAM 25. But THIS row does clear at GAM 42 and below, the
  upstream-documented 40.8 included -- status optimal, constrained growth 0.4101, biomass
  0.9754, CO2 0.9397, ethanol 0.9109. "No code could ever produce an ordering here" would
  therefore be false. Refitting a vendored constant so this gate passes is still tuning a
  model to a wanted answer, so the shipped 55.3 stands and the set stays refused.

What would resolve it, and it is not a bound moved here. A host model with an overflow or
  proteome-allocation mechanism that reaches the super-critical branch, or a GAM refit on
  the 9.x network published upstream with its own fitting data. Until one of those exists,
  biomass at D = 0.40 is a contradiction between model and measurement, and printing that is
  the result.
"""


def banner(text: str) -> None:
    print("\n" + "=" * 78 + f"\n{text}\n" + "=" * 78)


def _attainable_capacities(model, reactions, growth_fraction: float = 0.90):
    """Each reaction's own FVA maximum at a fraction of maximum growth, unregulated.

    A second, tighter reference for E-Flux than the model's implied ceiling, computed
    rather than asserted. It is still an over-estimate of what any one enzyme carries in a
    given flux distribution -- it is a maximum over all of them -- but it is bounded by
    the same LP the product range is read off, so a scale against it is comparable.
    """
    with model as scoped:
        scoped.objective = BIOMASS
        maximum = growth_or_none(scoped)
        if maximum is None:
            raise RuntimeError(
                f"{model} cannot grow, so there is no optimum to take "
                "a fraction of. Reading this with `slim_optimize` set the biomass lower "
                "bound to nan instead, which FVA accepts and answers around")
        scoped.reactions.get_by_id(BIOMASS).lower_bound = maximum * growth_fraction
        fva = cobra.flux_analysis.flux_variability_analysis(
            scoped, reaction_list=list(reactions), fraction_of_optimum=0.0,
            processes=FVA_PROCESSES,
        )
    return {rid: max(float(fva.loc[rid, "maximum"]), 0.0) for rid in reactions}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Condition-matched regulation and complete measured-task audit")
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--tasks-only", action="store_true")
    parser.add_argument("--ec-model", type=pathlib.Path)
    parser.add_argument("--plain-model", type=pathlib.Path)
    parser.add_argument("--strain-transfer", help="explicit reason for comparing the model host with DS28911")
    parser.add_argument("--uptake-mode", choices=("fixed", "cap"), default="fixed")
    return parser.parse_args(argv)


def model_identity(model):
    source = Source("simulated", "data/gem/MANIFEST.md and the loaded SBML model id/taxonomy",
                    "model provenance, not a measured strain-specific validation")
    if model.id == "yeastGEM_v9__46__0__46__2":
        return ModelIdentity(model.id, "Saccharomyces cerevisiae", "S288C", "yeast-GEM 9.0.2",
                             source, taxonomy_id="559292")
    if model.id == "M_ecYeastGEM_batch_v8__46__3__46__4":
        return ModelIdentity(model.id, "Saccharomyces cerevisiae", "not annotated", "yeast-GEM 8.3.4", source)
    raise ValueError("unknown model identity; provide a provenance-bound identity through run_state_allocation.py")


def reference_environment(reference=REFERENCE_AEROBIC_BATCH):
    return Environment(
        "van-hoek-1998-table1-D0.40", "Saccharomyces cerevisiae", "DS28911",
        {"glucose_uptake": reference.glucose_uptake, "oxygen_uptake": reference.oxygen_uptake,
         "growth_rate": reference.growth_rate},
        {"glucose_uptake": "mmol/gDW/h", "oxygen_uptake": "mmol/gDW/h", "growth_rate": "1/h"},
        Source("measured", reference.source, "one condition-matched chemostat row"),
    )


def measured_tasks(reference=REFERENCE_AEROBIC_BATCH):
    env = reference_environment(reference)
    metadata = {"condition_id": env.condition_id, "host": env.organism}
    return [
        MetabolicTask("biomass", BIOMASS, reference.growth_rate, reference.source, unit="1/h", **metadata),
        MetabolicTask("ethanol", ETHANOL, reference.ethanol_secretion, reference.source,
                      direction="secretion", **metadata),
        MetabolicTask("CO2", CO2, reference.co2_secretion, reference.source,
                      direction="secretion", **metadata),
    ]


def measured_task_report(model, *, reference=REFERENCE_AEROBIC_BATCH, uptake_mode="fixed"):
    """Audit the entire declared measurement set; host validation belongs to the paired caller."""
    tasks = measured_tasks(reference)
    env = reference_environment(reference)
    result = {"model_id": model.id, "condition_id": env.condition_id, "source": reference.source,
              "tasks": [task.name for task in tasks], "order": [], "uptake_mode": uptake_mode,
              "protein_budget": "native; not calibrated", "biological_validation": False}
    with model as scoped:
        constraints = [constrain_uptake(scoped, exchange, value, mode=uptake_mode, source=reference.source)
                       for exchange, value in (("r_1714", reference.glucose_uptake),
                                               ("r_1992", reference.oxygen_uptake))]
        result["uptake_constraints"] = [row.report() if row is not None else None for row in constraints]
        scoped.objective = BIOMASS
        scoped.objective_direction = "max"
        result["baseline_growth"] = growth_or_none(scoped)
        result["baseline_status"] = str(scoped.solver.status)
        try:
            order = task_priority_order(scoped, tasks, condition_id=env.condition_id, host=env.organism)
            applied = apply_task_bounds(scoped, order, biomass_reaction=BIOMASS)
            result.update(status=applied.status, order=[asdict(row) for row in order],
                          constrained_growth=applied.max_growth, refusal=None)
        except InfeasibleRegulation as exc:
            result.update(status=exc.status, constrained_growth=None, refusal=str(exc))
        except (ValueError, KeyError) as exc:
            result.update(status="refused", constrained_growth=None, refusal=str(exc))
    return result


def matched_task_comparison(plain, ec, *, strain_transfer=None, uptake_mode="fixed"):
    env = reference_environment()
    result = {"condition_id": env.condition_id, "organism": env.organism, "measured_strain": env.strain,
              "strain_transfer": strain_transfer, "uptake_mode": uptake_mode,
              "tasks": [task.name for task in measured_tasks()], "models": {},
              "protein_budget_calibration": "none; never automatically tuned", "biological_validation": False}
    try:
        identities = [model_identity(model) for model in (plain, ec)]
        require_matched_models(*identities)
        if "prot_pool_exchange" in plain.reactions or "prot_pool_exchange" not in ec.reactions:
            raise ValueError("comparison needs a plain model and an enzyme-constrained model with an explicit protein pool")
        for model, identity in zip((plain, ec), identities):
            validate_host(model, identity, env, strain_transfer=strain_transfer)
        result["identities"] = [asdict(identity) for identity in identities]
    except ValueError as exc:
        result.update(status="incompatible", refusal=str(exc))
        return result
    for label, model in (("plain", plain), ("ec", ec)):
        result["models"][label] = measured_task_report(model, uptake_mode=uptake_mode)
    result["status"] = "completed"
    result["interpretation"] = "matched condition/base network; a refused full task set remains refused, with biomass retained"
    return result


def _print_relative_width_summary(frame: pd.DataFrame) -> None:
    defined = (np.isfinite(frame.relative_width_regulated)
               & np.isfinite(frame.relative_width_unregulated))
    narrowed = (frame.loc[defined, "relative_width_regulated"]
                < frame.loc[defined, "relative_width_unregulated"] - 1e-9)
    print(f"\n  rows where regulation narrowed the relative width: {int(narrowed.sum())} "
          f"of {int(defined.sum())} defined comparisons; {int((~defined).sum())} undefined")


def main(argv=None) -> None:
    global OUT
    args = parse_args(argv)
    default_ec = (paths.ec_yeast_gem() if os.environ.get("YSTWIN_EC_YEAST_GEM") else
                  paths.data_dir() / "gem" / "ecYeastGEM_yeast902.xml.gz")
    selected_ec = args.ec_model if args.ec_model is not None else default_ec
    ec_path = paths.resolve_or_exit(
        selected_ec if selected_ec is not None and selected_ec.is_file() else None,
        "a matched GECKO model", "YSTWIN_EC_YEAST_GEM")
    plain_path = paths.resolve_or_exit(args.plain_model or paths.yeast_gem(), "the plain yeast model", "YSTWIN_YEAST_GEM")
    ec, settings = load_model(ec_path)
    plain, plain_settings = load_model(plain_path)
    OUT = args.output_dir.resolve()
    if OUT.exists() and (not OUT.is_dir() or any(OUT.iterdir())):
        raise SystemExit("output directory must be empty: refusing to overwrite or mix prior artifacts")
    OUT.mkdir(parents=True, exist_ok=True)
    comparison = matched_task_comparison(plain, ec, strain_transfer=args.strain_transfer,
                                         uptake_mode=args.uptake_mode)
    model_paths = {"plain": plain_path, "ec": ec_path}
    comparison["model_paths"] = {name: paths.display_path(path) for name, path in model_paths.items()}
    comparison["model_sha256"] = {}
    for name, path in model_paths.items():
        with path.open("rb") as stream:
            comparison["model_sha256"][name] = hashlib.file_digest(stream, "sha256").hexdigest()
    comparison["solvers"] = {"plain": str(plain_settings), "ec": str(settings)}
    (OUT / "regulation_matched_tasks.json").write_text(json.dumps(comparison, indent=2, allow_nan=False) + "\n")
    print(json.dumps(comparison, indent=2, allow_nan=False))
    if args.tasks_only or comparison["status"] == "incompatible":
        print(f"wrote {OUT / 'regulation_matched_tasks.json'}")
        return
    print(f"model: {ec_path.name}  solver: {settings}")
    model = add_beta_carotene_pathway(ec)
    for exchange, value in (("r_1714", REFERENCE_AEROBIC_BATCH.glucose_uptake),
                            ("r_1992", REFERENCE_AEROBIC_BATCH.oxygen_uptake)):
        constrain_uptake(model, exchange, value, mode=args.uptake_mode, source=REFERENCE_AEROBIC_BATCH.source)
    regulons = load_regulons()

    # ----------------------------------------------------------------- 1
    banner("D3.1  What the module -> gene -> reaction mapping actually reaches")
    print("module -> gene: SGD curated regulation records, restricted to MacIsaac 2006")
    print("               conserved binding motifs (PMID 16522208). Binding evidence only,")
    print("               so the map is independent of any expression dataset.")
    print("gene -> reaction: the GEM's own gene_reaction_rule. Not curated here.")
    present = {g.id for g in model.genes}
    rows = []
    for module in MODULES:
        regulon = regulons.get(module)
        if regulon is None:
            rows.append({"module": module, "factors": "", "regulon_genes": 0,
                         "genes_in_model": 0, "reactions": 0, "enzyme_draws": 0,
                         "note": "metabolite pool: no transcription factor"})
            continue
        genes = sorted(set(regulon.genes) & present)
        reactions = set()
        for gene in genes:
            reactions |= {r.id for r in model.genes.get_by_id(gene).reactions}
        draws = {r for r in reactions if r.startswith("draw_prot")}
        if genes:
            note = ""
        elif regulon.genes:
            note = f"{len(regulon.genes)}-gene regulon, none of them in this model"
        else:
            note = "factor has no target genes at the MacIsaac evidence level"
        rows.append({
            "module": module, "factors": "/".join(regulon.factors),
            "regulon_genes": len(regulon.genes), "genes_in_model": len(genes),
            "reactions": len(reactions), "enzyme_draws": len(draws), "note": note,
        })
    reach = pd.DataFrame(rows)
    print(reach.to_string(index=False))
    reach.to_csv(OUT / "regulation_module_reach.csv", index=False)
    covered = reach[reach.reactions > 0]
    print(f"\n  {len(covered)} of {len(MODULES)} modules reach at least one reaction; "
          f"{len(reach[reach.regulon_genes == 0])} have no regulon at all.")

    # ----------------------------------------------------------------- 2
    banner("D3.2  Which direction can this panel push? (signed weights, not absolute)")
    print("A negative weight in STRESSORS[...].targets means down-regulation, and")
    print("module_response carries the sign through. So this is a property of the panel:")
    signs = []
    for stressor, dose in PROBES:
        activity = module_response(stressor, dose)
        for module, value in activity.items():
            if value == 0.0:
                continue
            signs.append({
                "stressor": stressor, "module": module, "activity": round(value, 4),
                "direction": "induced" if value > 0 else "REPRESSED",
                "has_regulon": bool(regulons.get(module) and regulons[module].genes),
            })
    sign_table = pd.DataFrame(signs)
    print(sign_table.to_string(index=False))
    repressive = sign_table[(sign_table.direction == "REPRESSED") & sign_table.has_regulon]
    print(f"\n  repressed module-arms WITH a regulon: {len(repressive)}")
    if repressive.empty:
        print("  -> Every module this panel represses is a metabolite pool, which has no")
        print("     transcription factor and therefore no regulon. On today's panel the")
        print("     regulation layer can only RELAX capacity bounds, never tighten them.")
        print("     That is a finding about the panel, not a limitation of the mechanism.")

    # ----------------------------------------------------------------- 3
    banner(f"D3.3  E-Flux (PMID {COLIJN_PMID}) vs no regulation, at D1's growth fractions")
    print(f"  E-Flux lifts a lower bound off zero: {EFLUX_LIFTS_LOWER_BOUND}")
    print("  (it bounds |v|; a capacity constraint never requires flux)")
    results = []
    for stressor, dose in PROBES:
        activity = module_response(stressor, dose)
        layer = eflux_layer(model, activity, regulons, allow_unmapped=UNMAPPED)
        print(f"\n{stressor} @ {dose}: {layer.summary()}")
        if layer.modules_excluded:
            print(f"  excluded by declaration: {list(layer.modules_excluded)}")
        print(f"{'growth %':>9} {'min_off':>11} {'max_off':>11} {'relw_off':>9}"
              f" {'min_on':>11} {'max_on':>11} {'relw_on':>9} {'ceiling':>9}")
        for fraction in GROWTH_FRACTIONS:
            off = product_flux_range(model, PRODUCT_DEMAND_ID, glucose_uptake=None,
                                     growth_fraction=fraction)
            with model as scoped:
                report = apply_eflux(scoped, layer, biomass_reaction=BIOMASS,
                                    acknowledge_infeasible=True)
                on = product_flux_range(scoped, PRODUCT_DEMAND_ID, glucose_uptake=None,
                                        growth_fraction=fraction)
            # At 100% of maximum growth the ceiling is zero to solver tolerance, so the
            # ratio of two values at 1e-15 is noise and is printed as such rather than as
            # a 27-fold narrowing.
            degenerate = abs(off.maximum) < 1e-9
            shift = float("nan") if degenerate else on.maximum / off.maximum
            tag = "  (zero)" if degenerate else f"{shift:>8.3f}x"
            print(f"{fraction:>8.0%} {off.minimum:>11.3e} {off.maximum:>11.3e} "
                  f"{off.relative_width:>9.3f} {on.minimum:>11.3e} {on.maximum:>11.3e} "
                  f"{on.relative_width:>9.3f} {tag}")
            results.append({
                "stressor": stressor, "dose": dose, "growth_fraction": fraction,
                "bounds_changed": report.bounds_changed,
                "bounds_tightened": report.bounds_tightened,
                "min_unregulated": off.minimum, "max_unregulated": off.maximum,
                "relative_width_unregulated": off.relative_width,
                "min_regulated": on.minimum, "max_regulated": on.maximum,
                "relative_width_regulated": on.relative_width,
                "ceiling_ratio": shift,
                "max_growth_unregulated": off.max_growth,
                "max_growth_regulated": on.max_growth,
            })
    frame = pd.DataFrame(results)
    frame.to_csv(OUT / "regulation_product_range.csv", index=False)
    _print_relative_width_summary(frame)
    print(f"  every regulated floor at zero: "
          f"{bool((frame.min_regulated.abs() < 1e-9).all())}")

    # ----------------------------------------------------------------- 4
    banner("D3.4  COUNTERFACTUAL: what the same layer does when a module IS repressed")
    print("Not a result about this panel. The sign of every module activity is flipped, so")
    print("each regulon is repressed instead of induced, to show whether the mechanism")
    print("binds at all -- i.e. whether section 3's null is the panel or the layer.")
    print("\nRun against TWO reference capacities, because the choice of reference is the")
    print("load-bearing decision and a null under only the loose one proves nothing:")
    print("  model-implied  = prot_pool budget / enzyme MW, the model's own hard ceiling")
    print("                   on one enzyme holding the entire protein pool")
    print("  attainable     = each reaction's own FVA maximum at 90% growth, unregulated,")
    print("                   which is tighter by however much the pool is actually shared")
    base = product_flux_range(model, PRODUCT_DEMAND_ID, glucose_uptake=None,
                              growth_fraction=0.90)
    print(f"\n  unregulated reference: {base.summary()}")

    counter = []
    for stressor, dose in PROBES:
        flipped = {m: -v for m, v in module_response(stressor, dose).items()}
        loose = eflux_layer(model, flipped, regulons, allow_unmapped=UNMAPPED)
        attainable = _attainable_capacities(model, [s.reaction for s in loose.scales])
        tight = eflux_layer(model, flipped, regulons, allow_unmapped=UNMAPPED,
                           capacities=attainable)
        for layer in (loose, tight):
            with model as scoped:
                report = apply_eflux(scoped, layer, biomass_reaction=BIOMASS,
                                    acknowledge_infeasible=True)
                feasible = report.feasible and report.max_growth is not None and report.max_growth > 0
                rng = (product_flux_range(scoped, PRODUCT_DEMAND_ID, glucose_uptake=None,
                                          growth_fraction=0.90) if feasible else None)
            row = {
                "stressor": stressor, "reference": layer.reference,
                "tightened": len(layer.binding),
                "min_scale": min((s.scale for s in layer.scales), default=float("nan")),
                "max_growth": report.max_growth, "solver_status": report.status,
                "growth_ratio": report.max_growth / base.max_growth if report.max_growth is not None else float("nan"),
                "ceiling": rng.maximum if rng else float("nan"),
                "ceiling_ratio": (rng.maximum / base.maximum) if rng else float("nan"),
                "floor": rng.minimum if rng else float("nan"),
                "relative_width": rng.relative_width if rng else float("nan"),
            }
            counter.append(row)
            print(f"  {stressor:12s} {layer.reference:26s} tightened "
                  f"{row['tightened']:>4d}, min scale {row['min_scale']:.3f} -> growth "
                  f"{row['growth_ratio']:.4f}x, ceiling {row['ceiling_ratio']:.4f}x, "
                  f"floor {row['floor']:.3e}, rel width {row['relative_width']:.3f}")
    pd.DataFrame(counter).to_csv(OUT / "regulation_counterfactual.csv", index=False)

    print("\n  How deep would uniform repression of these reactions have to go to bite?")
    print("  Every reaction the H2O2 regulons reach, scaled by one factor, attainable")
    print("  reference. This turns 'no effect' into a number: the depth at which it starts.")
    flipped = {m: -v for m, v in module_response("H2O2", 0.15).items()}
    probe = eflux_layer(model, flipped, regulons, allow_unmapped=UNMAPPED)
    attainable = _attainable_capacities(model, [s.reaction for s in probe.scales])
    depth = []
    print(f"  {'scale':>7} {'growth':>10} {'growth x':>9} {'ceiling':>11} "
          f"{'ceiling x':>10} {'floor':>11} {'rel width':>10}")
    for scale in (1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.0):
        with model as scoped:
            for rid, capacity in attainable.items():
                scoped.reactions.get_by_id(rid).upper_bound = scale * capacity
            scoped.objective = BIOMASS
            growth = growth_or_none(scoped)  # nan > 1e-9 is False, but not by intent
            status = str(scoped.solver.status)
            ok = status == "optimal" and growth is not None and growth > 1e-9
            rng = (product_flux_range(scoped, PRODUCT_DEMAND_ID, glucose_uptake=None,
                                      growth_fraction=0.90) if ok else None)
        row = {
            "scale": scale, "solver_status": status,
            "max_growth": growth if growth is not None else float("nan"),
            "growth_ratio": (growth / base.max_growth) if growth is not None else float("nan"),
            "ceiling": rng.maximum if rng else float("nan"),
            "ceiling_ratio": (rng.maximum / base.maximum) if rng else float("nan"),
            "floor": rng.minimum if rng else float("nan"),
            "relative_width": rng.relative_width if rng else float("nan"),
        }
        depth.append(row)
        print(f"  {scale:>7.2f} {row['max_growth']:>10.5f} {row['growth_ratio']:>9.4f} "
              f"{row['ceiling']:>11.4e} {row['ceiling_ratio']:>10.4f} "
              f"{row['floor']:>11.3e} {row['relative_width']:>10.3f}")
    pd.DataFrame(depth).to_csv(OUT / "regulation_depth_sweep.csv", index=False)
    print("\n  Read solver_status before the flux columns: infeasible rows have no interval.")
    print("  In feasible rows here, the product floor stays at zero and relative width at")
    print("  1.000. A capacity limit by itself does not require product flux; infeasibility")
    print("  is not a zero-growth or zero-production observation.")

    # ----------------------------------------------------------------- 5
    banner(f"D3.5  Measured task-efficiency lower bounds (COSMIC-inspired, PMID {COSMIC_PMID})")
    print("Task efficiency = measured flux / individually maximised FBA flux, with both")
    print("nutrient rates from the SAME measured condition. Task floors leave upper bounds")
    print("unchanged. This is consistency with measurements, not a forecast.\n")
    print("COSMIC section 4.3 and the supplied CHO task tables now verify the source procedure;")
    print("the exact CHO model and complete fitted execution configuration remain unavailable.")
    print("See data/cosmic_sources.json; no exact-paper reproduction is claimed.")
    reference = REFERENCE_AEROBIC_BATCH
    tasks = measured_tasks(reference)
    measured_result = comparison.get("models", {}).get("ec")
    if measured_result is None or measured_result["status"] != "optimal":
        refusal = comparison.get("refusal") if measured_result is None else measured_result["refusal"]
        print(f"  REFUSED on the full measured set: {refusal}")
        print("  Biomass, ethanol and CO2 remain in the set; no fallback order is produced.")
        for line in COSMIC_SCOPED_OUT.splitlines():
            print("  " + line if line else "")
    else:
        with model as scoped:
            for exchange, value in (("r_1714", reference.glucose_uptake), ("r_1992", reference.oxygen_uptake)):
                constrain_uptake(scoped, exchange, value, mode=args.uptake_mode, source=reference.source)
            order = task_priority_order(scoped, tasks, condition_id=reference_environment().condition_id,
                                         host=reference_environment().organism)
            for priority in order:
                print("  " + priority.summary())
            pd.DataFrame([{**vars(priority), "condition_id": comparison["condition_id"],
                           "model_role": "enzyme_constrained", "model_id": ec.id,
                           "model_sha256": comparison["model_sha256"]["ec"]}
                          for priority in order]).to_csv(
                OUT / "regulation_matched_task_priority.csv", index=False)
            report = apply_task_bounds(scoped, order, biomass_reaction=BIOMASS)
            print("\n  " + report.summary())
            rng = product_flux_range(scoped, PRODUCT_DEMAND_ID, glucose_uptake=None, growth_fraction=0.90)
            print(f"  beta-carotene with those task lower bounds: {rng.summary()}")
            print("  Other-task floors can couple to a product only through this model's stoichiometry;")
            print("  they are not measurements of beta-carotene production.")

    print("\n  And the product itself:")
    product = MetabolicTask("beta-carotene", PRODUCT_DEMAND_ID)
    try:
        task_priority_order(model, [*tasks, product])
    except MissingMeasurement as exc:
        print(f"  REFUSED, as it must -- {exc}")

    print("\n  NOT A RESULT -- arithmetic showing an explicitly ASSUMED product floor.")
    print("  A lower bound is not an equality and does not force a collapsed interval:")
    for q_p in (1.0e-3, 5.0e-3, 1.0e-2):
        with model as scoped:
            for exchange, value in (("r_1714", reference.glucose_uptake), ("r_1992", reference.oxygen_uptake)):
                constrain_uptake(scoped, exchange, value, mode=args.uptake_mode, source=reference.source)
            scoped.reactions.get_by_id(PRODUCT_DEMAND_ID).lower_bound = q_p
            rng = product_flux_range(scoped, PRODUCT_DEMAND_ID, glucose_uptake=None,
                                     growth_fraction=0.90)
        print(f"    q_p >= {q_p:.1e} -> [{rng.minimum:.4e}, {rng.maximum:.4e}], "
              f"relative width {rng.relative_width:.3f}")
    print("  This floor was supplied, not predicted. Remaining width is real model freedom.")
    print("  State-specific allocation and training-only environmental classification are")
    print("  available through run_state_allocation.py with explicit sourced policies/data.")
    print("  For this project's G-e it means the same thing D1 meant: the product number")
    print("  has to come from outside the stoichiometry.")
    print("\n  THE MISSING MEASUREMENT, named precisely: the specific beta-carotene")
    print("  production rate q_p in mmol/gDW/h for a strain carrying crtE/crtYB/crtI,")
    print("  from a carbon and product balance over a known interval (content in mg/gDCW")
    print("  and the specific growth rate over the same interval suffice, via")
    print("  q_p = content/MW * mu at steady content; add its accumulation rate otherwise).")
    print("  No strain of OURS carries the pathway, so our numerator does not exist, and")
    print("  substituting a simulated flux would set the efficiency to 1.0 by construction")
    print("  and make the resulting floor an artefact.")
    print("  A PUBLISHED numerator does exist, and docs/PARKED.md says so in its own words:")
    print("  Elizondo 2025 ran three strains at two dilution rates with crt in them, and")
    print("  q_betacarotene for all six steady states is vendored here at")
    print("  data/carotenoid/elizondo2025_steady_states.tsv, in exactly the q_p = content/MW")
    print("  * mu form derived above. It is not used as the measured task in this section")
    print("  because it reports no q_O2, so it cannot meet this section's own requirement")
    print("  that both nutrient rates come from the same measured condition, and because it")
    print("  is this project's own pathway-flux calibration set (scripts/fit_pathway_flux.py)")
    print("  and so is not independent of what it would be used to check.")
    print(f"\nwrote regulation tables to {OUT}; task priorities exist only if the full set solved")


if __name__ == "__main__":
    main()
