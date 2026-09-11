"""Train the stress state at scale, and save the model that gets applied to real plates.

Everything upstream asks whether a general stress state can be learned. This fits one and
writes it out. Two models are trained: a wide one over every promoter fusion, which says
what the landscape supports, and the five-channel build that can actually be made, which
is the one a plate will be read with.

Scoring is leave-one-stressor-out in module units -- trained without a stressor, does the
model report the right biology when it arrives? Both the latent basis and the readout are
blind to the held-out stressor.

Usage: python scripts/run_training.py [--quick]
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.experiment_design import RECOMMENDED_DESIGN
from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD, THREE_SENSOR_BUILD
from ystwin.analysis.stress_model import (
    module_transfer,
    select_width_for_modules,
    train_stress_model,
)
from ystwin.analysis.transfer import select_dimension_by_transfer
from ystwin.generator.panel_experiment import (
    MEASURED_ACTIVITY_CV,
    MEASURED_GROWTH_RATE_SE,
    panel_dataset,
    pool_replicates,
    with_growth_channel,
)
from ystwin.generator.context import CultureContext
from ystwin.generator.panel_gates import accept_panel
from ystwin.generator.stress_panel import MODULES, REPORTERS, STRESSORS, transcriptional_reporters

OUT = paths.outputs_dir()
WIDE = transcriptional_reporters() + ["roGFP2-Grx1"]
MEASURED = OUT / "sensor_characterisation.csv"
CONTEXTS = [
    CultureContext(),
    CultureContext(growth_phase="diauxic"),
    CultureContext(growth_phase="stationary"),
    CultureContext(carbon_source="galactose"),
    CultureContext(carbon_source="ethanol"),
    CultureContext(oxygen=0.03),
    CultureContext(temperature_c=37.0),
]
"""The conditions a plate is actually read in. Everything before this trained in one corner
-- exponential glucose at 30 degrees -- where a raised general stress response always means
a stressor, which on a late plate it does not."""


def reference_readings():
    """Real dilution-corrected activity, for the posterior-predictive gate to compare with."""
    if not MEASURED.exists():
        return None
    real = pd.read_csv(MEASURED)
    return real[real.construct.isin(["UPRE1", "UPRE2"])].pivot_table(
        index=["plate", "dose_mM"], columns="construct",
        values="activity_late").dropna().to_numpy()


def gate(data, reference, tag):
    """Refuse to train on a dataset that could not have come off a plate.

    The posterior-predictive check runs on a matched slice, not the whole panel. The
    reference is one promoter under one agent; a twenty-channel panel across twenty-five
    stressors legitimately spans far more than that, and holding it to the spread of a
    single ladder would refuse it for covering more biology.
    """
    matched = panel_dataset(
        reporters=["UPRE-ER"], stressors=["DTT"], replicates=8,
        noise_cv=MEASURED_ACTIVITY_CV, growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=0)
    physical = accept_panel(data.readings, reporters=data.reporters)
    predictive = accept_panel(matched.readings, reference=reference,
                              reporters=matched.reporters)
    for report, what in ((physical, "physical"), (predictive, "matched slice")):
        verdict = "accepted" if report.accepted else f"REFUSED by {report.refused_by}"
        print(f"  gates on {tag} ({what}): {verdict}")
        if not report.accepted:
            raise SystemExit(f"training data refused: {report.reason}")


def deployable(data):
    """Turn a plate into what the model is actually fitted on.

    Two things a plate gives for free and the fit was throwing away. Replicates have to be
    pooled first, because the latent model gives every well its own state and extra wells
    otherwise buy parameters instead of averaging noise. And growth rate comes from optical
    density, which separates culture context from stress at no cost in fluorophores.
    """
    return pool_replicates(with_growth_channel(data))


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


DATASET_SEEDS = (0, 1, 2)
"""Dataset seeds to replicate over.

The seed here draws a different simulated experiment, not just different simulation
noise, so a result that moves across these is a result about one draw. Everything
below reports the spread alongside the value.
"""


def build(reporters, doses, replicates, seed=0, contexts=CONTEXTS):
    return panel_dataset(reporters=reporters, stressors=list(STRESSORS), doses=doses,
                         replicates=replicates, noise_cv=MEASURED_ACTIVITY_CV,
                         growth_rate_se=MEASURED_GROWTH_RATE_SE, seed=seed,
                         combinations=RECOMMENDED_DESIGN, contexts=contexts)


def train_and_score(data, reporters, label, reference, doses, replicates, seed=0):
    """Gate the raw plate, then fit what the plate supports."""
    gate(data, reference, label)
    fitted = deployable(data)
    k = select_width_for_modules(fitted, max_states=fitted.readings.shape[1])
    model = train_stress_model(fitted, n_states=k)
    reached = sum(v > 0.25 for v in model.recovery.values())
    print(f"  latent width {k} of {fitted.readings.shape[1]} channels "
          f"({len(reporters)} reporters + growth); recovers {reached}/{len(MODULES)} in-sample")
    return model, transfer_table(fitted, k, "leave-one-stressor-out")


def replicated_transfer(reporters, doses, replicates, label, reference,
                        seeds=DATASET_SEEDS):
    """Train and score once per dataset seed, and report the spread.

    Returns the per-seed tables and a summary. A single seed cannot distinguish a
    real margin from a lucky draw, which is the whole point.
    """
    tables, widths = [], []
    for seed in seeds:
        data = build(reporters, doses, replicates, seed=seed)
        gate(data, reference, f"{label} (seed {seed})")
        fitted = deployable(data)
        k = select_width_for_modules(fitted, max_states=fitted.readings.shape[1])
        widths.append(k)
        table = transfer_table(fitted, k, f"{label} seed {seed}", quiet=True)
        table["seed"] = seed
        tables.append(table)

    combined = pd.concat(tables, ignore_index=True)
    per_seed = combined.groupby("seed").modules_recovered.median()
    scores = combined.groupby("seed").mean_over_driven.mean()
    print(f"  {label}: latent width {widths}")
    print(f"    modules recovered, median per seed : {list(per_seed.round(1))}"
          f"  -> {per_seed.median():.1f} [{per_seed.min():.0f}, {per_seed.max():.0f}]")
    print(f"    mean score on driven modules       : {list(scores.round(3))}"
          f"  -> {scores.median():.3f} [{scores.min():.3f}, {scores.max():.3f}]")
    if len(set(widths)) > 1:
        print(f"    !! latent width is not stable across seeds {widths}; "
              "the chosen width is a property of the draw")
    return combined


def transfer_table(data, n_states, tag, quiet=False):
    rows = []
    for name in sorted(STRESSORS):
        scored = module_transfer(data, name, n_states)
        reached = [m for m, v in scored.items() if v > 0.25]
        rows.append({"held_out": name, "modules_recovered": len(reached),
                     "mean_over_driven": float(np.mean([scored[m] for m in scored
                                                        if m in STRESSORS[name].targets]))})
    frame = pd.DataFrame(rows)
    if not quiet:
        print(f"  {tag}: median {frame['modules_recovered'].median():.0f} modules recovered "
              f"per held-out stressor, mean score on the modules it drives "
              f"{frame['mean_over_driven'].mean():.3f}")
    return frame


def main(quick: bool) -> None:
    doses = (0.25, 0.5, 1.0, 2.0) if quick else (0.125, 0.25, 0.5, 1.0, 2.0, 4.0)
    replicates = 3 if quick else 6

    rule("training set")
    reference = reference_readings()
    wide = build(WIDE, doses, replicates)
    gate(wide, reference, "the wide panel")
    print(f"  {len(wide.labels)} wells, {len(set(wide.labels))} treatments, "
          f"{len(WIDE)} channels, {len(STRESSORS)} stressors, {len(MODULES)} state variables")
    print(f"  {len(set(wide.contexts))} culture contexts: {sorted(set(wide.contexts))}")

    fitted_wide = deployable(wide)
    k = select_width_for_modules(fitted_wide, max_states=8)
    channel_k = select_dimension_by_transfer(fitted_wide, max_states=6, observed=tuple(range(10)))
    print(f"  latent width chosen on module recovery: {k}")
    print(f"  (choosing on channel prediction instead would say {channel_k}, a different task)")

    rule("wide model: every promoter fusion plus one ratiometric sensor")
    model = train_stress_model(fitted_wide, n_states=k)
    recovered = {m: v for m, v in model.recovery.items() if v > 0.25}
    print(f"  recovers {len(recovered)}/{len(MODULES)} state variables in-sample")
    print("  " + ", ".join(f"{m} {v:.2f}" for m, v in
                           sorted(recovered.items(), key=lambda kv: -kv[1])[:12]))
    wide_scores = transfer_table(fitted_wide, k, "leave-one-stressor-out")

    rule("the five-channel build that can actually be made")
    reporters = RECOMMENDED_BUILD.stress_reporters
    kinds = ", ".join(f"{r} ({REPORTERS[r].kind.value[:5]})" for r in reporters)
    print(f"  {kinds}")
    print(f"  + {RECOMMENDED_BUILD.reference_fluorophore} reference, buildable: "
          f"{RECOMMENDED_BUILD.buildable}")
    deployed, built_scores = train_and_score(
        build(reporters, doses, replicates), reporters, "the five-channel build",
        reference, doses, replicates)
    # Printed, not written: the per-seed spread has no column in any tracked table.
    replicated_transfer(reporters, doses, replicates,
                        "five-channel across seeds", reference)

    rule("the three-sensor build: one channel per stress axis")
    kinds = ", ".join(f"{r} -> {REPORTERS[r].module}" for r in THREE_SENSOR_BUILD)
    print(f"  {kinds}")
    per_strain = replicates * len(RECOMMENDED_BUILD.stress_reporters) // len(THREE_SENSOR_BUILD)
    print(f"  {per_strain} replicates, not {replicates}: a plate holds a fixed number of "
          f"wells, so fewer strains means more wells each")
    minimal, minimal_scores = train_and_score(
        build(THREE_SENSOR_BUILD, doses, per_strain), THREE_SENSOR_BUILD,
        "the three-sensor build", reference, doses, per_strain)
    replicated_transfer(THREE_SENSOR_BUILD, doses, per_strain,
                        "three-sensor across seeds", reference)
    direct = {REPORTERS[r].module for r in THREE_SENSOR_BUILD}
    indirect = sorted(m for m, v in minimal.recovery.items() if v > 0.25 and m not in direct)
    print(f"  reads {len(direct)} modules, recovers {len(indirect)} more it never sees: "
          f"{', '.join(indirect)}")

    rule("does it read an unstressed late plate as unstressed?")
    for phase in ("exponential", "diauxic", "stationary"):
        quiet = panel_dataset(reporters=WIDE, stressors=["DTT"], doses=(0.0,), replicates=6,
                              noise_cv=MEASURED_ACTIVITY_CV, seed=7,
                              contexts=[CultureContext(growth_phase=phase)])
        said = float(model.predict_modules(deployable(quiet).readings)["ESR"].mean())
        real = float(quiet.modules_frame()["ESR"].mean())
        print(f"  {phase:<12} true ESR {real:.2f}   reported {said:5.2f}   error {abs(said-real):.2f}")

    rule("what the deployable model reports, module by module")
    table = pd.DataFrame({
        "module": list(MODULES),
        "wide_in_sample": [model.recovery[m] for m in MODULES],
        "built_in_sample": [deployed.recovery[m] for m in MODULES],
        "three_sensor_in_sample": [minimal.recovery[m] for m in MODULES],
    })
    print(table[table[["wide_in_sample", "built_in_sample",
                       "three_sensor_in_sample"]].max(axis=1) > 0.01]
          .round(3).to_string(index=False))

    model.save(OUT / "stress_model_wide.npz")
    deployed.save(OUT / "stress_model_build.npz")
    minimal.save(OUT / "stress_model_three_sensor.npz")
    # Stamp the run mode into every table. --quick uses four doses and three replicates
    # against six and six, which moves the numbers enough to change conclusions -- and a
    # table with no mode recorded was quoted as a full run once already.
    mode = "quick" if quick else "full"
    for frame, name in ((wide_scores, "training_transfer_wide"),
                        (built_scores, "training_transfer_build"),
                        (minimal_scores, "training_transfer_three_sensor"),
                        (table, "training_recovery")):
        stamped = frame.copy()
        stamped["run_mode"] = mode
        stamped["n_doses"] = len(doses)
        stamped["replicates"] = replicates
        stamped.to_csv(OUT / f"{name}.csv", index=False)
    print(f"\n  tables written in {mode} mode "
          f"({len(doses)} doses, {replicates} replicates)")
    print(f"\nwrote models and tables to {OUT}")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
