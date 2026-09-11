"""Generate a dataset spanning the environment, not just the dose ladder.

Every simulated plate in this repository so far was one environment: exponential-phase glucose,
30 degrees, air, one pH. A stress model fitted there has never seen the growth rate move for a
reason that is not stress, and the two things a plate reader measures -- optical density and
reporter -- both move with growth. So the model has no way to tell "slowed" from "stressed", and
will read a cold incubator as a stressor.

This sweeps the environment axes the generator supports, crosses them with the dose ladder and
the nutrient axis, simulates the replicate wells of every cell, and writes **one row per well
carrying every condition that produced it**. A downstream fit can then be conditioned on the
environment instead of having to infer it.

Three things it deliberately does not do.

It does not extrapolate physiology. Carbon and oxygen fluxes come from van Hoek 1998's measured
chemostat curve (PMID 9797269), run between 0.025 and 0.40 /h. A cell whose growth rate falls
outside that raises; the refusal is caught and written into the ``physiology_refusal`` column,
and the flux columns are left empty. How many cells are refused is a result -- it says how much
of the environment space this generator can actually ground.

It does not simulate a plate nobody would run. An environment whose unstressed growth rate is
zero -- an anoxic ethanol culture, an incubator above the maximum growth temperature -- is
reported and skipped rather than integrated. A well of cells that never divide is not a weak
measurement, it is no measurement, and giving it a reporter trace would put a number where the
wet lab has none.

It does not quietly reproduce a nicer figure than the plates show. Well-to-well spread defaults
to ``analysis.power.DEFAULT_WELL_CV``, 0.052 measured from the zero-dose wells of the real
NewProtocol plates, rather than to the 0.04 that ``generator/plate.py`` asserts.

Usage: python scripts/run_environment_sweep.py   (--help lists every axis)
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from dataclasses import replace

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.analysis.power import DEFAULT_WELL_CV
from ystwin.generator.context import citrine_ph_response, context_growth_rate
from ystwin.generator.culture import (
    CRITICAL_GROWTH_RATE_PER_H,
    simulate_culture,
    substrate_limited_capacity,
)
from ystwin.generator.design import context_grid, environment_grid
from ystwin.generator.literature import parameters_from_literature
from ystwin.observation import ReporterOptics

DEFAULT_CONSTRUCTS = ("UPRE1", "UPRE2", "NativeYap1", "AlteredYap1")

OPTICS = ReporterOptics(
    gain=2.6e6, background=400.0, autofluorescence=900.0, inner_filter_coeff=None
)
"""The same optics ``generator/plate.py`` uses, so the sweep is comparable to those plates.

Every field is **asserted**. ``docs/CLAIM_BOUNDARY.md`` records autofluorescence as never having
been measured on any plate -- both reporter-free BY4741 controls were read in OD600 only -- and
it is 900 against a background of 400, so it is not a small term. It is also held constant per
gram here, which is a further assumption: autofluorescence in yeast rises as a culture leaves
exponential growth, so the constant is least right in exactly the contexts this sweep adds."""

GDCW_PER_OD = 0.42
"""Dry weight per OD unit. **Asserted**, carried over from ``generator/plate.py``.

Strain- and medium-specific; protocol P1 in ``docs/PROTOCOLS.md`` says to measure it by
filtering and drying a known volume. Until then it scales every biomass here by an unknown
constant, which is harmless for a growth *rate* and not harmless for a yield."""

OD_BLANK = 0.098
"""Medium-only absorbance, from the real plates. Held fixed across the sweep, which is a
simplification: a real blank drifts upward with evaporation and differs by well position."""


def _parse() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--carbon", nargs="+", default=["glucose", "galactose", "ethanol"])
    parser.add_argument("--phase", nargs="+", default=["exponential", "diauxic", "stationary"])
    parser.add_argument("--temperature", nargs="+", type=float, default=[25.0, 30.0, 37.0])
    parser.add_argument("--oxygen", nargs="+", type=float, default=[0.21, 0.05])
    parser.add_argument("--ph-cytosolic", nargs="+", type=float, default=[7.0, 6.5])
    parser.add_argument("--ph-medium", nargs="+", type=float, default=[5.5])
    parser.add_argument("--glucose", nargs="+", type=float, default=[20.0])
    parser.add_argument("--supplements", nargs="+", type=int, default=[0],
                        help="1 where the medium carries ergosterol and Tween 80, which is "
                             "what anaerobic growth requires at all (PMID 13034889)")
    parser.add_argument("--dose", nargs="+", type=float,
                        default=[0.0, 0.1, 0.25, 0.5, 1.0, 2.0])
    parser.add_argument("--nutrient", nargs="+", type=float, default=[1.0, 0.6])
    parser.add_argument("--construct", nargs="+", default=list(DEFAULT_CONSTRUCTS))
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--duration-h", type=float, default=4.14)
    parser.add_argument("--timepoints", type=int, default=25)
    parser.add_argument("--well-cv", type=float, default=DEFAULT_WELL_CV)
    parser.add_argument("--reader-cv", type=float, default=0.01)
    parser.add_argument("--k-mat", type=float, default=None,
                        help="chromophore maturation rate, 1/h. Unset means instantaneous, "
                             "which is an assumption and not a measurement; see "
                             "CultureParameters.k_mat")
    parser.add_argument("--inoculum-od", type=float, default=0.16)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=pathlib.Path, default=None)
    return parser.parse_args()


def _wells(condition, construct: str, params, args, rng) -> list[dict]:
    """Simulate the replicate wells of one design cell.

    The strain's ``mu_max`` is replaced by the environment's unstressed rate before integrating.
    That is not a fudge: ``simulate_culture`` applies the dose as ``mu_max * (1 - hill(dose))``,
    so substituting the context rate for ``mu_max`` makes the dose act as the same *fraction* of
    whatever the environment allows, which is exactly the composition ``design.py`` documents and
    marks asserted. The promoter parameters are untouched, so the reporter still sees the real
    dose.

    The carrying capacity is derived from the glucose in the well and the **measured** biomass
    yield at this growth rate wherever that rate is inside van Hoek's range, and falls back to
    the strain's asserted capacity where it is not. Which of the two was used is written into
    the row: they differ by more than a factor of two across the Crabtree switch, and a reader
    cannot tell from the numbers alone.
    """
    times = np.linspace(0.0, args.duration_h, args.timepoints)
    base = condition.to_row()
    ph_factor = citrine_ph_response(condition.context.ph_cytosolic)
    unstressed = context_growth_rate(condition.context)
    rows = []
    for replicate in range(args.replicates):
        start_od = args.inoculum_od * float(np.exp(rng.normal(0.0, args.well_cv)))
        biomass0 = start_od * GDCW_PER_OD
        if condition.physiology is None:
            capacity = params.carrying_capacity
            capacity_source = "asserted per strain"
        else:
            capacity = substrate_limited_capacity(
                biomass0, condition.context.glucose_g_per_L, condition.growth_rate)
            capacity_source = "measured yield, van Hoek 1998 PMID 9797269"

        strain = replace(params, mu_max=unstressed, carrying_capacity=capacity,
                         k_mat=args.k_mat)
        out = simulate_culture(times, condition.dose_mM, strain, biomass0,
                               nutrient_factor=condition.nutrient_factor)

        true_od = out["biomass"] / GDCW_PER_OD + OD_BLANK
        per_cell = OPTICS.autofluorescence + OPTICS.gain * out["reporter"] * ph_factor
        true_rfu = OPTICS.background + per_cell * out["biomass"]
        od = true_od * (1.0 + rng.normal(0.0, args.reader_cv, times.size))
        rfu = true_rfu * (1.0 + rng.normal(0.0, args.reader_cv, times.size))

        rows.append({
            **base,
            "construct": construct,
            "replicate": replicate,
            "inoculum_od": float(start_od),
            "carrying_capacity_g_per_L": float(capacity),
            "carrying_capacity_source": capacity_source,
            "citrine_ph_factor": float(ph_factor),
            "k_mat": args.k_mat,
            "well_cv": float(args.well_cv),
            "reader_cv": float(args.reader_cv),
            "above_critical_growth_rate": bool(
                condition.growth_rate > CRITICAL_GROWTH_RATE_PER_H),
            "od_final": float(od[-1]),
            "od_max": float(np.max(od)),
            "rfu_final": float(rfu[-1]),
            "specific_fluorescence_final": float(
                (rfu[-1] - OPTICS.background) / max(od[-1] - OD_BLANK, 1e-9)),
            "reporter_final": float(out["reporter"][-1]),
            "realised_growth_rate": float(out["growth_rate"][0]),
        })
    return rows


def main() -> None:
    args = _parse()
    rng = np.random.default_rng(args.seed)
    contexts = context_grid(
        carbon_sources=args.carbon,
        growth_phases=args.phase,
        temperatures_c=args.temperature,
        oxygen_fractions=args.oxygen,
        ph_medium=args.ph_medium,
        ph_cytosolic=args.ph_cytosolic,
        glucose_g_per_L=args.glucose,
        anaerobic_supplements=[bool(flag) for flag in args.supplements],
    )
    alive = [c for c in contexts if context_growth_rate(c) > 0.0]
    dead = len(contexts) - len(alive)

    rows: list[dict] = []
    for construct in args.construct:
        params = parameters_from_literature(construct)
        for condition in environment_grid(params, alive, args.dose, args.nutrient):
            rows.extend(_wells(condition, construct, params, args, rng))

    frame = pd.DataFrame(rows)
    out = args.out or (paths.outputs_dir() / "environment_sweep.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)

    refused = int((frame.physiology_refusal.astype(str) != "").sum())
    print(f"{len(frame)} wells: {len(alive)} environments x {len(args.dose)} doses x "
          f"{len(args.nutrient)} nutrient levels x {len(args.construct)} constructs x "
          f"{args.replicates} replicates")
    if dead:
        print(f"{dead} of {len(contexts)} environments support no growth at all and were not "
              "simulated")
    print(f"physiology grounded on {len(frame) - refused} rows, refused on {refused} "
          f"({refused / max(len(frame), 1):.0%}) whose growth rate is outside the measured "
          "chemostat range 0.025-0.40 /h")
    grounded = frame.fermentative.dropna()
    if len(grounded):
        print(f"of the grounded rows {float(np.mean(grounded.astype(bool))):.0%} are "
              f"respiro-fermentative, i.e. above {CRITICAL_GROWTH_RATE_PER_H} /h")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
