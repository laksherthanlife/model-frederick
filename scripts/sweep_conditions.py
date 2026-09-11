"""Sweep dilution and expression for a named fixed empirical product comparison.

    python3 scripts/sweep_conditions.py
    python3 scripts/sweep_conditions.py --pathway beta_carotene --feed 20

Every cell calls ``predict_product`` with the shipped CrtE rate law and branch kinetics,
in the default culture environment. This exploratory fixed-gene comparison is not a causal
cassette-dosage law or nested validation of choosing CrtE and the model family.
Unsupported growth rates or pathway/calibration combinations retain an explicit refusal.

``--feed`` is used ONLY in the conditional biomass calculation, ``biomass_yield * feed``.
It is not passed as a changed product-model environment: that empirical API refuses feed
changes for which it has no calibrated route. Thus a feed sweep here rescales titre and
productivity under a stated biomass assumption; it does not predict how feed changes
expression, entry flux or pathway capacity.

``score_targets.py`` compares independently refitted target laws. Its ranking is not an
error penalty for converting a single rate prediction: relative and log errors remain
identical when the same measured growth, uptake and biomass factors are applied to both
prediction and observation. Those factors cancel exactly. Here biomass is estimated,
not measured, so that identity does not quantify biomass-model error or culture transfer.

Content bounds map the fixed entry law's held-out log spread through the monotone pathway
solver at fixed growth and kinetics. The endpoint mapping is exact for this conditional
calculation, not a calibrated prediction interval: it omits gene/model selection, input
measurement error, kinetic uncertainty and model misspecification. A fitted capacity is
not a universal genotype-independent ceiling.

Biomass yield and glucose uptake come from van Hoek 1998's glucose-limited chemostats,
whereas the Elizondo calibration strains were glucose-excess. Read titre as content scaled
by biomass assumed for a different culture, not as an independently validated flask assay.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.generator.context import CultureContext
from ystwin.generator.culture import chemostat_physiology
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, SetpointUnreachable, predict_product

#: Mean CrtE transcript levels across the two Elizondo 2025 steady states per strain.
#: These define fixed-expression comparisons, not expression predictions at each dilution.
STRAINS = {"b-car2": 0.25806685149999997,
           "b-car3": 0.58856938,
           "b-car4": 0.8220245785}

#: Glucose, g per mmol -- 180.16 g/mol.
_GLUCOSE_G_PER_MMOL = 0.18016

DILUTION_RATES = (0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.30)


def sweep(spec, feed_g_per_L: float) -> pd.DataFrame:
    """One fixed-model row per (dilution rate, strain), with conditional biomass scaling."""
    if not math.isfinite(feed_g_per_L) or feed_g_per_L <= 0:
        raise ValueError("feed_g_per_L must be positive and finite for the biomass calculation")
    context = CultureContext()
    rows = []
    for rate in DILUTION_RATES:
        for strain, expression in STRAINS.items():
            row = {"dilution_rate_per_h": rate, "strain": strain,
                   "entry_expression": expression,
                   "model": f"fixed_{BETA_CAROTENE_FLUX.entry_enzyme}_rate_law",
                   "prediction_scope": "exploratory fixed-gene empirical comparison; not nested validation",
                   "biomass_feed_g_per_L": feed_g_per_L,
                   "biomass_source": "van Hoek 1998 glucose-limited chemostats",
                   "biomass_scope": "conditional scaling only; no feed effect on product model"}
            try:
                prediction = predict_product(
                    spec, Genotype(entry_expression=expression),
                    Environment(context=context, growth_rate_setpoint_per_h=rate),
                    BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS)
            except SetpointUnreachable as exc:
                rows.append({**row, "answered": False, "refusal": f"washed out: {exc}"})
                continue
            except ValueError as exc:
                # First sentence, not first "." -- the message quotes growth rates, so
                # splitting on the period truncated it at "growth rate 0".
                message = " ".join(str(exc).split())
                rows.append({**row, "answered": False,
                             "refusal": message.split(". ")[0]})
                continue
            physiology = chemostat_physiology(rate)
            # Glucose-limited steady state with near-complete conversion: X = Y_sx * S_feed.
            biomass = physiology.biomass_yield * feed_g_per_L
            content = prediction.content_mg_per_gdcw
            # The BAND, not just the point. `content_interval` re-solves at the ends of the
            # flux calibration's leave-one-strain-out spread, which is exact here because the
            # solve is monotone in the entry flux. A titre quoted without it is a number a
            # reader has no way to discount, and this column existed unused until 2026-09-01.
            molar_mass = spec.nodes[-1].molar_mass_g_per_mol or 1.0
            low, high = prediction.content_interval(spec, BETA_CAROTENE_KINETICS)
            content_low, content_high = low * molar_mass, high * molar_mass
            rows.append({
                **row,
                "answered": True,
                "content_mg_per_gdcw": content,
                "content_low_mg_per_gdcw": content_low,
                "content_high_mg_per_gdcw": content_high,
                "biomass_g_per_L": biomass,
                "titre_mg_per_L": content * biomass,
                "titre_low_mg_per_L": content_low * biomass,
                "titre_high_mg_per_L": content_high * biomass,
                # What a fermenter is actually rated on: mass out per litre per hour.
                "productivity_mg_per_L_h": content * biomass * rate,
                "product_rate_mmol_per_gdcw_h": prediction.rate_mmol_per_gdcw_h,
                # Yield on substrate: mg product per gram of glucose taken up. Second-best
                # target in `score_targets.py`, and the one a process cares about after
                # titre -- it is what says whether the carbon is going where you want.
                # q_glucose is negative by convention (uptake), hence the abs().
                "yield_mg_per_g_glucose": (
                    prediction.rate_mmol_per_gdcw_h * molar_mass
                    / (abs(physiology.glucose_uptake) * _GLUCOSE_G_PER_MMOL)
                    if physiology.glucose_uptake else float("nan")),
                "biomass_yield_g_per_g": physiology.biomass_yield,
                "refusal": "",
            })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathway", default="beta_carotene")
    parser.add_argument("--feed", type=float, default=20.0,
                        help="glucose feed in g/L for conditional biomass scaling, not a product-model input")
    args = parser.parse_args()

    spec = load_pathway(args.pathway)
    table = sweep(spec, args.feed)
    answered = table[table.answered]

    print(f"\n{args.pathway} -- {len(answered)} of {len(table)} fixed-model comparisons answered "
          f"(biomass-scaling feed {args.feed:g} g/L glucose)\n")
    print("Exploratory fixed-CrtE rate law, not nested validation or a causal dosage law.")
    print("Feed rescales assumed biomass only; product-model culture inputs stay at their defaults.\n")
    if not answered.empty:
        show = answered.assign(
            content=lambda d: [f"{v:.3f} [{lo:.3f}, {hi:.3f}]" for v, lo, hi in zip(
                d.content_mg_per_gdcw, d.content_low_mg_per_gdcw, d.content_high_mg_per_gdcw)],
            titre=lambda d: [f"{v:.2f} [{lo:.2f}, {hi:.2f}]" for v, lo, hi in zip(
                d.titre_mg_per_L, d.titre_low_mg_per_L, d.titre_high_mg_per_L)],
        )[["dilution_rate_per_h", "strain", "content", "biomass_g_per_L", "titre",
           "productivity_mg_per_L_h"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    refused = table[~table.answered]
    if not refused.empty:
        print(f"\n{len(refused)} refused by the empirical model's scope or fitted range:")
        for rate, group in refused.groupby("dilution_rate_per_h"):
            print(f"  D = {rate:.2f} /h -- {group.iloc[0].refusal[:96]}")

    if not answered.empty:
        best = answered.loc[answered.productivity_mg_per_L_h.idxmax()]
        print(f"\nHighest conditional volumetric productivity in this grid: {best.strain} at D = "
              f"{best.dilution_rate_per_h:.2f} /h -> "
              f"{best.productivity_mg_per_L_h:.2f} mg/L/h "
              f"({best.titre_mg_per_L:.2f} mg/L at {best.biomass_g_per_L:.1f} gDCW/L)")
        print("\nThe fixed rate-law comparison makes content fall as growth dilutes the pool. "
              "This is a model consequence, not validation of a growth-rate mechanism.\n"
              "Titre and biomass use glucose-limited physiology for glucose-excess calibration "
              "strains; the bands exclude that transfer uncertainty and kinetic uncertainty.")

    # Each column's own provenance, derived rather than remembered. `rate` and `content`
    # rest only on constants fitted here; `yield` and `titre` each lean on van Hoek's
    # glucose-limited table, measured on a culture that is not this one.
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "parameter_provenance",
            pathlib.Path(__file__).resolve().parent / "parameter_provenance.py")
        provenance = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(provenance)
        print("\nwhat each column above rests on:")
        for quantity in ("rate", "content", "yield", "titre"):
            print(f"  {provenance.provenance_summary(quantity)}")
    except Exception as exc:                      # never let a footnote break a prediction
        print(f"\n(provenance summary unavailable: {type(exc).__name__})")

    out = paths.outputs_dir() / "condition_sweep.csv"
    table.to_csv(out, index=False)
    print(f"\n-> {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
