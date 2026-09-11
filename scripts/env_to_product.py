"""Retrospective PHB scoring and fixed-gene beta-carotene environment comparisons.

    python3 scripts/env_to_product.py

Writes ``env_to_product.csv`` under ``YSTWIN_OUTPUTS`` (default ``outputs``). Uses the
vendored Kocharin states, beta-carotene pathway definition and physiology table; this
empirical comparison does not request a GEM audit.

PHB compares mu-only and environment-response regressions leave-one-state-out on all eleven
Kocharin states. This is retrospective scoring within one source cohort, not independent
external validation. The shipped environment factor is reported separately, not scored as
though it were fitted without the held-out state.

BETA-CAROTENE uses a fitted flux law at fixed entry expression. Context audits reachability
of each requested held growth rate; the direct carbon-source response remains unmeasured.
Only reachable, in-range states return a mu-only content, with the predictor's environment
layer status. Known scientific refusals retain their requested rows with missing predictions
and explicit reasons. Equality among answered cases is not environment validation, and a
missing carbon-source case cannot establish invariance.

``mu_per_h`` is the requested setpoint for beta-carotene and the source-state rate for PHB.
``growth_rate_per_h`` is the returned beta-carotene rate, missing on refusals. Content columns
retain the existing mg/gDW dry-biomass convention; the predictor names that basis gDCW.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.generator.context import CultureContext
from ystwin.pathway import calibrations
from ystwin.pathway.environment_flux import carbon_descriptor, environment_factor
from ystwin.pathway.solve import ImplausibleContent
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, LayerState, SetpointUnreachable, predict_product

STATES = pathlib.Path(__file__).resolve().parents[1] / "data" / "phb" / \
    "kocharin2013_chemostat_states.tsv"


def _phb_rows() -> pd.DataFrame:
    """Retrospective leave-one-state-out PHB scoring within one Kocharin source cohort."""
    frame = pd.read_csv(STATES, sep="\t")
    frame["z"] = [carbon_descriptor(g, e) for g, e
                  in zip(frame.feed_glucose_g_per_l, frame.feed_ethanol_g_per_l)]
    frame["log_content"] = np.log(frame.phb_mg_per_gdw)
    frame["log_mu"] = np.log(frame.mu_per_h)
    frame["log_mu_z"] = frame.log_mu * frame.z

    def loo(columns):
        design = np.column_stack([np.ones(len(frame))] + [frame[c].values for c in columns])
        target = frame.log_content.values
        predicted = np.empty(len(frame))
        for index in range(len(frame)):
            mask = np.ones(len(frame), bool)
            mask[index] = False
            beta, *_ = np.linalg.lstsq(design[mask], target[mask], rcond=None)
            predicted[index] = design[index] @ beta
        return np.exp(predicted)

    frame["predicted_mu_only"] = loo(["log_mu"])
    frame["predicted_with_environment"] = loo(["log_mu", "z", "log_mu_z"])
    # The shipped law, applied forward rather than refitted: what the module would say.
    frame["shipped_law_fold_vs_glucose"] = [
        round(environment_factor("phb", z, mu), 4)
        for z, mu in zip(frame.z, frame.mu_per_h)]
    frame["prediction_scope"] = (
        "retrospective leave-one-state-out within one Kocharin 2013 source cohort; "
        "not independent external validation")
    return frame


def _carotene_rows() -> pd.DataFrame:
    """Retain every fixed-gene request, with a prediction or a typed scientific refusal."""
    spec = load_pathway("beta_carotene")
    genotype = Genotype(entry_expression=1.0)
    rows = []
    for carbon in ("glucose", "ethanol"):
        for growth_rate in (0.101, 0.15, 0.2543):
            environment = Environment(
                context=CultureContext(carbon_source=carbon),
                growth_rate_setpoint_per_h=growth_rate)
            row = {"product": "beta_carotene", "carbon_source": carbon,
                   "mu_per_h": growth_rate, "prediction_mode": "empirical",
                   "prediction_scope": "fixed-gene comparison, not environment validation"}
            try:
                got = predict_product(spec, genotype, environment,
                                      calibrations.BETA_CAROTENE_FLUX,
                                      calibrations.BETA_CAROTENE_KINETICS, mode="empirical")
            except (SetpointUnreachable, ImplausibleContent) as refusal:
                row.update({"refused": type(refusal).__name__,
                            "refusal_reason": str(refusal),
                            "growth_rate_per_h": None, "content_mg_per_gdw": None,
                            "environment_layer_state": LayerState.NOT_RUN,
                            "environment_sets_the_number": None})
            else:
                layer = got.layer("environment")
                row.update({"refused": "", "refusal_reason": "",
                            "growth_rate_per_h": got.growth_rate_per_h,
                            "content_mg_per_gdw": round(got.content_mg_per_gdcw, 6),
                            "environment_layer_state": layer.state,
                            "environment_sets_the_number": layer.sets_the_number})
            rows.append(row)
    return pd.DataFrame(rows).astype({"growth_rate_per_h": "float64",
                                      "content_mg_per_gdw": "float64",
                                      "environment_sets_the_number": "boolean"})


def main() -> int:
    phb = _phb_rows()
    carotene = _carotene_rows()

    out = paths.outputs_dir() / "env_to_product.csv"
    combined = pd.concat([
        phb.assign(product="phb")[[
            "product", "state_id", "carbon_source", "mu_per_h", "z",
            "phb_mg_per_gdw", "predicted_mu_only", "predicted_with_environment",
            "shipped_law_fold_vs_glucose", "prediction_scope"]].rename(
                columns={"phb_mg_per_gdw": "measured_mg_per_gdw"}),
        carotene], ignore_index=True)
    combined.to_csv(out, index=False)

    def fold_error(predicted, measured):
        return float(np.median(np.exp(np.abs(np.log(predicted / measured)))))

    mu_only = fold_error(phb.predicted_mu_only, phb.phb_mg_per_gdw)
    with_env = fold_error(phb.predicted_with_environment, phb.phb_mg_per_gdw)

    print(f"PHB -- {phb.prediction_scope.iloc[0]} ({len(phb)} states):")
    print(f"    mu only            {mu_only:.3f}x median fold error")
    print(f"    mu + environment   {with_env:.3f}x")
    print("    within this cohort, adding environment terms reduces median fold error "
          f"in excess of 1 by {100 * (mu_only - with_env) / (mu_only - 1):.0f}%\n")
    print(phb[["state_id", "mu_per_h", "carbon_source", "phb_mg_per_gdw",
               "predicted_mu_only", "predicted_with_environment"]]
          .round(3).to_string(index=False))

    print("\nbeta-carotene -- fixed-gene empirical comparison, not environment validation:")
    print(carotene[["carbon_source", "mu_per_h", "growth_rate_per_h", "content_mg_per_gdw",
                    "environment_layer_state", "refused"]].to_string(index=False))
    answered = carotene[carotene.refused.eq("")]
    refusals = carotene[carotene.refused.ne("")]
    print(f"    {len(answered)} answered, {len(refusals)} refused out of "
          f"{len(carotene)} requested cases")
    for rate, attempted in carotene.groupby("mu_per_h"):
        part = attempted[attempted.refused.eq("")]
        refused_count = len(attempted) - len(part)
        comparison = ("matched-carbon comparison unavailable" if refused_count else
                      f"{part.content_mg_per_gdw.nunique()} distinct answered content value(s)")
        print(f"    mu = {rate:g} /h: {len(part)} answered, {refused_count} refused; {comparison}")
    print("    Context audits held-setpoint reachability; the direct carbon response remains "
          "unmeasured. Equal answered contents do not establish biological invariance.")
    print("    Missing refused cases cannot establish carbon-source invariance.")
    for row in refusals.itertuples(index=False):
        print(f"    refusal ({row.carbon_source}, requested mu = {row.mu_per_h:g} /h): "
              f"{row.refused}: {row.refusal_reason}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
