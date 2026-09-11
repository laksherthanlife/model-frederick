"""Commit the GEM's environment predictions for every declared pathway, before the bench test.

    python3 scripts/gem_environment_predictions.py

Writes ``outputs/gem_environment_predictions.csv``. Needs the Yeast9 GSMM and nothing else --
no wet-lab data, no fitted parameter, no calibration. That is the point: the descriptor is
stoichiometric, so it can be computed for a product with no data at all, which makes it a
PREDICTION rather than a fit.

WHAT IS BEING COMMITTED. `pathway/environment_flux.py` carries a MEASURED environment channel
for PHB: carbon source moves content 3.82x at an identical growth rate, 64 standard deviations,
one genotype and one vessel. Whether that transfers to a second product is the open question
that decides the project's generalisation claim, and a matched-growth-rate bench test is
planned to settle it. This script writes down what the GEM says the answer will be, first.

The prediction is that beta-carotene shares PHB's SIGN: ethanol beats glucose for GGPP as it
does for acetyl-CoA. At mu = 0.05 /h both are STABLE across the plausible carbon supply
(4-25 mmol C/gDCW/h): beta-carotene 1.236-1.950x, PHB 1.314-2.347x, neither crossing 1.0.

**A THIRD CLAIM WAS WITHDRAWN ON 2026-09-04: THE PREDICTIONS HOLD AT mu = 0.05 AND NOWHERE
ELSE ON THIS GRID.** The mu = 0.10 and mu = 0.15 rows shipped with ``is_a_prediction=True``
and should not have. `sign_is_stable` swept the carbon supply and ``continue``d past every
supply the model could not hold at that growth rate, so those verdicts were reached on 5 and
4 of 6 requested points -- and the surviving points are the ones nearest the feasibility
boundary, where both maxima collapse toward zero and their ratio diverges. PHB reads 3.0139x
at mu = 0.15 against 1.3584x at mu = 0.05 for that reason and not for a chemical one; its
glucose maximum has fallen from 1.537 to 0.147 mmol/gDCW/h. `SignStability` now carries
``points_requested``, ``points_feasible`` and ``min_growth_headroom``, all three of which are
written into the CSV, and `is_a_prediction` requires a complete sweep every point of which
sits at least `MIN_GROWTH_HEADROOM` below its own maximum growth rate. This is the SECOND
instance of the cf34701 finding in this module -- see `pathway/gem_environment.py`'s docstring
for why returning the sweep's provenance, and not just sweeping one more parameter, is the
repair.

**A SECOND CLAIM WAS WITHDRAWN ON 2026-09-03 AND SHOULD NOT COME BACK.** This script used to
add that "the transfer is NOT uniform -- glycogen and gadusol flip at low growth rate". That
was an artefact of the matched carbon supply shipping at 20.0 mmol C/gDCW/h, which is 4.7x
Kocharin's own measured uptake at mu = 0.05 (3.99, 4.52, 4.23). At the measured supply both
reverse, 0.86 to 1.63, and the pattern holds only for roughly 8 <= S <= 25 -- a window never
swept in the committed work. Those two rows are now reported as FORMULATION-DEPENDENT rather
than as predictions, and they agree with each other to 0.15% so they were one row, not two.

WHAT IT DOES NOT PREDICT, restated here because a committed number invites over-reading: the
MAGNITUDE (the GEM reaches about 30% of the measured log effect) and the REVERSAL with growth
rate (Kocharin's ethanol advantage shrinks and flips at mu = 0.1718 /h; the GEM's grows). The
second is the boundary between the two halves of the vision -- the GEM says what is possible,
and how much of it a cell takes is the part a learned regulatory state must supply.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.pathway.gem_environment import (
    DEFAULT_FORMULATION,
    MIN_GROWTH_HEADROOM,
    Formulation,
    InfeasibleState,
    sign_is_stable,
    yield_ratio,
)
from ystwin.pathway.spec import available_pathways, load_pathway

warnings.filterwarnings("ignore")

#: The four growth rates Kocharin measured, so the prediction is on the same grid as the one
#: measurement that can test it.
GROWTH_RATES = (0.05, 0.10, 0.15, 0.20)


#: The two matched-carbon supplies the module works at: the prediction default (Kocharin's own
#: measured median uptake) and the model-comparison supply that admits all eleven states.
FORMULATION_SUPPLIES_MMOL_C = (11.0, 20.0)


def _formulation_sensitivity(model, out_dir) -> None:
    """Both formulations of the same question, side by side, into a committed table.

    ADDED 2026-09-04, and it exists because `pathway/gem_environment.py` has quoted a
    formulation sensitivity since the day it was written -- "a bare precursor sink gives
    1.09-1.25 for PHB while a full PHB reaction carrying its NADPH cost gives 1.40-2.00" --
    that NOTHING IN THIS REPOSITORY COULD PRODUCE. `Formulation.precursor_sink_only` was
    declared, printed in `label()`, copied through `sign_is_stable` and used in the
    `comparable_with` equality contract, and never read by the computation, so both arms
    returned bit-identical numbers under different labels. The sink arm's half then went stale
    when the precursor-carrier fix landed, and the full arm's half could not go stale because
    it had never been computed.

    A pathway that has not declared its net cofactor balance is SKIPPED here and refused by
    `yield_ratio`, which is why only PHB has rows: beta-carotene's CrtI cofactor chemistry is
    not written down anywhere in this repository, and inventing it to fill a table is the
    failure this whole module is built against.
    """
    rows = []
    for name in available_pathways():
        spec = load_pathway(name)
        if not spec.full_pathway_stoichiometry:
            continue
        for supply in FORMULATION_SUPPLIES_MMOL_C:
            for growth_rate in GROWTH_RATES:
                measured = {}
                for sink_only in (True, False):
                    formulation = Formulation(
                        matched_carbon_cmol_per_gdcw_h=supply * 1e-3,
                        precursor_sink_only=sink_only)
                    try:
                        measured[sink_only] = yield_ratio(
                            model, spec, growth_rate, formulation).ratio
                    except InfeasibleState:
                        measured[sink_only] = None
                if measured[True] is None or measured[False] is None:
                    continue
                rows.append({
                    "product": name,
                    "supply_mmol_c_per_gdcw_h": supply,
                    "growth_rate_per_h": growth_rate,
                    "ratio_precursor_sink": round(measured[True], 4),
                    "ratio_full_pathway": round(measured[False], 4),
                    "full_over_sink": round(measured[False] / measured[True], 4)})
    if not rows:
        return
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "gem_environment_formulation.csv", index=False)
    print("\nformulation sensitivity: the same question asked two ways")
    print(frame.to_string(index=False))
    print("  a pathway with no declared net cofactor balance is absent from this table by "
          "refusal, not by omission -- see gem_environment.FullPathwayUndeclared")


def _score_against_kocharin(model, out_dir) -> None:
    """Leave-one-CARBON-SOURCE-out: can each descriptor predict a feed it has never seen?

    The bar that separates a stoichiometric descriptor from a fitted one. An empirical `z`
    fitted on two feeds must EXTRAPOLATE to a third; the GEM computes the third from
    stoichiometry and needs no fit at all. Both are scored here so the comparison is a number
    rather than an argument.
    """
    import math

    import numpy as np

    from ystwin.pathway.environment_flux import carbon_descriptor

    states = pathlib.Path(__file__).resolve().parents[1] / "data" / "phb" / \
        "kocharin2013_chemostat_states.tsv"
    if not states.exists():
        return
    frame = pd.read_csv(states, sep="\t")
    frame["z"] = [carbon_descriptor(g, e) for g, e
                  in zip(frame.feed_glucose_g_per_l, frame.feed_ethanol_g_per_l)]
    frame["log_content"] = np.log(frame.phb_mg_per_gdw)
    frame["log_mu"] = np.log(frame.mu_per_h)

    spec = load_pathway("phb")
    # SCORED AT A STATED SUPPLY THAT ADMITS ALL ELEVEN STATES, which is deliberately NOT the
    # prediction default. The two ask different questions: a PREDICTION should use the
    # measured uptake (11.0) and be gated on sign stability, while a MODEL COMPARISON must
    # keep every state or the laws are scored on different data and the numbers stop being
    # comparable. At 11.0 the mu = 0.20 states go infeasible and z-only appears to degrade
    # from 1.344x to 2.176x purely from losing them.
    scoring = Formulation(matched_carbon_cmol_per_gdcw_h=20.0e-3)
    gem = []
    for growth_rate, share in zip(frame.mu_per_h, frame.z):
        try:
            ratio = yield_ratio(model, spec, growth_rate, scoring).ratio
        except InfeasibleState:
            gem.append(np.nan)
        else:
            # Linear in the ethanol carbon share between the two pure feeds, which is the
            # interpolation `environment_descriptor` documents.
            gem.append(share * math.log(ratio))
    frame["gem"] = gem
    # The ABSOLUTE maximum as well as the contrast, because the comparison between them is
    # the finding: the absolute is 92% collinear with growth rate, so its apparent win over
    # mu-only is mostly mu. Scoring only the flattering one would hide that.
    absolute = []
    for growth_rate, share in zip(frame.mu_per_h, frame.z):
        try:
            result = yield_ratio(model, spec, growth_rate, scoring)
        except InfeasibleState:
            absolute.append(np.nan)
        else:
            absolute.append(math.log(
                result.glucose_max_mmol_per_gdcw_h * (1 - share)
                + result.ethanol_max_mmol_per_gdcw_h * share))
    frame["gem_absolute"] = absolute
    frame = frame.dropna(subset=["gem", "gem_absolute"])

    def leave_one_source_out(columns):
        design = np.column_stack([np.ones(len(frame))] + [frame[c].values for c in columns])
        target = frame.log_content.values
        groups = frame.carbon_source.values
        errors = []
        for held in np.unique(groups):
            mask = groups != held
            if mask.sum() < design.shape[1]:
                continue
            beta, *_ = np.linalg.lstsq(design[mask], target[mask], rcond=None)
            errors += list(abs(target[~mask] - design[~mask] @ beta))
        return math.exp(float(np.median(errors)))

    laws = [("mu only", ["log_mu"], 0),
            ("log GEM max", ["gem_absolute"], 0),
            ("GEM carbon contrast", ["gem"], 0),
            ("z empirical", ["z"], 1),
            ("mu + GEM", ["log_mu", "gem"], 0)]
    scored = pd.DataFrame([
        {"law": name, "fitted_environment_parameters": fitted,
         "leave_one_carbon_source_out_fold_error": round(leave_one_source_out(columns), 4)}
        for name, columns, fitted in laws])
    scored["corr_with_log_content"] = [
        round(float(np.corrcoef(frame[c[0]] if len(c) == 1 else frame[c[0]],
                                frame.log_content)[0, 1]), 4) if len(c) == 1 else None
        for _, c, _ in laws]
    scored["scoring_supply_mmol_c_per_gdcw_h"] = 20.0
    scored["corr_with_log_mu"] = [
        round(float(np.corrcoef(frame[c[0]], frame.log_mu)[0, 1]), 4) if len(c) == 1 else None
        for _, c, _ in laws]
    scored.to_csv(out_dir / "gem_environment_scores.csv", index=False)
    print("\nleave-one-CARBON-SOURCE-out on the Kocharin states "
          "(predicting a feed never seen):")
    print(scored.to_string(index=False))


def main() -> int:
    path = paths.yeast_gem()
    if path is None:
        print("yeast-GEM not present; set YSTWIN_YEAST_GEM", file=sys.stderr)
        return 1

    import cobra

    model = cobra.io.read_sbml_model(str(path))
    rows = []
    for name in available_pathways():
        spec = load_pathway(name)
        for growth_rate in GROWTH_RATES:
            row = {"product": name,
                   "precursor_metabolite": spec.precursor_metabolite,
                   "precursor_stoichiometry": spec.precursor_stoichiometry,
                   "growth_rate_per_h": growth_rate,
                   "formulation": DEFAULT_FORMULATION.label()}
            try:
                result = yield_ratio(model, spec, growth_rate, DEFAULT_FORMULATION)
            except InfeasibleState:
                row.update({"glucose_max_mmol_per_gdcw_h": None,
                            "ethanol_max_mmol_per_gdcw_h": None,
                            "ethanol_over_glucose": None, "favours_ethanol": None})
            else:
                sweep = sign_is_stable(model, spec, growth_rate)
                row.update({
                    "glucose_max_mmol_per_gdcw_h": round(
                        result.glucose_max_mmol_per_gdcw_h, 6),
                    "ethanol_max_mmol_per_gdcw_h": round(
                        result.ethanol_max_mmol_per_gdcw_h, 6),
                    "ethanol_over_glucose": round(result.ratio, 4),
                    "favours_ethanol": bool(result.favours_ethanol),
                    # A sign that reverses inside the plausible supply range is a property of
                    # the parameter, not of the chemistry, and is not reported as a prediction.
                    "sign_stable_over_supply_range": bool(sweep.stable),
                    "ratio_low": round(sweep.ratio_low, 4),
                    "ratio_high": round(sweep.ratio_high, 4),
                    # THE SWEEP'S OWN PROVENANCE, written into the artefact so the gate is
                    # visible in it rather than only in the code. Both columns were added on
                    # 2026-09-04 and both change committed verdicts -- see the module
                    # docstring's second cf34701 note.
                    "points_requested": sweep.points_requested,
                    "points_feasible": sweep.points_feasible,
                    "min_growth_headroom": round(sweep.min_growth_headroom, 4),
                    "min_growth_headroom_required": MIN_GROWTH_HEADROOM,
                    "is_a_prediction": bool(sweep.is_a_prediction())})
            rows.append(row)

    frame = pd.DataFrame(rows)
    out = paths.outputs_dir() / "gem_environment_predictions.csv"
    frame.to_csv(out, index=False)
    _formulation_sensitivity(model, out.parent)
    _score_against_kocharin(model, out.parent)

    pivot = frame.pivot(index="product", columns="growth_rate_per_h",
                        values="ethanol_over_glucose")
    print(f"ethanol / glucose precursor yield  [{DEFAULT_FORMULATION.label()}]")
    print("above 1.0 = ethanol is the better feed for that product\n")
    print(pivot.to_string())
    predicted = frame[frame.is_a_prediction.eq(True)]
    refused = frame[frame.is_a_prediction.eq(False)]

    def _rows(subset):
        return ", ".join(f"{r.product}@mu={r.growth_rate_per_h:g}"
                         for r in subset.itertuples()) or "none"

    print(f"\nCOMMITTED PREDICTIONS (sign stable over the WHOLE 4-25 mmol C/gDCW/h sweep, "
          f"every point at least {MIN_GROWTH_HEADROOM:.0%} below its own mu_max): "
          f"{_rows(predicted)}")
    print(f"NOT PREDICTIONS: {_rows(refused)}")
    print("  -- a row is refused if the sign reverses, if the sweep could not evaluate every "
          "requested supply, or if any point sat too close to the feasibility boundary, where "
          "the ratio diverges. The three reasons are separate columns.")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
