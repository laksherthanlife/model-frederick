"""What the maintenance constraint costs the GEM, swept over OXYGEN as well as glucose.

    python3 scripts/maintenance_scale.py

Writes ``outputs/maintenance_scale.csv``. Needs the Yeast9 GSMM and nothing else -- no
wet-lab data, no fitted parameter.

`bridge/latent_bridge.py` turns a stress reporter into an ATP maintenance bound. Its scale,
``_MAX_STRESS_MAINTENANCE = 6.5`` mmol ATP/gDW/h, is an envelope its own comment disowns.
`bridge/maintenance_calibration.py` carries the measurement it stood in for (PMID 40181231:
0.066 against 0.082 mmol glucose/gDW/h, GSR intact against absent). This script measures what
that difference DOES to growth, because none of it was checkable while it lived in prose.

WHY IT WAS REWRITTEN, and it is a correction rather than an extension. The first version
swept GLUCOSE ONLY and wrote ``oxygen_lower_bound = -1000.0`` into every row. Its headline
finding -- one slope,
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.004642 /h per mmol ATP/gDW/h, unchanged from glucose -1.0 to
-10.0 -- is true and was reported as though it were general. It is not. Every row of that
table secreted
<!-- audit:value table=outputs/maintenance_scale.csv column=ethanol_at_zero_ngam_mmol_per_gdcw_h row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.000 mmol/gDCW/h of ethanol: the invariance was measured across
the carbon supply while the cell stayed fully respiratory, which is the one axis that does
not move it. Oxygen is the axis that does.

WHAT THE SWEEP FINDS. At this repository's own validated operating point --
`fba/physiology.REFERENCE_AEROBIC_BATCH`, glucose 11.1 against oxygen 3.7 mmol/gDCW/h, where
overflow is ON and the model secretes
<!-- audit:value table=outputs/maintenance_scale.csv column=ethanol_at_zero_ngam_mmol_per_gdcw_h row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->15.923 mmol/gDCW/h -- a unit of NGAM costs
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.010778 /h per mmol ATP/gDW/h, which is
<!-- audit:value table=outputs/maintenance_scale.csv column=slope_ratio_to_oxygen_unlimited row="glucose_lower_bound_mmol_per_gdcw_h=-11.1;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->2.322x the
respiratory slope both bridge modules quoted. A fermenting cell buys ATP at a worse exchange
rate, so maintenance costs it more growth.

AND THE REGIME IS NOT SET BY THE BOUNDS ALONE -- THE NGAM BOUND IS PART OF IT. Raising
maintenance on a respiring cell eventually forces overflow, and the slope changes at exactly
that point. At glucose -1.0 against oxygen -3.7 the model secretes no ethanol at all up to
NGAM 7.2 and slopes at the respiratory 0.004642; by NGAM 10.0 it ferments and the secant has
risen to
<!-- audit:value table=outputs/maintenance_scale.csv column=slope_at_top_probe_per_h_per_mmol_atp row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.006105. That is the whole content of ``ngam_slope_spread_relative``, which reaches
<!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_spread_relative row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.510426 in the worst row, and it is why the table carries
``regime_at_top_probe`` beside ``regime``: a row where those two disagree is a row where a
single slope is describing two different metabolisms.

THREE CLAIMS THE OLD TABLE MADE THAT THIS ONE NARROWS.

1.  "Growth is EXACTLY linear in the NGAM bound." True wherever the regime does not change,
    where the secant spread is
    <!-- audit:value table=outputs/maintenance_scale.csv column=ngam_slope_spread_relative row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->0.0. Where it does change, one slope is wrong by more than
    the growth assay can resolve: ``slope_linearity_growth_error_per_h`` reaches
    <!-- audit:value table=outputs/maintenance_scale.csv column=slope_linearity_growth_error_per_h row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.048399 /h against a
    <!-- audit:value table=outputs/maintenance_scale.csv column=growth_assay_floor_per_h row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->0.0117 /h floor. Both columns are in
    the table so a reader can see which rows a linear reduction is licensed on.
2.  "The penalty ratio is condition-independent." It is
    <!-- audit:value table=outputs/maintenance_scale.csv column=penalty_ratio_asserted_over_measured row="glucose_lower_bound_mmol_per_gdcw_h=-1.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->4.39 wherever the response is
    linear, which is a property of the two NGAM values and not of the metabolism, and it
    reaches
    <!-- audit:value table=outputs/maintenance_scale.csv column=penalty_ratio_asserted_over_measured row="glucose_lower_bound_mmol_per_gdcw_h=-1.5;oxygen_lower_bound_mmol_per_gdcw_h=-3.7" -->8.491 where it is not. It is carried per row rather than once.
3.  "The no-growth threshold is a property of the carbon supply." Of the carbon supply AND
    the oxygen supply, and the second moves it further:
    <!-- audit:value table=outputs/maintenance_scale.csv column=no_growth_ngam_mmol_atp_per_gdw_h row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-1000.0" -->191.92 at glucose -10.0 with oxygen
    free, <!-- audit:value table=outputs/maintenance_scale.csv column=no_growth_ngam_mmol_atp_per_gdw_h row="glucose_lower_bound_mmol_per_gdcw_h=-10.0;oxygen_lower_bound_mmol_per_gdcw_h=-2.0" -->25.73 at the same glucose under an oxygen cap of 2.0.

THE SLOPE IS MODEL-DERIVED / IN_SILICO, not MEASURED. It is the gradient of a linear program
solved on `data/gem/yeast-GEM.xml.gz` under bounds this script chose. No culture was read to
obtain it, and the two modules that quote it now say so.

Every condition is written into the table beside its result, because a GEM number without its
glucose and oxygen bounds is not comparable to any other GEM number in this repository -- and
the previous version of this table is the proof of it, having carried the bound it never
varied in a column nobody read.
"""
from __future__ import annotations

import pathlib
import sys
import warnings

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from ystwin import paths
from ystwin.bridge.latent_bridge import _ACTIVITY_FULL_SCALE, _MAX_STRESS_MAINTENANCE
from ystwin.bridge.latent_bridge import _RESTING_MAINTENANCE
from ystwin.bridge.maintenance_calibration import MEASURED_MAINTENANCE, regime_label
from ystwin.bridge.stress_energetics import maintenance_from_reporter
from ystwin.fba.physiology import REFERENCE_AEROBIC_BATCH
from ystwin.generator.panel_experiment import MEASURED_GROWTH_RATE_SE

NGAM_REACTION = "r_4046"
GLUCOSE_EXCHANGE = "r_1714"
OXYGEN_EXCHANGE = "r_1992"
ETHANOL_EXCHANGE = "r_1761"

#: Glucose uptake caps, mmol/gDCW/h, as lower bounds. The first three are the sweep this
#: table has always run; the fourth is the repository's own validated operating point.
GLUCOSE_BOUNDS = (-1.0, -1.5, -10.0, -REFERENCE_AEROBIC_BATCH.glucose_uptake)

#: Oxygen uptake caps. ``-1000.0`` is "unlimited", which is what every row of the previous
#: version of this table silently held; the rest walk the cell into overflow. The second is
#: `REFERENCE_AEROBIC_BATCH`'s own measured oxygen uptake, so one row of the cross product is
#: exactly the phenotype `fba/physiology.py` validates against.
OXYGEN_BOUNDS = (-1000.0, -REFERENCE_AEROBIC_BATCH.oxygen_uptake, -2.0, -1.0)

#: NGAM values the secant slopes are taken at, spanning the measured full-stress total and
#: the asserted one so a change of regime between them cannot hide.
SLOPE_PROBES = (1.0, 2.0, 5.0, 7.2, 10.0)


def _solve(model, ngam):
    """Growth and ethanol flux at one NGAM bound; ``(None, None)`` where it is infeasible.

    The warning filter is narrow and deliberate: probing past the no-growth threshold is what
    the bisection below is FOR, so cobra's infeasible-status warning is an expected result
    here, not a defect.
    """
    with model:
        model.reactions.get_by_id(NGAM_REACTION).bounds = (ngam, 1000.0)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "Solver status is", UserWarning)
            solution = model.optimize()
        if solution.status != "optimal" or solution.objective_value is None:
            return None, None
        return solution.objective_value, float(solution.fluxes[ETHANOL_EXCHANGE])


def _growth(model, ngam):
    return _solve(model, ngam)[0]


def _no_growth_threshold(model):
    low, high = 0.0, 400.0
    while high - low > 0.005:
        mid = (low + high) / 2.0
        value = _growth(model, mid)
        low, high = (low, mid) if value is None or value < 1e-9 else (mid, high)
    return low


def _round(value, digits):
    """``round`` that passes ``None`` through, so a refused cell stays empty in the CSV."""
    return None if value is None else round(value, digits)


def _penalty_percent(model, base, ngam, threshold):
    """Growth lost at one NGAM bound, as a percentage, or ``None`` past the threshold.

    ``None`` rather than 100: past the threshold the LP is INFEASIBLE, which is a different
    statement from zero growth, and `audit_claims.py` refuses to resolve a marker against an
    empty cell. Writing 100 there would launder an infeasibility into a number.
    """
    if ngam >= threshold:
        return None
    return 100.0 * (base - _growth(model, ngam)) / base


def _row(model, glucose_bound, oxygen_bound, measured_total, asserted_total):
    """One (glucose, oxygen) operating point, solved end to end.

    The no-growth threshold is found FIRST because it decides which of the fixed probes are
    reachable here: at glucose -1.0 under a tight oxygen cap the cell cannot pay 10 mmol
    ATP/gDW/h at all, and a secant taken across an infeasible point is not a slope.
    """
    threshold = _no_growth_threshold(model)
    base, ethanol = _solve(model, 0.0)

    probes = tuple(probe for probe in SLOPE_PROBES if probe < threshold)
    solved = {probe: _solve(model, probe) for probe in probes}
    slopes = [(base - solved[probe][0]) / probe for probe in probes]
    slope = slopes[0] if slopes else None

    # What a SINGLE slope costs, in the units the growth assay reports. A relative spread
    # says the response is bent; this says whether the bend is visible.
    linearity_error = (max(abs((base - slope * probe) - solved[probe][0]) for probe in probes)
                       if slope is not None else None)
    spread = ((max(slopes) - min(slopes)) / max(slopes)) if len(slopes) > 1 else None

    top = probes[-1] if probes else None
    top_ethanol = solved[top][1] if top is not None else None
    top_slope = slopes[-1] if slopes else None

    penalty_measured = _penalty_percent(model, base, measured_total, threshold)
    penalty_asserted = _penalty_percent(model, base, asserted_total, threshold)
    ratio = (penalty_asserted / penalty_measured
             if None not in (penalty_measured, penalty_asserted) and penalty_measured
             else None)

    return {
        "glucose_lower_bound_mmol_per_gdcw_h": glucose_bound,
        "oxygen_lower_bound_mmol_per_gdcw_h": oxygen_bound,
        "regime": regime_label(ethanol),
        "growth_at_zero_ngam_per_h": round(base, 6),
        "ethanol_at_zero_ngam_mmol_per_gdcw_h": round(ethanol, 3),
        "ngam_slope_per_h_per_mmol_atp": _round(slope, 6),
        "ngam_slope_spread_relative": _round(spread, 6),
        "slope_linearity_growth_error_per_h": _round(linearity_error, 6),
        "growth_assay_floor_per_h": MEASURED_GROWTH_RATE_SE,
        "slope_probes_used": len(probes),
        "top_slope_probe_mmol_atp_per_gdw_h": top,
        "ethanol_at_top_probe_mmol_per_gdcw_h": _round(top_ethanol, 3),
        "regime_at_top_probe": None if top_ethanol is None else regime_label(top_ethanol),
        "slope_at_top_probe_per_h_per_mmol_atp": _round(top_slope, 6),
        "measured_total_ngam_mmol_atp_per_gdw_h": round(measured_total, 4),
        "asserted_total_ngam_mmol_atp_per_gdw_h": round(asserted_total, 4),
        "growth_penalty_measured_percent": _round(penalty_measured, 2),
        "growth_penalty_asserted_percent": _round(penalty_asserted, 2),
        "penalty_ratio_asserted_over_measured": _round(ratio, 3),
        "no_growth_ngam_mmol_atp_per_gdw_h": round(threshold, 2),
    }


def _add_ablation(frame, measured_total, asserted_total):
    """What using the oxygen-unlimited slope everywhere costs, per row, in /h.

    The ablation this rewrite exists for. Both modules quoted ONE slope, so this asks what
    that slope gets wrong at each operating point and compares it against the growth assay's
    own floor rather than against nothing.
    """
    reference = {row.glucose_lower_bound_mmol_per_gdcw_h: row.ngam_slope_per_h_per_mmol_atp
                 for row in frame.itertuples()
                 if row.oxygen_lower_bound_mmol_per_gdcw_h == -1000.0}
    respiratory = frame["glucose_lower_bound_mmol_per_gdcw_h"].map(reference)
    difference = (frame["ngam_slope_per_h_per_mmol_atp"] - respiratory).abs()

    frame["slope_ratio_to_oxygen_unlimited"] = (
        frame["ngam_slope_per_h_per_mmol_atp"] / respiratory).round(3)
    frame["ablation_growth_error_at_measured_per_h"] = (
        difference * measured_total).round(6)
    frame["ablation_growth_error_at_asserted_per_h"] = (
        difference * asserted_total).round(6)
    frame["ablation_over_growth_floor_at_asserted"] = (
        difference * asserted_total / MEASURED_GROWTH_RATE_SE).round(3)
    return frame


def main() -> int:
    path = paths.yeast_gem()
    if path is None:
        print("yeast-GEM not present; set YSTWIN_YEAST_GEM", file=sys.stderr)
        return 1

    import cobra

    model = cobra.io.read_sbml_model(str(path))
    glucose = model.reactions.get_by_id(GLUCOSE_EXCHANGE)
    oxygen = model.reactions.get_by_id(OXYGEN_EXCHANGE)

    measured_total = maintenance_from_reporter(
        _ACTIVITY_FULL_SCALE, _ACTIVITY_FULL_SCALE).ngam_high
    asserted_total = _RESTING_MAINTENANCE + _MAX_STRESS_MAINTENANCE

    rows = []
    for glucose_bound in GLUCOSE_BOUNDS:
        for oxygen_bound in OXYGEN_BOUNDS:
            glucose.lower_bound, oxygen.lower_bound = glucose_bound, oxygen_bound
            rows.append(_row(model, glucose_bound, oxygen_bound,
                             measured_total, asserted_total))

    frame = _add_ablation(pd.DataFrame(rows), measured_total, asserted_total)
    out = paths.outputs_dir() / "maintenance_scale.csv"
    frame.to_csv(out, index=False)

    intact = MEASURED_MAINTENANCE["gsr_intact"]
    print(f"measured resting  {intact.atp_interval()[0]:.3f}-{intact.atp_interval()[1]:.3f} "
          f"mmol ATP/gDW/h; measured full-stress total {measured_total:.3f}; "
          f"asserted {asserted_total:.1f}")
    print(frame.to_string(index=False))
    print(f"\nwrote {paths.display_path(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
