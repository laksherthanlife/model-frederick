"""Write a forward prediction down BEFORE the measurement that would score it.

    python3 scripts/register_prediction.py            # the standing D = 0.18 /h entry
    python3 scripts/register_prediction.py --rate 0.2 # register another one

`docs/CLAIM_BOUNDARY.md` opens with "Nothing in this repository is Tier 3", and the reason
is not that the predictions are bad. It is that every score in the repository was computed
after the measurement it scores. A prediction registered in advance can be wrong in public,
which is the only thing that separates Tier 3 from a well-argued Tier 1.

**Why 0.18 /h in particular.** The cyclase capacity is fitted as proportional to growth
rate over [0.101, 0.254] /h, and 0.18 sits near the middle of that interval with no
measurement in it.

Not because the law is unpinned -- an earlier version of this docstring said the exponent
was unidentifiable at two dilution rates and that was wrong. Three strains at those two
rates carry three different pathway fluxes, and a Michaelis-Menten cyclase responds to flux
with curvature, so the exponent fits freely to 1.04 and survives leave-one-strain-out. See
tests/test_growth_exponent_identifiability.py.

What 0.18 tests is the law's SHAPE BETWEEN the anchors rather than its exponent: every state
behind the fit sits at one end or the other, so the interior is interpolation that nothing
has checked. That is a smaller claim than the one this file used to make, and it is the one
a measurement here would actually settle.

**What goes in, and what used to.** The environment, and ONE GENOTYPE NUMBER: the strain's
relative CrtE expression. Nothing measured about the product enters.

That is a change. This script previously took the pathway flux as the mean of the strain's
own measured `q_lycopene + q_betacarotene`, which meant the registered prediction rested on
a measurement of the very quantity being predicted. It was registered in advance and it was
still circular. Both pools now come out of the chain and neither goes in.

The numbers therefore differ from the first registration. That is recorded rather than
hidden: the earlier file predicted from measured flux and this one predicts from genotype,
and only the second is a prediction in the sense the claim ladder means.

Re-running this reproduces the committed table exactly, which is the point: a registered
prediction nobody can regenerate is a number in a file rather than a commitment.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.pathway.calibrations import BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS
from ystwin.pathway.spec import load_pathway
from ystwin.predict import Environment, Genotype, predict_product

CALIBRATION = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"
REGISTERED_RATE = 0.18


def strain_expression(directory: pathlib.Path) -> pd.Series:
    """Relative CrtE expression per strain -- the genotype, and the only strain input.

    Averaged over the strain's states because expression from a constitutive cassette is a
    property of the construct. Where the two states disagree, that disagreement is inside
    the leave-one-strain-out error the calibration reports.
    """
    mrna = pd.read_csv(directory / "elizondo2025_relative_mrna.tsv", sep="\t")
    entry = mrna[mrna.gene == BETA_CAROTENE_FLUX.entry_enzyme]
    states = pd.read_csv(directory / "elizondo2025_steady_states.tsv", sep="\t")
    joined = states[["condition", "strain"]].merge(
        entry[["condition", "rel_expression"]], on="condition")
    return joined.groupby("strain").rel_expression.mean()


def register(rate: float) -> pd.DataFrame:
    directory = paths.data_dir() / "carotenoid"
    if not CALIBRATION.exists():
        raise SystemExit(f"{CALIBRATION} is missing; it is tracked, so this is a bad checkout")
    spec = load_pathway("beta_carotene")
    rows = []
    for strain, expression in strain_expression(directory).items():
        got = predict_product(spec, Genotype(float(expression), strain),
                              Environment(growth_rate_setpoint_per_h=rate),
                              BETA_CAROTENE_FLUX, BETA_CAROTENE_KINETICS,
                              mode="empirical")
        low, high = got.flux.interval()
        rows.append({
            "strain": strain,
            # THE REGISTERED SCHEMA IS FROZEN, INCLUDING THIS COLUMN NAME. The vessel this
            # prediction was registered for was a chemostat, so the column is
            # `dilution_rate` and stays `dilution_rate`. `Environment`'s field was renamed
            # to `growth_rate_setpoint_per_h` on 2026-09-02 when `fba/fedbatch.py` replaced
            # the chemostat, and renaming the column to match would edit a sealed record --
            # a pre-registration says what was predicted AND under what plan, and rewriting
            # it after the fact is the failure the registration exists to prevent.
            #
            # The prediction remains RUNNABLE under the new vessel, which is the thing worth
            # checking: 0.18 /h sits 29% below this strain's mu_max of 0.254 and opens at
            # q_S(0) = 0.18/0.098 = 1.84 g/gDCW/h against an uptake ceiling of 2.59, so both
            # of `design_fedbatch`'s refusals pass. See docs/PROTOCOLS.md P3.
            "dilution_rate": rate,
            "entry_gene": BETA_CAROTENE_FLUX.entry_enzyme,
            "entry_expression": float(expression),
            "predicted_pathway_flux": got.flux.flux_mmol_per_gdcw_h,
            "flux_low": low,
            "flux_high": high,
            "predicted_bcar_rate": got.rate_mmol_per_gdcw_h,
            "predicted_lycopene_content": got.intermediates()["lycopene"],
            "predicted_bcar_content_mg_per_gdcw": got.content_mg_per_gdcw,
        })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rate", type=float, default=REGISTERED_RATE,
                        help=f"dilution rate to register at (default {REGISTERED_RATE})")
    args = parser.parse_args()

    table = register(args.rate)
    # ``:g`` rather than a fixed width, so 0.18 gives D018 and not D0180 -- the name a
    # registered prediction is cited by must not change when the format string does.
    #
    # Written as one f-string ending in .csv rather than assembled with with_suffix(),
    # because scripts/audit_reproducibility.py reads destinations out of the syntax tree
    # and a name built in two steps resolves to no filename at all. This file then reports
    # as having no writer, which is exactly the complaint that produced this script.
    digits = f"{args.rate:g}".replace(".", "")
    out = paths.outputs_dir() / f"registered_prediction_D{digits}.csv"
    table.to_csv(out, index=False)

    print("=" * 88)
    print(f"Registered at D = {args.rate:g} /h, before any measurement exists there")
    print("=" * 88)
    print(f"  {'strain':8s} {'CrtE in':>9s} {'flux out':>11s} {'band':>26s} "
          f"{'rate out':>11s} {'content':>10s}")
    for r in table.itertuples():
        print(f"  {r.strain:8s} {r.entry_expression:9.3f} {r.predicted_pathway_flux:11.4e} "
              f"[{r.flux_low:10.4e}, {r.flux_high:10.4e}] {r.predicted_bcar_rate:11.4e} "
              f"{r.predicted_bcar_content_mg_per_gdcw:8.3f} mg/g")
    print(f"\n  The band is the flux law's leave-one-strain-out error, "
          f"{BETA_CAROTENE_FLUX.typical_fold_error:.2f}x, which is")
    print("  the spread a user of this number actually faces rather than an interval on a")
    print("  fitted constant. Nothing about the product of these strains went in.")
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
