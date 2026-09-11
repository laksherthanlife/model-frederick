"""Score the environment layer against Kocharin's PHB chemostats. It fails, and cleanly.

    python3 scripts/score_phb_environment.py

Kocharin & Nielsen 2013 (PMID 23514405, PMC3610212, CC-BY) ran ELEVEN aerobic chemostat
steady states on ONE genotype: four dilution rates crossed with three carbon feeds, all
matched at 0.666 Cmol/L. One genotype across eleven environments is the design that makes
"flux is a strain constant" falsifiable instead of fitted, and it is the first dataset this
repository has that can test the environment layer at all.

WHAT THIS CAN AND CANNOT TEST. PHB carries no expression measurement of any kind, so the
flux law -- ``flux = alpha * entry_enzyme_expression`` -- cannot be run or scored here at
all. And the solver's prediction for a passthrough chain is ``content = flux / mu``, which
is the identity Kocharin used to report their own numbers, so reproducing it tests
arithmetic rather than biology. What is left is the one thing worth testing: with the
genotype held fixed, does anything in the environment predict the flux?

THE ANSWER IS NO, and the shape of the failure is more useful than the number:

  - Carbon source moves the flux 4.2x AT IDENTICAL DILUTION RATE. The model has no term
    that can express that, because environment reaches the product only through mu.
  - The growth exponent is not one law. Glucose gives q ~ mu^1.20, ethanol mu^0.04,
    the 1:2 mix mu^0.03. A single exponent cannot describe all three.
  - The measured flux itself moves <!-- audit:value table=outputs/phb_environment_score.csv column=q_span_over_all_states row="state_id=glc_D005" -->6.2309x
    across the eleven states, on ONE genotype. Whatever explains that is environmental, and
    this architecture has nowhere to put it.

**A GEM CLAIM WAS DELETED FROM THIS DOCSTRING ON 2026-09-04 RATHER THAN CORRECTED.** It read:
"Pinning growth at each D and maximising a demand on cytosolic acetyl-CoA gives a ceiling that
moves 1.3x across all eleven states while the measured flux moves 5.2x. Pearson against the
measurement is -0.016." Three numbers, and this script imports no cobra and computes none of
them -- it never has. Two things are wrong with them and only one is fixable here.

The fixable one: "the measured flux moves 5.2x across all eleven states" is not a span over
eleven states. Over all eleven, q_PHB moves 6.2309x; 5.1963x is the span with `etoh_D010` and
`mix_D010` removed -- a partial state set labelled as the whole. That number IS computable
from the vendored file, so it is now computed, written to
`outputs/phb_environment_score.csv`, and marked above.

The unfixable one, and the reason the rest was deleted instead of recomputed: 1.3x and -0.016
came from a CEILING formulation -- "how much product could the host supply?" -- that
`pathway/gem_environment.py` has since measured to be the wrong question and abandoned. Its
docstring records why in terms this paragraph cannot improve on: the chain's predicted flux
sits 12.8-195x below the ceiling so it never binds, and scaling by each state's measured
uptake gives a captured fraction spanning 216x. Recomputing 1.3x under a superseded
formulation would be reviving a retracted method to keep a sentence; the successor result --
the GEM asked for a YIELD instead, scored leave-one-carbon-source-out on these same eleven
states -- lives there and has cells.

That verdict is the same wall as everywhere else in this architecture: FBA returns CEILINGS,
the biology runs far below them, and a ceiling therefore carries no information about the
flux. It is why `fba/audit.py` audits instead of supplying.

READ TOGETHER WITH KOCHARIN 2012 (PMID 23009357), which supplies the other half. All four
strains there carry the same phaA on the same PGK1 promoter -- entry-enzyme expression is
1.0 across the panel by construction -- and PHB specific productivity still moves 16.5-fold
on host acetyl-CoA supply alone. So for this pathway the flux law's input does not vary and
the flux does. The beta-carotene result does not transfer here unexamined.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin import paths

STATES = "kocharin2013_chemostat_states.tsv"


def load() -> pd.DataFrame:
    path = paths.data_dir() / "phb" / STATES
    if not path.exists():
        raise SystemExit(f"{path} is tracked and missing; bad checkout")
    return pd.read_csv(path, sep="\t")


def rule(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def growth_exponents(states: pd.DataFrame) -> pd.DataFrame:
    """Fit q ~ mu^b within each feed. One law would give one b."""
    rows = []
    for feed, group in states.groupby("carbon_source"):
        slope, intercept = np.polyfit(np.log(group.dilution_rate_per_h),
                                      np.log(group.q_phb_mmol_per_gdcw_h), 1)
        rows.append({"carbon_source": feed, "exponent": slope,
                     "n_states": len(group)})
    return pd.DataFrame(rows)


def variance_explained(states: pd.DataFrame) -> dict[str, float]:
    """Fraction of log-flux variance a per-group mean removes. Not a model, a bound on one."""
    log_flux = np.log(states.q_phb_mmol_per_gdcw_h)
    total = float(np.var(log_flux))
    out = {}
    for key in ("carbon_source", "dilution_rate_per_h"):
        fitted = states.groupby(key).q_phb_mmol_per_gdcw_h.transform(
            lambda s: np.log(s).mean())
        out[key] = 1.0 - float(np.var(log_flux - fitted)) / total
    return out


def main() -> int:
    states = load()

    rule(f"Kocharin 2013: {len(states)} chemostat steady states, ONE genotype (SCKK006)")
    table = states.pivot_table(index="dilution_rate_per_h", columns="carbon_source",
                               values="q_phb_mmol_per_gdcw_h")
    print(table.to_string(float_format=lambda v: f"{v:.5f}"))
    print("\n  q_PHB, mmol/gDCW/h. Rows are dilution rate; columns are the feed.")

    rule("Carbon source moves the flux at IDENTICAL growth rate")
    low = table.loc[0.05]
    ratio = low.max() / low.min()
    print(f"  at D = 0.05 /h: {low.min():.5f} ({low.idxmin()})  ->  "
          f"{low.max():.5f} ({low.idxmax()})   = {ratio:.2f}x")
    print("\n  The model has no term that can express this. Environment reaches the product")
    print("  only through mu, and mu is identical across that row.")

    rule("And the growth-rate dependence is not one law")
    exponents = growth_exponents(states)
    for r in exponents.itertuples():
        print(f"  {r.carbon_source:<22} q ~ mu^{r.exponent:+.2f}   ({r.n_states} states)")
    span = exponents.exponent.max() - exponents.exponent.min()
    print(f"\n  spread {span:.2f} across three feeds. A single fitted exponent describes none")
    print("  of them, and the sign is not even stable: glucose rises, the other two are flat.")

    rule("How much of it is explainable at all")
    explained = variance_explained(states)
    for key, value in explained.items():
        print(f"  {key:<22} removes {value:>6.1%} of the variance in log(flux)")
    print("\n  Neither dominates, and together they are two grouping variables on eleven")
    print("  points -- which is a description of the data, not a model of it.")

    rule("The measured flux span, which is what this file can actually check")
    q = states.q_phb_mmol_per_gdcw_h
    span = float(q.max() / q.min())
    print(f"  q_PHB over all {len(states)} states: {q.min():.5f} -> {q.max():.5f}"
          f"  = {span:.4f}x")
    print("  This docstring said 5.2x until 2026-09-04. That is the span with etoh_D010 and")
    print("  mix_D010 removed -- a partial state set labelled as the whole -- and it had no")
    print("  cell to be checked against because nothing here computed it.")

    out = paths.outputs_dir() / "phb_environment_score.csv"
    merged = states.merge(exponents, on="carbon_source")
    # Written onto every row rather than one, because a marker addresses this file by
    # (column, row) and a scalar on one arbitrary row would make every other row's cell a
    # blank that reads as "not measured". It is a property of the whole table.
    merged["q_span_over_all_states"] = round(span, 4)
    merged.to_csv(out, index=False)
    print(f"\n  wrote {out}")

    rule("The verdict this script exists to record")
    print("  With the genotype held fixed and the environment varied -- which is the")
    print("  question the whole model is pitched at -- nothing in this architecture")
    print("  predicts the flux. That is a property of the model, not of the data:")
    print("  Kocharin's eleven states are exactly the experiment that would identify an")
    print("  environment term, and the architecture has nowhere to put one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
