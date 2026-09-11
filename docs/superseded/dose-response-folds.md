# Superseded: naive vs corrected dose-response folds

**Current result:** [`../DATA_INVENTORY.md`](../DATA_INVENTORY.md) — the plate-paired
column. All eight rows remain INCONCLUSIVE at n=2.

## v1 — ratio of pooled means

    AlteredYap1 @ 1.0 mM H2O2   1.43

The fold was computed as a ratio of means pooled across plates. A plate-level shift that
multiplies a plate's dosed and control wells alike **cancels in a per-plate ratio and does
not cancel in a ratio of means across plates** — so pooling imported batch variation into a
quantity meant to be free of it, and moved the point estimate rather than just its spread.

## v2 — plate-paired, per-cell inversion

    AlteredYap1 @ 1.0 mM H2O2   1.31

Folds computed within each plate and combined geometrically, matching the unit of analysis
`analysis/power.py` already used. Seven of the eight reported rows moved by at most 0.02;
this one moved by 0.12.

## Why v2 moved again

The activity itself came from the per-cell inversion, which needs a growth rate estimated
from a noisy optical trace. The total-signal form cancels `mu` algebraically, so that error
is removed rather than propagated:

    AlteredYap1 @ 1.0 mM H2O2   1.55

Seven rows again move by at most 0.02. The eighth — the one genuine induction in the table —
gets **larger**, which is the expected direction: the per-cell route was charging it for
noise in a quantity the total route never needs.

## What has not changed across any version

Every row is INCONCLUSIVE. `analysis/uncertainty.py::fold_change` refuses an interval below
three biological replicates and there are two. The point estimates have moved three times;
the reason none of them can be quoted as a result has not moved at all.
