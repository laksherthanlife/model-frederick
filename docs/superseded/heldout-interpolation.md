# Superseded: held-out interpolation skill

**Current result:** [`../WHY_INTERPOLATION_WORKS.md`](../WHY_INTERPOLATION_WORKS.md) — the
fit beats the baseline at every interior dose at or below 0.5 mM and loses above it, sign
separating on `dose / EC50` with the crossover between 3.7 and 7.4.

Three versions preceded it. Each was reported as a finding, and each was wrong in a
different way.

---

## v1 — the model loses on both splits

    extrapolation   model 504.2   nearest-dose 323.1   skill -0.561
    interpolation   model 374.0   nearest-dose 363.3   skill -0.029

**Claimed:** a fitted dose-response does not beat carrying the nearest measured dose
forward, on either split.

**Why it was wrong:** the model was a **log-linear** fit chosen for simplicity. The
measured ladders are biphasic, and `generator/panel_calibration.py` exists precisely
because a saturating Hill could not describe the fall at high dose. The right functional
form was already in the repository. Scoring the wrong one and reporting the loss as a
property of the problem was the error.

**What survived:** that "use the nearest dose you measured" is a strong baseline. It beat a
wrong functional form everywhere.

---

## v2 — the model wins on interpolation, skill +0.301

    interpolation   biphasic 294.5   nearest-dose 363.3   skill +0.189 (all fits)
                                                               +0.301 (identifiable only)

**Claimed:** the project's first predicted-then-measured win.

**Why it was wrong:** one partition, one seed. Seed 0 turned out to be the best of three.
Reporting a single partition as a result is the error the whole held-out exercise exists
to prevent, committed while building the thing that prevents it.

**What survived:** the biphasic form does beat the baseline somewhere, and filtering on
`DoseFit.identifiable` matters — the flag tracks the outcome exactly.

---

## v3 — skill +0.236 over three seeds, range spanning zero

    interpolation   skill +0.236   range [-0.073, +0.301]   3 partitions, one losing

**Claimed:** a real signal but not an established win, with the range crossing zero.

**Why it was wrong:** **confounded.** Those three seeded splits withheld two, three and
four doses respectively, so the partitions differed in *training size* as well as in which
dose was held out. A spread across them mixed the two effects and could not attribute
either.

**What survived:** the direction — that a single partition had overstated the result.

---

## Why v4 is different

There are only five interior rungs, so the partition space is **exhausted** rather than
sampled, one dose at a time. That makes the partitions comparable (equal training size) and
the result complete rather than a draw from it. And the structure that emerged — sign
separating on `dose / EC50` — has a mechanism behind it that predicts the direction before
the data is seen, rather than being read off it.

The pattern across all four versions is the same mistake in three costumes: **a number
reported before the thing that would have qualified it was in place.** First the functional
form, then replication, then comparability.

---

## v4's magnitudes were themselves at n=2

The exhaustive-partition result above was correct in structure and optimistic in size.
It was computed on two biological replicates, because a third plate was being dropped by
a parsing fault -- the export writes the plate down the page and the reader stopped after
21 wells of 87. With that plate read correctly:

| held out | v4 (n=2) | current (n=3) |
| ---: | ---: | ---: |
| 0.1 mM | +0.164 | +0.086 |
| 0.2 mM | +0.429 | +0.368 |
| 0.5 mM | +0.236 | +0.101 |
| 1.0 mM | −0.807 | −0.683 |
| 2.0 mM | −0.178 | −0.178 |

Nothing about the mechanism or the crossover changed. The lesson is narrower and worth
keeping: **an effect size measured at the smallest n the estimator will accept should be
expected to shrink**, and this one halved at two of the three doses where the model wins.
Every version in this file has moved in that direction.
