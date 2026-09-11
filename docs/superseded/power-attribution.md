# Superseded: attribution power by design

**Current result:** [`../FINDINGS.md`](../FINDINGS.md), "What the design has to be" —
7.5% / 16.7% / 52.5% / 96.7% for the dose-only design at 1.2x / 1.5x / 2.0x / 3.0x, pinned
cell by cell to `outputs/power_analysis.csv`.

> This line read *1% / 5% / 44% / 96%* until 2026-09-09, and `FINDINGS.md` carried the same
> four numbers in a second, unpinned table directly below its pinned one. Neither matched the
> committed artefact. That is the failure this directory exists to prevent — a superseded
> version kept alongside the live result eventually gets read as the live result — so the
> version table below is now the only copy, and the live document links here instead of
> repeating it.

## v1 — no within-plate variance at all

    dose only   0%  4%  46%  98%

`analysis/power.py` modelled a plate-level batch effect and nothing between wells on one
plate. That omission was real: the plate term is shared by every condition and largely
divides out of a within-plate dose slope, while the well term lands on it directly.

## v2 — well-level sigma guessed at 0.10

    dose only   1%  8%  32%  85%

Reported as a correction, with the assertion that published well-to-well spread in yeast
plate assays "is generally larger, not smaller" — **written without a reference, and
backwards.**

## Why v2 was wrong

Measured from the zero-dose wells of the two usable plates, three wells per construct per
plate differing only in position, the lognormal sigma of recovered activity has a median of
**0.052** over eight construct-by-plate groups (range 0.006-0.083). Kensy et al. 2009
report under 5% between replicate wells of one clone, which agrees with the measurement and
not with the guess.

At 0.052 the table lands within two points of where it started. So the missing variance
component was real and **the alarm about its size was not** — the intermediate report that
a threefold effect had fallen from 98% to 85% was an artefact of a constant roughly twice
too large. Replicates needed at 1.2x go from 6 to 7, not to 8-9.

The lesson is narrower than "don't guess": the guess was *stated as a guess* and still
propagated, because a guessed constant with an unreferenced justification beside it reads
like a sourced one.
