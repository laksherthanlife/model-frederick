# Superseded results

Every number here was reported, then replaced. It is kept because **how a result was
wrong is usually more informative than the result** — a reader deciding whether to trust
the current figure wants to know what the last three were and why each moved.

Keeping them inline was making the live documents unreadable: `NULL_RESULTS.md` had grown
three nested retractions of one claim. Keeping them in git alone is not enough either,
since nobody reads a commit log to check whether a number has a history.

## The rule

A live document states **only the current result** and links here. When a result is
replaced, the superseded version moves into a file below with three things: what was
claimed, why it was wrong, and what replaced it.

Nothing here is deleted. A superseded result that turns out to have been right is a
finding of its own, and it cannot be one if the record is gone.

## What is here

| file | claim | why it moved |
| --- | --- | --- |
| [`heldout-interpolation.md`](heldout-interpolation.md) | the held-out prediction skill | three versions: one lucky seed, one confounded spread, then exhaustive |
| [`power-attribution.md`](power-attribution.md) | attribution power by design | a guessed within-plate variance, twice too large |
| [`dose-response-folds.md`](dose-response-folds.md) | naive vs corrected dose response | pooled-ratio estimator, then a noisier inversion |
| [`module-transfer-inflation.md`](module-transfer-inflation.md) | how many modules transfer to an unseen stressor | a scoring bug in `module_transfer`; every count it produced was too high |
| [`pool-before-regulon.md`](pool-before-regulon.md) | pools respond before the regulons they share an agent with | no general ordering exists — H2O2's two pool arms go opposite ways |
| [`maintenance-coefficient.md`](maintenance-coefficient.md) | the stress-to-ATP coefficient | `300.0` was dimensionally incoherent and wrong by 10⁴; now a 6.5 mmol ATP/gDW/h envelope |
| [`reference-phenotype.md`](reference-phenotype.md) | the measured phenotype the GSMMs are validated against | the numbers were attributed to van Hoek 1998 and are not in it; two findings reverse |
| [`autofluorescence-first-pass.md`](autofluorescence-first-pass.md) | the reporter-free plate cannot give a value for the `a` sweep | `FINDINGS.md` measured the slope on 2026-08-30: the median `RFU/OD` it used is mostly background, and a = +18.2 RFU/OD closes the sweep |

## All eight are now reachable from the claim they overturned

This section used to record that `dose-response-folds.md` and `power-attribution.md` were
named by **no** document, test or script outside this directory, while the other five were
linked from between one and four places. Both are now linked: `DATA_INVENTORY.md` points here
for the n=2 fold table and the two estimator changes, and `FINDINGS.md` points here for the
two earlier well-level sigmas.

`autofluorescence-first-pass.md` arrived already linked: it was moved here from
`docs/research/` on 2026-09-11 with **zero** inbound references, and `FINDINGS.md` and
`DATA_INVENTORY.md` — the two live documents that overturned it — now name it.

Fixing the link found what an unreachable archive costs. Both live documents had kept a copy
of the superseded version inline rather than linking, and one of those copies had gone stale
without anyone noticing: `FINDINGS.md` carried *1% / 5% / 44% / 96%* in an unpinned table
directly below its own pinned table reading 7.5% / 16.7% / 52.5% / 96.7%, and this file
repeated the same four wrong numbers as the "current result". A superseded result kept beside
the live one is eventually read as the live one — which is the whole reason for the rule that
a live document states only the current result and links here.

## What this is not

Not a changelog: mechanical changes, refactors and new work do not belong here. Only
**a reported number that was replaced by a different number**.
