# The ATP sensor plates contain a response the t = 10 h comparison hid

Three constructs, all reported as not working: an ICL-UAS build tested 2026-08-05, and an
ACS-UAS and second ICL-UAS build tested 2026-08-09. Each was read across four
galactose:glucose ratios at four total sugar concentrations (2%, 1.5%, 1%, 0.5%), 18 h,
OD600 and mCitrine at Ex480/Em530.

Reproduce with `python scripts/analyse_atp_sensor.py`. Tests: `tests/test_atp_sensor.py`.

## There is a large, ordered, well-timed response

Each well's RFU/OD600 is flat through exponential growth and then rises to a peak:

| total sugar | glu 10:0 | gal 3:7 | gal 6:4 | gal 9:1 |
| --- | --- | --- | --- | --- |
| 2% | **2.84× @ h12.6** | 1.93× @ h10.8 | 1.38× (no peak) | 1.26× (no peak) |
| 1.5% | 2.12× @ h11.5 | 1.46× @ h9.1 | 1.26× (no peak) | 1.29× (no peak) |
| 1% | 1.48× @ h9.0 | 1.15× (no peak) | — | — |
| 0.5% | 1.14× (no peak) | — | — | — |

Two orderings, both clean. Induction needs glucose and scales with it. And **the peak
arrives later the more glucose there was** — h9.0, h11.5, h12.6 for 1%, 1.5%, 2%. A dose
response in time is what running out of a sugar looks like, and it is the strongest evidence
in the dataset.

Inverting dilution properly (`k = d(F/OD)/dt + mu·(F/OD)`) puts synthesis in the rise window
at 1.64× and 1.37× its exponential-phase value for 2% and 1.5% glucose.

## Why comparing at t = 10 h hid it

At h10 the conditions are at different points of their own transients:

| | RFU/OD at h10 | its own peak | |
| --- | --- | --- | --- |
| 2% glucose | 19381 | 24026 @ h12.6 | still rising |
| 1.5% glucose | 15736 | 17324 @ h11.5 | still rising |
| 1% glucose | 11913 | 12379 @ h9.0 | already declining |
| 0.5% glucose | 7412 | 9884 @ h0.1 | never peaked |

A bar chart at one timepoint therefore reports how far through its own transient each
condition happens to be, mixed with how much signal it makes. The same effect drives the
2.9–4.7× spread at t = 12 h across the whole plate, which collapses to 1.08–1.21× when
conditions are matched on optical density instead — they sit up to 11 h apart at the same
density.

## The expectation was also wrong for this promoter

The slides expect normalized fluorescence to rise monotonically as growth conditions worsen,
as an ATP sensor would. A CSRE element is not one. It answers to Cat8 and Adr1, which answer
to running out of *fermentable* sugar — ACS1 derepresses several hundredfold on ethanol or
under sugar limitation (Kratzer & Schuller 1997, PMID 9427394); ICL1 is a Cat8-driven
glyoxylate-cycle gene. The correct prediction is a transient at the diauxic shift, timed by
how much sugar there was, which is what the plates show. Galactose is another fermentable
sugar, which is why the galactose-rich columns never peak.

The panel already routes these constructs to the **carbon** regulon (Snf1 → Adr1/Cat8,
marker genes ADH2/ICL1) rather than to ATP. QUEEN is the only ATP reporter in the library
because no promoter fusion reads the adenylate pool.

## What still cannot be concluded

The response is real and ordered. Whether it is *the promoter* is not settled here:

1. **No promoterless strain on either plate.** At the diauxic shift cells build
   mitochondria, and flavin autofluorescence at 480/530 rises with them — scaling with
   glucose consumed, exactly as the observed signal does. This alternative fits every number
   above and only a promoterless control excludes it.
2. **Raw RFU falls 9–28% after its peak**, and the loss scales with peak height. Nothing
   stable does that: mCitrine is not degraded and cells do not lose it. Some of the signal
   is being quenched or destroyed in a way none of the kinetics here models.
3. **No constitutive-promoter strain**, so growth-rate coupling cannot be subtracted.

## What to run next

- **A promoterless strain and a constitutive strain on the same plate.** This is the single
  highest-value addition and it settles point 1.
- **Ethanol or glycerol against glucose** — the axis these promoters actually read, where
  the expected effect is hundredfold rather than threefold.
- **Compare on density or on each well's own peak, never on one clock time.**
  `ystwin.analysis.matched_density` does the first.

## What this changed in the model

Nothing. `_CARBON["galactose"]["carbon"]` is **0.05** and has been since the baseline commit
`f562d8f`; `git log -S` finds no change to that line.

This paragraph previously said the value "was briefly lowered to 0.10 ... and has been put
back to 0.45". None of that happened: it was never 0.10, never 0.45, and never edited. The
conclusion it reached was right for the wrong reason and still holds — the galactose columns
are flat, but until a promoterless control says the reporter reads anything, a flat response
calibrates nothing, so nothing here has moved a constant. Corrected 2026-08-30.
