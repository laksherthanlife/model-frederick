# Superseded: the reporter-free plate's autofluorescence read

**Current result:** [`../FINDINGS.md`](../FINDINGS.md), "Autofluorescence measured at last, and
it is too small to matter" — **a = +18.2 RFU/OD**, 95% CI **[-6.1, +42.4]**, measured
2026-08-30 by `scripts/measure_autofluorescence.py` as one slope of RFU on OD per well and a t
interval over the thirty wells that grew. The sweep this document concluded had to stay is
closed: 18.2 RFU/OD is 0.42% of a reporter well with an upper bound of 0.97%, and the lowest
fraction at which any verdict in `outputs/autofluorescence_sensitivity.csv` changes is 5%.

**Why it moved.** The plate reading below is right and the method on top of it is the error.
Taking a median of `RFU/OD` measures `a + background/OD`, which at a four-hour plate's optical
densities is mostly background; a per-well slope removes the background term, which is why the
ratio returns hundreds of RFU/OD and the slope returns 18. So "Cannot: a value for the `a`
sweep" and "So the sweep stays" below are both overturned — by a reanalysis of this same plate,
not by a new one.

**One number here appears nowhere else.** The **239 RFU/OD** median in the table below is this
document's own: it is taken over the 28 wells it counts as grown, blank-corrected against the
60 cell-free wells of this plate. `docs/FINDINGS.md`, `docs/DATA_INVENTORY.md`,
`docs/research/XPT_INVENTORY.md` §2.2, `src/ystwin/readings.py` and
`scripts/measure_autofluorescence.py` all carry **249 RFU/OD** for the same median, taken over
thirty wells with each plate's own floor removed. Neither is a measurement of `a`; both are the
discarded ratio form.

**What survives.** Two findings below are not superseded and are current elsewhere: the
mCitrine channel *is* present in `20260807_BY4741_ER_Oxidative.xpt` and the `.xlsx` export
dropped it (`docs/DATA_INVENTORY.md` carries that correction), and `AFL` is the auto feedback
loop circuit, not autofluorescence (`../research/AUTO_FEEDBACK_LOOP.md`).

This file was `docs/research/AUTOFLUORESCENCE_FOUND.md` until 2026-09-11. Nothing below it
is edited but the one relative link the move would otherwise have broken; it stands as
written on 2026-08-29.

---

# The reporter-free plate WAS read in the fluorescence channel — and still does not measure autofluorescence

**2026-08-29.** `docs/DATA_INVENTORY.md` states, flatly, that autofluorescence "has never been
measured, on any plate", and that the two reporter-free BY4741 plates were "**never read in the
mCitrine channel**". README section 3 builds on that: one plate of the plain strain read in the
glow channel "settles three of our claims", and `outputs/autofluorescence_sensitivity.csv` exists
only because the analysis has to *sweep* assumed values instead of using a measured one.

**The first half of that is false.** The raw instrument file exists and it has the channel.

## What the file contains

`20260807_BY4741_ER_Oxidative.xpt`, read with `plate/gen5.py`:

| | |
| --- | --- |
| protocol | `4h-10min_mCitrine_OD600.prt` — the same file the sensor plates ran |
| channels | `OD600:600` **and `mCitrine:480,530`** |
| shape | 25 timepoints × 96 wells, both channels |

So the mCitrine read happened. The `.xlsx` export that reached this repository dropped it; the
instrument file did not. That sentence in `DATA_INVENTORY.md` needs correcting, and this file is
the correction.

## Why it still does not give a number

Read off the plate at the last timepoint, the layout is unambiguous:

    columns 1-6, rows A-F   cells      OD 0.09-0.70   RFU 168-217
    everything else         no cells   OD ~0.094      RFU ~56

That looks like a measurement of autofluorescence and it is not one, for a reason visible in the
numbers. Blank-corrected against the 60 cell-free wells, on the 28 wells that actually grew:

| quantity | value |
| --- | ---: |
| regression of excess RFU on excess OD | **slope −50 RFU/OD**, intercept **+150 RFU** |
| median per-well (excess RFU / excess OD) | 239 RFU/OD |
| max per-well ratio | 764 RFU/OD |

**The slope is negative.** More biomass, less excess fluorescence. Autofluorescence is by
definition a per-biomass term — `observation.py` writes the model as
`RFU = background + [autofluorescence + gain·R]·biomass` — so a term that does not rise with
biomass is not it. Nearly all the excess is the constant +150 RFU offset.

The likely cause is in the file's own name: this is `BY4741_ER_Oxidative`, so columns 1-6 carry
DTT and H2O2 ladders and columns 7-12 do not. The cell wells and the cell-free wells therefore
differ by **two** things — cells and stressor — and the design cannot separate them. DTT and
peroxide at millimolar concentrations are not optically neutral at 480/530.

## What can and cannot be concluded

**Can:** the measurement was taken, the channel is there, and `DATA_INVENTORY.md`'s "never read"
is wrong.

**Cannot:** a value for the `a` sweep. Taking the median per-well ratio at face value gives
12.7%–27.1% of the four constructs' own 0 mM per-cell signal (1876, 1859, 1078, 880 RFU/OD) —
which lands *inside* the 20% and 30% columns of
`outputs/autofluorescence_sensitivity.csv` rather than excluding them. An earlier reading of this
plate claimed those columns were excluded by measurement; that claim came from comparing against
a different set of blank wells and does not survive using the cell-free wells on the same plate.

**So the sweep stays.** What changes is that we now know the experiment exists, that it was run
with the stressor confounded into the comparison, and exactly what a clean version costs:

> One plate of BY4741 in **plain medium**, no stressor, read in `mCitrine:480,530` against
> cell-free wells of the **same** medium. Four hours, the existing protocol. That isolates the
> per-biomass term and settles the three claims README section 3 names.

The instrument file for it may already exist — `docs/research/DATA_SURVEY_2026_08.md` lists the
other candidates in the same directory.

**`AFL` does not stand for autofluorescence.** It is the team's **auto feedback loop** circuit —
a LacI autorepressor, which is why the glycerol-stock box lists `AFL`, `AFL control` and
`No-LacI` side by side. The four `AFL*.xpt` files are a separate experimental system and are
described in [`AUTO_FEEDBACK_LOOP.md`](../research/AUTO_FEEDBACK_LOOP.md); none of them is an
autofluorescence run, and the name is close enough to mislead, which it did.
