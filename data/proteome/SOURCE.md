# PaxDb, *S. cerevisiae* whole-organism integrated abundances

`paxdb_scerevisiae_integrated.tsv` — 6,351 proteins, abundance in **molar ppm** (parts per
million of protein MOLECULES, not of protein mass).

## Provenance

Fetched 2026-09-02 from
`https://pax-db.org/downloads/latest/datasets/4932/4932-WHOLE_ORGANISM-integrated.txt`.
The file's own header records what it is:

    #name: S.cerevisiae - Whole organism (Integrated)
    #score: 27.5
    #coverage: 96
    #description: integrated dataset: weighted average of all S.cerevisiae WHOLE_ORGANISM
    #organ: WHOLE_ORGANISM
    #integrated: true
    #publication_year: 2026

It is a weighted average over 26 constituent datasets, the heaviest being
`Saccharomyces_cerevisiae_SC_biomart_17916_O` (0.122), `Yeast_2009-0522` (0.110) and
`Saccharomyces-cerevisiae_PA_2013-3` (0.093). The upstream header carries the full weight
list; it was dropped from the vendored TSV to keep the file one table, and is reproduced in
the fetch script's output.

Upstream source bytes: sha256 `6d1b614fa95d4740…`, 153,753 bytes.

## The one check worth running on it

The abundances sum to **1,000,009 ppm**. That is the arithmetic that says the column really
is a molar fraction of the proteome and not an intensity in arbitrary units — a distinction
that decides whether it can be converted to mmol/gDCW at all. If a future refresh does not
sum to ~1e6, the column has changed meaning and nothing downstream should use it.

## What it is for, and the conversion it feeds

`src/ystwin/pathway/proteome.py` turns a ppm into an enzyme content:

    [E] mmol/gDCW  =  ppm/1e6  x  (total protein g/gDCW / average protein MW)  x  1000

which then feeds `pathway/enzyme_capacity.py`'s `vmax = kcat * [E]`. **The average protein
molar mass in that conversion is an assumption and it scales the answer linearly** — it is
declared in the module rather than hidden here.

## What it does NOT contain, and this is the limitation that matters

**No heterologous protein.** PaxDb is a measurement of the wild-type proteome, so CrtE,
CrtYB, CrtI, EEVS and every other cassette gene is absent by construction. For those the
abundance has to come from somewhere else: targeted proteomics on the producing strain (the
right answer), or an upper bound anchored on the promoter's own native protein — TDH3 sits
at 23,256 ppm, PGK1 at 13,767, TEF1 at 8,658. An anchor of that kind is a CEILING and must
be labelled one, because a heterologous protein folds and is degraded differently from the
native protein whose promoter it borrowed.

## How far the conversion was trusted, measured rather than assumed

Checked on 2026-09-02 against native enzymes where kcat, abundance and flux are each known
independently — kcat from the vendored `ecYeastGEM_batch`, abundance from this file, flux
from FBA at a fixed growth rate. Reported as capacity/flux, where 1.0 means the enzyme runs
exactly at capacity:

| enzyme | mu 0.10 | mu 0.25 | mu 0.40 |
| --- | ---: | ---: | ---: |
| ERG9 squalene synthase | 2.66 | **1.06** | 0.66 |
| ERG13 HMG-CoA synthase | 0.58 | 0.23 | 0.14 |
| ERG1 squalene epoxidase | 0.17 | 0.07 | 0.04 |

ERG9 lands within 6% at mu 0.25. ERG13 and ERG1 come out BELOW 1, which is infeasible — an
enzyme cannot carry more flux than its capacity — so for those the in vitro kcat understates
the in vivo rate. That is the known in vitro / in vivo apparent-kcat gap, here as a
measurement rather than a citation, and ERG1 being the worst is consistent with squalene
epoxidase being the classical rate-limiting step of sterol synthesis.

**So the honest error bar on a capacity derived this way is about an order of magnitude, and
it is one-sided: the derived capacity tends to be too LOW.** That is enough to catch a claim
that is wrong by 63x and not enough to set a headline number, which is exactly why
`pathway/enzyme_capacity.py` audits and refuses.
