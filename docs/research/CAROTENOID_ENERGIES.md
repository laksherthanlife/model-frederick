# Why the thermodynamic gate said `cannot_say` about beta-carotene, and why it no longer does

> **RESOLVED 2026-08-30, and neither reason was a fact about carotenoid chemistry.** The first
> was a naming gap: the C40 species were in the tables under ModelSEED ids nobody had
> declared. The second, found the same day, was that `fba/carotenoid.py` wrote CrtI with four
> free FAD as terminal acceptor — **+166.2 ± 12.8 kJ/mol, thermodynamically impossible**, on a
> step Verwaal 2007 (PMID 17496128) measured running to completion. That impossible number is
> what made the two estimators straddle zero, which is what this file recorded as irreducible
> disagreement. Written as the flavin oxidase it is (`+ 4 O₂ → + 4 H₂O₂`, −292.7 ± 17.4), the
> sign disagreement is gone and all three steps gate as **`runs`**.
>
> What survives is the magnitude disagreement below, 22 to 124 kJ/mol — and it is **not
> specific to C40 polyenes**, which is what this file used to argue: `CRTE`, which touches no
> C40 species at all, disagrees by the largest ratio of the four at 8.6×. So the refutation now
> refuses a NUMBER and not a VERDICT. The analysis below is kept as written because its
> conclusion — do not wire one of these energies in as truth — is still right.

**2026-08-29, revised 2026-08-30.** `pathway/thermo_gate.py` is wired into `predict.py` and
works: on `phb` it returns `cannot_run` at 10 µM cytosolic acetyl-CoA — and, after the
kcal/kJ correction described below, `cannot_run` at 425 µM as well. That correction refuted
the thiolase lead; `scripts/thiolase_threshold.py` is the record. On `beta_carotene` all three
declared steps return `cannot_say`. This file records why, and why that is the right answer to
leave in place.

**The unit error, because it is the same defect this file was opened about.** An earlier
analysis here was refuted three votes to zero for reading a kcal/mol table as kJ/mol. The
refutation was right about the table and did not go far enough: `bridge/thermodynamic.py`
itself had the same error, reading the vendored ModelSEED database through `pytfa` under its
default `thermo_unit='kJ/mol'`, so the pH and ionic-strength transform was computed in kJ on
top of a tabulated energy in kcal. Water came out at +24.85 kJ/mol. Corrected it is −155.70
against Alberty's −155.66, and over the 1,098 single-compartment yeast-GEM reactions shared
with eQuilibrator the agreement moves from RMSE 320.1 / r 0.435 to RMSE 116.6 / r 0.942.

The estimator disagreement below shrank with it — 22, 225 and 124 kJ/mol where this file
first recorded 91, 275 and 252 — because part of the old gap was the unit and not the
chemistry. It did not shrink enough to change the conclusion, and `CRTI` still disagrees in
sign: −59 kJ/mol against +166.

## The vendored tables do carry the energies

`data/thermo/compounds.tsv` has all three C40 species:

| ModelSEED id | species | formula | `deltag` | `deltagerr` |
| --- | --- | --- | ---: | ---: |
| `cpd03205` | phytoene | C40H64 | 844.55 | 4.78 |
| `cpd03217` | lycopene | C40H56 | 840.60 | 5.45 |
| `cpd01420` | beta-carotene | C40H56 | 771.52 | 5.35 |

So the earlier statement in `data/pathways/beta_carotene.toml` — that the tables "carry no
formation energy for a C40 carotenoid" — is **wrong as written**. They carry values.

## But the loader cannot reach them, and the reason is structural

`ThermodynamicData.load().dgf()` returns `None` for all three ids. The gap is not coverage,
it is naming: `fba/carotenoid.py` invents the metabolite ids `phytoene_c`, `lycopene_c` and
`betacarotene_c` when it installs the pathway, and `data/thermo/aliases.tsv` maps ModelSEED
ids onto ids from *published* models (AlgaGEM, BiGG, AraCyc and so on). An id this repository
made up an hour ago is in no alias table anywhere, and never will be.

Closing that is a small, real piece of work: have the spec or `carotenoid.py` declare each new
metabolite's ModelSEED id, and have `thermodynamic.py` resolve through it. **It is deliberately
not done here**, because of the next section.

## Supplying an energy would make the gate refute a strain that works

Two estimators are available and they do not agree.

**The units are a trap, and the first attempt fell into it.** `data/thermo/thermo_data.thermodb`
is in **kcal/mol** — its own units field says so, and `cpd00001` (water) reads −56.687, which is
−237.2 kJ/mol, the standard value. `pytfa`'s `MetaboliteThermo` defaults to `thermo_unit='kJ/mol'`,
so constructing it without passing the unit inflates every energy by 4.184×. An adversarial pass
caught this three votes to zero on a first version of this analysis; the numbers below are the
corrected ones.

Corrected, the ModelSEED and eQuilibrator estimates for the three pathway steps still disagree by
**91, 101 and 252 kJ/mol** while quoting uncertainties of **5, 13 and 15** — 8 to 18 sigma apart.
Both are extrapolations: a C40 polyene with eleven conjugated double bonds is outside the group
decomposition either method was trained on.

And the consequence is not neutral. Supplying either estimate turns the CrtI desaturase step into
`cannot_run`, so `require_feasible` would refuse **a strain Elizondo 2025 measured making
beta-carotene**. A gate that refutes a strain known to work is worse than a gate that abstains.

## The recommendation, and what is in the tree

**Do not wire an energy in.** `cannot_say` is an accurate report of the state of knowledge;
a number 15 sigma wide is not. Nothing was wired in, and the module that carried the first,
unit-erroneous attempt was deleted rather than corrected — a refuted analysis is not a
capability.

What would change this is a *measured* formation energy for a carotenoid, not a better
extrapolation from the same groups.
