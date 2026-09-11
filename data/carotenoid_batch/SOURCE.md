# `data/carotenoid_batch/` — published *S. cerevisiae* carotenoid strains

One file, `published_batch_titres.tsv`: 29 published strains with their medium, cultivation
mode, titre and content, and — for the nine rows worked up in detail — their cassette
dosage, promoter and integration site.

This document exists because the file did not have one, and more importantly because the file
had **no code reader**. Five module docstrings made categorical claims about it
(`pathway/capacity.py`, `pathway/enzyme_capacity.py`, `pathway/environment_flux.py`,
`pathway/solve.py`, `predict.py`) and `grep -rn published_batch_titres` returned two prose
mentions and nothing else. Every one of those claims had gone false, and this table is what
falsified them. It now has a reader — `src/ystwin/pathway/published_cassettes.py` — which
validates every column against a declared vocabulary and RAISES on anything else.

## What each column means, and its domain

The reader enforces the domain. A cell outside it is a load-time `CassetteDomainError`
naming the row, not a silent pass.

| column | domain |
| --- | --- |
| `pmid`, `reference`, `journal`, `strain` | free text, provenance |
| `carbon_source` | free text, parsed against `published_cassettes.CARBON_SPECIES`; a cell naming no known species and not in `UNSTATED` raises |
| `mode`, `hours`, `bioreactor`, `measurand` | free text, cultivation description |
| `content_mg_per_gdcw`, `titre_mg_per_l`, `biomass_g_per_l` | a number or blank |
| `dfba_ready` | exactly `yes` or `no` |
| `copies_crtE`, `copies_crtYB`, `copies_crtI` | an integer, or `2u` (episomal, dose unstated), or `high-copy`, or `unrecorded`, or blank |
| `promoter` | a promoter name, or the sentinel `unrecorded`, or blank — **and see the population rule below** |
| `integration`, `copy_evidence` | free text; `copy_evidence` names WHERE in the source the dosage was read |

### The population rule, which was a real convention written down nowhere

**`promoter`, `integration` and the `copies_*` columns are populated exactly on the
`dfba_ready = yes` rows, and blank on every `dfba_ready = no` row.** Nine of twenty-nine are
`yes`, and the correspondence is exact in both directions. That is a statement about how much
work has been done on a row, not about what its source reports: a `no` row's paper may state
its promoter perfectly clearly, and the blank means nobody has transcribed it.

The reader enforces this in both directions, so filling a promoter on a `no` row fails until
someone changes the rule here first.

**A `no` row whose cassette IS now known keeps its blanks and puts the genotype in `notes`.**
Bu 2020 (PMID 33062054) is the worked example: its Table 1 prints one genomic copy each of
crtE, crtYB and crtI with their promoters, and that transcription is in the row's `notes`
rather than in `copies_*`, because the row is blocked on biomass and not on genotype. The rule
above is what makes `copies_*` mean "worked up end to end" instead of "someone happened to
read this one"; breaking it silently would cost more than the columns are worth.

`unrecorded` is the sentinel for the other case: a `yes` row whose source states the value and
whose cell has not been filled. It exists because the Bubphasawan 2025 row held the free-text
string `not recorded here` in a column otherwise holding promoter names — a non-value no
reader could tell from a promoter, and one contradicted by its own row's `copy_evidence`,
which says the genotype is "countable from the printed construct". **That to-do is now done:
on 2026-09-09 the full text (PMC12433382) resolved it from the printed genotype —
`1021b::TPRM9-SpCrtE(opt)-PGAL1-PGAL10-SpCrtYB(opt)-TGAL10-PGAL7-SpCrtI(opt)-TCPS1`, so crtE
is on PGAL1, crtYB on PGAL10 and crtI on PGAL7.** No cell holds `unrecorded` any more. The
sentinel stays declared, because the next partially-worked row needs it and rediscovering it
is how `not recorded here` got written the first time.

## `2u` and `high-copy` are not numbers, deliberately

An episomal 2-micron plasmid has a copy number that is a distribution over a population, not
an integer, and it varies with medium and selection. Converting either to a number would
manufacture the very quantitative dosage the survey is being asked whether it has. They are
distinct declared values so that `states_a_numeric_crtyb_dosage` can be false for them.

## What the table currently supports, and what it refutes

Computed by the reader rather than asserted here — run
`published_cassettes.crtyb_dosage_unfittable_reason()` for the live sentence.

* **Seven of twenty-nine rows state a crtYB dosage as a number**, and they still hold only
  **two** distinct values (1 and 2). Every pair of them differs in host, promoter, medium,
  cultivation mode and assay as well as in dosage, so the dosage → capacity coefficient is
  unfittable through **confounding**, not through absent numbers. `pathway/capacity.py`'s
  refusal rests on that, and `tests/test_capacity.py` fails the day a third comparable dosage
  arrives — which is exactly when the refusal should be revisited. The 2026-09-09 literature
  pass added two more numeric dosages (Fathi 2021 at 1, Bu 2022 at 2) and **moved neither the
  distinct count nor the refusal**, which is the sharpest evidence yet that the obstacle is
  confounding rather than sample size.
* **Not one row has ethanol in the feed.** The species seen are glucose, galactose, sucrose,
  xylose, fatty acid (olive oil, oleic acid), fruit juice, and undeclared rich media.
  `pathway/environment_flux.py` fits its descriptor on the feed's *ethanol carbon fraction*,
  so `z = 0` on every row here and no carotenoid row could fit that coefficient. That reason
  is robust to a new xylose row; the sentence it replaced — "every carotenoid row is glucose
  or galactose" — was not, and was already false.
  **`acetate` joined `ethanol` as a vocabulary entry matching nothing** on 2026-09-09: the
  only acetate row was Bu 2020, whose `carbon_source` had carried "+ 10 g/L acetate" from a
  *different strain's* experiment than its own 67.8 mg/L titre. Two real feeds the vocabulary
  still cannot name are Arhar 2024's 5 g/L lactic acid and 2 g/L acetic acid, so that row's
  `carbon_species` reads `{glucose}` alone; extending `CARBON_SPECIES` is a live to-do.
* **Six rows report an HPLC β-carotene content above 1.2483 mg/gDCW, of which four have no
  storage engineering**: Lange 2011 (3.90), Verwaal 2007 (5.90), Xie 2014 (7.41) and Arhar
  2024 (79.00). See `docs/research/KINETIC_FIT.md` §6 and `pathway/capacity.py` for what that
  does to the "membrane-holding limit" reading of the fitted ceiling — and for the caveat that
  travels with it: all four are batch or stationary harvests, where μ → 0 removes the growth
  dilution the steady-state law depends on.
  **The two rows added on 2026-09-09 do NOT extend that list, and the distinction is the whole
  point.** Fathi 2021 (46.50) engineers lipid accumulation with *Y. lipolytica* lipases and
  Bu 2022 (11.40) engineers lipid-droplet TAG metabolism directly — both are storage
  engineering, so both are consistent with a membrane-holding limit rather than against it.
  Bubphasawan 2025 (2.15) is HPLC-specific and storage-naive but only 1.7× the ceiling, well
  inside what a batch harvest can produce. The refutation still rests on the same four rows.

## `dfba_ready` and `bioreactor` are NOT conflated — the 7/22 match was a coincidence

Both columns split 7/22 before the 2026-09-09 pass, which looks like one column copied into
the other. **It is not.** The two sets are very nearly *disjoint*: they shared exactly one
row, Bubphasawan 2025 (Jaccard 0.067). A conflation would make them identical; a coincidence
of counts makes them complementary, which is what the file shows. Recomputing the two sets is
one call and settles it:

```
dfba_ready = yes  (9): Arhar 2024, Bu 2022, Bubphasawan 2025, Fathi 2021, Lange 2011,
                       Li 2013, Ukibe 2009, Verwaal 2007, Xie 2014
bioreactor = yes  (7): Bubphasawan 2025, Lin 2025, López 2019, Olson 2016, Su 2020,
                       Sun 2020, Xie 2015
```

The 2026-09-09 pass broke the tie by promoting two flask rows, so the counts are 9/7 and the
coincidence cannot be misread again.

**The report that "Verwaal 2007 is a shake flask yet carries `bioreactor = yes`" is false.**
Verwaal 2007 carries `bioreactor = no` in the working tree and in `git show HEAD:` alike, and
`no` is correct — the full text (PMC1932764) says only "Cells were grown at 225 rpm and 30 °C
in a shaking incubator". What Verwaal 2007 *does* carry is `dfba_ready = yes`, which is the
column the two are easy to swap, and that is almost certainly where the report came from.

The real reason the two columns look correlated is causal, not clerical: a row earns
`dfba_ready` by having a **β-carotene-specific assay and a content in mg/gDCW**, and the
bioreactor papers in this survey overwhelmingly quantify by 453 nm absorbance instead
(López 2019, Olson 2016, Sun 2020, Xie 2015 — all `total_carotenoid`). The columns are
anti-correlated *through the assay*, so treating one as a proxy for the other would be wrong
in both directions.

## What is UNRESOLVED, and why — the register

Recorded per row in `notes`; collected here so a reader can see the shape of what is missing.
**A refusal is a result.** Nothing below was upgraded on an abstract.

| blocked on | rows |
| --- | --- |
| **Paywalled, no accessible full text** | Yamano 1994, Yan 2012 ×2 (PDF read; assay simply not described), Li 2013, Reyes 2014, Wang 2014, Xie 2015, Yamada 2018, Zhao 2021, Yamada 2024, Lin 2025 |
| **Biomass — a specific β-carotene titre with no gDCW/L and no OD→DCW factor** | Bu 2020 (67.8 mg/L, HPLC 450 nm, OD600 only), Fan 2024 (166.79 mg/L, HPLC with authentic standards, the string "DCW" does not occur in the paper) |
| **Absolute number is in a figure or supplement only** | Rabeharindranto 2019 (HPLC-PDA with authentic standards; prose gives only "two times more β-carotene") |
| **Assay is a 453 nm absorbance sum, so no upgrade is possible** | López 2019, Cheng 2020, Sun 2020, Olson 2016 (OD454), Yamada 2022 (colorimetry) |
| **Culture volume and flask size not stated** | Verwaal 2007, Lange 2011, Fathi 2021, Rabeharindranto 2019, Arhar 2024 |
| **Temperature and rpm not stated** | Fan 2024 |
| **Cultivation duration not stated** | Cheng 2020 ("harvested at the end of fermentation") |

Two of these are worth naming because they cost a promotion each. **Bu 2020 and Fan 2024 both
have a clean, specific HPLC β-carotene assay and are blocked purely on biomass** — neither
paper reports a dry cell weight or a conversion factor, so their titres cannot become
contents. They are the two rows a single sentence in a supplement would unblock, and they are
the best targets if anyone revisits this survey.

**The Yamada lineage is a deliberate non-promotion.** Yamada 2022 (PMC8727735) quantifies by
"colorimetry [8]" against a β-carotene-in-hexane standard, uses the same YPH499/Mo3Crt79
strain as Yamada 2018, and its reference [8] *is* Yamada 2018 (PMID 30138874). That is strong
evidence that Yamada 2018 and Yamada 2024 are absorbance sums — but it is a sentence in the
2022 paper, not in theirs. Both stay `unverified` with the evidence in `notes`. Setting a
measurand from a sibling paper's citation is the same move as reading one out of an abstract,
which is what the Yamano 1994 audit note exists to stop.

## Provenance

Each row carries its own `pmid` and `reference`; `copy_evidence` names the specific figure,
table or sentence a dosage was read from. There is no single upstream export — the file was
assembled paper by paper — so the per-row citation IS the provenance, and a row without one
is a row that has not been checked.
