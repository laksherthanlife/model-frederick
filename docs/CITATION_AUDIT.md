# Citation audit of `docs/`

Every literature identifier in `docs/` resolved against PubMed, PubMed Central, Crossref,
DataCite and arXiv on 2026-08-26, and each resolved title read against the claim beside it.
The guard holding this in place is `citations.py`, `scripts/refresh_citations.py` and
`tests/test_citations.py`, with four vendored tables under `data/citations/`.

Counts below are a **snapshot at 13:10 on 2026-08-26** and exclude this file. `docs/` was
being edited concurrently throughout — one document added, seven modified — so the totals
will already have moved; the additions were audited as they landed, see §7. Re-run
`python3 scripts/refresh_citations.py --check` for the current position.

---

## 1. Verdict

**357 distinct identifiers, across 608 citations, in 18 of the 34 documents under
`docs/`:** 181 PMIDs, 30 PMCIDs, 115 DOIs, 31 arXiv IDs. **All 357 resolve.** 356 resolve
to a bibliographic record; one — a PLoS supplementary-file DOI — has a live handle and no
record at any registry, and is now written down as exactly that. **282 of the 357 were
cited only from `docs/`** and sat entirely outside the old guard, which walked `src/`,
`scripts/` and `tests/` and understood only PMIDs.

**No identifier in `docs/` points at the wrong paper.** Every one of the 15 that resolves to
an unrelated paper is one `docs/SOURCE_AUDIT.md` deliberately records *as* an error, and
the vendored tables independently confirm each of its descriptions.

**The nine defects found are all of the opposite kind: prose that has gone stale around
correct identifiers.** Two are HIGH — two documents assert that fourteen citations are
still broken when they have all been fixed, and one of them carries a same-day
verification stamp. Four are MEDIUM: a stale cross-file warning, and two load-bearing
attributions carrying no identifier at all.

---

## 2. The wrong ones

Rows citing `docs/CLEANUP.md` point at a document that no longer exists: it was removed
in `7ec4106` once its open items were either resolved or moved into this audit. The line
numbers are kept because they locate the claim in that commit's version of the file, which
is still reachable with `git show 7ec4106^:docs/CLEANUP.md`.

| # | File : line | What the doc claims | What it actually is | Severity |
|---|---|---|---|---|
| 1 | `docs/research/PARAMETER_SOURCES.md`:239–247 | "**Still outstanding from `SOURCE_AUDIT.md` §3.1** (verified again on **2026-08-26** against the current file, and *not* yet fixed): the 14 broken module PMIDs are still present — `stress_panel.py:210` (10361303 → Tovar 1999), `:213` (9663389 → Galant 1998), `:250` (32160542 …), `:253` (24815987 …), `:259` (21982710 …). The QUEEN-2m *reporter* entry was repaired; the `MODULES` entries were not." | **False, and stamped with today's date.** All fourteen identifiers return **zero** hits in `src/ystwin/generator/stress_panel.py`. Every `MODULES` entry now carries the exact replacement §3.1 prescribed: cell_wall 10594829, dna_damage 9741624, calcium 9407035, carbon 12676948, nitrogen 10604478, hypoxia 9510529, copper 3043194, zinc 9271382, sulfur 9409150, retrograde 8422683, peroxide 32130885, atp 25283467, nadh 21982714, and alkaline_ph now 11698381 + 11050096. The five line numbers point at unrelated code: `:210` is blank, `:213` is a docstring, and `:250`/`:253`/`:259` are the tails of the *corrected* source strings. The trailing "one further slip" about the QUEEN-2m source citing "Takaine 2021" is stale too — see finding 3 | **HIGH** |
| 2 | `docs/CLEANUP.md`:723 | "~~`docs/research/PARAMETER_SOURCES.md:240`~~ — ~~'…and **not** yet fixed: the 14 broken module PMIDs are still present'~~ — **RESOLVED** — `abcad97`. **The assertion is gone from the file**" | **The assertion is not gone from the file.** It is at `docs/research/PARAMETER_SOURCES.md`:240, verbatim, at exactly the line this row names. The two findings compound: a reader who checks the resolution log believes the problem closed, and a reader who reads `PARAMETER_SOURCES.md` believes fourteen citations broken. Neither is right, and each document's error is what makes the other's survive | **HIGH** |
| 3 | `docs/CLEANUP.md`:385–389 | "### N1d (low) — **OPEN** — the vendored table and a `source` string disagree on a year … `stress_panel.py:433` cites *'Takaine 2021, PMID 33654827'*" | **Fixed, and the same document says so 307 lines earlier.** `docs/CLEANUP.md`:78 records N1d as **FIXED**. The code at `src/ystwin/generator/stress_panel.py`:558–563 now reads "the Bio-protocol imaging method in PMID 33654827 — also 2019 in PubMed, though its PMID looks 2021 because Bio-protocol indexed it late". The string "Takaine 2021" appears nowhere in the file, and line 433 is a NADH reporter's `target_ec50` | MEDIUM |
| 4 | `docs/research/UPR_ANCHOR.md`:735 and :808 | "⚠ `G4_STATISTICS.md` cites PMC3326325 — that identifier belongs to an unrelated ribosome paper", repeated in the closing list as "ValidPrime is PMC3326333, not PMC3326325" | **The defect is fixed and the warning is stale.** `docs/research/G4_STATISTICS.md`:110 cites PMC3326333 and documents the old error in place. UPR_ANCHOR sends a reader to check a file that is already correct — and its own §7 entry is the authority the fix was made from | MEDIUM |
| 5 | `docs/research/PREDICTION.md`:418–423 | "Bracher et al. describe the components as *sharpness* and *calibration*, and the word 'dispersion' appears exactly once in their paper — in the PIT-histogram section, not as a component name." Used to instruct the reader **not** to attribute the three-way WIS decomposition to Bracher | **A precise, falsifiable, load-bearing claim about a paper that is never identified.** No DOI, PMID or arXiv ID anywhere in the document. It is Bracher, Ray, Gneiting & Reich, "Evaluating epidemic forecasts in an interval format", *PLoS Comput Biol* 17(2):e1008618, 2021 — doi:10.1371/journal.pcbi.1008618, PMID 33577550. A guard built on identifiers cannot see a citation that has none | MEDIUM |
| 6 | `docs/research/PREDICTION.md`:145–150 | The capability-ladder correction: the 0–5 ladder is *not* in Rasheed, San & Kvamsdal's *IEEE Access* 2020 paper, which "lists eight 'values' of a digital twin and no capability levels at all"; it is in San, Rasheed & Kvamsdal, *GAMM-Mitteilungen* 44(4):e202100007, 2021 | **The correction is itself uncitable.** Neither the paper being ruled out nor the one being credited carries an identifier, and the claim is a *negative* about a named paper's contents — the hardest kind for a reader to check and the one most needing a resolvable handle. This is the attribution the document itself says is "commonly got wrong" | MEDIUM |
| 7 | `docs/research/PREDICTION.md`, §8 | Annotated bibliography, §8.1 → §8.5, then §8.7 | **§8.6 does not exist.** Finding 5 seen from the other side: the forecast-hub and WIS references, Bracher included, have no section to live in, which is why the most checkable claim in §2.6 has no entry | LOW |
| 8 | `docs/SOURCE_AUDIT.md`:97 and :393 | PMC3326325 is "Sharifulin **2011**, *NAR* 40:3056" | The record is *Nucleic Acids Res* volume **40, 2012** (e-pub 2011 Dec 16). Meanwhile `docs/research/G4_STATISTICS.md`:110 calls it an unrelated paper "in the same journal and **year**" as the 2012 ValidPrime paper. Both cannot hold; the volume settles it as 2012 | LOW |
| 9 | `docs/research/PREDICTION.md`:1263 | "DOI **prefix** `10.1175/MWR-D-14-00269.1` is stated on the AMS landing page" | Not a prefix. `10.1175` is the prefix; that string is the complete DOI and it resolves cleanly to Scheuerer & Hamill, "Variogram-Based Proper Scoring Rules for Probabilistic Forecasts of Multivariate Quantities", *Monthly Weather Review* 143, 2015. The hedge understates what the document already has | LOW |

Findings 1–4 share a shape worth naming: **the citations were fixed and the prose about
them was not.** That is the failure mode a title-resolving guard creates if nothing else
changes — it makes identifiers trustworthy and leaves the commentary about identifiers
unchecked. Findings 1 and 2 are the expensive version, because they tell a reader that the
repository's most-imported module is carrying fourteen fabricated-looking citations when it
is carrying none.

### One live wrong identifier, and it is not in `docs/`

`src/ystwin/generator/context.py`:35 still cites **PMID 8710508** for "Schuller & Entian
1996 — ICL1 moves >200-fold only on transfer to a non-fermentative source". That PMID is
*"The complete nucleotide sequence of bacteriophage HP1 DNA"*. Not a new finding —
`docs/SOURCE_AUDIT.md` §3.1 records it at MEDIUM and notes that no Schüller/Entian 1996
paper exists to substitute — and it was already in the old vendored table. It is the only
one of the fifteen documented-wrong identifiers still doing work in code rather than being
quoted as evidence, and the only one the prose is right about.

### Two flags that do not survive checking

**The reviewer's Takaine flag is misattributed.** It was reported that
`docs/research/UPR_ANCHOR.md` gives PMID 33654827 as Takaine **2021** when it is **2019**.
`UPR_ANCHOR.md` never mentions Takaine and never cites 33654827. The "2021" was in
`src/ystwin/generator/stress_panel.py`, inside code, already covered by the old guard —
and it has since been fixed there too, so the error now exists in no file at all. What
remains is the stale prose about it, finding 3. The record: PMID 33654827, *Bio Protoc*
9(15):e3320, pubdate 2019 Aug 5. Note the repository cites two Takaine papers a year
apart, correctly: 30858198 (*J Cell Sci* 132, 2019) is the yeast demonstration, 33654827
is the protocol.

**And one I nearly filed myself.** I briefly "found" a spliced reference entry in
`UPR_ANCHOR.md` — Rüegsegger's authors and title over Teste's journal, DOI and PMID — which
turned out to be two correct entries read across a seam in my own `sed` output. Both are
right: Rüegsegger, Leber & Walter, *Cell* 107(1), 2001, PMID 11595189 at line 609; Teste
et al., *BMC Mol Biol* 10:99, 2009, PMID 19874630 at line 711. Recorded because the
vendored table is what disproved it in seconds, which is the only argument for the table
that matters.

---

## 3. Identifiers that resolve to unrelated papers **on purpose**

`docs/SOURCE_AUDIT.md` is an audit document. It tabulates the fourteen broken PMIDs from
`src/ystwin/generator/stress_panel.py` against what each really points at, so it *contains*
the wrong identifiers as evidence. **These are not errors in `docs/` and must not be
reported as such.** All fifteen resolve exactly as the audit describes — the vendored table
is independent confirmation of the audit's own claims, including both "off by one in the
same issue" cases:

| PMID cited (wrong) | Resolves to | Correct PMID |
|---|---|---|
| 10361303 | Tovar 1999 — the mitosome in *Entamoeba histolytica* | 10594829 |
| 9663389 | Galant 1998 — butterfly achaete-scute homolog, wing scales | 9741624 |
| 9420334 | Sellers 1998 — retinoblastoma protein / E2F | 9407035 |
| 12482873 | Despouy 2003 — cyclin D3 / retinoic acid receptors | 12676948 |
| 10604474 | Liang 1999 — beclin 1 and autophagy | 10604478 |
| 9515917 | Maupin-Furlow 1998 — 20S proteasome of *Methanosarcina* | 9510529 |
| 3062374 | Robinson 1988 — vacuolar protein sorting mutants | 3043194 |
| 9223284 | Minvielle-Sebastia 1997 — poly(A)-binding protein | 9271382 |
| 9409149 | Brown 1997 — "Archaea and the prokaryote-to-eukaryote transition" | 9409150 — which resolves to Thomas & Surdin-Kerjan, "Metabolism of sulfur amino acids in *Saccharomyces cerevisiae*", so the audit's "off by one" is exact |
| 8500171 | Fütterer 1993 — ribosome migration on cauliflower mosaic virus RNA | 8422683 |
| 12588995 | Kren 2003 — War1p and weak-acid stress | 12509465 |
| 32160542 | Laver 2020 — RNA-binding protein Rasputin/G3BP | 32130885 |
| 24815987 | Vasefi 2014 — polarization-sensitive hyperspectral **dermoscopy** | 25283467 |
| 21982710 | Krahmer 2011 — phosphatidylcholine synthesis / lipid droplets | 21982714 |
| PMC3326325 (a PMCID, tabulated with these) | Sharifulin — ribosomal protein S26 | PMC3326333 |

Fourteen of the fifteen are now cited only from the audit documents; 8710508 is the
exception, above.

---

## 4. Deliberately unverified, and kept

Seventy markers — 20 `UNVERIFIED`, 50 `⚠` — across eleven documents. **Every one is honest
and every one stays.** They are the difference between a bibliography and a bibliography
you can trust, and removing a marker to make a page look finished would be the worst edit
available here. A marker records that the *identifier* is right and the *claim* was not
confirmed from primary text — precisely the distinction the guard cannot make.

| File | Markers | What they cover |
|---|---|---|
| `docs/research/UPR_ANCHOR.md` | 23 | Its own §7 preamble says so: "⚠ marks something the reader should check." Paywalled or unretrievable primary text (Travers 2000; Cell blocks automation), a scanned-image-only PMC record (PMID 8313910 / PMC394856), an abstract-only PMC record (PMID 11102521 / PMC15070), graphical-only magnitudes, and two bibliographic traps it corrects rather than propagates: the re-issued Elsevier suffix on PMID 9382810, and the `rrp6Δ`/`cbc1Δ` ratio misread as a wild-type fold change behind the PDI1 "2–2.5-fold" figure |
| `docs/SOURCE_AUDIT.md` | 17 | Claims it could not confirm behind a paywall or a failed fetch: the Alberty magnesium constants at I = 0.25 M, the Hahn 2006 MMS 0.02 % dose, the Görner 1998 *hog1Δ* detail (correct paper, unread detail), PMC7646510 and PMC3847955 not re-read, HyPer7 in *S. cerevisiae* found only for *Komagataella phaffii* |
| `docs/ARCHITECTURE.md` | 16 | Design cautions rather than citations |
| `docs/research/CIRCULARITY.md` | 3 | Preprint-only status called out per entry — "arXiv only, never published", "preprint, no venue" — plus two genuine traps: the corrupt DOI field on the arXiv record for Frazier, Robert & Rousseau, and the published TMLR title of Hermans et al. differing from the arXiv "A Trust Crisis…" (arXiv:2110.06581) |
| `docs/research/IDENTIFIABILITY.md` | 3 | "metadata verified, full text unread" on Shapiro 1985, Bekker & ten Berge 1997, Eriksson & Koivunen 2004; and Anderson & Rubin (1956) listed under **Unverified** as not indexed in Crossref and not retrieved |
| `docs/research/PARAMETER_SOURCES.md` | 3 | The Hahn 2006 dose, the YRE-versus-HSE contribution, the HyPer7 yeast demonstration |
| `docs/CLEANUP.md`, `docs/research/EXTERNAL_DATA.md`, `docs/research/G4_STATISTICS.md`, `docs/research/PREDICTION.md`, `docs/research/TOOLING.md` | 1 each | Single scoped caveats, including PREDICTION's own statement of the convention |

Two of these are models rather than caveats. `docs/research/CIRCULARITY.md` flags an
arXiv/published **title** divergence in prose; the arXiv table now carries the preprint
title in a column, so that divergence is visible in a diff as well as in a sentence. And
`docs/research/G4_STATISTICS.md`:110 keeps the wrong PMCID in the text while explaining it
— which is why both PMC3326325 and PMC3326333 are legitimately cited, and why
`tests/test_citations.py` pins which is which so the explanation cannot quietly invert.

---

## 5. What the guard now covers

`citations.py` is one scanner, imported by both `scripts/refresh_citations.py` and
`tests/test_citations.py`. Each file previously carried its own copy of the pattern and of
a `_logical_text` helper; two scanners can drift so that the refresh vendors an identifier
the test never asks about, or the reverse. There is now one definition.

**Surfaces:** `src/`, `scripts/`, `tests/`, `docs/` — `.py` and `.md`.

| Table | Rows | Why it is its own file |
|---|---|---|
| `data/citations/pmid_titles.csv` | 182 (was 53) | PMIDs from code *and* prose stay in one table: same resolver, same schema, and a PMID cited from both surfaces is one paper that should appear once with every citer listed. Splitting by citing surface would duplicate rows and break the no-duplicates invariant |
| `data/citations/pmcid_titles.csv` | 30 | A different NCBI database, and it carries the **PMC → PMID mapping**, which is the actual check. The one PMCID error here was an identifier one digit off; what exposes that fastest is that its PMID is not the PMID named beside it |
| `data/citations/doi_titles.csv` | 117 | Crossref, then DataCite, then the DOI handle system, with a `registry` column recording which answered. Three tiers because a DOI's registration agency is not knowable from the string: 8 dataset DOIs here (Zenodo, Edinburgh DataShare) are DataCite-only and absent from Crossref, and calling them unresolvable would have been wrong |
| `data/citations/arxiv_titles.csv` | 31 | Two extra columns that would be dead weight elsewhere and are the whole point here: `version`, because a result can move between preprint versions so an unversioned identifier pins nothing — 30 of 31 citations are unpinned, and the one that is pinned, `arXiv:2012.05841v2` in `docs/research/PREDICTION.md`, is correct practice; and `published_title`, because an arXiv title differing from the published title lets a citation name one document and point at another while reading as consistent |

### The markdown analogue of the string-concatenation hole

The Python surface hides a citation through implicit string concatenation, and rejoining
adjacent literals fixes it. Markdown's version is the **soft line break**: a reference entry
wraps so that `PMID` ends one line and the digits begin the next.
`docs/research/UPR_ANCHOR.md` does this twice — PMID 34047633 at line 653, PMID 12184808 at
line 722 — and the old pattern was `PMID[: ]+(\d{6,8})`, in which `[: ]` does not match a
newline. **Both citations were invisible to the guard while perfectly legible to a reader.**

The tempting fix is to unwrap paragraphs before matching, by analogy with the Python fix.
**That is the wrong tool, and I tested it: it finds zero additional identifiers.** Markdown
renders a soft break *as whitespace*, so a citation split mid-token is broken in the
rendered document too, not merely hidden, and unwrapping cannot repair one. The correct fix
is that the separator between a label and its identifier is whitespace **including
newlines**, which is what markdown says it is. No unwrap step was shipped.

Three further markdown-only forms, all real here and all previously invisible:

- **Emphasis between label and digits** — `PMID **12676948**`
  (`docs/research/PARAMETER_SOURCES.md`:73).
- **A plural label introducing a list** — `PMIDs 11102521, 11179418, 9712873`
  (`PARAMETER_SOURCES.md`:89, 93, 94, 96, 153). `PMIDs` does not match `PMID[: ]` at all,
  and the second and third identifiers carry no label of their own.
- **A table column whose header names PMID**, with no inline label anywhere on the row.
  `docs/SOURCE_AUDIT.md` §3.1, §4.1 and §4.3 are built this way, and those rows are the
  repository's record of which identifiers were checked — the most load-bearing citations
  in it. **31 identifiers are reachable only this way**, including the entire 34-PMID
  summary cell at `SOURCE_AUDIT.md`:387. Header detection is structural: a row is a header
  if and only if the `---` rule is the next row, so a data row mentioning PMID cannot
  promote a column of page numbers into a column of citations. A `PMC`-prefixed number in
  such a column is masked before the digits are read, because a PMCID is never a PMID —
  without that, an audit row tabulating a *PMCID* error under a PMID heading mints a PMID
  that does not exist. That bug existed in this scanner for an afternoon and is pinned now.

Each form has a regression pin in `tests/test_citations.py`, so narrowing the separator back
to a literal space and colon fails the suite rather than silently shrinking coverage. So
does `test_documentation_is_actually_covered`, which fails if the walk stops reaching
`docs/` — the way the previous hole survived was that nothing failed when nothing was
looked at.

---

## 6. What the guard still cannot do

**It cannot judge whether a paper supports a sentence,** and that is not a temporary
limitation. Every finding in §2 was found by reading. Findings 5 and 6 are invisible to
*any* identifier-based check because the defect is the absence of an identifier, and
findings 1–4 are invisible because the identifiers involved are all correct — what is wrong
is the prose describing them. **A guard that resolves identifiers makes identifiers
trustworthy and leaves commentary about identifiers entirely unchecked**, which is where
this audit's two HIGH findings live.

Specifically:

- **An uncited claim is uncovered.** Findings 5 and 6 are precise assertions about named
  papers, and nothing automated can want an identifier that was never written.
- **A right identifier on a wrong claim passes.** The tables pin identifier → title.
  Whether *"Ledermann's bound gives r ≤ …"* is what Ledermann 1937 proves is a judgement,
  and 356 of these titles were read against their claims by a person, not a test.
- **Prose about the state of the code is unchecked.** `scripts/audit_claims.py` verifies
  that a document does not name a module that does not exist; nothing verifies that a
  document's claim about what is *inside* a module is still true. Findings 1, 2 and 3 are
  all of that shape, and a check is buildable: a document asserting "PMID X is still cited
  at `file:line`" is mechanically falsifiable. It is not built here.
- **Verbatim quotes are only spot-checked.** `docs/SOURCE_AUDIT.md` marks several claims
  "verified verbatim". Two spot-checks confirmed the quoted string exactly in the abstract
  (PMID 10411744, "The two sets of proteins overlapped only slightly"; PMID 18627600, "a
  positive transcriptional loop"). Two others quote figures not in the abstract (−306 mV
  from PMID 22705944, pH 7.5 from PMID 23762325) — which is *not* evidence against them,
  since both would live in the full text, and reporting them as findings would be exactly
  the sloppy move this audit exists to prevent.
- **Author and year are not pinned.** The tables carry `year`, and nothing asserts that a
  document's "Takaine 2021" matches it. Adding a first-author column would strengthen the
  diff; I did not, to keep the schema the tests already depend on stable.
- **A PMID list interleaved with author names is partly invisible.** The label scan stops at
  the first letter, so in `**34 PMIDs** — 8898194 Sidrauski 1996; 9323131 Sidrauski 1997; …`
  only the table-column path reaches the rest. It does reach it here, so that cell is fully
  covered — but a list of that shape *outside* a PMID-headed column would not be.
- **URL-only citations are uncovered by design.** Ojala & Garriga 2010 (JMLR) and
  DNVGL-RP-A204 have no DOI to resolve. A property of the sources, not a gap to paper over.
- **The tables are a snapshot.** Titles were resolved on 2026-08-26. A publisher retitling
  a record shows up only on the next refresh, which is the intent: the refresh is
  deliberate and the diff is the review.

---

## 7. The corpus moved during the audit

`docs/superseded/pool-before-regulon.md` was added, and seven documents modified, by
concurrent work while this was running. That document's citations grew from ten to thirteen
mid-audit and **all thirteen are correct** — author, year and title match the record in
every case: PMIDs 10468205 (Winterbourn & Metodiewa 1999), 17210445 (Ogusucu 2007),
31389622 (Calabrese 2019), 12437921 (Delaunay 2002), 11013218 (Delaunay 2000), 28418333
(Goulev 2017), 23762325 (Ayer 2013), 22561702 (Dardalhon 2012), 14556853 (Azevedo 2003),
26944189 (Tomalin 2016), 34118234 (Kritsiligkou 2021), 23242256 (Morgan 2013) and 28540740
(Deponte 2017). Deponte's title is quoted verbatim in the prose and matches the record
exactly. The two identifiers added to `docs/CLEANUP.md` in the same window, PMIDs 11698381
(Xu & Mitchell 2001) and 11050096 (Lamb 2001), are correct as well and are the
`alkaline_ph` repair.

That a live document accreted thirteen citations during a single afternoon's audit is the
argument for the guard being a test rather than a document: this file is already out of
date, and `tests/test_citations.py` is not.

---

## 8. Reproducing this

```bash
python3 scripts/refresh_citations.py                              # all four tables
python3 scripts/refresh_citations.py --check --kind pmid          # drift only, no write
python3 -m pytest tests/test_citations.py -q -p no:cacheprovider  # 35 passed
python3 scripts/audit_claims.py                                   # all checks pass
```

`--check` reports drift and exits non-zero without writing, so it can gate a release.
`--kind` narrows the refresh, which matters because the DOI pass makes 117 Crossref
requests at a polite rate.
