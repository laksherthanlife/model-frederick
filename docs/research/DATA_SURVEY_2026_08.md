# What wet-lab data exists on this machine, and what the repository knows about it

Survey of `~/Desktop/igem-shit/`, run 2026-08-29 against the repository at commit `a8932f5`
(branch `hardening`). Nothing was copied into the repository. Every claim below was read off
the files themselves; where a number was extracted by parsing a binary the extraction is
described so it can be repeated.

**Headline, stated first because it is a negative:** *no carotenoid measurement of any kind
exists on this machine.* Neither M8 (phytoene) nor M5 (β-carotene past 1.25 mg/gDCW) is
unblocked, and the reason is structural rather than accidental — this lab never ran a
carotenoid strain. What the survey did find is a different, cheaper unblock: the
**reporter-free autofluorescence measurement that `docs/DATA_INVENTORY.md` says "has never
been measured, on any plate" is present, in full, inside an instrument file the repository
has never opened.**

---

## 1. The directory

`~/Desktop/igem-shit/` — 236 files, three data trees plus unrelated media.

| Tree | Contents |
| --- | --- |
| `iGEM 2026/Experiments/` | **22 `.xpt` Gen5 instrument files** — the plate reader's own format, one per run, 2026-06-12 to 2026-08-08 |
| `iGEM 2026/Protocols/` | 3 `.prt` Gen5 protocol definitions |
| `iGEM 2026/Results (excel sheets)/` | 25 `.xlsx` exports |
| `Result & Analysis/Plate Reader Result/` | 20 `.xlsx` RAW/ANALYSED pairs at top level, plus `Biosensor Testing/` (5) and `Characterising Plasmid Loss/` (7) |
| `Result & Analysis/rt-qPCR results/` | 4 run families, ~60 Bio-Rad CFX exports + 2 `.pcrd` + 3 zips |
| `Result & Analysis/Gel Doc/` | 47 gel images (`.tif`/`.png`), colony PCR and construction checks |
| `Result & Analysis/Sequencing Result/` | 3 zips of Sanger `.ab1`/`.seq` — GFP-RFP, URA, 1xUPRE1/4xUPRE2 |
| `Result & Analysis/Wetlab results.pptx` | 90-slide summary deck |
| everything else | promo `.mp4`s, two textbook PDFs, a contacts CSV, a timetable — not data |

`Result & Analysis-20260821T102706Z-1-001.zip` is a Google-Drive export of the
`Result & Analysis/` tree beside it; spot-checked entries are the same files.

---

## 2. No carotenoid data exists. M8 and M5 stay blocked.

Searched by filename and by Spotlight full-text over the whole home directory for `HPLC`,
`carotene`, `phytoene`, `lycopene`, `chromatogram`, and for the chromatography file
extensions `.d`, `.cdf`, `.ch`.

- **Zero hits inside `~/Desktop/igem-shit/`.** No chromatogram, no integration report, no
  extraction record, no absorbance scan at 285 nm or 450 nm.
- Every `HPLC` / `carotene` / `phytoene` hit elsewhere on the machine is **this repository's
  own prose and code**, sibling software projects (`Models/`, `dente/`, `syn-vision/`,
  `Validation-Notebooks/`, `some-shit-.../studio`), literature PDFs in `~/Downloads/iGEM/papers/`,
  or a wiki illustration named *"HPLC and Gas Chromatography machines.png"*. Not one is a
  measurement.
- `Wetlab results.pptx` contains **no** occurrence of `caroten`, `HPLC`, `lycopene`,
  `phytoene`, `chromato`, `LC-MS` or `GC-MS` across all 90 slides.
- Slide 87 photographs the −80 °C glycerol-stock box and lists its contents:
  `AFL`, `AFL control`, `No-LacI`, `UPRE1`, `UPRE2`, `NY1`, `AY1`, `ACS1`. **There is no
  carotenoid strain in the collection**, high-producing or otherwise.

So:

- **M8 (phytoene, run order 1) — not unblocked.** `docs/MEASUREMENTS_NEEDED.md` says
  "re-integrate HPLC runs that already exist". Those runs exist in *Elizondo 2025*, in
  someone else's freezer. They do not exist here, and no local re-integration is possible.
- **M5 (β-carotene ceiling, run order 2) — not unblocked, and worse than unmeasured.** The
  spec's precondition is "one chemostat run, **if a high producer already exists**". No
  β-carotene-producing strain exists in this lab at all. M5 needs a strain before it needs a
  chemostat.

This is consistent with `docs/SECOND_DATASET_HUNT.md` and the memory note
`igem-second-dataset-negative`; it extends them from "no second published chemostat series"
to "no local carotenoid material either".

---

## 3. What *is* unblocked: the autofluorescence measurement

`docs/DATA_INVENTORY.md` states, under **"Autofluorescence has never been measured, on any
plate"**:

> `20260807_BY4741_ER_Oxidative.xlsx` — 33 wells, **OD600 channel only** […] Neither was
> read in the mCitrine channel. The reporter-free control exists and the quantity it was
> made to measure was not measured.

**The xlsx export is OD-only. The instrument file is not.**

`~/Desktop/igem-shit/iGEM 2026/Experiments/20260807_BY4741_ER_Oxidative.xpt` is an OLE2
compound document (Gen5 experiment). Its `SUBSETS/1/DATA` stream holds **4800 readings in
24-byte records** — two channels × 96 wells × 25 timepoints:

| block | n | range | identity |
| --- | ---: | --- | --- |
| A | 2400 | 0.0874 – 0.7014 | OD600 |
| B | 2400 | **44 – 222** | **mCitrine, 480/530** |

All 2400 mCitrine values are non-zero. The `Contents` stream names the protocol the run
used: `4h-10min_mCitrine_OD600.prt` — **the same protocol file as
`20260804_ER_Oxidative_Replicate4.xpt`**, i.e. the reporter-free strain was read in the
mCitrine channel *at the same settings as the sensor plates*, which is exactly the condition
`observation.py`'s docstring imposes ("on an isogenic reporter-free strain, not by
subtracting a guess").

For scale: on `20260804` the reporter constructs read A1 = 378 → 2564 RFU over the same 4 h.
BY4741 sits at 44 – 222 RFU. That is the `autofluorescence` term, and it has been sitting in
the instrument file since 7 August.

**This is item 1 of `DATA_INVENTORY.md`'s own "what to do, in cost order"** — described there
as "the cheapest item here", costing "one plate". It costs no plate. It costs one export.
Three of the reported dose-response readings — UPRE1 at 2 mM and 5 mM DTT, UPRE2 at 5 mM
DTT, the rows whose verdict is most sensitive to the assumed background — depend on it.
(Where they land in the sweep is what `a` decides; whether they show no response is not,
and that reading has since been withdrawn from the registry.)

### How the extraction was validated

The record layout was not guessed. On `20260804_ER_Oxidative_Replicate4.xpt` the same parse
yields 4800 records; `604` of the `1340` distinct mCitrine values in the committed
`data/plates/20260804_ERandOxidativeStress_Replicate4__mCitrine.csv` appear verbatim in the
`.xpt` block, and `141/408` of the OD values do. The *dis*agreement is the confirmation:
the `.xpt` mCitrine block runs **46 – 4236** while the committed csv runs **−15 – 3898**, and
the csv's H1–H3 blank wells read −3, −1, 4. The export is blank-subtracted; the instrument
file is not. Which is precisely what `DATA_INVENTORY.md` says, and what makes the `.xpt` the
only route to the background.

Parse recipe, for whoever writes the extractor: `olefile` → `SUBSETS/1/DATA`; the 8-byte
sequence occurring exactly 4800 times is the record delimiter; the reading is the
little-endian `double` starting 9 bytes before each delimiter; records are OD600 block then
mCitrine block, 96 wells × 25 timepoints each. **A Gen5 re-export should still be preferred**
— it carries well labels and timestamps, which this parse infers from position.

### The same file settles the other open item

`DATA_INVENTORY.md`: *"20260804 stays out […] The raw data exists only in the .xpt instrument
file, so **re-exporting that one file takes the panel to n=4** — which is the cheapest
remaining gain in the project."* That file is
`~/Desktop/igem-shit/iGEM 2026/Experiments/20260804_ER_Oxidative_Replicate4.xpt`. It exists,
it holds all 96 wells (not the 87 exported, and not the 21 of the superseded table), and its
mCitrine block is un-subtracted. The claim is correct and the file is now located.

All 22 `.xpt` files are raw in this sense. `README.md:156` and `DATA_INVENTORY.md:57` are the
only two places the repository mentions `.xpt` at all, and neither records where they live.

---

## 4. Cross-reference: what the repository already knows

| Found | Repository status |
| --- | --- |
| 4 NewProtocol biosensor xlsx (`20260722/28`, `0803`, `0804`) | **known** — committed as text under `data/plates/`, with sha256, via `manifest.csv` |
| `20260807_BY4741`, `20260814_BY4741` xlsx | **known** — `data/plates/20260814_..._OD600.csv`; discussed in `DATA_INVENTORY.md` |
| ATP-sensor xlsx (`20260805`, `20260809`) | **known** — vendored at `data/atp_sensor/`, `docs/ATP_SENSOR.md` |
| qPCR: 24 Jul, 11 Aug, 13 Aug ER&Ox replicates | **known** — all three analysed in `docs/research/G4_ANCHOR.md` §"2.41 (24 Jul), 3.37 (11 Aug), 4.08 (13 Aug)" |
| `ER&OX-summary.xlsx` crosstalk workbook | **known, deliberately uncommitted** — `data/crosstalk/SOURCE.md`. Note: it is not in `igem-shit/`, it is at `~/Downloads/ER&OX-summary.xlsx` |
| **22 `.xpt` instrument files** | **not known** — mentioned twice in prose, never located, never parsed |
| **3 `.prt` protocol files** | **not known** — the gain/filter definitions behind every plate |
| AFL / **auto feedback loop** plates — a LacI autorepressor, with `AFL control` and `No-LacI` (20 xlsx across the two trees + deck slides 1–37, and four 24-hour `.xpt` runs) | **not known** — no mention anywhere in `docs/` or `README.md`. Described in [`AUTO_FEEDBACK_LOOP.md`](AUTO_FEEDBACK_LOOP.md) |
| Plasmid-loss plates (7 xlsx) + plasmid-loss qPCR (2 runs) | **not known** |
| 47 gel-doc images | **not known** |
| 3 Sanger sequencing zips (incl. `1xUPRE1_and_4xUPRE2`) | **not known** — the construct-identity evidence for the two ER sensors |
| `Wetlab results.pptx`, 90 slides | **not known** — includes the written protocols (slides 84–86) and the strain inventory (87–89) |

The AFL, plasmid-loss, gel and sequencing material is genuinely outside the model's current
scope; it is listed so the gap is deliberate rather than accidental.

---

## 5. A path bug this survey found

`src/ystwin/paths.py` resolves two datasets by convention, and **both conventions are now
stale**:

```python
biosensor_plates():  IGEM_ROOT / "plates" / "Biosensor Testing"
                     ~/Desktop/Result & Analysis/Plate Reader Result/Biosensor Testing
                     ~/Desktop/iGEM 2026/Results (excel sheets)
qpcr_raw_dir():      IGEM_ROOT / "qpcr" / "ER and Oxidative Stress Biosensors"
                     ~/Desktop/Result & Analysis/rt-qPCR results/ER and Oxidative Stress Biosensors
```

Checked: `IGEM/plates` and `IGEM/qpcr` **do not exist**; `~/Desktop/Result & Analysis` and
`~/Desktop/iGEM 2026` **do not exist**. Every one of those directories now lives one level
deeper, under `~/Desktop/igem-shit/`. So on this machine both resolvers return `None` unless
`YSTWIN_PLATES` / `YSTWIN_QPCR_RAW` are set.

This is not currently a wrong-answer bug — the module is built to return `None` rather than
guess, and the biosensor consumers fall through to the committed text in `data/plates/`.
But the `~/Desktop/...` fallbacks are dead weight that read as live, and they are the same
class of thing `audit_reproducibility.py` already removed from `igem_results()` and
`crosstalk_workbook()`: **a path found by convention that resolves on one machine — and this
one no longer resolves even there.** The consistent fix is to delete both Desktop fallbacks
and let the environment variables be the only route, exactly as `crosstalk_workbook()` does.

---

## 6. Privacy

62 Office files under `igem-shit/` were inspected for `dc:creator` / `cp:lastModifiedBy`.
**23 carry a non-generic personal name.** They are named in this report only by filename, no
name is transcribed, and none was copied. `~/Downloads/ER&OX-summary.xlsx` carries both
fields populated, confirming `data/crosstalk/SOURCE.md`'s stated reason for keeping it out.

The Bio-Rad qPCR `Run Information` exports record `Created By User,admin` and instrument
serial `795BR04057` — no personal name. The `.xpt` files embed the Windows path
`C:\Users\NUS\Desktop\iGEM 2026\Protocols\…`, a shared lab account rather than an individual.
**The `.xpt` and `.prt` files are, on this evidence, free of personal names** — which matters,
because unlike the workbooks they could be vendored.

---

## 7. What to do

1. **Export the mCitrine channel of `20260807_BY4741_ER_Oxidative.xpt`.** No bench time.
   It is `DATA_INVENTORY.md`'s own cheapest item and it decides three of eight reported rows.
2. **Export the raw mCitrine channel of `20260804_ER_Oxidative_Replicate4.xpt`** — the n=4
   panel the repository already identified and could not reach.
3. **Repoint or delete the stale `~/Desktop/...` fallbacks in `paths.py`** (§5).
4. **Record in `MEASUREMENTS_NEEDED.md` that M8 and M5 have no local material** — M8's
   chromatograms are Elizondo's, and M5 has no strain, not merely no measurement.
