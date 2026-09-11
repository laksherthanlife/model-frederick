# Talemi 2016 -- the yeast osmo-stat, and the Slt2/CWI arm inside it

Talemi SR, Tiger CF, Andersson M, Babazadeh R, Welkenhuysen N, Klipp E, Hohmann S, Schaber J
(2016). "Systems Level Analysis of the Yeast Osmo-Stat." *Scientific Reports* **6**:30950.
doi:10.1038/srep30950 -- PMID 27515486 -- PMC4981887.

`src/ystwin/mech/cell_wall_talemi2016.py` is a hand transcription of the SBML deposit below.

## Licence

Two licences, and they are different files.

* **The encoded model is CC0 1.0.** `MODEL1606100000_url.xml` carries, in its own RDF block,
  `http://creativecommons.org/publicdomain/zero/1.0/` with the BioModels dedication: "To the
  extent possible under law, all copyright and related or neighbouring rights to this encoded
  model have been dedicated to the public domain worldwide." Redistribution is unconditional.
* **The article and its supplementary files are CC BY 4.0**, declared in the article's own
  `<permissions>` block: "This work is licensed under a Creative Commons Attribution 4.0
  International License." Redistribution is permitted with the attribution above. This covers
  `PMC4981887.xml`, `srep30950-s2-OsmoStat_Final.xml` and the two text extracts.

This is a stronger position than `jalihal2021/`, which carries no licence, is custody-only and
is gitignored. Nothing here is gitignored.

## Files

| file | sha256 | what it is |
|---|---|---|
| `MODEL1606100000_url.xml` | `c4937f14c2d8e36b13a449b5ac9771d4e17426543298d7a2c9ca9c4df835e5dc` | BioModels MODEL1606100000, "Talemi2016 - Yeast osmo-homoestasis", SBML Level 2 Version 4. **The transcription source.** 11 species, 12 reactions, 100 global parameters, 8 function definitions, 29 assignment rules and 3 rate rules |
| `srep30950-s2-OsmoStat_Final.xml` | `6d18c0ddaf42c9ea3d5d18d42e14157b833be5284ed639fd8d9d85e7a9edfde3` | the authors' own SBML, Supplementary file S2, `Model_COPSI_SBML/OsmoStat_Talemi et al. 2016_Final.xml`. Kept because it establishes the deposit's identity -- see below |
| `PMC4981887.xml` | `ff98279a6e3d443c761b5395dbb746f793a57ae5d0c27171d62676a03f7272c6` | article full text (NCBI efetch, `db=pmc&id=4981887&retmode=xml`) |
| `srep30950-s1-p21-33-tables-S4-S8.txt` | `15d4fc1d551ba8eadab9a53be86f06ecd27b557245c90f5765941fe0e96e2175` | **derived, ours**: `pdftotext -layout -f 21 -l 33` over the supplementary PDF. Tables S4 (the ODE system), S5 (the rate laws), S6 (initial conditions), S7 (auxiliary variables) and S8 (the parameter table) |
| `srep30950-s1-p06-07-calcofluor-module.txt` | `93e1b248685239cf6a8d2ca5351a29bdd55e22b2caef540dc9c5daf5fdfb88b7` | **derived, ours**: `pdftotext -layout -f 6 -l 7` over the same PDF. The "Calcofluor mediated Slt2 activating module", which is the only cell-wall-damage input the source has and which the deposit does not contain |

Retrieval, in the order it was done:

```
curl -L -o talemi.zip https://www.ebi.ac.uk/biomodels/model/download/MODEL1606100000   # then unzip
curl -o PMC4981887.xml 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=4981887&retmode=xml'
curl -o supp.zip https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4981887/supplementaryFiles
```

## The supplementary PDF is pinned and NOT vendored

`srep30950-s1.pdf`, sha256
`97e1c32b05f447767a40f8a51a00684b5f6885c03b7647ea049baadd85030a91`, 12,491,506 bytes, is inside
the `supplementaryFiles` archive above. It is **not** copied here: at 12.5 MB it is larger than
this whole directory tree (`data/native_reference_models` is 2.8 MB) and would become the single
largest tracked file in the repository. The two text extracts carry the five tables the
transcription is scored against; they were produced with poppler `pdftotext` 26.04.0 and the
page ranges above, so anyone with the pinned PDF can regenerate and diff them.

The extracts are **ours**, not publisher bytes, and are labelled so. They are corroboration
rather than the load-bearing artifact: every Table S6 and Table S8 row the module quotes is
independently confirmed by recomputing the deposit's own stored value from it -- a
misread formula does not reproduce a 15-significant-figure initial condition by accident.

## Why the authors' own SBML is kept as well

**The accession the paper prints does not exist.** The Supplementary Information's "Simulation
Instructions" says the selected model is in BioModels under "access identifier
MODEL1604100004". `https://www.ebi.ac.uk/biomodels/MODEL1604100004` answers HTTP 404. The model
is actually MODEL1606100000, which no page of the article names.

So the deposit's identity was established by comparison instead, against the authors' own
Supplementary S2 SBML, and it holds exactly:

* the same 4 compartments, the same 11 species with the same `initialAmount`,
  `boundaryCondition` and `hasOnlySubstanceUnits`, and the same 100 global parameters;
* 13 of those 100 parameter *strings* differ and **no value does** -- every one is a
  scientific-notation reformat (`1E-9` against `1e-09`, `6.022E23` against `6.022e+23`,
  `5.02026122077738E5` against `502026.122077738`).

Both files are kept so that this comparison stays runnable rather than being a claim in a
document.

## Scope, stated here because the directory name does not say it

The transcription is registered as axis **`cell_wall_slt2`**, not `cell_wall`. The deposit
contains **no Rlm1, no RLM1 box, no transcription of any kind, and no cell-wall-damage input**:
`grep -ci` over `MODEL1606100000_url.xml` returns 0 for each of Rlm1, congo, caffeine,
calcofluor, caspofungin, zymolyase, Pkc1, Bck1, Wsc1, Mid2 and "transcription". What it
contains is Slt2 activation driven by **cell volume**, and the direction is the one that reads
backwards: the paper's CWI arm answers a volume **increase** (hypo-osmotic adaptation) while
HOG answers a volume decrease. `mech/cell_wall_talemi2016.py` refuses the missing input and the
missing output by name, each as a `Param.refused`.
