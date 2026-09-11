# `data/phb/` — index

**Two different papers put files in this directory.** They are the same product and the same
`data/pathways/phb.toml` spec, but different experimental axes, and each has its own
provenance document. Read the one that matches the file you are using.

| dataset | files | provenance | axis |
| --- | --- | --- | --- |
| **Kocharin 2012**, AMB Express 2:52, PMID `23009357` | `kocharin2012_states.tsv`, `kocharin2012_text_ratios.tsv`, `kocharin2012_derived.tsv` | **`SOURCE_kocharin2012_genotype.md`** | genotype: four strains, batch flask and bioreactor |
| **Kocharin 2013**, PMID `23514405` | `kocharin2013_chemostat_states.tsv` | *written by whoever added that file — not covered by `SOURCE_kocharin2012_genotype.md`* | environment: carbon source × dilution rate, chemostat |

This index file is deliberately thin so that either contributor can rewrite it without
destroying the other's provenance. **The substance lives in the per-paper documents.**

## One thing worth recording across both

Both datasets express PHB per **3-hydroxybutyrate repeat unit at 86.09 g/mol**, which is the
basis `data/pathways/phb.toml` declares. In `kocharin2013_chemostat_states.tsv` the ratio
`phb_mg_per_gdw / phb_mmol_per_gdcw` is 86.09 on every row, arrived at independently. The
polymer itself has no molar mass — KEGG `C06143` gives `(C4H6O2)n` — so agreeing on the
repeat unit is what makes the two tables comparable at all.
