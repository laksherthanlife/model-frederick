# Replacing the ER stress anchor

Research note. Companion to `G4_STATISTICS.md` (which settles that the existing ER data cannot be
rescued) and `G4_ANCHOR.md` (which settles what the existing data can and cannot say). This note
answers only the forward question: **what to measure instead, and how.**

Every primer sequence quoted below was located in the S288C reference sequence computationally
before being recommended — coordinates, amplicon sizes and uniqueness counts are given so the work
can be repeated. Every citation was resolved against PubMed/PMC. Claims I could not verify against
primary text are marked ⚠.

---

## 1. Verdict

**Switch the ER anchor to `KAR2`**, normalised to the geometric mean of **`TAF10`, `ALG9` and
`TFC1`**, and run **`HAC1i` exon-junction qPCR on the same cDNA** as a genomic-DNA-proof witness
that the pathway fired.

Total `HAC1` cannot report the UPR: Ire1 converts the unspliced message *directly* into the spliced
one, so the sum barely moves. Its intron is a poor genomic-DNA guard, because `HAC1u` is a stable
cytoplasmic mRNA identical to the locus. `KAR2` is the native analogue of a UPRE reporter, is
abundant, and its heat-shock element is inert under DTT at 30 °C.

The largest single change, though, is the dose ladder: five doses to 2 mM, sampled at 45 min, is
worth about **20× the replicates**.

---

## 2. Why `HAC1` fails, and whether any `HAC1` assay is designable

### 2.1 The splicing biology, from the primary literature

`HAC1` encodes the bZIP activator of the UPR (Cox & Walter 1996). Its message carries a **252-nt
intron** that is *not* removed by the spliceosome: Ire1's cytoplasmic endoribonuclease domain
cleaves both splice sites, and tRNA ligase Trl1/Rlg1 joins the exons (Sidrauski, Cox & Walter 1996;
Sidrauski & Walter 1997; Kawahara et al. 1997, 1998). I confirmed the intron length independently
from SGD: `YFL031W` is 969 nt of CDS+intron, CDS 1–661, **intron 662–913 (252 nt)**, CDS 914–969
(exon 2 is only **56 nt** of coding sequence).

The intron is a *translational* switch, not a maturation step. Long-range base pairing between the
5′ UTR and the intron blocks translation of `HAC1u` (Rüegsegger, Leber & Walter 2001; Chapman &
Walter 1997), and unspliced message that escapes that block is silenced again post-translationally
(Di Santo, Aboulhouda & Weinberg 2016). Consequently `HAC1u` is not a transient nuclear pre-mRNA —
**it is an abundant, stable, cytoplasmic mRNA in unstressed cells.** That single fact is what
destroys the standard qPCR design (§2.3).

### 2.2 Is total `HAC1` mRNA constant during the UPR? Yes — quantified, twice

The strongest primary statement is Kawahara, Yanagi, Yura & Mori (1997), *Mol Biol Cell* 8:1845,
tracking the 1.4 kb (unspliced) and 1.2 kb (spliced) species by Northern blot after tunicamycin:

> "The sum of the two mRNA species remained almost constant, suggesting that 1.4-kb mRNA was
> directly converted to 1.2-kb mRNA."

Leber, Bernales & Walter (2004), *PLoS Biol* 2:e235, checked the same thing with 6 mM DTT and
1 µg/mL tunicamycin and found "rapid and efficient splicing" with no boost in abundance — their
section heading is "**ER-Distal Secretory Stress** Boosts `HAC1` mRNA Abundance", and the 3–4×
abundance increase they report occurs only in *sec12-1*, *sec14-1* and *sec1-1* strains shifted to
the non-permissive temperature. Under DTT — the team's stressor — total `HAC1` does not move.

So the repo's reasoning was right, and now it has numbers behind it: a total-`HAC1` assay is
measuring a quantity the pathway holds constant by construction. **Even a perfect, gDNA-free
total-`HAC1` assay would have returned a flat anchor.** The 91% gDNA problem and the wrong-transcript
problem are genuinely independent faults, exactly as `G4_ANCHOR.md` states.

### 2.3 The crux: does non-spliceosomal splicing break intron-spanning primer design?

**Yes for the textbook design, and for a sharper reason than "the intron is unusual".**

The textbook gDNA fix places primers in two exons flanking an intron, so that genomic template gives
a longer product (or none). It works because the only intron-containing species in the cell is a
short-lived nuclear pre-mRNA. For `HAC1` that assumption is false: `HAC1u` is a *bona fide, abundant,
cytoplasmic* mRNA with the intron in place. An exon-flanking pair therefore yields **the same long
amplicon from genomic DNA and from the dominant mRNA isoform** — size cannot separate them. And the
usual escape hatch (the genomic product is too long to amplify) is closed: 252 nt is trivially
amplifiable. This is exactly what the Kimata lab's competitive RT-PCR exploits — one primer pair,
two band sizes — and it is why that assay is run on gel, not as SYBR qPCR.

Concretely, using the Kimata-lab competitive pair (Geronimo et al. 2025, Table 1), located in the
reference sequence:

| Template | Amplicon |
|---|---|
| spliced `HAC1i` | **396 bp** |
| unspliced `HAC1u` | **648 bp** |
| **genomic DNA** | **648 bp — indistinguishable from `HAC1u`** |

**But no for the junction design, and this is the answer that matters.** The exon1–exon2 junction
sequence created by Trl1 does not exist anywhere in the genome. A primer that crosses it is
genome-absent, and therefore gDNA-insensitive by the ValidPrime paper's own criterion:

> "qPCR assays can be designed to be gDNA insensitive, such as those designed to target exons
> flanking a long intron **or with primers that cross exon–exon junctions**, [whereas] qPCR assays
> for single-exon genes will readily amplify contaminating gDNA."
> — Laurell et al. 2012, *Nucleic Acids Research* 40(7):e51

### 2.4 Has anyone published working splice-specific `HAC1` qPCR primers? Yes — three labs

| Assay | Primers (5′→3′) | Where it binds | gDNA? |
|---|---|---|---|
| **`HAC1i`** (Hata, Ishiwata-Kimata & Kimata 2022) | F `ACCTGCCGTAGACAACAACA` / R `ACCTGACTGCGCTTCTGGAT` | F exon 1 @530–549; R spans junction (14 nt exon 2 + 6 nt exon 1), found **only** in spliced sequence @655–674. Amplicon **145 bp** | **immune** |
| **`HAC1i`** (Gast et al. 2021) | F `GCGTCGGACCAAGAGACTTC` / R `CTGACTGCGCTTCTGGATTAC` | F exon 1 @459–478; R spans junction, spliced-only @652–672. Amplicon **214 bp** | **immune** |
| **`HAC1u`** (Gast et al. 2021) | F `CAATTGGCGTAATCCAGCCG` / R `AGCTGGGGCTAGTGTTCTTG` | F exon1→intron @644–663; R inside intron @694–713. Amplicon **70 bp** | amplifies gDNA |
| **total `HAC1`** (Hata et al. 2022) | F `GCGTCGGACCAAGAGACTT` / R `TCGTCGACTCTGGTACATTTTC` | both inside exon 1, @459–477 / @504–525. Amplicon **67 bp** | amplifies gDNA |
| **`HAC1u` vs `HAC1i`, gel** (Geronimo et al. 2025) | F `TACAGGGATTTCCAGAGCACG` / R `TGAAGTGATGAAGAAATCATTCAATTC` | F exon 1 @318; R exon 2. 396 bp spliced / 648 bp unspliced | 648 bp band = gDNA |

Independent verification I ran: searching the whole SGD ORF set (6039 ORFs, both strands), each of
the eight non-junction primers above occurs **exactly once**, while both junction primers
(`ACCTGACTGCGCTTCTGGAT`, `CTGACTGCGCTTCTGGATTAC`) occur **zero times**. That zero is the point — it
is what "gDNA-proof" means, demonstrated rather than asserted.

A fourth lab (Matabishi-Bibi et al. 2022, *Nat Commun* 13:6331) runs the same three-assay scheme
(`HAC1u` / `HAC1i` / `HAC1Tot`, with one primer spanning the exon–intron junction and one the
exon–exon junction) and, importantly, **published the specificity control you should copy**: they
validated the junction primers on defined mixtures of two plasmids carrying the unspliced and
spliced forms at u:s ratios of 1:0, 0.75:0.25, 0.5:0.5, 0.25:0.75, 0:1. That titration is the only
honest way to show a junction primer is not cross-priming on the large `HAC1u` excess, and it costs
two plasmids and one plate row.

### 2.5 So why is `HAC1i` not my anchor?

Three reasons, and none of them is "it can't be measured".

1. **Shape mismatch breaks the gate's own criterion.** `HAC1i` is the *switch upstream* of the
   reporter. Splicing efficiency runs from ≈0 in unstressed cells to 0.955–1.0 under stress (Tehfe
   et al. 2021, RNA-seq/ribosome-profiling re-analysis, Table 1) — near-binary. `anchor_agreement`
   scores a **Pearson r across doses**, and a saturating anchor beside a graded reporter gives a poor
   r *even when everything is working correctly*. `KAR2` sits downstream of Hac1p on the same UPRE
   step as the reporter, so their dose–response shapes match by construction.
2. **Undefined baseline.** With splicing ≈0 at 0 mM, the control `Cq` is at or past the limit of
   detection, and `delta_delta_cq` divides by it. A splicing *percentage* is well defined at zero; a
   *fold change* is not.
3. **Low copy number.** `HAC1` is a lowly expressed transcript; `G4_STATISTICS.md` already notes that
   technical Cq SD rises from 0.1–0.3 to 0.5–2 cycles above Cq ≈ 30–35, where Poisson sampling
   dominates. `KAR2`/BiP is one of the most abundant ER proteins and its message amplifies early.

**`HAC1i` is nonetheless mandatory as a control**, because it answers a question `KAR2` cannot: *did
Ire1 fire at this dose at all?* It is the only yeast assay in this whole experiment that is immune to
genomic DNA by construction, so it also serves as a DNase-independent cross-check — if `HAC1i` moves
and `KAR2` does not, the fault is downstream (see §7.6), not in the RNA.

---

## 3. Anchor candidates

Intron status was computed from `SGD_features.tab` (retrieved 2026-08-26) rather than looked up:
**282 of 5783 verified ORFs carry an intron — 4.9%.** Not one ER or oxidative candidate below has
one. Intron-spanning design is simply unavailable, which is §4's whole problem.

| Gene | Systematic | Intron | Published induction (verified) | Independence risk | Verdict |
|---|---|---|---|---|---|
| **KAR2** (BiP) | YJL034W | none | Ire1-dependent induction by 1.0 µg/mL tunicamycin, 30 min, confirmed by RT-qPCR in *IRE1⁺* but not *ire1Δ* (Geronimo et al. 2025, Fig. 3C). Kohno et al. 1993: BiP mRNA "synthesized at a high basal rate and further induced by … unfolded proteins", and "induced severalfold by heat shock"; UPRE + HSE mapped to −245…−9. ⚠ magnitudes are shown graphically, not stated numerically, in both papers | Carries **both** a UPRE and an HSE → heat/ethanol-inducible independently of Ire1 (16% ethanol induces `KAR2` equally in *IRE1⁺* and *ire1Δ*, Geronimo 2025). **Not a confound here**: DTT, tunicamycin and inositol depletion do *not* induce `HSP104` (ibid., Fig. S5), i.e. they do not fire the HSR at 30 °C | **RANK 1 — use it** |
| **SIL1** | YOL031C | none | Ire1-dependent Tm induction; **not** induced by ethanol in either *IRE1⁺* or *ire1Δ* (Geronimo 2025, Fig. 3E) | Cleanest UPR-specificity in the tested panel — no HSE reported | **RANK 2** — the choice if you ever compare stressors |
| **JEM1** | YJL073W | none | as SIL1 (ibid., Fig. 3F) | as SIL1 | RANK 3 |
| **ERO1** | YML130C | none | Ire1-dependent Tm induction (ibid., Fig. 3D) | Also carries an HSE (Kohno 1993; ⚠ Takemori 2006, cited by Geronimo 2025, not independently retrieved) — same caveat as KAR2 with none of KAR2's abundance advantage | RANK 4 |
| **EUG1** | YDR518W | none | Member of the canonical UPR target set; qPCR-confirmed UPR target (Cui et al. 2019) ⚠ magnitude graphical only | Low basal expression → high Cq | RANK 5 |
| **LHS1** | YKL073W | none | as EUG1 | ER Hsp70; also induced in `gas1Δ` basally | RANK 5 |
| **PDI1** | YCL043C | none | **⚠ bibliographic trap.** The "2–2.5-fold" figure in `generator/stress_panel.py` traces to Sarkar, Paira & Das 2018 *NAR* 46:1139, whose sentence is: "*Examination of the steady-state levels of two ER chaperone messages, BiP/KAR2 and PDI, revealed a consistent 2- to 2.5-fold upregulation in the presence of tunicamycin **in the exosome/DRN-deficient rrp6-Δ and cbc1-Δ strains***". That is a **mutant-vs-wild-type** ratio, not a wild-type stress induction. Do not quote it as the latter | — | RANK 6 — and fix the annotation |
| **FKB2** | YDR519W | none | canonical set ⚠ | — | RANK 7 |
| **SCJ1** | YMR214W | none | Hac1-induced (Ishiwata-Kimata et al. 2025, gene list) | — | RANK 8 |
| **HRD1** | YOL013C | none | ERAD UPR target (Travers et al. 2000) | ERAD arm induction is modest | RANK 9 |
| **INO1** | YJL153C | none | Largest dynamic range of the canonical set | **Inositol-repressed.** SD/SC media contain inositol, so basal and induced levels are set by medium composition, not by ER stress | **Reject** |
| **HAC1 total** | YFL031W | 252 nt, non-spliceosomal | Sum of isoforms "almost constant" (Kawahara 1997); no abundance boost under DTT/Tm (Leber 2004) | — | **Reject — this is the current failure** |
| **HAC1i** (junction) | YFL031W | — | Splicing efficiency ≈0 unstressed → 0.955–1.0 induced (Tehfe 2021); ~2× rise on a 30→39 °C shift, sustained ≥12 h (Hata 2022) | Upstream of the reporter; switch-like | **Adopt as pathway control, not as the correlation anchor** (§2.5) |

⚠ **On Travers et al. 2000.** Cell blocks automated retrieval and I could not parse the full text or
supplement. I therefore cite it for what secondary sources agree on and for its design — DTT and
tunicamycin, 60 min, fold change measured against untreated, with *ire1Δ* and *hac1Δ* strains used to
subtract UPR-independent changes — but **the frequently-quoted "canonical seven" set (`KAR2`,
`EUG1`, `PDI1`, `LHS1`, `FKB2`, `ERO1`, `INO1`) is marked ⚠ here**: I saw it only in secondary
descriptions, never in the paper's own words. If the team needs that list in a report, retrieve the
PDF by hand first.

### Why `KAR2` wins

- **It is the endogenous version of the thing being tested.** The biosensor is a synthetic
  UPRE→mCitrine cassette. `KAR2` is a native UPRE gene whose UPRE was mapped by mutagenesis (Kohno
  et al. 1993). Asking "does the synthetic UPRE track the native UPRE?" is the sharpest possible
  form of gate G4, and it makes the two dose-response *shapes* commensurable.
- **Abundance buys precision for free.** `G4_STATISTICS.md` shows the precision requirement is
  σ-limited: at σ = 0.3 cycles a ≤20% relative SD is unachievable at any margin. Cq noise is a
  function of copy number, and BiP is abundant. The failed `Hac1 Exon` pair amplified 3.72 cycles
  later than `UBC` on identical −RT template (13× worse); a `KAR2` assay should land at or before
  the reference.
- **The HSE objection does not bite under DTT.** The repo's own note ("KAR2 carries both a UPRE and
  an HSE, so it is heat-inducible") is correct and primary-sourced. But the confound requires the
  HSR to fire, and Geronimo et al. 2025 show directly that DTT, tunicamycin and inositol depletion
  leave `HSP104` uninduced while ethanol induces it drastically. At constant 30 °C with DTT, `KAR2`'s
  HSE is inert. **If the team ever runs the heat arm of `stress_panel.py`, switch to `SIL1`.**
- **Primers already exist and are published** (§7.2).

---

## 4. Killing genomic DNA in a yeast RNA prep

### 4.1 The structural problem, quantified

**4.9% of verified yeast ORFs have an intron** (282/5783, computed above). Every gene in §3 is
single-exon. Yeast RNA prep also *starts* with more genomic DNA than most systems, because lysis is
mechanical (bead-beating / dismembrator) or enzymatic (zymolyase), both of which shear chromatin
into the lysate. So in yeast the textbook fix is unavailable exactly where it is most needed. Four
things to do instead, in order of how much they buy.

### 4.2 Use a real DNase step, and use it in solution

Teste et al. 2009 — the reference-gene paper this note also relies on for §5 — used RNeasy
(Qiagen) **plus an additional on-column RNase-free DNase digestion**, and that is the realistic
iGEM baseline. Two refinements are worth the effort:

- **In-solution beats on-column for high-gDNA input.** On-column digestion is time- and
  enzyme-limited by the column format; an in-solution digest (e.g. Ambion `DNA-free`, which is what
  Laurell et al. 2012 used) lets you extend the incubation and add enzyme. ⚠ I could not retrieve a
  head-to-head yeast benchmark with residual-gDNA numbers; treat this as mechanism-based, and
  measure it on your own samples with the −RT wells you are already running.
- **Do not heat-inactivate DNase I in the presence of Mg²⁺.** DNase I is thermostable when
  magnesium is present, and heating RNA with Mg²⁺ degrades it. Wiame et al. 2000 (*BioTechniques*
  29:252) give the fix: chelate with EGTA, then heat — irreversible inactivation without RNA
  degradation. Inactivation-reagent kits (`DNA-free`) do the same thing by removing the enzyme.
- **Carry-over inhibits RT.** Laurell et al. cap it explicitly: "the volume of DNase-treated RNA did
  not exceed 25% of the total volume during RT." Worth copying verbatim.
- A single-tube alternative used by the yeast UPR labs themselves: **ReverTra Ace qPCR RT Master
  Mix with gDNA Remover** (Toyobo), which is what Geronimo et al. 2025 used for the exact `KAR2` /
  `SIL1` / `JEM1` / `TAF10` assays recommended here.

### 4.3 Design the assay so gDNA cannot be amplified

Only one yeast option exists in this experiment, and it is `HAC1i` (§2.4) — a genuine exon–exon
junction assay in a genome where almost nothing else offers one. That is a second, independent
reason to run it. UTR-spanning primers do **not** help: a 3′-UTR amplicon is as present in genomic
DNA as a coding one. Poly(A) selection reduces but does not eliminate gDNA and costs RNA yield;
Đermić et al. explicitly note that "any procedure implemented to reduce the concentration of DNA in
the sample certainly also causes the reduction of initial concentration of the RNA itself."

### 4.4 The structural fix if DNase is not enough: Đermić et al. 2023

Đermić, Ljubić, Matulić, Procino, Feliciello, Ugarković & Feliciello (*Sci Rep* 2023;13:11470) solve
precisely this class of problem — a target with no intron where gDNA and cDNA are sequence-identical.
The method:

1. Reverse-transcribe with a **gene-specific RT primer carrying four mismatches, alternating with
   matched bases, starting at the 3′ end** ("primer-specific mismatch", PSM). The mismatches become
   part of the cDNA.
2. Use the **same modified primer** in the PCR, at a discriminating annealing temperature (~60 °C).
   It still primes on cDNA (perfect match to the sequence it created) but dissociates from genomic
   DNA.

They report that "no recA sequence amplification was observed using our method unless cDNA was
created by reverse transcription", and that persistent human α-satellite gDNA contamination — a
target with no introns and no length discrimination, i.e. the hardest case — disappeared from −RT
reactions.

Costs, stated plainly: it needs a **gene-specific RT** per target (no random hexamers, no oligo-dT),
the mismatch count and discriminating temperature must be optimised per assay, and it is a 2023
method without much independent replication. **Recommendation: hold it in reserve.** If, after §4.2,
the `KAR2` −RT margin is still under ~3 cycles, this is the principled next move rather than another
round of DNase.

### 4.5 What −RT margin to expect

The number the field actually uses, from the ValidPrime paper itself:

> "A difference of at least five quantification cycles (Cq) between RT(+) and RT(−) reactions
> indicates that <3% of the total signal originates from gDNA, and is commonly used as limit to
> ensure accurate estimation of GOI expression. **Smaller differences typically call for DNase
> treatment of samples.**"

Note what that sentence does and does not say. It is a *statement about what 5 cycles means*, and a
trigger for DNase — not a guideline mandate, consistent with `G4_STATISTICS.md`'s finding that MIQE
never mandated a numeric margin. Practical targets for this experiment:

| Margin | %gDNA | Read as |
|---:|---:|---|
| ≥ 5 cycles | ≤ 3.1% | comfortable; what a good DNase step should deliver on an abundant target like `KAR2` |
| 3–5 cycles | 3–12.5% | acceptable; the repo's existing `_MIN_RT_MINUS_MARGIN` bar |
| < 3 cycles | > 12.5% | re-treat with DNase, or escalate to §4.4 |
| < 0.74 cycles | > 60% | past the ValidPrime ceiling; not a measurement |

The team's `TRX2` assay already achieves 8.09–8.6 cycles (0.3–0.4% gDNA) on RNA from the *same*
preps. **That is proof the RNA prep can clear the bar** — the ER failure was the primer pair and the
target, not the extraction. Expect `KAR2` to behave like `TRX2`.

One diagnostic to keep from `G4_ANCHOR.md` §2.7 and apply to every new pair: run it on the −RT
wells. All these loci are single-copy in a haploid, so on identical genomic template a healthy pair
lands within ~1 cycle of `UBC`/`TAF10`. The old `Hac1 Exon` pair was 3.72 cycles late — 13-fold less
efficient on DNA. Any replacement that fails this test should be redesigned before use.

---

## 5. Normalisers: use three, and not the ones you are using

### 5.1 Why the geometric mean of several

Vandesompele et al. 2002 (*Genome Biology* 3(7):RESEARCH0034) is the origin of the convention and of
geNorm: single reference genes carry their own biological regulation, and the geometric mean of
several validated genes averages that away while resisting outliers (an arithmetic mean in Cq space
*is* a geometric mean in quantity space, which is why the convention is stated that way).

The repo has independently rediscovered why this matters. From `G4_ANCHOR.md` §R4: for
`NativeYap1`, `Var(ΔCq_UBC) = 0.859` against `Var(ΔCq_TRX2) = 0.202` — **the reference gene
contributes four times more variance than the anchor**, and normalising to `UBC` makes the
measurement noisier than not normalising at all. `UBC` swings up to 2.23 cycles across doses within
one construct-replicate. And the *sign* of the pooled oxidative result flips depending on the
normalisation choice. No statistical treatment fixes that; only a validated multi-gene normaliser
does.

Expected gain, from the same section: with independent reference errors `s_R → s_R/√k`, going from
one reference to three cuts SD(ΔΔCq) by ~1.4×, **worth roughly twice the replicates for the price of
two extra plate rows.**

### 5.2 What the yeast literature validated

Teste, Duquenne, François & Parrou (2009), *BMC Molecular Biology* 10:99, is the definitive yeast
study: 13 candidates, geNorm ranking across sample sets spanning growth phases, carbon sources and
strain backgrounds, cross-checked against ~30 SGD microarray experiments. Their conclusions:

- **Stable:** `ALG9`, `TAF10`, `TFC1`, `UBC6` (plus `KRE11`, `FRP2` in some sets).
- **Invalidated, explicitly:** `ACT1`, `PDA1`, `TDH3`, `IPP1`, `RDN18` — "we invalidated and
  discourage the use of `ACT1` as well as other commonly used reference genes … as internal
  controls".
- Their working normalisation factor is `NF(UBC6, TAF10, ALG9)`; for glucose growth,
  `TAF10 + UBC6` alone sufficed (pairwise variation V2/3 = 0.098, well under the 0.15 convention).

Vaudano, Noti, Costantini & Garcia-Moruno (2011), *Biotechnology Letters* 33:1593, independently
evaluated reference genes in *S. cerevisiae* under fermentation with three algorithms and reached the
same headline conclusion — a validated set beats `ACT1`, and using an unreliable single reference
risks "misleading interpretation of expression data". ⚠ I verified the abstract but not the gene
list; do not quote specific Vaudano genes without checking the table.

### 5.3 Recommendation

**Use `TAF10` + `ALG9` + `TFC1`, geometric mean.** Drop `UBC6` for this experiment specifically:

1. Its measured PCR efficiency is the lowest of Teste's stable set (**84%** vs `TAF10` 96%, `ALG9`
   93%, `TFC1` 91%), and `G4_STATISTICS.md` shows efficiency error propagates worse than Cq error.
2. `UBC6` is an **ER-membrane-anchored ubiquitin-conjugating enzyme functioning in ER-associated
   degradation** — Teste's own functional annotation reads "ER-associated protein catabolic
   process". Travers et al. 2000's central finding is that the UPR and ERAD are co-regulated. A
   normaliser sitting inside the pathway under study is a bad normaliser *in principle*, even if
   Teste found it stable across growth conditions (which did not include ER stress). ⚠ I found no
   paper explicitly measuring `UBC6` induction by DTT — this is a mechanism-based objection, not a
   measured one. It costs nothing to avoid.

**Resolve what "UBC" currently means.** `qpcr.py` sets `REFERENCE_TARGET = "UBC"`, which is not a
yeast gene name — *S. cerevisiae* has `UBC1`–`UBC13`. If the logbook's primers are Teste's `UBC6`
pair, objection 2 applies directly. If they are `UBC4` (`YBR082C`), note that `UBC4` is one of the
few yeast genes that **does** carry an intron (95 nt) and could have been made gDNA-proof. Either
way the identity must be read off the ordered oligo sequences before the next run, because it
changes the interpretation of the existing `UBC` instability.

Validated primer pairs, verified in the reference sequence (positions in ORF, both primers unique
across all 6039 ORFs):

| Gene | Systematic | Forward | Reverse | Amplicon | Eff. (Teste) |
|---|---|---|---|---|---|
| `TAF10` | YDR167W | `ATATTCCAGGATCAGGTCTTCCGTAGC` | `GTAGTCTTCTCATTCTGTTGATGTTGTTGTTG` | 141 bp | 96% |
| `ALG9` | YNL219C | `CACGGATAGTGGCTTTGGTGAACAATTAC` | `TATGATTATCTGGCAGCAGGAAAGAACTTGGG` | 162 bp | 93% |
| `TFC1` | YBR123C | `GCTGGCACTCATATCTTATCGTTTCACAATGG` | `GAACCTGCTGTCAATACCGCCTGGAG` | 223 bp | 91% |
| (`UBC6`) | YER100W | `GATACTTGGAATCCTGGCTGGTCTGTCTC` | `AAAGGGTCTTCTGTTTCATCACCTGTATTTGC` | 272 bp | 84% |

A convergence worth noting: the `TAF10` pair Geronimo et al. 2025 use for their UPR qPCR is
*character-for-character* Teste's pair. The yeast UPR field has already adopted the recommendation.

**Implementation.** `qpcr.py::delta_delta_cq` takes `reference_target: str`. It needs
`reference_targets: Sequence[str]` normalising to the arithmetic mean of their Cq — this is already
written up as `G4_ANCHOR.md` §R4 and is not re-derived here.

---

## 6. The oxidative anchor: keep `TRX2`, move the dose

### 6.1 The effect-size argument, and why it does not favour a gene swap

`g4_anchor.py` computes `reps_needed = n × (target_snr / snr)²`. The quadratic is the whole story:
doubling the effect is worth four times the replicates. So the instinct is to hunt for a bigger
gene. Here that instinct is wrong, for three reasons.

**`TRX2` is not broken.** Its −RT margin is 8.09–8.6 cycles (0.3–0.4% gDNA) and it passes QC on 100%
of readings. It is the one assay in this project that already works. Discarding a validated assay to
chase fold change trades a known quantity for an unknown one.

**The obvious big-fold candidates are disqualified by the gate's own logic.**

| Gene | Systematic | Regulation | Verdict |
|---|---|---|---|
| **CTT1** | YGR088W | Msn2/4 via STRE. Huge fold: 30× at 45 min after 2 mM H₂O₂ in early exponential phase, and ~1000× on entry into stationary phase (Vázquez et al. 2017, Fig. 3) | **Reject.** Msn2/4 activity is coupled to *reduced growth rate* — that is the definition of the ESR (Gasch et al. 2000; Brauer et al. 2008 show ESR expression tracks growth rate across nutrient limitations). G4's `growth_confound_r2` test exists to exclude exactly this. An anchor whose dominant regulator is growth rate is the worst possible choice: it would maximise, not resolve, the confound. Its 1000× swing with growth phase makes the point |
| **GRE2** | YOL151W | Msn2/4 + Yap1 + Hog1 | Reject — same growth coupling, plus multi-pathway |
| **TSA1** | YML028W | Yap1/Skn7; very high basal abundance | Weak candidate: high basal compresses fold change |
| **SOD1** | YJR104C | Largely post-translational regulation | Reject — small transcriptional response |
| **GPX2** | YBR244W | Yap1 + Skn7 (Tsuzi et al. 2004a) **and** calcineurin/Crz1 Ca²⁺ signalling (Tsuzi et al. 2004b) | Reject — second, unrelated input |
| **AHP1** | YLR109W | Yap1 target | Plausible second target; no advantage established |
| **TRX2** | YGR209C (ORF only **315 bp**) | Yap1-dependent (Kuge & Jones 1994); Yap1/Skn7 regulon (Lee et al. 1999) | **Keep** |

**And the diagnosis was never the gene.** `G4_ANCHOR.md` §5.1 shows the reporter's H₂O₂ response is
1.68–1.89× at 0.5 mM and the anchor RNA was taken at 0, 0.2, 0.5 mM — the bottom 40% of the
reporter's own range. The anchor was asked to confirm a signal at doses where the *reporter* has
barely moved.

### 6.2 Does 1.0 mM recover more SNR than extra replicates? Yes, by 4.5×

`SE(slope) = σ / √(R · Sxx)` with `Sxx = Σ(x − x̄)²`. Computing `Sxx` for the candidate ladders:

| H₂O₂ ladder | `Sxx` | SE(slope) relative to current | equivalent replicate multiplier |
|---|---:|---:|---:|
| 0, 0.2, 0.5 (**current**) | 0.127 | 1.00 | 1.0× |
| 0, 0.2, 0.5, **1.0** | 0.567 | **0.47** | **4.5×** |
| 0, 0.25, 0.5, 1.0 | 0.547 | 0.48 | 4.3× |

**Adding a single dose row at 1.0 mM is worth about 4.5 times as many biological replicates.** At
n = 3 with the extended ladder, SE(slope) = 0.345, which beats n = 6 on the current ladder (0.516).
Against `replicates_needed` of 80 for `NativeYap1`, this is the difference between impossible and
routine — the requirement falls to roughly 18, and combining it with the free pooling across the two
Yap1 constructs (`G4_ANCHOR.md` §R2, λ = 0.12, near-complete pooling justified) brings it within
reach.

Caveat, and it is a real one: at 1 mM H₂O₂ the repo records growth at 23% of control and the
reporter response goes *negative* at 2 and 4 mM. So **1.0 mM is the ceiling, sampled early**, and
1.0 mM must be paired with a per-well growth-rate covariate — which the paired design supplies.

### 6.3 If a second oxidative target is wanted

Add `AHP1` (YLR109W) or `SRX1` (YKL086W) alongside `TRX2` rather than instead of it — both are
Yap1/Skn7-regulon members with lower basal expression than `TSA1`, so a larger fold is plausible.
⚠ I did not verify published induction magnitudes for either against primary text; treat this as a
hypothesis to test on the same plate, not as a recommendation with numbers behind it. The cost is
two plate rows, and the decision rule is simple: whichever of the three shows the larger
effect/spread ratio on the pilot becomes the anchor.

---

## 7. The confirmatory protocol

Hand this to a bench member.

### 7.1 What it must fix

Four things, in descending order of how much each is worth:

1. **The dose ladder** — worth ~20× the replicates (§7.3).
2. **The timepoint** — the current design samples the anchor after the UPR has switched itself off
   (§7.4). This is probably the second-largest error and nobody has flagged it.
3. **The normaliser** — three genes instead of one unstable one; ~2× the replicates (§5).
4. **The target and the gDNA** — `KAR2` instead of total `HAC1`, with DNase (§3, §4).

### 7.2 Assays

All primers verified in the S288C reference; each occurs exactly once across all 6039 ORFs.

| Role | Gene | Forward | Reverse | Amplicon | Source |
|---|---|---|---|---|---|
| **Anchor** | `KAR2` | `TCTGAAGGTGTCTGCCACAG` | `TTAGTGATGGTGATAGATTCGGATT` | 60 bp | Geronimo et al. 2025 |
| Pathway control | `HAC1i` | `ACCTGCCGTAGACAACAACA` | `ACCTGACTGCGCTTCTGGAT` | 145 bp | Hata et al. 2022 |
| Pathway control | `HAC1` total | `GCGTCGGACCAAGAGACTT` | `TCGTCGACTCTGGTACATTTTC` | 67 bp | Hata et al. 2022 |
| Reference 1 | `TAF10` | `ATATTCCAGGATCAGGTCTTCCGTAGC` | `GTAGTCTTCTCATTCTGTTGATGTTGTTGTTG` | 141 bp | Teste et al. 2009 |
| Reference 2 | `ALG9` | `CACGGATAGTGGCTTTGGTGAACAATTAC` | `TATGATTATCTGGCAGCAGGAAAGAACTTGGG` | 162 bp | Teste et al. 2009 |
| Reference 3 | `TFC1` | `GCTGGCACTCATATCTTATCGTTTCACAATGG` | `GAACCTGCTGTCAATACCGCCTGGAG` | 223 bp | Teste et al. 2009 |
| (optional) ER-specific | `SIL1` | `AGCCAGGCAATCTAACTTGG` | `TCAAGTGTCTGCTGTCGATAAGT` | 92 bp | Geronimo et al. 2025 |
| (optional) ER-specific | `JEM1` | `CCAAGATAACGGCCTCTCAG` | `GATGTCGTTGGTGTTGTTGC` | 158 bp | Geronimo et al. 2025 |
| (optional) HSR check | `HSP104` | `AAGGACGACGCTGCTAACAT` | `CACTTGGTTCAGCGACTTCA` | — | Geronimo et al. 2025 |

Report `HAC1i` as **splicing efficiency**, `2^(Cq_total − Cq_HAC1i)`, exactly as Hata et al. do —
self-normalising, needing no reference gene. Report `KAR2` as ΔΔCq against the geometric mean of the
three references.

### 7.3 Doses

**DTT: 0, 0.25, 0.5, 1.0, 2.0 mM.**

| DTT ladder | `Sxx` | SE(slope) relative | equivalent replicate multiplier |
|---|---:|---:|---:|
| 0, 0.2, 0.5 (**current**) | 0.127 | 1.00 | 1.0× |
| 0, 0.5, 1.0, 2.0 | 2.188 | 0.24 | **17.3×** |
| **0, 0.25, 0.5, 1.0, 2.0** | 2.500 | **0.23** | **19.7×** |

Why 2.0 mM is safe here even though the repo measures a lethal dose of 1.55 mM: that figure is for a
multi-hour plate run. At a 45-minute exposure the literature is unambiguous — Pincus et al. 2010
report that cells in 2.2 mM DTT "continue to grow" while only 5 mM arrests division. The dose fix and
the timing fix are the same fix.

Where this sits on the published dose–response also explains the original failure. In Pincus et al.'s
system, **DTT doses of 1.5 mM and below produce less than 10% of the maximal UPR reporter response**;
2.2 mM is about half-maximal, 3.3 mM about 75%, 5 mM near-saturating. The anchor was sampled at 0.2
and 0.5 mM — deep in the sub-10% region of the pathway's own transfer function. (These specific
percentages are from the model fitted to their measured titration; the measured Northern data behind
it show clear activation at 1.5 and 2.2 mM.)

### 7.4 Timepoint — the finding that is not in any existing note

**Sample at 45 min, not at the end of the plate read.** At low DTT the yeast UPR is a *pulse*, not a
step. Pincus et al. 2010 measured `HAC1` splicing by Northern blot across a DTT titration: at 5 mM,
splicing is maximal and sustained; **at 2.2 mM and 1.5 mM it activates and then deactivates within
4 h and 2 h respectively.** The same adaptation appears with tunicamycin — Geronimo et al. 2025 had
to raise the dose to 2.5 µg/mL "because HAC1 mRNA splicing subsides when cells are sustainedly
treated with 1.0 µg/mL tunicamycin".

The consequence for this experiment is severe and specific. mCitrine is a stable protein: the plate
reader records the **time-integral** of promoter activity. qPCR records the **instantaneous**
transcript pool. Sampling RNA hours after dosing at 0.2–0.5 mM catches the pathway *after it has
switched itself off*, while the reporter still carries the accumulated signal. A flat anchor beside a
risen reporter is then the *expected* result of a correctly functioning system — which is precisely
the failure mode `g4_anchor.py`'s docstring warns about, arriving by a route nobody had modelled.

Two corollaries:

- Harvest RNA at **45 min** (and, if budget allows, a second set at 90 min to see the decay).
- Compare `KAR2` not against accumulated RFU but against the reporter's **instantaneous promoter
  activity at 45 min**, which the repo can already compute from RFU/OD given `k_deg` (see
  `PROTOCOLS.md` P2 — this is a second, independent reason to finish that measurement).

A third caveat from the same literature: do not push DTT so high that translation fails. Geronimo et
al. 2025's central result is that when global protein synthesis is compromised, `HAC1` is spliced but
**UPR target genes are not induced** — cycloheximide added 10 min after DTT abolished `KAR2`, `ERO1`,
`SIL1` and `JEM1` induction despite successful splicing. Running `HAC1i` alongside `KAR2` detects
this directly: splicing up, `KAR2` flat means translational shutdown, not a broken reporter.

### 7.5 Replicates, plate design, and preserving the pairing

- **n = 4 biological replicates**, not 3. From `G4_ANCHOR.md` §R6: `t(0.975, 2) = 4.303` versus
  `t(0.975, 3) = 3.182` — **the fourth replicate shrinks every interval by 26% before a single extra
  data point is considered**, and it is worth more than the third or the fifth.
- **Run `UPRE1` and `UPRE2` together.** They share the anchor, the reference set, the stressor and
  the ladder, so pooling is licensed and free: it halves the SE, taking the minimum detectable
  effect from 2.77-fold to 1.56-fold (`G4_ANCHOR.md` §R2). If only one can be run, run **`UPRE2`** —
  it has the largest reporter response (2.18× at 2 mM) and it is the construct the 3-dose growth test
  would wrongly have REFUTED (growth R² 0.98 on 3 doses, 0.78 on 7).
- **Preserve and extend the pairing.** The existing protocol's real strength is that RNA came from
  the *same wells* as the fluorescence read. Keep that, and fix the timing conflict this way:
  designate **sacrifice wells** — one per dose per construct — and harvest them at 45 min while the
  rest of the plate keeps reading. By 45 min the reader has already taken 4–5 points on each
  sacrificed well, which is exactly what is needed to estimate its instantaneous promoter activity
  and growth rate at the harvest moment. The pairing is then at the level of the individual well,
  not the condition: `Var(paired difference) = 2σ²(1 − ρ)`, so at ρ = 0.5 the effective n doubles and
  at ρ = 0.7 it triples.
- **Break the replicate/run confound.** In the existing data, biological replicate is *perfectly
  confounded with qPCR run* (one culture set per plate). No model can separate them. Fix it in the
  design: put an aliquot of one pooled cDNA on **every** qPCR plate as an inter-run calibrator, and
  where possible split one biological replicate across two runs.

### 7.6 Controls

| Control | Purpose | Failure reading |
|---|---|---|
| **−RT, every sample × every gDNA-sensitive target** | contamination | margin < 3 cycles → re-DNase; < 0.74 → discard, do not correct |
| **NTC, every target** | primer-dimer | any amplification |
| **gDNA sensitivity grading** — each new pair run on a genomic DNA dilution series spanning ≥3 log₁₀ | ValidPrime Fig. 1A procedure; grades an assay `A+` if it does not amplify gDNA | `HAC1i` should be `A+`; nothing else will be |
| **−RT cross-gene check** | assay quality on DNA (§4.5) | a pair landing >1 cycle from `TAF10` on identical −RT template is a bad pair |
| **Efficiency**: 5-point, 4-fold dilution of pooled cDNA, per assay | `G4_STATISTICS.md` §4 requires one measured, assay-specific efficiency with its SE, propagated — never per-curve fitting | E outside ~0.9–1.05 → redesign |
| **Melt curve**, every well | specificity | secondary peaks |
| **Positive control**: 2 µg/mL tunicamycin, 60 min, one well per construct | the pathway's ceiling; validates every primer pair in one shot | no `KAR2` induction here means the assay, not the biology, failed |
| **`HAC1i` splicing** on every sample | did Ire1 fire? | splicing up + `KAR2` flat = translational shutdown (§7.4), not a dead reporter |
| **`HSP104`** on the top dose only | is the HSR firing, i.e. does `KAR2`'s HSE matter? | `HSP104` induced → switch the anchor to `SIL1` |
| **Reporter-free parent (BY4741)** wells | autofluorescence, already standard per `PROTOCOLS.md` P3 | — |
| **Junction-primer specificity titration** | Matabishi-Bibi et al.'s plasmid mixtures at u:s = 1:0 … 0:1 | non-linear response → `HAC1i` primer is cross-priming on `HAC1u` |

### 7.7 Analysis, and what counts as decisive

1. QC gate on `min(anchor pass rate, reference pass rate)`, not the anchor alone — `run_g4.py`
   currently gates only on the anchor, and the reference gene was half the failures
   (`G4_ANCHOR.md` §2.6).
2. Require **`n_doses ≥ 5`** before the correlation and growth tests run at all. At 3 doses Fisher's
   z has `SE = 1/√(n−3)` — division by zero — so `min_correlation = 0.7` has no confidence interval
   whatsoever, and `growth_confound_r2` can refute a construct on an R² with one degree of freedom.
   The 5-dose ladder in §7.3 is chosen partly for this.
3. Regress **anchor on reporter**, not the reverse, replicate-centred — the reporter is ~25× the
   more precise measurement, so it is the predictor (`G4_ANCHOR.md` §R3).
4. Report the equivalence test (TOST) against a **pre-registered** SESOI of 2-fold, fixed before
   looking at the data.

**Decisive means:** `KAR2` moves ≥2-fold from 0 to 2.0 mM DTT with SNR ≥ 2 across four replicates,
the `HAC1i` splicing fraction rises monotonically over the same ladder, the reporter tracks `KAR2`
with r ≥ 0.7 over five doses, and growth R² over those five doses is below 0.9. Any of those failing
is informative in a way the current INCONCLUSIVE is not.

---

## 8. Annotated bibliography

All entries resolved against PubMed/PMC on 2026-08-26. ⚠ marks something the reader should check.

### `HAC1` splicing and the UPR switch

- **Cox JS, Walter P.** "A novel mechanism for regulating activity of a transcription factor that
  controls the unfolded protein response." *Cell* 1996;87(3):391–404.
  doi:[10.1016/S0092-8674(00)81360-4](https://doi.org/10.1016/S0092-8674(00)81360-4). PMID 8898193.
  — `HAC1` identified; splicing regulates Hac1p level in an Ire1-dependent way.
- **Sidrauski C, Cox JS, Walter P.** "tRNA ligase is required for regulated mRNA splicing in the
  unfolded protein response." *Cell* 1996;87(3):405–413.
  doi:[10.1016/S0092-8674(00)81361-6](https://doi.org/10.1016/S0092-8674(00)81361-6). PMID 8898194.
  — Trl1/Rlg1 is the ligase. The companion paper to the above; note both are in the same issue and
  are frequently conflated in citation lists.
- **Sidrauski C, Walter P.** "The transmembrane kinase Ire1p is a site-specific endonuclease that
  initiates mRNA splicing in the unfolded protein response." *Cell* 1997;90(6):1031–1039.
  doi:[10.1016/S0092-8674(00)80369-4](https://doi.org/10.1016/S0092-8674(00)80369-4). PMID 9323131.
- **Kawahara T, Yanagi H, Yura T, Mori K.** "Endoplasmic reticulum stress-induced mRNA splicing
  permits synthesis of transcription factor Hac1p/Ern4p that activates the unfolded protein
  response." *Mol Biol Cell* 1997;8(10):1845–1862.
  doi:[10.1091/mbc.8.10.1845](https://doi.org/10.1091/mbc.8.10.1845). PMID 9348528. PMC25627.
  — **The source for "total `HAC1` is constant"**: the sum of the 1.4 kb and 1.2 kb species "remained
  almost constant", i.e. direct conversion, not induction. Also gives the 252-nt intron.
- **Kawahara T, Yanagi H, Yura T, Mori K.** "Unconventional splicing of `HAC1`/`ERN4` mRNA required
  for the unfolded protein response: sequence-specific and non-sequential cleavage of the splice
  sites." *J Biol Chem* 1998;273(3):1802–1807.
  doi:[10.1074/jbc.273.3.1802](https://doi.org/10.1074/jbc.273.3.1802). PMID 9430730.
- **Chapman RE, Walter P.** "Translational attenuation mediated by an mRNA intron." *Curr Biol*
  1997;7(11):850–859. doi:[10.1016/S0960-9822(06)00373-3](https://doi.org/10.1016/S0960-9822(06)00373-3).
  PMID 9382810. ⚠ **Bibliographic trap:** the DOI suffix `S0960-9822(06)00373-3` is a re-issued
  identifier, and several databases list an incorrect one (`…00038-8` belongs to a different *Current
  Biology* dispatch on the same topic). Check before citing.
- **Rüegsegger U, Leber JH, Walter P.** "Block of `HAC1` mRNA translation by long-range base pairing
  is released by cytoplasmic splicing upon induction of the unfolded protein response." *Cell*
  2001;107(1):103–114. doi:[10.1016/S0092-8674(01)00505-0](https://doi.org/10.1016/S0092-8674(01)00505-0).
  PMID 11595189. — Why `HAC1u` is stable and cytoplasmic rather than a transient pre-mRNA; the fact
  that breaks intron-flanking primer design.
- **Di Santo R, Aboulhouda S, Weinberg DE.** "The fail-safe mechanism of post-transcriptional
  silencing of unspliced `HAC1` mRNA." *eLife* 2016;5:e20069.
  doi:[10.7554/eLife.20069](https://doi.org/10.7554/eLife.20069). PMID 27692069. PMC5114014.
- **Leber JH, Bernales S, Walter P.** "IRE1-independent gain control of the unfolded protein
  response." *PLoS Biol* 2004;2(8):e235.
  doi:[10.1371/journal.pbio.0020235](https://doi.org/10.1371/journal.pbio.0020235). PMID 15314654.
  PMC509300. — The one paper that says `HAC1` abundance *can* rise (3–4×), and it is specific to
  ER-**distal** secretory stress in *sec* mutants, not to DTT or tunicamycin. Quote it precisely or
  it becomes an argument against §2.2 rather than for it.
- **Xia X.** "Translation control of `HAC1` by regulation of splicing in *Saccharomyces cerevisiae*."
  *Int J Mol Sci* 2019;20(12):2860. doi:[10.3390/ijms20122860](https://doi.org/10.3390/ijms20122860).
  PMC6627864. — Review; source for the 252-nt intron and independent-order cleavage.
- **Tehfe A, Roseshter T, Wei Y, Xia X.** "Does *Saccharomyces cerevisiae* require specific
  post-translational silencing against leaky translation of Hac1up?" *Microorganisms*
  2021;9(3):620. doi:[10.3390/microorganisms9030620](https://doi.org/10.3390/microorganisms9030620).
  PMC8002603. — Splicing efficiency 0.955–1.0 in UPR-induced cells versus essentially nil in
  unstressed; the source for "near-binary".

### Splice-specific qPCR assays (the crux)

- **Hata T, Ishiwata-Kimata Y, Kimata Y.** "Induction of the unfolded protein response at high
  temperature in *Saccharomyces cerevisiae*." *Int J Mol Sci* 2022;23(3):1669.
  doi:[10.3390/ijms23031669](https://doi.org/10.3390/ijms23031669). PMID 35163590. PMC8836091.
  — **The primary reference for a working `HAC1i` exon-junction real-time RT-PCR**, with exact
  primer sequences (Forward-1/Reverse-5 for `HAC1i`, Forward-2/Reverse-7 for total) and the
  self-normalising readout `2^(Cq_total − Cq_HAC1i)`. Also the 30→39 °C, ~2× splicing result.
- **Matabishi-Bibi L, Challal D, Barucco M, et al.** "Termination of the unfolded protein response is
  guided by ER stress-induced `HAC1` mRNA nuclear retention." *Nat Commun* 2022;13:6331.
  doi:[10.1038/s41467-022-34133-8](https://doi.org/10.1038/s41467-022-34133-8). PMID 36284099.
  PMC9596429. — A second, independent `HAC1u`/`HAC1i`/`HAC1Tot` qPCR scheme, and the **plasmid-ratio
  specificity titration** that should be copied.
- **Matabishi-Bibi L, Goncalves C, Babour A.** "RNA exosome-driven RNA processing instructs the
  duration of the unfolded protein response." *Nucleic Acids Research* 2025;gkaf088.
  doi:[10.1093/nar/gkaf088](https://doi.org/10.1093/nar/gkaf088). PMID 39995043. PMC11850225.
  — Same assays applied to a tunicamycin activation/deactivation time course.
- **Gast V, Campbell K, Picazo C, Engqvist M, Siewers V, Molin M.** "The yeast eIF2 kinase Gcn2
  facilitates H₂O₂-mediated feedback inhibition of both protein synthesis and endoplasmic reticulum
  oxidative folding during recombinant protein production." *Appl Environ Microbiol*
  2021;87(14):e0030121. doi:[10.1128/AEM.00301-21](https://doi.org/10.1128/AEM.00301-21). PMID
  34047633. PMC8276805. — A third `HAC1` spliced/unspliced primer set, plus `KAR2`/`PDI1`/`ERO1`/
  `EUG1`/`SCJ1`/`LHS1`/`JEM1` and `TRX2`/`TSA1`/`SRX1`/`CTT1` pairs in one table. Uses `ACT1` as
  reference — do not copy that part.

### UPR target genes and their regulation

- **Travers KJ, Patil CK, Wodicka L, Lockhart DJ, Weissman JS, Walter P.** "Functional and genomic
  analyses reveal an essential coordination between the unfolded protein response and ER-associated
  degradation." *Cell* 2000;101(3):249–258.
  doi:[10.1016/S0092-8674(00)80835-1](https://doi.org/10.1016/S0092-8674(00)80835-1). PMID 10847680.
  — The canonical genome-wide UPR target set. ⚠ **Full text not retrievable** (Cell blocks automated
  access, including the "bronze OA" PDF link Semantic Scholar advertises). Everything I state about
  its internals — the 60-min DTT/tunicamycin design, the *ire1Δ*/*hac1Δ* subtraction, and the
  seven-gene "canonical profile" — comes from secondary descriptions and is marked ⚠ in §3.
- **Kimata Y, Ishiwata-Kimata Y, Yamada S, Kohno K.** "Yeast unfolded protein response pathway
  regulates expression of genes for anti-oxidative stress and for cell surface proteins." *Genes
  Cells* 2006;11(1):59–69. doi:[10.1111/j.1365-2443.2005.00921.x](https://doi.org/10.1111/j.1365-2443.2005.00921.x).
  PMID 16371132. — 90 upregulated genes; UPR-dependence established by comparing constitutive
  `HAC1i` against *ire1Δ HAC1u* and confirming with tunicamycin. ⚠ Per-gene fold values are in
  Supplementary Table S1, which I could not retrieve.
- **Kohno K, Normington K, Sambrook J, Gething MJ, Mori K.** "The promoter region of the yeast
  `KAR2` (BiP) gene contains a regulatory domain that responds to the presence of unfolded proteins
  in the endoplasmic reticulum." *Mol Cell Biol* 1993;13(2):877–890.
  doi:[10.1128/mcb.13.2.877-890.1993](https://doi.org/10.1128/mcb.13.2.877-890.1993). PMID 8423809.
  PMC358971. — **The primary source for the repo's own note** that `KAR2` carries both a UPRE and an
  HSE: the two elements are "functionally independent of each other but work additively", and yeast
  BiP "is also induced severalfold by heat shock, albeit in a transient fashion".
- **Geronimo RAC, Ishiwata-Kimata Y, Funahashi Y, Izawa S, Kimata Y.** "Impairment in global protein
  synthesis uncouples UPR gene induction from `HAC1` mRNA splicing in *Saccharomyces cerevisiae*."
  *Front Microbiol* 2025;16:1629132. doi:[10.3389/fmicb.2025.1629132](https://doi.org/10.3389/fmicb.2025.1629132).
  PMID 40959222. PMC12435717. — The most useful single paper for this protocol. Supplies working
  primers for `KAR2`, `ERO1`, `SIL1`, `JEM1`, `HSP104` and `TAF10`; shows `SIL1`/`JEM1` are
  UPR-specific while `KAR2`/`ERO1` are HSR-inducible; shows `HSP104` is *not* induced by DTT,
  tunicamycin or inositol depletion; and demonstrates that splicing without translation gives no
  target-gene induction.
- **Ishiwata-Kimata Y, Nguyen PTM, Sugimoto M, Kimata Y.** "Potential of a constitutive-UPR and
  histone deacetylase A-deficient *Saccharomyces cerevisiae* strain for biomolecule production."
  *Appl Environ Microbiol* 2025;91(9):e0064425. doi:[10.1128/aem.00644-25](https://doi.org/10.1128/aem.00644-25).
  PMID 40772766. PMC12442391. — Current Hac1-induced gene list.
- **Sarkar D, Paira S, Das B.** "Nuclear mRNA degradation tunes the gain of the unfolded protein
  response in *Saccharomyces cerevisiae*." *Nucleic Acids Research* 2018;46(3):1139–1156.
  doi:[10.1093/nar/gkx1160](https://doi.org/10.1093/nar/gkx1160). PMID 29165698. PMC5814838.
  — ⚠ **The trap corrected in §3.** Its "2- to 2.5-fold" is a `rrp6Δ`/`cbc1Δ`-versus-wild-type ratio
  under tunicamycin, not a wild-type induction. Uses `SCR1` (a Pol III transcript) as loading
  control, which is not transferable to a poly(A)-primed RT.
- **Pincus D, Chevalier MW, Aragón T, van Anken E, Vidal SE, El-Samad H, Walter P.** "BiP binding to
  the ER-stress sensor Ire1 tunes the homeostatic behavior of the unfolded protein response." *PLoS
  Biol* 2010;8(7):e1000415. doi:[10.1371/journal.pbio.1000415](https://doi.org/10.1371/journal.pbio.1000415).
  PMID 20625545. PMC2897766. — **The DTT dose–response and the adaptation kinetics.** ≤1.5 mM gives
  <10% of maximal response; 2.2 mM ≈ half-maximal; 5 mM saturating and growth-arresting; splicing
  deactivates within 2 h at 1.5 mM and 4 h at 2.2 mM. The basis for §7.3 and §7.4.
- **Rubio C, Pincus D, Korennykh A, Schuck S, El-Samad H, Walter P.** "Homeostatic adaptation to
  endoplasmic reticulum stress depends on Ire1 kinase activity." *J Cell Biol* 2011;193(1):171–184.
  doi:[10.1083/jcb.201007077](https://doi.org/10.1083/jcb.201007077). PMID 21444684. PMC3082176.
  — Standard practice is 2 mM DTT for 40–60 min; useful corroboration for the dose window.

### Reference genes

- **Teste MA, Duquenne M, François JM, Parrou JL.** "Validation of reference genes for quantitative
  expression analysis by real-time RT-PCR in *Saccharomyces cerevisiae*." *BMC Molecular Biology*
  2009;10:99. doi:[10.1186/1471-2199-10-99](https://doi.org/10.1186/1471-2199-10-99). PMID 19874630.
  PMC2776018. — **The definitive yeast study.** `ALG9`, `TAF10`, `TFC1`, `UBC6` stable; `ACT1`,
  `PDA1`, `TDH3`, `IPP1`, `RDN18` explicitly invalidated. Table 1 carries the primer sequences and
  measured efficiencies reproduced in §5.3. ⚠ Note the volume/issue: it is **10:99**, not the
  "2009;10(1)" form some citation managers produce.
- **Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F.** "Accurate
  normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal
  control genes." *Genome Biology* 2002;3(7):RESEARCH0034.
  doi:[10.1186/gb-2002-3-7-research0034](https://doi.org/10.1186/gb-2002-3-7-research0034). PMID
  12184808. PMC126239. — geNorm; the origin of the multi-gene geometric-mean convention.
- **Vaudano E, Noti O, Costantini A, Garcia-Moruno E.** "Identification of reference genes suitable
  for normalization of RT-qPCR expression data in *Saccharomyces cerevisiae* during alcoholic
  fermentation." *Biotechnology Letters* 2011;33(8):1593–1599.
  doi:[10.1007/s10529-011-0603-y](https://doi.org/10.1007/s10529-011-0603-y). PMID 21452013.
  — Independent confirmation that a validated set beats `ACT1`. ⚠ Gene list verified only via
  abstract.

### Genomic DNA

- **Laurell H, Iacovoni JS, Abot A, Svec D, Maoret JJ, Arnal JF, Kubista M.** "Correction of RT-qPCR
  data for genomic DNA-derived signals with ValidPrime." *Nucleic Acids Research* 2012;40(7):e51.
  doi:[10.1093/nar/gkr1259](https://doi.org/10.1093/nar/gkr1259). PMID 22228834. **PMC3326333**.
  ⚠ `G4_STATISTICS.md` cites PMC3326325 — that identifier belongs to an unrelated ribosome paper.
  Correct it. Source for the 60% ceiling, the "5 cycles ⇒ <3%" reading, the exon–exon-junction
  gDNA-insensitivity criterion, and the 25%-of-RT-volume rule.
- **Đermić D, Ljubić S, Matulić M, Procino A, Feliciello MC, Ugarković Đ, Feliciello I.** "Reverse
  transcription-quantitative PCR (RT-qPCR) without the need for prior removal of DNA." *Scientific
  Reports* 2023;13(1):11470. doi:[10.1038/s41598-023-38383-4](https://doi.org/10.1038/s41598-023-38383-4).
  PMID 37454173. PMC10349872. — The mismatched-RT-primer method; four alternating 3′ mismatches for
  20–26 nt primers, discriminating anneal at ~60 °C. Demonstrated on intronless bacterial genes and
  on satellite DNA, i.e. exactly the hard case yeast presents.
- **Wiame I, Remy S, Swennen R, Sági L.** "Irreversible heat inactivation of DNase I without RNA
  degradation." *BioTechniques* 2000;29(2):252–256.
  doi:[10.2144/00292bm11](https://doi.org/10.2144/00292bm11). PMID 10948426. — EGTA before heat.
- **Bickler SW, Heinrich MC, Bagby GC.** "Magnesium-dependent thermostability of DNase I."
  *BioTechniques* 1992;13(1):64–66. PMID 1503777. — Why the previous entry is necessary.

### Oxidative arm

- **Kuge S, Jones N.** "`YAP1`-dependent activation of `TRX2` is essential for the response of
  *Saccharomyces cerevisiae* to oxidative stress by hydroperoxides." *EMBO J* 1994;13(3):655–664.
  doi:[10.1002/j.1460-2075.1994.tb06304.x](https://doi.org/10.1002/j.1460-2075.1994.tb06304.x).
  PMID 8313910. PMC394856. ⚠ PMC holds a scanned image only; no machine-readable text.
- **Lee J, Godon C, Lagniel G, Spector D, Garin J, Labarre J, Toledano MB.** "Yap1 and Skn7 control
  two specialized oxidative stress response regulons in yeast." *J Biol Chem* 1999;274(23):16040–16046.
  doi:[10.1074/jbc.274.23.16040](https://doi.org/10.1074/jbc.274.23.16040). PMID 10347154.
- **Godon C, Lagniel G, Lee J, et al.** "The H₂O₂ stimulon in *Saccharomyces cerevisiae*." *J Biol
  Chem* 1998;273(35):22480–22489. doi:[10.1074/jbc.273.35.22480](https://doi.org/10.1074/jbc.273.35.22480).
  PMID 9712873.
- **Tsuzi D, Maeta K, Takatsume Y, Izawa S, Inoue Y.** "Regulation of the yeast phospholipid
  hydroperoxide glutathione peroxidase `GPX2` by oxidative stress is mediated by Yap1 and Skn7."
  *FEBS Lett* 2004;565(1–3):148–154. doi:[10.1016/j.febslet.2004.03.091](https://doi.org/10.1016/j.febslet.2004.03.091).
  PMID 15135069. — and the companion showing calcineurin/Crz1 control:
  *FEBS Lett* 2004;569(1–3):301–306. doi:[10.1016/j.febslet.2004.05.077](https://doi.org/10.1016/j.febslet.2004.05.077).
  PMID 15225652. Together, the reason `GPX2` is rejected as a "clean" Yap1 anchor.
- **Vázquez J, González B, Sempere V, Mas A, Torija MJ, Beltran G.** "Melatonin reduces oxidative
  stress damage induced by hydrogen peroxide in *Saccharomyces cerevisiae*." *Front Microbiol*
  2017;8:1066. doi:[10.3389/fmicb.2017.01066](https://doi.org/10.3389/fmicb.2017.01066). PMC5471302.
  — `CTT1` up ~30× at 45 min after 2 mM H₂O₂ in early exponential phase, and ~1000× on entry into
  stationary phase. The second number is why `CTT1` is rejected: an anchor that moves 1000-fold with
  growth phase cannot arbitrate a growth confound.
- **Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D, Brown PO.**
  "Genomic expression programs in the response of yeast cells to environmental changes." *Mol Biol
  Cell* 2000;11(12):4241–4257. doi:[10.1091/mbc.11.12.4241](https://doi.org/10.1091/mbc.11.12.4241).
  PMID 11102521. PMC15070. ⚠ PMC holds abstract only; per-gene values require the supplement, whose
  original Stanford host is dead.
- **Brauer MJ, Huttenhower C, Airoldi EM, et al.** "Coordination of growth rate, cell cycle, stress
  response, and metabolic activity in yeast." *Mol Biol Cell* 2008;19(1):352–367.
  doi:[10.1091/mbc.e07-08-0779](https://doi.org/10.1091/mbc.e07-08-0779). PMID 17959824. PMC2174172.
  — Chemostat evidence that ESR expression is a function of growth rate; the formal basis for
  rejecting Msn2/4-driven genes as anchors in a gate that tests for a growth confound.

### Data sources used computationally

- **SGD** `SGD_features.tab`, chromosomal feature file, retrieved 2026-08-26 — intron counts
  (282/5783 verified ORFs = 4.9%) and per-gene intron status.
- **SGD** `orf_genomic.fasta` (6039 ORFs), retrieved 2026-08-26 — primer localisation, amplicon
  sizes, and uniqueness counts.
- **SGD** locus API, `YFL031W` sequence details — `HAC1` exon/intron coordinates used to verify every
  junction primer in §2.4.

---

## What this changes in the repo

Not made here (this note owns only itself), but implied:

- `qpcr.py::TARGET_FOR_CONSTRUCT` — `"Hac1"` → `"KAR2"` for `UPRE1`/`UPRE2`.
- `qpcr.py::REFERENCE_TARGET` — single string → sequence of three, normalising on the mean Cq
  (already specified in `G4_ANCHOR.md` §R4).
- `qpcr.py::QPCR_PLATE` — the recorded layout has 3 dose rows per gene half; the new ladder needs 5,
  and the plate must carry three reference genes plus `KAR2` plus two `HAC1` assays.
- `gates/g4_anchor.py` — require `n_doses ≥ 5` before the correlation and growth tests run.
- `generator/stress_panel.py` — the `tunicamycin` annotation "KAR2/PDI 2-2.5 fold, NAR 46:1139"
  misattributes a mutant-vs-wild-type ratio; see §3.
- `docs/research/G4_STATISTICS.md` — ValidPrime is PMC3326333, not PMC3326325.
