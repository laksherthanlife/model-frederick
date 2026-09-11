# Can the existing G4 replicates be rescued? No, and here is the arithmetic

## Verdict

**The ER anchor is unsalvageable analytically. Do not attempt a gDNA correction on it.** At UPRE1's observed +0.13-cycle margin, genomic DNA contributes **91.4%** of the +RT signal and transcript only 8.6%. The only published gDNA-correction method is validated up to **60% gDNA** — a margin of about 0.74 cycles — so the ER data sits far outside the correctable regime. At literature-typical Cq noise the correction returns a *negative* transcript quantity about a third of the time. UPRE2, at +0.60 cycles and 66% gDNA, is also past the ceiling, though only just.

**The oxidative anchor is clean and merely underpowered, and that is a different and much better problem.** TRX2's 8.6-cycle margin is 0.3% gDNA. What it needs is either more replicates or a larger effect size — and a larger effect size is worth far more, because the requirement scales with the square of the inverse SNR.

One genuinely useful correction to received wisdom: **MIQE never mandated a 5-cycle RT− margin.** That number is a software default whose real content is "≤3% background", and it originated as a *no-template-control* rule. MIQE 2009, 2010 and 2.0 (2025) all require only that the +RT/−RT comparison be *reported* and that the experimenter *state their own* tolerance. MIQE 2.0 goes further and explicitly endorses correcting rather than discarding — which is why the arithmetic below matters rather than the convention.

---

## 1. The right unit is %gDNA, not cycles

`qpcr.py::rt_minus_margin` currently gates on `_MIN_RT_MINUS_MARGIN = 3.0` cycles. Three cycles is a more defensible threshold than the community's five, and the docstring's reasoning is right. But cycles are a log-scale proxy for the quantity that actually decides whether a measurement means anything:

```
%gDNA of the +RT signal  =  2^(−ΔCq) × 100,   ΔCq = Cq(−RT) − Cq(+RT)
```

Everything below is computed from that identity. I verified each row independently.

| ΔCq margin | %gDNA | %cDNA | variance inflation | P(negative estimate) at σ=0.2 | status |
|---:|---:|---:|---:|---:|---|
| 8.60 | 0.3 | 99.7 | 1.00× | 0.0% | **TRX2 observed — clean** |
| 5.00 | 3.1 | 96.9 | 1.03× | 0.0% | the community convention |
| 3.00 | 12.5 | 87.5 | 1.15× | 0.0% | this repo's current gate |
| 2.00 | 25.0 | 75.0 | 1.37× | 0.0% | |
| 1.00 | 50.0 | 50.0 | 2.24× | 0.0% | |
| **0.74** | **59.9** | 40.1 | 2.90× | 0.4% | **ValidPrime correctability ceiling** |
| **0.60** | **66.0** | 34.0 | 3.52× | 1.7% | **UPRE2 observed — past the ceiling** |
| **0.13** | **91.4** | **8.6** | **15.72×** | **32.3%** | **UPRE1 observed — hopeless** |

Variance inflation is the delta-method factor for the subtraction relative to the uncorrected quantity, `√(1+r²)/(1−r)` with `r` the gDNA fraction — note it is **independent of σ**, so it cannot be beaten by better pipetting. The final column is the probability that the observed margin goes to zero or below and the subtraction returns a negative transcript quantity, taking SD(observed ΔCq) = σ√2.

The precision requirement runs the same way. For a target of ≤20% relative SD on the corrected quantity:

- at σ = 0.1 cycles you need ΔCq ≥ **0.75**
- at σ = 0.2 cycles you need ΔCq ≥ **1.84**
- at σ = 0.3 cycles it is **unachievable at any margin** — the `ln2·σ` noise floor is already 20.8%

Published technical SD of Cq for well-behaved replicates is **0.1–0.3 cycles** (Nordgård et al. 2006, Table 1: SD 0.21–0.26 across four transcripts; Bustin et al. 2025 assume a fixed ΔCq SD of 0.2 for CI modelling), rising to 0.5–2 cycles above Cq ≈ 30–35 where Poisson sampling dominates. So the σ = 0.2 row is the realistic one, and it demands a margin fourteen times larger than UPRE1's.

For scale: qbase+'s own flagging arithmetic notes that **0.14 cycles of difference between duplicates corresponds to about 0.1 SD units.** UPRE1's entire margin is at the resolution limit of the software convention it would be judged against.

## 2. The subtraction exists, was written down, and was deliberately rejected

The formula is not novel. Laurell et al. (*Nucleic Acids Research* 2012;40(7):e51 — the ValidPrime paper) write it as their Equation (3):

```
Cq_RNA = −log₂( 2^(−Cq(RT+)) − 2^(−Cq(RT−)) )
```

and then decline to build on it:

> "Traditionally, determination of the RNA component using RT(−) controls would be achieved using Equation (3). **However… low reproducibility and other factors detract from the accuracy of this approach.**"

and, in the Discussion, the sentence that settles the question:

> "It is possible that **the lack of accuracy and low reproducibility generally observed in RT(−) reactions has previously restrained the development of a correction-based model similar to that proposed in Equation (3).**"

Their reasons are specific and all apply here: the gene-of-interest assay was designed to amplify transcript, not genomic DNA, so its efficiency *on gDNA* is unknown; and low template copy number makes RT− replicates dominated by stochastic effects.

The method they built instead estimates the gDNA contribution from a separate non-transcribed-locus assay plus a gDNA standard, and its published validity limit is explicit and repeated three times in the paper: correction works **"as long as the DNA contribution to the total signal is <60%"**, above which their software emits `HIGHDNA` and *"correction is not recommended"*.

**No publication derives the variance, bias or identifiability limit of the RT− subtraction.** ValidPrime gives only the empirical 60% ceiling. The variance-inflation and sign-flip columns in §1 are therefore my own derivations from the published σ values, using the delta method of Nordgård et al. and Hellemans et al. — presented as derived, not quoted.

Two related notes worth carrying. Software in this space **discards rather than subtracts**: qbase+ auto-*excludes* data points within its negative-control threshold and contains no subtraction formula anywhere. And where subtraction has been tried in an adjacent setting — Li et al. (*BMC Genomics* 2022;23:554) subtracting gDNA-derived FPKM in RNA-seq — it **backfired when the background estimate was poor**, increasing the differential-expression count rather than reducing it.

## 3. What to do with the oxidative anchor instead

TRX2 passes RT− QC on 100% of readings with an 8.6-cycle margin. Its problem is spread between biological replicates: NativeYap1 at SNR 0.39 needs roughly 80 replicates, AlteredYap1 at SNR 1.19 needs about 9.

`gates/g4_anchor.py` derives that count as `n_reps × (target_snr / snr)²` — correct, and the quadratic is the whole story. **Doubling the effect size is worth four times as many replicates.** So the productive move on the oxidative side is not 80 plates; it is a better-behaved anchor gene, or a dose closer to the responsive window. The dose-response analysis already identified that window as roughly 0.5–1 mM, and the anchor RNA was taken at 0, 0.2 and 0.5 mM — so two of three anchor points sit at or below the bottom of the responsive range, which is a plausible and cheap explanation for the small effect. Sampling the anchor at 1.0 mM as well may recover more SNR than any number of extra replicates.

## 4. Three statistical improvements that are worth making regardless

**Express the gate in %gDNA.** The physically meaningful quantity is the contamination fraction; cycles are a log proxy. Reporting both, and tying the threshold to the ValidPrime ceiling and to the achievable-precision requirement, replaces a convention with a derivation. This follows the repo's own rule about using the authoritative signal rather than the convenient one.

**Distinguish REFUTED from INCONCLUSIVE statistically, not just verbally.** G4 already draws the three-way distinction, and the statistical form of it is an equivalence test — two one-sided tests against a pre-named margin, as now implemented in `analysis/uncertainty.py::equivalence`. "The anchor effect is positively small" and "we could not measure the anchor effect" are different findings and the gate should be able to say which it has.

**Bootstrap, never the delta method, for anything near a boundary.** Bilgrau et al. (*BMC Bioinformatics* 2016;17:159) benchmark delta method against Monte Carlo and bootstrap for efficiency-adjusted ΔΔCq and note that division by a quantity near zero "can increase the variance dramatic[ally]". At a marginal ΔCq the delta method will hand back a confidently narrow, symmetric interval on a quantity that is a third likely to be negative. Their other finding is worth acting on independently: the common practice of treating amplification efficiency as known **systematically underestimates the standard error** and inflates false positives.

One efficiency note, since the repo's ddCq assumes a value: Ruijter et al. (*Clinical Chemistry* 2021;67(6):829–842) show that with efficiencies in the plausible range 1.93–2.05, a reported ΔCq consistent with an 8-fold ratio actually admits ratios from **1.45 to 42.6**. But do *not* fix this by fitting efficiency per amplification curve — Nordgård et al. found that "sample-specific amplification efficiencies determined from individual amplification curves primarily increase the random error… and should be avoided", with CVs up to 250%, because "errors in the base number (amplification efficiency) are more severe than errors in the exponent (threshold cycle)". Use one measured, assay-specific, dilution-series efficiency with its standard error, and propagate that error.

## 5. What this means for the four verdicts

| Construct | Anchor | Current | After this analysis |
|---|---|---|---|
| UPRE1 | HAC1 | INCONCLUSIVE | **REFUTED as a measurement** — 91.4% gDNA; the assay did not measure transcript. Not a statement about the biosensor. |
| UPRE2 | HAC1 | INCONCLUSIVE | **REFUTED as a measurement** — 66% gDNA, past the correctable ceiling. |
| NativeYap1 | TRX2 | INCONCLUSIVE | INCONCLUSIVE, genuinely — clean anchor, underpowered. Needs a bigger effect, not 80 plates. |
| AlteredYap1 | TRX2 | INCONCLUSIVE | INCONCLUSIVE, genuinely — closest to decidable; ~9 replicates, or fewer with a better dose. |

The distinction in the right-hand column matters for how the team talks about this. For the ER pair, **nothing has been learned about the biosensors** — the anchor measurement failed, and that is a protocol finding, not a biosensor finding. For the oxidative pair the measurement worked and the answer is simply not yet resolved. Only the second pair is evidence of anything about the constructs.

The replacement anchor for the ER pair is the subject of `UPR_ANCHOR.md`.

---

## Bibliography

**Guidelines**
- Bustin et al. "The MIQE Guidelines." *Clinical Chemistry* 2009;55(4):611–622. doi:10.1373/clinchem.2008.112797. PMID 19246619 — requires reporting the +RT/−RT comparison; **mandates no numeric margin**
- Bustin et al. "MIQE précis." *BMC Molecular Biology* 2010;11:74. doi:10.1186/1471-2199-11-74
- Bustin et al. "MIQE 2.0." *Clinical Chemistry* 2025;71(6):634–651. doi:10.1093/clinchem/hvaf043. PMID 40272429 — explicitly endorses gDNA correction, citing ValidPrime
- dMIQE Group. "The Digital MIQE Guidelines Update." *Clinical Chemistry* 2020;66(8):1012–1029. doi:10.1093/clinchem/hvaa125

**gDNA correction**
- Laurell et al. "Correction of RT–qPCR data for genomic DNA-derived signals with ValidPrime." *Nucleic Acids Research* 2012;40(7):e51. doi:10.1093/nar/gkr1259. PMC3326333 — **the primary reference; the 60% ceiling**. An earlier version of this document gave PMC3326325, which is an unrelated ribosomal-protein paper in the same journal and year; adjacent identifiers in one issue are an easy way to be wrong while looking right, and the DOI is the safer handle.
- Li, Zhang, Wang & Yu. "Genes expressed at low levels raise false discovery rates in RNA samples contaminated with genomic DNA." *BMC Genomics* 2022;23:554. doi:10.1186/s12864-022-08785-1
- Đermić et al. "Reverse transcription-quantitative PCR without the need for prior removal of DNA." *Sci Rep* 2023;13:11470. doi:10.1038/s41598-023-38383-4
- Padhi et al. "A PCR-based approach to assess genomic DNA contamination in RNA." *Anal Biochem* 2016;494:49–51. doi:10.1016/j.ab.2015.10.012

**Error propagation**
- Hellemans, Mortier, De Paepe, Speleman & Vandesompele. "qBase relative quantification framework." *Genome Biology* 2007;8(2):R19. doi:10.1186/gb-2007-8-2-r19 — delta method; formulas 5–6, 12, 14, 16. Contains **no** −RT term
- Nordgård, Kvaløy, Farmen & Heikkilä. "Error propagation in relative real-time RT-PCR quantification models." *Analytical Biochemistry* 2006;356(2):182–193. doi:10.1016/j.ab.2006.06.020. PMID 16899212
- Bilgrau et al. "Unaccounted uncertainty from qPCR efficiency estimates entails uncontrolled false positive rates." *BMC Bioinformatics* 2016;17:159. doi:10.1186/s12859-016-0997-6
- Bustin, Kirvell, Nolan, Mueller & Shipley. "When Two-Fold Is Not Enough: Quantifying Uncertainty in Low-Copy qPCR." *IJMS* 2025;26(16):7796. doi:10.3390/ijms26167796
- Ruiz-Villalba, Ruijter & van den Hoff. *Life* 2021;11:496 — technical replicate tolerances

**Efficiency**
- Ruijter et al. "Efficiency Correction Is Required for Accurate Quantitative PCR Analysis and Reporting." *Clinical Chemistry* 2021;67(6):829–842. doi:10.1093/clinchem/hvab052 — ⚠️ quotes obtained via the publisher page, PDF not independently parsed
- Ruijter et al. "Amplification efficiency: linking baseline and bias." *Nucleic Acids Research* 2009;37(6):e45. doi:10.1093/nar/gkp045
- Nolan, Hands & Bustin. "Quantification of mRNA using real-time RT-PCR." *Nature Protocols* 2006;1(3):1559–1582. doi:10.1038/nprot.2006.236 — the 5-cycle rule appears here, for **NTC**, not −RT

**Unverified** — no primary Bio-Rad or Thermo Fisher technical note stating a numeric −RT margin could be retrieved; the "≥5–7 cycles" line circulates only in secondary sources. Do not attribute it to a vendor.
