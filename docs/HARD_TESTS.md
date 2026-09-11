# Hard tests: what happened when the model was attacked

Three independent passes were pointed at the architecture with instructions to make it
produce a confidently wrong answer, to find the controls it would pass while picking up
nothing real, and to design the hardest held-out split the data supports.

**They found eight things. Six were real, and two of the six were false claims in this
repository's own documentation.** Everything below was measured, not argued, and each entry
says what it cost and what was done.

---

## 1. No pool had a physical bound — the model would say a cell is 36× its own weight

At an entry expression of 1e4, the lycopene pool returned **36,093 mg/gDCW — 3,609% of dry
weight** — as an ordinary float, with no note and no refusal.

The steady-state solve has exactly one positive root and no opinion about whether that root
describes anything. `solve_pathway` now raises `ImplausibleContent` above 1 g per gram of
dry weight. The bound is deliberately absurd: anything tighter would be a contested claim
about how much product a yeast can hold, and that does not belong hidden in a solver.

## 2. A ceiling nobody knew the model was asserting

**β-carotene content can never exceed 1.2483 mg/gDCW.** For a terminal node fed by a
saturating step, `content = flux_out/μ` and `flux_out ≤ vmax_per_growth · μ`, so
`content ≤ vmax_per_growth` and **μ cancels exactly**. No genotype at any growth rate.

| entry expression | content mg/gDCW |
| ---: | ---: |
| 1 | 1.1065 |
| 10 | 1.2369 |
| 100 | 1.2476 |
| **ceiling** | **1.2483** |

A hundredfold increase in cassette dosage buys 0.06%, and
`test_more_entry_expression_gives_more_product` passed throughout because it never went
high enough to notice.

It is **not** a fitting artefact — leave-one-out moves it between 1.11 and 1.41 mg/gDCW,
and one calibration state sits at 96.7% of it, so it is set by data. That makes it the
sharpest falsifiable claim in the package, and it is **17× below** López 2019's 21 mg/gDCW.

### THE CEILING IS REFUTED (2026-08-28)

**It was the sharpest falsifiable claim in the package, and it has been falsified — by a
factor of about 63.** Stated here rather than in a footnote, because a repository that
publishes a falsifiable claim and then buries its falsification has not been doing the thing
it says it does.

**Arhar et al. 2024, PMID 39215465, *J Appl Microbiol* 135(9):lxae224 — 79 mg/gDCW.** The
abstract states it outright: *"a strain with a β-carotene content of 79 mg g-1 cell dry
weight, corresponding to a 76-fold improvement over the starting strain."* It survives both
objections that dispose of the absorbance tier below:

- **Chromatographically resolved.** Waters e2695 HPLC, C18 reverse phase, β-carotene detected
  at 448 nm against a Sigma authentic standard, phytoene quantified separately at 286 nm
  against its own standard. A lycopene standard was run and *"biological samples never
  contained lycopene above the detection limit"* — so the one interferent that matters here
  is affirmatively excluded, not assumed away.
- **Gravimetric dry weight.** *"filtration of a defined culture volume through 0.45 um
  nitrocellulose filters, followed by overnight drying at 65C."* No OD conversion.
- **Internally consistent.** 672 mg/L ÷ 79 mg/gCDW = 8.5 g/L, a plausible shake-flask biomass.

**Fathi et al. 2021, PMID 33332529 — 46.5 mg/gDCW** is a second, independent confirmation:
different lab, different column chemistry (Discovery HS F5), different medium, β-carotene at
RT 7.6 min against a Sigma standard, DCW dried in a pre-weighed tube at 60 °C for 96 h. It
reports no lycopene at all, so it does not itself *demonstrate* separation from lycopene —
Arhar alone carries the refutation and Fathi corroborates it.

**The universal ceiling failed; a dosage-dependent replacement remains unidentified.**
`MEASUREMENTS_NEEDED.md` M5 called for *"a capacity that reads the cyclase's own dosage"*.
The cyclase fit uses a shared capacity across Elizondo's six states; it does not establish a
crtYB-dosage scaling law. `Genotype.cassette` now represents per-gene dosage, but
`capacity_for` refuses to scale the anchor at another crtYB dosage. The survey reader
`pathway/published_cassettes.py` finds seven rows with numeric crtYB dosages but still only
two distinct values, 1 and 2, across different study contexts. The 2026-09-09 literature pass
added two of those seven and moved neither the distinct count nor the refusal. These comparisons do not isolate
a dosage or localisation effect.

**What survives and what does not:**

| | Status |
| --- | --- |
| The fit to Elizondo's six states (contents 0.5–1.0 mg/gDCW) | **survives** — inside its calibration range |
| `content <= 1.2483 mg/gDCW` as a claim about *any* genotype at *any* growth rate | **dead** |
| `NodeKinetics.vmax_per_growth` as a *universal constant* | **dead as a universal bound** — `Genotype.cassette` carries dosage, but no validated dosage-to-capacity scaling is implemented |
| The μ-independence of the ceiling algebra | **untested** — no paper in the set reports a growth rate, so nothing here bears on it |

The missing representation now exists in `Genotype.cassette`. `capacity_for` retains the
anchor at the calibration dosage and refuses other crtYB dosages; changing the constant
would not identify the missing dosage-to-capacity relation.

### The López row, separately: this paragraph was right

**That gap is now closed, and this paragraph was right.** PMID 31380362's methods quantify
"the absorbance at **453 nm** of the hexane extracts" against a β-carotene standard, with no
chromatography anywhere, so the 21 mg/gDCW is a summed absorbance of every coloured
carotenoid — β-carotene **and lycopene** — not β-carotene. `EXTERNAL_CAROTENOID_BOUND.md`
had recorded all four López rows as analyte `beta-car`; that was the error, and it is
corrected there. The distinction decides the comparison rather than softening it: on the
model's own D = 0.18 prediction, an A453 assay reads **1.42× and 1.98× above the ceiling**
for β-car3 and β-car4 while β-carotene itself stays below it, so an absorbance number above
1.2483 is the *expected* reading and not evidence against the bound. Predictions within 10% of it now say
so, because two strains differing by 0.4% at an asymptote read as a ranking when both are
just reporting the fitted capacity.

## 3. The FBA audit was never in the chain, and this repository said it was

```
grep -n '^from' src/ystwin/predict.py   →   no ystwin.fba
```

A commit message said "the FBA layer is now called" and `WHAT_IS_LEFT.md` said "reached
from `predict.py`". `fba/audit.py` is reached only from `scripts/predict_product.py
--audit`. Opt-in is defensible — a GSMM costs an 18-second SBML parse — but describing it
as part of the chain was not. The docstring now reads `NOT CALLED` and says why.

## 4. The growth-rate guard was an accident of one node's rate law

The range check lived inside the SATURATING and PROPORTIONAL arms of the walk, so the
terminal and passthrough arms never reached it. `phb` and `glycogen`, both passthrough end
to end, had **no growth-rate guard at all**. It now runs once over every node before
anything is solved, so a refusal cannot arrive after a half-solve.

## 5. Carbon source was accepted and discarded in silence

Two calls differing only in `carbon_source` returned bit-identical numbers **and identical
notes**, so nothing distinguished *"the model says the carbon source does not matter"* from
*"the model cannot see the carbon source"*.

The number is still identical — inventing a term is what
[`EXTERNAL_PRODUCT_VALIDATION.md`](EXTERNAL_PRODUCT_VALIDATION.md) §3 shows the data cannot
support — but every prediction now names what it discarded, and carries the 4.24× the feed
moves the flux by on Kocharin's series.

## 6. The physiology stamped on every prediction describes a different culture

`chemostat_physiology` reads van Hoek 1998's **glucose-limited** chemostat. Elizondo's six
calibration states are **glucose-excess**:

| state | μ | q_glucose model | measured | ethanol model | measured |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2D01 | 0.101 | 1.11 | 6.68 | 0.00 | 6.40 |
| 3D025 | 0.254 | 2.89 | 11.72 | 0.02 | 15.84 |

4.1–6.4× on uptake, and every state ferments at every rate including μ = 0.101 where a
glucose-limited culture makes none. So `ProductPrediction.fermentative` reports `False` for
every calibration state while those cultures were visibly fermenting. The flag is kept — it
is right about van Hoek's cultures and two consumers depend on it — and every prediction now
says which culture it describes.

## 7. The extrapolation window was two literals unrelated to the data

`flux.py` hard-coded `0.2 <= expression <= 5.0` and the note described it as "roughly
0.5–3×". The measured CrtE range is **[0.245, 1.000]** — a 4.08-fold span whose maximum sits
**5× below the literal ceiling**, so an expression of 2.0, twice anything ever measured,
passed silently. The window is now the calibration's own recorded range, in both directions,
with the magnitude of the excursion in the note.

## 8. Phytoene's `passthrough` is now known-wrong rather than unverified

`beta_carotene.toml` declared phytoene as holding no pool and said the assumption was
untested. It has since been tested by someone else: Chen et al. 2016 (PMID 27329233)
quantify phytoene at **3.99%** and neurosporene at **4.87%** of total carotenoid, and state
that CrtI conversion is rate-limiting.

It stays `passthrough`, because changing it to `saturating` needs a CrtI vmax and Km that do
not exist for these strains — substituting invented kinetics for a declared-wrong assumption
trades a visible error for a hidden one. The note now carries the refutation. What the flux
law predicts is the **post-CrtI** flux, and any claim about total pathway flux is low by an
unmeasured amount.

---

## What the hard tests did not find

No arithmetic error. The generic solver reproduces the hand-written one bitwise, carbon
closes to 1e-19, and the steady-state identity holds in Kocharin's own numbers to 1e-6. The
leave-one-strain-out score of 14.2% survived every attack on it.

**The failures were all at the boundaries** — what the model does outside its calibration,
what it refuses, and what it claims about itself. That is where a model that is honest inside
its range is most likely to be dishonest, and it is where none of the existing 2,400 tests
were looking.
