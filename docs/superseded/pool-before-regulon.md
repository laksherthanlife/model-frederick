# Superseded: pools respond before the regulons they share an agent with

**Current result:** there is no general pool-before-regulon ordering, and the reason is
sharper than "the rule is wrong." **H2O2's two pool arms go opposite ways.** Cytosolic
H2O2 leads its regulon threefold; the glutathione pool lags it more than tenfold. One
agent, one insult, two orderings — which is exactly what a single multiplier cannot
encode.

The old rule got `peroxide` roughly right by coincidence and `redox` wrong by more than an
order of magnitude, and neither number came from a measurement of a pool.

---

## What was asserted

`generator/stress_panel.py` carried a blanket constant:

```python
_POOL_POTENCY = 0.35
_POOLS = frozenset({"redox", "peroxide", "atp", "ph", "nadh"})
```

Any pool arm with no explicit override got `agent_ec50 x 0.35`, justified in
`module_ec50`'s docstring as *"Pools lead because their chemistry has no threshold to
cross"* and in the test module as *"peroxide oxidises glutathione stoichiometrically, with
no threshold to cross."* Two tests asserted it:

```python
def test_a_pool_responds_below_the_regulon_it_shares_an_agent_with(self):
    """H2O2 oxidises the glutathione pool before Yap1 has fired."""
    assert module_ec50("H2O2", "peroxide") < module_ec50("H2O2", "oxidative")
```

## Why it is wrong

**The two arms are different observables and were treated as one.** In this panel
`peroxide` is *cytosolic H2O2*, read by HyPer7; `redox` is the *glutathione pool*, read by
roGFP2-Grx1. They share a `_POOLS` membership and nothing else. H2O2 reaches the first by
diffusion, with no enzyme in between, and the second only through a peroxidase.

**So peroxide genuinely leads.** Kritsiligkou 2021 (**PMID 34118234**), in *S. cerevisiae*
BY4742 with this module's own sensor: *"The lowest concentration triggering a detectable
response was ≈ 5 μM for roGFP2-Tsa2ΔC**R** and ≈ 20 μM for HyPer7."* Against 100 µM for the
first detectable oxidised Yap1. The test asserting `peroxide < oxidative` was **right**,
for a module its own docstring was not describing.

**And glutathione genuinely lags. The chemistry there is wrong, not approximate.** H2O2 does not oxidise glutathione directly at
any relevant rate. The rate constant for H2O2 with the glutathione *thiolate* is
18–26 M⁻¹s⁻¹, and only about 2% of GSH is thiolate at cytosolic pH, so the effective
constant is below 1 M⁻¹s⁻¹ (Winterbourn & Metodiewa 1999, **PMID 10468205**). The same
H2O2 reacts with the yeast peroxiredoxin Tsa1 at ~10⁷ M⁻¹s⁻¹ (Ogusucu 2007,
**PMID 17210445**). That is six to seven orders of magnitude, and it means GSSG formation
in vivo is enzyme-catalysed and downstream — Calabrese 2019 (**PMID 31389622**) shows the
mitochondrial E_GSH response to H2O2 is abolished in a Prx1 active-site mutant. No
peroxidase, no glutathione oxidation.

**And Yap1 is wired to the fast sensor on purpose.** Yap1 is not oxidised by H2O2 at all;
it is relayed through the thiol peroxidase Gpx3/Orp1, which is in Tsa1's kinetic class
(Delaunay 2002, **PMID 12437921**). The relay is catalytic and amplifying. So the
transcription factor fires *before* the bulk pools move — that is the point of the
architecture, and it is why the antioxidant response is useful at all.

**The measured ordering, opposite to the claim:**

| observable | dose | delivery | citation |
| --- | ---: | --- | --- |
| cytosolic H2O2 (HyPer7) — lowest detectable response | **20 µM** | bolus | Kritsiligkou 2021, PMID 34118234 |
| cytosolic H2O2 (roGFP2-Tsa2ΔC<sup>R</sup>) — lowest detectable | 5 µM | bolus | Kritsiligkou 2021 |
| oxidised Yap1 first detectable (25 and 50 µM below detection) | 100 µM | bolus | Delaunay 2000, PMID 11013218 |
| Yap1-GFP partial nuclear entry, no growth effect | 0.1 mM | perfusion | Goulev 2017, PMID 28418333 |
| Yap1-GFP saturated | >0.2 mM | perfusion | Goulev 2017 |
| cytosolic E_GSH — **no shift at all** | 0.2 mM | bolus, 60 min | Ayer 2013, PMID 23762325 |
| cytosolic rxYFP ~30% oxidised | 0.4 mM | bolus | Dardalhon 2012, PMID 22561702 |
| cytosolic E_GSH first significant shift, 40–50 mV | 1 mM | bolus | Ayer 2013 |
| cytosolic roGFP2 **half** oxidised | 2 mM | bolus | Ayer 2013 |

So the glutathione pool lags the regulon by roughly 13× half-maximal to half-maximal
(2 mM against 0.15 mM), while cytosolic H2O2 leads it by 3×.

Morgan 2013 (**PMID 23242256**) is the strongest single source, and it is the authors' own
section heading: *"Supplementary Figure 1 … **The cytosolic glutathione pool is robustly
resistant to perturbation.**"* Their H2O2 series starts at **0.5 mM** — they did not test
lower — and their working figure for the total pool is 10 mM (Deponte 2017,
**PMID 28540740**, puts the estimate at 13 mM, and titles his review *"The Incomplete
Glutathione Puzzle: Just Guessing at Numbers and Figures?"*, which is the right level of
confidence for it).

**One over-reach corrected in the course of this.** An earlier draft of this note cited
Delaunay 2002 for "Yap1 target-gene induction at 50 µM." That is wrong: the 50–800 µM
range in that paper establishes that *Gpx3 is required* for normal TRX2 induction across
it, measured in Δgpx3 cells. It is not a wild-type activation threshold. The defensible
floor is Delaunay 2000's 100 µM, by direct redox western, with 25 and 50 µM tested and
below detection.

**A second error, in the test rather than the panel.** The test docstring says *"H2O2
oxidises the glutathione pool before Yap1 has fired"* but asserts on `peroxide` — and in
this panel `peroxide` is **cytosolic H2O2** read by HyPer7, while `redox` is the
glutathione pool read by roGFP2-Grx1. The test justified itself with one module and
measured another.

**A third, in the H2O2 spec.** Its source note said the pool arms were *"pinned at the
value they already had"* because *"Goulev measured Yap1 nuclear entry, not glutathione
chemistry, so the correction has no claim on them."* But the value they already had was
`0.5 × 0.35 = 0.175` — the product of the very EC50 being corrected and a default that had
the ordering backwards. Pinning preserved the product of two errors and presented it as an
independent datum. The argument that immunised 0.175 against the correction is the same
argument that would have prevented 0.175 from existing.

The pin also silently inverted the invariant: with the agent EC50 at 0.15 and the pool at
0.175, the pool led nothing. The two tests had been failing for that reason.

## How far it spread

Eight of twenty-five agents have pool arms — twelve arms in total. Auditing each against
the literature:

| verdict | arms |
| --- | --- |
| **refuted** | H2O2 → redox (lags 13×, not leads), menadione → redox/peroxide |
| **right by coincidence** | H2O2 → peroxide — 0.5 × 0.35 = 52 µM sits between the measured 20 µM floor and the 150 µM regulon EC50, but nothing in the derivation knew that |
| **unestablished** | diamide → redox, DTT → redox, antimycin_A → nadh, menadione → nadh |
| **defensible** | acetic_acid → ph, sodium_hydroxide → ph, glucose_starvation → atp, antimycin_A → atp |

Three of eight agents survive on the blanket rule, one of them by accident. For menadione the nearest evidence
argues actively against a lead: 2 mM paraquat, the same redox-cycling class, produces no
change in cytosolic E_GSH while still inducing SOD2 (Ayer 2013), and menadione activates
Yap1 partly by direct cysteine modification, a route that never touches the pool
(Azevedo 2003, **PMID 14556853**). For diamide the mechanism *does* favour a lead — it
oxidises glutathione directly — and an earlier draft of this note inferred from the 8–20 mM
**saturating calibration** doses that the pool must therefore lag. **That inference was
withdrawn**: a saturating calibration point says nothing about a threshold. Morgan 2013
measured cytosolic E_GSH under **1.5 mM** diamide, exactly this panel's EC50, so diamide
does engage the pool at its own EC50. No dose-response, EC50 or Hill for cytosolic E_GSH
under diamide in *S. cerevisiae* exists that could be found, so the arm is marked unsourced
and no direction is asserted.

## What replaced it

`_POOL_POTENCY` is gone. In its place:

- `_DIRECT_POOLS = {"atp", "ph"}` with `_DIRECT_POOL_LEAD = 0.35`. A weak acid entering
  the cytosol or a blocked respiratory chain draining ATP needs no enzyme, so it moves at
  the concentration of the insult. The mechanism is sound and the factor is asserted.
- `redox`, `peroxide` and `nadh` get **no default**. `module_ec50` raises on an arm with
  no value in `target_ec50`, so a future pool arm cannot silently inherit an ordering that
  the same agent's other arm contradicts.
- H2O2 is measured, in both directions: `peroxide` **0.05 mM**, a bracket between
  Kritsiligkou's 20 µM detection floor and the 150 µM regulon EC50 — labelled a bracket
  because no HyPer7 half-max exists for yeast; `redox` **2.0 mM**, which is a true
  half-activation (Ayer: 2 mM "led to oxidation of half of the roGFP2 probe").
- The unestablished arms assert **no separation** — pool EC50 equal to the agent's own —
  with the reason recorded in each source note. Asserting no ordering is the honest
  position when the sign is not known; inventing a factor is not.

## Two consequences worth keeping

**The glutathione arm of the H2O2 panel is unreachable.** `redox` at 2.0 mM is *twice*
this agent's measured lethal dose of 1.0 mM, so `healthy_ladder` — capped at 0.55 mM —
never comes near it. roGFP2-Grx1 will read flat across every well the panel can afford to
spend. That is a fact about the design, not a number to shrink until it looks usable, and
it is the same trap the ESR arm fell into at 1.25 mM before that was caught.

It also shows up downstream, which is the useful part: the three-sensor build recovers
three of H2O2's four target modules instead of four, and the one it loses is `redox`.
`test_and_h2o2_loses_its_glutathione_arm_to_its_own_lethal_dose` asserts that, and asserts
that menadione — whose pool arms sit at its own EC50 — keeps its `redox` arm. The panel
cannot report what it cannot dose for.

**The dose axis is unflattened in both directions now.** `peroxide` at 0.05, `oxidative`
at 0.15, `ESR` at 0.375 and `redox` at 2.0 are four separated rungs on one agent, each
with its own source, where before they were one number and a multiplier.

## The lesson

A default that fills in every gap will be inherited by cases nobody checked. The blanket
factor read as a modelling convenience and functioned as an unsourced empirical claim about
eight agents — and it survived because nothing forced any single arm to name its evidence.
Refusing is better than defaulting where the sign of the effect is genuinely unknown.

Both tests were failing before this was investigated, in the permissive direction, and the
first instinct was that the constant was stale and needed re-deriving as `0.15 × 0.35`.
That would have restored a false ordering with more precision.
