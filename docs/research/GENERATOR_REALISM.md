# How close the generator is to a wet-lab plate

Lane C audit of `src/ystwin/generator/`. What it models, what it does not, what changed, and
the `od_linear_max` verdict.

The owner wants to sweep many environmental conditions, so the generator's realism is what
every simulated result rests on. That framing sets the priority order below: a gap is serious
in proportion to how many downstream numbers move when it is wrong, not to how visibly wrong
it looks.

---

## 1. The prioritised gap list

Twelve gaps, most-serious first. "Fixed" means fixed in this pass; "reported" means it is
outside this lane's files and is written down here instead.

| # | Gap | Was | Now |
| --- | --- | --- | --- |
| 1 | **Growth-rate-dependent physiology and the Crabtree switch** | Absent. Growth rate only diluted the reporter. Biomass yield was a per-strain constant, so a slowed culture ended up *sparser*. | **Fixed.** `culture.chemostat_physiology` reads van Hoek 1998 Table 1; `substrate_limited_capacity` puts the measured yield into the carrying capacity. §2. |
| 2 | **OD non-linearity at high density** | Absent. `generator/plate.py` computes `true_od = biomass / gdcw_per_od + od_blank`, exactly linear at every density, even though `observation.observe_od` carries a `saturation_k` term that would model it. | **Reported.** Not this lane's file. The consequence is stated in §4, and it is the same defect as the `od_linear_max` placeholder seen from the generator's side. |
| 3 | **Reader and well noise** | Three different numbers in three places, only one of them measured. | **Reported and used.** §3. |
| 4 | **Temperature** | `1 - 0.04 * abs(T - 30)`. Symmetric and linear; yeast is neither. | **Fixed.** Rosso's cardinal temperature model at Salvadó's measured cardinals. §5. |
| 5 | **Oxygen limitation** | `1 - respiring * (1 - O2/0.21)`, linear, and it put an unsupplemented anaerobic glucose culture at 0.28 /h. | **Endpoint fixed, middle still asserted.** §6. |
| 6 | **Reporter maturation** | `reporter.py` has modelled it since it was written and `culture.py` never used it. Every simulated well assumed a reporter fluorescent the instant it is translated. | **Fixed.** `CultureParameters.k_mat`, default `None` so nothing moves. §7. |
| 7 | **pH** | Not on the context at all. Not on growth, not on the reporter — even though `redox.py` is an entire module about pH confounding a sensor. | **Recorded; reporter effect added; growth effect refused.** §8. |
| 8 | **Autofluorescence** | 900 RFU/gDCW, constant, never measured on any plate. | **Reported.** §9. |
| 9 | **Carbon source and diauxie** | A three-row static table. `context.py`'s own docstring says a batch glucose well ferments, exhausts, then respires its own ethanol — and no code does that. | **Reported.** §9. |
| 10 | **Plate position** | `plate_layout.edge_multiplier` exists; `panel_experiment` takes `edge_effect` and defaults it to 0.0; `generator/plate.py` does not model position at all. | **Reported.** §9. |
| 11 | **Blank handling** | One fixed `od_blank = 0.098` for the whole plate and the whole run. | **Reported.** §9. |
| 12 | **pH → growth** | — | **Deliberately refused.** §8. |

---

## 2. Growth physiology, grounded in the measured chemostat curve

### The problem in one sentence

A stressor in this generator did exactly one thing to physiology — it lowered the growth rate —
and in a real aerobic glucose culture lowering the growth rate is not one thing, it is a switch.

### What the data says

`data/physiology/chemostatData_VanHoek1998.tsv`, from van Hoek, van Dijken & Pronk 1998,
**PMID 9797269**, "Effect of specific growth rate on fermentative capacity of baker's yeast",
Table 1. Verified: the title above is the string PubMed returns for that identifier, and
`docs/research/CHEMOSTAT_REFERENCE.md` records the arithmetic check of the transcription
against the paper's own carbon balance. Strain DS28911, aerobic glucose-limited chemostat,
30 °C, dissolved oxygen above 60 % of air saturation. Ten dilution rates, 0.025 to 0.40 /h.

Three facts, and the generator contradicted all three:

- **Ethanol is undetectable up to D = 0.25 /h, appears at 0.28 and reaches 13.9 mmol/gDW/h at
  0.40.** There is a threshold in the growth rate, not a slope.
- **Biomass yield is flat at 0.45–0.49 g/g below the threshold and collapses to 0.20 above it.**
  The *faster* culture makes less biomass per gram of sugar — a factor of 2.4.
- **Oxygen uptake is not monotonic.** It peaks at 7.4 at the threshold and falls to 3.7 at
  0.40 /h. Code that assumes oxygen demand rises with growth rate is wrong on the four rows
  that matter.

### Why this is the biggest gap and not a detail

The generator's untreated glucose culture grows at **0.40 /h**, which is the top of van Hoek's
range and well above the 0.28 /h threshold. A mid-range stress dose takes it below. So the
switch sits in the middle of every dose ladder the project runs — it is not a corner case
reached only by an exotic sweep.

### What changed

`src/ystwin/generator/culture.py` gains:

- `CRITICAL_GROWTH_RATE_PER_H = 0.28` — named because it is a threshold *in the growth rate*,
  and the growth rate is what every stressor moves.
- `chemostat_physiology(growth_rate)` → `ChemostatPhysiology`, carrying glucose uptake, O₂,
  CO₂, ethanol, glycerol, acetate, pyruvate and biomass yield, interpolated between measured
  rows. Also `fermentative` (the switch as a boolean) and `respiratory_quotient`.
- `substrate_limited_capacity(initial_biomass, glucose_g_per_L, growth_rate)` — the carrying
  capacity as `inoculum + glucose × Y(µ)` rather than an asserted per-strain constant.

The numbers are **read from the TSV, not transcribed into source**. A second copy would be one
edit away from disagreeing with the first, and the checked copy is the file.

### It refuses rather than extrapolating

`chemostat_physiology` raises outside [0.025, 0.40] /h. The tempting bad fix is to clamp to the
nearest measured row, and it is silent: a stressed culture in this generator goes to zero
growth, so the low edge is hit constantly. Clamping would report a 0.45 g/g biomass yield for a
poisoned well on the authority of a chemostat held at 0.025 /h — a completely different
physiological state that happens to share a number with nothing.

The sweep catches that refusal and writes it into a column rather than stopping. **On the
default sweep, 46 % of rows are refused.** That is a result: it says how much of the
environment space this generator can currently ground.

### Is the Crabtree transition reproduced? Yes

One construct, one environment (exponential glucose, 30 °C, air, full nutrients), dose ladder,
from `scripts/run_environment_sweep.py`:

| dose (mM) | µ (1/h) | fermentative | Y (g/g) | q_ethanol | q_O₂ | capacity (g/L) |
| ---: | ---: | :---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.400 | yes | 0.200 | 13.90 | 3.70 | 4.06 |
| 0.10 | 0.397 | yes | 0.202 | 13.65 | 3.78 | 4.10 |
| 0.25 | 0.388 | yes | 0.207 | 12.85 | 4.03 | 4.21 |
| 0.50 | 0.366 | yes | 0.221 | 10.89 | 4.66 | 4.48 |
| 1.00 | 0.312 | yes | 0.337 | 3.99 | 5.87 | 6.81 |
| 2.00 | 0.215 | **no** | 0.480 | 0.00 | 5.82 | 9.67 |

The dose ladder walks the culture across the switch. Ethanol goes to zero, the yield doubles,
oxygen uptake **rises** as growth falls, and the biomass the well can reach goes up 2.4×.

**One honest caveat, and it is load-bearing.** On the 4.14 h read the real plates use, the
final OD still *falls* with dose, because 4 h never approaches the carrying capacity. The
yield effect is invisible on a short read and inverts the density ordering on a run to
substrate exhaustion. So this fix changes nothing about the existing 4 h plates and changes a
great deal about any sweep that runs to saturation — which is exactly the kind of sweep the
generator is now being asked to produce.

**A second caveat about regime.** van Hoek is glucose-*limited*; a microplate well is
glucose-*excess*. In glucose excess a Crabtree-positive yeast ferments at any growth rate
because glucose repression is on, so the true batch yield at 0.40 /h is if anything *below*
0.20 g/g. Applying this curve to a batch well is therefore a regime extrapolation in the
conservative direction on yield, and the *threshold* it locates is a chemostat property. This
is stated rather than corrected because no equivalent batch curve was found, and inventing one
would be worse than naming the limitation.

**Tier.** The curve itself is a measurement, faithfully read. Its application to a batch
microplate well is **Tier 0 / asserted** in the sense of `docs/CLAIM_BOUNDARY.md`: an adopted
model, not a fitted or validated one.

---

## 3. Where the measured well CV lives, and where it does not

The measured number is `analysis/power.py::DEFAULT_WELL_CV = 0.052` — the lognormal sigma of
recovered activity across the zero-dose wells of the two NewProtocol plates that carry a usable
blank, median over eight construct-by-plate groups, range 0.006 to 0.083. It is genuinely
measured and its docstring says how.

**The generator does not use it.** Three numbers, three places:

| where | constant | value | status |
| --- | --- | ---: | --- |
| `analysis/power.py` | `DEFAULT_WELL_CV` | 0.052 | **measured** on the real plates |
| `generator/plate.py` | `PlateConditions.well_cv` | 0.040 | asserted; 23 % below the measured value |
| `generator/plate.py` | `PlateConditions.reader_cv` | 0.010 | asserted |
| `generator/panel_experiment.py` | `panel_dataset(noise_cv=...)` | 0.050 | asserted default, against a measured total plate-to-plate CV of 0.146 in the same file |

`generator/plate.py` is not this lane's file, so the constant is not changed there. What this
lane did instead: **`scripts/run_environment_sweep.py` defaults `--well-cv` to
`analysis.power.DEFAULT_WELL_CV`**, so there is one definition and the sweep inherits the
measured value rather than a second copy of an asserted one.

The recommended fix for another lane is a one-line import in `plate.py`, not a new number.
`generator/calibrate.py` already estimates both CVs from a real plate, so the honest end state
is that `PlateConditions` defaults come from a calibration rather than from a literal.

---

## 4. `od_linear_max` — the verdict

**Keep the placeholder. The literature bounds the answer and cannot set it.**

### What the literature does establish

Stevenson, McVey, Clark, Swain & Pilizota 2016, *Scientific Reports* 6:38828,
**PMID 27958314**, doi 10.1038/srep38828, **"General calibration of microbial growth in
microplate readers"**. Verified: that title is what PubMed returns for 27958314.

This is the right paper for the question and not a generic one. It measured **yeast** — S288C
derivatives BY4741, BY4743 and a `cln3` deletion — at 30 °C, in a **Costar flat-bottom 96-well
plate with lid, 200 µL per well**, read at **600 nm**. That is the strain, format, fill volume
and channel this repository uses.

Two quotes carry the answer:

> "In the single scattering regime (small N and OD600 ≲ 0.2) … where the Beer-Lambert law
> (OD ~ N) approximately holds, exact solutions to the scattering problem exist"

> "In this regime the Beer-Lambert law is no longer a suitable approximation, and OD is
> expected to have a parabolic dependency on N"

So: **OD600 ≈ 0.2 is where proportionality between optical density and cell number ends.**

### Why that number still cannot be `od_linear_max`

Three reasons, and each on its own is sufficient.

**Units.** `OpticalQualityGate` compares `max_raw_od` — the raw plate-reader number with the
blank still in it. Stevenson's ≲ 0.2 is in path-corrected, cuvette-equivalent units: their
methods state that main-text values are "the measured OD multiplied by 1.0560 for 300 µl,
1.5848 for 200 µl and 6.3694 for 100 µl". Two things follow. At 200 µL the raw threshold would
be ≈ 0.13. And the 100 µL factor is **four times** the 200 µL one, not two times, so the factor
is not a path length and cannot be derived from well geometry — it is an instrument-and-plate
measurement. No literature number converts to this repository's reader and fill volume.

**The scatterer is not fixed.** Stevenson's central result is that the OD↔concentration
calibration depends on the size and shape of the cell, demonstrated by comparing haploid yeast,
diploid yeast and a large-cell `cln3` mutant. Stress changes cell size. So the "linear range"
is not a property of the reader alone, and a single number gates a moving target.

**Two different non-linearities are being conflated.** Multiple scattering is a property of the
*sample*; detector saturation and stray light are properties of the *instrument*. The literature
figure is about the first. `calib/od_linearity.py`'s dilution series measures their combination
in the medium and format actually used, which is the quantity the gate actually wants.

### What happens if the literature number is taken at face value

The real plates read raw OD from 0.09 (blank) to 0.69. A 200 µL raw threshold of ≈ 0.13 is
below where most of these wells *start* once the blank is included, and far below where they
end. **On the strict Beer–Lambert reading, essentially the whole of every plate in this
repository is in the multiple-scattering regime.**

That has a direction, and stating it is the point of doing this rather than shrugging. In the
parabolic regime OD ≈ aN + bN² with b > 0, so `d(ln OD)/dt < d(ln N)/dt`: **the growth rate
estimated from OD is biased low, and increasingly so as the culture gets denser.** A fast
control reaches high density sooner than a slowed treated well, so the control is under-read
more, the apparent growth *difference* between them shrinks, and the dilution correction
`k = dR/dt + µR` then under-corrects the control relative to the treated well. That is a bias
with dose structure running in the same direction as the confound this project already tracks.

### What the current default does

`OpticalQualityGate.od_linear_max` defaults to **1.0**, which is the cuvette figure. Nothing in
the literature supports 1.0 for a 96-well raw reading, and the real plates top out at 0.69, so
**the linearity check never fires.** The gate is currently inert on this axis. That is the
honest description of the placeholder and it is worse than "unmeasured": it reads as a check.
`gates/g1_optical.py` is not this lane's file and was not changed.

### What measurement would settle it

Protocol P1/P4's dilution series, and `calib/od_linearity.py::fit_od_linearity` already turns
it into the number. Stevenson adds two requirements P1 does not currently state:

1. **The cells must come from the same growth state as the wells being gated**, not just the
   same medium — because the calibration moves with cell size, and P1 as written grows a
   culture to OD 2.5–3 and dilutes it, which is a *stationary* scatterer used to gate
   exponential wells.
2. **Report the path-correction factor for the exact fill volume**, since without it the number
   is not comparable to any published figure, including Stevenson's own.

A bead standard would make the result comparable across instruments: Beal et al. 2020,
*Communications Biology* 3:512, **PMID 32943734**, doi 10.1038/s42003-020-01127-5, **"Robust
estimation of bacterial cell count from optical density"** — the iGEM Interlab study, which
exists precisely because plate-reader OD is not comparable between labs without one.

---

## 5. Temperature

**Was:** `rate *= max(0.0, 1.0 - 0.04 * abs(T - 30))`. Symmetric about 30 °C and linear.

**Now:** the cardinal temperature model with inflection of Rosso, Lobry & Flandrois 1993,
**PMID 8412234**, doi 10.1006/jtbi.1993.1099, "An unexpected correlation between cardinal
temperatures of microbial growth highlighted by a new model", evaluated at cardinals measured
for *S. cerevisiae* and **normalised to 1.0 at 30 °C**.

Cardinals from Salvadó, Arroyo-López, Guillamón, Salazar, Querol & Barrio 2011,
*Applied and Environmental Microbiology* 77:2292, **PMID 21317255**, doi 10.1128/AEM.01861-10,
"Temperature adaptation markedly determines evolution within the genus Saccharomyces", Table 2.
Ten *S. cerevisiae* strains from ten independent origins, YNB + 20 g/L glucose, **96-well
microplates read at OD600**, eleven temperatures from 4 to 46 °C — the same medium, format and
channel this generator simulates.

`T_opt = 32.3` and `T_max = 45.4` are the species means the paper states in its own abstract.
`T_min = 2.8` is the mean of the same ten Table 2 rows; the abstract states no species mean for
that column. Averaging the other two columns the same way reproduces 32.27 and 45.40 exactly,
which is what licenses reading `T_min` off the table by the same arithmetic.

The normalisation matters. The carbon-source growth rates in `context.py` were read from
30 °C cultures, so the temperature term has to be a *shape* that leaves 30 °C alone. Using the
model's own µ_opt would double-count the rate and silently rescale every existing result.

What it changes:

| T (°C) | old linear | measured CTMI |
| ---: | ---: | ---: |
| 20 | 0.600 | 0.609 |
| 25 | 0.800 | 0.847 |
| **30** | **1.000** | **1.000** |
| 32.3 | 0.908 | 1.020 |
| **37** | **0.720** | **0.916** |
| 40 | 0.600 | 0.716 |
| **45** | **0.400** | **0.072** |
| 46 | 0.360 | 0.000 |

The two disagree most exactly where a heat-shock experiment is run. At 37 °C — the standard
heat-stress condition — the old form was 21 % too harsh. At 45 °C it was **5.6× too generous**,
because it had no idea a maximum growth temperature existed; the measured maximum is 45.4 °C.

---

## 6. Oxygen

**The endpoint is now a measurement. The middle is still asserted, and says so.**

**Was:** `rate *= 1 - respiring * (1 - O2/0.21)`, linear in the gas fraction. At O₂ = 0 on
glucose that gives `0.40 × (1 − 0.3) = 0.28 /h`. That number is wrong for a reason no amount of
tuning fixes: **sterol and unsaturated fatty acid synthesis both consume molecular oxygen**, so
a defined medium without them supports no anaerobic growth at all — Andreasen & Stier 1953,
*J Cell Comp Physiol* 41:23, **PMID 13034889**, doi 10.1002/jcp.1030410103, "Anaerobic nutrition
of Saccharomyces cerevisiae. I. Ergosterol requirement for growth in a defined medium".

**Now:**

- `CultureContext.anaerobic_supplements` (default `False`) records whether the medium carries
  ergosterol and an unsaturated fatty acid.
- At O₂ = 0: zero on a non-fermentable substrate; zero on a fermentable one without supplements;
  otherwise capped at `ANAEROBIC_MU_MAX_PER_H = 0.31`.
- Between air and anoxia: **unchanged and explicitly asserted.**

0.31 /h is Verduyn, Postma, Scheffers & van Dijken 1990, *J Gen Microbiol* 136:395,
**PMID 1975265**, doi 10.1099/00221287-136-3-395, "Physiology of Saccharomyces cerevisiae in
anaerobic glucose-limited chemostat cultures": CBS 8066 in mineral medium supplemented with
ergosterol and Tween 80, µ_max 0.31 /h, maximal biomass yield 0.10 g/g at D = 0.10. It is
applied as a **ceiling** rather than a value, because the measurement is glucose and the carbon
table also holds galactose, for which no anaerobic rate was found.

That paper carries a second fact the generator still does not represent: **the yield collapses
long before the rate does.** 0.31 against 0.42 /h is a 26 % loss of rate; 0.10 g/g against van
Hoek's 0.48 g/g at the same D = 0.10 is a **4.8-fold** loss of yield. Any model that represents
oxygen limitation as a growth-rate multiplier alone gets the small effect right and the large
one wrong. Closing that would mean an anaerobic branch in `substrate_limited_capacity`, and it
is left open because a two-point yield curve is not a curve.

**Why the middle was not replaced.** The linear shortfall is very probably too gradual — yeast's
affinity for oxygen is high enough that respiration runs almost to anoxia, so the true curve
should stay flat and then fall off a cliff. No verified half-saturation constant for
*S. cerevisiae* oxygen uptake was found in this pass, and swapping one invented shape for
another buys nothing. The two regimes therefore **do not meet**, and the discontinuity at O₂ = 0
is left visible: the asserted curve carries an unsupplemented anaerobic glucose culture to
0.28 /h and the answer is zero. A measured half-saturation constant would close it properly.

---

## 7. Reporter maturation

`reporter.py` has modelled chromophore maturation since it was written — `ReporterKinetics`
takes `k_mat`, `simulate_reporter` integrates a two-species system, and `promoter_activity`
inverts it. **`culture.py` never passed it.** Every simulated well in this repository has
therefore assumed a reporter that is fluorescent the instant it is translated.

`CultureParameters.k_mat` now exists and is passed through. It defaults to `None`
(instantaneous), so **no existing result moves** — the point is that the assumption is now
reachable and its cost is measurable.

No default value is supplied, and that is deliberate. Nagai, Ibata, Park, Kubota, Mikoshiba &
Miyawaki 2002, *Nat Biotechnol* 20:87, **PMID 11753368**, doi 10.1038/nbt0102-87, "A variant of
yellow fluorescent protein with fast and efficient maturation for cell-biological applications",
identifies chromophore oxidation as the rate-limiting step of YFP maturation and shows a single
F46L substitution "greatly accelerat[ing]" it. The rate belongs to the exact variant, so there
is no such thing as "the YFP maturation rate" to borrow. It has never been measured for the
mCitrine constructs on these plates, and that is a protocol-shaped hole: it matters most on the
short reads this project runs, because an unmatured reporter is invisible and a 4.14 h read
gives it little time to arrive.

---

## 8. pH

### Recorded

`CultureContext` gains `ph_medium` and `ph_cytosolic`. They are different numbers and the
distinction is not pedantry — yeast holds a near-neutral cytosol against a medium at 4, and
passing a medium pH to a reporter model predicts a dark cell that is not dark.

### The effect that is measured is on the reporter, not on growth

`citrine_ph_response(ph)` implements a single-protonation titration at Citrine's measured pKa.
Griesbeck, Baird, Campbell, Zacharias & Tsien 2001, *J Biol Chem* 276:29188, **PMID 11387331**,
doi 10.1074/jbc.M102815200, "Reducing the environmental sensitivity of yellow fluorescent
protein. Mechanism and applications": the Q69M substitution "confers a much lower pK(a) (5.7)
than for previous YFPs".

The pKa is measured; the single-site Henderson–Hasselbalch *form* is **asserted** — Griesbeck
reports a pKa, not a Hill coefficient, so the steepness of the curve is not pinned by that
paper.

The size of the effect is why this belongs in the generator at all. A full unit of cytosolic
acidification from 7.0 to 6.0 removes **30 %** of the fluorescence with the promoter untouched.
The induction folds this project measures are around 1.5. So the artifact is a third of the
signal, and it has dose structure, because the stressors that acidify are the ones being dosed.
This is the same argument `redox.py` already makes for the glutathione sensor, applied to the
reporter channel.

### The effect on growth is refused

`ph_growth_factor(ph, cardinal_ph)` implements the cardinal pH model of Rosso, Lobry, Bajard &
Flandrois 1995, *Appl Environ Microbiol* 61:610, **PMID 16534932**, doi
10.1128/aem.61.2.610-616.1995, "Convenient Model To Describe the Combined Effects of Temperature
and pH on Microbial Growth" — and **raises unless the caller supplies measured cardinal pH
values.**

Three reasons, and they compound:

- **No cardinal pH values for *S. cerevisiae* were found in this pass.**
- **A relevant negative result argues against guessing them.** Arroyo-López, Orlić, Querol &
  Barrio 2009, *Int J Food Microbiol* 131:120, **PMID 19246112**, doi
  10.1016/j.ijfoodmicro.2009.01.035, fit response surfaces to *S. cerevisiae* T73,
  *S. kudriavzevii* and their hybrid; pH was significant for the hybrid and for
  *S. kudriavzevii* and **not** for *S. cerevisiae*, where temperature dominated.
- **A bare pH knob models the wrong variable.** Below about pH 5 what damages yeast is the
  *undissociated weak acid* that low pH creates. Verduyn 1990 measured exactly this: acetate or
  propionate added to anaerobic chemostats lowered the biomass yield and raised specific ethanol
  production through proton import, an effect "independent of the dilution rate at a given acid
  concentration". Medium pH and weak-acid concentration are not separable knobs.

Returning 1.0 instead of raising would quietly assert that pH does not matter, in every sweep
that varies it. The refusal is the honest answer and the test suite pins it.

One structural thing the same Rosso 1995 paper does license: it builds its combined model "using
the hypothesis that the temperature and pH effects on the µ_max are independent". That is the
citation for multiplying the environment terms together rather than fitting interactions, which
is what `context_growth_rate` does.

---

## 9. Gaps left open, and what each would cost

**Autofluorescence.** `ReporterOptics.autofluorescence = 900` RFU/gDCW against a background of
400 — not a small term, and `docs/CLAIM_BOUNDARY.md` records it as never measured on any plate
(both reporter-free BY4741 controls were read in OD600 only). It is also modelled as *constant
per gram*, which is a second assumption: yeast autofluorescence rises as a culture leaves
exponential growth, so the constant is least right in exactly the contexts this sweep adds.
Fixing it needs a plate, not a paper: read the existing reporter-free controls in the mCitrine
channel.

**OD non-linearity in the generator.** §4 establishes that the real plates sit in the
multiple-scattering regime. `generator/plate.py` nevertheless produces perfectly linear OD, so a
pipeline validated on synthetic plates has never met the bias it will meet on real ones — and
the bias has dose structure. `observation.observe_od` already carries the `saturation_k` term
that would model it. This is the highest-value remaining fix and it is one line in another
lane's file, gated on a measured `saturation_k` from protocol P1.

**Carbon source and diauxie.** `_CARBON` is three static rows. `context.py`'s own docstring
states the problem: "A well started on glucose ferments it, then switches to the ethanol it
made, and arrives at the strongest derepressing state on the plate — which is why a late reading
of a glucose well can exceed a galactose one." No code does that. `_PHASE` labels the shift
rather than producing it. With the yield curve now in place the pieces exist to make glucose
exhaustion a *consequence* — glucose charge, uptake rate and yield are all present — but the
second growth phase on ethanol would need its own rate, and this lane did not build it.

**Plate position.** `plate_layout.edge_multiplier` models the outer ring, `panel_experiment`
exposes it as `edge_effect` and defaults it to 0.0, and `generator/plate.py` ignores position
entirely. Its own docstring makes the argument for why this cannot be folded into noise:
replicates in the same ring inherit the same bias, so averaging does nothing, and a generator
with noise but no geometry overstates what replication buys. The default of 0.0 is the part
that matters — the machinery exists and is off.

**Blank handling.** One `od_blank = 0.098` for the whole plate and the whole run, with reader
noise on top. Real blanks drift upward with evaporation and differ by position, and the same
edge wells that evaporate most are the ones whose blank is furthest off. Since every
blank-corrected OD divides into a growth rate, a drifting blank is a systematic error in µ that
grows through the run.

---

## 10. Every constant added or changed, and where it comes from

| Constant / function | Where | Value | Provenance |
| --- | --- | --- | --- |
| `CRITICAL_GROWTH_RATE_PER_H` | `culture.py` | 0.28 /h | **Measured.** van Hoek 1998, PMID 9797269, Table 1 |
| `chemostat_physiology` curve | `culture.py` | 10 rows × 8 fluxes | **Measured.** Same table, read from the vendored TSV |
| `substrate_limited_capacity` | `culture.py` | `X₀ + S·Y(µ)` | **Measured** yield; the mass-balance form is arithmetic |
| `CultureParameters.k_mat` | `culture.py` | `None` | **No default supplied.** Nagai 2002, PMID 11753368, for why one cannot be borrowed |
| `CARDINAL_TEMPERATURES_C` | `context.py` | 2.8 / 32.3 / 45.4 °C | **Measured.** Salvadó 2011, PMID 21317255, Table 2, *S. cerevisiae* means |
| CTMI functional form | `context.py` | — | **Adopted model.** Rosso 1993, PMID 8412234 |
| normalisation at 30 °C | `context.py` | ÷ CTMI(30) | **Arithmetic**, to keep the existing carbon-source rates intact |
| `ANAEROBIC_MU_MAX_PER_H` | `context.py` | 0.31 /h | **Measured.** Verduyn 1990, PMID 1975265. Applied as a ceiling |
| unsupplemented anoxic growth = 0 | `context.py` | 0.0 | **Measured requirement.** Andreasen & Stier 1953, PMID 13034889 |
| sub-air oxygen shortfall | `context.py` | linear | **Asserted**, unchanged, and now marked as such |
| `CITRINE_PKA` | `context.py` | 5.7 | **Measured.** Griesbeck 2001, PMID 11387331 |
| single-site titration form | `context.py` | — | **Asserted shape**, measured pKa |
| `ph_growth_factor` | `context.py` | refuses | **Refusal.** Rosso 1995, PMID 16534932, is the model; no *S. cerevisiae* cardinals found; Arroyo-López 2009, PMID 19246112, found pH non-significant for *S. cerevisiae* |
| multiplicative environment terms | `context.py` | — | **Adopted hypothesis.** Rosso 1995 independence assumption |
| dose as a *fraction* of the context rate | `design.py` | — | **Asserted.** Consistent with `context.py`'s stated position; nothing measures whether a stressor's fractional inhibition is temperature-independent |
| sweep `--well-cv` default | `run_environment_sweep.py` | 0.052 | **Measured**, imported from `analysis.power.DEFAULT_WELL_CV` rather than copied |
| sweep optics, `GDCW_PER_OD`, `OD_BLANK` | `run_environment_sweep.py` | 2.6e6 / 400 / 900, 0.42, 0.098 | **Asserted**, carried over from `generator/plate.py` and marked in the module |

Nothing was tuned toward a nicer-looking result. The one place a threshold was chosen rather
than derived is a test tolerance in `tests/test_environment_sweep.py`, and it is stated there
with the measured collinearity figures beside it.

---

## 11. The sweep entry point

`scripts/run_environment_sweep.py`. Full factorial over carbon source, growth phase,
temperature, aeration, anaerobic supplements, glucose charge, medium pH and cytosolic pH,
crossed with the dose ladder, the nutrient axis and the constructs; replicate wells simulated in
every cell; **one row per well carrying every condition that produced it.**

Two behaviours worth knowing before reading its output.

**Refusals are data.** A cell whose growth rate falls outside van Hoek's measured range is
written with a `physiology_refusal` string and empty flux columns rather than raising. The
default sweep grounds 54 % of rows and refuses 46 %.

**Environments that support no growth are skipped and counted.** An anoxic ethanol culture or an
incubator above 45.4 °C produces no data in a wet lab either, and giving it a reporter trace
would put a number where the measurement has none.

The design axes are also available without simulating anything, through
`generator.design.context_grid` and `generator.design.environment_grid`, so a power or
identifiability calculation can be run over the same space at no cost.

`collinearity` over that space is the reason for all of it, and the size of the effect is not
automatic:

| design | \|corr(stress, growth)\| |
| --- | ---: |
| one environment, dose ladder only | 0.895 |
| 3 temperatures × 2 aeration levels | 0.705 |
| 2 carbon × 5 temperatures × 3 aeration × 2 nutrient | 0.344 |

A narrow environment grid buys about 0.19 and is nowhere near enough. It takes a grid that
spans several axes together to get where a latent stress state is separable from growth.

---

## 12. Two existing tests now fail, and they should

`tests/test_three_sensor_build.py` goes from 30 passed to 28 passed / 2 failed. **The cause is
the temperature model and nothing else** — restoring the old
`1 - 0.04*|T-30|` in place and re-running gives 30 passed again, so this is attributed, not
guessed. That file is not this lane's and was not edited.

The mechanism is one line of its own fixture. `CONTEXTS` includes
`CultureContext(temperature_c=37.0)`, and a 37 °C culture's growth rate moves from
`0.40 × 0.72 = 0.288` /h to `0.40 × 0.916 = 0.366` /h. Growth rate sets reporter dilution in
every well of that context, so the fitted latent basis moves with it.

| test | asserts | now | how marginal |
| --- | --- | --- | --- |
| `test_pooling_the_same_wells_lifts_recovery_from_five_modules_to_nine` | 9 modules recovered | 10 | `heat` recovers at **0.2653** against a threshold of 0.25 |
| `test_a_stressor_holding_up_its_own_axis_does_not[antimycin_A]` | none of antimycin A's own targets recover | `{'ESR'}` recovers | ESR is a weight-0.35 target of antimycin A; its weight-0.9 `retrograde` and weight-−0.85 `atp` targets still do not recover |

Both are **Tier 1** figures in the sense of `docs/CLAIM_BOUNDARY.md` — counts fitted to this
generator and scored on the same generator — and both are exactly what that file says will move
when the generator's *structure* is perturbed rather than its noise. The direction is the one
the change predicts: making 37 °C a live experiment instead of a nearly-dead one gives the heat
module enough signal to cross a threshold it was sitting just under.

Neither failure contradicts the argument its test was written to defend. The pooling test's
claim is that pooling lifts recovery, and it still does, by more. The transfer test's stated
claim, in its own docstring, is that a stressor's **ATP** axis does not come back when it is
held out — and `atp` still does not; what crosses is an incidental weight-0.35 ESR target, so
the assertion is stricter than the sentence it defends.

The decision belongs to whoever owns that file, and there are only two honest options: re-pin
the counts with the new generator and say in the docstring that they moved and why, or narrow
the transfer assertion to the ATP claim the docstring actually makes. Loosening the temperature
model to keep the numbers is the one option that is not available.
