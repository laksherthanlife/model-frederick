# Two measurements the model needs and cannot infer

Both feed code that is already written and tested. Each section states what to run,
and the exact function that consumes the result.

---

## P1 — OD600 linear range on this reader and plate format

**Why.** `OpticalQualityGate.od_linear_max` currently carries a guess. It decides how
much of every plate is usable: at the old inoculation density it rejected 81% of one
experiment. In a 96-well plate the optical path is ~0.3–0.5 cm, so linearity ends well
below the cuvette figure of 1.0, but *where* is specific to this reader, this plate and
this fill volume.

**Run.**
1. Grow a culture of the working strain to OD600 ≈ 2.5–3 in the normal medium.
2. Prepare a two-fold dilution series in that same medium: 12 steps, stock to 1/2048.
3. Load at the **exact fill volume, plate type and lid/seal used in real runs** — the
   path length is set by volume, so a calibration at 200 µL does not transfer to 150 µL.
4. Include ≥3 media-only wells as blanks.
5. Read OD600 with the standard protocol, three technical replicates per dilution.

**Analyse.**
```python
from ystwin.calib.od import fit_od_calibration
calib = fit_od_calibration(dilution_factors, measured_od, blank=0.09, gdcw_per_od=None)
print(calib.summary())
gate = calib.optical_gate(tolerance=0.05)      # drops straight into G1
```

**What you get.** The ceiling is exactly `tolerance / k`, so a reader that saturates
twice as hard halves the usable range. Report `k` — it is a property of the instrument
and belongs in the methods section.

**Worth adding while the plate is set up.** Filter and dry a known volume of the stock
to get `gdcw_per_od`. It is the only thing standing between OD and a biomass state in
real units, and it is strain- and medium-specific, so a literature value will not do.

---

## P2 — Reporter loss rate (`k_deg`) for the actual construct

**Why.** `k_deg` separates a reporter that reports from one that reports `1/mu`.
Assuming zero is the textbook position for a stable YFP, and on the newer plates that
assumption holds. On the older, badly-conditioned plates it did not. Measuring it
settles the matter instead of inferring it from data whose OD channel is in doubt.

**Run — chase.**
1. Grow the reporter strain to mid-exponential under the standard protocol.
2. Split into three: (a) chase, (b) no-chase control, (c) reporter-free parent.
3. Add cycloheximide to (a) at 100 µg/mL to block translation. Vehicle only to (b), (c).
4. Read OD600 and mCitrine every 10 min for 6 h on the standard schedule.

**Analyse.**
```python
from ystwin.calib.kdeg import fit_decay_rate, fit_kdeg_from_chase
from ystwin.readings import RawOD, RawRFU
bleach = fit_decay_rate(t, control_rfu_per_od).rate        # from the no-chase control
# RAW readings, blank included: this function subtracts od_blank itself. Handing it a
# corrected trace is a TypeError rather than a number -- see ystwin.readings.
fit = fit_kdeg_from_chase(t, rfu=RawRFU(chase_rfu), optical_density=RawOD(chase_od),
                          od_blank=0.09, bleaching_rate=bleach)
print(fit.summary())
```
Residual growth and photobleaching are subtracted; a **negative** `k_deg` is reported
rather than clamped, and means translation was not actually blocked.

**Run — separating photobleaching from degradation (only if the chase shows real loss).**
Degradation runs with the clock; bleaching runs with the number of reads. They separate
only across read intervals. Run two identical plates, one read every 10 min, one every
30 min, then:
```python
from ystwin.calib.kdeg import partition_loss
split = partition_loss([(10/60, loss_fast), (30/60, loss_slow)])
```
One read interval gives only the total, however many plates are run.

**Free lower bound, no experiment.** `minimum_consistent_kdeg(t, RawOD(od), RawRFU(rfu))` returns the
smallest loss rate under which inferred promoter activity never goes negative. Use it
as the prior until the chase is done.

---

## P3 — Two controls worth adding to every plate

Neither is a separate experiment; both are wells.

**Reporter-free parent (BY4741).** Gives cellular autofluorescence as a function of
biomass, which is `ReporterOptics.autofluorescence`. Subtracting a guessed constant
instead biases the reporter estimate low at high density. The 2026-08-14 and 2026-08-07
BY4741 plates already do this — keep it standard.

**Non-producing isogenic strain — only once beta-carotene is in.** Not applicable to the
current biosensor strains, which carry no carotenoid pathway; the reporter channel is
clean of product absorbance today. Listed here so it is not forgotten later. Beta-carotene absorbs at
450–478 nm, on mCitrine's excitation. So the product dims the reporter: a path from the
quantity being predicted back into the channel predicting it. Shuffled-channel and
causal-null controls cannot detect it, because it is real signal.

To calibrate the coefficient, spike purified beta-carotene into wells of the
non-producer across the expected titre range and fit the attenuation. Then
`correct_inner_filter(rfu, carotenoid, optics)` removes it before any dilution
correction or state inference. Until it is measured, `ReporterOptics.inner_filter_coeff`
should stay `None` — correction is then refused rather than silently assumed absent.


## P4. Optical density linear range

`OpticalQualityGate.od_linear_max` is the one number the gate flags as most worth
measuring, and it is still a placeholder. Everything downstream divides by OD, so a growth
rate taken above the linear range is wrong by whatever the reader compresses.

1. Grow a culture to roughly OD 2 in the medium the experiments use, not in water --
   scattering is medium-dependent.
2. Prepare a two-fold series in the plate format the experiments use, at the same fill
   volume: eight or nine wells covering roughly OD 0.02 to 2. Path length is set by fill
   volume, so a different volume gives a different answer.
3. Include at least three medium-only blanks and subtract their mean.
4. Read on the same instrument, wavelength and settings as the experiments.
5. Record the nominal density of each well from its dilution factor and the blank-corrected
   reading, then:

```python
from ystwin.calib.od_linearity import fit_linear_range

fit = fit_linear_range(expected_od, measured_od, tolerance=0.05)
print(fit.summary())
```

6. Put `fit.od_linear_max` into `OpticalQualityGate` and re-run `scripts/run_gates.py`.

If `saturation_observed` is False the series never left the linear range, so the answer is
a lower bound -- extend the series upward rather than quoting the top well as the limit.

---

## P3 — Exponential fed-batch, the run that replaces a chemostat state

**Why.** Every steady-state number this repository produces divides by a growth rate, and
until 2026-09-02 that rate came from a chemostat pump. The chemostat is the one piece of
equipment in the design nobody has. An exponential fed-batch gives the same held growth rate
in a single vessel with no effluent line, no level control and no overflow weir — and gives
one thing the chemostat could not: **a growth rate that is measured rather than asserted.**

**The control law.** Feed

```
F(t) = F0 · exp(mu_set · t),      F0 = mu_set · X0 · V0 / (Y · Sf)
```

and nothing else. There is no feedback loop to tune, because the vessel supplies it: writing
`B = X·V` for total biomass, carbon limitation gives `q_S = F·Sf/B` and therefore
`d(ln q_S)/dt = mu_set − mu(q_S)`, whose fixed point is `mu = mu_set` for **any** monotone
uptake→growth relationship. The yield `Y` sets only where you start, not where you end up —
simulated with `Y` wrong by tenfold, `mu` still lands within 2% of setpoint
(`tests/test_fedbatch.py::TestTheSetpointIsAStableFixedPoint`).

**Plan the run.**
```bash
python3 scripts/design_fedbatch_run.py      # writes outputs/fedbatch_design.csv
```
```python
from ystwin.fba.fedbatch import design_fedbatch
design = design_fedbatch(
    0.101, 0.254,                       # setpoint, and the strain's measured mu_max
    initial_volume_l=1.0, max_volume_l=2.0,
    initial_biomass_g_per_l=0.05,       # LOW on purpose -- see below
    feed_substrate_g_per_l=500.0,
    assumed_yield_g_per_g=0.098,
    max_uptake_g_per_gdcw_h=2.59)
print(design.summary())
```

**Two ways it refuses, both before anything is inoculated.**

1. **The setpoint must sit below `mu_max` with margin.** The convergence rate is
   `q_S · dmu/dq_S`, which collapses as the uptake curve flattens, so at the ceiling the feed
   stops setting `mu` and past it glucose accumulates without bound.
2. **The feed must not open too rich:** `q_S(0) = mu_set / Y_assumed` has to sit under the
   strain's uptake ceiling. This one was found by simulation contradicting the algebra, and
   it is the sharper of the two — a culture whose feed opens above its ceiling starts
   *already* maxed out, grows at `mu_max` rather than the setpoint, and never enters the
   regime the feed controls. **A bigger vessel cannot rescue it**, because `F0` scales with
   `X0·V0` and `q_S(0)` divides that straight back out.

   The practical rule: **when the yield is uncertain, guess it HIGH.** A feed that starts too
   lean converges; one that starts too rich never begins.

> **Correction (2026-09-03).** An earlier version of this section called refusal 2 "the
> sharper of the two". It is not, at the numbers this protocol uses. Refusal 1 fires above
> `(1 − margin)·mu_max` = 0.2286 /h; refusal 2 fires above `Y·q_max` = 0.098 × 2.59 =
> 0.2538 /h. **Refusal 1 always fires first, by a factor of 1.110**, so refusal 2 is
> unreachable for these strains and the two are not independent guards in practice. Both are
> kept because a strain with a higher `mu_max` or a lower assumed yield reverses which binds,
> but at the Elizondo numbers only the margin does any work.

**Inoculate LOW.** The generation budget is
`log2(1 + (Vmax − V0)·Y·Sf/(X0·V0))`, and `X0` is the only term a protocol controls. On the
2 L plan above, dropping the inoculum from 0.5 to 0.05 g/L takes usable generations from
**1.99 to 5.29** — it is the single highest-value decision in the protocol.

**Sample, and re-anchor when you do.** A chemostat is sampled from its overflow for free;
here every sample removes biomass and pushes `mu` above setpoint. Ignoring that costs a few
percent RMSE. The fix is arithmetic, not equipment:
```python
from ystwin.fba.fedbatch import reanchor_feed
new_feed_l_per_h = reanchor_feed(design, measured_biomass_g_per_l, current_volume_l)
```

**Measure mu — do not assert it from the pump.** This is the step that makes the fed-batch
*better* than the vessel it replaces, not merely equal to it. A chemostat's `D = F/V` rests
on pump and level calibration and is never independently checked; here the biomass series you
are already taking gives `mu` with a standard error:
```python
from ystwin.fba.fedbatch import estimate_growth_rate
estimate = estimate_growth_rate(times_h, biomass_g_per_l, volumes_l)
print(estimate.summary())        # mu = 0.1010 +/- 0.0008 /h (0.76%) from 9 points over 24 h
```
Pass **total** biomass — the volume argument is not optional in spirit. The vessel is filling,
so concentration alone understates `mu` by exactly the volume's own growth.

That 0.76% is not a hopeful figure: at this repository's own OD600 replicate CV of 1.80%,
nine points over 24 h give `SE(mu) = CV/sqrt(Sxx) = 0.00077 /h`, and a 2,000-run simulation
of the estimator returns a median SE of 0.00074 /h against a spread of the estimate itself of
0.00077 /h — analytic and empirical agree to three digits. **This beats what a chemostat can
claim**, because `D = F/V` rests on pump and level calibration (typically 2–4%) and carries no
error bar at all. More points and a longer span both tighten it, and `Sxx` says by how much
before you run anything.

**What consumes the result.** The measured `mu` goes straight in as
`Environment(growth_rate_setpoint_per_h=...)`; `pathway/solve.py` divides by it exactly as it
divided by `D`, because the settled fed-batch residual substrate **is** the chemostat's Monod
residual at the same `mu`. Nothing downstream changed.

**It cannot reproduce EITHER Elizondo calibration state, and an earlier version of this
section wrongly said the substitution was exact.** The control law is derived from "under
carbon limitation every fed molecule is consumed, so `q_S = F·Sf/B`". Elizondo's six states
are glucose-**EXCESS**: at `mu = 0.101` the b-car2 culture takes up glucose at 6.68
mmol/gDCW/h while secreting ethanol at 6.40, a biomass yield of 0.084 g/g against the
0.47–0.50 g/g of a carbon-limited chemostat at the same rate — a 5.6–6.0× deficit, and
`docs/HARD_TESTS.md` §6 records the excess explicitly. **A glucose-excess culture has no
monotone `q_S → mu` map, so the fixed-point argument has nothing to stand on.** The upper
state additionally sits on `mu_max` with zero headroom. What the fed-batch gives is a
*carbon-limited* state at a chosen `mu`, which is a cleaner regime than the calibration set —
not the same one.

**And it is not the low-equipment option.** The shipped design reaches `X_final` = 24.5 g/L,
which at van Hoek 1998's own `qO2` = 2.5 mmol/gDCW/h at `D` = 0.10 demands **61.2 mmol
O₂/L/h**. A baffled shake flask delivers roughly 10–40. This protocol therefore needs a
sparged, stirred, DO-controlled vessel. **If that is not available, use §P4 instead** — the
balance does not require a held growth rate at all.

---

## P4 — The same measurement in a shake flask, with no chemostat and no fermenter

**Why this exists, and why it supersedes the equipment ask in P3.** Every recommendation in
this repository used to require a *held* growth rate, and therefore a chemostat — or the
fed-batch of §P3, which an audit on 2026-09-03 showed needs **61.2 mmol O₂/L/h** and so a
sparged, stirred, DO-controlled vessel. That premise was attacked and did not survive.

**The balance the model closes is per gram of biomass:**

```
d[X]/dt = v_in − v_out(X) − mu·[X]
```

**The vessel appears nowhere in it.** Integrating it forward at constant `mu` converges on
exactly the steady state `pathway/solve.py` returns, at a rate set by `mu` alone. So the cost
of having no chemostat is not a different answer — it is a **generation count**, and it is the
same count in either vessel (`scripts/batch_sufficiency.py` → `outputs/batch_sufficiency.csv`):

| within | generations | hours at mu = 0.101 | at mu = 0.2543 |
|---|---|---|---|
| 10% | 3.7–4.0 | 25.7 | 10.9 |
| 5% | 4.8–5.1 | 32.9 | 13.8 |
| 1% | 7.2–7.5 | 49.5 | 20.4 |

**`mu` must be KNOWN and roughly CONSTANT. It does not have to be HELD.**

**Run.**

1. **Two serial exponential precultures** in the same medium, never letting either leave
   exponential phase. This is not fussiness: starting cold, the entry enzyme's own dilution
   transient costs 6–9 generations on its own — more than a 20 g/L flask supplies. From a
   preculture already at the same `mu`, the deviation over eight generations is under 0.002%.
   Two flasks remove the problem entirely.
2. **Inoculate LOW.** This is the whole experiment. Generations of exponential phase available
   in 20 g/L glucose: `X0` = 0.5 g/L gives ~2.6, reaching only **−36%**; `X0` = 0.05 gives
   ~5.7, reaching **−5%**; `X0` = 0.005 gives ~9.0, reaching **−0.6%**.
3. **Read OD600 continuously** — a plate reader every 10 minutes is ideal, but manual flask
   reads work. Fit `mu` by regression on `ln(OD)`.
4. **Harvest for content at the end of exponential phase**, and take at least two harvest
   points so the approach can be seen rather than assumed.

**Certify steadiness on the OD trace, not on the content series.** The lag error tracks the
fractional **fall in `mu` per doubling** almost one-for-one: 1% drift per doubling → about
−1% bias, 5% → about −6%, 10% → about −15%. The OD trace has ~100 points at a 1.80% CV; the
content series has four at 5–10%. A flatness test on content **under-certifies by about
2×** — at 5% assay scatter it first declares "flat" around 4 generations, where the true bias
is still about −10%.

**What you get on `mu`, measured.** Using this repo's own `estimate_growth_rate` at its own
OD600 CV of 1.80%:

| method | points | SE(mu) | relative |
|---|---|---|---|
| §P3 fed-batch as written | 9 over 24 h | 0.00074 | 0.73% |
| **shake flask, manual reads** | 9 over 24 h | **0.00073** | **0.73%** |
| shake flask, short run | 6 over 10 h | 0.00197 | 1.95% |
| **plate reader, 10-min reads** | 145 over 24 h | **0.00021** | **0.21%** |

**A plate reader measures `mu` 3.6× more precisely than the fed-batch protocol does** — with
equipment this project already parses exports from.

**What consumes the result.** Exactly what P3's would: `Environment(growth_rate_setpoint_per_h=`
the measured `mu`. Nothing downstream distinguishes the vessel, which is the point.

**The one real limit.** `mu` in a flask is whatever the medium gives, not what you choose, and
the fitted window is `[0.101, 0.2543] /h`. Six of twelve carbon × temperature combinations land
inside it — glucose at 20 °C and 15 °C, galactose at 25 °C and 20 °C, ethanol at 30 °C and
25 °C. Temperature and carbon source are the two dials that move `mu` into range.
