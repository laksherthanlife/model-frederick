# Superseded: the stress-to-ATP coefficient

**Current value:** `_MAX_STRESS_MAINTENANCE = 6.5` mmol ATP/gDW/h at full scale, on a
normalised activity, clipped at the measured range. Range produced: 0.70 to 7.20.

The old value was not merely unfitted. It was dimensionally incoherent, and it made the
branch unrunnable on real data.

---

## What it was

```python
_MAINTENANCE_PER_ACTIVITY = 300.0  # mmol ATP/gDW/h per unit promoter activity
...
atp_maintenance = resting_maintenance + maintenance_per_activity * excess
```

`excess` is promoter activity above basal, in **RFU/OD/h**. The comment calls the
coefficient "mmol ATP/gDW/h per unit promoter activity", which is a real unit only if
"unit promoter activity" means a normalised quantity near 1. It never did: the estimator
supplies dilution-corrected activity in reporter units, and on the four measured
constructs the excess reaches ~1305.

## What that produced

| construct | max excess (RFU/OD/h) | implied NGAM |
| --- | ---: | ---: |
| AlteredYap1 | 1187.8 | 356,341 |
| NativeYap1 | 1304.7 | 391,410 |
| UPRE1 | 1055.7 | 316,710 |
| UPRE2 | 1127.0 | 338,096 |

Against a resting NGAM of **0.7**. Bisecting Yeast9 at its default glucose bound, the
model goes infeasible at **19.19** mmol ATP/gDW/h, and growth there is already zero. So
the constant sat roughly **2 x 10^4** past the point where the cell stops growing, and
about 5 x 10^4 past the plausible envelope.

The branch had therefore never been run on a real activity. It is gated behind
`acknowledge_unvalidated=True` and refused by G4, so nothing had ever fed it one --
which is exactly how an error of this size survives.

## How the replacement is bounded

Not by fitting, and not by reading a value off a figure. Lahtvee 2016, **PMID 27307591**,
is the right anchor for the mechanism: glucose-limited chemostats at a fixed dilution
rate, "in which specific growth rate-dependent changes are eliminated", concluding that
raised **maintenance ATP** underpins the general stress response -- which is the mechanism
this branch already assumed. Their design is ethanol 20/40/60 g/L, NaCl 0.2/0.4/0.6 M and
33/36/38 C at D = 0.1 /h in CEN.PK113-7D. Their additional maintenance energy is published
as a bar chart in mmol ATP/gDW/h, so the numbers were not transcribed.

Instead the design was reproduced in the model. At D = 0.1 /h a biomass yield near 0.5
g/g puts q_glucose at 1.1-1.5 mmol/gDW/h, and the NGAM headroom above resting at
q_glucose = 1.5 is **6.54**. Hence 6.5 at full scale.

**An independent consistency check, not engineered.** Applying the resulting coefficient
at full stress gives NGAM 7.20, at which Yeast9 grows at **0.1002 /h** -- essentially
exactly the dilution rate Lahtvee held. An envelope derived from feasibility alone lands
on the growth rate of the paper it is anchored to.

## What is still asserted

The scale is Tier 0. It is an envelope consistent with a measured design, not a measured
value. Extracting Lahtvee's Figure 1B numbers, or their supplementary tables, would make
it Tier 1 and is the obvious next step. And normalising by the *observed* activity range
means the coefficient inherits that range: a plate with a brighter reporter would need it
restated, which is a real limitation and the reason `activity_full_scale` is an argument
rather than a hidden constant.

G4 still refuses the mapping. Fixing the arithmetic does not validate it.

## The lesson

A constant whose units are stated in a comment and never checked against the quantity it
multiplies can be wrong by any factor at all. This one was wrong by 10^4 to 10^5 and the
suite was green throughout, because nothing had ever called it with real data. Where a
coefficient converts between two units, the test that matters is not its value but
whether the product lands somewhere the model can survive.
