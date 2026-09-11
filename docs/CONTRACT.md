# What this twin takes in, and what it gives back

A competing group can state its model's job in one sentence: *"fixed environment
`e = [I, O, S]` → complete product curve `P(t)`; the model never receives an environmental
time series."* This project could not, and that was a scope problem wearing a documentation
problem's clothes. This is the fix.

## The vocabulary

Borrowed rather than invented, from Kapteyn, Pretorius & Willcox, "A probabilistic
graphical model foundation for enabling predictive digital twins at scale", *Nature
Computational Science* **1**(5):337–347, 2021 (doi:10.1038/s43588-021-00069-0). Their
abstraction of an asset-twin system is six quantities. Filled in for this project:

| | Kapteyn's definition | This twin |
| --- | --- | --- |
| **S** physical state | "Parametrized state of the physical asset" | The culture's actual condition: activity of each stress regulon, specific growth rate, viability, intracellular redox and ATP |
| **D** digital state | "Parameters (model inputs) that define the computational models comprising the digital twin" | `biomass`, `reporter`, `promoter_activity`, `growth_rate` — the four states in `estimator.py`, plus the reporter kinetics `k_deg`, `k_mat` |
| **O** observational data | "Available information describing the state of the physical asset" | Synergy H1 kinetic export: OD600 and mCitrine per well per timepoint. Separately, Cq values from the qPCR anchor |
| **U** control inputs | "Actions or decisions that influence the physical asset" | Stressor and dose, nutrient level, inoculation density, run length, sampling interval — and, upstream, which reporters were built at all |
| **Q** quantities of interest | "Quantities describing the asset, estimated via model outputs" | Dilution-corrected promoter activity; module activities; fold change against the zero-dose control; and, design-ahead, the feasible product-flux range |
| **R** reward | "Quantifies overall performance of the asset-twin system" | **Not defined.** See below — this is a real gap, not an omission |

The distinction that earns its keep here is **S versus D**, and Kapteyn et al. state the
consequence directly: *"The digital state space will generally be only a subset of the
physical state space."*

That is not a caveat in this project; it is the central finding.
`docs/research/IDENTIFIABILITY.md` analyzes four fitted channels against the original
seven-programme design. Under the static factor-analysis assumptions, the Ledermann bound
permits one generic factor; it is not a universal ceiling for a supervised dynamical observer.
The broader catalogue is not a complete physical state either. The newer teacher/student
loop's three labelled synthetic aggregates and informative virtual reporters are a different
identification problem, not evidence that the biological state space has been recovered.

## The contract, at three ambition levels

### Level 1 — what can be committed to today

> **Given the control inputs `U` for one well and its observed OD600 and mCitrine time
> series `O`, return a posterior over that well's digital state `D` at every observed
> timepoint, and the dilution-corrected promoter activity `Q` that follows from it — with
> an interval, and with a refusal wherever the optical channel is not quantitative.**

Every component exists. `estimator.py` produces the posterior, `reporter.py` inverts the
dilution, `analysis/uncertainty.py` supplies the interval and refuses below three
biological replicates, `gates/g1_optical.py` supplies the refusal. What is missing is only
the wiring: the architecture map records that `estimator.py` — the repository's only
posterior — is imported by nothing but its own tests, and that every number in `outputs/`
is a point estimate.

**This is the contract to commit to.** It is honest, it is defensible, and it is one
integration away.

### Level 2 — the first scored prediction

> **Given `U` and the first `n` hours of `O`, return a predictive distribution over the
> remainder of the trace, scored against what actually happened.**

This is Level 1 plus a forward pass and a scoring rule. It needs no producing strain and no
new wet-lab work: the data is already on disk. It converts the twin from a describer into a
predictor, which is the single largest gap between this project and its comparators.

Two things make it real rather than nominal. It must be scored with a **proper scoring
rule** — CRPS or the weighted interval score — against a **baseline ladder**: a
mean-predictor, the naive uncorrected RFU/OD, and a growth-only predictor. And it must be
reported as a **skill score** relative to that baseline, so that doing worse than the
trivial prediction shows up as a negative number rather than hiding inside a plausible
RMSE. `analysis/nulls.py::skill_score` already implements that convention.

A caution carried from the verification pass: the weighted interval score's components are
**sharpness (interval width) plus over- and under-prediction penalties**. The three-way
"dispersion" naming in common use is downstream vocabulary from the `scoringutils` package,
not from Bracher et al., who use the word once and not as a component name.

### Level 3 — the ambition, declared as unbuilt

> **Given a fixed environment, return a distribution over the β-carotene titre curve
> `P(t)`.**

Blocked, and blocked for two independent reasons that should be stated together whenever
this is described. **No strain produces the product** — the constructs under
characterisation carry no carotenoid pathway, so nothing in the current data constrains this
layer and nothing in this layer constrains the current data. And **FBA cannot supply it**:
D1 established that the feasible product-flux range is `[0, ceiling]` at every growth level
tested, relative width 1.000, because a capacity bound moves the top of the interval and
never lifts the floor off zero.

Level 3 is a design target. `kinetic/carotenoid.py` exists so the product layer is settled
before the pathway lands rather than after. It should never appear in a results table.

## What the twin is never given

Worth stating explicitly, because it is what makes the contract falsifiable:

- **The true module activities.** They are latent, and per the identifiability work most of
  them are not recoverable from four channels regardless of design.
- **The plate map**, in the general case. It was *recovered* by matching against the team's
  own derived sheets and then confirmed against the logbook — which is a result, not an input.
- **Its own future observations.** `estimator.py::forecast` propagates from the last update
  only, and `update()` refuses an observation that predates the filter's position. That
  firewall is deliberate and must survive any smoother added for retrospective work.
- **Any measurement of autofluorescence.** It has never been made on any plate — see
  `docs/DATA_INVENTORY.md` — so the measurement model in `observation.py` carries a term
  that no data constrains.

## The reward, and why it is empty

`R` is undefined here, and that is worth naming rather than quietly dropping.

In Kapteyn's framework the reward is what turns a twin from a description into a
decision-support tool; the segment from now to the prediction horizon *"can be viewed as a
partially observable Markov decision process"*, and without a reward there is no policy and
no decision. This project has the machinery that *would* feed one — `analysis/design.py`
scores sensor sets on Fisher information, `analysis/power.py` costs a design in replicates,
`experiment_design.py` scores combination designs on transfer — but nothing composes them
into a single objective over which an experiment could be chosen.

That is the natural next layer and it is not yet built. Until it is, this is a **predictive**
twin and not a prescriptive one. On the DNV capability scale as rendered by San, Rasheed &
Kvamsdal (*GAMM-Mitteilungen* 44(4):e202100007, 2021 — note the scale originates with
DNVGL-RP-A204, not with the more commonly cited Rasheed 2020 paper), that places it at
**level 3, predictive**, reaching for level 4, prescriptive — and honestly at level 2,
diagnostic, until Level 1 above is actually wired up.

## One sentence

For the top of the README, and for anyone who asks what the model does:

> **Given a dose and a well's OD and fluorescence time series, this twin returns a
> posterior over that well's growth rate and promoter activity — dilution removed,
> uncertainty attached, and refused outright where the optical channel is not
> quantitative.**

Everything else here is either a diagnostic that supports that sentence, or design work for
a strain that does not exist yet.
