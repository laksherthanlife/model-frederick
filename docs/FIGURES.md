# The figures

Four figures carry the argument in this repository. They are not a gallery: each one
settles a question that a table alone does not, and two of the four settle it negatively.
Everything else worth drawing is either a diagnostic (the report set, below) or a number
better read from the CSV it lives in.

Regenerate all of them, and the report set, with one command:

```bash
python3 scripts/make_figures.py
```

That writes the story set as SVG and 200-dpi PNG into `outputs/`, alongside
`outputs/figures_manifest.csv`, which maps every figure to the committed tables it reads
and to the script that writes each of those. Nothing in a figure is transcribed from prose
in `docs/`: if a number is not in a committed CSV, it is either derived from one by a
builder in `viz/figures.py` or it is not on the figure.

The conventions are enforced in the drawing code rather than left to whoever renders next.
Every categorical distinction is doubled by marker shape, marker fill or line style, so a
greyscale print keeps every distinction the colour made; the palette is Okabe-Ito and
refuses to cycle rather than giving two categories one colour. Axis labels carry units.
Each figure's own caption block, along the bottom of the canvas, names its source tables
and states what a reader should conclude — including when the conclusion is that nothing
can be concluded.

---

## fig01 — The dilution correction, with intervals

**Question.** After correcting for growth dilution, which apparent inductions are still
there, and which of those can be told apart from no change at all?

**Source.** `outputs/sensor_characterisation.csv` — `python3
scripts/run_sensor_characterisation.py`. Folds and intervals via
`viz/figures.py::dose_response_table`, which calls
`analysis/uncertainty.py::fold_change`.

**What it shows.** Three things at once, per reporter. The open dashed series is the raw
per-cell readout, which rises to about 2-fold at the top of both ladders. The filled solid
series is the same data after the reporter ODE is inverted for growth dilution, and it
rises far less. The whisker is the 95% cluster-bootstrap interval with the plate as the
resampling unit, which is what decides whether either number means anything.

**Conclude.** Six of twenty corrected calls clear 1.0 — UPRE1 and UPRE2 at 0.2, 0.5 and
1.0 mM DTT — and every other corrected call is inconclusive, drawn hollow. At 2 and 5 mM
DTT, where the naive readout is at its highest, the corrected fold collapses to about 1.0:
those apparent inductions are the growth slowdown, not the sensor. At 2 and 4 mM H2O2 the
recovered activity goes negative and there is no ratio to take at all; those doses are
marked as gaps rather than plotted.

This figure changed shape when a third biological replicate was recovered. At two plates
`fold_change` refuses a cluster bootstrap outright, so the same figure previously carried
no intervals and said so on its face. It now carries them, and the note switches between
the two statements from the data rather than from a hardcoded sentence.

## fig02 — The held-out prediction, and where it fails

**Question.** Over what dose range does the fitted dose response beat carrying the nearest
measured dose forward, and why does it stop?

**Source.** `outputs/heldout_interpolation_by_dose.csv` — `python3
scripts/run_heldout_score.py`. The EC50 on the x axis is refitted from
`outputs/sensor_characterisation.csv` by `viz/figures.py::identifiable_ec50`, because no
committed table carries it: `outputs/panel_calibration.csv` leaves `fitted_ec50` empty,
since that module refuses a fit it cannot identify.

**What it shows.** Skill against the nearest-dose baseline, per held-out interior rung,
plotted against dose in multiples of the fitted EC50 rather than against dose in mM — that
is the axis on which the sign separates. The zero line is the point of the figure: below it
the model is worse than persistence. The lower panel is the mechanism and not more data:
the induction term is Michaelis-Menten, so its sensitivity to dose relative to its value at
`d = K` is `4 / (1 + d/K)²`, a property of the functional form with no free parameter.

**Conclude.** The fit beats the baseline at 0.1, 0.2 and 0.5 mM and loses at 1.0 and 2.0 mM,
with the crossover between 3.5 and 7.1 multiples of the fitted EC50. It is a **modest**
win, not a decisive one: +0.086, +0.368 and +0.101. A prediction from this model outside
its EC50 range should not be used, because deep in saturation the model cannot move its
prediction with dose while the measured activity still varies.

**On staleness.** This table was first computed at two biological replicates and its skill
values were roughly twice as high. `heldout_skill_table` therefore checks freshness before
the figure is drawn, and states the verdict in the caption: it multiplies the vouched fit
count the scorer recorded by the wells per construct in the current characterisation table,
and compares that with the test rows actually scored. If the characterisation table has
grown since the scorer last ran, the figure says STALE rather than plotting old numbers
quietly.

## fig03 — What transfer does and does not recover

**Question.** What does the latent stress state actually recover on an unseen stressor, and
does the fitted basis beat a random subspace of the same rank?

**Source.** `outputs/external_redundancy.csv`, `outputs/external_nulls.csv` and
`outputs/external_nulls_independent.csv` — `python3 scripts/run_external_validation.py`;
`outputs/cross_family_transfer.csv` — `python3 scripts/run_cross_family.py`.

**What it shows.** Panel A is recovery against the number of *other* stressors in the panel
that drive the same module, under both peer definitions, on the Gasch 2000 arrays. Panel B
is every headline claim beside the null it was tested against, on the axis
`observed − null median`; zero means the null did as well. Failures are grouped at the top.

**Conclude.** Both results are negative, and both are the honest description of what the
transfer numbers mean.

- A module is recovered when some other stressor in the panel drives it, and not otherwise.
  What looks like generalisation to an unseen stressor is interpolation across redundant
  coverage — a fact about the design of the panel, not about the latent state.
- The transfer score does not beat a random rotation of the same rank. Not in simulation
  (+0.0499 against +0.0358, p = 0.333), not on the eight-stressor real panel (+0.6876
  against +0.6953, p = 0.990), and not on the six-stressor panel with nothing traceable to
  Gasch (+0.6048 against +0.6077, p = 0.726). It does clear the much weaker
  `matched_marginals` null, which says only that the channels are correlated at all.

The one claim in the panel that beats its null is the readout through the literature
loading matrix, which is *given* rather than fitted. That asymmetry is the finding: what is
known works, and what is estimated does not.

**Not on this figure.** The simulation-side peer-count table — 0 peers recovering 0 of 6,
6+ recovering 84%, Spearman +0.591 — lives only in
[`superseded/module-transfer-inflation.md`](superseded/module-transfer-inflation.md) and no
committed CSV carries it, so panel A plots the real-data version instead of transcribing
the simulated one.

## fig04 — Do the resolved calls survive the analysis choices

**Question.** Are the six resolved calls properties of the plates, or of where the late
window was cut and how much autofluorescence was assumed?

**Source.** `outputs/late_window_sensitivity.csv` and
`outputs/autofluorescence_sensitivity.csv` — `python3
scripts/run_sensor_characterisation.py`. Which calls to track comes from the same fold
table fig01 draws, passed in rather than recomputed, so the two figures always describe the
same claims.

**What it shows.** The *lower* bound of each interval — not the point estimate — against
each knob, on a shared y axis so that a knob which moves nothing does not look like one
that moves a lot. A point estimate that stays above 1.0 while its lower bound crosses has
not survived anything, which is why the lower bound is the quantity plotted.

**Conclude.** Four of the six calls keep an interval clear of 1.0 everywhere swept. All six
survive the assumed autofluorescence across 0-30% of the 0 mM reporter signal, and that
sweep barely moves the bounds at all. Two lose their lower bound at an extreme of the late
window: UPRE1 at 0.2 mM when the late fraction is 0.9, and UPRE2 at 0.5 mM at 0.8 and 0.9.
The point estimates stay above 1.0 throughout — but the interval is the claim, so the
honest statement is four robust calls and two that depend on the window.

---

## The report set

`figures/01`-`05` are a different set for a different purpose: `scripts/build_report.py`
inlines them into one self-contained HTML page, one figure per section, including the two
that are simulated (the decoupling grid and the attribution power) and the one that exists
to show a refusal (the G4 anchor contamination). They are regenerated by the same command.
Only their SVGs are tracked; the PNGs and the HTML are not, because they are derived from
the SVGs and inflate every diff.

`figures/02_dose_response_correction` and `outputs/fig01_dose_response_correction` are the
same drawing from the same builder, rendered into both sets, because the report needs it in
its reading order and the story set needs it as the headline.
