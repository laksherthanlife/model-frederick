# Held-out splits: which claims this dataset can support, and which it cannot

Every held-out number in this repository was produced by a split written inline in
whichever script needed one. No two of them partition the data the same way, none was
recorded, and none can be re-created by a reader. A held-out R² quoted against an unnamed
partition is not a claim about generalisation — it is a claim about an arrangement of rows
that no longer exists.

`src/ystwin/analysis/splits.py` makes the partition an artefact. `scripts/make_splits.py`
writes it to `outputs/split_manifest.csv` — one row per group, carrying its assignment, the
seed, the group key and a SHA-256 hash of the whole partition — with counts in
`outputs/split_manifest_summary.csv`. A result cites a row. A reader re-runs the script and
compares hashes.

The short version of what follows: **on the real biosensor data, three of the five
applicable split kinds are available today and two are blocked, and the one that would
license a genuinely forward-looking claim is blocked for every construct.**

---

## 1. The two constraints that decide the split kinds

**The unit of replication is the plate.** Wells on one plate share an inoculum, a medium
batch, a reader and a position in an incubator. A split that puts wells from one plate on
both sides has trained on the batch effect it is about to be scored against.
`analysis/uncertainty.py` refuses that mistake on the analysis side and `analysis/power.py`
on the design side; `splits.py` is the matching refusal for evaluation, and it is
structural rather than conventional: the caller declares the group key, and a split whose
groups straddle the axis being held out is refused, not permitted.

```python
make_split(frame, "heldout_replicate", group_key=("construct", "dose_mM"))
# SplitNotPossible: the heldout_replicate split holds out whole values of 'plate', but
# the group key ['construct', 'dose_mM'] produces group(s) containing more than one …
# Assigning such a group puts one plate on both sides of the split.
```

**The two sensor pairs were challenged with different agents.** UPRE1 and UPRE2 saw DTT,
NativeYap1 and AlteredYap1 saw H₂O₂, and their dose axes are not comparable.
`qpcr.assert_single_stressor` already refuses to pool them; the construct-level split
delegates to it, so the refusal and its wording come from the one place that records which
agent challenged which construct. The frame's own labelling is checked against
`qpcr.STRESSOR_FOR_CONSTRUCT` rather than trusted, so a mislabelled inventory is caught
before it becomes a committed manifest.

---

## 2. The split kinds

| Kind | Held out | Stratum | What a score licenses |
| --- | --- | --- | --- |
| `interpolation` | one **interior** dose rung | construct (real) / stressor (panel) | the dose response is smooth between the rungs that were run |
| `extrapolation` | the **top** dose | construct / stressor | the dose response continues **above** the top rung that was run |
| `heldout_replicate` | one whole plate | construct | the result would have held on a plate that was never run |
| `heldout_construct` | one construct | stressor | the model is about the pathway, not the promoter it was fitted to |
| `heldout_stressor` | one stressor, **and every combination containing it** | — | a latent stress state reaches a stressor it never saw |
| `heldout_combination` | one co-applied pair, components kept in training | — | the interaction is predicted from the parts |

### `interpolation` vs `extrapolation`, and why they are two kinds

These are different claims of very different difficulty, and conflating them is the
specific error this module exists to prevent. Interpolating between rungs that were run is
a smoothness claim. Predicting above the top rung is the claim the dose-response work in
this repository has the most reason to doubt: the top of every ladder is where growth
arrest, dilution artefacts and the unmeasured autofluorescence term all concentrate — see
`docs/DATA_INVENTORY.md`, where the 2–5 mM DTT readings are exactly the ones that move under
a plausible autofluorescence share while the 0.1–0.2 mM readings do not. Neither group is
resolved by that stability: the registry withdrew both the "entirely dilution" and the "no
induction" readings for want of predeclared equivalence margins.

So they get separate names, separate partitions and separate hashes, and `extrapolation`
is **deterministic**: the top dose is the top dose, no seed moves it, and asking for a
different one is refused rather than granted.

```python
make_split(frame, "extrapolation", group_key=..., hold_out=0.5)
# SplitNotPossible: dose_mM 0.5 cannot be the held-out value of the extrapolation split
# for construct=UPRE1; that split holds out the top rung, which is ['5']
```

The top dose is taken **per construct**, because DTT tops out at 5.0 mM and H₂O₂ at 4.0.
A pooled split would hold out only DTT's top rung and none of H₂O₂'s.

`interpolation` never holds out the lowest rung. Predicting below the training range is
extrapolation downward wearing the wrong name, and on the real ladder the lowest rung is
the 0 mM control that is the denominator of every fold change — withholding it leaves a
training set that cannot express the quantity being predicted. It also requires four rungs:
holding out the middle of three leaves two, and two points cannot distinguish a curve from
a line, so the split would not test the thing it is named for.

### `heldout_replicate`, the only forward-looking claim

Leave one plate out entirely. This is the only split here whose score answers "would this
have worked on a plate we had not run", which is the question anyone reading the results
actually has.

The training side must keep at least `uncertainty.MIN_PLATES_FOR_INTERVAL` = 3 plates. That
number is not invented for this module: it is this package's own answer to how few clusters
admit a statement at all, and below it `fold_change` returns `INCONCLUSIVE` rather than a
tight-looking number. A claim about generalising across batches rests on having measured
batch variation, so a training side on which the package refuses to state an interval
cannot support one. **Four plates is therefore the minimum.**

### `heldout_construct`: one training construct is enough, one training plate is not

Train on UPRE1, predict UPRE2 — same stressor, different promoter. One training construct
is allowed here, unlike one training plate, and the distinction is not arbitrary: this is a
named one-to-one transfer claim whose scope is stated by naming both constructs, whereas a
single training plate leaves the very quantity the replicate split is about — batch spread
— unmeasured.

### `heldout_stressor` takes the combinations with it

Holding out H₂O₂ while `DTT+H2O2` stays in training does not hold out H₂O₂: the model has
seen it, at a dose, alongside another agent. The leak is invisible in the counts, because
the held-out label really was withheld. Every treatment containing the agent goes with it,
which is why the panel's `heldout_stressor` test side has **two** groups (`H2O2` and
`DTT+H2O2`) and not one.

### Deliberately absent: a plain random split

A random well-level train/test split leaks the plate and leaks the dose group, scores
near-perfectly for that reason, and licenses nothing. It is not offered.

`validation` is an **assignment**, not a kind: further units moved out of training by the
same grouping discipline, so that model selection never touches the held-out set. For
`extrapolation` the validation units are the next rungs down, so tuning happens on the same
question as scoring.

---

## 3. What is available today

Regenerate with `python3 scripts/make_splits.py`. Numbers below are seed 0.

### Real biosensor plates — `dataset = real_biosensor`

Coverage read off the exports, not off the protocol. Replicates 2 and 4 fluoresced columns
1–3 only, which under the recorded layout is UPRE1:

| Construct | Plates read in mCitrine | Plates with a blank in mCitrine |
| --- | ---: | ---: |
| UPRE1 | 4 | 2 |
| UPRE2 | 2 | 2 |
| NativeYap1 | 2 | 2 |
| AlteredYap1 | 2 | 2 |

| Kind | Status | Train / test groups | Train / test rows |
| --- | --- | ---: | ---: |
| `interpolation` | **available** | 60 / 10 | 180 / 30 |
| `extrapolation` | **available** | 60 / 10 | 180 / 30 |
| `heldout_construct` | **available** | 42 / 28 | 126 / 84 |
| `heldout_replicate` | **blocked for all four constructs**; UPRE1 only under a caveat | 21 / 7 | 63 / 21 |
| `heldout_stressor` | **impossible by design** | — | — |

Held out at seed 0: interpolation 0.2 mM (UPRE1, AlteredYap1) and 0.5 mM (UPRE2,
NativeYap1); extrapolation 5.0 mM (DTT constructs) and 4.0 mM (H₂O₂ constructs);
`heldout_construct` UPRE2 within DTT and AlteredYap1 within H₂O₂.

### Simulated panel — `dataset = simulated_panel`

25 stressors plus the `DTT+H2O2` combination, 5 geometric rungs each, 3 replicate wells.

| Kind | Status | Train / validation / test groups | Rows |
| --- | --- | ---: | ---: |
| `interpolation` | available | 104 / 0 / 26 | 312 / 0 / 78 |
| `extrapolation` | available | 104 / 0 / 26 | 312 / 0 / 78 |
| `heldout_stressor` | available | 23 / 1 / 2 | 345 / 15 / 30 |
| `heldout_combination` | available | 25 / 0 / 1 | 375 / 0 / 15 |
| `heldout_replicate` | not applicable — no plate axis | — | — |
| `heldout_construct` | not applicable — no promoter axis | — | — |

The two datasets support **disjoint** claims. The panel has the stressor coverage and no
batch structure; the real plates have the batch structure and two stressors that cannot be
separated from their reporters. Nothing scored on the panel is evidence about a plate, and
nothing scored on the plates is evidence about an unseen agent.

### Partition hashes, seed 0

| Dataset | Kind | Hash |
| --- | --- | --- |
| real_biosensor | interpolation | `90c3f17dfe07d3296f7259785b6e75d6ac9d1b6943520583bc451c720b886ef0` |
| real_biosensor | extrapolation | `1558856c1ff26c946eabb1a189f1668957609daa91ddc7ba8e22687302f7bf41` |
| real_biosensor | heldout_construct | `b242b5ba797c84a751bc9c959818248bd947bdc4a465289c3121671c5af8a534` |
| real_biosensor | heldout_replicate (UPRE1 only) | `19bdc70068fe7e691c5703b04bd2b84ace703eb6a4708671301b36e4228122a0` |
| simulated_panel | interpolation | `e82a651624efc2b1612985060ebd2f0a4e9a412400874caf4065fb09e994a423` |
| simulated_panel | extrapolation | `aa0c6e84ee136633b9415a505e00cdc4bf9185bcd61dc502ac20961d08068c46` |
| simulated_panel | heldout_stressor | `bb5e875c8ad6ce747eebfaa0930a13c2ae16286112471ad4cfc84f868d543d67` |
| simulated_panel | heldout_combination | `4974190efd4b1d65e5e66d59548acb63e0bad41cebd2634564ca663aad81a582` |

The hash is taken over the sorted group → assignment mapping, so it does not move with row
order, dict iteration order, pandas version or machine. Two results that quote the same
hash were scored on the same partition. Quote the first 12 characters in prose.

---

## 4. What is currently impossible, and why

### `heldout_replicate` — blocked for every construct

This is the important one, because it is the only split that licenses a forward-looking
claim, and it is exactly the one the dataset cannot support.

At **n = 2 biological replicates**, holding one out leaves one plate to train on. That is
not a small training set; it is a training set on which this package's own machinery
refuses to state anything (`MIN_PLATES_FOR_INTERVAL = 3`), so the split would produce a
held-out score whose training side has no measured batch variation at all. The module
refuses rather than returning it:

```
the heldout_replicate split is not available on this dataset:
  2 distinct plate value(s) for construct=AlteredYap1; the heldout_replicate split needs at least 4;
  2 distinct plate value(s) for construct=NativeYap1;  the heldout_replicate split needs at least 4;
  2 distinct plate value(s) for construct=UPRE2;       the heldout_replicate split needs at least 4.
That split would have licensed: the result would have held on a plate that was never run.
```

Per construct, on the plates read in the reporter channel:

| Construct | n plates | Replicate split |
| --- | ---: | --- |
| UPRE1 | 4 | available — **but see the caveat below** |
| UPRE2 | 2 | **blocked** |
| NativeYap1 | 2 | **blocked** |
| AlteredYap1 | 2 | **blocked** |

**The UPRE1 caveat.** UPRE1 reaches four plates only by counting 20260728 and 20260804,
neither of which carries a blank in the reporter channel. The background term is unmeasured
on both, so a score on that partition is licensed only under a stated background bound —
`analysis/uncertainty.py` supports that shape of statement and it is item 5 of the cost-ordered
list in `docs/DATA_INVENTORY.md`. On the plates whose blanks were also read — the set every
current result uses — **UPRE1 is at n = 2 as well, and the replicate split is blocked for
all four constructs.**

The manifest carries the UPRE1 partition anyway, because it is the only forward-looking one
the real data offers and hiding it would not make it less true. Anything scored on it must
name the caveat.

**What would lift it:** two more full biological replicates read in all twelve columns,
with blanks inside every channel that is read. That is items 2–4 of `docs/DATA_INVENTORY.md`
and it takes the whole four-construct comparison to n = 4 — enough for a replicate split
with a three-plate training side, on every construct.

### `heldout_stressor` — impossible by design, not by sample size

More plates will not fix this one. Every construct in the real design sits under exactly
one stressor: UPRE1 and UPRE2 under DTT, NativeYap1 and AlteredYap1 under H₂O₂. Holding out
H₂O₂ therefore also holds out both Yap1 reporters entirely, and nothing in the resulting
score separates "unseen agent" from "unseen reporter". The refusal is detected from the
frame's own crossing rather than from the construct names, so it fires for any dataset
built the same way:

```
every construct in this frame appears under exactly one stressor
({'AlteredYap1': 'H2O2', 'NativeYap1': 'H2O2', 'UPRE1': 'DTT', 'UPRE2': 'DTT'}), so
holding out a stressor also holds out its constructs entirely and the score cannot
distinguish an unseen agent from an unseen reporter; use heldout_construct within a
stressor instead
```

**What would lift it:** running at least one construct under both agents, so that stressor
and reporter stop being the same column. Failing that, `heldout_stressor` is a
simulated-panel claim only, and the leave-one-stressor-out numbers in
`outputs/transfer_by_stressor.csv` should be read as statements about the generator rather
than about the bench.

### `heldout_combination` — nothing to hold out

No real plate here co-applied two agents, so there is no interaction on the bench to
predict. Simulation only.

---

## 5. Using a manifest

```python
import pandas as pd
from ystwin.analysis.splits import make_split

split = make_split(inventory, "extrapolation", group_key=("plate", "construct", "dose_mM"),
                   seed=0, dataset="real_biosensor")
side = split.assign_rows(inventory)          # a Series of train / validation / test
train, test = inventory[side == "train"], inventory[side == "test"]
print(split.hash)                            # cite this next to the score
```

`assign_rows` refuses a frame the manifest does not describe, so a partition cannot be
silently applied to different data.

`feasibility(frame, kind, group_key=…)` returns the same checks as a table instead of
raising, which is where the counts in section 4 come from — they are recomputed on every
run rather than typed in once.

## Column reference: `outputs/split_manifest.csv`

| Column | Meaning |
| --- | --- |
| `dataset` | `real_biosensor` or `simulated_panel` |
| `split_kind`, `seed` | which kind, drawn with which seed |
| `group_key` | the columns whose combinations are indivisible, `\|`-joined |
| `within` | the stratum columns; each stratum was split independently |
| `axis` | the column whose values were held out |
| `stratum`, `group_id` | which stratum this group is in, and its id |
| `construct`, `stressor`, `plate`, `dose_mM` | the group key exploded, blank where not part of it |
| `axis_value` | this group's value on the held-out axis |
| `assignment` | `train`, `validation` or `test` |
| `n_rows` | wells in this group |
| `held_out` | the axis values withheld **in this row's own stratum** |
| `partition_hash` | SHA-256 of the whole sorted group → assignment mapping |
