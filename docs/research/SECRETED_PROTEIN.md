# A secreted protein is the product this lab's instrument was built for

**2026-08-30.** Every product this repository has modelled is a small molecule made by a
linear enzyme chain: beta-carotene, PHB, glycogen. A **secreted recombinant protein** is a
different kind of thing, and it is the one where the biosensors this team actually built stop
being a side project and become the measurement.

## The connection, and it is mechanistic rather than thematic

`UPRE1` and `UPRE2` report the **unfolded protein response** through the Hac1/Ire1 pathway
(Cox & Walter, PMID 8898193, cited in `generator/stress_panel.py`). The ladder they were
characterised over is **DTT**, and DTT induces that response by **reducing disulfide bonds**.

Human insulin carries **three disulfide bonds** — two interchain, one intrachain within the A
chain. Forming each one in the ER runs through **Ero1**, a flavoprotein that takes the
electrons from protein disulfide isomerase and passes them to molecular oxygen, producing
**hydrogen peroxide**. That is the same chemistry as CrtI, corrected in `fba/carotenoid.py`
the same day this file was written, and for the same reason: a flavoprotein reduces O2 by two
electrons at one site.

So secreting a disulfide-rich protein loads **both** of this lab's sensor families at once,
and by two different routes:

| load | route | sensor | calibrated over |
| --- | --- | --- | --- |
| unfolded protein in the ER | folding demand exceeds capacity | `UPRE1`, `UPRE2` | DTT, 0–5 mM |
| hydrogen peroxide | Ero1 makes one per disulfide formed | `NativeYap1`, `AlteredYap1` | H2O2, 0–4 mM |

No other product in this repository does that. Carotenoid and PHB are read by neither.

## And the metabolic model cannot represent either of them

Checked against the vendored yeast-GEM v9.0.2, not assumed:

- **There is no protein-disulfide-formation reaction.** Searching reaction names for
  disulfide, thiol oxidation, Ero1 or oxidoreductin returns four hits and all four are
  cofactor recycling — glutathione oxidoreductase twice, a glutathione disulfide exchange,
  and a methionine/thioredoxin oxidoreductase. None forms a disulfide in a substrate protein.
- **There is no `ERO1` or `PDI1` gene** in the model at all.
- **There is no H2O2 species in the `er` compartment.** Hydrogen peroxide exists as `s_0837`
  (cytosol), `s_0838` (mitochondrion), `s_0839` (nucleus) and `s_0840` (peroxisome). The ER
  carries 135 metabolites and peroxide is not one of them.

So the oxidative cost of secretory folding — the thing that makes a disulfide-rich protein
expensive, and the thing this lab's oxidative sensors are built to read — was **absent from
the model**, not merely unparameterised.

## So it was built

`src/ystwin/fba/secretion.py`, the same day this was found. It adds what the GEM lacks rather
than reporting that the GEM lacks it:

- **`h2o2_er`**, hydrogen peroxide in the ER compartment, which did not exist.
- **`H2O2ter`**, its export to the cytosol — where the model's catalase (`r_0255`) is, and
  where `NativeYap1` and `AlteredYap1` report from. Without a sink the folding reaction
  optimises to zero and reads as a capacity limit rather than a missing reaction.
- **`folding_reaction(model, protein, disulfides)`**, which charges one O2 and one H2O2 per
  bond and gives the folded species a formula **two hydrogens lighter per bond** —
  `2 R-SH + O2 -> R-S-S-R + H2O2` moves them onto the peroxide. The first version did not,
  and `check_mass_balance` reported it unbalanced by exactly 2n hydrogen.

It is generalisable by construction: the peroxide pool and its export install once per model,
and any secreted protein attaches its own disulfide count to them. Insulin's three bonds are
not special-cased anywhere.

**What is deliberately not modelled** is the protein thiols themselves. A cobra model has no
species for "reduced cysteine pair in a folding polypeptide", and inventing one would mean
inventing its pool size and its turnover. What is installed is the NET oxygen and peroxide
accounting, which is what a flux balance can use and what a Yap1 sensor can read. It is not a
mechanism and should not be read as one.

`tests/test_secretory_folding.py` pins both directions: the folding reaction carries flux
aerobically, and returns exactly zero when oxygen uptake is closed.

## What that means for a claim about insulin

An FBA ceiling on insulin can be computed and it will be honest about amino acids and ATP. It
will be **silent about the cost that matters**, because the reaction generating that cost is
not in the model. Anyone quoting such a ceiling as a yield should be told that.

This is not a reason to skip the calculation. It is the reason to say what the calculation
leaves out, which is what `src/ystwin/fba/insulin.py` records.

## What would close it

In cost order, and the first is free:

1. **Nothing new.** The DTT and H2O2 ladders already measured are the dose-response of the two
   loads a secreted protein imposes. What is missing is not data but a strain that expresses
   one — and the sensors would then be read against a load the lab created rather than one it
   pipetted.
2. ~~**An `ERO1` reaction and an ER peroxide species** in the GEM.~~ **Done** —
   `src/ystwin/fba/secretion.py`, above.
3. A secreted-protein strain carrying one of the existing sensors. That is the experiment this
   whole instrument has been pointing at, and this repository has never said so.
