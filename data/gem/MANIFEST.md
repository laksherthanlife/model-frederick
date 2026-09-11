# Genome-scale models, tracked

Both are published models, vendored gzipped so every FBA result here is reproducible from
a clone. cobrapy reads `.xml.gz` directly; `ystwin.paths` prefers these over any
sibling checkout.

| file | model id | reactions | metabolites | sha256 of the uncompressed SBML |
| --- | --- | ---: | ---: | --- |
| `yeast-GEM.xml.gz` | yeastGEM_v9.0.2 | 4131 | 2806 | `9fd2c572cace73c2ea835205617313554d1bf89f4ef077f49defbdfa219a4ad7` |
| `ecYeastGEM_batch.xml.gz` | GECKO batch, on yeast-GEM 8.3.4 | see below | | `3a9d97a68e2408ea243114794c027b11d70f253ecf6a758f4bcd5f2de50170b8` |
| `ecYeastGEM_yeast902.xml.gz` | GECKO 3.2.5, on yeast-GEM **9.0.2** | 8269 | 3969 | `0a725a63cf0395a2330eff3134713255ce4704af903e59ad36dcf4b62ca509de` |

Sources: yeast-GEM v9.0.2 from SysBioChalmers/yeast-GEM releases; ecYeastGEM_batch from
SysBioChalmers/ecModels. Both are distributed by the Systems Biology group at Chalmers
under the terms stated in their repositories.

**The version mismatch is real and load-bearing.** The ec model is built on yeast-GEM
8.3.4, not 9.0.2: its metabolite ids are compartment-suffixed and its charges are all
zero. `paths.ec_yeast_gem` says so, and `fba/physiology.py` records why the ec model is
the one used for the product ceiling -- plain Yeast9 handed the measured uptakes still
overpredicts growth by 73% (0.69 against a measured 0.40 /h), where the ec model is
within 5.8% from its protein pool alone.

**The ec model is irreversibly split.** `r_1714` is pinned at (0, 0) and uptake flows
through `r_1714_REV`, which is why `fba/physiology.py::cap_uptake` reads the split from
the model rather than trusting a reaction id.

## Upgrades: what was checked on 2026-09-03, and what is blocked

**yeast-GEM v9.1.1 exists (released 2026-08-30) and is mechanically safe to adopt, but is
NOT adopted.** `scripts/gem_version_diff.py` measures the change and writes
`outputs/gem_version_diff.csv`; rerun it against any candidate before adopting one.

    reactions   4131 -> 4105   (31 removed, 5 added)
    metabolites 2806 -> 2748
    genes       1161 -> 1143
    every identifier this repository pins by name survives -- see PINNED in that script
    growth falls a UNIFORM 5.53% at glucose -1.0, -1.5 and -10.0 with oxygen unlimited

It is held for two reasons, neither of which is difficulty. First, **nothing is fixed by it**:
the change is a uniform scaling, so the overprediction against Elizondo's six measured states
improves only from +309..+526% to +287..+491% and stays useless there -- plain Yeast9 has no
Crabtree mechanism and those cultures are glucose-excess and fermenting, which is why the ec
model carries the ceiling. Second, adopting it is a **regeneration event**: every marked GEM
number in `outputs/maintenance_scale.csv` and `outputs/gem_environment_*.csv` moves, along
with every docstring marker pinned to those cells. That is a deliberate, self-contained piece
of work and it should not ride along with anything else.

**UPDATE, later on 2026-09-03: the rebuild was done, and it is both newer and better.**
`ecYeastGEM_yeast902.xml.gz` is a GECKO 3.2.5 ecModel on this repository's own yeast-GEM
9.0.2, built by `scripts/gecko/build_ecyeastgem.m`. Measured against the adapter's
experimental growth rate of 0.41 /h on unlimited glucose:

    vendored ecYeastGEM_batch, on 8.3.4     0.3768 /h    -8.1%
    ecYeastGEM_yeast902, on 9.0.2           0.3821 /h    -6.8%

**It is NOT the default, and the reason changed once it was measured.** The swap was made,
scored, and reverted the same day. Against van Hoek 1998's aerobic batch phenotype
(growth 0.40 /h, glucose 11.1, ethanol 13.9 mmol/gDCW/h):

    model    growth          glucose          ethanol
    8.3.4    0.377  (+5.8%)  17.93  (+61.5%)  29.57  (+113%)
    9.0.2    0.382  (+4.5%)   4.39  (-60.5%)   0.00  (-100%)

The 9.0.2 build is FULLY RESPIRATORY -- same growth on 4.4 glucose and 10.0 oxygen, no
ethanol at all -- where the 8.3.4 model ferments. Every culture this repository is
calibrated against ferments: Elizondo's six states secrete 2.9-16.1 mmol/gDCW/h at BOTH
dilution rates. A model that cannot produce overflow is qualitatively wrong for them, and a
percentage point of growth accuracy does not buy that back. Overflow emerges in an ecModel
when the protein pool is tight enough to trade protein-expensive respiration for
protein-cheap fermentation; the DLKcat-derived kcats here are evidently generous enough that
respiration stays affordable, which is a sigma/pool calibration exercise rather than a
rebuild. `tests/test_ec_model_choice.py` pins all of this so the swap cannot happen quietly.

One build artefact was fixed before scoring: the exported model shipped `r_1714` at -1.0,
a leftover bound that starved it to 0.086 /h and made the first comparison meaningless. It
now ships unlimited carbon and oxygen, as the 8.3.4 batch model does, so the protein pool is
what binds. `fba/physiology.py::cap_uptake` already handles both, because it
reads the irreversible split from the model rather than trusting a reaction id: GECKO 2 pins
`r_1714` at (0,0) and supplies through `r_1714_REV`, GECKO 3 leaves `r_1714` unsplit. The
protein pool differs too -- `(0, 0.1037)` grams against `(-125, 0)` milligrams.

Three things that made the build non-obvious, recorded so it need not be rediscovered:
the conventional GEM must be the **.yml and not the SBML**, because DLKcat needs substrate
SMILES and `metSmiles` does not survive an SBML import -- the input file is written EMPTY
with no error (yeast-GEM 9.0.2's yml carries 1800 SMILES; GECKO's bundled 8.6.2 yml carries
none, so the tutorial cannot be copied verbatim); **DLKcat is not optional**, since fuzzy
BRENDA matching alone reaches only 0.2602 /h after tuning, worse than the model it would
replace; and the **Optimization Toolbox is not needed**, because RAVEN bundles GLPK and
libSBML and RAVEN requires only one LP solver.

**What was originally recorded here, and what remains true of it.** Three findings, checked
2026-09-03 before the build was attempted:

  * The **published** `ecYeastGEM_batch` in `SysBioChalmers/ecModels` still carries the model
    id `M_ecYeastGEM_batch_v8__46__3__46__4` -- byte-identical in version to the copy vendored
    here. Upstream has not rebuilt on any yeast-GEM 9.x.
  * `ecModels/config.ini` pins the build chain to **MATLAB R2020b**, COBRA Toolbox v3.1,
    RAVEN v2.4.1 and the libSBML MATLAB interface. RAVEN describes itself as "a software suite
    for MATLAB"; neither it nor GECKO claims Octave support.
  * **Octave 11.3.0 is installed here and does not help**, because it is not in that chain.
    Installing it was worth doing and it settles the question in the other direction: the
    blocker is MATLAB-only tooling, not a missing interpreter.

A rebuild needs a MATLAB licence, and that is what changed: one was available, so the build
above was run under MATLAB R2026a. The first two findings stand -- nobody publishes an
ecYeastGEM on 9.x, and the chain is MATLAB-only. The third is now moot: Octave is still not
usable, but it is no longer the blocker.

Verify a copy against the table above with:

```bash
gunzip -c data/gem/yeast-GEM.xml.gz | shasum -a 256
```
