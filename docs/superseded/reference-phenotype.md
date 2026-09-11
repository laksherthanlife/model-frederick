# Superseded: the measured phenotype the models were validated against

**Current values:** van Hoek 1998's own Table 1 at D = 0.40 /h — glucose 11.1, oxygen 3.7,
ethanol 13.9, CO2 18.9 mmol/gDW/h, growth 0.40 /h, strain DS28911, aerobic glucose-limited
chemostat.

The previous values were attributed to the same paper and are not in it. Two published
findings rested on them and both reverse.

---

## What was there

```python
REFERENCE_AEROBIC_BATCH = PhysiologyReference(
    name="aerobic glucose batch",
    growth_rate=0.40, glucose_uptake=21.3, oxygen_uptake=7.8,
    ethanol_secretion=27.4, co2_secretion=20.4,
    source="van Hoek, van Dijken & Pronk 1998, Appl Environ Microbiol 64:4226",
)
```

## What the paper contains

Checked directly against the full text of **PMID 9797269** (PMC106631):

| string | occurrences |
| --- | ---: |
| `21.3` | **0** |
| `27.4` | **0** |
| `20.4` | **0** |
| `CEN.PK113-7D` — the strain `tests/test_physiology.py` claimed | **0** |
| `DS28911` — the paper's actual strain | 10 |
| `11.1` — its glucose uptake at D = 0.40 | 1 |

The paper ran **no batch culture**. Every steady state in it is a glucose-limited
chemostat. The other van Hoek 1998 record, PMID 9603825, does use CEN.PK113-7D but reports
µ_max and enzyme activities only, no extracellular fluxes, and does not contain those
numbers either. Those are the only two van Hoek 1998 records in PubMed.

The old fluxes are close to **twice** the paper's at the same growth rate, which reads
like a genuine glucose-excess batch phenotype from a source nobody has found rather than
transcription drift. They were not preserved: a number that cannot be sourced is not a
measurement, and there is no honest way to keep it while saying so.

## The two findings that reverse

**"Plain Yeast9 overpredicts growth by 73%."** It does not. Given the paper's own uptakes
it grows at 0.346 against a measured 0.40 — a **13.5% underprediction**. Against the
corrected reference both models pass at 35% tolerance:

| | growth | ethanol | CO2 |
| --- | ---: | ---: | ---: |
| plain Yeast9 | −13.5% | +15.3% | +8.1% |
| GECKO ec | −12.9% | +13.2% | +2.7% |

**"The ec model misses oxygen by 65%."** It is 27%. Two artefacts inflated it: the
reference oxygen was 7.8 against a measured 3.7, and `fba/physiology.py::cap_uptake` was a
no-op on the ec model because that model is irreversibly split — `r_1992` was never
actually constrained. Both were fixed the same day.

## What survives, and what it does not license

The ec model's growth agreement is real and is the reason to prefer it: left unconstrained
its protein pool alone puts growth within **5.8%** of the measured rate. But it reaches
that growth on **61% more glucose** and **113% more ethanol** than the culture it is
compared against — a glucose-excess optimum measured against a glucose-limited chemostat.
Right growth, wrong route. Oxygen-, OUR- and respiration-linked claims remain unsupported.

## A second mis-citation, same paper

`bridge/physiology_bridge.py` carried `_GLUCOSE_QMAX = 22.0 # mmol/gDW/h, van Hoek et al.
1998 aerobic batch`. That number **is** in the paper, and it is a different quantity:

> "the fermentative capacity showed only a small further increase, up to **22.0 mmol of
> ethanol** · g of dry yeast biomass⁻¹ · h⁻¹ at D = 0.40 h⁻¹"

An offline assay in which chemostat-grown cells are moved to **anaerobic** conditions
under CO2 with 2% glucose. So an anaerobic ethanol production capacity was read as an
aerobic glucose uptake ceiling — wrong metabolite, wrong gas regime. The highest in-situ
glucose uptake the paper reports is 11.1.

## The lesson

The first error is a citation that was never checked against its source. The second is
worse and more interesting: the number *was* in the cited paper, so any check that
verifies "does this identifier resolve" or even "does this figure appear in that paper"
passes it. Only reading the sentence around the number catches it.

`tests/test_citations.py` says as much in its own docstring — *"whether a paper supports
the sentence next to it is a judgement, and no test makes it"*. This is what that residue
looks like when it bites.
