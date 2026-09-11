"""Which enzyme sets a pathway's capacity, and why this repository cannot yet say by how much.

**The claim this module exists to make honest.** For a terminal node fed by a saturating
step, `pathway/solve.py` gives ``content = flux_out / mu`` with
``flux_out <= vmax_per_growth * mu``, so **mu cancels exactly** and

    content <= vmax_per_growth = 2.3252e-3 mmol/gDCW x 536.87 g/mol = 1.2483 mg/gDCW

for every genotype at every growth rate. Arhar 2024 (PMID 39215465) measured **79 mg/gDCW**
of beta-carotene by HPLC on gravimetric dry weight with lycopene below detection. **The model
is refuted by about 63x**, and this is the sharpest falsifiable claim the package makes.

**Raising the one genotype input the model accepts does not help**, and that is the whole
diagnosis. `predict.py`'s `entry_expression` feeds `pathway/flux.py`'s ``flux = alpha *
expression``, which sets the flux INTO the chain. The ceiling is set by ``vmax_per_growth``,
a property of the saturating step at the END of it -- the crtYB cyclase, not the crtE the
flux law is fitted on. Measured on the shipped calibration:

    entry expression      1      10     100      ceiling
    content mg/gDCW  1.1065  1.2369  1.2476       1.2483

A hundredfold increase buys 0.06%. The input the model exposes is not the input the answer
depends on.

**What is missing is one coefficient, and the reason it is missing was RESTATED on
2026-09-04.** This paragraph read: "every cassette is qualitative -- 'crtI x2', '2 copies
crtY', 'high-copy episomal', 'integrated'. A dosage that cannot be read as a number cannot be
regressed against a content." That is false, and the repository's own survey is what falsifies
it: `data/carotenoid_batch/published_batch_titres.tsv` carries ``copies_crtE``,
``copies_crtYB``, ``copies_crtI`` and a ``copy_evidence`` column, and five rows state a crtYB
dosage as a NUMBER with named evidence -- Verwaal 2007 ("confirmed by Southern blotting"),
Li 2013, Xie 2014 ("Table I genotype names every locus"), Arhar 2024 and Bubphasawan 2025. The
claim went stale three days after it was written, because it was a claim ABOUT A FILE that
nothing recomputed.

The refusal below is UNCHANGED, and the true reason is sharper.
:func:`ystwin.pathway.published_cassettes.crtyb_dosage_unfittable_reason` computes it from the
table at the point of use: those five rows hold only TWO distinct numeric dosages, 1 and 2,
and every pair of them differs in host, promoter, medium, cultivation mode and assay as well
as in dosage. **The coefficient is unfittable through CONFOUNDING, not through absent
numbers** -- and unlike the old reason, this one stops being true the day a third comparable
dosage arrives, which is exactly when the refusal should be revisited.
`tests/test_capacity.py` asserts that condition rather than grepping a document for the words
"high-copy" and "integrated", which is what it did before.

So this module takes the cassette and **refuses**, rather than inventing a scaling that would
make the refuted ceiling move for a reason nobody measured.

That refusal is the point. A model that cannot express the question is worse than one that
expresses it and says the answer is unmeasured, because only the second names an experiment.
"""
from __future__ import annotations

from collections.abc import Mapping

__all__ = [
    "CAPACITY_SETTING_ENZYME",
    "CapacityUnmeasured",
    "capacity_for",
]

#: The enzyme whose dosage sets the ceiling, as distinct from the one the flux law is fitted
#: on. For `beta_carotene` that is the crtYB cyclase: `kinetic/carotenoid.py`'s
#: ``vmax_per_growth`` is a property of the lycopene -> beta-carotene step, and `flux.py`'s
#: ``alpha`` is fitted on crtE. Naming them separately is what makes the confusion visible.
CAPACITY_SETTING_ENZYME = {"beta_carotene": "crtYB"}


class CapacityUnmeasured(RuntimeError):
    """Raised when a caller asks for a capacity at a dosage nothing has measured.

    Deliberately not a ``ValueError``: the caller did nothing wrong. The question is
    well-formed and the repository cannot answer it, which is a different thing from a bad
    argument and should not be caught by the same ``except``.
    """


def capacity_for(product: str, anchor_mmol_per_gdcw: float,
                 cassette: Mapping[str, float] | None) -> float:
    """The capacity for this strain's cassette, or a refusal naming what would settle it.

    Args:
        product: Pathway name, used to look up which enzyme sets the ceiling.
        anchor_mmol_per_gdcw: The measured capacity of the CALIBRATION strain -- for
            `beta_carotene`, ``kinetic/carotenoid.ELIZONDO2025.capacity_mmol_per_gdcw``.
        cassette: Relative gene dosages on the calibration strain's scale, or ``None``.

    Returns:
        The anchor, when the strain is the one the anchor was measured on -- either because
        no cassette was stated, or because the capacity-setting gene is at dosage 1.0.
        Returned UNCHANGED rather than rebuilt, so a
        :class:`~ystwin.pathway.solve.ContentCeiling` handed in comes back still carrying the
        measurement that refutes it. This used to read ``float(anchor_mmol_per_gdcw)``, which
        is the one line in this module that could quietly drop that.

    Raises:
        CapacityUnmeasured: when the capacity-setting gene is at any other dosage. **This is
            the honest answer and not a gap to route around.** Extrapolating would mean
            choosing a functional form -- linear in copy number? saturating? -- that no
            measurement in this repository distinguishes, and that the published survey
            cannot distinguish either: it holds two distinct numeric crtYB dosages, and two
            points cannot separate a line from a curve even before the confounding.
    """
    enzyme = CAPACITY_SETTING_ENZYME.get(product)
    if enzyme is None or not cassette or enzyme not in cassette:
        return _anchor(anchor_mmol_per_gdcw)

    dosage = float(cassette[enzyme])
    if dosage == 1.0:
        return _anchor(anchor_mmol_per_gdcw)

    raise CapacityUnmeasured(
        f"{product}: this strain carries {enzyme} at dosage {dosage:g} against the "
        f"calibration strain's 1.0, and NOTHING MEASURED how capacity scales with it. The "
        f"anchor {anchor_mmol_per_gdcw:.4g} mmol/gDCW is the capacity OF THE CALIBRATION "
        f"STRAIN, and returning it here would report a ceiling for a genotype it was not "
        f"measured on.\n\n"
        f"This matters more than an ordinary gap: the ceiling it sets, 1.2483 mg/gDCW for "
        f"beta-carotene, is refuted about 63-fold by Arhar 2024 (PMID 39215465, 79 mg/gDCW "
        f"by HPLC). The repair is exactly this coefficient.\n\n"
        f"WHAT WOULD SETTLE IT: content measured on three or more strains differing in "
        f"{enzyme} dosage ALONE, with the dosage stated as a number, and with the same "
        f"analyte and the same assay. Why the published literature does not already settle "
        f"it, computed from data/carotenoid_batch/published_batch_titres.tsv rather than "
        f"asserted: {_literature_reason()}.")


def _anchor(anchor_mmol_per_gdcw: float) -> float:
    """Normalise the anchor to a float WITHOUT stripping a marked one.

    ``float(x)`` on a :class:`~ystwin.pathway.solve.ContentCeiling` returns a bare float and
    discards the refutation it carries, which is exactly the loss that class exists to
    prevent -- and this function is the only place in the capacity path where the number
    changes hands. A ``ContentCeiling`` IS a float, so the identity branch covers it without
    this module importing the solver, and an int or a numpy scalar still converts.
    """
    return (anchor_mmol_per_gdcw if isinstance(anchor_mmol_per_gdcw, float)
            else float(anchor_mmol_per_gdcw))


def _literature_reason() -> str:
    """Why the survey cannot fit the coefficient, derived from the survey.

    Imported lazily and degraded to a short sentence if the table is absent, because a
    refusal message must not itself fail. Until 2026-09-04 this message carried its own copy
    of the reason -- "28 published strains and not one reports its cassette quantitatively" --
    which was false by the time anyone read it. Five rows do.
    """
    try:
        from .published_cassettes import crtyb_dosage_unfittable_reason

        return crtyb_dosage_unfittable_reason()
    except Exception as failure:                     # the survey is optional at runtime
        return (f"the published-strain survey could not be read ({failure}), so the reason "
                f"cannot be computed here; see data/carotenoid_batch/SOURCE.md")
