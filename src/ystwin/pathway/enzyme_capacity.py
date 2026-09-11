"""Capacity from an enzyme instead of from a fit: ``vmax = kcat * [E]``.

The pathway solver's ceiling, ``NodeKinetics.vmax_per_growth``, is FITTED on one product's
own steady states. That is why nothing here generalises: a new product needs a new
calibration before it can be answered at all. This module computes the same quantity from
two things that are properties of a PROTEIN rather than of a product --

    vmax = kcat * [E]        kcat  turnover number, 1/s, a property of the enzyme
                             [E]   enzyme content, mmol/gDCW, a property of the strain

-- and **refuses to let the answer be used**. It is an AUDIT, in the same sense as
:mod:`ystwin.fba.audit` and :mod:`ystwin.pathway.thermo_gate`: it can contradict the fitted
capacity and it can never replace it. What promotes it is stated at the bottom of this
docstring and is not a matter of taste.

WHY AN AUDIT AND NOT A REPLACEMENT. Two independent reasons, and the second is the one that
would bite a caller.

1. It is not identified. The dominant term is the cassette's proteome mass fraction, which
   is not established by the published dosage metadata. The corpus contains both numeric
   copy counts and qualitative expression constructs; neither measures active enzyme mass.
   Over the plausible window the derived capacity moves by the full width of that window,
   so that spread is a diagnostic, not an identified parameter.
2. Promoting it would make the headline number settable by a caller. The fitted capacity is
   scored leave-one-strain-out against six measured states; a derived one is scored against
   nothing. A path where a caller supplies a kcat and the product content moves is exactly
   what `tests/test_everything_is_wired.py` exists to prevent, and its docstring says so:
   "This is what keeps an unvalidated branch from laundering itself into a headline number."

WHAT IT IS FOR, and it has already earned its place once. The shipped beta-carotene ceiling
is 1.2483 mg/gDCW and **Arhar 2024 (PMID 39215465) measured 79 mg/gDCW**, so the fitted
ceiling is low by about 63x -- see :mod:`ystwin.pathway.capacity`. That refutation arrived
from a paper. This module reaches it from the enzyme side instead, and the claim is
DIRECTIONAL over a two-parameter box rather than a value at any point in it.

THE CLAIM, STATED OVER ITS WHOLE BOX. This is a correction of 2026-09-04. The sentence here
used to read "at ANY turnover in the host model's own isoprenoid range the derived capacity
lands on the measurement" -- quantifying over kcat while holding the enzyme content at
whatever the tests happened to pass, which was a single hard-coded 0.5% of protein. The
second axis is the one this module's own next paragraph calls DOMINANT, so a claim
quantified over the first while the second sat pinned at one corner was the same shape of
error commit cf34701 had already corrected once for a supply parameter. The box has two
axes and both are declared:

    kcat            0.011 to 29900 /s   :data:`ISOPRENOID_KCAT_RANGE_PER_S`, a min/max over
                                        :data:`ISOPRENOID_BRANCH_GENES`
    [E], as a MASS  0.1% to 5% of       the window a heterologous cassette enzyme plausibly
    fraction        total protein       occupies, from a weak single integration up to a
                                        strong high-copy one. UNMEASURED for every strain
                                        here; the promotion criterion at the bottom of this
                                        docstring is the experiment that would close it.

Over that whole box, at mu = 0.101 and CrtYB's 74736 g/mol, the derived capacity is between
1.13x and 1.5e8x the fitted 2.3252e-3 mmol/gDCW -- **above it everywhere, at every corner,
including the most conservative one**. That is the property the tests pin and it is all this
module may assert. What it may NOT assert is a magnitude: the same box spans 0.018x to
2.4e6x Arhar's 0.14715 mmol/gDCW, so it BRACKETS the measurement rather than landing on it.
**A first-principles check that would have flagged the ceiling before the paper did** is
worth having even though it cannot set the number.

THE PRIOR ON kcat IS TAKEN FROM THE HOST MODEL, NOT INVENTED. This is a correction of a
mistake made on 2026-09-02 and it is recorded because the mistake was invisible. The first
pass at this asserted "carotenoid enzymes are slow, kcat ~ 0.01-0.1 /s" from intuition, and
every downstream conclusion inherited it. The vendored `data/gem/ecYeastGEM_batch.xml.gz`
carries a turnover number for 90.6% of its GPR reactions, encoded as 1/kcat coefficients on
the protein pool, and its NATIVE ISOPRENOID BRANCH -- the closest measurable neighbours of
CrtE, CrtI and CrtYB -- puts ERG1, the squalene epoxidase, at 0.076 /s. The asserted prior
was low by one to three orders of magnitude. Anything this module says about a carotenoid
enzyme is anchored on that branch and on nothing else, and a caller who supplies their own
kcat is supplying a measurement, not a guess.

AND THE CORRECTION HAD A SECOND HALF, made on 2026-09-04. The interval that came out of that
first repair, ``(0.076, 46.0)``, was itself a derived number frozen as a literal: it was two
NAMED enzymes -- ERG1 and HMG1 -- not a min/max, and the gene set and direction rule that
would have made it one were never written down. Over the branch as this module now declares
it, :data:`ISOPRENOID_BRANCH_GENES`, the true interval is 0.011 to 29900 /s: ERG6 is slower
than ERG1 and ERG10 and IDI1, both core MVA steps, are ~20x and ~650x faster than 46. The
constant is now a min/max over a DECLARED gene set under a DECLARED direction rule, and
`tests/test_enzyme_capacity.py` reopens the vendored file and recomputes it rather than
retyping it. See the constant's own comment for what that width costs the claim.

THIS MODULE NO LONGER OWNS A ppm -> mmol/gDCW CONVERSION, and that is a correction of
2026-09-04. Until then `derived_capacity` took a `proteome_fraction` and computed
``total_protein_g / fraction / MW_enzyme`` -- a MASS map -- while `pathway/proteome.py`
computed ``ppm/1e6 x total_protein_moles`` -- a MOLAR map, which is the correct one for
PaxDb's molar ppm. Both were reachable, `proteome.proteome_fraction`'s docstring advertised
the molar number as being "in the form `enzyme_capacity` takes", and nothing in either
signature said which convention it carried. Feeding a molar fraction down the mass path is
wrong by exactly ``AVERAGE_PROTEIN_MW / MW_enzyme``: 0.6690x for CrtYB at 74736 g/mol, a
1.49x understatement of [E] -- comparable to the whole 1.75x width of the fitted capacity's
own CI. `TOTAL_PROTEIN_G_PER_GDCW` was also defined independently in both modules.

The repair was to remove the ambiguous argument rather than to pick a side for it.
:func:`derived_capacity` now takes the physical quantity ``vmax = kcat * [E]`` actually
needs -- ``enzyme_mmol_per_gdcw`` -- and :mod:`ystwin.pathway.proteome` owns both maps, each
named for its own convention: :func:`~ystwin.pathway.proteome.enzyme_content_mmol_per_gdcw`
for a molar ppm and :func:`~ystwin.pathway.proteome.enzyme_content_from_mass_fraction` for
the mass fraction a targeted PRM assay reports. A caller now states its convention by
choosing a function name, which a docstring cannot be relied on to do.

WHAT WOULD PROMOTE IT to setting the number, stated in advance so it cannot be argued
afterwards. Two things, both measurements:

    1. The cassette's proteome mass fraction, per strain. Targeted PRM/SRM with a heavy
       peptide per Crt enzyme on the three Elizondo strains at both dilution rates -- six
       samples. That collapses the dominant uncertainty from a tenfold window to an interval.
    2. A held-out product content the derived capacity predicts BETTER than the fitted one.
       Not "as well": better, on a state neither was fitted on.

Until both, this module reports and refuses.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "ISOPRENOID_BRANCH_GENES",
    "ISOPRENOID_KCAT_RANGE_PER_S",
    "DerivedCapacity",
    "EnzymeCapacityUnmeasured",
    "derived_capacity",
]

#: The genes that ARE "the host's own isoprenoid branch", declared as data because the
#: interval below is a min/max over them and an unrecorded gene set is what let a wrong
#: interval ship. This is the ergosterol pathway end to end -- acetyl-CoA to ergosterol --
#: which is the branch CrtE, CrtI and CrtYB tap into and compete with for the same precursor
#: pool: the MVA arm (ERG10, ERG13, HMG1, HMG2, ERG12, ERG8, MVD1, IDI1), the prenyl
#: transferases whose products the cassette diverts (ERG20, BTS1, the CrtE homologue), and
#: the post-squalene sterol arm (ERG9, ERG1, ERG7, ERG11, ERG24-27, ERG6, ERG2, ERG3, ERG5,
#: ERG4). The three enzymes nearest the cassette are in it by structure, not by preference:
#: BTS1 is CrtE's yeast orthologue, ERG9 shares CrtYB's prenyltransferase fold, and ERG1 is
#: the FAD desaturase closest to CrtI.
#:
#: DIRECTION RULE, applied when the interval is recomputed: forward arms only. Reaction ids
#: containing ``_REV`` are excluded, and so are GECKO's two classes of pseudo-reaction,
#: ``draw_prot_*`` (protein draw from the pool) and ``arm_*`` (isoenzyme pooling), neither of
#: which carries a turnover at all.
ISOPRENOID_BRANCH_GENES = (
    "ERG10", "ERG13", "HMG1", "HMG2", "ERG12", "ERG8", "MVD1", "IDI1",
    "ERG20", "BTS1", "ERG9", "ERG1", "ERG7", "ERG11", "ERG24", "ERG25",
    "ERG26", "ERG27", "ERG6", "ERG2", "ERG3", "ERG5", "ERG4",
)

#: Turnover range of :data:`ISOPRENOID_BRANCH_GENES`, 1/s, read out of
#: `data/gem/ecYeastGEM_batch.xml.gz` and pinned by a test that reopens the file and
#: recomputes it -- see `tests/test_enzyme_capacity.py`. ERG6 (0.0109998 /s) is the slow end
#: and IDI1's forward arm (29900.086 /s) the fast one. Used ONLY to refuse a kcat outside
#: anything the host's own branch does -- never as a default, because a default here would
#: be the Tier 0 constant this whole module exists to avoid.
#:
#: The endpoints are the branch's own extremes ROUNDED OUTWARD to six significant figures,
#: which is about the precision the file's stoichiometries carry. Outward and not to-nearest
#: because a to-nearest round puts the two enzymes that DEFINE the interval outside it: the
#: first attempt shipped 0.011 and 29900, and the layer in `predict.py` then reported IDI1 --
#: a core MVA step of this very branch -- as "an enzyme unlike the ones the vendored ecModel
#: carries". The test asserts both properties: the endpoints match the recomputed min and max
#: to three significant figures, and every branch turnover lies inside them.
#:
#: CORRECTED 2026-09-04. This shipped as ``(0.076, 46.0)`` with the comment "ERG1 is the slow
#: end and is itself the slowest enzyme in that branch". Both halves were false. The pair was
#: two named enzymes -- ERG1 and HMG1/HMG2 -- and not a min/max over any declared set: five
#: branch enzymes are slower than ERG1 (ERG6 0.011, ERG5 0.015, ERG3 0.015, IDI1's reverse
#: arm 0.061, ERG24 0.065) and five are faster than 46 (ERG10 900, ERG25 924, ERG26 11920,
#: ERG2 16700, IDI1 29900), two of them -- ERG10 and IDI1 -- core MVA steps. The stated upper
#: bound was ~650x low. Nothing could catch it: the gene set lived in a sentence, the audit
#: ledger exempted the constant from marker-checking on the ground that "a vendored model
#: file is not an output table", and the only test asserted `0 < low < high` and
#: `high/low > 100`, which almost any pair satisfies.
#:
#: THE WIDTH IS ITSELF THE FINDING, and it weakens what this constant can do. A single
#: metabolic branch whose ecModel turnovers span 2.7 MILLION fold is not a prior; the range
#: rejects a typo and nothing finer. Read the module's disagreement with the fitted ceiling
#: at ERG1's own 0.076 /s -- a named enzyme -- rather than at "any kcat in the range".
ISOPRENOID_KCAT_RANGE_PER_S = (0.0109997, 29900.1)


class EnzymeCapacityUnmeasured(RuntimeError):
    """The derived capacity needs a number nobody has measured.

    Raised rather than defaulted, for the reason the module docstring gives: the enzyme
    content moves the answer by the width of its own window, so a default would not be a
    conservative choice, it would be the whole result.
    """


@dataclass(frozen=True)
class DerivedCapacity:
    """``vmax = kcat * [E]`` for one enzyme, with everything it rests on attached."""

    enzyme: str
    kcat_per_s: float
    enzyme_mmol_per_gdcw: float
    vmax_mmol_per_gdcw_h: float

    def capacity_mmol_per_gdcw(self, growth_rate: float) -> float:
        """The quantity `NodeKinetics.vmax_per_growth` holds, for comparison with it.

        The fitted parameter is a capacity PER UNIT GROWTH, so the derived absolute vmax is
        divided by mu to be comparable. That the fitted one carries a mu and this one does
        not is a real difference between the two, not a units nuisance: `vmax_per_growth`
        asserts the capacity scales with growth, and an enzyme's kcat does not know about
        growth at all.
        """
        if growth_rate <= 0:
            raise ValueError(f"growth_rate must be positive, got {growth_rate}")
        return self.vmax_mmol_per_gdcw_h / growth_rate

    def disagreement_with(self, fitted_capacity: float, growth_rate: float) -> float:
        """Derived over fitted, as a fold. 1.0 is agreement; large is the interesting case."""
        return self.capacity_mmol_per_gdcw(growth_rate) / fitted_capacity

    def summary(self, fitted_capacity: float | None = None,
                growth_rate: float | None = None) -> str:
        head = (f"{self.enzyme}: kcat {self.kcat_per_s:.3g}/s at "
                f"{self.enzyme_mmol_per_gdcw:.4g} mmol/gDCW of enzyme gives vmax "
                f"{self.vmax_mmol_per_gdcw_h:.4g} mmol/gDCW/h")
        if fitted_capacity is None or growth_rate is None:
            return head
        fold = self.disagreement_with(fitted_capacity, growth_rate)
        direction = "above" if fold >= 1.0 else "below"
        return (f"{head}; capacity {self.capacity_mmol_per_gdcw(growth_rate):.4g} vs fitted "
                f"{fitted_capacity:.4g} mmol/gDCW, {fold:.3g}x {direction} the fit")


def derived_capacity(enzyme: str,
                     kcat_per_s: float | None,
                     enzyme_mmol_per_gdcw: float | None,
                     ) -> DerivedCapacity:
    """``vmax = kcat * [E]``, or a refusal naming the measurement that is missing.

    Args:
        enzyme: Gene name, for the message. Not looked up anywhere -- this module holds no
            table of enzymes, deliberately, because a table would become a set of defaults.
        kcat_per_s: Turnover number, 1/s. **Required.** No default: see the module docstring
            on the 2026-09-02 mistake, where an asserted turnover was wrong by up to three
            orders of magnitude and silently carried every conclusion with it.
        enzyme_mmol_per_gdcw: The enzyme's content, mmol/gDCW. **Required**, and it is the
            term that dominates the answer. This module does NOT convert an abundance into
            a content -- :mod:`ystwin.pathway.proteome` owns that map, in two functions
            named for their conventions, because the 2026-09-04 correction in the module
            docstring is what happens when a fraction's convention lives in prose.

    Raises:
        EnzymeCapacityUnmeasured: when any required input is absent, naming it.
        ValueError: on a non-finite or out-of-range value.
    """
    import math

    missing = [name for name, value in (("kcat_per_s", kcat_per_s),
                                        ("enzyme_mmol_per_gdcw", enzyme_mmol_per_gdcw))
               if value is None]
    if missing:
        raise EnzymeCapacityUnmeasured(
            f"cannot derive a capacity for {enzyme!r}: {', '.join(missing)} not supplied. "
            f"This module holds no table of enzyme defaults, because a default here would "
            f"BE the result -- the enzyme content alone moves the derived capacity by the "
            f"width of its own plausible window. kcat comes from a turnover measurement or "
            f"from a sequence-based predictor with its own error stated; "
            f"enzyme_mmol_per_gdcw comes from targeted proteomics on the strain being "
            f"predicted, which is one heavy peptide per enzyme, converted by "
            f"pathway/proteome.py's enzyme_content_from_mass_fraction. Neither is a number "
            f"this package may invent")

    for name, value in (("kcat_per_s", kcat_per_s),
                        ("enzyme_mmol_per_gdcw", enzyme_mmol_per_gdcw)):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive, got {value}")

    vmax = kcat_per_s * 3600.0 * enzyme_mmol_per_gdcw
    return DerivedCapacity(enzyme=enzyme, kcat_per_s=float(kcat_per_s),
                           enzyme_mmol_per_gdcw=float(enzyme_mmol_per_gdcw),
                           vmax_mmol_per_gdcw_h=vmax)
