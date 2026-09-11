"""Within-paper matched contrasts against the flat genotype-independent beta-carotene ceiling.

    python3 scripts/matched_contrast.py [--out outputs/matched_contrast.csv]

**What this is, named precisely, because the name is the whole honesty question.** A
within-paper before/after is a MATCHED CONTRAST, not a held-out prediction. The model issues
no number for either arm of either contrast: `predict.py` takes `entry_expression`, the
published table states no expression for any of its rows, and no genotype-to-content law
exists in this package. Calling either of these a "prediction" would be a lie about which
direction the arrow points.

**What a matched contrast CAN test, and why it is worth running.** `pathway/solve.py`
asserts one thing that needs no genotype input at all:

    content <= vmax_per_growth,  because mu cancels exactly

for EVERY genotype at EVERY growth rate. That is the one model output a table with no
expression column can still falsify, and it is why this script exists. The refutation on the
books -- Arhar 2024 at 79 mg/gDCW against 1.2483 -- is a BETWEEN-paper comparison, and every
between-paper comparison carries the confounds `capacity.py` names: host, promoter, medium,
cultivation, assay and biomass basis all move together with the cassette. A within-paper
contrast removes all of them at once. It is the design `capacity.py`'s refusal says the
literature lacks, and the literature turns out to hold two of them.

**THE TWO CONTRASTS, and they are not equivalent.**

*Verwaal 2007* moves ONE crtI copy, Southern-confirmed, in one host, one TDH3p, one
YNB/2% glucose, one 72 h 225 rpm shake flask, one HPLC run reported in one table: 501 ->
5,918 ug/gDCW, 11.81x. The control sits INSIDE the ceiling at 0.40x and the derivative sits
4.74x OUTSIDE it. The contrast STRADDLES the bound, so no appeal to a different assay class,
a different biomass basis or a different host is available -- both numbers are the same
assay, the same balance and the same strain background.

And the model's answer here is not merely unsupported, it is a NUMBER, and the number is
wrong. `data/pathways/beta_carotene.toml` declares the phytoene node -- the node CrtI drains
-- as `passthrough`: flux in equals flux out, no capacity, no parameter. `BETA_CAROTENE_KINETICS`
holds kinetics for the lycopene node ALONE. So the shipped model predicts a fold of exactly
1.000 for any change in crtI dosage, at any dosage, and Verwaal measured 11.81. The spec's
own header cites Verwaal 2007 for its chemistry; its own citation refutes its lumping.

*Ukibe 2009* moves BTS1, the *S. cerevisiae* GGPP synthase, in one INVSc1 host, one ADH1p,
one SC/2% glucose, one 20 mL-in-100 mL flask at 20 C and 250 rpm, one HPLC: ~18 -> ~390
ug/gDCW, 21.7x, against the paper's own printed "a prominent 22-fold increase". **Both arms
sit inside the ceiling** -- 0.014x and 0.31x -- so this contrast does not test the bound at
all, and this script says so rather than counting it. BTS1 supplies the chain's precursor,
so its route into the model is the entry flux, which is the one input the model exposes and
the table does not state. Ukibe is a null for the ceiling and an unmet input requirement for
the flux law, and reporting it as a second refutation would be double-counting.

**WHAT NEITHER CONTRAST CAN ESTABLISH.** Neither is a held-out prediction. Neither carries a
dispersion: Verwaal's Table 4 is "averages of two independent cultures, except for
YB/I/E+tHMG1+I, which are triplicates" with no standard deviations printed, and Ukibe's
Figure 2b shows no error bars, so no interval attaches to either fold. In BOTH papers the
derivative arm also gains a vector and its selection marker -- YIplac128 (LEU2) over a
leu2,3-112 host in Verwaal, pHV-BTS1 (HIS3) in Ukibe -- so each pair differs in one
complemented auxotrophy as well as in the cassette. That is a real confound, it is stated on
the contrast rather than in a footnote, and it is why `also_changed` is a field.

**PROVENANCE, and the label this run introduces.** The control arms are not in any vendored
file. They were transcribed HERE from the open-access full text on the date each contrast
records, and they are OUR artifact, not a checksummed upstream input: `Arm.provenance` says
so on every one of them. What makes them checkable rather than asserted is that each
contrast also declares which of its two arms the vendored survey holds, and
:func:`survey_agreement` binds the contrast to
`data/carotenoid_batch/published_batch_titres.tsv` by PMID and refuses to run if the two
transcriptions disagree. A silent edit to the survey breaks this script instead of moving a
result.

**THE THIRD CONTRAST THIS SCRIPT DECLINES.** Xie 2014's row records an internal four-rung
copy ladder in `copy_evidence`. It is declined, with the reason computed from the row rather
than remembered, because the rungs were never transcribed as rows and the row's measurand is
`both_reported` -- so which assay each rung is in is unrecorded, and a ladder scored across
two assay classes is the exact mistake `EXTERNAL_CAROTENOID_BOUND.md` records Lopez 2019
falling into.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
from dataclasses import dataclass

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS
from ystwin.pathway.capacity import (
    CAPACITY_SETTING_ENZYME,
    CapacityUnmeasured,
    capacity_for,
)
from ystwin.pathway.flux import CAROTENOID_GENES, score_product_validation
from ystwin.pathway.published_cassettes import PublishedCassette, published_cassettes
from ystwin.pathway.solve import ContentCeiling, content_ceiling
from ystwin.pathway.spec import PathwaySpec, RateLaw, load_pathway

PRODUCT = "beta_carotene"

#: What an arm's number IS, as a level distinct from a vendored checksummed input. A value
#: transcribed here from a paper is this repository's artifact and must never be recorded as
#: a verified local copy of an upstream file. Every arm below carries it.
TRANSCRIBED_HERE = "transcribed_here_from_open_access_full_text"

#: The survey stores each paper's printed content rounded (Verwaal's 5,918 ug/g as 5.9), so
#: agreement is checked at a relative tolerance rather than by equality.
SURVEY_ROUNDING_TOLERANCE = 0.01

#: A growth rate handed to the nested-validation reproduction below. The published table
#: states none, so this is a FABRICATION and is named one: it exists only to prove that the
#: test still cannot run after being granted the most favourable value available.
GENEROUS_MU_PER_H = 0.1777


class ContrastProvenanceError(ValueError):
    """A declared contrast disagrees with the vendored survey, or names something absent.

    Raised rather than warned. A matched contrast whose derivative no longer matches the
    survey row it claims to come from is a transcription that has drifted from its source,
    and a drifted transcription reported as a result is the failure this whole run exists to
    avoid.
    """


@dataclass(frozen=True)
class Arm:
    """One side of a within-paper contrast, with what the source does and does not report."""

    label: str
    content_mg_per_gdcw: float
    #: Independent cultures behind the value, or ``None`` where the source states none.
    n_cultures: int | None
    #: Whether the source prints a standard deviation or error bar for THIS value. Where it
    #: is False no interval can be attached to the fold, and the report says so.
    dispersion_reported: bool
    provenance: str = TRANSCRIBED_HERE


@dataclass(frozen=True)
class MatchedContrast:
    """A before/after pair from ONE paper, holding everything but the cassette constant."""

    pmid: str
    reference: str
    control: Arm
    derivative: Arm
    #: The gene whose dosage moves. Checked against the spec by :func:`model_response`.
    changed_gene: str
    #: The spec node the changed gene DRAINS, or ``None`` when the gene acts upstream of the
    #: chain's first node and can only reach the answer through the entry flux.
    changed_node: str | None
    changed_step: str
    #: Which arm the vendored survey row holds. Verwaal's row is the derivative; Ukibe's row
    #: is the CONTROL, and getting this backwards would check the wrong number.
    survey_arm: str
    held_constant: tuple[str, ...]
    #: What moves BESIDES the cassette. Empty is a claim, not a default.
    also_changed: tuple[str, ...]
    #: A fold the paper itself prints, for reproduction. ``None`` where it prints none, which
    #: is a fact about the source and not a licence to invent one.
    paper_stated_fold: float | None
    locator: str
    accessed: str

    @property
    def measured_fold(self) -> float:
        return self.derivative.content_mg_per_gdcw / self.control.content_mg_per_gdcw


@dataclass(frozen=True)
class DeclinedContrast:
    """A within-paper contrast this script can see and refuses to score."""

    pmid: str
    reference: str
    what_is_recorded: str


@dataclass(frozen=True)
class ModelResponse:
    """What the shipped model says the contrast's fold is, and by which route."""

    #: The predicted fold, or ``None`` where the model issues no number.
    fold: float | None
    route: str
    reason: str


@dataclass(frozen=True)
class ContrastVerdict:
    contrast: MatchedContrast
    control_fold_of_ceiling: float
    derivative_fold_of_ceiling: float
    measured_fold: float
    model: ModelResponse
    verdict: str
    survey_agreement: str
    limitations: tuple[str, ...]

    @property
    def discrepancy(self) -> float | None:
        """Measured fold over the model's fold, where the model commits to one."""
        return None if self.model.fold is None else self.measured_fold / self.model.fold


CONTRASTS: tuple[MatchedContrast, ...] = (
    MatchedContrast(
        pmid="17496128",
        reference="Verwaal 2007",
        control=Arm("CEN.PK113-6B + YIplac211 YB/I/E* + YIplac204 tHMG1",
                    0.501, n_cultures=2, dispersion_reported=False),
        derivative=Arm("the same, + YIplac128 crtI (one added crtI copy)",
                       5.918, n_cultures=3, dispersion_reported=False),
        changed_gene="CrtI",
        changed_node="phytoene",
        changed_step="phytoene -> lycopene, four desaturations",
        survey_arm="derivative",
        held_constant=("host CEN.PK113-6B", "TDH3p on every crt gene", "YNB + 2% glucose",
                       "72 h shake flask, 30 C, 225 rpm", "one HPLC diode-array method",
                       "one table in one paper"),
        also_changed=("the YIplac128 vector and its LEU2 marker, which complements the "
                      "leu2,3-112 host the control arm retains",),
        paper_stated_fold=None,
        locator="Appl Environ Microbiol 73:4342, Table 4 (501 and 5,918 ug/g dry weight); "
                "single-copy integration confirmed by Southern blotting; PMC1932764",
        accessed="2026-09-09"),
    MatchedContrast(
        pmid="19801484",
        reference="Ukibe 2009",
        control=Arm("INVSc1 + pTV-crtI + pUV-crtYB",
                    0.018, n_cultures=None, dispersion_reported=False),
        derivative=Arm("the same, + pHV-BTS1 (S. cerevisiae GGPP synthase)",
                       0.390, n_cultures=None, dispersion_reported=False),
        changed_gene="BTS1",
        changed_node=None,
        changed_step="FPP -> GGPP, upstream of the chain's first node",
        survey_arm="control",
        held_constant=("host INVSc1", "ADH1p", "SC + 2% glucose",
                       "20 mL in a 100 mL flask, 20 C, 250 rpm, 5 d",
                       "one HPLC method separating beta-carotene at 25.0 min",
                       "one figure in one paper"),
        also_changed=("the pHV-BTS1 plasmid and its HIS3 marker, so the two arms need "
                      "different SC dropout supplementation",),
        paper_stated_fold=22.0,
        locator="Appl Environ Microbiol 75:7205, Figure 2b and Results "
                "('approximately 18 ug per g [dry weight]', 'approximately 390 ug per g "
                "[dry weight]', 'a prominent 22-fold increase'); PMC2786542",
        accessed="2026-09-09"),
)

DECLINED: tuple[DeclinedContrast, ...] = (
    DeclinedContrast(
        pmid="23860829",
        reference="Xie 2014",
        what_is_recorded="a four-rung internal copy ladder in the row's copy_evidence"),
)


def spec() -> PathwaySpec:
    return load_pathway(PRODUCT)


def ceiling() -> ContentCeiling:
    """The bound itself, from the shipped calibration rather than from a literal.

    Every fold in this script divides by this, so a refit moves the whole report instead of
    leaving a stale 1.2483 behind in a comment.
    """
    bound = content_ceiling(spec(), BETA_CAROTENE_KINETICS)
    if bound is None:
        raise ContrastProvenanceError(
            f"{PRODUCT} has no content ceiling under the shipped calibration, so there is "
            f"nothing for a matched contrast to bracket. See pathway/solve.py::content_ceiling")
    return bound


def survey_row(pmid: str) -> PublishedCassette:
    """The vendored survey row for a PMID, or a refusal naming the file."""
    for cassette in published_cassettes():
        if cassette.pmid == pmid:
            return cassette
    raise ContrastProvenanceError(
        f"PMID {pmid} is not in data/carotenoid_batch/published_batch_titres.tsv. A contrast "
        f"that cannot be bound to a survey row has no cross-check on its transcription")


def survey_agreement(contrast: MatchedContrast) -> str:
    """Cross-check the contrast's transcription against the independently transcribed survey.

    The survey holds ONE arm of each pair -- which one is declared on the contrast -- and its
    cell is the paper's printed value rounded. Agreement at
    :data:`SURVEY_ROUNDING_TOLERANCE` is therefore the check, and disagreement raises.
    """
    if contrast.survey_arm not in ("control", "derivative"):
        raise ContrastProvenanceError(
            f"{contrast.reference}: survey_arm is {contrast.survey_arm!r}, expected "
            f"'control' or 'derivative'")
    row = survey_row(contrast.pmid)
    arm = getattr(contrast, contrast.survey_arm)
    if row.content_mg_per_gdcw is None:
        raise ContrastProvenanceError(
            f"{contrast.reference}: the survey row states no content, so the transcription "
            f"of its {contrast.survey_arm} arm has nothing to be checked against")
    relative = abs(row.content_mg_per_gdcw - arm.content_mg_per_gdcw) / arm.content_mg_per_gdcw
    if relative > SURVEY_ROUNDING_TOLERANCE:
        raise ContrastProvenanceError(
            f"{contrast.reference}: the survey holds {row.content_mg_per_gdcw:g} mg/gDCW for "
            f"PMID {contrast.pmid} and this contrast's {contrast.survey_arm} arm is "
            f"{arm.content_mg_per_gdcw:g}, a {relative:.1%} disagreement above the "
            f"{SURVEY_ROUNDING_TOLERANCE:.0%} rounding tolerance. One of the two "
            f"transcriptions is wrong and this script will not report either")
    return (f"survey row {contrast.pmid} holds {row.content_mg_per_gdcw:g} mg/gDCW, this "
            f"contrast's {contrast.survey_arm} arm holds {arm.content_mg_per_gdcw:g} "
            f"({relative:.2%} apart, both from {row.reference})")


def reproduces_paper_fold(contrast: MatchedContrast,
                          tolerance: float = 0.05) -> tuple[bool | None, str]:
    """Does the fold derived from the two arms reproduce the fold the paper itself prints?

    Returns ``(None, reason)`` where the paper prints none. That is the honest outcome for
    Verwaal, whose Table 4 gives both contents and no ratio, and inventing a stated fold to
    reproduce would make this check circular.
    """
    if contrast.paper_stated_fold is None:
        return None, (f"{contrast.reference} prints no fold for this step, so there is "
                      f"nothing to reproduce; the arms themselves are the reproduction "
                      f"target and are cross-checked against the survey")
    relative = abs(contrast.measured_fold - contrast.paper_stated_fold) / contrast.paper_stated_fold
    return relative <= tolerance, (
        f"{contrast.reference} prints {contrast.paper_stated_fold:g}x; the two arms give "
        f"{contrast.measured_fold:.2f}x, {relative:.1%} apart")


def _node(chain: PathwaySpec, name: str):
    for node in chain.nodes:
        if node.name == name:
            return node
    raise ContrastProvenanceError(
        f"{name!r} is not a node of the {chain.product} chain "
        f"({', '.join(node.name for node in chain.nodes)})")


def model_response(contrast: MatchedContrast, chain: PathwaySpec) -> ModelResponse:
    """What the shipped model says this contrast's fold is, derived from the spec.

    Three routes exist by which a cassette change could reach a content, and the answer is
    which one the changed gene takes:

    * the capacity-setting enzyme, where :func:`capacity_for` REFUSES rather than extrapolate;
    * a node the gene drains, where the node's rate law decides -- a `passthrough` node has
      no capacity and no parameter, so the model's fold is exactly 1.000 at any dosage;
    * upstream of the chain, where the only route is the entry flux and the model needs an
      expression the published table does not state.

    Read off the spec rather than looked up per gene, so a spec edit changes the verdict.
    """
    capacity_gene = CAPACITY_SETTING_ENZYME.get(chain.product)
    if capacity_gene is not None and contrast.changed_gene.casefold() == capacity_gene.casefold():
        # The anchor is the feeder node's capacity, read off the chain the way
        # `content_ceiling` reads it, so this does not hard-code which node that is.
        anchor = BETA_CAROTENE_KINETICS[chain.nodes[-2].name].vmax_per_growth
        try:
            capacity_for(chain.product, anchor, {capacity_gene: 2.0})
        except CapacityUnmeasured as refusal:
            return ModelResponse(None, f"capacity ({capacity_gene})", str(refusal))
        return ModelResponse(None, f"capacity ({capacity_gene})",
                             "capacity_for returned a number for a dosage it has not "
                             "measured; that refusal has been weakened and this verdict "
                             "must be re-derived before it is reported")

    named = {node.enzyme.casefold(): node for node in chain.nodes if node.enzyme}
    if contrast.changed_node is None:
        if contrast.changed_gene.casefold() in named:
            raise ContrastProvenanceError(
                f"{contrast.reference}: changed_node is None but {contrast.changed_gene!r} "
                f"drains node {named[contrast.changed_gene.casefold()].name!r}")
        if contrast.changed_gene.casefold() == chain.entry_enzyme.casefold():
            route = f"entry expression ({chain.entry_enzyme})"
        else:
            route = (f"entry flux, upstream of {chain.nodes[0].name!r} "
                     f"(precursor {chain.precursor_metabolite})")
        return ModelResponse(None, route, (
            f"{contrast.changed_gene} acts on {contrast.changed_step}, so its only route "
            f"into this model is the entry flux, and the entry flux reads a relative "
            f"{chain.entry_enzyme} expression that the published table states for no row. "
            f"The model issues NO number for either arm"))

    node = _node(chain, contrast.changed_node)
    if not node.enzyme or node.enzyme.casefold() != contrast.changed_gene.casefold():
        raise ContrastProvenanceError(
            f"{contrast.reference}: node {node.name!r} is drained by {node.enzyme!r}, not by "
            f"the declared {contrast.changed_gene!r}")
    if node.rate_law == RateLaw.PASSTHROUGH:
        return ModelResponse(1.0, f"passthrough node {node.name!r}", (
            f"{contrast.changed_gene} drains the {node.name!r} node, which "
            f"data/pathways/{chain.product}.toml declares `passthrough`: flux in equals flux "
            f"out, with no capacity and no fitted parameter -- BETA_CAROTENE_KINETICS holds "
            f"an entry for {', '.join(sorted(BETA_CAROTENE_KINETICS))} and nothing else. The "
            f"model therefore predicts a fold of exactly 1.000 at ANY {contrast.changed_gene} "
            f"dosage"))
    return ModelResponse(None, f"{node.rate_law} node {node.name!r}", (
        f"{contrast.changed_gene} drains the {node.name!r} node under a "
        f"{node.rate_law} law, whose parameters are fitted on the calibration strain "
        f"and carry no dosage coefficient"))


def limitations(contrast: MatchedContrast, model: ModelResponse) -> tuple[str, ...]:
    """What this contrast cannot show, assembled from its declared fields.

    Written as a computation over the record rather than as prose beside it, so a contrast
    that later gains replicates or loses a marker confound stops carrying the caveat.
    """
    found = [
        "a within-paper before/after is a matched contrast, not a held-out prediction: "
        "nothing was withheld and nothing was forecast",
    ]
    if model.fold is None:
        found.append(f"the model issues no number for either arm ({model.route}), so the "
                     f"measured fold cannot be scored against a prediction")
    silent = [arm for arm in (contrast.control, contrast.derivative)
              if not arm.dispersion_reported]
    if silent:
        cultures = [str(arm.n_cultures) if arm.n_cultures else "unstated"
                    for arm in (contrast.control, contrast.derivative)]
        found.append(f"the source prints no dispersion for "
                     f"{'either arm' if len(silent) == 2 else 'one arm'} "
                     f"(independent cultures: {' and '.join(cultures)}), so no interval "
                     f"attaches to the {contrast.measured_fold:.2f}x")
    for confound in contrast.also_changed:
        found.append(f"the arms are not isogenic: {confound}")
    found.append("one paper is one laboratory; a matched contrast removes between-paper "
                 "confounds and does not replace independent replication")
    return tuple(found)


def verdict_name(control_fold_of_ceiling: float, derivative_fold_of_ceiling: float) -> str:
    """Name the contrast from where its two arms sit relative to the bound.

    Separate from :func:`verdict_for` because the naming rule is the claim: refuting an
    upper bound and STRADDLING it are different statements, and a pair whose control is
    already outside cannot say the added copy carried the strain across.
    """
    if derivative_fold_of_ceiling <= 1.0:
        return "not_a_ceiling_test"
    if control_fold_of_ceiling <= 1.0:
        return "refutes_the_flat_ceiling_within_one_paper"
    return "refutes_the_ceiling_but_both_arms_are_already_outside"


def verdict_for(contrast: MatchedContrast, chain: PathwaySpec,
                bound: ContentCeiling) -> ContrastVerdict:
    """Score one contrast against the ceiling, one-sided.

    The bracket is one-sided by construction: a content ABOVE the bound refutes it, and a
    content below is compatible with the bound and with every smaller value, so it confirms
    nothing. A contrast whose two arms both sit below therefore is not a ceiling test, and
    saying otherwise would be counting a null as evidence.
    """
    control = contrast.control.content_mg_per_gdcw / bound.mg_per_gdcw
    derivative = contrast.derivative.content_mg_per_gdcw / bound.mg_per_gdcw
    model = model_response(contrast, chain)
    return ContrastVerdict(
        contrast=contrast, control_fold_of_ceiling=control,
        derivative_fold_of_ceiling=derivative, measured_fold=contrast.measured_fold,
        model=model, verdict=verdict_name(control, derivative),
        survey_agreement=survey_agreement(contrast),
        limitations=limitations(contrast, model))


def ceiling_bracket(bound: ContentCeiling | None = None) -> pd.DataFrame:
    """The one-sided bracket over every ``dfba_ready`` row of the published survey.

    One row per published strain, with the fold of the ceiling its content reaches. `inside`
    means only "does not refute": the bound is an upper limit, so nothing below it is
    evidence for it.
    """
    bound = bound or ceiling()
    rows = []
    for cassette in published_cassettes():
        if not cassette.dfba_ready or cassette.content_mg_per_gdcw is None:
            continue
        fold = cassette.content_mg_per_gdcw / bound.mg_per_gdcw
        rows.append({
            "pmid": cassette.pmid,
            "reference": cassette.reference,
            "content_mg_per_gdcw": cassette.content_mg_per_gdcw,
            "ceiling_mg_per_gdcw": bound.mg_per_gdcw,
            "fold_of_ceiling": fold,
            "position": "inside" if fold <= 1.0 else "outside",
            "measurand": cassette.measurand,
            "copies_crtE": cassette.copies.get("crtE", ""),
            "copies_crtYB": cassette.copies.get("crtYB", ""),
            "copies_crtI": cassette.copies.get("crtI", ""),
            "promoter": cassette.promoter,
        })
    return pd.DataFrame(rows).sort_values("fold_of_ceiling").reset_index(drop=True)


def dosage_identifiability() -> dict[str, float | int | str]:
    """Whether crtE dosage is identifiable BETWEEN papers on the dfba_ready rows.

    Two numbers, because they answer different questions and the difference matters. Treating
    dosage as a FACTOR gives the between-level variance fraction, which is the largest share
    a dosage effect could possibly explain here. Treating it as a linear covariate gives
    less, because the highest-content strain carries the lowest dosage.
    """
    rows = [row for row in published_cassettes()
            if row.dfba_ready and row.content_mg_per_gdcw is not None
            and isinstance(row.copies.get("crtE"), int)]
    dosage = np.array([float(row.copies["crtE"]) for row in rows])
    log_content = np.log10([row.content_mg_per_gdcw for row in rows])
    design = np.column_stack([np.ones_like(dosage), dosage])
    fitted = design @ np.linalg.lstsq(design, log_content, rcond=None)[0]
    total = float(np.sum((log_content - log_content.mean()) ** 2))
    within = float(sum(np.sum((log_content[dosage == level] - log_content[dosage == level].mean()) ** 2)
                       for level in np.unique(dosage)))
    top = max(rows, key=lambda row: row.content_mg_per_gdcw)
    return {
        "n_rows": len(rows),
        "n_levels": int(np.unique(dosage).size),
        "r2_linear_in_dosage": 1.0 - float(np.sum((log_content - fitted) ** 2)) / total,
        "r2_dosage_as_factor": 1.0 - within / total,
        "within_level_variance_fraction": within / total,
        "highest_content_reference": top.reference,
        "highest_content_mg_per_gdcw": top.content_mg_per_gdcw,
        "highest_content_copies_crtE": top.copies.get("crtE", ""),
    }


def nested_validation_on_published_rows() -> pd.DataFrame:
    """Run the real nested validator on the most generous frame the survey can produce.

    Not a test of the model: a test of the TABLE. Every missing quantity is granted its most
    favourable available value -- one distinct strain per paper, a fabricated mu at
    :data:`GENEROUS_MU_PER_H`, and the content read as a rate -- and the run still fails,
    because the 14 expression channels and the lycopene channel simply do not exist in this
    table and no generosity can conjure them. The failure is the finding.
    """
    rows = [row for row in published_cassettes()
            if row.dfba_ready and row.content_mg_per_gdcw is not None]
    terminal = spec().nodes[-1]
    frame = pd.DataFrame([{
        "condition": row.pmid,
        "strain": row.reference,
        "mu_per_h": GENEROUS_MU_PER_H,
        "q_betacarotene": row.content_mg_per_gdcw / terminal.molar_mass_g_per_mol * GENEROUS_MU_PER_H,
        "q_lycopene": np.nan,
    } for row in rows])
    for gene in CAROTENOID_GENES:
        frame[gene] = np.nan
    return score_product_validation(spec(), frame, draws=0)


def _print_bracket(bracket: pd.DataFrame, bound: ContentCeiling) -> None:
    inside = bracket[bracket.position == "inside"]
    print(f"\nONE-SIDED CEILING BRACKET -- {len(bracket)} dfba_ready published strains "
          f"against {bound.mg_per_gdcw:.4f} mg/gDCW")
    print(f"  the bound is set at {bound.capacity_enzyme} on the {bound.capacity_node!r} node "
          f"and bounds {bound.bounded_node!r}; the flux law reads {spec().entry_enzyme}")
    print(f"\n  {'reference':<20} {'mg/gDCW':>9} {'x ceiling':>10}  {'position':<8} measurand")
    for _, row in bracket.iterrows():
        print(f"  {row.reference:<20} {row.content_mg_per_gdcw:9.3f} {row.fold_of_ceiling:10.3f}"
              f"  {row.position:<8} {row.measurand}")
    print(f"\n  {len(inside)} inside, {len(bracket) - len(inside)} outside, "
          f"worst {bracket.fold_of_ceiling.max():.1f}x "
          f"({bracket.loc[bracket.fold_of_ceiling.idxmax(), 'reference']})")
    print("  'inside' means only that the row does not refute an UPPER bound. Nothing below "
          "a ceiling is evidence for it.")


def _print_contrast(scored: ContrastVerdict) -> None:
    contrast = scored.contrast
    reproduced, note = reproduces_paper_fold(contrast)
    print(f"\n{'=' * 78}\n{contrast.reference} (PMID {contrast.pmid}) -- {contrast.changed_step}")
    print(f"  control    {contrast.control.content_mg_per_gdcw:9.3f} mg/gDCW  "
          f"{scored.control_fold_of_ceiling:7.3f}x ceiling   {contrast.control.label}")
    print(f"  derivative {contrast.derivative.content_mg_per_gdcw:9.3f} mg/gDCW  "
          f"{scored.derivative_fold_of_ceiling:7.3f}x ceiling   {contrast.derivative.label}")
    print(f"  measured fold {scored.measured_fold:.2f}x on {contrast.changed_gene}")
    print(f"  held constant: {'; '.join(contrast.held_constant)}")
    print(f"  provenance: {contrast.control.provenance}, {contrast.locator}, "
          f"read {contrast.accessed}")
    print(f"  cross-check: {scored.survey_agreement}")
    print(f"  paper fold reproduced: {reproduced} -- {note}")
    predicted = "no number" if scored.model.fold is None else f"{scored.model.fold:.3f}x"
    print(f"  MODEL says {predicted} via {scored.model.route}")
    print(f"    {scored.model.reason}")
    if scored.discrepancy is not None:
        print(f"    measured / predicted = {scored.discrepancy:.2f}x")
    print(f"  VERDICT: {scored.verdict}")
    for limitation in scored.limitations:
        print(f"    cannot show: {limitation}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=pathlib.Path, default=None,
                        help="optional CSV of the ceiling bracket and the scored contrasts")
    args = parser.parse_args(argv)

    chain, bound = spec(), ceiling()
    bracket = ceiling_bracket(bound)
    _print_bracket(bracket, bound)

    scored = [verdict_for(contrast, chain, bound) for contrast in CONTRASTS]
    for one in scored:
        _print_contrast(one)

    print(f"\n{'=' * 78}\nDECLINED, with the reason computed from the survey row:")
    for declined in DECLINED:
        row = survey_row(declined.pmid)
        print(f"  {declined.reference}: {declined.what_is_recorded}, but the rungs were never "
              f"transcribed as rows and this row's measurand is {row.measurand!r}, so which "
              f"assay each rung is in is unrecorded. Scoring it would cross assay classes.")

    identifiability = dosage_identifiability()
    print(f"\n{'=' * 78}\nWHY A BETWEEN-PAPER DOSAGE REGRESSION CANNOT REPLACE THESE:")
    print(f"  {identifiability['n_rows']} rows state a numeric crtE dosage over "
          f"{identifiability['n_levels']} levels; log-content variance explained is "
          f"{identifiability['r2_dosage_as_factor']:.3f} with dosage as a factor and "
          f"{identifiability['r2_linear_in_dosage']:.3f} linear in dosage, so "
          f"{identifiability['within_level_variance_fraction']:.1%} sits WITHIN levels")
    print(f"  the highest content, {identifiability['highest_content_reference']} at "
          f"{identifiability['highest_content_mg_per_gdcw']:g} mg/gDCW, carries "
          f"crtE = {identifiability['highest_content_copies_crtE']}")

    validation = nested_validation_on_published_rows()
    failed = validation[validation.status == "failed"]
    print(f"\n{'=' * 78}\nAND WHY THIS IS NOT A LEAVE-ONE-PUBLISHED-STRAIN-OUT TEST:")
    print(f"  score_product_validation on the most generous frame: {len(failed)}/"
          f"{len(validation)} outer folds failed, "
          f"{int(validation.n_selection_candidates.iloc[0])} candidates offered")
    print(f"  reason: {failed.failure_reason.iloc[0]}")
    print(f"  mu was FABRICATED at {GENEROUS_MU_PER_H} /h and the table still cannot run it.")

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        contrasts = pd.DataFrame([{
            "pmid": one.contrast.pmid,
            "reference": one.contrast.reference,
            "changed_gene": one.contrast.changed_gene,
            "control_mg_per_gdcw": one.contrast.control.content_mg_per_gdcw,
            "derivative_mg_per_gdcw": one.contrast.derivative.content_mg_per_gdcw,
            "measured_fold": one.measured_fold,
            "control_fold_of_ceiling": one.control_fold_of_ceiling,
            "derivative_fold_of_ceiling": one.derivative_fold_of_ceiling,
            "model_fold": one.model.fold,
            "model_route": one.model.route,
            "verdict": one.verdict,
            "provenance": one.contrast.control.provenance,
            "locator": one.contrast.locator,
        } for one in scored])
        merged = pd.concat([bracket.assign(record="ceiling_bracket"),
                            contrasts.assign(record="matched_contrast")], ignore_index=True)
        merged.to_csv(args.out, index=False)
        print(f"\n  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
