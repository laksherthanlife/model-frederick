"""An A453 assay reads the sum of the coloured carotenoids, and the sum breaks the ceiling.

The model asserts `content <= 1.2483 mg/gDCW` for **beta-carotene**, and mu cancels exactly.
Someone will try to test that with a plate reader, because a plate reader is what a teaching
lab has and because a beta-carotene standard makes the assay feel specific. It is not.

Single-wavelength absorbance of a crude hexane extract sums every COLOURED carotenoid --
beta-carotene and **lycopene**, plus gamma-carotene and neurosporene. Phytoene is colourless
(lambda-max about 286 nm) and is the one intermediate that does not interfere. A pure
standard fixes the extinction coefficient; it confers no specificity.

That matters here because **lycopene is the direct precursor and the model predicts it
accumulating past the product**. So the sum can sit above the ceiling while the ceiling is
being obeyed -- which is not a hypothetical: it is exactly what Lopez 2019 (PMID 31380362)
reported, and this repository then recorded their absorbance sum as analyte `beta-car` and
called it a 17x refutation.

These tests pin the arithmetic behind `MEASUREMENTS_NEEDED.md` M5 step 6 and
`EXTERNAL_CAROTENOID_BOUND.md`'s correction, so the tables there cannot go stale in silence.
"""

from __future__ import annotations

import csv
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]

# Both are C40H56, so one molar mass serves. From data/pathways/beta_carotene.toml.
MOLAR_MASS_G_PER_MOL = 536.87

# MEASUREMENTS_NEEDED.md M5: content <= vmax_per_growth, with mu cancelling.
CEILING_MG_PER_GDCW = 1.2483


@pytest.fixture(scope="module")
def registered():
    path = REPO / "outputs" / "registered_prediction_D018.csv"
    if not path.is_file():
        pytest.skip("registered prediction not present; run scripts/register_prediction.py")
    # Read to a string first: handing csv an open file leaves the handle for the garbage
    # collector, and pytest turns that into an unraisable-exception error.
    return list(csv.DictReader(path.read_text().splitlines()))


def _rows(registered):
    for r in registered:
        beta = float(r["predicted_bcar_content_mg_per_gdcw"])
        lycopene = float(r["predicted_lycopene_content"]) * MOLAR_MASS_G_PER_MOL
        yield r["strain"], beta, lycopene, beta + lycopene


def test_beta_carotene_alone_obeys_the_ceiling_in_every_registered_strain(registered):
    """The claim under test. If this fails the ceiling is refuted by the model itself."""
    for strain, beta, _, _ in _rows(registered):
        assert beta <= CEILING_MG_PER_GDCW, (
            f"{strain} predicts {beta:.4f} mg/gDCW of beta-carotene, above the "
            f"{CEILING_MG_PER_GDCW} ceiling the model asserts")


def test_but_the_absorbance_sum_breaks_it_in_two_of_the_three(registered):
    """The point of the file: the assay disagrees with the claim while the claim holds."""
    over = [(s, total) for s, _, _, total in _rows(registered)
            if total > CEILING_MG_PER_GDCW]
    assert len(over) == 2, (
        f"expected exactly two strains whose A453 sum exceeds the ceiling, got "
        f"{[(s, round(t, 4)) for s, t in over]}. The M5 and EXTERNAL_CAROTENOID_BOUND tables "
        "quote two; if this number moved, both are stale")
    assert {s for s, _ in over} == {"b-car3", "b-car4"}


def test_lycopene_exceeds_the_product_in_the_highest_strain(registered):
    """Why the interference is not a rounding concern -- the precursor is the larger pool."""
    by_strain = {s: (beta, lyc) for s, beta, lyc, _ in _rows(registered)}
    beta, lycopene = by_strain["b-car4"]
    assert lycopene > beta, (
        f"b-car4 predicts lycopene {lycopene:.4f} vs beta-carotene {beta:.4f} mg/gDCW. If the "
        "precursor no longer dominates, the absorbance argument weakens and the docs should "
        "say so")


@pytest.mark.parametrize("strain,expected_total", [
    ("b-car2", 0.7762), ("b-car3", 1.7703), ("b-car4", 2.4726),
])
def test_the_published_sums_are_what_the_docs_quote(registered, strain, expected_total):
    """Pins the exact three numbers in M5 step 6 and the EXTERNAL_CAROTENOID_BOUND table."""
    totals = {s: total for s, _, _, total in _rows(registered)}
    assert totals[strain] == pytest.approx(expected_total, abs=5e-4)


def test_phytoene_is_not_part_of_the_interference():
    """The one intermediate that does NOT interfere, and the docs say so.

    Phytoene has three conjugated double bonds and absorbs near 286 nm, so it is invisible at
    453. Pinned as a documentation claim rather than a computation because getting it
    backwards would overstate the problem, and overstating it is its own error.
    """
    for name in ("MEASUREMENTS_NEEDED.md", "EXTERNAL_CAROTENOID_BOUND.md"):
        text = (REPO / "docs" / name).read_text()
        assert "colourless" in text and "286" in text, (
            f"docs/{name} should record that phytoene is colourless (~286 nm) and does not "
            "interfere at 453 nm")


class TestTheCeilingIsRefuted:
    """The claim was the sharpest in the package and it is false. Pinned so it stays said.

    A repository that publishes a falsifiable claim and then quietly stops mentioning its
    falsification has not been doing the thing it says it does. These are documentation
    tests on purpose: the arithmetic is not in dispute, the honesty of the record is.
    """

    def test_the_refutation_is_recorded_with_its_evidence(self):
        text = (REPO / "docs" / "HARD_TESTS.md").read_text()
        assert "39215465" in text, "Arhar 2024 must be cited by PMID"
        assert "79 mg" in text or "79.00" in text or "79 mg/gDCW" in text
        for token in ("gravimetric", "lycopene", "REFUTED"):
            assert token.lower() in text.lower(), f"HARD_TESTS.md should address {token!r}"

    def test_the_prediction_layer_says_so_where_a_caller_will_see_it(self):
        """The near-ceiling note must not still claim the bound is unreachable in reality."""
        source = (REPO / "src" / "ystwin" / "predict.py").read_text()
        assert "39215465" in source, (
            "predict.py's near-ceiling note should name the paper that refuted the ceiling, "
            "so a caller reading the note is not told a false thing about the organism")
        assert "reachable by NO genotype" not in source, (
            "the old note asserted the ceiling is reachable by no genotype at any growth "
            "rate. Arhar 2024 reached 63x it")

    def test_the_withdrawn_paper_is_not_cited_as_evidence(self):
        """Bubphasawan 2024 (PMID 38710418) was withdrawn by Elsevier; Crossref says so.

        **Scoped to one file's table rows until 2026-08-30, and that is how it stayed green
        while `docs/FINDINGS.md` cited the retracted contents as live evidence in prose for
        two days.** A guard against a citation has to look wherever a citation can appear:
        every Markdown file in the repository, and every line of it, not the rows of the one
        table somebody remembered to strike. Found by the 2026-08-30 audit.
        """
        offenders = []
        for path in sorted([REPO / "README.md", *(REPO / "docs").rglob("*.md")]):
            for number, line in enumerate(path.read_text().splitlines(), 1):
                if "38710418" not in line:
                    continue
                # "retracted" is the publisher's own word for it and was missing from this
                # list, which failed the first line written to record the retraction.
                if any(mark in line.lower() for mark in
                       ("withdrawn", "struck", "retract", "~~")):
                    continue
                offenders.append(f"{path.relative_to(REPO)}:{number}: {line.strip()[:120]}")

        assert not offenders, (
            "the withdrawn PMID 38710418 is cited without being marked withdrawn in:\n  "
            + "\n  ".join(offenders))
