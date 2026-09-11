"""Which enzyme sets the ceiling, and the refusal that replaces a coefficient nobody measured.

The model's sharpest claim is that no genotype at any growth rate exceeds 1.2483 mg/gDCW of
beta-carotene, because for a terminal node fed by a saturating step `mu` cancels and content
is capped at `vmax_per_growth`. Arhar 2024 (PMID 39215465) measured 79 mg/gDCW. The claim is
refuted about 63-fold.

What makes it a STRUCTURAL failure rather than a wrong constant: the one genotype input the
model accepted, `entry_expression`, feeds the flux INTO the chain, while the ceiling is set by
the saturating step at the END of it. Raising the input a hundredfold moves the answer by
0.06%. The model could not express the question.
"""

from __future__ import annotations

import pytest

from ystwin.kinetic.carotenoid import ELIZONDO2025
from ystwin.pathway.capacity import (
    CAPACITY_SETTING_ENZYME,
    CapacityUnmeasured,
    capacity_for,
)
from ystwin.predict import Genotype

ANCHOR = ELIZONDO2025.capacity_mmol_per_gdcw


class TestTheGenotypeCanNowExpressACassette:
    def test_a_cassette_is_optional(self):
        """`None` means "not stated", which is not the same as one copy of everything."""
        assert Genotype(1.0, "x").cassette is None

    def test_and_is_carried_when_given(self):
        genotype = Genotype(1.0, "x", cassette={"crtE": 1.0, "crtYB": 2.0})

        assert genotype.cassette["crtYB"] == 2.0

    def test_a_zero_dosage_is_refused(self):
        """A gene the strain does not carry is absent from the mapping, not zero."""
        with pytest.raises(ValueError, match="must be positive"):
            Genotype(1.0, "x", cassette={"crtYB": 0.0})


class TestTheCeilingIsSetByADifferentEnzymeThanTheFluxLaw:
    def test_the_two_are_named_separately(self):
        """`flux.py`'s alpha is fitted on crtE; the capacity is the crtYB cyclase step.
        Conflating them is why raising the model's only input does nothing."""
        assert CAPACITY_SETTING_ENZYME["beta_carotene"] == "crtYB"

    def test_raising_the_flux_enzyme_does_not_move_the_ceiling(self):
        """crtE at ten times dosage returns the anchor unchanged, and that is correct."""
        assert capacity_for("beta_carotene", ANCHOR,
                            {"crtE": 10.0, "crtYB": 1.0}) == ANCHOR

    def test_the_anchor_is_the_refuted_ceiling(self):
        """1.2483 mg/gDCW, against Arhar 2024's measured 79."""
        assert ANCHOR * 536.87 == pytest.approx(1.2483, abs=1e-3)


class TestItRefusesRatherThanExtrapolating:
    def test_a_different_capacity_gene_dosage_is_refused(self):
        with pytest.raises(CapacityUnmeasured):
            capacity_for("beta_carotene", ANCHOR, {"crtYB": 4.0})

    def test_the_refusal_names_the_experiment(self):
        with pytest.raises(CapacityUnmeasured) as raised:
            capacity_for("beta_carotene", ANCHOR, {"crtYB": 4.0})
        message = str(raised.value)

        assert "WHAT WOULD SETTLE IT" in message
        assert "dosage ALONE" in message
        assert "39215465" in message

    def test_it_is_not_a_value_error(self):
        """The caller did nothing wrong. The question is well formed and unanswerable, and
        that must not be caught by the same `except` as a bad argument."""
        with pytest.raises(CapacityUnmeasured) as raised:
            capacity_for("beta_carotene", ANCHOR, {"crtYB": 2.0})

        assert not isinstance(raised.value, ValueError)

    def test_the_calibration_strain_itself_is_answerable(self):
        assert capacity_for("beta_carotene", ANCHOR, {"crtYB": 1.0}) == ANCHOR

    def test_an_unstated_cassette_is_answerable(self):
        """Every existing caller passes no cassette and must keep working unchanged."""
        assert capacity_for("beta_carotene", ANCHOR, None) == ANCHOR

    def test_a_product_with_no_known_capacity_enzyme_is_answerable(self):
        """Only beta-carotene has been diagnosed. Refusing everything else would be a guess
        in the other direction."""
        assert capacity_for("phb", ANCHOR, {"phaA": 4.0}) == ANCHOR


class TestWhyTheCoefficientCannotBeFittedFromTheLiterature:
    """The condition the refusal actually rests on, asserted against the survey table.

    REPLACES `test_the_survey_reports_cassettes_qualitatively`, which read a MARKDOWN file
    and asserted that the substrings "high-copy" and "integrated" appeared somewhere in it.
    That passes whatever the strains are, and it went on passing after
    `data/carotenoid_batch/published_batch_titres.tsv` gained ``copies_crtE``,
    ``copies_crtYB``, ``copies_crtI`` and ``copy_evidence`` columns -- which made the claim it
    was guarding ("28 published strains and not one states its dosage as a number") false
    while the test stayed green.

    The refusal rests on CONFOUNDING, not on absent numbers: too few distinct comparable
    dosages to fit anything. This asserts that, so it FAILS the day a third arrives -- which
    is exactly when the refusal should be revisited.
    """

    def test_fewer_than_three_distinct_numeric_crtyb_dosages_are_published(self):
        from ystwin.pathway.published_cassettes import numeric_crtyb_dosages

        dosages = numeric_crtyb_dosages()

        assert len(dosages) < 3, (
            f"{len(dosages)} distinct crtYB dosages are now published ({sorted(dosages)}); "
            f"the capacity refusal in pathway/capacity.py rests on there being too few to "
            f"fit, so it is due for revisiting rather than for a looser test")

    def test_some_papers_DO_state_a_number_which_is_why_the_old_reason_was_wrong(self):
        """The half that falsified the previous claim. Recorded as an assertion so that
        "every cassette is qualitative" cannot be written again."""
        from ystwin.pathway.published_cassettes import (
            numeric_crtyb_dosages,
            published_cassettes,
        )

        stating = [c for c in published_cassettes() if c.states_a_numeric_crtyb_dosage]

        assert len(stating) >= 4
        assert all(c.copy_evidence for c in stating), "a stated dosage needs named evidence"
        assert numeric_crtyb_dosages()

    def test_the_refusal_message_quotes_the_derived_reason(self):
        """`capacity.py` must not carry its own copy of the reason. That is how the last one
        went stale three days after it was written."""
        from ystwin.pathway.published_cassettes import crtyb_dosage_unfittable_reason

        with pytest.raises(CapacityUnmeasured) as raised:
            capacity_for("beta_carotene", ANCHOR, {"crtYB": 2.0})

        assert crtyb_dosage_unfittable_reason() in str(raised.value)
