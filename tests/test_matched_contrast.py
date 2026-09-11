"""The within-paper matched contrasts, and the boundary of what they are allowed to claim.

Two things are pinned here and they pull in opposite directions, which is the point.

**That the test is sharp.** Verwaal 2007 moves one Southern-confirmed crtI copy inside one
paper and the content crosses the ceiling: control 0.401x, derivative 4.741x. The shipped
model's answer for that step is not "unsupported" -- it is a NUMBER, 1.000x, because
`data/pathways/beta_carotene.toml` declares the phytoene node `passthrough` and
`BETA_CAROTENE_KINETICS` holds kinetics for the lycopene node alone. The model is out by
11.81x on a step whose chemistry the spec cites this very paper for.

**That the test is not more than it is.** Ukibe 2009 moves 21.7x with BOTH arms inside the
ceiling, so it must come back `not_a_ceiling_test` and must NOT be counted as a second
refutation. Neither contrast is a held-out prediction, neither carries a dispersion, and in
both papers the derivative arm also gains a selection marker. Every one of those is asserted
below, because a matched contrast reported as a prediction is exactly the mislabelling this
run exists to avoid.

The transcriptions are OUR artifact and never a verified upstream copy. What keeps them
honest is that each contrast declares which arm the vendored survey holds and is refused if
the two transcriptions disagree -- so the tests below feed it a disagreeing survey and check
that it raises rather than reports.
"""

from __future__ import annotations

import dataclasses
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

import matched_contrast as module

from ystwin.pathway.capacity import CAPACITY_SETTING_ENZYME
from ystwin.pathway.published_cassettes import PublishedCassette
from ystwin.pathway.spec import Fate, Node, PathwaySpec, RateLaw

VERWAAL = "17496128"
UKIBE = "19801484"


@pytest.fixture(scope="module")
def chain():
    return module.spec()


@pytest.fixture(scope="module")
def bound():
    return module.ceiling()


@pytest.fixture(scope="module")
def scored(chain, bound):
    return {one.pmid: module.verdict_for(one, chain, bound) for one in module.CONTRASTS}


def _contrast(pmid):
    return next(one for one in module.CONTRASTS if one.pmid == pmid)


class TestTheOneSidedCeilingBracket:
    """The bracket a prior pass reported as 2 inside, 7 outside, up to 63.3x."""

    def test_every_dfba_ready_row_with_a_content_is_bracketed(self):
        table = module.ceiling_bracket()

        assert len(table) == 9
        assert set(table.position) == {"inside", "outside"}

    def test_two_inside_and_seven_outside(self):
        table = module.ceiling_bracket()

        assert (table.position == "inside").sum() == 2
        assert (table.position == "outside").sum() == 7
        assert sorted(table.loc[table.position == "inside", "reference"]) == [
            "Li 2013", "Ukibe 2009"]

    def test_the_worst_row_is_arhar_at_about_63x(self):
        table = module.ceiling_bracket()
        worst = table.loc[table.fold_of_ceiling.idxmax()]

        assert worst.reference == "Arhar 2024"
        assert worst.fold_of_ceiling == pytest.approx(63.3, rel=0.01)

    def test_the_fold_comes_from_the_shipped_calibration_not_a_literal(self, bound):
        """A refit must move this table. `refuted_fold` is computed in `solve.py` from the
        same constant, so the two must agree without either being written down."""
        table = module.ceiling_bracket()
        worst = table.fold_of_ceiling.max()

        assert bound.refuted_fold == pytest.approx(worst, rel=1e-9)
        assert set(table.ceiling_mg_per_gdcw) == {bound.mg_per_gdcw}

    def test_inside_is_not_evidence_for_the_bound(self, bound):
        """Both inside rows are two to three orders below it. A one-sided bracket can only
        refute, and a test that read these as support would be reading a null as a result."""
        table = module.ceiling_bracket()
        inside = table[table.position == "inside"]

        assert (inside.fold_of_ceiling < 0.5).all()


class TestTheTranscriptionsAgreeWithTheVendoredSurvey:
    def test_each_contrast_binds_to_a_survey_row_by_pmid(self):
        for contrast in module.CONTRASTS:
            assert module.survey_row(contrast.pmid).reference == contrast.reference

    def test_the_survey_holds_the_derivative_for_verwaal_and_the_control_for_ukibe(self):
        """Getting this backwards checks the wrong number and would pass silently."""
        assert _contrast(VERWAAL).survey_arm == "derivative"
        assert _contrast(UKIBE).survey_arm == "control"

    def test_agreement_is_within_the_survey_rounding(self, scored):
        """5.9 in the survey against 5,918 ug/g in Table 4 is 0.3% apart, not equal."""
        assert "0.30% apart" in scored[VERWAAL].survey_agreement
        assert "0.00% apart" in scored[UKIBE].survey_agreement

    def test_a_drifted_survey_refuses_rather_than_reports(self, monkeypatch):
        drifted = (PublishedCassette(
            pmid=VERWAAL, reference="Verwaal 2007", strain="s", carbon_source="2% glucose",
            carbon_species=frozenset({"glucose"}), dfba_ready=True, copies={}, promoter="TDH3p",
            integration="genomic", copy_evidence="", content_mg_per_gdcw=59.0,
            measurand="beta_carotene"),)
        monkeypatch.setattr(module, "published_cassettes", lambda: drifted)

        with pytest.raises(module.ContrastProvenanceError, match="disagreement above"):
            module.survey_agreement(_contrast(VERWAAL))

    def test_a_pmid_absent_from_the_survey_refuses(self, monkeypatch):
        monkeypatch.setattr(module, "published_cassettes", tuple)

        with pytest.raises(module.ContrastProvenanceError, match="not in data/carotenoid_batch"):
            module.survey_row(VERWAAL)

    def test_every_arm_is_labelled_as_this_repository_s_transcription(self):
        """Never a verified copy of an upstream file. These numbers were read off a paper
        here, and the level they carry has to say which kind they are."""
        for contrast in module.CONTRASTS:
            for arm in (contrast.control, contrast.derivative):
                assert arm.provenance == module.TRANSCRIBED_HERE
            assert contrast.accessed and "PMC" in contrast.locator


class TestReproducingTheSourcePapersOwnNumber:
    def test_ukibe_reproduces_the_fold_the_paper_itself_prints(self):
        """The paper's own words are 'a prominent 22-fold increase'. Deriving 21.67x from
        its two printed contents is the reproduction, and it is what earns this row a place."""
        reproduced, note = module.reproduces_paper_fold(_contrast(UKIBE))

        assert reproduced is True
        assert "22x" in note and "21.67x" in note

    def test_verwaal_prints_no_fold_and_none_is_invented(self):
        """`None` is the honest outcome, not a failure. Table 4 gives both contents and no
        ratio; manufacturing a stated fold to reproduce would make the check circular."""
        reproduced, note = module.reproduces_paper_fold(_contrast(VERWAAL))

        assert reproduced is None
        assert "prints no fold" in note

    def test_a_fold_that_missed_the_paper_would_fail_rather_than_be_tuned(self):
        wrong = dataclasses.replace(_contrast(UKIBE), paper_stated_fold=5.0)

        assert module.reproduces_paper_fold(wrong)[0] is False


class TestWhatTheModelSaysAndWhyItIsAStructuralAnswer:
    def test_verwaal_gets_a_number_and_the_number_is_one(self, scored):
        """Not "unsupported". The phytoene node is `passthrough`, so crtI dosage cannot move
        the answer at ANY dosage, and that is a prediction rather than a gap."""
        response = scored[VERWAAL].model

        assert response.fold == 1.0
        assert "passthrough node 'phytoene'" in response.route

    def test_and_it_is_read_off_the_spec_not_looked_up_per_gene(self, chain):
        node = next(node for node in chain.nodes if node.enzyme == "CrtI")

        assert node.name == "phytoene"
        assert node.rate_law == RateLaw.PASSTHROUGH

    def test_the_discrepancy_is_the_whole_measured_fold(self, scored):
        assert scored[VERWAAL].measured_fold == pytest.approx(11.81, abs=0.01)
        assert scored[VERWAAL].discrepancy == pytest.approx(11.81, abs=0.01)

    def test_ukibe_gets_no_number_because_its_gene_is_upstream_of_the_chain(self, scored, chain):
        response = scored[UKIBE].model

        assert response.fold is None
        assert chain.precursor_metabolite in response.route
        assert chain.entry_enzyme in response.reason

    def test_a_capacity_gene_change_is_refused_and_the_refusal_names_the_experiment(self, chain):
        """`capacity_for` must not quietly return the anchor for a dosage nothing measured,
        so a contrast moving crtYB comes back as a refusal carrying that message."""
        crtyb = dataclasses.replace(
            _contrast(VERWAAL), changed_gene=CAPACITY_SETTING_ENZYME["beta_carotene"],
            changed_node=None)
        response = module.model_response(crtyb, chain)

        assert response.fold is None
        assert "WHAT WOULD SETTLE IT" in response.reason

    def test_a_node_that_is_not_in_the_spec_refuses(self, chain):
        missing = dataclasses.replace(_contrast(VERWAAL), changed_node="carotene")

        with pytest.raises(module.ContrastProvenanceError, match="is not a node"):
            module.model_response(missing, chain)

    def test_a_node_drained_by_a_different_enzyme_refuses(self, chain):
        mismatched = dataclasses.replace(_contrast(VERWAAL), changed_node="lycopene")

        with pytest.raises(module.ContrastProvenanceError, match="is drained by"):
            module.model_response(mismatched, chain)

    def test_declaring_a_gene_upstream_when_it_drains_a_node_refuses(self, chain):
        """The two branches must not both be reachable for the same gene, or the verdict
        would depend on which field the author happened to fill in."""
        lying = dataclasses.replace(_contrast(VERWAAL), changed_node=None)

        with pytest.raises(module.ContrastProvenanceError, match="drains node"):
            module.model_response(lying, chain)

    def test_a_saturating_node_gets_no_number_either(self):
        """Only a passthrough node licenses the exact 1.000. A fitted saturating step has
        parameters, just no dosage coefficient, and must not be reported as a prediction."""
        chain = PathwaySpec(
            product="prod", organism="S. cerevisiae", entry_enzyme="E0", precursor_metabolite="p",
            nodes=(Node("n0", RateLaw.SATURATING, Fate.DILUTED, "CrtI"),
                   Node("prod", RateLaw.PASSTHROUGH, Fate.DILUTED, "", molar_mass_g_per_mol=536.87)))
        contrast = dataclasses.replace(_contrast(VERWAAL), changed_node="n0")
        response = module.model_response(contrast, chain)

        assert response.fold is None
        assert "saturating node 'n0'" in response.route


class TestTheVerdictsAndTheirBoundary:
    def test_verwaal_straddles_the_ceiling_within_one_paper(self, scored):
        one = scored[VERWAAL]

        assert one.control_fold_of_ceiling == pytest.approx(0.401, abs=0.001)
        assert one.derivative_fold_of_ceiling == pytest.approx(4.741, abs=0.001)
        assert one.verdict == "refutes_the_flat_ceiling_within_one_paper"

    def test_ukibe_is_not_a_ceiling_test_and_is_not_counted_as_one(self, scored):
        """21.7x, and it refutes nothing here: both arms are already inside the bound. A
        one-sided bracket cannot be refuted from below."""
        one = scored[UKIBE]

        assert one.measured_fold == pytest.approx(21.67, abs=0.01)
        assert one.derivative_fold_of_ceiling < 1.0
        assert one.verdict == "not_a_ceiling_test"

    def test_only_one_of_the_two_contrasts_refutes(self, scored):
        assert sum(one.verdict.startswith("refutes") for one in scored.values()) == 1

    def test_a_pair_already_outside_the_bound_is_named_differently(self):
        """Refuting the ceiling and STRADDLING it are different claims. A pair whose control
        is already outside cannot say the added copy carried it across."""
        assert module.verdict_name(1.6, 16.0) == "refutes_the_ceiling_but_both_arms_are_already_outside"
        assert module.verdict_name(0.4, 4.7) == "refutes_the_flat_ceiling_within_one_paper"
        assert module.verdict_name(0.014, 0.31) == "not_a_ceiling_test"

    def test_a_derivative_exactly_at_the_bound_does_not_refute_it(self):
        """The bound is `content <= vmax_per_growth`, so equality satisfies it. An off-by-one
        here would turn the one row sitting closest to the bound into a refutation."""
        assert module.verdict_name(0.4, 1.0) == "not_a_ceiling_test"


class TestItStatesWhatItCannotShow:
    def test_neither_contrast_is_called_a_held_out_prediction(self, scored):
        for one in scored.values():
            assert any("not a held-out prediction" in limit for limit in one.limitations)

    def test_the_absent_dispersion_is_carried_on_both(self, scored):
        for one in scored.values():
            assert any("no dispersion" in limit for limit in one.limitations)
            assert not one.contrast.control.dispersion_reported

    def test_the_marker_confound_is_carried_on_both(self, scored):
        assert any("LEU2" in limit for limit in scored[VERWAAL].limitations)
        assert any("HIS3" in limit for limit in scored[UKIBE].limitations)

    def test_a_contrast_with_dispersion_and_no_confound_sheds_those_caveats(self, chain):
        """The caveats are computed from the record, so a better-reported pair must stop
        carrying them. Otherwise they are prose that happens to live in a tuple."""
        clean = dataclasses.replace(
            _contrast(VERWAAL), control=module.Arm("c", 0.5, 3, True),
            derivative=module.Arm("d", 5.9, 3, True), also_changed=())
        found = module.limitations(clean, module.model_response(clean, chain))

        assert not any("no dispersion" in limit for limit in found)
        assert not any("not isogenic" in limit for limit in found)
        assert any("not a held-out prediction" in limit for limit in found)


class TestTheThingsThisTableCannotSupport:
    def test_the_published_rows_cannot_run_a_leave_one_strain_out_validation(self):
        """The prior read-only pass reported 7/7. On all nine dfba_ready rows it is 9/9, and
        the reason is the same: no candidate has a complete finite inner score."""
        result = module.nested_validation_on_published_rows()

        assert len(result) == 9
        assert (result.status == "failed").all()
        assert set(result.failure_reason) == {
            "no candidate has a complete finite inner validation score"}

    def test_all_fifty_eight_candidates_were_offered_and_none_scored(self):
        result = module.nested_validation_on_published_rows()

        assert set(result.n_selection_candidates) == {58}
        assert result.selected_gene.isna().all()
        assert np.isinf(result.inner_rmse_log).all()

    def test_the_growth_rate_it_was_granted_is_a_fabrication_and_is_named_one(self):
        """Nine rows, none of which states a mu. The frame invents one so that the failure
        cannot be blamed on a missing growth rate."""
        result = module.nested_validation_on_published_rows()

        assert set(result.mu_per_h) == {module.GENEROUS_MU_PER_H}
        assert not hasattr(module.survey_row(VERWAAL), "mu_per_h")

    def test_between_paper_crte_dosage_explains_a_fifth_of_the_log_variance(self):
        found = module.dosage_identifiability()

        assert found["n_rows"] == 8 and found["n_levels"] == 3
        assert found["r2_dosage_as_factor"] == pytest.approx(0.205, abs=0.001)
        assert found["r2_linear_in_dosage"] == pytest.approx(0.173, abs=0.001)
        assert found["within_level_variance_fraction"] == pytest.approx(0.795, abs=0.001)

    def test_the_highest_published_content_carries_no_crte_at_all(self):
        """Which is why a dosage regression across papers cannot stand in for a matched
        contrast: the strongest strain is at the lowest level of the covariate."""
        found = module.dosage_identifiability()

        assert found["highest_content_reference"] == "Arhar 2024"
        assert found["highest_content_copies_crtE"] == 0

    def test_the_xie_ladder_is_declined_and_the_reason_is_the_rows_own_measurand(self):
        declined = next(one for one in module.DECLINED if one.reference == "Xie 2014")

        assert module.survey_row(declined.pmid).measurand == "both_reported"
        assert "copy_evidence" in declined.what_is_recorded
        assert declined.pmid not in {one.pmid for one in module.CONTRASTS}


class TestTheRunner:
    def test_it_runs_and_writes_both_records_when_asked(self, tmp_path, capsys):
        out = tmp_path / "matched_contrast.csv"

        assert module.main(["--out", str(out)]) == 0

        written = pd.read_csv(out)
        assert set(written.record) == {"ceiling_bracket", "matched_contrast"}
        assert len(written[written.record == "ceiling_bracket"]) == 9
        assert len(written[written.record == "matched_contrast"]) == 2

    def test_it_prints_the_verdicts_and_the_refusal(self, capsys):
        module.main([])
        printed = capsys.readouterr().out

        assert "refutes_the_flat_ceiling_within_one_paper" in printed
        assert "not_a_ceiling_test" in printed
        assert "9/9 outer folds failed" in printed
