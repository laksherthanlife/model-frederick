"""Grounded parameters and the published results the generator must reproduce.

Every value carries its source. The expectations are separate from the parameters on
purpose: a generator that reproduces a published observation it was not fitted to is
evidence; one tuned until it matches is decoration.
"""

import pytest

from ystwin.generator.literature import (
    EXPECTATIONS,
    LIBRARY,
    check_expectations,
    parameters_from_literature,
)


class TestTheLibraryIsCitable:
    def test_every_value_names_a_source(self):
        for name, entry in LIBRARY.items():
            assert entry.source, name

    def test_every_value_carries_units(self):
        for name, entry in LIBRARY.items():
            assert entry.units, name

    def test_values_taken_on_trust_are_marked_as_such(self):
        """Anything not read in full must say so, so it can be checked later."""
        unverified = [n for n, e in LIBRARY.items() if not e.verified]
        assert unverified, "the unverified flag exists precisely because some are"
        for name in unverified:
            assert LIBRARY[name].note, name


class TestPublishedExpectations:
    def test_total_hac1_transcript_is_expected_not_to_move(self):
        assert EXPECTATIONS["hac1_total_fold"].value == pytest.approx(1.0, abs=0.2)

    def test_kar2_is_expected_to_move(self):
        assert EXPECTATIONS["kar2_fold"].value > 1.5

    def test_low_peroxide_is_expected_to_be_sublethal_in_wild_type(self):
        assert EXPECTATIONS["h2o2_sublethal_mM"].value >= 0.5

    def test_the_reference_phenotype_is_present(self):
        for key in ("batch_growth_rate", "batch_glucose_uptake", "batch_ethanol"):
            assert key in EXPECTATIONS


class TestTheGeneratorReproducesThem:
    def test_a_construct_can_be_built_from_the_library(self):
        params = parameters_from_literature("UPRE1")

        assert params.mu_max > 0
        assert params.promoter_peak > params.promoter_basal

    def test_peroxide_at_the_sublethal_dose_barely_slows_growth(self):
        """The published claim, checked against the model rather than assumed."""
        params = parameters_from_literature("NativeYap1")
        dose = EXPECTATIONS["h2o2_sublethal_mM"].value

        assert params.growth_rate_at(dose) > 0.75 * params.mu_max

    def test_dtt_at_the_published_stress_dose_clearly_slows_growth(self):
        params = parameters_from_literature("UPRE2")
        dose = EXPECTATIONS["dtt_stress_mM"].value

        assert params.growth_rate_at(dose) < 0.7 * params.mu_max

    def test_the_sublethal_dose_still_induces_the_reporter(self):
        """A window where stress rises and growth barely moves is what decoupling needs."""
        params = parameters_from_literature("NativeYap1")
        dose = EXPECTATIONS["h2o2_sublethal_mM"].value

        induction = params.promoter_activity_at(dose) / params.promoter_basal
        assert induction > 1.2

    def test_all_four_constructs_are_available(self):
        for construct in ("UPRE1", "UPRE2", "NativeYap1", "AlteredYap1"):
            assert parameters_from_literature(construct) is not None

    def test_an_unknown_construct_is_refused(self):
        with pytest.raises(KeyError, match="mScarletSensor"):
            parameters_from_literature("mScarletSensor")


def test_the_expectation_check_reports_pass_and_fail_per_item():
    report = check_expectations()

    assert set(report.columns) >= {"expectation", "expected", "observed", "passed"}
    assert len(report) >= 3


def test_the_expectation_check_is_not_vacuous():
    """It must be capable of failing, or it is testing nothing."""
    report = check_expectations()

    assert report.observed.notna().any()
