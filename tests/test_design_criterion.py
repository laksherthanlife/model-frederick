"""Design on Fisher information, replacing three ad-hoc criteria with one statistical one.

The earlier score deduplicated rows, normalised them and summed log(1+sigma); a module
counted as identifiable if its projection weight cleared 0.35; and the general-stress
reporter had to be forced in by hand. All three were symptoms of the same gap: with no
noise model, "identifiable" had no statistical meaning.

With one, everything falls out of the same object. Observations are y = L x + e with
e ~ N(0, S), so the information is L' S^-1 L, a module's standard error is the root of
the corresponding diagonal of its inverse, and a target state is served by weighting that
target rather than by overriding the result.
"""

import numpy as np
import pytest

from ystwin.analysis.design import (
    NoiseModel,
    fisher_information,
    module_standard_errors,
    design_score,
    identifiable_by_precision,
)
from ystwin.generator.stress_panel import MODULES, reporter_loadings

QUIET = NoiseModel(relative_cv=0.02)
"""Deliberately quieter than any real reader, and named so.

Most tests here are about the shape of the criterion -- symmetry, rank, scale invariance,
whether a penalty lands on one channel -- and those are true at any noise level. Keeping an
explicit unrealistic value in them stops a reader mistaking the number for a claim about a
plate. The claims about a plate are in `TestTheDefaultNoiseIsTheMeasuredNoise`.
"""


class TestFisherInformation:
    def test_it_is_square_in_the_number_of_modules(self):
        info = fisher_information(reporter_loadings(["UPRE-ER", "TRX2-oxidative"]), QUIET)

        assert info.shape == (len(MODULES), len(MODULES))

    def test_it_is_symmetric_and_positive_semidefinite(self):
        info = fisher_information(reporter_loadings(["UPRE-ER", "FeRE-iron"]), QUIET)

        assert info == pytest.approx(info.T)
        assert np.linalg.eigvalsh(info).min() > -1e-9

    def test_a_noisier_reader_carries_less_information(self):
        loadings = reporter_loadings(["UPRE-ER", "FeRE-iron"])
        quiet = fisher_information(loadings, NoiseModel(relative_cv=0.01))
        noisy = fisher_information(loadings, NoiseModel(relative_cv=0.20))

        assert np.trace(quiet) > np.trace(noisy)

    def test_a_repeated_reporter_does_add_information(self):
        """Unlike the old score, a replicate legitimately counts -- it averages noise down."""
        once = fisher_information(reporter_loadings(["UPRE-ER"]), QUIET)
        twice = fisher_information(reporter_loadings(["UPRE-ER", "UPRE-ER"]), QUIET)

        assert np.trace(twice) > np.trace(once)


class TestStandardErrorsReplaceTheMagicThreshold:
    def test_every_module_gets_a_standard_error(self):
        errors = module_standard_errors(
            reporter_loadings(list(__import__(
                "ystwin.generator.stress_panel", fromlist=["REPORTERS"]).REPORTERS)), QUIET)

        assert set(errors) == set(MODULES)

    def test_an_unobserved_module_has_infinite_error(self):
        errors = module_standard_errors(reporter_loadings(["FeRE-iron"]), QUIET)

        assert np.isinf(errors["UPR"])
        assert np.isfinite(errors["iron"])

    def test_a_module_seen_only_through_crosstalk_is_estimated_poorly(self):
        errors = module_standard_errors(reporter_loadings(["UPRE-ER", "HSE-heat"]), QUIET)

        assert errors["ESR"] > errors["UPR"]

    def test_more_reporters_shrink_the_errors(self):
        few = module_standard_errors(reporter_loadings(["UPRE-ER", "FeRE-iron"]), QUIET)
        many = module_standard_errors(
            reporter_loadings(["UPRE-ER", "FeRE-iron", "HSE-heat", "STRE-general"]), QUIET)

        assert many["ESR"] < few["ESR"]

    def test_identifiability_is_a_precision_claim_not_a_projection_weight(self):
        loadings = reporter_loadings(["STRE-general", "UPRE-ER", "HSE-heat", "FeRE-iron"])

        strict = identifiable_by_precision(loadings, QUIET, min_effect=0.05)
        lenient = identifiable_by_precision(loadings, QUIET, min_effect=1.0)

        assert set(strict) <= set(lenient)

    def test_a_noisier_reader_makes_fewer_modules_identifiable(self):
        loadings = reporter_loadings(["STRE-general", "UPRE-ER", "HSE-heat", "FeRE-iron"])

        clean = identifiable_by_precision(loadings, NoiseModel(relative_cv=0.01), min_effect=0.2)
        dirty = identifiable_by_precision(loadings, NoiseModel(relative_cv=0.50), min_effect=0.2)

        assert len(dirty) <= len(clean)


class TestTargetingReplacesForcingAReporterIn:
    def test_an_unweighted_score_is_d_optimal(self):
        """log det over the estimable directions, which for a full-rank set is log det."""
        loadings = reporter_loadings(list(__import__(
            "ystwin.generator.stress_panel", fromlist=["REPORTERS"]).REPORTERS))
        info = fisher_information(loadings, QUIET)
        sign, logdet = np.linalg.slogdet(info)

        assert sign > 0
        assert design_score(loadings, QUIET) == pytest.approx(float(logdet), rel=1e-6)

    def test_weighting_a_target_module_changes_which_set_wins(self):
        general = reporter_loadings(["STRE-general", "UPRE-ER", "HSE-heat", "FeRE-iron"])
        clean = reporter_loadings(["TRX2-oxidative", "UPRE-ER", "HSE-heat", "FeRE-iron"])

        unweighted = design_score(general, QUIET) - design_score(clean, QUIET)
        targeted = (design_score(general, QUIET, priority={"ESR": 20.0})
                    - design_score(clean, QUIET, priority={"ESR": 20.0}))

        assert targeted > unweighted

    def test_a_priority_on_a_module_no_set_reads_changes_nothing(self):
        loadings = reporter_loadings(["UPRE-ER", "FeRE-iron"])

        plain = design_score(loadings, QUIET)
        weighted = design_score(loadings, QUIET, priority={"osmotic": 50.0})

        assert weighted == pytest.approx(plain, rel=1e-9)

    def test_an_unknown_priority_module_is_refused(self):
        with pytest.raises(KeyError, match="vibes"):
            design_score(reporter_loadings(["UPRE-ER"]), QUIET, priority={"vibes": 2.0})


class TestTheScoreDoesNotLeanOnARidge:
    """A rank-deficient set has determinant zero, so log det needs a regulariser -- and the
    size of that regulariser silently decides how much an extra estimable module is worth.
    The pseudo-determinant ignores the null space instead, so nothing arbitrary is priced in.
    """

    def test_the_score_is_finite_for_a_rank_deficient_set(self):
        assert np.isfinite(design_score(reporter_loadings(["UPRE-ER"]), QUIET))

    def test_the_score_ignores_the_null_space(self):
        """Two reporters span two directions; the score must be their information alone."""
        loadings = reporter_loadings(["UPRE-ER", "FeRE-iron"])
        info = fisher_information(loadings, QUIET)
        nonzero = np.linalg.eigvalsh(info)[-2:]

        assert design_score(loadings, QUIET) == pytest.approx(float(np.sum(np.log(nonzero))), rel=1e-6)

    def test_it_does_not_reward_a_vanishing_eigenvalue(self):
        """Under a ridge, adding a near-collinear reporter can raise the score by ~log(1/ridge)."""
        two = design_score(reporter_loadings(["STRE-general", "UPRE-ER"]), QUIET)
        plus_duplicate = design_score(
            reporter_loadings(["STRE-general", "UPRE-ER", "STRE-general"]), QUIET)

        assert plus_duplicate - two < 5.0


class TestTheDefaultNoiseIsTheMeasuredNoise:
    """A default is what a caller who has not thought about the reader silently gets.

    It was 0.02, against 0.146 measured on this project's own plates -- so any precision
    claim made without an explicit noise argument was seven times too confident, and nothing
    in the code said so. The fix is not to forbid the default but to make it the measured
    number, which is the only value that cannot mislead.
    """

    def test_it_is_the_cv_measured_on_this_projects_plates(self):
        from ystwin.generator.panel_experiment import OBSERVED_ACTIVITY_CV

        assert NoiseModel().relative_cv == OBSERVED_ACTIVITY_CV

    def test_it_is_not_an_optimistic_round_number(self):
        """The regression pin. 0.02 and 0.05 were both in the tree; both were guesses."""
        assert NoiseModel().relative_cv > 0.10

    def test_standard_errors_scale_linearly_in_it(self):
        """Which is why the correction moves every precision claim and no ranking.

        Multiplicative noise enters the information as 1/cv^2 uniformly across channels, so
        the whole covariance scales by cv^2 and every standard error by cv. Two consequences
        follow, and they are opposite: the ordering of sensor sets is invariant, and the
        smallest resolvable effect is not.
        """
        loadings = reporter_loadings(["STRE-general", "UPRE-ER", "HSE-heat", "FeRE-iron"])
        low = module_standard_errors(loadings, NoiseModel(relative_cv=0.05))
        high = module_standard_errors(loadings, NoiseModel(relative_cv=0.146))

        for name, error in low.items():
            assert high[name] == pytest.approx(error * 0.146 / 0.05, rel=1e-9)

    def test_the_correction_costs_real_verdicts_at_a_fine_min_effect(self):
        """The linearity is not a licence to ignore the change.

        Four modules the old default called resolvable at an effect of 0.1 are not, and
        saying so is the deliverable: the ranking survived, the claims did not.
        """
        loadings = reporter_loadings(["STRE-general", "UPRE-ER", "HSE-heat", "FeRE-iron"])
        optimistic = identifiable_by_precision(loadings, QUIET, min_effect=0.1)
        measured = identifiable_by_precision(loadings, NoiseModel(), min_effect=0.1)

        assert len(optimistic) == 4
        assert measured == []


class TestSelectingOnWhatCanActuallyBeMeasured:
    def test_the_value_counts_identifiable_modules_first(self):
        from ystwin.analysis.design import design_value

        one = design_value(reporter_loadings(["UPRE-ER"]), QUIET, min_effect=0.5)
        three = design_value(
            reporter_loadings(["UPRE-ER", "FeRE-iron", "HSE-heat"]), QUIET, min_effect=0.5)

        assert three[0] > one[0]

    def test_information_only_breaks_ties(self):
        from ystwin.analysis.design import design_value

        value = design_value(reporter_loadings(["UPRE-ER", "FeRE-iron"]), QUIET, min_effect=0.5)

        assert value[0] == 2
        assert value[1] == pytest.approx(design_score(reporter_loadings(["UPRE-ER", "FeRE-iron"]), QUIET))

    def test_a_priority_makes_one_module_worth_more_than_another(self):
        from ystwin.analysis.design import design_value

        plain = design_value(reporter_loadings(["UPRE-ER", "FeRE-iron"]), QUIET, min_effect=0.5)
        weighted = design_value(
            reporter_loadings(["UPRE-ER", "FeRE-iron"]), QUIET, min_effect=0.5,
            priority={"iron": 9.0})

        assert weighted[0] == pytest.approx(plain[0] + 8.0)


class TestNoiseScalesWithSignal:
    """Plate noise is multiplicative, so a brighter reporter is not a better one.

    Treating noise as a constant per channel makes information grow with the size of a
    reporter's loadings, so a promoter that picks up crosstalk outscores a clean one purely
    for reading more loudly. That is backwards: crosstalk is what makes a module hard to
    separate. With relative noise the two effects cancel and only orthogonality is left,
    which is why row-normalising was the right instinct in the earlier ad-hoc score even
    though it had no justification there. Here it falls out of the noise model.
    """

    def test_scaling_every_loading_leaves_the_score_unchanged(self):
        loadings = reporter_loadings(["UPRE-ER", "FeRE-iron"])

        assert design_score(2.5 * loadings, QUIET) == pytest.approx(
            design_score(loadings, QUIET), rel=1e-6)

    def test_a_clean_sensor_is_not_beaten_by_a_louder_crosstalking_one(self):
        """roGFP2 reads one pool exactly; HSE-heat reads its own module plus an STRE."""
        clean = fisher_information(reporter_loadings(["roGFP2-Grx1"]), QUIET)
        crosstalking = fisher_information(reporter_loadings(["HSE-heat"]), QUIET)

        assert np.trace(clean) >= np.trace(crosstalking)

    def test_a_noise_floor_still_penalises_a_dim_reporter(self):
        """With an additive floor, absolute brightness matters again -- as it should."""
        floored = NoiseModel(relative_cv=0.02, floor=0.5)
        loadings = reporter_loadings(["FeRE-iron"])

        assert design_score(3.0 * loadings, floored) > design_score(0.1 * loadings, floored)

    def test_orthogonal_sensors_still_beat_overlapping_ones(self):
        orthogonal = design_score(reporter_loadings(["FeRE-iron", "ZRE-zinc"]), QUIET)
        overlapping = design_score(reporter_loadings(["STRE-general", "STRE-osmotic"]), QUIET)

        assert orthogonal > overlapping

    def test_a_row_of_zeros_carries_no_information_and_does_not_divide_by_zero(self):
        loadings = np.zeros((1, len(MODULES)))

        assert np.all(np.isfinite(fisher_information(loadings, QUIET)))
        assert np.trace(fisher_information(loadings, QUIET)) == pytest.approx(0.0)


class TestTheGrowthConfoundBelongsInTheNoiseModel:
    """A promoter fusion pays a variance penalty a ratiometric sensor does not.

    Reporter protein is diluted by growth, so turning fluorescence into promoter activity
    means dividing out a growth rate that was itself estimated from noisy optical density.
    That correction has variance, and it lands on every transcriptional channel. A
    ratiometric sensor reports an equilibrium between two forms of one molecule, so the
    ratio cancels concentration and growth never enters.

    Without this the criterion cannot tell the two kinds apart, because after normalising
    for multiplicative noise a clean fusion and a clean biochemical sensor look identical.
    They are not: one of them needs a growth measurement to be interpreted at all.
    """

    def test_a_growth_penalty_lowers_the_information_a_channel_carries(self):
        loadings = reporter_loadings(["STRE-general"])
        clean = fisher_information(loadings, QUIET)
        penalised = fisher_information(loadings, QUIET, extra_cv=np.array([0.30]))

        assert np.trace(penalised) < np.trace(clean)

    def test_it_applies_per_channel_not_to_the_whole_set(self):
        loadings = reporter_loadings(["STRE-general", "roGFP2-Grx1"])
        both = fisher_information(loadings, QUIET, extra_cv=np.array([0.30, 0.30]))
        only_first = fisher_information(loadings, QUIET, extra_cv=np.array([0.30, 0.0]))

        assert np.trace(only_first) > np.trace(both)

    def test_no_penalty_reproduces_the_plain_information(self):
        loadings = reporter_loadings(["UPRE-ER", "FeRE-iron"])

        assert fisher_information(loadings, QUIET, extra_cv=np.zeros(2)) == pytest.approx(
            fisher_information(loadings, QUIET))

    def test_a_wrong_length_penalty_is_refused(self):
        with pytest.raises(ValueError, match="one entry per channel"):
            fisher_information(reporter_loadings(["UPRE-ER"]), QUIET, extra_cv=np.zeros(3))


class TestSelectionPrefersSensorsThatDoNotNeedGrowthDividedOut:
    def test_a_ratiometric_sensor_is_chosen_once_the_penalty_is_real(self):
        from ystwin.analysis.sensor_selection import select_sensors
        from ystwin.generator.stress_panel import REPORTERS, Kind

        chosen = select_sensors(n_channels=5, noise=NoiseModel(relative_cv=0.05),
                                min_effect=0.5, growth_cv=0.30)

        assert any(REPORTERS[r].kind is Kind.RATIOMETRIC for r in chosen)

    def test_without_a_penalty_the_two_kinds_are_interchangeable(self):
        from ystwin.analysis.design import design_score

        fusion = design_score(reporter_loadings(["FeRE-iron"]), QUIET)
        sensor = design_score(reporter_loadings(["roGFP2-Grx1"]), QUIET)

        assert fusion == pytest.approx(sensor, rel=1e-9)
