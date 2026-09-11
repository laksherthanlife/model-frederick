"""Choosing reporters under a channel cap, scored by Fisher information.

Selection used to maximise a spanning score with no noise model, which forced two
workarounds: a bare 0.35 threshold to call a module "identifiable", and hand-forcing the
general-stress reporter into the answer. Both are gone. Precision defines identifiability
and a priority weight expresses what the build is for, so the selector is asked the right
question rather than corrected after it answers.
"""

import pytest

from ystwin.analysis.design import NoiseModel, design_value, identifiable_by_precision
from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD, select_sensors
from ystwin.generator.stress_panel import MODULES, REPORTERS, reporter_loadings

PLATE = NoiseModel(relative_cv=0.05)


def modules_of(names, noise=PLATE, min_effect=0.5):
    return identifiable_by_precision(reporter_loadings(names), noise, min_effect=min_effect)


class TestSelection:
    def test_it_returns_the_requested_number_of_sensors(self):
        assert len(select_sensors(n_channels=4, noise=PLATE)) == 4

    def test_it_only_returns_reporters_that_exist(self):
        assert set(select_sensors(n_channels=4, noise=PLATE)) <= set(REPORTERS)

    def test_it_beats_an_arbitrary_set_of_the_same_size(self):
        chosen = design_value(
            reporter_loadings(select_sensors(n_channels=4, noise=PLATE)), PLATE, 0.5)
        arbitrary = design_value(
            reporter_loadings(["STRE-general", "STRE-osmotic", "UPRE-ER", "HSE-heat"]), PLATE, 0.5)

        assert chosen >= arbitrary

    def test_it_reserves_reference_channels(self):
        assert len(select_sensors(n_channels=5, n_reference=1, noise=PLATE)) == 4

    def test_it_honours_a_required_reporter(self):
        assert "UPRE-ER" in select_sensors(n_channels=3, required=["UPRE-ER"], noise=PLATE)

    def test_it_refuses_a_budget_the_library_cannot_fill(self):
        with pytest.raises(ValueError, match="library has only"):
            select_sensors(n_channels=99, noise=PLATE)

    def test_it_refuses_a_required_reporter_that_does_not_exist(self):
        with pytest.raises(KeyError, match="not in the library"):
            select_sensors(n_channels=3, required=["GFP-vibes"], noise=PLATE)

    def test_it_never_returns_a_duplicate(self):
        chosen = select_sensors(n_channels=4, noise=PLATE)

        assert len(set(chosen)) == len(chosen)


class TestNoiseChangesTheAnswer:
    def test_a_noisier_reader_identifies_no_more_modules(self):
        clean = modules_of(select_sensors(n_channels=4, noise=NoiseModel(relative_cv=0.01)),
                           NoiseModel(relative_cv=0.01))
        dirty = modules_of(select_sensors(n_channels=4, noise=NoiseModel(relative_cv=0.30)),
                           NoiseModel(relative_cv=0.30))

        assert len(dirty) <= len(clean)

    def test_more_channels_identify_at_least_as_many_modules(self):
        four = modules_of(select_sensors(n_channels=4, noise=PLATE))
        six = modules_of(select_sensors(n_channels=6, noise=PLATE))

        assert len(six) >= len(four)

    def test_the_whole_panel_identifies_more_than_half_the_landscape(self):
        assert len(modules_of(list(REPORTERS))) > len(MODULES) / 2


class TestPriorityReplacesForcing:
    def test_prioritising_the_general_state_makes_it_identifiable(self):
        """The old code forced STRE-general in by hand; now the objective asks for ESR."""
        targeted = select_sensors(n_channels=4, noise=PLATE, priority={"ESR": 20.0})

        assert "ESR" in modules_of(targeted)

    def test_the_targeted_build_measures_as_many_modules_as_the_free_one(self):
        """Prioritising ESR should redirect the set, not shrink what it can resolve."""
        free = modules_of(select_sensors(n_channels=4, noise=PLATE))
        targeted = modules_of(select_sensors(n_channels=4, noise=PLATE, priority={"ESR": 20.0}))

        assert len(targeted) >= len(free)

    def test_a_priority_is_not_a_guarantee_of_a_particular_reporter(self):
        """It weights the goal, so any set that serves it may win -- that is the point."""
        chosen = select_sensors(n_channels=4, noise=PLATE, priority={"iron": 20.0})

        assert "iron" in modules_of(chosen)


class TestRecommendedBuild:
    def test_it_fits_five_channels_with_one_reference(self):
        assert len(RECOMMENDED_BUILD.stress_reporters) == 4
        assert RECOMMENDED_BUILD.n_channels == 5

    def test_it_names_a_slow_and_a_fast_fluorophore_for_the_growth_pair(self):
        assert RECOMMENDED_BUILD.reference_fluorophore == "mOrange2"

    def test_it_reaches_the_general_stress_state(self):
        assert "ESR" in modules_of(RECOMMENDED_BUILD.stress_reporters)

    def test_it_records_the_noise_it_was_designed_against(self):
        assert RECOMMENDED_BUILD.noise.relative_cv > 0

    def test_every_module_it_claims_really_is_identifiable(self):
        claimed = set(RECOMMENDED_BUILD.identifiable)
        actual = set(modules_of(RECOMMENDED_BUILD.stress_reporters,
                                RECOMMENDED_BUILD.noise, RECOMMENDED_BUILD.min_effect))

        assert claimed == actual

    def test_it_is_what_the_selector_returns_for_its_own_settings(self):
        assert set(RECOMMENDED_BUILD.stress_reporters) == set(select_sensors(
            n_channels=RECOMMENDED_BUILD.n_channels, n_reference=1,
            noise=RECOMMENDED_BUILD.noise, priority=RECOMMENDED_BUILD.priority,
            min_effect=RECOMMENDED_BUILD.min_effect, growth_cv=RECOMMENDED_BUILD.growth_cv))


def test_an_explicit_empty_library_does_not_fall_back_to_the_entire_panel():
    with pytest.raises(ValueError, match="library has only 0"):
        select_sensors(1, library=[])


def test_required_reporters_cannot_overrun_or_repeat_channels():
    with pytest.raises(ValueError, match="exceed"):
        select_sensors(1, required=["UPRE-ER", "TRX2-oxidative"])
    with pytest.raises(ValueError, match="unique"):
        select_sensors(2, required=["UPRE-ER", "UPRE-ER"])
    with pytest.raises(ValueError, match="multiple library entries"):
        select_sensors(2, library=["UPRE-ER", "UPRE-ER"])


def test_no_interpretable_build_is_a_named_refusal():
    name = next(name for name, reporter in REPORTERS.items()
                if reporter.ph_sensitive and reporter.module != "ph")
    with pytest.raises(ValueError, match="no spectrally compatible, interpretable"):
        select_sensors(1, library=[name])


def test_growth_correction_noise_cannot_be_negative_or_nonfinite():
    from ystwin.analysis.sensor_selection import growth_penalty_cv

    for error in (-0.1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="standard error"):
            growth_penalty_cv(error, 0.2)
