"""The environment channel, and the two limits that keep it honest.

`scripts/product_environment_sweep.py` measured that content was a strict function of mu:
62 environments at a held growth rate returned ONE content value. `pathway/environment_flux.py`
carries the measurement that says yeast disagrees -- carbon source moves PHB content 3.82x at
an identical growth rate, 64 standard deviations, one genotype and one vessel.

These tests pin three things:
  * the channel EXISTS, from the raw data, with no model in the way;
  * the effect REVERSES with growth rate, so a single coefficient is the wrong shape --
    this is the finding an audit caught after the fitted-scalar form had already been
    proposed, and it is the reason the module fits an interaction;
  * the coefficients are NOT established and the module refuses to export them to a product
    that has not measured its own.
"""

from __future__ import annotations

import math
import pathlib

import pandas as pd
import pytest

from ystwin.pathway.environment_flux import (
    KOCHARIN_PHB,
    MEASURED_ENVIRONMENT_RESPONSES,
    EnvironmentUnmeasured,
    carbon_descriptor,
    environment_factor,
)

_STATES = pathlib.Path(__file__).resolve().parents[1] / "data" / "phb" / \
    "kocharin2013_chemostat_states.tsv"


@pytest.fixture(scope="module")
def states():
    import pandas as pd
    if not _STATES.exists():
        pytest.skip("Kocharin chemostat states not present")
    return pd.read_csv(_STATES, sep="\t")


class TestTheChannelExistsInTheRawData:
    """No fit, no model. Same genotype, same vessel, carbon matched, mu equal by row."""

    def test_one_genotype_and_one_vessel(self, states):
        """The confound that would explain it away, ruled out by the design itself."""
        assert states.feed_carbon_cmol_per_l.nunique() == 1
        assert len(states) == 11

    def test_carbon_source_moves_content_at_identical_growth_rate(self, states):
        low = states[states.mu_per_h == 0.05]
        glucose = low[low.carbon_source == "glucose"].phb_mg_per_gdw.iloc[0]
        ethanol = low[low.carbon_source == "ethanol"].phb_mg_per_gdw.iloc[0]

        assert glucose == pytest.approx(4.33)
        assert ethanol == pytest.approx(16.55)
        assert ethanol / glucose > 3.8

    def test_the_separation_is_far_beyond_the_measurement_error(self, states):
        """64 standard deviations. This is why the CHANNEL needs no significance argument,
        even though its functional form cannot get one."""
        low = states[states.mu_per_h == 0.05]
        g = low[low.carbon_source == "glucose"].iloc[0]
        e = low[low.carbon_source == "ethanol"].iloc[0]
        pooled = math.hypot(g.phb_sd_mg_per_gdw, e.phb_sd_mg_per_gdw)

        assert (e.phb_mg_per_gdw - g.phb_mg_per_gdw) / pooled > 50.0

    def test_the_effect_beats_the_calibration_noise_floor(self, states):
        """22.2% is the leave-one-strain-out error of the entry-flux law. Any channel below
        it cannot be checked; this one clears it by more than an order of magnitude."""
        from ystwin.pathway import calibrations

        low = states[states.mu_per_h == 0.05]
        fold = (low[low.carbon_source == "ethanol"].phb_mg_per_gdw.iloc[0]
                / low[low.carbon_source == "glucose"].phb_mg_per_gdw.iloc[0])
        floor = calibrations.BETA_CAROTENE_FLUX.typical_fold_error

        assert (fold - 1.0) > 10 * (floor - 1.0)


class TestTheEffectReversesWithGrowthRate:
    """The audit finding. A single scalar coefficient would be wrong above the crossover, and
    the crossover sits INSIDE the measured window rather than out at an extrapolated edge."""

    def test_ethanol_helps_at_low_growth_and_hurts_at_high(self, states):
        def fold(mu):
            block = states[states.mu_per_h == mu]
            return (block[block.carbon_source == "ethanol"].phb_mg_per_gdw.iloc[0]
                    / block[block.carbon_source == "glucose"].phb_mg_per_gdw.iloc[0])

        assert fold(0.05) > 3.0
        assert fold(0.15) < 1.0

    def test_the_fitted_law_reproduces_that_reversal(self):
        assert environment_factor("phb", 1.0, 0.05) > 3.0
        assert environment_factor("phb", 1.0, 0.20) < 1.0

    def test_the_sign_flip_is_inside_the_measured_window(self):
        low, high = KOCHARIN_PHB.growth_rate_range
        flip = KOCHARIN_PHB.sign_flip_growth_rate_per_h

        assert low < flip < high
        assert flip == pytest.approx(0.1718, abs=0.001)

    def test_a_constant_coefficient_would_be_wrong_above_the_flip(self):
        """Pins why the interaction term is not optional. Without it the law carries one
        direction everywhere, and the data changes direction at 0.1718 /h."""
        below = environment_factor("phb", 1.0, 0.10)
        above = environment_factor("phb", 1.0, 0.20)

        assert (below - 1.0) * (above - 1.0) < 0

    def test_no_effect_when_the_feed_is_all_glucose(self):
        for mu in (0.05, 0.10, 0.20):
            assert environment_factor("phb", 0.0, mu) == pytest.approx(1.0)


class TestItRefusesToExportTheCoefficients:
    def test_beta_carotene_has_no_measured_response(self):
        assert "beta_carotene" not in MEASURED_ENVIRONMENT_RESPONSES

    def test_and_asking_for_one_is_refused(self):
        with pytest.raises(EnvironmentUnmeasured, match="no measured environment response"):
            environment_factor("beta_carotene", 1.0, 0.15)

    def test_the_refusal_explains_why_transfer_is_not_safe(self):
        """PHB is acetyl-CoA-limited and ethanol feeds acetyl-CoA directly through ACS --
        a pathway mechanism, not a host one. The refusal has to say that, or a reader will
        reasonably assume the coefficient is host-level and reuse it."""
        with pytest.raises(EnvironmentUnmeasured, match="pathway mechanism and not a host"):
            environment_factor("gadusol", 1.0, 0.15)

    def test_the_refusal_names_what_would_lift_it(self):
        with pytest.raises(EnvironmentUnmeasured, match="SAME growth rate by"):
            environment_factor("glycogen", 1.0, 0.15)


class TestTheCoefficientsAreNotEstablished:
    """Three carbon sources means the smallest attainable permutation p is 1/3! = 0.1667.
    The module carries that number so anyone quoting the law also has the reason it is not
    significant -- the same trap the entry-flux ranking claim fell into at 1/3! for strains."""

    def test_the_permutation_floor_is_reported(self):
        assert KOCHARIN_PHB.permutation_floor == pytest.approx(1.0 / 6.0)

    def test_it_cannot_reach_conventional_significance(self):
        assert KOCHARIN_PHB.permutation_floor > 0.05

    def test_the_response_names_its_source_and_its_sample_size(self):
        assert "Kocharin" in KOCHARIN_PHB.source
        assert KOCHARIN_PHB.n_states == 11
        assert KOCHARIN_PHB.n_carbon_sources == 3

    def test_the_summary_carries_the_floor(self):
        assert "permutation floor" in KOCHARIN_PHB.summary()


class TestTheCarbonDescriptor:
    def test_it_is_a_carbon_fraction_not_a_mass_fraction(self):
        """The feeds were matched on CARBON, which is what makes this a carbon-quality
        contrast rather than a carbon-quantity one. A mass-based descriptor would confound
        the two and would not give 0.677 for Kocharin's 1:2 mix."""
        assert carbon_descriptor(6.35, 10.21) == pytest.approx(0.677, abs=0.001)

    @pytest.mark.parametrize("glucose,ethanol,expected", [
        (20.0, 0.0, 0.0),
        (0.0, 15.32, 1.0),
    ])
    def test_the_pure_feeds_sit_at_the_ends(self, glucose, ethanol, expected):
        assert carbon_descriptor(glucose, ethanol) == pytest.approx(expected)

    def test_it_matches_the_committed_feed_compositions(self, states):
        for _, row in states.iterrows():
            z = carbon_descriptor(row.feed_glucose_g_per_l, row.feed_ethanol_g_per_l)
            if row.carbon_source == "glucose":
                assert z == pytest.approx(0.0)
            elif row.carbon_source == "ethanol":
                assert z == pytest.approx(1.0)
            else:
                assert z == pytest.approx(0.677, abs=0.001)

    def test_a_carbonless_feed_is_refused(self):
        with pytest.raises(ValueError, match="no carbon"):
            carbon_descriptor(0.0, 0.0)

    @pytest.mark.parametrize("fraction", [-0.1, 1.1])
    def test_an_impossible_fraction_is_refused(self, fraction):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            environment_factor("phb", fraction, 0.1)

    def test_a_non_positive_growth_rate_is_refused(self):
        with pytest.raises(ValueError, match="must be positive"):
            environment_factor("phb", 1.0, 0.0)


class TestTheLayerReportsAndDoesNotSetTheNumber:
    """`predict.py::_environment_layer`. The point of the layer is that the chain stops being
    SILENT about an axis where the truth is known to move by 3.82x."""

    @pytest.fixture(scope="class")
    def carotene(self):
        from ystwin.pathway.spec import load_pathway
        return load_pathway("beta_carotene")

    @staticmethod
    def _layer(spec, carbon_source, growth_rate):
        from ystwin.generator.context import CultureContext
        from ystwin.predict import Environment, _environment_layer
        return _environment_layer(
            spec, Environment(context=CultureContext(carbon_source=carbon_source)),
            growth_rate)

    def test_it_reports_rather_than_going_inert_for_an_unmeasured_product(self, carotene):
        """INERT means 'provably could not have changed anything'. That is false here: the
        channel is real, it is unmeasured for THIS product. Calling it INERT would record the
        absence of a coefficient as the absence of an effect."""
        from ystwin.predict import LayerState

        status = self._layer(carotene, "glucose", 0.15)

        assert status.state == LayerState.REPORTED
        assert not status.sets_the_number

    def test_the_report_says_the_answer_is_mu_only_on_a_moving_axis(self, carotene):
        status = self._layer(carotene, "ethanol", 0.15)

        assert "mu-only answer" in status.detail

    def test_an_unmapped_carbon_source_is_not_run_rather_than_interpolated(self, carotene):
        """galactose sits between glucose and ethanol on every physiological axis, which is
        exactly why interpolating it would be inventing a number."""
        from ystwin.predict import LayerState

        status = self._layer(carotene, "galactose", 0.15)

        assert status.state == LayerState.NOT_RUN

    def test_it_reports_rather_than_setting_the_number_even_where_it_was_measured(self):
        """Corrected 2026-09-04. This asserted ``IN_CHAIN`` and ``sets_the_number``, and
        both were false: the factor was computed, formatted into ``detail``, and never
        applied. The test passed anyway because it only ever looked at the LABEL --
        which is how a layer came to advertise "4.679x on content" over a content that
        did not move. See ``test_the_content_does_not_move`` below, which is the check
        this class was missing.
        """
        from ystwin.pathway.spec import load_pathway
        from ystwin.predict import LayerState

        status = self._layer(load_pathway("phb"), "ethanol", 0.05)

        assert status.state == LayerState.REPORTED
        assert not status.sets_the_number

    def test_the_report_says_it_was_not_applied_and_why(self):
        """A REPORTED layer that prints a factor without saying it is unapplied reads
        exactly like an applied one."""
        from ystwin.pathway.spec import load_pathway

        detail = self._layer(load_pathway("phb"), "ethanol", 0.05).detail

        assert "NOT applied" in detail
        assert "permutation floor" in detail

    def test_the_content_does_not_move_with_carbon_source(self):
        """THE regression, and the one no test made. ``sets_the_number`` is a claim about
        ``predict_product``'s return value, so it has to be checked there and not on the
        status object. If a future change wires the factor in, this fails and the label,
        the detail string and this test must all move together."""
        from ystwin.generator.context import CultureContext
        from ystwin.pathway.flux import FluxCalibration
        from ystwin.pathway.spec import load_pathway
        from ystwin.predict import Environment, Genotype, predict_product

        spec = load_pathway("phb")
        calibration = FluxCalibration(
            alpha=1e-3, entry_enzyme=spec.entry_enzyme, loso_rmse_log=0.0,
            loso_skill=0.0, n_states=2, expression_range=(1.0, 1.0),
            source="synthetic layer fixture, not a PHB calibration")
        contents = [
            predict_product(
                spec, Genotype(1.0, "synthetic-phb"),
                Environment(context=CultureContext(carbon_source=carbon),
                            growth_rate_setpoint_per_h=0.05),
                calibration, {}, mode="legacy").calculation.content_mg_per_gdcw
            for carbon in ("glucose", "ethanol")
        ]

        assert contents[0] == contents[1], (
            "the environment layer moved the content; if that is intended, it must return "
            "IN_CHAIN and sets_the_number=True, and the factor must be scored")

    def test_and_it_carries_the_reversal_through_the_layer(self):
        from ystwin.pathway.spec import load_pathway

        phb = load_pathway("phb")
        low = self._layer(phb, "ethanol", 0.05).detail
        high = self._layer(phb, "ethanol", 0.20).detail

        assert "4.679x" in low
        assert "0.827x" in high

    def test_the_layer_runs_on_every_prediction(self, carotene):
        from ystwin.generator.context import CultureContext
        from ystwin.pathway import calibrations
        from ystwin.predict import Environment, Genotype, predict_product

        got = predict_product(
            carotene, Genotype(entry_expression=1.0),
            Environment(context=CultureContext(carbon_source="glucose"),
                        growth_rate_setpoint_per_h=0.15),
            calibrations.BETA_CAROTENE_FLUX, calibrations.BETA_CAROTENE_KINETICS)

        assert [layer.name for layer in got.layers].count("environment") == 1

    def test_it_does_not_change_the_beta_carotene_number(self, carotene):
        """The layer reports; it must not set a number for a product with no coefficient.
        If this ever fails, PHB's acetyl-CoA coefficient has leaked into an isoprenoid."""
        from ystwin.generator.context import CultureContext
        from ystwin.pathway import calibrations
        from ystwin.predict import Environment, Genotype, predict_product

        def content(carbon_source):
            return predict_product(
                carotene, Genotype(entry_expression=1.0),
                Environment(context=CultureContext(carbon_source=carbon_source),
                            growth_rate_setpoint_per_h=0.12),
                calibrations.BETA_CAROTENE_FLUX,
                calibrations.BETA_CAROTENE_KINETICS,
                mode="empirical").content_mg_per_gdcw

        assert content("glucose") == content("ethanol")

    def test_the_layer_does_not_turn_unreachable_ethanol_into_a_prediction(self, carotene):
        from ystwin.generator.context import CultureContext
        from ystwin.pathway import calibrations
        from ystwin.predict import Environment, Genotype, SetpointUnreachable, predict_product

        with pytest.raises(SetpointUnreachable, match="at or below the setpoint"):
            predict_product(
                carotene, Genotype(entry_expression=1.0),
                Environment(context=CultureContext(carbon_source="ethanol"),
                            growth_rate_setpoint_per_h=0.15),
                calibrations.BETA_CAROTENE_FLUX,
                calibrations.BETA_CAROTENE_KINETICS, mode="empirical")


class TestTheLawMustNotBeExtrapolatedToANewFeed:
    """Measured, and it inverts the ranking. The interaction earns its place on a held-out
    STATE and is the worst of the three on a held-out FEED, where its worst case more than
    doubles. Both are true; the module documents the split rather than picking one.
    """

    @pytest.fixture(scope="class")
    def scores(self, states):
        import math

        import numpy as np

        frame = states.copy()
        frame["z"] = [carbon_descriptor(g, e) for g, e
                      in zip(frame.feed_glucose_g_per_l, frame.feed_ethanol_g_per_l)]
        frame["logc"] = np.log(frame.phb_mg_per_gdw)
        frame["logmu"] = np.log(frame.mu_per_h)
        frame["mz"] = frame.logmu * frame.z

        def loo(columns, group):
            design = np.column_stack(
                [np.ones(len(frame))] + [frame[c].values for c in columns])
            target = frame.logc.values
            groups = frame[group].values if group else np.arange(len(frame))
            errors = []
            for value in np.unique(groups):
                mask = groups != value
                if mask.sum() < design.shape[1]:
                    continue
                beta, *_ = np.linalg.lstsq(design[mask], target[mask], rcond=None)
                errors += list(abs(target[~mask] - design[~mask] @ beta))
            errors = np.array(errors)
            return math.exp(float(np.median(errors))), math.exp(errors.max())

        return {"interaction_state": loo(["logmu", "z", "mz"], None),
                "interaction_feed": loo(["logmu", "z", "mz"], "carbon_source"),
                "z_only_state": loo(["z"], None),
                "z_only_feed": loo(["z"], "carbon_source")}

    def test_the_interaction_wins_on_a_held_out_state(self, scores):
        assert scores["interaction_state"][0] < scores["z_only_state"][0]

    def test_and_loses_on_a_held_out_feed(self, scores):
        assert scores["interaction_feed"][0] > scores["z_only_feed"][0]

    def test_its_worst_case_more_than_doubles_on_an_unseen_feed(self, scores):
        """2.32x for z-only against 5.61x for the interaction. This is the number that says
        do not extrapolate."""
        assert scores["interaction_feed"][1] > 2 * scores["z_only_feed"][1]

    def test_the_module_says_so_where_it_will_be_read(self):
        import ystwin.pathway.environment_flux as module

        assert "DO NOT EXTRAPOLATE IT TO A FEED IT WAS NOT FITTED ON" in module.__doc__


class TestTheDescriptorCoversExactlyTheFeedsTheResponsesClaim:
    """A prescription is bound to the thing that has to realise it.

    The module docstring used to prescribe "a glycerol or acetate feed ... at which point the
    coefficients become earnable for the first time". `carbon_descriptor`'s entire signature
    is ``(glucose_g_per_l, ethanol_g_per_l)``, so that feed cannot be expressed at all and
    adding those rows could not earn anything for THIS law -- it needs a new descriptor, not
    new rows. The sentence was written from the descriptor's docstring rather than from its
    signature, and these tests make the signature the thing that has to agree.
    """

    DECLARED = {"glucose", "ethanol"}

    def test_the_signature_names_exactly_the_declared_carbon_species(self):
        """The check that fails when the prose widens and the code does not."""
        import inspect

        from ystwin.pathway.environment_flux import _CARBON_ATOMS, _MOLAR_MASS_G_PER_MOL

        parameters = set(inspect.signature(carbon_descriptor).parameters)
        from_signature = {name.removesuffix("_g_per_l") for name in parameters}

        assert from_signature == self.DECLARED
        assert set(_CARBON_ATOMS) == self.DECLARED
        assert set(_MOLAR_MASS_G_PER_MOL) == self.DECLARED

    def test_every_measured_response_is_fitted_on_feeds_the_descriptor_can_express(self,
                                                                                   states):
        """The other direction. A response claiming n carbon sources must be describable by
        the same n distinct descriptor values, or its coefficient is fitted on a contrast the
        descriptor cannot see."""
        from ystwin.pathway.environment_flux import MEASURED_ENVIRONMENT_RESPONSES

        response = MEASURED_ENVIRONMENT_RESPONSES["phb"]
        values = {round(carbon_descriptor(g, e), 6) for g, e
                  in zip(states.feed_glucose_g_per_l, states.feed_ethanol_g_per_l)}

        assert len(values) == response.n_carbon_sources

    def test_a_carbon_species_the_descriptor_does_not_carry_cannot_be_passed_at_all(self):
        """Structural refusal, not a silent zero -- there is no parameter to put it in."""
        with pytest.raises(TypeError):
            carbon_descriptor(glucose_g_per_l=0.0, glycerol_g_per_l=15.0)

    def test_and_routing_it_through_an_existing_parameter_would_be_silently_wrong(self):
        """Why the prescription mattered. If a glycerol feed were entered as anything the
        descriptor does carry, z would be a number about the WRONG species -- and entered as
        "not ethanol" it is exactly 0.0, so the law predicts a glycerol state as a glucose
        state and the environment factor comes back exactly 1.0."""
        from ystwin.pathway.environment_flux import environment_factor

        assert carbon_descriptor(20.0, 0.0) == 0.0
        assert environment_factor("phb", 0.0, 0.05) == pytest.approx(1.0)
        assert environment_factor("phb", 1.0, 0.05) > 4.0


class TestNoQuantityIsBroadcastAcrossLawsThatCannotShareIt:
    """One statistic must not be written onto rows that never produced it.

    `scripts/environment_channel.py` assigned `b_z`, `b_mu_z`, `sign_flip_growth_rate_per_h`
    and `permutation_p` as SCALARS to the whole frame, so five rows of
    `outputs/environment_law_scores.csv` carried the z-only law's relabelling p and the
    interaction fit's coefficients -- including the rows "content constant" and "mu only",
    which contain no z term at all. Nothing quoted the wrong cells because every existing
    marker happened to target row="law=z only"; but a marker addresses this file by
    (column, row), so a future marker for row="law=mu only" column="permutation_p" would have
    pinned prose to a number that law never computed.

    This is the same disease as a stale literal, one scale down: a value frozen and carried to
    rows that never produced it.
    """

    @pytest.fixture(scope="class")
    def scores(self):
        import pathlib as _pathlib

        import pandas as pd

        path = _pathlib.Path(__file__).resolve().parents[1] / "outputs" / \
            "environment_law_scores.csv"
        if not path.exists():
            pytest.skip("environment_law_scores.csv not committed")
        return pd.read_csv(path).set_index("law")

    def test_a_law_with_no_z_term_carries_no_z_coefficient(self, scores):
        for law in ("content constant", "mu only"):
            assert pd.isna(scores.loc[law, "b_z"]), law
            assert pd.isna(scores.loc[law, "b_mu_z"]), law

    def test_only_the_interaction_law_defines_a_sign_flip(self, scores):
        """The crossover needs both a main effect and an interaction to be solved for."""
        defined = scores.sign_flip_growth_rate_per_h.notna()

        assert list(scores.index[defined]) == ["mu + z + mu*z"]

    def test_a_law_that_cannot_read_z_gets_the_honest_permutation_p_of_one(self, scores):
        """Permuting the z labels cannot change a prediction that does not read z, so every
        relabelling ties and the p-value is 1.0 -- not the z-only law's 0.1667."""
        for law in ("content constant", "mu only"):
            assert scores.loc[law, "permutation_p"] == 1.0, law

    def test_the_permutation_p_is_not_constant_across_laws_with_different_terms(self, scores):
        """The general form of the check, so a re-broadcast fails immediately whichever
        column it lands in."""
        assert scores.permutation_p.nunique() > 1

    def test_the_z_only_row_still_holds_the_number_the_module_quotes(self, scores):
        """The one cell every existing marker points at. It must not have moved."""
        assert scores.loc["z only", "permutation_p"] == pytest.approx(0.1667, abs=1e-4)
        assert scores.loc["z only", "permutation_floor"] == pytest.approx(0.1667, abs=1e-4)

