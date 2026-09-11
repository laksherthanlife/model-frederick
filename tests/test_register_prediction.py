"""The one prediction in this repository written down before its measurement.

`docs/CLAIM_BOUNDARY.md` opens with "Nothing in this repository is Tier 3", and the reason
was never that the predictions were bad: it is that every score here was computed after the
measurement it scores. `outputs/registered_prediction_D018.csv` is the first thing that can
be wrong in public.

That only means something if the file can be regenerated. A registered prediction nobody
can re-derive is a number in a file, not a commitment — and until this script existed the
reproducibility audit correctly reported that nothing in `scripts/` wrote it.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pandas as pd
import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "register_prediction.py"
_COMMITTED = _REPO / "outputs" / "registered_prediction_D018.csv"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("register_prediction", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def registered(script):
    return script.register(script.REGISTERED_RATE)


class TestItReproducesTheRegisteredFile:
    def test_the_committed_table_is_still_there(self):
        assert _COMMITTED.exists()

    def test_every_number_matches_the_committed_one(self, registered):
        """Exactly, not approximately. A registered prediction that drifts is not a
        registered prediction; it is a moving target with a filename."""
        committed = pd.read_csv(_COMMITTED)

        pd.testing.assert_frame_equal(
            registered.sort_values("strain").reset_index(drop=True),
            committed.sort_values("strain").reset_index(drop=True))


class TestTheRegisteredSchemaIsFrozen:
    """A pre-registration records what was predicted AND under what plan.

    On 2026-09-02 `fba/fedbatch.py` replaced the chemostat and `Environment.dilution_rate`
    became `growth_rate_setpoint_per_h`. A blanket rename swept this file's column with it,
    which broke the reproduction -- correctly, because the committed table is a SEALED
    RECORD and the rename would have edited it after the fact. The column stays
    `dilution_rate` because the vessel it was registered for was a chemostat.

    These tests exist so the next rename cannot quietly do the same thing.
    """

    def test_the_committed_column_is_still_the_registered_one(self):
        committed = pd.read_csv(_COMMITTED)

        assert "dilution_rate" in committed.columns
        assert "growth_rate_setpoint_per_h" not in committed.columns

    def test_the_script_still_writes_the_registered_column(self, registered):
        assert "dilution_rate" in registered.columns

    def test_the_registered_rate_is_still_feasible_under_the_new_vessel(self, script):
        """The substantive question the rename raises. A sealed record is worth keeping only
        if it can still be run -- and 0.18 /h passes BOTH of the fed-batch refusals, where
        the upper Elizondo state at 0.254 /h passes neither."""
        from ystwin.fba.fedbatch import design_fedbatch

        design = design_fedbatch(
            script.REGISTERED_RATE, 0.254,
            initial_volume_l=1.0, max_volume_l=2.0, initial_biomass_g_per_l=0.05,
            feed_substrate_g_per_l=500.0, assumed_yield_g_per_g=0.098,
            max_uptake_g_per_gdcw_h=2.59)

        assert design.capacity_margin > 0.25
        assert design.initial_uptake_g_per_gdcw_h < 2.59


class TestOnlyTheEnvironmentAndOneGENOTYPENumberGoIn:
    """This class used to be called ...AndOneStrainConstantGoIn, and the constant was the
    strain's own measured ``q_lycopene + q_betacarotene``. A prediction registered in
    advance that rests on a measurement of the quantity being predicted is still circular,
    and registering it early does not fix that. The input is now relative CrtE expression.
    """

    def test_the_input_is_expression_and_nothing_about_the_product(self, registered):
        assert "entry_expression" in registered.columns
        assert not any("q_" in c or "measured" in c for c in registered.columns)

    def test_it_is_the_strain_s_own_entry_gene_expression(self, script, registered):
        directory = script.paths.data_dir() / "carotenoid"

        expected = script.strain_expression(directory)

        for row in registered.itertuples():
            assert row.entry_expression == pytest.approx(expected[row.strain])

    def test_the_gene_it_used_is_recorded_beside_every_row(self, registered, script):
        """A relative expression number means nothing without saying relative to what."""
        assert set(registered.entry_gene) == {
            script.BETA_CAROTENE_FLUX.entry_enzyme}

    def test_the_genotype_does_not_depend_on_the_dilution_rate(self, script):
        """The whole reason a genotype number can be carried to an unmeasured rate. If it
        were re-derived per rate there would be nothing to register."""
        one = script.register(0.15).set_index("strain").entry_expression
        other = script.register(0.22).set_index("strain").entry_expression

        pd.testing.assert_series_equal(one, other)

    def test_the_flux_is_predicted_rather_than_supplied(self, registered):
        assert (registered.predicted_pathway_flux > 0).all()

    def test_and_so_is_the_intermediate(self, registered):
        """Lycopene comes OUT. If it were supplied the prediction would turn on the one
        number nobody measured."""
        assert (registered.predicted_lycopene_content > 0).all()


class TestTheRateItRegistersAt:
    def test_it_is_inside_the_calibrated_range_and_is_neither_endpoint(self, script):
        """A third rate is only informative between the two that were fitted. Outside them
        the prediction is an extrapolation and cannot separate `vmax ~ mu` from `vmax ~
        mu^2` any better than the endpoints already do."""
        from ystwin.kinetic.carotenoid import ELIZONDO2025

        low, high = ELIZONDO2025.growth_rate_range

        assert low < script.REGISTERED_RATE < high
        assert script.REGISTERED_RATE not in (low, high)

    def test_a_rate_outside_the_calibration_is_refused_rather_than_registered(self, script):
        from ystwin.generator.culture import measured_growth_rate_range

        low, _ = measured_growth_rate_range()
        with pytest.raises(ValueError, match="outside the measured chemostat curve"):
            script.register(low / 2)

    def test_a_reachable_rate_outside_the_kinetic_fit_is_refused(self, script):
        with pytest.raises(ValueError, match="outside the range its kinetics were fitted"):
            script.register(0.3)

    def test_an_unreachable_rate_is_refused_before_extrapolation(self, script):
        from ystwin.predict import SetpointUnreachable

        with pytest.raises(SetpointUnreachable, match="at or below the setpoint"):
            script.register(0.9)

    def test_the_filename_a_prediction_is_cited_by_does_not_move(self, script, tmp_path,
                                                                 monkeypatch):
        """``D018``, not ``D0180``. A ``:.3f`` in the stem would have renamed the file the
        claim boundary cites, and a registered prediction that changes address is one
        nobody can go back to."""
        monkeypatch.setenv("YSTWIN_OUTPUTS", str(tmp_path))
        monkeypatch.setattr(script.sys, "argv", ["register_prediction.py"])

        script.main()

        assert (tmp_path / "registered_prediction_D018.csv").exists()


class TestTheBandIsReportedAndNotHidden:
    def test_every_row_carries_a_band_that_brackets_its_flux(self, registered):
        assert (registered.flux_low <= registered.predicted_pathway_flux).all()
        assert (registered.predicted_pathway_flux <= registered.flux_high).all()

    def test_the_band_is_the_held_out_error_and_not_a_parameter_interval(self, registered,
                                                                        script):
        """What a user of this number faces is the spread of held-out predictions, not a
        confidence interval on a fitted constant. Quoting the second where the first is
        meant is how [0.85, 1.20] got shipped for a parameter whose real interval was
        [0.53, 1.78] -- see tests/test_growth_exponent_identifiability.py."""
        expected = script.BETA_CAROTENE_FLUX.typical_fold_error
        spans = registered.flux_high / registered.flux_low

        assert spans.max() == pytest.approx(expected ** 2, rel=1e-6)
