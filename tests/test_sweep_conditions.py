"""The condition sweep: the project's own question asked in bulk.

`scripts/sweep_conditions.py` exists because "given these conditions, how much do I get"
is what the vision promises, and until it was written the only way to ask was to compose
`predict_product` by hand. These tests pin the two things that make the answer honest: the
numbers move in the direction the biology requires, and the inputs that reach nothing are
visibly not reaching anything.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[1]
_SCRIPT = _REPO / "scripts" / "sweep_conditions.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("sweep_conditions", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def table(script):
    from ystwin.pathway.spec import load_pathway
    return script.sweep(load_pathway("beta_carotene"), feed_g_per_L=20.0)


class TestTheSweepAnswersAndRefusesInTheRightPlaces:
    def test_most_of_the_grid_is_answered(self, table):
        assert table.answered.sum() >= 15

    def test_it_refuses_outside_the_fitted_growth_range(self, table):
        """The kinetics were fitted over [0.101, 0.254] /h. A sweep that answered at 0.30
        would be extrapolating a law nobody checked there, which is the failure the solver's
        range guard was added to prevent."""
        refused = table[~table.answered]

        assert set(refused.dilution_rate_per_h) == {0.10, 0.30}
        assert all("outside the range" in r for r in refused.refusal)

    def test_every_refusal_says_why(self, table):
        assert all(r.strip() for r in table[~table.answered].refusal)


class TestTheNumbersMoveTheWayTheBiologyRequires:
    def test_content_falls_as_the_pump_runs_faster(self, table):
        """Growth dilutes the pool, so content per gDCW must fall monotonically with D at
        fixed genotype. This is the one qualitative claim the pathway solve makes."""
        for strain, group in table[table.answered].groupby("strain"):
            ordered = group.sort_values("dilution_rate_per_h").content_mg_per_gdcw

            assert (ordered.diff().dropna() < 0).all(), f"{strain} content did not fall"

    def test_content_rises_with_entry_enzyme_expression(self, table):
        for rate, group in table[table.answered].groupby("dilution_rate_per_h"):
            ordered = group.sort_values("entry_expression").content_mg_per_gdcw

            assert (ordered.diff().dropna() > 0).all(), f"D={rate} content did not rise"

    def test_productivity_rises_with_dilution_where_content_falls(self, table):
        """The trade-off a fermenter is actually run on: less per cell, more cells per hour."""
        best = table[table.answered].sort_values("productivity_mg_per_L_h").iloc[-1]

        assert best.dilution_rate_per_h == 0.25


class TestWhatDoesNotReachTheAnswerIsVisiblyNotReachingIt:
    def test_the_feed_scales_titre_but_not_content(self, script):
        """Doubling the feed doubles the cells and therefore the titre. It is arithmetic,
        not metabolism: content per gDCW is untouched, because nothing in the chain reads
        the carbon supply. Kocharin's series moves the product flux 4.24x at fixed mu on
        feed alone, so the real effect exists and this is not it."""
        from ystwin.pathway.spec import load_pathway
        spec = load_pathway("beta_carotene")

        lean = script.sweep(spec, feed_g_per_L=10.0)
        rich = script.sweep(spec, feed_g_per_L=20.0)

        lean_ok, rich_ok = lean[lean.answered], rich[rich.answered]
        assert lean_ok.content_mg_per_gdcw.tolist() == pytest.approx(
            rich_ok.content_mg_per_gdcw.tolist())
        assert rich_ok.titre_mg_per_L.tolist() == pytest.approx(
            (2.0 * lean_ok.titre_mg_per_L).tolist())
