"""Every product scored, or the exact column that stops it -- as a test, not a paragraph.

`scripts/score_all_products.py` answers "does this work on more than one compound". The
answer has a positive half and a negative half, and the negative half is the one that rots:
a blocker recorded in prose gets quietly forgotten, and somebody scores PHB on a `q` column
that was computed from the content it is being checked against.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pytest

_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "score_all_products.py"


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("score_all_products", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheTwoProductsThatAreScored:
    def test_both_compounds_in_the_chain_are_scored(self, script):
        scored = {r["product"] for r in script.scored_products()}

        assert scored == {"beta_carotene", "lycopene"}

    def test_the_intermediate_is_the_harder_one(self, script):
        """Lycopene is the node before the product and fits about twice as badly. Pinned
        because the README quoted only beta-carotene's 14.2% until 2026-08-31, and a chain
        scored on one molecule and reported on one molecule is one nobody has seen fail."""
        by_product = {r["product"]: r for r in script.scored_products()}

        assert by_product["lycopene"]["median_abs_error_pct"] > \
               by_product["beta_carotene"]["median_abs_error_pct"]

    def test_neither_score_is_suspiciously_perfect(self, script):
        """A product fitted on its own measurement would come back near zero error."""
        for row in script.scored_products():
            assert row["median_abs_error_pct"] > 1.0


class TestPHBCannotBeScoredAndTheReasonIsArithmetic:
    """The trap. Eleven states, three carbon sources, content and q and mu all present --
    and `q` is `content x mu`, so `content = q/mu` recovers its own input. A fit that
    cannot fail is the most convincing wrong answer this repository has available."""

    def test_phb_is_reported_as_unscoreable(self, script):
        assert script.phb_is_circular()["status"] == "CANNOT SCORE"

    def test_the_circularity_is_measured_not_asserted(self, script):
        """If Kocharin ever publishes an independently measured q, this test fails and the
        refusal should be revisited -- which is the point of computing it rather than
        writing it down."""
        import pandas as pd

        repo = pathlib.Path(__file__).resolve().parents[1]
        frame = pd.read_csv(repo / "data" / "phb" / "kocharin2013_chemostat_states.tsv",
                            sep="\t")
        frame = frame[frame.phb_mmol_per_gdcw.notna() & frame.q_phb_mmol_per_gdcw_h.notna()]
        recovered = frame.q_phb_mmol_per_gdcw_h / frame.mu_per_h

        assert np.allclose(recovered, frame.phb_mmol_per_gdcw, rtol=1e-4)


class TestTheNegativesStayNamed:
    def test_glycogen_loads_and_is_refused_at_solve(self, script):
        """It moved on 2026-09-01 when `solve_pathway` gained a degradation outlet. The
        guarantee is unchanged -- no number without the rate -- and the refusal now names
        the missing measurement instead of calling the pathway out of scope."""
        rows = {r["product"]: r for r in script.spec_refusals()}

        assert rows["glycogen"]["status"] == "REFUSED AT SOLVE"
        assert "degradation_rate_per_h" in rows["glycogen"]["detail"]

    def test_every_investigated_candidate_names_its_blocking_column(self, script):
        assert len(script.NOT_VENDORED) >= 5
        for name, where, why in script.NOT_VENDORED:
            assert name and where and len(why) > 40, f"{name} has no stated blocker"

    def test_the_two_chemostat_near_misses_are_recorded(self, script):
        """alpha-santalene and resveratrol have real q and mu and no expression. They are
        the cheapest third product and must not drop off the list."""
        names = {name for name, _, _ in script.NOT_VENDORED}

        assert "alpha-santalene" in names
        assert "resveratrol" in names
