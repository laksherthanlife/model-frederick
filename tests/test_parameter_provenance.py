"""The inventory that makes "build in silico, calibrate later" operable.

Every refusal in this package already names its own missing measurement -- `calibrated_kinetics`
names the phytoene peak, `calibrate.py` names the dry weight, `solve_pathway` names the
degradation rate. Scattered across a dozen docstrings that is discipline nobody can act on.
Gathered into one table it is a work plan.

These tests keep it a work plan: every constant carries a grade, a source and the experiment
that would move it, and the grades stay distinguishable from each other.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

_SCRIPT = (pathlib.Path(__file__).resolve().parents[1]
           / "scripts" / "parameter_provenance.py")


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("parameter_provenance", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def table(script):
    return script.inventory()


class TestEveryConstantIsAccountedFor:
    def test_nothing_is_ungraded(self, script, table):
        grades = {script.MEASURED, script.BOUNDED, script.BORROWED,
                  script.ASSERTED, script.REFUSED}

        assert set(table.grade) <= grades
        assert not table.grade.isna().any()

    def test_every_row_names_where_it_came_from(self, table):
        for row in table.itertuples():
            assert len(row.source) > 20, row.parameter

    def test_every_row_names_what_would_move_it(self, table):
        """A parameter with no stated route to improvement is one nobody can act on, which
        is the state this table exists to end."""
        for row in table.itertuples():
            assert len(row.upgrade) > 10, row.parameter

    def test_all_five_grades_are_used(self, table):
        """If a grade never appears the distinction is decorative."""
        assert table.grade.nunique() == 5


class TestTheGradesSayWhatTheyClaim:
    def test_the_fitted_scalar_is_measured_and_carries_its_n(self, script, table):
        row = table[table.parameter.str.startswith("alpha")].iloc[0]

        assert row.grade == script.MEASURED
        assert "6 chemostat states" in row.source

    def test_the_borrowed_physiology_says_whose_culture_it_is(self, script, table):
        borrowed = table[table.grade == script.BORROWED]

        assert len(borrowed) >= 2
        for row in borrowed.itertuples():
            assert "van Hoek" in row.source

    def test_the_asserted_dry_weight_factor_names_its_confound(self, script, table):
        """It is not merely unmeasured: `calibrate.py` fits it as one product with the
        optics gain and reports the pair as unidentifiable."""
        row = table[table.parameter == "gdcw_per_od"].iloc[0]

        assert row.grade == script.ASSERTED
        assert "unidentifiable" in row.source

    def test_the_refused_entries_have_no_value(self, script, table):
        for row in table[table.grade == script.REFUSED].itertuples():
            assert row.value == "none", row.parameter


class TestTheCeilingIsFlaggedDespiteBeingMeasured:
    """The row that proves a grade is not a quality score. It is MEASURED -- fitted, on real
    data -- and it is wrong by 63x, which is worse than any asserted constant here."""

    def test_the_ceiling_is_measured(self, script, table):
        row = table[table.parameter == "content ceiling"].iloc[0]

        assert row.grade == script.MEASURED

    def test_and_its_note_says_it_is_refuted(self, table):
        row = table[table.parameter == "content ceiling"].iloc[0]

        assert "REFUTED" in row.upgrade
        assert "63x" in row.upgrade

    def test_its_value_matches_the_live_calibration(self, table):
        """Read from the code, not typed in, so it cannot drift from what the model uses."""
        from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS

        expected = BETA_CAROTENE_KINETICS["lycopene"].vmax_per_growth * 536.87
        row = table[table.parameter == "content ceiling"].iloc[0]

        assert float(row.value) == pytest.approx(expected, abs=5e-4)


class TestTheInventoryIsConnectedToTheNumbers:
    """The piece that was missing when this was checked on 2026-09-01: the inventory existed
    and the prediction existed, and nothing joined them. `ProductPrediction.layers` reports
    which LAYERS ran and says nothing about which CONSTANTS were measured — so the van Hoek
    caveat reached a caller only because somebody typed it at that call site, and a new
    borrowed constant would have arrived silently.
    """

    def test_every_reported_quantity_has_a_dependency_list(self, script):
        assert set(script.DEPENDS_ON) == {"rate", "content", "yield", "titre"}

    def test_an_unmapped_quantity_raises_rather_than_returning_nothing(self, script):
        """An empty frame would read as 'this rests on nothing', which is the opposite of
        true and exactly the wrong direction for this table to fail in."""
        with pytest.raises(KeyError, match="no dependency list"):
            script.provenance_for("productivity")

    def test_every_named_dependency_exists_in_the_inventory(self, script, table):
        """A dependency list naming a constant the table does not carry is a lie that makes
        the summary look better than it is."""
        known = set(table.parameter)
        for quantity, names in script.DEPENDS_ON.items():
            missing = [n for n in names if n not in known]
            assert missing == [], f"{quantity} names unknown constants: {missing}"

    def test_rate_and_content_rest_only_on_measured_constants(self, script):
        """Which is why they are the two quantities worth quoting."""
        for quantity in ("rate", "content"):
            grades = set(script.provenance_for(quantity).grade)
            assert grades == {script.MEASURED}, f"{quantity} -> {grades}"

    def test_titre_and_yield_each_carry_a_borrowed_constant(self, script):
        """van Hoek's glucose-limited table, on a culture that is not this one. The summary
        has to surface that without anybody remembering to write it down."""
        for quantity in ("titre", "yield"):
            summary = script.provenance_summary(quantity)
            assert "borrowed" in summary
            assert "weakest link BORROWED" in summary

    def test_the_summary_names_the_weakest_link_not_the_average(self, script):
        """A number resting on six measured constants and one asserted one is as good as the
        asserted one. Averaging would hide exactly the constant a reader needs."""
        assert script.provenance_summary("rate").endswith("MEASURED")
        assert script.provenance_summary("titre").endswith("BORROWED")
