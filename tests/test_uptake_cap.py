"""Capping uptake has to follow the model's own structure, not a reaction id.

An enzyme-constrained GSMM is irreversibly split: the exchange keeps its id but is fixed
at ``(0, 0)`` in the export direction, and a generated ``<id>_REV`` carries the supply.
Setting ``lower_bound = -magnitude`` on the export then caps nothing -- and worse, opens a
route that creates extracellular substrate. On ecYeastGEM_batch that made
``aerobic_batch_constraints`` a silent no-op: growth identical before and after, glucose
flowing uncapped at 17.93 against a requested 21.3.

These build the split by hand so the rule is tested without the GSMM files, which are not
in the repository.
"""

import cobra
import pytest

from ystwin.fba.physiology import cap_uptake


def _unsplit():
    """A plain reversible exchange, as Yeast9 carries it."""
    model = cobra.Model("unsplit")
    s = cobra.Metabolite("s_0565", compartment="e")
    ex = cobra.Reaction("r_1714", lower_bound=-1.0, upper_bound=1000.0)
    ex.add_metabolites({s: -1.0})
    sink = cobra.Reaction("consume", lower_bound=0.0, upper_bound=1000.0)
    sink.add_metabolites({s: -1.0})
    model.add_reactions([ex, sink])
    model.objective = sink
    return model


def _split():
    """The irreversible split an ec model carries: export pinned, _REV supplies."""
    model = cobra.Model("split")
    s = cobra.Metabolite("s_0565", compartment="e")
    out = cobra.Reaction("r_1714", lower_bound=0.0, upper_bound=0.0)
    out.add_metabolites({s: -1.0})
    rev = cobra.Reaction("r_1714_REV", lower_bound=0.0, upper_bound=1000.0)
    rev.add_metabolites({s: 1.0})
    sink = cobra.Reaction("consume", lower_bound=0.0, upper_bound=1000.0)
    sink.add_metabolites({s: -1.0})
    model.add_reactions([out, rev, sink])
    model.objective = sink
    return model


class TestItPicksTheReactionThatCarriesSupply:
    def test_an_unsplit_model_is_capped_on_the_exchange(self):
        model = _unsplit()

        assert cap_uptake(model, "r_1714", 21.3) == "r_1714"
        assert model.reactions.get_by_id("r_1714").lower_bound == pytest.approx(-21.3)

    def test_a_split_model_is_capped_on_the_reverse_reaction(self):
        """The regression. Constraining `r_1714` here caps nothing."""
        model = _split()

        assert cap_uptake(model, "r_1714", 21.3) == "r_1714_REV"
        assert model.reactions.get_by_id("r_1714_REV").upper_bound == pytest.approx(21.3)

    def test_it_reads_the_split_from_the_model_not_from_a_name(self):
        """A model named like an ec model but not split must still be capped correctly,
        and the reverse. Filenames are not a signal."""
        plain = _unsplit()
        plain.id = "ecYeastGEM_batch"

        assert cap_uptake(plain, "r_1714", 5.0) == "r_1714"


class TestTheCapActuallyBinds:
    def test_uptake_is_limited_on_a_split_model(self):
        """What the bug cost: the objective was free to consume without limit."""
        model = _split()
        assert model.slim_optimize() > 900.0, "unbounded before the cap"

        cap_uptake(model, "r_1714", 21.3)

        assert model.slim_optimize() == pytest.approx(21.3)

    def test_uptake_is_limited_on_an_unsplit_model(self):
        model = _unsplit()
        cap_uptake(model, "r_1714", 21.3)

        assert model.slim_optimize() == pytest.approx(21.3)

    def test_the_old_approach_would_have_left_it_unbounded(self):
        """States the defect as a fact rather than as prose, so the fix cannot silently
        revert to the id-based form."""
        model = _split()
        model.reactions.get_by_id("r_1714").lower_bound = -21.3

        assert model.slim_optimize() > 900.0

    def test_and_would_have_opened_a_route_that_creates_substrate(self):
        """The export runs backwards under a negative lower bound, so the constraint
        intended to restrict supply becomes a second source of it."""
        model = _split()
        model.reactions.get_by_id("r_1714_REV").upper_bound = 0.0
        assert model.slim_optimize() == pytest.approx(0.0), "no supply route left"

        model.reactions.get_by_id("r_1714").lower_bound = -21.3

        assert model.slim_optimize() == pytest.approx(21.3)


def test_a_magnitude_is_taken_as_a_magnitude():
    """Callers pass either sign; uptake is the same either way."""
    a, b = _split(), _split()
    cap_uptake(a, "r_1714", 21.3)
    cap_uptake(b, "r_1714", -21.3)

    assert a.reactions.get_by_id("r_1714_REV").upper_bound == pytest.approx(
        b.reactions.get_by_id("r_1714_REV").upper_bound)


@pytest.mark.parametrize("factory", [_unsplit, _split])
def test_explicit_conditions_can_reopen_supply_and_replace_a_previous_medium_bound(factory):
    model = factory()

    for magnitude in (0.0, 2.0, 21.3, 0.0, 5.0):
        cap_uptake(model, "r_1714", magnitude)
        solution = model.optimize(raise_error=True)
        assert solution.status == "optimal"
        assert solution.objective_value == pytest.approx(magnitude)
