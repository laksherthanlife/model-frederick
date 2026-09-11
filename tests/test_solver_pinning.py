"""A flux number is a property of the solver as well as of the model.

cobrapy picks its solver from whatever optlang finds installed. This package pins
``swiglpk``, so a clean install gets GLPK -- but a machine that also has Gurobi or CPLEX
silently gets a different simplex, which breaks ties among alternate optima differently.
Nothing in the code would notice.

That matters here specifically because `fba/fva.py` measured the product flux interval
at relative width 1.000: the feasible range runs from zero to the ceiling. A wide
interval is exactly where the particular optimum returned is arbitrary, so an unpinned
solver would make every single-solve flux partly a machine artefact.
"""

import cobra
import pytest

from ystwin.fba import solver


class TestTheSettingsAreActuallyApplied:
    def test_configure_sets_the_pinned_solver_and_tolerance(self):
        model = cobra.Model("probe")

        settings = solver.configure(model)

        assert solver.PINNED_SOLVER in type(model.solver).__module__
        assert model.tolerance == pytest.approx(solver.PINNED_TOLERANCE)
        assert settings.solver == solver.PINNED_SOLVER

    def test_the_tolerance_is_pinned_at_the_default_not_tightened(self):
        """Deliberate. The goal is that two machines agree, not that the numbers move --
        tightening it would silently restate every FBA result in the repository."""
        assert solver.PINNED_TOLERANCE == pytest.approx(cobra.Model("probe").tolerance)

    def test_an_explicit_tolerance_overrides_the_default(self):
        model = cobra.Model("probe")

        settings = solver.configure(model, tolerance=1e-9)

        assert model.tolerance == pytest.approx(1e-9)
        assert settings.tolerance == pytest.approx(1e-9)

    def test_fva_workers_are_pinned_to_one(self):
        """FVA forks a worker pool by default. The count cannot change the answer in
        principle, and does change whether a run reproduces."""
        assert solver.FVA_PROCESSES == 1


class TestTheSettingsCanBeRecorded:
    def test_they_carry_the_cobra_version(self):
        """Every table here records the conditions it was produced under. A flux table
        without its solver is missing one of them."""
        settings = solver.configure(cobra.Model("probe"))

        assert settings.cobra_version == cobra.__version__

    def test_they_render_as_one_line_for_a_table_cell(self):
        rendered = str(solver.configure(cobra.Model("probe")))

        assert solver.PINNED_SOLVER in rendered
        assert cobra.__version__ in rendered


class TestItRefusesRatherThanFallingBack:
    def test_an_unavailable_solver_raises(self, monkeypatch):
        """A silent fallback would reintroduce the machine dependence while looking
        configured, which is worse than not pinning at all."""
        monkeypatch.setattr(solver, "PINNED_SOLVER", "no_such_solver")

        with pytest.raises(RuntimeError, match="unavailable"):
            solver.configure(cobra.Model("probe"))

    def test_the_refusal_says_why_it_matters(self, monkeypatch):
        monkeypatch.setattr(solver, "PINNED_SOLVER", "no_such_solver")

        with pytest.raises(RuntimeError, match="alternate optima"):
            solver.configure(cobra.Model("probe"))


class TestTheRangeItProtects:
    def test_a_degenerate_toy_range_is_reported_as_full_width(self):
        """The property that makes pinning load-bearing, on a model small enough to
        reason about: two routes of equal cost to the same product, so which one carries
        flux is arbitrary and the interval is the whole of it."""
        from ystwin.fba.fva import product_flux_range

        model = cobra.Model("toy")
        a, b, p = (cobra.Metabolite(i, compartment="c") for i in ("a", "b", "p"))
        feed = cobra.Reaction("feed", lower_bound=0, upper_bound=10)
        feed.add_metabolites({a: 1})
        r1 = cobra.Reaction("route1", lower_bound=0, upper_bound=10)
        r1.add_metabolites({a: -1, p: 1})
        r2 = cobra.Reaction("route2", lower_bound=0, upper_bound=10)
        r2.add_metabolites({a: -1, p: 1})
        grow = cobra.Reaction("growth", lower_bound=0, upper_bound=10)
        grow.add_metabolites({a: -1, b: 1})
        sink_b = cobra.Reaction("sink_b", lower_bound=0, upper_bound=10)
        sink_b.add_metabolites({b: -1})
        out = cobra.Reaction("DM_p", lower_bound=0, upper_bound=10)
        out.add_metabolites({p: -1})
        model.add_reactions([feed, r1, r2, grow, sink_b, out])

        rng = product_flux_range(model, "DM_p", glucose_uptake=None,
                                 growth_fraction=0.5, biomass_reaction="growth")

        assert rng.minimum == pytest.approx(0.0)
        assert rng.maximum > 0.0
        assert rng.relative_width == pytest.approx(1.0)
