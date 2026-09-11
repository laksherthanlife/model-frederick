"""Replacing hand-written chemistry with a generic walker, without moving a digit.

``kinetic/carotenoid.py::solve_branch_from_flux`` is beta-carotene written out by hand: one
quadratic, two named pools, three published fits hanging off it. ``pathway/solve.py`` does
the same arithmetic for any declared chain, which is only an improvement if it gives the
same answers -- a generic solver that agrees to three decimals has quietly changed every
number the repository has published.

So the first class here is a regression guard and nothing else: the same flux and growth
rate through both code paths, to machine precision, on the shipped calibration.

The rest pin the four closed forms the module claims. Each is a mass balance rather than a
fit, which is why they can be asserted exactly:

    saturating    mu*X^2 + (mu*km + vmax - v_in)*X - v_in*km = 0
    proportional  X = v_in/(k + mu)
    passthrough   X = 0
    terminal      X = v_in/mu

``carbon_closes`` is the same statement from the other end -- everything entering the chain
leaves as somebody's dilution term -- so it is zero by identity, not by tolerance, and an
arithmetic slip anywhere in the walk breaks it.
"""

from __future__ import annotations

import pytest

from ystwin.kinetic.carotenoid import ELIZONDO2025, solve_branch_from_flux
from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS
from ystwin.pathway.solve import NodeKinetics, solve_pathway
from ystwin.pathway.spec import Node, PathwaySpec, RateLaw, load_pathway

# Flux in mmol/gDCW/h, growth in /h. The growth rates stay inside the calibrated window of
# the hand-written branch, which refuses outside [0.101, 0.254].
CASES = (
    (1.0085e-3, 0.100987987),
    (5.0e-4, 0.15),
    (2.0e-3, 0.20),
    (7.5e-4, 0.254320182),
)


@pytest.fixture(scope="module")
def beta_carotene():
    return load_pathway("beta_carotene")


def _one_node(rate_law=RateLaw.PASSTHROUGH):
    """A chain that is nothing but its product, so the terminal form stands alone."""
    return PathwaySpec(
        product="squalene", organism="Saccharomyces cerevisiae",
        nodes=(Node("squalene", rate_law, molar_mass_g_per_mol=410.72),),
        entry_enzyme="ERG9", precursor_metabolite="s_0190")


def _two_nodes(rate_law):
    """One intermediate under ``rate_law``, then the product."""
    return PathwaySpec(
        product="squalene", organism="Saccharomyces cerevisiae",
        nodes=(Node("presqualene", rate_law, enzyme="ERG9"),
               Node("squalene", RateLaw.PASSTHROUGH, molar_mass_g_per_mol=410.72)),
        entry_enzyme="ERG9", precursor_metabolite="s_0190")


class TestItReproducesTheHandWrittenBranch:
    def test_the_lycopene_pool_matches_the_hand_written_quadratic(self, beta_carotene):
        """The intermediate is the number the generic walk could most easily get wrong.

        It comes out of a quadratic whose coefficients are assembled in a different file,
        from a calibration stored in a different form -- ``vmax_per_growth`` here against
        ``capacity_mmol_per_gdcw`` there. Agreement to machine precision says the two are
        the same law and not two roundings of it.
        """
        for flux, mu in CASES:
            generic = solve_pathway(beta_carotene, flux, mu, BETA_CAROTENE_KINETICS)
            hand = solve_branch_from_flux(flux, mu)

            assert generic.node("lycopene").content_mmol_per_gdcw == pytest.approx(
                hand.lycopene_content, rel=1e-15), f"flux={flux} mu={mu}"

    def test_the_beta_carotene_pool_matches_the_hand_written_branch(self, beta_carotene):
        for flux, mu in CASES:
            generic = solve_pathway(beta_carotene, flux, mu, BETA_CAROTENE_KINETICS)
            hand = solve_branch_from_flux(flux, mu)

            assert generic.terminal.content_mmol_per_gdcw == pytest.approx(
                hand.beta_carotene_content, rel=1e-15), f"flux={flux} mu={mu}"

    def test_the_cyclase_rate_matches_the_hand_written_branch(self, beta_carotene):
        """The rate, not just the pool: these differ by a factor of mu and both are quoted.

        ``beta_carotene_rate`` is what an HPLC on a chemostat pellet measures, and it is the
        flux leaving the lycopene node here. Checking only the contents would let a
        growth-rate factor go missing in one path and not the other.
        """
        for flux, mu in CASES:
            generic = solve_pathway(beta_carotene, flux, mu, BETA_CAROTENE_KINETICS)
            hand = solve_branch_from_flux(flux, mu)

            assert generic.node("lycopene").flux_out == pytest.approx(
                hand.beta_carotene_rate, rel=1e-15), f"flux={flux} mu={mu}"

    def test_the_calibrated_capacity_is_the_same_number_in_both_files(self):
        """The two modules hold the cyclase fit separately, and nothing joins them.

        ``calibrations.py`` states ``vmax_per_growth`` and ``ELIZONDO2025`` states
        ``capacity_mmol_per_gdcw``; they are the same fitted constant written twice. A
        refit that updates one and not the other would make the reproduction tests above
        fail with no indication of which copy moved.
        """
        assert (BETA_CAROTENE_KINETICS["lycopene"].vmax_per_growth
                == pytest.approx(ELIZONDO2025.capacity_mmol_per_gdcw, rel=1e-12))

    def test_the_calibrated_km_is_the_same_number_in_both_files(self):
        assert (BETA_CAROTENE_KINETICS["lycopene"].km
                == pytest.approx(ELIZONDO2025.km_mmol_per_gdcw, rel=1e-12))


class TestTheChainClosesOnCarbon:
    def test_nothing_is_lost_walking_the_shipped_pathway(self, beta_carotene):
        """An identity, so the tolerance is machine epsilon and not a fitted margin.

        At steady state every mole entering the chain leaves as somebody's ``mu*[X]``. If
        the walk ever dropped a node's dilution, or carried a flux forward that a step did
        not actually produce, this is the number that would move -- and it would move by
        the size of the mistake, not by a rounding.
        """
        for flux, mu in CASES:
            solution = solve_pathway(beta_carotene, flux, mu, BETA_CAROTENE_KINETICS)

            assert abs(solution.carbon_closes) < 1e-15, f"flux={flux} mu={mu}"

    def test_nothing_is_lost_through_a_proportional_step(self):
        for flux, mu in CASES:
            solution = solve_pathway(
                _two_nodes(RateLaw.PROPORTIONAL), flux, mu,
                {"presqualene": NodeKinetics(rate_constant=0.4)})

            assert abs(solution.carbon_closes) < 1e-15, f"flux={flux} mu={mu}"

    def test_a_zero_flux_leaves_every_pool_empty(self, beta_carotene):
        solution = solve_pathway(beta_carotene, 0.0, 0.15, BETA_CAROTENE_KINETICS)

        assert all(node.content_mmol_per_gdcw == 0.0 for node in solution.nodes)


class TestEachRateLawHasItsClosedForm:
    def test_a_terminal_pool_is_exactly_the_flux_over_the_growth_rate(self):
        """Nothing consumes the product, so washout is the only outlet there is.

        This one line is why the module is intracellular-only: a secreted product leaves
        through a transporter instead, and ``v_in/mu`` would then be a number about a
        culture that does not exist.
        """
        solution = solve_pathway(_one_node(), 1.0085e-3, 0.15)

        assert solution.terminal.content_mmol_per_gdcw == 1.0085e-3 / 0.15

    def test_a_passthrough_intermediate_holds_exactly_nothing(self, beta_carotene):
        """Declared of phytoene, and it is an assumption the spec makes visible.

        Exactly zero rather than nearly zero, because ``passthrough`` is a statement that
        the node has no pool at all -- if it acquires one, the declaration was wrong and
        the right fix is a rate law, not a small number.
        """
        solution = solve_pathway(beta_carotene, 1.0085e-3, 0.15, BETA_CAROTENE_KINETICS)

        assert solution.node("phytoene").content_mmol_per_gdcw == 0.0

    def test_a_passthrough_intermediate_passes_its_flux_on_undiminished(self,
                                                                       beta_carotene):
        solution = solve_pathway(beta_carotene, 1.0085e-3, 0.15, BETA_CAROTENE_KINETICS)

        assert solution.node("phytoene").flux_out == 1.0085e-3

    def test_a_proportional_pool_is_the_flux_over_the_rate_constant_plus_growth(self):
        """First order, so the pool has no ceiling and cannot saturate.

        That is a different claim from a Michaelis-Menten step with a large km, which is
        why the module refuses to reach one law by taking a limit of the other.
        """
        solution = solve_pathway(_two_nodes(RateLaw.PROPORTIONAL), 2e-3, 0.2,
                                 {"presqualene": NodeKinetics(rate_constant=0.4)})

        assert solution.node("presqualene").content_mmol_per_gdcw == 2e-3 / (0.4 + 0.2)

    def test_a_proportional_step_carries_its_rate_constant_times_its_pool(self):
        solution = solve_pathway(_two_nodes(RateLaw.PROPORTIONAL), 2e-3, 0.2,
                                 {"presqualene": NodeKinetics(rate_constant=0.4)})
        node = solution.node("presqualene")

        assert node.flux_out == pytest.approx(0.4 * node.content_mmol_per_gdcw, rel=1e-15)

    def test_a_saturating_step_never_carries_more_than_its_vmax(self):
        """The ceiling is the point of declaring the law, and it is where the root matters.

        The quadratic has one positive root; taking the other would give a negative pool
        and a flux above vmax, which is the arithmetic slip this bound catches.
        """
        vmax_per_growth, mu = 0.0023252, 0.15
        solution = solve_pathway(
            _two_nodes(RateLaw.SATURATING), 1.0, mu,
            {"presqualene": NodeKinetics(vmax_per_growth=vmax_per_growth, km=5.974e-4)})

        assert solution.node("presqualene").flux_out < vmax_per_growth * mu

    def test_a_saturating_pool_is_positive_even_under_a_flood(self):
        solution = solve_pathway(
            _two_nodes(RateLaw.SATURATING), 1.0, 0.15,
            {"presqualene": NodeKinetics(vmax_per_growth=0.0023252, km=5.974e-4)})

        assert solution.node("presqualene").content_mmol_per_gdcw > 0.0

    def test_content_is_reported_in_the_unit_the_literature_uses(self):
        """mg/gDCW is what a paper prints, and the conversion belongs beside the mass.

        ``None`` where no molar mass was declared, rather than a zero that would read as an
        empty pool.
        """
        solution = solve_pathway(_one_node(), 1.0085e-3, 0.15)
        terminal = solution.terminal

        assert terminal.content_mg_per_gdcw == pytest.approx(
            terminal.content_mmol_per_gdcw * 410.72, rel=1e-15)


class TestGrowthWashesTheProductOut:
    def test_a_faster_culture_holds_less_product_at_the_same_flux(self, beta_carotene):
        """The washout term, which is the only thing removing the product.

        A twin that reported the same content at every dilution rate would be describing a
        pool nothing empties. Every measured state in the calibration set is the reverse of
        that, and so is the arithmetic.
        """
        slow = solve_pathway(beta_carotene, 1.0085e-3, 0.11, BETA_CAROTENE_KINETICS)
        fast = solve_pathway(beta_carotene, 1.0085e-3, 0.25, BETA_CAROTENE_KINETICS)

        assert (fast.terminal.content_mmol_per_gdcw
                < slow.terminal.content_mmol_per_gdcw)

    def test_a_faster_culture_holds_less_of_the_intermediate_too(self, beta_carotene):
        slow = solve_pathway(beta_carotene, 1.0085e-3, 0.11, BETA_CAROTENE_KINETICS)
        fast = solve_pathway(beta_carotene, 1.0085e-3, 0.25, BETA_CAROTENE_KINETICS)

        assert (fast.node("lycopene").content_mmol_per_gdcw
                < slow.node("lycopene").content_mmol_per_gdcw)

    def test_more_flux_gives_more_product_at_a_fixed_growth_rate(self, beta_carotene):
        low = solve_pathway(beta_carotene, 5e-4, 0.15, BETA_CAROTENE_KINETICS)
        high = solve_pathway(beta_carotene, 2e-3, 0.15, BETA_CAROTENE_KINETICS)

        assert (high.terminal.content_mmol_per_gdcw
                > low.terminal.content_mmol_per_gdcw)


class TestWhatTheSolverRefuses:
    def test_a_non_growing_culture_is_refused(self, beta_carotene):
        """At mu = 0 every pool here diverges, because washout is the only outlet.

        Returning the limit -- or a large number -- would read as a prediction about a
        stationary-phase culture, which is precisely the regime this arithmetic does not
        describe.
        """
        with pytest.raises(ValueError, match="no steady state"):
            solve_pathway(beta_carotene, 1e-3, 0.0, BETA_CAROTENE_KINETICS)

    def test_a_negative_growth_rate_is_refused(self, beta_carotene):
        with pytest.raises(ValueError, match="growth_rate must be positive"):
            solve_pathway(beta_carotene, 1e-3, -0.1, BETA_CAROTENE_KINETICS)

    def test_a_negative_flux_is_refused(self, beta_carotene):
        with pytest.raises(ValueError, match="entry_flux must be non-negative"):
            solve_pathway(beta_carotene, -1e-6, 0.15, BETA_CAROTENE_KINETICS)

    def test_a_saturating_node_with_no_parameters_names_the_node(self, beta_carotene):
        """A refusal that does not name the node sends the caller reading the whole spec.

        The kinetics dictionary is keyed by node name and defaults to empty, so this is
        also what a renamed node looks like from the solver's side.
        """
        with pytest.raises(ValueError, match="lycopene"):
            solve_pathway(beta_carotene, 1e-3, 0.15, kinetics={})

    def test_a_saturating_node_missing_only_its_km_is_still_refused(self):
        """Half a Michaelis-Menten pair is not a rate law, and the half that is present
        would otherwise define the whole step through the dataclass default.
        """
        with pytest.raises(ValueError, match="needs both"):
            solve_pathway(_two_nodes(RateLaw.SATURATING), 1e-3, 0.15,
                          {"presqualene": NodeKinetics(vmax_per_growth=0.0023252)})

    def test_a_non_positive_vmax_is_refused(self):
        with pytest.raises(ValueError, match="must be positive"):
            solve_pathway(_two_nodes(RateLaw.SATURATING), 1e-3, 0.15,
                          {"presqualene": NodeKinetics(vmax_per_growth=0.0, km=5.974e-4)})

    def test_a_proportional_node_with_no_rate_constant_names_the_node(self):
        with pytest.raises(ValueError, match="presqualene"):
            solve_pathway(_two_nodes(RateLaw.PROPORTIONAL), 1e-3, 0.15, kinetics={})

    def test_an_unknown_node_name_fails_by_name(self, beta_carotene):
        solution = solve_pathway(beta_carotene, 1e-3, 0.15, BETA_CAROTENE_KINETICS)

        with pytest.raises(KeyError, match="phytofluene"):
            solution.node("phytofluene")


class TestItRefusesOutsideTheRangeItsKineticsWereFittedOver:
    """A guard that existed, was dropped, and had to come back.

    The hand-written carotenoid solver refused outside [0.101, 0.254] /h. Replacing it with
    this generic one dropped that, because `NodeKinetics` had no notion of a fitted range --
    and the chain then answered at **D = 0.05**, half the lower bound, with no warning of any
    kind. That is the rate a third chemostat point would most usefully sit at, so the silent
    answer was waiting exactly where somebody would go looking.

    Refusing rather than warning is deliberate. A warning attached to a returned number is
    read by whoever is looking and by nobody else, and the number reaches the table either
    way. The rate a user most wants -- batch, mu around 0.4 -- is the one furthest outside
    the fit, which is precisely when a soft signal is least likely to stop anyone.
    """

    @staticmethod
    def _kinetics(low=0.101, high=0.254):
        return {"lycopene": NodeKinetics(vmax_per_growth=2.3252e-3, km=5.974e-4,
                                         growth_rate_range=(low, high))}

    def test_a_rate_below_the_fitted_window_is_refused(self, beta_carotene):
        with pytest.raises(ValueError, match="outside the range"):
            solve_pathway(beta_carotene, 6.0e-4, 0.05, self._kinetics())

    def test_a_rate_above_it_is_refused_too(self, beta_carotene):
        """Batch. The single most likely thing a user asks for."""
        with pytest.raises(ValueError, match="outside the range"):
            solve_pathway(beta_carotene, 6.0e-4, 0.40, self._kinetics())

    def test_a_rate_inside_it_is_answered(self, beta_carotene):
        solved = solve_pathway(beta_carotene, 6.0e-4, 0.18, self._kinetics())

        assert solved.terminal.content_mmol_per_gdcw > 0

    def test_the_boundaries_themselves_are_inside(self, beta_carotene):
        """Half-open would make the fitted endpoints unanswerable, which would be absurd --
        they are the two rates every parameter was measured at."""
        for rate in (0.101, 0.254):
            assert solve_pathway(beta_carotene, 6.0e-4, rate, self._kinetics()).entry_flux > 0

    def test_the_refusal_names_the_range_and_the_rate(self, beta_carotene):
        """A refusal a reader cannot act on is an obstacle. Which rate, and against what."""
        with pytest.raises(ValueError, match=r"0\.05.*\[0\.101, 0\.254\]"):
            solve_pathway(beta_carotene, 6.0e-4, 0.05, self._kinetics())

    def test_kinetics_with_no_declared_range_still_answer(self, beta_carotene):
        """A rate law fitted over no stated range is a real situation, and inventing a fake
        one to satisfy the guard would be worse than recording the absence."""
        unranged = {"lycopene": NodeKinetics(vmax_per_growth=2.3252e-3, km=5.974e-4)}

        assert solve_pathway(beta_carotene, 6.0e-4, 0.90, unranged).entry_flux > 0

    def test_the_shipped_calibration_declares_its_range(self, beta_carotene):
        """The whole point. An unranged shipped constant would pass every test above and
        still answer at 0.05."""
        from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS

        assert BETA_CAROTENE_KINETICS["lycopene"].growth_rate_range is not None

    def test_and_it_is_the_range_the_fit_was_made_over(self, beta_carotene):
        from ystwin.kinetic.carotenoid import ELIZONDO2025
        from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS

        assert (BETA_CAROTENE_KINETICS["lycopene"].growth_rate_range
                == ELIZONDO2025.growth_rate_range)
