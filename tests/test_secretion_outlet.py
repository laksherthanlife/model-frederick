"""The export outlet: a terminal node that leaves through a transporter.

`Fate.SECRETED` used to be refused by `PathwaySpec` at load, which made every secreted
product unrepresentable rather than under-parameterised. That was the same mistake the
`DEGRADED` refusal made and it is corrected the same way: a pathway whose chemistry this
walk cannot describe is refused at LOAD, and a pathway one measured number short is refused
at SOLVE, by name. A secreted TERMINAL node is the second kind. A secreted INTERMEDIATE is
still the first, and `tests/test_pathway_spec.py` keeps that half.

WHY SATURATING AND NOT FIRST ORDER, which is the one scientific decision in this file.
`NodeKinetics.degradation_rate_per_h` is deliberately first order, and its docstring argues
that inventing two constants to avoid approximating with one is the worse trade. That
argument does not carry over, and the reason is a measurement rather than a preference.

Kastberg 2025 (PMID 39971732) ran five K. phaffii production strains in a glucose-limited
chemostat at D = 0.1 /h and measured the intracellular proteome and the secretome ON THE
SAME SAMPLES. The human insulin precursor rises 4.12 log2FC -- 17.4x -- inside the cell
between I1G and I1S while the secretome does not move at all; Mambalgin-1 rises 4.15 log2FC
inside against 1.8 log2FC secreted. A first-order outlet secretes in proportion to the pool
BY CONSTRUCTION, so no value of `k` reproduces a 17x pool beside a flat supernatant. That
pair is a direct refutation of `v_sec = k*X`, and `test_no_first_order_outlet_can_reproduce
_the_kastberg_pair` makes it executable.

WHAT THIS DOES NOT DO. Not one shipped pathway gains a number. No `data/pathways/*.toml`
declares a secreted node and none honestly can: Kastberg itself contains no absolute
quantity of product anywhere -- every product value in it is a relative label-free log2 fold
change -- so it motivates the shape and cannot calibrate it. The outlet ships as a capacity
and a refusal that names the experiment, which is what `test_no_secretion_constant_is
_shipped_anywhere` pins.
"""

from __future__ import annotations

import math

import pytest

from ystwin.pathway.solve import (
    NodeKinetics,
    content_ceiling,
    solve_pathway,
)
from ystwin.pathway.spec import Fate, Node, PathwaySpec, RateLaw

MOLAR_MASS = 50.0
MU = 0.1
WINDOW = (0.05, 0.20)

# The operating point that reproduces Kastberg's measured pair. Chosen so both states sit
# above the half-saturation content, which is what "the exporter is at capacity" means.
#
# THESE ARE ILLUSTRATIVE AND THE ONLY REAL PAIRED DATA SITS SOMEWHERE ELSE ENTIRELY.
# Recorded 2026-09-02, when Adelantado 2017 (PMID 28526017) was fitted for the first time --
# four glucose-limited chemostat states of K. phaffii Fab 2F5 at D ~ 0.095, the one published
# dataset carrying paired intracellular content and supernatant rate in absolute mass. Least
# squares on v = vmax*X/(km + X) over its (X = cytosolic + membrane, v = qFab) pairs returns
# vmax = 0.153 mg/gDCW/h and km = 0.149 mg/gDCW, with standard errors of 489% and 619% and
# R^2 = 0.29. Three things follow, and the third contradicts the sentence above:
#
#   1. The pairing is NOT MONOTONIC. X = 0.038 secretes 0.049 while the larger X = 0.051
#      secretes 0.033 (Pearson r = 0.53). No monotone saturating law fits all four points, so
#      vmax would have to move with oxygenation -- which confounds it with km.
#   2. The constants are ESTIMABLE IN PRINCIPLE AND UNIDENTIFIED IN FACT. The pool spans 1.8x,
#      nowhere near enough curvature to separate the two. Any fit must ship the ridge.
#   3. The fitted km sits 2.9x ABOVE the largest measured pool, so every real state is on the
#      LINEAR limb -- the opposite of the regime these constants illustrate. That means this
#      dataset CANNOT distinguish a saturating outlet from a first-order one, which is the one
#      structural choice `NodeKinetics.secretion_km_mmol_per_gdcw` says it made on evidence.
#      The Kastberg pair still refutes first-order EXPORT, because a 17x pool beside a flat
#      supernatant is not something k*X can produce -- but that is a refutation of the
#      alternative, not a confirmation of this one, and the difference should not be blurred.
#
# In this module's units (Fab ~ 48 kDa) the measured pair is vmax ~ 3.2e-3 mmol/gDCW/h and
# km ~ 3.1e-3 mmol/gDCW. The values below are 125x and 0.3x those, and are kept because the
# tests exercise ARITHMETIC and limits, not a calibration. Nothing here is a fitted constant.
VMAX_SEC = 0.4
KM_SEC = 1e-3


def _chain(fate: str, law: str = RateLaw.PASSTHROUGH) -> PathwaySpec:
    """Two nodes, so the terminal arm is exercised with a feeder above it."""
    return PathwaySpec(
        product="prod", organism="S. cerevisiae", entry_enzyme="E0",
        precursor_metabolite="p",
        nodes=(Node("n0", law, Fate.DILUTED, "E0"),
               Node("prod", RateLaw.PASSTHROUGH, fate, "",
                    molar_mass_g_per_mol=MOLAR_MASS)))


def _export(vmax: float = VMAX_SEC, km: float = KM_SEC,
            window: tuple[float, float] | None = WINDOW) -> dict[str, NodeKinetics]:
    return {"prod": NodeKinetics(secretion_vmax_mmol_per_gdcw_h=vmax,
                                 secretion_km_mmol_per_gdcw=km,
                                 growth_rate_range=window)}


class TestTheRefusalMovedRatherThanLifted:
    """The load-time refusal is gone for a terminal node; the guarantee it protected is not."""

    def test_a_secreted_terminal_node_now_loads(self):
        spec = _chain(Fate.SECRETED)

        assert spec.nodes[-1].fate == Fate.SECRETED

    def test_solving_without_an_export_capacity_is_refused(self):
        with pytest.raises(ValueError, match="secretion_vmax_mmol_per_gdcw_h"):
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU)

    def test_the_refusal_still_names_the_dilution_term_that_is_missing(self):
        """The exact string the load-time refusal carried, and the reason is unchanged.

        `tests/test_pathway_spec.py` and `tests/test_tal_pathway.py` both matched on
        `mu*[X]` because the refusal is arithmetic and not policy. Moving where it is raised
        must not cost the reader the sentence that says so.
        """
        with pytest.raises(ValueError) as raised:
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU)

        assert "mu*[X]" in str(raised.value)

    def test_the_refusal_names_the_paired_experiment_and_not_just_the_constants(self):
        """A refusal naming a missing number but not how to get it sends the reader nowhere.

        Two constants cannot be had from one steady state -- it fixes only the ratio
        `vmax/(km + X)` -- so the message has to say TWO OR MORE, or a wet lab will run the
        single-state experiment and come back with an unidentifiable pair.
        """
        with pytest.raises(ValueError) as raised:
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU)
        message = str(raised.value)

        assert "supernatant accumulation rate" in message
        assert "TWO OR MORE" in message

    def test_the_refusal_says_how_large_the_omitted_outlet_can_be(self):
        """A refusal without a magnitude reads as pedantry. Pfeffer 2011 supplies one."""
        with pytest.raises(ValueError, match=r"5\.5x the dilution flux"):
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU)

    @pytest.mark.parametrize("kinetics", [
        NodeKinetics(secretion_vmax_mmol_per_gdcw_h=VMAX_SEC, growth_rate_range=WINDOW),
        NodeKinetics(secretion_km_mmol_per_gdcw=KM_SEC, growth_rate_range=WINDOW),
    ])
    def test_one_of_the_two_constants_is_not_enough(self, kinetics):
        """Half a parameterisation must refuse, because the ridge is the whole problem."""
        with pytest.raises(ValueError, match="secretion_"):
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU, {"prod": kinetics})

    def test_an_export_capacity_without_a_growth_rate_range_is_refused(self):
        """Stricter than the degradation outlet, and the asymmetry is dimensional.

        `degradation_rate_per_h` is first order and carries its own `mu` implicitly through
        the balance. An export vmax is an ABSOLUTE flux: it does not scale with growth, so
        applying one at a dilution rate it was not measured at is a category error rather
        than an extrapolation. Rebnegger 2024 measured qp falling 34-fold across the range.
        """
        with pytest.raises(ValueError, match="growth_rate_range"):
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU, _export(window=None))

    @pytest.mark.parametrize("bad", [math.inf, -math.inf, math.nan, 0.0, -1.0])
    def test_a_non_finite_or_non_positive_capacity_is_refused(self, bad):
        with pytest.raises(ValueError, match="finite and positive"):
            solve_pathway(_chain(Fate.SECRETED), 0.5, MU, _export(vmax=bad))

    def test_export_constants_on_a_diluted_node_are_refused_rather_than_ignored(self):
        """A field the solver reads and drops is a claim the reader believes.

        `NodeKinetics` is frozen with no `__post_init__`, so the point of use is the only
        validation site there is. Silently returning the diluted answer while the caller
        believes an export was modelled is the worst of the available outcomes.
        """
        with pytest.raises(ValueError, match="only a SECRETED node reads"):
            solve_pathway(_chain(Fate.DILUTED), 0.5, MU, _export())


class TestTheDilutedAnswerIsUntouched:
    """The inertness claim, asserted with `==` rather than `approx` throughout.

    The claim is exactness, not agreement: the terminal arm branches on `None` before any
    arithmetic rather than passing a zero capacity through the quadratic, precisely so that
    an existing result cannot move in its last bits.
    """

    def test_a_diluted_terminal_is_bit_for_bit_flux_over_mu(self):
        solved = solve_pathway(_chain(Fate.DILUTED), 0.5, MU)

        assert solved.terminal.content_mmol_per_gdcw == 0.5 / MU

    def test_a_diluted_node_carries_a_literal_zero_loss_flux(self):
        """So the `carbon_closes` sum is exact rather than merely small."""
        solved = solve_pathway(_chain(Fate.DILUTED), 0.5, MU)

        assert all(n.loss_flux_mmol_per_gdcw_h == 0.0 for n in solved.nodes)

    def test_a_degraded_terminal_is_bit_for_bit_flux_over_the_summed_loss(self):
        """The other pre-existing outlet, which also runs through the new branch."""
        solved = solve_pathway(_chain(Fate.DEGRADED), 0.5, MU,
                               {"prod": NodeKinetics(degradation_rate_per_h=0.4)})

        assert solved.terminal.content_mmol_per_gdcw == 0.5 / (MU + 0.4)

    def test_the_shipped_beta_carotene_ceiling_is_unchanged(self, beta_carotene=None):
        """`content_ceiling` gained a fate check; the shipped spec must not notice."""
        from ystwin.pathway.calibrations import BETA_CAROTENE_KINETICS
        from ystwin.pathway.spec import load_pathway

        spec = load_pathway("beta_carotene")

        assert content_ceiling(spec, BETA_CAROTENE_KINETICS) == 2.3252e-3


class TestTheArithmeticOfTheOutlet:
    def test_the_solved_root_satisfies_the_balance_it_was_derived_from(self):
        """The quadratic is a rearrangement, so the original balance is the real check."""
        for flux in (1e-3, 1e-2, 0.1, 0.5, 1.0):
            for mu in (0.05, 0.1, 0.2):
                solved = solve_pathway(_chain(Fate.SECRETED), flux, mu, _export())
                x = solved.terminal.content_mmol_per_gdcw
                residual = flux - mu * x - VMAX_SEC * x / (KM_SEC + x)

                assert abs(residual) <= 1e-12 * flux, f"flux={flux} mu={mu}"

    def test_carbon_closes_on_a_secreting_chain(self):
        """The identity, on the case that would have broken it."""
        solved = solve_pathway(_chain(Fate.SECRETED), 0.5, MU, _export())

        assert solved.terminal.loss_flux_mmol_per_gdcw_h > 0.0
        assert abs(solved.carbon_closes) <= 1e-12 * solved.entry_flux

    def test_carbon_closes_on_a_degrading_chain(self):
        """THE MISSING TEST, and it failed before this change.

        `carbon_closes` subtracted only `mu*[X]` while its docstring called the result an
        identity that must be zero. From the day the DEGRADED outlet began solving it
        returned that outlet's flux instead: 0.4 here, which is 80% of the entry flux,
        as an ordinary float with no flag. Every existing test of the identity ran on
        DILUTED chains, where the omitted term is zero.
        """
        solved = solve_pathway(_chain(Fate.DEGRADED), 0.5, MU,
                               {"prod": NodeKinetics(degradation_rate_per_h=0.4)})

        assert solved.carbon_closes == 0.0

    def test_a_vanishing_capacity_degenerates_to_the_diluted_answer(self):
        """The outlet is a strict generalisation, checked at the limit rather than argued."""
        solved = solve_pathway(_chain(Fate.SECRETED), 0.5, MU, _export(vmax=1e-300))

        assert solved.terminal.content_mmol_per_gdcw == pytest.approx(0.5 / MU, rel=1e-12)

    def test_a_free_transporter_empties_the_cell(self):
        """The other limit: infinite capacity secretes everything and holds nothing."""
        solved = solve_pathway(_chain(Fate.SECRETED), 0.5, MU, _export(vmax=1e9))

        assert solved.terminal.content_mmol_per_gdcw < 1e-8
        assert solved.terminal.loss_flux_mmol_per_gdcw_h == pytest.approx(0.5, rel=1e-9)

    def test_the_stable_branch_beats_the_naive_root_where_it_matters(self):
        """Why this arm does not reuse the SATURATING arm's expression verbatim.

        `(-b + sqrt(D))/(2a)` subtracts nearly equal numbers whenever `b > 0`, and for this
        outlet `b > 0` is the UNSATURATED regime -- a strain below its export capacity --
        which is the common case and the one a caller most wants. The SATURATING arm is
        deliberately left alone: it has the same weakness, but at the shipped calibration
        the two branches agree to 8.2e-16 relative, and rewriting it would move committed
        numbers to fix nothing that is biting.
        """
        mu, vmax, km, flux = 0.15, 1e-2, 1e-18, 1e-3
        a = mu
        b = mu * km + vmax - flux
        c = -flux * km
        naive = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
        solved = solve_pathway(_chain(Fate.SECRETED), flux, mu,
                               _export(vmax=vmax, km=km)).terminal

        assert naive == 0.0
        assert solved.content_mmol_per_gdcw > 0.0
        residual = flux - mu * solved.content_mmol_per_gdcw - (
            vmax * solved.content_mmol_per_gdcw / (km + solved.content_mmol_per_gdcw))
        assert abs(residual) <= 1e-12 * flux


class TestTheKastbergSignature:
    """The measurement that chose the rate law, as an executable claim.

    Kastberg 2025 measured a 4.12 log2FC intracellular rise in the human insulin precursor
    against a secretome that did not move. These two tests assert that a saturating outlet
    produces that shape and that no first-order outlet can.
    """

    LOW_CONTENT = 0.05
    HIGH_CONTENT = 0.87

    def _flux_for(self, content: float) -> float:
        return MU * content + VMAX_SEC * content / (KM_SEC + content)

    def test_past_the_capacity_the_pool_climbs_and_the_export_does_not(self):
        low = solve_pathway(_chain(Fate.SECRETED), self._flux_for(self.LOW_CONTENT),
                            MU, _export()).terminal
        high = solve_pathway(_chain(Fate.SECRETED), self._flux_for(self.HIGH_CONTENT),
                             MU, _export()).terminal

        content_log2 = math.log2(high.content_mmol_per_gdcw / low.content_mmol_per_gdcw)
        export_log2 = math.log2(high.loss_flux_mmol_per_gdcw_h
                                / low.loss_flux_mmol_per_gdcw_h)

        # Kastberg's hIP pair: +4.12 log2FC intracellular, secretome not significant.
        assert content_log2 == pytest.approx(4.12, abs=0.01)
        assert export_log2 < 0.05

    def test_no_first_order_outlet_can_reproduce_that_pair(self):
        """The falsification of the alternative design, rather than an argument for this one.

        For `v_sec = k*X` the export ratio IS the content ratio, for every `k` and at every
        operating point. A 17.4x pool beside a flat supernatant is therefore not a fit that
        a first-order law fits badly; it is one the law cannot express.
        """
        ratio = self.HIGH_CONTENT / self.LOW_CONTENT

        for k in (1e-6, 1e-3, 0.1, 1.0, 1e3, 1e6):
            assert (k * self.HIGH_CONTENT) / (k * self.LOW_CONTENT) == pytest.approx(ratio)

        assert ratio == pytest.approx(17.4, abs=0.05)


class TestTheChangeIsCapabilityOnly:
    def test_no_shipped_pathway_declares_a_secreted_node(self):
        """A spec must not arrive with the code. Kastberg carries no absolute titre at all,
        so nothing in reach could calibrate one honestly."""
        from ystwin.pathway.spec import available_pathways, load_pathway

        for name in available_pathways():
            spec = load_pathway(name)

            assert all(n.fate != Fate.SECRETED for n in spec.nodes), name

    def test_no_secretion_constant_is_shipped_anywhere(self):
        """The Tier 0 guard: a default smuggled into a calibration would be a measurement
        nobody made."""
        from ystwin.pathway import calibrations

        for name in dir(calibrations):
            value = getattr(calibrations, name)
            if not isinstance(value, dict):
                continue
            for parameters in value.values():
                if isinstance(parameters, NodeKinetics):
                    assert parameters.secretion_vmax_mmol_per_gdcw_h is None
                    assert parameters.secretion_km_mmol_per_gdcw is None


class TestTheCeilingDeclinesRatherThanOverstating:
    def test_the_ceiling_is_none_for_a_secreting_terminal(self):
        """`mu` cancels in that bound only when dilution is the sole outlet.

        With a third outlet `X <= vmax_per_growth * mu/(mu + k_loss)`, which is strictly
        increasing in `mu` and reaches `vmax_per_growth` only in the limit. The old return
        stayed a valid upper bound and stopped being the asymptote the docstring calls it,
        which is the kind of half-true a caller acts on.
        """
        spec = _chain(Fate.SECRETED, law=RateLaw.SATURATING)
        kinetics = {"n0": NodeKinetics(vmax_per_growth=2.3252e-3, km=5.974e-4)}

        assert content_ceiling(spec, kinetics) is None

    def test_the_ceiling_is_none_for_a_degrading_terminal_too(self):
        spec = _chain(Fate.DEGRADED, law=RateLaw.SATURATING)
        kinetics = {"n0": NodeKinetics(vmax_per_growth=2.3252e-3, km=5.974e-4)}

        assert content_ceiling(spec, kinetics) is None


class TestTheProductionRateRefusesRatherThanUnderReporting:
    def test_content_times_growth_is_refused_for_a_lossy_terminal(self):
        """Two defensible answers differ by the degraded fraction, so the scalar refuses.

        On Pfeffer 2011's chemostat the washout term is 6.4% of what the cell makes, so a
        property that returned `content * mu` for a secreted product would be sixteen times
        too small -- and three scripts write that number into committed CSVs.
        """
        from ystwin.predict import ProductPrediction

        assert hasattr(ProductPrediction, "rate_mmol_per_gdcw_h")
