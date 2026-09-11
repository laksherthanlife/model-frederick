"""The GEM asked for a yield, not a ceiling -- and the reason that distinction was needed.

Two workflows concluded the GEM could not carry the environment signal. Both asked it for a
CEILING, and that failure is real: the predicted flux sits 12.8-195x below the ceiling, and
scaling by each state's measured uptake gives a captured fraction spanning 216x. The
conclusion drawn from it was not.

Asked for a YIELD at matched carbon and matched growth rate, Yeast9 ranks the carbon sources
correctly on the one axis where a measurement exists. These tests pin that, pin the formulation
sensitivity that made the earlier answer unstable, and pin the committed prediction for
beta-carotene so it can be falsified by the planned bench test rather than quietly revised.

Marked `integration` throughout: every test here solves the real GSMM.
"""

from __future__ import annotations

import math
import pathlib

import pytest

from ystwin.pathway.gem_environment import (
    DEFAULT_FORMULATION,
    MIN_GROWTH_HEADROOM,
    Formulation,
    InfeasibleState,
    environment_descriptor,
    sign_is_stable,
    yield_ratio,
)
from ystwin.pathway.spec import load_pathway

pytestmark = pytest.mark.integration

_OUTPUTS = pathlib.Path(__file__).resolve().parents[1] / "outputs"


def _cell(column: str, law: str) -> float:
    """One committed number from `outputs/gem_environment_scores.csv`.

    Tests in this class quoted precise values in their docstrings under assertions loose
    enough to admit an 80% move -- `< 0.60`, `> 0.0`, `< 2.0` -- so two of them went stale
    across the precursor-carrier fix and stayed green. Comparing against the cell instead is
    what makes the docstring and the assertion the same claim.
    """
    import pandas as pd

    frame = pd.read_csv(_OUTPUTS / "gem_environment_scores.csv")
    return float(frame.loc[frame.law == law, column].iloc[0])


def _prediction_cell(column: str, product: str, growth_rate: float) -> float:
    """The same, from `outputs/gem_environment_predictions.csv`."""
    import pandas as pd

    frame = pd.read_csv(_OUTPUTS / "gem_environment_predictions.csv")
    row = frame[(frame["product"] == product)
                & (frame.growth_rate_per_h == growth_rate)]
    return float(row[column].iloc[0])


@pytest.fixture(scope="module")
def phb():
    return load_pathway("phb")


@pytest.fixture(scope="module")
def carotene():
    return load_pathway("beta_carotene")


class TestTheGemRanksTheMeasuredAxisCorrectly:
    """The finding that reopened the question. Kocharin measured ethanol 3.82x glucose for
    PHB content at mu = 0.05 /h; the GEM has to at least agree on the sign."""

    def test_ethanol_beats_glucose_for_phb(self, yeast_gem, phb):
        assert yield_ratio(yeast_gem, phb, 0.05).favours_ethanol

    @pytest.mark.parametrize("growth_rate", [0.05, 0.10, 0.15])
    def test_at_every_growth_rate_measured(self, yeast_gem, phb, growth_rate):
        assert yield_ratio(yeast_gem, phb, growth_rate).ratio > 1.0

    def test_the_ratio_is_a_real_yield_difference_not_a_rounding_artefact(self, yeast_gem,
                                                                         phb):
        result = yield_ratio(yeast_gem, phb, 0.05)

        assert result.glucose_max_mmol_per_gdcw_h > 0
        assert result.ethanol_max_mmol_per_gdcw_h > result.glucose_max_mmol_per_gdcw_h * 1.05

    def test_it_understates_the_measured_effect_rather_than_matching_it(self, yeast_gem, phb):
        """Honesty check, and it must not silently become false. The GEM reaches about 30% of
        the measured log effect; a version that suddenly matched 3.82x would mean the
        formulation had drifted into something that is fitting rather than predicting."""
        ratio = yield_ratio(yeast_gem, phb, 0.05).ratio

        assert 1.0 < ratio < 2.0


class TestTheGemDoesNotReproduceTheReversal:
    """The known limit, pinned so nobody claims it. Kocharin's ethanol advantage SHRINKS with
    growth rate and flips at mu = 0.1718 /h. The GEM's grows. That gap is the boundary between
    the two halves of the vision -- what is possible versus how much a cell takes."""

    def test_the_gem_advantage_grows_with_growth_rate(self, yeast_gem, phb):
        low = yield_ratio(yeast_gem, phb, 0.05).ratio
        high = yield_ratio(yeast_gem, phb, 0.15).ratio

        assert high > low

    def test_which_is_the_opposite_of_the_measurement(self):
        from ystwin.pathway.environment_flux import environment_factor

        measured_low = environment_factor("phb", 1.0, 0.05)
        measured_high = environment_factor("phb", 1.0, 0.20)

        assert measured_high < measured_low

    def test_the_gem_never_flips_sign_inside_the_measured_window(self, yeast_gem, phb):
        """So it cannot supply the reversal, and a caller must not expect it to."""
        for growth_rate in (0.05, 0.10, 0.15):
            assert yield_ratio(yeast_gem, phb, growth_rate).ratio > 1.0


class TestTheCommittedPredictionForBetaCarotene:
    """Stoichiometric, so it needs no data and can be wrong. That is what makes it a
    prediction. The planned matched-growth-rate bench test is what settles it."""

    @pytest.mark.parametrize("growth_rate", [0.05, 0.10, 0.15])
    def test_beta_carotene_is_predicted_to_share_phbs_sign(self, yeast_gem, carotene,
                                                           growth_rate):
        assert yield_ratio(yeast_gem, carotene, growth_rate).favours_ethanol

    def test_the_prediction_needs_no_calibration(self, yeast_gem, carotene):
        """beta_carotene has an entry-flux calibration but this does not read it; phb has
        none at all and still gets a descriptor. That is the transfer the vision needs."""
        result = yield_ratio(yeast_gem, load_pathway("gadusol"), 0.10)

        assert result.product == "gadusol"
        assert result.ratio > 0

    def test_only_the_two_products_with_data_have_a_stable_sign(self, yeast_gem, carotene):
        """The claim that replaced a withdrawn one. This module used to advertise that
        "gadusol and glycogen flip at low growth rate", which was an artefact of a matched
        carbon supply shipping at 4.7x Kocharin's own measured uptake. Swept across the
        plausible range those two cross 1.0 and are not predictions; beta-carotene and PHB
        do not and are.

        AT mu = 0.05, which the 2026-09-04 fix made explicit rather than incidental."""
        from ystwin.pathway.gem_environment import sign_is_stable

        for product in ("beta_carotene", "phb"):
            sweep = sign_is_stable(yeast_gem, load_pathway(product), 0.05)
            assert sweep.stable, f"{product} sign should be stable"
            assert sweep.is_a_prediction(), sweep.summary()
            assert sweep.ratio_low > 1.0

        for product in ("gadusol", "glycogen"):
            sweep = sign_is_stable(yeast_gem, load_pathway(product), 0.05)
            assert not sweep.stable, f"{product} sign should NOT be a prediction"
            assert not sweep.is_a_prediction()
            assert sweep.ratio_low < 1.0 < sweep.ratio_high


class TestASweepReportsWhatItActuallyCovered:
    """`sign_is_stable` returns its own provenance, added 2026-09-04.

    The second instance of the cf34701 finding in this module. `sign_is_stable` was written
    to repair a claim quantified over an unswept SUPPLY, and reproduced the same shape one
    level down: it swept supply, left the GROWTH RATE wherever the caller put it, returned a
    bare `(stable, low, high)` that recorded nothing about coverage, and `continue`d past
    every infeasible point. A verdict reached on four of six points -- the four nearest the
    feasibility boundary -- was written out as a committed prediction.

    These tests pin the provenance rather than today's two demoted rows, because the next
    forgotten parameter will not be the growth rate.
    """

    def test_a_complete_sweep_says_so(self, yeast_gem, phb):
        sweep = sign_is_stable(yeast_gem, phb, 0.05)

        assert sweep.points_feasible == sweep.points_requested == 6
        assert sweep.swept_the_whole_range
        assert sweep.min_growth_headroom > 0.0

    def test_an_incomplete_sweep_is_refused_as_a_prediction_even_when_stable(
            self, yeast_gem, phb):
        """The exact shipped defect. At mu = 0.15 every ratio PHB can still be evaluated at
        is above 1.0, so `stable` is True -- and the sweep only reached four of the six
        supplies it asked for, so it is not a prediction."""
        sweep = sign_is_stable(yeast_gem, phb, 0.15)

        assert sweep.stable
        assert sweep.points_feasible < sweep.points_requested
        assert not sweep.is_a_prediction()
        assert "NOT a prediction" in sweep.summary()

    def test_a_growth_rate_near_the_feasibility_boundary_is_refused(self, yeast_gem, phb):
        """The other half of the gate, isolated from the coverage half. Swept over a narrow
        supply range where every point IS feasible, a growth rate close to the model's own
        maximum still fails, because the ratio there is set by the boundary and not by the
        chemistry."""
        from ystwin.pathway.gem_environment import Formulation, max_growth_rate

        supply = 11.0e-3
        ceiling = max_growth_rate(
            yeast_gem, "glucose", Formulation(matched_carbon_cmol_per_gdcw_h=supply))
        near = sign_is_stable(yeast_gem, phb, 0.97 * ceiling,
                              supply_range_cmol_per_gdcw_h=(supply, 1.02 * supply))
        far = sign_is_stable(yeast_gem, phb, 0.30 * ceiling,
                             supply_range_cmol_per_gdcw_h=(supply, 1.02 * supply))

        assert near.swept_the_whole_range and far.swept_the_whole_range
        assert near.min_growth_headroom < MIN_GROWTH_HEADROOM <= far.min_growth_headroom
        assert not near.is_a_prediction()
        assert far.is_a_prediction()
        assert near.ratio_high > 1.5 * far.ratio_high, (
            "the near-boundary ratio should be inflated, which is why it is gated")

    def test_a_never_feasible_product_is_not_reported_as_a_sign_reversal(
            self, yeast_gem, phb):
        """The conflation `InfeasibleState`'s own docstring forbids. This used to return
        `(False, nan, nan)`, which the caller wrote out as `is_a_prediction=False` under the
        banner "sign reverses inside that range" -- recording "the vessel cannot run" as
        "this feed is bad for this product"."""
        with pytest.raises(InfeasibleState, match="not an unstable sign"):
            sign_is_stable(yeast_gem, phb, 5.0)

    def test_the_growth_headroom_is_exactly_the_feasibility_boundary(self, yeast_gem, phb):
        """`max_growth_rate` and `yield_ratio` must agree about where the boundary is, or the
        gate is measuring a different constraint from the one that bites."""
        from ystwin.pathway.gem_environment import Formulation, max_growth_rate

        formulation = Formulation(matched_carbon_cmol_per_gdcw_h=11.0e-3)
        ceiling = max_growth_rate(yeast_gem, "glucose", formulation)

        yield_ratio(yeast_gem, phb, 0.98 * ceiling, formulation)
        with pytest.raises(InfeasibleState):
            yield_ratio(yeast_gem, phb, 1.02 * ceiling, formulation)

    def test_gadusol_and_glycogen_were_never_two_predictions(self, yeast_gem):
        """They agree to well under a percent, because the descriptor collapses to which
        side of pyruvate dehydrogenase the precursor sits on."""
        a = yield_ratio(yeast_gem, load_pathway("gadusol"), 0.05).ratio
        b = yield_ratio(yeast_gem, load_pathway("glycogen"), 0.05).ratio

        assert abs(a - b) / a < 0.02


class TestTheSecondFormulationIsActuallyComputed:
    """`precursor_sink_only=False` must produce a DIFFERENT number, or the flag is a label.

    THE ONE ASSERTION THIS SUITE LACKED, and it is what would have caught the defect on the
    day the field was written. `Formulation.precursor_sink_only` was declared, printed in
    `label()`, copied through `sign_is_stable` and used in the `comparable_with` equality
    contract -- and never read by `_max_precursor`. So two results carrying "precursor sink"
    and "full pathway" were bit-identical, and `comparable_with` refused to compare two
    results that were the same computation. The module docstring meanwhile quoted a
    formulation sensitivity (1.09-1.25 against 1.40-2.00) that nothing in the repository could
    produce, and it therefore could not go stale, because it had never been computed.
    """

    @pytest.fixture(scope="class")
    def full(self):
        return Formulation(precursor_sink_only=False)

    def test_the_two_formulations_disagree_on_a_product_with_a_declared_cofactor_cost(
            self, yeast_gem, phb, full):
        sink = yield_ratio(yeast_gem, phb, 0.05).ratio
        pathway = yield_ratio(yeast_gem, phb, 0.05, full).ratio

        assert pathway != pytest.approx(sink, rel=1e-6), (
            "the flag changes the label and must change the computation")

    def test_the_full_pathway_consumes_the_nadph_the_sink_cannot_see(self, yeast_gem, phb):
        """Which is WHY they differ. PhaB is NADPH-dependent; the bare precursor pull has no
        way to express that, and the spec's declared net balance does."""
        assert phb.full_pathway_stoichiometry["s_1212"] == -1.0
        assert phb.full_pathway_stoichiometry["s_0373"] == -phb.precursor_stoichiometry

    def test_a_spec_that_has_not_declared_its_balance_is_refused(self, yeast_gem, carotene,
                                                                 full):
        """Refused rather than falling back to the sink, which is the same refusal shape as
        `PrecursorCarrierUnknown`: a cofactor cost omitted is not a cofactor cost of zero, and
        a result labelled 'full pathway' that silently held the sink's number is the defect."""
        from ystwin.pathway.gem_environment import FullPathwayUndeclared

        assert carotene.full_pathway_stoichiometry is None
        with pytest.raises(FullPathwayUndeclared, match="no full_pathway_stoichiometry"):
            yield_ratio(yeast_gem, carotene, 0.05, full)

    def test_results_from_the_two_formulations_are_refused_as_incomparable(self, yeast_gem,
                                                                          phb, full):
        """`comparable_with` now refuses a comparison that is genuinely incomparable, where
        before it refused one between two identical computations."""
        sink = yield_ratio(yeast_gem, phb, 0.05)
        pathway = yield_ratio(yeast_gem, phb, 0.05, full)

        assert not sink.comparable_with(pathway)
        assert "precursor sink" in sink.formulation.label()
        assert "full pathway" in pathway.formulation.label()


class TestTheFormulationIsPinned:
    """The sensitivity that made the earlier answer unstable. An unpinned formulation has no
    defined result, so every number carries the choices that produced it."""

    def test_a_result_carries_its_formulation(self, yeast_gem, phb):
        assert yield_ratio(yeast_gem, phb, 0.10).formulation == DEFAULT_FORMULATION

    def test_results_from_different_formulations_are_not_comparable(self, yeast_gem, phb):
        wide = Formulation(matched_carbon_cmol_per_gdcw_h=40.0e-3)

        assert not yield_ratio(yeast_gem, phb, 0.10).comparable_with(
            yield_ratio(yeast_gem, phb, 0.10, wide))

    def test_the_matched_carbon_supply_really_does_move_the_answer(self, yeast_gem, phb):
        """Pins that the sensitivity is real rather than a story: if this ever stops being
        true, the formulation field has become decoration."""
        lean = yield_ratio(yeast_gem, phb, 0.10,
                           Formulation(matched_carbon_cmol_per_gdcw_h=12.0e-3)).ratio
        rich = yield_ratio(yeast_gem, phb, 0.10,
                           Formulation(matched_carbon_cmol_per_gdcw_h=40.0e-3)).ratio

        assert lean != rich

    def test_the_label_names_every_choice_that_matters(self):
        label = DEFAULT_FORMULATION.label()

        assert "matched" in label and "O2" in label and "sink" in label


class TestTheDescriptor:
    def test_pure_glucose_is_zero_by_construction(self, yeast_gem, carotene):
        assert environment_descriptor(yeast_gem, carotene, 0.0, 0.10) == 0.0

    def test_it_is_the_log_ratio_at_full_ethanol(self, yeast_gem_factory, carotene):
        descriptor = environment_descriptor(yeast_gem_factory(), carotene, 1.0, 0.10)
        ratio = yield_ratio(yeast_gem_factory(), carotene, 0.10).ratio

        assert descriptor == pytest.approx(math.log(ratio))

    def test_a_blend_sits_between_the_pure_feeds(self, yeast_gem, carotene):
        blend = environment_descriptor(yeast_gem, carotene, 0.677, 0.10)
        full = environment_descriptor(yeast_gem, carotene, 1.0, 0.10)

        assert 0.0 < blend < full

    @pytest.mark.parametrize("fraction", [-0.1, 1.1])
    def test_an_impossible_fraction_is_refused(self, yeast_gem, carotene, fraction):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            environment_descriptor(yeast_gem, carotene, fraction, 0.10)


class TestItRefusesRatherThanReturningZero:
    def test_an_unholdable_growth_rate_is_refused(self, yeast_gem, phb):
        """Zero is a legitimate descriptor value -- a product the host cannot make -- and
        infeasibility is a different claim. Conflating them would report 'this vessel cannot
        run' as 'this feed is bad for this product'."""
        starved = Formulation(matched_carbon_cmol_per_gdcw_h=1.0e-3)

        with pytest.raises(InfeasibleState, match="infeasible"):
            yield_ratio(yeast_gem, phb, 0.40, starved)

    def test_the_refusal_says_it_is_not_a_zero_yield(self, yeast_gem, phb):
        with pytest.raises(InfeasibleState, match="not read this as a zero yield"):
            yield_ratio(yeast_gem, phb, 0.40,
                        Formulation(matched_carbon_cmol_per_gdcw_h=1.0e-3))

    def test_a_non_positive_growth_rate_is_refused(self, yeast_gem, phb):
        with pytest.raises(ValueError, match="must be positive"):
            yield_ratio(yeast_gem, phb, 0.0)

    def test_the_model_is_not_mutated(self, yeast_gem, phb):
        """Every change happens inside `with model:`. A leaked bound would silently corrupt
        every later solve in the session -- the fixture is module-scoped."""
        before = yeast_gem.reactions.get_by_id("r_1714").lower_bound
        yield_ratio(yeast_gem, phb, 0.10)

        assert yeast_gem.reactions.get_by_id("r_1714").lower_bound == before
        assert "YSTWIN_PRECURSOR_SINK" not in [r.id for r in yeast_gem.reactions]


class TestTheMagnitudeIsNotPredicted:
    """The correction. An earlier version of this module's docstring claimed the GEM
    descriptor beat mu-only with zero fitted parameters. It does not, and the reason is
    collinearity: the absolute maximum precursor flux is 92% correlated with growth rate, so
    it beats mu by BEING mu. These tests pin the honest version so the claim cannot come back.
    """

    @pytest.fixture(scope="class")
    def _kocharin_template(self, yeast_gem_factory):
        yeast_gem = yeast_gem_factory()
        import numpy as np
        import pandas as pd

        from ystwin.pathway.environment_flux import carbon_descriptor

        path = pytest.importorskip("pathlib").Path(
            "data/phb/kocharin2013_chemostat_states.tsv")
        if not path.exists():
            pytest.skip("Kocharin states not present")
        frame = pd.read_csv(path, sep="\t")
        frame["z"] = [carbon_descriptor(g, e) for g, e
                      in zip(frame.feed_glucose_g_per_l, frame.feed_ethanol_g_per_l)]
        frame["log_content"] = np.log(frame.phb_mg_per_gdw)
        frame["log_mu"] = np.log(frame.mu_per_h)
        spec = load_pathway("phb")
        # THE SCORING FORMULATION, which is deliberately not the PREDICTION default, and this
        # fixture used the default until 2026-09-04. `scripts/gem_environment_predictions.py`
        # scores the model comparison at 20 mmol C/gDCW/h because that supply admits all
        # eleven states -- a comparison that drops states scores its laws on different data.
        # Scoring here at the 11.0 default silently held out the mu = 0.20 rows, so this
        # fixture reproduced NONE of the cells the tests below quote: -0.8519 against the
        # committed -0.8190, and -0.0059 against +0.3656. Nothing caught it because the
        # assertions were `> 0.60` and `< 0.60`. Same formulation, same numbers, now.
        scoring = Formulation(matched_carbon_cmol_per_gdcw_h=20.0e-3)
        absolute, contrast = [], []
        for growth_rate, share in zip(frame.mu_per_h, frame.z):
            try:
                result = yield_ratio(yeast_gem, spec, growth_rate, scoring)
            except InfeasibleState:
                absolute.append(float("nan"))
                contrast.append(float("nan"))
                continue
            blended = (result.glucose_max_mmol_per_gdcw_h * (1 - share)
                       + result.ethanol_max_mmol_per_gdcw_h * share)
            absolute.append(math.log(blended))
            contrast.append(share * math.log(result.ratio))
        frame["gem_absolute"] = absolute
        frame["gem_contrast"] = contrast
        return frame.dropna(subset=["gem_absolute", "gem_contrast"])

    @pytest.fixture
    def kocharin(self, _kocharin_template):
        return _kocharin_template.copy(deep=True)

    def test_the_absolute_maximum_is_almost_entirely_growth_rate(self, kocharin):
        """-0.8190 -- 82% collinear with growth rate, so any 'win' over mu is mostly mu.

        The docstring said "0.92 in magnitude" until 2026-09-04, which was the value from
        before the precursor-carrier fix (385b6b8), and the assertion was `> 0.60` -- loose
        enough to admit an 80% move, so the stale number sat under a test that could not see
        it. Both are now pinned to the committed cell.
        """
        import numpy as np

        correlation = np.corrcoef(kocharin.gem_absolute, kocharin.log_mu)[0, 1]

        assert correlation == pytest.approx(_cell("corr_with_log_mu", "log GEM max"),
                                            abs=0.01)

    def test_the_pure_carbon_contrast_is_a_weak_predictor(self, kocharin):
        """+0.3656 against the empirical descriptor's +0.5801, after the carrier fix roughly
        doubled it. Still a sign and not a magnitude: it scores 2.1942x
        leave-one-carbon-source-out against mu-only's 1.7237x, so a better signal is not yet a
        better predictor.

        These four numbers were CURRENT -- the carrier fix updated this docstring correctly --
        and the assertions were `gem < 0.60` and `empirical > 0.0`, which would have passed on
        almost any pair. A precise number under a vacuous assertion is what let two stale
        docstrings survive in this same class, so both now read the cell.
        """
        import numpy as np

        gem = np.corrcoef(kocharin.gem_contrast, kocharin.log_content)[0, 1]
        empirical = np.corrcoef(kocharin.z, kocharin.log_content)[0, 1]

        assert gem == pytest.approx(
            _cell("corr_with_log_content", "GEM carbon contrast"), abs=0.01)
        assert empirical == pytest.approx(
            _cell("corr_with_log_content", "z empirical"), abs=0.01)
        assert gem < empirical

    def test_the_gem_understates_the_measured_effect_by_a_lot(self, yeast_gem, phb):
        """The magnitude gap, stated as a ratio so it cannot be read as agreement.
        The fitted law gives 4.6794x at mu = 0.05 (raw measured contrast 3.822x); the GEM
        says 1.3584x.

        The docstring said "about 1.09x" until 2026-09-04. That was the value from before the
        precursor-carrier fix (385b6b8), and the function returns 1.3584 -- but the assertion
        was `predicted < 2.0`, which the stale number and the true one both satisfy, so
        nothing could see it. Pinned to the committed prediction cell now.
        """
        from ystwin.pathway.environment_flux import environment_factor

        predicted = yield_ratio(yeast_gem, phb, 0.05).ratio
        measured = environment_factor("phb", 1.0, 0.05)

        assert predicted == pytest.approx(
            _prediction_cell("ethanol_over_glucose", "phb", 0.05), abs=1e-3)
        # The FITTED law at z = 1, which is 4.679x -- not the raw contrast, which is 3.822x.
        # The old assertion `measured > 3.0` admitted either and the docstring named the raw
        # one, so the two were interchangeable here; they are different quantities and the
        # gap this test measures is against the law the chain would actually apply.
        assert measured == pytest.approx(4.6794, abs=0.01)
        assert measured / predicted > 2.0


class TestThePrecursorCarrierIsReturned:
    """The bug that deflated one precursor family fifteenfold and left the other untouched.

    A pathway pulling acetyl-CoA consumes the ACETYL and hands back the CoA -- thiolase
    condenses two acetyl-CoA into acetoacetyl-CoA plus free CoA. Draining it without
    returning the carrier forces de novo coenzyme A synthesis for every molecule, so the
    number reported is a CoA-BIOSYNTHESIS ceiling and not a carbon yield. GGPP carries no
    group and was unaffected, so the defect destroyed exactly the between-product contrast
    this module computes.
    """

    def test_a_precursor_with_no_declared_carrier_is_refused(self, yeast_gem):
        """Refused rather than defaulted, because "no carrier" is right for GGPP and wrong
        by 15.5x for acetyl-CoA, and nothing in a metabolite id says which."""
        from ystwin.pathway.gem_environment import PrecursorCarrierUnknown
        from ystwin.pathway.spec import load_pathway

        spec = load_pathway("phb")
        # The net cofactor balance goes too: it is written per mole of product AROUND a named
        # precursor, so a spec that changes the precursor and keeps the balance is refused by
        # `PathwaySpec` itself before it can reach this refusal.
        invented = type(spec)(**{**spec.__dict__, "precursor_metabolite": "s_0001",
                                 "full_pathway_stoichiometry": None}) \
            if hasattr(spec, "__dict__") else None
        if invented is None:
            pytest.skip("spec is not reconstructable in this build")

        with pytest.raises(PrecursorCarrierUnknown, match="carrier chemistry"):
            yield_ratio(yeast_gem, invented, 0.10)

    def test_returning_coa_raises_the_acetyl_coa_yield_about_fifteenfold(self, yeast_gem):
        """The measurement that found it: 0.1749 without the carrier, 2.7130 with it, at
        20 mmol C/gDCW/h on glucose at mu = 0.10."""
        import cobra

        def supply(with_carrier):
            model = yeast_gem.copy()
            sink = cobra.Reaction("PROBE")
            sink.lower_bound, sink.upper_bound = 0.0, 1000.0
            mets = {model.metabolites.get_by_id("s_0373"): -2.0}
            if with_carrier:
                mets[model.metabolites.get_by_id("s_0529")] = 2.0
            sink.add_metabolites(mets)
            model.add_reactions([sink])
            model.reactions.get_by_id("r_1714").lower_bound = -20.0 / 6.0
            model.reactions.get_by_id("r_1761").lower_bound = 0.0
            model.reactions.get_by_id("r_1992").lower_bound = -1000.0
            model.reactions.get_by_id("r_2111").bounds = (0.10, 0.10)
            model.objective = "PROBE"
            return model.optimize().objective_value

        assert supply(True) / supply(False) > 10.0

    def test_the_beta_carotene_prediction_is_unaffected(self, yeast_gem, carotene):
        """GGPP carries no group, so the committed falsifiable prediction survived the fix
        unchanged -- which is why it is still a prediction and not a retrofit."""
        for growth_rate in (0.05, 0.10, 0.15):
            assert yield_ratio(yeast_gem, carotene, growth_rate).ratio > 1.2

    def test_phb_now_predicts_a_larger_ethanol_advantage(self, yeast_gem, phb):
        """And it moves TOWARD the measured 3.82x rather than away: 1.09-1.25x before the
        fix, 1.32-1.61x after."""
        assert yield_ratio(yeast_gem, phb, 0.05).ratio > 1.3
