"""D1: can constraint-based modelling predict beta-carotene flux at all?

If the feasible range of the product flux at a realistic growth rate spans most of
its own maximum, FBA is not a product predictor -- it bounds the product rather
than forecasting it, and the carotenoid branch needs kinetics instead.
"""

import numpy as np
import pytest

from ystwin.fba.carotenoid import PRODUCT_DEMAND_ID
from ystwin.fba.fva import capacity_sweep, product_flux_range

pytestmark = pytest.mark.integration


def test_forcing_maximum_growth_leaves_essentially_no_room_for_product(carotenoid_model):
    rng = product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=1.0)

    assert rng.width < 1e-6


def test_relaxing_the_growth_requirement_widens_the_product_range(carotenoid_model):
    tight = product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=0.99)
    loose = product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=0.50)

    assert loose.width > tight.width


def test_the_minimum_is_zero_whenever_the_pathway_is_merely_optional(carotenoid_model):
    rng = product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=0.90)

    assert rng.minimum == pytest.approx(0.0, abs=1e-9)
    assert rng.relative_width == pytest.approx(1.0, abs=1e-6)


def test_an_uncapped_pathway_leaves_the_product_flux_completely_undetermined(carotenoid_model):
    """The decisive check: without a capacity bound the model predicts 'anything'."""
    rng = product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=0.90)

    assert rng.maximum > 1e-3
    assert rng.relative_width > 0.99


def test_an_effective_capacity_bound_is_what_actually_narrows_the_prediction(carotenoid_model):
    free = product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=0.90)
    capped = product_flux_range(
        carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=0.90,
        pathway_capacity=free.maximum * 0.1,
    )

    assert capped.maximum == pytest.approx(free.maximum * 0.1, rel=1e-3)
    assert capped.width < free.width


def test_capacity_sweep_reports_one_row_per_requested_capacity(carotenoid_model):
    sweep = capacity_sweep(
        carotenoid_model, PRODUCT_DEMAND_ID,
        capacities=[1e-4, 1e-3, 1e-2], growth_fraction=0.90,
    )

    assert list(sweep.capacity) == [1e-4, 1e-3, 1e-2]
    assert np.all(np.diff(sweep.maximum) > 0)


def test_growth_fraction_outside_the_unit_interval_is_rejected(carotenoid_model):
    with pytest.raises(ValueError, match="growth_fraction"):
        product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, growth_fraction=1.5)


def test_the_model_is_not_mutated_by_running_the_analysis(carotenoid_model):
    before = carotenoid_model.reactions.get_by_id("r_1714").bounds
    product_flux_range(carotenoid_model, PRODUCT_DEMAND_ID, glucose_uptake=15.0, growth_fraction=0.9)

    assert carotenoid_model.reactions.get_by_id("r_1714").bounds == before
