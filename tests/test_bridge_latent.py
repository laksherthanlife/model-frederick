"""Branch B: latent stress state to FBA constraints. Research only.

This is the half of g_psi the plan describes and the data do not yet support. G4
returned INCONCLUSIVE for all four constructs, so no module is admissible, and §9's
rule is that a latent-to-constraint mapping is kept only once it predicts an excluded
anchor.

It is built anyway, because seeing what it *would* add is worth having. But it is
built so that its output cannot be mistaken for Branch A's: the branch label says
"latent", the validation status travels with every prediction, and calling it without
acknowledging that status raises.
"""

import pytest

from ystwin.bridge.latent_bridge import (
    LatentState,
    latent_constraints,
    compare_branches,
)
from ystwin.bridge.physiology_bridge import MeasuredState, physiology_constraints

pytestmark = pytest.mark.integration


def _measured(**over):
    base = dict(biomass_gl=0.35, growth_rate=0.28, glucose_mM=20.0, time_h=2.0)
    base.update(over)
    return MeasuredState(**base)


def _latent(activity=135.2, construct="UPRE1"):
    """Activities are RFU/OD/h on the measured axis, where full scale is
    ``_ACTIVITY_FULL_SCALE`` (1465.0, the largest excess on the four characterised
    constructs).

    They were 1.0e-3 to 6.0e-3 until 2026-09-04, and those were NOT small readings -- they
    were left over from baseline commit f562d8f, when the coefficient was
    ``_MAINTENANCE_PER_ACTIVITY = 300.0`` and activity was MULTIPLIED by it: 6.0e-3 x 300
    gave NGAM 2.5 mmol ATP/gDW/h, a real stress load. The unit fix inverted the operation to
    a division by full scale and the inputs here were never rescaled, so from then on they
    meant 4e-6 of full scale -- an NGAM increment of 2.7e-5 over resting. The tests below
    kept passing on solver noise, and `test_the_latent_branch_changes_the_flux_it_claims_to
    _change` only revealed it when a 12% change in the denominator tipped its 1e-6 relative
    tolerance. Each value here is the old one carried onto the new axis at equal NGAM."""
    return LatentState(promoter_activity=activity, construct=construct, time_h=2.0)


class TestItCannotBeMistakenForValidatedWork:
    def test_calling_it_without_acknowledging_the_status_raises(self):
        with pytest.raises(ValueError, match="not validated"):
            latent_constraints(_latent(), physiology_constraints(_measured()))

    def test_acknowledging_it_returns_constraints_labelled_latent(self):
        constraints = latent_constraints(
            _latent(), physiology_constraints(_measured()), acknowledge_unvalidated=True
        )

        assert constraints.branch == "latent"

    def test_the_validation_status_travels_with_the_constraints(self):
        constraints = latent_constraints(
            _latent(), physiology_constraints(_measured()), acknowledge_unvalidated=True
        )

        assert constraints.validation.startswith("REFUSED")
        assert "G4" in constraints.validation

    def test_the_two_branches_are_never_labelled_the_same(self):
        physical = physiology_constraints(_measured())
        latent = latent_constraints(_latent(), physical, acknowledge_unvalidated=True)

        assert physical.branch != latent.branch


class TestTheMappingItself:
    def test_more_stress_raises_the_maintenance_demand(self):
        low = latent_constraints(
            _latent(activity=67.6), physiology_constraints(_measured()),
            acknowledge_unvalidated=True,
        )
        high = latent_constraints(
            _latent(activity=270.5), physiology_constraints(_measured()),
            acknowledge_unvalidated=True,
        )

        assert high.atp_maintenance > low.atp_maintenance

    def test_basal_activity_leaves_maintenance_at_its_resting_value(self):
        constraints = latent_constraints(
            _latent(activity=67.6), physiology_constraints(_measured()),
            acknowledge_unvalidated=True, basal_activity=67.6,
        )

        assert constraints.atp_maintenance == pytest.approx(constraints.resting_maintenance)

    def test_it_inherits_the_physical_bounds_it_was_given(self):
        physical = physiology_constraints(_measured())
        latent = latent_constraints(_latent(), physical, acknowledge_unvalidated=True)

        assert latent.glucose_uptake == pytest.approx(physical.glucose_uptake)
        assert latent.growth_lower_bound == pytest.approx(physical.growth_lower_bound)

    def test_the_effect_is_named_not_a_bare_multiplier(self):
        """§9 forbids post-solve product multipliers; the effect must be a constraint."""
        constraints = latent_constraints(
            _latent(), physiology_constraints(_measured()), acknowledge_unvalidated=True
        )

        assert constraints.named_effect == "ATP maintenance"


def test_comparing_branches_reports_both_and_says_which_is_load_bearing(yeast_gem):
    physical = physiology_constraints(_measured())
    latent = latent_constraints(_latent(activity=270.5), physical, acknowledge_unvalidated=True)

    table = compare_branches(yeast_gem, physical, latent, reactions=["r_2111", "r_1761"])

    assert set(table.branch) == {"physiology", "latent"}
    assert table.set_index("branch").loc["physiology", "load_bearing"]
    assert not table.set_index("branch").loc["latent", "load_bearing"]


def test_the_latent_branch_changes_the_flux_it_claims_to_change(yeast_gem):
    """If the constraint made no difference, the module would be doing nothing."""
    physical = physiology_constraints(_measured())
    latent = latent_constraints(_latent(activity=405.7), physical, acknowledge_unvalidated=True)

    table = compare_branches(yeast_gem, physical, latent, reactions=["r_1992"])
    oxygen = table.set_index("branch")["flux_r_1992"]

    assert oxygen["latent"] != pytest.approx(oxygen["physiology"], rel=1e-6)


def test_a_higher_maintenance_demand_costs_more_substrate_at_the_same_growth(yeast_gem):
    """The direction the mapping predicts: stress is paid for in carbon and oxygen."""
    physical = physiology_constraints(_measured())
    latent = latent_constraints(_latent(activity=405.7), physical, acknowledge_unvalidated=True)

    table = compare_branches(
        yeast_gem, physical, latent, reactions=["r_2111", "r_1714", "r_1992"]
    ).set_index("branch")

    assert table.loc["latent", "flux_r_2111"] == pytest.approx(table.loc["physiology", "flux_r_2111"])
    assert abs(table.loc["latent", "flux_r_1714"]) > abs(table.loc["physiology", "flux_r_1714"])
    assert abs(table.loc["latent", "flux_r_1992"]) > abs(table.loc["physiology", "flux_r_1992"])


def test_the_known_respiratory_limitation_shows_up_here_too(yeast_gem):
    """With oxygen unconstrained the model fully respires and makes no ethanol.

    Documented rather than worked around: it is the same limitation D1 found, and it
    is why no oxygen- or overflow-linked claim rests on this solve.
    """
    physical = physiology_constraints(_measured())
    latent = latent_constraints(_latent(), physical, acknowledge_unvalidated=True)

    table = compare_branches(yeast_gem, physical, latent, reactions=["r_1761"])

    assert (table["flux_r_1761"].abs() < 1e-6).all()


def test_the_activities_these_tests_use_are_a_meaningful_fraction_of_full_scale():
    """The guard for the failure mode found on 2026-09-04.

    Every assertion in this module is about what the latent branch DOES to a flux, so it is
    only meaningful if the activity fed in represents real stress. For eight months it did
    not: the inputs were 4e-6 of full scale, the NGAM increment was 2.7e-5 mmol ATP/gDW/h,
    and the tests passed on floating-point noise rather than on the effect they name.

    A unit change that inverts an operation -- here multiply-by-300 becoming
    divide-by-full-scale -- rescales every input by ~10^5 and breaks no test, which is
    exactly why nothing caught it. This one would.
    """
    from ystwin.bridge.latent_bridge import _ACTIVITY_FULL_SCALE

    used = [67.6, 135.2, 270.5, 405.7]

    for activity in used:
        fraction = activity / _ACTIVITY_FULL_SCALE
        assert 0.01 < fraction < 1.0, (
            f"activity {activity} is {fraction:.2e} of full scale "
            f"{_ACTIVITY_FULL_SCALE}; at that level the constraint moves nothing and any "
            f"test asserting it does is reading solver noise")


def test_the_flux_change_is_far_larger_than_the_tolerance_that_checks_it():
    """The specific weakness that let the above go unnoticed: an assertion of
    `!= approx(rel=1e-6)` passes on noise. Pin real headroom instead."""
    from ystwin.bridge.latent_bridge import (
        _ACTIVITY_FULL_SCALE,
        _MAX_STRESS_MAINTENANCE,
        _RESTING_MAINTENANCE,
    )

    increment = (405.7 / _ACTIVITY_FULL_SCALE) * _MAX_STRESS_MAINTENANCE

    assert increment > 0.1 * _RESTING_MAINTENANCE
