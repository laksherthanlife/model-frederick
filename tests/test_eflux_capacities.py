"""The `capacities` override, which used to silently make every reaction irreversible.

`eflux_layer(..., capacities=...)` is the repair its own module docstring recommends: the
model-implied reference capacity is a very loose ceiling, so a caller with a computed
alternative passes it in. It returned ``(value, 0.0)`` -- forward from the caller, reverse
**zero** -- and zero reverse capacity is not "unspecified", it is the assertion that the
reaction cannot run backwards. On yeast-GEM 9.0.2, handing the layer back its own forward
capacities (the identity, which must be a no-op) locked 92 reversible reactions forward-only
and took max growth to 0.00000 /h at a module activity of 1e-12, i.e. with no regulation at
all. That is the regression pinned below.

The second half of this file pins the other Phase-0 finding about this layer: it is
`LayerState.INERT`, permanently, and the count behind that label is zero.
`tests/test_fba_audit.py::TestRegulationIsWiredAndProvablyInert` already pins that for DTT at
one dose on one model; what is added here is the panel-wide version -- all 25 stressors at 3x
their own EC50, on both shipped GEMs -- and the correction that the reason usually given for
it ("every module activity is non-negative") is false as written.
"""
from __future__ import annotations

import cobra
import pytest

from ystwin.bridge.regulation import (
    EFLUX_BOUNDS_TIGHTENED_AT_3X_EC50,
    EFLUX_LAYER_STATE,
    MODULE_FACTORS,
    Regulon,
    apply_eflux,
    eflux_layer,
    load_regulons,
    override_capacity,
)
from ystwin.generator.stress_panel import (
    METABOLITE_POOLS,
    MODULES,
    STRESSORS,
    module_response,
)

#: What `gene_scales` refuses rather than guesses, taken from the panel's own declaration
#: instead of re-derived here: the five metabolite pools have no transcription factor, and
#: Crz1 has no target gene at MacIsaac's evidence level. Both are facts about the mapping.
UNMAPPED = frozenset(METABOLITE_POOLS | {"calcium"})


def _reaction(model, rid, stoichiometry, gpr="", bounds=(0.0, 1000.0)):
    rxn = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([rxn])
    rxn.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoichiometry.items()})
    rxn.gene_reaction_rule = gpr
    return rxn


@pytest.fixture
def toy():
    """Two routes from B to C: one irreversible, one reversible, both on gene ``g1``.

    ``R_rev`` at (-10, 10) is the shape the bug lives on -- a reversible reaction whose
    reverse capacity the model states and the override used to throw away.
    """
    model = cobra.Model("toy")
    model.add_metabolites([cobra.Metabolite(m) for m in ("A_e", "A", "B", "C")])
    _reaction(model, "EX_A", {"A_e": -1}, bounds=(-10.0, 1000.0))
    _reaction(model, "R_in", {"A_e": -1, "A": 1})
    _reaction(model, "R1", {"A": -1, "B": 1}, gpr="g1")
    _reaction(model, "R_rev", {"B": -1, "C": 1}, gpr="g1", bounds=(-10.0, 10.0))
    biomass = _reaction(model, "BIO", {"C": -1})
    model.objective = biomass
    return model


@pytest.fixture
def regulons():
    return {"m1": Regulon("m1", ("F1",), ("g1",), "test")}


# --------------------------------------------------------------------------- #
# the bug
# --------------------------------------------------------------------------- #


def test_a_supplied_forward_capacity_leaves_the_reverse_direction_where_the_model_put_it(
    toy, regulons
):
    """The regression, at unit scale. ``(value, 0.0)`` is what this used to return.

    The caller said one thing -- how much forward capacity ``R_rev`` has -- and the layer
    silently answered a second question it was not asked, by asserting the reaction is
    irreversible. Overriding a capacity may change a reaction's magnitude; it must not
    change its direction.
    """
    layer = eflux_layer(toy, {"m1": 0.5}, regulons, capacities={"R_rev": 4.0})
    row = next(r for r in layer.scales if r.reaction == "R_rev")

    assert row.reference_capacity == 4.0
    assert row.reverse_capacity == 10.0
    assert (row.lower_bound, row.upper_bound) == pytest.approx((-15.0, 6.0))


def test_the_layers_own_capacities_round_trip_to_exactly_the_same_bounds(toy, regulons):
    """The identity. Handing the layer back what it computed must change nothing.

    This is the property that failed on the real GEM, and it is the cheapest possible check
    of it: a round-trip through the override is a no-op or the override is lying about what
    it overrides.
    """
    implied = eflux_layer(toy, {"m1": -0.5}, regulons)
    for capacities in (
        {r.reaction: r.reference_capacity for r in implied.scales},
        {r.reaction: (r.reference_capacity, r.reverse_capacity) for r in implied.scales},
    ):
        supplied = eflux_layer(toy, {"m1": -0.5}, regulons, capacities=capacities)

        assert [(r.reaction, r.lower_bound, r.upper_bound) for r in supplied.scales] == [
            (r.reaction, r.lower_bound, r.upper_bound) for r in implied.scales]


def test_a_pair_states_both_directions_and_the_caller_owns_both(toy, regulons):
    """FVA under a regime gives a minimum and a maximum, so a caller can state both."""
    layer = eflux_layer(toy, {"m1": -0.5}, regulons, capacities={"R_rev": (8.0, 2.0)})
    row = next(r for r in layer.scales if r.reaction == "R_rev")

    assert (row.lower_bound, row.upper_bound) == pytest.approx((-1.0, 4.0))


def test_an_irreversible_reaction_is_unchanged_by_the_repair(toy, regulons):
    """``R1`` is (0, 1000), so its reverse capacity was, and stays, genuinely zero.

    The old behaviour was correct on exactly this case, which is why the bug survived: every
    reaction in an enzyme-constrained GECKO export is irreversible by construction.
    """
    layer = eflux_layer(toy, {"m1": -0.5}, regulons, capacities={"R1": 2.0})
    row = next(r for r in layer.scales if r.reaction == "R1")

    assert (row.lower_bound, row.upper_bound) == pytest.approx((0.0, 1.0))
    assert layer.reference == "caller-supplied capacities"


def test_a_negative_capacity_is_refused_because_it_is_not_a_quantity(toy):
    with pytest.raises(ValueError, match="non-negative magnitude"):
        override_capacity(toy.reactions.get_by_id("R_rev"), -1.0)
    with pytest.raises(ValueError, match="non-negative magnitude"):
        override_capacity(toy.reactions.get_by_id("R_rev"), (1.0, -1.0))


def test_an_infinite_reverse_bound_is_refused_rather_than_read_as_zero(toy):
    """No finite reference is not the same as no reverse capacity, and this is the seam
    where the original bug's reasoning would reappear."""
    toy.reactions.get_by_id("R_rev").lower_bound = -float("inf")

    with pytest.raises(ValueError, match="infinite reverse bound"):
        override_capacity(toy.reactions.get_by_id("R_rev"), 4.0)

    assert override_capacity(
        toy.reactions.get_by_id("R_rev"), (4.0, 10.0)) == (4.0, 10.0)


# --------------------------------------------------------------------------- #
# the layer is INERT, permanently
# --------------------------------------------------------------------------- #


def test_the_layer_carries_its_own_inert_label_and_prints_it(toy, regulons):
    """A report should read the state off the layer, not paraphrase a docstring."""
    layer = eflux_layer(toy, {"m1": 0.5}, regulons)

    assert EFLUX_LAYER_STATE == "inert"
    assert layer.layer_state == "inert"
    assert "[E-Flux, INERT]" in layer.summary()


def test_no_panel_activity_can_tighten_a_bound_because_the_scale_never_falls_below_one(
    toy, regulons
):
    """The structural reason, at the level of one reaction, so it is not just a count.

    On the range the test below measures ``module_response`` to occupy for a module that has
    a regulon, ``scale = max(0, 1 + a) >= 1``, and an upper bound is only ever raised.
    """
    for activity in (0.0, 1e-12, 0.5, 3.0):
        layer = eflux_layer(toy, {"m1": activity}, regulons)

        assert layer.binding == ()
        assert all(row.scale >= 1.0 for row in layer.scales)


def test_only_the_modules_with_no_regulon_ever_go_negative():
    """The premise above, measured -- and the architecture document has it wrong.

    "Every module activity is non-negative" is false: `atp` reaches -0.381, `ph` -0.448,
    `nadh` -0.150 and `redox` -0.109 across the panel at 3x EC50. What is true, and is what
    the inertness actually rests on, is that those four are metabolite POOLS with no
    transcription factor -- so they have no regulon, E-Flux cannot see them, and every module
    that does carry one is non-negative at every dose from 0 to 3x its EC50.
    """
    doses = [f / 10.0 for f in range(31)]
    with_regulon = [value
                    for name, spec in STRESSORS.items() for f in doses
                    for module, value in module_response(name, f * spec.ec50).items()
                    if module in MODULE_FACTORS]
    without = {module: value
               for name, spec in STRESSORS.items()
               for module, value in module_response(name, 3.0 * spec.ec50).items()
               if module not in MODULE_FACTORS and value < 0.0}

    assert min(with_regulon) == 0.0
    assert sorted(without) == ["atp", "nadh", "ph", "redox"]
    assert min(without.values()) == pytest.approx(-0.4478, abs=1e-4)


# --------------------------------------------------------------------------- #
# against the real model
# --------------------------------------------------------------------------- #


@pytest.mark.integration
def test_the_identity_override_no_longer_kills_growth_on_the_real_gem(yeast_gem):
    """The measured regression: 0.08584 /h -> 0.00000 /h with no regulation applied.

    Activity 1e-12 is as close to unregulated as a float gets, so every bound this layer
    writes is the reaction's own capacity back again and growth must not move. Before the
    repair, 92 reversible reactions lost their reverse capacity here and the LP died.
    """
    regulons = load_regulons()
    activity = {module: 1e-12 for module in MODULES}
    implied = eflux_layer(yeast_gem, activity, regulons, allow_unmapped=UNMAPPED)
    capacities = {row.reaction: row.reference_capacity for row in implied.scales}
    supplied = eflux_layer(yeast_gem, activity, regulons, allow_unmapped=UNMAPPED,
                           capacities=capacities)

    at_risk = [row.reaction for row in implied.scales if row.reverse_capacity > 0.0]
    assert len(at_risk) == 92
    assert [(r.reaction, r.lower_bound, r.upper_bound) for r in supplied.scales] == [
        (r.reaction, r.lower_bound, r.upper_bound) for r in implied.scales]

    with yeast_gem:
        implied_growth = apply_eflux(yeast_gem, implied,
                                     acknowledge_infeasible=True).max_growth
    with yeast_gem:
        supplied_growth = apply_eflux(yeast_gem, supplied,
                                      acknowledge_infeasible=True).max_growth

    assert implied_growth == pytest.approx(0.08584, abs=1e-5)
    assert supplied_growth == pytest.approx(implied_growth, rel=1e-12)


@pytest.mark.integration
def test_gecko_never_had_the_bug_because_it_has_no_reversible_reaction_to_lose(ec_yeast_gem):
    """Why this survived: the enzyme-constrained model, which most callers use, is immune.

    GECKO splits every reversible reaction into two irreversible ones, so no row here
    carries a reverse capacity and ``(value, 0.0)`` was accidentally right.
    """
    layer = eflux_layer(ec_yeast_gem, {module: 1e-12 for module in MODULES},
                        load_regulons(), allow_unmapped=UNMAPPED)

    assert layer.scales
    assert all(row.reverse_capacity == 0.0 for row in layer.scales)


@pytest.mark.integration
def test_zero_bounds_are_tightened_across_all_twenty_five_stressors_at_three_times_ec50(
    yeast_gem, ec_yeast_gem
):
    """The measurement behind the INERT label, on both shipped GEMs.

    This is the whole case for retiring E-Flux as a live regulation route: driven as hard as
    the panel goes, on the models this project actually solves, it changes nothing. Kept as
    the null arm, not offered as a mechanism.
    """
    regulons = load_regulons()
    for model, fewest in ((yeast_gem, 54), (ec_yeast_gem, 16)):
        layers = [eflux_layer(model, module_response(name, 3.0 * spec.ec50), regulons,
                              allow_unmapped=UNMAPPED)
                  for name, spec in STRESSORS.items()]

        assert len(STRESSORS) == 25
        assert sum(len(layer.binding) for layer in layers) == 0
        assert EFLUX_BOUNDS_TIGHTENED_AT_3X_EC50 == 0
        # Zero out of scales that computed, not zero because nothing ran. The second
        # would be a wiring fault; only the first is a property of the panel.
        assert min(len(layer.scales) for layer in layers) == fewest
