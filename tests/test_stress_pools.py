"""Retained osmolyte pools: a demand carrying mu*[X], priced against exporting the same flux.

Unless a test says otherwise every number here was measured on Yeast9 v9.0.2 at glucose <= 10
mmol/gDCW/h with oxygen free, where maximum growth is 0.887685 /h. The export ladder this file
pins (-3.10% / -6.20% / -31.02% at 0.5 / 1.0 / 5.0 mmol/gDCW/h out of r_1808) is the
four-auditor GEM audit's, reproduced here so the retained-pool costs sit beside a number that
predates them. Two tests repeat the comparison at the repo's REFERENCE_AEROBIC_BATCH.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.integration

GLUCOSE_UPTAKE = 10.0
UNSTRESSED_GROWTH = 0.887685
CYTOSOLIC_VOLUME_ML_PER_GDCW = 2.0


@pytest.fixture(scope="module")
def _pooled_template(yeast_gem_factory):
    from ystwin.fba.stress_pools import add_retained_pool_sinks

    return add_retained_pool_sinks(yeast_gem_factory())


@pytest.fixture
def pooled(_pooled_template, model_copy):
    return model_copy(_pooled_template)


def _growth(model, glucose_uptake=GLUCOSE_UPTAKE):
    with model as m:
        m.reactions.get_by_id("r_1714").lower_bound = -abs(glucose_uptake)
        return m.slim_optimize()


# --------------------------------------------------------------------------- #
# the gap this module fills
# --------------------------------------------------------------------------- #


def test_the_model_cannot_hold_glycerol_before_this_module_touches_it(yeast_gem):
    """The gap: every glycerol reaction moves it somewhere, none lets it sit."""
    glycerol = yeast_gem.metabolites.get_by_id("s_0765")

    assert not any(rxn.boundary for rxn in glycerol.reactions)
    assert yeast_gem.reactions.get_by_id("r_1172").reaction == "s_0765 --> s_0766"


def test_the_paired_measurement_is_retention_and_not_export():
    """One set of samples: inside rises 18.2x, the supernatant 1.47x, the control falls."""
    from ystwin.fba.stress_pools import FPS1_CLOSURE_MOL_PER_L, fps1_closure_folds

    folds = fps1_closure_folds()

    assert FPS1_CLOSURE_MOL_PER_L["internal_peak"] == pytest.approx(0.99685, rel=1e-4)
    assert folds["internal"] == pytest.approx(18.18, rel=1e-3)
    assert folds["external"] == pytest.approx(1.470, rel=1e-3)
    assert folds["unstressed"] == pytest.approx(0.7042, rel=1e-3)
    assert folds["internal"] > 12 * folds["external"]


# --------------------------------------------------------------------------- #
# what installing the pools does to the model
# --------------------------------------------------------------------------- #


def test_each_pool_is_a_boundary_demand_on_the_cytosolic_osmolyte(pooled):
    from ystwin.fba.stress_pools import GLYCEROL_C, TREHALOSE_C, retained_demand_id

    for base_id, formula in ((GLYCEROL_C, "C3H8O3"), (TREHALOSE_C, "C12H22O11")):
        demand = pooled.reactions.get_by_id(retained_demand_id(base_id))
        metabolite, coefficient = next(iter(demand.metabolites.items()))

        assert demand.boundary and demand.bounds == (0, 1000.0)
        assert metabolite.id == base_id and coefficient == -1
        assert metabolite.formula == formula and metabolite.charge == 0
        assert metabolite.compartment == "c"


def test_a_demand_is_unbalanced_by_construction_and_says_so(pooled):
    """A drain is not a reaction. Its imbalance is exactly the neutral species it drains."""
    from ystwin.fba.stress_pools import GLYCEROL_C, TREHALOSE_C, retained_demand_id

    glycerol = pooled.reactions.get_by_id(retained_demand_id(GLYCEROL_C))
    trehalose = pooled.reactions.get_by_id(retained_demand_id(TREHALOSE_C))

    assert glycerol.check_mass_balance() == {"C": -3, "H": -8, "O": -3}
    assert trehalose.check_mass_balance() == {"C": -12, "H": -22, "O": -11}
    assert "charge" not in glycerol.check_mass_balance()
    assert "charge" not in trehalose.check_mass_balance()


def test_installing_pools_adds_no_unbalanced_internal_reaction(yeast_gem, pooled):
    added = set(r.id for r in pooled.reactions) - set(r.id for r in yeast_gem.reactions)

    assert added == {"DM_s_0765_retained", "DM_s_1520_retained"}
    assert all(pooled.reactions.get_by_id(rid).boundary for rid in added)


def test_the_source_model_is_left_untouched(yeast_gem, pooled):
    assert not yeast_gem.reactions.has_id("DM_s_0765_retained")
    assert pooled.reactions.has_id("DM_s_0765_retained")


def test_installing_pools_does_not_move_growth(yeast_gem, pooled):
    """An open drain the objective does not want sits at zero, so no existing answer moves."""
    assert _growth(pooled) == pytest.approx(_growth(yeast_gem), rel=1e-9)
    assert _growth(yeast_gem) == pytest.approx(UNSTRESSED_GROWTH, rel=1e-5)
    assert pooled.slim_optimize() == pytest.approx(yeast_gem.slim_optimize(), rel=1e-6)


def test_both_pools_can_carry_flux(pooled):
    from cobra.flux_analysis import flux_variability_analysis

    from ystwin.fba.stress_pools import GLYCEROL_C, TREHALOSE_C, retained_demand_id

    with pooled as m:
        m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE_UPTAKE
        fva = flux_variability_analysis(
            m, [retained_demand_id(GLYCEROL_C), retained_demand_id(TREHALOSE_C)],
            fraction_of_optimum=0.0, processes=1)

    assert fva.loc[retained_demand_id(GLYCEROL_C), "maximum"] == pytest.approx(16.1204, rel=1e-3)
    assert fva.loc[retained_demand_id(TREHALOSE_C), "maximum"] == pytest.approx(4.6581, rel=1e-3)
    assert (fva["minimum"] == 0.0).all()


def test_installing_the_same_pool_twice_is_refused(pooled):
    from ystwin.fba.stress_pools import add_retained_pool_sinks

    with pytest.raises(ValueError, match="already installed"):
        add_retained_pool_sinks(pooled)


# --------------------------------------------------------------------------- #
# the measured pools, and the units they were measured in
# --------------------------------------------------------------------------- #


def test_a_molar_pool_refuses_to_become_a_content_without_a_volume():
    from ystwin.fba.stress_pools import MEASURED_POOLS

    with pytest.raises(ValueError, match="cytosolic volume"):
        MEASURED_POOLS["glycerol_osmotic"].content()


def test_a_per_dry_weight_pool_refuses_a_volume_it_never_used():
    from ystwin.fba.stress_pools import MEASURED_POOLS

    with pytest.raises(ValueError, match="needs no cytosolic volume"):
        MEASURED_POOLS["trehalose_chemostat"].content(2.0)


def test_every_pool_carries_its_citation_and_its_condition():
    from ystwin.fba.stress_pools import MEASURED_POOLS

    for name, pool in MEASURED_POOLS.items():
        assert pool.source.strip(), f"{name} has no citation"
        assert pool.condition.strip(), f"{name} has no measured condition"
        assert (pool.concentration_m is None) != (pool.content_mmol_per_gdcw is None)


def test_the_measured_osmotic_pool_is_the_one_petelenz_kurdziel_reports():
    """0.9969 M internal glycerol, 90 min after 0.4 M NaCl, against 0.0548 M before it."""
    from ystwin.fba.stress_pools import MEASURED_POOLS

    stressed = MEASURED_POOLS["glycerol_osmotic"]
    basal = MEASURED_POOLS["glycerol_unstressed"]

    assert stressed.concentration_m == pytest.approx(0.99685, rel=1e-4)
    assert stressed.concentration_m / basal.concentration_m == pytest.approx(18.18, rel=1e-3)
    assert stressed.content(2.0) == pytest.approx(1.99371, rel=1e-4)
    assert stressed.content_mg_per_gdcw(2.7) == pytest.approx(247.9, rel=1e-3)


def test_the_two_independent_glycerol_measurements_bracket_each_other():
    """They agree at 2.62 mL/gDCW, inside thermo_gate's 1.0-2.7 range. Not a calibration."""
    from ystwin.fba.stress_pools import MEASURED_POOLS

    per_dry_weight = MEASURED_POOLS["glycerol_osmotic_dry_weight"].content()
    molar = MEASURED_POOLS["glycerol_osmotic"]

    assert per_dry_weight / molar.concentration_m == pytest.approx(2.624, rel=1e-3)
    assert molar.content(1.0) < per_dry_weight < molar.content(2.7)


def test_retention_flux_is_the_diluted_balance_pathway_solve_uses():
    """v = mu*[X] is the inverse of solve.py's terminal DILUTED line, X = v_in/mu."""
    from ystwin.fba.stress_pools import retention_flux

    assert retention_flux(1.99371, 0.8) == pytest.approx(1.594968)
    assert retention_flux(1.99371, 0.8) / 0.8 == pytest.approx(1.99371)


# --------------------------------------------------------------------------- #
# what holding a pool costs, beside what exporting it costs
# --------------------------------------------------------------------------- #


def test_holding_the_measured_osmotic_pool_costs_ten_percent_of_growth(pooled):
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, retention_versus_export

    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)
    cost = retention_versus_export(pooled, GLYCEROL_C, content, GLUCOSE_UPTAKE)

    assert cost.growth_free == pytest.approx(0.887685, rel=1e-5)
    assert cost.growth_held == pytest.approx(0.799871, rel=1e-4)
    assert cost.cost_of_holding == pytest.approx(0.098925, rel=1e-3)
    assert cost.retention_flux == pytest.approx(1.5947, rel=1e-3)


def test_holding_a_flux_and_exporting_it_cost_the_same(pooled):
    """The GEM prices the SYNTHESIS. r_1172 is a free channel, so the destination is free."""
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, retention_versus_export

    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)
    cost = retention_versus_export(pooled, GLYCEROL_C, content, GLUCOSE_UPTAKE)

    assert cost.growth_held == pytest.approx(cost.growth_exported, abs=1e-9)


def test_the_export_ladder_the_audit_measured_reproduces(pooled):
    from ystwin.fba.stress_pools import cost_of_exporting

    free = _growth(pooled)
    costs = {flux: 1 - cost_of_exporting(pooled, flux, GLUCOSE_UPTAKE) / free
             for flux in (0.5, 1.0, 5.0)}

    assert costs[0.5] == pytest.approx(0.0310, abs=1e-4)
    assert costs[1.0] == pytest.approx(0.0620, abs=1e-4)
    assert costs[5.0] == pytest.approx(0.3102, abs=1e-4)


def test_an_open_fps1_costs_more_than_a_shut_one_at_the_same_pool(pooled):
    """Closing the channel is worth 5.6 points of growth at the measured osmotic pool."""
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, retention_versus_export

    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)
    leaking = retention_versus_export(
        pooled, GLYCEROL_C, content, GLUCOSE_UPTAKE, leak_flux_mmol_per_gdcw_h=1.0)

    assert leaking.cost_of_leaking == pytest.approx(0.1548, abs=1e-3)
    assert leaking.cost_of_leaking - leaking.cost_of_holding == pytest.approx(0.0559, abs=1e-3)


def test_a_retained_pool_is_a_yield_tax_and_a_forced_export_is_not(pooled):
    """The point of the pool. mu*[X] costs the same FRACTION at every uptake; 1.0 out does not."""
    from ystwin.fba.stress_pools import (
        GLYCEROL_C,
        MEASURED_POOLS,
        cost_of_exporting,
        retention_versus_export,
    )

    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)
    rich = retention_versus_export(pooled, GLYCEROL_C, content, 10.0)
    lean = retention_versus_export(pooled, GLYCEROL_C, content, 2.0)
    export_rich = 1 - cost_of_exporting(pooled, 1.0, 10.0) / rich.growth_free
    export_lean = 1 - cost_of_exporting(pooled, 1.0, 2.0) / lean.growth_free

    assert lean.cost_of_holding == pytest.approx(rich.cost_of_holding, rel=1e-6)
    assert export_rich == pytest.approx(0.0620, abs=1e-4)
    assert export_lean == pytest.approx(0.3148, abs=1e-4)


def test_the_same_comparison_in_the_repos_own_reference_regime(pooled):
    """Glucose 11.1 / oxygen 3.7, fba/physiology.py: the pool costs 3.79%, a 1.0 leak 5.71%."""
    from ystwin.fba.physiology import REFERENCE_AEROBIC_BATCH
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, retention_versus_export

    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)
    cost = retention_versus_export(
        pooled, GLYCEROL_C, content,
        REFERENCE_AEROBIC_BATCH.glucose_uptake, REFERENCE_AEROBIC_BATCH.oxygen_uptake,
        leak_flux_mmol_per_gdcw_h=1.0)

    assert cost.growth_free == pytest.approx(0.345984, rel=1e-5)
    assert cost.cost_of_holding == pytest.approx(0.037901, abs=1e-4)
    assert cost.cost_of_holding == pytest.approx(cost.cost_of_exporting, abs=1e-9)
    assert cost.cost_of_leaking - cost.cost_of_holding == pytest.approx(0.0549, abs=1e-3)


def test_the_volume_free_glycerol_pool_costs_thirteen_percent(pooled):
    """A salt-stressed pool measured per DRY WEIGHT, so no cytosolic volume enters the answer."""
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, retention_versus_export

    pool = MEASURED_POOLS["glycerol_osmotic_dry_weight"]
    cost = retention_versus_export(pooled, GLYCEROL_C, pool.content(), GLUCOSE_UPTAKE)

    assert pool.content_mg_per_gdcw() == pytest.approx(240.87, rel=1e-9)
    assert cost.growth_held == pytest.approx(0.775928, rel=1e-4)
    assert cost.cost_of_holding == pytest.approx(0.12590, abs=1e-4)
    assert cost.retention_flux == pytest.approx(2.0295, rel=1e-3)


def test_the_osmotic_trehalose_pool_is_smaller_than_the_biomass_already_pins(pooled):
    """r_4048 pins trehalose at 0.13655 mmol/gDCW, 5x the measured osmotic pool."""
    from ystwin.fba.stress_pools import MEASURED_POOLS, TREHALOSE_C, retention_versus_export

    content = MEASURED_POOLS["trehalose_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)
    cost = retention_versus_export(
        pooled, TREHALOSE_C, content, GLUCOSE_UPTAKE, exchange_reaction="r_1650")

    assert content == pytest.approx(0.027516, rel=1e-4)
    assert content < 0.13655
    assert cost.cost_of_holding == pytest.approx(0.00522, abs=1e-4)


def test_the_chemostat_trehalose_pool_needs_no_volume_at_all(pooled):
    from ystwin.fba.stress_pools import MEASURED_POOLS, TREHALOSE_C, retention_versus_export

    pool = MEASURED_POOLS["trehalose_chemostat"]
    cost = retention_versus_export(
        pooled, TREHALOSE_C, pool.content(), GLUCOSE_UPTAKE, exchange_reaction="r_1650")

    assert pool.content_mg_per_gdcw() == pytest.approx(19.4, rel=1e-6)
    assert cost.cost_of_holding == pytest.approx(0.01069, abs=1e-4)


# --------------------------------------------------------------------------- #
# what it changes downstream, and what it refuses
# --------------------------------------------------------------------------- #


def test_holding_the_pool_lowers_the_product_ceiling_the_audit_checks_against(yeast_gem):
    """What this changes: fba/audit.py's feasible_max falls 14.2% at mu = 0.5 /h."""
    from ystwin.fba.carotenoid import add_beta_carotene_pathway
    from ystwin.fba.stress_pools import (
        GLYCEROL_C,
        MEASURED_POOLS,
        add_retained_pool_sinks,
        hold_pool,
    )

    model = add_retained_pool_sinks(add_beta_carotene_pathway(yeast_gem))
    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)

    def ceiling(hold):
        with model as m:
            m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE_UPTAKE
            m.reactions.get_by_id("r_2111").bounds = (0.5, 0.5)
            if hold:
                hold_pool(m, GLYCEROL_C, content)
            m.objective = "DM_betacarotene_c"
            return m.slim_optimize()

    free, held = ceiling(False), ceiling(True)

    assert free == pytest.approx(0.427996, rel=1e-4)
    assert held == pytest.approx(0.367396, rel=1e-4)
    assert 1 - held / free == pytest.approx(0.1416, abs=1e-3)


def test_the_demand_adds_nothing_until_fps1_shuts_and_then_it_is_the_only_route(pooled, yeast_gem):
    """Criterion (a). With Fps1 open the sink is redundant; shut, it is the whole module."""
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, hold_pool

    content = MEASURED_POOLS["glycerol_osmotic"].content(CYTOSOLIC_VOLUME_ML_PER_GDCW)

    def held_via_exchange(model, fps1_open):
        with model as m:
            m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE_UPTAKE
            if not fps1_open:
                m.reactions.get_by_id("r_1172").bounds = (0.0, 0.0)
            export, biomass = m.reactions.get_by_id("r_1808"), m.reactions.get_by_id("r_2111")
            m.add_cons_vars([m.problem.Constraint(
                export.flux_expression - content * biomass.flux_expression,
                lb=0.0, ub=0.0, name="retain_via_export")])
            return m.slim_optimize()

    def held_via_demand(fps1_open):
        with pooled as m:
            m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE_UPTAKE
            if not fps1_open:
                m.reactions.get_by_id("r_1172").bounds = (0.0, 0.0)
            hold_pool(m, GLYCEROL_C, content)
            return m.slim_optimize()

    # Fps1 open: the pre-existing exchange reproduces the sink exactly, on a model with no
    # sink installed at all, so `add_retained_pool_sinks` moved no answer.
    assert held_via_exchange(yeast_gem, True) == pytest.approx(0.799871, rel=1e-4)
    assert held_via_demand(True) == pytest.approx(held_via_exchange(yeast_gem, True), abs=1e-9)

    # Fps1 shut: glycerol cannot reach s_0766, so the export route pins growth at zero.
    assert held_via_exchange(pooled, False) == pytest.approx(0.0, abs=1e-9)
    assert held_via_demand(False) == pytest.approx(0.799871, rel=1e-4)


def test_the_cost_is_a_band_across_the_volumes_nobody_has_measured(pooled):
    """9.89% is one volume's answer. thermo_gate's range makes it 5.20-12.91%, a 2.48x band."""
    from ystwin.pathway.thermo_gate import CYTOSOLIC_VOLUMES_ML_PER_GDCW
    from ystwin.fba.stress_pools import GLYCEROL_C, MEASURED_POOLS, retention_versus_export

    pool = MEASURED_POOLS["glycerol_osmotic"]
    costs = [retention_versus_export(pooled, GLYCEROL_C, pool.content(v), GLUCOSE_UPTAKE)
             .cost_of_holding for v in CYTOSOLIC_VOLUMES_ML_PER_GDCW]

    assert CYTOSOLIC_VOLUMES_ML_PER_GDCW == (1.0, 2.0, 2.7)
    assert costs == pytest.approx([0.052036, 0.098925, 0.129080], abs=1e-4)
    assert costs[2] / costs[0] == pytest.approx(2.48, rel=1e-2)


def test_a_negative_export_would_have_been_an_uptake_and_is_refused(pooled):
    """The heterologous module's bug in another dress: -5.0 out of r_1808 grew 1.1398 /h."""
    from ystwin.fba.stress_pools import GLYCEROL_C, cost_of_exporting, retention_versus_export

    with pytest.raises(ValueError, match="export flux must be non-negative"):
        cost_of_exporting(pooled, -5.0, GLUCOSE_UPTAKE)
    with pytest.raises(ValueError, match="leak flux must be non-negative"):
        retention_versus_export(
            pooled, GLYCEROL_C, 1.0, GLUCOSE_UPTAKE, leak_flux_mmol_per_gdcw_h=-5.0)


def test_a_regime_with_no_growth_is_refused_rather_than_divided_by(pooled):
    from ystwin.fba.stress_pools import GLYCEROL_C, retention_versus_export

    with pytest.raises(ValueError, match="does not grow"):
        retention_versus_export(pooled, GLYCEROL_C, 1.0, glucose_uptake=0.0, oxygen_uptake=0.0)


def test_an_osmolyte_the_model_does_not_carry_is_named_not_guessed(pooled):
    from ystwin.fba.stress_pools import add_retained_pool_sinks

    with pytest.raises(KeyError, match="s_9999"):
        add_retained_pool_sinks(pooled, ("s_9999",))


def test_the_coupling_constraint_leaves_the_model_as_it_found_it(pooled):
    from ystwin.fba.stress_pools import GLYCEROL_C, hold_pool

    before = _growth(pooled)
    with pooled as m:
        m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE_UPTAKE
        hold_pool(m, GLYCEROL_C, 2.0)
        assert m.slim_optimize() < before

    assert _growth(pooled) == pytest.approx(before, rel=1e-9)
