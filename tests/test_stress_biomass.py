"""Biomass composition as a parameter, and what the GEM charges for moving it.

Yeast9 pins the six carbohydrate components of r_4048, so no adaptation can change composition and
every one of them reads as an exact -100% lethal. These tests pin the shipped numbers, check that
each renormalisation conserves the gram it claims to, and record the growth costs the repair makes
measurable. Growth is at glucose 10 mmol/gDCW/h with oxygen free, base 0.887685 /h, unless the test
says REFERENCE_AEROBIC_BATCH -- glucose 11.1, oxygen 3.7, base 0.345984 /h.
"""

import cobra
import pytest

from ystwin.fba import stress_biomass as sb

pytestmark = pytest.mark.integration

GLUCOSE, OXYGEN = "r_1714", "r_1992"
TREHALOSE_SYNTHASE, GLYCOGEN_SYNTHASE = "r_1051", "r_0510"
BASE_GROWTH = 0.8876853585594559
REFERENCE_GROWTH = 0.3459844408121987


def grow(model, glucose=10.0, oxygen=1000.0):
    with model as m:
        m.reactions.get_by_id(GLUCOSE).lower_bound = -abs(glucose)
        m.reactions.get_by_id(OXYGEN).lower_bound = -abs(oxygen)
        return float(m.slim_optimize())


def grow_without(model, reaction_id, glucose=10.0, oxygen=1000.0):
    with model as m:
        m.reactions.get_by_id(GLUCOSE).lower_bound = -abs(glucose)
        m.reactions.get_by_id(OXYGEN).lower_bound = -abs(oxygen)
        m.reactions.get_by_id(reaction_id).knock_out()
        return float(m.slim_optimize())


@pytest.fixture(scope="module")
def _deletant_model_template(yeast_gem_factory):
    """The measured tps1 gsy1 gsy2 composition: no trehalose, no glycogen."""
    return sb.apply_measured(yeast_gem_factory(), "sillje1999_scu10", renormalise="whole_cell")


@pytest.fixture
def deletant_model(_deletant_model_template, model_copy):
    return model_copy(_deletant_model_template)


@pytest.fixture(scope="module")
def _thick_wall_model_template(yeast_gem_factory):
    return sb.apply_measured(yeast_gem_factory(), "cwi_gas1_ram1998", renormalise="whole_cell")


@pytest.fixture
def thick_wall_model(_thick_wall_model_template, model_copy):
    return model_copy(_thick_wall_model_template)


@pytest.fixture(scope="module")
def _heat_shock_model_template(yeast_gem_factory):
    """Hottiger 1987's 1 g trehalose per g protein, as the whole-cell fixed point."""
    return sb.apply_measured(yeast_gem_factory(), "hottiger1987_heat_shock", renormalise="whole_cell")


@pytest.fixture
def heat_shock_model(_heat_shock_model_template, model_copy):
    return model_copy(_heat_shock_model_template)


@pytest.fixture(scope="module")
def _zeroed_models_template(yeast_gem_factory):
    """Each component driven out of r_4048 in turn -- the fixture the ablation study needed."""
    yeast_gem = yeast_gem_factory()
    return {
        name: sb.set_component_masses(yeast_gem, {name: 0.0}, renormalise="whole_cell")
        for name in sb.CARBOHYDRATE_COMPONENTS
    }


@pytest.fixture
def zeroed_models(_zeroed_models_template, model_copy):
    return {name: model_copy(model) for name, model in _zeroed_models_template.items()}


# --- what Yeast9 actually ships ---------------------------------------------------------------


def test_the_shipped_carbohydrate_coefficients_are_the_audited_constants(yeast_gem):
    coefficients = {
        name: -yeast_gem.reactions.get_by_id("r_4048").get_coefficient(mid)
        for name, mid in sb.CARBOHYDRATE_COMPONENTS.items()
    }

    assert coefficients == pytest.approx(
        {
            "beta_1_3_glucan": 0.73914,
            "beta_1_6_glucan": 0.24696,
            "chitin": 0.02361,
            "glycogen": 0.35689,
            "mannan": 0.70204,
            "trehalose": 0.13655,
        },
        rel=1e-9,
    )


def test_the_carbohydrate_pool_weighs_383_mg_of_which_the_wall_is_278(yeast_gem):
    masses = sb.component_masses(yeast_gem)

    assert sum(masses.values()) == pytest.approx(383.1204, rel=1e-6)
    assert sum(masses[k] for k in sb.WALL_COMPONENTS) == pytest.approx(278.5134, rel=1e-6)
    assert masses["trehalose"] == pytest.approx(46.7406, rel=1e-6)


def test_the_biomass_reaction_sums_to_957_mg_not_the_1000_it_is_named_for(yeast_gem):
    """Worth pinning: renormalising to an assumed 1 gDCW would silently move every coefficient."""
    assert sb.total_biomass_mass(yeast_gem) == pytest.approx(956.9167, rel=1e-6)
    assert sb.pool_masses(yeast_gem)["protein"] == pytest.approx(464.0244, rel=1e-6)


def test_the_shipped_pseudoreaction_is_not_element_balanced_because_the_pool_carries_no_formula(
    yeast_gem,
):
    imbalance = yeast_gem.reactions.get_by_id("r_4048").check_mass_balance()

    assert set(imbalance) == {"C", "H", "O", "N"}
    assert "charge" not in imbalance
    assert sb.reaction_mass(yeast_gem, "r_4048") == pytest.approx(383.1204, rel=1e-6)


def test_the_shipped_protein_and_rna_pools_are_already_charge_imbalanced(yeast_gem):
    """Not something this module introduces -- the check below is that it does not make it worse."""
    assert yeast_gem.reactions.get_by_id("r_4047").check_mass_balance()["charge"] == pytest.approx(
        -4.05850098, rel=1e-6)
    assert yeast_gem.reactions.get_by_id("r_4049").check_mass_balance()["charge"] == pytest.approx(
        0.38067600, rel=1e-6)
    assert yeast_gem.reactions.get_by_id("r_4041").check_mass_balance() == {}


# --- the measured states carry their provenance -------------------------------------------------


def test_every_measured_state_names_a_source_a_regime_and_a_derivation():
    for name, record in sb.MEASURED_STATES.items():
        assert record.source and record.regime and record.derivation, name
        given = [record.masses, record.factors, record.protein_ratios]
        assert sum(field is not None for field in given) == 1, name


def test_a_state_that_declares_two_kinds_of_number_at_once_is_refused():
    with pytest.raises(ValueError, match="exactly one of"):
        sb.MeasuredComposition(
            name="both", regime="r", source="s", derivation="d",
            masses={"trehalose": 1.0}, factors={"trehalose": 2.0},
        )


def test_the_sillje_masses_are_the_published_per_cell_numbers_divided_by_dry_weight():
    slow = sb.MEASURED_STATES["sillje1999_slow"].masses

    assert slow["trehalose"] == pytest.approx(6.2 * (342.29648 / 2) / 17.4, rel=1e-12)
    assert slow["glycogen"] == pytest.approx(11.0 * 162.1406 / 17.4, rel=1e-12)
    assert sb.MEASURED_STATES["sillje1999_fast"].masses["glycogen"] == pytest.approx(
        2.0 * 162.1406 / 10.6, rel=1e-12
    )


def test_the_cwi_factors_are_the_published_ram_1998_ratios():
    gas1 = sb.MEASURED_STATES["cwi_gas1_ram1998"].factors

    assert gas1["chitin"] == pytest.approx(4.7 / 2.8, rel=1e-12)
    assert gas1["mannan"] == pytest.approx(167.5 / 120.7, rel=1e-12)
    assert gas1["beta_1_3_glucan"] == pytest.approx(65.3 / 76.6, rel=1e-12)


def test_the_hottiger_ratios_are_the_published_grams_of_trehalose_per_gram_of_protein():
    assert sb.MEASURED_STATES["hottiger1987_27C"].protein_ratios == {"trehalose": 0.01}
    assert sb.MEASURED_STATES["hottiger1987_heat_shock"].protein_ratios == {"trehalose": 1.0}


def test_the_lumped_glucan_split_is_declared_asserted_not_passed_off_as_measured():
    for name in ("cwi_gas1_ram1998", "cwi_fks1_ram1998"):
        assert sb.MEASURED_STATES[name].asserted, name


def test_there_is_no_default_stressed_composition_only_cited_states():
    with pytest.raises(KeyError, match="no measured state"):
        sb.apply_measured(None, "heat_shock")


def test_no_futile_cycle_rate_is_supplied_because_none_was_measured():
    assert sb.TREHALOSE_CYCLE_RATE is None


def test_the_three_modelling_closures_stay_declared_as_asserted():
    """The proportional residual, the untouched GAM and the 957 mg total are conventions."""
    assert len(sb.MODELLING_CLOSURES) == 3
    assert all(text.startswith("ASSERTED:") for text in sb.MODELLING_CLOSURES)


# --- renormalisation conserves what it says it conserves ----------------------------------------


def test_pool_mode_conserves_the_carbohydrate_pool_and_touches_nothing_else(yeast_gem):
    shipped_pool = sum(sb.component_masses(yeast_gem).values())
    shipped_protein = sb.pool_masses(yeast_gem)["protein"]

    out = sb.set_component_masses(yeast_gem, {"trehalose": 200.0}, renormalise="pool")

    assert sum(sb.component_masses(out).values()) == pytest.approx(shipped_pool, rel=1e-12)
    assert sb.total_biomass_mass(out) == pytest.approx(sb.total_biomass_mass(yeast_gem), rel=1e-12)
    assert sb.pool_masses(out)["protein"] == pytest.approx(shipped_protein, rel=1e-12)
    assert sb.component_masses(out)["trehalose"] == pytest.approx(200.0, rel=1e-12)


def test_whole_cell_mode_conserves_the_total_and_lets_the_carbohydrate_fraction_move(yeast_gem):
    shipped_pool = sum(sb.component_masses(yeast_gem).values())
    shipped_protein = sb.pool_masses(yeast_gem)["protein"]

    out = sb.set_component_masses(yeast_gem, {"trehalose": 200.0}, renormalise="whole_cell")

    assert sb.total_biomass_mass(out) == pytest.approx(sb.total_biomass_mass(yeast_gem), rel=1e-12)
    assert sum(sb.component_masses(out).values()) > shipped_pool
    assert sb.pool_masses(out)["protein"] < shipped_protein


def test_the_unnamed_components_keep_their_shipped_proportions_to_each_other(yeast_gem):
    shipped = sb.component_masses(yeast_gem)

    out = sb.set_component_masses(yeast_gem, {"trehalose": 200.0}, renormalise="whole_cell")

    moved = sb.component_masses(out)
    assert moved["mannan"] / moved["glycogen"] == pytest.approx(
        shipped["mannan"] / shipped["glycogen"], rel=1e-12)


def test_rescaling_a_pool_scales_its_charge_imbalance_by_exactly_the_mass_factor(
    yeast_gem, thick_wall_model
):
    """The pools ship charge-imbalanced; a proportional rewrite must not make that worse."""
    for rid in ("r_4047", "r_4049", "r_4063"):
        before = yeast_gem.reactions.get_by_id(rid).check_mass_balance()["charge"]
        after = thick_wall_model.reactions.get_by_id(rid).check_mass_balance()["charge"]
        mass_factor = sb.reaction_mass(thick_wall_model, rid) / sb.reaction_mass(yeast_gem, rid)
        assert after / before == pytest.approx(mass_factor, rel=1e-12), rid


def test_the_input_model_is_never_modified(yeast_gem):
    before = sb.component_masses(yeast_gem)

    sb.set_component_masses(yeast_gem, {"trehalose": 300.0}, renormalise="whole_cell")
    sb.set_pool_masses(yeast_gem, {"protein": 300.0})
    sb.pin_trehalose_cycle(yeast_gem, 1.0)

    assert sb.component_masses(yeast_gem) == before
    assert yeast_gem.reactions.get_by_id("r_0194").lower_bound == 0.0


def test_the_rewritten_pseudoreaction_stays_charge_neutral_and_carries_flux(thick_wall_model):
    from cobra.flux_analysis import flux_variability_analysis

    reaction = thick_wall_model.reactions.get_by_id("r_4048")
    assert "charge" not in reaction.check_mass_balance()
    assert sb.reaction_mass(thick_wall_model, "r_4048") == pytest.approx(
        sum(sb.component_masses(thick_wall_model).values()), rel=1e-12
    )

    with thick_wall_model as m:
        m.reactions.get_by_id(GLUCOSE).lower_bound = -10.0
        fva = flux_variability_analysis(
            m, ["r_4048", "r_2111"], fraction_of_optimum=0.0, processes=1)

    assert fva.loc["r_4048", "maximum"] == pytest.approx(0.892827, rel=1e-4)
    assert fva.loc["r_2111", "maximum"] == pytest.approx(fva.loc["r_4048", "maximum"], rel=1e-9)


def test_a_request_that_leaves_no_mass_for_the_rest_is_refused_not_clamped(yeast_gem):
    with pytest.raises(ValueError, match="only"):
        sb.set_component_masses(yeast_gem, {"trehalose": 400.0}, renormalise="pool")


def test_negative_masses_unknown_components_and_unknown_modes_are_refused(yeast_gem):
    with pytest.raises(ValueError, match="negative mass"):
        sb.set_component_masses(yeast_gem, {"trehalose": -1.0})
    with pytest.raises(KeyError, match="unknown carbohydrate components"):
        sb.set_component_masses(yeast_gem, {"mannoprotein": 10.0})
    with pytest.raises(ValueError, match="renormalise must be"):
        sb.set_component_masses(yeast_gem, {"trehalose": 10.0}, renormalise="none")


# --- pools are parameters too, and have no measured default --------------------------------------


def test_a_named_pool_holds_its_mass_and_the_rest_absorb_the_difference(yeast_gem):
    shipped = sb.pool_masses(yeast_gem)

    out = sb.set_pool_masses(yeast_gem, {"protein": shipped["protein"] / 2})

    after = sb.pool_masses(out)
    assert after["protein"] == pytest.approx(shipped["protein"] / 2, rel=1e-12)
    assert sb.total_biomass_mass(out) == pytest.approx(sb.total_biomass_mass(yeast_gem), rel=1e-12)
    assert after["RNA"] / after["DNA"] == pytest.approx(shipped["RNA"] / shipped["DNA"], rel=1e-12)


def test_rescaling_a_pool_leaves_the_components_inside_it_in_proportion(yeast_gem):
    shipped = sb.component_masses(yeast_gem)

    out = sb.set_pool_masses(yeast_gem, {"protein": 232.0122189})

    moved = sb.component_masses(out)
    assert moved["trehalose"] / moved["mannan"] == pytest.approx(
        shipped["trehalose"] / shipped["mannan"], rel=1e-12)
    assert grow(out) == pytest.approx(0.895744, rel=1e-5)


def test_halving_protein_is_worth_less_than_one_percent_of_growth(yeast_gem):
    """The whole 464 mg protein pool is 48.5% of the cell and moving it barely moves growth."""
    out = sb.set_pool_masses(yeast_gem, {"protein": 232.0122189})

    assert (grow(out) - BASE_GROWTH) / BASE_GROWTH == pytest.approx(0.00908, abs=1e-4)


def test_unknown_pools_absent_pools_and_negative_pool_masses_are_refused(yeast_gem):
    with pytest.raises(KeyError, match="unknown biomass pools"):
        sb.set_pool_masses(yeast_gem, {"mannoprotein": 10.0})
    with pytest.raises(ValueError, match="negative mass"):
        sb.set_pool_masses(yeast_gem, {"protein": -1.0})
    with pytest.raises(ValueError, match="only"):
        sb.set_pool_masses(yeast_gem, {"protein": 2000.0})


def test_moving_one_lipid_pool_without_the_other_is_refused_not_silently_lethal(yeast_gem):
    """The bug this guard exists for: an unequal move used to hand back a model at growth 0."""
    shipped = sb.pool_masses(yeast_gem)

    for name in sb.SLIME_COUPLED_POOLS:
        with pytest.raises(ValueError, match="SLIME"):
            sb.set_pool_masses(yeast_gem, {name: shipped[name] * 1.2})


def test_the_two_lipid_pools_move_together_and_then_the_cost_is_finite(yeast_gem):
    """Named in their shipped ratio they are a real parameter: doubling lipid costs 8%."""
    shipped = sb.pool_masses(yeast_gem)
    doubled = sb.set_pool_masses(yeast_gem, {name: shipped[name] * 2.0
                                             for name in sb.SLIME_COUPLED_POOLS})

    assert grow(doubled) == pytest.approx(0.816709, rel=1e-5)
    assert grow(sb.set_pool_masses(yeast_gem, {name: shipped[name] * 0.5
                                               for name in sb.SLIME_COUPLED_POOLS})
               ) == pytest.approx(0.927913, rel=1e-5)


def test_the_unequal_lipid_move_the_guard_refuses_really_would_have_been_lethal(yeast_gem):
    """Why the refusal is not pedantry: a 1% mismatch is enough to zero the objective."""
    mismatched = yeast_gem.copy()
    sb._scale_pool(mismatched, "r_4065", 1.01)

    assert grow(mismatched) == pytest.approx(0.0, abs=1e-9)


def test_a_component_rewrite_keeps_the_two_lipid_pools_locked_together(yeast_gem):
    """set_component_masses uses one factor for every pool, so it cannot break the coupling."""
    shipped = sb.pool_masses(yeast_gem)

    out = sb.set_pool_masses(yeast_gem, {"protein": 232.0122189})
    moved = sb.pool_masses(out)

    assert moved["lipid_backbone"] / shipped["lipid_backbone"] == pytest.approx(
        moved["lipid_chain"] / shipped["lipid_chain"], rel=1e-9)
    assert grow(out) > 0.8


def test_reaction_mass_cannot_see_the_gram_denominated_SLIME_species(yeast_gem):
    """KNOWN DEFECT, recorded not hidden: the lipid pools are read in the wrong unit.

    yeast-GEM's SLIME chain and backbone species are counted in GRAMS -- the chain coefficient in
    r_3975 is palmitate's MW in g/mmol -- so multiplying them by a formula weight is a unit error,
    and s_0694 carries no formula at all. The conserved total is understated by about 33 mg.
    """
    chain = yeast_gem.reactions.get_by_id("r_4065")
    grams = sum(-c for m, c in chain.metabolites.items() if c < 0) * 1000.0
    backbone_gap = -yeast_gem.reactions.get_by_id("r_4063").get_coefficient("s_0694") * 1000.0

    assert sb.reaction_mass(yeast_gem, "r_4065") == pytest.approx(11.2146, rel=1e-4)
    assert grams == pytest.approx(42.7890, rel=1e-4)
    assert backbone_gap == pytest.approx(1.4801, rel=1e-4)
    corrected = sb.total_biomass_mass(yeast_gem) - 11.2146 + grams + backbone_gap
    assert corrected == pytest.approx(989.9713, rel=1e-5)


def test_no_measured_protein_or_lipid_state_is_shipped_because_none_was_found():
    """set_pool_masses is the swept parameter the brief asks for; it gets no default."""
    for record in sb.MEASURED_STATES.values():
        named = set(record.masses or record.factors or record.protein_ratios)
        assert named <= set(sb.CARBOHYDRATE_COMPONENTS), record.name


# --- the protein-ratio fixed point ---------------------------------------------------------------


def test_the_trehalose_to_protein_fixed_point_holds_after_renormalisation(yeast_gem):
    """Hottiger 1987 reports trehalose per gram of protein, and protein itself moves."""
    target = sb.masses_for_protein_ratios(yeast_gem, {"trehalose": 1.0})["trehalose"]
    out = sb.set_component_masses(yeast_gem, {"trehalose": target}, renormalise="whole_cell")

    assert target == pytest.approx(323.121, rel=1e-5)
    assert sb.component_masses(out)["trehalose"] / sb.pool_masses(out)["protein"] == pytest.approx(
        1.0, rel=1e-9
    )


def test_at_a_fixed_pool_the_ratio_is_just_ratio_times_shipped_protein(yeast_gem):
    protein = sb.pool_masses(yeast_gem)["protein"]

    at_pool = sb.masses_for_protein_ratios(yeast_gem, {"trehalose": 0.01}, renormalise="pool")

    assert at_pool["trehalose"] == pytest.approx(0.01 * protein, rel=1e-12)


def test_the_heat_shock_load_does_not_fit_inside_the_carbohydrate_pool_at_all(yeast_gem):
    """464 mg/gDCW of trehalose against a 383 mg pool: at fixed pool the state is impossible."""
    with pytest.raises(ValueError, match="only 383"):
        sb.apply_measured(yeast_gem, "hottiger1987_heat_shock", renormalise="pool")


def test_the_heat_shock_state_puts_trehalose_and_protein_at_the_same_mass(
    yeast_gem, heat_shock_model
):
    masses = sb.component_masses(heat_shock_model)

    assert masses["trehalose"] == pytest.approx(323.1208, rel=1e-5)
    assert sb.pool_masses(heat_shock_model)["protein"] == pytest.approx(323.1208, rel=1e-5)
    assert sb.total_biomass_mass(heat_shock_model) == pytest.approx(
        sb.total_biomass_mass(yeast_gem), rel=1e-12)


# --- (a) the false lethal becomes a finite cost --------------------------------------------------


def test_the_shipped_model_calls_every_carbohydrate_synthase_absolutely_essential(yeast_gem):
    """All six, not just trehalose: -100% each, and every one of them is the coefficient talking."""
    assert grow(yeast_gem) == pytest.approx(BASE_GROWTH, rel=1e-5)

    for name, reaction_id in sb.COMPONENT_SYNTHASES.items():
        assert grow_without(yeast_gem, reaction_id) == pytest.approx(0.0, abs=1e-9), name


def test_driving_a_component_out_of_the_composition_makes_its_synthase_free(zeroed_models):
    """The false lethal in full: every one of the six costs exactly nothing once it is not demanded."""
    for name, model in zeroed_models.items():
        unblocked = grow(model)
        assert grow_without(model, sb.COMPONENT_SYNTHASES[name]) == pytest.approx(
            unblocked, rel=1e-9), name


def test_the_composition_costs_themselves_range_from_0_05_to_2_4_percent(zeroed_models):
    """What each component is actually worth, in place of six identical -100% readings."""
    costs = {
        name: (grow(model) - BASE_GROWTH) / BASE_GROWTH for name, model in zeroed_models.items()
    }

    assert costs == pytest.approx(
        {
            "beta_1_3_glucan": -0.0240082,
            "beta_1_6_glucan": -0.0074477,
            "chitin": -0.0004855,
            "glycogen": -0.0109376,
            "mannan": -0.0226714,
            "trehalose": -0.0106964,
        },
        abs=1e-5,
    )


def test_mannan_has_no_synthase_at_all_so_its_ablation_is_of_protein_mannosylation(yeast_gem):
    """Honest limit: s_1107 is filled only by transport from what r_0362 makes in the ER."""
    producers = [
        r.id for r in yeast_gem.metabolites.get_by_id("s_1107").reactions
        if r.get_coefficient("s_1107") > 0
    ]

    assert producers == []
    assert sb.COMPONENT_SYNTHASES["mannan"] == "r_0362"
    assert "mannosyltransferase" in yeast_gem.reactions.get_by_id("r_0362").name


def test_at_the_measured_deletant_composition_neither_knockout_costs_anything(deletant_model):
    """Silljé 1999: tps1 gsy1 gsy2 grew at the wild-type rate. -100% was the coefficient talking."""
    unblocked = grow(deletant_model)

    assert grow_without(deletant_model, TREHALOSE_SYNTHASE) == pytest.approx(unblocked, rel=1e-9)
    assert grow_without(deletant_model, GLYCOGEN_SYNTHASE) == pytest.approx(unblocked, rel=1e-9)


def test_the_deletant_composition_itself_costs_2_27_percent_not_100(deletant_model):
    growth = grow(deletant_model)

    assert growth == pytest.approx(0.867549, rel=1e-5)
    assert (growth - BASE_GROWTH) / BASE_GROWTH == pytest.approx(-0.02268, abs=1e-4)


def test_the_lethal_returns_the_moment_the_composition_demands_trehalose_again(yeast_gem):
    """Honest limit: Yeast9 has one trehalose source, so a nonzero coefficient is still a cliff."""
    slow = sb.apply_measured(yeast_gem, "sillje1999_slow", renormalise="whole_cell")

    assert grow(slow) > 0.8
    assert grow_without(slow, TREHALOSE_SYNTHASE) == 0.0


def test_accumulating_trehalose_is_priced_as_a_gain_under_both_closures(yeast_gem):
    """6.29 mmol glucose/g against 6.74 for the hexose polymers: trehalose is the cheapest gram."""
    at_pool = [
        grow(sb.set_component_masses(yeast_gem, {"trehalose": t}, renormalise="pool"))
        for t in (0.0, 46.7406, 100.0, 200.0)
    ]
    at_cell = [
        grow(sb.set_component_masses(yeast_gem, {"trehalose": t}, renormalise="whole_cell"))
        for t in (0.0, 46.7406, 100.0, 200.0)
    ]

    assert at_pool == sorted(at_pool) and at_cell == sorted(at_cell)
    assert at_pool == pytest.approx([0.885963, 0.887685, 0.889656, 0.893381], rel=1e-5)
    assert at_cell == pytest.approx([0.878190, 0.887685, 0.898758, 0.920312], rel=1e-5)


def test_the_measured_heat_shock_load_is_a_seven_percent_growth_GAIN_not_a_cost(heat_shock_model):
    """The composition arm cannot charge for trehalose accumulation. This is the refutation."""
    growth = grow(heat_shock_model)

    assert growth == pytest.approx(0.948313, rel=1e-5)
    assert (growth - BASE_GROWTH) / BASE_GROWTH == pytest.approx(0.06830, abs=1e-4)


def test_the_cost_of_heat_shock_trehalose_is_the_turnover_and_the_gem_prices_that(yeast_gem):
    """Hottiger's own title is futile cycling; forcing trehalase makes synthesis follow."""
    costs = {
        rate: (grow(sb.pin_trehalose_cycle(yeast_gem, rate)) - BASE_GROWTH) / BASE_GROWTH
        for rate in (0.0, 0.1, 1.0, 5.0)
    }

    assert costs == pytest.approx(
        {0.0: 0.0, 0.1: -0.0013946, 1.0: -0.0139455, 5.0: -0.0697277}, abs=1e-6)


def test_the_forced_cycle_is_a_cycle_and_not_an_import(yeast_gem):
    cycled = sb.pin_trehalose_cycle(yeast_gem, 1.0)

    with cycled as m:
        m.reactions.get_by_id(GLUCOSE).lower_bound = -10.0
        solution = m.optimize()

    assert solution.fluxes["r_0194"] == pytest.approx(1.0, rel=1e-6)
    assert solution.fluxes["r_1051"] == pytest.approx(
        1.0 + 0.13655 * solution.objective_value, rel=1e-6)
    assert solution.fluxes["r_1650"] == pytest.approx(0.0, abs=1e-9)


def test_pinning_the_cycle_refuses_a_negative_rate_and_a_model_that_lacks_the_reaction(yeast_gem):
    with pytest.raises(ValueError, match="non-negative"):
        sb.pin_trehalose_cycle(yeast_gem, -1.0)
    with pytest.raises(KeyError, match="r_0194"):
        sb.pin_trehalose_cycle(cobra.Model("empty"), 1.0)


def test_pinning_the_cycle_refuses_a_medium_that_could_feed_it(yeast_gem):
    """Without this the forced trehalase flux is an import, and the ATP cost never appears."""
    fed = yeast_gem.copy()
    fed.reactions.get_by_id("r_1650").lower_bound = -10.0

    with pytest.raises(ValueError, match="allows uptake"):
        sb.pin_trehalose_cycle(fed, 1.0)


# --- (b) the growth cost of a thickened wall -----------------------------------------------------


def test_the_measured_cwi_response_thickens_the_wall_by_8_6_percent(yeast_gem, thick_wall_model):
    shipped = sum(sb.component_masses(yeast_gem)[k] for k in sb.WALL_COMPONENTS)
    thickened = sum(sb.component_masses(thick_wall_model)[k] for k in sb.WALL_COMPONENTS)

    assert thickened / shipped == pytest.approx(1.0855, rel=1e-3)
    assert sb.component_masses(thick_wall_model)["chitin"] == pytest.approx(8.0527, rel=1e-4)


def test_the_cwi_wall_costs_0_07_percent_against_the_carbohydrate_pool_and_gains_0_58_against_the_cell(
    yeast_gem, thick_wall_model
):
    at_fixed_pool = grow(sb.apply_measured(yeast_gem, "cwi_gas1_ram1998", renormalise="pool"))

    assert (at_fixed_pool - BASE_GROWTH) / BASE_GROWTH == pytest.approx(-0.00069, abs=5e-5)
    assert (grow(thick_wall_model) - BASE_GROWTH) / BASE_GROWTH == pytest.approx(0.00579, abs=5e-5)


def test_the_fks1_wall_is_thinner_and_costs_0_40_percent_against_the_cell(yeast_gem):
    """The other measured CWI lesion: less glucan, more mannan, and the sign flips."""
    fks1 = sb.apply_measured(yeast_gem, "cwi_fks1_ram1998", renormalise="whole_cell")
    wall = sum(sb.component_masses(fks1)[k] for k in sb.WALL_COMPONENTS)
    shipped = sum(sb.component_masses(yeast_gem)[k] for k in sb.WALL_COMPONENTS)

    assert wall / shipped == pytest.approx(0.9472, rel=1e-3)
    assert (grow(fks1) - BASE_GROWTH) / BASE_GROWTH == pytest.approx(-0.00404, abs=5e-5)


def test_chitin_is_the_one_carbohydrate_the_cell_pays_extra_for(yeast_gem):
    """Chitin is 13% dearer per gram, so at a fixed pool a chitin-rich wall is a real cost."""
    chitin_rich = sb.scale_components(yeast_gem, {"chitin": 5.0}, renormalise="pool")

    assert (grow(chitin_rich) - BASE_GROWTH) / BASE_GROWTH == pytest.approx(-0.00161, abs=5e-5)
    assert sb.component_masses(chitin_rich)["chitin"] == pytest.approx(5 * 4.7974, rel=1e-4)


# --- the regime is load-bearing, so every number above is repeated in the repo's own -------------


def test_the_same_states_at_the_reference_aerobic_batch_regime(yeast_gem, thick_wall_model,
                                                               heat_shock_model):
    """glucose 11.1 / oxygen 3.7 (van Hoek 1998). Composition costs barely move; the cycle's 6x."""
    def at_reference(model):
        return (grow(model, 11.1, 3.7) - REFERENCE_GROWTH) / REFERENCE_GROWTH

    assert grow(yeast_gem, 11.1, 3.7) == pytest.approx(REFERENCE_GROWTH, rel=1e-6)
    assert at_reference(thick_wall_model) == pytest.approx(0.00565, abs=5e-5)
    assert at_reference(heat_shock_model) == pytest.approx(0.06573, abs=1e-4)
    assert at_reference(sb.pin_trehalose_cycle(yeast_gem, 1.0)) == pytest.approx(-0.08307, abs=1e-4)


# --- the other vendored model --------------------------------------------------------------------


def test_the_composition_reads_on_gecko_too_and_that_model_has_no_chitin(ec_yeast_gem):
    """ecYeastGEM_batch is built on yeast-GEM 8.3.4: different coefficients, chitin absent."""
    masses = sb.component_masses(ec_yeast_gem)

    assert masses["chitin"] == 0.0
    assert masses["trehalose"] == pytest.approx(43.2854, rel=1e-4)
    assert sb.total_biomass_mass(ec_yeast_gem) == pytest.approx(979.2411, rel=1e-6)


def test_the_futile_cycle_refuses_gecko_rather_than_guessing_which_split_reaction_to_pin(
    ec_yeast_gem,
):
    with pytest.raises(KeyError, match="r_0194"):
        sb.pin_trehalose_cycle(ec_yeast_gem, 1.0)
