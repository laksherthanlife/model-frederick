"""The proton budget has to bind, and a weak acid has to be billed to Pma1 at 1 ATP.

Model-touching tests are marked integration: the SBML parse costs about 18 s and is done once,
through the session-scoped ``yeast_gem`` fixture in conftest, which resolves the model with
:mod:`ystwin.paths`. An earlier draft of this file walked up the directory tree looking for
``dente/data/yeast-GEM.xml`` in a sibling checkout. That is exactly the resolve-by-convention
that `audit_reproducibility.py` refuses, and it silently read a DIFFERENT model from the one
`data/gem/yeast-GEM.xml.gz` pins.
"""
from __future__ import annotations

import cobra
import pytest

from ystwin.fba import stress_ph as sp
from ystwin.fba.physiology import aerobic_batch_constraints
from ystwin.fba.solver import FVA_PROCESSES, configure

GLUCOSE, OXYGEN = "r_1714", "r_1992"


@pytest.fixture(scope="module")
def _gem_template(yeast_gem_factory):
    """Yeast9 with the pinned solver. Every test scopes its edits with ``with model:``."""
    yeast_gem = yeast_gem_factory()
    configure(yeast_gem)
    for rid in (sp.PYROPHOSPHATE_EXCHANGE, sp.PROTON_EXCHANGE, sp.PMA1_REACTION):
        if not yeast_gem.reactions.has_id(rid):
            pytest.skip(f"this model does not carry {rid}")
    return yeast_gem


@pytest.fixture
def gem(_gem_template, model_copy):
    return model_copy(_gem_template)


@pytest.fixture(scope="module")
def _acid_models_template(_gem_template, model_copy):
    return {key: sp.add_weak_acid_uncoupling(model_copy(_gem_template), key)
            for key in sp.WEAK_ACIDS}


@pytest.fixture
def acid_models(_acid_models_template, model_copy):
    return {key: model_copy(model) for key, model in _acid_models_template.items()}


def _glc10(model):
    model.reactions.get_by_id(GLUCOSE).lower_bound = -10.0
    model.reactions.get_by_id(OXYGEN).lower_bound = -1000.0
    return model


def _anaerobic(model):
    return sp.anaerobic_constraints(model, glucose_uptake=10.0)


#: The four regimes every proton-budget number in this module is reported in.
REGIMES = {
    "shipped default": lambda m: m,
    "glucose<=10, O2 free": _glc10,
    "REFERENCE_AEROBIC_BATCH": aerobic_batch_constraints,
    "anaerobic glucose 10": _anaerobic,
}


# ---------------------------------------------------------------- chemistry, no model needed

def test_henderson_hasselbalch_is_half_at_the_pka():
    assert sp.undissociated_fraction(4.757, 4.757) == pytest.approx(0.5)


def test_undissociated_fraction_falls_monotonically_with_ph():
    fractions = [sp.undissociated_fraction(ph, 4.757) for ph in (3.0, 4.0, 5.0, 6.0, 7.0)]
    assert fractions == sorted(fractions, reverse=True)
    assert fractions[0] > 0.98 and fractions[-1] < 0.01


def test_cytosol_is_far_enough_above_every_pka_that_dissociation_is_one_way():
    for acid in sp.WEAK_ACIDS.values():
        assert sp.undissociated_fraction(sp.CYTOSOLIC_PH, acid.pka) < 0.01


def test_influx_is_linear_in_concentration_and_scales_with_the_undissociated_fraction():
    at_ph4 = sp.weak_acid_influx(50.0, 4.0, 4.757, permeability=0.1)
    assert sp.weak_acid_influx(100.0, 4.0, 4.757, permeability=0.1) == pytest.approx(2 * at_ph4)
    assert at_ph4 == pytest.approx(0.1 * 50.0 * sp.undissociated_fraction(4.0, 4.757))


def test_registry_carries_the_cited_pka_values():
    assert sp.WEAK_ACIDS["acetic"].pka == 4.757
    assert sp.WEAK_ACIDS["propionic"].pka == 4.874
    assert sp.WEAK_ACIDS["lactic"].pka == 3.86
    for acid in sp.WEAK_ACIDS.values():
        assert acid.source


def test_deprotonation_removes_exactly_one_hydrogen():
    assert sp._deprotonate("C2H4O2") == "C2H3O2"
    assert sp._deprotonate("C3H6O3") == "C3H5O3"
    assert sp._deprotonate("CH2O2") == "CHO2"
    with pytest.raises(ValueError):
        sp._deprotonate("CO2")


# ---------------------------------------------------------------- the additions themselves

@pytest.mark.integration
def test_every_added_internal_reaction_is_mass_and_charge_balanced(acid_models):
    for key, model in acid_models.items():
        for rid in sp.acid_reaction_ids(sp.WEAK_ACIDS[key]):
            rxn = model.reactions.get_by_id(rid)
            if rxn.boundary:
                continue
            assert rxn.check_mass_balance() == {}, f"{rid} is unbalanced"


@pytest.mark.integration
def test_every_added_reaction_can_carry_flux(acid_models):
    """FVA at fraction_of_optimum=0 with the acid exchange merely OPENED, not pinned.

    Pinning the exchange at the load makes its own range [-load, -load] trivially non-zero,
    which is a test of set_acid_load rather than of whether the route can carry anything.
    """
    from cobra.flux_analysis import flux_variability_analysis

    for key, model in acid_models.items():
        spec = sp.WEAK_ACIDS[key]
        ids = list(sp.acid_reaction_ids(spec))
        with model as m:
            _glc10(m)
            m.reactions.get_by_id(f"EX_{spec.acid_e_id}").lower_bound = -1000.0
            ranges = flux_variability_analysis(m, reaction_list=ids, fraction_of_optimum=0.0,
                                               processes=FVA_PROCESSES)
        for rid in ids:
            span = max(abs(ranges.loc[rid, "minimum"]), abs(ranges.loc[rid, "maximum"]))
            assert span > 1e-6, f"{rid} cannot carry flux"


@pytest.mark.integration
def test_an_installed_acid_that_is_not_loaded_is_a_structurally_dead_route(acid_models):
    """FVA, not a growth comparison, is the authoritative check that nothing new can flow.

    The entry route is [0, 0] with the exchange shut. Lactate's anion efflux is NOT, and
    should not be: it moves the native s_0063/s_0064 pair, so the model could already do it.
    """
    from cobra.flux_analysis import flux_variability_analysis

    for key, model in acid_models.items():
        spec = sp.WEAK_ACIDS[key]
        entry = [f"EX_{spec.acid_e_id}", f"{key.upper()}_DIFF", f"{key.upper()}_DISS"]
        with model as m:
            _glc10(m)
            ranges = flux_variability_analysis(m, reaction_list=entry, fraction_of_optimum=0.0,
                                               processes=FVA_PROCESSES)
        for rid in entry:
            assert ranges.loc[rid, "minimum"] == pytest.approx(0.0, abs=1e-9), rid
            assert ranges.loc[rid, "maximum"] == pytest.approx(0.0, abs=1e-9), rid


@pytest.mark.integration
def test_installing_an_acid_that_is_not_loaded_does_not_move_growth(gem, acid_models):
    """Tolerance is 1e-5, not 1e-9, and the reason is the solver rather than the model.

    The route is proved dead above by FVA. Adding three to five columns to a degenerate LP
    still changes which optimal vertex GLPK lands on, by about 1e-6 relative here.
    """
    for name, setup in REGIMES.items():
        with gem as m:
            setup(m)
            base = m.slim_optimize()
        for key, model in acid_models.items():
            with model as m:
                setup(m)
                assert m.slim_optimize() == pytest.approx(base, rel=1e-5), f"{key} in {name}"


@pytest.mark.integration
def test_the_source_model_is_left_untouched(gem, acid_models):
    assert not gem.metabolites.has_id("acetic_ha_e")
    assert gem.reactions.get_by_id(sp.PYROPHOSPHATE_EXCHANGE).bounds == (0.0, 1000.0)
    assert gem.reactions.get_by_id(sp.PROTON_EXCHANGE).bounds == (-1000.0, 1000.0)
    assert gem.reactions.get_by_id(sp.PMA1_REACTION).gene_reaction_rule.endswith("YER005W")


@pytest.mark.integration
def test_installing_twice_is_refused(gem):
    with pytest.raises(ValueError):
        sp.add_weak_acid_uncoupling(sp.add_weak_acid_uncoupling(gem, "acetic"), "acetic")


@pytest.mark.integration
def test_the_pump_the_bill_is_charged_to_already_moves_one_proton_per_atp(gem):
    """r_0227 deposits the hydrolysis proton OUTSIDE, so it is a 1:1 pump as it ships."""
    rxn = gem.reactions.get_by_id(sp.PMA1_REACTION)
    assert rxn.check_mass_balance() == {}
    assert rxn.metabolites[gem.metabolites.get_by_id("s_0434")] == -1   # ATP
    assert rxn.metabolites[gem.metabolites.get_by_id("s_0796")] == 1    # extracellular H+
    assert gem.metabolites.get_by_id("s_0794") not in rxn.metabolites   # no cytosolic H+ term


# ---------------------------------------------------------------- the proton budget

@pytest.mark.integration
def test_shipped_yeast9_under_prices_a_cytosolic_proton_in_every_regime(gem):
    for name, setup in REGIMES.items():
        with gem as m:
            setup(m)
            price = sp.proton_price(m)
        assert price < 0.95, f"{name} already prices a proton at {price}"


@pytest.mark.integration
def test_closing_the_budget_prices_a_proton_at_exactly_one_atp_in_every_regime(gem):
    for name, setup in REGIMES.items():
        with gem as m:
            setup(m)
            sp.constrain_proton_budget(m)
            assert sp.proton_price(m) == pytest.approx(1.0, abs=1e-3), name


@pytest.mark.integration
def test_shipped_yeast9_excretes_pyrophosphate_instead_of_running_its_proton_pump(gem):
    with gem as m:
        _glc10(m)
        sol = m.optimize()
        assert sol.fluxes[sp.PMA1_REACTION] == pytest.approx(0.0, abs=1e-6)
        assert sol.fluxes[sp.PYROPHOSPHATE_EXCHANGE] > 1.0
        sp.constrain_proton_budget(m)
        sol = m.optimize()
        assert sol.fluxes[sp.PYROPHOSPHATE_EXCHANGE] == pytest.approx(0.0, abs=1e-6)
        assert sol.fluxes[sp.PMA1_REACTION] > 1.0


@pytest.mark.integration
def test_closing_the_budget_costs_under_two_percent_of_growth_in_every_regime(gem):
    """The measured costs are 0.726 / 0.726 / 1.669 / 1.527 percent, in REGIMES order."""
    for name, setup in REGIMES.items():
        with gem as m:
            setup(m)
            report = sp.proton_budget_cost(m)
        assert 0.0 < report["relative_growth_cost"] < 0.02, f"{name}: {report}"
        assert report["pma1_before"] == pytest.approx(0.0, abs=1e-6)
        assert report["pma1_after"] > 0.1


@pytest.mark.integration
def test_the_escape_route_diagnostic_names_pyrophosphate_then_the_proton_itself(gem):
    with gem as m:
        _glc10(m)
        leaky = sp.proton_escape_routes(m, 1.0)
        assert sp.PYROPHOSPHATE_EXCHANGE in {rid for rid, _, _ in leaky[:3]}
        sp.constrain_proton_budget(m)
        tight = sp.proton_escape_routes(m, 1.0)
    assert tight[0][0] == sp.PROTON_EXCHANGE
    assert tight[0][2] == pytest.approx(1.0, abs=0.05)


@pytest.mark.integration
def test_the_gpr_fix_makes_the_pump_answer_to_pma1(gem):
    with gem as m:
        _glc10(m)
        sp.constrain_proton_budget(m)
        for gene in ("YGL008C", "YPL036W"):
            m.genes.get_by_id(gene).knock_out()
        assert m.reactions.get_by_id(sp.PMA1_REACTION).bounds == (0, 0)
    with gem as m:
        _glc10(m)
        for gene in ("YGL008C", "YPL036W"):
            m.genes.get_by_id(gene).knock_out()
        assert m.reactions.get_by_id(sp.PMA1_REACTION).upper_bound > 0


@pytest.mark.integration
def test_the_report_says_what_it_changed(gem):
    with gem as m:
        report = sp.constrain_proton_budget(m)
    assert sp.PYROPHOSPHATE_EXCHANGE in report.changed
    assert sp.PROTON_EXCHANGE in report.changed
    assert report.missing == ()
    assert sp.PROTON_EXCHANGE in report.summary()


# ---------------------------------------------------------------- Verduyn's mechanism

@pytest.mark.integration
def test_a_pure_uncoupler_costs_exactly_one_atp_per_mol_once_the_budget_is_closed(acid_models):
    for key in ("propionic", "lactic"):
        with acid_models[key] as m:
            _anaerobic(m)
            sp.constrain_proton_budget(m)
            assert sp.atp_per_acid(m, key, 2.0, 5.0) == pytest.approx(1.0, abs=0.02)


@pytest.mark.integration
def test_the_shipped_budget_under_charges_that_acid_threefold(acid_models):
    with acid_models["propionic"] as m:
        _glc10(m)
        loose = sp.atp_per_acid(m, "propionic", 2.0, 5.0)
    with acid_models["propionic"] as m:
        _glc10(m)
        sp.constrain_proton_budget(m)
        tight = sp.atp_per_acid(m, "propionic", 2.0, 5.0)
    assert loose == pytest.approx(1.0 / 3.0, abs=0.02)
    assert tight == pytest.approx(1.0, abs=0.02)
    assert tight / loose == pytest.approx(3.0, abs=0.1)


@pytest.mark.integration
def test_aerobically_acetate_and_lactate_pay_for_themselves_and_propionate_does_not(acid_models):
    """The honest limit of the claim: two of the three acids are carbon sources with O2 free."""
    prices = {}
    for key in sp.WEAK_ACIDS:
        with acid_models[key] as m:
            _glc10(m)
            sp.constrain_proton_budget(m)
            prices[key] = sp.atp_per_acid(m, key, 2.0, 5.0)
    assert prices["acetic"] < -1.0 and prices["lactic"] < -1.0
    assert prices["propionic"] == pytest.approx(1.0, abs=0.02)


@pytest.mark.integration
def test_the_pma1_bill_rises_one_for_one_with_the_acid_load(acid_models):
    """Verduyn 1990 measured this linearity between acid added and the energy to expel it."""
    with acid_models["acetic"] as m:
        _anaerobic(m)
        sp.constrain_proton_budget(m)
        rows = sp.acid_load_response(m, "acetic", [2.0, 3.0, 4.0, 5.0])
    steps = [rows[i + 1]["pma1_flux"] - rows[i]["pma1_flux"] for i in range(3)]
    for step in steps:
        assert step == pytest.approx(1.0, abs=0.01)
    growth = [r["growth"] for r in rows]
    assert growth == sorted(growth, reverse=True)


@pytest.mark.integration
def test_the_anaerobic_regime_reproduces_the_measured_biomass_yield(gem):
    """Verduyn 1990 (PMID 1975265) measured 0.10 g biomass per g glucose at D = 0.10 /h.

    Growth is pinned at their dilution rate and glucose minimised, which is the quantity they
    measured; the model answers 0.1068 g/g, +6.8%.
    """
    with gem as m:
        sp.anaerobic_constraints(m, glucose_uptake=1000.0)
        sp.constrain_proton_budget(m)
        report = sp.verduyn_anaerobic_yield(m)
    assert abs(report["relative_error"]) < 0.15
    assert report["ethanol"] > 1.0, "an anaerobic culture that makes no ethanol is not one"


@pytest.mark.integration
def test_the_models_marginal_atp_yield_is_below_the_measured_yatp_and_says_so(gem):
    """Verduyn's YATP is about 16 g/mol; this model's marginal value is 10.3, so every ATP
    cost here is a 0.65x under-statement of the growth penalty they measured."""
    with gem as m:
        _anaerobic(m)
        sp.constrain_proton_budget(m)
        g_per_mol = sp.marginal_atp_price(m) * 1000.0
    assert 8.0 < g_per_mol < 13.0
    assert g_per_mol < sp.VERDUYN_YATP
    assert "10.32 g biomass/mol ATP" in sp.PROVENANCE[
        "weak-acid uncoupling costs ATP, linearly in acid added"]


@pytest.mark.integration
def test_the_anaerobic_vitamin_bound_is_not_what_limits_growth(gem):
    with gem as m:
        sp.anaerobic_constraints(m, glucose_uptake=10.0, vitamin_uptake=0.01)
        tight = m.slim_optimize()
    with gem as m:
        sp.anaerobic_constraints(m, glucose_uptake=10.0, vitamin_uptake=1.0)
        loose = m.slim_optimize()
    assert loose == pytest.approx(tight, rel=1e-6)


# ---------------------------------------------------------------- the curve and its depth

@pytest.mark.integration
def test_the_ph_curve_costs_growth_as_the_medium_acidifies(acid_models):
    with acid_models["propionic"] as m:
        _anaerobic(m)
        sp.constrain_proton_budget(m)
        rows = sp.ph_response_curve(m, "propionic", 50.0, [3.0, 4.0, 5.0, 6.0, 7.0],
                                    permeability=0.1)
    fractions = [r["undissociated_fraction"] for r in rows]
    assert fractions == sorted(fractions, reverse=True)
    assert [r["growth"] for r in rows] == sorted(r["growth"] for r in rows)
    assert rows[0]["pma1_flux"] > 5 * rows[-1]["pma1_flux"]


@pytest.mark.integration
def test_the_ph_curve_reports_an_atp_cost_that_matches_the_pump_it_is_charged_to(acid_models):
    """One ATP per proton, so the ATP burden read off growth IS the acid that came in.

    The Pma1 bill is that plus a growth-dependent baseline -- biomass synthesis makes protons
    too -- so the pump flux is above the load and does not track it one-for-one across a pH
    curve, only within one.
    """
    with acid_models["propionic"] as m:
        _anaerobic(m)
        sp.constrain_proton_budget(m)
        rows = sp.ph_response_curve(m, "propionic", 50.0, [3.0, 5.0, 7.0], permeability=0.1)
    for row in rows:
        assert row["atp_cost"] == row["pma1_flux"]
        assert row["atp_burden"] == pytest.approx(row["acid_influx"], rel=2e-3)
        assert row["pma1_flux"] > row["acid_influx"]
        assert row["proton_export"] > 0.0


@pytest.mark.integration
def test_closing_the_budget_changes_the_ph_answer(acid_models):
    """The whole point: the same pH now moves growth by a different amount."""
    curves = {}
    for fix in (False, True):
        with acid_models["propionic"] as m:
            _anaerobic(m)
            if fix:
                sp.constrain_proton_budget(m)
            curves[fix] = sp.ph_response_curve(m, "propionic", 50.0, [3.0, 5.5, 7.0],
                                               permeability=0.1)
    for i, ph in enumerate((3.0, 5.5)):
        loose = 1.0 - curves[False][i]["relative_growth"]
        tight = 1.0 - curves[True][i]["relative_growth"]
        assert tight > loose * 1.15, f"pH {ph}: {loose} -> {tight}"
    # At pH 5.5 the shipped model calls the load free; the closed budget charges for it.
    assert curves[False][1]["pma1_flux"] == pytest.approx(0.0, abs=1e-6)
    assert curves[True][1]["pma1_flux"] > 1.0


@pytest.mark.integration
def test_one_acid_load_reaches_hundreds_of_reactions_across_the_compartments(acid_models):
    model = acid_models["propionic"]
    with model as m:
        _anaerobic(m)
        sp.constrain_proton_budget(m)
        moved = sp.proton_load_depth(m, "propionic", 3.0)
    ids = {rid for rid, _, _ in moved}
    assert len(ids) > 200
    installed = set(sp.acid_reaction_ids(sp.WEAK_ACIDS["propionic"]))
    assert installed <= ids
    assert {sp.PMA1_REACTION, sp.PROTON_EXCHANGE} <= ids
    # Depth: the pump is not the end of it -- fermentation and nitrogen assimilation move too.
    assert {"r_1761", "r_0959", "r_0569", "r_1654"} <= ids
    compartments = {met.compartment for rid in ids
                    for met in model.reactions.get_by_id(rid).metabolites}
    assert len(compartments) >= 6


@pytest.mark.integration
def test_the_depth_is_ordered_with_the_installed_route_first(acid_models):
    with acid_models["propionic"] as m:
        _anaerobic(m)
        sp.constrain_proton_budget(m)
        moved = sp.proton_load_depth(m, "propionic", 3.0)
    top = [rid for rid, _, _ in moved[:7]]
    assert set(sp.acid_reaction_ids(sp.WEAK_ACIDS["propionic"])) <= set(top)
    assert sp.PMA1_REACTION in top and sp.PROTON_EXCHANGE in top


@pytest.mark.integration
def test_the_proton_price_does_not_measure_the_pumps_stoichiometry(gem):
    """It returns 1.0 for a 2:1 and a 3:1 pump, and with the pump deleted. Not a coupling ratio.

    Its denominator is an ATP-hydrolysis probe that itself releases one cytosolic proton, so
    once that proton dominates the probe's cost the ratio is 1.0 whatever r_0227 does.
    """
    for n in (1, 2, 3):
        with gem as m:
            _glc10(m)
            sp.constrain_proton_budget(m)
            m.reactions.get_by_id(sp.PMA1_REACTION).add_metabolites(
                {m.metabolites.get_by_id("s_0796"): n - 1})
            assert sp.proton_price(m) == pytest.approx(1.0, abs=1e-3), f"{n} H+/ATP"
    with gem as m:
        _glc10(m)
        sp.constrain_proton_budget(m)
        m.reactions.get_by_id(sp.PMA1_REACTION).bounds = (0.0, 0.0)
        assert sp.proton_price(m) == pytest.approx(1.0, abs=1e-3)


@pytest.mark.integration
def test_deleting_pma1_changes_neither_growth_nor_the_price_of_the_acid(acid_models):
    """The pump is a correlate at the optimum, not the payer, so atp_cost must not be quoted."""
    for setup in (_glc10, _anaerobic):
        with acid_models["propionic"] as m:
            setup(m)
            sp.constrain_proton_budget(m)
            with_pump, growth_on = sp.atp_per_acid(m, "propionic", 2.0, 5.0), m.slim_optimize()
        with acid_models["propionic"] as m:
            setup(m)
            sp.constrain_proton_budget(m)
            m.reactions.get_by_id(sp.PMA1_REACTION).bounds = (0.0, 0.0)
            no_pump, growth_off = sp.atp_per_acid(m, "propionic", 2.0, 5.0), m.slim_optimize()
            rows = sp.acid_load_response(m, "propionic", [3.0])
        assert growth_off == pytest.approx(growth_on, rel=1e-9)
        assert no_pump == pytest.approx(with_pump, abs=1e-3) == pytest.approx(1.0, abs=0.02)
        assert rows[0]["pma1_flux"] == pytest.approx(0.0, abs=1e-9)


@pytest.mark.integration
def test_atp_per_acid_counts_protons_delivered_not_atp_spent(acid_models):
    """Doctor the dissociation step to release 2 or 3 H+ and the answer is 2 or 3."""
    for n in (1, 2, 3):
        with acid_models["propionic"] as m:
            _anaerobic(m)
            sp.constrain_proton_budget(m)
            m.reactions.get_by_id("PROPIONIC_DISS").add_metabolites(
                {m.metabolites.get_by_id("s_0794"): n - 1})
            assert sp.atp_per_acid(m, "propionic", 2.0, 5.0) == pytest.approx(n, abs=0.02)


@pytest.mark.integration
def test_nothing_installed_or_constrained_creates_free_atp_or_free_growth(gem, acid_models):
    """The classic failure when adding exchanges: shut every boundary, ask for ATP from nothing."""
    maintenance = [r.id for r in gem.reactions if r.lower_bound > 0 and not r.boundary]

    def free_atp(model):
        with model as m:
            for r in m.reactions:
                if r.boundary:
                    r.bounds = (0.0, 0.0)
            for rid in maintenance:
                m.reactions.get_by_id(rid).lower_bound = 0.0
            drain = cobra.Reaction("_egc", lower_bound=0.0, upper_bound=1000.0)
            m.add_reactions([drain])
            drain.add_metabolites({m.metabolites.get_by_id(i): v for i, v in (
                ("s_0434", -1), ("s_0803", -1), ("s_0394", 1), ("s_1322", 1), ("s_0794", 1))})
            m.objective = drain
            atp = m.slim_optimize()
            m.objective = m.reactions.get_by_id(sp.BIOMASS_REACTION)
            return float(atp), float(m.slim_optimize())

    assert free_atp(gem) == (pytest.approx(0.0, abs=1e-6), pytest.approx(0.0, abs=1e-6))
    for key, model in acid_models.items():
        assert free_atp(model) == (pytest.approx(0.0, abs=1e-6), pytest.approx(0.0, abs=1e-6)), key
        with model as m:
            m.reactions.get_by_id(f"EX_{sp.WEAK_ACIDS[key].acid_e_id}").lower_bound = -1000.0
            assert free_atp(m) == (pytest.approx(0.0, abs=1e-6), pytest.approx(0.0, abs=1e-6)), key


def test_every_number_in_the_module_has_a_stated_provenance():
    text = " ".join(sp.PROVENANCE.values())
    assert "PLACEHOLDER" in sp.PROVENANCE["DEFAULT_PERMEABILITY = 1.0"]
    for author in ("Verduyn", "Malpartida", "Perlin", "Venema", "Orij", "Casal", "Martell"):
        assert author in text, author
    for key in ("growth cost of closing the budget, by regime",
                "one pH input moves 441 reactions",
                "ATP charged per mol of acid, by acid and regime"):
        assert sp.PROVENANCE[key].startswith("MEASURED")
    for acid in sp.WEAK_ACIDS.values():
        for rid, side, reason in acid.escapes:
            assert side in ("upper", "lower") and reason
