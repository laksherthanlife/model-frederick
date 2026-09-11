"""Oxidative stress becomes dosable, and the price list changes.

Every number here was measured on yeast-GEM v9.0.2 through GLPK with cobra 0.32.1, at
glucose <= 10 with oxygen free unless a test says otherwise. The claim under test
throughout is one thing: before this module, ablating a peroxide defence was free and
forcing catalase priced the MANUFACTURE of peroxide; after it, a dose arrives from the
medium and the defences cost growth for clearing it.
"""

import cobra
import pytest

from ystwin.fba.solver import FVA_PROCESSES, PINNED_TOLERANCE, configure
from ystwin.fba.stress_ros import (
    CITATIONS,
    CLEARANCE_REACTIONS,
    CLEARANCE_ROUTES,
    MITOCHONDRIAL_LEAK,
    NADPH_OXIDASE,
    PEROXIDE_EXCHANGE,
    PEROXIDE_TRANSPORT,
    REDOX_THROUGHPUT_CEILING,
    ROS_REACTION_IDS,
    THIOREDOXIN_PEROXIDASE,
    TSA1_ORF,
    TSA2_ORF,
    YNO1_ORF,
    add_ros_module,
    clearance_cost,
    clearance_ladder,
    closed_medium_ranges,
    dose_peroxide,
    energy_from_peroxide,
    install_tsa1,
    peroxide_dose_response,
    superoxide_leak_response,
)

pytestmark = pytest.mark.integration

# glucose <= 10 with oxygen free, the regime every cost below was measured in.
GLUCOSE, BASE_GROWTH = 10.0, 0.8876853585594526
CYTOSOLIC_CATALASE_ORF, PEROXISOMAL_CATALASE_ORF = "YGR088W", "YDR256C"
CATALASES = (CYTOSOLIC_CATALASE_ORF, PEROXISOMAL_CATALASE_ORF)
GLUTATHIONE_PEROXIDASES = ("YNL229C", "YBR244W", "YCL035C", "YDR513W", "YIR037W", "YKL026C")
SOD1, SOD2 = "r_4270", "r_4190"
SUPEROXIDE_C, SUPEROXIDE_M = "s_3931", "s_3813"
# Numerically zero. The solver is pinned at 1e-7 and an LP vertex lands a few 1e-12 off it.
ZERO = 1e-8


@pytest.fixture(scope="module")
def _gem_template(yeast_gem_factory):
    model = yeast_gem_factory()
    configure(model)
    return model


@pytest.fixture
def gem(_gem_template, model_copy):
    return model_copy(_gem_template)


@pytest.fixture(scope="module")
def _ros_template(_gem_template, model_copy):
    model = add_ros_module(model_copy(_gem_template))
    configure(model)
    return model


@pytest.fixture
def ros(_ros_template, model_copy):
    return model_copy(_ros_template)


def _growth(model, dose=0.0, knockouts=(), glucose=GLUCOSE, oxygen=None):
    """Growth under a regime. ``dose=None`` leaves the exchange at whatever bounds it has."""
    with model as m:
        m.reactions.get_by_id("r_1714").lower_bound = -abs(glucose)
        if oxygen is not None:
            m.reactions.get_by_id("r_1992").lower_bound = -abs(oxygen)
        for orf in knockouts:
            m.genes.get_by_id(orf).knock_out()
        if dose is None:
            pass
        elif m.reactions.has_id(PEROXIDE_EXCHANGE):
            dose_peroxide(m, dose)
        elif dose:
            raise KeyError(f"{model.id!r} has no ROS module, so it cannot be dosed")
        value = m.slim_optimize()
        optimal = m.solver.status == "optimal"
    return float(value) if optimal else 0.0


def _fva(model, reaction_ids, glucose=GLUCOSE, open_ros=False):
    """FVA at a regime. ``open_ros`` puts the model in the state a caller's dose puts it in."""
    with model as m:
        m.reactions.get_by_id("r_1714").lower_bound = -abs(glucose)
        if open_ros and m.reactions.has_id(PEROXIDE_EXCHANGE):
            m.reactions.get_by_id(PEROXIDE_EXCHANGE).bounds = (-1000.0, 1000.0)
            for rid in (MITOCHONDRIAL_LEAK, NADPH_OXIDASE):
                m.reactions.get_by_id(rid).bounds = (0.0, 1000.0)
        return cobra.flux_analysis.flux_variability_analysis(
            m, reaction_list=list(reaction_ids), fraction_of_optimum=0.0,
            processes=FVA_PROCESSES)


class TestTheHoleThisModuleFills:
    def test_the_shipped_model_has_no_way_to_make_superoxide(self, gem):
        """The finding this module exists to fix: both SODs are dead reactions."""
        for mid in (SUPEROXIDE_C, SUPEROXIDE_M):
            met = gem.metabolites.get_by_id(mid)
            producers = [r.id for r in met.reactions if r.metabolites[met] > 0]

            assert producers == [], f"{mid} already has a source: {producers}"

    def test_the_shipped_model_has_no_extracellular_peroxide(self, gem):
        extracellular = [m.id for m in gem.metabolites
                         if m.formula == "H2O2" and m.compartment == "e"]

        assert extracellular == []

    def test_tsa1_and_yno1_are_absent_from_the_shipped_gene_list(self, gem):
        """1161 genes and neither the dominant cytosolic peroxiredoxin nor the NADPH oxidase."""
        assert not gem.genes.has_id(TSA1_ORF)
        assert not gem.genes.has_id(YNO1_ORF)
        assert gem.genes.has_id(TSA2_ORF)


class TestTheAdditionsThemselves:
    def test_installing_the_module_does_not_move_growth(self, gem, ros):
        """Nothing is forced, so the price list gains rows without any price changing."""
        with gem as m:
            m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE
            before = m.slim_optimize()

        assert before == pytest.approx(BASE_GROWTH, rel=1e-9)
        assert _growth(ros) == pytest.approx(before, rel=1e-9)

    def test_installing_the_module_does_not_move_the_shipped_default_either(self, gem, ros):
        """Glucose -1.0 as shipped. The two differ by 6.8e-8, below the 1e-7 tolerance."""
        assert ros.slim_optimize() == pytest.approx(gem.slim_optimize(), abs=1e-6)

    def test_every_added_reaction_is_mass_and_charge_balanced(self, ros):
        """check_mass_balance() reports a 'charge' key too, so an empty dict covers both."""
        for rid in ROS_REACTION_IDS:
            rxn = ros.reactions.get_by_id(rid)
            if rxn.boundary:
                continue

            assert rxn.check_mass_balance() == {}, f"{rid}: {rxn.check_mass_balance()}"

    def test_the_exchange_is_imbalanced_by_exactly_one_peroxide(self, ros):
        """A boundary reaction is meant to be unbalanced; this pins what it exchanges."""
        imbalance = ros.reactions.get_by_id(PEROXIDE_EXCHANGE).check_mass_balance()

        assert imbalance == {"H": -2, "O": -2}

    def test_every_added_reaction_can_carry_flux(self, ros):
        """Measured with the ROS routes open, because shut they pin the transport too."""
        fva = _fva(ros, ROS_REACTION_IDS, open_ros=True)

        for rid in ROS_REACTION_IDS:
            assert fva.loc[rid, "maximum"] - fva.loc[rid, "minimum"] > 1.0, f"{rid} is dead"

    def test_the_transport_is_ungated_because_diffusion_has_no_gene(self, ros):
        """Yeast9's own H2O2 transport r_1839 carries no GPR either; this follows it."""
        assert ros.reactions.get_by_id(PEROXIDE_TRANSPORT).gene_reaction_rule == ""
        assert ros.reactions.get_by_id("r_1839").gene_reaction_rule == ""

    def test_every_ros_route_is_shut_until_a_caller_asks_for_one(self, ros):
        """All three. An open route is a free sink or a free bypass, not an absent dose."""
        for rid in (PEROXIDE_EXCHANGE, MITOCHONDRIAL_LEAK, NADPH_OXIDASE):
            assert ros.reactions.get_by_id(rid).bounds == (0.0, 0.0), rid

    def test_every_addition_carries_a_citation_or_a_stated_assertion(self, ros):
        for key in list(ROS_REACTION_IDS) + [TSA1_ORF]:
            assert CITATIONS[key].strip(), f"{key} has no provenance"


class TestNoFreeLunch:
    def test_nothing_the_module_added_is_unbounded(self, ros):
        """The classic way a new exchange goes wrong. Measured with the routes opened."""
        fva = _fva(ros, [MITOCHONDRIAL_LEAK, NADPH_OXIDASE], open_ros=True)

        assert fva["maximum"].max() < 999.0
        assert fva["minimum"].min() > -999.0
        assert fva.loc[MITOCHONDRIAL_LEAK, "maximum"] == pytest.approx(60.0, rel=1e-6)
        assert fva.loc[NADPH_OXIDASE, "maximum"] == pytest.approx(114.0533, rel=1e-5)

    def test_an_open_leak_is_an_alternative_oxidase_the_zero_bound_is_what_stops(self, gem,
                                                                                ros):
        """The bug the (0, 0) default exists to close, pinned so it cannot come back.

        Leak, then SOD2, then cytochrome c peroxidase r_0437 reoxidises cytochrome c with no
        complex IV -- a terminal oxidase S. cerevisiae does not have. Complex III's GPR does
        not stop it: complex III is present and needed, and a GPR constrains no flux here.
        """
        cox1 = "Q0045"

        assert _growth(ros, dose=None, knockouts=(cox1,)) == pytest.approx(
            _growth(gem, dose=None, knockouts=(cox1,)), rel=1e-9)
        with ros as m:
            m.reactions.get_by_id(MITOCHONDRIAL_LEAK).bounds = (0.0, 1000.0)

            assert _growth(m, dose=None, knockouts=(cox1,)) == pytest.approx(0.371940,
                                                                             rel=1e-4)

    def test_installing_the_module_cannot_raise_growth_on_any_medium(self, gem, ros):
        """Every exchange open: the leak was worth 4.16% here before the bound was zeroed."""
        def wide_open(model):
            with model as m:
                for rxn in m.exchanges:
                    rxn.lower_bound = -1000.0
                return m.slim_optimize()

        assert wide_open(ros) == pytest.approx(wide_open(gem), rel=1e-9)

    def test_the_default_exchange_is_not_a_free_peroxide_sink(self, gem, ros):
        """The bug this default was changed to fix, pinned so it cannot come back.

        At (0, 1000) the exchange secretes without limit, which is a free disposal route for
        the peroxide the native oxidases make. The mutant lacking all nine cytosolic
        clearance genes is dead in the shipped model and went to full wild-type growth with
        the module installed and no dose asked for -- an existing answer changed by
        installing alone. Shut both ways it is dead in both models.
        """
        orfs = tuple(orf for _label, genes, _rid in CLEARANCE_ROUTES for orf in genes)
        shipped = tuple(o for o in orfs if gem.genes.has_id(o))

        assert _growth(gem, dose=None, knockouts=shipped) < 1e-9
        assert _growth(ros, dose=None, knockouts=orfs) < 1e-9
        with ros as m:
            m.reactions.get_by_id(PEROXIDE_EXCHANGE).bounds = (0.0, 1000.0)

            assert _growth(m, dose=None, knockouts=orfs) == pytest.approx(BASE_GROWTH,
                                                                          rel=1e-9)

    def test_no_addition_can_carry_flux_on_a_closed_medium(self, ros):
        """Every exchange shut both ways: a non-zero range would be an infeasible cycle."""
        ranges = closed_medium_ranges(
            ros, list(ROS_REACTION_IDS) + [SOD1, SOD2, "r_2111", "r_4046"]
            + list(CLEARANCE_REACTIONS))

        assert ranges.abs().max().max() < ZERO, ranges

    def test_peroxide_is_not_an_energy_source(self, ros):
        """Maximum ATP turnover on peroxide and nothing else."""
        assert energy_from_peroxide(ros, dose=10.0) == pytest.approx(0.0, abs=1e-6)

    def test_peroxide_is_not_a_growth_substrate(self, ros):
        with ros as m:
            for rxn in m.exchanges:
                rxn.lower_bound = 0.0
            m.reactions.get_by_id(PEROXIDE_EXCHANGE).bounds = (-100.0, 0.0)
            m.slim_optimize()

            assert m.solver.status != "optimal"


class TestSuperoxideNowHasASource:
    def test_both_superoxide_dismutases_carry_flux_once_there_is_a_leak(self, ros):
        """r_4190 (SOD2) and r_4270 (SOD1) were structurally dead; they are not now."""
        with ros as m:
            m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE
            m.reactions.get_by_id(MITOCHONDRIAL_LEAK).bounds = (1.0, 1.0)
            m.reactions.get_by_id(NADPH_OXIDASE).bounds = (1.0, 1.0)
            m.slim_optimize()

            assert m.reactions.get_by_id(SOD2).flux == pytest.approx(1.0, rel=1e-6)
            assert m.reactions.get_by_id(SOD1).flux == pytest.approx(1.0, rel=1e-6)

    def test_the_mitochondrial_leak_price_list(self, ros):
        """One ubiquinol diverted per two superoxide: 0.0117542 /h per mmol/gDCW/h."""
        frame = superoxide_leak_response(ros, [0.0, 1.0, 40.0, 60.0, 61.0])
        growth = dict(zip(frame["leak"], frame["growth"]))

        assert growth[0.0] - growth[1.0] == pytest.approx(0.0117542, rel=1e-4)
        assert 1.0 - growth[40.0] / BASE_GROWTH == pytest.approx(0.529652, rel=1e-4)
        assert growth[60.0] == pytest.approx(0.0, abs=1e-9)
        assert not frame["feasible"].iloc[-1]

    def test_the_nadph_oxidase_price_list(self, ros):
        """One NADPH per two superoxide: 0.0074244 /h per mmol/gDCW/h, the peroxide slope.

        Linear to 50 and no further -- past that something else binds and the leak gets
        dearer, which is why the module docstring gives a range rather than a slope alone.
        """
        frame = superoxide_leak_response(ros, [0.0, 1.0, 50.0, 100.0, 115.0],
                                         source=NADPH_OXIDASE)
        growth = dict(zip(frame["leak"], frame["growth"]))
        slope = growth[0.0] - growth[1.0]

        assert slope == pytest.approx(0.0074244, rel=1e-4)
        assert (growth[0.0] - growth[50.0]) / 50.0 == pytest.approx(slope, rel=1e-6)
        assert (growth[0.0] - growth[100.0]) / 100.0 > slope * 1.04
        assert 1.0 - growth[100.0] / BASE_GROWTH == pytest.approx(0.872165, rel=1e-4)
        assert not frame["feasible"].iloc[-1]


class TestADoseIsNowPayable:
    def test_a_dose_through_catalase_is_free(self, ros):
        """The honest, uncomfortable result: catalase needs no cofactor, so clearing is free.

        Peroxide toxicity is Fenton chemistry and damage repair, and stoichiometry has none
        of it. Anything the ODE layer charges for a dose comes from somewhere else.
        """
        frame = peroxide_dose_response(ros, [0.0, 5.0, 50.0], glucose_uptake=GLUCOSE)

        assert frame["growth_cost_fraction"].abs().max() < 1e-9
        assert frame.loc[frame["dose"] == 50.0, "r_0255"].iloc[0] == pytest.approx(25.0,
                                                                                   abs=1e-3)

    def test_a_dose_without_catalase_costs_growth_linearly(self, ros):
        """The number the ODE layer prices its ROS flux against: -7.4244e-3 /h per mmol."""
        frame = peroxide_dose_response(ros, [0.0, 1.0, 10.0, 50.0], glucose_uptake=GLUCOSE,
                                       knockouts=CATALASES)
        growth = dict(zip(frame["dose"], frame["growth"]))
        slope_1 = (growth[1.0] - growth[0.0]) / 1.0

        assert slope_1 == pytest.approx(-7.424447e-3, rel=1e-4)
        assert (growth[10.0] - growth[0.0]) / 10.0 == pytest.approx(slope_1, rel=1e-6)
        assert (growth[50.0] - growth[0.0]) / 50.0 == pytest.approx(slope_1, rel=1e-6)

    def test_the_thiol_route_is_what_pays_when_catalase_is_gone(self, ros):
        """Flux equals the dose, through r_0550, the reaction TSA1 was missing from."""
        frame = peroxide_dose_response(ros, [5.0], glucose_uptake=GLUCOSE,
                                       knockouts=CATALASES)

        assert frame["r_0255"].iloc[0] == pytest.approx(0.0, abs=1e-9)
        assert frame["r_0550"].iloc[0] == pytest.approx(5.0, abs=1e-3)

    def test_a_dose_above_the_redox_throughput_ceiling_is_infeasible(self, ros):
        """114.053 mmol/gDCW/h is all the reducing power glucose 10 can raise.

        The same number is the FVA ceiling of the exchange and of the NADPH oxidase,
        because all three spend the one currency.
        """
        frame = peroxide_dose_response(ros, [114.0, 114.1], glucose_uptake=GLUCOSE,
                                       knockouts=CATALASES)

        assert REDOX_THROUGHPUT_CEILING == pytest.approx(114.0533, abs=1e-3)
        assert list(frame["feasible"]) == [True, False]
        assert frame["growth"].iloc[0] == pytest.approx(0.000431, abs=1e-5)

    def test_deleting_catalase_goes_from_free_to_costly_once_there_is_a_dose(self, ros):
        """The whole point, in one number: 1.4e-6 undosed, 4.2e-2 at a dose of 5.

        The GEM audit found ablation of nearly every stress adaptation moves growth by less
        than 1e-5. That is still true undosed. It is not true of a dosed model.
        """
        cost = clearance_cost(ros, CATALASES, dose=5.0, glucose_uptake=GLUCOSE)

        assert cost.undosed_cost < 1e-5
        assert cost.dosed_cost == pytest.approx(0.041819, rel=1e-3)
        assert cost.dosed_cost / cost.undosed_cost > 1e4

    def test_forcing_catalase_used_to_price_manufacture_and_a_dose_halves_it(self, ros):
        """Forced at 2.0 catalase costs 3.345% undosed, 1.673% with half its H2O2 supplied.

        Undosed the model has to MAKE the peroxide through polyamine oxidase, which is what
        the old cost measured. Supplying half the substrate removes half of it.
        """
        def forced(dose):
            with ros as m:
                m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE
                dose_peroxide(m, dose)
                m.reactions.get_by_id("r_0255").bounds = (2.0, 2.0)
                return 1.0 - m.slim_optimize() / BASE_GROWTH

        assert forced(0.0) == pytest.approx(0.033454, rel=1e-3)
        assert forced(2.0) == pytest.approx(0.016726, rel=1e-3)


class TestWhichRoutePaysAndWhatItCosts:
    def test_the_ladder_prices_every_cytosolic_route_at_a_dose_of_five(self, ros):
        """Catalase free, either thiol route 4.182%, none of them infeasible."""
        ladder = clearance_ladder(ros, dose=5.0, glucose_uptake=GLUCOSE)

        assert list(ladder["routes_deleted"]) == [0, 1, 2, 3]
        assert ladder["growth_cost_fraction"].iloc[0] == pytest.approx(0.0, abs=1e-9)
        assert ladder["growth_cost_fraction"].iloc[1] == pytest.approx(0.041819, rel=1e-3)
        assert ladder["r_0255"].iloc[0] == pytest.approx(2.5, abs=1e-3)
        assert ladder["r_0550"].iloc[1] == pytest.approx(5.0, abs=1e-3)

    def test_the_glutathione_and_thioredoxin_routes_are_priced_identically(self, ros):
        """Both spend one NADPH per H2O2, so deleting either alone is exactly free.

        This is why the audit kept finding stress ablations cost nothing: with a dose
        imposed the ROUTE still costs 4.182%, it is the CHOICE of route that is free.
        """
        ladder = clearance_ladder(ros, dose=5.0, glucose_uptake=GLUCOSE)

        assert ladder["growth"].iloc[2] == pytest.approx(ladder["growth"].iloc[1], rel=1e-9)

    def test_with_every_cytosolic_route_gone_a_dose_is_infeasible_not_merely_costly(
            self, ros):
        """Clearance is work the model cannot decline -- and it cannot decline it undosed
        either, because the native oxidases make peroxide that also has to go somewhere."""
        ladder = clearance_ladder(ros, dose=5.0, glucose_uptake=GLUCOSE)

        # The solver is pinned at 1e-7, so 1e-9 was tighter than "zero" can mean;
        # it flaked here at 5.2e-9.
        assert ladder["undosed_growth"].iloc[3] == pytest.approx(
            0.0, abs=PINNED_TOLERANCE)
        assert ladder["growth"].isna().iloc[3]

    def test_the_peroxisomal_catalase_cannot_reach_a_cytosolic_dose(self, ros):
        """Yeast9 has no cytosol<->peroxisome peroxide transport, so CTA1 is not a route.

        r_0256 acts on s_0840 and nothing moves s_0837 there. Deleting CTA1 alone therefore
        changes nothing at a dose, which is why CLEARANCE_ROUTES names only three routes.
        """
        with ros as m:
            m.reactions.get_by_id("r_1714").lower_bound = -GLUCOSE
            dose_peroxide(m, 5.0)
            m.slim_optimize()

            assert m.reactions.get_by_id("r_0256").flux == pytest.approx(0.0, abs=1e-9)
        assert _growth(ros, dose=5.0, knockouts=(PEROXISOMAL_CATALASE_ORF,)) == pytest.approx(
            BASE_GROWTH, rel=1e-9)

    def test_the_nuclear_peroxiredoxin_route_is_structurally_blocked(self, ros):
        """r_1037 is named in CLEARANCE_REACTIONS and can never carry flux.

        Its nuclear thioredoxin has no reductase, so r_1839 and r_1037 are both dead. Named
        rather than dropped, because a route that exists and cannot run is worth knowing.
        """
        fva = _fva(ros, ["r_1037", "r_1839"])

        assert fva.abs().max().max() < ZERO

    def test_the_clearance_routes_named_are_the_ones_that_exist(self, ros):
        for rid, description in CLEARANCE_REACTIONS.items():
            assert ros.reactions.has_id(rid), f"{rid} ({description}) is not in the model"
        for _label, orfs, rid in CLEARANCE_ROUTES:
            assert ros.reactions.has_id(rid)
            for orf in orfs:
                assert ros.genes.has_id(orf), f"{orf} is not a gene in this model"


class TestTheOxygenTrap:
    def test_an_oxygen_limited_dose_raises_growth_instead_of_lowering_it(self, ros):
        """Pinned so it cannot change quietly: catalase makes O2 the bound withheld."""
        limited = peroxide_dose_response(ros, [0.0, 5.0], glucose_uptake=11.1,
                                         oxygen_uptake=3.7)
        free = peroxide_dose_response(ros, [0.0, 5.0], glucose_uptake=11.1,
                                      oxygen_uptake=None)

        assert limited["growth_cost_fraction"].iloc[1] == pytest.approx(-0.223145, rel=1e-3)
        assert free["growth_cost_fraction"].iloc[1] == pytest.approx(0.0, abs=1e-9)

    def test_the_trap_saturates_and_the_cell_ends_up_secreting_oxygen(self, ros):
        """At glucose 11.1 the dose lifts growth 0.3460 -> 0.9857 and then stops.

        The mechanism is visible in the oxygen exchange, which turns positive: past the
        saturation point catalase supplies more O2 than the cell can use.
        """
        frame = peroxide_dose_response(ros, [50.0, 100.0], glucose_uptake=11.1,
                                       oxygen_uptake=3.7)

        assert frame["growth"].iloc[0] == pytest.approx(0.985688, rel=1e-4)
        assert frame["growth"].iloc[1] == pytest.approx(frame["growth"].iloc[0], rel=1e-9)
        with ros as m:
            m.reactions.get_by_id("r_1714").lower_bound = -11.1
            m.reactions.get_by_id("r_1992").lower_bound = -3.7
            dose_peroxide(m, 100.0)
            m.slim_optimize()

            assert m.reactions.get_by_id("r_1992").flux == pytest.approx(25.585, rel=1e-3)


class TestTsa1:
    def test_tsa1_is_on_the_cytosolic_thioredoxin_peroxidase_beside_tsa2(self, ros):
        rule = ros.reactions.get_by_id(THIOREDOXIN_PEROXIDASE).gene_reaction_rule

        assert ros.genes.has_id(TSA1_ORF)
        assert ros.genes.get_by_id(TSA1_ORF).name == "TSA1"
        assert [r.id for r in ros.genes.get_by_id(TSA1_ORF).reactions] == [
            THIOREDOXIN_PEROXIDASE]
        assert rule.count(TSA1_ORF) == rule.count(TSA2_ORF) == 2

    def test_tsa1_changes_an_answer_the_shipped_gene_list_gets_wrong(self, gem, ros):
        """Without TSA1, deleting TSA2 on top of the catalases and the GPx pool is LETHAL.

        r_0550 is the last route a cytosolic dose has and the shipped GPR gates every clause
        of it on TSA2. With TSA1 present the same octuple mutant grows at 0.8877 /h undosed
        and 0.8506 /h at a dose of 5. Not a bound moving by a percent: 0 against 0.89.
        """
        knockouts = CATALASES + GLUTATHIONE_PEROXIDASES + (TSA2_ORF,)

        assert _growth(gem, dose=0.0, knockouts=knockouts) < 1e-9
        assert _growth(ros, dose=0.0, knockouts=knockouts) == pytest.approx(BASE_GROWTH,
                                                                           rel=1e-5)
        assert _growth(ros, dose=5.0, knockouts=knockouts) == pytest.approx(0.850562,
                                                                           rel=1e-4)


class TestRefusals:
    def test_install_tsa1_refuses_a_model_that_already_has_it(self, ros):
        with pytest.raises(ValueError, match=TSA1_ORF):
            install_tsa1(ros.copy())

    def test_add_ros_module_refuses_a_second_installation(self, ros):
        with pytest.raises(ValueError, match="already installed"):
            add_ros_module(ros)

    def test_add_ros_module_refuses_an_enzyme_constrained_model(self, gem):
        """A GECKO reaction without a protein draw is free, the one thing GECKO is for."""
        with gem as m:
            m.add_reactions([cobra.Reaction("prot_pool_exchange")])

            with pytest.raises(NotImplementedError, match="enzyme-constrained"):
                add_ros_module(m)

    def test_install_tsa1_refuses_a_gecko_split_reaction(self, ec_yeast_gem):
        """ecYeastGEM splits r_0550 into r_0550No1/No2, each drawing on the protein pool.

        What is published for Tsa1 is a second-order rate constant, not a kcat, so there is
        nothing honest to put in the draw. The refusal names both split reactions.
        """
        assert ec_yeast_gem.reactions.has_id("r_0550No1")
        assert not ec_yeast_gem.reactions.has_id("r_0550")

        with pytest.raises(KeyError, match="kcat"):
            install_tsa1(ec_yeast_gem.copy())

    def test_dose_peroxide_refuses_a_negative_dose(self, ros):
        with ros as m, pytest.raises(ValueError, match="non-negative"):
            dose_peroxide(m, -1.0)

    def test_dosing_a_model_without_the_module_names_the_function_that_installs_it(self,
                                                                                   gem):
        with pytest.raises(KeyError, match="add_ros_module"):
            peroxide_dose_response(gem, [1.0])
